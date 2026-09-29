// Gaussian ("blobby") molecular surface: a full-screen ray march through
// the density volume baked by gaussian_density.wgsl (Blinn 1982;
// Krone, Stone, Ertl & Schulten, EuroVis 2012). One filtered
// texture fetch per step, and whole 8^3 bricks with nothing at the
// isovalue are skipped, so frame time follows the visible surface rather
// than atom count (8GLV, 4M atoms: see benchmarks/README.md).
//
// Marches in world space: one ray per pixel is moved into world space
// via `params.view_inv`, the hit goes back to view space once for depth
// and shading.

struct Camera {
    view_proj: mat4x4<f32>,
    view: mat4x4<f32>,
    proj: mat4x4<f32>,
    planes: array<vec4<f32>, 6>,
    viewport: vec2<f32>,
    quad_px_threshold: f32,
    proj_scale: f32,
    hiz_mip_count: u32,
    occlusion: u32,
    near: f32,
    // 0 = perspective, 1 = orthographic -- see `vv_render::camera::
    // Projection`. `fs_gaussian_surface` branches its ray generation on
    // this (perspective casts from a single eye point; orthographic
    // needs parallel rays -- see that function's comment).
    projection: u32,
    // View-space clip plane: `dot(xyz, p) + w < 0` is cut away; (0, 0,
    // 0, 1) keeps everything. `vv_render::camera::CameraUniform::clip`.
    clip: vec4<f32>,
    // `vv_render::camera::CameraUniform::cut`, in `clip`'s form.
    cut: vec4<f32>,
};
@group(0) @binding(0) var<uniform> cam: Camera;


// Same layout as `Params` in gaussian_density.wgsl and
// `vv_render::scene::GaussianSurfaceParams`.
struct GaussianSurfaceParams {
    view_inv: mat4x4<f32>,
    vol_origin: vec3<f32>,
    voxel: f32,
    dims: vec3<u32>,
    isovalue: f32,
    brick_dims: vec3<u32>,
    blob_factor: f32,
    material0: vec4<f32>,
    material1: vec4<f32>,
    material2: vec4<f32>,
};
@group(1) @binding(0) var<uniform> params: GaussianSurfaceParams;
// rgb = color, a = density; baked by gaussian_density.wgsl.
@group(1) @binding(1) var volume: texture_3d<f32>;
@group(1) @binding(2) var volume_sampler: sampler;
// Nonzero for every 8^3 brick that holds a voxel at or above the
// isovalue; every other brick is skipped whole.
@group(1) @binding(3) var<storage, read> bricks: array<u32>;

const BRICK: f32 = 8.0;
const MAX_BRICK_VISITS: u32 = 2048u;
const MAX_SAMPLES: u32 = 4096u;
const REFINE_STEPS: u32 = 8u;

// Voxel i sits at world `vol_origin + i * voxel`; texel centers are at
// (i + 0.5) / dims in texture coordinates.
fn sample_volume(p: vec3<f32>) -> vec4<f32> {
    let uvw = ((p - params.vol_origin) / params.voxel + 0.5) / vec3<f32>(params.dims);
    return textureSampleLevel(volume, volume_sampler, uvw, 0.0);
}

fn density(p: vec3<f32>) -> f32 {
    return sample_volume(p).a;
}

fn brick_occupied(b: vec3<i32>) -> bool {
    let d = vec3<i32>(params.brick_dims);
    return bricks[u32(b.x + b.y * d.x + b.z * d.x * d.y)] != 0u;
}

struct VsOut {
    @builtin(position) clip: vec4<f32>,
    @location(0) ndc: vec2<f32>,
};

// Full-screen triangle, no vertex buffer (as in blit.wgsl).
@vertex
fn vs_gaussian_surface(@builtin(vertex_index) vertex_index: u32) -> VsOut {
    let ndc = vec2<f32>(
        select(-1.0, 3.0, vertex_index == 1u),
        select(-1.0, 3.0, vertex_index == 2u),
    );
    var out: VsOut;
    out.clip = vec4<f32>(ndc, 0.0, 1.0);
    out.ndc = ndc;
    return out;
}

struct FsOut {
    @location(0) color: vec4<f32>,
    @location(1) normal: vec4<f32>,
    @location(2) id: u32,
    @builtin(frag_depth) depth: f32,
};

// Where a pixel's ray first meets the surface (or discards).
struct GaussianHit {
    depth: f32,
    view_depth: f32,
    normal_view: vec3<f32>,
    dir_view: vec3<f32>,
    base: vec3<f32>,
    // On the clip plane's cross-section (shaded matte).
    cap: bool,
};

fn gaussian_hit(in: VsOut) -> GaussianHit {
    // World-space ray. `ndc / proj[i][i]` inverts the projection for both
    // perspective (a direction) and orthographic (an offset); see
    // `vv_render::camera::Projection`.
    let offset_view = vec3<f32>(in.ndc.x / cam.proj[0][0], in.ndc.y / cam.proj[1][1], 0.0);
    var ro: vec3<f32>;
    var rd: vec3<f32>;
    if (cam.projection == 1u) {
        ro = (params.view_inv * vec4<f32>(offset_view, 1.0)).xyz;
        rd = normalize((params.view_inv * vec4<f32>(0.0, 0.0, -1.0, 0.0)).xyz);
    } else {
        ro = (params.view_inv * vec4<f32>(0.0, 0.0, 0.0, 1.0)).xyz;
        rd = normalize((params.view_inv * vec4<f32>(offset_view.x, offset_view.y, -1.0, 0.0)).xyz);
    }

    // Clip the ray to the volume's box.
    let box_min = params.vol_origin;
    let box_max = params.vol_origin + vec3<f32>(params.dims - 1u) * params.voxel;
    let inv = 1.0 / select(rd, vec3<f32>(1e-12), abs(rd) < vec3<f32>(1e-12));
    let ta = (box_min - ro) * inv;
    let tb = (box_max - ro) * inv;
    let t_near = max(max(min(ta.x, tb.x), min(ta.y, tb.y)), min(ta.z, tb.z));
    let t_far = min(min(max(ta.x, tb.x), max(ta.y, tb.y)), max(ta.z, tb.z));
    // Keep the part of the ray the cut planes keep, and cap where it
    // starts inside the surface.
    let kept = kept_span((cam.view * vec4<f32>(ro, 1.0)).xyz, (cam.view * vec4<f32>(rd, 0.0)).xyz, max(t_near, 0.0), t_far);
    if (kept.end <= kept.start) {
        discard;
    }
    let capped = kept.start > max(t_near, 0.0);
    let t0 = kept.start;
    let t_far_kept = kept.end;
    if (capped && density(ro + rd * t0) >= params.isovalue) {
        let hit_view = (cam.view * vec4<f32>(ro + rd * t0, 1.0)).xyz;
        let clip = cam.proj * vec4<f32>(hit_view, 1.0);
        return GaussianHit(
            clip.z / clip.w,
            -hit_view.z,
            kept.normal,
            normalize((cam.view * vec4<f32>(rd, 0.0)).xyz),
            sample_volume(ro + rd * t0).rgb,
            true,
        );
    }

    // Brick DDA (Amanatides & Woo): visit bricks front to back, march
    // only inside occupied ones.
    let brick_size = BRICK * params.voxel;
    let bdims = vec3<i32>(params.brick_dims);
    let start = ro + rd * t0;
    var b = clamp(vec3<i32>(floor((start - box_min) / brick_size)), vec3<i32>(0), bdims - 1);
    let dir_step = vec3<i32>(sign(rd));
    let next_edge = box_min + (vec3<f32>(b) + select(vec3<f32>(0.0), vec3<f32>(1.0), rd > vec3<f32>(0.0))) * brick_size;
    var t_max = select((next_edge - ro) * inv, vec3<f32>(3.4e38), dir_step == vec3<i32>(0));
    let t_delta = abs(brick_size * inv);

    let h = 0.5 * params.voxel;
    var t = t0;
    var prev_t = t0;
    var prev_d = density(start);
    var hit_t = -1.0;
    var samples = 0u;
    for (var visit = 0u; visit < MAX_BRICK_VISITS; visit++) {
        let t_exit = min(min(min(t_max.x, t_max.y), t_max.z), t_far_kept);
        if (brick_occupied(b)) {
            loop {
                t = min(t + h, t_exit);
                let d = density(ro + rd * t);
                samples++;
                if (prev_d < params.isovalue && d >= params.isovalue) {
                    var lo = prev_t;
                    var hi = t;
                    for (var r = 0u; r < REFINE_STEPS; r++) {
                        let mid = 0.5 * (lo + hi);
                        if (density(ro + rd * mid) < params.isovalue) {
                            lo = mid;
                        } else {
                            hi = mid;
                        }
                    }
                    hit_t = hi;
                    break;
                }
                prev_t = t;
                prev_d = d;
                if (t >= t_exit || samples >= MAX_SAMPLES) {
                    break;
                }
            }
            if (hit_t >= 0.0 || samples >= MAX_SAMPLES) {
                break;
            }
        } else {
            // Nothing in this brick reaches the isovalue.
            t = t_exit;
            prev_t = t_exit;
            prev_d = 0.0;
        }
        if (t_exit >= t_far_kept) {
            break;
        }
        if (t_max.x < t_max.y && t_max.x < t_max.z) {
            b.x += dir_step.x;
            t_max.x += t_delta.x;
        } else if (t_max.y < t_max.z) {
            b.y += dir_step.y;
            t_max.y += t_delta.y;
        } else {
            b.z += dir_step.z;
            t_max.z += t_delta.z;
        }
        if (any(b < vec3<i32>(0)) || any(b >= bdims)) {
            break;
        }
    }
    if (hit_t < 0.0) {
        discard;
    }

    let hit_world = ro + rd * hit_t;
    let e = params.voxel;
    let grad = vec3<f32>(
        density(hit_world + vec3<f32>(e, 0.0, 0.0)) - density(hit_world - vec3<f32>(e, 0.0, 0.0)),
        density(hit_world + vec3<f32>(0.0, e, 0.0)) - density(hit_world - vec3<f32>(0.0, e, 0.0)),
        density(hit_world + vec3<f32>(0.0, 0.0, e)) - density(hit_world - vec3<f32>(0.0, 0.0, e)),
    );
    let normal_world = normalize(-grad);

    let hit_view = (cam.view * vec4<f32>(hit_world, 1.0)).xyz;
    let dir_view = normalize((cam.view * vec4<f32>(rd, 0.0)).xyz);
    let normal_view = normalize((cam.view * vec4<f32>(normal_world, 0.0)).xyz);
    let clip = cam.proj * vec4<f32>(hit_view, 1.0);
    return GaussianHit(clip.z / clip.w, -hit_view.z, normal_view, dir_view, sample_volume(hit_world).rgb, false);
}

@fragment
fn fs_gaussian_surface(in: VsOut) -> FsOut {
    let h = gaussian_hit(in);
    var m = material(params.material0, params.material1, params.material2);
    if (h.cap) {
        m.specular = 0.0;
    }
    var out: FsOut;
    out.depth = h.depth;
    out.color = vec4<f32>(shade(h.base, h.normal_view, h.dir_view, m), 1.0);
    out.normal = vec4<f32>(h.normal_view * 0.5 + 0.5, 1.0);
    out.id = 0u; // Not pickable yet.
    return out;
}

// The march stops at the first crossing into the surface, which always
// faces the eye: one layer per pixel.
@fragment
fn fs_gaussian_surface_glass(in: VsOut) -> GlassOut {
    let h = gaussian_hit(in);
    var m = material(params.material0, params.material1, params.material2);
    if (h.cap) {
        m.specular = 0.0;
    }
    return glass(h.base, h.normal_view, h.dir_view, m, h.view_depth, h.depth);
}
