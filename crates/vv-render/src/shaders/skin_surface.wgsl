// GPU ray-caster for the skin surface: one screen-aligned billboard per
// patch, sized from the patch's bounding sphere (`vs_skin_surface`, the
// same construction as `draw.wgsl`'s `vs_sphere`), whose fragments
// ray-cast only that patch (skin_patch.wgsl), depth-tested against every
// other patch and every other representation by the shared depth buffer.

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
    // 0 = perspective, 1 = orthographic -- see `vv_render::camera::Projection`.
    projection: u32,
    // View-space clip plane: `dot(xyz, p) + w < 0` is cut away; (0, 0,
    // 0, 1) keeps everything. `vv_render::camera::CameraUniform::clip`.
    clip: vec4<f32>,
    // `vv_render::camera::CameraUniform::cut`, in `clip`'s form.
    cut: vec4<f32>,
};
@group(0) @binding(0) var<uniform> cam: Camera;

// xyz = position (world space), w = weight (unshrunk, `r^2 / shrink`).
@group(1) @binding(0) var<storage, read> atoms: array<vec4<f32>>;
@group(1) @binding(1) var<storage, read> colors: array<u32>;

// Layout must match `vv_render::scene::SkinSurfaceParams`.
struct SkinSurfaceParams {
    view_inv: mat4x4<f32>,
    shrink: f32,
    atom_count: u32,
    patch_count: u32,
    id_base: u32,
    material0: vec4<f32>,
    material1: vec4<f32>,
    material2: vec4<f32>,
};
@group(1) @binding(2) var<uniform> params: SkinSurfaceParams;

@group(1) @binding(3) var<storage, read> patches: array<Patch>;

// CSR over `patches`: patch i's competitors are
// `competitors[competitor_starts[i]..competitor_starts[i + 1]]`.
@group(1) @binding(4) var<storage, read> competitor_starts: array<u32>;
@group(1) @binding(5) var<storage, read> competitors: array<u32>;

fn patch_bound(p: Patch) -> vec3<f32> {
    return vec3<f32>(p.bx, p.by, p.bz);
}

fn unpack_color(c: u32) -> vec3<f32> {
    return vec3<f32>(
        f32(c & 0xffu),
        f32((c >> 8u) & 0xffu),
        f32((c >> 16u) & 0xffu),
    ) / 255.0;
}

@vertex
fn vs_skin_surface(@builtin(vertex_index) vertex: u32) -> VsOut {
    // Six consecutive vertices make one patch's quad.
    let idx = vertex / 6u;
    let p = patches[idx];
    let center = (cam.view * vec4<f32>(patch_bound(p), 1.0)).xyz;
    var out: VsOut;
    out.patch_index = idx;
    if (wholly_cut(center, p.br)) {
        // Past the cut: it would only discard, and being the nearest
        // thing on its own rays, nothing else can early-Z reject it --
        // skip rasterizing it at all.
        out.clip = vec4<f32>(2.0, 2.0, 2.0, 1.0);
        out.view_pos = center;
        return out;
    }
    let pos = billboard(center, p.br, vertex);
    out.clip = cam.proj * vec4<f32>(pos, 1.0);
    out.view_pos = pos;
    return out;
}

struct FsOut {
    @location(0) color: vec4<f32>,
    @location(1) normal: vec4<f32>,
    @location(2) id: u32,
    @builtin(frag_depth) depth: f32,
};


// Where a pixel's ray meets this patch's piece of the surface (or
// discards).
struct SkinHit {
    depth: f32,
    view_depth: f32,
    normal_view: vec3<f32>,
    dir_view: vec3<f32>,
    atom: u32,
};

// This patch's own surface, cut open where a plane crosses it: the solid
// cross-section is drawn separately, by every atom's own disc
// (`vs_atom_cap`/`fs_atom_cap` below) -- a patch only ever covers the
// solvent-exposed shell, so a buried atom has none; recasting just this
// one patch past the cut would leave buried interior pixels an
// inconsistent patchwork of hits and holes.
fn skin_hit(in: VsOut) -> SkinHit {
    // View-space ray exactly as `fs_sphere` builds it, then into world
    // space (the patches and the grid live there) via the inverse view.
    var ro_v: vec3<f32>;
    var rd_v: vec3<f32>;
    if (cam.projection == 1u) {
        ro_v = in.view_pos;
        rd_v = vec3<f32>(0.0, 0.0, -1.0);
    } else {
        ro_v = vec3<f32>(0.0, 0.0, 0.0);
        rd_v = normalize(in.view_pos);
    }
    let ro = (params.view_inv * vec4<f32>(ro_v, 1.0)).xyz;
    let rd = normalize((params.view_inv * vec4<f32>(rd_v, 0.0)).xyz);

    let piece = patches[in.patch_index];
    let s = skin_cast(piece, in.patch_index, ro, rd, 1e-4, params.shrink);
    if (s.t < 0.0) {
        discard;
    }
    let hit_view = (cam.view * vec4<f32>(ro + rd * s.t, 1.0)).xyz;
    if (clip_distance(hit_view) < 0.0) {
        discard;
    }
    let dir_view = normalize((cam.view * vec4<f32>(rd, 0.0)).xyz);
    let normal_view = normalize((cam.view * vec4<f32>(s.normal, 0.0)).xyz);
    let clip = cam.proj * vec4<f32>(hit_view, 1.0);
    return SkinHit(clip.z / clip.w, -hit_view.z, normal_view, dir_view, s.atom);
}

@fragment
fn fs_skin_surface(in: VsOut) -> FsOut {
    let h = skin_hit(in);
    let m = material(params.material0, params.material1, params.material2);
    var out: FsOut;
    out.depth = h.depth;
    out.color = vec4<f32>(shade(srgb_to_linear(unpack_color(colors[h.atom])), h.normal_view, h.dir_view, m), 1.0);
    out.normal = vec4<f32>(h.normal_view * 0.5 + 0.5, 1.0);
    out.id = params.id_base + h.atom + 1u;
    return out;
}

// Glass prepass: which patch (of which surface, by `id_base`) owns the
// nearest front-facing hit per pixel. The glass pass keeps only that
// patch's fragment, so the glass is one layer: no doubled seams where
// patches meet, no cavity walls buried inside the molecule. An id, not
// a depth-equality test, so it holds on every driver.
struct OwnerOut {
    @location(0) owner: vec2<u32>,
    @builtin(frag_depth) depth: f32,
};

@fragment
fn fs_skin_surface_owner(in: VsOut) -> OwnerOut {
    let h = skin_hit(in);
    if (dot(h.normal_view, h.dir_view) > 0.0) {
        discard;
    }
    return OwnerOut(vec2<u32>(params.id_base, in.patch_index + 1u), h.depth);
}

@group(2) @binding(0) var glass_owner: texture_2d<u32>;

@fragment
fn fs_skin_surface_glass(in: VsOut) -> GlassOut {
    let owner = textureLoad(glass_owner, vec2<i32>(in.clip.xy), 0).xy;
    if (any(owner != vec2<u32>(params.id_base, in.patch_index + 1u))) {
        discard;
    }
    let h = skin_hit(in);
    let m = material(params.material0, params.material1, params.material2);
    return glass(srgb_to_linear(unpack_color(colors[h.atom])), h.normal_view, h.dir_view, m, h.view_depth, h.depth);
}

// The solid interior a clip plane exposes, one atom at a time -- see
// `ses_surface.wgsl`'s `vs_atom_cap`/`fs_atom_cap`, which this mirrors.
// `atoms[i].w` is the skin weight `r^2 / shrink` (`skin_hit`'s doc), not a
// radius: `weight_for_radius` picks it so `sqrt(shrink * w) == r`, the
// atom's true van der Waals radius, at every shrink factor -- the same
// sphere an isolated atom's own vertex patch sits on.
struct AtomCapVs {
    @builtin(position) clip: vec4<f32>,
    @location(0) view_pos: vec3<f32>,
    @location(1) @interpolate(flat) center: vec3<f32>,
    @location(2) @interpolate(flat) radius: f32,
    @location(3) @interpolate(flat) atom: u32,
};

@vertex
fn vs_atom_cap(@builtin(vertex_index) vertex: u32) -> AtomCapVs {
    let idx = vertex / 6u;
    let a = atoms[idx];
    let r = sqrt(params.shrink * a.w);
    let center = (cam.view * vec4<f32>(a.xyz, 1.0)).xyz;
    var out: AtomCapVs;
    out.center = center;
    out.radius = r;
    out.atom = idx;
    if (no_cap_possible(center, r)) {
        out.clip = vec4<f32>(2.0, 2.0, 2.0, 1.0); // outside the frustum
        out.view_pos = center;
        return out;
    }
    out.view_pos = billboard(center, r, vertex);
    out.clip = cam.proj * vec4<f32>(out.view_pos, 1.0);
    return out;
}

@fragment
fn fs_atom_cap(in: AtomCapVs) -> FsOut {
    var ro: vec3<f32>;
    var rd: vec3<f32>;
    if (cam.projection == 1u) {
        ro = in.view_pos;
        rd = vec3<f32>(0.0, 0.0, -1.0);
    } else {
        ro = vec3<f32>(0.0, 0.0, 0.0);
        rd = normalize(in.view_pos);
    }
    let c = in.center;
    let r = in.radius;
    let oc = ro - c;
    let b = dot(rd, oc);
    let disc = b * b - dot(oc, oc) + r * r;
    if (disc < 0.0) {
        discard;
    }
    let t_near = -b - sqrt(disc);
    let kept = kept_span(ro, rd, t_near, -b + sqrt(disc));
    if (kept.start > kept.end || kept.start <= t_near) {
        discard;
    }
    var t = kept.start;
    let q = ro + rd * t - c;
    t -= cap_nudge(r * r - dot(q, q), t, dot(in.view_pos - ro, rd));
    if (t <= 0.0) {
        discard;
    }
    let p = ro + rd * t;
    let clip = cam.proj * vec4<f32>(p, 1.0);
    var m = material(params.material0, params.material1, params.material2);
    m.specular = 0.0;
    var out: FsOut;
    out.depth = clip.z / clip.w;
    out.color = vec4<f32>(shade(srgb_to_linear(unpack_color(colors[in.atom])), kept.normal, rd, m), 1.0);
    out.normal = vec4<f32>(kept.normal * 0.5 + 0.5, 1.0);
    out.id = params.id_base + in.atom + 1u;
    return out;
}
