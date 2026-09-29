// Screen-space ambient occlusion over the opaque pass's depth and normal
// targets, so every representation gets it (impostors, LOD points,
// cartoon, surfaces). Alchemy-style (McGuire et al., HPG 2011): sample
// points on a spiral around the pixel, occlusion from how far each one
// rises above the pixel's tangent plane, fading out past `ao_radius`.
// Runs at half resolution. The spiral is rotated per texel by a 2x2
// pattern that the composite pass (outline.wgsl) averages away with its
// 4x4 full-resolution (2x2 half-resolution) blur. The vertex shader
// comes from fullscreen.wgsl.
//
// The same pass writes screen-space shadows to the green channel: a
// short ray from the pixel toward the key light, marched through the
// depth buffer; the light is blocked where the ray passes just behind
// a nearer surface. Only what is on screen casts shadows.
//
// With an occlusion volume (`vv_render::OcclusionVolume`, the whole
// scene's occupancy in world space) both also get a long-range term:
// cones traced through the volume's mips for occlusion, and a march
// toward the key light for shadows cast by anything, on screen or not.
// They start a few voxels off the surface; the screen-space terms cover
// the contact range.

// Layout must match `vv_render::renderer::PostUniform` (and `Post` in
// outline.wgsl).
struct Post {
    inv_proj: mat4x4<f32>,
    outline_enabled: u32,
    outline_width: u32,
    near: f32,
    depth_threshold: f32,
    normal_threshold: f32,
    outline_strength: f32,
    ao_strength: f32,
    ao_radius: f32,
    fog_strength: f32,
    fog_near: f32,
    fog_far: f32,
    proj_y: f32,
    background: vec4<f32>,
    orthographic: u32,
    shadow_strength: f32,
    _pad0: u32,
    _pad1: u32,
    background_top: vec4<f32>,
    // View-space direction toward the key light (xyz).
    light_dir: vec4<f32>,
    view_inv: mat4x4<f32>,
    // Occlusion volume: corner (xyz) and voxel size (w); size in world
    // units (xyz) and 1 when there is a volume (w).
    vol_origin: vec4<f32>,
    vol_extent: vec4<f32>,
    // World-space clip plane; the volume is empty on its cut side.
    clip_world: vec4<f32>,
    // The camera's near cut, as `clip_world`.
    cut_world: vec4<f32>,
};
@group(1) @binding(0) var volume: texture_3d<f32>;
@group(1) @binding(1) var volume_sampler: sampler;

@group(0) @binding(0) var normal: texture_2d<f32>;
@group(0) @binding(1) var depth: texture_depth_2d;
@group(0) @binding(2) var<uniform> post: Post;

const SAMPLES: u32 = 12u;
const TURNS: f32 = 7.0;
// Radius limits on screen: small enough to stay cheap, large enough that
// a whole-structure view (atoms ~1 px) still gets large-scale shape.
const MIN_RADIUS_PX: f32 = 10.0;
const MAX_RADIUS_PX: f32 = 256.0;

fn view_pos(p: vec2<i32>, d: f32, dims: vec2<f32>) -> vec3<f32> {
    let uv = (vec2<f32>(p) + 0.5) / dims;
    let ndc = vec2<f32>(uv.x * 2.0 - 1.0, 1.0 - uv.y * 2.0);
    let v = post.inv_proj * vec4<f32>(ndc, d, 1.0);
    return v.xyz / v.w;
}

fn bayer2(p: vec2<i32>) -> f32 {
    let m = array<f32, 4>(0.0, 2.0, 3.0, 1.0);
    let q = vec2<u32>(p) % 2u;
    return (m[q.y * 2u + q.x] + 0.5) / 4.0;
}

// Pixel of a view-space point (both projections are centred).
fn to_pixel(v: vec3<f32>, dims: vec2<f32>) -> vec2<f32> {
    let w = select(-v.z, 1.0, post.orthographic == 1u);
    let proj_x = post.proj_y * dims.y / dims.x;
    let ndc = vec2<f32>(proj_x * v.x, post.proj_y * v.y) / w;
    return vec2<f32>(ndc.x * 0.5 + 0.5, 0.5 - ndc.y * 0.5) * dims;
}

const SHADOW_STEPS: u32 = 16u;
const SHADOW_LENGTH: f32 = 20.0;
const SHADOW_MIN_PX: f32 = 12.0;
const SHADOW_MAX_PX: f32 = 120.0;

// 1 = lit, 0 = the key light is blocked.
fn light_visibility(pos: vec3<f32>, n: vec3<f32>, texel: vec2<i32>, px_per_unit: f32, dims_i: vec2<i32>) -> f32 {
    let l = post.light_dir.xyz;
    if (dot(n, l) <= 0.0) {
        return 1.0;
    }
    let dims = vec2<f32>(dims_i);
    let length_px = clamp(SHADOW_LENGTH * px_per_unit, SHADOW_MIN_PX, SHADOW_MAX_PX);
    let len = length_px / px_per_unit;
    // An occluder at most this thick (about an atom): a surface far in
    // front of the ray is something else, not what blocks it.
    let thickness = 0.5 * len;
    let start = pos + n * (0.02 * len);
    let bias = 0.02 * len;
    let jitter = bayer2(texel);
    var blocked = 0.0;
    for (var i = 0u; i < SHADOW_STEPS; i++) {
        let along = (f32(i) + jitter) / f32(SHADOW_STEPS);
        let r = start + l * (along * len);
        let q = vec2<i32>(to_pixel(r, dims));
        if (any(q < vec2<i32>(0)) || any(q >= dims_i)) {
            break;
        }
        let dq = textureLoad(depth, q, 0);
        if (dq <= 0.0) {
            continue;
        }
        // How far the surface seen at q lies in front of the ray point;
        // soft edges at both ends, and a blocker far along the ray
        // shadows less (penumbra).
        let behind = view_pos(q, dq, dims).z - r.z;
        let inside = smoothstep(bias, 4.0 * bias, behind)
            * (1.0 - smoothstep(0.5 * thickness, thickness, behind));
        blocked = max(blocked, inside * (1.0 - along));
    }
    // Faces turned away from the light are already dark from shading.
    return 1.0 - blocked * smoothstep(0.0, 0.2, dot(n, l));
}

// Mean occupancy of a sphere of `radius` around world point `p`.
fn occupancy(p: vec3<f32>, radius: f32) -> f32 {
    if (min(dot(post.clip_world.xyz, p) + post.clip_world.w, dot(post.cut_world.xyz, p) + post.cut_world.w) < 0.0) {
        return 0.0;
    }
    let uvw = (p - post.vol_origin.xyz) / post.vol_extent.xyz;
    let lod = log2(max(radius / post.vol_origin.w, 1.0));
    return textureSampleLevel(volume, volume_sampler, uvw, lod).r;
}

const CONE_STEPS: u32 = 4u;

// Fraction of the sky above `n` the volume leaves open: three cones (the
// normal, and two opposite ones tilted 50 degrees and spun per pixel by
// the same 2x2 pattern the screen-space term uses, so the composite blur
// averages eight directions),
// each marched with doubling steps from three voxels out, with a
// footprint growing with distance. The normal cone counts double.
fn volume_ao(p: vec3<f32>, n: vec3<f32>, spin: f32) -> f32 {
    let voxel = post.vol_origin.w;
    let helper = select(vec3<f32>(1.0, 0.0, 0.0), vec3<f32>(0.0, 1.0, 0.0), abs(n.x) > 0.7);
    let t = normalize(cross(n, helper));
    let b = cross(n, t);
    var open = 0.0;
    for (var c = 0u; c < 3u; c++) {
        var d = n;
        if (c > 0u) {
            let angle = spin + f32(c) * 3.1415927;
            d = normalize(n * 0.64 + (t * cos(angle) + b * sin(angle)) * 0.77);
        }
        var vis = 1.0;
        var dist = 3.0 * voxel;
        for (var s = 0u; s < CONE_STEPS; s++) {
            let a = occupancy(p + d * dist, 0.5 * dist);
            vis *= 1.0 - clamp(a, 0.0, 1.0);
            dist *= 2.0;
        }
        open += select(1.0, 2.0, c == 0u) * vis;
    }
    return open / 4.0;
}

const VOLUME_SHADOW_STEPS: u32 = 10u;

// Light reaching `p` through the volume from the key light `l`.
fn volume_light(p: vec3<f32>, n: vec3<f32>, l: vec3<f32>) -> f32 {
    if (dot(n, l) <= 0.0) {
        return 1.0;
    }
    let voxel = post.vol_origin.w;
    var vis = 1.0;
    var dist = 2.0 * voxel;
    for (var s = 0u; s < VOLUME_SHADOW_STEPS; s++) {
        let a = occupancy(p + l * dist, max(0.1 * dist, 0.5 * voxel));
        vis *= 1.0 - clamp(3.0 * a, 0.0, 1.0);
        dist *= 1.5;
    }
    return vis;
}

@fragment
fn fs_ao(in: Vs) -> @location(0) vec4<f32> {
    let dims_i = vec2<i32>(textureDimensions(depth));
    let dims = vec2<f32>(dims_i);
    let texel = vec2<i32>(in.clip.xy);
    let p = min(texel * 2, dims_i - 1);
    let d = textureLoad(depth, p, 0);
    if (d <= 0.0) {
        return vec4<f32>(1.0);
    }
    let pos = view_pos(p, d, dims);
    let n = normalize(textureLoad(normal, p, 0).xyz * 2.0 - 1.0);

    // World units per pixel at this depth.
    var px_per_unit = 0.5 * dims.y * post.proj_y;
    if (post.orthographic == 0u) {
        px_per_unit = px_per_unit / max(-pos.z, 1e-3);
    }
    let has_volume = post.vol_extent.w > 0.5;
    let p_world = (post.view_inv * vec4<f32>(pos, 1.0)).xyz;
    let n_world = normalize((post.view_inv * vec4<f32>(n, 0.0)).xyz);
    var lit = 1.0;
    if (post.shadow_strength > 0.0) {
        lit = light_visibility(pos, n, texel, px_per_unit, dims_i);
        if (has_volume) {
            let l_world = normalize((post.view_inv * vec4<f32>(post.light_dir.xyz, 0.0)).xyz);
            lit *= volume_light(p_world, n_world, l_world);
        }
    }
    if (post.ao_strength <= 0.0) {
        return vec4<f32>(1.0, lit, 1.0, 1.0);
    }
    var far_ao = 1.0;
    if (has_volume) {
        far_ao = pow(volume_ao(p_world, n_world, bayer2(texel) * 6.2831853), post.ao_strength);
    }
    let radius_px = clamp(post.ao_radius * px_per_unit, MIN_RADIUS_PX, MAX_RADIUS_PX);
    let radius = radius_px / px_per_unit;

    let spin = bayer2(texel) * 6.2831853;
    var occlusion = 0.0;
    for (var i = 0u; i < SAMPLES; i++) {
        let alpha = (f32(i) + 0.5) / f32(SAMPLES);
        let angle = alpha * TURNS * 6.2831853 + spin;
        let offset = vec2<f32>(cos(angle), sin(angle)) * (alpha * radius_px);
        let q = p + vec2<i32>(round(offset));
        if (any(q < vec2<i32>(0)) || any(q >= dims_i)) {
            continue;
        }
        let dq = textureLoad(depth, q, 0);
        if (dq <= 0.0) {
            continue;
        }
        let v = view_pos(q, dq, dims) - pos;
        let dist = length(v);
        if (dist < 1e-4) {
            continue;
        }
        let rise = max(dot(n, v / dist) - 0.1, 0.0);
        let falloff = max(1.0 - (dist * dist) / (radius * radius), 0.0);
        occlusion += rise * falloff;
    }
    let ao = clamp(1.0 - post.ao_strength * 3.5 * occlusion / f32(SAMPLES), 0.0, 1.0);
    return vec4<f32>(ao * far_ao, lit, 1.0, 1.0);
}
