// Shared declarations for the cull, sphere, and point shaders.

// `projection`: 0 = perspective, 1 = orthographic -- see
// `vv_render::camera::Projection`. Everywhere this shader (and cull.wgsl,
// which shares this declaration) casts a ray from a single eye point or
// divides an on-screen size by view depth, it must branch on this: both
// assume perspective and are simply wrong for parallel-ray orthographic.
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
    projection: u32,
    // View-space clip plane: `dot(xyz, p) + w < 0` is cut away; (0, 0,
    // 0, 1) keeps everything. `vv_render::camera::CameraUniform::clip`.
    clip: vec4<f32>,
    // `vv_render::camera::CameraUniform::cut`, in `clip`'s form.
    cut: vec4<f32>,
};

@group(0) @binding(0) var<uniform> cam: Camera;

// Per-page parameters; layout must match `vv_render::scene::PageParams`.
// Per page rather than per frame so several structures, each with its own
// representation, can share one frame and one picking buffer.
struct PageParams {
    radius_scale: f32,
    bond_radius: f32,
    atom_id_base: u32,
    bond_id_base: u32,
    // Added to every drawn radius (SAS: the probe; licorice: the stick).
    radius_offset: f32,
    // 1: every bond as a 1-px line (the lines style).
    lines: u32,
    min_quad_radius: f32,
    // World-space view-depth margin for glass self-occlusion culling
    // (cull.wgsl's `glass_hidden`); 0 disables it. `vv_render::scene::
    // glass_cull_margin`.
    glass_cull_margin: f32,
    material0: vec4<f32>,
    material1: vec4<f32>,
    material2: vec4<f32>,
};

// xyz = position, w = van der Waals radius (scaled by params.radius_scale on use).
@group(1) @binding(0) var<storage, read> atoms: array<vec4<f32>>;
@group(1) @binding(1) var<storage, read> colors: array<u32>;
@group(1) @binding(7) var<uniform> params: PageParams;

fn atom_radius(a: vec4<f32>) -> f32 {
    return a.w * params.radius_scale + params.radius_offset;
}

fn unpack_color(c: u32) -> vec4<f32> {
    return vec4<f32>(
        f32(c & 0xffu),
        f32((c >> 8u) & 0xffu),
        f32((c >> 16u) & 0xffu),
        f32((c >> 24u) & 0xffu),
    ) / 255.0;
}
