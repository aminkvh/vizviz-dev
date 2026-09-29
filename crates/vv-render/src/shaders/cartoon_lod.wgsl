// Cartoon level of detail, once per frame before drawing: each residue
// step (span of cross-sections) outside the frustum is dropped; the rest
// emit just enough section pairs for their on-screen length (one per
// TARGET_PX, up to every join). Spans share end sections, so dropping
// sections inside one never opens a crack.

// Must match `atoms.wgsl`'s `Camera`.
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

// Layout must match `vv_render::scene::CartoonSectionGpu`.
struct Section {
    center: vec3<f32>,
    half_width: f32,
    across: vec3<f32>,
    half_thickness: f32,
    up: vec3<f32>,
    roundness: f32,
    color: u32,
    _pad0: u32,
    _pad1: u32,
    _pad2: u32,
};
@group(1) @binding(0) var<storage, read> sections: array<Section>;
@group(1) @binding(1) var<storage, read> spans: array<vec2<u32>>;
@group(1) @binding(2) var<storage, read_write> pairs: array<vec2<u32>>;
// DrawIndexedIndirectArgs: index count, instance count, ...
@group(1) @binding(3) var<storage, read_write> indirect: array<atomic<u32>, 5>;

const TARGET_PX: f32 = 3.0;
// Slack for the ribbon bending away from the chord between the span's
// ends (a residue step is ~3.8 A).
const BEND_MARGIN: f32 = 2.0;

fn in_frustum(c: vec3<f32>, r: f32) -> bool {
    for (var k = 0u; k < 6u; k++) {
        let pl = cam.planes[k];
        if (dot(pl.xyz, c) + pl.w < -r) {
            return false;
        }
    }
    return true;
}

@compute @workgroup_size(64)
fn lod(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    if (i >= arrayLength(&spans)) {
        return;
    }
    let span = spans[i];
    let a = sections[span.x];
    let b = sections[span.y];
    let length = distance(a.center, b.center);
    let reach = max(max(a.half_width, b.half_width), max(a.half_thickness, b.half_thickness));
    let center = (a.center + b.center) * 0.5;
    if (!in_frustum(center, length * 0.5 + reach + BEND_MARGIN)) {
        return;
    }
    var px = length * cam.proj_scale;
    if (cam.projection == 0u) {
        px = px / max(-(cam.view * vec4<f32>(center, 1.0)).z, 1e-3);
    }
    let joins = span.y - span.x;
    let count = clamp(u32(ceil(px / TARGET_PX)), 1u, joins);
    let base = atomicAdd(&indirect[1], count);
    for (var j = 0u; j < count; j++) {
        pairs[base + j] = vec2<u32>(span.x + j * joins / count, span.x + (j + 1u) * joins / count);
    }
}
