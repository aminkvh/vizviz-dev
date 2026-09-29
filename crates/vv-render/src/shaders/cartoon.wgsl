// Cartoon ribbons, arrows and coil tubes (`vv_core::cartoon`). Only the
// cross-sections are stored; each vertex of the tube between two
// sections is built here from `vertex_index`, so there are no vertex or
// index buffers. Ordinary rasterized triangles with hardware depth.

// Must match `atoms.wgsl`'s `Camera` and `draw.wgsl`'s `Style` layouts -
// same bind group (group 0), so this shares `cam`/`style` with every
// other draw pipeline without a second uniform upload.
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

// Per-mesh: this cartoon's slice of the frame-wide pick id space and its
// material (`vv_render::scene::CartoonParams`, rewritten every frame).
// Scalar padding, not `vec3<u32>` or `array<u32, 3>`: both change the
// layout in the uniform address space.
struct CartoonUniform {
    id_base: u32,
    _pad0: u32,
    _pad1: u32,
    _pad2: u32,
    material0: vec4<f32>,
    material1: vec4<f32>,
    material2: vec4<f32>,
};
@group(1) @binding(0) var<uniform> cartoon: CartoonUniform;

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
@group(1) @binding(1) var<storage, read> sections: array<Section>;
// Section pairs to join this frame (shaders/cartoon_lod.wgsl).
@group(1) @binding(2) var<storage, read> pairs: array<vec2<u32>>;

// Must match `vv_core::cartoon::RING`.
const RING: u32 = 8u;

fn unpack_color(c: u32) -> vec4<f32> {
    return vec4<f32>(
        f32(c & 0xffu),
        f32((c >> 8u) & 0xffu),
        f32((c >> 16u) & 0xffu),
        f32((c >> 24u) & 0xffu),
    ) / 255.0;
}

fn signed_pow(x: f32, e: f32) -> f32 {
    return sign(x) * pow(abs(x), e);
}

struct VsOut {
    @builtin(position) clip: vec4<f32>,
    @location(0) view_pos: vec3<f32>,
    @location(1) view_normal: vec3<f32>,
    @location(2) color: vec4<f32>,
    // Frame-wide pick id, index + base + 1 (0 stays free for "nothing
    // here" - the same convention `atoms.wgsl`'s draw shaders use).
    @location(3) @interpolate(flat) id: u32,
};

// cos and sin of 2 pi i / RING.
const RING_COS = array<f32, 8>(1.0, 0.70710678, 0.0, -0.70710678, -1.0, -0.70710678, 0.0, 0.70710678);
const RING_SIN = array<f32, 8>(0.0, 0.70710678, 1.0, 0.70710678, 0.0, -0.70710678, -1.0, -0.70710678);

// Twin of `vv_core::cartoon::CartoonMesh::vertex`. One instance per
// section pair, 2 * RING vertices: vertex v is ring vertex v % RING of
// the pair's first section (v < RING) or second (the shared index
// pattern `scene::cartoon_join_indices` makes the triangles, so each
// vertex is shaded once per pair).
@vertex
fn vs_cartoon(@builtin(vertex_index) v: u32, @builtin(instance_index) pair: u32) -> VsOut {
    let k = select(pairs[pair].x, pairs[pair].y, v >= RING);
    let s = sections[k];
    let c = RING_COS[v % RING];
    let sn = RING_SIN[v % RING];
    let x = s.half_width * signed_pow(c, s.roundness);
    let y = s.half_thickness * signed_pow(sn, s.roundness);
    let nx = signed_pow(c, 2.0 - s.roundness) / max(s.half_width, 1e-4);
    let ny = signed_pow(sn, 2.0 - s.roundness) / max(s.half_thickness, 1e-4);
    let position = s.center + s.across * x + s.up * y;
    let normal = normalize(s.across * nx + s.up * ny);

    let view_pos = (cam.view * vec4<f32>(position, 1.0)).xyz;
    var out: VsOut;
    out.clip = cam.proj * vec4<f32>(view_pos, 1.0);
    out.view_pos = view_pos;
    out.view_normal = (cam.view * vec4<f32>(normal, 0.0)).xyz;
    out.color = linear_rgba(unpack_color(s.color));
    out.id = cartoon.id_base + k + 1u;
    return out;
}

struct FsOut {
    @location(0) color: vec4<f32>,
    @location(1) normal: vec4<f32>,
    @location(2) id: u32,
};

// Back faces are drawn only as caps: the closed tube seen from inside,
// where a plane cut it open. `normal.a` 0 marks a cap for the pass that
// moves its depth onto the plane (cap_depth.wgsl).
@fragment
fn fs_cartoon(in: VsOut, @builtin(front_facing) front: bool) -> FsOut {
    if (clip_distance(in.view_pos) < 0.0) {
        discard;
    }
    let rd = normalize(in.view_pos);
    var n = normalize(in.view_normal);
    var m = material(cartoon.material0, cartoon.material1, cartoon.material2);
    if (!front) {
        let k = cap_on_ray(in.view_pos, cam.projection == 1u);
        if (k.start < 0.0) {
            discard;
        }
        n = k.normal;
        m.specular = 0.0;
    }
    var out: FsOut;
    out.color = vec4<f32>(shade(in.color.rgb, n, rd, m), in.color.a);
    out.normal = vec4<f32>(n * 0.5 + 0.5, select(0.0, 1.0, front));
    out.id = in.id;
    return out;
}

// Transparent: weighted-blended OIT (shading.wgsl's `glass_lit`),
// rasterized depth (no ray-cast gap to correct, unlike a sphere/cylinder
// impostor).
@fragment
fn fs_cartoon_glass(in: VsOut, @builtin(front_facing) front: bool) -> GlassColor {
    if (clip_distance(in.view_pos) < 0.0) {
        discard;
    }
    let rd = normalize(in.view_pos);
    var n = normalize(in.view_normal);
    var m = material(cartoon.material0, cartoon.material1, cartoon.material2);
    if (!front) {
        let k = cap_on_ray(in.view_pos, cam.projection == 1u);
        if (k.start < 0.0) {
            discard;
        }
        n = k.normal;
        m.specular = 0.0;
    }
    return glass_lit(in.color.rgb, n, rd, m, -in.view_pos.z);
}

@fragment
fn fs_cartoon_pick(in: VsOut, @builtin(front_facing) front: bool) -> @location(0) u32 {
    if (clip_distance(in.view_pos) < 0.0) {
        discard;
    }
    if (!front && cap_on_ray(in.view_pos, cam.projection == 1u).start < 0.0) {
        discard;
    }
    return in.id;
}
