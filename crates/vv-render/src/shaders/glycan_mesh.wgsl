// 3D-SNFG glycan glyphs (`vv_core::glycan::mesh`): a plain triangle soup
// rebuilt on the CPU every frame (a residue's shape only needs its ring
// centroid and two reference points, not a spline) and uploaded as one
// vertex-buffer-free storage array, unlike `cartoon.wgsl`'s GPU-built
// ring sections. Ordinary rasterized triangles with hardware depth.

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
    clip: vec4<f32>,
    cut: vec4<f32>,
};
@group(0) @binding(0) var<uniform> cam: Camera;

// Per-mesh: this glycan's slice of the frame-wide pick id space and its
// material. Layout must match `vv_render::scene::CartoonParams`, reused
// as-is (a glycan has no properties of its own beyond size, already
// baked into the mesh's vertex positions on the CPU).
struct GlycanUniform {
    id_base: u32,
    _pad0: u32,
    _pad1: u32,
    _pad2: u32,
    material0: vec4<f32>,
    material1: vec4<f32>,
    material2: vec4<f32>,
};
@group(1) @binding(0) var<uniform> glycan: GlycanUniform;

// Layout must match `vv_render::scene::GlycanVertexGpu`.
struct Vertex {
    position: vec3<f32>,
    _pad0: f32,
    normal: vec3<f32>,
    color: u32,
};
@group(1) @binding(1) var<storage, read> verts: array<Vertex>;

fn unpack_color(c: u32) -> vec4<f32> {
    return vec4<f32>(
        f32(c & 0xffu),
        f32((c >> 8u) & 0xffu),
        f32((c >> 16u) & 0xffu),
        f32((c >> 24u) & 0xffu),
    ) / 255.0;
}

struct VsOut {
    @builtin(position) clip: vec4<f32>,
    @location(0) view_pos: vec3<f32>,
    @location(1) view_normal: vec3<f32>,
    @location(2) color: vec4<f32>,
    // Frame-wide pick id: `vertex_index` names one vertex of the flat-
    // shaded triangle it belongs to, so any of its 3 vertices resolves
    // to the same source atom (`GlycanGpu::source`) — the same "index +
    // base + 1" convention `atoms.wgsl`'s draw shaders use.
    @location(3) @interpolate(flat) id: u32,
};

@vertex
fn vs_glycan(@builtin(vertex_index) v: u32) -> VsOut {
    let vert = verts[v];
    let view_pos = (cam.view * vec4<f32>(vert.position, 1.0)).xyz;
    var out: VsOut;
    out.clip = cam.proj * vec4<f32>(view_pos, 1.0);
    out.view_pos = view_pos;
    out.view_normal = (cam.view * vec4<f32>(vert.normal, 0.0)).xyz;
    out.color = linear_rgba(unpack_color(vert.color));
    out.id = glycan.id_base + v + 1u;
    return out;
}

struct FsOut {
    @location(0) color: vec4<f32>,
    @location(1) normal: vec4<f32>,
    @location(2) id: u32,
};

@fragment
fn fs_glycan(in: VsOut, @builtin(front_facing) front: bool) -> FsOut {
    if (clip_distance(in.view_pos) < 0.0) {
        discard;
    }
    let rd = normalize(in.view_pos);
    // A clip plane cuts a glyph away whole rather than capping it open
    // (unlike `fs_cartoon`'s `cap_on_ray`): a solid few-Angstrom glyph
    // has no meaningful interior cross-section to show.
    var n = normalize(in.view_normal);
    if (!front) {
        n = -n;
    }
    let m = material(glycan.material0, glycan.material1, glycan.material2);
    var out: FsOut;
    out.color = vec4<f32>(shade(in.color.rgb, n, rd, m), in.color.a);
    out.normal = vec4<f32>(n * 0.5 + 0.5, 1.0);
    out.id = in.id;
    return out;
}

// Transparent: weighted-blended OIT (shading.wgsl's `glass_lit`),
// rasterized depth (no ray-cast gap, unlike a sphere/cylinder impostor).
@fragment
fn fs_glycan_glass(in: VsOut, @builtin(front_facing) front: bool) -> GlassColor {
    if (clip_distance(in.view_pos) < 0.0) {
        discard;
    }
    let rd = normalize(in.view_pos);
    var n = normalize(in.view_normal);
    if (!front) {
        n = -n;
    }
    let m = material(glycan.material0, glycan.material1, glycan.material2);
    return glass_lit(in.color.rgb, n, rd, m, -in.view_pos.z);
}

@fragment
fn fs_glycan_pick(in: VsOut) -> @location(0) u32 {
    if (clip_distance(in.view_pos) < 0.0) {
        discard;
    }
    return in.id;
}
