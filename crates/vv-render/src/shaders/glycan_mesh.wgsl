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

// Per-mesh: this glycan's slice of the frame-wide pick id space, its
// material, and how it is drawn (0 filled; 1 only the triangles' edges,
// `display_width` pixels wide -- the surfaces' mesh display). Layout must
// match `vv_render::scene::CartoonParams`.
struct GlycanUniform {
    id_base: u32,
    display: u32,
    display_width: f32,
    _pad: u32,
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
    // 1 at this vertex of its triangle, 0 at the other two: the distance
    // to an edge, in pixels, is a barycentric coordinate over its screen
    // gradient (`off_edges`).
    @location(4) bary: vec3<f32>,
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
    let corner = v % 3u;
    out.bary = vec3<f32>(f32(corner == 0u), f32(corner == 1u), f32(corner == 2u));
    return out;
}

// Screen-space length of each barycentric coordinate's gradient: pixels
// per unit. Taken before any `discard`, where derivatives are defined.
fn bary_gradient(bary: vec3<f32>) -> vec3<f32> {
    return vec3<f32>(
        length(vec2<f32>(dpdx(bary.x), dpdy(bary.x))),
        length(vec2<f32>(dpdx(bary.y), dpdy(bary.y))),
        length(vec2<f32>(dpdx(bary.z), dpdy(bary.z))),
    );
}

// In the mesh display, whether this pixel is further than half a line width
// from every edge of its triangle.
fn off_edges(bary: vec3<f32>, gradient: vec3<f32>) -> bool {
    if (glycan.display == 0u) {
        return false;
    }
    let edge_px = bary / max(gradient, vec3<f32>(1e-6));
    return min(edge_px.x, min(edge_px.y, edge_px.z)) > 0.5 * glycan.display_width;
}

struct FsOut {
    @location(0) color: vec4<f32>,
    @location(1) normal: vec4<f32>,
    @location(2) id: u32,
};

@fragment
fn fs_glycan(in: VsOut, @builtin(front_facing) front: bool) -> FsOut {
    let gradient = bary_gradient(in.bary);
    if (clip_distance(in.view_pos) < 0.0 || off_edges(in.bary, gradient)) {
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
    let gradient = bary_gradient(in.bary);
    if (clip_distance(in.view_pos) < 0.0 || off_edges(in.bary, gradient)) {
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
    let gradient = bary_gradient(in.bary);
    if (clip_distance(in.view_pos) < 0.0 || off_edges(in.bary, gradient)) {
        discard;
    }
    return in.id;
}
