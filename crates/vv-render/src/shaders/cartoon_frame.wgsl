// A cartoon's cross-sections at a new frame: the WGSL twin of
// `vv_core::cartoon::section`. Each section's place (`Recipe`) is fixed;
// each residue slot's control point and guide (`spline`) come with the
// frame, so a playing trajectory uploads a few vectors per residue
// instead of every section.

struct Recipe {
    slot: u32,
    t: f32,
    first: u32,
    end: u32,
};

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

@group(0) @binding(0) var<storage, read> recipes: array<Recipe>;
// Slot s: control point at 2s, guide at 2s + 1.
@group(0) @binding(1) var<storage, read> spline: array<vec4<f32>>;
@group(0) @binding(2) var<storage, read_write> sections: array<Section>;

fn control(k: i32, first: i32, last: i32) -> vec3<f32> {
    return spline[2 * clamp(k, first, last)].xyz;
}

// `vv_core::backbone::catmull_rom`.
fn catmull_rom(p0: vec3<f32>, p1: vec3<f32>, p2: vec3<f32>, p3: vec3<f32>, t: f32) -> vec3<f32> {
    let t2 = t * t;
    let t3 = t2 * t;
    return 0.5 * ((2.0 * p1)
        + (-p0 + p2) * t
        + (2.0 * p0 - 5.0 * p1 + 4.0 * p2 - p3) * t2
        + (-p0 + 3.0 * p1 - 3.0 * p2 + p3) * t3);
}

// `v` without its component along unit `axis`, normalized; zero when
// parallel (`vv_core::cartoon::perpendicular`).
fn perpendicular(v: vec3<f32>, axis: vec3<f32>) -> vec3<f32> {
    let p = v - axis * dot(v, axis);
    let len = length(p);
    return select(vec3<f32>(0.0), p / len, len > 1e-12);
}

fn any_perpendicular(axis: vec3<f32>) -> vec3<f32> {
    let helper = select(vec3<f32>(0.0, 1.0, 0.0), vec3<f32>(1.0, 0.0, 0.0), abs(axis.x) < 0.9);
    return perpendicular(helper, axis);
}

fn normalize_or_zero(v: vec3<f32>) -> vec3<f32> {
    let len = length(v);
    return select(vec3<f32>(0.0), v / len, len > 1e-12);
}

@compute @workgroup_size(64)
fn cartoon_frame(@builtin(global_invocation_id) gid: vec3<u32>) {
    let k = gid.x;
    if (k >= arrayLength(&recipes)) {
        return;
    }
    let r = recipes[k];
    let i = i32(r.slot);
    let first = i32(r.first);
    let last = i32(r.end) - 1;
    let p0 = control(i - 1, first, last);
    let p1 = control(i, first, last);
    let p2 = control(i + 1, first, last);
    let p3 = control(i + 2, first, last);
    let t = r.t;
    let center = catmull_rom(p0, p1, p2, p3, t);
    let d = 1e-3;
    var tangent = normalize_or_zero(
        catmull_rom(p0, p1, p2, p3, min(t + d, 1.0)) - catmull_rom(p0, p1, p2, p3, max(t - d, 0.0)),
    );
    if (all(tangent == vec3<f32>(0.0))) {
        tangent = normalize_or_zero(p2 - p1);
    }
    let guide = mix(spline[2 * i + 1].xyz, spline[2 * min(i + 1, last) + 1].xyz, t);
    var across = perpendicular(guide, tangent);
    if (all(across == vec3<f32>(0.0))) {
        across = any_perpendicular(tangent);
    }
    sections[k].center = center;
    sections[k].across = across;
    sections[k].up = cross(tangent, across);
}
