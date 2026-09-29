// Path-traced spheres, bond cylinders, triangles and SES and skin
// patches (`vv_render::path_trace`; the patches are hit by
// ses_patch.wgsl and skin_patch.wgsl, appended after this file): per
// sample, a jittered camera ray finds the nearest surface through a BVH;
// the surface is shaded as the viewport shades it (shading.wgsl's
// `light_surface_seen`, appended after this file), each light counted
// only if a shadow ray toward a point of it (a narrow cone: soft
// shadows) gets out, and the ambient light only if an occlusion ray to a
// random direction of the sky above the surface does. Transparent
// materials let a camera ray through with probability 1 - opacity.

struct TraceParams {
    view_inv: mat4x4<f32>,
    view: mat4x4<f32>,
    background: vec4<f32>,
    background_top: vec4<f32>,
    // proj[0][0], proj[1][1], orthographic, transparent background.
    lens: vec4<f32>,
    size: vec2<u32>,
    tile_origin: vec2<u32>,
    first_sample: u32,
    sample_count: u32,
    unshadowed: u32,
    light_spread: f32,
    scene_radius: f32,
    probe_radius: f32,
    shrink: f32,
    // Scale the sky and key+fill terms (independently of `lighting`'s own
    // weights) before they reach `light_surface_seen`, so a render can
    // lean on occlusion without touching the viewport's `Lighting`
    // uniform. 1.0 both is parity with the viewport.
    ao_scale: f32,
    direct_scale: f32,
    // Three scalars, not `vec3<f32>`/`array<f32,3>`: a uniform-address
    // struct pads either to 16-byte alignment, which would desync this
    // hand-packed layout from `TraceParams` in `path_trace.rs`.
    _pad0: f32,
    _pad1: f32,
    _pad2: f32,
};

// What shading.wgsl reads of the camera: the view-space clip plane.
struct Camera {
    clip: vec4<f32>,
    // `vv_render::camera::CameraUniform::cut`, in `clip`'s form.
    cut: vec4<f32>,
};

struct Node {
    min: vec3<f32>,
    first: u32,
    max: vec3<f32>,
    count: u32,
};

struct Cylinder {
    // xyz = end a, w = radius.
    a: vec4<f32>,
    b: vec4<f32>,
};

struct MaterialGpu {
    m0: vec4<f32>,
    m1: vec4<f32>,
    m2: vec4<f32>,
};

@group(0) @binding(0) var<uniform> params: TraceParams;
@group(0) @binding(1) var<uniform> cam: Camera;
// @group(0) @binding(2): `lighting`, declared in shading.wgsl.
@group(0) @binding(3) var<storage, read> nodes: array<Node>;
// Top two bits: sphere 0, cylinder 1, triangle 2, patch 3; the rest
// indexes them.
@group(0) @binding(4) var<storage, read> refs: array<u32>;
// xyz = centre, w = radius.
@group(0) @binding(5) var<storage, read> spheres: array<vec4<f32>>;
// (colour, material).
@group(0) @binding(6) var<storage, read> sphere_looks: array<vec2<u32>>;
@group(0) @binding(7) var<storage, read> cylinders: array<Cylinder>;
// (colour at a, colour at b, material, -).
@group(0) @binding(8) var<storage, read> cylinder_looks: array<vec4<u32>>;
@group(0) @binding(9) var<storage, read> materials: array<MaterialGpu>;
// Per tile pixel: linear colour sum (premultiplied by coverage), coverage sum.
@group(0) @binding(10) var<storage, read_write> accum: array<vec4<f32>>;
// Per tile pixel: RGBA8, straight alpha, sRGB.
@group(0) @binding(11) var<storage, read_write> pixels: array<u32>;
// Mesh vertex v at 2v: (position, colour bits), (normal, material bits).
@group(0) @binding(12) var<storage, read> vertices: array<vec4<f32>>;
// Three vertex indices, and one unused.
@group(0) @binding(13) var<storage, read> triangles: array<vec4<u32>>;
// SES and skin patches (`vv_render::trace_surfaces`): each record, and
// what records point into, as ses_patch.wgsl and skin_patch.wgsl read
// them. Top two bits of a record: SES convex 0, torus 1, concave 2
// (`ses_patch.wgsl`), or skin 3 with its index in `skin_patches`.
@group(0) @binding(14) var<storage, read> patches: array<vec4<u32>>;
// xyz = position; w = van der Waals radius (SES) or weight (skin).
@group(0) @binding(15) var<storage, read> atoms: array<vec4<f32>>;
// (colour, material, -, -).
@group(0) @binding(16) var<storage, read> atom_looks: array<vec4<u32>>;
@group(0) @binding(17) var<storage, read> caps: array<u32>;
@group(0) @binding(18) var<storage, read> probes: array<Probe>;
@group(0) @binding(19) var<storage, read> probe_neighbors: array<u32>;
@group(0) @binding(20) var<storage, read> skin_patches: array<Patch>;
@group(0) @binding(21) var<storage, read> competitor_starts: array<u32>;
@group(0) @binding(22) var<storage, read> competitors: array<u32>;

const TILE: u32 = 256u;
const CYLINDER: u32 = 0x40000000u;
const TRIANGLE: u32 = 0x80000000u;
const PATCH: u32 = 0xc0000000u;
const KIND_CONCAVE: u32 = 2u;
const INDEX: u32 = 0x3fffffffu;
const NO_HIT: u32 = 0xffffffffu;
const FAR: f32 = 3.4e38;
const PI: f32 = 3.14159265;
const MAX_LAYERS: u32 = 8u;

// PCG hash (Jarzynski & Olano 2020).
fn pcg(v: u32) -> u32 {
    let state = v * 747796405u + 2891336453u;
    let word = ((state >> ((state >> 28u) + 4u)) ^ state) * 277803737u;
    return (word >> 22u) ^ word;
}

var<private> rng: u32;
// This pixel's hash and the sample index within it, for `stratified2`:
// set once per sample in `trace_tile`, read by every stratified draw
// `shade_sample` makes for that sample (jitter, shadow cones, AO).
var<private> pixel_hash: u32;
var<private> sample_index: u32;

fn random() -> f32 {
    rng = pcg(rng);
    return f32(rng >> 8u) / 16777216.0;
}

// Roberts (2018) R2 sequence: low-discrepancy, so successive samples
// fill [0,1)^2 more evenly than i.i.d. random pairs -- less noise at
// equal sample count. `salt` gives each dimension (pixel jitter, key and
// fill shadow cones, the AO hemisphere) its own Cranley-Patterson
// rotation, so they don't share a pattern with each other or across
// pixels.
const R2_A: vec2<f32> = vec2<f32>(0.7548776662, 0.5698402909);

fn stratified2(salt: u32) -> vec2<f32> {
    let h = pcg(pixel_hash + salt);
    let rotation = vec2<f32>(f32(h & 0xffffu), f32(h >> 16u)) / 65536.0;
    return fract(R2_A * f32(sample_index) + rotation);
}

fn unpack(c: u32) -> vec3<f32> {
    let s = vec3<f32>(f32(c & 0xffu), f32((c >> 8u) & 0xffu), f32((c >> 16u) & 0xffu)) / 255.0;
    return srgb_to_linear(s);
}

// Kept by the clip plane (world point).
fn kept(p: vec3<f32>) -> bool {
    return clip_distance((params.view * vec4<f32>(p, 1.0)).xyz) >= 0.0;
}

// A hit with an optional flat cap: `cap` marks a plane cross-section (the
// solid interior a cut exposes), `n` its outward normal in world space,
// valid only when `cap` is set -- `surface()` uses it in place of the
// primitive's own geometric normal, exactly as `draw.wgsl`'s `fs_sphere`
// caps the same primitive on screen.
struct RefHit {
    t: f32,
    cap: bool,
    n: vec3<f32>,
};

const NO_HIT_T: RefHit = RefHit(FAR, false, vec3<f32>(0.0));

// Nearest kept point of a sphere in (t_min, t_max), or the plane cap where
// the near side is cut away (view-space `ro_v`/`rd_v`, for `kept_span`).
fn hit_sphere(ro: vec3<f32>, rd: vec3<f32>, ro_v: vec3<f32>, rd_v: vec3<f32>, s: vec4<f32>, t_min: f32, t_max: f32) -> RefHit {
    let oc = ro - s.xyz;
    let b = dot(oc, rd);
    let disc = b * b - dot(oc, oc) + s.w * s.w;
    if (disc < 0.0) {
        return NO_HIT_T;
    }
    let root = sqrt(disc);
    let near = -b - root;
    let far = -b + root;
    // Unclamped by (t_min, t_max): a cap is a real plane crossing of the
    // sphere itself, not an artifact of where this search happens to start.
    let k = kept_span(ro_v, rd_v, near, far);
    if (k.start > k.end || k.start <= t_min || k.start >= t_max) {
        return NO_HIT_T;
    }
    let cap = k.start > near + 1e-4;
    var n = vec3<f32>(0.0);
    if (cap) {
        n = normalize((params.view_inv * vec4<f32>(k.normal, 0.0)).xyz);
    }
    return RefHit(k.start, cap, n);
}

// Nearest root of an open cylinder in (t_min, t_max) whose point is kept.
// Uncapped where cut (a known gap: bonds are cut open in `render`, capped
// in the viewport by `draw.wgsl`'s `fs_cylinder`).
fn hit_cylinder(ro: vec3<f32>, rd: vec3<f32>, c: Cylinder, t_min: f32, t_max: f32) -> f32 {
    let axis = c.b.xyz - c.a.xyz;
    let len2 = dot(axis, axis);
    if (len2 <= 0.0) {
        return FAR;
    }
    let oc = ro - c.a.xyz;
    let dd = dot(rd, axis);
    let od = dot(oc, axis);
    let qa = len2 - dd * dd;
    let qb = len2 * dot(oc, rd) - od * dd;
    let qc = len2 * dot(oc, oc) - od * od - c.a.w * c.a.w * len2;
    let h = qb * qb - qa * qc;
    if (qa <= 1e-12 || h < 0.0) {
        return FAR;
    }
    let root = sqrt(h);
    for (var k = 0; k < 2; k++) {
        let t = (select(-qb - root, -qb + root, k == 1)) / qa;
        let along = od + t * dd;
        if (t > t_min && t < t_max && along >= 0.0 && along <= len2 && kept(ro + rd * t)) {
            return t;
        }
    }
    return FAR;
}

// Distance and barycentrics (u, v of the second and third vertex) of a
// two-sided hit on triangle `k` (Moller & Trumbore 1997); t is FAR on a
// miss.
fn triangle_hit(ro: vec3<f32>, rd: vec3<f32>, k: u32) -> vec3<f32> {
    let t = triangles[k];
    let a = vertices[2u * t.x].xyz;
    let e1 = vertices[2u * t.y].xyz - a;
    let e2 = vertices[2u * t.z].xyz - a;
    let p = cross(rd, e2);
    let det = dot(e1, p);
    if (abs(det) < 1e-12) {
        return vec3<f32>(FAR, 0.0, 0.0);
    }
    let inv = 1.0 / det;
    let s = ro - a;
    let u = dot(s, p) * inv;
    let q = cross(s, e1);
    let v = dot(rd, q) * inv;
    if (u < 0.0 || v < 0.0 || u + v > 1.0) {
        return vec3<f32>(FAR, 0.0, 0.0);
    }
    return vec3<f32>(dot(e2, q) * inv, u, v);
}

// Cut open, uncapped: a mesh's own back faces would be needed for a solid
// cap (as `cartoon.wgsl` draws on screen), which the tracer's flat vertex
// buffer does not carry -- a known gap, see docs/COMMANDS.md.
fn hit_triangle(ro: vec3<f32>, rd: vec3<f32>, k: u32, t_min: f32, t_max: f32) -> RefHit {
    let t = triangle_hit(ro, rd, k).x;
    if (t > t_min && t < t_max && kept(ro + rd * t)) {
        return RefHit(t, false, vec3<f32>(0.0));
    }
    return NO_HIT_T;
}

fn probe_radius() -> f32 {
    return params.probe_radius;
}

fn atom_color(a: u32) -> vec3<f32> {
    return unpack(atom_looks[a].x);
}

// Where the ray first meets patch record `k` past `t_min`, cast exactly
// as the viewport casts it (`ok` false on a miss).
fn patch_surf(ro: vec3<f32>, rd: vec3<f32>, k: u32, t_min: f32) -> Surf {
    let p = patches[k];
    switch (kind_of(p)) {
        case KIND_CONVEX: {
            return hit_convex(p, ro, rd, t_min, 0.0, false);
        }
        case KIND_TORUS: {
            return hit_torus(p, ro, rd, t_min);
        }
        case KIND_CONCAVE: {
            return hit_concave(p, ro, rd, t_min);
        }
        default: {
            let i = first_of(p);
            let s = skin_cast(skin_patches[i], i, ro, rd, t_min, params.shrink);
            return Surf(s.t >= 0.0, s.t, s.normal, atom_color(s.atom), s.atom);
        }
    }
}

// Cut open, uncapped: a per-patch cap would need the plane point to fall
// on this same patch's own bounded piece of the surface, which leaves
// buried atoms (no patch of their own) unfilled -- see `ses_surface.wgsl`
// and `skin_surface.wgsl`'s own `vs_atom_cap`/`fs_atom_cap` for why the
// viewport uses a separate atom-disc layer instead, and docs/COMMANDS.md
// for this being a known viewport/tracer gap.
fn hit_patch(ro: vec3<f32>, rd: vec3<f32>, k: u32, t_min: f32, t_max: f32) -> RefHit {
    let s = patch_surf(ro, rd, k, t_min);
    if (s.ok && s.t < t_max && kept(ro + rd * s.t)) {
        return RefHit(s.t, false, vec3<f32>(0.0));
    }
    return NO_HIT_T;
}

fn hit_ref(ro: vec3<f32>, rd: vec3<f32>, ro_v: vec3<f32>, rd_v: vec3<f32>, r: u32, t_min: f32, t_max: f32) -> RefHit {
    switch (r & ~INDEX) {
        case PATCH: {
            return hit_patch(ro, rd, r & INDEX, t_min, t_max);
        }
        case CYLINDER: {
            let t = hit_cylinder(ro, rd, cylinders[r & INDEX], t_min, t_max);
            return RefHit(t, false, vec3<f32>(0.0));
        }
        case TRIANGLE: {
            return hit_triangle(ro, rd, r & INDEX, t_min, t_max);
        }
        default: {
            return hit_sphere(ro, rd, ro_v, rd_v, spheres[r], t_min, t_max);
        }
    }
}

// Entry distance of the ray into a node's box, or FAR.
fn enter_box(ro: vec3<f32>, inv: vec3<f32>, n: Node, t_max: f32) -> f32 {
    let t0 = (n.min - ro) * inv;
    let t1 = (n.max - ro) * inv;
    let near = max(max(min(t0.x, t1.x), min(t0.y, t1.y)), min(t0.z, t1.z));
    let far = min(min(max(t0.x, t1.x), max(t0.y, t1.y)), max(t0.z, t1.z));
    if (far < max(near, 0.0) || near > t_max) {
        return FAR;
    }
    return near;
}

struct Hit {
    t: f32,
    r: u32,
    cap: bool,
    // World-space flat cap normal; valid only when `cap` is set.
    cap_normal: vec3<f32>,
};

// The nearest hit in (t_min, t_max); with `first_found`, any hit.
fn trace(ro: vec3<f32>, rd: vec3<f32>, t_min: f32, t_max: f32, first_found: bool) -> Hit {
    let ro_v = (params.view * vec4<f32>(ro, 1.0)).xyz;
    let rd_v = (params.view * vec4<f32>(rd, 0.0)).xyz;
    let inv = 1.0 / select(rd, vec3<f32>(1e-20), abs(rd) < vec3<f32>(1e-20));
    var best = Hit(t_max, NO_HIT, false, vec3<f32>(0.0));
    var stack: array<u32, 48>;
    var top = 0u;
    if (enter_box(ro, inv, nodes[0], best.t) < FAR) {
        stack[0] = 0u;
        top = 1u;
    }
    while (top > 0u) {
        top -= 1u;
        let node = nodes[stack[top]];
        if (enter_box(ro, inv, node, best.t) >= FAR) {
            continue;
        }
        if (node.count > 0u) {
            for (var k = node.first; k < node.first + node.count; k++) {
                let r = refs[k];
                let h = hit_ref(ro, rd, ro_v, rd_v, r, t_min, best.t);
                if (h.t < best.t) {
                    best = Hit(h.t, r, h.cap, h.n);
                    if (first_found) {
                        return best;
                    }
                }
            }
            continue;
        }
        // Nearer child on top, so it is searched first.
        let a = enter_box(ro, inv, nodes[node.first], best.t);
        let b = enter_box(ro, inv, nodes[node.first + 1u], best.t);
        let near_first = a <= b;
        let near = select(node.first + 1u, node.first, near_first);
        let far = select(node.first, node.first + 1u, near_first);
        if (max(a, b) < FAR && top < 47u) {
            stack[top] = far;
            top += 1u;
        }
        if (min(a, b) < FAR && top < 48u) {
            stack[top] = near;
            top += 1u;
        }
    }
    return best;
}

struct Surface {
    // Shading normal (interpolated on a mesh) and the surface's own,
    // which says which side the ray hit.
    normal: vec3<f32>,
    geometric: vec3<f32>,
    base: vec3<f32>,
    m: Material,
};

fn surface(ro: vec3<f32>, rd: vec3<f32>, hit: Hit) -> Surface {
    let p = ro + rd * hit.t;
    if ((hit.r & ~INDEX) == TRIANGLE) {
        // Interpolated as the raster interpolates its vertices.
        let k = hit.r & INDEX;
        let b = triangle_hit(ro, rd, k);
        let w = vec3<f32>(1.0 - b.y - b.z, b.y, b.z);
        let t = triangles[k];
        let i = vec3<u32>(2u * t.x, 2u * t.y, 2u * t.z);
        let n = normalize(vertices[i.x + 1u].xyz * w.x + vertices[i.y + 1u].xyz * w.y + vertices[i.z + 1u].xyz * w.z);
        let base = unpack(bitcast<u32>(vertices[i.x].w)) * w.x
            + unpack(bitcast<u32>(vertices[i.y].w)) * w.y
            + unpack(bitcast<u32>(vertices[i.z].w)) * w.z;
        let mg = materials[bitcast<u32>(vertices[i.x + 1u].w)];
        let face = cross(vertices[i.y].xyz - vertices[i.x].xyz, vertices[i.z].xyz - vertices[i.x].xyz);
        // Wound counter-clockwise from outside, like the shading normals.
        return Surface(n, normalize(face), base, material(mg.m0, mg.m1, mg.m2));
    }
    if ((hit.r & ~INDEX) == PATCH) {
        // Found again from just before the hit, for its normal and colour.
        let s = patch_surf(ro, rd, hit.r & INDEX, hit.t - 1e-3 * max(1.0, hit.t * 1e-3));
        let mg = materials[atom_looks[s.atom].y];
        return Surface(s.normal, s.normal, s.color, material(mg.m0, mg.m1, mg.m2));
    }
    if ((hit.r & ~INDEX) == CYLINDER) {
        let k = hit.r & INDEX;
        let c = cylinders[k];
        let look = cylinder_looks[k];
        let axis = c.b.xyz - c.a.xyz;
        let along = dot(p - c.a.xyz, axis) / dot(axis, axis);
        let n = normalize(p - (c.a.xyz + axis * along));
        let mg = materials[look.z];
        return Surface(n, n, unpack(select(look.y, look.x, along < 0.5)), material(mg.m0, mg.m1, mg.m2));
    }
    let s = spheres[hit.r];
    let look = sphere_looks[hit.r];
    let mg = materials[look.y];
    var m = material(mg.m0, mg.m1, mg.m2);
    var n = normalize(p - s.xyz);
    if (hit.cap) {
        // A flat cut facing the light would mirror it: keep caps matte,
        // as `draw.wgsl`'s `fs_sphere` does for the same cross-section.
        n = hit.cap_normal;
        m.specular = 0.0;
    }
    return Surface(n, n, unpack(look.x), m);
}

// Two unit vectors perpendicular to unit `n`.
fn tangents(n: vec3<f32>) -> mat2x3<f32> {
    let helper = select(vec3<f32>(1.0, 0.0, 0.0), vec3<f32>(0.0, 1.0, 0.0), abs(n.x) > 0.9);
    let t = normalize(cross(n, helper));
    return mat2x3<f32>(t, cross(n, t));
}

// A direction within `spread` radians of unit `d`, uniform over the cap;
// `u` in [0,1)^2 (`stratified2`) picks the sample.
fn in_cone(d: vec3<f32>, spread: f32, u: vec2<f32>) -> vec3<f32> {
    let cos_t = 1.0 - u.x * (1.0 - cos(spread));
    let sin_t = sqrt(max(1.0 - cos_t * cos_t, 0.0));
    let phi = 2.0 * PI * u.y;
    let tb = tangents(d);
    return normalize(d * cos_t + (tb[0] * cos(phi) + tb[1] * sin(phi)) * sin_t);
}

// A cosine-weighted direction of the hemisphere about unit `n`; `u` in
// [0,1)^2 (`stratified2`) picks the sample.
fn cosine_hemisphere(n: vec3<f32>, u: vec2<f32>) -> vec3<f32> {
    let r = sqrt(u.x);
    let phi = 2.0 * PI * u.y;
    let tb = tangents(n);
    return normalize(n * sqrt(max(1.0 - r * r, 0.0)) + (tb[0] * cos(phi) + tb[1] * sin(phi)) * r);
}

// 1 when a ray from `p` toward `d` reaches the sky.
fn open_sky(p: vec3<f32>, d: vec3<f32>, eps: f32) -> f32 {
    if (params.unshadowed != 0u) {
        return 1.0;
    }
    return select(0.0, 1.0, trace(p, d, eps, 4.0 * params.scene_radius, true).r == NO_HIT);
}

fn background(py: f32) -> vec3<f32> {
    let f = py / f32(params.size.y);
    return mix(params.background_top.rgb, params.background.rgb, f);
}

// One sample of pixel `px`: linear colour premultiplied by coverage, and
// coverage.
fn shade_sample(px: vec2<f32>) -> vec4<f32> {
    let ndc = vec2<f32>(
        px.x / f32(params.size.x) * 2.0 - 1.0,
        1.0 - px.y / f32(params.size.y) * 2.0,
    );
    let offset = vec2<f32>(ndc.x / params.lens.x, ndc.y / params.lens.y);
    var ro: vec3<f32>;
    var rd: vec3<f32>;
    if (params.lens.z > 0.5) {
        ro = (params.view_inv * vec4<f32>(offset, 0.0, 1.0)).xyz;
        rd = normalize((params.view_inv * vec4<f32>(0.0, 0.0, -1.0, 0.0)).xyz);
    } else {
        ro = (params.view_inv * vec4<f32>(0.0, 0.0, 0.0, 1.0)).xyz;
        rd = normalize((params.view_inv * vec4<f32>(offset, -1.0, 0.0)).xyz);
    }
    let view_dir = normalize((params.view * vec4<f32>(rd, 0.0)).xyz);
    var t_min = 0.0;
    for (var layer = 0u; layer < MAX_LAYERS; layer++) {
        let hit = trace(ro, rd, t_min, FAR, false);
        if (hit.r == NO_HIT) {
            break;
        }
        let s = surface(ro, rd, hit);
        let n_view = normalize((params.view * vec4<f32>(s.normal, 0.0)).xyz);
        // Inside a primitive (the near side clipped away): seen from within.
        let inside = dot(s.geometric, rd) > 0.0;
        let n = select(s.normal, -s.normal, inside);
        let nv = select(n_view, -n_view, inside);
        if (random() >= material_alpha(nv, view_dir, s.m)) {
            t_min = hit.t * (1.0 + 1e-5) + 1e-4;
            continue;
        }
        let p = ro + rd * hit.t;
        let eps = 1e-3 * max(1.0, length(p) * 1e-3);
        let origin = p + select(s.geometric, -s.geometric, inside) * eps;
        // One shadow ray per live light, toward its world-space direction,
        // scaled by `direct_scale`; the sky/ambient ray by `ao_scale` --
        // see `TraceParams`.
        // Each ray draws from `stratified2` with its own salt (1..4 for up
        // to 4 lights, 5 for the sky) so none of them share a pattern.
        var seen_lights = array<f32, 4>(0.0, 0.0, 0.0, 0.0);
        let count = u32(lighting.count);
        for (var i = 0u; i < count; i++) {
            let dir = normalize((params.view_inv * vec4<f32>(lighting.lights[i].dir, 0.0)).xyz);
            if (dot(n, dir) > 0.0) {
                seen_lights[i] = params.direct_scale * open_sky(origin, in_cone(dir, params.light_spread, stratified2(1u + i)), eps);
            }
        }
        let seen_sky = params.ao_scale * open_sky(origin, cosine_hemisphere(n, stratified2(5u)), eps);
        let lit = light_surface_seen(s.base, nv, view_dir, s.m, seen_lights, seen_sky);
        return vec4<f32>(tonemap(lit.diffuse + vec3<f32>(lit.specular)), 1.0);
    }
    if (params.lens.w > 0.5) {
        return vec4<f32>(0.0);
    }
    return vec4<f32>(background(px.y), 1.0);
}

@compute @workgroup_size(8, 8, 1)
fn trace_tile(@builtin(global_invocation_id) gid: vec3<u32>) {
    let pixel = params.tile_origin + gid.xy;
    let slot = gid.y * TILE + gid.x;
    if (any(pixel >= params.size)) {
        return;
    }
    var sum = select(accum[slot], vec4<f32>(0.0), params.first_sample == 0u);
    pixel_hash = pixel.x * 2654435761u + pixel.y * 2246822519u;
    for (var s = params.first_sample; s < params.first_sample + params.sample_count; s++) {
        rng = pcg(pixel.x + pcg(pixel.y + pcg(s)));
        sample_index = s;
        // The first sample through the pixel centre, as the raster samples.
        let jitter = select(stratified2(0u), vec2<f32>(0.5), s == 0u);
        sum += shade_sample(vec2<f32>(pixel) + jitter);
    }
    accum[slot] = sum;
}

fn linear_to_srgb(c: vec3<f32>) -> vec3<f32> {
    let x = clamp(c, vec3<f32>(0.0), vec3<f32>(1.0));
    return select(1.055 * pow(x, vec3<f32>(1.0 / 2.4)) - 0.055, x * 12.92, x <= vec3<f32>(0.0031308));
}

@compute @workgroup_size(8, 8, 1)
fn resolve(@builtin(global_invocation_id) gid: vec3<u32>) {
    let pixel = params.tile_origin + gid.xy;
    let slot = gid.y * TILE + gid.x;
    if (any(pixel >= params.size)) {
        return;
    }
    let mean = accum[slot] / f32(params.first_sample + params.sample_count);
    // Straight alpha: the colour of what covers the pixel.
    let rgb = linear_to_srgb(select(mean.rgb, mean.rgb / mean.a, mean.a > 0.0));
    let byte = vec4<u32>(round(vec4<f32>(rgb, clamp(mean.a, 0.0, 1.0)) * 255.0));
    pixels[slot] = byte.x | (byte.y << 8u) | (byte.z << 16u) | (byte.w << 24u);
}
