// Bakes the solvent-excluded surface into the volume the Gaussian
// surface's march draws (gaussian_surface.wgsl), as a distance field:
// per voxel, D = distance to the nearest place a probe centre can be,
// which exceeds the probe radius exactly inside the excluded volume.
// Stored as a = 1 + D - probe (isovalue 1), rgb = colour of the atom the
// nearest boundary piece belongs to. `vv_core::ses::Ses::probe_distance`
// is the CPU reference: D is 0 outside the SAS; inside, the distance to
// the SAS boundary -- an atom's SAS sphere where no other sphere covers
// its radial point, or a free arc of a circle where two meet (from
// `vv_core::ses::build`). Past `clamp` it is clamped: only voxels near the
// surface need exact values.

struct Params {
    vol_origin: vec3<f32>,
    voxel: f32,
    dims: vec3<u32>,
    probe: f32,
    brick_dims: vec3<u32>,
    // Largest D stored exactly (probe radius plus two voxels).
    clamp: f32,
    grid_origin: vec3<f32>,
    // At least the largest SAS radius plus `clamp`, so a 3x3x3 cell
    // neighbourhood holds every sphere and arc within `clamp` of a voxel,
    // and every sphere that could cover a radial point that close.
    cell_size: f32,
    grid_dims: vec3<i32>,
    _pad: u32,
};

// xyz = position, w = van der Waals radius.
@group(0) @binding(0) var<storage, read> atoms: array<vec4<f32>>;
@group(0) @binding(1) var<storage, read> colors: array<u32>;
@group(0) @binding(2) var<uniform> params: Params;
@group(0) @binding(3) var<storage, read> cell_starts: array<u32>;
@group(0) @binding(4) var<storage, read> cell_atoms: array<u32>;
@group(0) @binding(5) var volume: texture_storage_3d<rgba16float, write>;
@group(0) @binding(6) var<storage, read_write> bricks: array<atomic<u32>>;
// CSR by first atom: atom i's tori are tori[torus_starts[i]..[i + 1]],
// each (other atom, first probe, last probe, -); probes FULL_CIRCLE for a
// free circle.
@group(0) @binding(7) var<storage, read> torus_starts: array<u32>;
@group(0) @binding(8) var<storage, read> tori: array<vec4<u32>>;
@group(0) @binding(9) var<storage, read> probes: array<vec4<f32>>;

// The z-slab this dispatch fills (one per submit, so no single submit
// runs long enough for the OS to reset the GPU).
struct Slab {
    z0: u32,
};
@group(1) @binding(0) var<uniform> slab: Slab;

const BRICK: u32 = 8u;
const FULL_CIRCLE: u32 = 0xffffffffu;
const TAU: f32 = 6.2831853;
// `vv_core::ses::INSIDE_EPS`.
const INSIDE_EPS: f32 = 1e-3;

fn unpack_color(c: u32) -> vec3<f32> {
    let s = vec3<f32>(f32(c & 0xffu), f32((c >> 8u) & 0xffu), f32((c >> 16u) & 0xffu)) / 255.0;
    return select(pow((s + 0.055) / 1.055, vec3<f32>(2.4)), s / 12.92, s <= vec3<f32>(0.04045));
}

fn flag_brick(b: vec3<u32>) {
    let i = b.x + b.y * params.brick_dims.x + b.z * params.brick_dims.x * params.brick_dims.y;
    atomicOr(&bricks[i], 1u);
}

fn cell_index(c: vec3<i32>) -> i32 {
    if (any(c < vec3<i32>(0)) || any(c >= params.grid_dims)) {
        return -1;
    }
    return c.x + c.y * params.grid_dims.x + c.z * params.grid_dims.x * params.grid_dims.y;
}

// `q` is in no SAS sphere but atom `skip`'s.
fn exposed(q: vec3<f32>, base: vec3<i32>, skip: u32) -> bool {
    for (var dz = -1; dz <= 1; dz++) {
        for (var dy = -1; dy <= 1; dy++) {
            for (var dx = -1; dx <= 1; dx++) {
                let ci = cell_index(base + vec3<i32>(dx, dy, dz));
                if (ci < 0) {
                    continue;
                }
                for (var k = cell_starts[ci]; k < cell_starts[ci + 1]; k++) {
                    let j = cell_atoms[k];
                    let a = atoms[j];
                    let big = a.w + params.probe;
                    let d = q - a.xyz;
                    if (j != skip && dot(d, d) < big * big - INSIDE_EPS) {
                        return false;
                    }
                }
            }
        }
    }
    return true;
}

// Angle of `v` about `n`, counter-clockwise from `start`, in [0, 2pi)
// (`vv_core::ses::angle_about`).
fn angle_about(n: vec3<f32>, start: vec3<f32>, v: vec3<f32>) -> f32 {
    let a = atan2(dot(n, cross(start, v)), dot(start, v));
    return select(a, a + TAU, a < 0.0);
}

// Distance from `p` to the free arc of the circle where atoms i's and
// j's SAS spheres meet, between probes a and b (`Ses::arc_distance`);
// a huge value when the spheres don't meet.
fn arc_distance(p: vec3<f32>, ci: vec3<f32>, ri: f32, cj: vec3<f32>, rj: f32, a: u32, b: u32) -> f32 {
    let v = cj - ci;
    let d = length(v);
    if (d <= 0.0 || d >= ri + rj || d + rj <= ri || d + ri <= rj) {
        return 3.4e38;
    }
    let n = v / d;
    let offset = (d * d + ri * ri - rj * rj) / (2.0 * d);
    let r2 = ri * ri - offset * offset;
    if (r2 <= 0.0) {
        return 3.4e38;
    }
    let center = ci + n * offset;
    let radius = sqrt(r2);
    let w = p - center;
    let z = dot(w, n);
    let radial = w - n * z;
    let s = length(radial);
    let on_circle = sqrt((s - radius) * (s - radius) + z * z);
    if (a == FULL_CIRCLE || s <= 1e-6) {
        return on_circle;
    }
    let pa = probes[a].xyz;
    let pb = probes[b].xyz;
    let start = pa - center;
    let span = select(angle_about(n, start, pb - center), TAU, a == b);
    if (angle_about(n, start, radial) <= span) {
        return on_circle;
    }
    return min(distance(p, pa), distance(p, pb));
}

@compute @workgroup_size(4, 4, 4)
fn bake(@builtin(global_invocation_id) local_gid: vec3<u32>) {
    let gid = local_gid + vec3<u32>(0u, 0u, slab.z0);
    if (any(gid >= params.dims)) {
        return;
    }
    let p = params.vol_origin + vec3<f32>(gid) * params.voxel;
    let base = vec3<i32>(floor((p - params.grid_origin) / params.cell_size));

    // Inside the SAS, how deep at least, and the nearest atom (the colour
    // of voxels no boundary piece reaches).
    var inside = false;
    var deep = -3.4e38;
    var nearest = 3.4e38;
    var owner = 0u;
    for (var dz = -1; dz <= 1; dz++) {
        for (var dy = -1; dy <= 1; dy++) {
            for (var dx = -1; dx <= 1; dx++) {
                let ci = cell_index(base + vec3<i32>(dx, dy, dz));
                if (ci < 0) {
                    continue;
                }
                for (var k = cell_starts[ci]; k < cell_starts[ci + 1]; k++) {
                    let i = cell_atoms[k];
                    let a = atoms[i];
                    let r = distance(p, a.xyz);
                    let big = a.w + params.probe;
                    inside = inside || r < big;
                    deep = max(deep, big - r);
                    if (r - a.w < nearest) {
                        nearest = r - a.w;
                        owner = i;
                    }
                }
            }
        }
    }

    var depth = 0.0;
    if (inside) {
        var best = params.clamp;
        if (deep < best) {
            for (var dz = -1; dz <= 1; dz++) {
                for (var dy = -1; dy <= 1; dy++) {
                    for (var dx = -1; dx <= 1; dx++) {
                        let ci = cell_index(base + vec3<i32>(dx, dy, dz));
                        if (ci < 0) {
                            continue;
                        }
                        for (var k = cell_starts[ci]; k < cell_starts[ci + 1]; k++) {
                            let i = cell_atoms[k];
                            let a = atoms[i];
                            let big = a.w + params.probe;
                            let r = distance(p, a.xyz);
                            // Its arcs lie on its sphere: none is nearer.
                            let on_sphere = abs(big - r);
                            if (on_sphere >= best) {
                                continue;
                            }
                            if (r > 1e-6) {
                                let q = a.xyz + (p - a.xyz) * (big / r);
                                if (exposed(q, base, i)) {
                                    best = on_sphere;
                                    owner = i;
                                }
                            }
                            for (var t = torus_starts[i]; t < torus_starts[i + 1u]; t++) {
                                let tor = tori[t];
                                let b = atoms[tor.x];
                                let arc = arc_distance(p, a.xyz, big, b.xyz, b.w + params.probe, tor.y, tor.z);
                                if (arc < best) {
                                    best = arc;
                                    owner = select(i, tor.x, distance(p, b.xyz) - b.w < r - a.w);
                                }
                            }
                        }
                    }
                }
            }
        }
        depth = best;
    }
    let value = 1.0 + depth - params.probe;
    textureStore(volume, gid, vec4<f32>(unpack_color(colors[owner]), value));

    // As in gaussian_density.wgsl: a voxel at 8b + 0 also belongs to brick
    // b - 1 along that axis.
    if (value >= 1.0) {
        let own = gid / BRICK;
        let edge = (gid % BRICK == vec3<u32>(0u)) & (gid > vec3<u32>(0u));
        let hi_b = min(own, params.brick_dims - 1u);
        let lo_b = min(select(own, own - 1u, edge), hi_b);
        for (var z = lo_b.z; z <= hi_b.z; z++) {
            for (var y = lo_b.y; y <= hi_b.y; y++) {
                for (var x = lo_b.x; x <= hi_b.x; x++) {
                    flag_brick(vec3<u32>(x, y, z));
                }
            }
        }
    }
}
