// The skin surface's patches (Edelsbrunner 1999), hit by a ray: the WGSL
// twin of `vv_cpu::skin_surface`, same patch quadric and the same exact
// mixed-cell membership test (`vv_core::skin_surface` -- read its module
// doc first; every geometric decision lives there, verified against the
// primary source, not here). The membership test's Voronoi half checks
// the patch's own competitor list (`SkinComplex::owned_by_competitors`
// explains why that short list is exact). Shared by the viewport
// (`skin_surface.wgsl`) and the path tracer (`path_trace.wgsl`), which
// each declare `atoms` (xyz = position, w = weight), `competitor_starts`
// and `competitors`.

// Layout must match `vv_render::scene::SkinPatchGpu` field-for-field
// (sixteen 4-byte scalars, natural stride 64).
struct Patch {
    cx: f32,
    cy: f32,
    cz: f32,
    s_weight: f32,
    ax: f32,
    ay: f32,
    az: f32,
    kind: u32,
    bx: f32,
    by: f32,
    bz: f32,
    br: f32,
    a0: u32,
    a1: u32,
    a2: u32,
    a3: u32,
};

const KIND_VERTEX: u32 = 0u;
const KIND_EDGE: u32 = 1u;
const KIND_TRIANGLE: u32 = 2u;
const KIND_TET: u32 = 3u;
const NO_ATOM: u32 = 0xffffffffu;
// `vv_core::skin_surface::MEMBERSHIP_EPS`.
const MEMBERSHIP_EPS: f32 = 2e-2;

fn patch_center(p: Patch) -> vec3<f32> {
    return vec3<f32>(p.cx, p.cy, p.cz);
}

fn patch_axis(p: Patch) -> vec3<f32> {
    return vec3<f32>(p.ax, p.ay, p.az);
}

fn power_distance(p: vec3<f32>, atom_index: u32) -> f32 {
    let a = atoms[atom_index];
    let d = p - a.xyz;
    return dot(d, d) - a.w;
}

// The `(alpha, beta)` of `F = alpha |u|^2 + beta (u . axis)^2 + s w_z`
// -- exact port of `vv_core::skin_surface::coefficients`.
fn coefficients(kind: u32, shrink: f32) -> vec2<f32> {
    let k = shrink / (1.0 - shrink);
    switch (kind) {
        case KIND_VERTEX: { return vec2<f32>(1.0, 0.0); }
        case KIND_EDGE: { return vec2<f32>(1.0, -(1.0 + k)); }
        case KIND_TRIANGLE: { return vec2<f32>(-k, 1.0 + k); }
        default: { return vec2<f32>(-k, 0.0); }
    }
}

// Exact port of `vv_core::skin_surface::SkinComplex::owned_by_competitors`:
// no competitor of the simplex is nearer to `q` in power distance than
// the simplex's own atoms (which all tie there).
fn voronoi_owned(p: Patch, patch_index: u32, q: vec3<f32>) -> bool {
    let bound = power_distance(q, p.a0) - MEMBERSHIP_EPS;
    let lo = competitor_starts[patch_index];
    let hi = competitor_starts[patch_index + 1u];
    for (var m = lo; m < hi; m = m + 1u) {
        if (power_distance(q, competitors[m]) < bound) {
            return false;
        }
    }
    return true;
}

// Exact port of `Patch::delaunay_contains`: is `pt` (already in the
// simplex's affine hull) inside the simplex, by barycentric coordinates.
fn delaunay_contains(p: Patch, pt: vec3<f32>) -> bool {
    let eps = MEMBERSHIP_EPS;
    let a = atoms[p.a0].xyz;
    switch (p.kind) {
        case KIND_VERTEX: {
            return true;
        }
        case KIND_EDGE: {
            let e = atoms[p.a1].xyz - a;
            let t = dot(pt - a, e) / dot(e, e);
            return t >= -eps && t <= 1.0 + eps;
        }
        case KIND_TRIANGLE: {
            let v0 = atoms[p.a1].xyz - a;
            let v1 = atoms[p.a2].xyz - a;
            let v2 = pt - a;
            let d00 = dot(v0, v0);
            let d01 = dot(v0, v1);
            let d11 = dot(v1, v1);
            let d20 = dot(v2, v0);
            let d21 = dot(v2, v1);
            let denom = d00 * d11 - d01 * d01;
            let v = (d11 * d20 - d01 * d21) / denom;
            let w = (d00 * d21 - d01 * d20) / denom;
            let u = 1.0 - v - w;
            return u >= -eps && v >= -eps && w >= -eps;
        }
        default: {
            let e1 = atoms[p.a1].xyz - a;
            let e2 = atoms[p.a2].xyz - a;
            let e3 = atoms[p.a3].xyz - a;
            let r = pt - a;
            let det = dot(e1, cross(e2, e3));
            let l1 = dot(r, cross(e2, e3)) / det;
            let l2 = dot(e1, cross(r, e3)) / det;
            let l3 = dot(e1, cross(e2, r)) / det;
            let l0 = 1.0 - l1 - l2 - l3;
            return l0 >= -eps && l1 >= -eps && l2 >= -eps && l3 >= -eps;
        }
    }
}

// Exact port of `Patch::contains`: decompose `x - z` into its component
// along the simplex (scaled by 1/(1-s), must land in the simplex) and
// across it (scaled by 1/s, must land in the simplex's power cell).
fn contains(p: Patch, patch_index: u32, x: vec3<f32>, shrink: f32) -> bool {
    let z = patch_center(p);
    let axis = patch_axis(p);
    let u = x - z;
    let d = dot(u, axis);
    var u_par: vec3<f32>;
    var u_perp: vec3<f32>;
    switch (p.kind) {
        case KIND_VERTEX: { u_par = vec3<f32>(0.0); u_perp = u; }
        case KIND_EDGE: { u_par = axis * d; u_perp = u - axis * d; }
        case KIND_TRIANGLE: { u_par = u - axis * d; u_perp = axis * d; }
        default: { u_par = u; u_perp = vec3<f32>(0.0); }
    }
    let pt = z + u_par / (1.0 - shrink);
    if (!delaunay_contains(p, pt)) {
        return false;
    }
    let q = z + u_perp / shrink;
    return voronoi_owned(p, patch_index, q);
}

// Exact port of `Patch::nearest_atom`.
fn nearest_member(p: Patch, x: vec3<f32>) -> u32 {
    var best = p.a0;
    var best_d = power_distance(x, p.a0);
    if (p.a1 != NO_ATOM) {
        let d = power_distance(x, p.a1);
        if (d < best_d) { best = p.a1; best_d = d; }
    }
    if (p.a2 != NO_ATOM) {
        let d = power_distance(x, p.a2);
        if (d < best_d) { best = p.a2; best_d = d; }
    }
    if (p.a3 != NO_ATOM) {
        let d = power_distance(x, p.a3);
        if (d < best_d) { best = p.a3; best_d = d; }
    }
    return best;
}

// Ascending roots of `a t^2 + b t + c = 0` that lie past `t_min`, as
// (count, t1, t2) -- exact port of `vv_core::skin_surface::
// quadratic_roots` (the cancellation-free q-form; see its doc for why
// the textbook formula drops hits near a hyperboloid's asymptotic cone)
// plus `Patch::ray_roots`' ordering/filtering. `shift` is added to each
// root first (the caller solves from the ray's closest approach).
fn positive_roots(a: f32, b: f32, c: f32, shift: f32, t_min: f32) -> vec3<f32> {
    var r1 = 0.0;
    var r2 = 0.0;
    var n = 0u;
    if (abs(a) < 1e-12) {
        if (abs(b) < 1e-12) {
            return vec3<f32>(0.0, 0.0, 0.0);
        }
        r1 = -c / b;
        n = 1u;
    } else {
        let disc = b * b - 4.0 * a * c;
        if (disc < 0.0) {
            return vec3<f32>(0.0, 0.0, 0.0);
        }
        // Not sign(b): WGSL's sign(0) is 0 (Rust's signum is 1), and b is
        // exactly 0 for axis-aligned orthographic rays.
        let q = -0.5 * (b + select(-1.0, 1.0, b >= 0.0) * sqrt(disc));
        if (abs(q) < 1e-20) {
            r1 = 0.0;
            n = 1u;
        } else {
            r1 = q / a;
            r2 = c / q;
            n = 2u;
        }
    }
    var out = vec3<f32>(0.0, 0.0, 0.0);
    let t1 = shift + r1;
    let t2 = shift + r2;
    if (n == 2u) {
        let lo = min(t1, t2);
        let hi = max(t1, t2);
        if (lo > t_min) {
            return vec3<f32>(2.0, lo, hi);
        }
        if (hi > t_min) {
            return vec3<f32>(1.0, hi, 0.0);
        }
        return out;
    }
    if (t1 > t_min) {
        out = vec3<f32>(1.0, t1, 0.0);
    }
    return out;
}

// Where a ray meets patch `p` (number `index` in the competitor lists)
// past `t_min`, `t < 0` when it misses: exact port of `Patch::ray_roots`
// (solved from the ray's closest approach to the centre, see there for
// why), the nearer root inside the patch's mixed cell, and
// `Patch::normal_at` there.
struct SkinCast {
    t: f32,
    normal: vec3<f32>,
    atom: u32,
};

fn skin_cast(p: Patch, index: u32, ro: vec3<f32>, rd: vec3<f32>, t_min: f32, shrink: f32) -> SkinCast {
    let z = patch_center(p);
    let axis = patch_axis(p);
    let ab = coefficients(p.kind, shrink);
    let alpha = ab.x;
    let beta = ab.y;
    let t0 = dot(z - ro, rd);
    let u0 = ro + rd * t0 - z;
    let d0 = dot(u0, axis);
    let dd = dot(rd, axis);
    let qa = alpha + beta * dd * dd;
    let qb = 2.0 * (alpha * dot(u0, rd) + beta * d0 * dd);
    let qc = alpha * dot(u0, u0) + beta * d0 * d0 + p.s_weight;
    let roots = positive_roots(qa, qb, qc, t0, t_min);
    let count = u32(roots.x);
    var t = -1.0;
    if (count > 0u && contains(p, index, ro + rd * roots.y, shrink)) {
        t = roots.y;
    } else if (count == 2u && contains(p, index, ro + rd * roots.z, shrink)) {
        t = roots.z;
    }
    if (t < 0.0) {
        return SkinCast(-1.0, vec3<f32>(0.0), 0u);
    }
    let hit = ro + rd * t;
    let u = hit - z;
    let n = normalize(alpha * u + beta * dot(u, axis) * axis);
    return SkinCast(t, n, nearest_member(p, hit));
}
