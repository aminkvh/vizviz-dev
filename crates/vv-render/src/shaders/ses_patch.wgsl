// The solvent-excluded surface's patches, hit by a ray: the WGSL twin of
// `vv_core::ses::Ses::cast` (read `vv_core::ses` first; the geometry and
// its tests live there). Shared by the viewport (`ses_surface.wgsl`) and
// the path tracer (`path_trace.wgsl`), which each declare what it reads:
// `atoms` (xyz = position, w = van der Waals radius), `caps`, `probes`,
// `probe_neighbors`, `probe_radius()` and `atom_color(atom)` (linear).

// `vv_render::ses_surface::ProbeGpu`.
struct Probe {
    cx: f32,
    cy: f32,
    cz: f32,
    a0: u32,
    a1: u32,
    a2: u32,
    lo: u32,
    hi: u32,
};

const KIND_CONVEX: u32 = 0u;
const KIND_TORUS: u32 = 1u;
const FULL_CIRCLE: u32 = 0xffffffffu;
const TAU: f32 = 6.28318530718;
// `vv_core::ses::SEAM`: patches overlap slightly at their borders.
const SEAM: f32 = 1e-3;

fn kind_of(p: vec4<u32>) -> u32 {
    return p.x >> 30u;
}

fn first_of(p: vec4<u32>) -> u32 {
    return p.x & 0x3fffffffu;
}

fn probe_center(p: Probe) -> vec3<f32> {
    return vec3<f32>(p.cx, p.cy, p.cz);
}

// `vv_core::ses::circle` for atoms i, j's SAS spheres: center, normal,
// radius, offset (from atom i), separation.
struct Circle {
    center: vec3<f32>,
    normal: vec3<f32>,
    radius: f32,
    offset: f32,
    separation: f32,
};

fn circle(i: u32, j: u32) -> Circle {
    let rp = probe_radius();
    let ci = atoms[i].xyz;
    let ri = atoms[i].w + rp;
    let rj = atoms[j].w + rp;
    let v = atoms[j].xyz - ci;
    let d = length(v);
    let n = v / d;
    let offset = (d * d + ri * ri - rj * rj) / (2.0 * d);
    return Circle(ci + n * offset, n, sqrt(max(ri * ri - offset * offset, 0.0)), offset, d);
}

// `vv_core::ses::Ses::torus_arc`: where torus patch p's arc of probe
// centres starts (from its circle's centre) and how far it turns; a full
// circle turns TAU.
fn torus_arc(p: vec4<u32>, c: Circle) -> vec4<f32> {
    if (p.z == FULL_CIRCLE) {
        return vec4<f32>(0.0, 0.0, 0.0, TAU);
    }
    let start = probe_center(probes[p.z]) - c.center;
    var span = TAU;
    if (p.z != p.w) {
        span = angle_about(c.normal, start, probe_center(probes[p.w]) - c.center);
    }
    return vec4<f32>(start, span);
}

// `vv_core::ses::torus_bound`: a sphere (xyz, r) holding the torus
// patch between atoms i and j: between its touch points' heights, in the
// annular sector its arc sweeps.
fn torus_bound(c: Circle, i: u32, j: u32, arc: vec4<f32>) -> vec4<f32> {
    let rp = probe_radius();
    let ri = atoms[i].w;
    let rj = atoms[j].w;
    let reach = c.radius * max(ri / (ri + rp), rj / (rj + rp));
    let inner = max(c.radius - rp, 0.0);
    let zi = -c.offset * rp / (ri + rp);
    let zj = (c.separation - c.offset) * rp / (rj + rp);
    let lo = min(zi, zj);
    let half = 0.5 * (max(zi, zj) - lo);
    var center = c.center + c.normal * (lo + half);
    var planar = reach;
    if (arc.w < 3.14159265) {
        let s = sin(0.5 * arc.w);
        let co = cos(0.5 * arc.w);
        let e0 = normalize(arc.xyz);
        let bisector = e0 * co + cross(c.normal, e0) * s;
        let rho = 0.5 * (reach + inner * co);
        let outer_corner = length(vec2<f32>(reach * co - rho, reach * s));
        let inner_corner = length(vec2<f32>(inner * co - rho, inner * s));
        center = center + bisector * rho;
        planar = max(max(outer_corner, inner_corner), reach - rho);
    }
    return vec4<f32>(center, sqrt(planar * planar + half * half) + 1e-3);
}

// Bounding sphere (xyz, r) of patch `p`.
fn ses_bound(p: vec4<u32>) -> vec4<f32> {
    let rp = probe_radius();
    switch (kind_of(p)) {
        case KIND_CONVEX: {
            return atoms[first_of(p)];
        }
        case KIND_TORUS: {
            let c = circle(first_of(p), p.y);
            return torus_bound(c, first_of(p), p.y, torus_arc(p, c));
        }
        default: {
            return vec4<f32>(probe_center(probes[first_of(p)]), rp);
        }
    }
}

// Both roots of |o + t d - c|^2 = r^2 (unit d), ascending; x < 0 when none.
fn sphere_roots(o: vec3<f32>, d: vec3<f32>, c: vec3<f32>, r: f32) -> vec3<f32> {
    let oc = o - c;
    let b = dot(oc, d);
    let h = b * b - (dot(oc, oc) - r * r);
    if (h < 0.0) {
        return vec3<f32>(-1.0, 0.0, 0.0);
    }
    let s = sqrt(h);
    return vec3<f32>(1.0, -b - s, -b + s);
}

// Polynomials of degree <= 4: c.x + c.y t + c.z t^2 + c.w t^3 + c4 t^4.
fn peval(c: vec4<f32>, c4: f32, t: f32) -> f32 {
    return (((c4 * t + c.w) * t + c.z) * t + c.y) * t + c.x;
}

// `vv_core::ses::bracket_root`: the root in [a, b], where the polynomial
// is monotonic, as (found, t).
fn bracket(c: vec4<f32>, c4: f32, a_in: f32, b_in: f32) -> vec2<f32> {
    var a = a_in;
    var b = b_in;
    let fa = peval(c, c4, a);
    let fb = peval(c, c4, b);
    if (fa == 0.0) {
        return vec2<f32>(1.0, a);
    }
    if (fa * fb > 0.0) {
        return vec2<f32>(0.0, 0.0);
    }
    let rising = fb > fa;
    let dc = vec4<f32>(c.y, 2.0 * c.z, 3.0 * c.w, 4.0 * c4);
    var t = 0.5 * (a + b);
    for (var k = 0; k < 24; k = k + 1) {
        let f = peval(c, c4, t);
        if ((f > 0.0) == rising) {
            b = t;
        } else {
            a = t;
        }
        let slope = peval(dc, 0.0, t);
        let newton = t - f / slope;
        var next = 0.5 * (a + b);
        if (slope != 0.0 && newton > a && newton < b) {
            next = newton;
        }
        let step = abs(next - t);
        t = next;
        if (step < 1e-6) {
            break;
        }
    }
    return vec2<f32>(1.0, t);
}

// Roots of the monic quartic t^4 + a.w t^3 + a.z t^2 + a.y t + a.x in
// [lo, hi], ascending, as (count, roots): the quadratic q'' in closed
// form splits q' into monotonic pieces, q''s roots split q.
struct Roots {
    count: u32,
    t: array<f32, 4>,
};

fn quartic_roots(a: vec4<f32>, lo: f32, hi: f32) -> Roots {
    // q' = a.y + 2 a.z t + 3 a.w t^2 + 4 t^3; q'' = 2 a.z + 6 a.w t + 12 t^2.
    let d1 = vec4<f32>(a.y, 2.0 * a.z, 3.0 * a.w, 4.0);
    var cuts2 = array<f32, 4>(lo, hi, hi, hi);
    var n2 = 1u;
    let disc = 36.0 * a.w * a.w - 96.0 * a.z;
    if (disc > 0.0) {
        let s = sqrt(disc);
        let r1 = (-6.0 * a.w - s) / 24.0;
        let r2 = (-6.0 * a.w + s) / 24.0;
        if (r1 > lo && r1 < hi) {
            cuts2[n2] = r1;
            n2 = n2 + 1u;
        }
        if (r2 > lo && r2 < hi) {
            cuts2[n2] = r2;
            n2 = n2 + 1u;
        }
    }
    cuts2[n2] = hi;
    var cuts1 = array<f32, 5>(lo, hi, hi, hi, hi);
    var n1 = 1u;
    for (var k = 0u; k < n2; k = k + 1u) {
        let r = bracket(d1, 0.0, cuts2[k], cuts2[k + 1u]);
        if (r.x > 0.0) {
            cuts1[n1] = r.y;
            n1 = n1 + 1u;
        }
    }
    cuts1[n1] = hi;
    var out: Roots;
    out.count = 0u;
    for (var k = 0u; k < n1; k = k + 1u) {
        let r = bracket(a, 1.0, cuts1[k], cuts1[k + 1u]);
        if (r.x > 0.0 && out.count < 4u) {
            out.t[out.count] = r.y;
            out.count = out.count + 1u;
        }
    }
    return out;
}

// Angle of v about n, counter-clockwise from `start`, in [0, 2pi).
fn angle_about(n: vec3<f32>, start: vec3<f32>, v: vec3<f32>) -> f32 {
    let a = atan2(dot(n, cross(start, v)), dot(start, v));
    return select(a, a + TAU, a < 0.0);
}

struct Surf {
    ok: bool,
    t: f32,
    normal: vec3<f32>,
    color: vec3<f32>,
    atom: u32,
};

fn miss() -> Surf {
    return Surf(false, 0.0, vec3<f32>(0.0), vec3<f32>(0.0), 0u);
}

// The bound radius of the whole circle's torus between atoms i and j: the
// one size the viewport judges that pair's fillets by.
fn fillet_size(i: u32, j: u32) -> f32 {
    return torus_bound(circle(i, j), i, j, vec4<f32>(0.0, 0.0, 0.0, TAU)).w;
}

// Nearest hits past `t_min`. Caps whose fillets are smaller than
// `min_fillet` are ignored, and all of them for a `whole` patch: the
// viewport skips those fillets, so the sphere under them shows instead of
// a hole. The path tracer passes (0, false).
fn hit_convex(p: vec4<u32>, o: vec3<f32>, d: vec3<f32>, t_min: f32, min_fillet: f32, whole: bool) -> Surf {
    let i = first_of(p);
    let c = atoms[i].xyz;
    let r = atoms[i].w;
    let ri = r + probe_radius();
    let roots = sphere_roots(o, d, c, r);
    if (roots.x < 0.0) {
        return miss();
    }
    for (var k = 0; k < 2; k = k + 1) {
        let t = select(roots.z, roots.y, k == 0);
        if (t < t_min) {
            continue;
        }
        let u = (o + d * t - c) / r;
        var free = true;
        for (var m = p.y; m < select(p.z, p.y, whole); m = m + 1u) {
            // `vv_core::ses::cap`.
            let j = caps[m];
            let v = atoms[j].xyz - c;
            let dist = length(v);
            let rj = atoms[j].w + probe_radius();
            let in_cap = dot(u, v) / dist > (dist * dist + ri * ri - rj * rj) / (2.0 * dist * ri) + SEAM;
            if (in_cap && !(min_fillet > 0.0 && fillet_size(i, j) < min_fillet)) {
                free = false;
                break;
            }
        }
        if (free) {
            return Surf(true, t, u, atom_color(i), i);
        }
    }
    return miss();
}

fn hit_torus(p: vec4<u32>, o: vec3<f32>, d: vec3<f32>, t_min: f32) -> Surf {
    let rp = probe_radius();
    let i = first_of(p);
    let j = p.y;
    let c = circle(i, j);
    let scale = c.radius + rp;
    let arc = torus_arc(p, c);
    let bound = torus_bound(c, i, j, arc);
    let chord = sphere_roots(o, d, bound.xyz, bound.w);
    if (chord.x < 0.0) {
        return miss();
    }
    // `vv_core::ses::torus_roots`: solve on the patch bound's chord,
    // shifted to its midpoint and scaled by the torus's size.
    let mid = 0.5 * (chord.y + chord.z);
    let oo = (o + d * mid - c.center) / scale;
    let r = c.radius / scale;
    let q = rp / scale;
    let b = dot(oo, d);
    let oz = dot(oo, c.normal);
    let dz = dot(d, c.normal);
    let o2 = dot(oo, oo);
    let e = o2 + r * r - q * q;
    let r4 = 4.0 * r * r;
    let a = vec4<f32>(
        e * e - r4 * (o2 - oz * oz),
        4.0 * b * e - 2.0 * r4 * (b - oz * dz),
        4.0 * b * b + 2.0 * e - r4 * (1.0 - dz * dz),
        4.0 * b,
    );
    let half = 0.5 * (chord.z - chord.y) / scale;
    let roots = quartic_roots(a, -half, half);
    for (var k = 0u; k < roots.count; k = k + 1u) {
        let t = mid + roots.t[k] * scale;
        if (t < t_min) {
            continue;
        }
        let x = o + d * t;
        let v = x - c.center;
        let z = dot(v, c.normal);
        let radial = v - c.normal * z;
        let s = length(radial);
        var out_dir = radial / max(s, 1e-6);
        if (s <= 1e-6) {
            let helper = select(vec3<f32>(1.0, 0.0, 0.0), vec3<f32>(0.0, 1.0, 0.0), abs(c.normal.x) >= 0.9);
            out_dir = normalize(cross(c.normal, helper));
        }
        let u = vec2<f32>(s - c.radius, z);
        // The spindle torus's inner lemon: not surface.
        if (abs(length(u) - rp) > 1e-3 * max(rp, 1.0) + 1e-3) {
            continue;
        }
        let to_i = vec2<f32>(-c.radius, -c.offset);
        let to_j = vec2<f32>(-c.radius, c.separation - c.offset);
        if (to_i.x * u.y - to_i.y * u.x > SEAM * rp * length(to_i)
            || u.x * to_j.y - u.y * to_j.x > SEAM * rp * length(to_j)) {
            continue;
        }
        if (p.z != FULL_CIRCLE) {
            let turn = angle_about(c.normal, arc.xyz, out_dir);
            if (turn > arc.w + SEAM && turn < TAU - SEAM) {
                continue;
            }
        }
        let probe_at = c.center + out_dir * c.radius;
        // Colour runs from atom i to atom j across the saddle.
        let w = clamp((z + c.offset) / c.separation, 0.0, 1.0);
        let color = mix(atom_color(i), atom_color(j), w);
        return Surf(true, t, (probe_at - x) / rp, color, select(i, j, w > 0.5));
    }
    return miss();
}

fn hit_concave(p: vec4<u32>, o: vec3<f32>, d: vec3<f32>, t_min: f32) -> Surf {
    let rp = probe_radius();
    let pr = probes[first_of(p)];
    let x = probe_center(pr);
    let a = atoms[pr.a0].xyz - x;
    let b = atoms[pr.a1].xyz - x;
    let c = atoms[pr.a2].xyz - x;
    let det = dot(a, cross(b, c));
    let roots = sphere_roots(o, d, x, rp);
    if (roots.x < 0.0) {
        return miss();
    }
    for (var k = 0; k < 2; k = k + 1) {
        let t = select(roots.z, roots.y, k == 0);
        if (t < t_min) {
            continue;
        }
        let hit = o + d * t;
        let u = hit - x;
        let w = vec3<f32>(dot(u, cross(b, c)), dot(a, cross(u, c)), dot(a, cross(b, u))) / det;
        if (min(w.x, min(w.y, w.z)) < -SEAM) {
            continue;
        }
        var covered = false;
        for (var m = pr.lo; m < pr.hi; m = m + 1u) {
            let y = probe_center(probes[probe_neighbors[m]]);
            if (dot(hit - y, hit - y) < rp * rp - 1e-4) {
                covered = true;
                break;
            }
        }
        if (covered) {
            continue;
        }
        let wn = max(w, vec3<f32>(0.0)) / max(w.x + w.y + w.z, 1e-6);
        let color = atom_color(pr.a0) * wn.x + atom_color(pr.a1) * wn.y + atom_color(pr.a2) * wn.z;
        var atom = pr.a0;
        if (wn.y > wn.x && wn.y >= wn.z) {
            atom = pr.a1;
        } else if (wn.z > wn.x && wn.z > wn.y) {
            atom = pr.a2;
        }
        return Surf(true, t, (x - hit) / rp, color, atom);
    }
    return miss();
}
