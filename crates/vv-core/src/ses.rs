//! The solvent-excluded surface (Richards' SES): the boundary of what a
//! probe sphere rolling over the van der Waals spheres cannot reach.
//!
//! Built analytically as three kinds of patch, after Quan & Stamm (2016)
//! and Plateau-Holleville et al. (TVCG 2024), from the geometry of the
//! solvent-accessible surface (SAS: each atom's sphere grown by the probe
//! radius, `R = r + probe`):
//!
//! - **convex**: part of an atom's van der Waals sphere, minus the caps its
//!   neighbours' SAS spheres cover (stored as cones: an axis and an angle);
//! - **toroidal**: where the probe touches two atoms, swept along the free
//!   arcs of the circle where their SAS spheres meet; points only count
//!   inside the triangle (probe centre, both atom centres), which is what
//!   removes the self-intersecting "spindle" part without special cases;
//! - **concave**: where the probe touches three atoms, a piece of the probe
//!   sphere inside the cone to the three centres, minus overlapping probe
//!   spheres.
//!
//! Everything is complete (cavities included) and is ray-cast directly --
//! [`Ses::cast`] is the CPU twin of the renderer's shader and the test
//! oracle's counterpart. Radii passed in are van der Waals radii.

use glam::{Vec2, Vec3, Vec4};
use rayon::prelude::*;

use crate::spatial::Grid;

/// Water, in Å: the usual probe radius.
pub const WATER_PROBE: f32 = 1.4;

/// `Torus::probes` for a circle free all the way round.
pub const FULL_CIRCLE: u32 = u32::MAX;

/// Squared-distance slack for "strictly inside another SAS sphere", in
/// Å²: a probe position touching a fourth sphere (co-spherical atoms) is
/// kept, not lost to rounding.
const INSIDE_EPS: f32 = 1e-3;

/// How far each patch's border test reaches past the border (as a cosine,
/// a barycentric weight, radians, or a sine), so neighbouring patches
/// overlap by a few thousandths of an Å instead of leaving f32 pinholes.
const SEAM: f32 = 1e-3;

/// A probe position touching three atoms' SAS spheres and no other.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Probe {
    pub center: Vec3,
    /// Ascending.
    pub atoms: [u32; 3],
}

/// A toroidal patch between atoms `atoms[0] < atoms[1]`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Torus {
    pub atoms: [u32; 2],
    /// The arc of probe centres, counter-clockwise about
    /// `c[atoms[1]] - c[atoms[0]]` from `probes[0]` to `probes[1]` (equal
    /// for a circle with one free gap); both [`FULL_CIRCLE`] when free all
    /// the way round.
    pub probes: [u32; 2],
}

/// A convex patch: atom `atom`'s sphere outside the caps its neighbours
/// `Ses::caps[start..end]` cut off.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Convex {
    pub atom: u32,
    pub caps: [u32; 2],
}

#[derive(Clone, Debug, Default)]
pub struct Ses {
    pub probe_radius: f32,
    pub convex: Vec<Convex>,
    /// Neighbour atoms whose SAS spheres cut a convex patch, widest cap
    /// first per patch; the cap's cone is recomputed from the two atoms
    /// ([`cap`]).
    pub caps: Vec<u32>,
    pub tori: Vec<Torus>,
    pub probes: Vec<Probe>,
    /// CSR: probe `p`'s overlapping probes (centres closer than two probe
    /// radii) are `probe_neighbors[probe_neighbor_starts[p]..[p + 1]]`.
    pub probe_neighbor_starts: Vec<u32>,
    pub probe_neighbors: Vec<u32>,
}

/// Where two SAS spheres meet, seen from the first atom.
#[derive(Clone, Copy, Debug)]
pub struct Circle {
    pub center: Vec3,
    /// Unit, from the first atom toward the second.
    pub normal: Vec3,
    pub radius: f32,
    /// Signed distance from the first atom's centre to the circle's plane.
    pub offset: f32,
    /// Distance between the two centres.
    pub separation: f32,
}

/// The circle where spheres `(ci, ri)` and `(cj, rj)` meet, or `None` when
/// they miss or one holds the other.
pub fn circle(ci: Vec3, ri: f32, cj: Vec3, rj: f32) -> Option<Circle> {
    let v = cj - ci;
    let d = v.length();
    if d <= 0.0 || d >= ri + rj || d + rj <= ri || d + ri <= rj {
        return None;
    }
    let normal = v / d;
    let offset = (d * d + ri * ri - rj * rj) / (2.0 * d);
    let r2 = ri * ri - offset * offset;
    (r2 > 0.0).then(|| Circle {
        center: ci + normal * offset,
        normal,
        radius: r2.sqrt(),
        offset,
        separation: d,
    })
}

/// The cap sphere `(cj, rj)` cuts off sphere `(ci, ri)`, as a cone from
/// `ci`: `xyz` the unit axis, `w` the cosine of its half-angle (a
/// direction `u` is inside when `u . axis > w`).
pub fn cap(ci: Vec3, ri: f32, cj: Vec3, rj: f32) -> Vec4 {
    let v = cj - ci;
    let d = v.length();
    (v / d).extend((d * d + ri * ri - rj * rj) / (2.0 * d * ri))
}

/// Two unit vectors spanning the plane normal to `n`.
fn basis(n: Vec3) -> (Vec3, Vec3) {
    let helper = if n.x.abs() < 0.9 { Vec3::X } else { Vec3::Y };
    let e1 = n.cross(helper).normalize();
    (e1, n.cross(e1))
}

/// Angle of `v` about `n`, counter-clockwise from `from`, in `[0, 2pi)`.
fn angle_about(n: Vec3, from: Vec3, v: Vec3) -> f32 {
    let a = n.dot(from.cross(v)).atan2(from.dot(v));
    if a < 0.0 {
        a + std::f32::consts::TAU
    } else {
        a
    }
}

struct Neighbor {
    atom: u32,
    center: Vec3,
    sas: f32,
}

/// Atom `i`'s SAS neighbours (spheres that cut or touch its own), and
/// whether its SAS sphere lies inside another's (then it adds nothing).
fn neighbors_of(
    i: usize,
    positions: &[Vec3],
    sas: &[f32],
    sas_max: f32,
    grid: &Grid,
    out: &mut Vec<Neighbor>,
) -> bool {
    out.clear();
    let (ci, ri) = (positions[i], sas[i]);
    let mut covered = false;
    grid.for_each_within(positions, ci, ri + sas_max, |j, d2| {
        let j = j as usize;
        if j == i {
            return;
        }
        let rj = sas[j];
        let d = d2.sqrt();
        if d >= ri + rj {
            return;
        }
        // Inside another sphere (ties go to the lower index, so of two
        // identical atoms one survives).
        if d + ri < rj || (d + ri <= rj && j < i) {
            covered = true;
        }
        out.push(Neighbor {
            atom: j as u32,
            center: positions[j],
            sas: rj,
        });
    });
    covered
}

fn inside_any(p: Vec3, neighbors: &[Neighbor], skip: [u32; 2], eps: f32) -> bool {
    neighbors.iter().any(|n| {
        n.atom != skip[0] && n.atom != skip[1] && p.distance_squared(n.center) < n.sas * n.sas - eps
    })
}

#[derive(Default)]
struct AtomPass {
    covered: bool,
    has_circle: bool,
    full: Vec<u32>,
    probes: Vec<Probe>,
}

/// Builds the SES of spheres at `positions` with van der Waals `radii`.
pub fn build(positions: &[Vec3], radii: &[f32], probe_radius: f32) -> Ses {
    assert_eq!(positions.len(), radii.len());
    assert!(probe_radius > 0.0);
    let n = positions.len();
    if n == 0 {
        return Ses {
            probe_radius,
            ..Ses::default()
        };
    }
    let sas: Vec<f32> = radii.iter().map(|r| r + probe_radius).collect();
    let sas_max = sas.iter().copied().fold(0.0, f32::max);
    let indices: Vec<u32> = (0..n as u32).collect();
    let grid = Grid::build(positions, &indices, 2.0 * sas_max);

    // Pass 1, per atom: classify its circles with higher-indexed atoms,
    // and find the probe positions it is the lowest atom of.
    let pass: Vec<AtomPass> = (0..n)
        .into_par_iter()
        .map_init(Vec::new, |nb, i| {
            let mut out = AtomPass {
                covered: neighbors_of(i, positions, &sas, sas_max, &grid, nb),
                ..AtomPass::default()
            };
            if out.covered {
                return out;
            }
            let (ci, ri) = (positions[i], sas[i]);
            nb.sort_unstable_by(|a, b| {
                (a.center.distance_squared(ci) - a.sas * a.sas)
                    .total_cmp(&(b.center.distance_squared(ci) - b.sas * b.sas))
            });
            for j in nb.iter() {
                let Some(c) = circle(ci, ri, j.center, j.sas) else {
                    continue;
                };
                out.has_circle = true;
                if j.atom < i as u32 {
                    continue;
                }
                let mut buried = false;
                let mut cut = false;
                for k in nb.iter() {
                    if k.atom == j.atom {
                        continue;
                    }
                    let w = k.center - c.center;
                    let h = w.dot(c.normal);
                    let rho = (w.length_squared() - h * h).max(0.0).sqrt();
                    let rk2 = k.sas * k.sas;
                    if h * h + (rho + c.radius) * (rho + c.radius) < rk2 {
                        buried = true;
                        break;
                    }
                    cut |= h * h + (rho - c.radius) * (rho - c.radius) < rk2;
                }
                if buried {
                    continue;
                }
                if !cut {
                    out.full.push(j.atom);
                    continue;
                }
                // Where circle (i, j) crosses sphere k, for k > j.
                for k in nb.iter() {
                    if k.atom <= j.atom || k.center.distance(j.center) >= j.sas + k.sas {
                        continue;
                    }
                    let w = k.center - c.center;
                    let h = w.dot(c.normal);
                    let q = w - c.normal * h;
                    let rho = q.length();
                    if rho <= 0.0 {
                        continue;
                    }
                    let cos_phi = (h * h + c.radius * c.radius + rho * rho - k.sas * k.sas)
                        / (2.0 * c.radius * rho);
                    if cos_phi.abs() >= 1.0 {
                        continue;
                    }
                    let sin_phi = (1.0 - cos_phi * cos_phi).sqrt();
                    let qh = q / rho;
                    let side = c.normal.cross(qh);
                    for s in [-1.0, 1.0] {
                        let x = c.center + (qh * cos_phi + side * (s * sin_phi)) * c.radius;
                        if !inside_any(x, nb, [j.atom, k.atom], INSIDE_EPS) {
                            out.probes.push(Probe {
                                center: x,
                                atoms: [i as u32, j.atom, k.atom],
                            });
                        }
                    }
                }
            }
            out
        })
        .collect();

    let mut tori: Vec<Torus> = Vec::new();
    let mut probes: Vec<Probe> = Vec::new();
    let mut covered = vec![false; n];
    let mut has_circle = vec![false; n];
    for (i, p) in pass.into_iter().enumerate() {
        covered[i] = p.covered;
        has_circle[i] = p.has_circle;
        tori.extend(p.full.into_iter().map(|j| Torus {
            atoms: [i as u32, j],
            probes: [FULL_CIRCLE; 2],
        }));
        probes.extend(p.probes);
    }

    // Pass 2: each cut circle's probes, in angular order, bound its arcs;
    // an arc is free when its midpoint is.
    let mut refs: Vec<(u32, u32, u32)> = Vec::with_capacity(probes.len() * 3);
    for (index, p) in probes.iter().enumerate() {
        let [a, b, c] = p.atoms;
        let index = index as u32;
        refs.extend([(a, b, index), (a, c, index), (b, c, index)]);
    }
    refs.par_sort_unstable();
    let mut groups: Vec<std::ops::Range<usize>> = Vec::new();
    let mut start = 0;
    for k in 1..=refs.len() {
        if k == refs.len() || refs[k].0 != refs[start].0 {
            groups.push(start..k);
            start = k;
        }
    }
    let segments: Vec<Vec<Torus>> = groups
        .into_par_iter()
        .map_init(Vec::new, |nb, range| {
            let group = &refs[range];
            let i = group[0].0 as usize;
            neighbors_of(i, positions, &sas, sas_max, &grid, nb);
            let mut out = Vec::new();
            let mut angles: Vec<(f32, u32)> = Vec::new();
            for run in group.chunk_by(|a, b| a.1 == b.1) {
                let j = run[0].1 as usize;
                let Some(c) = circle(positions[i], sas[i], positions[j], sas[j]) else {
                    continue;
                };
                let from = probes[run[0].2 as usize].center - c.center;
                let e1 = from.normalize();
                let e2 = c.normal.cross(e1);
                angles.clear();
                angles.extend(run.iter().map(|&(_, _, p)| {
                    (
                        angle_about(c.normal, from, probes[p as usize].center - c.center),
                        p,
                    )
                }));
                angles.sort_by(|a, b| a.0.total_cmp(&b.0));
                for (s, &(a0, p0)) in angles.iter().enumerate() {
                    let (a1, p1) = angles[(s + 1) % angles.len()];
                    let a1 = if s + 1 == angles.len() {
                        a1 + std::f32::consts::TAU
                    } else {
                        a1
                    };
                    let mid = 0.5 * (a0 + a1);
                    let x = c.center + (e1 * mid.cos() + e2 * mid.sin()) * c.radius;
                    if a1 > a0 && !inside_any(x, nb, [j as u32, j as u32], 0.0) {
                        out.push(Torus {
                            atoms: [i as u32, j as u32],
                            probes: [p0, p1],
                        });
                    }
                }
            }
            out
        })
        .collect();
    tori.extend(segments.into_iter().flatten());
    tori.sort_unstable_by_key(|t| (t.atoms, t.probes));

    // Pass 3: a convex patch for every atom with some SAS left -- one on a
    // torus, or one no other sphere touches -- cut by the caps of the
    // circles bounding it. Caps of circles with no free arc are left out:
    // what they alone cover is sealed inside the excluded volume, drawn
    // only through glass or a clip plane, and they are most of the caps.
    let mut pairs: Vec<(u32, u32, f32)> = tori
        .iter()
        .flat_map(|t| [(t.atoms[0], t.atoms[1]), (t.atoms[1], t.atoms[0])])
        .map(|(i, j)| {
            let (i, j) = (i as usize, j as usize);
            let offset =
                circle(positions[i], sas[i], positions[j], sas[j]).map_or(f32::MAX, |c| c.offset);
            (i as u32, j as u32, offset)
        })
        .collect();
    // Widest cap (smallest offset) first, so covered points fail early.
    pairs.par_sort_unstable_by(|a, b| (a.0, a.2, a.1).partial_cmp(&(b.0, b.2, b.1)).unwrap());
    pairs.dedup_by_key(|p| (p.0, p.1));
    let mut convex = Vec::new();
    let mut caps = Vec::with_capacity(pairs.len());
    let mut next = pairs.iter().peekable();
    for i in 0..n as u32 {
        let start = caps.len() as u32;
        while let Some(&(_, j, _)) = next.next_if(|p| p.0 == i) {
            caps.push(j);
        }
        let isolated = !covered[i as usize] && !has_circle[i as usize];
        if caps.len() as u32 > start || isolated {
            convex.push(Convex {
                atom: i,
                caps: [start, caps.len() as u32],
            });
        }
    }

    // Pass 4: overlapping probe spheres, which cut each other's patches --
    // only where a probe sphere dips through its atoms' plane can another
    // reach into its cone.
    let centers: Vec<Vec3> = probes.iter().map(|p| p.center).collect();
    let probe_index: Vec<u32> = (0..centers.len() as u32).collect();
    let (probe_neighbor_starts, probe_neighbors) = if centers.is_empty() {
        (vec![0], Vec::new())
    } else {
        let probe_grid = Grid::build(&centers, &probe_index, 2.0 * probe_radius);
        let lists: Vec<Vec<u32>> = probes
            .par_iter()
            .enumerate()
            .map(|(p, probe)| {
                let [a, b, c] = probe.atoms.map(|k| positions[k as usize]);
                let plane = (b - a).cross(c - a).normalize();
                let mut list = Vec::new();
                if (probe.center - a).dot(plane).abs() < probe_radius {
                    probe_grid.for_each_within(
                        &centers,
                        probe.center,
                        2.0 * probe_radius,
                        |q, d2| {
                            if q as usize != p && d2 < 4.0 * probe_radius * probe_radius {
                                list.push(q);
                            }
                        },
                    );
                    list.sort_unstable();
                }
                list
            })
            .collect();
        let mut starts = Vec::with_capacity(lists.len() + 1);
        let mut flat = Vec::new();
        starts.push(0);
        for list in lists {
            flat.extend(list);
            starts.push(flat.len() as u32);
        }
        (starts, flat)
    };

    Ses {
        probe_radius,
        convex,
        caps,
        tori,
        probes,
        probe_neighbor_starts,
        probe_neighbors,
    }
}

/// What a ray hit.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum PatchRef {
    Convex(u32),
    Torus(u32),
    Concave(u32),
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Hit {
    pub t: f32,
    /// Unit, pointing out of the excluded volume (into the solvent).
    pub normal: Vec3,
    pub patch: PatchRef,
}

/// Both roots of `|o + t d - c|^2 = r^2` for unit `d`, ascending.
fn sphere_roots(o: Vec3, d: Vec3, c: Vec3, r: f32) -> Option<[f32; 2]> {
    let oc = o - c;
    let b = oc.dot(d);
    let h = b * b - (oc.length_squared() - r * r);
    (h >= 0.0).then(|| {
        let h = h.sqrt();
        [-b - h, -b + h]
    })
}

impl Ses {
    /// First patch hit at `t >= t_min` along the ray `origin + t dir`
    /// (`dir` unit), testing every patch: the reference the shader and
    /// tests are checked against, not a fast path.
    pub fn cast(
        &self,
        positions: &[Vec3],
        radii: &[f32],
        origin: Vec3,
        dir: Vec3,
        t_min: f32,
    ) -> Option<Hit> {
        let mut best: Option<Hit> = None;
        let mut take = |t: f32, normal: Vec3, patch: PatchRef| {
            if best.is_none_or(|b| t < b.t) {
                best = Some(Hit { t, normal, patch });
            }
        };
        for (k, c) in self.convex.iter().enumerate() {
            if let Some((t, n)) = self.hit_convex(c, positions, radii, origin, dir, t_min) {
                take(t, n, PatchRef::Convex(k as u32));
            }
        }
        for (k, torus) in self.tori.iter().enumerate() {
            if let Some((t, n)) = self.hit_torus(torus, positions, radii, origin, dir, t_min) {
                take(t, n, PatchRef::Torus(k as u32));
            }
        }
        for k in 0..self.probes.len() {
            if let Some((t, n)) = self.hit_concave(k, positions, origin, dir, t_min) {
                take(t, n, PatchRef::Concave(k as u32));
            }
        }
        best
    }

    pub fn hit_convex(
        &self,
        patch: &Convex,
        positions: &[Vec3],
        radii: &[f32],
        o: Vec3,
        d: Vec3,
        t_min: f32,
    ) -> Option<(f32, Vec3)> {
        let i = patch.atom as usize;
        let (c, r) = (positions[i], radii[i]);
        let rp = self.probe_radius;
        let cones: Vec<Vec4> = self.caps[patch.caps[0] as usize..patch.caps[1] as usize]
            .iter()
            .map(|&j| cap(c, r + rp, positions[j as usize], radii[j as usize] + rp))
            .collect();
        sphere_roots(o, d, c, r)?
            .into_iter()
            .filter(|&t| t >= t_min)
            .find_map(|t| {
                let u = (o + d * t - c) / r;
                cones
                    .iter()
                    .all(|s| u.dot(s.truncate()) <= s.w + SEAM)
                    .then_some((t, u))
            })
    }

    /// Where `torus`'s arc of probe centres starts (from the centre of its
    /// circle `c`) and how far it turns about `c.normal`; `None` for a full
    /// circle.
    pub fn torus_arc(&self, torus: &Torus, c: &Circle) -> Option<(Vec3, f32)> {
        (torus.probes[0] != FULL_CIRCLE).then(|| {
            let a = self.probes[torus.probes[0] as usize].center - c.center;
            let b = self.probes[torus.probes[1] as usize].center - c.center;
            let span = if torus.probes[0] == torus.probes[1] {
                std::f32::consts::TAU
            } else {
                angle_about(c.normal, a, b)
            };
            (a, span)
        })
    }

    pub fn hit_torus(
        &self,
        torus: &Torus,
        positions: &[Vec3],
        radii: &[f32],
        o: Vec3,
        d: Vec3,
        t_min: f32,
    ) -> Option<(f32, Vec3)> {
        let rp = self.probe_radius;
        let [i, j] = torus.atoms.map(|a| a as usize);
        let c = circle(positions[i], radii[i] + rp, positions[j], radii[j] + rp)?;
        let arc = self.torus_arc(torus, &c);
        let (center, radius) = torus_bound(&c, radii[i], radii[j], rp, arc);
        let window = sphere_roots(o, d, center, radius)?;
        torus_roots(o, d, c.center, c.normal, c.radius, rp, window)
            .into_iter()
            .flatten()
            .filter(|&t| t >= t_min)
            .find_map(|t| {
                let p = o + d * t;
                let v = p - c.center;
                let z = v.dot(c.normal);
                let radial = v - c.normal * z;
                let s = radial.length();
                let out = if s > 1e-6 {
                    radial / s
                } else {
                    basis(c.normal).0
                };
                // Inside the triangle (probe centre, both atom centres),
                // in the half-plane through the axis and `p`.
                let u = Vec2::new(s - c.radius, z);
                // A spindle torus's quartic also vanishes on the inner
                // "lemon", one probe radius from the far side of the circle
                // and inside the near side's probe sphere: not surface.
                if (u.length() - rp).abs() > 1e-3 * rp.max(1.0) + 1e-3 {
                    return None;
                }
                let to_i = Vec2::new(-c.radius, -c.offset);
                let to_j = Vec2::new(-c.radius, c.separation - c.offset);
                if to_i.perp_dot(u) > SEAM * rp * to_i.length()
                    || u.perp_dot(to_j) > SEAM * rp * to_j.length()
                {
                    return None;
                }
                if let Some((from, span)) = arc {
                    let turn = angle_about(c.normal, from, out);
                    if turn > span + SEAM && turn < std::f32::consts::TAU - SEAM {
                        return None;
                    }
                }
                let probe_center = c.center + out * c.radius;
                Some((t, (probe_center - p) / rp))
            })
    }

    /// Where probe `probe`'s sphere bounds the excluded volume: inside the
    /// cone from its centre through its three atoms, and in no overlapping
    /// probe's sphere.
    pub fn hit_concave(
        &self,
        probe: usize,
        positions: &[Vec3],
        o: Vec3,
        d: Vec3,
        t_min: f32,
    ) -> Option<(f32, Vec3)> {
        let rp = self.probe_radius;
        let x = self.probes[probe].center;
        let others = &self.probe_neighbors[self.probe_neighbor_starts[probe] as usize
            ..self.probe_neighbor_starts[probe + 1] as usize];
        sphere_roots(o, d, x, rp)?
            .into_iter()
            .filter(|&t| t >= t_min)
            .find_map(|t| {
                let p = o + d * t;
                let covered = others
                    .iter()
                    .any(|&q| p.distance_squared(self.probes[q as usize].center) < rp * rp - 1e-4);
                (!covered && self.in_cone(probe, positions, p)).then_some((t, (x - p) / rp))
            })
    }

    /// How far `p` is from the nearest place a probe centre can be (outside
    /// every SAS sphere): zero outside the SAS, and more than the probe
    /// radius exactly inside the excluded volume, so the SES is where it
    /// equals the probe radius. Inside the SAS the nearest such place lies
    /// on the SAS boundary: on a sphere where no other sphere covers it, or on a
    /// free arc of a circle where two spheres meet (whose ends are the
    /// probe positions). Brute force over every patch: the reference the
    /// renderer's volume bake (`shaders/ses_volume.wgsl`) follows.
    pub fn probe_distance(&self, positions: &[Vec3], radii: &[f32], p: Vec3) -> f32 {
        let rp = self.probe_radius;
        let in_sas = positions
            .iter()
            .zip(radii)
            .any(|(c, r)| p.distance_squared(*c) < (r + rp) * (r + rp));
        if !in_sas {
            return 0.0;
        }
        // Every sphere, tested against every other: convex patches' caps
        // leave out what is sealed inside the excluded volume, which is not
        // SAS boundary.
        let mut best = f32::INFINITY;
        for i in 0..positions.len() {
            let (c, big) = (positions[i], radii[i] + rp);
            let r = p.distance(c);
            if r < 1e-6 || (big - r).abs() >= best {
                continue;
            }
            let q = c + (p - c) * (big / r);
            let open = positions
                .iter()
                .zip(radii)
                .enumerate()
                .all(|(j, (cj, rj))| {
                    j == i || q.distance_squared(*cj) >= (rj + rp) * (rj + rp) - INSIDE_EPS
                });
            if open {
                best = (big - r).abs();
            }
        }
        for torus in &self.tori {
            let [i, j] = torus.atoms.map(|a| a as usize);
            if let Some(c) = circle(positions[i], radii[i] + rp, positions[j], radii[j] + rp) {
                best = best.min(self.arc_distance(torus, &c, p));
            }
        }
        best
    }

    /// Distance from `p` to `torus`'s free arc of probe centres on its
    /// circle `c`.
    fn arc_distance(&self, torus: &Torus, c: &Circle, p: Vec3) -> f32 {
        let v = p - c.center;
        let z = v.dot(c.normal);
        let radial = v - c.normal * z;
        let s = radial.length();
        let on_circle = ((s - c.radius) * (s - c.radius) + z * z).sqrt();
        match self.torus_arc(torus, c) {
            Some((from, span)) if s > 1e-6 && angle_about(c.normal, from, radial) > span => torus
                .probes
                .map(|k| p.distance(self.probes[k as usize].center))
                .into_iter()
                .fold(f32::INFINITY, f32::min),
            _ => on_circle,
        }
    }

    /// `p` inside the cone from probe `probe`'s centre through its atoms.
    fn in_cone(&self, probe: usize, positions: &[Vec3], p: Vec3) -> bool {
        let x = self.probes[probe].center;
        let [a, b, c] = self.probes[probe].atoms.map(|k| positions[k as usize] - x);
        let u = p - x;
        let w =
            Vec3::new(u.dot(b.cross(c)), a.dot(u.cross(c)), a.dot(b.cross(u))) / a.dot(b.cross(c));
        w.min_element() >= -SEAM
    }
}

/// A sphere holding torus patch `c`'s surface (the circle of SAS spheres
/// of van der Waals radii `ri`, `rj`; `arc` as [`Ses::torus_arc`]): the
/// patch runs from its touch point on one atom to the other, bulging
/// toward the axis, so it lies between the touch points' heights and, in
/// the plane, in the annular sector the arc sweeps between `c.radius - rp`
/// and the touch points' distance from the axis -- far tighter than the
/// whole torus, which is what the billboards and the cull are sized from.
pub fn torus_bound(c: &Circle, ri: f32, rj: f32, rp: f32, arc: Option<(Vec3, f32)>) -> (Vec3, f32) {
    let (sas_i, sas_j) = (ri + rp, rj + rp);
    let reach = c.radius * (ri / sas_i).max(rj / sas_j);
    let inner = (c.radius - rp).max(0.0);
    let zi = -c.offset * rp / sas_i;
    let zj = (c.separation - c.offset) * rp / sas_j;
    let (lo, hi) = (zi.min(zj), zi.max(zj));
    let half = 0.5 * (hi - lo);
    let axial = c.center + c.normal * (lo + half);
    let (shift, planar) = match arc {
        Some((from, span)) if span < std::f32::consts::PI => {
            // On the sector's bisector, between its inner corners and its
            // outer arc; the farthest points are then the corners and the
            // arc's middle.
            let (sin, cos) = (0.5 * span).sin_cos();
            let e0 = from.normalize();
            let bisector = e0 * cos + c.normal.cross(e0) * sin;
            let rho = 0.5 * (reach + inner * cos);
            let outer_corner = Vec2::new(reach * cos - rho, reach * sin).length();
            let inner_corner = Vec2::new(inner * cos - rho, inner * sin).length();
            (
                bisector * rho,
                outer_corner.max(inner_corner).max(reach - rho),
            )
        }
        _ => (Vec3::ZERO, reach),
    };
    (axial + shift, (planar * planar + half * half).sqrt() + 1e-3)
}

/// Up to four ascending roots of the ray `o + t d` (unit `d`) against the
/// torus of centre `c`, axis `n`, major radius `major`, tube `minor`, for
/// `t` in `window` (a chord of a sphere around what is wanted).
///
/// The quartic is solved on that chord, shifted to its midpoint and
/// scaled by the torus's size (so f32 holds up far from the origin), by
/// bracketing: the quadratic's closed-form roots split the cubic
/// derivative into monotonic pieces, the cubic's roots split the quartic,
/// and each sign change is refined with safeguarded Newton.
pub fn torus_roots(
    o: Vec3,
    d: Vec3,
    c: Vec3,
    n: Vec3,
    major: f32,
    minor: f32,
    window: [f32; 2],
) -> [Option<f32>; 4] {
    let scale = major + minor;
    let [t0, t1] = window;
    let mid = 0.5 * (t0 + t1);
    let oo = (o + d * mid - c) / scale;
    let (r, rp) = (major / scale, minor / scale);
    let b = oo.dot(d);
    let oz = oo.dot(n);
    let dz = d.dot(n);
    let o2 = oo.length_squared();
    let e = o2 + r * r - rp * rp;
    let r4 = 4.0 * r * r;
    // t^4 + a3 t^3 + a2 t^2 + a1 t + a0
    let a3 = 4.0 * b;
    let a2 = 4.0 * b * b + 2.0 * e - r4 * (1.0 - dz * dz);
    let a1 = 4.0 * b * e - 2.0 * r4 * (b - oz * dz);
    let a0 = e * e - r4 * (o2 - oz * oz);
    let quartic = [a0, a1, a2, a3, 1.0];
    let half = 0.5 * (t1 - t0) / scale;
    let roots = poly_roots(&quartic, -half, half);
    roots.map(|r| r.map(|t| mid + t * scale))
}

fn eval(p: &[f32], t: f32) -> f32 {
    p.iter().rev().fold(0.0, |acc, &k| acc * t + k)
}

fn derivative(p: &[f32]) -> Vec<f32> {
    p.iter()
        .enumerate()
        .skip(1)
        .map(|(k, &c)| c * k as f32)
        .collect()
}

/// Real roots of `p` (coefficients low to high, degree <= 4) in `[lo, hi]`,
/// ascending.
fn poly_roots(p: &[f32], lo: f32, hi: f32) -> [Option<f32>; 4] {
    let mut out = [None; 4];
    let mut count = 0;
    let dp = derivative(p);
    let mut cuts = vec![lo];
    if dp.len() >= 2 {
        cuts.extend(poly_roots(&dp, lo, hi).into_iter().flatten());
    }
    cuts.push(hi);
    for w in cuts.windows(2) {
        if let Some(t) = bracket_root(p, &dp, w[0], w[1]) {
            if count < 4 {
                out[count] = Some(t);
                count += 1;
            }
        }
    }
    out
}

/// The root of `p` in `[a, b]`, where it is monotonic, if it changes sign.
fn bracket_root(p: &[f32], dp: &[f32], mut a: f32, mut b: f32) -> Option<f32> {
    let (fa, fb) = (eval(p, a), eval(p, b));
    if fa == 0.0 {
        return Some(a);
    }
    if fa * fb > 0.0 {
        return None;
    }
    let rising = fb > fa;
    let mut t = 0.5 * (a + b);
    for _ in 0..24 {
        let f = eval(p, t);
        if (f > 0.0) == rising {
            b = t;
        } else {
            a = t;
        }
        let slope = eval(dp, t);
        let newton = t - f / slope;
        let next = if slope != 0.0 && newton > a && newton < b {
            newton
        } else {
            0.5 * (a + b)
        };
        // Converged (the interval is scaled to about [-1, 1]).
        let step = (next - t).abs();
        t = next;
        if step < 1e-6 {
            break;
        }
    }
    Some(t)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::element::Element;

    /// Crambin (1CRN) from the fixtures: 327 atoms.
    fn crambin() -> (Vec<Vec3>, Vec<f32>) {
        let path =
            std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/small/1CRN.pdb");
        let text = std::fs::read_to_string(path).unwrap();
        text.lines()
            .filter(|l| l.starts_with("ATOM") || l.starts_with("HETATM"))
            .map(|l| {
                let f = |a: usize, b: usize| l[a..b].trim().parse::<f32>().unwrap();
                let element = Element::from_symbol(l[76..78].trim().as_bytes());
                (
                    Vec3::new(f(30, 38), f(38, 46), f(46, 54)),
                    element.vdw_radius(),
                )
            })
            .unzip()
    }

    /// Points on the SAS boundary (sphere points no other SAS sphere
    /// holds), spaced ~0.1 Å: the SES is exactly the points one probe
    /// radius from this set, measured from inside the SAS.
    struct Oracle {
        points: Vec<Vec3>,
        grid: Grid,
        centers: Vec<Vec3>,
        sas: Vec<f32>,
    }

    impl Oracle {
        fn new(positions: &[Vec3], radii: &[f32], probe: f32) -> Oracle {
            let sas: Vec<f32> = radii.iter().map(|r| r + probe).collect();
            let per_sphere = 12_000;
            let golden = std::f32::consts::PI * (3.0 - 5f32.sqrt());
            let unit: Vec<Vec3> = (0..per_sphere)
                .map(|k| {
                    let y = 1.0 - 2.0 * (k as f32 + 0.5) / per_sphere as f32;
                    let s = (1.0 - y * y).sqrt();
                    let a = golden * k as f32;
                    Vec3::new(s * a.cos(), y, s * a.sin())
                })
                .collect();
            let sas_max = sas.iter().copied().fold(0.0, f32::max);
            let all: Vec<u32> = (0..positions.len() as u32).collect();
            let atom_grid = Grid::build(positions, &all, sas_max);
            let points: Vec<Vec3> = (0..positions.len())
                .into_par_iter()
                .flat_map_iter(|i| {
                    let c = positions[i];
                    let r = sas[i];
                    let unit = &unit;
                    let atom_grid = &atom_grid;
                    let sas = &sas;
                    unit.iter().map(move |u| c + *u * r).filter(move |&p| {
                        let mut inside = false;
                        atom_grid.for_each_within(positions, p, sas_max, |j, d2| {
                            inside |= j as usize != i && d2 < sas[j as usize] * sas[j as usize];
                        });
                        !inside
                    })
                })
                .collect();
            let index: Vec<u32> = (0..points.len() as u32).collect();
            let grid = Grid::build(&points, &index, 0.5);
            Oracle {
                points,
                grid,
                centers: positions.to_vec(),
                sas,
            }
        }

        fn in_sas(&self, p: Vec3) -> bool {
            self.centers
                .iter()
                .zip(&self.sas)
                .any(|(c, r)| p.distance_squared(*c) < r * r)
        }

        /// Distance to the SAS boundary, capped at `cap`.
        fn boundary_distance(&self, p: Vec3, cap: f32) -> f32 {
            let mut best = cap * cap;
            self.grid
                .for_each_within(&self.points, p, cap, |_, d2| best = best.min(d2));
            best.sqrt()
        }

        /// Inside the excluded volume by at least `margin`: in the SAS,
        /// and no probe centre (a point outside the SAS) within one probe
        /// radius plus `margin`.
        fn excluded(&self, p: Vec3, probe: f32, margin: f32) -> bool {
            self.in_sas(p) && self.boundary_distance(p, probe + 0.5) >= probe + margin
        }
    }

    fn view_rays(positions: &[Vec3], n: usize) -> Vec<(Vec3, Vec3)> {
        let lo = positions
            .iter()
            .copied()
            .fold(Vec3::splat(f32::MAX), Vec3::min);
        let hi = positions
            .iter()
            .copied()
            .fold(Vec3::splat(f32::MIN), Vec3::max);
        let centre = 0.5 * (lo + hi);
        let half = 0.5 * (hi - lo).max_element() + 4.0;
        let dir = Vec3::new(0.3, -0.2, -1.0).normalize();
        let (e1, e2) = basis(dir);
        (0..n * n)
            .map(|k| {
                let (x, y) = ((k % n) as f32, (k / n) as f32);
                let sx = (x + 0.5) / n as f32 * 2.0 - 1.0;
                let sy = (y + 0.5) / n as f32 * 2.0 - 1.0;
                (centre - dir * 60.0 + (e1 * sx + e2 * sy) * half, dir)
            })
            .collect()
    }

    #[test]
    fn two_atoms_make_two_caps_and_one_full_torus() {
        let positions = [Vec3::ZERO, Vec3::new(3.0, 0.0, 0.0)];
        let radii = [1.7, 1.7];
        let ses = build(&positions, &radii, 1.4);
        assert_eq!(ses.convex.len(), 2);
        assert_eq!(ses.caps, [1, 0]);
        assert!(ses.probes.is_empty());
        assert_eq!(
            ses.tori,
            [Torus {
                atoms: [0, 1],
                probes: [FULL_CIRCLE; 2]
            }]
        );
        // Down the bond's perpendicular bisector, the ray meets the
        // torus's saddle: one probe radius from the probe circle.
        let hit = ses
            .cast(
                &positions,
                &radii,
                Vec3::new(1.5, 10.0, 0.0),
                Vec3::NEG_Y,
                0.0,
            )
            .unwrap();
        let c = circle(positions[0], 3.1, positions[1], 3.1).unwrap();
        assert!(matches!(hit.patch, PatchRef::Torus(0)));
        assert!((hit.t - (10.0 - (c.radius - 1.4))).abs() < 1e-3, "{hit:?}");
        assert!((hit.normal - Vec3::Y).length() < 1e-3);
    }

    #[test]
    fn torus_roots_are_on_the_torus_and_ascending() {
        let (c, n) = (
            Vec3::new(100.0, -50.0, 20.0),
            Vec3::new(1.0, 2.0, 2.0).normalize(),
        );
        let (major, minor) = (1.2, 1.4); // a spindle torus
        for k in 0..200 {
            let a = k as f32 * 0.37;
            let o = c + Vec3::new(a.cos() * 9.0, (a * 1.3).sin() * 9.0, (a * 0.7).cos() * 9.0);
            let target = c + Vec3::new((a * 2.1).sin(), (a * 1.7).cos(), a.sin()) * 1.5;
            let d = (target - o).normalize();
            let roots = sphere_roots(o, d, c, major + minor)
                .map_or([None; 4], |w| torus_roots(o, d, c, n, major, minor, w));
            let mut last = f32::MIN;
            for t in roots.into_iter().flatten() {
                let v = o + d * t - c;
                let z = v.dot(n);
                let s = (v - n * z).length();
                // Spindle tori also have roots on the inner lemon, a tube
                // radius from the circle's far side.
                let near = Vec2::new(s - major, z).length() - minor;
                let far = Vec2::new(s + major, z).length() - minor;
                let off = near.abs().min(far.abs());
                assert!(off < 2e-3, "ray {k}: off the torus by {off}");
                assert!(t >= last);
                last = t;
            }
        }
    }

    #[test]
    fn probe_distance_is_the_probe_radius_on_the_surface() {
        let (positions, radii) = crambin();
        let probe = 1.4;
        let ses = build(&positions, &radii, probe);
        let rays = view_rays(&positions, 30);
        let hits: Vec<(Vec3, Vec3)> = rays
            .par_iter()
            .filter_map(|&(o, d)| {
                ses.cast(&positions, &radii, o, d, 0.0)
                    .map(|h| (o + d * h.t, d))
            })
            .collect();
        assert!(hits.len() > 300, "only {} hits", hits.len());
        for &(p, d) in &hits {
            let on = ses.probe_distance(&positions, &radii, p);
            assert!((on - probe).abs() < 2e-3, "{on} at a hit {p}");
            assert!(ses.probe_distance(&positions, &radii, p - d * 0.05) < probe);
            assert!(ses.probe_distance(&positions, &radii, p + d * 0.05) > probe);
        }
    }

    #[test]
    fn every_hit_on_crambin_is_one_probe_radius_from_the_sas() {
        let (positions, radii) = crambin();
        let probe = 1.4;
        let ses = build(&positions, &radii, probe);
        let oracle = Oracle::new(&positions, &radii, probe);
        let rays = view_rays(&positions, 40);
        let hits: Vec<Hit> = rays
            .par_iter()
            .filter_map(|&(o, d)| ses.cast(&positions, &radii, o, d, 0.0).map(|h| (o, d, h)))
            .map(|(o, d, h)| {
                let p = o + d * h.t;
                let dist = oracle.boundary_distance(p, probe + 0.5);
                assert!(
                    (dist - probe).abs() < 0.08,
                    "{:?} at {p}: {dist} from the SAS",
                    h.patch
                );
                h
            })
            .collect();
        assert!(hits.len() > 500, "only {} hits", hits.len());
        for kind in 0..3 {
            assert!(
                hits.iter().any(|h| match h.patch {
                    PatchRef::Convex(_) => kind == 0,
                    PatchRef::Torus(_) => kind == 1,
                    PatchRef::Concave(_) => kind == 2,
                }),
                "no hit of patch kind {kind}"
            );
        }
    }

    /// Casts against the patches and marches the oracle's excluded volume
    /// along the same rays; returns every disagreement over 0.1 Å
    /// (infinite when one side misses) and the ray count. A ray counts as
    /// hitting when marching finds the volume at all, and as missing only
    /// when it finds no point clear of the sampling error -- so slivers
    /// where probe spheres all but meet go either way.
    fn disagreements(probe: f32) -> (Vec<f32>, usize) {
        let (positions, radii) = crambin();
        let ses = build(&positions, &radii, probe);
        let oracle = Oracle::new(&positions, &radii, probe);
        let rays = view_rays(&positions, 32);
        let march = |o: Vec3, d: Vec3, margin: f32| {
            let mut t = positions
                .iter()
                .zip(&radii)
                .filter_map(|(&c, &r)| sphere_roots(o, d, c, r + probe))
                .map(|[t, _]| t)
                .fold(f32::MAX, f32::min);
            while t < 120.0 {
                if oracle.excluded(o + d * t, probe, margin) {
                    return Some(t);
                }
                t += 0.02;
            }
            None
        };
        let off = rays
            .par_iter()
            .filter_map(|&(o, d)| {
                let cast = ses.cast(&positions, &radii, o, d, 0.0).map(|h| h.t);
                match (cast, march(o, d, 0.0)) {
                    (Some(a), Some(b)) => {
                        let clear = march(o, d, 0.02).map_or(f32::INFINITY, |c| (a - c).abs());
                        let e = (a - b).abs().min(clear);
                        (e >= 0.1).then_some(e)
                    }
                    (None, Some(_)) => march(o, d, 0.02).map(|_| f32::INFINITY),
                    (Some(_), None) => Some(f32::INFINITY),
                    (None, None) => None,
                }
            })
            .collect();
        (off, rays.len())
    }

    #[test]
    fn crambin_hits_agree_with_marching_the_excluded_volume() {
        // A 3 Å probe makes spindle tori and overlapping probes common.
        // Near-misses come from the oracle's ~0.15 Å sampling on grazing
        // rays; a hole or a stray patch would be off by more than 1 Å.
        for probe in [1.4, 3.0] {
            let (off, total) = disagreements(probe);
            eprintln!("probe {probe}: {off:?} of {total} rays");
            assert!(off.iter().all(|&e| e < 1.0), "probe {probe}: {off:?}");
            assert!(off.len() * 100 <= total, "probe {probe}: {off:?}");
        }
    }
}
