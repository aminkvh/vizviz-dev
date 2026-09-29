//! Skin surface patches (Edelsbrunner, "Deformable Smooth Surface
//! Design," Discrete & Computational Geometry 21:87-115, 1999) -- the
//! mixed complex over the local power cells of [`crate::feasible_cell`],
//! one quadric patch per Delaunay simplex, plus the exact test for which
//! patch owns which point.
//!
//! **Source.** Checked against Cheng, Dey, Edelsbrunner & Sullivan,
//! "Dynamic Skin Triangulation," DCG 25:525-568, 2001, Section 2 (which
//! restates the 1999 construction), not re-derived from memory -- this
//! module's earlier life re-derived its edge patch twice and got it wrong
//! both times. Their definitions, in this module's notation:
//!
//! - A weighted point is a sphere `S_i = (z_i, w_i)` (`w = radius^2`,
//!   possibly negative), the zero set of `f_i(x) = |x - z_i|^2 - w_i`.
//!   Affine combinations `sum(gamma_i * f_i)`, `sum(gamma_i) = 1`, are
//!   again such functions, so they define spheres.
//! - Shrinking by `s` maps `(z, w)` to `(z, s*w)`; the paper's fixed
//!   `s = 1/2` is [`DEFAULT_SHRINK`]. The skin is the envelope of the
//!   shrunk spheres of the *convex* hull of the `S_i`.
//! - The **mixed cell** of a Delaunay simplex `X` (vertex, edge,
//!   triangle, tetrahedron of the regular triangulation) is the
//!   Minkowski sum `mu_X = (1-s)*delta_X (+) s*nu_X` of the simplex and
//!   its dual power-diagram cell `nu_X` (a polyhedron, polygon, edge,
//!   vertex respectively). The mixed cells tile space.
//! - **Inside `mu_X` the skin equals the envelope of the shrunk spheres
//!   of the *affine* hull of `X`'s spheres** ("Skin patches" in their
//!   Section 2): a sphere for a vertex or tetrahedron, a hyperboloid of
//!   revolution for an edge or triangle, centered at the point `z_X`
//!   common to the affine hulls of `delta_X` and `nu_X` -- the
//!   orthocenter of `X` (its orthosphere `(z_X, w_z)` has
//!   `|z_X - z_i|^2 - w_i = w_z` for every `i` in `X`).
//!
//! **The one formula.** Minimizing the shrunk affine family's distance
//! function at `x` over the family (its center `c` ranges over
//! `aff delta_X`; every member has `f_c(z_X) = w_z`, so `f_c(x) = |x -
//! c|^2 - (|z_X - c|^2 - w_z)`) is a convex quadratic in `c` with
//! minimizer `c = z_X + u_par / (1 - s)`, `u = x - z_X`, `u_par` its
//! component along `aff delta_X`. Substituting:
//!
//! ```text
//! F_X(x) = |u_perp|^2 - K |u_par|^2 + s w_z,    K = s / (1 - s)
//! ```
//!
//! negative inside the body, zero on the patch. `u_perp` lies in the
//! dual cell's directions. Vertex: `|x - a|^2 = s w_a` (the shrunk
//! sphere). Edge: `rho^2 - K t^2 = -s w_z`, the hyperboloid this
//! module's edge patch already had (waist on the power bisector,
//! one-sheeted when `w_z < 0`, i.e. overlapping atoms; two-sheeted caps
//! otherwise). Triangle: `d^2 - K rho^2 = -s w_z` with `d` the offset
//! along the triangle normal -- a two-sheeted saddle blend when `w_z <
//! 0`, a one-sheeted tunnel through the triangle otherwise. Tetrahedron:
//! `|x - z|^2 = (1 - s) w_z`, a spherical void, present only when `w_z >
//! 0`. At `s = 1/2` every hyperboloid's asymptotic cone has a right
//! opening angle, as the paper says it must. The minimizer `c` is also
//! why the patch only counts inside `mu_X`: there, and only there, the
//! affine family's minimizer is a *convex* combination.
//!
//! **Membership, exactly.** `x` is in `mu_X` iff `x = (1-s) p + s q` with
//! `p` in `delta_X` and `q` in `nu_X`; since the two live in orthogonal
//! affine subspaces through `z_X`, the decomposition is unique:
//! `p = z_X + u_par / (1 - s)` ([`Patch::delaunay_point`], tested by barycentric
//! coordinates) and `q = z_X + u_perp / s` ([`Patch::voronoi_point`],
//! tested by "no atom outside `X` has a smaller power distance to `q`
//! than `X`'s own atoms" -- an ordinary nearest-power-neighbor query,
//! [`voronoi_owned_by`]). This replaces the earlier heuristic (an axial
//! window around the waist plus a nearest-atom check at `x` itself),
//! which drew vertex spheres into regions that belong to edge patches
//! and had no triangle patches at all -- the visible ridges where three
//! necks met.
//!
//! **What is still local.** The simplices come from
//! [`crate::feasible_cell`]'s per-atom power cells (edges = surviving
//! faces, triangles = two faces sharing an edge, tetrahedra = three faces
//! sharing a vertex), so they are the regular triangulation only as far
//! as each atom's window ([`CELL_BOUND_RADIUS`]) reaches; a simplex the
//! window misses leaves its mixed cell undrawn. The membership tests
//! themselves are exact.

use glam::Vec3;
use rayon::prelude::*;

use crate::feasible_cell;
use crate::spatial::Grid;

/// The paper's own construction (`sqrt(S) = (z, r / sqrt 2)`), i.e.
/// weights halved; CGAL's `Skin_surface_3` exposes the same generalized
/// `s in (0, 1)` this module takes.
pub const DEFAULT_SHRINK: f32 = 0.5;

/// Radius of the neighbor search that builds each atom's power cell:
/// generous against bonded/non-bonded spacing (1-4 Å) so a cell's real
/// faces are all found.
pub const NEIGHBOR_SEARCH_RADIUS: f32 = 6.0;

/// Half-extent of each atom's feasible-cell window
/// (`feasible_cell::build`'s `bound_radius`) -- several times the real
/// neighbor spacing, per that module's own measured guidance.
pub const CELL_BOUND_RADIUS: f32 = 5.0;

/// Cell edge of the grid the ownership queries ([`voronoi_owned_by`])
/// scan. Each query derives its own exact radius (`sqrt(pi_X(q) +
/// w_max)`, see there), typically 3-4 Å for van-der-Waals-scale weights,
/// so a 4 Å cell answers most queries from the 3x3x3 block and larger
/// radii just widen the block ([`Grid::for_each_within`]).
pub const OWNERSHIP_CELL: f32 = 4.0;

/// Tolerance for the membership tests, in Å^2 of power distance (and,
/// for the barycentric half, as a fraction of the simplex): patches of
/// adjacent mixed cells meet tangent-continuously, so accepting a little
/// past a cell boundary draws the same surface twice to within a
/// pixel, while rejecting a little short opens a pinhole where each
/// side's quadric root lands a hair inside the other's cell.
const MEMBERSHIP_EPS: f32 = 2e-2;

/// Slack for the build-time visibility test ([`build_complex`]): a patch
/// is dropped only when its quadric provably stays clear of its mixed
/// cell by more than this.
const VISIBILITY_EPS: f32 = 5e-2;

/// Chavent, Levy & Maigret 2008 (J. Mol. Graph. Model. 27:209-216) eq. 1:
/// weight `r^2 / shrink`, so that `shrink * w == r^2` and an isolated
/// atom's own patch sits on its true van der Waals sphere at every
/// shrink factor, the same size every other representation draws it.
pub fn weight_for_radius(radius: f32, shrink: f32) -> f32 {
    radius * radius / shrink
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(u8)]
pub enum PatchKind {
    Vertex = 0,
    Edge = 1,
    Triangle = 2,
    Tetrahedron = 3,
}

impl PatchKind {
    pub fn simplex_size(self) -> usize {
        self as usize + 1
    }
}

/// One mixed cell's quadric: everything a renderer needs to ray-cast
/// it and decide which hits are really in the cell.
#[derive(Clone, Copy, Debug)]
pub struct Patch {
    pub kind: PatchKind,
    /// The simplex's atoms; unused slots are `u32::MAX`.
    pub atoms: [u32; 4],
    /// The orthocenter `z_X`.
    pub center: Vec3,
    /// Unit edge direction (edge) or unit triangle normal (triangle);
    /// zero for a vertex or tetrahedron, whose quadric is isotropic.
    pub axis: Vec3,
    /// The orthosphere weight `w_z`.
    pub weight: f32,
    /// A sphere containing the whole patch, or `bound_radius <= 0` when
    /// the quadric never enters the cell (a tetrahedron with no void, a
    /// triangle whose tunnel/saddle lies outside its window).
    pub bound_center: Vec3,
    pub bound_radius: f32,
}

/// The `(alpha, beta)` of `F = alpha |u|^2 + beta (u . axis)^2 + s w_z`,
/// the module doc's formula regrouped so one quadratic serves all four
/// kinds (`|u_perp|^2 = |u|^2 - d^2` for an edge, `d^2` for a triangle).
fn coefficients(kind: PatchKind, shrink: f32) -> (f32, f32) {
    let k = shrink / (1.0 - shrink);
    match kind {
        PatchKind::Vertex => (1.0, 0.0),
        PatchKind::Edge => (1.0, -(1.0 + k)),
        PatchKind::Triangle => (-k, 1.0 + k),
        PatchKind::Tetrahedron => (-k, 0.0),
    }
}

impl Patch {
    /// The module doc's `F_X(x)`: negative inside the body, zero on the
    /// patch, positive outside -- meaningful only inside the mixed cell
    /// ([`Self::contains`]).
    pub fn field(&self, x: Vec3, shrink: f32) -> f32 {
        let (alpha, beta) = coefficients(self.kind, shrink);
        let u = x - self.center;
        let d = u.dot(self.axis);
        alpha * u.length_squared() + beta * d * d + shrink * self.weight
    }

    /// The outward unit normal at a point on the patch: `F`'s gradient,
    /// which increases out of the body for every kind (for a void it
    /// points into the cavity, which *is* outward for the body).
    pub fn normal_at(&self, x: Vec3, shrink: f32) -> Vec3 {
        let (alpha, beta) = coefficients(self.kind, shrink);
        let u = x - self.center;
        let d = u.dot(self.axis);
        (alpha * u + beta * d * self.axis).normalize_or_zero()
    }

    /// Both ray parameters where `origin + t*dir` (`dir` unit length)
    /// crosses the quadric, ascending, `None` where there is no such
    /// positive root. A renderer must check *both* against
    /// [`Self::contains`]: the nearer crossing routinely lies outside the
    /// cell while the farther one is the real hit (looking down the axis
    /// of a neck, for instance).
    pub fn ray_roots(&self, origin: Vec3, dir: Vec3, shrink: f32) -> [Option<f32>; 2] {
        let (alpha, beta) = coefficients(self.kind, shrink);
        // Solve from the ray's closest approach to the center, not from a
        // camera tens of Å away: keeps every coefficient at the patch's
        // own scale, so the quadratic's cancellation-prone case (`a -> 0`,
        // a ray along a hyperboloid's asymptotic cone -- 45 degrees off
        // its axis at shrink 1/2, hardly rare) stays accurate in f32.
        let t0 = (self.center - origin).dot(dir);
        let u0 = origin + dir * t0 - self.center;
        let d0 = u0.dot(self.axis);
        let dd = dir.dot(self.axis);
        let a = alpha + beta * dd * dd;
        let b = 2.0 * (alpha * u0.dot(dir) + beta * d0 * dd);
        let c = alpha * u0.length_squared() + beta * d0 * d0 + shrink * self.weight;
        let [r1, r2] = quadratic_roots(a, b, c);
        let mut out = [None, None];
        let mut n = 0;
        for r in [r1, r2].into_iter().flatten() {
            let t = t0 + r;
            if t > RAY_EPS {
                out[n] = Some(t);
                n += 1;
            }
        }
        if let [Some(x), Some(y)] = out {
            if y < x {
                out = [Some(y), Some(x)];
            }
        }
        out
    }

    /// The Delaunay component of `x`'s decomposition (`p` in the module
    /// doc): the point of `aff delta_X` whose scaled offset reproduces
    /// `x`'s component along the simplex.
    pub fn delaunay_point(&self, x: Vec3, shrink: f32) -> Vec3 {
        let u = x - self.center;
        let u_par = match self.kind {
            PatchKind::Vertex => Vec3::ZERO,
            PatchKind::Edge => self.axis * u.dot(self.axis),
            PatchKind::Triangle => u - self.axis * u.dot(self.axis),
            PatchKind::Tetrahedron => u,
        };
        self.center + u_par / (1.0 - shrink)
    }

    /// The power-diagram component of `x`'s decomposition (`q` in the
    /// module doc): the point of `aff nu_X` whose scaled offset
    /// reproduces `x`'s component across the simplex.
    pub fn voronoi_point(&self, x: Vec3, shrink: f32) -> Vec3 {
        let u = x - self.center;
        let u_perp = match self.kind {
            PatchKind::Vertex => u,
            PatchKind::Edge => u - self.axis * u.dot(self.axis),
            PatchKind::Triangle => self.axis * u.dot(self.axis),
            PatchKind::Tetrahedron => Vec3::ZERO,
        };
        self.center + u_perp / shrink
    }

    /// Is `p` (already in `aff delta_X`) inside the simplex? Barycentric
    /// coordinates, with [`MEMBERSHIP_EPS`] of slack.
    pub fn delaunay_contains(&self, p: Vec3, positions: &[Vec3]) -> bool {
        let a = positions[self.atoms[0] as usize];
        let eps = MEMBERSHIP_EPS;
        match self.kind {
            PatchKind::Vertex => true,
            PatchKind::Edge => {
                let b = positions[self.atoms[1] as usize];
                let e = b - a;
                let t = (p - a).dot(e) / e.length_squared();
                (-eps..=1.0 + eps).contains(&t)
            }
            PatchKind::Triangle => {
                let b = positions[self.atoms[1] as usize];
                let c = positions[self.atoms[2] as usize];
                let (v0, v1, v2) = (b - a, c - a, p - a);
                let (d00, d01, d11) = (v0.dot(v0), v0.dot(v1), v1.dot(v1));
                let (d20, d21) = (v2.dot(v0), v2.dot(v1));
                let denom = d00 * d11 - d01 * d01;
                let v = (d11 * d20 - d01 * d21) / denom;
                let w = (d00 * d21 - d01 * d20) / denom;
                let u = 1.0 - v - w;
                u >= -eps && v >= -eps && w >= -eps
            }
            PatchKind::Tetrahedron => {
                let b = positions[self.atoms[1] as usize];
                let c = positions[self.atoms[2] as usize];
                let d = positions[self.atoms[3] as usize];
                let (e1, e2, e3, r) = (b - a, c - a, d - a, p - a);
                let det = e1.dot(e2.cross(e3));
                let l1 = r.dot(e2.cross(e3)) / det;
                let l2 = e1.dot(r.cross(e3)) / det;
                let l3 = e1.dot(e2.cross(r)) / det;
                let l0 = 1.0 - l1 - l2 - l3;
                l0 >= -eps && l1 >= -eps && l2 >= -eps && l3 >= -eps
            }
        }
    }

    /// Is `x` in this patch's mixed cell? `voronoi_owned` answers
    /// "does `X` own this point of the power diagram" for the `q` this
    /// computes -- [`voronoi_owned_by`] over a [`Grid`] for real
    /// renderers, or a brute-force scan in tests.
    pub fn contains(
        &self,
        x: Vec3,
        shrink: f32,
        positions: &[Vec3],
        voronoi_owned: impl FnOnce(Vec3) -> bool,
    ) -> bool {
        self.delaunay_contains(self.delaunay_point(x, shrink), positions)
            && voronoi_owned(self.voronoi_point(x, shrink))
    }

    /// Which of the simplex's atoms `x` is nearest to in power distance
    /// -- the atom whose color a hit takes.
    pub fn nearest_atom(&self, x: Vec3, positions: &[Vec3], weights: &[f32]) -> u32 {
        self.members()
            .iter()
            .copied()
            .min_by(|&i, &j| {
                power_distance(x, positions[i as usize], weights[i as usize]).total_cmp(
                    &power_distance(x, positions[j as usize], weights[j as usize]),
                )
            })
            .expect("a patch has at least one atom")
    }

    pub fn members(&self) -> &[u32] {
        &self.atoms[..self.kind.simplex_size()]
    }

    pub fn has_surface(&self) -> bool {
        self.bound_radius > 0.0
    }
}

pub fn power_distance(x: Vec3, center: Vec3, weight: f32) -> f32 {
    (x - center).length_squared() - weight
}

/// `X`'s power-diagram cell contains `q` iff no atom outside `X` is
/// nearer to `q` in power distance than `X`'s own atoms are (all of
/// which tie at any point of `aff nu_X`). An atom `m` can only win if
/// `|q - z_m|^2 < w_m + pi_X(q)`, so the scan radius is exactly
/// `sqrt(pi_X(q) + weight_max)` (`weight_max` the largest weight among
/// all atoms) -- no fixed neighborhood to get wrong, and nothing at all
/// to scan when that is negative.
pub fn voronoi_owned_by(
    grid: &Grid,
    positions: &[Vec3],
    weights: &[f32],
    weight_max: f32,
    q: Vec3,
    members: &[u32],
) -> bool {
    let first = members[0] as usize;
    let pi_x = power_distance(q, positions[first], weights[first]);
    let reach2 = pi_x + weight_max;
    if reach2 <= 0.0 {
        return true;
    }
    let bound = pi_x - MEMBERSHIP_EPS;
    let mut owned = true;
    grid.for_each_within(positions, q, reach2.sqrt(), |m, _| {
        if owned
            && !members.contains(&m)
            && power_distance(q, positions[m as usize], weights[m as usize]) < bound
        {
            owned = false;
        }
    });
    owned
}

/// A hit must be this far in front of the ray origin to count.
const RAY_EPS: f32 = 1e-4;

/// Real roots of `a t^2 + b t + c = 0` in no particular order, by the
/// cancellation-free form (`q = -(b + sign(b) sqrt(disc)) / 2`, roots
/// `q / a` and `c / q`): the textbook `(-b +- sqrt) / 2a` loses the
/// small root entirely once `a` is small against `b`, which for these
/// quadrics is every ray near a hyperboloid's asymptotic cone. `a` and
/// `q` near zero fall back to the linear/other root.
fn quadratic_roots(a: f32, b: f32, c: f32) -> [Option<f32>; 2] {
    if a.abs() < 1e-12 {
        if b.abs() < 1e-12 {
            return [None, None];
        }
        return [Some(-c / b), None];
    }
    let disc = b * b - 4.0 * a * c;
    if disc < 0.0 {
        return [None, None];
    }
    let q = -0.5 * (b + b.signum() * disc.sqrt());
    if q.abs() < 1e-20 {
        // b == 0 and c == 0 (or disc == 0 with b == 0): a double root at 0.
        return [Some(0.0), None];
    }
    [Some(q / a), Some(c / q)]
}

/// Every patch of the local mixed complex over `positions`/`weights`.
pub struct SkinComplex {
    /// One patch per simplex; those whose quadric never enters their
    /// own mixed cell (buried atoms and the simplices between them --
    /// most of a protein's interior) have `bound_radius == 0` and
    /// [`Patch::has_surface`] false, so renderers skip them while the
    /// membership tests still see every cell.
    pub patches: Vec<Patch>,
    pub shrink: f32,
    /// The grid the ownership queries run against (cell edge
    /// [`OWNERSHIP_CELL`]), shared so renderers need not build a second
    /// one over the same atoms.
    pub grid: Grid,
    /// The largest weight of any atom, the ownership queries' reach.
    pub weight_max: f32,
    /// Per patch (same order as `patches`, CSR: patch `i`'s list is
    /// `competitors[competitor_starts[i]..competitor_starts[i + 1]]`),
    /// the only atoms that can ever take a point away from the patch's
    /// simplex -- see [`SkinComplex::owned_by_competitors`]. Empty for
    /// patches without a surface.
    pub competitor_starts: Vec<u32>,
    pub competitors: Vec<u32>,
}

impl SkinComplex {
    /// Patch counts by kind (vertices, edges, triangles, tetrahedra),
    /// counting only patches with a surface.
    pub fn counts(&self) -> [usize; 4] {
        let mut counts = [0; 4];
        for p in self.patches.iter().filter(|p| p.has_surface()) {
            counts[p.kind as usize] += 1;
        }
        counts
    }

    /// [`voronoi_owned_by`] over this complex's own grid.
    pub fn voronoi_owned_by(
        &self,
        positions: &[Vec3],
        weights: &[f32],
        q: Vec3,
        members: &[u32],
    ) -> bool {
        voronoi_owned_by(&self.grid, positions, weights, self.weight_max, q, members)
    }

    /// The atoms that can undercut patch `index`'s simplex anywhere.
    pub fn competitors(&self, index: usize) -> &[u32] {
        &self.competitors
            [self.competitor_starts[index] as usize..self.competitor_starts[index + 1] as usize]
    }

    /// The Voronoi half of the membership test for patch `index`, over
    /// its competitor list instead of the grid. Exact, by a short
    /// argument: `q` lies outside atom `i`'s power cell only if some atom
    /// `m` has `pi_m(q) < pi_i(q)`; walking the segment from a point of
    /// the cell to `q`, the first bisector crossed belongs to a cell
    /// that shares a face with `i`'s, and since `pi_m - pi_i` is linear
    /// along the segment and already negative there, it is still
    /// negative at `q`. So the facet neighbors of `X`'s atoms -- their
    /// feasible cells' surviving neighbors -- suffice, and a patch's
    /// fragments test a dozen or two atoms instead of a grid block.
    pub fn owned_by_competitors(
        &self,
        index: usize,
        positions: &[Vec3],
        weights: &[f32],
        q: Vec3,
    ) -> bool {
        let patch = &self.patches[index];
        let first = patch.atoms[0] as usize;
        let bound = power_distance(q, positions[first], weights[first]) - MEMBERSHIP_EPS;
        self.competitors(index)
            .iter()
            .all(|&m| power_distance(q, positions[m as usize], weights[m as usize]) >= bound)
    }

    /// [`Patch::contains`] with the competitor-list ownership test.
    pub fn contains(&self, index: usize, x: Vec3, positions: &[Vec3], weights: &[f32]) -> bool {
        self.patches[index].contains(x, self.shrink, positions, |q| {
            self.owned_by_competitors(index, positions, weights, q)
        })
    }
}

/// What one atom's cell says about one simplex it belongs to: the patch,
/// plus the extents of the simplex's dual cell `nu_X` (as seen through
/// this atom's window) that decide whether the patch is visible at all
/// and how tightly it can be bounded, plus the atom's own facet
/// neighbors (its share of the simplex's competitor list).
struct SimplexRecord {
    patch: Patch,
    /// `min |q - z|^2` over `nu_X`.
    nu_min2: f32,
    /// `max |q - z|^2` over `nu_X`, or `INFINITY` when `nu_X` reaches the
    /// window (the true cell is unbounded there).
    nu_max2: f32,
    /// A sphere around the mixed cell `(1-s) delta_X (+) s nu_X` itself,
    /// when `nu_X` is bounded -- usually far tighter than the quadric's
    /// own bound for slivers, whose orthocenters sit far away.
    mixed_bound: Option<(Vec3, f32)>,
    /// The atoms whose cells reported this simplex (always members of
    /// it, so at most four). Its competitors are the union of their
    /// surviving neighbors, looked up once per atom at the end rather
    /// than copied into every record: copying made the merge move ~40
    /// ids per record and dominated build time (8GLV).
    sources: [u32; 4],
    source_count: u8,
}

impl SimplexRecord {
    fn merge(&mut self, other: SimplexRecord) {
        self.nu_min2 = self.nu_min2.min(other.nu_min2);
        self.nu_max2 = self.nu_max2.max(other.nu_max2);
        self.mixed_bound = match (self.mixed_bound, other.mixed_bound) {
            (Some(a), Some(b)) => Some(if a.1 >= b.1 { a } else { b }),
            _ => None,
        };
        for &src in &other.sources[..other.source_count as usize] {
            if !self.sources[..self.source_count as usize].contains(&src) {
                self.sources[self.source_count as usize] = src;
                self.source_count += 1;
            }
        }
    }

    fn from_source(
        patch: Patch,
        nu_min2: f32,
        nu_max2: f32,
        mixed_bound: Option<(Vec3, f32)>,
        source: u32,
    ) -> Self {
        SimplexRecord {
            patch,
            nu_min2,
            nu_max2,
            mixed_bound,
            sources: [source, u32::MAX, u32::MAX, u32::MAX],
            source_count: 1,
        }
    }
}

/// Builds the local mixed complex: each atom's power cell
/// ([`feasible_cell::build`], in parallel), the simplices read off the
/// cells' combinatorics, and one [`Patch`] per simplex. `weights` are
/// the unshrunk `r^2 / shrink` values ([`weight_for_radius`]).
///
/// **Visibility.** Writing a point of the cell as `x = (1-s) p + s q`
/// (`p` in the simplex, `q` in its dual cell), the module doc's field is
/// `F(x) = s * G(p, q)` with `G = s |q - z|^2 - (1-s) |p - z|^2 + w_z`,
/// so the quadric enters the cell iff `G` takes both signs over
/// `delta_X x nu_X`: `min G = s * min|q-z|^2 - (1-s) * max|p-z|^2 + w_z`
/// and `max G = s * max|q-z|^2 - (1-s) * min|p-z|^2 + w_z` (the extremes
/// over `p` are at the simplex's own atoms, `|z_i - z|^2 = w_i + w_z`,
/// and at the point of the simplex nearest `z`; over `q`, at the dual
/// cell's polytope vertices and its point nearest `z`). A patch failing
/// this can never be hit, so it gets no bounding sphere and no
/// billboard -- which is what makes the drawn patch count grow with a
/// structure's surface area rather than its volume: a buried atom's
/// mixed cell is a fraction of an Ångström across and its sphere never
/// reaches it.
pub fn build_complex(positions: &[Vec3], weights: &[f32], shrink: f32) -> SkinComplex {
    assert_eq!(positions.len(), weights.len());
    assert!(shrink > 0.0 && shrink < 1.0, "shrink must be in (0, 1)");
    let indices: Vec<u32> = (0..positions.len() as u32).collect();
    let grid = Grid::build(positions, &indices, OWNERSHIP_CELL);
    let weight_max = weights.iter().copied().fold(f32::NEG_INFINITY, f32::max);

    let per_atom: Vec<AtomRecords> = (0..positions.len())
        .into_par_iter()
        .map(|i| atom_records(i, positions, weights, shrink, &grid))
        .collect();
    let (per_atom_records, neighbors): (Vec<_>, Vec<_>) = per_atom.into_iter().unzip();
    let mut flat: Vec<([u32; 4], SimplexRecord)> = per_atom_records.into_iter().flatten().collect();
    flat.par_sort_unstable_by_key(|(key, _)| *key);

    // Each simplex arrives once per member atom whose cell saw it; merge
    // the runs of equal keys.
    let mut merged: Vec<SimplexRecord> = Vec::with_capacity(flat.len() / 2);
    let mut last_key = None;
    for (key, record) in flat {
        if last_key == Some(key) {
            merged
                .last_mut()
                .expect("run has a first record")
                .merge(record);
        } else {
            merged.push(record);
            last_key = Some(key);
        }
    }

    let mut finished: Vec<(Patch, Vec<u32>)> = merged
        .into_par_iter()
        .filter(|rec| {
            // A tetrahedron's dual cell is the single point `z`; if some
            // fifth atom undercuts the four there, the tetrahedron is not
            // in the regular triangulation at all (a degenerate shared
            // vertex of the local cells).
            rec.patch.kind != PatchKind::Tetrahedron
                || voronoi_owned_by(
                    &grid,
                    positions,
                    weights,
                    weight_max,
                    rec.patch.center,
                    &rec.patch.atoms,
                )
        })
        .map(|rec| finish_patch(rec, &neighbors, positions, weights, shrink))
        .collect();
    // Deterministic order (the hash map's is not), vertices first.
    finished.sort_by_key(|(p, _)| (p.kind as u8, p.atoms));
    let mut patches = Vec::with_capacity(finished.len());
    let mut competitor_starts = Vec::with_capacity(finished.len() + 1);
    let mut competitors = Vec::new();
    competitor_starts.push(0u32);
    for (patch, list) in finished {
        competitors.extend(list);
        competitor_starts.push(competitors.len() as u32);
        patches.push(patch);
    }
    SkinComplex {
        patches,
        shrink,
        grid,
        weight_max,
        competitor_starts,
        competitors,
    }
}

/// Atom `i`'s cell and every simplex record it contributes.
/// One atom's simplex records, keyed by sorted member atoms, and its
/// cell's surviving neighbors.
type AtomRecords = (Vec<([u32; 4], SimplexRecord)>, Vec<u32>);

/// Also returns the cell's surviving neighbors: every record's
/// competitors (see [`SimplexRecord::sources`]).
fn atom_records(
    i: usize,
    positions: &[Vec3],
    weights: &[f32],
    shrink: f32,
    grid: &Grid,
) -> AtomRecords {
    let me = i as u32;
    let origin = positions[i];
    let mut neighbors = Vec::new();
    grid.for_each_within(positions, origin, NEIGHBOR_SEARCH_RADIUS, |j, _| {
        if j as usize != i {
            neighbors.push((j, positions[j as usize] - origin, weights[j as usize]));
        }
    });
    let Some(cell) = feasible_cell::build(weights[i], neighbors, CELL_BOUND_RADIUS) else {
        // Dominated: an empty power cell, so nothing can ever be drawn
        // for it, but it keeps a (surfaceless) patch so `patches` has
        // one vertex entry per atom.
        let mut patch = vertex_patch(me, positions, weights, shrink);
        patch.bound_radius = 0.0;
        return (
            vec![(
                [me, u32::MAX, u32::MAX, u32::MAX],
                SimplexRecord::from_source(patch, 0.0, 0.0, None, me),
            )],
            Vec::new(),
        );
    };
    let world = |v: u32| origin + cell.vertices[v as usize];
    let touches = |v: u32| cell.vertex_touches_window(v);
    let s = shrink;
    let my_neighbors: Vec<u32> = cell.surviving_neighbors().collect();
    let mut out = Vec::new();

    // Vertex: nu is the whole cell.
    {
        let patch = vertex_patch(me, positions, weights, s);
        let verts: Vec<u32> = {
            let mut all: Vec<u32> = cell
                .faces
                .iter()
                .flat_map(|f| f.verts.iter().copied())
                .collect();
            all.sort_unstable();
            all.dedup();
            all
        };
        let unbounded = verts.iter().any(|&v| touches(v));
        let nu_max2 = if unbounded {
            f32::INFINITY
        } else {
            verts
                .iter()
                .map(|&v| (world(v) - patch.center).length_squared())
                .fold(0.0, f32::max)
        };
        let mixed_bound = (!unbounded)
            .then(|| bounding_sphere_of(verts.iter().map(|&v| origin * (1.0 - s) + world(v) * s)));
        out.push((
            [me, u32::MAX, u32::MAX, u32::MAX],
            SimplexRecord {
                patch,
                nu_min2: 0.0,
                nu_max2,
                mixed_bound,
                sources: [me, u32::MAX, u32::MAX, u32::MAX],
                source_count: 1,
            },
        ));
    }

    // Edges: nu is the face polygon.
    for face in cell.faces.iter().filter(|f| f.neighbor.is_some()) {
        let j = face.neighbor.unwrap();
        let key = sorted2([me, j]);
        let Some(patch) = edge_patch(key, positions, weights, s) else {
            continue;
        };
        let verts_w: Vec<Vec3> = face.verts.iter().map(|&v| world(v)).collect();
        let unbounded = face.verts.iter().any(|&v| touches(v));
        let nu_max2 = if unbounded {
            f32::INFINITY
        } else {
            verts_w
                .iter()
                .map(|&v| (v - patch.center).length_squared())
                .fold(0.0, f32::max)
        };
        let nu_min2 = distance2_to_convex_polygon(patch.center, &verts_w, face.normal);
        let mixed_bound = (!unbounded).then(|| {
            let ends = [origin, positions[j as usize]];
            bounding_sphere_of(
                ends.iter()
                    .flat_map(|&e| verts_w.iter().map(move |&v| e * (1.0 - s) + v * s)),
            )
        });
        out.push((
            [key[0], key[1], u32::MAX, u32::MAX],
            SimplexRecord {
                patch,
                nu_min2,
                nu_max2,
                mixed_bound,
                sources: [me, u32::MAX, u32::MAX, u32::MAX],
                source_count: 1,
            },
        ));
    }

    // Triangles: nu is the polygon edge two faces share.
    let real: Vec<&feasible_cell::Face> =
        cell.faces.iter().filter(|f| f.neighbor.is_some()).collect();
    for (x, fj) in real.iter().enumerate() {
        for fk in &real[x + 1..] {
            let shared: Vec<u32> = fj
                .verts
                .iter()
                .copied()
                .filter(|v| fk.verts.contains(v))
                .collect();
            if shared.len() < 2 {
                continue;
            }
            let (j, k) = (fj.neighbor.unwrap(), fk.neighbor.unwrap());
            let key = sorted3([me, j, k]);
            let Some(patch) = triangle_patch(key, positions, weights, s) else {
                continue;
            };
            let verts_w: Vec<Vec3> = shared.iter().map(|&v| world(v)).collect();
            let unbounded = shared.iter().any(|&v| touches(v));
            let nu_max2 = if unbounded {
                f32::INFINITY
            } else {
                verts_w
                    .iter()
                    .map(|&v| (v - patch.center).length_squared())
                    .fold(0.0, f32::max)
            };
            let mut nu_min2 = f32::INFINITY;
            for a in 0..verts_w.len() {
                for b in a + 1..verts_w.len() {
                    nu_min2 =
                        nu_min2.min(distance2_to_segment(patch.center, verts_w[a], verts_w[b]));
                }
            }
            let mixed_bound = (!unbounded).then(|| {
                let corners = [origin, positions[j as usize], positions[k as usize]];
                bounding_sphere_of(
                    corners
                        .iter()
                        .flat_map(|&c| verts_w.iter().map(move |&v| c * (1.0 - s) + v * s)),
                )
            });
            out.push((
                [key[0], key[1], key[2], u32::MAX],
                SimplexRecord {
                    patch,
                    nu_min2,
                    nu_max2,
                    mixed_bound,
                    sources: [me, u32::MAX, u32::MAX, u32::MAX],
                    source_count: 1,
                },
            ));
        }
    }

    // Tetrahedra: nu is a polytope vertex three faces share.
    let mut at_vertex: std::collections::HashMap<u32, Vec<u32>> = std::collections::HashMap::new();
    for f in &real {
        for &v in &f.verts {
            at_vertex.entry(v).or_default().push(f.neighbor.unwrap());
        }
    }
    for (v, ns) in at_vertex {
        if ns.len() < 3 || touches(v) {
            continue;
        }
        let v_w = world(v);
        for a in 0..ns.len() {
            for b in a + 1..ns.len() {
                for c in b + 1..ns.len() {
                    let key = sorted4([me, ns[a], ns[b], ns[c]]);
                    let Some(patch) = tetrahedron_patch(key, positions, weights, s) else {
                        continue;
                    };
                    let d2 = (v_w - patch.center).length_squared();
                    let corners = key.map(|i| positions[i as usize]);
                    let mixed_bound = Some(bounding_sphere_of(
                        corners.iter().map(|&c| c * (1.0 - s) + v_w * s),
                    ));
                    out.push((
                        key,
                        SimplexRecord {
                            patch,
                            nu_min2: d2,
                            nu_max2: d2,
                            mixed_bound,
                            sources: [me, u32::MAX, u32::MAX, u32::MAX],
                            source_count: 1,
                        },
                    ));
                }
            }
        }
    }
    (out, my_neighbors)
}

/// Applies the visibility test (see [`build_complex`]) and the tighter
/// of the two bounding spheres to a merged record; returns the patch
/// with its deduplicated competitor list (empty when it has no surface).
fn finish_patch(
    rec: SimplexRecord,
    neighbors: &[Vec<u32>],
    positions: &[Vec3],
    weights: &[f32],
    shrink: f32,
) -> (Patch, Vec<u32>) {
    let mut patch = rec.patch;
    if !patch.has_surface() {
        return (patch, Vec::new());
    }
    let z = patch.center;
    let w_max = patch
        .members()
        .iter()
        .map(|&i| weights[i as usize])
        .fold(f32::NEG_INFINITY, f32::max);
    let g_min = shrink * rec.nu_min2 - (1.0 - shrink) * (w_max + patch.weight) + patch.weight;
    let d2_simplex = match patch.kind {
        PatchKind::Vertex | PatchKind::Triangle | PatchKind::Tetrahedron => 0.0,
        PatchKind::Edge => distance2_to_segment(
            z,
            positions[patch.atoms[0] as usize],
            positions[patch.atoms[1] as usize],
        ),
    };
    let g_max = if rec.nu_max2.is_infinite() {
        f32::INFINITY
    } else {
        shrink * rec.nu_max2 - (1.0 - shrink) * d2_simplex + patch.weight
    };
    if g_min > VISIBILITY_EPS || g_max < -VISIBILITY_EPS {
        patch.bound_radius = 0.0;
        return (patch, Vec::new());
    }
    if let Some((center, radius)) = rec.mixed_bound {
        if radius < patch.bound_radius {
            patch.bound_center = center;
            patch.bound_radius = radius;
        }
    }
    let mut competitors: Vec<u32> = rec.sources[..rec.source_count as usize]
        .iter()
        .flat_map(|&src| neighbors[src as usize].iter().copied())
        .collect();
    competitors.sort_unstable();
    competitors.dedup();
    competitors.retain(|m| !patch.members().contains(m));
    (patch, competitors)
}

fn bounding_sphere_of(points: impl Iterator<Item = Vec3> + Clone) -> (Vec3, f32) {
    let (sum, n) = points
        .clone()
        .fold((Vec3::ZERO, 0usize), |(s, n), p| (s + p, n + 1));
    let center = sum / n.max(1) as f32;
    let radius = points.map(|p| (p - center).length()).fold(0.0, f32::max);
    (center, radius + 1e-3)
}

fn distance2_to_segment(z: Vec3, a: Vec3, b: Vec3) -> f32 {
    let ab = b - a;
    let l2 = ab.length_squared();
    if l2 < 1e-12 {
        return (z - a).length_squared();
    }
    let t = ((z - a).dot(ab) / l2).clamp(0.0, 1.0);
    (a + ab * t - z).length_squared()
}

/// Squared distance from `z` (in the polygon's plane) to a convex
/// polygon given by its cyclic corners: zero inside, else the nearest
/// edge.
fn distance2_to_convex_polygon(z: Vec3, verts: &[Vec3], normal: Vec3) -> f32 {
    let n = verts.len();
    if n < 3 {
        return verts
            .iter()
            .map(|&v| (v - z).length_squared())
            .fold(f32::INFINITY, f32::min);
    }
    let (mut all_pos, mut all_neg) = (true, true);
    let mut best = f32::INFINITY;
    for i in 0..n {
        let (a, b) = (verts[i], verts[(i + 1) % n]);
        let side = (b - a).cross(z - a).dot(normal);
        all_pos &= side >= 0.0;
        all_neg &= side <= 0.0;
        best = best.min(distance2_to_segment(z, a, b));
    }
    if all_pos || all_neg {
        0.0
    } else {
        best
    }
}

fn sorted2(mut v: [u32; 2]) -> [u32; 2] {
    v.sort_unstable();
    v
}

fn sorted3(mut v: [u32; 3]) -> [u32; 3] {
    v.sort_unstable();
    v
}

fn sorted4(mut v: [u32; 4]) -> [u32; 4] {
    v.sort_unstable();
    v
}

fn vertex_patch(a: u32, positions: &[Vec3], weights: &[f32], shrink: f32) -> Patch {
    let (p, w) = (positions[a as usize], weights[a as usize]);
    let r2 = shrink * w;
    Patch {
        kind: PatchKind::Vertex,
        atoms: [a, u32::MAX, u32::MAX, u32::MAX],
        center: p,
        axis: Vec3::ZERO,
        weight: -w,
        bound_center: p,
        bound_radius: if r2 > 0.0 { r2.sqrt() } else { 0.0 },
    }
}

/// `2 u . e_i = |e_i|^2 - w_i + w_0` for each edge `e_i = z_i - z_0` of
/// the simplex, the linear system the orthocenter `z_0 + u` solves
/// (`u` in the simplex's own span).
fn edge_patch(atoms: [u32; 2], positions: &[Vec3], weights: &[f32], shrink: f32) -> Option<Patch> {
    let [a, b] = atoms;
    let (pa, pb) = (positions[a as usize], positions[b as usize]);
    let (wa, wb) = (weights[a as usize], weights[b as usize]);
    let e = pb - pa;
    let len = e.length();
    if len < 1e-6 {
        return None;
    }
    let axis = e / len;
    let t_z = (len * len - wb + wa) / (2.0 * len);
    let center = pa + axis * t_z;
    let weight = t_z * t_z - wa;
    let k = shrink / (1.0 - shrink);
    // The window of axial offsets `d = (x - z) . axis` inside the cell:
    // the segment's own extent about `z`, scaled by `1 - s`.
    let (d_lo, d_hi) = (-(1.0 - shrink) * t_z, (1.0 - shrink) * (len - t_z));
    let d_max = d_lo.abs().max(d_hi.abs());
    let rho2 = k * d_max * d_max - shrink * weight;
    let bound_radius = if rho2 > 0.0 {
        let half = 0.5 * (d_hi - d_lo);
        (half * half + rho2).sqrt()
    } else {
        0.0
    };
    Some(Patch {
        kind: PatchKind::Edge,
        atoms: [a, b, u32::MAX, u32::MAX],
        center,
        axis,
        weight,
        bound_center: center + axis * (0.5 * (d_lo + d_hi)),
        bound_radius,
    })
}

fn triangle_patch(
    atoms: [u32; 3],
    positions: &[Vec3],
    weights: &[f32],
    shrink: f32,
) -> Option<Patch> {
    let [a, b, c] = atoms;
    let (pa, pb, pc) = (
        positions[a as usize],
        positions[b as usize],
        positions[c as usize],
    );
    let (wa, wb, wc) = (
        weights[a as usize],
        weights[b as usize],
        weights[c as usize],
    );
    let (e1, e2) = (pb - pa, pc - pa);
    let (g11, g12, g22) = (e1.dot(e1), e1.dot(e2), e2.dot(e2));
    let det = g11 * g22 - g12 * g12;
    if det <= 1e-8 * g11 * g22 {
        return None;
    }
    let (b1, b2) = (0.5 * (g11 - wb + wa), 0.5 * (g22 - wc + wa));
    let alpha = (b1 * g22 - b2 * g12) / det;
    let beta = (g11 * b2 - g12 * b1) / det;
    let u = e1 * alpha + e2 * beta;
    let center = pa + u;
    let weight = u.length_squared() - wa;
    let axis = e1.cross(e2).normalize();
    let k = shrink / (1.0 - shrink);
    let rho_max = (1.0 - shrink)
        * [pa, pb, pc]
            .iter()
            .map(|p| (*p - center).length())
            .fold(0.0f32, f32::max);
    let d2 = k * rho_max * rho_max - shrink * weight;
    let bound_radius = if d2 > 0.0 {
        (rho_max * rho_max + d2).sqrt()
    } else {
        0.0
    };
    Some(Patch {
        kind: PatchKind::Triangle,
        atoms: [a, b, c, u32::MAX],
        center,
        axis,
        weight,
        bound_center: center,
        bound_radius,
    })
}

fn tetrahedron_patch(
    atoms: [u32; 4],
    positions: &[Vec3],
    weights: &[f32],
    shrink: f32,
) -> Option<Patch> {
    let [a, b, c, d] = atoms;
    let pa = positions[a as usize];
    let wa = weights[a as usize];
    let es = [
        positions[b as usize] - pa,
        positions[c as usize] - pa,
        positions[d as usize] - pa,
    ];
    let rhs = [
        0.5 * (es[0].length_squared() - weights[b as usize] + wa),
        0.5 * (es[1].length_squared() - weights[c as usize] + wa),
        0.5 * (es[2].length_squared() - weights[d as usize] + wa),
    ];
    // Solve `E^T u = rhs` for `u` directly (rows of `E^T` are the edges).
    let det = es[0].dot(es[1].cross(es[2]));
    let scale = es[0].length() * es[1].length() * es[2].length();
    if det.abs() <= 1e-6 * scale {
        return None;
    }
    let u =
        (es[1].cross(es[2]) * rhs[0] + es[2].cross(es[0]) * rhs[1] + es[0].cross(es[1]) * rhs[2])
            / det;
    let center = pa + u;
    let weight = u.length_squared() - wa;
    let r2 = (1.0 - shrink) * weight;
    Some(Patch {
        kind: PatchKind::Tetrahedron,
        atoms: [a, b, c, d],
        center,
        axis: Vec3::ZERO,
        weight,
        bound_center: center,
        bound_radius: if r2 > 0.0 { r2.sqrt() } else { 0.0 },
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    struct SplitMix64(u64);
    impl SplitMix64 {
        fn next_u64(&mut self) -> u64 {
            self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut z = self.0;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            z ^ (z >> 31)
        }
        fn next_f32(&mut self) -> f32 {
            (self.next_u64() >> 40) as f32 / (1u64 << 24) as f32
        }
        fn vec3(&mut self, spread: f32) -> Vec3 {
            Vec3::new(
                (self.next_f32() * 2.0 - 1.0) * spread,
                (self.next_f32() * 2.0 - 1.0) * spread,
                (self.next_f32() * 2.0 - 1.0) * spread,
            )
        }
    }

    fn orthocenter_is_equidistant(patch: &Patch, positions: &[Vec3], weights: &[f32]) {
        let pd: Vec<f32> = patch
            .members()
            .iter()
            .map(|&i| power_distance(patch.center, positions[i as usize], weights[i as usize]))
            .collect();
        for &d in &pd {
            assert!(
                (d - patch.weight).abs() < 1e-3,
                "power distances {pd:?} vs weight {}",
                patch.weight
            );
        }
    }

    #[test]
    fn skin_weight_puts_an_isolated_atom_on_its_van_der_waals_sphere() {
        for shrink in [0.3, 0.5, 0.8] {
            let w = weight_for_radius(1.7, shrink);
            let patch = vertex_patch(0, &[Vec3::ZERO], &[w], shrink);
            assert!((patch.bound_radius - 1.7).abs() < 1e-5);
            // On the sphere, F == 0; inside negative; outside positive.
            assert!(patch.field(Vec3::X * 1.7, shrink).abs() < 1e-4);
            assert!(patch.field(Vec3::X * 1.0, shrink) < 0.0);
            assert!(patch.field(Vec3::X * 2.0, shrink) > 0.0);
            assert!(patch.normal_at(Vec3::X * 1.7, shrink).dot(Vec3::X) > 0.999);
        }
    }

    #[test]
    fn orthocenters_are_equidistant_in_power_distance_for_every_simplex_kind() {
        let mut rng = SplitMix64(7);
        for _ in 0..50 {
            let positions: Vec<Vec3> = (0..4).map(|_| rng.vec3(3.0)).collect();
            let weights: Vec<f32> = (0..4).map(|_| 1.0 + rng.next_f32() * 4.0).collect();
            if let Some(p) = edge_patch([0, 1], &positions, &weights, 0.5) {
                orthocenter_is_equidistant(&p, &positions, &weights);
                assert!(
                    (p.center - positions[0]).cross(p.axis).length() < 1e-4,
                    "edge center off the edge line"
                );
            }
            if let Some(p) = triangle_patch([0, 1, 2], &positions, &weights, 0.5) {
                orthocenter_is_equidistant(&p, &positions, &weights);
                assert!(
                    (p.center - positions[0]).dot(p.axis).abs() < 1e-3,
                    "triangle center off the plane"
                );
            }
            if let Some(p) = tetrahedron_patch([0, 1, 2, 3], &positions, &weights, 0.5) {
                orthocenter_is_equidistant(&p, &positions, &weights);
            }
        }
    }

    #[test]
    fn edge_waist_coincides_with_feasible_cells_own_bisector_plane() {
        // Cross-check between two independently-built modules: the
        // hyperboloid's symmetry point (its orthocenter) sits on the
        // same plane `feasible_cell::power_bisector` computes between
        // the same two atoms.
        let a = Vec3::new(1.0, -2.0, 0.5);
        let b = Vec3::new(5.0, 0.5, -1.0);
        let (wa, wb) = (2.0, 5.0);
        let patch = edge_patch([0, 1], &[a, b], &[wa, wb], 0.5).unwrap();
        let (_, bisector_offset) = crate::feasible_cell::power_bisector(wa, wb, b - a);
        let waist = (patch.center - a).dot(patch.axis);
        assert!(
            (waist - bisector_offset).abs() < 1e-4,
            "waist at {waist} vs bisector at {bisector_offset}"
        );
    }

    #[test]
    fn equal_weight_edge_patch_is_an_hourglass_symmetric_about_the_midpoint() {
        let a = Vec3::new(-2.0, 1.0, 0.0);
        let b = Vec3::new(4.0, 1.0, 0.0);
        let w = 12.0;
        let patch = edge_patch([0, 1], &[a, b], &[w, w], 0.5).unwrap();
        let rho_at = |t: f32| {
            let on_axis = a + patch.axis * t;
            let (mut lo, mut hi) = (0.0f32, 20.0f32);
            for _ in 0..60 {
                let mid = 0.5 * (lo + hi);
                if patch.field(on_axis + Vec3::Y * mid, 0.5) < 0.0 {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            0.5 * (lo + hi)
        };
        let waist = rho_at(3.0);
        for offset in [1.0f32, 2.0, 3.0] {
            let (toward_a, toward_b) = (rho_at(3.0 - offset), rho_at(3.0 + offset));
            assert!(
                toward_a > waist && toward_b > waist,
                "widening away from the waist"
            );
            assert!(
                (toward_a - toward_b).abs() < 1e-2,
                "symmetric about the midpoint"
            );
        }
    }

    /// The load-bearing check for the whole module, for all four kinds
    /// at once: the closed-form `F <= 0` must agree with "inside the
    /// union of the shrunk spheres of the simplex's *affine* family",
    /// sampled brute-force with the paper's own combination rule
    /// (`sum(gamma_i f_i)`, `sum(gamma_i) = 1`, gamma unrestricted in
    /// sign). Shared-premise mistakes are exactly what this module got
    /// wrong before, so the brute-force side uses nothing but the
    /// primary source's definitions and f64.
    #[test]
    fn closed_form_matches_the_brute_force_envelope_of_the_affine_family() {
        let shrink = 0.5f32;
        let mut rng = SplitMix64(11);
        let combined = |gamma: &[f64], z: &[Vec3], w: &[f32]| -> (Vec3, f64) {
            // f_gamma(x) = |x|^2 - 2 x.c + sum gamma_i (|z_i|^2 - w_i)
            //            = |x - c|^2 - W,  W = |c|^2 - sum gamma_i (|z_i|^2 - w_i)
            let mut c = glam::DVec3::ZERO;
            let mut k = 0.0f64;
            for (i, &g) in gamma.iter().enumerate() {
                let zi = z[i].as_dvec3();
                c += zi * g;
                k += g * (zi.length_squared() - w[i] as f64);
            }
            (c.as_vec3(), c.length_squared() - k)
        };
        // Affine (barycentric) coordinates of a point of `aff delta_X`
        // with respect to the simplex, so a probe can be skipped when the
        // minimizing family member would fall outside the sampled range.
        let affine_coords = |p: Vec3, positions: &[Vec3]| -> Vec<f64> {
            let n = positions.len();
            let a = positions[0].as_dvec3();
            let es: Vec<glam::DVec3> = positions[1..].iter().map(|z| z.as_dvec3() - a).collect();
            let r = p.as_dvec3() - a;
            // Least squares on the Gram system (exact for a point in the hull).
            let mut g = vec![vec![0.0f64; n - 1]; n - 1];
            let mut rhs = vec![0.0f64; n - 1];
            for i in 0..n - 1 {
                rhs[i] = r.dot(es[i]);
                for j in 0..n - 1 {
                    g[i][j] = es[i].dot(es[j]);
                }
            }
            // Gaussian elimination, tiny system.
            for i in 0..n - 1 {
                let piv = g[i][i];
                for j in i + 1..n - 1 {
                    let f = g[j][i] / piv;
                    let (above, below) = g.split_at_mut(j);
                    for (x, p) in below[0][i..n - 1].iter_mut().zip(&above[i][i..n - 1]) {
                        *x -= f * p;
                    }
                    rhs[j] -= f * rhs[i];
                }
            }
            let mut lam = vec![0.0f64; n - 1];
            for i in (0..n - 1).rev() {
                let mut s = rhs[i];
                for k in i + 1..n - 1 {
                    s -= g[i][k] * lam[k];
                }
                lam[i] = s / g[i][i];
            }
            let mut out = vec![1.0 - lam.iter().sum::<f64>()];
            out.extend(lam);
            out
        };
        let mut checked = 0;
        for kind in [
            PatchKind::Vertex,
            PatchKind::Edge,
            PatchKind::Triangle,
            PatchKind::Tetrahedron,
        ] {
            let n = kind.simplex_size();
            // Coarser gamma grids for the higher-dimensional families
            // (sample counts grow as steps^(n-1)); the probe-skipping
            // threshold below is set from the worst-case quantization
            // error of each (see the assert message).
            let (range, steps, skip_below): (f64, i32, f32) = match kind {
                PatchKind::Vertex => (1.0, 1, 0.02),
                PatchKind::Edge => (6.0, 120, 0.1),
                PatchKind::Triangle => (4.0, 64, 0.2),
                PatchKind::Tetrahedron => (3.0, 36, 0.3),
            };
            for _ in 0..10 {
                // Well-separated simplices (no edge shorter than 1 Å), so
                // the gamma grid's step is a small fraction of any edge.
                let positions: Vec<Vec3> = loop {
                    let p: Vec<Vec3> = (0..n).map(|_| rng.vec3(2.0)).collect();
                    let ok = (0..n).all(|i| (i + 1..n).all(|j| p[i].distance(p[j]) > 1.0));
                    if ok {
                        break p;
                    }
                };
                let weights: Vec<f32> = (0..n).map(|_| 2.0 + rng.next_f32() * 6.0).collect();
                let patch = match kind {
                    PatchKind::Vertex => Some(vertex_patch(0, &positions, &weights, shrink)),
                    PatchKind::Edge => edge_patch([0, 1], &positions, &weights, shrink),
                    PatchKind::Triangle => triangle_patch([0, 1, 2], &positions, &weights, shrink),
                    PatchKind::Tetrahedron => {
                        tetrahedron_patch([0, 1, 2, 3], &positions, &weights, shrink)
                    }
                };
                let Some(patch) = patch else { continue };
                // Sample gamma on a grid around the simplex (the affine
                // family, not just the convex hull); the last coordinate
                // is 1 - sum of the others.
                let mut family: Vec<(Vec3, f64)> = Vec::new();
                let mut idx = vec![0i32; n - 1];
                loop {
                    let mut gamma: Vec<f64> = idx
                        .iter()
                        .map(|&i| -range + 2.0 * range * i as f64 / steps as f64)
                        .collect();
                    let last = 1.0 - gamma.iter().sum::<f64>();
                    gamma.push(last);
                    family.push(combined(&gamma, &positions, &weights));
                    if n == 1 {
                        break;
                    }
                    let mut carry = 0;
                    while carry < n - 1 {
                        idx[carry] += 1;
                        if idx[carry] <= steps {
                            break;
                        }
                        idx[carry] = 0;
                        carry += 1;
                    }
                    if carry == n - 1 {
                        break;
                    }
                }
                for _ in 0..40 {
                    let x = patch.center + rng.vec3(2.0);
                    let f = patch.field(x, shrink);
                    if f.abs() < skip_below {
                        continue; // too close to the surface to resolve on the gamma grid
                    }
                    // Skip probes whose minimizing member lies outside the
                    // sampled gamma range (the last coordinate is implied).
                    let c = patch.delaunay_point(x, shrink);
                    let coords = affine_coords(c, &positions);
                    if n > 1 && coords[1..].iter().any(|g| g.abs() > range - 0.3) {
                        continue;
                    }
                    let inside = family.iter().any(|&(c, w)| {
                        let sw = shrink as f64 * w;
                        sw > 0.0 && (x - c).as_dvec3().length_squared() <= sw
                    });
                    assert_eq!(
                        f < 0.0,
                        inside,
                        "{kind:?}: F={f} at {x:?} but brute force says inside={inside}"
                    );
                    checked += 1;
                }
            }
        }
        assert!(checked > 400, "expected plenty of probes, got {checked}");
    }

    #[test]
    fn tetrahedron_void_has_the_predicted_radius_and_faces_inward() {
        // Four equal atoms at a regular tetrahedron's corners, far enough
        // apart that their orthocenter (the centroid) is uncovered: the
        // patch is a spherical void of radius sqrt((1 - s) w_z).
        let r = 2.0;
        let positions = vec![
            Vec3::new(1.0, 1.0, 1.0).normalize() * r,
            Vec3::new(1.0, -1.0, -1.0).normalize() * r,
            Vec3::new(-1.0, 1.0, -1.0).normalize() * r,
            Vec3::new(-1.0, -1.0, 1.0).normalize() * r,
        ];
        let weights = vec![1.0; 4];
        let patch = tetrahedron_patch([0, 1, 2, 3], &positions, &weights, 0.5).unwrap();
        assert!(
            patch.center.length() < 1e-4,
            "orthocenter should be the centroid"
        );
        assert!((patch.weight - (r * r - 1.0)).abs() < 1e-4);
        assert!((patch.bound_radius - (0.5 * patch.weight).sqrt()).abs() < 1e-4);
        let on = Vec3::X * patch.bound_radius;
        assert!(patch.field(on, 0.5).abs() < 1e-4);
        assert!(
            patch.field(Vec3::ZERO, 0.5) > 0.0,
            "the void's center is outside the body"
        );
        assert!(
            patch.normal_at(on, 0.5).dot(Vec3::X) < -0.999,
            "outward (into the void) normal"
        );
    }

    /// Rays along a hyperboloid's asymptotic cone (`a -> 0`) from a
    /// camera far away: the textbook formula loses the near root here;
    /// the q-form plus re-origining must still land on the zero set.
    #[test]
    fn ray_roots_survive_the_asymptotic_cone_direction_from_far_away() {
        let a = Vec3::new(0.0, 0.0, 0.0);
        let b = Vec3::new(3.0, 0.0, 0.0);
        let patch = edge_patch([0, 1], &[a, b], &[6.0, 6.0], 0.5).unwrap();
        // At shrink 1/2 the cone half-angle is 45 degrees: dir with
        // dd^2 = 1/2 (+ a hair) gives |a| ~ 1e-4..1e-6.
        let mut hits = 0;
        for i in 0..200 {
            let wobble = 1e-6 * i as f32;
            let dir = Vec3::new(
                std::f32::consts::FRAC_1_SQRT_2 + wobble,
                std::f32::consts::FRAC_1_SQRT_2 - wobble,
                0.0,
            )
            .normalize();
            // Aim through a point on the waist circle, from 60 Å back.
            let target = patch.center + Vec3::Y * 1.0;
            let origin = target - dir * 60.0;
            // The second root runs off along the asymptote (to ~1e7 Å
            // as `a -> 0`); only roots within a sane distance are
            // meaningful, and the one near the target must be found.
            for t in patch
                .ray_roots(origin, dir, 0.5)
                .into_iter()
                .flatten()
                .filter(|&t| t < 200.0)
            {
                let f = patch.field(origin + dir * t, 0.5);
                assert!(f.abs() < 2e-2, "root {t} has F={f} (dir wobble {wobble})");
                if (t - 60.0).abs() < 5.0 {
                    hits += 1;
                }
            }
        }
        assert!(
            hits >= 200,
            "expected the near hit for every cone-direction ray, got {hits}"
        );
    }

    #[test]
    fn ray_roots_land_on_the_zero_set_and_come_ascending() {
        let mut rng = SplitMix64(5);
        let positions: Vec<Vec3> = (0..3).map(|_| rng.vec3(2.0)).collect();
        let weights = vec![3.0, 4.0, 5.0];
        let patch = triangle_patch([0, 1, 2], &positions, &weights, 0.5).unwrap();
        let mut hits = 0;
        for _ in 0..200 {
            let origin = patch.center + rng.vec3(6.0);
            let dir = (patch.center + rng.vec3(1.0) - origin).normalize();
            let roots = patch.ray_roots(origin, dir, 0.5);
            if let [Some(t1), Some(t2)] = roots {
                assert!(t1 < t2);
            }
            for t in roots.into_iter().flatten() {
                let f = patch.field(origin + dir * t, 0.5);
                assert!(f.abs() < 1e-2, "root {t} has F={f}");
                hits += 1;
            }
        }
        assert!(hits > 50, "expected plenty of crossings, got {hits}");
    }

    /// The mixed cells tile space: with every simplex of a small cluster
    /// present (a window far larger than the cluster), each probe point
    /// is claimed by exactly one patch's cell. Uses brute-force
    /// ownership so this tests the membership math, not the grid.
    #[test]
    fn mixed_cells_tile_space_around_a_small_cluster() {
        let mut rng = SplitMix64(3);
        let positions: Vec<Vec3> = (0..14).map(|_| rng.vec3(2.5)).collect();
        let weights: Vec<f32> = (0..14).map(|_| 1.5 + rng.next_f32() * 2.0).collect();
        let shrink = 0.5;
        let complex = build_complex(&positions, &weights, shrink);
        let mut all = [0; 4];
        for p in &complex.patches {
            all[p.kind as usize] += 1;
        }
        assert!(
            all[1] > 20 && all[2] > 20 && all[3] > 5,
            "expected a real complex, got {all:?}"
        );
        let owned = |q: Vec3, members: &[u32]| {
            let bound = power_distance(
                q,
                positions[members[0] as usize],
                weights[members[0] as usize],
            ) - MEMBERSHIP_EPS;
            (0..positions.len() as u32).all(|m| {
                members.contains(&m)
                    || power_distance(q, positions[m as usize], weights[m as usize]) >= bound
            })
        };
        let (mut once, mut never, mut multi) = (0, 0, 0);
        for _ in 0..600 {
            let x = rng.vec3(3.0);
            let n = complex
                .patches
                .iter()
                .filter(|p| p.contains(x, shrink, &positions, |q| owned(q, p.members())))
                .count();
            match n {
                0 => never += 1,
                1 => once += 1,
                _ => multi += 1,
            }
        }
        println!("tiling: once={once} never={never} multi={multi}");
        assert_eq!(
            never, 0,
            "a point no mixed cell claims means a missing simplex"
        );
        // Double claims come from the deliberate membership tolerance at
        // cell boundaries (see `MEMBERSHIP_EPS`); a few percent of
        // uniformly random points is that band, not a tiling error.
        assert!(
            multi <= 45,
            "cells should only overlap within the membership tolerance, got {multi} multi-claims"
        );
    }

    /// The competitor-list ownership test must agree with the exhaustive
    /// one everywhere the cells are complete -- the check behind
    /// `SkinComplex::owned_by_competitors`'s neighbor-sufficiency
    /// argument, run on real data rather than trusted.
    #[test]
    fn competitor_lists_reproduce_exhaustive_ownership() {
        let mut rng = SplitMix64(17);
        let positions: Vec<Vec3> = (0..16).map(|_| rng.vec3(2.5)).collect();
        let weights: Vec<f32> = (0..16).map(|_| 1.5 + rng.next_f32() * 2.0).collect();
        let complex = build_complex(&positions, &weights, 0.5);
        let exhaustive = |q: Vec3, members: &[u32]| {
            let bound = power_distance(
                q,
                positions[members[0] as usize],
                weights[members[0] as usize],
            ) - MEMBERSHIP_EPS;
            (0..positions.len() as u32).all(|m| {
                members.contains(&m)
                    || power_distance(q, positions[m as usize], weights[m as usize]) >= bound
            })
        };
        let (mut checked, mut owned) = (0, 0);
        for (i, p) in complex.patches.iter().enumerate() {
            if !p.has_surface() {
                continue;
            }
            for _ in 0..30 {
                let q = p.voronoi_point(p.bound_center + rng.vec3(p.bound_radius), 0.5);
                let a = complex.owned_by_competitors(i, &positions, &weights, q);
                let b = exhaustive(q, p.members());
                assert_eq!(
                    a,
                    b,
                    "patch {i} ({:?}, atoms {:?}) at q={q:?}",
                    p.kind,
                    p.members()
                );
                checked += 1;
                owned += a as usize;
            }
        }
        assert!(
            checked > 300 && owned > 20 && owned < checked,
            "checked {checked}, owned {owned}"
        );
    }

    #[test]
    fn delaunay_and_voronoi_points_recompose_to_x() {
        let mut rng = SplitMix64(9);
        let positions: Vec<Vec3> = (0..4).map(|_| rng.vec3(2.0)).collect();
        let weights = vec![2.0, 3.0, 2.5, 4.0];
        let shrink = 0.4;
        let patches = [
            vertex_patch(0, &positions, &weights, shrink),
            edge_patch([0, 1], &positions, &weights, shrink).unwrap(),
            triangle_patch([0, 1, 2], &positions, &weights, shrink).unwrap(),
            tetrahedron_patch([0, 1, 2, 3], &positions, &weights, shrink).unwrap(),
        ];
        for patch in &patches {
            for _ in 0..20 {
                let x = patch.center + rng.vec3(3.0);
                let p = patch.delaunay_point(x, shrink);
                let q = patch.voronoi_point(x, shrink);
                let back = p * (1.0 - shrink) + q * shrink;
                assert!(
                    (back - x).length() < 1e-4,
                    "{:?}: {x:?} -> {back:?}",
                    patch.kind
                );
                // q lies in aff nu_X: equal power distance to every member.
                let pd: Vec<f32> = patch
                    .members()
                    .iter()
                    .map(|&i| power_distance(q, positions[i as usize], weights[i as usize]))
                    .collect();
                for d in &pd {
                    assert!(
                        (d - pd[0]).abs() < 1e-3,
                        "{:?}: q not equidistant: {pd:?}",
                        patch.kind
                    );
                }
            }
        }
    }

    /// A buried atom's sphere never reaches its own (tiny) mixed cell, so
    /// its patch is culled; every atom on the outside keeps one. Checked
    /// by brute force against the real definition (does any point of the
    /// sphere pass `contains`?), not by trusting the analytic test.
    #[test]
    fn visibility_culling_matches_a_brute_force_search_of_each_sphere() {
        // A dense cubic cluster: the center atom is fully buried.
        let mut positions = Vec::new();
        for x in -2..=2 {
            for y in -2..=2 {
                for z in -2..=2 {
                    positions.push(Vec3::new(x as f32, y as f32, z as f32) * 1.6);
                }
            }
        }
        let weights = vec![weight_for_radius(1.5, 0.5); positions.len()];
        let complex = build_complex(&positions, &weights, 0.5);
        let dirs: Vec<Vec3> = {
            let mut rng = SplitMix64(21);
            (0..400)
                .map(|_| rng.vec3(1.0).normalize_or(Vec3::X))
                .collect()
        };
        let (mut culled, mut kept) = (0, 0);
        for p in complex
            .patches
            .iter()
            .filter(|p| p.kind == PatchKind::Vertex)
        {
            let r = (0.5 * weights[p.atoms[0] as usize]).sqrt();
            let any_visible = dirs.iter().any(|&d| {
                let x = p.center + d * r;
                p.contains(x, 0.5, &positions, |q| {
                    complex.voronoi_owned_by(&positions, &weights, q, p.members())
                })
            });
            if p.has_surface() {
                kept += 1;
            } else {
                culled += 1;
                assert!(
                    !any_visible,
                    "atom {} was culled but a point of its sphere is in its cell",
                    p.atoms[0]
                );
            }
        }
        let center = complex
            .patches
            .iter()
            .find(|p| p.kind == PatchKind::Vertex && p.atoms[0] == 62)
            .unwrap();
        assert!(!center.has_surface(), "the central atom is buried");
        assert!(culled >= 1 && kept >= 8 * 8, "culled {culled}, kept {kept}");
    }

    #[test]
    fn a_dominated_atom_still_gets_a_vertex_patch_but_owns_no_point() {
        // A tiny atom inside a huge neighbor: its power cell is empty, so
        // no point's Voronoi test can ever pick it, and its vertex patch
        // is simply never drawn.
        let positions = vec![Vec3::ZERO, Vec3::new(0.5, 0.0, 0.0)];
        let weights = vec![100.0, 0.5];
        let complex = build_complex(&positions, &weights, 0.5);
        let small = complex
            .patches
            .iter()
            .find(|p| p.kind == PatchKind::Vertex && p.atoms[0] == 1)
            .unwrap();
        for probe in [
            Vec3::new(0.5, 0.3, 0.0),
            Vec3::new(0.7, 0.0, 0.0),
            Vec3::new(0.5, 0.0, -0.4),
        ] {
            let q = small.voronoi_point(probe, 0.5);
            assert!(!complex.voronoi_owned_by(&positions, &weights, q, small.members()));
        }
    }
}
