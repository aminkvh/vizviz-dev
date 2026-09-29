//! Per-shape local meshes and their per-frame world-space assembly.
//!
//! Each SNFG glyph is built in the local orthonormal frame `(u, v, w)`
//! `super::frame` returns (`u` toward the attachment, `v` from the ring
//! oxygen, `w = u x v`) and placed at the residue's ring centroid. Sizes
//! below are proportioned from the source script's per-shape constants
//! (`cube_size = size*0.806`, `diamond_size = size*1.3`, ...; see
//! `Shape::size_factor`) but not literally its vertex-by-vertex
//! construction, which the source builds with Tcl's `trans angle`
//! rotations for a fixed-function immediate-mode renderer; an explicit
//! local-space triangle list is the equivalent for an instanced mesh.
//!
//! Two colors: recognized residues pair two different colors on the
//! cube (hexosamines: white/blue, white/green, ...), the diamond
//! (uronic acids: blue/white, white/brown, ...) and the cone
//! (deoxyhexNAc: white/blue, white/red, ...) — every other shape's
//! recognized residues use one color for both. The cube's split
//! reproduces the source script's own per-face diagonal exactly (its
//! `linked_carb.oriented.cube`, six triangles each color, one diagonal
//! per cube face); the diamond and cone use a plain half/half split,
//! simpler than the source's alternating-face pattern but the same
//! silhouette and color placement.

use glam::Vec3;

use super::{GlycanFrame, GlycanPlan, Shape};

/// A frame's glycan glyphs as one triangle soup: flat-shaded (each
/// triangle's 3 vertices repeat its own normal), so hard SNFG edges and
/// two-tone faces need no shared-vertex blending.
#[derive(Default)]
pub struct PolytopeMesh {
    pub positions: Vec<Vec3>,
    pub normals: Vec<Vec3>,
    pub colors: Vec<[u8; 3]>,
    /// The residue's anomeric carbon, one entry per vertex (matches
    /// `positions`): what a pick on this triangle names.
    pub source_atom: Vec<u32>,
}

impl PolytopeMesh {
    fn clear(&mut self) {
        self.positions.clear();
        self.normals.clear();
        self.colors.clear();
        self.source_atom.clear();
    }

    pub fn triangle_count(&self) -> usize {
        self.positions.len() / 3
    }

    fn tri(&mut self, a: Vec3, b: Vec3, c: Vec3, color: [u8; 3], atom: u32) {
        let n = (b - a).cross(c - a).normalize_or_zero();
        for p in [a, b, c] {
            self.positions.push(p);
            self.normals.push(n);
            self.colors.push(color);
            self.source_atom.push(atom);
        }
    }
}

/// Maps local `(u, v, w)`-frame coordinates into world space.
struct Frame {
    center: Vec3,
    u: Vec3,
    v: Vec3,
    w: Vec3,
}

impl Frame {
    fn at(&self, lu: f32, lv: f32, lw: f32) -> Vec3 {
        self.center + self.u * lu + self.v * lv + self.w * lw
    }
}

/// Rebuilds `out` for every residue `keep` accepts, at shape size `size`
/// (the `Glycan` rep's `size` option, matching the source script's
/// global `size` variable). Reuses `out`'s buffers: after the first
/// call sizes them to the frame's triangle count, replaying a
/// trajectory allocates only when that count changes.
pub fn build_mesh(
    plan: &GlycanPlan,
    frame: &GlycanFrame,
    size: f32,
    keep: impl Fn(u32) -> bool,
    out: &mut PolytopeMesh,
) {
    out.clear();
    for (i, res) in plan.residues.iter().enumerate() {
        if !keep(res.residue) {
            continue;
        }
        let center = frame.centroid[i];
        let (u, v, w) = super::frame(center, frame.attach_point[i], frame.ring_oxygen[i]);
        let f = Frame { center, u, v, w };
        let s = size * res.shape.size_factor();
        let (c1, c2, atom) = (res.color1.rgb(), res.color2.rgb(), res.anomeric());
        match res.shape {
            Shape::Sphere => icosahedron(out, &f, s, c1, atom),
            Shape::Cube => box_shape(out, &f, Vec3::splat(s / 2.0), c1, c2, atom),
            // Elongated, not cubic: `v` (the ring-pucker axis) longest,
            // `u` (toward the linkage) thinnest so it reads apart from a
            // connecting cylinder.
            Shape::Rectangle => box_shape(out, &f, Vec3::new(0.15, 0.45, 0.25) * s, c1, c1, atom),
            Shape::Diamond => diamond(out, &f, s, c1, c2, atom),
            Shape::FlatDiamond => flat_diamond(out, &f, s, c1, atom),
            Shape::Cone => cone(out, &f, s, c1, c2, atom),
            Shape::Star => star(out, &f, s, c1, atom),
            Shape::Hexagon => flat_polygon(out, &f, 6, 0.5 * s, 0.25 * s, c1, atom),
            Shape::Pentagon => flat_polygon(out, &f, 5, 0.5 * s, 0.25 * s, c1, atom),
        }
    }
}

/// The gray glycosidic-linkage cylinders and protein-attachment markers
/// (the source script's default `cylinder_radius`, `redfac = 0`: a
/// plain connector, no reduced-size cap at a free terminus). Appends to
/// `out` *without* clearing it, so a caller runs this right after
/// [`build_mesh`] into the same buffer -- one mesh, one upload, one draw
/// call for both the shapes and their connections.
pub fn build_linkage_mesh(
    plan: &GlycanPlan,
    frame: &GlycanFrame,
    radius: f32,
    keep: impl Fn(u32) -> bool,
    out: &mut PolytopeMesh,
) {
    const GRAY: [u8; 3] = [140, 140, 140];
    if radius <= 0.0 {
        return;
    }
    for (i, res) in plan.residues.iter().enumerate() {
        if !keep(res.residue) {
            continue;
        }
        let linked = matches!(
            res.attachment,
            super::Attachment::Residue(_) | super::Attachment::ProteinCa(_)
        );
        if !linked {
            continue;
        }
        let atom = res.anomeric();
        cylinder(
            out,
            frame.centroid[i],
            frame.attach_point[i],
            radius,
            GRAY,
            atom,
        );
        if matches!(res.attachment, super::Attachment::ProteinCa(_)) {
            let anchor = Frame {
                center: frame.attach_point[i],
                u: Vec3::X,
                v: Vec3::Y,
                w: Vec3::Z,
            };
            icosahedron(out, &anchor, radius, GRAY, atom);
        }
    }
}

/// A capped-free cylinder from `p0` to `p1`, `n`-sided, flat-shaded.
fn cylinder(out: &mut PolytopeMesh, p0: Vec3, p1: Vec3, radius: f32, color: [u8; 3], atom: u32) {
    const SIDES: usize = 8;
    let axis = p1 - p0;
    let len = axis.length();
    if len < 1e-6 {
        return;
    }
    let dir = axis / len;
    let seed = if dir.x.abs() < 0.9 { Vec3::X } else { Vec3::Y };
    let right = (seed - dir * dir.dot(seed)).normalize() * radius;
    let up = dir.cross(right);
    let ring = |center: Vec3| -> Vec<Vec3> {
        (0..SIDES)
            .map(|i| {
                let a = i as f32 / SIDES as f32 * std::f32::consts::TAU;
                center + right * a.cos() + up * a.sin()
            })
            .collect()
    };
    let (bottom, top) = (ring(p0), ring(p1));
    for i in 0..SIDES {
        let j = (i + 1) % SIDES;
        out.tri(bottom[i], bottom[j], top[j], color, atom);
        out.tri(bottom[i], top[j], top[i], color, atom);
    }
}

/// A regular icosahedron of radius `r` (Glc/Man/Gal/... spheres): 20
/// flat facets read as a ball at typical viewing distance, far cheaper
/// than a subdivided sphere for a shape this small.
fn icosahedron(out: &mut PolytopeMesh, f: &Frame, r: f32, color: [u8; 3], atom: u32) {
    let t = (1.0 + 5f32.sqrt()) / 2.0;
    // Three golden rectangles, one per axis pair; each raw vertex is
    // normalized to `r` in the local `(u, v, w)` frame.
    let verts: Vec<Vec3> = [
        (-1.0, t, 0.0),
        (1.0, t, 0.0),
        (-1.0, -t, 0.0),
        (1.0, -t, 0.0),
        (0.0, -1.0, t),
        (0.0, 1.0, t),
        (0.0, -1.0, -t),
        (0.0, 1.0, -t),
        (t, 0.0, -1.0),
        (t, 0.0, 1.0),
        (-t, 0.0, -1.0),
        (-t, 0.0, 1.0),
    ]
    .into_iter()
    .map(|(x, y, z)| {
        let local = Vec3::new(x, y, z).normalize() * r;
        f.at(local.x, local.y, local.z)
    })
    .collect();
    const FACES: [(usize, usize, usize); 20] = [
        (0, 11, 5),
        (0, 5, 1),
        (0, 1, 7),
        (0, 7, 10),
        (0, 10, 11),
        (1, 5, 9),
        (5, 11, 4),
        (11, 10, 2),
        (10, 7, 6),
        (7, 1, 8),
        (3, 9, 4),
        (3, 4, 2),
        (3, 2, 6),
        (3, 6, 8),
        (3, 8, 9),
        (4, 9, 5),
        (2, 4, 11),
        (6, 2, 10),
        (8, 6, 7),
        (9, 8, 1),
    ];
    for &(a, b, c) in &FACES {
        out.tri(verts[a], verts[b], verts[c], color, atom);
    }
}

/// A box of half-extents `half` (along `u, v, w`), each face split by
/// one diagonal between `color1` and `color2`: a cube (equal extents)
/// reproduces the source script's own crossed-cube split exactly
/// (`linked_carb.oriented.cube`'s `s1..s8` corners); unequal extents
/// give the flat rectangle its shape, monochrome (`color1 == color2`).
fn box_shape(
    out: &mut PolytopeMesh,
    f: &Frame,
    half: Vec3,
    color1: [u8; 3],
    color2: [u8; 3],
    atom: u32,
) {
    let p = |lu: f32, lv: f32, lw: f32| f.at(lu * half.x, lv * half.y, lw * half.z);
    let (s1, s2, s3, s4) = (
        p(1., 1., 1.),
        p(-1., 1., 1.),
        p(1., -1., 1.),
        p(-1., -1., 1.),
    );
    let (s5, s6, s7, s8) = (
        p(1., 1., -1.),
        p(-1., 1., -1.),
        p(1., -1., -1.),
        p(-1., -1., -1.),
    );
    for (a, b, c) in [
        (s2, s4, s3),
        (s1, s6, s2),
        (s4, s8, s7),
        (s5, s8, s6),
        (s2, s8, s4),
        (s1, s7, s5),
    ] {
        out.tri(a, b, c, color1, atom);
    }
    for (a, b, c) in [
        (s1, s2, s3),
        (s3, s4, s7),
        (s1, s5, s6),
        (s5, s7, s8),
        (s2, s6, s8),
        (s1, s3, s7),
    ] {
        out.tri(a, b, c, color2, atom);
    }
}

/// A bipyramid (octahedron): `color1` for the 4 faces to `+w`, `color2`
/// for the 4 to `-w`. The source's divided diamond instead alternates
/// color by equatorial edge on both pyramids (a pinwheel); this is a
/// plain front/back split with the same two colors and silhouette.
fn diamond(
    out: &mut PolytopeMesh,
    f: &Frame,
    size: f32,
    color1: [u8; 3],
    color2: [u8; 3],
    atom: u32,
) {
    let r = size / 2.0;
    let (e1, e2, e3, e4) = (
        f.at(r, 0., 0.),
        f.at(0., r, 0.),
        f.at(-r, 0., 0.),
        f.at(0., -r, 0.),
    );
    let top = f.at(0., 0., r);
    let bottom = f.at(0., 0., -r);
    for (a, b) in [(e1, e2), (e2, e3), (e3, e4), (e4, e1)] {
        out.tri(a, b, top, color1, atom);
    }
    for (a, b) in [(e1, e2), (e2, e3), (e3, e4), (e4, e1)] {
        out.tri(b, a, bottom, color2, atom);
    }
}

/// The di-deoxynonulosonates' "Flattened Diamond" (Pse, Leg, Aci,
/// 4eLeg): the same bipyramid as [`diamond`] but squashed along `w`, so
/// it reads apart from the sialic acids' full `Diamond` at a glance, as
/// the flat 2D symbol does from the filled one. Always monochrome: the
/// spec pairs no divided variant with this shape.
fn flat_diamond(out: &mut PolytopeMesh, f: &Frame, size: f32, color: [u8; 3], atom: u32) {
    const FLATTEN: f32 = 0.4;
    let r = size / 2.0;
    let (e1, e2, e3, e4) = (
        f.at(r, 0., 0.),
        f.at(0., r, 0.),
        f.at(-r, 0., 0.),
        f.at(0., -r, 0.),
    );
    let top = f.at(0., 0., r * FLATTEN);
    let bottom = f.at(0., 0., -r * FLATTEN);
    for (a, b) in [(e1, e2), (e2, e3), (e3, e4), (e4, e1)] {
        out.tri(a, b, top, color, atom);
        out.tri(b, a, bottom, color, atom);
    }
}

/// A cone along `u`: base (radius `0.5*size`) toward the attachment,
/// apex away from it, as the source script's does (its base sits at
/// `0.66*half_length` forward, its apex `1.33*half_length` back).
/// `color1`/`color2` split the rim at the halfway angle, the "Divided
/// Triangle" (deoxyhexNAc) shape's two-tone equivalent; monochrome
/// (`color1 == color2`) draws the plain "Filled Triangle" (deoxyhexose).
fn cone(out: &mut PolytopeMesh, f: &Frame, size: f32, color1: [u8; 3], color2: [u8; 3], atom: u32) {
    const SEGMENTS: usize = 12;
    let radius = 0.5 * size;
    let base_u = 0.3 * size;
    let apex = f.at(-0.7 * size, 0., 0.);
    let rim: Vec<Vec3> = (0..SEGMENTS)
        .map(|i| {
            let a = i as f32 / SEGMENTS as f32 * std::f32::consts::TAU;
            f.at(base_u, a.cos() * radius, a.sin() * radius)
        })
        .collect();
    let base_center = f.at(base_u, 0., 0.);
    for i in 0..SEGMENTS {
        let j = (i + 1) % SEGMENTS;
        let color = if i < SEGMENTS / 2 { color1 } else { color2 };
        out.tri(base_center, rim[i], rim[j], color, atom); // cap, outward +u
        out.tri(apex, rim[j], rim[i], color, atom); // side, outward from axis
    }
}

/// A regular `n`-gon extruded by `thickness` along `w`, lying in the
/// `u`/`v` plane (hexagon, pentagon).
fn flat_polygon(
    out: &mut PolytopeMesh,
    f: &Frame,
    n: usize,
    radius: f32,
    thickness: f32,
    color: [u8; 3],
    atom: u32,
) {
    let ring = |t: f32| -> Vec<Vec3> {
        (0..n)
            .map(|i| {
                let a = i as f32 / n as f32 * std::f32::consts::TAU;
                f.at(a.cos() * radius, a.sin() * radius, t)
            })
            .collect()
    };
    let front = ring(thickness);
    let back = ring(-thickness);
    let front_center = f.at(0., 0., thickness);
    let back_center = f.at(0., 0., -thickness);
    for i in 0..n {
        let j = (i + 1) % n;
        out.tri(front_center, front[i], front[j], color, atom);
        out.tri(back_center, back[j], back[i], color, atom);
        out.tri(front[i], back[i], back[j], color, atom);
        out.tri(front[i], back[j], front[j], color, atom);
    }
}

/// A 5-pointed star extruded by `thickness` along `w`, lying in the
/// `u`/`v` plane.
fn star(out: &mut PolytopeMesh, f: &Frame, size: f32, color: [u8; 3], atom: u32) {
    const POINTS: usize = 5;
    let thickness = 0.3 * size;
    let outer = 0.75 * size;
    let inner = 0.38 * size;
    let point = |i: usize, r: f32, t: f32| -> Vec3 {
        let a = i as f32 / (POINTS * 2) as f32 * std::f32::consts::TAU;
        f.at(a.cos() * r, a.sin() * r, t)
    };
    let ring = |t: f32| -> Vec<Vec3> {
        (0..POINTS * 2)
            .map(|i| point(i, if i % 2 == 0 { outer } else { inner }, t))
            .collect()
    };
    let front = ring(thickness);
    let back = ring(-thickness);
    let front_center = f.at(0., 0., thickness);
    let back_center = f.at(0., 0., -thickness);
    let n = POINTS * 2;
    for i in 0..n {
        let j = (i + 1) % n;
        out.tri(front_center, front[i], front[j], color, atom);
        out.tri(back_center, back[j], back[i], color, atom);
        out.tri(front[i], back[i], back[j], color, atom);
        out.tri(front[i], back[j], front[j], color, atom);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::glycan::{Attachment, GlycanResidue, SnfgColor};

    fn one_residue(
        shape: Shape,
        color1: SnfgColor,
        color2: SnfgColor,
    ) -> (GlycanPlan, GlycanFrame) {
        let residue = GlycanResidue {
            residue: 0,
            ring_atoms: [0, 1, 2, 3, 4, 5],
            shape,
            color1,
            color2,
            label: "test",
            attachment: Attachment::None,
        };
        let plan = GlycanPlan {
            residues: vec![residue],
        };
        let frame = GlycanFrame {
            centroid: vec![Vec3::ZERO],
            ring_oxygen: vec![Vec3::new(0.0, 0.0, 1.0)],
            attach_point: vec![Vec3::new(1.0, 0.0, 0.0)],
        };
        (plan, frame)
    }

    fn triangle_count(shape: Shape) -> usize {
        let (plan, frame) = one_residue(shape, SnfgColor::Blue, SnfgColor::Blue);
        let mut out = PolytopeMesh::default();
        build_mesh(&plan, &frame, 4.0, |_| true, &mut out);
        out.triangle_count()
    }

    #[test]
    fn every_shape_produces_a_closed_nonempty_mesh() {
        for shape in [
            Shape::Sphere,
            Shape::Cube,
            Shape::Diamond,
            Shape::FlatDiamond,
            Shape::Cone,
            Shape::Rectangle,
            Shape::Star,
            Shape::Hexagon,
            Shape::Pentagon,
        ] {
            let n = triangle_count(shape);
            assert!(n > 0, "{shape:?} produced no triangles");
            assert!(
                out_len_matches(shape),
                "{shape:?}: positions/normals/colors/source_atom length mismatch"
            );
        }
    }

    fn out_len_matches(shape: Shape) -> bool {
        let (plan, frame) = one_residue(shape, SnfgColor::White, SnfgColor::Blue);
        let mut out = PolytopeMesh::default();
        build_mesh(&plan, &frame, 4.0, |_| true, &mut out);
        let n = out.positions.len();
        n == out.normals.len() && n == out.colors.len() && n == out.source_atom.len() && n % 3 == 0
    }

    #[test]
    fn cube_face_diagonals_are_evenly_split_between_the_two_colors() {
        let (plan, frame) = one_residue(Shape::Cube, SnfgColor::White, SnfgColor::Blue);
        let mut out = PolytopeMesh::default();
        build_mesh(&plan, &frame, 4.0, |_| true, &mut out);
        assert_eq!(out.triangle_count(), 12);
        let white = out
            .colors
            .iter()
            .step_by(3)
            .filter(|&&c| c == SnfgColor::White.rgb())
            .count();
        let blue = out
            .colors
            .iter()
            .step_by(3)
            .filter(|&&c| c == SnfgColor::Blue.rgb())
            .count();
        assert_eq!((white, blue), (6, 6));
    }

    #[test]
    fn cone_rim_is_evenly_split_between_the_two_colors() {
        // The "Divided Triangle" (deoxyhexNAc, e.g. QuiNAc: white/blue).
        let (plan, frame) = one_residue(Shape::Cone, SnfgColor::White, SnfgColor::Blue);
        let mut out = PolytopeMesh::default();
        build_mesh(&plan, &frame, 4.0, |_| true, &mut out);
        let white = out
            .colors
            .iter()
            .step_by(3)
            .filter(|&&c| c == SnfgColor::White.rgb())
            .count();
        let blue = out
            .colors
            .iter()
            .step_by(3)
            .filter(|&&c| c == SnfgColor::Blue.rgb())
            .count();
        assert_eq!((white, blue), (12, 12), "12 segments, cap+side each");
    }

    #[test]
    fn keep_filters_by_residue_index() {
        let (plan, frame) = one_residue(Shape::Sphere, SnfgColor::Green, SnfgColor::Green);
        let mut out = PolytopeMesh::default();
        build_mesh(&plan, &frame, 4.0, |_| false, &mut out);
        assert_eq!(out.triangle_count(), 0);
    }

    #[test]
    fn convex_shapes_have_outward_facing_normals() {
        // Every shape but the star is a convex solid centered on the
        // residue centroid (the origin here): a lit face's normal must
        // point away from the center it bulges from, or shading (and
        // backface culling, if the pipeline uses it) would turn it
        // inside-out. The star's points fold the surface concave, so it
        // is exempt.
        for shape in [
            Shape::Sphere,
            Shape::Cube,
            Shape::Diamond,
            Shape::FlatDiamond,
            Shape::Cone,
            Shape::Rectangle,
            Shape::Hexagon,
            Shape::Pentagon,
        ] {
            let (plan, frame) = one_residue(shape, SnfgColor::Blue, SnfgColor::White);
            let mut out = PolytopeMesh::default();
            build_mesh(&plan, &frame, 4.0, |_| true, &mut out);
            for tri in out.positions.chunks_exact(3) {
                let centroid = (tri[0] + tri[1] + tri[2]) / 3.0;
                let normal = (tri[1] - tri[0]).cross(tri[2] - tri[0]);
                assert!(
                    normal.dot(centroid) > 0.0,
                    "{shape:?}: inward-facing triangle at {centroid:?}"
                );
            }
        }
    }

    #[test]
    fn vertices_stay_within_a_size_scaled_bound() {
        // A loose bound: every shape's local construction above should
        // keep vertices within a couple of `size` of the centroid.
        let (plan, frame) = one_residue(Shape::Star, SnfgColor::Orange, SnfgColor::Orange);
        let mut out = PolytopeMesh::default();
        build_mesh(&plan, &frame, 4.0, |_| true, &mut out);
        for p in &out.positions {
            assert!(
                p.length() < 8.0,
                "{p:?} too far from the residue centroid at origin"
            );
        }
    }

    #[test]
    fn linkage_cylinder_has_outward_facing_normals_and_is_gray() {
        let mut out = PolytopeMesh::default();
        cylinder(
            &mut out,
            Vec3::ZERO,
            Vec3::new(0.0, 0.0, 10.0),
            1.0,
            [140, 140, 140],
            7,
        );
        assert_eq!(out.triangle_count(), 16); // 8 sides * 2 triangles
        for tri in out.positions.chunks_exact(3) {
            let centroid = (tri[0] + tri[1] + tri[2]) / 3.0;
            // Distance from the axis (z), not from the origin: a
            // cylinder is convex around its axis line, not a point.
            let radial = Vec3::new(centroid.x, centroid.y, 0.0);
            let normal = (tri[1] - tri[0]).cross(tri[2] - tri[0]);
            assert!(
                normal.dot(radial) > 0.0,
                "inward-facing cylinder triangle at {centroid:?}"
            );
        }
        assert!(out.colors.iter().all(|&c| c == [140, 140, 140]));
        assert!(out.source_atom.iter().all(|&a| a == 7));
    }

    #[test]
    fn build_linkage_mesh_draws_only_resolved_links() {
        // Residue 0: linked to residue 1 (a `Residue` attachment).
        // Residue 1: `None` (a free reducing end) -- no cylinder for it.
        let residues = vec![
            GlycanResidue {
                residue: 0,
                ring_atoms: [0, 1, 2, 3, 4, 5],
                shape: Shape::Cube,
                color1: SnfgColor::Blue,
                color2: SnfgColor::Blue,
                label: "a",
                attachment: Attachment::Residue(1),
            },
            GlycanResidue {
                residue: 1,
                ring_atoms: [6, 7, 8, 9, 10, 11],
                shape: Shape::Sphere,
                color1: SnfgColor::Green,
                color2: SnfgColor::Green,
                label: "b",
                attachment: Attachment::None,
            },
        ];
        let plan = GlycanPlan { residues };
        let frame = GlycanFrame {
            centroid: vec![Vec3::ZERO, Vec3::new(5.0, 0.0, 0.0)],
            ring_oxygen: vec![Vec3::new(0.0, 0.0, 1.0), Vec3::new(5.0, 0.0, 1.0)],
            attach_point: vec![Vec3::new(5.0, 0.0, 0.0), Vec3::new(6.0, 0.0, 0.0)],
        };
        let mut out = PolytopeMesh::default();
        build_linkage_mesh(&plan, &frame, 0.5, |_| true, &mut out);
        assert_eq!(
            out.triangle_count(),
            16,
            "exactly one cylinder, no marker sphere"
        );

        // A zero radius (the rep's `radius` option at its minimum) draws
        // nothing at all, matching the source script's icon-mode default.
        let mut zeroed = PolytopeMesh::default();
        build_linkage_mesh(&plan, &frame, 0.0, |_| true, &mut zeroed);
        assert_eq!(zeroed.triangle_count(), 0);
    }
}
