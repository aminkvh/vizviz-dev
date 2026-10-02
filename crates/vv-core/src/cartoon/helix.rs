//! Helices as straight cylinders along their axes.
//!
//! The axis comes from the local-axis construction on four consecutive
//! CAs (Sugeta & Miyazawa 1967, Biopolymers 5:673; Kahn 1989, Comput.
//! Chem. 13:185; the form used by Bansal, Kumar & Velavan 2000, J. Biomol.
//! Struct. Dyn. 17:811): the second difference `r(i-1) - 2 r(i) + r(i+1)`
//! of a helix points from CA `i` at the axis, and its length is
//! `2 R (1 - cos w)` for helix radius `R` and twist `w` per residue, with
//! `w` read from two neighbouring second differences. A least-squares
//! line (the principal axis) through those local axis points is the
//! cylinder's axis. A helix whose local axis points stray from their line
//! by more than [`KINK`] is cut where they stray most and each part fitted
//! again, so a kinked or curved helix draws as a few straight cylinders.
//!
//! Each CA's spline control point is moved onto its cylinder's axis, so the
//! Catmull-Rom spline through collinear points is the straight cylinder and
//! the shared mesh, picking and colouring paths need no cylinder case.

use std::cmp::Ordering;

use glam::{Mat3, Vec3};

use super::{CartoonFrame, CartoonPlan, Shape};

/// Fewest residues in a cylinder: four CAs give two local axis points, the
/// least that define a line.
const MIN_RESIDUES: usize = 4;
/// How far, in Angstroms, local axis points may stray from their fitted
/// line before the helix is cut there.
const KINK: f32 = 1.0;
/// Local radius estimates beyond this come from a nearly straight chain,
/// where the twist is too small to divide by.
const MAX_RADIUS: f32 = 3.5;
const POWER_STEPS: usize = 24;

struct Line {
    origin: Vec3,
    dir: Vec3,
}

impl Line {
    fn project(&self, p: Vec3) -> Vec3 {
        self.origin + self.dir * (p - self.origin).dot(self.dir)
    }

    fn deviation(&self, p: Vec3) -> f32 {
        p.distance(self.project(p))
    }
}

/// The local helix axis point beside each of `ca[1..len - 1]`, for four or
/// more CAs.
fn local_axis(ca: &[Vec3]) -> Vec<Vec3> {
    let inward: Vec<Vec3> = ca.windows(3).map(|w| w[0] - w[1] * 2.0 + w[2]).collect();
    (0..inward.len())
        .map(|k| {
            let pair = k.min(inward.len() - 2);
            let cos_twist = inward[pair]
                .normalize_or_zero()
                .dot(inward[pair + 1].normalize_or_zero());
            let radius = (inward[k].length() / (2.0 * (1.0 - cos_twist)).max(1e-3)).min(MAX_RADIUS);
            ca[k + 1] + inward[k].normalize_or_zero() * radius
        })
        .collect()
}

/// The least-squares line through `points`: through their centroid along
/// the covariance matrix's largest eigenvector, found by power iteration
/// and signed to run from the first point toward the last.
fn fit_line(points: &[Vec3]) -> Line {
    let centroid = points.iter().copied().sum::<Vec3>() / points.len() as f32;
    let mut covariance = Mat3::ZERO;
    for &p in points {
        let d = p - centroid;
        covariance += Mat3::from_cols(d * d.x, d * d.y, d * d.z);
    }
    let along = points[points.len() - 1] - points[0];
    let mut dir = along.normalize_or_zero();
    if dir == Vec3::ZERO {
        dir = Vec3::X;
    }
    for _ in 0..POWER_STEPS {
        let next = (covariance * dir).normalize_or_zero();
        if next == Vec3::ZERO {
            break;
        }
        dir = next;
    }
    if dir.dot(along) < 0.0 {
        dir = -dir;
    }
    Line {
        origin: centroid,
        dir,
    }
}

/// `[from, to)` ranges of `ca` (offset by `at`) that each fit one straight
/// axis, cut where the local axis points leave their line.
fn straight_runs(ca: &[Vec3], at: usize, out: &mut Vec<[usize; 2]>) {
    let n = ca.len();
    if n >= 2 * MIN_RESIDUES {
        let axis = local_axis(ca);
        let line = fit_line(&axis);
        let worst = axis
            .iter()
            .enumerate()
            .map(|(k, &p)| (k, line.deviation(p)))
            .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(Ordering::Equal));
        if let Some((k, _)) = worst.filter(|&(_, d)| d > KINK) {
            // Axis point `k` sits beside residue `k + 1`.
            let cut = (k + 1).clamp(MIN_RESIDUES, n - MIN_RESIDUES);
            straight_runs(&ca[..cut], at, out);
            straight_runs(&ca[cut..], at + cut, out);
            return;
        }
    }
    out.push([at, at + n]);
}

/// Moves each CA of one cylinder (`controls`, every one a helix CA) onto
/// the cylinder's axis.
pub(super) fn straighten(controls: &mut [Vec3]) {
    if controls.len() < MIN_RESIDUES {
        return;
    }
    let line = fit_line(&local_axis(controls));
    for c in controls {
        *c = line.project(*c);
    }
}

/// `[first, end)` slots of each maximal run of helix residues in a segment.
fn helix_runs(shapes: &[Shape], segments: &[[u32; 2]]) -> Vec<[usize; 2]> {
    let mut runs = Vec::new();
    for &[first, end] in segments {
        let (first, end) = (first as usize, end as usize);
        let mut i = first;
        while i < end {
            if shapes[i] != Shape::Helix {
                i += 1;
                continue;
            }
            let start = i;
            while i < end && shapes[i] == Shape::Helix {
                i += 1;
            }
            runs.push([start, i]);
        }
    }
    runs
}

impl CartoonPlan {
    /// This plan with every helix drawn as straight cylinders along its
    /// axis (helices of fewer than four residues become coil), its
    /// sections `samples_per_residue` to a residue step. `frame` is any
    /// frame of the plan's slots, which decides where kinked helices are
    /// cut; use the result's own [`CartoonPlan::frame`] to draw it, since
    /// that is what puts the controls on the axes.
    pub fn with_cylinder_helices(
        &self,
        frame: &CartoonFrame,
        samples_per_residue: usize,
    ) -> CartoonPlan {
        let mut shapes = self.shapes.clone();
        let mut cylinders = Vec::new();
        for [first, end] in helix_runs(&self.shapes, &self.segments) {
            if end - first < MIN_RESIDUES {
                shapes[first..end].fill(Shape::Coil);
                continue;
            }
            shapes[first..end].fill(Shape::Cylinder);
            let mut runs = Vec::new();
            straight_runs(&frame.controls[first..end], first, &mut runs);
            cylinders.extend(runs.into_iter().map(|[a, b]| [a as u32, b as u32]));
        }
        let mut out = CartoonPlan {
            shapes,
            cylinders,
            ..self.without_sections()
        };
        for &[first, end] in &self.segments {
            out.plan_segment(first as usize, end as usize, samples_per_residue);
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::builder::{AtomRow, TopologyBuilder};
    use crate::dssp::DsspCode;

    const RADIUS: f32 = 2.3;
    const RISE: f32 = 1.5;
    const TWIST: f32 = 100.0 * std::f32::consts::PI / 180.0;

    /// CAs of an ideal alpha helix along +z.
    fn helix(n: usize) -> Vec<Vec3> {
        (0..n)
            .map(|i| {
                let a = i as f32 * TWIST;
                Vec3::new(RADIUS * a.cos(), RADIUS * a.sin(), RISE * i as f32)
            })
            .collect()
    }

    fn angle_between(a: Vec3, b: Vec3) -> f32 {
        a.normalize()
            .dot(b.normalize())
            .clamp(-1.0, 1.0)
            .acos()
            .to_degrees()
    }

    #[test]
    fn local_axis_points_of_an_ideal_helix_lie_on_its_axis() {
        for p in local_axis(&helix(14)) {
            assert!(p.truncate().length() < 0.05, "{p}");
        }
    }

    #[test]
    fn the_fitted_axis_of_an_ideal_helix_is_z_within_two_degrees() {
        let line = fit_line(&local_axis(&helix(14)));
        assert!(angle_between(line.dir, Vec3::Z) < 2.0, "{}", line.dir);
    }

    #[test]
    fn straightening_puts_every_ca_on_one_line_keeping_its_height() {
        let mut ca = helix(12);
        straighten(&mut ca);
        for (i, p) in ca.iter().enumerate() {
            assert!(p.truncate().length() < 0.1, "{p}");
            assert!((p.z - RISE * i as f32).abs() < 0.1, "{p}");
        }
    }

    /// Two ideal helices of `n` residues, the second turned `degrees`
    /// about the x axis through the first's last axis point.
    fn kinked(n: usize, degrees: f32) -> Vec<Vec3> {
        let first = helix(n);
        let pivot = Vec3::new(0.0, 0.0, RISE * (n - 1) as f32 + RISE);
        let turn = glam::Quat::from_rotation_x(degrees.to_radians());
        let second = helix(n).into_iter().map(|p| pivot + turn * p);
        first.into_iter().chain(second).collect()
    }

    #[test]
    fn a_straight_helix_is_one_run_and_a_kinked_one_is_cut_at_the_kink() {
        let mut runs = Vec::new();
        straight_runs(&helix(20), 0, &mut runs);
        assert_eq!(runs, [[0, 20]]);

        let mut runs = Vec::new();
        straight_runs(&kinked(12, 35.0), 0, &mut runs);
        assert_eq!(runs.len(), 2, "{runs:?}");
        assert!(runs[0][1].abs_diff(12) <= 2, "{runs:?}");
        assert_eq!(runs[1][1], 24);
    }

    fn chain(points: &[Vec3]) -> (crate::Topology, Vec<Vec3>) {
        let mut b = TopologyBuilder::new();
        for (i, &position) in points.iter().enumerate() {
            let mut name = [b' '; 4];
            name[..2].copy_from_slice(b"CA");
            b.push(&AtomRow {
                element: crate::Element::CARBON,
                name,
                serial: i as u32 + 1,
                alt_loc: 0,
                comp: "ALA",
                asym: "A",
                auth_asym: "A",
                seq_id: i as i32 + 1,
                auth_seq_id: i as i32 + 1,
                ins_code: 0,
                entity: 1,
                position,
                occupancy: 1.0,
                b_factor: 0.0,
                charge: 0,
                hetero: false,
            });
        }
        let structure = b.finish().unwrap();
        let positions = structure.frame(0).positions().to_vec();
        (
            std::sync::Arc::try_unwrap(structure.topology).unwrap(),
            positions,
        )
    }

    /// Coil, a helix of ten residues, coil.
    fn helix_in_coil() -> (CartoonPlan, Vec<Vec3>) {
        let mut points = vec![Vec3::new(8.0, 0.0, -3.0), Vec3::new(5.0, 1.0, -1.5)];
        points.extend(helix(10));
        points.extend([Vec3::new(5.0, -1.0, 17.0), Vec3::new(8.0, 0.0, 19.0)]);
        let codes: Vec<DsspCode> = (0..points.len())
            .map(|i| {
                if (2..12).contains(&i) {
                    DsspCode::AlphaHelix
                } else {
                    DsspCode::Coil
                }
            })
            .collect();
        let (topology, positions) = chain(&points);
        (super::super::plan(&topology, &positions, &codes), positions)
    }

    #[test]
    fn a_cylinder_plan_draws_the_helix_as_a_collinear_round_tube() {
        let (plan, positions) = helix_in_coil();
        let ribbon_frame = plan.frame(&positions);
        let styled = plan.with_cylinder_helices(&ribbon_frame, 6);
        assert_eq!(styled.cylinders, [[2, 12]]);

        let mesh = styled.mesh(&styled.frame(&positions));
        let cylinder_sections: Vec<_> = mesh
            .sections
            .iter()
            .filter(|s| s.half_width > RADIUS * 0.99)
            .collect();
        assert!(
            cylinder_sections.len() >= 9 * 6,
            "{}",
            cylinder_sections.len()
        );
        for s in &cylinder_sections {
            assert!((s.half_thickness - s.half_width).abs() < 1e-5);
            assert!((s.roundness - 1.0).abs() < 1e-5);
            // The end steps lean toward the neighbouring coil's control
            // point; everything between is exactly on the axis.
            let off_axis = s.center.truncate().length();
            let inside = (2.0 * RISE..8.0 * RISE).contains(&s.center.z);
            assert!(off_axis < if inside { 0.01 } else { 0.4 }, "{}", s.center);
        }
    }

    #[test]
    fn a_cylinder_ends_in_a_step_not_a_taper() {
        let (plan, positions) = helix_in_coil();
        let styled = plan.with_cylinder_helices(&plan.frame(&positions), 6);
        let widths: Vec<f32> = styled.recipes.iter().map(|r| r.half_width).collect();
        let wide = widths.iter().position(|&w| w > 2.0).unwrap();
        assert!(widths[wide - 1] < 0.5, "a step in: {:?}", &widths[..=wide]);
        let last = widths.iter().rposition(|&w| w > 2.0).unwrap();
        assert!(widths[last + 1] < 0.5, "a step out: {:?}", &widths[last..]);
    }

    #[test]
    fn a_short_helix_stays_coil_and_the_ribbon_plan_is_untouched() {
        let points: Vec<Vec3> = helix(3)
            .into_iter()
            .chain(helix(8).into_iter().map(|p| p + Vec3::new(9.0, 0.0, 0.0)))
            .collect();
        let codes: Vec<DsspCode> = (0..points.len())
            .map(|i| {
                if i < 3 {
                    DsspCode::AlphaHelix
                } else {
                    DsspCode::Coil
                }
            })
            .collect();
        let (topology, positions) = chain(&points);
        let plan = super::super::plan(&topology, &positions, &codes);
        let styled = plan.with_cylinder_helices(&plan.frame(&positions), 6);
        assert!(styled.cylinders.is_empty());
        assert_eq!(styled.recipes.len(), plan.recipes.len());
    }
}
