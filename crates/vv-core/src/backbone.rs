//! Backbone traces: the geometry behind the "tube" representation
//! (cartoon-lite). No secondary-structure shapes; a real cartoon (helices
//! as ribbons, strands as arrows) is `crate::cartoon`, which needs
//! validated assignment (docs/VALIDATION.md).
//!
//! One trace atom per polymer residue (`CA` for amino acids, `P` for
//! nucleotides, `C4'` when a 5' terminal nucleotide has no phosphate), in
//! chain order, split into segments wherever consecutive trace atoms are
//! too far apart to be bonded through the backbone (a chain break, a
//! missing loop). The Tube representation draws it as a round swept mesh
//! (`crate::cartoon::tube_plan`); `tube_from_trace` samples the same
//! Catmull-Rom spline into points.

use glam::Vec3;

use crate::residue_class::ResidueClass;
use crate::topology::Topology;

/// Tube radius in Angstroms. Thinner than a bond stick's 0.15 would read
/// as wire; 0.3 is the conventional "tube" thickness.
pub const TUBE_RADIUS: f32 = 0.3;
/// Spline samples per residue. Six keeps a helix (about 1.5 A rise per
/// residue) visibly round; more only costs draw calls.
pub const SAMPLES_PER_RESIDUE: usize = 6;

/// Longest CA-CA distance still treated as one continuous chain. Trans
/// peptides are 3.8 A; 4.5 A allows cis peptides and sloppy models but
/// not a missing residue (>= 7 A).
const PROTEIN_BREAK: f32 = 4.5;
/// Longest P-P distance in a continuous nucleic acid strand (~6-7 A).
const NUCLEIC_BREAK: f32 = 8.5;

/// Trace atoms of one polymer, grouped into continuous segments.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Trace {
    /// Atom indices, in chain order; every segment has at least two.
    pub segments: Vec<Vec<u32>>,
}

/// Sampled tube geometry, ready to upload as spheres plus cylinders.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Tube {
    pub positions: Vec<Vec3>,
    /// The trace atom each sample stands for (nearest control point).
    pub source: Vec<u32>,
    /// Consecutive sample pairs within a segment.
    pub bonds: Vec<[u32; 2]>,
}

fn trace_atom(topology: &Topology, residue: usize) -> Option<(u32, ResidueClass)> {
    let rec = &topology.residues[residue];
    let kind = topology.residue_class(residue);
    let wanted: &[&str] = match kind {
        ResidueClass::Protein => &["CA"],
        ResidueClass::Nucleic => &["P", "C4'"],
        _ => return None,
    };
    for name in wanted {
        if let Some(atom) = rec
            .atoms
            .clone()
            .find(|&a| topology.atom_name(a as usize) == *name)
        {
            return Some((atom, kind));
        }
    }
    None
}

/// A cartoon's nucleic-acid ladder: for every nucleotide with a trace
/// atom, the pair (trace atom, the base's Watson-Crick atom) -- N1 for a
/// purine (it has N9), N3 for a pyrimidine -- so the rungs of paired
/// bases meet in the middle.
pub fn nucleic_ladder(topology: &Topology) -> Vec<[u32; 2]> {
    (0..topology.residues.len())
        .filter_map(|r| {
            let (trace, kind) = trace_atom(topology, r)?;
            if kind != ResidueClass::Nucleic {
                return None;
            }
            let find = |name: &str| {
                topology.residues[r]
                    .atoms
                    .clone()
                    .find(|&a| topology.atom_name(a as usize) == name)
            };
            let base = if find("N9").is_some() {
                find("N1")
            } else {
                find("N3")
            }?;
            Some([trace, base])
        })
        .collect()
}

/// The backbone trace of every chain, split at breaks. Alternate
/// conformations are not special-cased: the first atom named `CA` (or
/// `P`) in the residue wins, which is the primary alt-loc in practice.
pub fn trace(topology: &Topology, positions: &[Vec3]) -> Trace {
    let mut segments = Vec::new();
    for chain in &topology.chains {
        let mut current: Vec<u32> = Vec::new();
        let mut last: Option<(Vec3, ResidueClass)> = None;
        let flush = |current: &mut Vec<u32>, segments: &mut Vec<Vec<u32>>| {
            if current.len() >= 2 {
                segments.push(std::mem::take(current));
            } else {
                current.clear();
            }
        };
        for residue in chain.residues.clone() {
            let Some((atom, kind)) = trace_atom(topology, residue as usize) else {
                // A ligand or water between polymer residues does not
                // break the trace by itself; distance decides below.
                continue;
            };
            let p = positions[atom as usize];
            if let Some((prev, prev_kind)) = last {
                let limit = match (prev_kind, kind) {
                    (ResidueClass::Protein, ResidueClass::Protein) => PROTEIN_BREAK,
                    (ResidueClass::Nucleic, ResidueClass::Nucleic) => NUCLEIC_BREAK,
                    _ => 0.0, // protein -> nucleic is always a new segment
                };
                if prev.distance(p) > limit {
                    flush(&mut current, &mut segments);
                }
            }
            current.push(atom);
            last = Some((p, kind));
        }
        flush(&mut current, &mut segments);
    }
    Trace { segments }
}

/// Centripetal-free (uniform) Catmull-Rom between `p1` and `p2` at `t`.
pub(crate) fn catmull_rom(p0: Vec3, p1: Vec3, p2: Vec3, p3: Vec3, t: f32) -> Vec3 {
    let t2 = t * t;
    let t3 = t2 * t;
    0.5 * ((2.0 * p1)
        + (-p0 + p2) * t
        + (2.0 * p0 - 5.0 * p1 + 4.0 * p2 - p3) * t2
        + (-p0 + 3.0 * p1 - 3.0 * p2 + p3) * t3)
}

/// Samples a smooth tube through `trace`. Passes exactly through every
/// trace atom (Catmull-Rom interpolates its control points).
pub fn tube_from_trace(trace: &Trace, positions: &[Vec3]) -> Tube {
    let mut tube = Tube::default();
    for segment in &trace.segments {
        let pts: Vec<Vec3> = segment.iter().map(|&a| positions[a as usize]).collect();
        let n = pts.len();
        let at = |i: isize| pts[i.clamp(0, n as isize - 1) as usize];
        let first_sample = tube.positions.len() as u32;
        for i in 0..n - 1 {
            let (p0, p1, p2, p3) = (at(i as isize - 1), pts[i], pts[i + 1], at(i as isize + 2));
            for s in 0..SAMPLES_PER_RESIDUE {
                let t = s as f32 / SAMPLES_PER_RESIDUE as f32;
                tube.positions.push(catmull_rom(p0, p1, p2, p3, t));
                tube.source
                    .push(if t < 0.5 { segment[i] } else { segment[i + 1] });
            }
        }
        tube.positions.push(pts[n - 1]);
        tube.source.push(segment[n - 1]);
        let last_sample = tube.positions.len() as u32 - 1;
        for k in first_sample..last_sample {
            tube.bonds.push([k, k + 1]);
        }
    }
    tube
}

/// `trace` + `tube_from_trace` in one call.
pub fn tube(topology: &Topology, positions: &[Vec3]) -> Tube {
    tube_from_trace(&trace(topology, positions), positions)
}

/// Putty tube radius, by B-factor:
/// linear from `radius_range.0` at the lowest B-factor of the drawn atoms
/// to `radius_range.1` at the highest, clamped to that span. A degenerate
/// range (every drawn atom has the same B-factor) draws at the range's
/// midpoint rather than dividing by zero.
pub fn putty_radius(b_factor: f32, b_range: (f32, f32), radius_range: (f32, f32)) -> f32 {
    let (b_lo, b_hi) = b_range;
    let (r_lo, r_hi) = radius_range;
    let span = b_hi - b_lo;
    let t = if span > f32::EPSILON {
        ((b_factor - b_lo) / span).clamp(0.0, 1.0)
    } else {
        0.5
    };
    r_lo + t * (r_hi - r_lo)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::builder::{AtomRow, TopologyBuilder};

    /// A straight "protein" of `n` alanines with only CA atoms, `step` A
    /// apart, optionally with a 20 A gap before residue `gap_before`.
    fn chain(n: usize, step: f32, gap_before: Option<usize>) -> (Topology, Vec<Vec3>) {
        let mut b = TopologyBuilder::new();
        let mut x = 0.0;
        for i in 0..n {
            if gap_before == Some(i) {
                x += 20.0;
            }
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
                position: Vec3::new(x, 0.0, 0.0),
                occupancy: 1.0,
                b_factor: 0.0,
                charge: 0,
                hetero: false,
            });
            x += step;
        }
        let structure = b.finish().unwrap();
        let positions = structure.frame(0).positions().to_vec();
        (
            std::sync::Arc::try_unwrap(structure.topology).unwrap(),
            positions,
        )
    }

    /// One purine (DG) and one pyrimidine (DC), a few named atoms each.
    #[test]
    fn nucleic_ladder_runs_from_the_trace_atom_to_the_base_pairing_atom() {
        let mut b = TopologyBuilder::new();
        let rows: [(&str, i32, &str); 8] = [
            ("DG", 1, "P"),
            ("DG", 1, "C4'"),
            ("DG", 1, "N9"),
            ("DG", 1, "N1"),
            ("DC", 2, "P"),
            ("DC", 2, "N1"),
            ("DC", 2, "N3"),
            ("HOH", 3, "O"),
        ];
        for (i, (comp, seq, atom)) in rows.iter().enumerate() {
            let mut name = [b' '; 4];
            name[..atom.len()].copy_from_slice(atom.as_bytes());
            b.push(&AtomRow {
                element: crate::Element::CARBON,
                name,
                serial: i as u32 + 1,
                alt_loc: 0,
                comp,
                asym: "A",
                auth_asym: "A",
                seq_id: *seq,
                auth_seq_id: *seq,
                ins_code: 0,
                entity: 1,
                position: Vec3::new(i as f32, 0.0, 0.0),
                occupancy: 1.0,
                b_factor: 0.0,
                charge: 0,
                hetero: false,
            });
        }
        let structure = b.finish().unwrap();
        let t = &structure.topology;
        let named = |i: u32| t.atom_name(i as usize).to_string();
        let ladder = nucleic_ladder(t);
        assert_eq!(ladder.len(), 2);
        assert_eq!(
            (named(ladder[0][0]), named(ladder[0][1])),
            ("P".into(), "N1".into())
        );
        assert_eq!(
            (named(ladder[1][0]), named(ladder[1][1])),
            ("P".into(), "N3".into())
        );
    }

    #[test]
    fn continuous_chain_is_one_segment_and_the_tube_passes_through_every_ca() {
        let (t, p) = chain(5, 3.8, None);
        let tr = trace(&t, &p);
        assert_eq!(tr.segments, vec![vec![0, 1, 2, 3, 4]]);
        let tube = tube_from_trace(&tr, &p);
        assert_eq!(tube.positions.len(), 4 * SAMPLES_PER_RESIDUE + 1);
        assert_eq!(tube.bonds.len(), tube.positions.len() - 1);
        for (i, &atom) in tr.segments[0].iter().enumerate() {
            let sample = tube.positions[i * SAMPLES_PER_RESIDUE];
            assert!(sample.distance(p[atom as usize]) < 1e-5, "knot {i}");
        }
        // Every sample maps to a real atom of the segment.
        assert!(tube.source.iter().all(|&a| (a as usize) < 5));
        assert_eq!(tube.source[0], 0);
        assert_eq!(*tube.source.last().unwrap(), 4);
    }

    #[test]
    fn a_gap_splits_the_trace_and_a_lone_residue_is_dropped() {
        let (t, p) = chain(6, 3.8, Some(3));
        let tr = trace(&t, &p);
        assert_eq!(tr.segments, vec![vec![0, 1, 2], vec![3, 4, 5]]);
        let tube = tube_from_trace(&tr, &p);
        // No bond crosses the gap.
        let gap = tube
            .bonds
            .iter()
            .filter(|[a, b]| {
                tube.positions[*a as usize].distance(tube.positions[*b as usize]) > 5.0
            })
            .count();
        assert_eq!(gap, 0);

        let (t, p) = chain(3, 3.8, Some(2));
        assert_eq!(
            trace(&t, &p).segments,
            vec![vec![0, 1]],
            "a single residue is no tube"
        );
    }

    #[test]
    fn spline_stays_near_a_gentle_curve() {
        // Points on a circle of radius 10: the interpolant between knots
        // should stay within a small fraction of an Angstrom of the arc.
        let (t, _) = chain(8, 3.8, None);
        let p: Vec<Vec3> = (0..8)
            .map(|i| {
                let a = i as f32 * 0.38;
                Vec3::new(10.0 * a.cos(), 10.0 * a.sin(), 0.0)
            })
            .collect();
        let tube = tube(&t, &p);
        for q in &tube.positions {
            assert!((q.length() - 10.0).abs() < 0.15, "{q}");
        }
    }

    fn close(a: f32, b: f32) -> bool {
        (a - b).abs() < 1e-6
    }

    #[test]
    fn putty_radius_spans_the_radius_range_over_the_b_factor_range() {
        assert!(close(putty_radius(10.0, (10.0, 50.0), (0.1, 1.0)), 0.1));
        assert!(close(putty_radius(50.0, (10.0, 50.0), (0.1, 1.0)), 1.0));
        assert!(close(putty_radius(30.0, (10.0, 50.0), (0.1, 1.0)), 0.55));
    }

    #[test]
    fn putty_radius_clamps_outside_the_b_factor_range() {
        assert!(close(putty_radius(0.0, (10.0, 50.0), (0.1, 1.0)), 0.1));
        assert!(close(putty_radius(100.0, (10.0, 50.0), (0.1, 1.0)), 1.0));
    }

    #[test]
    fn putty_radius_is_the_midpoint_when_every_drawn_atom_has_the_same_b_factor() {
        assert!(close(putty_radius(20.0, (20.0, 20.0), (0.1, 1.0)), 0.55));
    }
}
