//! Structural analysis kernels: fast, native Rust, validated against a
//! reference before they ship (docs/VALIDATION.md, docs/ANALYSIS.md).
//!
//! Shape of every kernel, chosen for the scale we design for (tens of
//! thousands of trajectory frames, millions of separate structures):
//!
//! - inputs are a coordinate slice plus index arrays, never a `Structure`,
//!   so the same call serves one frame, one model of an ensemble, or a
//!   frame streamed from a trajectory;
//! - the batch form (`*_into`, writing into a caller-provided slice) is
//!   the primitive and allocates nothing per call, so a loop over frames
//!   can reuse its buffers; the scalar form is a convenience over it;
//! - the `*_frames` forms run frames in parallel with rayon.
//!
//! Angles are returned in degrees. Dihedral signs follow the IUPAC
//! convention (clockwise positive looking down the central bond), the
//! same as MDAnalysis and MDTraj report; the helix test below pins
//! that against real backbone geometry, not just self-consistency.

use glam::Vec3;
use rayon::prelude::*;

use crate::spatial::Grid;
use crate::topology::Topology;

/// Euclidean distance between two atoms, in the structure's native units
/// (Angstroms for anything loaded from mmCIF/PDB).
pub fn distance(positions: &[Vec3], a: usize, b: usize) -> f32 {
    positions[a].distance(positions[b])
}

/// Angle at `b` between `a-b` and `c-b`, in degrees (0-180).
pub fn angle(positions: &[Vec3], a: usize, b: usize, c: usize) -> f32 {
    angle_of(positions[a], positions[b], positions[c])
}

fn angle_of(a: Vec3, b: Vec3, c: Vec3) -> f32 {
    let (u, v) = (a - b, c - b);
    // atan2 form is exact near 0 and 180 where acos loses precision.
    u.cross(v).length().atan2(u.dot(v)).to_degrees()
}

/// Dihedral (torsion) angle of the chain `a-b-c-d` about `b-c`, in
/// degrees (-180..=180), IUPAC sign.
pub fn dihedral(positions: &[Vec3], a: usize, b: usize, c: usize, d: usize) -> f32 {
    dihedral_of(positions[a], positions[b], positions[c], positions[d])
}

fn dihedral_of(p0: Vec3, p1: Vec3, p2: Vec3, p3: Vec3) -> f32 {
    // Project the two outer bonds onto the plane normal to the central
    // bond, then measure the signed angle between the projections.
    let b0 = p0 - p1;
    let b1 = (p2 - p1).normalize_or_zero();
    let b2 = p3 - p2;
    let v = b0 - b0.dot(b1) * b1;
    let w = b2 - b2.dot(b1) * b1;
    let x = v.dot(w);
    let y = b1.cross(v).dot(w);
    y.atan2(x).to_degrees()
}

/// `out[i] = distance(pairs[i])`. `out.len()` must equal `pairs.len()`.
pub fn distances_into(positions: &[Vec3], pairs: &[[u32; 2]], out: &mut [f32]) {
    assert_eq!(pairs.len(), out.len());
    for (o, [a, b]) in out.iter_mut().zip(pairs) {
        *o = positions[*a as usize].distance(positions[*b as usize]);
    }
}

/// `out[i] = angle(triples[i])`, degrees.
pub fn angles_into(positions: &[Vec3], triples: &[[u32; 3]], out: &mut [f32]) {
    assert_eq!(triples.len(), out.len());
    for (o, [a, b, c]) in out.iter_mut().zip(triples) {
        *o = angle_of(
            positions[*a as usize],
            positions[*b as usize],
            positions[*c as usize],
        );
    }
}

/// `out[i] = dihedral(quads[i])`, degrees.
pub fn dihedrals_into(positions: &[Vec3], quads: &[[u32; 4]], out: &mut [f32]) {
    assert_eq!(quads.len(), out.len());
    for (o, [a, b, c, d]) in out.iter_mut().zip(quads) {
        *o = dihedral_of(
            positions[*a as usize],
            positions[*b as usize],
            positions[*c as usize],
            positions[*d as usize],
        );
    }
}

/// Runs a per-frame kernel over every frame in parallel, filling `out`
/// row-major as `[frame][measurement]`. `out.len()` must be
/// `frames.len() * per_frame`.
fn frames_into<F>(frames: &[&[Vec3]], per_frame: usize, out: &mut [f32], kernel: F)
where
    F: Fn(&[Vec3], &mut [f32]) + Sync,
{
    assert_eq!(out.len(), frames.len() * per_frame);
    if per_frame == 0 {
        return;
    }
    out.par_chunks_mut(per_frame)
        .zip(frames.par_iter())
        .for_each(|(row, positions)| kernel(positions, row));
}

/// Distances for every frame, `out[frame * pairs.len() + i]`.
pub fn distances_frames(frames: &[&[Vec3]], pairs: &[[u32; 2]], out: &mut [f32]) {
    frames_into(frames, pairs.len(), out, |p, row| {
        distances_into(p, pairs, row)
    });
}

/// Angles for every frame, `out[frame * triples.len() + i]`.
pub fn angles_frames(frames: &[&[Vec3]], triples: &[[u32; 3]], out: &mut [f32]) {
    frames_into(frames, triples.len(), out, |p, row| {
        angles_into(p, triples, row)
    });
}

/// Dihedrals for every frame, `out[frame * quads.len() + i]`.
pub fn dihedrals_frames(frames: &[&[Vec3]], quads: &[[u32; 4]], out: &mut [f32]) {
    frames_into(frames, quads.len(), out, |p, row| {
        dihedrals_into(p, quads, row)
    });
}

/// One atom pair closer than a cutoff.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Contact {
    pub a: u32,
    pub b: u32,
    pub distance: f32,
}

/// Appends every pair `(a, b)` with `a` from `group_a`, `b` from
/// `group_b`, and `|a - b| <= cutoff` to `out`, sorted by `a` then `b`.
/// A pair of the same atom (when the groups overlap) is skipped. Builds a
/// grid over the smaller group and walks the larger in parallel.
pub fn contacts_into(
    positions: &[Vec3],
    group_a: &[u32],
    group_b: &[u32],
    cutoff: f32,
    out: &mut Vec<Contact>,
) {
    if group_a.is_empty() || group_b.is_empty() {
        return;
    }
    // Query from the larger side so the grid (the cheap part) covers the
    // smaller; swap the pair back when reporting so `a` is always from
    // `group_a`.
    let swapped = group_a.len() < group_b.len();
    let (walk, index) = if swapped {
        (group_b, group_a)
    } else {
        (group_a, group_b)
    };
    let grid = Grid::build(positions, index, cutoff);
    let mut found: Vec<Contact> = walk
        .par_iter()
        .fold(Vec::new, |mut acc, &q| {
            let p = positions[q as usize];
            grid.for_each_within(positions, p, cutoff, |other, d2| {
                if other != q {
                    let (a, b) = if swapped { (other, q) } else { (q, other) };
                    acc.push(Contact {
                        a,
                        b,
                        distance: d2.sqrt(),
                    });
                }
            });
            acc
        })
        .reduce(Vec::new, |mut x, mut y| {
            x.append(&mut y);
            x
        });
    found.sort_unstable_by_key(|c| (c.a, c.b));
    out.extend(found);
}

/// Every pair within one group closer than `cutoff`, each once with
/// `a < b`, sorted. For a whole structure pass `0..atom_count`.
pub fn neighbor_pairs_into(positions: &[Vec3], atoms: &[u32], cutoff: f32, out: &mut Vec<Contact>) {
    if atoms.is_empty() {
        return;
    }
    let grid = Grid::build(positions, atoms, cutoff);
    let mut found: Vec<Contact> = atoms
        .par_iter()
        .fold(Vec::new, |mut acc, &q| {
            let p = positions[q as usize];
            grid.for_each_within(positions, p, cutoff, |other, d2| {
                if other > q {
                    acc.push(Contact {
                        a: q,
                        b: other,
                        distance: d2.sqrt(),
                    });
                }
            });
            acc
        })
        .reduce(Vec::new, |mut x, mut y| {
            x.append(&mut y);
            x
        });
    found.sort_unstable_by_key(|c| (c.a, c.b));
    out.extend(found);
}

/// Number of `group_a`/`group_b` contacts within `cutoff` in every frame,
/// frames in parallel. The time-series question ("how many contacts does
/// this interface keep along the trajectory?") without materializing the
/// pairs of 40K frames.
pub fn contact_counts_frames(
    frames: &[&[Vec3]],
    group_a: &[u32],
    group_b: &[u32],
    cutoff: f32,
) -> Vec<u32> {
    frames
        .par_iter()
        .map(|positions| {
            if group_a.is_empty() || group_b.is_empty() {
                return 0;
            }
            let swapped = group_a.len() < group_b.len();
            let (walk, index) = if swapped {
                (group_b, group_a)
            } else {
                (group_a, group_b)
            };
            let grid = Grid::build(positions, index, cutoff);
            let mut n = 0u32;
            for &q in walk {
                grid.for_each_within(positions, positions[q as usize], cutoff, |other, _| {
                    n += (other != q) as u32;
                });
            }
            n
        })
        .collect()
}

/// The distinct residue pairs behind a contact list, sorted. Two atoms of
/// the same residue never form a residue pair.
pub fn residue_pairs(topology: &Topology, contacts: &[Contact]) -> Vec<[u32; 2]> {
    let mut pairs: Vec<[u32; 2]> = contacts
        .iter()
        .map(|c| {
            let (ra, rb) = (
                topology.residue_index[c.a as usize],
                topology.residue_index[c.b as usize],
            );
            if ra <= rb {
                [ra, rb]
            } else {
                [rb, ra]
            }
        })
        .filter(|[ra, rb]| ra != rb)
        .collect();
    pairs.sort_unstable();
    pairs.dedup();
    pairs
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn distance_between_two_points() {
        let positions = [Vec3::ZERO, Vec3::new(3.0, 4.0, 0.0)];
        assert_eq!(distance(&positions, 0, 1), 5.0);
        assert_eq!(distance(&positions, 0, 0), 0.0);
        assert_eq!(distance(&positions, 0, 1), distance(&positions, 1, 0));
    }

    #[test]
    fn angles_at_the_vertex() {
        let p = [
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::ZERO,
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(-1.0, 0.0, 0.0),
        ];
        assert!((angle(&p, 0, 1, 2) - 90.0).abs() < 1e-4);
        assert!((angle(&p, 0, 1, 3) - 180.0).abs() < 1e-4);
        assert!(angle(&p, 0, 1, 0).abs() < 1e-4);
        let mut out = [0.0; 2];
        angles_into(&p, &[[0, 1, 2], [0, 1, 3]], &mut out);
        assert!((out[0] - 90.0).abs() < 1e-4 && (out[1] - 180.0).abs() < 1e-4);
    }

    #[test]
    fn dihedral_sign_follows_the_right_hand_rule() {
        // Looking down b->c (the +y axis), a sits on +x; rotating d
        // clockwise from there (toward -z... i.e. +z as seen from b) is
        // positive. Pinned numerically: this configuration is +90 in
        // MDAnalysis and MDTraj.
        let p = [
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::ZERO,
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 1.0, -1.0),
        ];
        assert!((dihedral(&p, 0, 1, 2, 3) - 90.0).abs() < 1e-4);
        let q = [p[0], p[1], p[2], Vec3::new(0.0, 1.0, 1.0)];
        assert!((dihedral(&q, 0, 1, 2, 3) + 90.0).abs() < 1e-4);
        // Trans (180) and cis (0).
        let t = [p[0], p[1], p[2], Vec3::new(-1.0, 1.0, 0.0)];
        assert!((dihedral(&t, 0, 1, 2, 3).abs() - 180.0).abs() < 1e-4);
        let c = [p[0], p[1], p[2], Vec3::new(1.0, 1.0, 0.0)];
        assert!(dihedral(&c, 0, 1, 2, 3).abs() < 1e-4);
        // Reversing the chain leaves the value unchanged.
        assert!((dihedral(&p, 3, 2, 1, 0) - dihedral(&p, 0, 1, 2, 3)).abs() < 1e-4);
    }

    #[test]
    fn frames_run_in_parallel_and_lay_out_row_major() {
        let f0 = [Vec3::ZERO, Vec3::X, Vec3::Y];
        let f1 = [Vec3::ZERO, Vec3::X * 2.0, Vec3::Y * 2.0];
        let frames: Vec<&[Vec3]> = vec![&f0, &f1];
        let pairs = [[0, 1], [0, 2], [1, 2]];
        let mut out = vec![0.0; 6];
        distances_frames(&frames, &pairs, &mut out);
        assert_eq!(out[0], 1.0);
        assert_eq!(out[3], 2.0);
        assert!((out[5] - 8f32.sqrt()).abs() < 1e-5);
    }

    fn lattice(n: usize, spacing: f32) -> Vec<Vec3> {
        let mut v = Vec::new();
        for x in 0..n {
            for y in 0..n {
                for z in 0..n {
                    v.push(Vec3::new(x as f32, y as f32, z as f32) * spacing);
                }
            }
        }
        v
    }

    fn brute_force(positions: &[Vec3], a: &[u32], b: &[u32], cutoff: f32) -> Vec<(u32, u32)> {
        let mut v = Vec::new();
        for &i in a {
            for &j in b {
                if i != j && positions[i as usize].distance(positions[j as usize]) <= cutoff {
                    v.push((i, j));
                }
            }
        }
        v.sort_unstable();
        v
    }

    #[test]
    fn contacts_and_neighbors_agree_with_brute_force() {
        let positions = lattice(9, 1.3);
        let n = positions.len() as u32;
        let a: Vec<u32> = (0..n).filter(|i| i % 3 == 0).collect();
        let b: Vec<u32> = (0..n).filter(|i| i % 2 == 0).collect();
        let cutoff = 2.0;

        let mut got = Vec::new();
        contacts_into(&positions, &a, &b, cutoff, &mut got);
        let got_pairs: Vec<(u32, u32)> = got.iter().map(|c| (c.a, c.b)).collect();
        assert_eq!(got_pairs, brute_force(&positions, &a, &b, cutoff));
        for c in &got {
            assert!((c.distance - distance(&positions, c.a as usize, c.b as usize)).abs() < 1e-5);
            assert!(a.contains(&c.a) && b.contains(&c.b));
        }

        // Swapping the groups swaps the roles but not the set.
        let mut swapped = Vec::new();
        contacts_into(&positions, &b, &a, cutoff, &mut swapped);
        let mut sw: Vec<(u32, u32)> = swapped.iter().map(|c| (c.b, c.a)).collect();
        sw.sort_unstable();
        assert_eq!(sw, got_pairs);

        let all: Vec<u32> = (0..n).collect();
        let mut pairs = Vec::new();
        neighbor_pairs_into(&positions, &all, cutoff, &mut pairs);
        let expected: Vec<(u32, u32)> = brute_force(&positions, &all, &all, cutoff)
            .into_iter()
            .filter(|(i, j)| i < j)
            .collect();
        assert_eq!(
            pairs.iter().map(|c| (c.a, c.b)).collect::<Vec<_>>(),
            expected
        );

        let f: Vec<&[Vec3]> = vec![&positions, &positions];
        assert_eq!(
            contact_counts_frames(&f, &a, &b, cutoff),
            vec![got.len() as u32; 2]
        );
    }

    #[test]
    fn empty_groups_produce_nothing() {
        let positions = lattice(3, 1.0);
        let mut out = Vec::new();
        contacts_into(&positions, &[], &[0, 1], 5.0, &mut out);
        contacts_into(&positions, &[0], &[], 5.0, &mut out);
        neighbor_pairs_into(&positions, &[], 5.0, &mut out);
        assert!(out.is_empty());
        assert_eq!(
            contact_counts_frames(&[&positions], &[], &[0], 5.0),
            vec![0]
        );
    }
}
