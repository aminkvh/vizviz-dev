//! Residue-level contacts between two groups of residues: what a
//! binding-site or interface annotation is made of.

use std::collections::HashMap;

use glam::Vec3;

use crate::{sasa, Grid, Topology};

/// Residue `a` of the first group touches residue `b` of the second; the
/// closest heavy-atom distance in Å.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ResidueContact {
    pub a: u32,
    pub b: u32,
    pub distance: f32,
}

/// Contacts within `cutoff` Å between the heavy atoms of residues chosen
/// by `in_a` and those chosen by `in_b`, one entry per residue pair,
/// sorted by `(a, b)`. `allowed(a, b)` vetoes pairs, e.g. to keep a chain
/// from contacting itself. Only each residue's first conformer counts.
pub fn residue_contacts(
    top: &Topology,
    positions: &[Vec3],
    cutoff: f32,
    in_a: impl Fn(usize) -> bool,
    in_b: impl Fn(usize) -> bool,
    allowed: impl Fn(u32, u32) -> bool,
) -> Vec<ResidueContact> {
    let conformer = sasa::first_conformer(top);
    let heavy = |atom: usize| conformer[atom] && !top.element[atom].is_hydrogen();
    let members = |pick: &dyn Fn(usize) -> bool| -> Vec<u32> {
        (0..top.atom_count())
            .filter(|&a| heavy(a) && pick(top.residue_index[a] as usize))
            .map(|a| a as u32)
            .collect()
    };
    let (atoms_a, atoms_b) = (members(&in_a), members(&in_b));
    if atoms_a.is_empty() || atoms_b.is_empty() {
        return Vec::new();
    }
    let grid = Grid::build(positions, &atoms_b, cutoff.max(1.0));
    let mut closest: HashMap<(u32, u32), f32> = HashMap::new();
    for &atom in &atoms_a {
        let ra = top.residue_index[atom as usize];
        grid.for_each_within(positions, positions[atom as usize], cutoff, |other, d2| {
            let rb = top.residue_index[other as usize];
            if allowed(ra, rb) {
                let d = closest.entry((ra, rb)).or_insert(f32::MAX);
                *d = d.min(d2);
            }
        });
    }
    let mut out: Vec<ResidueContact> = closest
        .into_iter()
        .map(|((a, b), d2)| ResidueContact {
            a,
            b,
            distance: d2.sqrt(),
        })
        .collect();
    out.sort_unstable_by_key(|c| (c.a, c.b));
    out
}
