//! Covalent links between a protein residue and something else: disulfide
//! bridges and attached glycans, read off the bond table.

use crate::{BondTable, Element, ResidueClass, Topology};

fn class_of(top: &Topology, atom: u32) -> ResidueClass {
    top.residue_class(top.residue_index[atom as usize] as usize)
}

/// Disulfide bridges as residue index pairs `[a, b]` with `a < b`, sorted.
/// Inter-chain bridges count.
pub fn disulfides(top: &Topology, bonds: &BondTable) -> Vec<[u32; 2]> {
    let mut pairs: Vec<[u32; 2]> = bonds
        .pairs
        .iter()
        .filter(|&&[a, b]| {
            top.element[a as usize] == Element::SULFUR
                && top.element[b as usize] == Element::SULFUR
                && class_of(top, a) == ResidueClass::Protein
                && class_of(top, b) == ResidueClass::Protein
        })
        .map(|&[a, b]| [top.residue_index[a as usize], top.residue_index[b as usize]])
        .filter(|[ra, rb]| ra != rb)
        .map(|[ra, rb]| [ra.min(rb), ra.max(rb)])
        .collect();
    pairs.sort_unstable();
    pairs.dedup();
    pairs
}

/// A protein residue with a glycan residue bonded to it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Glycosylation {
    pub residue: u32,
    pub glycan: u32,
    /// The protein atom the glycan hangs from (`ND2`, `OG`, ...).
    pub atom: u32,
}

/// Every protein-to-glycan bond, sorted by protein residue.
pub fn glycosylated(top: &Topology, bonds: &BondTable) -> Vec<Glycosylation> {
    let mut found: Vec<Glycosylation> = bonds
        .pairs
        .iter()
        .filter_map(|&[a, b]| {
            let (protein, sugar) = match (class_of(top, a), class_of(top, b)) {
                (ResidueClass::Protein, ResidueClass::Glycan) => (a, b),
                (ResidueClass::Glycan, ResidueClass::Protein) => (b, a),
                _ => return None,
            };
            Some(Glycosylation {
                residue: top.residue_index[protein as usize],
                glycan: top.residue_index[sugar as usize],
                atom: protein,
            })
        })
        .collect();
    found.sort_unstable_by_key(|g| (g.residue, g.glycan));
    found.dedup_by_key(|g| (g.residue, g.glycan));
    found
}
