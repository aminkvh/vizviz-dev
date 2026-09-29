//! Roles: what a residue does in this structure, layered on its chemical
//! class ([`ResidueClass`]). A residue can hold several, so a fatty acid
//! in a protein pocket is class `lipid` with role `ligand`.
//!
//! - `ligand`: a non-polymer, non-solvent molecule that is not a
//!   crystallization additive: every small molecule, every lipid that is
//!   not part of a membrane, every glycan not covalently attached to a
//!   protein or nucleic acid. A ligand bonded to the polymer through a
//!   non-peptide bond is still a ligand.
//! - `membrane`: a lipid in a contact network of at least
//!   [`role_contacts::MEMBRANE_MIN_LIPIDS`] lipids (bilayer, micelle).
//! - `additive` and `cofactor`: small molecules named in [`role_names`].
//!
//! Decided once, right after the classes; cost is a name test per small
//! molecule plus one grid pass over lipid atoms and one over glycan atoms.

use std::ops::BitOr;

use glam::Vec3;
use rayon::prelude::*;

use super::{role_contacts, role_names, ResidueClass};
use crate::Topology;

/// A set of roles, one byte per residue.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Roles(u8);

impl Roles {
    pub const NONE: Roles = Roles(0);
    pub const LIGAND: Roles = Roles(1);
    pub const MEMBRANE: Roles = Roles(2);
    pub const ADDITIVE: Roles = Roles(4);
    pub const COFACTOR: Roles = Roles(8);

    pub const fn contains(self, other: Roles) -> bool {
        self.0 & other.0 == other.0
    }
}

impl BitOr for Roles {
    type Output = Roles;

    fn bitor(self, rhs: Roles) -> Roles {
        Roles(self.0 | rhs.0)
    }
}

/// Roles a small molecule's name decides; none for any other class.
fn small_molecule_roles(class: ResidueClass, name: &str) -> Roles {
    if class != ResidueClass::SmallMolecule {
        Roles::NONE
    } else if role_names::is_additive(name) {
        Roles::ADDITIVE
    } else if role_names::is_cofactor(name) {
        Roles::LIGAND | Roles::COFACTOR
    } else {
        Roles::LIGAND
    }
}

/// A residue's roles from its class and name alone, for a topology whose
/// roles were never assigned (no coordinates to test contacts with): every
/// lipid counts as membrane and no glycan as a ligand.
pub(crate) fn alone(class: ResidueClass, name: &str) -> Roles {
    match class {
        ResidueClass::Lipid => Roles::MEMBRANE,
        other => small_molecule_roles(other, name),
    }
}

/// One [`Roles`] per residue. `positions` are frame-0 coordinates.
pub(crate) fn assign(t: &Topology, classes: &[ResidueClass], positions: &[Vec3]) -> Vec<Roles> {
    let mut roles: Vec<Roles> = classes
        .par_iter()
        .enumerate()
        .map(|(r, &class)| small_molecule_roles(class, t.residue_name(r)))
        .collect();
    let lipids = residues_of(classes, ResidueClass::Lipid);
    let membrane = role_contacts::membrane(t, positions, &lipids);
    for (&r, in_membrane) in lipids.iter().zip(membrane) {
        roles[r as usize] = if in_membrane {
            Roles::MEMBRANE
        } else {
            Roles::LIGAND
        };
    }
    let glycans = residues_of(classes, ResidueClass::Glycan);
    let attached = role_contacts::attached_glycans(t, classes, positions, &glycans);
    for (&r, attached) in glycans.iter().zip(attached) {
        if !attached {
            roles[r as usize] = Roles::LIGAND;
        }
    }
    roles
}

fn residues_of(classes: &[ResidueClass], class: ResidueClass) -> Vec<u32> {
    classes
        .par_iter()
        .enumerate()
        .filter_map(|(r, &c)| (c == class).then_some(r as u32))
        .collect()
}
