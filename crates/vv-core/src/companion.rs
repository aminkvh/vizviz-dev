//! What a cartoon or tube draws beside the polymer: the non-polymer
//! content a reader expects to see with it, decided per residue from its
//! class and roles (`crate::residue_class`).

use crate::residue_class::{ResidueClass, Roles};
use crate::topology::Topology;

/// A kind of non-polymer content, each with its own default treatment
/// (docs/RENDERING.md).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Companion {
    /// A bound small molecule or cofactor.
    Ligand,
    /// A crystallization additive or buffer component.
    Additive,
    Ion,
    /// A monosaccharide, attached to the polymer or free.
    Glycan,
    /// A lipid that is not part of a membrane.
    BoundLipid,
    MembraneLipid,
    Water,
}

impl Companion {
    /// The companion a residue of `class` with `roles` is, if any:
    /// polymers and placeholder residues are not.
    pub fn of(class: ResidueClass, roles: Roles) -> Option<Companion> {
        match class {
            ResidueClass::Protein | ResidueClass::Nucleic | ResidueClass::Other => None,
            ResidueClass::Water => Some(Companion::Water),
            ResidueClass::Ion => Some(Companion::Ion),
            ResidueClass::Glycan => Some(Companion::Glycan),
            ResidueClass::Lipid if roles.contains(Roles::MEMBRANE) => {
                Some(Companion::MembraneLipid)
            }
            ResidueClass::Lipid => Some(Companion::BoundLipid),
            ResidueClass::SmallMolecule if roles.contains(Roles::ADDITIVE) => {
                Some(Companion::Additive)
            }
            ResidueClass::SmallMolecule => Some(Companion::Ligand),
        }
    }
}

/// Each residue's companion kind, `None` for polymer residues.
pub fn companions(topology: &Topology) -> Vec<Option<Companion>> {
    (0..topology.residues.len())
        .map(|r| Companion::of(topology.residue_class(r), topology.residue_roles(r)))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn each_class_and_role_maps_to_its_treatment() {
        let of = Companion::of;
        assert_eq!(of(ResidueClass::Protein, Roles::NONE), None);
        assert_eq!(of(ResidueClass::Nucleic, Roles::NONE), None);
        assert_eq!(
            of(ResidueClass::SmallMolecule, Roles::LIGAND | Roles::COFACTOR),
            Some(Companion::Ligand)
        );
        assert_eq!(
            of(ResidueClass::SmallMolecule, Roles::ADDITIVE),
            Some(Companion::Additive)
        );
        assert_eq!(of(ResidueClass::Ion, Roles::NONE), Some(Companion::Ion));
        assert_eq!(of(ResidueClass::Water, Roles::NONE), Some(Companion::Water));
        assert_eq!(
            of(ResidueClass::Glycan, Roles::NONE),
            Some(Companion::Glycan)
        );
        assert_eq!(
            of(ResidueClass::Lipid, Roles::MEMBRANE),
            Some(Companion::MembraneLipid)
        );
        assert_eq!(
            of(ResidueClass::Lipid, Roles::LIGAND),
            Some(Companion::BoundLipid)
        );
    }
}
