//! Per-residue molecule class (protein, nucleic, lipid, glycan, water,
//! ion, small molecule, other) and the roles layered on it.
//!
//! Decided once per structure ([`classify`], called from `Structure`
//! construction) in two layers, then cached on the topology as one byte
//! per residue (`Topology::residue_class`):
//!
//! 1. **Names** ([`names`] tables, `crate::glycan`, `crate::element`'s ion
//!    list): one lookup per distinct residue name, then a parallel table
//!    read per residue.
//! 2. **Structure** ([`fallback`]) for names no table knows: element rules
//!    (a lone metal or halide atom is an ion, one O with up to two H is
//!    water), then the bond graph of just those residues and their chain
//!    neighbours: a residue peptide-bonded to protein is protein, one
//!    phosphodiester-bonded to nucleic is nucleic, a long alkyl chain with
//!    an acyl/phosphate head (or a sterol ring system) is lipid; anything
//!    else is a small molecule.
//!
//! A class is chemistry; [`Roles`] say what a residue does in this
//! structure (bound ligand, membrane lipid, additive, cofactor) and
//! overlap: a lipid in a protein pocket is a `lipid` and a `ligand`.
//! Roles are decided after the classes, in [`roles`].

use std::collections::HashMap;
use std::sync::OnceLock;

use glam::Vec3;
use rayon::prelude::*;

use crate::element::ion_element;
use crate::Topology;

mod fallback;
#[cfg(test)]
mod fallback_tests;
mod lipid_shape;
mod names;
mod role_contacts;
mod role_names;
#[cfg(test)]
mod role_tests;
pub mod roles;

pub use roles::Roles;

/// The molecule class of one residue. `repr(u8)`: the topology stores one
/// per residue.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum ResidueClass {
    Protein,
    Nucleic,
    Lipid,
    /// A monosaccharide `crate::glycan` recognizes.
    Glycan,
    Water,
    /// A monatomic ion.
    Ion,
    /// A non-polymer molecule in no other class: a cofactor, a drug, a
    /// buffer component. Being a bound ligand is a role ([`Roles::LIGAND`]),
    /// not a class.
    #[default]
    SmallMolecule,
    /// Placeholders and residues with no known element (dummy atoms).
    Other,
}

impl ResidueClass {
    pub const COUNT: usize = 8;
    pub const ALL: [ResidueClass; Self::COUNT] = [
        Self::Protein,
        Self::Nucleic,
        Self::Lipid,
        Self::Glycan,
        Self::Water,
        Self::Ion,
        Self::SmallMolecule,
        Self::Other,
    ];

    pub fn index(self) -> usize {
        self as usize
    }

    pub fn label(self) -> &'static str {
        match self {
            Self::Protein => "protein",
            Self::Nucleic => "nucleic",
            Self::Lipid => "lipid",
            Self::Glycan => "glycan",
            Self::Water => "water",
            Self::Ion => "ion",
            Self::SmallMolecule => "small molecule",
            Self::Other => "other",
        }
    }
}

/// Residues per class, indexed by [`ResidueClass::index`].
pub type ClassCounts = [u32; ResidueClass::COUNT];

pub fn counts(classes: &[ResidueClass]) -> ClassCounts {
    let mut counts = [0u32; ResidueClass::COUNT];
    for c in classes {
        counts[c.index()] += 1;
    }
    counts
}

/// `"protein 574 · lipid 312 · water 12000"`: the non-empty classes in
/// [`ResidueClass::ALL`] order, or `None` for an unclassified topology.
pub fn describe(counts: &ClassCounts) -> Option<String> {
    let parts: Vec<String> = ResidueClass::ALL
        .iter()
        .filter(|c| counts[c.index()] > 0)
        .map(|c| format!("{} {}", c.label(), counts[c.index()]))
        .collect();
    (!parts.is_empty()).then(|| parts.join(" \u{b7} "))
}

fn name_index() -> &'static HashMap<String, ResidueClass> {
    static INDEX: OnceLock<HashMap<String, ResidueClass>> = OnceLock::new();
    INDEX.get_or_init(|| {
        let mut m = HashMap::new();
        for &(class, groups) in names::GROUPS {
            for name in groups.iter().flat_map(|g| g.iter()) {
                m.entry(name.to_string()).or_insert(class);
            }
        }
        m
    })
}

/// Amber/GROMACS terminal templates prefix `N`/`C` onto a protein code
/// (`NALA`, `CCYX`): a 4-character name whose last 3 characters are one.
/// Tried only after the plain tables miss, so it cannot shadow a real
/// 4-letter code.
fn is_terminal_protein(upper: &str) -> bool {
    upper.len() == 4
        && upper.is_char_boundary(1)
        && matches!(upper.as_bytes()[0], b'N' | b'C')
        && name_index().get(&upper[1..]) == Some(&ResidueClass::Protein)
}

fn glycan_or_lipid(upper: &str, atom_count: usize) -> ResidueClass {
    let lipid = atom_count >= names::LIPID_OVER_GLYCAN_MIN_ATOMS
        && names::LIPID_OVER_GLYCAN.contains(&upper);
    if lipid {
        ResidueClass::Lipid
    } else {
        ResidueClass::Glycan
    }
}

/// Smallest atom count of each size bucket [`of_name`]'s answer can differ
/// across: a single atom, a few, enough to be a lipid tail.
const SIZE_BUCKETS: [usize; 3] = [1, 2, names::LIPID_OVER_GLYCAN_MIN_ATOMS];

fn size_bucket(atom_count: usize) -> usize {
    SIZE_BUCKETS
        .iter()
        .rposition(|&min| atom_count >= min)
        .unwrap_or(0)
}

/// The class a residue name alone decides, or `None` when only structure
/// can. `atom_count` matters where a name has two readings: a monatomic
/// ion code is an ion only for a one-atom residue (`CO` may be carbon
/// monoxide), Amber `LA`/`AR`/`PA` are lipid only above one atom
/// (lanthanum, argon), and `DHA` is a lipid tail only when large.
pub fn of_name(name: &str, atom_count: usize) -> Option<ResidueClass> {
    let upper = name.to_ascii_uppercase();
    if let Some(&class) = name_index().get(&upper) {
        let single_atom_code = class == ResidueClass::Lipid
            && atom_count == 1
            && names::MULTI_ATOM_ONLY.contains(&upper.as_str());
        if !single_atom_code {
            return Some(class);
        }
    }
    if is_terminal_protein(&upper) {
        Some(ResidueClass::Protein)
    } else if atom_count == 1 && ion_element(name).is_some() {
        Some(ResidueClass::Ion)
    } else if crate::glycan::is_glycan_name(name) {
        Some(glycan_or_lipid(&upper, atom_count))
    } else {
        None
    }
}

/// A residue's class from its name and elements only, with no bond graph:
/// what `Topology::residue_class` answers before the topology is
/// classified.
pub(crate) fn classify_alone(t: &Topology, residue: usize) -> ResidueClass {
    let rec = &t.residues[residue];
    of_name(t.residue_name(residue), rec.atoms.len())
        .or_else(|| fallback::by_elements(t, rec))
        .unwrap_or_default()
}

/// One class per residue. `positions` (frame 0) lets the structural layer
/// perceive bonds for the residues no name table knows; an MD topology's
/// own bond list is used instead when it has one, and with neither the
/// unknown residues stay ligands.
pub fn classify(t: &Topology, positions: Option<&[Vec3]>) -> Vec<ResidueClass> {
    let by_name: Vec<[Option<ResidueClass>; SIZE_BUCKETS.len()]> = (0..t.names.len())
        .map(|i| {
            let name = t.names.get(crate::InternId(i as u32));
            SIZE_BUCKETS.map(|atoms| of_name(name, atoms))
        })
        .collect();
    let first: Vec<Option<ResidueClass>> = t
        .residues
        .par_iter()
        .map(|r| by_name[r.comp.0 as usize][size_bucket(r.atoms.len())])
        .collect();
    fallback::resolve(t, positions, first)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Codes the glycan table shares with a polymer residue (statine, D-Glu,
    /// the methyladenosines); the polymer reading wins in `of_name`.
    const GLYCAN_CODES_OF_POLYMER_RESIDUES: &[&str] = &["STA", "DGL", "6MA", "1MA"];

    #[test]
    fn tables_are_pairwise_disjoint() {
        let mut seen: HashMap<String, ResidueClass> = HashMap::new();
        let mut clashes = Vec::new();
        for &(class, groups) in names::GROUPS {
            for name in groups.iter().flat_map(|g| g.iter()) {
                let key = name.to_ascii_uppercase();
                if let Some(prev) = seen.insert(key.clone(), class) {
                    clashes.push(format!("{key}: {prev:?} and {class:?}"));
                }
            }
        }
        for name in seen.keys() {
            if ion_element(name).is_some() {
                clashes.push(format!("{name}: also an ion"));
            }
            if crate::glycan::is_glycan_name(name)
                && !GLYCAN_CODES_OF_POLYMER_RESIDUES.contains(&name.as_str())
            {
                clashes.push(format!("{name}: also a glycan"));
            }
        }
        assert!(clashes.is_empty(), "{clashes:?}");
    }

    #[test]
    fn names_decide_each_class() {
        let class = |n: &str| of_name(n, 10);
        assert_eq!(class("ALA"), Some(ResidueClass::Protein));
        assert_eq!(class("hse"), Some(ResidueClass::Protein));
        assert_eq!(class("NALA"), Some(ResidueClass::Protein));
        assert_eq!(class("CCYX"), Some(ResidueClass::Protein));
        assert_eq!(class("DA5"), Some(ResidueClass::Nucleic));
        assert_eq!(class("TIP3"), Some(ResidueClass::Water));
        assert_eq!(class("W"), Some(ResidueClass::Water));
        assert_eq!(of_name("SOD", 1), Some(ResidueClass::Ion));
        assert_eq!(of_name("Cl-", 1), Some(ResidueClass::Ion));
        assert_eq!(of_name("CO", 2), None, "a multi-atom CO is not cobalt");
        assert_eq!(class("NAG"), Some(ResidueClass::Glycan));
        assert_eq!(class("POPC"), Some(ResidueClass::Lipid));
        assert_eq!(class("CHL1"), Some(ResidueClass::Lipid));
        assert_eq!(class("OLC"), Some(ResidueClass::Lipid));
        assert_eq!(class("HEM"), None);
        assert_eq!(class("NME"), None);
        assert_eq!(class("ZZZ9"), None);
    }

    #[test]
    fn amber_split_lipid_pieces_all_classify() {
        for piece in [
            "PC", "PE", "PS", "PGR", "PA", "OL", "ST", "MY", "LA", "AR", "DHA",
        ] {
            assert_eq!(of_name(piece, 40), Some(ResidueClass::Lipid), "{piece}");
        }
    }

    #[test]
    fn one_atom_element_codes_are_not_lipid_tails() {
        assert_eq!(of_name("PA", 30), Some(ResidueClass::Lipid));
        assert_eq!(of_name("PA", 1), None);
        assert_eq!(of_name("AR", 1), None);
        assert_eq!(of_name("LA", 1), None);
    }

    #[test]
    fn describe_lists_non_empty_classes_in_order() {
        let mut c = [0u32; ResidueClass::COUNT];
        assert_eq!(describe(&c), None);
        c[ResidueClass::Water.index()] = 12;
        c[ResidueClass::Protein.index()] = 3;
        assert_eq!(describe(&c).unwrap(), "protein 3 \u{b7} water 12");
    }

    #[test]
    fn counts_add_up() {
        let c = counts(&[
            ResidueClass::Protein,
            ResidueClass::Protein,
            ResidueClass::Water,
        ]);
        assert_eq!(c[ResidueClass::Protein.index()], 2);
        assert_eq!(c[ResidueClass::Water.index()], 1);
        assert_eq!(c.iter().sum::<u32>(), 3);
    }
}
