//! Intra-residue bond templates: which atom-name pairs are bonded inside
//! one standard residue, independent of geometry.
//!
//! Source: wwPDB Chemical Component Dictionary (CCD),
//! <https://www.wwpdb.org/data/ccd>, `_chem_comp_bond` records fetched
//! from `https://files.rcsb.org/ligands/download/<ID>.cif` and
//! converted by `gen_templates.py` in this directory (not part of the
//! build; rerun it by hand to refresh the table).
//!
//! Force-field protonation-state names (HSD, HID, CYM, ASH, GLH, LYN, ...)
//! are not registered in the CCD -- those 3-letter codes collide with
//! unrelated ligands there (HSD is a cyclohexene-triol, CYM is
//! S-methylcysteine, checked by fetch). Each instead aliases its parent
//! standard residue: the parent's own CCD entry is already drawn as a
//! neutral, fully protonated free amino acid and so already carries every
//! hydrogen a protonation variant needs (HIS has both ND1-HD1 and
//! NE2-HE2; ASP has OD2-HD2; GLU has OE2-HE2; LYS has HZ1-3; CYS has
//! SG-HG) -- whichever ones are actually present in a given residue
//! still match by name.
//!
//! Terminal-atom aliases (`H1`-`H3` for an N-terminus, `OT1`/`OT2` for a
//! C-terminus) are not part of a free amino acid's CCD entry either (it
//! has no chain neighbor to be a terminus of); [`gen_templates.py`] appends
//! them by hand to every amino acid template.

use std::collections::HashMap;
use std::sync::OnceLock;

use super::BondOrder;

pub mod generated {
    include!("templates_generated.rs");
}

/// A residue's template, indexed once (not per residue -- `match_templates`
/// runs this over every residue in a structure, so a linear scan or a
/// map built per call would dominate at scale; see its own doc).
pub struct Template {
    pub bonds: &'static [(&'static str, &'static str, BondOrder)],
    /// Non-hydrogen atom names the template mentions, sorted and
    /// deduplicated, for a cheap coverage pre-check: an actual bond also
    /// needs the right distance, so this is only an upper bound, but a
    /// residue that fails even this can never pass, and skipping it
    /// avoids the expensive pairing pass. Checked by `binary_search`, not
    /// a `HashSet`: hashing a short string costs more than a handful of
    /// `Ord` comparisons against this (at most ~40 names).
    pub heavy_names: Vec<&'static str>,
}

static INDEX: OnceLock<HashMap<&'static str, Template>> = OnceLock::new();

fn index() -> &'static HashMap<&'static str, Template> {
    INDEX.get_or_init(|| {
        generated::TEMPLATES
            .iter()
            .map(|&(name, bonds)| {
                let mut heavy_names: Vec<&str> = bonds
                    .iter()
                    .flat_map(|&(a, b, _)| [a, b])
                    .filter(|n| !n.starts_with('H'))
                    .collect();
                heavy_names.sort_unstable();
                heavy_names.dedup();
                (name, Template { bonds, heavy_names })
            })
            .collect()
    })
}

/// This residue's template, or `None` if it has no template (ligands and
/// other unrecognized residues fall back to geometric perception).
pub fn get(residue_name: &str) -> Option<&'static Template> {
    index().get(residue_name)
}

/// This residue's intra-residue bonds by atom name, or `None` if it has
/// no template.
pub fn lookup(residue_name: &str) -> Option<&'static [(&'static str, &'static str, BondOrder)]> {
    get(residue_name).map(|t| t.bonds)
}

/// Common glycan residue names: their anomeric carbon (`C1`, or `C2` for
/// the sialic acids) forms the glycosidic bond to the next residue.
pub const GLYCANS: &[&str] = &["NAG", "MAN", "BMA", "GAL", "FUC", "SIA", "NAN", "GLC"];

/// The sialic acids link through their ketosidic C2, not C1 (IUPAC
/// carbohydrate nomenclature: C1 is their exocyclic carboxyl carbon).
pub const SIALIC_ACIDS: &[&str] = &["SIA", "NAN"];

/// Cofactor residues whose metal keeps its intra-residue coordination
/// bonds (Fe-N of the porphyrin) from the template; see module doc on
/// [`super`] for why free metal ions do not.
pub const METAL_COFACTORS: &[&str] = &["HEM", "HEC"];

#[cfg(test)]
mod tests {
    use super::*;

    /// `bonds` names `a`-`b` (either direction), regardless of order.
    fn has(bonds: &[(&str, &str, BondOrder)], a: &str, b: &str) -> bool {
        bonds
            .iter()
            .any(|&(x, y, _)| (x, y) == (a, b) || (x, y) == (b, a))
    }

    /// `bonds`' order for `a`-`b` (either direction), or `None` if absent.
    fn order_of(bonds: &[(&str, &str, BondOrder)], a: &str, b: &str) -> Option<BondOrder> {
        bonds
            .iter()
            .find(|&&(x, y, _)| (x, y) == (a, b) || (x, y) == (b, a))
            .map(|&(_, _, o)| o)
    }

    #[test]
    fn every_standard_amino_acid_has_a_backbone() {
        for name in [
            "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS",
            "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
        ] {
            let bonds = lookup(name).unwrap_or_else(|| panic!("no template for {name}"));
            assert!(has(bonds, "N", "CA"), "{name} N-CA");
            assert!(has(bonds, "CA", "C"), "{name} CA-C");
            assert!(has(bonds, "C", "O"), "{name} C-O");
        }
    }

    /// Every amino acid's backbone carbonyl is a double bond (wwPDB CCD
    /// `_chem_comp_bond.value_order`, e.g. `ALA C-O DOUB`).
    #[test]
    fn backbone_carbonyl_is_a_double_bond() {
        for name in ["ALA", "GLY", "PRO", "VAL"] {
            let bonds = lookup(name).unwrap();
            assert_eq!(
                order_of(bonds, "C", "O"),
                Some(BondOrder::Double),
                "{name} C=O"
            );
        }
    }

    /// PHE's ring is Kekulized in the CCD (`CG-CD1 DOUB`, `CG-CD2 SING`,
    /// ...): 3 alternating double bonds around the 6-membered ring, plus
    /// its backbone carbonyl -- 4 doubles total, none elsewhere.
    #[test]
    fn phenylalanine_ring_has_three_kekule_double_bonds() {
        let bonds = lookup("PHE").unwrap();
        let ring_doubles = [("CG", "CD1"), ("CD2", "CE2"), ("CE1", "CZ")];
        for (a, b) in ring_doubles {
            assert_eq!(
                order_of(bonds, a, b),
                Some(BondOrder::Double),
                "PHE {a}={b}"
            );
        }
        let count = bonds
            .iter()
            .filter(|&&(_, _, o)| o == BondOrder::Double)
            .count();
        assert_eq!(count, 4, "PHE: 3 ring doubles + the backbone C=O");
    }

    /// ASP/GLU's carboxylate and ARG's guanidinium each carry one double
    /// bond (the CCD's own Kekulé choice, not a delocalized/aromatic mark
    /// -- see `templates.rs` module doc and `BondOrder::Aromatic`).
    #[test]
    fn carboxylate_and_guanidinium_have_one_double_bond() {
        let asp = lookup("ASP").unwrap();
        assert_eq!(order_of(asp, "CG", "OD1"), Some(BondOrder::Double));
        assert_eq!(order_of(asp, "CG", "OD2"), Some(BondOrder::Single));

        let glu = lookup("GLU").unwrap();
        assert_eq!(order_of(glu, "CD", "OE1"), Some(BondOrder::Double));

        let arg = lookup("ARG").unwrap();
        assert_eq!(order_of(arg, "CZ", "NH1"), Some(BondOrder::Single));
        assert_eq!(order_of(arg, "CZ", "NH2"), Some(BondOrder::Double));
    }

    #[test]
    fn protonation_variants_alias_their_parent() {
        for (variant, parent) in [
            ("HSD", "HIS"),
            ("HSE", "HIS"),
            ("HSP", "HIS"),
            ("HID", "HIS"),
            ("HIE", "HIS"),
            ("HIP", "HIS"),
            ("CYX", "CYS"),
            ("CYM", "CYS"),
            ("ASH", "ASP"),
            ("GLH", "GLU"),
            ("LYN", "LYS"),
        ] {
            assert_eq!(lookup(variant), lookup(parent), "{variant} vs {parent}");
        }
    }

    #[test]
    fn nucleotides_and_water_and_glycans_are_present() {
        for name in ["A", "C", "G", "U", "DA", "DC", "DG", "DT", "HOH"] {
            assert!(lookup(name).is_some(), "{name}");
        }
        for name in GLYCANS {
            assert!(lookup(name).is_some(), "{name}");
        }
    }

    #[test]
    fn hem_and_hec_carry_the_four_porphyrin_fe_n_bonds() {
        for name in METAL_COFACTORS {
            let bonds = lookup(name).unwrap();
            for n in ["NA", "NB", "NC", "ND"] {
                assert!(has(bonds, "FE", n), "{name} FE-{n}");
            }
        }
    }

    #[test]
    fn unknown_residue_has_no_template() {
        assert!(lookup("LIG").is_none());
        assert!(lookup("XYZ").is_none());
    }
}
