//! Molecule classification on real and hand-built structures.

use std::path::PathBuf;

use vv_core::{select, ResidueClass, Structure};

fn fixture(path: &str) -> Structure {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../fixtures")
        .join(path);
    vv_io::load(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()))
}

fn classes_of(s: &Structure, resname: &str) -> Vec<ResidueClass> {
    let t = &s.topology;
    (0..t.residue_count())
        .filter(|&r| t.residue_name(r) == resname)
        .map(|r| t.residue_class(r))
        .collect()
}

fn all_are(s: &Structure, resname: &str, class: ResidueClass) {
    let found = classes_of(s, resname);
    assert!(!found.is_empty(), "no {resname} residue in the fixture");
    assert!(found.iter().all(|&c| c == class), "{resname}: {found:?}");
}

fn count(s: &Structure, expr: &str) -> usize {
    select(&s.topology, s.frame(0).positions(), expr)
        .unwrap_or_else(|e| panic!("{expr}: {e}"))
        .count_ones(..)
}

#[test]
fn hemoglobin_heme_is_a_ligand_not_protein() {
    let s = fixture("small/4HHB.cif");
    all_are(&s, "HEM", ResidueClass::SmallMolecule);
    all_are(&s, "PO4", ResidueClass::SmallMolecule);
    all_are(&s, "HOH", ResidueClass::Water);
    all_are(&s, "VAL", ResidueClass::Protein);
    let c = &s.topology.class_counts;
    assert_eq!(c[ResidueClass::SmallMolecule.index()], 6, "4 HEM and 2 PO4");
    assert_eq!(c.iter().sum::<u32>() as usize, s.topology.residue_count());
    assert_eq!(
        count(&s, "ligand or additive"),
        count(&s, "resname HEM PO4")
    );
}

#[test]
fn glycans_are_separate_from_the_protein_they_hang_off() {
    let s = fixture("glycan/6X3Z.pdb");
    all_are(&s, "NAG", ResidueClass::Glycan);
    all_are(&s, "MAN", ResidueClass::Glycan);
    all_are(&s, "ASN", ResidueClass::Protein);
    assert_eq!(s.topology.class_counts[ResidueClass::Glycan.index()], 21);
    assert_eq!(count(&s, "glycan"), count(&s, "resname NAG MAN BMA"));
}

#[test]
fn membrane_fixture_classes() {
    let s = fixture("membrane/membrane.pdb");
    all_are(&s, "POPC", ResidueClass::Lipid);
    all_are(&s, "CHL1", ResidueClass::Lipid);
    all_are(&s, "TIP3", ResidueClass::Water);
    all_are(&s, "SOD", ResidueClass::Ion);
    all_are(&s, "CLA", ResidueClass::Ion);
    for piece in ["PC", "PA", "OL"] {
        all_are(&s, piece, ResidueClass::Lipid);
    }
}

#[test]
fn unknown_names_are_decided_by_structure() {
    let s = fixture("membrane/membrane.pdb");
    all_are(&s, "XYZ", ResidueClass::Protein);
    all_are(&s, "LPX", ResidueClass::Lipid);
    all_are(&s, "SW9", ResidueClass::Water);
    all_are(&s, "PO4", ResidueClass::SmallMolecule);
    all_are(&s, "BZX", ResidueClass::SmallMolecule);
}

#[test]
fn keywords_agree_with_the_classes() {
    let s = fixture("membrane/membrane.pdb");
    let lipid = count(&s, "lipid");
    assert_eq!(
        lipid,
        count(&s, "resname POPC CHL1 PC PA OL LPX"),
        "lipid keyword covers every lipid residue and nothing else"
    );
    assert_eq!(count(&s, "solvent"), count(&s, "water or ion"));
    assert_eq!(count(&s, "polymer"), count(&s, "protein or nucleic"));
    assert_eq!(count(&s, "resname XYZ and protein"), 5);
    let total: usize = [
        "protein",
        "nucleic",
        "lipid",
        "glycan",
        "water",
        "ion",
        "(ligand or additive) and not (lipid or glycan)",
    ]
    .iter()
    .map(|k| count(&s, k))
    .sum();
    assert_eq!(total, count(&s, "all"), "the classes partition the atoms");
}
