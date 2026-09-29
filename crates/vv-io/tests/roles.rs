//! Roles (ligand, membrane, additive, cofactor) on real structures.
//!
//! The committed fixtures always run. The entries under `fixtures/large/`
//! (2RH1, 1E7H, 2HNX, 4DKL from the PDB) are downloaded by hand and
//! git-ignored; their tests skip when the file is absent.

use std::path::PathBuf;

use vv_core::{select, ResidueClass, Roles, Structure};

fn fixtures_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../fixtures")
}

fn load(path: &str) -> Structure {
    let path = fixtures_dir().join(path);
    vv_io::load(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()))
}

/// A downloaded entry, or `None` (test skips) when it is not there.
fn load_large(id: &str) -> Option<Structure> {
    let path = format!("large/{id}.pdb");
    fixtures_dir().join(&path).exists().then(|| load(&path))
}

fn named(s: &Structure, resname: &str) -> Vec<(ResidueClass, Roles)> {
    let t = &s.topology;
    let found: Vec<_> = (0..t.residue_count())
        .filter(|&r| t.residue_name(r) == resname)
        .map(|r| (t.residue_class(r), t.residue_roles(r)))
        .collect();
    assert!(!found.is_empty(), "no {resname} residue");
    found
}

fn assert_all(s: &Structure, resname: &str, class: ResidueClass, roles: Roles) {
    let found = named(s, resname);
    assert!(
        found.iter().all(|&(c, r)| c == class && r == roles),
        "{resname}: {found:?}"
    );
}

fn assert_bound_lipid(s: &Structure, resname: &str) {
    assert_all(s, resname, ResidueClass::Lipid, Roles::LIGAND);
}

fn count(s: &Structure, expr: &str) -> usize {
    select(&s.topology, s.frame(0).positions(), expr)
        .unwrap_or_else(|e| panic!("{expr}: {e}"))
        .count_ones(..)
}

#[test]
fn hemoglobin_heme_is_a_cofactor_ligand_and_phosphate_an_additive() {
    let s = load("small/4HHB.pdb");
    assert_all(
        &s,
        "HEM",
        ResidueClass::SmallMolecule,
        Roles::LIGAND | Roles::COFACTOR,
    );
    assert_all(&s, "PO4", ResidueClass::SmallMolecule, Roles::ADDITIVE);
    assert_eq!(count(&s, "resname PO4 and ligand"), 0);
    assert_eq!(count(&s, "cofactor"), count(&s, "resname HEM"));
}

#[test]
fn attached_glycans_are_not_ligands() {
    let s = load("glycan/6X3Z.pdb");
    assert!(count(&s, "glycan") > 0);
    assert_eq!(count(&s, "glycan and ligand"), 0);
}

#[test]
fn a_bilayer_is_membrane_and_a_lone_lipid_is_a_ligand() {
    let s = load("membrane/bilayer.pdb");
    assert_eq!(count(&s, "membrane"), count(&s, "chain L"));
    assert_eq!(count(&s, "membrane and ligand"), 0);
    assert_eq!(count(&s, "lipid and ligand"), count(&s, "chain B"));
    assert_eq!(count(&s, "lipid"), count(&s, "membrane or ligand"));
}

#[test]
fn a_handful_of_lipids_is_no_membrane() {
    let s = load("membrane/membrane.pdb");
    assert_eq!(count(&s, "membrane"), 0);
    assert_eq!(count(&s, "lipid and ligand"), count(&s, "lipid"));
    assert_all(&s, "TIP3", ResidueClass::Water, Roles::NONE);
}

#[test]
fn a_gpcr_crystal_has_ligand_lipids_and_additive_reagents() {
    let Some(s) = load_large("2RH1") else { return };
    for lipid in ["CLR", "PLM"] {
        assert_bound_lipid(&s, lipid);
    }
    assert_all(&s, "SO4", ResidueClass::SmallMolecule, Roles::ADDITIVE);
    assert_all(&s, "GLC", ResidueClass::Glycan, Roles::LIGAND);
    assert_eq!(count(&s, "membrane"), 0);
    assert_eq!(count(&s, "additive and ligand"), 0);
}

#[test]
fn opioid_receptor_crystal_keeps_lipids_and_drops_reagents() {
    let Some(s) = load_large("4DKL") else { return };
    assert_bound_lipid(&s, "CLR");
    // Monoacylglycerol: its two copies classify differently by shape, but
    // either way it is a bound ligand.
    assert!(named(&s, "MPG").iter().all(|&(_, r)| r == Roles::LIGAND));
    for reagent in ["SO4", "1PE"] {
        assert_all(&s, reagent, ResidueClass::SmallMolecule, Roles::ADDITIVE);
    }
}

#[test]
fn fatty_acids_in_binding_pockets_are_lipids_and_ligands() {
    for id in ["1E7H", "2HNX"] {
        let Some(s) = load_large(id) else { continue };
        assert_bound_lipid(&s, "PLM");
        assert_eq!(count(&s, "membrane"), 0, "{id}");
    }
}
