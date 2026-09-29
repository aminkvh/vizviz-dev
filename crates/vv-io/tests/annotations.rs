//! Header metadata survives parsing, and both formats land it under the
//! same mmCIF-style names.

use std::path::PathBuf;

use vv_core::{Annotations, Structure};

fn load(name: &str) -> Structure {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../fixtures/small")
        .join(name);
    let s = vv_io::load(path).unwrap_or_else(|e| panic!("{name}: {e}"));
    s.topology.validate().unwrap();
    s
}

fn lower(s: Option<&str>) -> String {
    s.unwrap_or("").to_ascii_lowercase()
}

/// What both 4HHB files must agree on, whichever record they came from.
fn check_hemoglobin(a: &Annotations, source: &str) {
    // The deposited title spells it the British way ("DEOXYHAEMOGLOBIN").
    assert!(
        lower(a.title()).contains("haemoglobin"),
        "{source} title: {:?}",
        a.title()
    );
    assert_eq!(a.method(), Some("X-RAY DIFFRACTION"), "{source}");
    let res = a
        .resolution()
        .unwrap_or_else(|| panic!("{source}: no resolution"));
    assert!((res - 1.74).abs() < 0.01, "{source} resolution {res}");
    assert_eq!(a.keywords(), Some("OXYGEN TRANSPORT"), "{source}");
    assert_eq!(lower(a.organism()), "homo sapiens", "{source}");

    let entities = a.entities();
    let globins = entities
        .iter()
        .filter(|(_, d)| {
            let d = d.to_ascii_lowercase();
            d.contains("hemoglobin") || d.contains("heme")
        })
        .count();
    assert!(globins >= 2, "{source} entities: {entities:?}");

    let unp = a.uniprot_accessions();
    for acc in ["P69905", "P68871"] {
        assert!(unp.iter().any(|u| u == acc), "{source} accessions: {unp:?}");
    }
    assert!(
        a.doi().is_some() || a.citation_title().is_some(),
        "{source}: no citation"
    );
}

#[test]
fn hemoglobin_mmcif_annotations() {
    let s = load("4HHB.cif");
    let a = &s.topology.annotations;
    check_hemoglobin(a, "4HHB.cif");
    assert_eq!(a.get("entry", "id"), Some("4HHB"));
    assert_eq!(a.deposition_date(), Some("1984-03-07"));
    assert_eq!(a.doi(), Some("10.1016/0022-2836(84)90472-8"));
    assert_eq!(a.entities().len(), 5);
    // The refine whitelist keeps the R factors but nothing else.
    let refine = a.category("refine").unwrap();
    assert_eq!(
        refine.items,
        vec!["ls_d_res_high", "ls_R_factor_R_work", "ls_R_factor_R_free"]
    );
    assert_eq!(
        refine.get("ls_R_factor_R_free", 0),
        None,
        "null reads absent"
    );
    assert!(a.get("cell", "length_a").is_some());
    assert!(a.category("pdbx_audit_revision_history").is_some());
    // Existing behaviour is untouched.
    assert_eq!(s.topology.id, "4HHB");
    assert_eq!(s.topology.title, a.title().unwrap());
    assert_eq!(s.atom_count(), 4779);
}

#[test]
fn hemoglobin_pdb_annotations_use_the_same_names() {
    let s = load("4HHB.pdb");
    let a = &s.topology.annotations;
    check_hemoglobin(a, "4HHB.pdb");
    assert_eq!(a.get("entry", "id"), Some("4HHB"));
    // Column split of the HEADER record: date is kept as written.
    assert_eq!(a.deposition_date(), Some("07-MAR-84"));
    assert_eq!(a.doi(), Some("10.1016/0022-2836(84)90472-8"));
    assert!(lower(a.citation_title()).contains("deoxyhaemoglobin"));
    assert_eq!(a.get("struct_keywords", "text"), Some("OXYGEN TRANSPORT"));
    assert_eq!(
        a.entities(),
        vec![
            ("1".to_string(), "HEMOGLOBIN SUBUNIT ALPHA".to_string()),
            ("2".to_string(), "HEMOGLOBIN SUBUNIT BETA".to_string()),
        ]
    );
    assert_eq!(s.topology.id, "4HHB");
    assert_eq!(s.topology.title, a.title().unwrap());
    assert_eq!(s.atom_count(), 4779);
}

#[test]
fn crambin_formats_agree_on_method_and_resolution() {
    let cif = load("1CRN.cif");
    let pdb = load("1CRN.pdb");
    let (c, p) = (&cif.topology.annotations, &pdb.topology.annotations);
    assert_eq!(c.method(), Some("X-RAY DIFFRACTION"));
    assert_eq!(p.method(), c.method());
    assert_eq!(c.resolution(), Some(1.5));
    assert_eq!(p.resolution(), c.resolution());
    assert_eq!(p.uniprot_accessions(), c.uniprot_accessions());
    assert_eq!(p.uniprot_accessions(), vec!["P01542"]);
    assert_eq!(cif.atom_count(), 327);
    assert_eq!(pdb.atom_count(), 327);
}

#[test]
fn every_fixture_keeps_its_atoms_and_gains_annotations() {
    for (name, atoms) in [("1UBQ", 660), ("1AKE", 3816)] {
        for ext in ["cif", "pdb"] {
            let s = load(&format!("{name}.{ext}"));
            assert_eq!(s.atom_count(), atoms, "{name}.{ext}");
            let a = &s.topology.annotations;
            assert!(!a.is_empty(), "{name}.{ext}");
            assert_eq!(a.get("entry", "id"), Some(name), "{name}.{ext}");
            assert!(a.resolution().is_some(), "{name}.{ext}");
        }
    }
}

#[test]
fn malformed_mmcif_header_blocks_are_ignored_not_errors() {
    // A tag with no value, a ragged citation loop, and an empty loop header
    // must not stop the atoms from loading.
    let src = "data_x
_entry.id BAD
_exptl.method
_struct.title 'Still parsed'
loop_
_citation.id
_citation.title
primary 'one' extra
loop_
_audit_author.name
loop_
_atom_site.group_PDB
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM N ALA A 0.0 0.0 0.0
ATOM CA ALA A 1.5 0.0 0.0
";
    let s = vv_io::parse(src.as_bytes(), vv_io::Format::Mmcif).unwrap();
    assert_eq!(s.atom_count(), 2);
    let a = &s.topology.annotations;
    assert_eq!(a.get("entry", "id"), Some("BAD"));
    assert_eq!(a.title(), Some("Still parsed"));
    assert_eq!(a.method(), None);
    assert_eq!(a.citation_title(), Some("one"));
}

#[test]
fn files_without_headers_have_empty_annotations() {
    let bare = b"ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n";
    let s = vv_io::parse(bare, vv_io::Format::Pdb).unwrap();
    assert!(s.topology.annotations.is_empty());
}
