//! Unobserved residues (`REMARK 465`, `_pdbx_unobs_or_zero_occ_residues`)
//! surface as one mmCIF-shaped annotation category in either format.

use std::path::PathBuf;

fn fixture(path: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../fixtures")
        .join(path)
}

#[test]
fn remark_465_rows_become_the_unobserved_residue_category() {
    let text = std::fs::read_to_string(fixture("glycan/6X3Z.pdb")).unwrap();
    let rows = text
        .lines()
        .skip_while(|l| !l.contains("M RES C SSSEQI"))
        .skip(1)
        .take_while(|l| l.starts_with("REMARK 465"))
        .filter(|l| !l[10..].trim().is_empty())
        .count();
    let structure = vv_io::load(fixture("glycan/6X3Z.pdb")).unwrap();
    let cat = structure
        .topology
        .annotations
        .category("pdbx_unobs_or_zero_occ_residues")
        .expect("category present");
    assert_eq!(cat.rows.len(), rows);
    assert_eq!(cat.get("auth_comp_id", 0), Some("GLN"));
    assert_eq!(cat.get("auth_asym_id", 0), Some("A"));
    assert_eq!(cat.get("auth_seq_id", 0), Some("1"));
}
