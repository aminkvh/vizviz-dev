//! The deposited sequence each reader attaches to the chain records.

use std::path::Path;

fn load(name: &str) -> std::sync::Arc<vv_core::Topology> {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../fixtures/small")
        .join(name);
    vv_io::load(path).unwrap().topology
}

fn first_chain_sequence(name: &str) -> String {
    let top = load(name);
    top.full_sequence[0].clone()
}

#[test]
fn pdb_and_mmcif_give_the_same_deposited_sequence() {
    let pdb = first_chain_sequence("4HHB.pdb");
    assert_eq!(pdb.len(), 141);
    assert!(pdb.starts_with("VLSPADKTNV"));
    assert_eq!(first_chain_sequence("4HHB.cif"), pdb);
}

#[test]
fn there_is_one_sequence_per_chain_record() {
    for name in ["4HHB.pdb", "4HHB.cif"] {
        let top = load(name);
        assert_eq!(top.full_sequence.len(), top.chains.len(), "{name}");
    }
}
