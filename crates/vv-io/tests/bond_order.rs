//! End-to-end bond-order resolution on a real structure: the file's own
//! `_chem_comp_bond` (source #1, highest priority) and the built-in
//! residue templates (source #2), merged by `vv_core::bonds::perceive`.

use std::path::PathBuf;

use vv_core::{BondOrder, Topology};

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../fixtures/small")
        .join(name)
}

/// The first atom named `atom_name` in the first residue named `resname`.
fn find_atom(topology: &Topology, resname: &str, atom_name: &str) -> u32 {
    let (_, res) = topology
        .residues
        .iter()
        .enumerate()
        .find(|(i, _)| topology.residue_name(*i) == resname)
        .unwrap_or_else(|| panic!("no {resname} residue"));
    res.atoms
        .clone()
        .find(|&a| topology.atom_name(a as usize) == atom_name)
        .unwrap_or_else(|| panic!("{resname} has no atom {atom_name}"))
}

fn order_between(bonds: &vv_core::BondTable, a: u32, b: u32) -> BondOrder {
    let pair = if a < b { [a, b] } else { [b, a] };
    let bond = bonds
        .pairs
        .binary_search(&pair)
        .unwrap_or_else(|_| panic!("{a}-{b} not bonded"));
    bonds.order_of(bond as u32)
}

#[test]
fn full_mmcif_chem_comp_bond_and_templates_agree_on_4hhb() {
    let structure = vv_io::load(fixture("4HHB.cif")).unwrap();
    let t = &structure.topology;
    let bonds = vv_core::bonds::perceive(t, structure.frame(0).positions());

    // ALA backbone carbonyl: the file's own `_chem_comp_bond` names it
    // (`ALA C O doub`), and so does the built-in template -- both agree.
    let ala_c = find_atom(t, "ALA", "C");
    let ala_o = find_atom(t, "ALA", "O");
    assert_eq!(
        order_between(&bonds, ala_c, ala_o),
        BondOrder::Double,
        "ALA C=O"
    );

    // HEM has no built-in template order beyond its Fe-N coordination
    // bonds (`METAL_COFACTORS`, always `Single`): its propionate's one
    // double bond comes entirely from the file's own chem_comp_bond
    // (`HEM CGD O1D doub`, `HEM CGD O2D sing`).
    let cgd = find_atom(t, "HEM", "CGD");
    let o1d = find_atom(t, "HEM", "O1D");
    let o2d = find_atom(t, "HEM", "O2D");
    assert_eq!(
        order_between(&bonds, cgd, o1d),
        BondOrder::Double,
        "HEM CGD=O1D"
    );
    assert_eq!(
        order_between(&bonds, cgd, o2d),
        BondOrder::Single,
        "HEM CGD-O2D"
    );
}

/// A structure whose file carries no `_chem_comp_bond` still gets orders
/// from the built-in residue templates alone (source #2): a plain legacy
/// PDB has no such category at all.
#[test]
fn pdb_load_with_no_chem_comp_bond_still_uses_the_built_in_template() {
    let structure = vv_io::load(fixture("1CRN.pdb")).unwrap();
    let t = &structure.topology;
    let bonds = vv_core::bonds::perceive(t, structure.frame(0).positions());
    let c = find_atom(t, "THR", "C");
    let o = find_atom(t, "THR", "O");
    assert_eq!(order_between(&bonds, c, o), BondOrder::Double, "THR C=O");
}
