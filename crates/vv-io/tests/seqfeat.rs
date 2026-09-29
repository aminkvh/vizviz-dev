//! `vv_core::seqfeat` kernels on real structures.

use std::path::PathBuf;

use vv_core::seqfeat;
use vv_core::{bonds, ResidueClass, Roles, Structure};

fn load(path: &str) -> Structure {
    vv_io::load(
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../fixtures")
            .join(path),
    )
    .unwrap()
}

fn bond_table(s: &Structure) -> vv_core::BondTable {
    bonds::perceive(&s.topology, s.frame(0).positions())
}

fn auth(s: &Structure, residue: u32) -> (String, i32) {
    let top = &s.topology;
    let rec = &top.residues[residue as usize];
    let chain = top.names.get(top.chains[rec.chain as usize].auth_asym);
    (chain.to_string(), rec.auth_seq_id)
}

#[test]
fn crambin_has_its_three_disulfides() {
    let s = load("small/1CRN.pdb");
    let pairs = seqfeat::disulfides(&s.topology, &bond_table(&s));
    let numbers: Vec<[i32; 2]> = pairs
        .iter()
        .map(|p| [auth(&s, p[0]).1, auth(&s, p[1]).1])
        .collect();
    assert_eq!(numbers, [[3, 40], [4, 32], [16, 26]]);
}

#[test]
fn hemoglobin_has_no_disulfides_and_contacts_its_hemes() {
    let s = load("small/4HHB.pdb");
    let top = &s.topology;
    assert!(seqfeat::disulfides(top, &bond_table(&s)).is_empty());
    let positions = s.frame(0).positions().to_vec();
    let heme = |r: usize| top.residue_roles(r).contains(Roles::LIGAND);
    let protein = |r: usize| top.residue_class(r) == ResidueClass::Protein;
    let contacts = seqfeat::residue_contacts(top, &positions, 4.0, protein, heme, |_, _| true);
    let touching: Vec<(String, i32)> = contacts.iter().map(|c| auth(&s, c.a)).collect();
    for chain in ["A", "B", "C", "D"] {
        assert!(
            touching.iter().any(|(c, _)| c == chain),
            "no heme contacts in chain {chain}"
        );
    }
    assert!(touching.contains(&("A".to_string(), 87)), "proximal His87");
    assert!(touching.contains(&("A".to_string(), 58)), "distal His58");
    assert!(contacts.iter().all(|c| c.distance <= 4.0));
}

#[test]
fn hemoglobin_interfaces_join_different_chains_only() {
    let s = load("small/4HHB.pdb");
    let top = &s.topology;
    let positions = s.frame(0).positions().to_vec();
    let protein = |r: usize| top.residue_class(r) == ResidueClass::Protein;
    let other_chain = |a: u32, b: u32| auth(&s, a).0 != auth(&s, b).0;
    let contacts = seqfeat::residue_contacts(top, &positions, 4.5, protein, protein, other_chain);
    assert!(contacts.len() > 100);
    assert!(contacts.iter().all(|c| auth(&s, c.a).0 != auth(&s, c.b).0));
}

#[test]
fn spike_glycans_hang_from_asparagines() {
    let s = load("glycan/6X3Z.pdb");
    let found = seqfeat::glycosylated(&s.topology, &bond_table(&s));
    assert!(!found.is_empty());
    assert!(found
        .iter()
        .all(|g| s.topology.residue_name(g.residue as usize) == "ASN"));
    assert!(found
        .iter()
        .all(|g| s.topology.atom_name(g.atom as usize) == "ND2"));
}

#[test]
fn unobserved_residues_are_placed_in_gaps() {
    let s = load("glycan/6X3Z.pdb");
    let top = &s.topology;
    let chain = top
        .chains
        .iter()
        .position(|c| top.names.get(c.auth_asym) == "A");
    let range = top.chains[chain.unwrap()].residues.clone();
    let gaps = seqfeat::unobserved(top, range.clone());
    assert!(!gaps.is_empty());
    assert!(gaps
        .iter()
        .all(|g| g.before >= range.start && g.before <= range.end));
    let listed: usize = gaps.iter().map(|g| g.residues.len()).sum();
    assert!(listed > 10, "{listed} unobserved residues in chain A");
    let starts_at = &gaps[0];
    assert_eq!(starts_at.residues[0].name, "GLN");
}

#[test]
fn a_fully_modeled_structure_has_no_gaps() {
    let s = load("small/1CRN.pdb");
    let range = 0..s.topology.residue_count() as u32;
    assert!(seqfeat::unobserved(&s.topology, range).is_empty());
}

#[test]
fn hydrophobic_residues_are_more_buried_than_charged_ones() {
    let s = load("small/1AKE.pdb");
    let top = &s.topology;
    let rel = seqfeat::relative_sasa(top, s.frame(0).positions(), 1.4);
    let mean = |names: &[&str]| {
        let v: Vec<f32> = (0..top.residue_count())
            .filter(|&r| names.contains(&top.residue_name(r)))
            .map(|r| rel[r])
            .collect();
        v.iter().sum::<f32>() / v.len() as f32
    };
    let buried = mean(&["LEU", "ILE", "VAL", "PHE"]);
    let exposed = mean(&["ARG", "LYS", "ASP", "GLU"]);
    assert!(buried < 0.2 && exposed > 0.3, "{buried} vs {exposed}");
}
