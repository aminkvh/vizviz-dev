//! Roles on small hand-built systems.

use glam::Vec3;

use super::role_contacts::MEMBRANE_MIN_LIPIDS;
use super::{ResidueClass, Roles};
use crate::builder::{AtomRow, TopologyBuilder};
use crate::{Element, Structure};

const C: Element = Element::CARBON;

/// `(residue name, element, atom name, position)`: one atom per residue.
type Entry = (&'static str, Element, &'static str, [f32; 3]);

/// One chain "A", one residue per entry.
fn structure(residues: &[Entry]) -> Structure {
    let mut b = TopologyBuilder::new();
    for (i, &(comp, element, name, p)) in residues.iter().enumerate() {
        let mut n = [b' '; 4];
        n[..name.len()].copy_from_slice(name.as_bytes());
        b.push(&AtomRow {
            element,
            name: n,
            serial: 0,
            alt_loc: 0,
            comp,
            asym: "A",
            auth_asym: "A",
            seq_id: i as i32 + 1,
            auth_seq_id: i as i32 + 1,
            ins_code: 0,
            entity: 1,
            position: Vec3::from(p),
            occupancy: 1.0,
            b_factor: 0.0,
            charge: 0,
            hetero: true,
        });
    }
    b.finish().unwrap()
}

fn roles(s: &Structure) -> Vec<Roles> {
    s.topology.residue_roles.clone()
}

/// `n` lipids in a row, `gap` A apart.
fn lipid_row(n: usize, gap: f32) -> Vec<Entry> {
    (0..n)
        .map(|i| ("POPC", C, "C1", [gap * i as f32, 0.0, 0.0]))
        .collect()
}

#[test]
fn lipids_in_a_dense_network_are_membrane_and_isolated_ones_are_ligands() {
    let s = structure(&lipid_row(20, 3.5));
    assert!(roles(&s).iter().all(|&r| r == Roles::MEMBRANE));
    let s = structure(&lipid_row(20, 9.0));
    assert!(roles(&s).iter().all(|&r| r == Roles::LIGAND));
}

#[test]
fn a_few_touching_lipids_are_still_ligands() {
    let s = structure(&lipid_row(MEMBRANE_MIN_LIPIDS - 1, 3.5));
    assert!(roles(&s).iter().all(|&r| r == Roles::LIGAND));
}

#[test]
fn a_bound_lipid_away_from_a_membrane_stays_a_ligand() {
    let mut residues = lipid_row(14, 3.5);
    residues.push(("POPC", C, "C1", [200.0, 0.0, 0.0]));
    let r = roles(&structure(&residues));
    assert_eq!(r[13], Roles::MEMBRANE);
    assert_eq!(r[14], Roles::LIGAND);
}

#[test]
fn additives_are_not_ligands_and_cofactors_are_both() {
    let s = structure(&[
        ("GOL", C, "C1", [0.0, 0.0, 0.0]),
        ("SO4", Element::from_symbol(b"S"), "S", [20.0, 0.0, 0.0]),
        ("HEM", C, "C1", [40.0, 0.0, 0.0]),
        ("QQ7", C, "C1", [60.0, 0.0, 0.0]),
    ]);
    let r = roles(&s);
    assert_eq!(r[0], Roles::ADDITIVE);
    assert_eq!(r[1], Roles::ADDITIVE);
    assert_eq!(r[2], Roles::LIGAND | Roles::COFACTOR);
    assert_eq!(r[3], Roles::LIGAND);
}

#[test]
fn only_glycans_free_of_the_polymer_are_ligands() {
    let s = structure(&[
        ("ASN", Element::NITROGEN, "ND2", [0.0, 0.0, 0.0]),
        ("NAG", C, "C1", [1.45, 0.0, 0.0]),
        ("NAG", C, "C1", [2.9, 0.0, 0.0]),
        ("NAG", C, "C1", [30.0, 0.0, 0.0]),
        ("NAG", C, "C1", [31.4, 0.0, 0.0]),
    ]);
    assert_eq!(s.topology.residue_class[1], ResidueClass::Glycan);
    let r = roles(&s);
    assert_eq!(
        r[..3],
        [Roles::NONE; 3],
        "attached, and the chain it starts"
    );
    assert_eq!(r[3..], [Roles::LIGAND; 2], "a free disaccharide");
}

#[test]
fn polymer_water_and_ions_have_no_role() {
    let s = structure(&[
        ("ALA", C, "CA", [0.0; 3]),
        ("SOD", Element::from_symbol(b"NA"), "SOD", [9.0, 0.0, 0.0]),
    ]);
    assert_eq!(roles(&s), [Roles::NONE; 2]);
}
