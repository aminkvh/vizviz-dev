//! The structural layer on small hand-built chains.

use glam::Vec3;

use super::{classify, ResidueClass};
use crate::builder::{AtomRow, TopologyBuilder};
use crate::{Element, Structure};

type Atom = (&'static str, Element, [f32; 3]);

const C: Element = Element::CARBON;
const N: Element = Element::NITROGEN;
const O: Element = Element::OXYGEN;

/// One chain "A", one residue per entry, numbered from 1.
fn structure(residues: &[(&str, &[Atom])]) -> Structure {
    let mut b = TopologyBuilder::new();
    for (i, (comp, atoms)) in residues.iter().enumerate() {
        for (name, element, p) in atoms.iter() {
            let mut n = [b' '; 4];
            n[..name.len()].copy_from_slice(name.as_bytes());
            b.push(&AtomRow {
                element: *element,
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
                position: Vec3::from(*p),
                occupancy: 1.0,
                b_factor: 0.0,
                charge: 0,
                hetero: false,
            });
        }
    }
    b.finish().unwrap()
}

fn classes(s: &Structure) -> Vec<ResidueClass> {
    s.topology.residue_class.clone()
}

/// A backbone stretch, peptide-linked to the one at `k - 1`.
fn backbone(k: f32) -> [Atom; 4] {
    let (x, y) = (3.42 * k, 0.74 * k);
    [
        ("N", N, [x, y, 0.0]),
        ("CA", C, [x + 1.46, y, 0.0]),
        ("C", C, [x + 2.22, y + 1.32, 0.0]),
        ("O", O, [x + 2.22, y + 2.55, 0.0]),
    ]
}

#[test]
fn lone_atoms_are_ions_ligands_or_other() {
    let s = structure(&[
        ("QQ1", &[("X", Element::from_symbol(b"ZN"), [0.0; 3])]),
        (
            "QQ2",
            &[("X", Element::from_symbol(b"CL"), [10.0, 0.0, 0.0])],
        ),
        ("QQ3", &[("X", C, [20.0, 0.0, 0.0])]),
        ("QQ4", &[("X", Element::UNKNOWN, [30.0, 0.0, 0.0])]),
    ]);
    use ResidueClass::*;
    assert_eq!(classes(&s), [Ion, Ion, SmallMolecule, Other]);
}

/// `n` carbons in a line, 1.5 A apart.
fn carbons(n: usize) -> Vec<Atom> {
    (0..n)
        .map(|i| ("C", C, [1.5 * i as f32, 0.0, 0.0]))
        .collect()
}

#[test]
fn a_name_with_two_readings_is_decided_by_residue_size() {
    let (big, small, pair) = (carbons(54), carbons(24), carbons(2));
    let s = structure(&[
        ("DHA", &big),
        ("DHA", &small),
        ("CO", &pair),
        (
            "CO",
            &[("CO", Element::from_symbol(b"CO"), [40.0, 0.0, 0.0])],
        ),
        ("PA", &big),
        (
            "PA",
            &[("PA", Element::from_symbol(b"PA"), [50.0, 0.0, 0.0])],
        ),
    ]);
    use ResidueClass::*;
    assert_eq!(
        classes(&s),
        [Lipid, Glycan, SmallMolecule, Ion, Lipid, SmallMolecule]
    );
}

#[test]
fn a_lone_oxygen_with_up_to_two_hydrogens_is_water() {
    let h = Element::HYDROGEN;
    let s = structure(&[
        (
            "QQ1",
            &[
                ("OW", O, [0.0; 3]),
                ("H1", h, [0.96, 0.0, 0.0]),
                ("H2", h, [-0.24, 0.93, 0.0]),
            ],
        ),
        ("QQ2", &[("O", O, [10.0, 0.0, 0.0])]),
        (
            "QQ3",
            &[("O", O, [20.0, 0.0, 0.0]), ("C", C, [21.4, 0.0, 0.0])],
        ),
    ]);
    use ResidueClass::*;
    assert_eq!(classes(&s), [Water, Water, SmallMolecule]);
}

#[test]
fn peptide_bonded_unknown_residues_are_protein_even_in_a_run() {
    let (a, b, c, d) = (backbone(0.0), backbone(1.0), backbone(2.0), backbone(3.0));
    let s = structure(&[("ALA", &a), ("QQ1", &b), ("QQ2", &c), ("GLY", &d)]);
    assert!(classes(&s).iter().all(|&c| c == ResidueClass::Protein));
}

#[test]
fn an_unlinked_residue_with_backbone_names_stays_a_ligand() {
    let (a, far) = (backbone(0.0), backbone(30.0));
    let s = structure(&[("ALA", &a), ("QQ1", &far)]);
    assert_eq!(
        classes(&s),
        [ResidueClass::Protein, ResidueClass::SmallMolecule]
    );
}

#[test]
fn phosphodiester_bonded_unknown_residues_are_nucleic() {
    let p = Element::from_symbol(b"P");
    let unit = |k: f32| -> [Atom; 2] {
        [
            ("P", p, [3.5 * k, 0.0, 0.0]),
            ("O3'", O, [3.5 * k + 1.9, 0.0, 0.0]),
        ]
    };
    let (a, b, c) = (unit(0.0), unit(1.0), unit(2.0));
    let s = structure(&[("DA", &a), ("QQ1", &b), ("DT", &c)]);
    assert!(classes(&s).iter().all(|&c| c == ResidueClass::Nucleic));
}

#[test]
fn md_bonds_stand_in_for_positions() {
    let (a, b) = (backbone(0.0), backbone(1.0));
    let s = structure(&[("ALA", &a), ("QQ1", &b)]);
    let mut t = s.topology.as_ref().clone();
    t.md_bonds = Some(vec![[2, 4]]); // C of residue 0, N of residue 1
    assert_eq!(
        classify(&t, None),
        [ResidueClass::Protein, ResidueClass::Protein]
    );
    t.md_bonds = None;
    assert_eq!(
        classify(&t, None),
        [ResidueClass::Protein, ResidueClass::SmallMolecule]
    );
}

/// `n` carbons on a 1.53 A zig-zag along z.
fn alkyl(n: usize) -> Vec<Atom> {
    (0..n)
        .map(|i| {
            let x = if i % 2 == 0 { 0.0 } else { 0.88 };
            ("C", C, [x, 0.0, 1.25 * i as f32])
        })
        .collect()
}

/// An acyl carbon one bond below the chain start, with its two oxygens.
fn ester_head() -> Vec<Atom> {
    vec![
        ("CA", C, [0.0, 0.0, -1.25]),
        ("O1", O, [-1.2, 0.0, -1.25]),
        ("O2", O, [0.5, 1.1, -1.95]),
    ]
}

#[test]
fn a_long_alkyl_chain_needs_a_polar_head_to_be_a_lipid() {
    let bare = alkyl(12);
    assert_eq!(
        classes(&structure(&[("QQ1", &bare)])),
        [ResidueClass::SmallMolecule]
    );

    let mut esterified = alkyl(12);
    esterified.extend(ester_head());
    assert_eq!(
        classes(&structure(&[("QQ1", &esterified)])),
        [ResidueClass::Lipid]
    );

    let mut short = alkyl(6);
    short.extend(ester_head());
    assert_eq!(
        classes(&structure(&[("QQ1", &short)])),
        [ResidueClass::SmallMolecule]
    );
}
