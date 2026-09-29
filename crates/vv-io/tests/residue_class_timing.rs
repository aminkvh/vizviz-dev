//! Classification and role time on a synthetic multi-million-atom system:
//! a membrane patch, waters, a protein block with attached glycans.
//! Timing only, so it is ignored by default:
//!
//!     cargo test --release -p vv-io --test residue_class_timing -- --ignored --nocapture

use std::time::Instant;

use vv_core::builder::{AtomRow, TopologyBuilder};
use vv_core::glam::Vec3;
use vv_core::{Element, ResidueClass, Roles, Topology};

const WATERS: usize = 2_000_000;
const LIPIDS: usize = 100_000;
const LIPID_ATOMS: usize = 30;
const UNKNOWN_LIGANDS: usize = 5_000;
const LIGAND_ATOMS: usize = 12;
const PROTEINS: usize = 1_000_000;
const PROTEIN_ATOMS: usize = 5;
const GLYCANS: usize = 2_000;

/// `atoms` carbons stacked along z from `origin`, 1.5 A apart.
fn push_residue(b: &mut TopologyBuilder, comp: &str, seq: i32, atoms: usize, origin: Vec3) {
    for i in 0..atoms {
        let mut name = *b"C   ";
        name[1] = b'A' + (i % 26) as u8;
        b.push(&AtomRow {
            element: Element::CARBON,
            name,
            serial: 0,
            alt_loc: 0,
            comp,
            asym: "A",
            auth_asym: "A",
            seq_id: seq,
            auth_seq_id: seq,
            ins_code: 0,
            entity: 1,
            position: origin + Vec3::new(0.0, 0.0, 1.5 * i as f32),
            occupancy: 1.0,
            b_factor: 0.0,
            charge: 0,
            hetero: true,
        });
    }
}

/// Lipids on a 317-wide square lattice (6 A pitch, all touching), the
/// unknown ligands in a block beside it, the waters on a cubic lattice
/// above.
fn lipid_origin(i: usize) -> Vec3 {
    Vec3::new(3.5 * (i % 317) as f32, 3.5 * (i / 317) as f32, 0.0)
}

fn ligand_origin(i: usize) -> Vec3 {
    Vec3::new(
        15.0 * (i % 100) as f32,
        -100.0 - 15.0 * (i / 100) as f32,
        0.0,
    )
}

/// Protein residues on a cubic lattice far from the rest.
fn protein_origin(i: usize) -> Vec3 {
    Vec3::new(
        4.0 * (i % 100) as f32,
        4.0 * ((i / 100) % 100) as f32,
        -1000.0 - 8.0 * (i / 10_000) as f32,
    )
}

fn water_origin(i: usize) -> Vec3 {
    Vec3::new(
        3.0 * (i % 126) as f32,
        3.0 * ((i / 126) % 126) as f32,
        100.0 + 3.0 * (i / 15_876) as f32,
    )
}

fn system(unknown: usize) -> TopologyBuilder {
    let mut b = TopologyBuilder::with_capacity(WATERS * 3 + LIPIDS * LIPID_ATOMS);
    let mut seq = 0;
    let mut next = || {
        seq += 1;
        seq
    };
    for i in 0..LIPIDS {
        push_residue(&mut b, "POPC", next(), LIPID_ATOMS, lipid_origin(i));
        if i < unknown {
            push_residue(&mut b, "QQ1", next(), LIGAND_ATOMS, ligand_origin(i));
        }
    }
    for i in 0..WATERS {
        push_residue(&mut b, "TIP3", next(), 3, water_origin(i));
    }
    for i in 0..PROTEINS {
        push_residue(&mut b, "ALA", next(), PROTEIN_ATOMS, protein_origin(i));
        if i % (PROTEINS / GLYCANS) == 0 {
            let bonded = protein_origin(i) + Vec3::new(1.45, 0.0, 0.0);
            push_residue(&mut b, "NAG", next(), 6, bonded);
        }
    }
    b
}

/// Milliseconds for the classes alone and for classes plus roles.
fn time_classification(unknown: usize) -> (f64, f64, Topology, usize) {
    let b = system(unknown);
    let atoms = b.positions.len();
    let (mut topology, positions) = (b.topology, b.positions);
    let started = Instant::now();
    vv_core::residue_class::classify(&topology, Some(&positions));
    let classes_ms = started.elapsed().as_secs_f64() * 1e3;
    let started = Instant::now();
    topology.assign_residue_classes(&positions);
    let total_ms = started.elapsed().as_secs_f64() * 1e3;
    (classes_ms, total_ms, topology, atoms)
}

#[test]
#[ignore = "timing"]
fn classification_time_at_millions_of_atoms() {
    for unknown in [0, UNKNOWN_LIGANDS] {
        let (classes_ms, total_ms, topology, atoms) = time_classification(unknown);
        let counts = topology.class_counts;
        println!(
            "{atoms} atoms, {} residues, {unknown} unknown-name residues: classes {classes_ms:.0} ms, classes + roles {total_ms:.0} ms",
            counts.iter().sum::<u32>()
        );
        assert_eq!(counts[ResidueClass::Water.index()] as usize, WATERS);
        assert_eq!(counts[ResidueClass::Lipid.index()] as usize, LIPIDS);
        assert_eq!(
            counts[ResidueClass::SmallMolecule.index()] as usize,
            unknown
        );
        assert_eq!(counts[ResidueClass::Glycan.index()] as usize, GLYCANS);
        let free_glycans = (0..topology.residue_count())
            .filter(|&r| topology.residue_class(r) == ResidueClass::Glycan)
            .filter(|&r| topology.residue_roles(r).contains(Roles::LIGAND))
            .count();
        assert_eq!(free_glycans, 0, "every glycan is bonded to the protein");
    }
}
