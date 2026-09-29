//! Round-trip tests for the structure writers: write, read back, compare
//! topology and coordinates against the original.

use std::path::PathBuf;

use vv_core::glam::Vec3;
use vv_core::{CoordSet, Structure, Topology};
use vv_io::synth::{protein_like, SynthParams};
use vv_io::{Format, SaveOptions};

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../fixtures/small")
        .join(name)
}

fn load(name: &str) -> Structure {
    let s = vv_io::load(fixture(name)).unwrap_or_else(|e| panic!("{name}: {e}"));
    s.topology.validate().unwrap();
    s
}

/// Every line of a fixed-column PDB record type is exactly 80 columns.
fn assert_fixed_width_lines(text: &str) {
    for line in text.lines() {
        let rec = &line[..line.len().min(6)];
        let fixed = rec.starts_with("ATOM")
            || rec.starts_with("HETATM")
            || rec.starts_with("TER")
            || rec.starts_with("CONECT")
            || rec.starts_with("MODEL")
            || rec.starts_with("ENDMDL")
            || rec.starts_with("CRYST1");
        if fixed {
            assert_eq!(line.len(), 80, "{rec} line not 80 columns: {line:?}");
        }
    }
}

/// The distinct bonded pairs, ignoring how many times the source file
/// happened to repeat one (some real files list a metal-coordination
/// bond on more than one base atom's `CONECT` line; that redundancy
/// isn't chemistry, so the writer normalizes it to one `CONECT` pair and
/// this comparison does too).
fn sorted_bond_pairs(t: &Topology) -> Vec<[u32; 2]> {
    let mut pairs: Vec<[u32; 2]> = t
        .explicit_bonds
        .iter()
        .map(|b| {
            let [a, c] = b.atoms;
            if a < c {
                [a, c]
            } else {
                [c, a]
            }
        })
        .collect();
    pairs.sort_unstable();
    pairs.dedup();
    pairs
}

fn assert_positions_close(a: &Structure, b: &Structure, frame_a: usize, frame_b: usize) {
    for (pa, pb) in a
        .frame(frame_a)
        .positions()
        .iter()
        .zip(b.frame(frame_b).positions())
    {
        assert!((*pa - *pb).abs().max_element() <= 0.001, "{pa} vs {pb}");
    }
}

fn assert_topology_matches(
    original: &Structure,
    parsed: &Structure,
    check_bonds: bool,
    check_chain_names: bool,
) {
    let (o, p) = (&original.topology, &parsed.topology);
    assert_eq!(p.atom_count(), o.atom_count());
    assert_eq!(p.residue_count(), o.residue_count());
    assert_eq!(p.chain_count(), o.chain_count());
    for a in 0..o.atom_count() {
        assert_eq!(p.element[a], o.element[a], "atom {a} element");
        assert_eq!(p.atom_name(a), o.atom_name(a), "atom {a} name");
        assert!(
            (p.occupancy.get(a).copied().unwrap_or(1.0)
                - o.occupancy.get(a).copied().unwrap_or(1.0))
            .abs()
                < 0.005,
            "atom {a} occupancy"
        );
        assert!(
            (p.b_factor.get(a).copied().unwrap_or(0.0) - o.b_factor.get(a).copied().unwrap_or(0.0))
                .abs()
                < 0.005,
            "atom {a} b_factor"
        );
    }
    for r in 0..o.residue_count() {
        assert_eq!(p.residue_name(r), o.residue_name(r), "residue {r} name");
        assert_eq!(
            p.residues[r].auth_seq_id, o.residues[r].auth_seq_id,
            "residue {r} auth_seq_id"
        );
    }
    if check_chain_names {
        for c in 0..o.chain_count() {
            assert_eq!(p.chain_name(c), o.chain_name(c), "chain {c} name");
        }
    }
    assert_positions_close(original, parsed, 0, 0);
    if check_bonds {
        assert_eq!(sorted_bond_pairs(p), sorted_bond_pairs(o), "bonds");
    }
}

const FIXTURES: [&str; 3] = ["1CRN", "4HHB", "1UBQ"];

#[test]
fn pdb_round_trips_through_the_pdb_writer() {
    for id in FIXTURES {
        let original = load(&format!("{id}.pdb"));
        let mut bytes = Vec::new();
        let warnings = vv_io::pdb_write::write(&original, None, &[0], &mut bytes).unwrap();
        // 4HHB's own file lists each protein chain's waters as a second,
        // non-contiguous run under the same letter (e.g. two chains named
        // "A"); one of the pair collides and is remapped with a warning.
        // That's an expected, documented loss of the literal id, not a
        // writer bug, so chain-name identity isn't checked when it fires.
        let text = String::from_utf8(bytes.clone()).unwrap();
        assert_fixed_width_lines(&text);
        let parsed = vv_io::parse(&bytes, Format::Pdb).unwrap_or_else(|e| panic!("{id}: {e}"));
        parsed.topology.validate().unwrap();
        assert_topology_matches(&original, &parsed, true, warnings.is_empty());
    }
}

#[test]
fn mmcif_round_trips_through_the_mmcif_writer() {
    for id in FIXTURES {
        let original = load(&format!("{id}.cif"));
        let mut bytes = Vec::new();
        let warnings = vv_io::mmcif_write::write(&original, None, &[0], &mut bytes).unwrap();
        assert!(warnings.is_empty(), "{id}: {warnings:?}");
        let parsed = vv_io::parse(&bytes, Format::Mmcif).unwrap_or_else(|e| panic!("{id}: {e}"));
        parsed.topology.validate().unwrap();
        assert_topology_matches(&original, &parsed, true, true);
        for (a, b) in parsed
            .topology
            .explicit_bonds
            .iter()
            .zip(&original.topology.explicit_bonds)
        {
            assert_eq!(
                std::mem::discriminant(&a.kind),
                std::mem::discriminant(&b.kind)
            );
        }
    }
}

#[test]
fn pdb_writer_emits_cryst1_when_the_cell_is_known() {
    let original = load("1CRN.cif");
    let mut bytes = Vec::new();
    vv_io::pdb_write::write(&original, None, &[0], &mut bytes).unwrap();
    let text = String::from_utf8(bytes).unwrap();
    let cryst1 = text.lines().find(|l| l.starts_with("CRYST1")).unwrap();
    assert_eq!(cryst1.len(), 80);
    assert!(cryst1.contains("40.960"), "{cryst1}");
    // 1CRN's `_cell.Z_PDB` is 2, not the fallback of 1.
    assert_eq!(cryst1[66..70].trim(), "2", "{cryst1}");
}

#[test]
fn selection_export_writes_only_the_selected_atoms() {
    let original = load("4HHB.cif");
    let mask = vv_core::select(
        &original.topology,
        original.frame(0).positions(),
        "chain A and name CA",
    )
    .unwrap();
    let expected = mask.count_ones(..);
    assert!(expected > 0);

    let opts = SaveOptions {
        atoms: Some(&mask),
        frames: &[0],
    };
    let dir = std::env::temp_dir().join("vizviz-tests");
    std::fs::create_dir_all(&dir).unwrap();

    for (path, format) in [
        (dir.join("4hhb_ca.pdb"), Format::Pdb),
        (dir.join("4hhb_ca.cif"), Format::Mmcif),
    ] {
        vv_io::save(&original, &path, &opts).unwrap();
        let parsed = vv_io::load(&path).unwrap();
        assert_eq!(parsed.atom_count(), expected, "{path:?}");
        for a in 0..parsed.atom_count() {
            assert_eq!(parsed.topology.atom_name(a), "CA");
        }
        let _ = format;
    }
}

#[test]
fn multi_model_round_trips_through_pdb_and_mmcif() {
    let original = load("1CRN.pdb");
    let frame0 = original.frame(0).positions().to_vec();
    let frame1: Vec<Vec3> = frame0
        .iter()
        .map(|p| *p + Vec3::new(1.0, 0.0, 0.0))
        .collect();
    let topology = (*original.topology).clone();
    let two_frames =
        Structure::with_frames(topology, vec![CoordSet::new(frame0), CoordSet::new(frame1)])
            .unwrap();

    let mut pdb_bytes = Vec::new();
    vv_io::pdb_write::write(&two_frames, None, &[0, 1], &mut pdb_bytes).unwrap();
    let pdb_parsed = vv_io::parse(&pdb_bytes, Format::Pdb).unwrap();
    assert_eq!(pdb_parsed.frame_count(), 2);
    assert_positions_close(&two_frames, &pdb_parsed, 0, 0);
    assert_positions_close(&two_frames, &pdb_parsed, 1, 1);

    let mut cif_bytes = Vec::new();
    vv_io::mmcif_write::write(&two_frames, None, &[0, 1], &mut cif_bytes).unwrap();
    let cif_parsed = vv_io::parse(&cif_bytes, Format::Mmcif).unwrap();
    assert_eq!(cif_parsed.frame_count(), 2);
    assert_positions_close(&two_frames, &cif_parsed, 0, 0);
    assert_positions_close(&two_frames, &cif_parsed, 1, 1);
}

#[test]
fn large_synthetic_structure_round_trips_with_hybrid36_pdb_serials() {
    let original = protein_like(&SynthParams {
        residues_per_chain: 2000,
        ..SynthParams::new(120_000)
    });
    assert!(original.atom_count() > 99_999);

    let mut cif_bytes = Vec::new();
    vv_io::mmcif_write::write(&original, None, &[0], &mut cif_bytes).unwrap();
    let cif_parsed = vv_io::parse(&cif_bytes, Format::Mmcif).unwrap();
    assert_eq!(cif_parsed.atom_count(), original.atom_count());
    assert_eq!(cif_parsed.topology.serial, original.topology.serial);
    assert_positions_close(&original, &cif_parsed, 0, 0);

    let mut pdb_bytes = Vec::new();
    let warnings = vv_io::pdb_write::write(&original, None, &[0], &mut pdb_bytes).unwrap();
    assert!(warnings.is_empty());
    let text = String::from_utf8(pdb_bytes.clone()).unwrap();
    assert_fixed_width_lines(&text);
    let pdb_parsed = vv_io::parse(&pdb_bytes, Format::Pdb).unwrap();
    assert_eq!(pdb_parsed.atom_count(), original.atom_count());
    assert!(
        pdb_parsed.topology.serial.iter().any(|&s| s > 99_999),
        "hybrid-36 serials should decode past 99999"
    );
    assert_eq!(pdb_parsed.topology.element, original.topology.element);
    assert_positions_close(&original, &pdb_parsed, 0, 0);
}

#[test]
fn more_than_62_chains_fall_back_to_shared_question_mark_with_warnings() {
    let original = protein_like(&SynthParams {
        residues_per_chain: 1,
        ..SynthParams::new(700)
    });
    let chain_count = original.topology.chain_count();
    assert!(chain_count > 62, "chain_count = {chain_count}");

    let mut bytes = Vec::new();
    let warnings = vv_io::pdb_write::write(&original, None, &[0], &mut bytes).unwrap();
    assert_eq!(warnings.len(), chain_count - 26, "{warnings:?}");
    // Still parses, even though chains past the 62-character pool share `?`.
    let parsed = vv_io::parse(&bytes, Format::Pdb).unwrap();
    parsed.topology.validate().unwrap();
}

#[test]
fn pdb_gz_round_trips() {
    let original = load("1CRN.pdb");
    let dir = std::env::temp_dir().join("vizviz-tests");
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("1crn_save.pdb.gz");
    vv_io::save(
        &original,
        &path,
        &SaveOptions {
            atoms: None,
            frames: &[0],
        },
    )
    .unwrap();
    let parsed = vv_io::load(&path).unwrap();
    assert_eq!(parsed.atom_count(), original.atom_count());
    assert_positions_close(&original, &parsed, 0, 0);
}

#[test]
fn xyz_writes_element_and_coordinates_per_frame() {
    let original = load("1CRN.pdb");
    let mut bytes = Vec::new();
    vv_io::xyz_write::write(&original, None, &[0], &mut bytes).unwrap();
    let text = String::from_utf8(bytes).unwrap();
    let mut lines = text.lines();
    let natoms: usize = lines.next().unwrap().trim().parse().unwrap();
    assert_eq!(natoms, original.atom_count());
    lines.next().unwrap(); // comment line
    let first: Vec<&str> = lines.next().unwrap().split_whitespace().collect();
    assert_eq!(first[0], original.topology.element[0].symbol());
    let x: f32 = first[1].parse().unwrap();
    assert!((x - original.frame(0).positions()[0].x).abs() < 0.001);
    assert_eq!(lines.count(), natoms - 1);
}

#[test]
fn pqr_writes_charge_and_radius_columns() {
    // 4HHB (4779 atoms, HETATM hemes and waters) exercises serials >= 10000
    // and residue numbers >= 1000, where fixed-width fields with no
    // literal separator between them would run together.
    let original = load("4HHB.pdb");
    let mut bytes = Vec::new();
    vv_io::pqr_write::write(&original, None, &[0], &mut bytes).unwrap();
    let text = String::from_utf8(bytes).unwrap();
    let lines: Vec<&str> = text.lines().collect();
    assert_eq!(lines.len(), original.atom_count());
    for line in &lines {
        // ATOM/HETATM serial name resName chain resSeq x y z charge radius
        let fields: Vec<&str> = line.split_whitespace().collect();
        assert_eq!(fields.len(), 11, "{line:?}");
        for &f in &fields[6..] {
            f.parse::<f32>()
                .unwrap_or_else(|e| panic!("{f:?} in {line:?}: {e}"));
        }
    }
    let radius: f32 = lines[0]
        .split_whitespace()
        .nth(10)
        .unwrap()
        .parse()
        .unwrap();
    assert!(radius > 0.0);
}

#[test]
fn gro_writes_nm_coordinates_and_a_box_line() {
    let original = load("1CRN.pdb");
    let mut bytes = Vec::new();
    vv_io::gro_write::write(&original, None, &[0], &mut bytes).unwrap();
    let text = String::from_utf8(bytes).unwrap();
    let lines: Vec<&str> = text.lines().collect();
    let natoms: usize = lines[1].trim().parse().unwrap();
    assert_eq!(natoms, original.atom_count());
    assert_eq!(lines.len(), 2 + natoms + 1, "title, count, atoms, box");
    let atom_line = lines[2];
    let x: f32 = atom_line[20..28].trim().parse().unwrap();
    assert!((x - original.frame(0).positions()[0].x / 10.0).abs() < 0.001);
    let box_line: Vec<f32> = lines
        .last()
        .unwrap()
        .split_whitespace()
        .map(|s| s.parse().unwrap())
        .collect();
    assert_eq!(box_line.len(), 3);
    assert!(box_line.iter().all(|&v| v > 0.0));
}
