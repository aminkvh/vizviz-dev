//! The antibody track on hand-made trastuzumab chains, on structures with
//! no antibody, and on a real Fab when its fixture is present.

use std::path::PathBuf;

use vv_core::antibody::{CdrDefinition, Scheme};
use vv_scene::LoadedStructure;

use super::Antibody;
use crate::sequence::rows::{chain_rows, one_letter};
use crate::sequence::tracks::{AntibodySettings, TrackContext, TrackData, TrackProvider};

const HEAVY: &str = "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS";
const LIGHT: &str = "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK";

const CODES: [(char, &str); 20] = [
    ('A', "ALA"),
    ('R', "ARG"),
    ('N', "ASN"),
    ('D', "ASP"),
    ('C', "CYS"),
    ('Q', "GLN"),
    ('E', "GLU"),
    ('G', "GLY"),
    ('H', "HIS"),
    ('I', "ILE"),
    ('L', "LEU"),
    ('K', "LYS"),
    ('M', "MET"),
    ('F', "PHE"),
    ('P', "PRO"),
    ('S', "SER"),
    ('T', "THR"),
    ('W', "TRP"),
    ('Y', "TYR"),
    ('V', "VAL"),
];

fn from_structure(structure: vv_core::Structure) -> LoadedStructure {
    let bonds = vv_core::bonds::perceive(&structure.topology, structure.frame(0).positions());
    LoadedStructure::new(structure, None, "test".into(), bonds)
}

fn fixture(path: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../fixtures")
        .join(path)
}

/// One CA-only chain per sequence, named A, B, ...
fn chains_of(seqs: &[&str]) -> LoadedStructure {
    let mut text = String::new();
    let mut serial = 0;
    for (chain, seq) in seqs.iter().enumerate() {
        for (i, letter) in seq.chars().enumerate() {
            serial += 1;
            let name = CODES.iter().find(|c| c.0 == letter).unwrap().1;
            text.push_str(&format!(
                "ATOM  {serial:>5}  CA  {name:<3} {}{:>4}    {:8.3}{:8.3}{:8.3}  1.00  0.00           C\n",
                (b'A' + chain as u8) as char,
                i + 1,
                serial as f32 * 3.8,
                0.0,
                0.0
            ));
        }
        text.push_str("TER\n");
    }
    from_structure(vv_io::pdb::parse(text.as_bytes()).unwrap())
}

fn run(l: &LoadedStructure, antibody: AntibodySettings) -> TrackData {
    let coords = l.structure.frame(0);
    let rows = chain_rows(&l.structure.topology);
    Antibody.compute(&TrackContext {
        loaded: l,
        positions: coords.positions(),
        rows: &rows,
        antibody,
        extras: Default::default(),
    })
}

fn residues(l: &LoadedStructure) -> u32 {
    l.structure.topology.residue_count() as u32
}

/// One-letter codes of the residues of chain row `chain` carrying `kind`.
fn letters_of_kind(l: &LoadedStructure, t: &TrackData, chain: usize, kind: u8) -> String {
    let top = &l.structure.topology;
    chain_rows(top)[chain]
        .1
        .clone()
        .filter(|&r| t.kind(r) == kind)
        .map(|r| one_letter(top.residue_name(r as usize)).unwrap())
        .collect()
}

#[test]
fn trastuzumab_cdr_bars_are_the_published_kabat_loops() {
    let l = chains_of(&[HEAVY, LIGHT]);
    let t = run(&l, AntibodySettings::default());
    for (chain, cdrs) in [
        (0, ["DTYIH", "RIYPTNGYTRYADSVKG", "WGGDGFYAMDY"]),
        (1, ["RASQDVNTAVA", "SASFLYS", "QQHYTTPPT"]),
    ] {
        for (n, seq) in cdrs.iter().enumerate() {
            assert_eq!(letters_of_kind(&l, &t, chain, n as u8 + 2), *seq);
        }
    }
}

#[test]
fn tooltip_names_the_number_region_and_definition() {
    let l = chains_of(&[HEAVY]);
    let t = run(&l, AntibodySettings::default());
    let notes: Vec<String> = (0..residues(&l)).filter_map(|r| t.describe(r)).collect();
    assert!(
        notes.contains(&"H52A \u{B7} CDR-H2 (Kabat)".to_string()),
        "{notes:?}"
    );
    assert!(notes.contains(&"H40 \u{B7} FR2 (Kabat)".to_string()));
}

#[test]
fn ticks_label_decades_and_insertions_and_the_domain_gets_a_badge() {
    let l = chains_of(&[HEAVY]);
    let t = run(&l, AntibodySettings::default());
    let ticks: Vec<&str> = (0..residues(&l)).filter_map(|r| t.tick(r)).collect();
    for wanted in ["10", "20", "52A", "82A"] {
        assert!(ticks.contains(&wanted), "{wanted} in {ticks:?}");
    }
    assert_eq!(t.badges().map(|b| b.1).collect::<Vec<_>>(), ["Heavy"]);
}

#[test]
fn scheme_and_definition_change_numbers_and_loops_independently() {
    let l = chains_of(&[HEAVY]);
    let chothia_loops = AntibodySettings {
        cdr: CdrDefinition::Chothia,
        ..Default::default()
    };
    let t = run(&l, chothia_loops);
    assert_eq!(letters_of_kind(&l, &t, 0, 2), "GFNIKDT");
    let imgt_numbers = AntibodySettings {
        scheme: Scheme::Imgt,
        ..Default::default()
    };
    let numbered = run(&l, imgt_numbers);
    assert_eq!(letters_of_kind(&l, &numbered, 0, 2), "DTYIH");
    let notes: String = (0..residues(&l))
        .filter_map(|r| numbered.describe(r))
        .collect();
    assert!(notes.contains("(IMGT numbering)"));
}

#[test]
fn non_antibody_structures_get_no_antibody_marks() {
    for path in ["small/4HHB.pdb", "small/1CRN.pdb"] {
        let l = from_structure(vv_io::load(fixture(path)).unwrap());
        let t = run(&l, AntibodySettings::default());
        assert!((0..residues(&l)).all(|r| t.kind(r) == 0), "{path}");
        assert_eq!(t.badges().count(), 0);
    }
}

#[test]
fn a_real_fab_gets_two_domains_when_the_fixture_is_present() {
    let path = fixture("real/1N8Z.cif");
    if !path.exists() {
        return;
    }
    let l = from_structure(vv_io::load(path).unwrap());
    let t = run(&l, AntibodySettings::default());
    let mut badges: Vec<&str> = t.badges().map(|b| b.1).collect();
    badges.sort_unstable();
    assert_eq!(badges, ["Heavy", "Kappa"]);
}
