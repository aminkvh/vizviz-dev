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
    let text = chains_pdb(seqs);
    from_structure(vv_io::pdb::parse(text.as_bytes()).unwrap())
}

fn chains_pdb(seqs: &[&str]) -> String {
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
    text
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

/// Loads the chains through a scene, as the application does.
fn scene_of(seqs: &[&str], tag: &str) -> (vv_scene::Scene, vv_scene::StructureId) {
    use vv_scene::{Command, CommandHistory, Scene};
    let path = std::env::temp_dir().join(format!("vizviz-{tag}-{}.pdb", std::process::id()));
    std::fs::write(&path, chains_pdb(seqs)).unwrap();
    let mut scene = Scene::new();
    let mut history = CommandHistory::new(10);
    let loaded = history.dispatch(&mut scene, Command::LoadStructure { path: path.clone() });
    std::fs::remove_file(path).unwrap();
    loaded.unwrap();
    let id = scene.structures().next().unwrap().0;
    (scene, id)
}

#[test]
fn anarci_without_its_executable_falls_back_to_native_and_says_why() {
    use crate::sequence::anarci::Backend;
    use crate::sequence::cache::Cache;
    use std::time::{Duration, Instant};

    let (scene, id) = scene_of(&[HEAVY, LIGHT], "anarci-missing");
    let loaded = scene.structure(id).unwrap();
    let mut cache = Cache::default();
    cache.set_anarci_exe(Some("no-such-anarci".into()));
    let settings = AntibodySettings {
        backend: Backend::Anarci,
        ..Default::default()
    };
    let env = cache.env(&scene, settings, false);
    let all = 0..residues(loaded);

    let first = cache.track(id, loaded, &Antibody, &env);
    assert!(!first.any_in(&all), "empty while ANARCI runs");
    assert_eq!(
        cache.anarci_notice(&scene).as_deref(),
        Some("ANARCI running\u{2026}")
    );

    let start = Instant::now();
    let shown = loop {
        let t = cache.track(id, loaded, &Antibody, &env);
        if t.any_in(&all) {
            break t;
        }
        assert!(start.elapsed() < Duration::from_secs(30), "never fell back");
        std::thread::sleep(Duration::from_millis(5));
    };
    let native = run(loaded, AntibodySettings::default());
    assert!((0..residues(loaded)).all(|r| shown.kind(r) == native.kind(r)));
    let notice = cache.anarci_notice(&scene).unwrap();
    assert!(
        notice.starts_with("ANARCI failed: `no-such-anarci` is not a file")
            && notice.ends_with("showing native numbering"),
        "{notice}"
    );
}

#[test]
fn native_backend_never_starts_anarci() {
    use crate::sequence::cache::Cache;

    let (scene, id) = scene_of(&[HEAVY], "anarci-native");
    let loaded = scene.structure(id).unwrap();
    let mut cache = Cache::default();
    cache.set_anarci_exe(Some("no-such-anarci".into()));
    let env = cache.env(&scene, AntibodySettings::default(), false);
    let t = cache.track(id, loaded, &Antibody, &env);
    assert!(t.any_in(&(0..residues(loaded))));
}

/// `(agreeing, total)` residue labels per scheme between native numbering
/// and the real ANARCI, over domains both found.
fn live_agreement(l: &LoadedStructure) -> Vec<(Scheme, usize, usize)> {
    use crate::sequence::anarci;
    use vv_core::antibody::find_in_residues;

    let top = &l.structure.topology;
    let rows = chain_rows(top);
    let chains = anarci::chain_sequences(top, &rows);
    let outcome = anarci::number(None, &chains).unwrap_or_else(|e| panic!("{e}"));
    Scheme::ALL
        .iter()
        .map(|&scheme| {
            let (mut agree, mut total) = (0, 0);
            for ((_, residues), external) in rows.iter().zip(&outcome) {
                for native in find_in_residues(top, residues.clone()) {
                    let Some(theirs) = external.iter().find(|e| e.chain == native.chain) else {
                        continue;
                    };
                    let Some(notes) = theirs.annotate(scheme, CdrDefinition::Kabat) else {
                        println!(
                            "{scheme:?}: no labels for {:?} {}..{}",
                            theirs.chain, theirs.start, theirs.end
                        );
                        continue;
                    };
                    for (index, label) in native.numbering(scheme) {
                        total += 1;
                        agree +=
                            usize::from(notes.iter().any(|a| a.index == index && a.label == label));
                    }
                }
            }
            (scheme, agree, total)
        })
        .collect()
}

#[test]
#[ignore = "needs ANARCI (VIZVIZ_ANARCI, PATH or WSL); prints agreement with native numbering"]
fn live_anarci_numbers_trastuzumab_like_native() {
    let l = chains_of(&[HEAVY, LIGHT]);
    for (scheme, agree, total) in live_agreement(&l) {
        println!("trastuzumab {scheme:?}: {agree}/{total} residues agree");
        assert!(total > 200 && agree * 100 >= total * 95, "{scheme:?}");
    }
}

#[test]
#[ignore = "needs ANARCI and $VIZVIZ_FIXTURES_REAL/1N8Z.cif; prints agreement with native numbering"]
fn live_anarci_numbers_1n8z_like_native() {
    let root = std::env::var_os("VIZVIZ_FIXTURES_REAL").map(PathBuf::from);
    let Some(path) = root.map(|r| r.join("1N8Z.cif")).filter(|p| p.exists()) else {
        return;
    };
    let l = from_structure(vv_io::load(path).unwrap());
    for (scheme, agree, total) in live_agreement(&l) {
        println!("1N8Z {scheme:?}: {agree}/{total} residues agree");
        assert!(total > 200);
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
