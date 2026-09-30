//! Checks of the strip's data on whole structures: the polymer-only rows
//! and ligand summary, unobserved residues placed by sequence alignment,
//! conservation and UniProt features. Entries under `fixtures/real` are
//! downloaded, not committed: a check that needs one is skipped without it.

use std::collections::{BTreeSet, HashMap};
use std::path::PathBuf;
use std::sync::Arc;

use vv_core::seqfeat::{unobserved, unobserved_in_entity, EntityChain, Gap};
use vv_core::AnnotationCategory;
use vv_scene::{Command, CommandHistory, Scene, StructureId};

use super::peers::Peers;
use super::providers::{Conservation, Uniprot};
use super::rows::{chain_rows, rows_of};
use super::tracks::{Extras, TrackContext, TrackData, TrackProvider};
use super::uniprot;

fn fixture(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../fixtures")
        .join(rel)
}

fn open(rel: &str) -> Option<(Scene, StructureId)> {
    let path = fixture(rel);
    if !path.exists() {
        eprintln!("skipped: {rel} is not downloaded");
        return None;
    }
    let mut scene = Scene::new();
    let mut history = CommandHistory::new(10);
    history
        .dispatch(&mut scene, Command::LoadStructure { path })
        .unwrap();
    let id = scene.structures().next().unwrap().0;
    Some((scene, id))
}

fn row_labels(scene: &Scene) -> Vec<String> {
    rows_of(scene, |_, l| {
        Arc::new(super::rows::ligand_groups(&l.structure.topology))
    })
    .into_iter()
    .map(|r| r.label)
    .collect()
}

fn ligand_texts(scene: &Scene) -> Vec<String> {
    rows_of(scene, |_, l| {
        Arc::new(super::rows::ligand_groups(&l.structure.topology))
    })
    .into_iter()
    .flat_map(|r| r.ligands.iter().map(|g| g.text()).collect::<Vec<_>>())
    .collect()
}

#[test]
fn hemoglobin_shows_four_chains_and_its_ligands() {
    let (scene, _) = open("small/4HHB.cif").unwrap();
    let labels = row_labels(&scene);
    assert_eq!(labels.len(), 5, "{labels:?}");
    assert!(labels[4].ends_with("ligands"));
    assert_eq!(ligand_texts(&scene), ["HEM \u{D7}4", "PO4 \u{D7}2"]);
}

#[test]
fn a_pdb_format_entry_gets_the_same_rows() {
    let (scene, _) = open("small/4HHB.pdb").unwrap();
    assert_eq!(row_labels(&scene).len(), 5);
    assert_eq!(ligand_texts(&scene), ["HEM \u{D7}4", "PO4 \u{D7}2"]);
}

#[test]
fn chains_of_only_sugars_are_not_sequence_rows() {
    let Some((scene, id)) = open("real/1N8Z.cif.gz") else {
        return;
    };
    let top = &scene.structure(id).unwrap().structure.topology;
    assert!(top.chain_count() > 3, "the entry has ligand-only chains");
    let labels = row_labels(&scene);
    assert_eq!(labels.len(), 4, "{labels:?}");
    assert!(labels[3].ends_with("ligands"));
    let ligands = ligand_texts(&scene);
    assert!(ligands.iter().any(|t| t.starts_with("NAG")), "{ligands:?}");
    assert!(!ligands.iter().any(|t| t.starts_with("HOH")));
}

/// Seven modeled residues numbered 11 12 13 1 2 3 5 of a ten-residue chain
/// whose numbering wraps, with 14 15 and 4 unobserved.
fn wrapped_chain() -> (vv_core::Topology, EntityChain) {
    let names = [
        "ALA", "LYS", "VAL", "LEU", "GLY", "ASP", "GLU", "ARG", "THR", "SER",
    ];
    let numbers = [11, 12, 13, 14, 15, 1, 2, 3, 4, 5];
    let unobserved_positions = [3, 4, 8];
    let mut text = String::from("REMARK 465   M RES C SSSEQI\n");
    for &p in &unobserved_positions {
        text += &format!("REMARK 465     {:<3} A{:>6}\n", names[p], numbers[p]);
    }
    let mut serial = 0;
    for p in (0..10).filter(|p| !unobserved_positions.contains(p)) {
        serial += 1;
        text += &format!(
            "ATOM  {serial:>5}  CA  {:<3} A{:>4}    {:8.3}{:8.3}{:8.3}  1.00  0.00           C\n",
            names[p],
            numbers[p],
            serial as f32 * 3.8,
            0.0,
            0.0
        );
    }
    let top = (*vv_io::pdb::parse(text.as_bytes()).unwrap().topology).clone();
    let entity = EntityChain {
        names: names.iter().map(|n| n.to_string()).collect(),
        numbers: numbers.iter().map(|&n| Some((n, 0))).collect(),
    };
    (top, entity)
}

fn gap_summary(gaps: &[Gap]) -> Vec<(u32, Vec<i32>)> {
    gaps.iter()
        .map(|g| {
            let seqs = g.residues.iter().filter_map(|u| u.number.map(|n| n.0));
            (g.before, seqs.collect())
        })
        .collect()
}

#[test]
fn non_monotonic_numbering_places_gaps_by_alignment_not_by_number() {
    let (top, entity) = wrapped_chain();
    let residues = 0..top.residue_count() as u32;
    let placed = unobserved_in_entity(&top, residues.clone(), &entity, 1);
    assert_eq!(gap_summary(&placed), [(3, vec![14, 15]), (6, vec![4])]);
    let by_number = unobserved(&top, residues);
    assert_ne!(gap_summary(&by_number), gap_summary(&placed));
}

#[test]
fn an_entity_without_numbers_borrows_them_from_the_unobserved_list() {
    let (top, mut entity) = wrapped_chain();
    entity.numbers.clear();
    let placed = unobserved_in_entity(&top, 0..top.residue_count() as u32, &entity, 1);
    assert_eq!(gap_summary(&placed), [(3, vec![14, 15]), (6, vec![4])]);
}

fn add_unobserved_rows(top: &mut vv_core::Topology, rows: &[(&str, &str)]) {
    let items = [
        "polymer_flag",
        "occupancy_flag",
        "auth_asym_id",
        "auth_comp_id",
        "auth_seq_id",
        "PDB_ins_code",
        "PDB_model_num",
    ];
    top.annotations.categories.push(AnnotationCategory {
        name: "pdbx_unobs_or_zero_occ_residues".into(),
        items: items.map(String::from).to_vec(),
        rows: rows
            .iter()
            .map(|(seq, model)| {
                ["Y", "1", "A", "LEU", seq, "", model]
                    .map(String::from)
                    .to_vec()
            })
            .collect(),
    });
}

#[test]
fn an_ensemble_uses_the_gap_list_of_the_shown_model() {
    let names = ["ALA", "LYS", "VAL", "LEU", "GLY", "ASP", "GLU", "ARG"];
    let text = (1..=6)
        .map(|i| {
            format!(
                "ATOM  {i:>5}  CA  {} A{i:>4}    {:8.3}{:8.3}{:8.3}  1.00  0.00           C\n",
                names[i - 1],
                i as f32 * 3.8,
                0.0,
                0.0
            )
        })
        .collect::<String>();
    let mut top = (*vv_io::pdb::parse(text.as_bytes()).unwrap().topology).clone();
    add_unobserved_rows(&mut top, &[("7", "1"), ("7", "2"), ("8", "2")]);
    let entity = EntityChain {
        names: names.map(String::from).to_vec(),
        numbers: (1..=8).map(|n| Some((n, 0))).collect(),
    };
    let listed = |model| {
        let gaps = unobserved_in_entity(&top, 0..6, &entity, model);
        gaps.iter().map(|g| g.residues.len()).sum::<usize>()
    };
    assert_eq!(listed(1), 1);
    assert_eq!(listed(2), 2);
}

fn numbers_of(gaps: &[Gap]) -> BTreeSet<(i32, u8)> {
    gaps.iter()
        .flat_map(|g| g.residues.iter().filter_map(|u| u.number))
        .collect()
}

/// One polymer chain: the numbers alignment finds missing, the numbers the
/// file lists that no modeled residue carries, and how many modeled
/// residues have an insertion code.
struct ChainMissing {
    name: String,
    placed: BTreeSet<(i32, u8)>,
    listed: BTreeSet<(i32, u8)>,
    inserted: usize,
}

fn missing_by_both_ways(scene: &Scene, id: StructureId) -> Vec<ChainMissing> {
    let loaded = scene.structure(id).unwrap();
    let top = &loaded.structure.topology;
    let entities = vv_io::seqdata::entity::read(loaded.path.as_deref().unwrap()).unwrap();
    let mut out = Vec::new();
    for (name, residues) in chain_rows(top) {
        let Some(entity) = entities.get(&name) else {
            continue;
        };
        let modeled: BTreeSet<(i32, u8)> = residues
            .clone()
            .map(|r| {
                (
                    top.residues[r as usize].auth_seq_id,
                    top.residues[r as usize].ins_code,
                )
            })
            .collect();
        let inserted = residues
            .clone()
            .filter(|&r| top.residues[r as usize].ins_code != 0)
            .count();
        let placed = numbers_of(&unobserved_in_entity(top, residues.clone(), entity, 1));
        let listed: BTreeSet<_> = numbers_of(&unobserved(top, residues))
            .difference(&modeled)
            .copied()
            .collect();
        out.push(ChainMissing {
            name,
            placed,
            listed,
            inserted,
        });
    }
    out
}

#[test]
fn alignment_finds_the_residues_the_file_lists_as_unobserved() {
    let mut chains = 0;
    for entry in ["small/4HHB.cif", "real/1N8Z.cif.gz", "real/1IGY.cif.gz"] {
        let Some((scene, id)) = open(entry) else {
            continue;
        };
        let top = &scene.structure(id).unwrap().structure.topology;
        let has_table = top
            .annotations
            .category("pdbx_unobs_or_zero_occ_residues")
            .is_some();
        for c in missing_by_both_ways(&scene, id) {
            if has_table {
                assert_eq!(c.placed, c.listed, "{entry} chain {}", c.name);
            }
            chains += 1;
        }
    }
    assert!(chains >= 3);
}

#[test]
fn the_her2_chain_of_1n8z_has_its_26_unobserved_residues_placed() {
    let Some((scene, id)) = open("real/1N8Z.cif.gz") else {
        return;
    };
    let chains = missing_by_both_ways(&scene, id);
    let her2 = chains.iter().find(|c| c.name == "C").unwrap();
    assert_eq!(her2.placed.len(), 26, "{:?}", her2.placed);
}

#[test]
fn an_antibody_with_insertion_codes_does_not_call_zero_occupancy_residues_missing() {
    let Some((scene, id)) = open("real/1IGY.cif.gz") else {
        return;
    };
    let chains = missing_by_both_ways(&scene, id);
    let inserted: usize = chains.iter().map(|c| c.inserted).sum();
    assert!(inserted > 5, "{inserted} insertion-coded residues");
    // The file lists four heavy-chain residues, which have coordinates at
    // zero occupancy: neither way calls them missing.
    assert!(chains[1].placed.is_empty() && chains[1].listed.is_empty());
}

fn run(provider: &dyn TrackProvider, scene: &Scene, id: StructureId, extras: Extras) -> TrackData {
    let loaded = scene.structure(id).unwrap();
    let coords = loaded.structure.frame(0);
    let rows = chain_rows(&loaded.structure.topology);
    provider.compute(&TrackContext {
        loaded,
        positions: coords.positions(),
        rows: &rows,
        antibody: Default::default(),
        extras,
    })
}

#[test]
fn hemoglobin_alpha_and_beta_conserve_about_a_third_of_their_columns() {
    let (scene, id) = open("small/4HHB.cif").unwrap();
    let peers = Peers::of(&scene);
    let track = run(
        &Conservation,
        &scene,
        id,
        Extras {
            peers: Some(&peers),
            structure: Some(id),
            ..Default::default()
        },
    );
    let top = &scene.structure(id).unwrap().structure.topology;
    let alpha = chain_rows(top)
        .into_iter()
        .find(|(n, _)| n == "A")
        .unwrap()
        .1;
    assert_eq!(alpha.len(), 141);
    let identical = alpha.clone().filter(|&r| track.kind(r) == 5).count();
    assert!(
        (40..100).contains(&identical),
        "{identical} of 141 columns identical in all four chains"
    );
    assert!(alpha.clone().all(|r| track.kind(r) != 0));
    let note = track.describe(alpha.start).unwrap();
    assert!(note.contains("of 4 chains"), "{note}");
}

fn sifts_and_features(entry: &str, accessions: &[&str], by_label: bool) -> uniprot::Data {
    let read =
        |name: &str| std::fs::read_to_string(fixture(&format!("small/uniprot/{name}"))).unwrap();
    let mut features = HashMap::new();
    for acc in accessions {
        features.insert(
            acc.to_string(),
            uniprot::parse_features(&read(&format!("{acc}.json"))),
        );
    }
    uniprot::Data {
        segments: uniprot::parse_sifts(&read(&format!("sifts_{entry}.json")), by_label),
        features,
    }
}

fn chain_start(scene: &Scene, id: StructureId, name: &str) -> u32 {
    let top = &scene.structure(id).unwrap().structure.topology;
    chain_rows(top)
        .into_iter()
        .find(|(n, _)| n == name)
        .unwrap()
        .1
        .start
}

#[test]
fn uniprot_features_mark_the_heme_histidines_of_hemoglobin() {
    let (scene, id) = open("small/4HHB.cif").unwrap();
    let data = sifts_and_features("4hhb", &["P69905", "P68871"], true);
    let track = run(
        &Uniprot,
        &scene,
        id,
        Extras {
            uniprot: Some(&data),
            ..Default::default()
        },
    );
    let top = &scene.structure(id).unwrap().structure.topology;
    let site = 4;
    for (chain, offset, ligand) in [("A", 87, "heme b"), ("B", 92, "heme b")] {
        let r = chain_start(&scene, id, chain) + offset - 1;
        assert_eq!(top.residue_name(r as usize), "HIS");
        assert_eq!(track.kind(r), site, "chain {chain}");
        let note = track.describe(r).unwrap();
        assert!(note.contains(ligand), "{note}");
    }
    let distal = chain_start(&scene, id, "A") + 57;
    assert_eq!(top.residue_name(distal as usize), "HIS");
    assert!(track.describe(distal).unwrap().contains("O2"));
}

#[test]
fn uniprot_features_reach_the_her2_chain_of_1n8z() {
    let Some((scene, id)) = open("real/1N8Z.cif.gz") else {
        return;
    };
    let data = sifts_and_features("1n8z", &["P04626"], true);
    let track = run(
        &Uniprot,
        &scene,
        id,
        Extras {
            uniprot: Some(&data),
            ..Default::default()
        },
    );
    let top = &scene.structure(id).unwrap().structure.topology;
    // UniProt Asn 187 carries an N-glycan.
    let range = chain_rows(top)
        .into_iter()
        .find(|(n, _)| n == "C")
        .unwrap()
        .1;
    let mapped = uniprot::map_residues(&data, "C", top, range);
    let asn = mapped.iter().find(|m| m.position == 187).unwrap().residue;
    assert_eq!(top.residue_name(asn as usize), "ASN");
    assert!(track.kind(asn) != 0);
    assert!(track.describe(asn).unwrap().contains("Glycosylation 187"));
    let rows = rows_of(&scene, |_, _| Default::default());
    let her2 = rows.iter().find(|r| r.label.ends_with(" C")).unwrap();
    assert!(track.any_in(&her2.residues));
}

/// Needs the network or a warm cache: `cargo test -- --ignored`.
#[test]
#[ignore]
fn the_live_fetch_reaches_hemoglobins_entries() {
    use super::cache::{Cache, Fetch};
    let (scene, id) = open("small/4HHB.cif").unwrap();
    let loaded = scene.structure(id).unwrap();
    let mut cache = Cache::default();
    let start = std::time::Instant::now();
    loop {
        match cache.uniprot_status(id, loaded) {
            Fetch::Pending => {}
            Fetch::Ready(_) => return,
            Fetch::Unavailable => panic!("unavailable"),
        }
        assert!(start.elapsed().as_secs() < 60);
        std::thread::sleep(std::time::Duration::from_millis(50));
    }
}
