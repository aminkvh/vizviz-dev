//! Entity descriptions, source organisms and branch descriptors survive a
//! read, write, read through both writers.

use std::collections::BTreeMap;
use std::path::PathBuf;

use vv_core::Structure;
use vv_io::SaveOptions;

fn fixture(dir: &str, name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../fixtures")
        .join(dir)
        .join(name)
}

fn round_trip(s: &Structure, ext: &str) -> (Structure, String) {
    static NEXT: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
    let n = NEXT.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    let path = std::env::temp_dir().join(format!("vv_entity_{}_{n}.{ext}", std::process::id()));
    let opts = SaveOptions {
        atoms: None,
        frames: &[0],
    };
    vv_io::save(s, &path, &opts).unwrap();
    let text = std::fs::read_to_string(&path).unwrap();
    let back = vv_io::load(&path).unwrap();
    std::fs::remove_file(&path).ok();
    (back, text)
}

/// Each entity's description with its source organism, keyed by description
/// so the renumbered ids do not matter.
fn described(s: &Structure) -> BTreeMap<String, Option<String>> {
    let a = &s.topology.annotations;
    a.entities()
        .into_iter()
        .map(|(id, description)| {
            let organism = ["entity_src_gen", "entity_src_nat"].iter().find_map(|c| {
                let cat = a.category(c)?;
                let row = (0..cat.rows.len()).find(|&r| cat.get("entity_id", r) == Some(&id))?;
                cat.get("pdbx_gene_src_scientific_name", row)
                    .or_else(|| cat.get("pdbx_organism_scientific", row))
                    .map(str::to_owned)
            });
            (description, organism)
        })
        .collect()
}

fn descriptors(s: &Structure) -> Vec<String> {
    let cat = s
        .topology
        .annotations
        .category("pdbx_entity_branch_descriptor")
        .expect("branch descriptors");
    let mut out: Vec<String> = (0..cat.rows.len())
        .filter_map(|r| cat.get("descriptor", r).map(str::to_owned))
        .collect();
    out.sort();
    out
}

#[test]
fn hemoglobin_mmcif_keeps_descriptions_and_organisms() {
    let s = vv_io::load(fixture("small", "4HHB.cif")).unwrap();
    let before = described(&s);
    assert!(before.len() >= 3, "{before:?}");
    assert!(before.values().any(Option::is_some));
    let (back, _) = round_trip(&s, "cif");
    assert_eq!(described(&back), before);
}

#[test]
fn hemoglobin_pdb_keeps_molecule_names_and_organisms() {
    let s = vv_io::load(fixture("small", "4HHB.pdb")).unwrap();
    let before = described(&s);
    let (back, text) = round_trip(&s, "pdb");
    assert!(text.contains("MOLECULE: HEMOGLOBIN SUBUNIT ALPHA;"));
    assert!(text.contains("ORGANISM_SCIENTIFIC: HOMO SAPIENS;"));
    assert_eq!(described(&back), before);
    let strands = back.topology.annotations.category("entity_poly").unwrap();
    let chains: Vec<&str> = (0..strands.rows.len())
        .filter_map(|r| strands.get("pdbx_strand_id", r))
        .collect();
    assert_eq!(chains, ["A,C", "B,D"]);
}

#[test]
fn hemoglobin_mmcif_to_pdb_and_back_keeps_molecule_names() {
    let s = vv_io::load(fixture("small", "4HHB.cif")).unwrap();
    let names: Vec<String> = described(&s).into_keys().collect();
    let (via_pdb, _) = round_trip(&s, "pdb");
    let polymers: Vec<String> = described(&via_pdb).into_keys().collect();
    assert_eq!(polymers.len(), 2, "PDB lists polymers only: {polymers:?}");
    assert!(polymers.iter().all(|p| names.contains(p)));
}

#[test]
fn glycosylated_receptor_keeps_branch_descriptors_and_sources() {
    let path = fixture("large", "6X3Z.cif");
    if !path.exists() {
        eprintln!("skipping: {} is not downloaded", path.display());
        return;
    }
    let s = vv_io::load(path).unwrap();
    let before = described(&s);
    let branch = descriptors(&s);
    assert!(!branch.is_empty());
    let (back, text) = round_trip(&s, "cif");
    assert!(text.contains("_pdbx_entity_branch_descriptor."));
    assert_eq!(described(&back), before);
    assert_eq!(descriptors(&back), branch);
}
