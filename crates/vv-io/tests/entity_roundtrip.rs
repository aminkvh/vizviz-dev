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
    assert_eq!(branched_entities(&back), branched_entities(&s));
}

/// 6X3Z's three branched entities keep their rows and descriptors (they
/// differ in composition; the synthetic test below covers linkage).
#[test]
fn same_composition_glycans_with_different_linkage_stay_separate() {
    let path = fixture("large", "6X3Z.cif");
    if !path.exists() {
        eprintln!("skipping: {} is not downloaded", path.display());
        return;
    }
    let s = vv_io::load(path).unwrap();
    assert_eq!(branched_entities(&s), 3);
    let (back, _) = round_trip(&s, "cif");
    assert_eq!(branched_entities(&back), 3);
    assert_eq!(descriptors(&back).len(), descriptors(&s).len());
}

fn branched_entities(s: &Structure) -> usize {
    s.topology
        .annotations
        .category("pdbx_entity_branch")
        .map_or(0, |c| c.rows.len())
}

/// A three-sugar chain `chain` with LINK records `links` as
/// `(donor atom, residue, acceptor residue)`: each residue has C1, O3, O4.
fn glycan_pdb(chain: char, links: &[(&str, u32, u32)], serial: &mut u32) -> String {
    let mut out = String::new();
    for seq in 1..=3u32 {
        for (k, (name, el)) in [("C1", "C"), ("O3", "O"), ("O4", "O")].iter().enumerate() {
            *serial += 1;
            let x = 3.0 * seq as f32 + k as f32;
            let y = if chain == 'B' { 0.0 } else { 20.0 };
            out += &format!(
                "HETATM{serial:5} {name:<4} NAG {chain}{seq:4}    {x:8.3}{y:8.3}{:8.3}  1.00  0.00          {el:>2}
",
                0.0
            );
        }
    }
    for &(atom, from, to) in links {
        let mut line = vec![b' '; 80];
        let mut put = |start: usize, text: &str| {
            line[start..start + text.len()].copy_from_slice(text.as_bytes())
        };
        put(0, "LINK");
        put(12, &format!("{atom:<4}"));
        put(17, "NAG");
        put(21, &chain.to_string());
        put(22, &format!("{from:>4}"));
        put(42, "C1");
        put(47, "NAG");
        put(51, &chain.to_string());
        put(52, &format!("{to:>4}"));
        out += &String::from_utf8(line).unwrap();
        out.push('\n');
    }
    out
}

fn branched_rows(text: &str) -> usize {
    text.lines()
        .filter(|l| l.ends_with(" branched") || l.contains(" branched "))
        .count()
}

#[test]
fn glycans_of_one_composition_linked_differently_are_two_entities() {
    let mut serial = 0;
    let mut pdb = glycan_pdb('B', &[("O4", 1, 2), ("O4", 2, 3)], &mut serial);
    pdb += &glycan_pdb('C', &[("O4", 1, 2), ("O3", 2, 3)], &mut serial);
    let same = glycan_pdb('D', &[("O4", 1, 2), ("O4", 2, 3)], &mut serial);
    let path = std::env::temp_dir().join(format!("vv_glycan_link_{}.pdb", std::process::id()));
    std::fs::write(
        &path,
        format!(
            "{pdb}{same}END
"
        ),
    )
    .unwrap();
    let s = vv_io::load(&path).unwrap();
    std::fs::remove_file(&path).ok();
    let (_, text) = round_trip(&s, "cif");
    assert_eq!(branched_rows(&text), 2, "{text}");
}
