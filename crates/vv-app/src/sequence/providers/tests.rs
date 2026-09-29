//! Each track on real structures, and on hand-made chains for motifs.

use std::path::PathBuf;

use vv_scene::LoadedStructure;

use super::*;
use crate::sequence::rows::chain_rows;
use crate::sequence::tracks::{TrackContext, TrackData, TrackProvider};

fn loaded(path: &str) -> LoadedStructure {
    let full = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../fixtures")
        .join(path);
    from_structure(vv_io::load(full).unwrap())
}

fn from_structure(structure: vv_core::Structure) -> LoadedStructure {
    let bonds = vv_core::bonds::perceive(&structure.topology, structure.frame(0).positions());
    LoadedStructure::new(structure, None, "test".into(), bonds)
}

fn run(provider: &dyn TrackProvider, loaded: &LoadedStructure) -> TrackData {
    let coords = loaded.structure.frame(0);
    let rows = chain_rows(&loaded.structure.topology);
    provider.compute(&TrackContext {
        loaded,
        positions: coords.positions(),
        rows: &rows,
    })
}

fn marked(t: &TrackData, n: usize) -> Vec<u32> {
    (0..n as u32).filter(|&r| t.kind(r) != 0).collect()
}

/// A chain of CA-only residues named by `names`, numbered from 1.
fn chain_of(names: &[&str]) -> LoadedStructure {
    let text: String = names
        .iter()
        .enumerate()
        .map(|(i, name)| {
            format!(
                "ATOM  {:>5}  CA  {name:<3} A{:>4}    {:8.3}{:8.3}{:8.3}  1.00  0.00           C\n",
                i + 1,
                i + 1,
                (i % 100) as f32 * 3.8,
                (i / 100) as f32 * 3.8,
                0.0
            )
        })
        .collect();
    from_structure(vv_io::pdb::parse(text.as_bytes()).unwrap())
}

fn residues(l: &LoadedStructure) -> usize {
    l.structure.topology.residue_count()
}

fn kind_at(t: &TrackData, r: u32) -> u8 {
    t.kind(r)
}

#[test]
fn crambin_disulfides_pair_six_cysteines() {
    let l = loaded("small/1CRN.pdb");
    let t = run(&Disulfides, &l);
    let top = &l.structure.topology;
    let cys = marked(&t, residues(&l));
    assert_eq!(cys.len(), 6);
    assert!(cys.iter().all(|&r| top.residue_name(r as usize) == "CYS"));
    let mut colors: Vec<u8> = cys.iter().map(|&r| t.kind(r)).collect();
    colors.sort_unstable();
    colors.dedup();
    assert_eq!(colors.len(), 3, "one color per bridge");
}

#[test]
fn hemoglobin_has_no_disulfides_and_marks_its_heme_pocket() {
    let l = loaded("small/4HHB.pdb");
    assert!(marked(&run(&Disulfides, &l), residues(&l)).is_empty());
    let site = run(&LigandSite, &l);
    let top = &l.structure.topology;
    let his87 = (0..residues(&l))
        .find(|&r| top.residues[r].auth_seq_id == 87 && top.residue_name(r) == "HIS")
        .unwrap();
    assert_eq!(site.kind(his87 as u32), 1);
    assert!(site.describe(his87 as u32).unwrap().contains("HEM"));
}

#[test]
fn hemoglobin_interface_and_secondary_structure_are_marked() {
    let l = loaded("small/4HHB.pdb");
    let iface = marked(&run(&Interface, &l), residues(&l));
    assert!(iface.len() > 50);
    let ss = run(&SecondaryStructure, &l);
    assert!(marked(&ss, residues(&l)).len() > 300);
    assert!((0..residues(&l) as u32).any(|r| ss.kind(r) == 1));
}

#[test]
fn spike_glycans_and_sequons_are_marked_apart() {
    let l = loaded("glycan/6X3Z.pdb");
    let t = run(&Glycans, &l);
    let n = residues(&l) as u32;
    let attached = (0..n).filter(|&r| t.kind(r) == 3).count();
    let open = (0..n).filter(|&r| t.kind(r) == 1).count();
    assert!(attached >= 1, "glycosylated Asn");
    assert!(open + attached >= 1);
    let r = (0..n).find(|&r| t.kind(r) == 3).unwrap();
    assert!(t.describe(r).unwrap().contains("N-glycosylated"));
}

#[test]
fn missing_residues_mark_the_edge_of_a_gap() {
    let l = loaded("glycan/6X3Z.pdb");
    let t = run(&Missing, &l);
    assert!(!marked(&t, residues(&l)).is_empty());
    let first = marked(&t, residues(&l))[0];
    assert!(t.describe(first).unwrap().contains("not modeled"));
    let full = loaded("small/1CRN.pdb");
    assert!(marked(&run(&Missing, &full), residues(&full)).is_empty());
}

#[test]
fn numbering_ticks_every_tenth_and_the_first() {
    let l = loaded("small/1CRN.pdb");
    let t = run(&Numbering, &l);
    let top = &l.structure.topology;
    let ticks: Vec<i32> = marked(&t, residues(&l))
        .iter()
        .map(|&r| top.residues[r as usize].auth_seq_id)
        .collect();
    assert_eq!(ticks, [1, 10, 20, 30, 40]);
}

#[test]
fn an_antibody_like_chain_shows_its_liabilities() {
    // Q at the N-terminus, NG, DG, DP, RGD, an unpaired C and a Met.
    let l = chain_of(&[
        "GLN", "VAL", "ASN", "GLY", "ALA", "ASP", "GLY", "ASP", "PRO", "ARG", "GLY", "ASP", "CYS",
        "MET", "ALA",
    ]);
    let t = run(&Liabilities, &l);
    let text = |r: u32| t.describe(r).unwrap_or_default();
    assert!(text(0).contains("N-term Gln"));
    assert!(text(2).contains("NG: Deamidation"));
    assert!(text(5).contains("DG: Isomerization"));
    assert!(text(7).contains("DP: Fragmentation"));
    assert!(text(10).contains("RGD: Integrin"));
    assert!(text(12).contains("Cys: Free cysteine"));
    assert!(text(13).contains("Met: Oxidation"));
    assert_eq!(kind_at(&t, 1), 0, "plain valine");
    assert_eq!(text(1), "");
}

#[test]
fn a_cysteine_in_a_disulfide_is_not_a_liability() {
    let l = loaded("small/1CRN.pdb");
    let t = run(&Liabilities, &l);
    let top = &l.structure.topology;
    let cys3 = (0..residues(&l))
        .find(|&r| top.residues[r].auth_seq_id == 3)
        .unwrap();
    assert_eq!(t.kind(cys3 as u32), 0);
}

#[test]
fn a_sequon_with_no_glycan_is_flagged_and_proline_blocks_it() {
    let l = chain_of(&[
        "ALA", "ASN", "VAL", "THR", "ALA", "ASN", "PRO", "SER", "ALA",
    ]);
    let t = run(&Glycans, &l);
    assert_eq!(t.kind(1), 1);
    assert!(t.describe(1).unwrap().contains("N-X-S/T"));
    assert_eq!(t.kind(5), 0);
}

#[test]
fn alternate_locations_are_marked() {
    let l = loaded("small/1AKE.pdb");
    let t = run(&AltLocs, &l);
    let n = marked(&t, residues(&l));
    assert!(!n.is_empty());
    assert!(t
        .describe(n[0])
        .unwrap()
        .starts_with("Alternate locations: A"));
}

#[test]
fn modified_residues_are_marked() {
    let l = chain_of(&["ALA", "MSE", "GLY"]);
    let t = run(&Modified, &l);
    assert_eq!(marked(&t, residues(&l)), [1]);
}

#[test]
fn a_nine_thousand_residue_chain_annotates_quickly() {
    let cycle = [
        "ALA", "ASN", "GLY", "SER", "ASP", "MET", "CYS", "LEU", "VAL", "THR",
    ];
    let names: Vec<&str> = cycle.iter().copied().cycle().take(9_000).collect();
    let l = chain_of(&names);
    let start = std::time::Instant::now();
    for p in crate::sequence::tracks::PROVIDERS {
        run(*p, &l);
    }
    assert!(start.elapsed().as_secs_f32() < 5.0, "{:?}", start.elapsed());
}
