//! Deuterium, long atom names, segment ids, polymer status from the file's
//! own records, and alternate-location display, through the readers and
//! writers. Real entries (downloaded, not committed) are `#[ignore]`d at
//! the bottom.

use std::path::PathBuf;

use vv_core::altloc::{visible_atoms, AltlocPolicy};
use vv_core::{flags, Element, ResidueClass, Structure};
use vv_io::{Format, SaveOptions};

fn small(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../fixtures/small")
        .join(name)
}

fn load(name: &str) -> Structure {
    let s = vv_io::load(small(name)).unwrap_or_else(|e| panic!("{name}: {e}"));
    s.topology.validate().unwrap();
    s
}

fn cif(text: &str) -> Structure {
    let s = vv_io::mmcif::parse(text.as_bytes()).unwrap();
    s.topology.validate().unwrap();
    s
}

static NEXT_FILE: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);

/// Writes `s` in `ext`'s format, reads it back, and returns it with the
/// writer's warnings.
fn round_trip(s: &Structure, ext: &str) -> (Structure, Vec<String>) {
    let frames: Vec<usize> = (0..s.frame_count()).collect();
    let n = NEXT_FILE.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    let path = std::env::temp_dir().join(format!("vv_identity_{}_{n}.{ext}", std::process::id()));
    let opts = SaveOptions {
        atoms: None,
        frames: &frames,
    };
    let warnings = vv_io::save(s, &path, &opts).unwrap();
    let back = vv_io::load(&path).unwrap();
    std::fs::remove_file(&path).ok();
    back.topology.validate().unwrap();
    (back, warnings)
}

fn pdb_text(s: &Structure) -> String {
    let mut out = Vec::new();
    vv_io::pdb_write::write(s, None, &[0], &mut out).unwrap();
    String::from_utf8(out).unwrap()
}

const ATOM_SITE: &str = "loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
";

/// One `_atom_site` row for the entity/name tests; `x` keeps atoms apart.
fn row(
    group: &str,
    id: u32,
    element: &str,
    name: &str,
    comp: &str,
    asym: &str,
    entity: u32,
) -> String {
    format!("{group} {id} {element} {name} . {comp} {asym} {entity} 1 {id}.0 0.0 0.0\n")
}

// ---- deuterium -------------------------------------------------------------

const HEAVY_WATER: &str =
    "HETATM    1  O   HOH A   1       0.000   0.000   0.000  1.00  0.00           O
HETATM    2  D1  HOH A   1       0.960   0.000   0.000  1.00  0.00           D
HETATM    3  H2  HOH A   1      -0.240   0.930   0.000  1.00  0.00           H
END
";

#[test]
fn pdb_element_d_is_kept_as_deuterium_but_is_hydrogen_everywhere_else() {
    let s = vv_io::pdb::parse(HEAVY_WATER.as_bytes()).unwrap();
    let t = &s.topology;
    assert!(t.is_deuterium(1) && !t.is_deuterium(0) && !t.is_deuterium(2));
    assert_eq!(t.element[1], Element::HYDROGEN);
    assert_eq!(t.flags[1], flags::HETERO | flags::DEUTERIUM);
    let positions = s.frame(0);
    let hydrogens = vv_core::select(t, positions.positions(), "hydrogen").unwrap();
    assert_eq!(hydrogens.ones().collect::<Vec<_>>(), vec![1, 2]);
    assert_eq!(t.element[1].vdw_radius(), Element::HYDROGEN.vdw_radius());
    let bonds = vv_core::bonds::perceive(t, positions.positions());
    assert_eq!(bonds.pairs.len(), 2, "both hydrogens bond to the oxygen");
}

#[test]
fn deuterium_survives_both_writers() {
    let s = vv_io::pdb::parse(HEAVY_WATER.as_bytes()).unwrap();
    assert!(pdb_text(&s)
        .lines()
        .any(|l| l.len() >= 78 && &l[76..78] == " D"));
    for ext in ["pdb", "cif"] {
        let (back, _) = round_trip(&s, ext);
        assert_eq!(back.topology.flags, s.topology.flags, "{ext}");
        assert_eq!(back.topology.element, s.topology.element, "{ext}");
    }
}

#[test]
fn mmcif_type_symbol_d_is_deuterium() {
    let text = format!(
        "data_x\n{ATOM_SITE}{}{}",
        row("HETATM", 1, "O", "O", "HOH", "A", 1),
        row("HETATM", 2, "D", "D1", "HOH", "A", 1)
    );
    let t = cif(&text).topology;
    assert!(t.is_deuterium(1) && !t.is_deuterium(0));
    assert_eq!(t.element[1], Element::HYDROGEN);
}

// ---- atom names longer than four characters ---------------------------------

fn long_name_ligand() -> Structure {
    let text = format!(
        "data_x\n{ATOM_SITE}{}{}{}{}",
        row("HETATM", 1, "C", "C1A1", "LIG", "A", 1),
        row("HETATM", 2, "C", "C1A1X", "LIG", "A", 1),
        row("HETATM", 3, "C", "C1A1Y", "LIG", "A", 1),
        row("HETATM", 4, "N", "N123456", "LIG", "A", 1),
    );
    cif(&text)
}

#[test]
fn mmcif_keeps_long_atom_names_whole_and_writes_them_back() {
    let s = long_name_ligand();
    let names: Vec<&str> = (0..4).map(|a| s.topology.atom_name(a)).collect();
    assert_eq!(names, ["C1A1", "C1A1X", "C1A1Y", "N123456"]);
    let (back, warnings) = round_trip(&s, "cif");
    assert!(warnings.is_empty());
    let back_names: Vec<&str> = (0..4).map(|a| back.topology.atom_name(a)).collect();
    assert_eq!(back_names, names);
}

#[test]
fn pdb_truncates_long_names_uniquely_and_says_so() {
    let (back, warnings) = round_trip(&long_name_ligand(), "pdb");
    assert_eq!(warnings.len(), 1, "{warnings:?}");
    assert!(warnings[0].contains("3 atom name(s)"), "{warnings:?}");
    let names: Vec<&str> = (0..4).map(|a| back.topology.atom_name(a)).collect();
    assert_eq!(names, ["C1A1", "C1A2", "C1A3", "N123"]);
}

// ---- segment ids ------------------------------------------------------------

#[test]
fn segid_splits_a_repeated_chain_letter_and_names_a_blank_chain() {
    let s = load("md_segments.pdb");
    let t = &s.topology;
    let chains: Vec<(&str, &str)> = (0..t.chain_count())
        .map(|c| (t.chain_name(c), t.segid(c)))
        .collect();
    assert_eq!(
        chains,
        [
            ("A", "PROA"),
            ("A", "PROB"),
            ("IONS", "IONS"),
            ("SOLV", "SOLV")
        ]
    );
    let index = t.chain_name_index();
    assert_eq!(index[0], index[1], "same letter, same colour group");
}

#[test]
fn segname_selects_by_segment_id() {
    let s = load("md_segments.pdb");
    let count = |expr: &str| {
        vv_core::select(&s.topology, s.frame(0).positions(), expr)
            .unwrap()
            .count_ones(..)
    };
    assert_eq!(count("segname PROB"), 3);
    assert_eq!(count("segname PROA PROB"), 11);
    assert_eq!(count("segid SOLV and hydrogen"), 4);
    assert_eq!(count("segname NOPE"), 0);
}

#[test]
fn segids_round_trip_through_pdb_and_mmcif() {
    let s = load("md_segments.pdb");
    for ext in ["pdb", "cif"] {
        let (back, _) = round_trip(&s, ext);
        let t = &back.topology;
        assert_eq!(t.chain_count(), s.topology.chain_count(), "{ext}");
        for c in 0..t.chain_count() {
            assert_eq!(t.segid(c), s.topology.segid(c), "{ext} chain {c}");
        }
    }
}

#[test]
fn a_structure_without_segids_writes_no_segid_column() {
    let (mut out, s) = (Vec::new(), load("1CRN.pdb"));
    vv_io::mmcif_write::write(&s, None, &[0], &mut out).unwrap();
    assert!(!String::from_utf8(out).unwrap().contains("vizviz_segid"));
}

// ---- polymer status from the file's own records -------------------------------

const ENTITIES: &str = "loop_
_entity.id
_entity.type
1 polymer
2 non-polymer
3 water
#
loop_
_entity_poly.entity_id
_entity_poly.type
1 'polypeptide(L)'
#
";

fn entity_structure(tables: &str) -> Structure {
    let mut text = format!("data_x\n{tables}{ATOM_SITE}");
    text += &row("ATOM", 1, "N", "N", "ALA", "A", 1);
    text += &row("ATOM", 2, "C", "CA", "ALA", "A", 1);
    text += &row("HETATM", 3, "N", "N", "ZZZ", "A", 1);
    text += &row("HETATM", 4, "C", "CA", "ZZZ", "A", 1);
    text += &row("HETATM", 5, "N", "N", "ALA", "B", 2);
    text += &row("HETATM", 6, "C", "CA", "ALA", "B", 2);
    text += &row("HETATM", 7, "O", "O", "HOH", "C", 3);
    cif(&text)
}

fn classes(s: &Structure) -> Vec<ResidueClass> {
    (0..s.topology.residue_count())
        .map(|r| s.topology.residue_class(r))
        .collect()
}

#[test]
fn entity_type_overrides_the_name_tables_where_they_disagree() {
    use ResidueClass::{Protein, SmallMolecule, Water};
    let with = entity_structure(ENTITIES);
    assert_eq!(classes(&with), [Protein, Protein, SmallMolecule, Water]);
    let without = entity_structure("");
    let by_name = classes(&without);
    assert_eq!(
        by_name[2], Protein,
        "the name table alone calls a free ALA protein"
    );
    assert_ne!(by_name[1], Protein, "and knows nothing of ZZZ");
}

#[test]
fn a_polymer_of_unknown_type_uses_each_monomers_chemical_component_type() {
    let tables = "loop_
_entity.id
_entity.type
1 polymer
2 non-polymer
3 water
loop_
_entity_poly.entity_id
_entity_poly.type
1 other
loop_
_chem_comp.id
_chem_comp.type
ZZZ 'L-peptide linking'
ALA 'L-peptide linking'
";
    let s = entity_structure(tables);
    assert_eq!(s.topology.residue_class(1), ResidueClass::Protein);
}

#[test]
fn real_4hhb_entities_and_seqres_agree_with_each_other() {
    let (c, p) = (load("4HHB.cif"), load("4HHB.pdb"));
    assert!(!c.topology.polymer_hint.is_empty() && !p.topology.polymer_hint.is_empty());
    assert_eq!(c.topology.class_counts, p.topology.class_counts);
    let hem = (0..c.topology.residue_count())
        .find(|&r| c.topology.residue_name(r) == "HEM")
        .unwrap();
    assert_eq!(c.topology.residue_class(hem), ResidueClass::SmallMolecule);
}

#[test]
fn seqres_marks_a_modified_residue_polymer_and_a_ligand_not() {
    let pdb = "SEQRES   1 A    3  ALA MSE GLY
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       1.000   0.000   0.000  1.00  0.00           C
HETATM    3  N   MSE A   2       2.000   0.000   0.000  1.00  0.00           N
HETATM    4  CA  MSE A   2       3.000   0.000   0.000  1.00  0.00           C
HETATM    5  N   ALA A 101       4.000   0.000   0.000  1.00  0.00           N
HETATM    6  CA  ALA A 101       5.000   0.000   0.000  1.00  0.00           C
END
";
    let s = vv_io::pdb::parse(pdb.as_bytes()).unwrap();
    assert_eq!(
        classes(&s),
        [
            ResidueClass::Protein,
            ResidueClass::Protein,
            ResidueClass::SmallMolecule
        ]
    );
}

// ---- alternate locations ------------------------------------------------------

#[test]
fn one_conformer_per_residue_is_displayed_by_default_and_all_stay_in_the_data() {
    let s = load("1AKE.pdb");
    let t = &s.topology;
    let alt_atoms = t.alt_loc.iter().filter(|&&a| a != 0).count();
    assert!(alt_atoms > 0);
    let shown = visible_atoms(t, AltlocPolicy::First).expect("1AKE has conformers");
    let hidden = t.atom_count() - shown.count_ones(..);
    assert!(hidden > 0 && hidden < alt_atoms);
    for res in &t.residues {
        let labels: std::collections::BTreeSet<u8> = res
            .atoms
            .clone()
            .filter(|&a| shown.contains(a as usize))
            .map(|a| t.alt_loc[a as usize])
            .filter(|&l| l != 0)
            .collect();
        assert!(labels.len() <= 1, "{labels:?}");
    }
    assert_eq!(visible_atoms(t, AltlocPolicy::All), None);
}

// ---- real entries (downloaded, not committed) -----------------------------------

fn real(name: &str) -> Option<PathBuf> {
    let p = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../fixtures/real")
        .join(name);
    p.exists().then_some(p)
}

fn load_real(name: &str) -> Structure {
    let path = real(name).unwrap_or_else(|| panic!("fixtures/real/{name} is not downloaded"));
    vv_io::load(path).unwrap()
}

#[test]
#[ignore = "needs fixtures/real/{3KCJ,5A93}.{pdb,cif,bcif}"]
fn neutron_entries_keep_deuterium_in_every_format_and_through_the_writers() {
    for id in ["3KCJ", "5A93"] {
        let (p, c, b) = (
            load_real(&format!("{id}.pdb")),
            load_real(&format!("{id}.cif")),
            load_real(&format!("{id}.bcif")),
        );
        let count = |s: &Structure| {
            (0..s.topology.atom_count())
                .filter(|&a| s.topology.is_deuterium(a))
                .count()
        };
        assert!(count(&p) > 1000, "{id}");
        assert_eq!((count(&p), count(&b)), (count(&c), count(&c)), "{id}");
        for ext in ["pdb", "cif"] {
            assert_eq!(count(&round_trip(&p, ext).0), count(&p), "{id} {ext}");
        }
    }
}

#[test]
#[ignore = "needs fixtures/real/3NIR.{pdb,cif,bcif}"]
fn real_alternate_locations_agree_across_formats_and_display_one_conformer() {
    let (p, c, b) = (
        load_real("3NIR.pdb"),
        load_real("3NIR.cif"),
        load_real("3NIR.bcif"),
    );
    assert_eq!(p.topology.alt_loc, c.topology.alt_loc);
    assert_eq!(b.topology.alt_loc, c.topology.alt_loc);
    let shown = visible_atoms(&p.topology, AltlocPolicy::First).unwrap();
    assert!(shown.count_ones(..) < p.topology.atom_count());
    let by_label = visible_atoms(&p.topology, AltlocPolicy::Label(b'B')).unwrap();
    assert_ne!(shown, by_label);
}

#[test]
#[ignore = "needs fixtures/real/*.bcif next to the .cif files"]
fn real_binary_cif_matches_mmcif_atom_for_atom() {
    for id in ["3KCJ", "5A93", "3NIR"] {
        let (b, c) = (
            load_real(&format!("{id}.bcif")),
            load_real(&format!("{id}.cif")),
        );
        let (tb, tc) = (&b.topology, &c.topology);
        assert_eq!(tb.atom_count(), tc.atom_count(), "{id}");
        assert_eq!(tb.name, tc.name, "{id}");
        assert_eq!(tb.element, tc.element, "{id}");
        assert_eq!(tb.flags, tc.flags, "{id}");
        assert_eq!(tb.alt_loc, tc.alt_loc, "{id}");
        assert_eq!(tb.chain_count(), tc.chain_count(), "{id}");
        for (pb, pc) in b.frame(0).positions().iter().zip(c.frame(0).positions()) {
            assert!((*pb - *pc).abs().max_element() < 1e-3, "{id}");
        }
    }
}

#[test]
fn format_is_picked_from_the_bcif_extension() {
    assert_eq!(Format::from_path(&small("x.bcif")), Some(Format::Bcif));
}
