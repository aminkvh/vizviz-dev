//! Polymer status through the writers: read, write, read again keeps each
//! residue's class and polymer hint (`_entity*`/`_chem_comp` in mmCIF,
//! `SEQRES` in PDB).

use std::path::PathBuf;

use vv_core::{PolymerHint, ResidueClass, Structure};
use vv_io::SaveOptions;

fn small(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../fixtures/small")
        .join(name)
}

fn load(name: &str) -> Structure {
    vv_io::load(small(name)).unwrap_or_else(|e| panic!("{name}: {e}"))
}

static NEXT: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);

fn round_trip(s: &Structure, ext: &str) -> Structure {
    let n = NEXT.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    let path = std::env::temp_dir().join(format!("vv_polymer_{}_{n}.{ext}", std::process::id()));
    let opts = SaveOptions {
        atoms: None,
        frames: &[0],
    };
    vv_io::save(s, &path, &opts).unwrap();
    let back = vv_io::load(&path).unwrap();
    std::fs::remove_file(&path).ok();
    back.topology.validate().unwrap();
    back
}

fn classes(s: &Structure) -> Vec<ResidueClass> {
    (0..s.topology.residue_count())
        .map(|r| s.topology.residue_class(r))
        .collect()
}

/// A hint with the distinctions the writers cannot keep folded together:
/// a water, a non-polymer and an unstated status all read "not polymer".
fn polymer_family(s: &Structure) -> Vec<PolymerHint> {
    (0..s.topology.residue_count())
        .map(|r| match s.topology.polymer_hint.get(r) {
            Some(&h @ (PolymerHint::Protein | PolymerHint::Nucleic)) => h,
            _ => PolymerHint::NonPolymer,
        })
        .collect()
}

#[test]
fn hints_and_classes_survive_both_writers() {
    for name in ["4HHB", "1AKE", "1CRN"] {
        for ext in ["pdb", "cif"] {
            let s = load(&format!("{name}.{ext}"));
            assert!(!s.topology.polymer_hint.is_empty(), "{name}.{ext}");
            for writer in ["pdb", "cif"] {
                let back = round_trip(&s, writer);
                let label = format!("{name}.{ext} via {writer}");
                assert_eq!(classes(&back), classes(&s), "{label}: classes");
                assert_eq!(polymer_family(&back), polymer_family(&s), "{label}: hints");
            }
        }
    }
}

#[test]
fn mmcif_hints_round_trip_exactly() {
    for name in ["4HHB", "1AKE", "1CRN"] {
        let s = load(&format!("{name}.cif"));
        let back = round_trip(&s, "cif");
        assert_eq!(
            back.topology.polymer_hint, s.topology.polymer_hint,
            "{name}"
        );
    }
}

const FREE_ALA_PDB: &str = "SEQRES   1 A    2  ALA GLY
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C
ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00  0.00           C
ATOM      4  N   GLY A   2       3.332   1.536   0.000  1.00  0.00           N
ATOM      5  CA  GLY A   2       4.000   2.800   0.000  1.00  0.00           C
TER       6      GLY A   2
HETATM    7  N   ALA A 101      10.000  10.000  10.000  1.00  0.00           N
HETATM    8  CA  ALA A 101      11.458  10.000  10.000  1.00  0.00           C
HETATM    9  C   ALA A 101      12.009  11.420  10.000  1.00  0.00           C
HETATM   10  O   HOH A 201      20.000  20.000  20.000  1.00  0.00           O
END
";

#[test]
fn a_free_amino_acid_ligand_stays_a_small_molecule() {
    use ResidueClass::{Protein, SmallMolecule, Water};
    let s = vv_io::pdb::parse(FREE_ALA_PDB.as_bytes()).unwrap();
    assert_eq!(classes(&s), [Protein, Protein, SmallMolecule, Water]);
    for ext in ["pdb", "cif"] {
        let back = round_trip(&s, ext);
        assert_eq!(classes(&back), classes(&s), "{ext}");
        assert_eq!(
            polymer_family(&back),
            [
                PolymerHint::Protein,
                PolymerHint::Protein,
                PolymerHint::NonPolymer,
                PolymerHint::NonPolymer
            ],
            "{ext}"
        );
    }
}

fn cif_text(s: &Structure) -> String {
    let mut out = Vec::new();
    vv_io::mmcif_write::write(s, None, &[0], &mut out).unwrap();
    String::from_utf8(out).unwrap()
}

fn pdb_text(s: &Structure) -> String {
    let mut out = Vec::new();
    vv_io::pdb_write::write(s, None, &[0], &mut out).unwrap();
    String::from_utf8(out).unwrap()
}

#[test]
fn identical_chains_share_an_entity_and_list_their_strands() {
    let text = cif_text(&load("4HHB.cif"));
    let types = "_entity.type\n1 polymer\n2 polymer\n3 non-polymer\n4 non-polymer\n5 water\n";
    assert!(text.contains(types));
    assert!(text.contains("_entity_poly.pdbx_strand_id\n1 'polypeptide(L)' VLSPADKTNV"));
    assert!(text.contains(" A,C\n2 'polypeptide(L)' VHLTPEEKSA"));
    assert!(text.contains(" B,D\n"));
    assert!(text.contains("3 HEM HEM\n4 PO4 PO4\n5 water HOH\n"));
    assert!(text.contains("HEM 'non-polymer'\n"));
    assert!(text.contains("VAL 'L-peptide linking'\n"));
}

#[test]
fn a_one_chain_entry_writes_its_canonical_sequence() {
    let text = cif_text(&load("1CRN.pdb"));
    assert!(text.contains("1 'polypeptide(L)' TTCCPSIVARSNFNVCRLPGTPEAICATYTGCIIIPGATCPGDYAN A\n"));
}

#[test]
fn pdb_seqres_lists_thirteen_residues_per_line_in_the_fixed_columns() {
    let text = pdb_text(&load("1CRN.cif"));
    let lines: Vec<&str> = text.lines().filter(|l| l.starts_with("SEQRES")).collect();
    assert_eq!(lines.len(), 4, "46 residues at 13 per line");
    assert_eq!(
        &lines[0][..19 + 13 * 4 - 1],
        "SEQRES   1 A   46  THR THR CYS CYS PRO SER ILE VAL ALA ARG SER ASN PHE"
    );
    assert_eq!(
        lines[3].trim_end(),
        "SEQRES   4 A   46  CYS PRO GLY ASP TYR ALA ASN"
    );
}

#[test]
fn a_chain_record_mixing_polymer_and_ligand_is_split_by_kind() {
    let s = vv_io::pdb::parse(FREE_ALA_PDB.as_bytes()).unwrap();
    let text = cif_text(&s);
    assert!(text.contains("_entity.type\n1 polymer\n2 non-polymer\n3 water\n"));
    let back = round_trip(&s, "cif");
    assert_eq!(back.topology.chain_count(), 3);
}

const FREE_ALA_OWN_CHAIN: &str = "SEQRES   1 A    2  ALA GLY
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C
ATOM      3  N   GLY A   2       3.332   1.536   0.000  1.00  0.00           N
TER       4      GLY A   2
HETATM    5  N   ALA B   1      10.000  10.000  10.000  1.00  0.00           N
HETATM    6  CA  ALA B   1      11.458  10.000  10.000  1.00  0.00           C
END
";

#[test]
fn a_free_amino_acid_in_a_chain_without_seqres_is_not_polymer() {
    use ResidueClass::{Protein, SmallMolecule};
    let s = vv_io::pdb::parse(FREE_ALA_OWN_CHAIN.as_bytes()).unwrap();
    assert_eq!(classes(&s), [Protein, Protein, SmallMolecule]);
    for ext in ["pdb", "cif"] {
        assert_eq!(classes(&round_trip(&s, ext)), classes(&s), "{ext}");
    }
}
