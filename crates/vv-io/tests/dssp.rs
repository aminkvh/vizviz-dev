//! `vv_core::dssp::assign` against real structures.
//!
//! This is a smoke test, not an accuracy gate: file `HELIX`/`SHEET`
//! records (PDB) and `_struct_conf` (mmCIF) are frequently
//! author-assigned rather than DSSP-derived, both parsers collapse
//! everything to three states (so H vs G vs I can never be checked this
//! way), and both assign secondary structure by residue *range*, so
//! end-cap residues disagree with any independent, H-bond-based
//! assignment as a matter of course. Real DSSP parity needs real DSSP
//! (`mkdssp`) output on these fixtures; that needs a C++20 toolchain this
//! environment does not have configured and is not attempted here.

use vv_core::dssp::{assign, DsspCode};
use vv_core::topology::SecondaryStructure as Ss;
use vv_core::Structure;

fn fixture(name: &str) -> Structure {
    let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../fixtures/small")
        .join(name);
    vv_io::load(path).unwrap()
}

fn agreement_with_file_records(s: &Structure) -> (usize, usize) {
    let codes = assign(&s.topology, s.frame(0).positions());
    let mut total = 0;
    let mut agree = 0;
    for (i, rec) in s.topology.residues.iter().enumerate() {
        if rec.ss == Ss::Unknown || codes[i] == DsspCode::None {
            continue;
        }
        total += 1;
        if codes[i].simplified() == rec.ss {
            agree += 1;
        }
    }
    (agree, total)
}

#[test]
fn agrees_broadly_with_file_provided_records() {
    for name in ["1CRN.cif", "1UBQ.cif", "4HHB.cif", "1AKE.cif"] {
        let s = fixture(name);
        let (agree, total) = agreement_with_file_records(&s);
        assert!(
            total > 0,
            "{name}: no residue has both a file SS record and a computed one"
        );
        let pct = 100.0 * agree as f32 / total as f32;
        assert!(
            pct > 70.0,
            "{name}: only {agree}/{total} ({pct:.0}%) residues agree with file SS records"
        );
        eprintln!("{name}: {agree}/{total} ({pct:.0}%) agree with file SS records");
    }
}

/// Crambin (1CRN) residues 7-19 (1-based) carry a HELIX record
/// (`analysis.rs`'s `crambin_helix_has_negative_phi_and_psi` pins the same
/// stretch by phi/psi). An independent, hydrogen-bond-based assignment
/// should find a run of alpha helix somewhere inside it - not necessarily
/// the exact same endpoints DSSP's own bookkeeping would give.
#[test]
fn crambins_helix_is_found_by_hydrogen_bonds() {
    let s = fixture("1CRN.cif");
    let codes = assign(&s.topology, s.frame(0).positions());
    let window = &codes[5..19];
    let run = window
        .windows(4)
        .any(|w| w.iter().all(|&c| c == DsspCode::AlphaHelix));
    assert!(run, "{window:?}");
}

/// Hemoglobin (4HHB) has heme groups and water: non-protein residues get
/// no DSSP code at all, not a default like `Coil`.
#[test]
fn non_protein_residues_get_no_code() {
    let s = fixture("4HHB.cif");
    let codes = assign(&s.topology, s.frame(0).positions());
    let any_none = s.topology.residues.iter().enumerate().any(|(i, _)| {
        s.topology.residue_class(i) != vv_core::ResidueClass::Protein && codes[i] == DsspCode::None
    });
    assert!(any_none, "4HHB has heme/water residues to exercise this");
}

/// Ubiquitin's beta-grasp fold has a prominent five-stranded sheet: this
/// reads out as `Strand`/`Bridge` somewhere, not just helix and coil - a
/// basic check that bridge detection fires at all on real geometry.
/// (Hemoglobin, an all-alpha globin fold, is not a useful check here even
/// though it loads correctly: it is not expected to have any strand.)
#[test]
fn strands_are_found_in_a_beta_sheet_containing_structure() {
    let s = fixture("1UBQ.cif");
    let codes = assign(&s.topology, s.frame(0).positions());
    let any_strand = codes
        .iter()
        .any(|&c| matches!(c, DsspCode::Strand | DsspCode::Bridge));
    assert!(any_strand, "{codes:?}");
}
