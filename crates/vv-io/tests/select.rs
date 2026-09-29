//! Selection language checked against a real entry: human hemoglobin 4HHB
//! (chains A-D, four HEM groups, waters).

use std::path::PathBuf;

use vv_core::fixedbitset::FixedBitSet;
use vv_core::{select, Element, Structure};

fn hemoglobin() -> Structure {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/small/4HHB.cif");
    let s = vv_io::load(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    s.topology.validate().unwrap();
    s
}

fn eval(s: &Structure, expr: &str) -> FixedBitSet {
    select(&s.topology, s.frame(0).positions(), expr).unwrap_or_else(|e| panic!("{expr}: {e}"))
}

fn count(s: &Structure, expr: &str) -> usize {
    eval(s, expr).count_ones(..)
}

/// Atoms in residues with the given name, counted straight from the topology.
fn atoms_named(s: &Structure, comp: &str) -> usize {
    let t = &s.topology;
    (0..t.residue_count())
        .filter(|&r| t.residue_name(r) == comp)
        .map(|r| t.residues[r].atoms.len())
        .sum()
}

#[test]
fn counts_match_the_entry() {
    let s = hemoglobin();
    assert_eq!(count(&s, "all"), 4779);
    assert_eq!(count(&s, "none"), 0);
    assert_eq!(count(&s, "chain A and name CA"), 141);
    assert_eq!(count(&s, "chain B and name CA"), 146);
    assert_eq!(count(&s, "resname HEM"), 172);
    assert_eq!(count(&s, "element FE"), 4);
    assert_eq!(count(&s, "index 0-9"), 10);
    assert_eq!(count(&s, "resid 1-10 and chain A and name CA"), 10);

    // Everything that is not a standard amino acid: 4 HEM, 221 HOH, 2 PO4.
    let hem = atoms_named(&s, "HEM");
    let water = atoms_named(&s, "HOH");
    let phosphate = atoms_named(&s, "PO4");
    assert_eq!((hem, water, phosphate), (172, 221, 2));
    assert_eq!(count(&s, "protein"), 4779 - hem - water - phosphate);

    let t = &s.topology;
    let water_residues = (0..t.residue_count())
        .filter(|&r| t.residue_name(r) == "HOH")
        .count();
    assert_eq!(count(&s, "water"), water_residues, "one atom per water");
    assert_eq!(eval(&s, "hetero"), eval(&s, "not protein"));
}

#[test]
fn within_finds_the_iron_coordination_sphere() {
    let s = hemoglobin();
    let sphere = eval(&s, "within 3 of element FE");
    assert!(sphere.count_ones(..) >= 20, "{}", sphere.count_ones(..));
    let iron = Element::from_symbol(b"Fe");
    for (i, e) in s.topology.element.iter().enumerate() {
        if *e == iron {
            assert!(sphere.contains(i), "Fe {i} is inside its own sphere");
        }
    }
    // Brute force on the first iron: every atom within 3 Å is selected.
    let coords = s.frame(0);
    let positions = coords.positions();
    let fe = s.topology.element.iter().position(|e| *e == iron).unwrap();
    for (i, p) in positions.iter().enumerate() {
        if p.distance(positions[fe]) <= 3.0 {
            assert!(sphere.contains(i), "atom {i}");
        }
    }
}

#[test]
fn byres_expands_the_heme_pocket_to_whole_residues() {
    let s = hemoglobin();
    let t = &s.topology;
    let pocket = eval(&s, "within 4 of resname HEM and protein");
    let expanded = eval(&s, "byres within 4 of resname HEM and protein");
    assert!(pocket.count_ones(..) > 0);
    assert!(pocket.is_subset(&expanded));
    assert!(expanded.count_ones(..) > pocket.count_ones(..));
    for r in &t.residues {
        let range = r.atoms.start as usize..r.atoms.end as usize;
        let selected = expanded.count_ones(range.clone());
        assert!(
            selected == 0 || selected == range.len(),
            "residue partially selected"
        );
    }
}

#[test]
fn errors_point_at_the_offending_token() {
    let e = vv_core::select::parse("chain A and resnam HEM").unwrap_err();
    assert_eq!(e.span, 12..18);
    assert!(e.message.contains("resnam"));
}
