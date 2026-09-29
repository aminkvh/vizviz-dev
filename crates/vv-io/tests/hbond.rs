//! `vv_core::hbond::hydrogen_bonds_into` against real structures.
//!
//! The real cross-check here is against `vv_core::dssp`: two
//! independently built kernels (a geometric Baker & Hubbard 1984 model
//! here, an electrostatic Kabsch & Sander 1983 model there) computing
//! different things should still substantially agree on the same real
//! backbone N-H...O=C pairs inside a helix. Substantial agreement, not
//! exact: DSSP's electrostatic threshold and this kernel's geometric one
//! are different criteria by design (module docs), so some helix turns
//! satisfy one and not the other.

use std::collections::HashSet;

use vv_core::dssp::{assign, DsspCode};
use vv_core::hbond::hydrogen_bonds_into;
use vv_core::{bonds, Structure};

fn fixture(name: &str) -> Structure {
    let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../fixtures/small")
        .join(name);
    vv_io::load(path).unwrap()
}

fn atom(s: &Structure, r: usize, name: &str) -> Option<u32> {
    s.topology.residues[r]
        .atoms
        .clone()
        .find(|&a| s.topology.atom_name(a as usize) == name)
}

#[test]
fn finds_a_reasonable_number_of_bonds_in_crambin() {
    let s = fixture("1CRN.cif");
    let coords = s.frame(0);
    let positions = coords.positions();
    let bond_table = bonds::perceive(&s.topology, positions);
    let mut out = Vec::new();
    hydrogen_bonds_into(&s.topology.element, &bond_table, positions, &mut out);
    // 46 residues, a real mix of helix, sheet and loop, plus waters: a
    // healthy double-digit count is expected, not zero and not the
    // thousands a broken angle check (finding everything) would give.
    assert!(out.len() > 15 && out.len() < 300, "{}", out.len());
}

#[test]
fn agrees_substantially_with_dssps_backbone_helix_turns() {
    for name in ["1CRN.cif", "1UBQ.cif", "4HHB.cif", "1AKE.cif"] {
        let s = fixture(name);
        let coords = s.frame(0);
        let positions = coords.positions();
        let codes = assign(&s.topology, positions);
        let bond_table = bonds::perceive(&s.topology, positions);
        let mut hbonds = Vec::new();
        hydrogen_bonds_into(&s.topology.element, &bond_table, positions, &mut hbonds);
        let pairs: HashSet<(u32, u32)> = hbonds
            .iter()
            .flat_map(|hb| [(hb.a, hb.b), (hb.b, hb.a)])
            .collect();

        // Every alpha-helix i, i+4 turn's backbone N(i+4)...O(i): DSSP
        // found an electrostatic bond there by construction. Does the
        // independent geometric kernel see the same pair?
        let mut checked = 0;
        let mut agreed = 0;
        for i in 0..codes.len().saturating_sub(4) {
            if codes[i] != DsspCode::AlphaHelix || codes[i + 4] != DsspCode::AlphaHelix {
                continue;
            }
            let (Some(o), Some(n)) = (atom(&s, i, "O"), atom(&s, i + 4, "N")) else {
                continue;
            };
            checked += 1;
            if pairs.contains(&(o, n)) {
                agreed += 1;
            }
        }
        assert!(
            checked > 5,
            "{name}: too few helix turns to mean anything: {checked}"
        );
        let pct = 100.0 * agreed as f32 / checked as f32;
        eprintln!("{name}: {agreed}/{checked} ({pct:.0}%) DSSP helix turns also found by the geometric kernel");
        assert!(
            pct > 60.0,
            "{name}: only {agreed}/{checked} ({pct:.0}%) DSSP helix turns also found by the geometric kernel"
        );
    }
}
