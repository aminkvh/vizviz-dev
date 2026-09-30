//! `vv_core::interactions` against real structures.

use vv_core::interactions::{hydrogen_bonds, metal_coordination, salt_bridges};
use vv_core::{bonds, Structure};

fn fixture(name: &str) -> Structure {
    let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../fixtures/small")
        .join(name);
    vv_io::load(path).unwrap()
}

fn label(s: &Structure, atom: u32) -> String {
    let t = &s.topology;
    let r = t.residue_index[atom as usize] as usize;
    format!(
        "{} {} {}",
        t.residue_name(r),
        t.residues[r].auth_seq_id,
        t.atom_name(atom as usize)
    )
}

#[test]
fn every_heme_iron_is_bound_to_its_proximal_histidine() {
    let s = fixture("4HHB.cif");
    let coords = s.frame(0);
    let contacts = metal_coordination(&s.topology, coords.positions(), &|_| true);
    let irons: Vec<&[u32; 2]> = contacts
        .iter()
        .filter(|[m, _]| s.topology.element[*m as usize].atomic_number() == 26)
        .collect();
    assert!(irons.len() >= 4, "{}", irons.len());
    let his: Vec<String> = irons.iter().map(|[_, d]| label(&s, *d)).collect();
    // Proximal (F8) histidine: 87 in the alpha chains, 92 in the beta.
    for expected in ["HIS 87 NE2", "HIS 92 NE2"] {
        assert_eq!(his.iter().filter(|h| **h == expected).count(), 2, "{his:?}");
    }
}

#[test]
fn crambins_backbone_hydrogen_bonds_are_in_the_dssp_range() {
    let s = fixture("1CRN.cif");
    let coords = s.frame(0);
    let positions = coords.positions();
    let table = bonds::perceive(&s.topology, positions);
    let found = hydrogen_bonds(&s.topology, &table, positions, &|_| true);
    let backbone = found
        .iter()
        .filter(|[a, b]| {
            let is = |x: u32, n: &str| s.topology.atom_name(x as usize) == n;
            (is(*a, "N") && is(*b, "O")) || (is(*a, "O") && is(*b, "N"))
        })
        .filter(|[a, b]| {
            let r = |x: u32| s.topology.residue_index[x as usize] as i64;
            (r(*a) - r(*b)).abs() >= 3
        })
        .count();
    // 41 here; DSSP finds a few dozen backbone N-H...O=C bonds in this 46-residue
    // protein (2 helices, a small sheet); the geometric criterion is looser.
    assert!((25..=70).contains(&backbone), "{backbone}");
}

#[test]
fn hemoglobin_has_salt_bridges_and_no_bridge_pairs_a_residue_with_itself() {
    let s = fixture("4HHB.cif");
    let coords = s.frame(0);
    let bridges = salt_bridges(&s.topology, coords.positions(), &|_| true);
    assert!(bridges.len() > 20, "{}", bridges.len());
    for [acid, base] in bridges {
        assert_ne!(
            s.topology.residue_index[acid as usize],
            s.topology.residue_index[base as usize]
        );
    }
}
