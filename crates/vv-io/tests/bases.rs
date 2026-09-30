//! `vv_core::bases` against real nucleic acids: slabs on real ring atoms
//! and base-pair detection.

use vv_core::bases::{bases, connectors, ladder, pairs, push_plate, HALF_THICKNESS};
use vv_core::glam::Vec3;
use vv_core::{PolytopeMesh, Structure};

fn fixture(name: &str) -> Structure {
    let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../fixtures/small")
        .join(name);
    vv_io::load(path).unwrap()
}

fn pair_count(name: &str) -> (usize, usize) {
    let s = fixture(name);
    let coords = s.frame(0);
    let found = bases(&s.topology, coords.positions(), &|_| true);
    let n = pairs(&found, coords.positions()).len();
    (found.len(), n)
}

#[test]
fn the_dickerson_dodecamer_has_twelve_base_pairs() {
    assert_eq!(pair_count("1BNA.cif"), (24, 12));
}

#[test]
fn trna_has_as_many_pairs_as_the_cloverleafs_four_stems() {
    let (bases, pairs) = pair_count("1EHZ.cif");
    assert_eq!(bases, 76);
    // 7+4+5+5 stem pairs; the count matches, the pairs are not compared one by one.
    assert_eq!(pairs, 21);
}

#[test]
fn a_ladder_has_a_stem_per_base_and_a_rung_per_pair() {
    let s = fixture("1BNA.cif");
    let coords = s.frame(0);
    let found = bases(&s.topology, coords.positions(), &|_| true);
    let p = pairs(&found, coords.positions());
    assert_eq!(ladder(&found, &p).len(), 24 + 12);
    assert_eq!(connectors(&found).len(), 24 + 24);
}

#[test]
fn purines_have_nine_ring_atoms_and_pyrimidines_six() {
    let s = fixture("1BNA.cif");
    let coords = s.frame(0);
    for b in bases(&s.topology, coords.positions(), &|_| true) {
        let name = s.topology.residue_name(b.residue as usize);
        let expected = if matches!(name, "DA" | "DG") { 9 } else { 6 };
        assert_eq!(b.ring.len(), expected, "{name}");
        assert_eq!(b.purine, expected == 9);
    }
}

#[test]
fn a_real_plates_vertices_lie_within_the_bases_own_thickness() {
    let s = fixture("1BNA.cif");
    let coords = s.frame(0);
    let positions = coords.positions();
    let all = bases(&s.topology, positions, &|_| true);
    let b = all.iter().find(|b| b.purine).expect("a purine");
    let mut mesh = PolytopeMesh::default();
    push_plate(b, positions, [9, 9, 9], &mut mesh);
    // Every vertex is a ring atom pushed out in-plane and half a
    // thickness along the normal, so it is that far from the plane through
    // the ring centroid; ring atoms of a real base are planar to ~0.1 A.
    let ring: Vec<Vec3> = b.ring.iter().map(|&a| positions[a as usize]).collect();
    let centroid = ring.iter().copied().sum::<Vec3>() / ring.len() as f32;
    let normal = mesh
        .normals
        .iter()
        .copied()
        .find(|n| n.length() > 0.5)
        .unwrap();
    let flat: Vec<f32> = mesh
        .positions
        .iter()
        .zip(&mesh.normals)
        .filter(|(_, n)| n.dot(normal).abs() > 0.99)
        .map(|(p, _)| (*p - centroid).dot(normal).abs())
        .collect();
    assert!(!flat.is_empty());
    for d in flat {
        assert!((d - HALF_THICKNESS).abs() < 0.25, "{d}");
    }
    for &atom in &mesh.source_atom {
        assert!(b.ring.contains(&atom) || atom == b.glycosidic);
    }
}
