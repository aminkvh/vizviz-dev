//! `vv_core::cartoon::build` against real structures, using real
//! `vv_core::dssp::assign` output (not hand-built codes, unlike
//! `vv-core`'s own unit tests): the geometry generator's actual input in
//! practice.

use vv_core::cartoon::{build, RIBBON_WIDTH, RING};
use vv_core::dssp::{assign, DsspCode};
use vv_core::Structure;

fn fixture(name: &str) -> Structure {
    let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../fixtures/small")
        .join(name);
    vv_io::load(path).unwrap()
}

#[test]
fn a_real_helix_and_a_real_strand_both_produce_well_formed_geometry() {
    for name in ["1CRN.cif", "1UBQ.cif", "4HHB.cif", "1AKE.cif"] {
        let s = fixture(name);
        let coords = s.frame(0);
        let positions = coords.positions();
        let codes = assign(&s.topology, positions);
        let mesh = build(&s.topology, positions, &codes).expand();

        assert!(!mesh.positions.is_empty(), "{name}: empty mesh");
        assert_eq!(mesh.positions.len(), mesh.normals.len());
        assert_eq!(mesh.positions.len(), mesh.source.len());
        assert_eq!(
            mesh.positions.len() % RING,
            0,
            "{name}: not whole cross-sections"
        );
        assert!(
            mesh.indices
                .iter()
                .all(|&i| (i as usize) < mesh.positions.len()),
            "{name}: an index points past the end of positions"
        );
        assert!(
            mesh.positions.iter().all(|p| p.is_finite()),
            "{name}: a non-finite vertex position"
        );
        assert!(
            mesh.normals.iter().all(|n| n.is_finite()),
            "{name}: a non-finite vertex normal"
        );
        assert!(
            codes.contains(&DsspCode::AlphaHelix)
                || codes
                    .iter()
                    .any(|&c| matches!(c, DsspCode::Strand | DsspCode::Bridge)),
            "{name}: has no helix or strand to exercise the ribbon path at all"
        );
    }
}

/// A structure with both real helix and real strand (crambin has both,
/// per its file's own HELIX/SHEET records): the wide ribbon shape should
/// actually show up wider than a coil region's, on real geometry, not
/// just the synthetic circle `vv-core`'s own unit test uses.
#[test]
fn a_real_helix_reads_wider_than_a_real_coil_region() {
    let s = fixture("1CRN.cif");
    let coords = s.frame(0);
    let positions = coords.positions();
    let codes = assign(&s.topology, positions);
    let mesh = build(&s.topology, positions, &codes).expand();

    // `source` gives each vertex's nearest trace (CA) atom. The width
    // ramps linearly across an SS boundary (by design - a hard jump would
    // look like a seam), so a cross-section right at the edge of a run is
    // partway between the two widths; only a residue *strictly inside* a
    // run (both neighbors the same code) is guaranteed to be at full
    // width on every cross-section that names it.
    let residue_of_atom = |atom: u32| s.topology.residue_index[atom as usize] as usize;
    let width_of = |base: usize| mesh.positions[base].distance(mesh.positions[base + RING / 2]);
    let is_coil_shaped =
        |c: DsspCode| matches!(c, DsspCode::Coil | DsspCode::Turn | DsspCode::Bend);
    let interior = |same: &dyn Fn(DsspCode) -> bool| {
        (1..codes.len() - 1).find(|&r| same(codes[r]) && same(codes[r - 1]) && same(codes[r + 1]))
    };
    let width_of_residue = |residue: usize| {
        (0..mesh.positions.len() / RING)
            .find(|&cs| residue_of_atom(mesh.source[cs * RING]) == residue)
            .map(|cs| width_of(cs * RING))
    };

    let helix_residue =
        interior(&|c| c == DsspCode::AlphaHelix).expect("1CRN has an interior helix residue");
    let coil_residue = interior(&is_coil_shaped).expect("1CRN has an interior coil-shaped residue");
    let (helix_width, coil_width) = (
        width_of_residue(helix_residue).expect("the helix residue has a cross-section"),
        width_of_residue(coil_residue).expect("the coil residue has a cross-section"),
    );
    assert!(
        (helix_width - RIBBON_WIDTH).abs() < 0.05,
        "helix width {helix_width}"
    );
    assert!(
        helix_width > coil_width * 2.0,
        "helix {helix_width} vs coil {coil_width}"
    );
}

/// Strands must lie flat in their sheet: each interior strand residue's
/// ribbon width should point at its neighbour strand (the nearest CA on
/// another strand, across the H-bonded ladder), not out of the sheet. A
/// curvature-based frame stands strands on edge and fails this.
#[test]
fn strands_lie_flat_in_their_sheet() {
    let s = fixture("1UBQ.cif");
    let coords = s.frame(0);
    let positions = coords.positions();
    let codes = assign(&s.topology, positions);
    let mesh = build(&s.topology, positions, &codes).expand();
    let residue_of = |atom: u32| s.topology.residue_index[atom as usize] as usize;
    let ca_of = |r: usize| {
        s.topology.residues[r]
            .atoms
            .clone()
            .find(|&a| s.topology.atom_name(a as usize) == "CA")
    };
    let strand = |r: usize| codes.get(r) == Some(&DsspCode::Strand);
    let (mut aligned, mut total) = (0, 0);
    for r in 1..codes.len().saturating_sub(1) {
        if !(strand(r - 1) && strand(r) && strand(r + 1)) {
            continue;
        }
        let Some(ca) = ca_of(r) else { continue };
        let here = positions[ca as usize];
        // Nearest CA on a strand residue at least 3 apart in sequence.
        let partner = (0..codes.len())
            .filter(|&q| strand(q) && q.abs_diff(r) >= 3)
            .filter_map(|q| ca_of(q).map(|a| positions[a as usize]))
            .min_by(|a, b| a.distance(here).total_cmp(&b.distance(here)));
        let Some(partner) = partner.filter(|p| p.distance(here) < 6.0) else {
            continue;
        };
        // The ring centred nearest this CA.
        let ring = (0..mesh.positions.len() / RING)
            .filter(|&k| residue_of(mesh.source[k * RING]) == r)
            .min_by(|&a, &b| {
                let c = |k: usize| {
                    mesh.positions[k * RING..(k + 1) * RING]
                        .iter()
                        .copied()
                        .sum::<vv_core::glam::Vec3>()
                        / RING as f32
                };
                c(a).distance(here).total_cmp(&c(b).distance(here))
            });
        let Some(k) = ring else { continue };
        let width_dir =
            (mesh.positions[k * RING] - mesh.positions[k * RING + RING / 2]).normalize();
        let to_partner = (partner - here).normalize();
        total += 1;
        if width_dir.dot(to_partner).abs() > 0.7 {
            aligned += 1;
        }
    }
    assert!(total >= 10, "too few strand residues checked: {total}");
    assert!(
        aligned * 10 >= total * 8,
        "only {aligned} of {total} strand residues have their ribbon in the sheet plane"
    );
}
