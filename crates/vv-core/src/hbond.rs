//! General hydrogen-bond detection from geometry alone, for any polar
//! pair - sidechain, ligand, water - not just the protein backbone:
//! E.N. Baker and R.E. Hubbard, "Hydrogen bonding in globular proteins,"
//! Progress in Biophysics and Molecular Biology 44(2):97-179, 1984
//! (doi:10.1016/0079-6107(84)90007-5). Their criterion for structures
//! with no observed hydrogen positions (the ordinary case for
//! crystallographic heavy-atom coordinates): a donor...acceptor distance
//! of at most 3.5 A, and the angle at each atom's own covalent neighbor -
//! antecedent-atom, the atom itself, its candidate partner - of at least
//! 90 degrees.
//!
//! This is deliberately separate from `vv_core::dssp`'s hydrogen-bond
//! step, not a generalization of it: DSSP's Kabsch & Sander 1983
//! electrostatic model is specific to the backbone amide/carbonyl dipole
//! pair and tuned for secondary-structure patterns; this kernel is
//! general-purpose and purely geometric. The two models will report
//! overlapping but not identical answers for the same backbone pair -
//! expected, since they are independently valid models measuring
//! different things, not the same computation twice.
//!
//! No hydrogen positions are assumed, so donor and acceptor roles are
//! not distinguished - Baker & Hubbard's no-hydrogen criterion is
//! symmetric in exactly this way, constraining the geometry around each
//! heavy atom's own covalent neighbor rather than a hydrogen nobody has
//! coordinates for. An atom with no covalent bond at all (a bare ion) is
//! never a candidate: a hydrogen bond needs a covalent antecedent to
//! measure the angle from, and metal coordination is a different
//! phenomenon with its own geometry.

use glam::Vec3;

use crate::analysis::{angle, neighbor_pairs_into};
use crate::bonds::BondTable;
use crate::element::Element;

/// Donor...acceptor distance beyond which no hydrogen bond is possible.
/// Baker & Hubbard 1984 (no observed hydrogens).
pub const MAX_DISTANCE: f32 = 3.5;
/// The angle at each atom's own covalent neighbor, through the atom,
/// toward its candidate partner, must be at least this. Baker & Hubbard
/// 1984.
pub const MIN_ANTECEDENT_ANGLE: f32 = 90.0;

/// One candidate hydrogen bond: `a` and `b` are polar atoms within range
/// with plausible antecedent geometry on both sides. Roles (donor/
/// acceptor) are not assigned - see module docs.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HydrogenBond {
    pub a: u32,
    pub b: u32,
    pub distance: f32,
}

/// Nitrogen, oxygen, or sulfur: the elements that plausibly donate or
/// accept a hydrogen bond in a biomolecule (sulfur mainly cysteine
/// thiols; rarer and weaker than N/O but a real, literature-recognized
/// case, not a stretch to include).
fn is_candidate(element: Element) -> bool {
    matches!(
        element,
        Element::NITROGEN | Element::OXYGEN | Element::SULFUR
    )
}

/// Every atom bonded to `atom`, from `bonds`. A linear scan of the whole
/// table per atom would be O(atoms x bonds); `by_atom` (built once by
/// [`hydrogen_bonds_into`]) makes this O(degree).
fn neighbors_of(by_atom: &[Vec<u32>], atom: u32) -> &[u32] {
    by_atom.get(atom as usize).map_or(&[], Vec::as_slice)
}

fn antecedent_angle_ok(by_atom: &[Vec<u32>], positions: &[Vec3], atom: u32, partner: u32) -> bool {
    let neighbors = neighbors_of(by_atom, atom);
    if neighbors.is_empty() {
        return false; // no covalent antecedent: not a real candidate
    }
    neighbors.iter().any(|&nb| {
        angle(positions, nb as usize, atom as usize, partner as usize) >= MIN_ANTECEDENT_ANGLE
    })
}

/// Every N/O/S...N/O/S pair meeting Baker & Hubbard's 1984 no-hydrogen
/// geometric criteria: distance at most [`MAX_DISTANCE`], and the angle
/// from each atom's own bonded neighbor through that atom to its partner
/// at least [`MIN_ANTECEDENT_ANGLE`] on *both* sides (module docs explain
/// why both, rather than picking one atom as "the" donor). `out` is
/// cleared first. `bonds` should already reflect `positions` (from
/// `vv_core::bonds::perceive` or a file's explicit bonds merged in);
/// stale bonds only misjudge which atoms are candidates, they cannot
/// crash this.
pub fn hydrogen_bonds_into(
    element: &[Element],
    bonds: &BondTable,
    positions: &[Vec3],
    out: &mut Vec<HydrogenBond>,
) {
    out.clear();
    assert_eq!(element.len(), positions.len());
    let candidates: Vec<u32> = (0..element.len() as u32)
        .filter(|&a| is_candidate(element[a as usize]))
        .collect();
    if candidates.is_empty() {
        return;
    }

    let mut by_atom = vec![Vec::new(); element.len()];
    for &[a, b] in &bonds.pairs {
        by_atom[a as usize].push(b);
        by_atom[b as usize].push(a);
    }

    let mut contacts = Vec::new();
    neighbor_pairs_into(positions, &candidates, MAX_DISTANCE, &mut contacts);
    out.extend(contacts.into_iter().filter_map(|c| {
        let ok = antecedent_angle_ok(&by_atom, positions, c.a, c.b)
            && antecedent_angle_ok(&by_atom, positions, c.b, c.a);
        ok.then_some(HydrogenBond {
            a: c.a,
            b: c.b,
            distance: c.distance,
        })
    }));
}

#[cfg(test)]
mod tests {
    use super::*;

    /// N=O with N's antecedent H1 straight across from O, and O's own
    /// antecedent C straight across from N: as textbook an H-bond
    /// geometry as a hand-built test gets.
    #[test]
    fn a_head_on_pair_with_good_antecedent_geometry_bonds() {
        let element = vec![
            Element::CARBON,
            Element::NITROGEN,
            Element::OXYGEN,
            Element::CARBON,
        ];
        // 0 (C) -- 1 (N) ... 2 (O) -- 3 (C), all colinear, N...O = 3.0 A.
        let positions = vec![
            Vec3::new(-1.5, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(3.0, 0.0, 0.0),
            Vec3::new(4.5, 0.0, 0.0),
        ];
        let bonds = BondTable {
            pairs: vec![[0, 1], [2, 3]],
            orders: Vec::new(),
        };
        let mut out = Vec::new();
        hydrogen_bonds_into(&element, &bonds, &positions, &mut out);
        assert_eq!(out.len(), 1);
        assert_eq!((out[0].a, out[0].b), (1, 2));
        assert!((out[0].distance - 3.0).abs() < 1e-4);
    }

    #[test]
    fn too_far_apart_does_not_bond() {
        let element = vec![
            Element::NITROGEN,
            Element::CARBON,
            Element::OXYGEN,
            Element::CARBON,
        ];
        let positions = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(-1.5, 0.0, 0.0),
            Vec3::new(4.0, 0.0, 0.0), // 4.0 A > MAX_DISTANCE
            Vec3::new(5.5, 0.0, 0.0),
        ];
        let bonds = BondTable {
            pairs: vec![[0, 1], [2, 3]],
            orders: Vec::new(),
        };
        let mut out = Vec::new();
        hydrogen_bonds_into(&element, &bonds, &positions, &mut out);
        assert!(out.is_empty());
    }

    #[test]
    fn a_grazing_angle_does_not_bond() {
        // Same N...O distance as the good-geometry case, but the
        // acceptor's antecedent puts the O-C bond nearly pointing back at
        // N instead of across from it: angle(C, O, N) is about 17 deg,
        // well under the 90 deg minimum.
        let element = vec![
            Element::CARBON,
            Element::NITROGEN,
            Element::OXYGEN,
            Element::CARBON,
        ];
        let positions = vec![
            Vec3::new(-1.5, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(3.0, 0.0, 0.0),
            Vec3::new(2.0, 0.3, 0.0),
        ];
        let bonds = BondTable {
            pairs: vec![[0, 1], [2, 3]],
            orders: Vec::new(),
        };
        let mut out = Vec::new();
        hydrogen_bonds_into(&element, &bonds, &positions, &mut out);
        assert!(out.is_empty(), "{out:?}");
    }

    #[test]
    fn an_atom_with_no_covalent_bond_is_never_a_candidate() {
        // An oxygen close enough to a carbonyl but with no perceived bond
        // of its own (an isolated ion, or a water oxygen bond perception
        // missed): a real N/O element, but nothing bonded to it to
        // measure an antecedent angle from.
        let element = vec![Element::OXYGEN, Element::CARBON, Element::OXYGEN];
        let positions = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(-1.2, 0.0, 0.0),
            Vec3::new(2.5, 0.0, 0.0),
        ];
        let bonds = BondTable {
            pairs: vec![[0, 1]],
            orders: Vec::new(),
        };
        let mut out = Vec::new();
        hydrogen_bonds_into(&element, &bonds, &positions, &mut out);
        assert!(out.is_empty(), "{out:?}");
    }

    #[test]
    fn real_geometry_perceived_bonds_and_elements_from_a_built_structure() {
        use crate::builder::{AtomRow, TopologyBuilder};

        // The same head-on geometry as the first test, but through the
        // real pipeline: `bonds::perceive` (geometric, from covalent
        // radii) instead of a hand-fed `BondTable`, and elements read
        // back off the built `Topology` instead of a literal slice. Real
        // structures (`vv-io/tests/`) are where a genuine fold belongs;
        // this only checks the plumbing between `TopologyBuilder`,
        // `bonds::perceive`, and this kernel.
        let atoms = [
            ("C1", Element::CARBON, [-1.5, 0.0, 0.0]),
            ("N", Element::NITROGEN, [0.0, 0.0, 0.0]),
            ("O", Element::OXYGEN, [3.0, 0.0, 0.0]),
            ("C2", Element::CARBON, [4.5, 0.0, 0.0]),
        ];
        let mut builder = TopologyBuilder::new();
        for (i, (name, element, pos)) in atoms.iter().enumerate() {
            let mut n = [b' '; 4];
            n[..name.len()].copy_from_slice(name.as_bytes());
            builder.push(&AtomRow {
                element: *element,
                name: n,
                serial: i as u32 + 1,
                alt_loc: 0,
                comp: "ALA",
                asym: "A",
                auth_asym: "A",
                seq_id: 1,
                auth_seq_id: 1,
                ins_code: 0,
                entity: 1,
                position: Vec3::from(*pos),
                occupancy: 1.0,
                b_factor: 0.0,
                charge: 0,
                hetero: false,
            });
        }
        let structure = builder.finish().unwrap();
        let bonds = crate::bonds::perceive(&structure.topology, structure.frame(0).positions());
        assert!(bonds.contains(0, 1), "C1-N should have perceived a bond");
        assert!(bonds.contains(2, 3), "O-C2 should have perceived a bond");
        let mut out = Vec::new();
        hydrogen_bonds_into(
            &structure.topology.element,
            &bonds,
            structure.frame(0).positions(),
            &mut out,
        );
        assert!(
            out.iter().any(|hb| (hb.a, hb.b) == (1, 2)),
            "expected N...O (indices 1, 2) among {out:?}"
        );
    }
}
