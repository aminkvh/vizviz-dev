//! Strand geometry for a double/triple/aromatic bond: where to place its
//! extra parallel cylinders relative to the one line a single bond draws.
//!
//! A multi-order bond draws as 2 (double, aromatic) or 3 (triple) strands
//! instead of the ordinary single cylinder, offset perpendicular to the
//! bond axis and, when a neighbouring bond gives one, in that bond's
//! plane -- so a ring's strands sit in the ring instead of at an
//! arbitrary angle to it.

use glam::Vec3;

use super::BondTable;

/// CSR adjacency over a [`BondTable`]'s `pairs` (assumed `a < b`, sorted
/// or not): `neighbors[starts[a]..starts[a + 1]]` is every atom bonded to
/// `a`. Built once per structure (O(bonds)), then queried in O(1) per
/// lookup -- the only way to find a multi-order bond's ring neighbour
/// without rescanning every bond for each one.
pub struct Adjacency {
    starts: Vec<u32>,
    neighbors: Vec<u32>,
}

impl Adjacency {
    pub fn build(pairs: &[[u32; 2]], atom_count: usize) -> Self {
        let mut starts = vec![0u32; atom_count + 1];
        for &[a, b] in pairs {
            starts[a as usize + 1] += 1;
            starts[b as usize + 1] += 1;
        }
        for i in 0..atom_count {
            starts[i + 1] += starts[i];
        }
        let mut cursor = starts.clone();
        let mut neighbors = vec![0u32; pairs.len() * 2];
        for &[a, b] in pairs {
            neighbors[cursor[a as usize] as usize] = b;
            cursor[a as usize] += 1;
            neighbors[cursor[b as usize] as usize] = a;
            cursor[b as usize] += 1;
        }
        Self { starts, neighbors }
    }

    fn neighbors_of(&self, atom: u32) -> &[u32] {
        let range = self.starts[atom as usize] as usize..self.starts[atom as usize + 1] as usize;
        &self.neighbors[range]
    }

    /// An atom bonded to `a` or `b`, other than the two of them, to orient
    /// a bond `a -> b`'s strands in its plane (a ring, a carbonyl's
    /// carbon also bonded to N and CA, ...). `None` when neither atom has
    /// another bond (an isolated diatomic).
    pub fn plane_neighbor(&self, a: u32, b: u32) -> Option<u32> {
        self.neighbors_of(a)
            .iter()
            .copied()
            .find(|&n| n != b)
            .or_else(|| self.neighbors_of(b).iter().copied().find(|&n| n != a))
    }
}

/// Shorthand for [`Adjacency::build`] over a [`BondTable`].
pub fn adjacency(bonds: &BondTable, atom_count: usize) -> Adjacency {
    Adjacency::build(&bonds.pairs, atom_count)
}

/// Perpendicular offset direction for a bond `a -> b`: in the plane of
/// `neighbor` (a third atom, typically from [`Adjacency::plane_neighbor`])
/// when given, so a ring's strands stay in the ring; otherwise a stable
/// perpendicular derived from whichever world axis is least parallel to
/// the bond, so it never degenerates. Always a unit vector (or `Vec3::
/// ZERO` if `a == b`, which a real bond never is).
pub fn strand_axis(a: Vec3, b: Vec3, neighbor: Option<Vec3>) -> Vec3 {
    let axis = (b - a).normalize_or_zero();
    if axis == Vec3::ZERO {
        return Vec3::ZERO;
    }
    let seed = neighbor.map(|n| n - a);
    let side = seed.map(|s| s - axis * s.dot(axis));
    side.and_then(|s| s.try_normalize())
        .unwrap_or_else(|| stable_perpendicular(axis))
}

/// A perpendicular to `axis` with no external reference to align to: the
/// world axis least parallel to `axis` can never be degenerate with it.
fn stable_perpendicular(axis: Vec3) -> Vec3 {
    let seed = if axis.x.abs() < 0.9 { Vec3::X } else { Vec3::Y };
    (seed - axis * seed.dot(axis)).normalize()
}

/// One parallel strand of a multi-order bond. `offset` is a signed
/// distance along [`strand_axis`]; `0.0` (a triple bond's centre strand)
/// draws on the true bond axis with no plane lookup needed.
#[derive(Clone, Copy, Debug)]
pub struct BondStrand {
    /// Index into the source [`BondTable::pairs`], for picking and colour.
    pub bond: u32,
    pub atoms: [u32; 2],
    pub plane_neighbor: Option<u32>,
    pub offset: f32,
}

/// Every strand `bonds`'s non-`Single` entries need, `separation` apart
/// (double: `±separation / 2`; triple: `0, ±separation`). Replaces, not
/// adds to, those bonds' ordinary single cylinder -- the caller excludes
/// any bond with `BondTable::order_of != Single` from its normal bond
/// list and draws these instead.
pub fn bond_strands(bonds: &BondTable, adjacency: &Adjacency, separation: f32) -> Vec<BondStrand> {
    let mut out = Vec::with_capacity(bonds.orders.len() * 2);
    for &(bond, order) in &bonds.orders {
        let [a, b] = bonds.pairs[bond as usize];
        let plane_neighbor = adjacency.plane_neighbor(a, b);
        let offsets: &[f32] = match order.strand_count() {
            2 => &[separation * 0.5, -separation * 0.5],
            3 => &[separation, 0.0, -separation],
            _ => &[0.0],
        };
        out.extend(offsets.iter().map(|&offset| BondStrand {
            bond,
            atoms: [a, b],
            plane_neighbor,
            offset,
        }));
    }
    out
}

/// World-space endpoints of `strand` at the given (possibly just-moved)
/// `positions`: recomputed from live atom positions rather than cached,
/// so a played trajectory's strands track their real bond and ring plane
/// every frame instead of drifting from a frame-0 offset.
pub fn bond_strand_endpoints(strand: &BondStrand, positions: &[Vec3]) -> (Vec3, Vec3) {
    let [a, b] = strand.atoms;
    let (pa, pb) = (positions[a as usize], positions[b as usize]);
    if strand.offset == 0.0 {
        return (pa, pb);
    }
    let neighbor_pos = strand.plane_neighbor.map(|n| positions[n as usize]);
    let axis = strand_axis(pa, pb, neighbor_pos) * strand.offset;
    (pa + axis, pb + axis)
}

#[cfg(test)]
mod tests {
    use super::super::BondOrder;
    use super::*;

    #[test]
    fn double_bond_strands_are_parallel_separated_and_symmetric() {
        let bonds = BondTable {
            pairs: vec![[0, 1]],
            orders: vec![(0, BondOrder::Double)],
        };
        let adj = Adjacency::build(&bonds.pairs, 2);
        let strands = bond_strands(&bonds, &adj, 0.3);
        assert_eq!(strands.len(), 2);
        let positions = [Vec3::ZERO, Vec3::new(1.5, 0.0, 0.0)];
        let (a0, b0) = bond_strand_endpoints(&strands[0], &positions);
        let (a1, b1) = bond_strand_endpoints(&strands[1], &positions);

        let bond_axis = (positions[1] - positions[0]).normalize();
        // Parallel: each strand runs along the same direction as the bond.
        assert!((a0 - b0).normalize().dot(bond_axis).abs() > 0.9999);
        assert!((a1 - b1).normalize().dot(bond_axis).abs() > 0.9999);
        // Separated by exactly `separation`.
        assert!((a0.distance(a1) - 0.3).abs() < 1e-5);
        assert!((b0.distance(b1) - 0.3).abs() < 1e-5);
        // Symmetric about the true bond axis: the two offsets cancel.
        let mid = (a0 - positions[0]) + (a1 - positions[0]);
        assert!(mid.length() < 1e-5);
        // Perpendicular to the bond, not along it.
        assert!((a0 - positions[0]).dot(bond_axis).abs() < 1e-5);
    }

    #[test]
    fn triple_bond_has_one_centre_and_two_offset_strands() {
        let bonds = BondTable {
            pairs: vec![[0, 1]],
            orders: vec![(0, BondOrder::Triple)],
        };
        let adj = Adjacency::build(&bonds.pairs, 2);
        let strands = bond_strands(&bonds, &adj, 0.3);
        assert_eq!(strands.len(), 3);
        let positions = [Vec3::ZERO, Vec3::new(1.2, 0.0, 0.0)];
        let centre = strands.iter().find(|s| s.offset == 0.0).unwrap();
        let (a, b) = bond_strand_endpoints(centre, &positions);
        assert_eq!((a, b), (positions[0], positions[1]));
    }

    #[test]
    fn strand_axis_prefers_the_neighbour_plane() {
        let a = Vec3::ZERO;
        let b = Vec3::new(1.0, 0.0, 0.0);
        // A neighbour in the XY plane: the strand offset must stay in it
        // (no Z component), unlike the arbitrary fallback perpendicular.
        let neighbor = Vec3::new(-0.5, 1.0, 0.0);
        let axis = strand_axis(a, b, Some(neighbor));
        assert!(axis.z.abs() < 1e-5);
        assert!(axis.dot((b - a).normalize()).abs() < 1e-5);
    }

    #[test]
    fn strand_axis_falls_back_when_no_neighbour_or_colinear() {
        let a = Vec3::ZERO;
        let b = Vec3::new(1.0, 0.0, 0.0);
        assert!(strand_axis(a, b, None).length() > 0.99);
        // Neighbour exactly on the bond line: no plane to derive.
        let colinear = Vec3::new(2.0, 0.0, 0.0);
        assert!(strand_axis(a, b, Some(colinear)).length() > 0.99);
    }

    #[test]
    fn adjacency_finds_a_ring_neighbor_on_either_side() {
        // A 3-atom ring: 0-1, 1-2, 2-0.
        let pairs = vec![[0, 1], [1, 2], [0, 2]];
        let adj = Adjacency::build(&pairs, 3);
        assert_eq!(adj.plane_neighbor(0, 1), Some(2));
        assert!(adj.plane_neighbor(1, 2).is_some());
    }
}
