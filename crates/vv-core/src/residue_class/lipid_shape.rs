//! Structural lipid test for a residue whose name no table knows.
//!
//! A residue is lipid-like when it has a long unbranched alkyl run **and**
//! a polar head made of an acyl group or a phosphate/sulfate; or when it
//! is a sterol (a carbon/oxygen skeleton of at least four fused rings and
//! 27 carbons, cholesterol's size). MD topologies carry no bond orders, so
//! "sp3" is judged from valence when the residue has hydrogens (a carbon
//! with four neighbours) and not at all when it has none.
//!
//! Limits: a drug with a decyl chain and an ester passes; a saturated
//! chain whose only polar group is an ether or alcohol (alkyl glycosides,
//! polyethylene glycol surfactants) does not; heavy-atom-only structures
//! cannot exclude a conjugated polyene chain.

use std::ops::Range;

use crate::Element;

/// Consecutive unbranched carbons that make a tail.
const CHAIN_MIN: usize = 8;
const STEROL_RINGS_MIN: usize = 4;
const STEROL_CARBONS_MIN: usize = 27;
const STEROL_OXYGENS_MAX: usize = 2;

/// One residue's bond graph, read out of a shared CSR adjacency.
pub(super) struct ResidueGraph<'a> {
    pub element: &'a [Element],
    pub offsets: &'a [u32],
    pub neighbors: &'a [u32],
    pub atoms: Range<u32>,
}

impl ResidueGraph<'_> {
    pub fn looks_like_lipid(&self) -> bool {
        (self.has_long_chain() && self.has_polar_head()) || self.is_sterol()
    }

    fn neighbors(&self, atom: u32) -> impl Iterator<Item = u32> + '_ {
        let (lo, hi) = (
            self.offsets[atom as usize] as usize,
            self.offsets[atom as usize + 1] as usize,
        );
        self.neighbors[lo..hi]
            .iter()
            .copied()
            .filter(|n| self.atoms.contains(n))
    }

    fn z(&self, atom: u32) -> u8 {
        self.element[atom as usize].atomic_number()
    }

    fn degree(&self, atom: u32) -> usize {
        self.neighbors(atom).count()
    }

    fn count_neighbors(&self, atom: u32, z: u8) -> usize {
        self.neighbors(atom).filter(|&n| self.z(n) == z).count()
    }

    fn has_hydrogen(&self) -> bool {
        self.atoms
            .clone()
            .any(|a| self.element[a as usize].is_hydrogen())
    }

    /// A carbon that can sit inside a saturated tail: bonded only to C and
    /// H, at most two of them carbon, and (with hydrogens present) sp3.
    fn is_chain_carbon(&self, atom: u32, check_valence: bool) -> bool {
        self.z(atom) == 6
            && self.neighbors(atom).all(|n| matches!(self.z(n), 1 | 6))
            && self.count_neighbors(atom, 6) <= 2
            && (!check_valence || self.degree(atom) == 4)
    }

    fn has_long_chain(&self) -> bool {
        let check_valence = self.has_hydrogen();
        let start = self.atoms.start;
        let mut seen = vec![false; self.atoms.len()];
        for seed in self.atoms.clone() {
            if seen[(seed - start) as usize] || !self.is_chain_carbon(seed, check_valence) {
                continue;
            }
            let (nodes, edges) = self.chain_component(seed, check_valence, &mut seen);
            if nodes >= CHAIN_MIN && edges + 1 == nodes {
                return true;
            }
        }
        false
    }

    /// Node and edge counts of the chain-carbon component around `seed`.
    fn chain_component(&self, seed: u32, check_valence: bool, seen: &mut [bool]) -> (usize, usize) {
        let start = self.atoms.start;
        let (mut nodes, mut edge_ends) = (0, 0);
        let mut stack = vec![seed];
        seen[(seed - start) as usize] = true;
        while let Some(a) = stack.pop() {
            nodes += 1;
            for n in self
                .neighbors(a)
                .filter(|&n| self.is_chain_carbon(n, check_valence))
            {
                edge_ends += 1;
                if !seen[(n - start) as usize] {
                    seen[(n - start) as usize] = true;
                    stack.push(n);
                }
            }
        }
        (nodes, edge_ends / 2)
    }

    /// Phosphate or sulfate (three or more O), or an acyl carbon: a carbon
    /// with a carbon neighbour, a terminal O and a second O or N (ester,
    /// acid, amide).
    fn has_polar_head(&self) -> bool {
        self.atoms.clone().any(|a| match self.z(a) {
            15 | 16 => self.count_neighbors(a, 8) >= 3,
            6 => self.is_acyl(a),
            _ => false,
        })
    }

    fn is_acyl(&self, atom: u32) -> bool {
        let terminal_oxygen = self
            .neighbors(atom)
            .any(|n| self.z(n) == 8 && self.degree(n) == 1);
        let hetero = self.count_neighbors(atom, 8) + self.count_neighbors(atom, 7);
        terminal_oxygen && hetero >= 2 && self.count_neighbors(atom, 6) >= 1
    }

    fn is_sterol(&self) -> bool {
        let heavy: Vec<u32> = self
            .atoms
            .clone()
            .filter(|&a| !self.element[a as usize].is_hydrogen())
            .collect();
        let count = |z| heavy.iter().filter(|&&a| self.z(a) == z).count();
        let (carbons, oxygens) = (count(6), count(8));
        if carbons < STEROL_CARBONS_MIN
            || oxygens > STEROL_OXYGENS_MAX
            || carbons + oxygens != heavy.len()
        {
            return false;
        }
        let edges: usize = heavy
            .iter()
            .map(|&a| {
                self.neighbors(a)
                    .filter(|&n| n > a && !self.element[n as usize].is_hydrogen())
                    .count()
            })
            .sum();
        edges + 1 >= heavy.len() + STEROL_RINGS_MIN
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Runs `check` on the graph with these atoms and undirected edges.
    fn with_graph<R>(
        elements: &[Element],
        edges: &[(u32, u32)],
        check: impl FnOnce(&ResidueGraph) -> R,
    ) -> R {
        let n = elements.len();
        let mut lists = vec![Vec::new(); n];
        for &(a, b) in edges {
            lists[a as usize].push(b);
            lists[b as usize].push(a);
        }
        let mut offsets = vec![0u32];
        let mut neighbors = Vec::new();
        for l in &lists {
            neighbors.extend(l);
            offsets.push(neighbors.len() as u32);
        }
        check(&ResidueGraph {
            element: elements,
            offsets: &offsets,
            neighbors: &neighbors,
            atoms: 0..n as u32,
        })
    }

    /// 27 carbons and one hydroxyl O in a chain, plus ring-closing bonds.
    fn sterol_like(closures: &[(u32, u32)]) -> bool {
        let mut elements = vec![Element::CARBON; 27];
        elements.push(Element::OXYGEN);
        let mut edges: Vec<(u32, u32)> = (0..27).map(|i| (i, i + 1)).collect();
        edges.extend_from_slice(closures);
        with_graph(&elements, &edges, |g| g.is_sterol())
    }

    #[test]
    fn four_rings_and_27_carbons_make_a_sterol() {
        assert!(sterol_like(&[(0, 5), (6, 11), (12, 17), (18, 23)]));
        assert!(!sterol_like(&[(0, 5), (6, 11), (12, 17)]));
    }

    #[test]
    fn a_ring_of_carbons_is_not_a_chain() {
        let elements = vec![Element::CARBON; 10];
        let mut edges: Vec<(u32, u32)> = (0..9).map(|i| (i, i + 1)).collect();
        assert!(with_graph(&elements, &edges, |g| g.has_long_chain()));
        edges.push((9, 0));
        assert!(!with_graph(&elements, &edges, |g| g.has_long_chain()));
    }

    #[test]
    fn a_branched_carbon_ends_the_run() {
        let elements = vec![Element::CARBON; 14];
        let mut edges: Vec<(u32, u32)> = (0..9).map(|i| (i, i + 1)).collect();
        edges.extend([(4, 10), (4, 11)]);
        assert!(!with_graph(&elements, &edges, |g| g.has_long_chain()));
    }
}
