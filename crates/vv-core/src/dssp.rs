//! Secondary structure assignment from backbone geometry, following Kabsch
//! & Sander's DSSP algorithm: W. Kabsch and C. Sander, "Dictionary of
//! protein secondary structure: pattern recognition of hydrogen-bonded and
//! geometrical features," Biopolymers 22(12):2577-2637, 1983
//! (doi:10.1002/bip.360221211). Cooperative structure is recognized as
//! repeats of two elementary hydrogen-bond patterns, "turn" and "bridge":
//! repeating turns are helices, repeating bridges are ladders, connected
//! ladders are sheets.
//!
//! This is an independent implementation from the published algorithm, not
//! a port of any existing DSSP codebase; the current, actively maintained
//! DSSP rewrite (`PDB-REDO/dssp`, BSD-2-Clause) was read to confirm
//! implementation details the 1983 paper leaves ambiguous (exact
//! hydrogen-position reconstruction, bridge partner geometry, and the
//! priority DSSP4 gives pi helices over alpha helices on overlap -
//! T.A.H. Touw et al., "A series of PDB-related databanks for everyday
//! needs," Nucleic Acids Research 43:D364-D368, 2015) - never its code.
//!
//! **Phase 1 scope**: helices (H, G, I),
//! bridges and strands (B, E), bends (S) and turns (T). Deliberately not
//! implemented yet: disulfide-bridge numbering, polyproline II, beta
//! bulges (a bridge one residue short of continuing a ladder is not
//! merged into it, so some real sheets will render as more, shorter
//! strands than DSSP proper would show), and DSSP's own sheet/ladder
//! lettering. Operates on one frame; a `*_frames` batch form (matching
//! `crate::analysis`'s convention) is not built yet.

use glam::Vec3;

use crate::analysis::{angle, neighbor_pairs_into};
use crate::residue_class::ResidueClass;
use crate::topology::{SecondaryStructure, Topology};

/// A hydrogen bond is real below this energy (kcal/mol); Kabsch & Sander
/// 1983.
const HBOND_ENERGY_THRESHOLD: f32 = -0.5;
/// `q1 * q2 * f`: partial charges 0.42e (C=O) and 0.20e (N-H), Coulomb
/// constant 332 kcal-A/mol/e^2. Kabsch & Sander 1983, eq. 1-2.
const COUPLING_CONSTANT: f32 = -27.888;
/// Below this atom-atom distance (A) the electrostatic term is a
/// numerical singularity, not real geometry; clamp to the strongest
/// possible bond instead of dividing by (near) zero.
const MIN_ATOM_DISTANCE: f32 = 0.5;
/// Energy clamp: nothing scores stronger than this, healthy geometry or
/// not.
const MIN_HBOND_ENERGY: f32 = -9.9;
/// Longer than this and residue `i+1`'s `N` is not really bonded to
/// residue `i`'s `C` - a chain break (a disordered loop, a new chain),
/// not a strained peptide.
const MAX_PEPTIDE_BOND: f32 = 2.5;
/// CA-CA distance beyond which two residues cannot hydrogen bond; the
/// candidate-pair prefilter that keeps this whole kernel far from O(n^2).
const CA_NEIGHBOR_CUTOFF: f32 = 9.0;
/// Kappa (the virtual bond angle CA(i-2)-CA(i)-CA(i+2)) above this many
/// degrees marks a bend.
const BEND_KAPPA_DEGREES: f32 = 70.0;

/// One residue's DSSP letter. `None` for anything that is not a protein
/// residue with a complete backbone (ligands, water, nucleic acids, and
/// protein residues missing N/CA/C/O, most often the ends of a disordered
/// loop).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
#[repr(u8)]
pub enum DsspCode {
    #[default]
    None = 0,
    /// '-': no assignment.
    Coil = 1,
    /// 'S': high local curvature, not otherwise classified.
    Bend = 2,
    /// 'T': hydrogen-bonded turn, not part of a full helix.
    Turn = 3,
    /// 'B': a single hydrogen-bonded bridge, not extended into a ladder.
    Bridge = 4,
    /// 'E': extended strand, part of a beta ladder of two or more bridges.
    Strand = 5,
    /// 'G': 3-10 helix (i, i+3 hydrogen bonds).
    Helix3_10 = 6,
    /// 'H': alpha helix (i, i+4 hydrogen bonds).
    AlphaHelix = 7,
    /// 'I': pi helix (i, i+5 hydrogen bonds).
    HelixPi = 8,
}

impl DsspCode {
    /// The character DSSP itself prints for this code.
    pub fn letter(self) -> char {
        match self {
            DsspCode::None => '?',
            DsspCode::Coil => '-',
            DsspCode::Bend => 'S',
            DsspCode::Turn => 'T',
            DsspCode::Bridge => 'B',
            DsspCode::Strand => 'E',
            DsspCode::Helix3_10 => 'G',
            DsspCode::AlphaHelix => 'H',
            DsspCode::HelixPi => 'I',
        }
    }

    /// The three-state scheme the rest of vizviz already understands
    /// (`vv_core::select`'s `helix`/`strand`/`coil` classes, the file
    /// parsers). Helix subtypes and bridges collapse; nothing here can be
    /// recovered from `SecondaryStructure` alone, which is why this type
    /// exists rather than widening that one.
    pub fn simplified(self) -> SecondaryStructure {
        match self {
            DsspCode::AlphaHelix | DsspCode::Helix3_10 | DsspCode::HelixPi => {
                SecondaryStructure::Helix
            }
            DsspCode::Strand | DsspCode::Bridge => SecondaryStructure::Strand,
            DsspCode::Bend | DsspCode::Turn | DsspCode::Coil => SecondaryStructure::Coil,
            DsspCode::None => SecondaryStructure::Unknown,
        }
    }
}

#[derive(Clone, Copy, Debug)]
struct Backbone {
    n: u32,
    ca: u32,
    c: u32,
    o: u32,
    chain: u32,
    proline: bool,
}

fn residue_backbone(topology: &Topology, residue: usize) -> Option<Backbone> {
    let rec = &topology.residues[residue];
    if topology.residue_class(residue) != ResidueClass::Protein {
        return None;
    }
    let find = |name: &str| {
        rec.atoms
            .clone()
            .find(|&a| topology.atom_name(a as usize) == name)
    };
    Some(Backbone {
        n: find("N")?,
        ca: find("CA")?,
        c: find("C")?,
        o: find("O")?,
        chain: rec.chain,
        proline: topology.residue_name(residue) == "PRO",
    })
}

/// The backbone atom indices of every residue in `topology`, `None` where
/// a residue is not a complete-backbone protein residue.
fn all_backbones(topology: &Topology) -> Vec<Option<Backbone>> {
    (0..topology.residues.len())
        .map(|r| residue_backbone(topology, r))
        .collect()
}

/// True when residue `i+1` does not continue the same polypeptide as `i`:
/// out of range, either lacks a backbone, on a different chain, or too far
/// for a real peptide bond (`MAX_PEPTIDE_BOND`, above).
fn chain_break(backbones: &[Option<Backbone>], positions: &[Vec3], i: usize) -> bool {
    let (Some(cur), Some(next)) = (
        backbones.get(i).copied().flatten(),
        backbones.get(i + 1).copied().flatten(),
    ) else {
        return true;
    };
    cur.chain != next.chain
        || positions[cur.c as usize].distance(positions[next.n as usize]) > MAX_PEPTIDE_BOND
}

/// The backbone amide hydrogen position DSSP itself uses: `N`, shifted one
/// Angstrom along the previous residue's `C=O` bond direction (trans
/// peptide geometry puts N-H roughly antiparallel to the preceding
/// carbonyl). Falls back to the bare `N` position at a chain start or for
/// proline, which has no amide hydrogen to begin with - `hbond_energy`
/// never uses it as a donor there regardless.
fn h_positions(backbones: &[Option<Backbone>], positions: &[Vec3]) -> Vec<Vec3> {
    backbones
        .iter()
        .enumerate()
        .map(|(i, bb)| {
            let Some(bb) = bb else { return Vec3::ZERO };
            let n = positions[bb.n as usize];
            if bb.proline || i == 0 {
                return n;
            }
            if chain_break(backbones, positions, i - 1) {
                return n;
            }
            let prev = backbones[i - 1].expect("no chain break implies a backbone");
            let (c, o) = (positions[prev.c as usize], positions[prev.o as usize]);
            n + (c - o).normalize_or_zero()
        })
        .collect()
}

/// The Kabsch & Sander electrostatic hydrogen-bond energy (kcal/mol)
/// between a donor's amide (`h`, `n`) and an acceptor's carbonyl (`c`,
/// `o`). Proline residues have no amide hydrogen and never donate;
/// callers skip those without calling this.
fn hbond_energy(h: Vec3, n: Vec3, c: Vec3, o: Vec3) -> f32 {
    let (r_ho, r_hc, r_nc, r_no) = (h.distance(o), h.distance(c), n.distance(c), n.distance(o));
    if [r_ho, r_hc, r_nc, r_no]
        .iter()
        .any(|&r| r < MIN_ATOM_DISTANCE)
    {
        return MIN_HBOND_ENERGY;
    }
    let e = COUPLING_CONSTANT * (1.0 / r_ho - 1.0 / r_hc + 1.0 / r_nc - 1.0 / r_no);
    e.max(MIN_HBOND_ENERGY)
}

/// A donor or acceptor's two strongest hydrogen bonds (DSSP allows each
/// backbone `N-H` and each `C=O` up to two partners at once).
#[derive(Clone, Copy)]
struct BestTwo {
    partner: [Option<u32>; 2],
    energy: [f32; 2],
}

impl Default for BestTwo {
    fn default() -> Self {
        Self {
            partner: [None, None],
            energy: [f32::INFINITY, f32::INFINITY],
        }
    }
}

impl BestTwo {
    fn offer(&mut self, partner: u32, energy: f32) {
        if energy < self.energy[0] {
            self.partner[1] = self.partner[0];
            self.energy[1] = self.energy[0];
            self.partner[0] = Some(partner);
            self.energy[0] = energy;
        } else if energy < self.energy[1] {
            self.partner[1] = Some(partner);
            self.energy[1] = energy;
        }
    }

    fn bonded_to(&self, partner: u32) -> bool {
        (0..2).any(|k| self.partner[k] == Some(partner) && self.energy[k] < HBOND_ENERGY_THRESHOLD)
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum BridgeKind {
    Parallel,
    AntiParallel,
}

/// Assigns DSSP-style secondary structure to every residue of `topology`
/// at `positions` (one frame; see module docs). `positions.len()` must be
/// `topology.atom_count()`.
pub fn assign(topology: &Topology, positions: &[Vec3]) -> Vec<DsspCode> {
    let backbones = all_backbones(topology);
    let n = backbones.len();
    let protein: Vec<usize> = (0..n).filter(|&i| backbones[i].is_some()).collect();
    let mut out = vec![DsspCode::None; n];
    for &i in &protein {
        out[i] = DsspCode::Coil;
    }

    let h = h_positions(&backbones, positions);

    // Candidate pairs via the existing CA neighbor-list kernel, so this
    // stays near-linear instead of testing every residue pair.
    let ca_atoms: Vec<u32> = protein.iter().map(|&i| backbones[i].unwrap().ca).collect();
    let mut contacts = Vec::new();
    neighbor_pairs_into(positions, &ca_atoms, CA_NEIGHBOR_CUTOFF, &mut contacts);

    let mut donors = vec![BestTwo::default(); n]; // this residue's N-H, by acceptor found
    for c in &contacts {
        let i = topology.residue_index[c.a as usize] as usize;
        let j = topology.residue_index[c.b as usize] as usize; // i < j: see below
        let (bi, bj) = (backbones[i].unwrap(), backbones[j].unwrap());
        // i donates to j (i's N-H .. j's C=O).
        if !bi.proline {
            let e = hbond_energy(
                h[i],
                positions[bi.n as usize],
                positions[bj.c as usize],
                positions[bj.o as usize],
            );
            donors[i].offer(j as u32, e);
        }
        // j donates to i, unless j is the very next residue after i: that
        // N-H and this C=O share the covalent peptide bond between them,
        // not a real hydrogen bond.
        let adjacent = j == i + 1 && !chain_break(&backbones, positions, i);
        if !bj.proline && !adjacent {
            let e = hbond_energy(
                h[j],
                positions[bj.n as usize],
                positions[bi.c as usize],
                positions[bi.o as usize],
            );
            donors[j].offer(i as u32, e);
        }
    }
    // `neighbor_pairs_into` reports `a < b` by atom index; CA atom index
    // increases with residue index (residues are chain-major, atom
    // contiguous - topology.rs), so `i < j` always holds above.
    let test_bond = |donor: usize, acceptor: usize| donors[donor].bonded_to(acceptor as u32);

    // A residue's own sequence neighbor, `None` across a chain break.
    let neighbor = |i: usize, delta: isize| -> Option<usize> {
        let j = usize::try_from(i as isize + delta).ok()?;
        if j >= n || backbones[j].is_none() {
            return None;
        }
        if chain_break(&backbones, positions, i.min(j)) {
            return None;
        }
        Some(j)
    };

    // n-turns: residue i+stride's N-H bonds to residue i's C=O. Three
    // strides (3, 4, 5) give 3-10, alpha and pi turns; two consecutive
    // turns of the same stride make a helix of `stride + 1` residues
    // (below).
    let mut turn_start = [vec![false; n], vec![false; n], vec![false; n]]; // stride 3, 4, 5
    for (k, &stride) in [3usize, 4, 5].iter().enumerate() {
        for &i in &protein {
            if i + stride >= n || backbones[i + stride].is_none() {
                continue;
            }
            if (i..i + stride).any(|j| chain_break(&backbones, positions, j)) {
                continue;
            }
            if test_bond(i + stride, i) {
                turn_start[k][i] = true;
            }
        }
    }

    // Alpha first (it wins any overlap it reaches), then 3-10 filling
    // whatever alpha left as coil, then pi - which DSSP4 lets override an
    // alpha assignment, unlike 3-10 (Touw et al. 2015; the 1983 paper
    // predates this nuance).
    for &i in &protein {
        if i == 0 || !turn_start[1][i] || !turn_start[1][i - 1] {
            continue;
        }
        out[i..=i + 3].fill(DsspCode::AlphaHelix);
    }
    for &i in &protein {
        if i == 0 || !turn_start[0][i] || !turn_start[0][i - 1] {
            continue;
        }
        if (i..=i + 2).all(|j| matches!(out[j], DsspCode::Coil | DsspCode::Helix3_10)) {
            out[i..=i + 2].fill(DsspCode::Helix3_10);
        }
    }
    for &i in &protein {
        if i == 0 || !turn_start[2][i] || !turn_start[2][i - 1] {
            continue;
        }
        if (i..=i + 4).all(|j| {
            matches!(
                out[j],
                DsspCode::Coil | DsspCode::HelixPi | DsspCode::AlphaHelix
            )
        }) {
            out[i..=i + 4].fill(DsspCode::HelixPi);
        }
    }

    // Bend: kappa, the virtual bond angle CA(i-2)-CA(i)-CA(i+2), over four
    // unbroken peptide bonds.
    let mut bend = vec![false; n];
    for &i in &protein {
        if i < 2 || i + 2 >= n || backbones[i - 2].is_none() || backbones[i + 2].is_none() {
            continue;
        }
        if (i - 2..i + 2).any(|j| chain_break(&backbones, positions, j)) {
            continue;
        }
        let (ca_prev, ca_cur, ca_next) = (
            backbones[i - 2].unwrap().ca,
            backbones[i].unwrap().ca,
            backbones[i + 2].unwrap().ca,
        );
        let kappa = angle(
            positions,
            ca_prev as usize,
            ca_cur as usize,
            ca_next as usize,
        );
        bend[i] = kappa > BEND_KAPPA_DEGREES;
    }

    // Bridges: the classic Kabsch-Sander hydrogen-bond patterns (1983,
    // fig. 4) between non-adjacent residues `i` (with sequence neighbors
    // `a`, `c`) and `j` (with sequence neighbors `d`, `f`).
    let mut bridge_partner: Vec<Vec<(usize, BridgeKind)>> = vec![Vec::new(); n];
    for c in &contacts {
        let i = topology.residue_index[c.a as usize] as usize;
        let j = topology.residue_index[c.b as usize] as usize;
        if j < i + 3 {
            continue; // this close, it is a turn, not a bridge
        }
        let Some(a) = neighbor(i, -1) else { continue };
        let Some(cc) = neighbor(i, 1) else { continue };
        let Some(d) = neighbor(j, -1) else { continue };
        let Some(f) = neighbor(j, 1) else { continue };
        let kind = if (test_bond(cc, j) && test_bond(j, a)) || (test_bond(f, i) && test_bond(i, d))
        {
            Some(BridgeKind::Parallel)
        } else if (test_bond(cc, d) && test_bond(f, a)) || (test_bond(j, i) && test_bond(i, j)) {
            Some(BridgeKind::AntiParallel)
        } else {
            None
        };
        if let Some(kind) = kind {
            bridge_partner[i].push((j, kind));
            bridge_partner[j].push((i, kind));
        }
    }

    // Any bridge makes a residue at least `Bridge` (B); extending it with
    // a sequence neighbor's bridge to the expected adjacent partner makes
    // the pair (and everything already `Strand` next to them) `Strand`
    // (E) - a ladder. `expected` follows the two diagrams above: a
    // parallel ladder's partner moves the same way as `i`, an
    // antiparallel one moves the opposite way.
    for &i in &protein {
        if !bridge_partner[i].is_empty()
            && !matches!(
                out[i],
                DsspCode::AlphaHelix | DsspCode::Helix3_10 | DsspCode::HelixPi
            )
        {
            out[i] = DsspCode::Bridge;
        }
    }
    loop {
        let mut changed = false;
        for &i in &protein {
            if out[i] != DsspCode::Bridge {
                continue;
            }
            let extends = bridge_partner[i].iter().any(|&(p, kind)| {
                [-1isize, 1].iter().any(|&d| {
                    let Some(ni) = neighbor(i, d) else {
                        return false;
                    };
                    let expected = match kind {
                        BridgeKind::Parallel => p as isize + d,
                        BridgeKind::AntiParallel => p as isize - d,
                    };
                    expected >= 0
                        && bridge_partner[ni]
                            .iter()
                            .any(|&(pj, k2)| k2 == kind && pj as isize == expected)
                        && matches!(out[ni], DsspCode::Bridge | DsspCode::Strand)
                })
            });
            if extends {
                out[i] = DsspCode::Strand;
                changed = true;
            }
        }
        if !changed {
            break;
        }
    }

    // Turn (T): within reach of an n-turn start that did not become a
    // full helix here; Bend (S): high local curvature; otherwise coil.
    for &i in &protein {
        if out[i] != DsspCode::Coil {
            continue;
        }
        let is_turn = [3usize, 4, 5]
            .iter()
            .enumerate()
            .any(|(k, &stride)| (1..stride).any(|back| back <= i && turn_start[k][i - back]));
        out[i] = if is_turn {
            DsspCode::Turn
        } else if bend[i] {
            DsspCode::Bend
        } else {
            DsspCode::Coil
        };
    }

    out
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Real-structure validation (agreement with file SS records, a
    /// textbook helix found by H-bond geometry) lives in
    /// `vv-io/tests/dssp.rs`: `vv-core` cannot depend on `vv-io` (which
    /// depends on it) even for tests, so real coordinates are not
    /// available here. These test the arithmetic in isolation instead.

    #[test]
    fn hbond_energy_matches_the_published_constant_at_ideal_geometry() {
        // A geometrically ideal N-H...O=C hydrogen bond: H...O = 2.0 A,
        // linear (N, H, O, C roughly colinear-ish), the textbook distance
        // an alpha-helix i, i+4 bond sits at. Kabsch & Sander's formula
        // should read out comfortably below the -0.5 kcal/mol threshold.
        let h = Vec3::new(0.0, 0.0, 0.0);
        let n = Vec3::new(0.0, 0.0, -1.0);
        let o = Vec3::new(0.0, 0.0, 2.0);
        let c = Vec3::new(0.0, 1.2, 2.6); // a plausible C=O geometry off O
        let e = hbond_energy(h, n, c, o);
        assert!(e < HBOND_ENERGY_THRESHOLD, "{e}");
        assert!(e >= MIN_HBOND_ENERGY, "{e}");
    }

    #[test]
    fn hbond_energy_is_clamped_for_degenerate_geometry() {
        let p = Vec3::ZERO;
        assert_eq!(hbond_energy(p, p, p, p), MIN_HBOND_ENERGY);
    }

    #[test]
    fn far_apart_atoms_do_not_hydrogen_bond() {
        let h = Vec3::new(0.0, 0.0, 0.0);
        let n = Vec3::new(0.0, 0.0, -1.0);
        let o = Vec3::new(0.0, 0.0, 20.0);
        let c = Vec3::new(0.0, 1.2, 20.6);
        assert!(hbond_energy(h, n, c, o) > HBOND_ENERGY_THRESHOLD);
    }

    #[test]
    fn dssp_code_simplifies_to_the_three_state_scheme() {
        assert_eq!(DsspCode::AlphaHelix.simplified(), SecondaryStructure::Helix);
        assert_eq!(DsspCode::Helix3_10.simplified(), SecondaryStructure::Helix);
        assert_eq!(DsspCode::HelixPi.simplified(), SecondaryStructure::Helix);
        assert_eq!(DsspCode::Strand.simplified(), SecondaryStructure::Strand);
        assert_eq!(DsspCode::Bridge.simplified(), SecondaryStructure::Strand);
        assert_eq!(DsspCode::Turn.simplified(), SecondaryStructure::Coil);
        assert_eq!(DsspCode::Bend.simplified(), SecondaryStructure::Coil);
        assert_eq!(DsspCode::Coil.simplified(), SecondaryStructure::Coil);
        assert_eq!(DsspCode::None.simplified(), SecondaryStructure::Unknown);
        for (code, letter) in [
            (DsspCode::AlphaHelix, 'H'),
            (DsspCode::Helix3_10, 'G'),
            (DsspCode::HelixPi, 'I'),
            (DsspCode::Strand, 'E'),
            (DsspCode::Bridge, 'B'),
            (DsspCode::Turn, 'T'),
            (DsspCode::Bend, 'S'),
            (DsspCode::Coil, '-'),
            (DsspCode::None, '?'),
        ] {
            assert_eq!(code.letter(), letter);
        }
    }
}
