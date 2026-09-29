//! The structural layer: classes for residues whose names no table knows.
//!
//! Only those residues, plus one chain neighbour on each side, are copied
//! into a small sub-topology whose bonds are perceived, so a system whose
//! names are all known costs nothing here and one with a few unknown
//! ligands costs a few small perceptions.

use std::ops::Range;

use glam::Vec3;
use rayon::prelude::*;

use super::lipid_shape::ResidueGraph;
use super::ResidueClass;
use crate::bonds::{self, BondTable};
use crate::{ChainRec, Element, ResidueRec, Topology};

/// Atoms and residues of a lone O with hydrogens, at most.
const WATER_ATOMS_MAX: usize = 3;

pub(super) fn resolve(
    t: &Topology,
    positions: Option<&[Vec3]>,
    mut classes: Vec<Option<ResidueClass>>,
) -> Vec<ResidueClass> {
    classes
        .par_iter_mut()
        .zip(t.residues.par_iter())
        .for_each(|(class, rec)| {
            if class.is_none() {
                *class = by_elements(t, rec);
            }
        });
    let unknown: Vec<u32> = classes
        .par_iter()
        .enumerate()
        .filter_map(|(i, c)| c.is_none().then_some(i as u32))
        .collect();
    if !unknown.is_empty() {
        for (residue, class) in resolve_by_bonds(t, positions, &unknown, &classes) {
            classes[residue as usize] = Some(class);
        }
    }
    classes.into_iter().map(Option::unwrap_or_default).collect()
}

/// Rules that need only a residue's own elements.
pub(super) fn by_elements(t: &Topology, rec: &ResidueRec) -> Option<ResidueClass> {
    let elements = &t.element[rec.atoms.start as usize..rec.atoms.end as usize];
    if elements.is_empty() || elements.iter().all(|e| e.is_unknown()) {
        return Some(ResidueClass::Other);
    }
    let heavy: Vec<Element> = elements
        .iter()
        .copied()
        .filter(|e| !e.is_hydrogen())
        .collect();
    match heavy[..] {
        [e] if elements.len() == 1 && is_ion_element(e) => Some(ResidueClass::Ion),
        [e] if e.atomic_number() == 8 && elements.len() <= WATER_ATOMS_MAX => {
            Some(ResidueClass::Water)
        }
        _ => None,
    }
}

/// Metals and the halides: what a lone atom of an unrecognized residue
/// name is taken to be an ion of.
fn is_ion_element(e: Element) -> bool {
    e.is_metal() || matches!(e.atomic_number(), 9 | 17 | 35 | 53)
}

/// Unknown residues (`unknown`, ascending) plus their chain neighbours,
/// as maximal residue ranges within a chain.
fn windows(t: &Topology, unknown: &[u32]) -> Vec<Range<u32>> {
    let mut out: Vec<Range<u32>> = Vec::new();
    for &u in unknown {
        let chain = &t.chains[t.residues[u as usize].chain as usize];
        let lo = u.saturating_sub(1).max(chain.residues.start);
        let hi = (u + 2).min(chain.residues.end);
        match out.last_mut() {
            Some(w) if lo <= w.end && same_chain(t, lo, w.end - 1) => w.end = w.end.max(hi),
            _ => out.push(lo..hi),
        }
    }
    out
}

fn same_chain(t: &Topology, a: u32, b: u32) -> bool {
    t.residues[a as usize].chain == t.residues[b as usize].chain
}

/// The windows' residues copied into one topology (one chain per window)
/// with their bonds and a CSR adjacency over the copied atoms.
struct SubSystem {
    topology: Topology,
    bonds: BondTable,
    offsets: Vec<u32>,
    neighbors: Vec<u32>,
    /// First sub-residue of each window.
    first_residue: Vec<u32>,
}

impl SubSystem {
    fn build(t: &Topology, positions: Option<&[Vec3]>, windows: &[Range<u32>]) -> Option<Self> {
        if positions.is_none() && t.md_bonds.is_none() {
            return None;
        }
        let mut sub = Topology {
            names: t.names.clone(),
            ..Topology::default()
        };
        let mut atom_map: Vec<u32> = Vec::new();
        let mut first_residue = Vec::with_capacity(windows.len());
        for w in windows {
            first_residue.push(sub.residues.len() as u32);
            copy_window(t, w, &mut sub, &mut atom_map);
        }
        let sub_positions: Vec<Vec3> = match positions {
            Some(p) => atom_map.iter().map(|&a| p[a as usize]).collect(),
            None => vec![Vec3::ZERO; atom_map.len()],
        };
        sub.md_bonds = t.md_bonds.as_ref().map(|b| remap_bonds(b, &atom_map));
        let bonds = bonds::perceive(&sub, &sub_positions);
        let (offsets, neighbors) = csr(&bonds, atom_map.len());
        Some(Self {
            topology: sub,
            bonds,
            offsets,
            neighbors,
            first_residue,
        })
    }

    fn residue_atoms(&self, residue: u32) -> Range<u32> {
        self.topology.residues[residue as usize].atoms.clone()
    }

    fn atom_named(&self, residue: u32, name: &str) -> Option<u32> {
        self.residue_atoms(residue)
            .find(|&a| self.topology.atom_name(a as usize) == name)
    }

    /// Whether `from`'s atom `from_name` is bonded to `to`'s atom `to_name`.
    fn linked(&self, from: u32, from_name: &str, to: u32, to_name: &str) -> bool {
        match (
            self.atom_named(from, from_name),
            self.atom_named(to, to_name),
        ) {
            (Some(a), Some(b)) => self.bonds.contains(a, b),
            _ => false,
        }
    }

    fn graph(&self, residue: u32) -> ResidueGraph<'_> {
        ResidueGraph {
            element: &self.topology.element,
            offsets: &self.offsets,
            neighbors: &self.neighbors,
            atoms: self.residue_atoms(residue),
        }
    }
}

fn copy_window(t: &Topology, w: &Range<u32>, sub: &mut Topology, atom_map: &mut Vec<u32>) {
    let chain = &t.chains[t.residues[w.start as usize].chain as usize];
    let first = sub.residues.len() as u32;
    for r in w.clone() {
        let rec = &t.residues[r as usize];
        let start = atom_map.len() as u32;
        for a in rec.atoms.clone() {
            atom_map.push(a);
            sub.element.push(t.element[a as usize]);
            if !t.name.is_empty() {
                sub.name.push(t.name[a as usize]);
            }
            if !t.alt_loc.is_empty() {
                sub.alt_loc.push(t.alt_loc[a as usize]);
            }
            sub.residue_index.push(sub.residues.len() as u32);
        }
        sub.residues.push(ResidueRec {
            atoms: start..atom_map.len() as u32,
            chain: sub.chains.len() as u32,
            ..rec.clone()
        });
    }
    sub.chains.push(ChainRec {
        residues: first..sub.residues.len() as u32,
        ..chain.clone()
    });
}

/// The MD bonds with both ends inside the sub-system, in sub indices.
fn remap_bonds(bonds: &[[u32; 2]], atom_map: &[u32]) -> Vec<[u32; 2]> {
    let find = |a: u32| atom_map.binary_search(&a).ok().map(|i| i as u32);
    bonds
        .iter()
        .filter_map(|&[a, b]| Some([find(a)?, find(b)?]))
        .collect()
}

fn csr(bonds: &BondTable, atoms: usize) -> (Vec<u32>, Vec<u32>) {
    let degree = bonds.degrees(atoms);
    let mut offsets = vec![0u32; atoms + 1];
    for (i, d) in degree.iter().enumerate() {
        offsets[i + 1] = offsets[i] + d;
    }
    let mut cursor = offsets.clone();
    let mut neighbors = vec![0u32; offsets[atoms] as usize];
    for &[a, b] in &bonds.pairs {
        neighbors[cursor[a as usize] as usize] = b;
        neighbors[cursor[b as usize] as usize] = a;
        cursor[a as usize] += 1;
        cursor[b as usize] += 1;
    }
    (offsets, neighbors)
}

fn resolve_by_bonds(
    t: &Topology,
    positions: Option<&[Vec3]>,
    unknown: &[u32],
    classes: &[Option<ResidueClass>],
) -> Vec<(u32, ResidueClass)> {
    let windows = windows(t, unknown);
    let Some(sub) = SubSystem::build(t, positions, &windows) else {
        return Vec::new();
    };
    windows
        .par_iter()
        .enumerate()
        .flat_map_iter(|(k, w)| resolve_window(&sub, sub.first_residue[k], w, classes))
        .collect()
}

/// Classes for the unknown residues of one window. Residues chained
/// through unknown ones (a modified stretch of a chain) resolve by
/// repeating the link test until nothing changes.
fn resolve_window(
    sub: &SubSystem,
    first: u32,
    w: &Range<u32>,
    classes: &[Option<ResidueClass>],
) -> Vec<(u32, ResidueClass)> {
    let mut local: Vec<Option<ResidueClass>> = classes[w.start as usize..w.end as usize].to_vec();
    let mut changed = true;
    while changed {
        changed = false;
        for i in 0..local.len() {
            if local[i].is_none() {
                if let Some(c) = linked_class(sub, first, &local, i) {
                    local[i] = Some(c);
                    changed = true;
                }
            }
        }
    }
    (0..local.len())
        .filter(|&i| classes[w.start as usize + i].is_none())
        .map(|i| {
            let class = local[i].unwrap_or_else(|| shape_class(sub, first + i as u32));
            (w.start + i as u32, class)
        })
        .collect()
}

/// Protein or nucleic when residue `i`'s backbone atoms bond to a
/// neighbour already of that class: peptide `C(i-1)-N(i)`, or the
/// phosphodiester `O3'(i-1)-P(i)`. Any other bond (a heme's Fe-N, a
/// glycan's link) does not count.
fn linked_class(
    sub: &SubSystem,
    first: u32,
    local: &[Option<ResidueClass>],
    i: usize,
) -> Option<ResidueClass> {
    let here = first + i as u32;
    let prev = i.checked_sub(1).map(|p| (first + p as u32, local[p]));
    let next = (i + 1 < local.len()).then(|| (first + i as u32 + 1, local[i + 1]));
    let joins = |class, from: &str, to: &str| {
        let before = prev.is_some_and(|(p, c)| c == Some(class) && sub.linked(p, from, here, to));
        let after = next.is_some_and(|(n, c)| c == Some(class) && sub.linked(here, from, n, to));
        before || after
    };
    if joins(ResidueClass::Protein, "C", "N") {
        Some(ResidueClass::Protein)
    } else if joins(ResidueClass::Nucleic, "O3'", "P") {
        Some(ResidueClass::Nucleic)
    } else {
        None
    }
}

fn shape_class(sub: &SubSystem, residue: u32) -> ResidueClass {
    if sub.graph(residue).looks_like_lipid() {
        ResidueClass::Lipid
    } else {
        ResidueClass::SmallMolecule
    }
}
