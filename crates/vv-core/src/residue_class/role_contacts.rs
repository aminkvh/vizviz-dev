//! The two spatial tests behind roles: which lipids form a membrane, and
//! which glycans are covalently attached to a polymer.

use glam::Vec3;
use rayon::prelude::*;

use super::ResidueClass;
use crate::{Grid, Topology};

/// Heavy atoms of two lipids closer than this are in contact.
const CONTACT: f32 = 4.0;

/// Contacts recorded per lipid; a lipid with more is inside a network
/// anyway, and a bilayer stays connected through any six of them.
const CONTACTS_MAX: usize = 6;

/// A contact network of at least this many lipid residues is a membrane
/// (or micelle): bound lipids form networks of a few, bilayer patches of
/// dozens or more.
pub(super) const MEMBRANE_MIN_LIPIDS: usize = 12;

/// Heavy atoms of different residues closer than this are bonded: the
/// longest glycosidic or glycan-protein bond (C-S, 1.82 A) plus slack, well
/// under the 2.6 A of a hydrogen bond.
const LINK_MAX: f32 = 1.9;

const NO_SLOT: u32 = u32::MAX;

/// Position of each listed residue within the list.
struct Slots(Vec<u32>);

impl Slots {
    fn new(residue_count: usize, listed: &[u32]) -> Self {
        let mut slots = vec![NO_SLOT; residue_count];
        for (k, &r) in listed.iter().enumerate() {
            slots[r as usize] = k as u32;
        }
        Slots(slots)
    }

    fn get(&self, residue: u32) -> Option<usize> {
        Some(self.0[residue as usize])
            .filter(|&s| s != NO_SLOT)
            .map(|s| s as usize)
    }
}

struct DisjointSets {
    parent: Vec<u32>,
}

impl DisjointSets {
    fn new(n: usize) -> Self {
        Self {
            parent: (0..n as u32).collect(),
        }
    }

    fn find(&mut self, mut x: usize) -> usize {
        while self.parent[x] as usize != x {
            self.parent[x] = self.parent[self.parent[x] as usize];
            x = self.parent[x] as usize;
        }
        x
    }

    fn union(&mut self, a: usize, b: usize) {
        let (ra, rb) = (self.find(a), self.find(b));
        self.parent[ra] = rb as u32;
    }

    fn sizes(&mut self) -> Vec<u32> {
        let mut sizes = vec![0u32; self.parent.len()];
        for x in 0..self.parent.len() {
            let root = self.find(x);
            sizes[root] += 1;
        }
        sizes
    }
}

fn heavy_atoms(t: &Topology, residues: &[u32]) -> Vec<u32> {
    residues
        .par_iter()
        .flat_map_iter(|&r| {
            t.residues[r as usize]
                .atoms
                .clone()
                .filter(|&a| !t.element[a as usize].is_hydrogen())
        })
        .collect()
}

/// For each of `lipids` (residue indices), whether it is in a contact
/// network of at least [`MEMBRANE_MIN_LIPIDS`] lipids.
pub(super) fn membrane(t: &Topology, positions: &[Vec3], lipids: &[u32]) -> Vec<bool> {
    if lipids.len() < MEMBRANE_MIN_LIPIDS {
        return vec![false; lipids.len()];
    }
    let slots = Slots::new(t.residues.len(), lipids);
    let grid = Grid::build(positions, &heavy_atoms(t, lipids), CONTACT);
    let edges: Vec<[usize; 2]> = (0..lipids.len())
        .into_par_iter()
        .flat_map_iter(|k| {
            let touching = touching_slots(t, positions, &grid, &slots, lipids[k], k);
            touching.into_iter().map(move |other| [k, other])
        })
        .collect();
    let mut sets = DisjointSets::new(lipids.len());
    for [a, b] in edges {
        sets.union(a, b);
    }
    let sizes = sets.sizes();
    (0..lipids.len())
        .map(|k| sizes[sets.find(k)] as usize >= MEMBRANE_MIN_LIPIDS)
        .collect()
}

/// Slots of the lipids touching `residue`, each once. Stops scanning at
/// [`CONTACTS_MAX`], so a lipid deep in a membrane costs a few queries.
fn touching_slots(
    t: &Topology,
    positions: &[Vec3],
    grid: &Grid,
    slots: &Slots,
    residue: u32,
    slot: usize,
) -> Vec<usize> {
    let mut touching: Vec<usize> = Vec::new();
    for a in t.residues[residue as usize].atoms.clone() {
        if t.element[a as usize].is_hydrogen() {
            continue;
        }
        grid.for_each_within(positions, positions[a as usize], CONTACT, |b, _| {
            let other = slots.get(t.residue_index[b as usize]);
            touching.extend(other.filter(|&o| o != slot));
        });
        touching.sort_unstable();
        touching.dedup();
        if touching.len() >= CONTACTS_MAX {
            break;
        }
    }
    touching
}

/// For each of `glycans` (residue indices), whether it belongs to a glycan
/// molecule with a covalent link to a protein or nucleic residue.
pub(super) fn attached_glycans(
    t: &Topology,
    classes: &[ResidueClass],
    positions: &[Vec3],
    glycans: &[u32],
) -> Vec<bool> {
    if glycans.is_empty() {
        return Vec::new();
    }
    let links = match &t.md_bonds {
        Some(bonds) => links_from_bonds(t, classes, bonds),
        None => links_from_geometry(t, classes, positions, glycans),
    };
    let slots = Slots::new(t.residues.len(), glycans);
    let mut sets = DisjointSets::new(glycans.len());
    let mut bound = vec![false; glycans.len()];
    for [a, b] in links {
        match (slots.get(a), slots.get(b)) {
            (Some(x), Some(y)) => sets.union(x, y),
            (Some(x), None) | (None, Some(x)) => bound[x] = true,
            (None, None) => {}
        }
    }
    let mut molecule_bound = vec![false; glycans.len()];
    for (x, &b) in bound.iter().enumerate() {
        let root = sets.find(x);
        molecule_bound[root] |= b;
    }
    (0..glycans.len())
        .map(|x| molecule_bound[sets.find(x)])
        .collect()
}

fn is_polymer_or_glycan(class: ResidueClass) -> bool {
    matches!(
        class,
        ResidueClass::Protein | ResidueClass::Nucleic | ResidueClass::Glycan
    )
}

/// Residue pairs joined by a bond of an MD topology, at least one a glycan.
fn links_from_bonds(t: &Topology, classes: &[ResidueClass], bonds: &[[u32; 2]]) -> Vec<[u32; 2]> {
    bonds
        .par_iter()
        .filter_map(|&[a, b]| {
            let (ra, rb) = (t.residue_index[a as usize], t.residue_index[b as usize]);
            let (ca, cb) = (classes[ra as usize], classes[rb as usize]);
            let glycan = ca == ResidueClass::Glycan || cb == ResidueClass::Glycan;
            (ra != rb && glycan && is_polymer_or_glycan(ca) && is_polymer_or_glycan(cb))
                .then_some([ra, rb])
        })
        .collect()
}

/// Residue pairs with heavy atoms within [`LINK_MAX`], at least one a
/// glycan: the grid holds only glycan atoms, so the many polymer atoms
/// far from any glycan cost one failed lookup each.
fn links_from_geometry(
    t: &Topology,
    classes: &[ResidueClass],
    positions: &[Vec3],
    glycans: &[u32],
) -> Vec<[u32; 2]> {
    let grid = Grid::build(positions, &heavy_atoms(t, glycans), LINK_MAX);
    classes
        .par_iter()
        .enumerate()
        .filter(|&(_, &c)| is_polymer_or_glycan(c))
        .flat_map_iter(|(r, _)| {
            let mut links = Vec::new();
            for a in t.residues[r].atoms.clone() {
                if t.element[a as usize].is_hydrogen() {
                    continue;
                }
                grid.for_each_within(positions, positions[a as usize], LINK_MAX, |b, _| {
                    let other = t.residue_index[b as usize];
                    if other as usize != r {
                        links.push([r as u32, other]);
                    }
                });
            }
            links
        })
        .collect()
}
