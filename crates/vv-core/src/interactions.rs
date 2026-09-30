//! Non-covalent contacts between atoms of one structure, for the overlay
//! that draws them as dashes: hydrogen bonds, metal coordination and salt
//! bridges. Each returns atom pairs `[a, b]`; `visible` limits both ends
//! to the shown conformer.

use glam::Vec3;

use crate::analysis::{neighbor_pairs_into, Contact};
use crate::bonds::BondTable;
use crate::element::Element;
use crate::residue_class::ResidueClass;
use crate::topology::Topology;

/// A kind of contact the overlay can draw.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum InteractionKind {
    Hbond,
    Metal,
    SaltBridge,
}

impl InteractionKind {
    pub const ALL: [InteractionKind; 3] = [Self::Hbond, Self::Metal, Self::SaltBridge];

    pub fn name(self) -> &'static str {
        match self {
            Self::Hbond => "hbond",
            Self::Metal => "metal",
            Self::SaltBridge => "saltbridge",
        }
    }

    pub fn parse(name: &str) -> Option<Self> {
        Self::ALL.into_iter().find(|k| k.name() == name)
    }
}

/// Which atoms a contact may end on.
pub type Visible<'a> = &'a dyn Fn(u32) -> bool;

/// Longest carboxylate O to basic N distance of a salt bridge: Barlow &
/// Thornton 1983, J Mol Biol 168:867.
pub const SALT_BRIDGE_DISTANCE: f32 = 4.0;

/// Metal-ligand distance statistics of high-resolution PDB entries, Zheng
/// et al. 2008, J Inorg Biochem 102:1765-1776, Table 3 (PDB-HR): metal
/// atomic number, ligand element, mean and SD in angstrom. Oxygen rows are
/// the carbonyl/carboxylate/hydroxyl ligand (waters are never drawn). Fe
/// and Cu take the wider of their two oxidation states.
const LIGAND_DISTANCES: &[(u8, u8, f32, f32)] = &[
    (11, 8, 2.43, 0.20),
    (12, 8, 2.21, 0.25),
    (19, 8, 2.76, 0.14),
    (20, 8, 2.37, 0.12),
    (25, 7, 2.20, 0.13),
    (25, 8, 2.15, 0.15),
    (26, 7, 2.16, 0.13),
    (26, 8, 2.14, 0.19),
    (26, 16, 2.29, 0.04),
    (27, 7, 2.07, 0.13),
    (27, 8, 2.09, 0.12),
    (27, 16, 2.33, 0.03),
    (28, 7, 1.99, 0.14),
    (28, 8, 2.17, 0.22),
    (28, 16, 2.24, 0.15),
    (29, 7, 2.04, 0.15),
    (29, 8, 2.24, 0.34),
    (29, 16, 2.36, 0.27),
    (30, 7, 2.07, 0.11),
    (30, 8, 2.08, 0.20),
    (30, 16, 2.32, 0.06),
];

/// Mean plus this many SD is the longest distance counted as coordination
/// (this project's choice; the source gives only the distributions).
const COORDINATION_SDS: f32 = 3.0;

/// The source's own contact radius, used for metals its tables omit.
const UNTABULATED_REACH: f32 = 3.0;

fn row_cutoff(&(_, _, mean, sd): &(u8, u8, f32, f32)) -> f32 {
    mean + COORDINATION_SDS * sd
}

/// The longest `metal`-`donor` distance counted as coordination. A pairing
/// the tables lack (Mn-S, Ca-N, ...) takes the metal's widest tabulated
/// cutoff: rare pairings still show, without the looser untabulated reach.
fn coordination_cutoff(metal: Element, donor: Element) -> f32 {
    let (m, d) = (metal.atomic_number(), donor.atomic_number());
    let rows = || LIGAND_DISTANCES.iter().filter(move |row| row.0 == m);
    match rows().find(|row| row.1 == d) {
        Some(row) => row_cutoff(row),
        None => rows()
            .map(row_cutoff)
            .reduce(f32::max)
            .unwrap_or(UNTABULATED_REACH),
    }
}

fn is_donor(element: Element) -> bool {
    matches!(
        element,
        Element::NITROGEN | Element::OXYGEN | Element::SULFUR
    )
}

fn residue_of(topology: &Topology, atom: u32) -> u32 {
    topology.residue_index[atom as usize]
}

fn in_water(topology: &Topology, atom: u32) -> bool {
    topology.residue_class(residue_of(topology, atom) as usize) == ResidueClass::Water
}

fn contacts_among(positions: &[Vec3], atoms: &[u32], cutoff: f32) -> Vec<Contact> {
    let mut out = Vec::new();
    neighbor_pairs_into(positions, atoms, cutoff, &mut out);
    out
}

/// Hydrogen bonds (`crate::hbond`, Baker & Hubbard 1984) between different
/// residues, leaving out water: a dash to an undrawn water is noise.
pub fn hydrogen_bonds(
    topology: &Topology,
    bonds: &BondTable,
    positions: &[Vec3],
    visible: Visible,
) -> Vec<[u32; 2]> {
    let mut found = Vec::new();
    crate::hbond::hydrogen_bonds_into(&topology.element, bonds, positions, &mut found);
    found
        .into_iter()
        .filter(|h| visible(h.a) && visible(h.b))
        .filter(|h| residue_of(topology, h.a) != residue_of(topology, h.b))
        .filter(|h| !in_water(topology, h.a) && !in_water(topology, h.b))
        .filter(|h| bonds.pairs.binary_search(&[h.a, h.b]).is_err())
        .map(|h| [h.a, h.b])
        .collect()
}

/// Each metal atom paired with every N, O or S within its element's
/// coordination cutoff, from other residues (a heme iron's own porphyrin
/// nitrogens are covalent to it in the model) and not water.
pub fn metal_coordination(
    topology: &Topology,
    positions: &[Vec3],
    visible: Visible,
) -> Vec<[u32; 2]> {
    let element = &topology.element;
    let atoms: Vec<u32> = (0..element.len() as u32)
        .filter(|&a| visible(a))
        .filter(|&a| element[a as usize].is_metal() || is_donor(element[a as usize]))
        .collect();
    let reach = LIGAND_DISTANCES
        .iter()
        .map(|&(_, _, mean, sd)| mean + COORDINATION_SDS * sd)
        .fold(UNTABULATED_REACH, f32::max);
    contacts_among(positions, &atoms, reach)
        .into_iter()
        .filter_map(|c| {
            let (m, d) = match (
                element[c.a as usize].is_metal(),
                element[c.b as usize].is_metal(),
            ) {
                (true, false) => (c.a, c.b),
                (false, true) => (c.b, c.a),
                _ => return None,
            };
            let near = c.distance <= coordination_cutoff(element[m as usize], element[d as usize]);
            let other = residue_of(topology, m) != residue_of(topology, d);
            (near && other && !in_water(topology, d)).then_some([m, d])
        })
        .collect()
}

/// Whether `atom` is a carboxylate oxygen of Asp/Glu (`Some(true)`) or a
/// charged nitrogen of Lys/Arg/His (`Some(false)`).
fn charge_side(topology: &Topology, atom: u32) -> Option<bool> {
    let residue = topology.residue_name(residue_of(topology, atom) as usize);
    let name = topology.atom_name(atom as usize);
    match (residue, name) {
        ("ASP", "OD1" | "OD2") | ("GLU", "OE1" | "OE2") => Some(true),
        ("LYS", "NZ") | ("ARG", "NE" | "NH1" | "NH2") | ("HIS", "ND1" | "NE2") => Some(false),
        _ => None,
    }
}

/// One line per acid/base residue pair: their closest carboxylate O to
/// basic N within [`SALT_BRIDGE_DISTANCE`].
pub fn salt_bridges(topology: &Topology, positions: &[Vec3], visible: Visible) -> Vec<[u32; 2]> {
    let atoms: Vec<u32> = (0..topology.atom_count() as u32)
        .filter(|&a| visible(a) && charge_side(topology, a).is_some())
        .collect();
    let mut best: std::collections::BTreeMap<(u32, u32), (f32, [u32; 2])> = Default::default();
    for c in contacts_among(positions, &atoms, SALT_BRIDGE_DISTANCE) {
        let (sa, sb) = (charge_side(topology, c.a), charge_side(topology, c.b));
        let (acid, base) = match (sa, sb) {
            (Some(true), Some(false)) => (c.a, c.b),
            (Some(false), Some(true)) => (c.b, c.a),
            _ => continue,
        };
        let key = (residue_of(topology, acid), residue_of(topology, base));
        let entry = best.entry(key).or_insert((f32::INFINITY, [acid, base]));
        if c.distance < entry.0 {
            *entry = (c.distance, [acid, base]);
        }
    }
    best.into_values().map(|(_, pair)| pair).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cutoff(metal: u8, donor: u8) -> f32 {
        let e = |z| Element::from_atomic_number(z).unwrap();
        coordination_cutoff(e(metal), e(donor))
    }

    #[test]
    fn an_untabulated_pairing_takes_the_metals_widest_cutoff() {
        let mn_n: f32 = 2.20 + 3.0 * 0.13;
        let mn_o = 2.15 + 3.0 * 0.15;
        assert!((cutoff(25, 16) - mn_n.max(mn_o)).abs() < 1e-5, "Mn-S");
        assert!((cutoff(20, 7) - (2.37 + 3.0 * 0.12)).abs() < 1e-5, "Ca-N");
        assert!((cutoff(30, 16) - (2.32 + 3.0 * 0.06)).abs() < 1e-5, "Zn-S");
        assert_eq!(cutoff(79, 16), UNTABULATED_REACH, "Au is untabulated");
    }
}
