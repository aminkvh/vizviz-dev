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

/// Longest metal-to-donor distance counted as coordination, by element:
/// the typical bond length for that metal (Harding 2006, Acta Cryst D62:
/// 678, whose tables give ~2.0-2.3 A for Zn, Fe and Cu, 2.1 for Mg, 2.4
/// for Ca and Na, 2.8 for K) plus about 0.5 A for coordinate error at
/// ordinary resolution. The cutoffs are this project's rounding of that
/// rule, not a table from the paper.
fn coordination_cutoff(metal: Element) -> f32 {
    match metal.atomic_number() {
        11 | 20 => 3.0,
        12 => 2.6,
        19 => 3.4,
        _ => 2.8,
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
    let reach = 3.4;
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
            let near = c.distance <= coordination_cutoff(element[m as usize]);
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
