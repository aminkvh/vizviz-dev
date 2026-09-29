//! Which alternate location (conformer) of each residue is drawn.
//!
//! All conformers stay in the data; a policy only picks the ones that
//! display, and take part in bonds and surfaces built from the displayed
//! atoms.

use std::fmt;
use std::str::FromStr;

use fixedbitset::FixedBitSet;

use crate::Topology;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum AltlocPolicy {
    /// Per residue, the conformer with the largest summed occupancy
    /// (ties go to the lowest label).
    #[default]
    First,
    All,
    /// This label; a residue without it falls back to `First`.
    Label(u8),
}

impl FromStr for AltlocPolicy {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, String> {
        match s.trim() {
            w if w.eq_ignore_ascii_case("first") => Ok(Self::First),
            w if w.eq_ignore_ascii_case("all") => Ok(Self::All),
            w if w.len() == 1 && w.as_bytes()[0].is_ascii_graphic() => {
                Ok(Self::Label(w.as_bytes()[0]))
            }
            w => Err(format!(
                "`{w}` is not an altloc policy; use `first`, `all` or a label such as `A`"
            )),
        }
    }
}

impl fmt::Display for AltlocPolicy {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::First => f.write_str("first"),
            Self::All => f.write_str("all"),
            Self::Label(l) => write!(f, "{}", *l as char),
        }
    }
}

/// The conformer label a residue shows under `First`, or `None` when it
/// has no alternate locations.
fn first_label(t: &Topology, atoms: std::ops::Range<usize>) -> Option<u8> {
    let mut weight: Vec<(u8, f32)> = Vec::new();
    for a in atoms {
        let label = t.alt_loc[a];
        if label == 0 {
            continue;
        }
        let occupancy = t.occupancy.get(a).copied().unwrap_or(1.0);
        match weight.iter_mut().find(|(l, _)| *l == label) {
            Some((_, w)) => *w += occupancy,
            None => weight.push((label, occupancy)),
        }
    }
    weight
        .into_iter()
        .min_by(|a, b| b.1.total_cmp(&a.1).then(a.0.cmp(&b.0)))
        .map(|(label, _)| label)
}

fn has_label(t: &Topology, atoms: std::ops::Range<usize>, label: u8) -> bool {
    atoms.into_iter().any(|a| t.alt_loc[a] == label)
}

/// Atoms that display under `policy`, or `None` when every atom does.
pub fn visible_atoms(t: &Topology, policy: AltlocPolicy) -> Option<FixedBitSet> {
    if policy == AltlocPolicy::All || t.alt_loc.iter().all(|&a| a == 0) {
        return None;
    }
    let mut keep = FixedBitSet::with_capacity(t.atom_count());
    keep.insert_range(..);
    for res in &t.residues {
        let atoms = res.atoms.start as usize..res.atoms.end as usize;
        let Some(first) = first_label(t, atoms.clone()) else {
            continue;
        };
        let shown = match policy {
            AltlocPolicy::Label(l) if has_label(t, atoms.clone(), l) => l,
            _ => first,
        };
        for a in atoms {
            if t.alt_loc[a] != 0 && t.alt_loc[a] != shown {
                keep.set(a, false);
            }
        }
    }
    Some(keep)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{AtomRow, Element, TopologyBuilder};
    use glam::Vec3;

    fn atom(seq: i32, name: &str, alt: u8, occupancy: f32) -> AtomRow<'_> {
        let mut n = [b' '; 4];
        n[..name.len()].copy_from_slice(name.as_bytes());
        AtomRow {
            element: Element::CARBON,
            name: n,
            serial: 1,
            alt_loc: alt,
            comp: "SER",
            asym: "A",
            auth_asym: "A",
            seq_id: seq,
            auth_seq_id: seq,
            ins_code: 0,
            entity: 1,
            position: Vec3::ZERO,
            occupancy,
            b_factor: 0.0,
            charge: 0,
            hetero: false,
        }
    }

    fn topology() -> Topology {
        let mut b = TopologyBuilder::new();
        b.push(&atom(1, "CA", 0, 1.0));
        b.push(&atom(1, "OG", b'A', 0.3));
        b.push(&atom(1, "OG", b'B', 0.7));
        b.push(&atom(2, "CA", 0, 1.0));
        b.push(&atom(2, "OG", b'A', 0.5));
        b.finish().unwrap().topology.as_ref().clone()
    }

    fn shown(t: &Topology, p: AltlocPolicy) -> Vec<usize> {
        visible_atoms(t, p).map_or_else(|| (0..t.atom_count()).collect(), |m| m.ones().collect())
    }

    #[test]
    fn first_keeps_the_highest_occupancy_conformer_per_residue() {
        assert_eq!(shown(&topology(), AltlocPolicy::First), vec![0, 2, 3, 4]);
    }

    #[test]
    fn a_label_falls_back_to_first_where_a_residue_lacks_it() {
        let t = topology();
        assert_eq!(shown(&t, AltlocPolicy::Label(b'A')), vec![0, 1, 3, 4]);
        assert_eq!(shown(&t, AltlocPolicy::Label(b'B')), vec![0, 2, 3, 4]);
    }

    #[test]
    fn all_and_altloc_free_structures_hide_nothing() {
        assert_eq!(visible_atoms(&topology(), AltlocPolicy::All), None);
        let mut b = TopologyBuilder::new();
        b.push(&atom(1, "CA", 0, 1.0));
        let t = b.finish().unwrap().topology.as_ref().clone();
        assert_eq!(visible_atoms(&t, AltlocPolicy::First), None);
    }

    #[test]
    fn policies_parse_and_print() {
        assert_eq!("first".parse(), Ok(AltlocPolicy::First));
        assert_eq!("ALL".parse(), Ok(AltlocPolicy::All));
        assert_eq!("B".parse(), Ok(AltlocPolicy::Label(b'B')));
        assert!("AB".parse::<AltlocPolicy>().is_err());
        assert_eq!(AltlocPolicy::Label(b'B').to_string(), "B");
    }
}
