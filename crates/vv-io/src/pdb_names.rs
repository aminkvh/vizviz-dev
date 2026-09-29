//! Atom names for the PDB format's four-character column.
//!
//! An mmCIF atom name may be longer. Such a name is cut to four
//! characters, and when the cut collides with another name of the same
//! residue and alternate location, its last character becomes a counter
//! (`1`-`9`, `A`-`Z`), so no two atoms of a residue ever share a name.

use std::collections::{BTreeMap, HashMap, HashSet};

use vv_core::fixedbitset::FixedBitSet;
use vv_core::Topology;

const COUNTERS: &str = "123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";

/// PDB names for every selected atom whose name exceeds four characters,
/// and one warning summarising what was cut or renamed.
pub(crate) fn shorten(
    topology: &Topology,
    atoms: Option<&FixedBitSet>,
) -> (HashMap<u32, String>, Vec<String>) {
    let selected = |a: u32| atoms.is_none_or(|m| m.contains(a as usize));
    let mut by_residue: BTreeMap<u32, Vec<u32>> = BTreeMap::new();
    for &(atom, _) in &topology.long_names {
        if selected(atom) {
            by_residue
                .entry(topology.residue_index[atom as usize])
                .or_default()
                .push(atom);
        }
    }
    let mut out = HashMap::new();
    let (mut cut, mut renamed) = (0usize, 0usize);
    let mut first: Option<(String, String)> = None;
    for (residue, long) in by_residue {
        let range = topology.residues[residue as usize].atoms.clone();
        let mut taken: HashSet<(String, u8)> = range
            .filter(|a| selected(*a) && !long.contains(a))
            .map(|a| {
                (
                    topology.atom_name(a as usize).to_string(),
                    topology.alt_loc[a as usize],
                )
            })
            .collect();
        for atom in long {
            let alt = topology.alt_loc[atom as usize];
            let full = topology.atom_name(atom as usize);
            let mut name: String = full.chars().take(4).collect();
            if taken.contains(&(name.clone(), alt)) {
                name = unused_counter_name(&name, alt, &taken);
                renamed += 1;
            }
            cut += 1;
            first.get_or_insert_with(|| (full.to_string(), name.clone()));
            taken.insert((name.clone(), alt));
            out.insert(atom, name);
        }
    }
    (out, warning(cut, renamed, first))
}

fn unused_counter_name(name: &str, alt: u8, taken: &HashSet<(String, u8)>) -> String {
    let stem: String = name.chars().take(3).collect();
    COUNTERS
        .chars()
        .map(|c| format!("{stem}{c}"))
        .find(|n| !taken.contains(&(n.clone(), alt)))
        .unwrap_or_else(|| name.to_string())
}

fn warning(cut: usize, renamed: usize, first: Option<(String, String)>) -> Vec<String> {
    let Some((from, to)) = first else {
        return Vec::new();
    };
    let mut message = format!(
        "{cut} atom name(s) longer than 4 characters were truncated for PDB (first: `{from}` written as `{to}`)"
    );
    if renamed > 0 {
        message.push_str(&format!(
            "; {renamed} renamed to stay unique within their residue"
        ));
    }
    vec![message]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn counter_replaces_the_last_character_until_free() {
        let mut taken = HashSet::new();
        taken.insert(("C1A1".to_string(), 0u8));
        taken.insert(("C1A2".to_string(), 0u8));
        assert_eq!(unused_counter_name("C1A1", 0, &taken), "C1A3");
        assert_eq!(unused_counter_name("C1A1", b'A', &taken), "C1A1");
    }

    #[test]
    fn no_warning_when_nothing_was_cut() {
        assert!(warning(0, 0, None).is_empty());
    }
}
