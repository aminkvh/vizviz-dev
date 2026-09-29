//! What the strip shows a row for: one `(chain name, residues)` per chain
//! name, and the one-letter sequence of a row.

use std::ops::Range;

use vv_core::residue_class::ResidueClass;
pub use vv_core::seqfeat::one_letter;
use vv_core::Topology;
use vv_scene::{Scene, StructureId};

pub struct Row {
    pub structure: StructureId,
    pub label: String,
    pub residues: Range<u32>,
}

pub fn rows_of(scene: &Scene) -> Vec<Row> {
    let mut rows = Vec::new();
    for (id, loaded) in scene.structures() {
        for (name, residues) in chain_rows(&loaded.structure.topology) {
            rows.push(Row {
                structure: id,
                label: format!("{} {name}", loaded.label),
                residues,
            });
        }
    }
    rows
}

/// One `(chain name, residues)` row per chain name. A PDB `TER` splits a
/// chain letter into polymer and water/ligand records: adjacent records
/// with one name join into one row, and a later non-polymer record that
/// repeats a name is left out rather than shown as a second "A".
pub fn chain_rows(top: &Topology) -> Vec<(String, Range<u32>)> {
    let mut rows: Vec<(String, Range<u32>)> = Vec::new();
    for (ci, chain) in top.chains.iter().enumerate() {
        let name = top.chain_name(ci);
        let same = rows.iter_mut().rev().find(|r| r.0 == name);
        match same {
            Some(row) if row.1.end == chain.residues.start => row.1.end = chain.residues.end,
            Some(_) if !is_polymer(top, chain.residues.start) => {}
            _ => rows.push((name.to_string(), chain.residues.clone())),
        }
    }
    rows
}

pub fn is_polymer(top: &Topology, residue: u32) -> bool {
    matches!(
        top.residue_class(residue as usize),
        ResidueClass::Protein | ResidueClass::Nucleic
    )
}

pub fn is_protein(top: &Topology, residue: u32) -> bool {
    top.residue_class(residue as usize) == ResidueClass::Protein
}

/// The letters of `residues` (`x` for anything without a code) and the
/// numbering breaks: `true` where the author numbering jumps forward, so
/// residues before and after it are not neighbours in the chain.
pub fn letters_and_breaks(top: &Topology, residues: Range<u32>) -> (Vec<u8>, Vec<bool>) {
    let letters = residues
        .clone()
        .map(|r| one_letter(top.residue_name(r as usize)).map_or(b'x', |c| c as u8))
        .collect();
    let mut prev: Option<i32> = None;
    let breaks = residues
        .map(|r| {
            let seq = top.residues[r as usize].auth_seq_id;
            let jump = prev.is_some_and(|p| seq - p > 1);
            prev = Some(seq);
            jump
        })
        .collect();
    (letters, breaks)
}

/// A residue's number as authors write it: `52A` for insertion code `A`.
pub fn residue_number(top: &Topology, residue: u32) -> String {
    let rec = &top.residues[residue as usize];
    match rec.ins_code {
        0 => rec.auth_seq_id.to_string(),
        c => format!("{}{}", rec.auth_seq_id, c as char),
    }
}

/// `NAG 1301 (B)`: name, author number and author chain.
pub fn describe_residue(top: &Topology, residue: u32) -> String {
    let rec = &top.residues[residue as usize];
    let chain = top.names.get(top.chains[rec.chain as usize].auth_asym);
    format!(
        "{} {} ({chain})",
        top.residue_name(residue as usize),
        residue_number(top, residue)
    )
}

/// Residues of the actively selected structure that contain at least one
/// selected atom. Recomputed every frame: it's one pass over the selected
/// atoms, not over the structure.
pub fn selected_residues(scene: &Scene) -> Option<(StructureId, Vec<bool>)> {
    let active = scene.active_selection()?;
    let top = &scene.structure(active.structure)?.structure.topology;
    let mut flags = vec![false; top.residue_count()];
    for atom in active.mask.ones() {
        flags[top.residue_index[atom] as usize] = true;
    }
    Some((active.structure, flags))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_ter_split_chain_letter_is_one_sequence_row() {
        let atom = |serial: u32, res: &str, chain: char, seq: u32| {
            format!(
                "{:<6}{serial:>5}  CA  {res:<4}{chain}{seq:>4}    {:8.3}{:8.3}{:8.3}  1.00  0.00           C\n",
                if res == "HOH" { "HETATM" } else { "ATOM" },
                seq as f32 * 3.8,
                0.0,
                0.0
            )
        };
        let text =
            atom(1, "ALA", 'A', 1) + &atom(2, "GLY", 'A', 2) + "TER\n" + &atom(3, "HOH", 'A', 3);
        let top = vv_io::pdb::parse(text.as_bytes()).unwrap().topology.clone();
        assert_eq!(top.chain_count(), 2);
        assert_eq!(chain_rows(&top), [("A".to_string(), 0..3)]);
    }

    #[test]
    fn numbering_jumps_are_breaks_and_insertion_codes_are_not() {
        let atom = |serial: u32, seq: u32, ins: char| {
            format!(
                "ATOM  {serial:>5}  CA  ALA A{seq:>4}{ins}   {:8.3}{:8.3}{:8.3}  1.00  0.00           C\n",
                serial as f32 * 3.8,
                0.0,
                0.0
            )
        };
        let text = atom(1, 1, ' ') + &atom(2, 1, 'A') + &atom(3, 2, ' ') + &atom(4, 5, ' ');
        let top = vv_io::pdb::parse(text.as_bytes()).unwrap().topology.clone();
        let (letters, breaks) = letters_and_breaks(&top, 0..4);
        assert_eq!(letters, b"AAAA");
        assert_eq!(breaks, [false, false, false, true]);
        assert_eq!(residue_number(&top, 1), "1A");
    }
}
