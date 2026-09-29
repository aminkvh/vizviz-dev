//! XYZ writer: `element x y z` per atom, one block per frame
//! (`natoms` / comment / atom lines), Angstroms as read. No topology
//! (residues, chains, bonds) survives this format.

use std::io::{self, Write};

use vv_core::fixedbitset::FixedBitSet;
use vv_core::Structure;

use crate::write::{read_frames, residue_atoms};

pub fn write(
    structure: &Structure,
    atoms: Option<&FixedBitSet>,
    frames: &[usize],
    out: &mut impl Write,
) -> io::Result<Vec<String>> {
    let t = &structure.topology;
    let selected: Vec<u32> = t
        .residues
        .iter()
        .flat_map(|r| residue_atoms(r, atoms))
        .collect();
    let id = if t.id.is_empty() {
        "vizviz"
    } else {
        t.id.as_str()
    };
    let held = read_frames(structure, frames)?;
    for (&frame, coords) in frames.iter().zip(&held) {
        let positions = coords.positions();
        writeln!(out, "{}", selected.len())?;
        writeln!(out, "{id} frame {frame}")?;
        for &a in &selected {
            let p = positions[a as usize];
            writeln!(
                out,
                "{:<2} {:>12.6} {:>12.6} {:>12.6}",
                t.element[a as usize].symbol(),
                p.x,
                p.y,
                p.z,
            )?;
        }
    }
    Ok(Vec::new())
}
