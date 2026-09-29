//! PQR writer: PDB-shaped ATOM/HETATM records with charge and radius in
//! place of occupancy/B-factor, whitespace-separated (PQR dropped fixed
//! columns so names and chain ids can run longer, unlike PDB).
//!
//! Neither a partial charge nor a per-atom radius is part of the data
//! model, so both are approximations, documented here rather than
//! silently guessed elsewhere: charge is the file's formal integer
//! charge (`Topology::charge`), not a force-field partial charge; radius
//! is each element's van der Waals radius (`Element::vdw_radius`, Bondi
//! 1964 / Rowland & Taylor 1996), not a per-atom-type radius from a
//! force field.

use std::io::{self, Write};

use vv_core::fixedbitset::FixedBitSet;
use vv_core::{flags, Structure};

use crate::write::{read_frames, residue_atoms};

pub fn write(
    structure: &Structure,
    atoms: Option<&FixedBitSet>,
    frames: &[usize],
    out: &mut impl Write,
) -> io::Result<Vec<String>> {
    let t = &structure.topology;
    let held = read_frames(structure, frames)?;
    for coords in &held {
        let positions = coords.positions();
        let mut serial = 1u32;
        for res in &t.residues {
            let chain = t.names.get(t.chains[res.chain as usize].auth_asym);
            let comp = t.names.get(res.comp);
            for a in residue_atoms(res, atoms) {
                let a = a as usize;
                let group = if t.flags.get(a).is_some_and(|f| f & flags::HETERO != 0) {
                    "HETATM"
                } else {
                    "ATOM"
                };
                let p = positions[a];
                // A literal space before every field: fixed widths alone
                // don't separate them (`ATOM` + a 5-digit serial, or two
                // 8-column coordinates that both fill their width, would
                // otherwise run together with nothing between).
                writeln!(
                    out,
                    "{group} {serial:>5} {} {comp} {chain} {} {:.3} {:.3} {:.3} {:.4} {:.4}",
                    t.atom_name(a),
                    res.auth_seq_id,
                    p.x,
                    p.y,
                    p.z,
                    t.charge.get(a).copied().unwrap_or(0) as f32,
                    t.element[a].vdw_radius(),
                )?;
                serial += 1;
            }
        }
    }
    Ok(Vec::new())
}
