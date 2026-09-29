//! GROMACS GRO writer: fixed columns `%5d%-5s%5s%5d%8.3f%8.3f%8.3f`
//! (residue number, residue name, atom name, atom number, x/y/z), nm
//! rather than the rest of this crate's Angstroms, one frame per block,
//! each ending in a box vector line.
//!
//! GRO's residue/atom number fields are 5 digits: real GROMACS tools
//! wrap them modulo 100000 rather than extending them (unlike PDB's
//! hybrid-36), so this writer does too. No unit cell is part of the data
//! model, so the box is the exported atoms' bounding box (in that frame)
//! padded by 1 nm on every side, not a real periodic cell.

use std::io::{self, Write};

use vv_core::fixedbitset::FixedBitSet;
use vv_core::glam::Vec3;
use vv_core::Structure;

use crate::write::{read_frames, residue_atoms};

const ANGSTROM_PER_NM: f32 = 10.0;
const BOX_PADDING_NM: f32 = 1.0;

fn wrap(n: i64) -> i64 {
    n.rem_euclid(100_000)
}

fn box_line(positions: &[Vec3], atoms: &[u32]) -> String {
    let (mut lo, mut hi) = (Vec3::splat(f32::MAX), Vec3::splat(f32::MIN));
    for &a in atoms {
        let p = positions[a as usize] / ANGSTROM_PER_NM;
        lo = lo.min(p);
        hi = hi.max(p);
    }
    if atoms.is_empty() {
        lo = Vec3::ZERO;
        hi = Vec3::ZERO;
    }
    let size = (hi - lo) + Vec3::splat(2.0 * BOX_PADDING_NM);
    format!("{:>10.5}{:>10.5}{:>10.5}", size.x, size.y, size.z)
}

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
        writeln!(out, "{id}, frame {frame}")?;
        writeln!(out, "{:>5}", selected.len())?;
        let mut atom_num = 1i64;
        for res in &t.residues {
            let resname = t.names.get(res.comp);
            let resname = &resname[..resname.len().min(5)];
            for &a in &residue_atoms(res, atoms) {
                let name = t.atom_name(a as usize);
                let name = &name[..name.len().min(5)];
                let p = positions[a as usize] / ANGSTROM_PER_NM;
                writeln!(
                    out,
                    "{:>5}{:<5}{:>5}{:>5}{:>8.3}{:>8.3}{:>8.3}",
                    wrap(res.auth_seq_id as i64),
                    resname,
                    name,
                    wrap(atom_num),
                    p.x,
                    p.y,
                    p.z,
                )?;
                atom_num += 1;
            }
        }
        writeln!(out, "{}", box_line(positions, &selected))?;
    }
    Ok(Vec::new())
}
