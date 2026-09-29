//! Shared helpers for the structure writers (`pdb_write`, `mmcif_write`,
//! `xyz_write`, `pqr_write`, `gro_write`).

use std::io;
use std::sync::Arc;

use vv_core::fixedbitset::FixedBitSet;
use vv_core::{CoordSet, SecondaryStructure, Structure, StructureError, Topology};

/// Atom indices of `residue` selected by `mask` (all of them when `mask`
/// is `None`), in ascending order.
pub(crate) fn residue_atoms(residue: &vv_core::ResidueRec, mask: Option<&FixedBitSet>) -> Vec<u32> {
    residue
        .atoms
        .clone()
        .filter(|&a| mask.is_none_or(|m| m.contains(a as usize)))
        .collect()
}

/// A maximal run of consecutive written residues in one chain sharing a
/// helix or strand assignment (residue indices).
pub(crate) struct SsRun {
    pub ss: SecondaryStructure,
    pub first: u32,
    pub last: u32,
}

pub(crate) fn ss_runs(topology: &Topology, mask: Option<&FixedBitSet>) -> Vec<SsRun> {
    let mut runs: Vec<SsRun> = Vec::new();
    for (i, res) in topology.residues.iter().enumerate() {
        let structured = matches!(
            res.ss,
            SecondaryStructure::Helix | SecondaryStructure::Strand
        );
        if !structured || residue_atoms(res, mask).is_empty() {
            continue;
        }
        let i = i as u32;
        let extends = runs.last().is_some_and(|r| {
            let prev = &topology.residues[r.last as usize];
            r.ss == res.ss && prev.chain == res.chain && r.last + 1 == i
        });
        if extends {
            runs.last_mut().expect("checked").last = i;
        } else {
            runs.push(SsRun {
                ss: res.ss,
                first: i,
                last: i,
            });
        }
    }
    runs
}

/// Reads every requested frame. A streamed frame that fails to read is an
/// I/O error, never silently frame 0 (unlike `Structure::frame`).
pub(crate) fn read_frames(
    structure: &Structure,
    frames: &[usize],
) -> io::Result<Vec<Arc<CoordSet>>> {
    frames
        .iter()
        .map(|&f| structure.try_frame(f).map_err(structure_err))
        .collect()
}

fn structure_err(e: StructureError) -> io::Error {
    io::Error::other(e.to_string())
}
