//! Secondary-structure ranges shared by the PDB (`HELIX`/`SHEET`) and mmCIF
//! (`_struct_conf`/`_struct_sheet_range`) readers.

use vv_core::{flags, SecondaryStructure, Topology};

/// A helix or strand: chain name, then `(seq_id, insertion code)` of the
/// first and last residue.
pub(crate) struct SsRange {
    pub chain: String,
    pub first: (i32, u8),
    pub last: (i32, u8),
}

/// Marks each range's residues by walking the chain in file order from the
/// first residue to the last, so insertion codes (`52`, `52A`, `53`) fall
/// inside the range as the formats intend. A range whose end residues are
/// missing from the file falls back to comparing residue numbers, skipping
/// HETATM residues (ligands and waters reuse low numbers).
pub(crate) fn apply(topology: &mut Topology, ranges: &[SsRange], ss: SecondaryStructure) {
    for range in ranges {
        for r in residues_in(topology, range) {
            topology.residues[r].ss = ss;
        }
    }
}

fn residues_in(topology: &Topology, range: &SsRange) -> Vec<usize> {
    let Some(asym) = topology.names.lookup(&range.chain) else {
        return Vec::new();
    };
    let residues: Vec<usize> = topology
        .chains
        .iter()
        .filter(|c| c.label_asym == asym)
        .flat_map(|c| c.residues.clone())
        .map(|r| r as usize)
        .collect();
    let key = |r: usize| (topology.residues[r].seq_id, topology.residues[r].ins_code);
    let start = residues.iter().position(|&r| key(r) == range.first);
    let end = start.and_then(|s| {
        residues[s..]
            .iter()
            .position(|&r| key(r) == range.last)
            .map(|e| s + e)
    });
    if let (Some(s), Some(e)) = (start, end) {
        return residues[s..=e].to_vec();
    }
    residues
        .into_iter()
        .filter(|&r| {
            let first_atom = topology.residues[r].atoms.start as usize;
            let hetero = topology.flags[first_atom] & flags::HETERO != 0;
            !hetero && key(r) >= range.first && key(r) <= range.last
        })
        .collect()
}
