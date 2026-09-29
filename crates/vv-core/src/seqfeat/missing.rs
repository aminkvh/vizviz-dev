//! Unobserved residues: those the file lists (`REMARK 465`, mmCIF
//! `_pdbx_unobs_or_zero_occ_residues`) but that have no coordinates,
//! placed against the residues that do.

use std::collections::HashSet;
use std::ops::Range;

use crate::{ResidueClass, Topology};

/// One unobserved residue, named as the file numbers it.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Unobserved {
    pub name: String,
    pub seq: i32,
    pub ins: u8,
}

/// Unobserved residues that sit just before modeled residue `before`
/// (a topology residue index), or after the last one when `before` is the
/// end of the chain range.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Gap {
    pub before: u32,
    pub residues: Vec<Unobserved>,
}

/// The polymer residues of `residues` (a chain's rows of one name) that
/// the file lists as unobserved, grouped into gaps between modeled ones.
/// Empty when the file carries no such table. Matches on the author chain
/// id and author numbering, the same numbering the strip shows.
pub fn unobserved(top: &Topology, residues: Range<u32>) -> Vec<Gap> {
    let Some(first) = residues.clone().next() else {
        return Vec::new();
    };
    if !matches!(
        top.residue_class(first as usize),
        ResidueClass::Protein | ResidueClass::Nucleic
    ) {
        return Vec::new();
    }
    let chain = top
        .names
        .get(top.chains[top.residues[first as usize].chain as usize].auth_asym);
    let listed = listed_for_chain(top, chain);
    place(top, residues, listed)
}

fn listed_for_chain(top: &Topology, chain: &str) -> Vec<Unobserved> {
    let Some(cat) = top.annotations.category("pdbx_unobs_or_zero_occ_residues") else {
        return Vec::new();
    };
    let mut seen: HashSet<(i32, u8)> = HashSet::new();
    (0..cat.rows.len())
        .filter(|&r| cat.get("auth_asym_id", r) == Some(chain))
        .filter(|&r| cat.get("polymer_flag", r).is_none_or(|f| f == "Y"))
        .filter(|&r| cat.get("occupancy_flag", r).is_none_or(|f| f == "0"))
        .filter(|&r| cat.get("PDB_model_num", r).is_none_or(|m| m == "1"))
        .filter_map(|r| {
            let seq = cat.get("auth_seq_id", r)?.parse().ok()?;
            let ins = cat.get("PDB_ins_code", r).map_or(0, |s| s.as_bytes()[0]);
            seen.insert((seq, ins)).then(|| Unobserved {
                name: cat.get("auth_comp_id", r).unwrap_or("?").to_string(),
                seq,
                ins,
            })
        })
        .collect()
}

/// Walks the modeled residues in order, assigning each listed residue to
/// the first modeled one that follows it in `(seq, ins)` order.
fn place(top: &Topology, residues: Range<u32>, mut listed: Vec<Unobserved>) -> Vec<Gap> {
    listed.sort_by_key(|u| (u.seq, u.ins));
    let key = |r: u32| {
        let rec = &top.residues[r as usize];
        (rec.auth_seq_id, rec.ins_code)
    };
    let mut at = residues.start;
    let mut gaps: Vec<Gap> = Vec::new();
    for u in listed {
        while at < residues.end && key(at) < (u.seq, u.ins) {
            at += 1;
        }
        match gaps.last_mut() {
            Some(g) if g.before == at => g.residues.push(u),
            _ => gaps.push(Gap {
                before: at,
                residues: vec![u],
            }),
        }
    }
    gaps
}
