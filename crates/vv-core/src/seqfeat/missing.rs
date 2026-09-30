//! Unobserved residues: those a chain's full sequence has but the model
//! lacks, placed against the residues that are modeled.
//!
//! With the chain's full sequence ([`EntityChain`]) the model is aligned to
//! it, so insertion codes and numbering that is not ascending do not
//! matter. Without one, the file's list of unobserved residues (`REMARK
//! 465`, mmCIF `_pdbx_unobs_or_zero_occ_residues`) is placed by comparing
//! author numbers, which assumes they ascend.

use std::collections::HashSet;
use std::ops::Range;

use super::align::align;
use super::letters::one_letter;
use crate::{ResidueClass, Topology};

/// One unobserved residue: its name, its 1-based position in the chain's
/// full sequence (0 when only a numbered list is known) and, when known,
/// its author number and insertion code.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Unobserved {
    pub name: String,
    pub position: u32,
    pub number: Option<(i32, u8)>,
}

/// Unobserved residues that sit just before modeled residue `before`
/// (a topology residue index), or after the last one when `before` is the
/// end of the chain range.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Gap {
    pub before: u32,
    pub residues: Vec<Unobserved>,
}

/// A chain's full sequence as its file states it, one entry per position.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct EntityChain {
    /// Residue names (`ALA`), in sequence order.
    pub names: Vec<String>,
    /// Author number and insertion code of each position, where the file
    /// gives them.
    pub numbers: Vec<Option<(i32, u8)>>,
}

fn is_polymer_chain(top: &Topology, residues: &Range<u32>) -> bool {
    residues.clone().next().is_some_and(|first| {
        matches!(
            top.residue_class(first as usize),
            ResidueClass::Protein | ResidueClass::Nucleic
        )
    })
}

fn author_chain(top: &Topology, first: u32) -> &str {
    top.names
        .get(top.chains[top.residues[first as usize].chain as usize].auth_asym)
}

/// The polymer residues of `residues` (a chain's rows of one name) that
/// the file lists as unobserved, grouped into gaps between modeled ones.
/// Empty when the file carries no such table. Matches on the author chain
/// id and author numbering, the same numbering the strip shows.
pub fn unobserved(top: &Topology, residues: Range<u32>) -> Vec<Gap> {
    if !is_polymer_chain(top, &residues) {
        return Vec::new();
    }
    let listed = listed_for_chain(top, author_chain(top, residues.start), 1);
    place(top, residues, listed)
}

/// `model` is the 1-based model whose table rows count.
fn listed_for_chain(top: &Topology, chain: &str, model: u32) -> Vec<Unobserved> {
    let Some(cat) = top.annotations.category("pdbx_unobs_or_zero_occ_residues") else {
        return Vec::new();
    };
    let model = model.to_string();
    // `occupancy_flag` is 1 for residues with no coordinates and 0 for
    // zero-occupancy ones, which are modeled. The table read from `REMARK
    // 465` has no model column and lists only residues with none.
    let from_remark = cat.column("PDB_model_num").is_none();
    let mut seen: HashSet<(i32, u8)> = HashSet::new();
    (0..cat.rows.len())
        .filter(|&r| cat.get("auth_asym_id", r) == Some(chain))
        .filter(|&r| cat.get("polymer_flag", r).is_none_or(|f| f == "Y"))
        .filter(|&r| from_remark || cat.get("occupancy_flag", r) == Some("1"))
        .filter(|&r| cat.get("PDB_model_num", r).is_none_or(|m| m == model))
        .filter_map(|r| {
            let seq = cat.get("auth_seq_id", r)?.parse().ok()?;
            let ins = cat.get("PDB_ins_code", r).map_or(0, |s| s.as_bytes()[0]);
            seen.insert((seq, ins)).then(|| Unobserved {
                name: cat.get("auth_comp_id", r).unwrap_or("?").to_string(),
                position: 0,
                number: Some((seq, ins)),
            })
        })
        .collect()
}

/// Walks the modeled residues in order, assigning each listed residue to
/// the first modeled one that follows it in `(seq, ins)` order.
fn place(top: &Topology, residues: Range<u32>, mut listed: Vec<Unobserved>) -> Vec<Gap> {
    listed.sort_by_key(|u| u.number);
    let key = |r: u32| {
        let rec = &top.residues[r as usize];
        Some((rec.auth_seq_id, rec.ins_code))
    };
    let mut at = residues.start;
    let mut gaps: Vec<Gap> = Vec::new();
    for u in listed {
        while at < residues.end && key(at) < u.number {
            at += 1;
        }
        push_into(&mut gaps, at, u);
    }
    gaps
}

fn push_into(gaps: &mut Vec<Gap>, before: u32, u: Unobserved) {
    match gaps.last_mut() {
        Some(g) if g.before == before => g.residues.push(u),
        _ => gaps.push(Gap {
            before,
            residues: vec![u],
        }),
    }
}

/// Like [`unobserved`], for a chain whose full sequence is known: the
/// residues of `entity` the modeled residues do not align to. `model` is
/// the shown 1-based model; when the file lists unobserved residues for
/// several models, only positions listed for that model count.
pub fn unobserved_in_entity(
    top: &Topology,
    residues: Range<u32>,
    entity: &EntityChain,
    model: u32,
) -> Vec<Gap> {
    if !is_polymer_chain(top, &residues) {
        return Vec::new();
    }
    let modeled: Vec<u8> = residues
        .clone()
        .map(|r| letter(top.residue_name(r as usize)))
        .collect();
    let full: Vec<u8> = entity.names.iter().map(|n| letter(n)).collect();
    let placed = align(&full, &modeled);
    let listed = model_listing(top, author_chain(top, residues.start), model);
    let numbers = numbered_positions(entity, &placed, top, &residues);
    let mut gaps: Vec<Gap> = Vec::new();
    let mut run: Vec<Unobserved> = Vec::new();
    for (i, at) in placed.iter().enumerate() {
        match at {
            None => {
                let number = numbers[i];
                let kept = listed
                    .as_ref()
                    .is_none_or(|l| number.is_some_and(|n| l.contains(&n)));
                if kept {
                    run.push(Unobserved {
                        name: entity.names[i].clone(),
                        position: i as u32 + 1,
                        number,
                    });
                }
            }
            Some(k) => {
                let before = residues.start + *k as u32;
                for u in run.drain(..) {
                    push_into(&mut gaps, before, u);
                }
            }
        }
    }
    for u in run {
        push_into(&mut gaps, residues.end, u);
    }
    gaps
}

fn letter(name: &str) -> u8 {
    one_letter(name).map_or(b'X', |c| c as u8)
}

/// Author numbers the shown model lists as unobserved, when the table has
/// rows for more than one model; `None` when there is nothing to filter by.
fn model_listing(top: &Topology, chain: &str, model: u32) -> Option<HashSet<(i32, u8)>> {
    let cat = top
        .annotations
        .category("pdbx_unobs_or_zero_occ_residues")?;
    let models: HashSet<&str> = (0..cat.rows.len())
        .filter_map(|r| cat.get("PDB_model_num", r))
        .collect();
    if models.len() < 2 {
        return None;
    }
    let shown = if models.contains(model.to_string().as_str()) {
        model
    } else {
        1
    };
    let listed = listed_for_chain(top, chain, shown);
    Some(listed.into_iter().filter_map(|u| u.number).collect())
}

/// The author number of each entity position: as the entity states it, or
/// else taken in order from the file's unobserved list, which is written in
/// sequence order, when it names exactly the positions left unaligned.
fn numbered_positions(
    entity: &EntityChain,
    placed: &[Option<usize>],
    top: &Topology,
    residues: &Range<u32>,
) -> Vec<Option<(i32, u8)>> {
    if entity.numbers.len() == entity.names.len() {
        return entity.numbers.clone();
    }
    let unaligned: Vec<usize> = (0..placed.len()).filter(|&i| placed[i].is_none()).collect();
    let listed = listed_for_chain(top, author_chain(top, residues.start), 1);
    let mut numbers = vec![None; placed.len()];
    if listed.len() == unaligned.len() {
        for (&i, u) in unaligned.iter().zip(&listed) {
            numbers[i] = u.number;
        }
    }
    numbers
}
