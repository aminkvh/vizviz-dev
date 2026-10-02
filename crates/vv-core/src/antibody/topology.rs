//! Antibody domains of the protein chains of a [`Topology`].

use std::fmt;
use std::ops::Range;
use std::sync::OnceLock;

use rayon::prelude::*;

use super::{find_domains, find_domains_in_chain, CdrDefinition, ChainType, Domain};
use crate::residue_class::ResidueClass;
use crate::seqfeat::protein_letter as name_letter;
use crate::Topology;

/// A residue inside a CDR.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CdrResidue {
    pub residue: u32,
    pub chain: ChainType,
    /// 1, 2 or 3.
    pub cdr: u8,
}

fn protein_letter(top: &Topology, residue: u32) -> char {
    let residue = residue as usize;
    match top.residue_class(residue) {
        ResidueClass::Protein => name_letter(top.residue_name(residue)),
        _ => 'X',
    }
}

/// Variable domains in the protein residues `residues`; `start` and `end`
/// of each are offsets from `residues.start`. Numbered from the chain's
/// deposited sequence when the file states one.
pub fn find_in_residues(top: &Topology, residues: Range<u32>) -> Vec<Domain> {
    match deposited(top, &residues) {
        Some(full) => find_against(top, residues, full),
        None => {
            let seq: String = residues.map(|r| protein_letter(top, r)).collect();
            find_domains(&seq)
        }
    }
}

/// The deposited sequence of the chain record holding `residues`.
fn deposited<'a>(top: &'a Topology, residues: &Range<u32>) -> Option<&'a str> {
    let first = top.residues.get(residues.start as usize)?;
    let full = top.full_sequence.get(first.chain as usize)?;
    (!full.is_empty()).then_some(full.as_str())
}

fn find_against(top: &Topology, residues: Range<u32>, full: &str) -> Vec<Domain> {
    let rows: Vec<u32> = residues
        .clone()
        .filter(|&r| top.residue_class(r as usize) == ResidueClass::Protein)
        .collect();
    let observed: String = rows.iter().map(|&r| protein_letter(top, r)).collect();
    find_domains_in_chain(full, &observed)
        .into_iter()
        .map(|d| d.reindexed(|i| (rows[i] - residues.start) as usize))
        .collect()
}

/// Variable domains of every chain of a [`Topology`], found on first use.
///
/// A clone starts empty, so a topology edited after cloning never sees a
/// stale answer; edit `chains` or `residues` in place only before the
/// first query, or call [`AntibodyCache::clear`].
#[derive(Default)]
pub struct AntibodyCache(OnceLock<Vec<Vec<Domain>>>);

impl AntibodyCache {
    pub fn clear(&mut self) {
        self.0 = OnceLock::new();
    }
}

impl Clone for AntibodyCache {
    fn clone(&self) -> Self {
        Self::default()
    }
}

impl fmt::Debug for AntibodyCache {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self.0.get() {
            Some(_) => "AntibodyCache(filled)",
            None => "AntibodyCache(empty)",
        })
    }
}

/// Domains of each chain, by chain index, computed once per topology (in
/// parallel) and shared by every later query.
pub fn chain_domains(top: &Topology) -> &[Vec<Domain>] {
    top.antibody.0.get_or_init(|| {
        top.chains
            .par_iter()
            .map(|chain| find_in_residues(top, chain.residues.clone()))
            .collect()
    })
}

/// Every CDR residue under `definition`, chain by chain.
pub fn cdr_residues(top: &Topology, definition: CdrDefinition) -> Vec<CdrResidue> {
    let scheme = definition.native_scheme();
    top.chains
        .iter()
        .zip(chain_domains(top))
        .flat_map(|(chain, domains)| {
            let first = chain.residues.start;
            domains.iter().flat_map(move |domain| {
                let kind = domain.chain;
                domain
                    .annotate(scheme, definition)
                    .into_iter()
                    .filter_map(move |a| {
                        Some(CdrResidue {
                            residue: first + a.index as u32,
                            chain: kind,
                            cdr: a.region.cdr()?,
                        })
                    })
            })
        })
        .collect()
}
