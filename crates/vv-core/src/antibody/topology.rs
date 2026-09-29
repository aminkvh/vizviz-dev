//! Antibody domains of the protein chains of a [`Topology`].

use std::ops::Range;

use rayon::prelude::*;

use super::{find_domains, CdrDefinition, ChainType, Domain};
use crate::residue_class::ResidueClass;
use crate::seqfeat::one_letter;
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
        ResidueClass::Protein => one_letter(top.residue_name(residue))
            .filter(char::is_ascii_uppercase)
            .unwrap_or('X'),
        _ => 'X',
    }
}

/// Variable domains in the protein residues `residues`; `start` and `end`
/// of each are offsets from `residues.start`.
pub fn find_in_residues(top: &Topology, residues: Range<u32>) -> Vec<Domain> {
    let seq: String = residues.map(|r| protein_letter(top, r)).collect();
    find_domains(&seq)
}

/// Every CDR residue under `definition`, chain by chain (in parallel).
pub fn cdr_residues(top: &Topology, definition: CdrDefinition) -> Vec<CdrResidue> {
    top.chains
        .par_iter()
        .flat_map_iter(|chain| {
            let first = chain.residues.start;
            find_in_residues(top, chain.residues.clone())
                .into_iter()
                .flat_map(move |domain| {
                    let notes = domain.annotate(definition.native_scheme(), definition);
                    let kind = domain.chain;
                    notes.into_iter().filter_map(move |a| {
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
