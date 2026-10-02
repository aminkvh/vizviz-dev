//! Numbering from a chain's full deposited sequence.
//!
//! Residues missing from the model still count: they take their numbers
//! from the full sequence and leave holes, and residues before the domain
//! (tags, signal peptides) get none.

use super::{find_domains, Domain};
use crate::seqfeat::align::align;

/// Variable domains of the chain whose deposited sequence is `full` and
/// whose modelled residues are `observed` (a subsequence of `full`, with
/// `X` for anything the one-letter code lacks). Indices in the result
/// refer to `observed`; `start..end` spans the modelled residues of each
/// domain, and a domain with none is dropped.
///
/// Falls back to [`find_domains`] on `observed` when `observed` cannot be
/// laid onto `full`.
pub fn find_domains_in_chain(full: &str, observed: &str) -> Vec<Domain> {
    let Some(placed) = place(full, observed) else {
        return find_domains(observed);
    };
    find_domains(full)
        .into_iter()
        .filter_map(|d| restrict(d, &placed))
        .collect()
}

/// For each position of `full`, the `observed` residue standing there.
/// `None` when fewer than half of the observed residues find a place.
fn place(full: &str, observed: &str) -> Option<Vec<Option<usize>>> {
    let placed = align(full.as_bytes(), observed.as_bytes());
    let found = placed.iter().flatten().count();
    (found * 2 >= observed.len() && found > 0).then_some(placed)
}

fn restrict(mut d: Domain, placed: &[Option<usize>]) -> Option<Domain> {
    d.at = placed[d.start..d.end].to_vec();
    let mut seen = d.at.iter().flatten().copied();
    let first = seen.next()?;
    let last = seen.last().unwrap_or(first);
    d.start = first;
    d.end = last + 1;
    Some(d)
}

impl Domain {
    /// The domain with every sequence index `i` replaced by `to(i)`, which
    /// must rise with `i`.
    pub(crate) fn reindexed(mut self, to: impl Fn(usize) -> usize) -> Domain {
        self.start = to(self.start);
        self.end = to(self.end - 1) + 1;
        for at in self.at.iter_mut().flatten() {
            *at = to(*at);
        }
        self
    }
}
