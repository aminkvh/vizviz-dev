//! Per-residue conservation of a chain against homologous chains: the
//! chain's sequence and its homologs are aligned pairwise to it, and each
//! of its residues gets the fraction of chains sharing it and a score.
//!
//! The score follows the properties Valdar (2002, Proteins 48:227) argues
//! a conservation score needs: it is a Shannon entropy, normalized to 0..1
//! (1 = one residue type, 0 = all twenty equally); sequences are weighted
//! by Henikoff and Henikoff (1994, J Mol Biol 243:574) so near-identical
//! copies do not dominate; and gaps lower the score in proportion to the
//! weight of the sequences that have one.

use std::collections::HashMap;

use super::align::{align, identical};

/// Pairwise identity, over the shorter sequence, from which two chains
/// count as homologs.
pub const MIN_IDENTITY: f32 = 0.30;
/// Shorter sequences are peptides, not comparable to a domain.
const MIN_LENGTH: usize = 10;
const GAP: u8 = b'-';

/// One residue of the chain and what the other chains have there.
#[derive(Clone, Debug, PartialEq)]
pub struct Column {
    pub residue: u8,
    /// Chains compared, the chain itself included.
    pub sequences: u32,
    /// Chains with the chain's own residue here.
    pub identical: u32,
    /// Residue letters found here and how many chains have each, most
    /// common first; gaps are not listed.
    pub counts: Vec<(u8, u32)>,
    /// Weighted normalized entropy with gap penalty: 1 is fully conserved.
    pub score: f32,
}

impl Column {
    pub fn identity(&self) -> f32 {
        self.identical as f32 / self.sequences as f32
    }
}

pub struct Conservation {
    pub columns: Vec<Column>,
    /// Homologous chains found, not counting the chain itself.
    pub homologs: usize,
}

/// The chain `target` against the sequences in `others`; `None` when none
/// is homologous.
pub fn conservation(target: &[u8], others: &[&[u8]]) -> Option<Conservation> {
    let rows = homolog_rows(target, others);
    if rows.is_empty() || target.is_empty() {
        return None;
    }
    let mut all: Vec<Vec<u8>> = vec![target.to_vec()];
    all.extend(rows);
    let weights = henikoff_weights(&all, target.len());
    let columns = (0..target.len())
        .map(|c| column(&all, &weights, c))
        .collect();
    Some(Conservation {
        columns,
        homologs: all.len() - 1,
    })
}

/// The residues each homolog has at each position of `target`, gaps as
/// `-`. A sequence appearing several times is aligned once.
fn homolog_rows(target: &[u8], others: &[&[u8]]) -> Vec<Vec<u8>> {
    let mut copies: HashMap<&[u8], usize> = HashMap::new();
    for other in others {
        *copies.entry(*other).or_default() += 1;
    }
    let mut rows = Vec::new();
    for (other, n) in copies {
        if let Some(row) = aligned_row(target, other) {
            rows.extend(std::iter::repeat_n(row, n));
        }
    }
    rows.sort();
    rows
}

fn aligned_row(target: &[u8], other: &[u8]) -> Option<Vec<u8>> {
    let shorter = target.len().min(other.len());
    if shorter < MIN_LENGTH {
        return None;
    }
    if target == other {
        return Some(target.to_vec());
    }
    let map = align(target, other);
    let same = identical(target, other, &map);
    if (same as f32) < MIN_IDENTITY * shorter as f32 {
        return None;
    }
    Some(map.iter().map(|m| m.map_or(GAP, |j| other[j])).collect())
}

/// Position-based weights, summing to 1: a sequence weighs more where it is
/// alone with its residue.
fn henikoff_weights(rows: &[Vec<u8>], width: usize) -> Vec<f64> {
    let mut weights = vec![0.0; rows.len()];
    for c in 0..width {
        let mut counts: HashMap<u8, u32> = HashMap::new();
        for row in rows {
            *counts.entry(row[c]).or_default() += 1;
        }
        let kinds = counts.len() as f64;
        for (w, row) in weights.iter_mut().zip(rows) {
            *w += 1.0 / (kinds * counts[&row[c]] as f64);
        }
    }
    let total: f64 = weights.iter().sum();
    weights.iter().map(|w| w / total).collect()
}

fn column(rows: &[Vec<u8>], weights: &[f64], c: usize) -> Column {
    let mut counts: HashMap<u8, u32> = HashMap::new();
    let mut weight_of: HashMap<u8, f64> = HashMap::new();
    for (row, w) in rows.iter().zip(weights) {
        if row[c] != GAP {
            *counts.entry(row[c]).or_default() += 1;
            *weight_of.entry(row[c]).or_default() += w;
        }
    }
    let residue = rows[0][c];
    let mut listed: Vec<(u8, u32)> = counts.into_iter().collect();
    listed.sort_by_key(|&(letter, n)| (std::cmp::Reverse(n), letter));
    Column {
        residue,
        sequences: rows.len() as u32,
        identical: rows.iter().filter(|r| r[c] == residue).count() as u32,
        counts: listed,
        score: entropy_score(&weight_of),
    }
}

/// `(1 - H / log2 20) * (1 - gap weight)`, from the summed sequence weight
/// of each residue type in the column.
fn entropy_score(weight_of: &HashMap<u8, f64>) -> f32 {
    let present: f64 = weight_of.values().sum();
    if present <= 0.0 {
        return 0.0;
    }
    let entropy: f64 = weight_of
        .values()
        .map(|w| w / present)
        .filter(|p| *p > 0.0)
        .map(|p| -p * p.log2())
        .sum();
    let conserved = 1.0 - entropy / 20f64.log2();
    (conserved.max(0.0) * present.min(1.0)) as f32
}

#[cfg(test)]
mod tests {
    use super::*;

    const CHAIN: &[u8] = b"MKVLAAGIVWWTTQQPSNDEHKRY";

    #[test]
    fn identical_chains_are_fully_conserved() {
        let c = conservation(CHAIN, &[CHAIN, CHAIN]).unwrap();
        assert_eq!(c.homologs, 2);
        for col in &c.columns {
            assert_eq!(col.identical, 3);
            assert!((col.score - 1.0).abs() < 1e-6, "{}", col.score);
        }
    }

    #[test]
    fn a_point_mutant_flags_only_its_column() {
        let mut mutant = CHAIN.to_vec();
        mutant[7] = b'D';
        let c = conservation(CHAIN, &[CHAIN, CHAIN, &mutant]).unwrap();
        for (i, col) in c.columns.iter().enumerate() {
            if i == 7 {
                assert_eq!(col.identical, 3);
                assert_eq!(col.sequences, 4);
                assert!(col.score < 0.9, "{}", col.score);
                assert_eq!(col.counts, [(b'I', 3), (b'D', 1)]);
            } else {
                assert!((col.score - 1.0).abs() < 1e-6);
            }
        }
    }

    #[test]
    fn an_unrelated_sequence_is_not_a_homolog() {
        let other = b"GGGGGGGGGGGGGGGGGGGGGGGG";
        assert!(conservation(CHAIN, &[other]).is_none());
    }

    #[test]
    fn gaps_lower_the_score() {
        let short = b"MKVLAAGIVWWTTQQ";
        let c = conservation(CHAIN, &[short]).unwrap();
        assert!((c.columns[0].score - 1.0).abs() < 1e-6);
        assert!(c.columns[20].score < 0.7);
        assert_eq!(c.columns[20].identical, 1);
    }
}
