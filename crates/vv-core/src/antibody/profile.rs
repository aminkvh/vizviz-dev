//! Framework-only position-specific scoring profiles, built once from the
//! seed sequences in `seeds.rs`.

use std::sync::OnceLock;

use super::seeds;
use super::ChainType;

/// One clustered framework example: `fr` is FR1..FR4 as one-letter strings,
/// gap columns already omitted, `-` marking a truncated FR4 tail.
pub(super) struct Seed {
    pub weight: u32,
    pub fr: [&'static str; 4],
}

/// IMGT framework column ranges (inclusive); the CDR columns between them
/// are never scored.
const FR_RANGES: [(u16, u16); 4] = [(1, 26), (39, 55), (66, 104), (118, 128)];

pub(super) const N_COLS: usize = 93;
/// Index of the first column after each CDR loop (FR2, FR3, FR4 starts).
pub(super) const LOOP_EXIT: [usize; 3] = [26, 43, 82];
pub(super) const C23: usize = 22;
pub(super) const C104: usize = 81;

/// Column of the conserved anchors used to sanity-check an alignment.
pub(super) const W41: usize = 26 + 2;

/// IMGT columns a germline leaves empty in FR1 and FR3.
fn empty_columns(chain: ChainType) -> &'static [u16] {
    match chain {
        ChainType::Heavy => &[10, 73],
        ChainType::Kappa => &[73, 81, 82],
        ChainType::Lambda => &[10, 73, 81, 82],
    }
}

pub(super) fn column_number(idx: usize) -> u16 {
    let (mut idx, mut out) = (idx as u16, 0);
    for (lo, hi) in FR_RANGES {
        let len = hi - lo + 1;
        if idx < len {
            out = lo + idx;
            break;
        }
        idx -= len;
    }
    out
}

/// Allowed CDR loop length and soft penalty, in residues.
pub(super) struct LoopPrior {
    pub min: usize,
    pub typical: (usize, usize),
    pub max: usize,
}

const SLOPE: f32 = 1.5;

impl LoopPrior {
    /// Cost per loop length `0..=max`; lengths below `min` are barred.
    fn costs(&self) -> Vec<f32> {
        (0..=self.max)
            .map(|len| {
                let outside =
                    self.typical.0.saturating_sub(len) + len.saturating_sub(self.typical.1);
                if len < self.min {
                    BARRED
                } else {
                    -SLOPE * outside as f32
                }
            })
            .collect()
    }
}

const BARRED: f32 = -1.0e8;

pub(super) struct Profile {
    pub chain: ChainType,
    /// Match log-odds in bits per column; index 20 is any non-standard letter.
    pub score: Vec<[f32; 21]>,
    /// Columns a germline leaves empty: deleting them is free.
    pub free_gap: Vec<bool>,
    /// Score of a CDR loop of each length (index), `NEG`-like when barred.
    pub loop_costs: [Vec<f32>; 3],
    /// Score of the per-column best residue; normalises a domain score.
    pub ideal: f32,
}

const AA: &[u8; 20] = b"ARNDCQEGHILKMFPSTWYV";

pub(super) fn aa_index(b: u8) -> usize {
    AA.iter()
        .position(|&a| a == b.to_ascii_uppercase())
        .unwrap_or(20)
}

/// Robinson & Robinson (1991) background frequencies, in `AA` order.
const BACKGROUND: [f32; 20] = [
    0.078, 0.051, 0.045, 0.054, 0.019, 0.043, 0.063, 0.074, 0.022, 0.051, 0.091, 0.057, 0.022,
    0.039, 0.052, 0.071, 0.058, 0.013, 0.032, 0.064,
];

/// BLOSUM62 (Henikoff & Henikoff 1992), in `AA` order.
#[rustfmt::skip]
const BLOSUM62: [[i8; 20]; 20] = [
    [ 4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0],
    [-1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3],
    [-2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3],
    [-2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3],
    [ 0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1],
    [-1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2],
    [-1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2],
    [ 0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3],
    [-2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3],
    [-1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3],
    [-1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1],
    [-1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2],
    [-1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1],
    [-2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1],
    [-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2],
    [ 1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2],
    [ 0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0],
    [-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3],
    [-2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1],
    [ 0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4],
];

/// Weight of the substitution-matrix pseudocounts against observed counts.
const PSEUDOCOUNT: f32 = 4.0;
/// Score of a residue placed in a column the germline leaves empty.
const EMPTY_COLUMN_SCORE: f32 = -1.0;

/// `target[b][a]`: probability of seeing `a` in a column where `b` is
/// conserved, from BLOSUM62 half-bit scores and the background.
fn substitution_targets() -> [[f32; 20]; 20] {
    let mut t = [[0.0; 20]; 20];
    for b in 0..20 {
        let mut sum = 0.0;
        for a in 0..20 {
            t[b][a] = BACKGROUND[a] * 2f32.powf(f32::from(BLOSUM62[a][b]) / 2.0);
            sum += t[b][a];
        }
        t[b].iter_mut().for_each(|x| *x /= sum);
    }
    t
}

fn column_counts(seeds: &[Seed], chain: ChainType) -> Vec<[f32; 20]> {
    let mut counts = vec![[0.0f32; 20]; N_COLS];
    let empty = empty_columns(chain);
    for seed in seeds {
        let w = 1.0 + (seed.weight as f32).ln();
        let mut idx = 0;
        for (block, &(lo, hi)) in FR_RANGES.iter().enumerate() {
            let mut letters = seed.fr[block].bytes();
            for col in lo..=hi {
                idx += 1;
                if empty.contains(&col) {
                    continue;
                }
                let a = letters.next().map_or(20, aa_index);
                if a < 20 {
                    counts[idx - 1][a] += w;
                }
            }
        }
    }
    counts
}

fn scores_from_counts(counts: &[[f32; 20]], free: &[bool]) -> Vec<[f32; 21]> {
    let target = substitution_targets();
    counts
        .iter()
        .zip(free)
        .map(|(c, &is_free)| {
            let n: f32 = c.iter().sum();
            let mut row = [0.0f32; 21];
            for a in 0..20 {
                let pseudo: f32 = (0..20).map(|b| c[b] / n.max(1e-6) * target[b][a]).sum();
                let p = (c[a] + PSEUDOCOUNT * pseudo) / (n + PSEUDOCOUNT);
                row[a] = if is_free {
                    EMPTY_COLUMN_SCORE
                } else {
                    (p / BACKGROUND[a]).log2()
                };
            }
            row
        })
        .collect()
}

fn loop_priors(chain: ChainType) -> [LoopPrior; 3] {
    let prior = |min, lo, hi, max| LoopPrior {
        min,
        typical: (lo, hi),
        max,
    };
    match chain {
        ChainType::Heavy => [
            prior(0, 7, 10, 22),
            prior(0, 5, 10, 30),
            prior(0, 5, 22, 50),
        ],
        ChainType::Kappa => [prior(0, 5, 12, 20), prior(0, 3, 3, 10), prior(0, 8, 11, 20)],
        ChainType::Lambda => [prior(0, 6, 11, 20), prior(0, 3, 7, 12), prior(0, 8, 13, 20)],
    }
}

fn build(chain: ChainType, seeds: &[Seed]) -> Profile {
    let empty = empty_columns(chain);
    let free_gap: Vec<bool> = (0..N_COLS)
        .map(|i| empty.contains(&column_number(i)))
        .collect();
    let score = scores_from_counts(&column_counts(seeds, chain), &free_gap);
    let ideal = score
        .iter()
        .zip(&free_gap)
        .filter(|(_, &free)| !free)
        .map(|(row, _)| row[..20].iter().copied().fold(f32::MIN, f32::max))
        .sum();
    let loop_costs = loop_priors(chain).map(|p| p.costs());
    Profile {
        chain,
        score,
        free_gap,
        loop_costs,
        ideal,
    }
}

pub(super) fn profiles() -> &'static [Profile; 3] {
    static PROFILES: OnceLock<[Profile; 3]> = OnceLock::new();
    PROFILES.get_or_init(|| {
        [
            build(ChainType::Heavy, seeds::HEAVY),
            build(ChainType::Kappa, seeds::KAPPA),
            build(ChainType::Lambda, seeds::LAMBDA),
        ]
    })
}
