//! Framework-only position-specific scoring profiles, built once from the
//! seed sequences in `seeds.rs` (antibodies) and `tcr_seeds.rs` (T-cell
//! receptors).

use std::collections::BTreeMap;
use std::sync::OnceLock;

use super::numbering::Label;
use super::ChainType;
use super::{seeds, tcr_seeds};

/// One clustered framework example: `fr` is FR1..FR4 as one-letter strings
/// with the family's empty columns omitted, `-` marking a column with no
/// residue (a truncated end, or disordered).
pub(super) struct Seed {
    pub weight: u32,
    pub fr: [&'static str; 4],
}

/// IMGT framework column ranges (inclusive); the CDR columns between them
/// are never scored.
const FR_RANGES: [(u16, u16); 4] = [(1, 26), (39, 55), (66, 104), (118, 128)];

/// Columns before the first CDR: the alignment may start anywhere up to the
/// first conserved Cys (IMGT 23).
pub(super) const C23: usize = 22;

/// Column of the conserved Trp (IMGT 41), used to sanity-check an alignment.
pub(super) const W41: usize = 26 + 2;

/// Chains of one type that leave the same framework columns empty share a
/// profile, the way the germline families of an antibody type do.
struct Family {
    chain: ChainType,
    /// IMGT columns the family leaves empty.
    empty: &'static [u16],
    /// Insertion columns the family fills (`(84, 'A')` is IMGT 84A).
    extras: &'static [(u16, char)],
    /// Insertion columns no seed fills: a residue there costs little, so a
    /// long framework does not spill into the neighbouring CDR.
    spare: &'static [(u16, char)],
    seeds: &'static [Seed],
}

const ALPHA_EXTRAS: &[(u16, char)] = &[(84, 'A'), (84, 'B'), (84, 'C')];
const LIGHT_FR3_SPARE: &[(u16, char)] = &[(82, 'A'), (82, 'B'), (82, 'C')];
const HEAVY_FR3_SPARE: &[(u16, char)] = &[(94, 'A'), (94, 'B'), (94, 'C')];

fn families() -> [Family; 7] {
    let family = |chain, empty, extras, seeds| Family {
        spare: &[],
        chain,
        empty,
        extras,
        seeds,
    };
    [
        Family {
            spare: HEAVY_FR3_SPARE,
            ..family(ChainType::Heavy, &[10, 73], &[], seeds::HEAVY)
        },
        Family {
            spare: LIGHT_FR3_SPARE,
            ..family(ChainType::Kappa, &[73, 81, 82], &[], seeds::KAPPA)
        },
        Family {
            spare: LIGHT_FR3_SPARE,
            ..family(ChainType::Lambda, &[10, 73, 81, 82], &[], seeds::LAMBDA)
        },
        family(
            ChainType::TcrAlpha,
            &[69, 70, 71, 72, 73],
            &[],
            tcr_seeds::ALPHA_69_73,
        ),
        family(
            ChainType::TcrAlpha,
            &[71, 72, 73, 74, 75, 76, 77],
            ALPHA_EXTRAS,
            tcr_seeds::ALPHA_71_77,
        ),
        family(ChainType::TcrBeta, &[73, 82], &[], tcr_seeds::BETA_73_82),
        family(ChainType::TcrBeta, &[82], &[], tcr_seeds::BETA_82),
    ]
}

/// IMGT label of every profile column, in column order.
fn column_labels(family: &Family) -> Vec<Label> {
    let mut out = Vec::new();
    for (lo, hi) in FR_RANGES {
        for number in lo..=hi {
            out.push(Label::new(number));
            for &(after, letter) in family
                .extras
                .iter()
                .chain(family.spare)
                .filter(|e| e.0 == number)
            {
                out.push(Label::with_insertion(after, letter));
            }
        }
    }
    out
}

/// Which of the four framework blocks a column number belongs to.
fn block_of(number: u16) -> usize {
    FR_RANGES
        .iter()
        .position(|&(_, hi)| number <= hi)
        .unwrap_or(3)
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
    /// IMGT label of each column.
    pub labels: Vec<Label>,
    /// Index of the first column after each CDR loop (FR2, FR3, FR4 starts).
    pub loop_exit: [usize; 3],
    /// Column of the second conserved Cys (IMGT 104).
    pub c104: usize,
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
/// Score of a residue placed in a column the germline leaves empty. Receptor
/// families are told apart by which columns they leave empty, so a residue
/// there costs more.
const EMPTY_COLUMN_SCORE: f32 = -1.0;
const RECEPTOR_EMPTY_SCORE: f32 = -4.0;

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

/// Residue counts per column label over one family's seeds; a seed counts
/// less than its multiplicity.
fn family_counts(family: &Family) -> BTreeMap<Label, [f32; 20]> {
    let labels = column_labels(family);
    let mut counts: BTreeMap<Label, [f32; 20]> = BTreeMap::new();
    for seed in family.seeds {
        let w = 1.0 + (seed.weight as f32).ln();
        let mut letters: Vec<_> = seed.fr.iter().map(|s| s.bytes()).collect();
        for label in &labels {
            if is_empty(family, label) {
                continue;
            }
            let a = letters[block_of(label.number)].next().map_or(20, aa_index);
            if a < 20 {
                counts.entry(*label).or_default()[a] += w;
            }
        }
    }
    counts
}

/// Counts pooled over every family of one chain type: families of a type
/// differ in which columns they fill, not in what the shared ones hold.
fn pooled_counts(chain: ChainType, families: &[Family]) -> BTreeMap<Label, [f32; 20]> {
    let mut pool: BTreeMap<Label, [f32; 20]> = BTreeMap::new();
    for family in families.iter().filter(|f| f.chain == chain) {
        for (label, counts) in family_counts(family) {
            let slot = pool.entry(label).or_default();
            slot.iter_mut().zip(counts).for_each(|(p, c)| *p += c);
        }
    }
    pool
}

fn is_empty(family: &Family, label: &Label) -> bool {
    match label.insertion() {
        None => family.empty.contains(&label.number),
        Some(letter) => family.spare.contains(&(label.number, letter)),
    }
}

fn scores_from_counts(counts: &[[f32; 20]], free: &[bool], empty_score: f32) -> Vec<[f32; 21]> {
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
                    empty_score
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
        ChainType::TcrAlpha => [prior(0, 5, 7, 14), prior(0, 4, 7, 14), prior(0, 8, 16, 30)],
        ChainType::TcrBeta => [prior(0, 5, 6, 12), prior(0, 6, 7, 12), prior(0, 8, 16, 30)],
    }
}

fn build(family: &Family, families: &[Family]) -> Profile {
    let chain = family.chain;
    let labels = column_labels(family);
    let free_gap: Vec<bool> = labels.iter().map(|l| is_empty(family, l)).collect();
    let pool = pooled_counts(chain, families);
    let counts: Vec<[f32; 20]> = labels
        .iter()
        .map(|l| pool.get(l).copied().unwrap_or([0.0; 20]))
        .collect();
    let empty_score = match chain.is_antibody() {
        true => EMPTY_COLUMN_SCORE,
        false => RECEPTOR_EMPTY_SCORE,
    };
    let score = scores_from_counts(&counts, &free_gap, empty_score);
    let column_of = |number: u16| {
        labels
            .iter()
            .position(|l| l.number == number && l.insertion().is_none())
            .unwrap_or(0)
    };
    let loop_exit = [39, 66, 118].map(column_of);
    let c104 = column_of(104);
    let ideal = score
        .iter()
        .zip(&free_gap)
        .filter(|(_, &free)| !free)
        .map(|(row, _)| row[..20].iter().copied().fold(f32::MIN, f32::max))
        .sum();
    let loop_costs = loop_priors(chain).map(|p| p.costs());
    Profile {
        chain,
        labels,
        loop_exit,
        c104,
        score,
        free_gap,
        loop_costs,
        ideal,
    }
}

pub(super) fn profiles() -> &'static [Profile] {
    static PROFILES: OnceLock<Vec<Profile>> = OnceLock::new();
    PROFILES.get_or_init(|| {
        let all = families();
        all.iter().map(|f| build(f, &all)).collect()
    })
}
