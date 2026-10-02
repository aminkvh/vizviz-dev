//! Places heavy-chain framework indels by similarity to a consensus.
//!
//! The label rules of `numbering` give one labeling per domain. Where the
//! reference numbering (the scheme authors' program) decides between
//! near-equivalent placements by the sequence, candidates around the rule's
//! labeling compete on the log-odds of each residue at its label under a
//! per-position consensus of that numbering (`consensus.rs`), minus a price
//! for each label or letter that differs from the rule. Two places are open:
//!
//! - FR2 (labels 36-49): the block of skipped labels may slide.
//! - The CDR-H2 tail to FR3 (labels 50-92): labels 50-65 may be skipped,
//!   and letters may move between the insertion sites of CDR-H2 and FR3.
//!
//! Every other label and letter stays where the rules put it: the reference
//! keeps the rules' site there whatever the sequence says, so a consensus
//! only adds errors (docs/ANTIBODY.md, "Placement by consensus").

use std::collections::BTreeMap;
use std::ops::{Range, RangeInclusive};
use std::sync::OnceLock;

use super::consensus::{Row, HEAVY};
use super::numbering::{insertion_sites, Label, Scheme};
use super::profile::BACKGROUND;
use super::ChainType;

/// Weight, in residues, of the background against observed counts.
const PSEUDOCOUNT: f32 = 5.0;
/// Bits per FR2 label used or skipped differently from the rule.
const FR2_PRICE: f32 = 4.0;
/// Bits per CDR-H2 tail label (50-65) used or skipped differently.
const TAIL_PRICE: f32 = 30.0;
/// Bits per insertion letter more or fewer than the rule's.
const LETTER_PRICE: f32 = 8.0;
/// Price of any change the boundary window does not open.
const SHUT: f32 = 1.0e3;
/// Letters beyond the rule's own count that a site may take.
const EXTRA_LETTERS: usize = 6;
/// A placement must beat the rule's by this much, so that float noise never
/// swaps equal labelings.
const TIE_MARGIN: f32 = 1.0e-3;
const NEG: f32 = -1.0e9;

const FR2: (u16, u16) = (36, 49);
const BOUNDARY: (u16, u16) = (50, 92);
const TAIL: RangeInclusive<u16> = 50..=65;

/// Log-odds of each residue per label.
struct Consensus {
    base: Vec<[f32; 20]>,
    inserted: BTreeMap<Label, [f32; 20]>,
}

impl Consensus {
    fn from_counts(rows: &[Row]) -> Self {
        let mut base = Vec::new();
        let mut inserted = BTreeMap::new();
        for &(number, insertion, counts) in rows {
            let scores = log_odds(&counts);
            if insertion == 0 {
                let at = usize::from(number);
                if base.len() <= at {
                    base.resize(at + 1, [0.0; 20]);
                }
                base[at] = scores;
            } else {
                inserted.insert(Label::with_insertion(number, insertion as char), scores);
            }
        }
        Consensus { base, inserted }
    }

    /// Zero for a label or residue the consensus has no data for.
    fn score(&self, label: Label, aa: u8) -> f32 {
        let row = match label.insertion() {
            None => self.base.get(usize::from(label.number)),
            Some(_) => self.inserted.get(&label),
        };
        row.and_then(|r| r.get(usize::from(aa)))
            .copied()
            .unwrap_or(0.0)
    }

    fn fit(&self, labels: &[Label], residues: &[u8]) -> f32 {
        labels
            .iter()
            .zip(residues)
            .map(|(l, &aa)| self.score(*l, aa))
            .sum()
    }
}

fn log_odds(counts: &[u16; 20]) -> [f32; 20] {
    let total: f32 = counts.iter().map(|&c| f32::from(c)).sum();
    std::array::from_fn(|a| {
        let p = (f32::from(counts[a]) + PSEUDOCOUNT * BACKGROUND[a]) / (total + PSEUDOCOUNT);
        (p / BACKGROUND[a]).log2()
    })
}

fn consensus(scheme: Scheme) -> Option<&'static Consensus> {
    static TABLES: OnceLock<Vec<Consensus>> = OnceLock::new();
    let tables = TABLES.get_or_init(|| {
        HEAVY
            .iter()
            .map(|rows| Consensus::from_counts(rows))
            .collect()
    });
    let at = match scheme {
        Scheme::Kabat => 0,
        Scheme::Chothia => 1,
        Scheme::Martin => 2,
        _ => return None,
    };
    Some(&tables[at])
}

/// Residue index range of the stretch whose rule labels lie in `lo..=hi`.
fn stretch(cur: &[Label], (lo, hi): (u16, u16)) -> Option<Range<usize>> {
    let inside = |l: &Label| (lo..=hi).contains(&l.number);
    let start = cur.iter().position(inside)?;
    let len = cur[start..].iter().take_while(|l| inside(l)).count();
    Some(start..start + len)
}

/// Slides the block of labels the rule skips in FR2 to where the residues
/// fit the consensus best; one block, so that the residues keep their order
/// around it.
fn slide_fr2(cur: &mut [Label], residues: &[u8], consensus: &Consensus) {
    let Some(range) =
        stretch(cur, FR2).filter(|r| cur[r.clone()].iter().all(|l| l.insertion().is_none()))
    else {
        return;
    };
    let slots = usize::from(FR2.1 - FR2.0) + 1;
    let Some(skipped) = slots.checked_sub(range.len()).filter(|&k| k > 0) else {
        return;
    };
    let rule: Vec<Label> = cur[range.clone()].to_vec();
    let rule_used = |s: usize| rule.iter().any(|l| l.number == FR2.0 + s as u16);
    let window = &residues[range.clone()];
    let mut best = (consensus.fit(&rule, window), rule.clone());
    for at in 0..=slots - skipped {
        let gone = at..at + skipped;
        let labels: Vec<Label> = (0..slots)
            .filter(|s| !gone.contains(s))
            .map(|s| Label::new(FR2.0 + s as u16))
            .collect();
        let moved = (0..slots)
            .filter(|s| gone.contains(s) == rule_used(*s))
            .count();
        let total = consensus.fit(&labels, window) - FR2_PRICE * moved as f32;
        if total > best.0 + TIE_MARGIN {
            best = (total, labels);
        }
    }
    cur[range].copy_from_slice(&best.1);
}

/// The boundary window's placement problem: residues of one stretch
/// against the base labels `lo..=hi` and their letters, relative to the
/// rule's labeling.
struct Problem<'a> {
    consensus: &'a Consensus,
    residues: &'a [u8],
    /// Base label numbers, ascending; the index is the track position.
    track: Vec<u16>,
    /// Whether the rule's labeling uses each base label.
    used: Vec<bool>,
    /// Letters the rule puts after each base label.
    letters: Vec<usize>,
    /// Letters each base label may take (0 where the scheme has no site).
    room: Vec<usize>,
}

impl<'a> Problem<'a> {
    fn new(
        consensus: &'a Consensus,
        residues: &'a [u8],
        sites: &[u16],
        (lo, hi): (u16, u16),
        rule: &[Label],
    ) -> Self {
        let track: Vec<u16> = (lo..=hi).collect();
        let mut used = vec![false; track.len()];
        let mut letters = vec![0usize; track.len()];
        for l in rule {
            let at = usize::from(l.number - lo);
            match l.insertion() {
                None => used[at] = true,
                Some(_) => letters[at] += 1,
            }
        }
        let room = (0..track.len())
            .map(|j| match sites.contains(&track[j]) || letters[j] > 0 {
                true => letters[j] + EXTRA_LETTERS,
                false => 0,
            })
            .collect();
        Problem {
            consensus,
            residues,
            track,
            used,
            letters,
            room,
        }
    }

    fn gap_price(&self, j: usize) -> f32 {
        match TAIL.contains(&self.track[j]) {
            true => TAIL_PRICE,
            false => SHUT,
        }
    }

    /// Price of skipping label `j` when the rule uses it.
    fn gap(&self, j: usize) -> f32 {
        if self.used[j] {
            -self.gap_price(j)
        } else {
            0.0
        }
    }

    /// Score of residue `r` at base label `j`, less the price of using a
    /// label the rule skips.
    fn fill(&self, j: usize, r: usize) -> f32 {
        let score = self
            .consensus
            .score(Label::new(self.track[j]), self.residues[r]);
        score - if self.used[j] { 0.0 } else { self.gap_price(j) }
    }

    /// Score of residue `r` as the `k`-th letter after `j`; letters past the
    /// rule's count cost, letters the rule has already cost are refunded.
    fn letter(&self, j: usize, k: usize, r: usize) -> f32 {
        let label = Label::nth_insertion(self.track[j], k - 1);
        let price = if k > self.letters[j] {
            -LETTER_PRICE
        } else {
            LETTER_PRICE
        };
        self.consensus.score(label, self.residues[r]) + price
    }

    /// Objective of a labeling, `None` when it is not a candidate.
    fn score(&self, labels: &[Label]) -> Option<f32> {
        if labels.len() != self.residues.len() {
            return None;
        }
        let lo = self.track[0];
        let mut used = vec![false; self.track.len()];
        let mut letters = vec![0usize; self.track.len()];
        let mut prev: Option<Label> = None;
        for &l in labels {
            let j = usize::from(l.number.checked_sub(lo)?);
            if j >= used.len() || prev.is_some_and(|p| p >= l) {
                return None;
            }
            prev = Some(l);
            match l.insertion() {
                None => used[j] = true,
                Some(_) => letters[j] += 1,
            }
            let spelled = Label::nth_insertion(l.number, letters[j].saturating_sub(1));
            if letters[j] > self.room[j] || (l.insertion().is_some() && l != spelled) {
                return None;
            }
        }
        let moved: f32 = (0..used.len())
            .filter(|&j| used[j] != self.used[j])
            .map(|j| self.gap_price(j))
            .sum();
        let lettered: usize = (0..used.len())
            .map(|j| letters[j].abs_diff(self.letters[j]))
            .sum();
        let fit = self.consensus.fit(labels, self.residues);
        Some(fit - moved - LETTER_PRICE * lettered as f32)
    }
}

/// Best scores after residue `r` is placed.
struct Step {
    /// Matched to base label `j`.
    at: Vec<f32>,
    /// Labels up to `j` consumed, `j` itself matched or skipped.
    closed: Vec<f32>,
    /// Labels up to `j` consumed and any letters after `j`.
    open: Vec<f32>,
    /// `lettered[j][k-1]`: `k` letters after `j`, this residue the last.
    lettered: Vec<Vec<f32>>,
}

fn best_index(values: &[f32]) -> Option<usize> {
    let mut best: Option<usize> = None;
    for (i, v) in values.iter().enumerate() {
        if best.is_none_or(|b| *v > values[b]) {
            best = Some(i);
        }
    }
    best
}

impl Problem<'_> {
    /// Before any residue: labels may be skipped from the start.
    fn start(&self) -> Step {
        let n = self.track.len();
        let mut closed = vec![NEG; n];
        let mut carried = 0.0;
        for (j, c) in closed.iter_mut().enumerate() {
            carried += self.gap(j);
            *c = carried;
        }
        Step {
            at: vec![NEG; n],
            open: closed.clone(),
            closed,
            lettered: self.room.iter().map(|&k| vec![NEG; k]).collect(),
        }
    }

    fn next(&self, prev: &Step, r: usize) -> Step {
        let n = self.track.len();
        let at: Vec<f32> = (0..n)
            .map(|j| {
                let before = match j {
                    0 if r == 0 => 0.0,
                    0 => NEG,
                    _ => prev.open[j - 1],
                };
                before + self.fill(j, r)
            })
            .collect();
        let lettered: Vec<Vec<f32>> = (0..n)
            .map(|j| {
                (1..=self.room[j])
                    .map(|k| {
                        let before = match k {
                            1 => prev.closed[j],
                            _ => prev.lettered[j][k - 2],
                        };
                        before + self.letter(j, k, r)
                    })
                    .collect()
            })
            .collect();
        let mut closed = vec![NEG; n];
        let mut open = vec![NEG; n];
        for j in 0..n {
            let skipped = if j == 0 {
                NEG
            } else {
                open[j - 1] + self.gap(j)
            };
            closed[j] = at[j].max(skipped);
            open[j] = lettered[j].iter().copied().fold(closed[j], f32::max);
        }
        Step {
            at,
            closed,
            open,
            lettered,
        }
    }

    /// Best labeling and its objective; `None` when no residues or no path.
    fn solve(&self) -> Option<(Vec<Label>, f32)> {
        let m = self.residues.len();
        if m == 0 {
            return None;
        }
        let mut steps = Vec::with_capacity(m + 1);
        steps.push(self.start());
        for r in 0..m {
            let step = self.next(&steps[r], r);
            steps.push(step);
        }
        // Sites the path never enters owe the price of the rule's letters
        // there, which the letter steps only charge where a path goes.
        let owed = LETTER_PRICE * self.letters.iter().sum::<usize>() as f32;
        let best = steps[m].open[self.track.len() - 1];
        (best > NEG / 2.0).then(|| (self.trace(&steps), best - owed))
    }

    /// Walks the table back from the last residue; ties prefer matching a
    /// label over skipping it and fewer letters over more.
    fn trace(&self, steps: &[Step]) -> Vec<Label> {
        enum At {
            Open(usize),
            Closed(usize),
            Match(usize),
            Letter(usize, usize),
        }
        let mut labels = vec![Label::new(0); self.residues.len()];
        let mut state = At::Open(self.track.len() - 1);
        let mut r = self.residues.len();
        while r > 0 {
            let step = &steps[r];
            state = match state {
                At::Open(j) => match best_index(&step.lettered[j]) {
                    Some(k) if step.lettered[j][k] > step.closed[j] => At::Letter(j, k + 1),
                    _ => At::Closed(j),
                },
                At::Closed(j) => match j > 0 && step.open[j - 1] + self.gap(j) > step.at[j] {
                    true => At::Open(j - 1),
                    false => At::Match(j),
                },
                At::Match(j) => {
                    labels[r - 1] = Label::new(self.track[j]);
                    r -= 1;
                    At::Open(j.saturating_sub(1))
                }
                At::Letter(j, k) => {
                    labels[r - 1] = Label::nth_insertion(self.track[j], k - 1);
                    r -= 1;
                    match k {
                        1 => At::Closed(j),
                        _ => At::Letter(j, k - 1),
                    }
                }
            };
        }
        labels
    }
}

/// Base label after which a scheme lets FR3 take letters.
fn fr3_site(scheme: Scheme) -> u16 {
    match scheme {
        Scheme::Martin => 72,
        _ => 82,
    }
}

/// Moves the FR3 letters between sites. The residues keep their order, so
/// the labels are sorted again.
fn move_fr3_letters(labels: &mut [Label], from: u16, to: u16) {
    for l in labels.iter_mut() {
        if let (Some(letter), true) = (l.insertion(), l.number == from) {
            *l = Label::with_insertion(to, letter);
        }
    }
    labels.sort_unstable();
}

/// Decides the CDR-H2 to FR3 split once, in Martin's frame, whichever
/// scheme is asked for: the split is a property of the residues, and the
/// schemes differ only in where FR3 letters sit.
fn place_boundary(scheme: Scheme, cur: &mut [Label], residues: &[u8]) {
    let (Some(range), Some(consensus)) = (stretch(cur, BOUNDARY), self::consensus(Scheme::Martin))
    else {
        return;
    };
    let mut rule = cur[range.clone()].to_vec();
    move_fr3_letters(&mut rule, fr3_site(scheme), fr3_site(Scheme::Martin));
    let sites = insertion_sites(Scheme::Martin, ChainType::Heavy);
    let window = &residues[range.clone()];
    let problem = Problem::new(consensus, window, &sites, BOUNDARY, &rule);
    let (Some(rule_score), Some((mut labels, best))) = (problem.score(&rule), problem.solve())
    else {
        return;
    };
    if best > rule_score + TIE_MARGIN {
        move_fr3_letters(&mut labels, fr3_site(Scheme::Martin), fr3_site(scheme));
        cur[range].copy_from_slice(&labels);
    }
}

/// `cur` moved to the best-scoring placement under the scheme's consensus;
/// unchanged for light chains, receptors and schemes without one.
pub(super) fn refine(
    scheme: Scheme,
    chain: ChainType,
    residues: &[u8],
    mut cur: Vec<Label>,
) -> Vec<Label> {
    let consensus = match (chain, consensus(scheme)) {
        (ChainType::Heavy, Some(c)) if residues.len() == cur.len() => c,
        _ => return cur,
    };
    slide_fr2(&mut cur, residues, consensus);
    place_boundary(scheme, &mut cur, residues);
    cur
}
