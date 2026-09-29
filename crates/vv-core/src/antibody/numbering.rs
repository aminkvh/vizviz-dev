//! Numbering schemes as label rules over the IMGT frame: fixed framework
//! offsets plus length-driven loops with a scheme-specific insertion site.

use std::fmt;

use super::align::Slot;
use super::profile::column_number;
use super::ChainType;

/// A residue number and optional insertion letter, e.g. `52A`.
///
/// For IMGT the letter stands for the numeric insertion index (`111A` is
/// IMGT `111.1`, `112B` is `112.2`).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Label {
    pub number: u16,
    insertion: u8,
}

impl Label {
    pub const fn new(number: u16) -> Self {
        Label {
            number,
            insertion: 0,
        }
    }

    pub const fn with_insertion(number: u16, letter: char) -> Self {
        Label {
            number,
            insertion: letter as u8,
        }
    }

    pub fn insertion(self) -> Option<char> {
        (self.insertion != 0).then_some(self.insertion as char)
    }

    /// `k`-th insertion letter: `A`..`Z`, then `a`..`z`.
    fn nth_insertion(number: u16, k: usize) -> Self {
        let letter = match k {
            0..=25 => b'A' + k as u8,
            26..=51 => b'a' + (k - 26) as u8,
            _ => b'~',
        };
        Label {
            number,
            insertion: letter,
        }
    }
}

impl fmt::Display for Label {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.number)?;
        self.insertion().map_or(Ok(()), |c| write!(f, "{c}"))
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Scheme {
    Imgt,
    Kabat,
    Chothia,
    Martin,
}

impl Scheme {
    pub const ALL: [Scheme; 4] = [Scheme::Kabat, Scheme::Chothia, Scheme::Imgt, Scheme::Martin];

    pub fn name(self) -> &'static str {
        match self {
            Scheme::Imgt => "IMGT",
            Scheme::Kabat => "Kabat",
            Scheme::Chothia => "Chothia",
            Scheme::Martin => "Martin",
        }
    }

    /// Case-insensitive [`Self::name`].
    pub fn parse(word: &str) -> Option<Scheme> {
        Self::ALL
            .into_iter()
            .find(|s| s.name().eq_ignore_ascii_case(word))
    }
}

/// `(first, last base position on the left of the axis, last)` of each
/// IMGT CDR loop; shorter or longer loops grow symmetrically around it.
const IMGT_LOOPS: [(u16, u16, u16); 3] = [(27, 32, 38), (56, 60, 65), (105, 111, 117)];

/// IMGT labels for a CDR loop of `n` residues: gaps and insertions are
/// placed symmetrically at the top of the loop, with the extra residue of
/// an odd length on the N-terminal side (Lefranc et al. 2003).
fn imgt_loop_labels(loop_id: usize, n: usize) -> Vec<Label> {
    let (first, axis, last) = IMGT_LOOPS[loop_id];
    let left = n.div_ceil(2);
    let left_base = usize::from(axis - first) + 1;
    let right_base = usize::from(last - axis);
    let mut out: Vec<Label> = (0..left)
        .map(|t| match t < left_base {
            true => Label::new(first + t as u16),
            false => Label::nth_insertion(axis, t - left_base),
        })
        .collect();
    let mut right: Vec<Label> = (0..n - left)
        .map(|t| match t < right_base {
            true => Label::new(last - t as u16),
            false => Label::nth_insertion(axis + 1, t - right_base),
        })
        .collect();
    right.reverse();
    out.extend(right);
    out
}

/// IMGT label of every aligned residue, in sequence order.
pub(super) fn imgt_labels(slots: &[Slot]) -> Vec<Label> {
    let mut out = Vec::with_capacity(slots.len());
    let mut i = 0;
    while i < slots.len() {
        match slots[i] {
            Slot::Col(c) => {
                out.push(Label::new(column_number(usize::from(c))));
                i += 1;
            }
            Slot::Loop(b) => {
                let run = slots[i..]
                    .iter()
                    .take_while(|s| **s == Slot::Loop(b))
                    .count();
                out.extend(imgt_loop_labels(usize::from(b), run));
                i += run;
            }
        }
    }
    out
}

/// A stretch of the IMGT frame and how a scheme labels its residues.
enum Span {
    /// One label per IMGT position, `number = imgt + offset`.
    Fixed { lo: u16, hi: u16, offset: i16 },
    /// Residues counted, not indexed: `base` labels, insertions lettered
    /// after `ins_after`, deletions taken in `delete` order.
    Var {
        lo: u16,
        hi: u16,
        base: (u16, u16),
        ins_after: u16,
        delete: &'static [u16],
    },
}

const fn fixed(lo: u16, hi: u16, offset: i16) -> Span {
    Span::Fixed { lo, hi, offset }
}

const fn var(lo: u16, hi: u16, base: (u16, u16), ins_after: u16, delete: &'static [u16]) -> Span {
    Span::Var {
        lo,
        hi,
        base,
        ins_after,
        delete,
    }
}

/// Insertion-only stretch (framework columns a germline leaves empty).
const fn extra(col: u16, hi: u16, ins_after: u16) -> Span {
    var(col, hi, (0, 0), ins_after, &[])
}

const H1_DELETE: &[u16] = &[32, 33, 31, 34, 30, 35];
const H2_DELETE: &[u16] = &[54, 53, 55, 52, 56, 51];
const H3_DELETE: &[u16] = &[100, 99, 98, 97, 96, 95];
const L1_DELETE: &[u16] = &[28, 29, 30, 27, 31];
const L2_DELETE: &[u16] = &[54, 53, 55, 52];
const L3_DELETE: &[u16] = &[95, 94, 96, 93];

fn heavy_spans(scheme: Scheme) -> Vec<Span> {
    let h1_ins = if scheme == Scheme::Kabat { 35 } else { 31 };
    vec![
        fixed(1, 9, 0),
        extra(10, 10, 9),
        fixed(11, 26, -1),
        var(27, 40, (26, 35), h1_ins, H1_DELETE),
        fixed(41, 54, -5),
        var(55, 74, (50, 65), 52, H2_DELETE),
        fixed(75, 91, -9),
        extra(92, 94, 82),
        fixed(95, 104, -12),
        var(105, 117, (93, 102), 100, H3_DELETE),
        fixed(118, 128, -15),
    ]
}

fn light_spans(scheme: Scheme, lambda: bool) -> Vec<Span> {
    let l1_ins = if scheme == Scheme::Kabat { 27 } else { 30 };
    let l2_ins = if scheme == Scheme::Martin { 52 } else { 54 };
    let mut spans = if lambda {
        vec![fixed(1, 9, 0), extra(10, 10, 9), fixed(11, 23, 0)]
    } else {
        vec![fixed(1, 23, 0)]
    };
    spans.extend([
        var(24, 40, (24, 34), l1_ins, L1_DELETE),
        fixed(41, 55, -6),
        var(56, 69, (50, 56), l2_ins, L2_DELETE),
        fixed(70, 72, -13),
        extra(73, 73, 59),
        fixed(74, 80, -14),
        extra(81, 82, 66),
        fixed(83, 104, -16),
        var(105, 117, (89, 97), 95, L3_DELETE),
    ]);
    spans.extend(if lambda {
        vec![
            fixed(118, 126, -20),
            extra(127, 127, 106),
            fixed(128, 128, -21),
        ]
    } else {
        vec![fixed(118, 128, -20)]
    });
    spans
}

fn spans(scheme: Scheme, chain: ChainType) -> Vec<Span> {
    match chain {
        ChainType::Heavy => heavy_spans(scheme),
        ChainType::Kappa => light_spans(scheme, false),
        ChainType::Lambda => light_spans(scheme, true),
    }
}

/// Labels for `m` residues filling `base`, with insertion letters after
/// `ins_after` or deletions in `delete` order.
fn var_labels(base: (u16, u16), ins_after: u16, delete: &[u16], m: usize) -> Vec<Label> {
    let mut positions: Vec<u16> = if base.1 == 0 {
        Vec::new()
    } else {
        (base.0..=base.1).collect()
    };
    if m < positions.len() {
        for gone in delete {
            if positions.len() == m {
                break;
            }
            positions.retain(|p| p != gone);
        }
        positions.truncate(m);
    }
    let extra_count = m.saturating_sub(positions.len());
    let letters = |n| (0..n).map(|k| Label::nth_insertion(ins_after, k));
    match positions.iter().position(|&p| p == ins_after) {
        Some(at) => {
            let mut out: Vec<Label> = positions[..=at].iter().map(|&p| Label::new(p)).collect();
            out.extend(letters(extra_count));
            out.extend(positions[at + 1..].iter().map(|&p| Label::new(p)));
            out
        }
        None => letters(extra_count)
            .chain(positions.iter().map(|&p| Label::new(p)))
            .collect(),
    }
}

fn span_of(spans: &[Span], number: u16) -> usize {
    spans
        .iter()
        .position(|s| match s {
            Span::Fixed { lo, hi, .. } | Span::Var { lo, hi, .. } => (*lo..=*hi).contains(&number),
        })
        .unwrap_or(spans.len() - 1)
}

/// Relabels residues (given by their IMGT labels) in a non-IMGT scheme.
pub(super) fn relabel(scheme: Scheme, chain: ChainType, imgt: &[Label]) -> Vec<Label> {
    if scheme == Scheme::Imgt {
        return imgt.to_vec();
    }
    let spans = spans(scheme, chain);
    let mut out = Vec::with_capacity(imgt.len());
    let mut i = 0;
    while i < imgt.len() {
        let s = span_of(&spans, imgt[i].number);
        let run = imgt[i..]
            .iter()
            .take_while(|l| span_of(&spans, l.number) == s)
            .count();
        match &spans[s] {
            Span::Fixed { offset, .. } => {
                out.extend(
                    imgt[i..i + run]
                        .iter()
                        .map(|l| Label::new((i32::from(l.number) + i32::from(*offset)) as u16)),
                );
            }
            Span::Var {
                base,
                ins_after,
                delete,
                ..
            } => {
                out.extend(var_labels(*base, *ins_after, delete, run));
            }
        }
        i += run;
    }
    out
}
