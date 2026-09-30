//! Numbering schemes as label rules over the IMGT frame: fixed framework
//! offsets plus length-driven loops with a scheme-specific insertion site.

use std::fmt;

use super::align::Slot;
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
    pub(super) fn nth_insertion(number: u16, k: usize) -> Self {
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
    Aho,
}

impl Scheme {
    pub const ALL: [Scheme; 5] = [
        Scheme::Kabat,
        Scheme::Chothia,
        Scheme::Imgt,
        Scheme::Martin,
        Scheme::Aho,
    ];

    pub fn name(self) -> &'static str {
        match self {
            Scheme::Imgt => "IMGT",
            Scheme::Kabat => "Kabat",
            Scheme::Chothia => "Chothia",
            Scheme::Martin => "Martin",
            Scheme::Aho => "AHo",
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

/// IMGT label of every aligned residue, in sequence order; `columns` labels
/// the profile columns the slots index.
pub(super) fn imgt_labels(slots: &[Slot], columns: &[Label]) -> Vec<Label> {
    let mut out = Vec::with_capacity(slots.len());
    let mut i = 0;
    while i < slots.len() {
        match slots[i] {
            Slot::Col(c) => {
                out.push(columns[usize::from(c)]);
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

/// Labels a loop or framework stretch gives up first when shorter than its
/// base run. Measured on structures numbered by each scheme's own authors
/// (docs/ANTIBODY.md, "Deletion order"): the orders run away from the
/// insertion site, except Kabat H1 and L1, which run from the loop's end.
const H1_KABAT: &[u16] = &[35, 34, 33, 32, 31, 30, 29, 28, 27, 26];
const H1_CHOTHIA: &[u16] = &[31, 30, 29, 28, 27, 26, 32, 33, 34, 35];
const H2_DELETE: &[u16] = &[53, 54, 55, 56, 57, 58, 52, 51, 50];
const H3_DELETE: &[u16] = &[100, 99, 98, 97, 96, 95, 94, 93];
const L1_KABAT: &[u16] = &[28, 29, 30, 31, 27, 26, 25, 24];
const L1_CHOTHIA: &[u16] = &[31, 32, 33, 34, 30, 29, 28, 27, 26, 25];
const L1_MARTIN: &[u16] = &[30, 29, 28, 27, 26, 25, 31, 32, 33, 34];
const L2_DELETE: &[u16] = &[54, 53, 55, 52, 56, 51, 50];
const L2_MARTIN: &[u16] = &[52, 51, 50, 53, 54, 55, 56];
const L3_DELETE: &[u16] = &[95, 94, 93, 92, 91, 90, 89];
/// Martin heavy FR2 (Abhinandan & Martin 2008 name H42; the depositions
/// measured here delete from H44).
const H_FR2_MARTIN: &[u16] = &[44, 43, 42, 41, 40, 39, 38, 37, 36];

fn heavy_spans(scheme: Scheme) -> Vec<Span> {
    if scheme == Scheme::Martin {
        return martin_heavy_spans();
    }
    let (h1_ins, h1_delete) = match scheme {
        Scheme::Kabat => (35, H1_KABAT),
        _ => (31, H1_CHOTHIA),
    };
    vec![
        fixed(1, 9, 0),
        extra(10, 10, 9),
        fixed(11, 26, -1),
        var(27, 40, (26, 35), h1_ins, h1_delete),
        fixed(41, 54, -5),
        var(55, 74, (50, 65), 52, H2_DELETE),
        fixed(75, 91, -9),
        extra(92, 94, 82),
        fixed(95, 104, -12),
        var(105, 117, (93, 102), 100, H3_DELETE),
        fixed(118, 128, -15),
    ]
}

/// Chothia plus indel sites at the framework positions the scheme moves:
/// H8, the FR2 deletion and the H72 insertion Kabat puts at H82.
fn martin_heavy_spans() -> Vec<Span> {
    vec![
        fixed(1, 7, 0),
        var(8, 10, (8, 9), 8, &[8]),
        fixed(11, 26, -1),
        var(27, 40, (26, 35), 31, H1_CHOTHIA),
        var(41, 54, (36, 49), 49, H_FR2_MARTIN),
        var(55, 74, (50, 65), 52, H2_DELETE),
        var(75, 104, (66, 92), 72, &[]),
        var(105, 117, (93, 102), 100, H3_DELETE),
        fixed(118, 128, -15),
    ]
}

fn light_spans(scheme: Scheme, lambda: bool) -> Vec<Span> {
    if scheme == Scheme::Martin {
        return martin_light_spans(lambda);
    }
    let (l1_ins, l1_delete) = match scheme {
        Scheme::Kabat => (27, L1_KABAT),
        _ => (30, L1_CHOTHIA),
    };
    let mut spans = if lambda {
        vec![fixed(1, 9, 0), extra(10, 10, 9), fixed(11, 23, 0)]
    } else {
        vec![fixed(1, 23, 0)]
    };
    spans.extend([
        var(24, 40, (24, 34), l1_ins, l1_delete),
        fixed(41, 55, -6),
        var(56, 69, (50, 56), 54, L2_DELETE),
        fixed(70, 72, -13),
        extra(73, 73, 59),
        fixed(74, 80, -14),
        extra(81, 82, 66),
        fixed(83, 104, -16),
        var(105, 117, (89, 97), 95, L3_DELETE),
    ]);
    spans.extend(lambda_fr4(lambda));
    spans
}

/// Chothia plus the indel sites at L7 (lambda FR1), L40A/L41, L68 and
/// the shifted L1 and L2 deletions.
fn martin_light_spans(lambda: bool) -> Vec<Span> {
    let mut spans = vec![
        fixed(1, 6, 0),
        var(7, 10, (7, 10), 10, &[7]),
        fixed(11, 23, 0),
        var(24, 40, (24, 34), 30, L1_MARTIN),
        var(41, 55, (35, 49), 40, &[41]),
        var(56, 69, (50, 56), 52, L2_MARTIN),
        var(70, 104, (57, 88), 68, &[68]),
        var(105, 117, (89, 97), 95, L3_DELETE),
    ];
    spans.extend(lambda_fr4(lambda));
    spans
}

/// Kabat's lambda FR4 carries 106A; kappa runs straight through.
fn lambda_fr4(lambda: bool) -> Vec<Span> {
    match lambda {
        true => vec![
            fixed(118, 126, -20),
            extra(127, 127, 106),
            fixed(128, 128, -21),
        ],
        false => vec![fixed(118, 128, -20)],
    }
}

fn spans(scheme: Scheme, chain: ChainType) -> Vec<Span> {
    match chain {
        ChainType::Heavy => heavy_spans(scheme),
        ChainType::Kappa => light_spans(scheme, false),
        ChainType::Lambda => light_spans(scheme, true),
        ChainType::TcrAlpha | ChainType::TcrBeta => Vec::new(),
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
    match scheme {
        Scheme::Aho if chain.is_antibody() => return super::aho::relabel(chain, imgt),
        Scheme::Imgt | Scheme::Aho => return imgt.to_vec(),
        _ if !chain.is_antibody() => return imgt.to_vec(),
        _ => {}
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
