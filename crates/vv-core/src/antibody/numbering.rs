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

    /// Case-insensitive [`Self::name`]; Martin's scheme is also known as
    /// enhanced Chothia (Abhinandan & Martin 2008).
    pub fn parse(word: &str) -> Option<Scheme> {
        let bare: String = word.chars().filter(|c| c.is_alphanumeric()).collect();
        if bare.eq_ignore_ascii_case("enhancedchothia") {
            return Some(Scheme::Martin);
        }
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
    /// after `ins_after`, deletions taken in `delete` order after any
    /// labels `gaps` takes from the aligner.
    Var {
        lo: u16,
        hi: u16,
        base: (u16, u16),
        ins_after: u16,
        delete: &'static [u16],
        gaps: Gaps,
    },
}

/// Which empty base columns of a counted stretch the aligner decides.
#[derive(Clone, Copy)]
enum Gaps {
    /// None: `delete` alone orders the deletions.
    Canonical,
    /// Every empty column of the table `(lo, hi, offset)`.
    Aligner(&'static [(u16, u16, i16)]),
    /// Only the columns before the first residue (an N-terminal truncation).
    Leading(&'static [(u16, u16, i16)]),
    /// The empty column of the table when exactly one is empty; with more,
    /// `delete` orders them.
    Single(&'static [(u16, u16, i16)]),
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
        gaps: Gaps::Canonical,
    }
}

/// Counted framework stretch: insertions after `ins_after`, deletions as
/// `gaps` and `delete` say.
const fn counted(
    lo: u16,
    hi: u16,
    base: (u16, u16),
    ins_after: u16,
    gaps: Gaps,
    delete: &'static [u16],
) -> Span {
    Span::Var {
        lo,
        hi,
        base,
        ins_after,
        delete,
        gaps,
    }
}

/// FR1 of either chain: the aligner fixes where the chain starts, `delete`
/// where the residues missing inside it were.
const fn framework1(
    hi: u16,
    base_hi: u16,
    columns: &'static [(u16, u16, i16)],
    ins_after: u16,
    delete: &'static [u16],
) -> Span {
    Span::Var {
        lo: 1,
        hi,
        base: (1, base_hi),
        ins_after,
        delete,
        gaps: Gaps::Leading(columns),
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
/// Heavy FR2 with two or more residues missing, in every scheme (the
/// published Martin site is H42; the reference program deletes from H44).
const H_FR2_DELETE: &[u16] = &[44, 43, 42, 41, 40, 39, 38, 37, 36];
/// Where FR1 loses residues inside the chain: light chains at 10 (Kabat,
/// Chothia) or 7 (Martin), heavy chains at 10 in all three; then spreading
/// away from that site.
const L_FR1_KABAT: &[u16] = &[10, 9, 8, 7, 6, 5, 4, 3, 2, 1];
const L_FR1_MARTIN: &[u16] = &[7, 6, 5, 4, 3, 2, 1];
const H_FR1_DELETE: &[u16] = &[10, 9, 11, 8, 12, 7, 6, 5, 4, 3, 2, 1];

/// Framework stretches of the IMGT frame that carry base labels, as
/// `(first IMGT column, last, offset)`. Heavy IMGT 10 has no Kabat label.
const LIGHT_FR1: &[(u16, u16, i16)] = &[(1, 23, 0)];
const HEAVY_FR1: &[(u16, u16, i16)] = &[(1, 9, 0), (11, 26, -1)];
const HEAVY_FR2: &[(u16, u16, i16)] = &[(41, 54, -5)];
const HEAVY_FR3: &[(u16, u16, i16)] = &[(75, 91, -9), (95, 104, -12)];
const LIGHT_FR3: &[(u16, u16, i16)] = &[(70, 72, -13), (74, 80, -14), (83, 104, -16)];

fn heavy_spans(scheme: Scheme) -> Vec<Span> {
    if scheme == Scheme::Martin {
        return martin_heavy_spans();
    }
    let (h1_ins, h1_delete) = match scheme {
        Scheme::Kabat => (35, H1_KABAT),
        _ => (31, H1_CHOTHIA),
    };
    vec![
        framework1(26, 25, HEAVY_FR1, 9, H_FR1_DELETE),
        var(27, 40, (26, 35), h1_ins, h1_delete),
        counted(41, 54, (36, 49), 49, Gaps::Single(HEAVY_FR2), H_FR2_DELETE),
        var(55, 74, (50, 65), 52, H2_DELETE),
        counted(75, 104, (66, 92), 82, Gaps::Aligner(HEAVY_FR3), &[]),
        var(105, 117, (93, 102), 100, H3_DELETE),
        var(118, 128, (103, 113), 113, &[]),
    ]
}

/// Chothia plus indel sites at the framework positions the scheme moves:
/// H8, the FR2 deletion and the H72 insertion Kabat puts at H82.
fn martin_heavy_spans() -> Vec<Span> {
    vec![
        framework1(26, 25, HEAVY_FR1, 8, H_FR1_DELETE),
        var(27, 40, (26, 35), 31, H1_CHOTHIA),
        counted(41, 54, (36, 49), 49, Gaps::Single(HEAVY_FR2), H_FR2_DELETE),
        var(55, 74, (50, 65), 52, H2_DELETE),
        counted(75, 104, (66, 92), 72, Gaps::Aligner(HEAVY_FR3), &[]),
        var(105, 117, (93, 102), 100, H3_DELETE),
        var(118, 128, (103, 113), 113, &[]),
    ]
}

fn light_spans(scheme: Scheme, lambda_j: bool) -> Vec<Span> {
    if scheme == Scheme::Martin {
        return martin_light_spans(lambda_j);
    }
    let (l1_ins, l1_delete) = match scheme {
        Scheme::Kabat => (27, L1_KABAT),
        _ => (30, L1_CHOTHIA),
    };
    let mut spans = vec![framework1(23, 23, LIGHT_FR1, 9, L_FR1_KABAT)];
    spans.extend([
        var(24, 40, (24, 34), l1_ins, l1_delete),
        fixed(41, 55, -6),
        var(56, 69, (50, 56), 54, L2_DELETE),
        counted(70, 104, (57, 88), 66, Gaps::Aligner(LIGHT_FR3), &[]),
        var(105, 117, (89, 97), 95, L3_DELETE),
    ]);
    spans.extend(lambda_fr4(lambda_j));
    spans
}

/// Chothia plus the indel sites at L7 (lambda FR1), L40A/L41, L68 and
/// the shifted L1 and L2 deletions.
fn martin_light_spans(lambda_j: bool) -> Vec<Span> {
    let mut spans = vec![
        framework1(23, 23, LIGHT_FR1, 10, L_FR1_MARTIN),
        var(24, 40, (24, 34), 30, L1_MARTIN),
        var(41, 55, (35, 49), 40, &[41]),
        var(56, 69, (50, 56), 52, L2_MARTIN),
        counted(70, 104, (57, 88), 68, Gaps::Aligner(LIGHT_FR3), &[68]),
        var(105, 117, (89, 97), 95, L3_DELETE),
    ];
    spans.extend(lambda_fr4(lambda_j));
    spans
}

/// Lambda J segments carry 106A; kappa chains and kappa-type J segments
/// (Lys or Arg at IMGT 127) run straight through.
fn lambda_fr4(lambda_j: bool) -> Vec<Span> {
    match lambda_j {
        true => vec![
            fixed(118, 126, -20),
            extra(127, 127, 106),
            fixed(128, 128, -21),
        ],
        false => vec![var(118, 128, (98, 108), 108, &[])],
    }
}

fn spans(scheme: Scheme, chain: ChainType, lambda_j: bool) -> Vec<Span> {
    match chain {
        ChainType::Heavy => heavy_spans(scheme),
        ChainType::Kappa => light_spans(scheme, false),
        ChainType::Lambda => light_spans(scheme, lambda_j),
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

/// Internal gaps up to this many move to the scheme's canonical site; a
/// longer one is a real deletion and stays where the aligner put it.
const CANONICAL_GAPS: usize = 2;

/// Base labels of `gaps` that no residue of `run` (IMGT labels) fills.
fn aligner_gaps(base: (u16, u16), gaps: Gaps, run: &[Label]) -> Vec<u16> {
    let (Gaps::Aligner(columns) | Gaps::Leading(columns) | Gaps::Single(columns)) = gaps else {
        return Vec::new();
    };
    let filled: Vec<u16> = run
        .iter()
        .filter(|l| l.insertion().is_none())
        .filter_map(|l| {
            let n = l.number;
            columns
                .iter()
                .find(|(lo, hi, _)| (*lo..=*hi).contains(&n))
                .map(|(_, _, off)| (i32::from(n) + i32::from(*off)) as u16)
        })
        .collect();
    let empty: Vec<u16> = (base.0..=base.1).filter(|p| !filled.contains(p)).collect();
    match gaps {
        Gaps::Single(_) if empty.len() != 1 => return Vec::new(),
        Gaps::Leading(_) => {}
        _ => return empty,
    }
    let first = filled.iter().copied().min().unwrap_or(base.1 + 1);
    let (leading, internal): (Vec<u16>, Vec<u16>) = empty.into_iter().partition(|p| *p < first);
    match internal.len() > CANONICAL_GAPS {
        true => [leading, internal].concat(),
        false => leading,
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
pub(super) fn relabel(
    scheme: Scheme,
    chain: ChainType,
    imgt: &[Label],
    lambda_j: bool,
) -> Vec<Label> {
    match scheme {
        Scheme::Aho if chain.is_antibody() => return super::aho::relabel(chain, imgt),
        Scheme::Imgt | Scheme::Aho => return imgt.to_vec(),
        _ if !chain.is_antibody() => return imgt.to_vec(),
        _ => {}
    }
    let spans = spans(scheme, chain, lambda_j);
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
                gaps,
                ..
            } => {
                let mut order = aligner_gaps(*base, *gaps, &imgt[i..i + run]);
                order.extend_from_slice(delete);
                out.extend(var_labels(*base, *ins_after, &order, run));
            }
        }
        i += run;
    }
    out
}

#[cfg(test)]
mod tests {
    use super::Scheme;

    #[test]
    fn enhanced_chothia_is_martins_scheme() {
        for word in [
            "martin",
            "Martin",
            "enhancedchothia",
            "enhanced-chothia",
            "Enhanced Chothia",
        ] {
            assert_eq!(Scheme::parse(word), Some(Scheme::Martin), "{word}");
        }
        assert_eq!(Scheme::parse("chothia"), Some(Scheme::Chothia));
    }
}
