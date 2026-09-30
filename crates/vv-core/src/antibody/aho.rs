//! AHo numbering (Honegger & Plückthun 2001, J Mol Biol 309:657) as a
//! relabelling of the IMGT frame.
//!
//! Every V domain gets positions 1 to 149. Framework columns are fixed
//! offsets from IMGT; each length-variable stretch has a fixed span of
//! positions and a gap placed around a centre the paper names (positions 8,
//! 28, 36, 63, 85-86, 123). Where the paper only draws the result (Figure 1),
//! the split of residues either side of a gap was read off that figure's
//! sequence alignment and is pinned by `tests/antibody_aho.rs`.

use super::numbering::Label;
use super::ChainType;

/// IMGT framework stretch and its AHo offset.
struct Shift {
    lo: u16,
    hi: u16,
    offset: i16,
}

const FIXED: [Shift; 7] = [
    Shift {
        lo: 1,
        hi: 7,
        offset: 0,
    },
    Shift {
        lo: 11,
        hi: 26,
        offset: 0,
    },
    Shift {
        lo: 39,
        hi: 55,
        offset: 2,
    },
    Shift {
        lo: 66,
        hi: 72,
        offset: 3,
    },
    Shift {
        lo: 74,
        hi: 80,
        offset: 2,
    },
    Shift {
        lo: 87,
        hi: 104,
        offset: 2,
    },
    Shift {
        lo: 118,
        hi: 128,
        offset: 21,
    },
];

/// Stretches whose residue count varies between chains; a residue run
/// there is placed by [`Stretch::place`].
#[derive(Clone, Copy)]
enum Stretch {
    /// IMGT 8-10: AHo 8-10, the first columns left empty (kappa fills all).
    NTerm,
    Cdr1,
    Cdr2,
    /// IMGT 73, which AHo has no column for.
    Extra73,
    /// IMGT 81-86: AHo 83-88 with a gap around 85-86.
    Hairpin,
    Cdr3,
}

const STRETCHES: [(u16, u16, Stretch); 6] = [
    (8, 10, Stretch::NTerm),
    (27, 38, Stretch::Cdr1),
    (56, 65, Stretch::Cdr2),
    (73, 73, Stretch::Extra73),
    (81, 86, Stretch::Hairpin),
    (105, 117, Stretch::Cdr3),
];

fn fixed_label(number: u16) -> Option<Label> {
    FIXED
        .iter()
        .find(|s| (s.lo..=s.hi).contains(&number))
        .map(|s| Label::new((i32::from(number) + i32::from(s.offset)) as u16))
}

fn stretch_of(number: u16) -> Option<Stretch> {
    STRETCHES
        .iter()
        .find(|(lo, hi, _)| (*lo..=*hi).contains(&number))
        .map(|s| s.2)
}

/// `n` residues in columns `lo..=hi`: `left` of them from the low end, the
/// rest from the high end. Past capacity, the surplus takes insertion
/// letters after `axis`.
fn spread(n: usize, (lo, hi): (u16, u16), axis: u16, left: usize) -> Vec<Label> {
    let capacity = usize::from(hi - lo) + 1;
    if n > capacity {
        let (before, after) = (usize::from(axis - lo), usize::from(hi - axis));
        let middle = n - before - after;
        let mut out: Vec<Label> = (lo..axis).map(Label::new).collect();
        out.push(Label::new(axis));
        out.extend((0..middle - 1).map(|k| Label::nth_insertion(axis, k)));
        out.extend((axis + 1..=hi).map(Label::new));
        return out;
    }
    let right = n - left;
    (lo..lo + left as u16)
        .chain(hi + 1 - right as u16..=hi)
        .map(Label::new)
        .collect()
}

/// First column of a `gap`-wide gap centred on `axis`, the extra column of
/// an even gap on the low side.
fn gap_start(axis: u16, gap: usize) -> u16 {
    axis - (gap / 2) as u16
}

/// Residues on the low side of a gap centred on `axis` in columns `lo..=hi`.
fn left_of_gap(n: usize, (lo, hi): (u16, u16), axis: u16) -> usize {
    let gap = usize::from(hi - lo) + 1 - n.min(usize::from(hi - lo) + 1);
    usize::from(gap_start(axis, gap) - lo)
}

/// Residues of CDR1 in the outer loop (columns 27-30, before the residue
/// fixed at 31): heavy and kappa have a fixed count, lambda germlines fall
/// into a short-loop and a long-loop group.
fn outer_count(chain: ChainType, cdr1_len: usize) -> usize {
    let n = match chain {
        ChainType::Heavy => 3,
        ChainType::Kappa => 2,
        ChainType::Lambda if cdr1_len <= 6 => 1,
        ChainType::Lambda | ChainType::TcrAlpha | ChainType::TcrBeta => 3,
    };
    n.min(cdr1_len.saturating_sub(1))
}

fn cdr1(chain: ChainType, n: usize) -> Vec<Label> {
    let outer = outer_count(chain, n);
    let mut out = spread(outer, (27, 30), 28, outer_left(outer));
    if n > outer {
        out.push(Label::new(31));
    }
    let inner = n.saturating_sub(outer + 1);
    out.extend(spread(
        inner,
        (32, 40),
        36,
        left_of_gap(inner, (32, 40), 36),
    ));
    out
}

/// Outer-loop residues fill the columns either side of the gap at 28.
fn outer_left(outer: usize) -> usize {
    usize::from(gap_start(28, 4 - outer).saturating_sub(27))
}

fn cdr2(n: usize) -> Vec<Label> {
    spread(n, (58, 68), 63, left_of_gap(n, (58, 68), 63))
}

/// The first two and the last residue of CDR3 sit in fixed columns
/// (107, 108, 138); the rest split evenly around 123, extra on the left.
fn cdr3(n: usize) -> Vec<Label> {
    let head = n.min(2);
    let tail = usize::from(n > 2);
    let middle = n - head - tail;
    let mut out: Vec<Label> = (107..107 + head as u16).map(Label::new).collect();
    out.extend(spread(middle, (109, 137), 123, middle.div_ceil(2)));
    out.extend((0..tail).map(|_| Label::new(138)));
    out
}

fn n_term(n: usize) -> Vec<Label> {
    (11 - n.min(3) as u16..=10).map(Label::new).collect()
}

fn hairpin(n: usize) -> Vec<Label> {
    let gap = 6usize.saturating_sub(n);
    let start = gap_start(86, gap);
    let mut out: Vec<Label> = (83..=88)
        .filter(|c| gap == 0 || !(start..start + gap as u16).contains(c))
        .map(Label::new)
        .collect();
    out.extend((0..n.saturating_sub(6)).map(|k| Label::nth_insertion(88, k)));
    out
}

impl Stretch {
    fn place(self, chain: ChainType, n: usize) -> Vec<Label> {
        match self {
            Stretch::NTerm => n_term(n),
            Stretch::Cdr1 => cdr1(chain, n),
            Stretch::Cdr2 => cdr2(n),
            Stretch::Extra73 => (0..n).map(|k| Label::nth_insertion(75, k)).collect(),
            Stretch::Hairpin => hairpin(n),
            Stretch::Cdr3 => cdr3(n),
        }
    }
}

/// AHo label of every residue, given its IMGT label, in sequence order.
pub(super) fn relabel(chain: ChainType, imgt: &[Label]) -> Vec<Label> {
    let mut out = Vec::with_capacity(imgt.len());
    let mut i = 0;
    while i < imgt.len() {
        let number = imgt[i].number;
        let Some(stretch) = stretch_of(number) else {
            out.push(fixed_label(number).unwrap_or(Label::new(number)));
            i += 1;
            continue;
        };
        let run = imgt[i..]
            .iter()
            .take_while(|l| stretch_of(l.number).is_some_and(|s| same(s, stretch)))
            .count();
        out.extend(stretch.place(chain, run));
        i += run;
    }
    out
}

fn same(a: Stretch, b: Stretch) -> bool {
    std::mem::discriminant(&a) == std::mem::discriminant(&b)
}
