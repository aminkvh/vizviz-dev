//! Global pairwise alignment of residue-letter sequences (Needleman and
//! Wunsch 1970) with affine gap costs (Gotoh 1982). Small and dependency
//! free: it places unobserved residues against a chain's full sequence and
//! lines homologous chains up for conservation.

const MATCH: i32 = 3;
const MISMATCH: i32 = -1;
const GAP_OPEN: i32 = -6;
const GAP_EXTEND: i32 = -1;
/// Stand-in for minus infinity that cannot overflow when a cost is added.
const NEG: i32 = i32::MIN / 4;
/// Largest `a.len() * b.len()` aligned: three score tables of this many
/// cells are held while aligning.
const MAX_CELLS: usize = 4_000_000;

fn is_unknown(c: u8) -> bool {
    c == b'x' || c == b'X'
}

/// An unknown residue (`x`, `X`) matches nothing and costs nothing.
fn substitution(a: u8, b: u8) -> i32 {
    if is_unknown(a) || is_unknown(b) {
        0
    } else if a == b {
        MATCH
    } else {
        MISMATCH
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum State {
    Pair,
    /// `a[i]` against a gap.
    GapInB,
    /// `b[j]` against a gap.
    GapInA,
}

struct Tables {
    cols: usize,
    pair: Vec<i32>,
    gap_in_b: Vec<i32>,
    gap_in_a: Vec<i32>,
}

impl Tables {
    fn at(&self, i: usize, j: usize) -> usize {
        i * self.cols + j
    }

    fn best(&self, i: usize, j: usize) -> i32 {
        let k = self.at(i, j);
        self.pair[k].max(self.gap_in_b[k]).max(self.gap_in_a[k])
    }
}

fn fill(a: &[u8], b: &[u8]) -> Tables {
    let (rows, cols) = (a.len() + 1, b.len() + 1);
    let mut t = Tables {
        cols,
        pair: vec![NEG; rows * cols],
        gap_in_b: vec![NEG; rows * cols],
        gap_in_a: vec![NEG; rows * cols],
    };
    t.pair[0] = 0;
    for i in 1..rows {
        t.gap_in_b[i * cols] = GAP_OPEN + GAP_EXTEND * (i as i32 - 1);
    }
    for j in 1..cols {
        t.gap_in_a[j] = GAP_OPEN + GAP_EXTEND * (j as i32 - 1);
    }
    for i in 1..rows {
        for j in 1..cols {
            let k = t.at(i, j);
            t.pair[k] = t.best(i - 1, j - 1) + substitution(a[i - 1], b[j - 1]);
            let up = t.at(i - 1, j);
            t.gap_in_b[k] = (t.best(i - 1, j) + GAP_OPEN).max(t.gap_in_b[up] + GAP_EXTEND);
            let left = t.at(i, j - 1);
            t.gap_in_a[k] = (t.best(i, j - 1) + GAP_OPEN).max(t.gap_in_a[left] + GAP_EXTEND);
        }
    }
    t
}

fn best_state(t: &Tables, i: usize, j: usize) -> State {
    let k = t.at(i, j);
    let best = t.best(i, j);
    if t.pair[k] == best {
        State::Pair
    } else if t.gap_in_b[k] == best {
        State::GapInB
    } else {
        State::GapInA
    }
}

/// The state the traceback was in before `state` at `(i, j)`.
fn previous_state(t: &Tables, state: State, i: usize, j: usize) -> State {
    let (pi, pj, table) = match state {
        State::Pair => return best_state(t, i - 1, j - 1),
        State::GapInB => (i - 1, j, &t.gap_in_b),
        State::GapInA => (i, j - 1, &t.gap_in_a),
    };
    if table[t.at(i, j)] == table[t.at(pi, pj)] + GAP_EXTEND {
        state
    } else {
        best_state(t, pi, pj)
    }
}

/// For each position of `a`, the position of `b` it is aligned to, or
/// `None` where it faces a gap. All `None` when the sequences are too long
/// to align in bounded memory.
pub fn align(a: &[u8], b: &[u8]) -> Vec<Option<usize>> {
    let mut map = vec![None; a.len()];
    if a.is_empty() || b.is_empty() || a.len().saturating_mul(b.len()) > MAX_CELLS {
        return map;
    }
    let t = fill(a, b);
    let (mut i, mut j) = (a.len(), b.len());
    let mut state = best_state(&t, i, j);
    while i > 0 && j > 0 {
        let next = previous_state(&t, state, i, j);
        match state {
            State::Pair => {
                map[i - 1] = Some(j - 1);
                i -= 1;
                j -= 1;
            }
            State::GapInB => i -= 1,
            State::GapInA => j -= 1,
        }
        state = next;
    }
    map
}

/// Identical aligned pairs, given `map` from [`align`] of `a` onto `b`.
pub fn identical(a: &[u8], b: &[u8], map: &[Option<usize>]) -> usize {
    map.iter()
        .enumerate()
        .filter(|&(i, m)| m.is_some_and(|j| a[i] == b[j] && !is_unknown(a[i])))
        .count()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn identical_sequences_align_position_by_position() {
        let s = b"MKVLAAGIV";
        let map = align(s, s);
        assert!(map.iter().enumerate().all(|(i, m)| *m == Some(i)));
    }

    #[test]
    fn a_missing_stretch_becomes_one_gap() {
        let map = align(b"MKVLAAGIVWWTTQQ", b"MKVLAAQQ");
        assert_eq!(map[5], Some(5));
        assert_eq!(map.iter().filter(|m| m.is_none()).count(), 7);
        assert_eq!(map[13], Some(6));
        assert_eq!(map[14], Some(7));
    }

    #[test]
    fn a_substitution_stays_aligned() {
        let map = align(b"MKVLAAG", b"MKVIAAG");
        assert!(map.iter().enumerate().all(|(i, m)| *m == Some(i)));
        assert_eq!(identical(b"MKVLAAG", b"MKVIAAG", &map), 6);
    }

    #[test]
    fn oversized_input_is_left_unaligned() {
        let big = vec![b'A'; 2100];
        assert!(align(&big, &big).iter().all(Option::is_none));
    }
}
