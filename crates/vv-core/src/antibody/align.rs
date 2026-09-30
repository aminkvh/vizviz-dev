//! Local alignment of a sequence to a framework profile: match/delete
//! states with affine gaps inside the framework, and free-length CDR loops
//! between framework blocks priced by a length prior.

use super::profile::{Profile, C23};

const NEG: f32 = -1.0e9;
const GAP_OPEN: f32 = 10.0;
const GAP_EXTEND: f32 = 2.0;

/// Where an aligned residue sits in the profile.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Slot {
    /// Framework column index into the profile.
    Col(u8),
    /// Residue of CDR loop 0, 1 or 2.
    Loop(u8),
}

pub(super) struct Alignment {
    pub score: f32,
    /// `(residue index, slot)` in sequence order; contiguous residues.
    pub slots: Vec<(usize, Slot)>,
}

/// Predecessor of a cell, for traceback.
#[derive(Clone, Copy)]
enum From {
    Start,
    Match,
    Delete,
    Loop,
}

/// One profile column's states, indexed by rows (residues consumed).
struct Column {
    m: Vec<f32>,
    d: Vec<f32>,
    m_from: Vec<From>,
    d_from: Vec<From>,
    /// Loop entering this column: best score per row and its source
    /// `(row, source was a deletion)`.
    entry: Vec<f32>,
    entry_src: Vec<(u32, bool)>,
}

impl Column {
    fn new(rows: usize) -> Self {
        Column {
            m: vec![NEG; rows],
            d: vec![NEG; rows],
            m_from: vec![From::Start; rows],
            d_from: vec![From::Start; rows],
            entry: Vec::new(),
            entry_src: Vec::new(),
        }
    }

    fn best(&self, row: usize) -> (f32, bool) {
        if self.d[row] > self.m[row] {
            (self.d[row], true)
        } else {
            (self.m[row], false)
        }
    }
}

fn gap_costs(profile: &Profile, col: usize) -> (f32, f32) {
    if profile.free_gap[col] {
        (0.0, 0.0)
    } else {
        (GAP_OPEN, GAP_EXTEND)
    }
}

/// Fills `col.entry`: best score of leaving the previous block after
/// `row` residues, having spent the loop's residues in between.
fn fill_loop_entry(col: &mut Column, prev: &Column, profile: &Profile, loop_id: usize) {
    let rows = prev.m.len();
    let costs = &profile.loop_costs[loop_id];
    let leave: Vec<(f32, bool)> = (0..rows).map(|r| prev.best(r)).collect();
    col.entry = vec![NEG; rows];
    col.entry_src = vec![(0, false); rows];
    for row in 0..rows {
        let (mut best, mut src) = (NEG, (0, false));
        for len in 0..=(costs.len() - 1).min(row) {
            let val = leave[row - len].0 + costs[len];
            if val > best {
                best = val;
                src = ((row - len) as u32, leave[row - len].1);
            }
        }
        col.entry[row] = best;
        col.entry_src[row] = src;
    }
}

fn fill_column(q: &[u8], profile: &Profile, idx: usize, prev: Option<&Column>, col: &mut Column) {
    let (open, ext) = gap_costs(profile, idx);
    let scores = &profile.score[idx];
    let entering_loop = !col.entry.is_empty();
    for row in 0..col.m.len() {
        if row > 0 {
            let (arrive, from) = match prev {
                None => (NEG, From::Start),
                Some(_) if entering_loop => (col.entry[row - 1], From::Loop),
                Some(p) => match p.best(row - 1) {
                    (v, true) => (v, From::Delete),
                    (v, false) => (v, From::Match),
                },
            };
            let free_start = if idx <= C23 { 0.0 } else { NEG };
            let (base, from) = if free_start >= arrive {
                (free_start, From::Start)
            } else {
                (arrive, from)
            };
            if base > NEG / 2.0 {
                col.m[row] = base + scores[usize::from(q[row - 1])];
                col.m_from[row] = from;
            }
        }
        let Some(p) = prev else { continue };
        let (via_m, from_m) = if entering_loop {
            (col.entry[row] - open, From::Loop)
        } else {
            (p.m[row] - open, From::Match)
        };
        let via_d = if entering_loop { NEG } else { p.d[row] - ext };
        let (val, from) = if via_d > via_m {
            (via_d, From::Delete)
        } else {
            (via_m, from_m)
        };
        if val > NEG / 2.0 {
            col.d[row] = val;
            col.d_from[row] = from;
        }
    }
}

fn fill(q: &[u8], profile: &Profile) -> Vec<Column> {
    let n_cols = profile.score.len();
    let mut cols: Vec<Column> = Vec::with_capacity(n_cols);
    for idx in 0..n_cols {
        let mut col = Column::new(q.len() + 1);
        if let Some(loop_id) = profile.loop_exit.iter().position(|&e| e == idx) {
            fill_loop_entry(&mut col, &cols[idx - 1], profile, loop_id);
        }
        fill_column(q, profile, idx, cols.last(), &mut col);
        cols.push(col);
    }
    cols
}

/// Best alignment end: a match in FR4 or later.
fn best_end(cols: &[Column], fr4_start: usize) -> Option<(usize, usize, f32)> {
    let mut best: Option<(usize, usize, f32)> = None;
    for (idx, col) in cols.iter().enumerate().skip(fr4_start) {
        for (row, &s) in col.m.iter().enumerate().skip(1) {
            if s > NEG / 2.0 && best.is_none_or(|(_, _, b)| s > b) {
                best = Some((row, idx, s));
            }
        }
    }
    best
}

fn trace(
    cols: &[Column],
    profile: &Profile,
    mut row: usize,
    mut idx: usize,
    score: f32,
) -> Alignment {
    let mut slots = Vec::new();
    let mut in_delete = false;
    loop {
        let col = &cols[idx];
        let from = if in_delete {
            col.d_from[row]
        } else {
            col.m_from[row]
        };
        if !in_delete {
            slots.push((row - 1, Slot::Col(idx as u8)));
        }
        let consumed = usize::from(!in_delete);
        match from {
            From::Start => break,
            From::Match | From::Delete => {
                row -= consumed;
                idx -= 1;
                in_delete = matches!(from, From::Delete);
            }
            From::Loop => {
                let loop_id = profile
                    .loop_exit
                    .iter()
                    .position(|&e| e == idx)
                    .expect("loop entry column");
                let loop_end = row - consumed;
                let (src_row, src_delete) = col.entry_src[loop_end];
                slots.extend(
                    (src_row as usize..loop_end)
                        .rev()
                        .map(|r| (r, Slot::Loop(loop_id as u8))),
                );
                row = src_row as usize;
                idx -= 1;
                in_delete = src_delete;
            }
        }
    }
    slots.reverse();
    Alignment { score, slots }
}

pub(super) fn align(q: &[u8], profile: &Profile) -> Option<Alignment> {
    let cols = fill(q, profile);
    let (row, idx, score) = best_end(&cols, profile.loop_exit[2])?;
    Some(trace(&cols, profile, row, idx, score))
}
