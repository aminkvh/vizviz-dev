//! Local alignment of a sequence to a framework profile: match/delete
//! states with affine gaps inside the framework, and free-length CDR loops
//! between framework blocks priced by a length prior.

use super::profile::Profile;

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

/// State of a cell: matched, deleted as part of a real gap, or deleted as
/// a free column (one the germline leaves empty).
#[derive(Clone, Copy, PartialEq, Eq)]
enum State {
    Match,
    Delete,
    Free,
}

/// Predecessor of a cell, for traceback.
#[derive(Clone, Copy)]
enum From {
    Start,
    Cell(State),
    Loop,
}

/// One profile column's states, indexed by rows (residues consumed).
struct Column {
    m: Vec<f32>,
    d: Vec<f32>,
    f: Vec<f32>,
    m_from: Vec<From>,
    d_from: Vec<From>,
    f_from: Vec<From>,
    /// Loop entering this column: best score per row and its source
    /// `(row, state left behind)`.
    entry: Vec<f32>,
    entry_src: Vec<(u32, State)>,
}

impl Column {
    fn new(rows: usize) -> Self {
        Column {
            m: vec![NEG; rows],
            d: vec![NEG; rows],
            f: vec![NEG; rows],
            m_from: vec![From::Start; rows],
            d_from: vec![From::Start; rows],
            f_from: vec![From::Start; rows],
            entry: Vec::new(),
            entry_src: Vec::new(),
        }
    }

    /// Best state at `row`; ties keep the earlier of match, delete, free.
    fn best(&self, row: usize) -> (f32, State) {
        let mut best = (self.m[row], State::Match);
        for (score, state) in [(self.d[row], State::Delete), (self.f[row], State::Free)] {
            if score > best.0 {
                best = (score, state);
            }
        }
        best
    }

    fn origin(&self, state: State, row: usize) -> From {
        match state {
            State::Match => self.m_from[row],
            State::Delete => self.d_from[row],
            State::Free => self.f_from[row],
        }
    }
}

/// Fills `col.entry`: best score of leaving the previous block after
/// `row` residues, having spent the loop's residues in between.
fn fill_loop_entry(col: &mut Column, prev: &Column, profile: &Profile, loop_id: usize) {
    let rows = prev.m.len();
    let costs = &profile.loop_costs[loop_id];
    let leave: Vec<(f32, State)> = (0..rows).map(|r| prev.best(r)).collect();
    col.entry = vec![NEG; rows];
    col.entry_src = vec![(0, State::Match); rows];
    for row in 0..rows {
        let (mut best, mut src) = (NEG, (0, State::Match));
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

/// Match state of `row` in a column: a residue scored against it.
fn fill_match(q: &[u8], profile: &Profile, idx: usize, prev: Option<&Column>, col: &mut Column) {
    let scores = &profile.score[idx];
    let entering_loop = !col.entry.is_empty();
    for row in 1..col.m.len() {
        let (arrive, from) = match prev {
            None => (NEG, From::Start),
            Some(_) if entering_loop => (col.entry[row - 1], From::Loop),
            Some(p) => {
                let (v, state) = p.best(row - 1);
                (v, From::Cell(state))
            }
        };
        let free_start = if idx <= profile.c23 { 0.0 } else { NEG };
        let (base, from) = match free_start >= arrive {
            true => (free_start, From::Start),
            false => (arrive, from),
        };
        if base > NEG / 2.0 {
            col.m[row] = base + scores[usize::from(q[row - 1])];
            col.m_from[row] = from;
        }
    }
}

/// Keeps `(val, from)` in cell `row` of `scores` when it beats `NEG`.
fn store(scores: &mut [f32], froms: &mut [From], row: usize, (val, from): (f32, From)) {
    if val > NEG / 2.0 {
        scores[row] = val;
        froms[row] = from;
    }
}

/// Deleting a real column: open from a match or free state, extend a gap.
fn fill_delete(prev: &Column, col: &mut Column, row: usize, entering_loop: bool) {
    let open = |v: f32| v - GAP_OPEN;
    let cand = if entering_loop {
        (open(col.entry[row]), From::Loop)
    } else {
        let mut best = (open(prev.m[row]), From::Cell(State::Match));
        for (val, from) in [
            (open(prev.f[row]), From::Cell(State::Free)),
            (prev.d[row] - GAP_EXTEND, From::Cell(State::Delete)),
        ] {
            if val > best.0 {
                best = (val, from);
            }
        }
        best
    };
    store(&mut col.d, &mut col.d_from, row, cand);
}

/// Deleting a free column costs nothing: it continues a real gap through
/// it, or starts a free run that does not make a real gap after it cheaper.
fn fill_free(prev: &Column, col: &mut Column, row: usize, entering_loop: bool) {
    if entering_loop {
        store(
            &mut col.f,
            &mut col.f_from,
            row,
            (col.entry[row], From::Loop),
        );
        return;
    }
    store(
        &mut col.d,
        &mut col.d_from,
        row,
        (prev.d[row], From::Cell(State::Delete)),
    );
    let (val, from) = match prev.m[row] >= prev.f[row] {
        true => (prev.m[row], From::Cell(State::Match)),
        false => (prev.f[row], From::Cell(State::Free)),
    };
    store(&mut col.f, &mut col.f_from, row, (val, from));
}

fn fill_column(q: &[u8], profile: &Profile, idx: usize, prev: Option<&Column>, col: &mut Column) {
    fill_match(q, profile, idx, prev, col);
    let Some(p) = prev else { return };
    let entering_loop = !col.entry.is_empty();
    for row in 0..col.m.len() {
        match profile.free_gap[idx] {
            true => fill_free(p, col, row, entering_loop),
            false => fill_delete(p, col, row, entering_loop),
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
    let mut state = State::Match;
    loop {
        let col = &cols[idx];
        let from = col.origin(state, row);
        if state == State::Match {
            slots.push((row - 1, Slot::Col(idx as u8)));
        }
        let consumed = usize::from(state == State::Match);
        match from {
            From::Start => break,
            From::Cell(previous) => {
                row -= consumed;
                idx -= 1;
                state = previous;
            }
            From::Loop => {
                let loop_id = profile
                    .loop_exit
                    .iter()
                    .position(|&e| e == idx)
                    .expect("loop entry column");
                let loop_end = row - consumed;
                let (src_row, src_state) = col.entry_src[loop_end];
                slots.extend(
                    (src_row as usize..loop_end)
                        .rev()
                        .map(|r| (r, Slot::Loop(loop_id as u8))),
                );
                row = src_row as usize;
                idx -= 1;
                state = src_state;
            }
        }
    }
    slots.reverse();
    Alignment { score, slots }
}

/// Longest stretch before the first conserved Cys that a chain can still
/// number in full (kappa FR1 has 22 positions before it).
const MAX_BEFORE_C23: usize = 22;

impl Alignment {
    /// Gives the residues in front of a late-starting alignment the
    /// framework columns just before its first, when the chain has no
    /// room for anything but framework there (an N-terminal stub the local
    /// score left out) and the columns exist for all of them.
    pub(super) fn claim_leading(&mut self, profile: &Profile) {
        let Some(&(first, Slot::Col(col))) = self.slots.first() else {
            return;
        };
        let cys = self
            .slots
            .iter()
            .find(|s| s.1 == Slot::Col(profile.c23 as u8))
            .map(|s| s.0);
        if first == 0 || cys.is_none_or(|i| i > MAX_BEFORE_C23) {
            return;
        }
        let free: Vec<usize> = (0..usize::from(col))
            .filter(|&c| !profile.free_gap[c])
            .collect();
        let Some(cols) = free.len().checked_sub(first).map(|at| &free[at..]) else {
            return;
        };
        let claimed = cols
            .iter()
            .enumerate()
            .map(|(i, &c)| (i, Slot::Col(c as u8)));
        self.slots.splice(0..0, claimed.collect::<Vec<_>>());
    }
}

pub(super) fn align(q: &[u8], profile: &Profile) -> Option<Alignment> {
    let cols = fill(q, profile);
    let (row, idx, score) = best_end(&cols, profile.loop_exit[2])?;
    Some(trace(&cols, profile, row, idx, score))
}
