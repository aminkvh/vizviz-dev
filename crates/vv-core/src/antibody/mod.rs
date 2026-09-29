//! Antibody variable-domain detection and numbering from sequence alone.
//!
//! A domain is found by aligning the chain to framework profiles (heavy,
//! kappa, lambda) built from public PDB sequences, with the CDRs treated
//! as free-length loops. The IMGT frame is canonical; Kabat, Chothia and
//! Martin numbers are derived from it by length rules. Method, data
//! provenance, licences and measured accuracy: `docs/ANTIBODY.md`.

mod align;
mod cdr;
mod numbering;
mod profile;
mod seeds;

pub use cdr::{CdrDefinition, Region};
pub use numbering::{Label, Scheme};

use align::{align, Alignment, Slot};
use profile::{aa_index, profiles, Profile, C104, C23, W41};

/// Variable-domain family of a chain.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum ChainType {
    Heavy,
    Kappa,
    Lambda,
}

/// A residue of a domain with its number and region.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Annotation {
    /// Index into the sequence given to [`find_domains`].
    pub index: usize,
    pub label: Label,
    pub region: Region,
}

/// A detected variable domain.
#[derive(Clone, Debug)]
pub struct Domain {
    pub chain: ChainType,
    /// Residue range `start..end` of the domain in the input sequence.
    pub start: usize,
    pub end: usize,
    /// Alignment score in bits.
    pub score: f32,
    /// Score relative to a perfect framework match, `0..=1`.
    pub confidence: f32,
    imgt: Vec<Label>,
}

/// Minimum `confidence` for a domain to be reported.
const MIN_CONFIDENCE: f32 = 0.30;
/// Shortest stretch worth aligning (a V domain is at least ~90 residues).
const MIN_STRETCH: usize = 80;
/// Distance range between the two conserved Cys of a V domain.
const CYS_SPAN: std::ops::RangeInclusive<usize> = 40..=110;
/// Residues kept before the first Cys (FR1) and after the last (CDR3+FR4).
const WINDOW_BEFORE: usize = 30;
const WINDOW_AFTER: usize = 70;
/// How far past the second Cys the FR4 `[WF]G.G` motif may start.
const FR4_REACH: usize = 55;

impl Domain {
    /// Label of each domain residue in `scheme`, as `(sequence index,
    /// label)` in sequence order. IMGT insertion letters stand for `.1`,
    /// `.2`, ... (`112A` is IMGT 112.1).
    pub fn numbering(&self, scheme: Scheme) -> Vec<(usize, Label)> {
        let labels = numbering::relabel(scheme, self.chain, &self.imgt);
        (self.start..self.end).zip(labels).collect()
    }

    /// Label in `scheme` and region under `definition` for each residue.
    pub fn annotate(&self, scheme: Scheme, definition: CdrDefinition) -> Vec<Annotation> {
        let native = numbering::relabel(definition.native_scheme(), self.chain, &self.imgt);
        let regions = definition.regions(self.chain, &native);
        let shown = numbering::relabel(scheme, self.chain, &self.imgt);
        (self.start..self.end)
            .zip(shown)
            .zip(regions)
            .map(|((index, label), region)| Annotation {
                index,
                label,
                region,
            })
            .collect()
    }
}

/// Loose FR4 test: `[WF]G`, `[WF]..G` or `.G.G` (germline `[WF]G.G` plus
/// the common variants `WAxG`, `RGQG`, `FSGG`).
fn has_fr4_motif(q: &[u8], from: usize) -> bool {
    let is = |i: usize, letters: &[u8]| {
        q.get(i)
            .is_some_and(|&a| letters.iter().any(|&l| a == aa_index(l) as u8))
    };
    (from..(from + FR4_REACH).min(q.len())).any(|i| {
        (is(i, b"WF") && (is(i + 1, b"G") || is(i + 3, b"G")))
            || (is(i + 1, b"G") && is(i + 3, b"G"))
    })
}

/// Span of `q` that can hold a domain: from before the first Cys that has
/// a partner one V-domain length away (and an FR4-like motif after it) to
/// after the last such partner.
fn candidate_window(q: &[u8]) -> Option<(usize, usize)> {
    let cys = aa_index(b'C') as u8;
    let sites: Vec<usize> = (0..q.len()).filter(|&i| q[i] == cys).collect();
    let pairs: Vec<(usize, usize)> = sites
        .iter()
        .flat_map(|&i| sites.iter().map(move |&j| (i, j)))
        .filter(|&(i, j)| j > i && CYS_SPAN.contains(&(j - i)) && has_fr4_motif(q, j + 3))
        .collect();
    let first = pairs.iter().map(|p| p.0).min()?;
    let last = pairs.iter().map(|p| p.1).max()?;
    Some((
        first.saturating_sub(WINDOW_BEFORE),
        (last + WINDOW_AFTER).min(q.len()),
    ))
}

fn anchors_hold(q: &[u8], slots: &[(usize, Slot)]) -> bool {
    let at = |col: usize| {
        slots
            .iter()
            .find(|(_, s)| *s == Slot::Col(col as u8))
            .map(|(i, _)| q[*i])
    };
    let is = |col, letter: u8| at(col) == Some(aa_index(letter) as u8);
    is(C23, b'C') && is(C104, b'C') && is(W41, b'W')
}

fn best_alignment(q: &[u8]) -> Option<(&'static Profile, Alignment)> {
    profiles()
        .iter()
        .filter_map(|p| align(q, p).map(|a| (p, a)))
        .max_by(|a, b| a.1.score.total_cmp(&b.1.score))
}

fn scan(q: &[u8], lo: usize, hi: usize, min_confidence: f32, out: &mut Vec<Domain>) {
    if hi - lo < MIN_STRETCH {
        return;
    }
    let Some((from, to)) = candidate_window(&q[lo..hi]) else {
        return;
    };
    let (from, to) = (lo + from, lo + to);
    let Some((profile, hit)) = best_alignment(&q[from..to]) else {
        return;
    };
    let confidence = hit.score / profile.ideal;
    if confidence < min_confidence || !anchors_hold(&q[from..to], &hit.slots) {
        return;
    }
    let start = from + hit.slots[0].0;
    let end = from + hit.slots.last().map_or(0, |s| s.0) + 1;
    let slots: Vec<Slot> = hit.slots.iter().map(|s| s.1).collect();
    out.push(Domain {
        chain: profile.chain,
        start,
        end,
        score: hit.score,
        confidence,
        imgt: numbering::imgt_labels(&slots),
    });
    scan(q, lo, start, min_confidence, out);
    scan(q, end, hi, min_confidence, out);
}

/// Finds every antibody variable domain (heavy, kappa, lambda) in a
/// one-letter protein sequence, in sequence order. Non-antibody chains
/// return an empty list.
pub fn find_domains(seq: &str) -> Vec<Domain> {
    find_domains_with(seq, MIN_CONFIDENCE)
}

/// [`find_domains`] with a caller-chosen minimum [`Domain::confidence`].
pub fn find_domains_with(seq: &str, min_confidence: f32) -> Vec<Domain> {
    let q: Vec<u8> = seq.bytes().map(|b| aa_index(b) as u8).collect();
    let mut out = Vec::new();
    scan(&q, 0, q.len(), min_confidence, &mut out);
    out.sort_by_key(|d| d.start);
    out
}
