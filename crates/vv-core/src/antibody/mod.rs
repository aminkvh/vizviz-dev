//! Antibody variable-domain detection and numbering from sequence alone.
//!
//! A domain is found by aligning the chain to framework profiles (heavy,
//! kappa, lambda, and T-cell receptor alpha and beta) built from public PDB
//! sequences, with the CDRs treated as free-length loops. The IMGT frame is
//! canonical; Kabat, Chothia, Martin and AHo numbers are derived from it by
//! length rules. Method, data provenance, licences and measured accuracy:
//! `docs/ANTIBODY.md`.

mod aho;
mod align;
mod cdr;
mod numbering;
mod profile;
mod seeds;
mod tcr_seeds;
mod topology;

pub use cdr::{CdrDefinition, Region};
pub use numbering::{Label, Scheme};
pub use topology::{cdr_residues, chain_domains, find_in_residues, AntibodyCache, CdrResidue};

use align::{align, Alignment, Slot};
use profile::{aa_index, profiles, Profile, C23, W41};

/// Variable-domain family of a chain: an antibody heavy or light chain, or
/// a T-cell receptor alpha or beta chain.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum ChainType {
    Heavy,
    Kappa,
    Lambda,
    TcrAlpha,
    TcrBeta,
}

impl ChainType {
    pub fn name(self) -> &'static str {
        match self {
            ChainType::Heavy => "Heavy",
            ChainType::Kappa => "Kappa",
            ChainType::Lambda => "Lambda",
            ChainType::TcrAlpha => "TCR alpha",
            ChainType::TcrBeta => "TCR beta",
        }
    }

    /// `H` for heavy, `L` for both light chains, `A` and `B` for the
    /// receptor chains.
    pub fn letter(self) -> char {
        match self {
            ChainType::Heavy => 'H',
            ChainType::Kappa | ChainType::Lambda => 'L',
            ChainType::TcrAlpha => 'A',
            ChainType::TcrBeta => 'B',
        }
    }

    pub fn is_antibody(self) -> bool {
        !matches!(self, ChainType::TcrAlpha | ChainType::TcrBeta)
    }
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
/// Antibody confidence above which receptor profiles are not consulted.
const ANTIBODY_SURE: f32 = 0.5;
/// Receptor profiles rest on fewer seed chains, so genuine receptors of a
/// gene family the seeds miss score lower; nothing else reaches that far.
const RECEPTOR_FLOOR_RATIO: f32 = 0.5;
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
    /// `.2`, ... (`112A` is IMGT 112.1). T-cell receptor domains are
    /// numbered in IMGT whatever `scheme` says.
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

fn anchors_hold(q: &[u8], slots: &[(usize, Slot)], c104: usize) -> bool {
    let at = |col: usize| {
        slots
            .iter()
            .find(|(_, s)| *s == Slot::Col(col as u8))
            .map(|(i, _)| q[*i])
    };
    let is = |col, letter: u8| at(col) == Some(aa_index(letter) as u8);
    is(C23, b'C') && is(c104, b'C') && is(W41, b'W')
}

/// Best-scoring profile among those selected by `wanted`.
fn best_of(q: &[u8], wanted: impl Fn(&Profile) -> bool) -> Option<(&'static Profile, Alignment)> {
    profiles()
        .iter()
        .filter(|p| wanted(p))
        .filter_map(|p| align(q, p).map(|a| (p, a)))
        .max_by(|a, b| a.1.score.total_cmp(&b.1.score))
}

/// Best profile overall. Receptor profiles are tried only when no antibody
/// profile fits convincingly: real antibodies clear `ANTIBODY_SURE` and
/// score near zero against receptors, so this saves their alignments.
fn best_alignment(q: &[u8]) -> Option<(&'static Profile, Alignment)> {
    let antibody = best_of(q, |p| p.chain.is_antibody());
    if antibody
        .as_ref()
        .is_some_and(|(p, a)| a.score / p.ideal >= ANTIBODY_SURE)
    {
        return antibody;
    }
    let receptor = best_of(q, |p| !p.chain.is_antibody());
    match (antibody, receptor) {
        (Some(a), Some(r)) => Some(if r.1.score > a.1.score { r } else { a }),
        (a, r) => a.or(r),
    }
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
    let floor = match profile.chain.is_antibody() {
        true => min_confidence,
        false => min_confidence * RECEPTOR_FLOOR_RATIO,
    };
    if confidence < floor || !anchors_hold(&q[from..to], &hit.slots, profile.c104) {
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
        imgt: numbering::imgt_labels(&slots, &profile.labels),
    });
    scan(q, lo, start, min_confidence, out);
    scan(q, end, hi, min_confidence, out);
}

/// Finds every antibody variable domain (heavy, kappa, lambda) in a
/// one-letter protein sequence, in sequence order. Non-antibody chains,
/// T-cell receptors included, return an empty list.
pub fn find_domains(seq: &str) -> Vec<Domain> {
    find_domains_with(seq, MIN_CONFIDENCE)
}

/// [`find_domains`] with a caller-chosen minimum [`Domain::confidence`].
pub fn find_domains_with(seq: &str, min_confidence: f32) -> Vec<Domain> {
    let mut all = find_variable_domains_with(seq, min_confidence);
    all.retain(|d| d.chain.is_antibody());
    all
}

/// Like [`find_domains`], but also reports T-cell receptor alpha and beta
/// variable domains. Each domain goes to the profile it fits best, so a
/// receptor is never reported as an antibody or the other way round.
pub fn find_variable_domains(seq: &str) -> Vec<Domain> {
    find_variable_domains_with(seq, MIN_CONFIDENCE)
}

/// [`find_variable_domains`] with a caller-chosen minimum confidence.
pub fn find_variable_domains_with(seq: &str, min_confidence: f32) -> Vec<Domain> {
    let q: Vec<u8> = seq.bytes().map(|b| aa_index(b) as u8).collect();
    let mut out = Vec::new();
    scan(&q, 0, q.len(), min_confidence, &mut out);
    out.sort_by_key(|d| d.start);
    out
}
