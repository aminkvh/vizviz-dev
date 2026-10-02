//! The antibody track: the numbering scheme's labels over framework and
//! CDR bars, for every variable domain found in a protein chain.

use rayon::prelude::*;
use vv_core::antibody::external::ExternalDomain;
use vv_core::antibody::{find_in_residues, Annotation, ChainType, Domain, Region};

use super::hex;
use crate::sequence::anarci::Backend;
use crate::sequence::tracks::{
    legend, AntibodySettings, External, Glyph, Inputs, TrackContext, TrackData, TrackProvider,
};

pub struct Antibody;

const FRAMEWORK: u8 = 1;
/// Fewest residues between two tick labels, so "100A" never overlaps its
/// neighbour.
const TICK_SPACING: usize = 3;

fn kind_of(region: Region) -> u8 {
    FRAMEWORK + region.cdr().unwrap_or(0)
}

/// Decade numbers and every insertion-coded position, thinned so labels
/// keep `TICK_SPACING` residues apart.
fn tick_indices(notes: &[Annotation]) -> Vec<usize> {
    let mut last: Option<usize> = None;
    let mut out = Vec::new();
    for (i, a) in notes.iter().enumerate() {
        let wanted = a.label.insertion().is_some() || a.label.number % 10 == 0;
        if wanted && last.is_none_or(|l| i >= l + TICK_SPACING) {
            out.push(i);
            last = Some(i);
        }
    }
    out
}

/// A numbered domain as the track paints it, whoever numbered it.
struct Painted {
    chain: ChainType,
    start: usize,
    notes: Vec<Annotation>,
}

impl Painted {
    fn native(domain: &Domain, s: AntibodySettings) -> Self {
        Self {
            chain: domain.chain,
            start: domain.start,
            notes: domain.annotate(s.scheme, s.cdr),
        }
    }

    /// `None` when ANARCI was not run for the scheme or definition.
    fn external(domain: &ExternalDomain, s: AntibodySettings) -> Option<Self> {
        Some(Self {
            chain: domain.chain,
            start: domain.start,
            notes: domain.annotate(s.scheme, s.cdr)?,
        })
    }
}

fn tooltip(chain: ChainType, note: &Annotation, s: AntibodySettings) -> String {
    let letter = chain.letter();
    let number = match s.scheme == s.cdr.native_scheme() {
        true => format!("{letter}{}", note.label),
        false => format!("{letter}{} ({} numbering)", note.label, s.scheme.name()),
    };
    let region = match note.region.cdr() {
        Some(n) => format!("CDR-{letter}{n}"),
        None => note.region.name().to_string(),
    };
    format!("{number} \u{B7} {region} ({})", s.cdr.name())
}

fn paint_domain(track: &mut TrackData, domain: &Painted, first: u32, s: AntibodySettings) {
    let notes = &domain.notes;
    for a in notes {
        let r = first + a.index as u32;
        track.mark(r, kind_of(a.region));
        track.note(r, tooltip(domain.chain, a, s));
    }
    for i in tick_indices(notes) {
        track.set_tick(first + notes[i].index as u32, notes[i].label.to_string());
    }
    track.set_badge(first + domain.start as u32, domain.chain.name());
}

/// Domains per chain row: ANARCI's when it is selected and has answered,
/// native otherwise (the track stays empty only while ANARCI runs).
fn painted_domains(ctx: &TrackContext) -> Vec<(u32, Vec<Painted>)> {
    let s = ctx.antibody;
    match (s.backend, ctx.extras.anarci) {
        (Backend::Anarci, External::Pending) => Vec::new(),
        (Backend::Anarci, External::Done(Ok(chains))) => ctx
            .rows
            .iter()
            .zip(chains)
            .map(|((_, rows), domains)| {
                let painted = domains.iter().filter_map(|d| Painted::external(d, s));
                (rows.start, painted.collect())
            })
            .collect(),
        _ => {
            let top = ctx.top();
            ctx.rows
                .par_iter()
                .map(|(_, rows)| {
                    let found = find_in_residues(top, rows.clone());
                    (
                        rows.start,
                        found.iter().map(|d| Painted::native(d, s)).collect(),
                    )
                })
                .collect()
        }
    }
}

impl TrackProvider for Antibody {
    fn id(&self) -> &'static str {
        "antibody"
    }

    fn label(&self) -> &'static str {
        "Antibody"
    }

    fn uses_antibody_settings(&self) -> bool {
        true
    }

    fn inputs(&self) -> Inputs {
        Inputs {
            anarci: true,
            ..Inputs::default()
        }
    }

    fn compute(&self, ctx: &TrackContext) -> TrackData {
        let found = painted_domains(ctx);
        let mut track = ctx.new_track(
            Glyph::LabeledBar,
            vec![
                legend("Framework", hex(0x9AA5B5)),
                legend("CDR1", hex(0xE8A33D)),
                legend("CDR2", hex(0x3FA796)),
                legend("CDR3", hex(0xC8628F)),
            ],
        );
        for (first, domains) in &found {
            for domain in domains {
                paint_domain(&mut track, domain, *first, ctx.antibody);
            }
        }
        track.finish()
    }
}
