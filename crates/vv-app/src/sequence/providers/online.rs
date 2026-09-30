//! The UniProt feature track: annotations fetched for the chains' UniProt
//! entries and drawn where SIFTS places them on the chain.

use std::collections::HashMap;

use super::hex;
use crate::sequence::rows::is_protein;
use crate::sequence::tracks::{legend, Glyph, Inputs, TrackContext, TrackData, TrackProvider};
use crate::sequence::uniprot::{map_residues, Data, Feature, Kind};

pub struct Uniprot;

/// Tooltip lines kept per residue; more are summarized.
const NOTES_PER_RESIDUE: u32 = 6;

fn kind_color(kind: Kind) -> u32 {
    hex(match kind {
        Kind::Region => 0xB4C4DC,
        Kind::Topology => 0x3AA88A,
        Kind::Variant => 0xA0A4AB,
        Kind::OtherSite => 0x8E7CC3,
        Kind::Modification => 0xE39B34,
        Kind::Site => 0xD1495B,
    })
}

/// Residue numbers a feature covers: a disulfide or cross-link joins its
/// two ends and covers nothing between them.
fn covered(f: &Feature) -> Vec<u32> {
    let joins_ends = matches!(f.label.as_str(), "Disulfide bond" | "Cross-link");
    if joins_ends && f.start != f.end {
        vec![f.start, f.end]
    } else {
        (f.start..=f.end).collect()
    }
}

fn note(f: &Feature, entry: &str) -> String {
    let span = if f.start == f.end {
        f.start.to_string()
    } else {
        format!("{}-{}", f.start, f.end)
    };
    if f.text.is_empty() {
        format!("{} {span} ({entry})", f.label)
    } else {
        format!("{} {span} ({entry}): {}", f.label, f.text)
    }
}

struct Painter<'a> {
    track: &'a mut TrackData,
    kinds: &'static [Kind],
    notes: HashMap<u32, u32>,
}

impl Painter<'_> {
    fn add(&mut self, residue: u32, kind: Kind, line: &str) {
        let Some(at) = self.kinds.iter().position(|&k| k == kind) else {
            return;
        };
        self.track.mark(residue, 1 + at as u8);
        let count = self.notes.entry(residue).or_default();
        *count += 1;
        match *count {
            n if n <= NOTES_PER_RESIDUE => self.track.note(residue, line),
            n if n == NOTES_PER_RESIDUE + 1 => self.track.note(residue, "and more"),
            _ => {}
        }
    }
}

/// What the feature track draws; variants are their own track, as there
/// are hundreds of them.
const FEATURE_KINDS: &[Kind] = &[
    Kind::Region,
    Kind::Topology,
    Kind::OtherSite,
    Kind::Modification,
    Kind::Site,
];

pub struct Variants;

impl TrackProvider for Variants {
    fn id(&self) -> &'static str {
        "variants"
    }

    fn label(&self) -> &'static str {
        "Variants"
    }

    fn online(&self) -> bool {
        true
    }

    fn inputs(&self) -> Inputs {
        Inputs {
            uniprot: true,
            ..Inputs::default()
        }
    }

    /// UniProt's natural variants; hover lists the substitution and its note.
    fn compute(&self, ctx: &TrackContext) -> TrackData {
        compute_kinds(ctx, &[Kind::Variant])
    }
}

fn compute_kinds(ctx: &TrackContext, kinds: &'static [Kind]) -> TrackData {
    let mut track = ctx.new_track(
        Glyph::Bar,
        kinds
            .iter()
            .map(|&k| legend(k.legend(), kind_color(k)))
            .collect(),
    );
    if let Some(data) = ctx.extras.uniprot {
        let mut painter = Painter {
            track: &mut track,
            kinds,
            notes: HashMap::new(),
        };
        for (name, residues) in ctx.rows {
            if !residues.is_empty() && is_protein(ctx.top(), residues.start) {
                paint_chain(&mut painter, ctx, data, name, residues.clone());
            }
        }
    }
    track.finish()
}

impl TrackProvider for Uniprot {
    fn id(&self) -> &'static str {
        "uniprot"
    }

    fn label(&self) -> &'static str {
        "UniProt"
    }

    fn online(&self) -> bool {
        true
    }

    fn inputs(&self) -> Inputs {
        Inputs {
            uniprot: true,
            ..Inputs::default()
        }
    }

    /// Features of the chain's UniProt entry: later kinds are painted over
    /// earlier ones (sites over modifications over regions), and every
    /// feature at a residue is listed in its tooltip.
    fn compute(&self, ctx: &TrackContext) -> TrackData {
        compute_kinds(ctx, FEATURE_KINDS)
    }
}

fn paint_chain(
    painter: &mut Painter,
    ctx: &TrackContext,
    data: &Data,
    name: &str,
    residues: std::ops::Range<u32>,
) {
    let mapped = map_residues(data, name, ctx.top(), residues);
    let at: HashMap<(&str, u32), u32> = mapped
        .iter()
        .map(|m| ((m.accession.as_str(), m.position), m.residue))
        .collect();
    let entries: HashMap<&str, &str> = data
        .segments
        .iter()
        .map(|s| (s.accession.as_str(), s.name.as_str()))
        .collect();
    let mut features: Vec<(&str, &Feature)> = data
        .features
        .iter()
        .filter(|(acc, _)| mapped.iter().any(|m| m.accession == **acc))
        .flat_map(|(acc, list)| list.iter().map(move |f| (acc.as_str(), f)))
        .collect();
    features.sort_by_key(|(acc, f)| (f.kind, *acc, f.start));
    for (acc, f) in features {
        let line = note(f, entries.get(acc).copied().unwrap_or(acc));
        for position in covered(f) {
            if let Some(&residue) = at.get(&(acc, position)) {
                painter.add(residue, f.kind, &line);
            }
        }
    }
}
