//! Annotation tracks: thin rows under a chain's sequence, one mark per
//! residue. A track is a [`TrackProvider`]; its result, [`TrackData`], is
//! computed once per structure (again per frame for geometry-dependent
//! ones) and cached (`cache.rs`).

use std::collections::BTreeMap;
use std::ops::Range;

use vv_core::glam::Vec3;
use vv_scene::LoadedStructure;

/// How a track's marks are painted.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Glyph {
    /// A filled bar across the residue, colored by mark kind.
    Bar,
    /// Secondary structure: kind 1 helix (thick bar), 2 strand (arrow),
    /// 3 turn (thin line).
    Structure,
    /// The residue's number as text.
    Ticks,
    /// A marker on the residue's edge: kind 1 left, kind 2 right.
    Edge,
}

pub struct LegendEntry {
    pub label: String,
    pub color: u32,
}

/// Marks for every residue of one structure. `kind` 0 is "no mark";
/// kind `k` paints with `colors[k - 1]`.
pub struct TrackData {
    pub glyph: Glyph,
    pub colors: Vec<u32>,
    pub legend: Vec<LegendEntry>,
    kind: Vec<u8>,
    notes: BTreeMap<u32, String>,
    marked_before: Vec<u32>,
}

impl TrackData {
    pub fn new(glyph: Glyph, residues: usize, legend: Vec<LegendEntry>) -> Self {
        let colors = legend.iter().map(|l| l.color).collect();
        Self {
            glyph,
            colors,
            legend,
            kind: vec![0; residues],
            notes: BTreeMap::new(),
            marked_before: Vec::new(),
        }
    }

    /// Paints residue `r` with `kind`, replacing any earlier mark.
    pub fn mark(&mut self, r: u32, kind: u8) {
        self.kind[r as usize] = kind;
    }

    /// Adds a line to residue `r`'s tooltip.
    pub fn note(&mut self, r: u32, line: impl AsRef<str>) {
        let entry = self.notes.entry(r).or_default();
        if !entry.is_empty() {
            entry.push('\n');
        }
        entry.push_str(line.as_ref());
    }

    /// Freezes the track and indexes it for [`TrackData::any_in`].
    pub fn finish(mut self) -> Self {
        let mut running = 0;
        self.marked_before = std::iter::once(0)
            .chain(self.kind.iter().map(|&k| {
                running += u32::from(k != 0);
                running
            }))
            .collect();
        self
    }

    pub fn kind(&self, r: u32) -> u8 {
        self.kind.get(r as usize).copied().unwrap_or(0)
    }

    pub fn color(&self, kind: u8) -> u32 {
        self.colors[kind as usize - 1]
    }

    /// What hovering residue `r` says about this track: its notes, or the
    /// legend label of its mark.
    pub fn describe(&self, r: u32) -> Option<String> {
        let kind = self.kind(r);
        if kind == 0 {
            return None;
        }
        match self.notes.get(&r) {
            Some(text) => Some(text.clone()),
            None => self.legend.get(kind as usize - 1).map(|l| l.label.clone()),
        }
    }

    /// Whether any residue in `range` is marked.
    pub fn any_in(&self, range: &Range<u32>) -> bool {
        let at = |i: u32| self.marked_before[i as usize];
        at(range.end) > at(range.start)
    }
}

/// Everything a provider may read about one structure at one frame.
pub struct TrackContext<'a> {
    pub loaded: &'a LoadedStructure,
    pub positions: &'a [Vec3],
    /// The strip's rows for this structure: chain name and residues.
    pub rows: &'a [(String, Range<u32>)],
}

impl TrackContext<'_> {
    pub fn top(&self) -> &vv_core::Topology {
        &self.loaded.structure.topology
    }

    pub fn new_track(&self, glyph: Glyph, legend: Vec<LegendEntry>) -> TrackData {
        TrackData::new(glyph, self.top().residue_count(), legend)
    }
}

/// One annotation track.
///
/// Hook for further tracks: implement this and add the unit struct to
/// [`PROVIDERS`]. The track appears in the header's Tracks menu and answers
/// `sequence track <id> on|off` with no other change. An antibody
/// numbering / CDR track registers here: it reads `ctx.rows` and
/// `ctx.top()` and returns a `TrackData` (`Glyph::Ticks` for a scheme's
/// numbers, `Glyph::Bar` with a legend for CDRs and frameworks).
pub trait TrackProvider: Sync {
    /// The name `sequence track` takes.
    fn id(&self) -> &'static str;
    /// The name in the Tracks menu and beside the track's row.
    fn label(&self) -> &'static str;
    /// Whether the result changes with the frame (geometry-dependent).
    fn per_frame(&self) -> bool {
        false
    }
    fn compute(&self, ctx: &TrackContext) -> TrackData;
}

pub fn provider(id: &str) -> Option<&'static dyn TrackProvider> {
    PROVIDERS.iter().copied().find(|p| p.id() == id)
}

/// The registry of tracks, in display order.
pub const PROVIDERS: &[&dyn TrackProvider] = &[
    &super::providers::SecondaryStructure,
    &super::providers::Numbering,
    &super::providers::Missing,
    &super::providers::Disulfides,
    &super::providers::Glycans,
    &super::providers::Liabilities,
    &super::providers::LigandSite,
    &super::providers::Interface,
    &super::providers::AltLocs,
    &super::providers::Modified,
];

pub fn legend(label: &str, color: u32) -> LegendEntry {
    LegendEntry {
        label: label.to_string(),
        color,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn any_in_counts_marks_inside_a_range_only() {
        let mut t = TrackData::new(Glyph::Bar, 6, vec![legend("x", 0xFF00_0000)]);
        t.mark(1, 1);
        t.mark(4, 1);
        let t = t.finish();
        assert!(t.any_in(&(0..2)));
        assert!(!t.any_in(&(2..4)));
        assert!(t.any_in(&(2..6)));
        assert!(!t.any_in(&(5..6)));
        assert!(!t.any_in(&(0..0)));
    }

    #[test]
    fn describe_prefers_notes_over_the_legend_label() {
        let mut t = TrackData::new(Glyph::Bar, 3, vec![legend("kind", 0)]);
        t.mark(0, 1);
        t.mark(1, 1);
        t.note(1, "first");
        t.note(1, "second");
        assert_eq!(t.describe(0).as_deref(), Some("kind"));
        assert_eq!(t.describe(1).as_deref(), Some("first\nsecond"));
        assert_eq!(t.describe(2), None);
    }
}
