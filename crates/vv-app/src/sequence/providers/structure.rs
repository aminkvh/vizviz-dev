//! Tracks about the chain itself: secondary structure, numbering, and the
//! residues the experiment did not place.

use vv_core::dssp::DsspCode;
use vv_render::color::by_secondary_structure;

use super::hex;
use crate::sequence::tracks::{legend, Glyph, TrackContext, TrackData, TrackProvider};

pub struct SecondaryStructure;

impl TrackProvider for SecondaryStructure {
    fn id(&self) -> &'static str {
        "ss"
    }

    fn label(&self) -> &'static str {
        "Secondary"
    }

    fn per_frame(&self) -> bool {
        true
    }

    /// The assignment the cartoon draws (file records refined by DSSP).
    fn compute(&self, ctx: &TrackContext) -> TrackData {
        let single = ctx.loaded.structure.frame_count() == 1;
        let codes = vv_core::cartoon::secondary_structure(ctx.top(), ctx.positions, single);
        let mut track = ctx.new_track(
            Glyph::Structure,
            vec![
                legend("Helix", by_secondary_structure(DsspCode::AlphaHelix)),
                legend("Strand", by_secondary_structure(DsspCode::Strand)),
                legend("Turn", by_secondary_structure(DsspCode::Turn)),
            ],
        );
        for (r, code) in codes.iter().enumerate() {
            match code {
                DsspCode::AlphaHelix | DsspCode::Helix3_10 | DsspCode::HelixPi => {
                    track.mark(r as u32, 1)
                }
                DsspCode::Strand | DsspCode::Bridge => track.mark(r as u32, 2),
                DsspCode::Turn => track.mark(r as u32, 3),
                _ => {}
            }
        }
        track.finish()
    }
}

pub struct Numbering;

impl Numbering {
    /// Whether `residue` gets a number: every tenth, plus the first of a
    /// row unless a tenth follows too closely to fit its label.
    fn ticks(top: &vv_core::Topology, residues: &std::ops::Range<u32>) -> Vec<u32> {
        let is_tenth = |r: u32| {
            let rec = &top.residues[r as usize];
            rec.auth_seq_id % 10 == 0 && rec.ins_code == 0
        };
        let mut ticks: Vec<u32> = residues.clone().filter(|&r| is_tenth(r)).collect();
        let first = residues.start;
        if ticks.first().is_none_or(|&t| t >= first + 3) {
            ticks.insert(0, first);
        }
        ticks
    }
}

impl TrackProvider for Numbering {
    fn id(&self) -> &'static str {
        "numbering"
    }

    fn label(&self) -> &'static str {
        "Numbering"
    }

    fn compute(&self, ctx: &TrackContext) -> TrackData {
        let mut track = ctx.new_track(
            Glyph::Ticks,
            vec![legend("Author numbering", hex(0x808080))],
        );
        for (_, residues) in ctx.rows {
            if residues.is_empty() {
                continue;
            }
            for r in Self::ticks(ctx.top(), residues) {
                track.mark(r, 1);
            }
        }
        track.finish()
    }
}

pub struct Missing;

const LISTED_NAMES: usize = 8;

fn gap_note(gap: &vv_core::seqfeat::Gap, after: bool) -> String {
    let shown: Vec<String> = gap
        .residues
        .iter()
        .take(LISTED_NAMES)
        .map(|u| match u.ins {
            0 => format!("{} {}", u.name, u.seq),
            c => format!("{} {}{}", u.name, u.seq, c as char),
        })
        .collect();
    let more = gap.residues.len().saturating_sub(LISTED_NAMES);
    let tail = if more > 0 {
        format!(" (+{more} more)")
    } else {
        String::new()
    };
    format!(
        "{} residue(s) not modeled {} this one: {}{tail}",
        gap.residues.len(),
        if after { "after" } else { "before" },
        shown.join(", ")
    )
}

impl TrackProvider for Missing {
    fn id(&self) -> &'static str {
        "missing"
    }

    fn label(&self) -> &'static str {
        "Missing"
    }

    fn compute(&self, ctx: &TrackContext) -> TrackData {
        let top = ctx.top();
        let mut track = ctx.new_track(
            Glyph::Edge,
            vec![
                legend(
                    "Residues in the sequence with no coordinates (gap before)",
                    hex(0xC0392B),
                ),
                legend("Residues with no coordinates (gap after)", hex(0xC0392B)),
            ],
        );
        for (_, residues) in ctx.rows {
            for gap in vv_core::seqfeat::unobserved(top, residues.clone()) {
                let after = gap.before >= residues.end;
                let at = if after { residues.end - 1 } else { gap.before };
                track.mark(at, if after { 2 } else { 1 });
                track.note(at, gap_note(&gap, after));
            }
        }
        track.finish()
    }
}
