//! Tracks about the chain itself: secondary structure, numbering, and the
//! residues the experiment did not place.

use std::ops::Range;

use vv_core::dssp::DsspCode;
use vv_core::seqfeat::{unobserved, unobserved_in_entity, Gap};
use vv_render::color::by_secondary_structure;

use super::hex;
use crate::sequence::tracks::{legend, Glyph, Inputs, TrackContext, TrackData, TrackProvider};

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
        .map(|u| match u.number {
            Some((seq, 0)) => format!("{} {seq}", u.name),
            Some((seq, ins)) => format!("{} {seq}{}", u.name, ins as char),
            None => format!("{} #{}", u.name, u.position),
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

    /// Only a file that lists unobserved residues for several models has
    /// a gap list that changes with the shown model.
    fn frame_key(&self, loaded: &vv_scene::LoadedStructure, frame: usize) -> usize {
        if lists_several_models(&loaded.structure.topology) {
            frame
        } else {
            0
        }
    }

    fn inputs(&self) -> Inputs {
        Inputs {
            entity: true,
            ..Inputs::default()
        }
    }

    /// Residues of the chain's full sequence that no modeled residue
    /// aligns to, when the file gives the sequence; else those the file
    /// lists, placed by author number.
    fn compute(&self, ctx: &TrackContext) -> TrackData {
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
        for (name, residues) in ctx.rows {
            for gap in gaps_of(ctx, name, residues) {
                let after = gap.before >= residues.end;
                let at = if after { residues.end - 1 } else { gap.before };
                track.mark(at, if after { 2 } else { 1 });
                track.note(at, gap_note(&gap, after));
            }
        }
        track.finish()
    }
}

fn gaps_of(ctx: &TrackContext, name: &str, residues: &Range<u32>) -> Vec<Gap> {
    let entity = ctx.extras.entity.and_then(|e| e.get(name));
    match entity {
        Some(entity) => unobserved_in_entity(ctx.top(), residues.clone(), entity, ctx.extras.model),
        None => unobserved(ctx.top(), residues.clone()),
    }
}

/// Whether the unobserved-residue table has rows for more than one model
/// (its rows are grouped by model).
fn lists_several_models(top: &vv_core::Topology) -> bool {
    let Some(cat) = top.annotations.category("pdbx_unobs_or_zero_occ_residues") else {
        return false;
    };
    let last = cat.rows.len().saturating_sub(1);
    cat.get("PDB_model_num", 0) != cat.get("PDB_model_num", last)
}
