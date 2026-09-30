//! Tracks computed from the structure's geometry or from comparing chains:
//! burial and conservation.

use vv_core::seqfeat::conservation;

use super::hex;
use crate::sequence::peers::Peers;
use crate::sequence::rows::{is_protein, letters_and_breaks};
use crate::sequence::tracks::{legend, Glyph, Inputs, TrackContext, TrackData, TrackProvider};

pub struct Burial;

/// Relative accessible area, as a fraction of the residue's maximum, below
/// which a residue falls in each class but the last.
const BURIAL_LIMITS: [f32; 3] = [0.10, 0.25, 0.50];

fn burial_kind(relative: f32) -> u8 {
    1 + BURIAL_LIMITS.iter().filter(|&&l| relative >= l).count() as u8
}

impl TrackProvider for Burial {
    fn id(&self) -> &'static str {
        "burial"
    }

    fn label(&self) -> &'static str {
        "Burial"
    }

    fn per_frame(&self) -> bool {
        true
    }

    fn inputs(&self) -> Inputs {
        Inputs {
            sasa: true,
            ..Inputs::default()
        }
    }

    /// Buried to exposed, by accessible area relative to the residue's
    /// maximum (Tien et al. 2013); a residue has no mark until the area
    /// is computed.
    fn compute(&self, ctx: &TrackContext) -> TrackData {
        let mut track = ctx.new_track(
            Glyph::Bar,
            vec![
                legend("Buried (under 10% exposed)", hex(0x1F3F8F)),
                legend("Mostly buried (10-25% exposed)", hex(0x4F7FCB)),
                legend("Partly exposed (25-50% exposed)", hex(0x9DBDE8)),
                legend("Exposed (50% or more)", hex(0xD9E5F6)),
            ],
        );
        for (r, &relative) in ctx.extras.rel_sasa.into_iter().flatten().enumerate() {
            if !relative.is_nan() {
                track.mark(r as u32, burial_kind(relative));
            }
        }
        track.finish()
    }
}

pub struct Conservation;

/// Most other chains a chain is compared with.
const MAX_COMPARED: usize = 150;

/// Score at or above which a residue is drawn as the most conserved class,
/// and the limits below which it falls in each lower one.
const CONSERVATION_LIMITS: [f32; 3] = [0.35, 0.60, 0.85];

fn conservation_kind(column: &conservation::Column) -> u8 {
    if column.identical == column.sequences {
        return 5;
    }
    1 + CONSERVATION_LIMITS
        .iter()
        .filter(|&&l| column.score >= l)
        .count() as u8
}

fn conservation_note(column: &conservation::Column) -> String {
    let letters: Vec<String> = column
        .counts
        .iter()
        .map(|&(letter, n)| format!("{} {n}", letter as char))
        .collect();
    format!(
        "{} of {} chains have {} ({:.0}%); score {:.2}\n{}",
        column.identical,
        column.sequences,
        column.residue as char,
        column.identity() * 100.0,
        column.score,
        letters.join(", ")
    )
}

impl TrackProvider for Conservation {
    fn id(&self) -> &'static str {
        "conservation"
    }

    fn label(&self) -> &'static str {
        "Conservation"
    }

    fn inputs(&self) -> Inputs {
        Inputs {
            peers: true,
            ..Inputs::default()
        }
    }

    /// Against every homologous protein chain of the loaded structures:
    /// the same entity, or at least 30% identical after alignment. A
    /// chain with no homolog has no track.
    fn compute(&self, ctx: &TrackContext) -> TrackData {
        let mut track = ctx.new_track(
            Glyph::Bar,
            vec![
                legend("Variable", hex(0xD8C9EA)),
                legend("Weakly conserved", hex(0xB79BD8)),
                legend("Conserved", hex(0x8E63BD)),
                legend("Highly conserved", hex(0x6437A0)),
                legend("Identical in every chain", hex(0x3B1F80)),
            ],
        );
        let Some(peers) = ctx.extras.peers else {
            return track.finish();
        };
        for (name, residues) in ctx.rows {
            let top = ctx.top();
            if residues.is_empty() || !is_protein(top, residues.start) {
                continue;
            }
            let (letters, _) = letters_and_breaks(top, residues.clone());
            let Some(result) = compare(peers, ctx, name, &letters) else {
                continue;
            };
            for (i, column) in result.columns.iter().enumerate() {
                let r = residues.start + i as u32;
                track.mark(r, conservation_kind(column));
                track.note(r, conservation_note(column));
            }
        }
        track.finish()
    }
}

fn compare(
    peers: &Peers,
    ctx: &TrackContext,
    name: &str,
    letters: &[u8],
) -> Option<conservation::Conservation> {
    let mut others = match ctx.extras.structure {
        Some(id) => peers.others(id, name),
        None => peers.chains.iter().map(|c| c.letters.as_slice()).collect(),
    };
    others.truncate(MAX_COMPARED);
    conservation::conservation(letters, &others)
}
