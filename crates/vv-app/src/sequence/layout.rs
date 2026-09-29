//! Where each chain's block sits in the strip: its sequence row plus the
//! thin track rows that have something to show for that chain.

use std::sync::Arc;

use vv_scene::Scene;

use super::cache::Cache;
use super::rows::{rows_of, Row};
use super::tracks::{AntibodySettings, Glyph, TrackData, TrackProvider, PROVIDERS};

pub const BAR_HEIGHT: f32 = 13.0;
pub const TICK_HEIGHT: f32 = 15.0;
const BLOCK_GAP: f32 = 4.0;

pub struct TrackLine {
    pub provider: &'static dyn TrackProvider,
    pub data: Arc<TrackData>,
    /// Distance below the block's top.
    pub offset: f32,
    pub height: f32,
}

pub struct Block {
    pub row: Row,
    pub top: f32,
    pub height: f32,
    pub tracks: Vec<TrackLine>,
}

fn line_height(glyph: Glyph) -> f32 {
    match glyph {
        Glyph::Ticks => TICK_HEIGHT,
        Glyph::LabeledBar => TICK_HEIGHT + BAR_HEIGHT,
        _ => BAR_HEIGHT,
    }
}

fn track_lines(
    row: &Row,
    scene: &Scene,
    cache: &mut Cache,
    enabled: &[&'static str],
    antibody: AntibodySettings,
    row_h: f32,
) -> Vec<TrackLine> {
    let Some(loaded) = scene.structure(row.structure) else {
        return Vec::new();
    };
    let mut offset = row_h;
    let mut lines = Vec::new();
    for provider in PROVIDERS.iter().filter(|p| enabled.contains(&p.id())) {
        let data = cache.track(row.structure, loaded, *provider, antibody);
        if data.any_in(&row.residues) {
            let height = line_height(data.glyph);
            lines.push(TrackLine {
                provider: *provider,
                data,
                offset,
                height,
            });
            offset += height;
        }
    }
    lines
}

/// One block per chain row, stacked from `y = 0`.
pub fn blocks(
    scene: &Scene,
    cache: &mut Cache,
    enabled: &[&'static str],
    antibody: AntibodySettings,
    row_h: f32,
) -> Vec<Block> {
    let mut top = 0.0;
    rows_of(scene)
        .into_iter()
        .map(|row| {
            let tracks = track_lines(&row, scene, cache, enabled, antibody, row_h);
            let height = tracks.last().map_or(row_h, |t| t.offset + t.height) + BLOCK_GAP;
            let block = Block {
                row,
                top,
                height,
                tracks,
            };
            top += height;
            block
        })
        .collect()
}

/// Index of the block containing vertical position `y`.
pub fn block_at(blocks: &[Block], y: f32) -> Option<usize> {
    let after = blocks.partition_point(|b| b.top <= y);
    let i = after.checked_sub(1)?;
    (y < blocks[i].top + blocks[i].height).then_some(i)
}
