//! The sequence strip: one row of one-letter residue codes per chain, for
//! every loaded structure, docked above the viewport by default. Click a
//! residue to select its atoms, Ctrl+click to add it to the selection,
//! double-click to zoom to it with its 5 Å neighbourhood (as the viewport
//! does), hover for its name, number and annotations. Drawn straight with
//! the painter and virtualized in both directions, so a 2.4M-atom
//! structure costs only the glyphs actually on screen.
//!
//! Residue letters can be colored by a scheme (`color.rs`) and each chain
//! can carry annotation tracks (`tracks.rs`, `providers/`); the header
//! (`header.rs`) chooses both and `sequence ...` commands (`command.rs`)
//! do the same.

mod background;
mod cache;
mod chain_props;
#[cfg(test)]
mod checks;
mod color;
mod command;
mod draw;
mod header;
mod layout;
mod peers;
mod prefs;
mod providers;
mod rows;
mod tooltip;
mod tracks;
mod uniprot;

use egui::{Sense, Ui};
use vv_scene::{Scene, StructureId};

use crate::ui::{AppUi, LayoutRequest};
use color::SeqColor;
use tracks::AntibodySettings;

/// Narrowest label column: room for "<structure> <chain>".
const MIN_LABEL_WIDTH: f32 = 96.0;
/// Panel height beyond the strip that a fit leaves free, for the
/// horizontal scroll bar.
const FIT_SLACK: f32 = 20.0;

/// What the strip shows, and what it has computed for it.
pub struct SequenceState {
    pub color: SeqColor,
    /// Ids of the enabled tracks (`tracks::PROVIDERS`).
    pub tracks: Vec<&'static str>,
    pub legend: bool,
    pub antibody: AntibodySettings,
    /// Strip height and free space the panel was last fitted at.
    fitted: Option<(i32, i32)>,
    cache: cache::Cache,
}

impl Default for SequenceState {
    fn default() -> Self {
        Self {
            color: SeqColor::None,
            tracks: Vec::new(),
            legend: false,
            antibody: AntibodySettings::default(),
            fitted: None,
            cache: cache::Cache::default(),
        }
    }
}

impl SequenceState {
    pub fn track_on(&self, id: &str) -> bool {
        self.tracks.contains(&id)
    }

    pub fn set_track(&mut self, id: &'static str, on: bool) {
        self.tracks.retain(|t| *t != id);
        if on {
            self.tracks.push(id);
        }
    }
}

/// Wide enough for every chain label and every shown track label.
fn label_width(
    ui: &Ui,
    blocks: &[layout::Block],
    font: &egui::FontId,
    small: &egui::FontId,
) -> f32 {
    let width = |text: &str, font: &egui::FontId| {
        ui.fonts_mut(|f| {
            f.layout_no_wrap(text.to_string(), font.clone(), egui::Color32::WHITE)
                .size()
                .x
        })
    };
    let mut wide = MIN_LABEL_WIDTH;
    for b in blocks {
        wide = wide.max(width(&b.row.label, font) + 12.0);
        for line in &b.tracks {
            let label = width(line.provider.label(), small);
            wide = wide.max(label + draw::TRACK_INDENT + 12.0);
        }
    }
    wide
}

/// A click on the strip: the residues it selects (one, or every residue
/// of a ligand chip), and whether they add to the selection.
struct Pick {
    structure: StructureId,
    residues: Vec<u32>,
    add: bool,
}

/// The structure and residues a click on `target` acts on: one residue,
/// or all of a ligand chip's.
fn clicked_residues(
    blocks: &[layout::Block],
    target: layout::Hit,
) -> Option<(StructureId, Vec<u32>)> {
    match target {
        layout::Hit::Residue { block, residue } => {
            Some((blocks[block].row.structure, vec![residue]))
        }
        layout::Hit::Ligands { block, group } => {
            let row = &blocks[block].row;
            Some((row.structure, row.ligands[group].residues.clone()))
        }
        layout::Hit::Label { .. } => None,
    }
}

impl AppUi<'_> {
    pub(crate) fn sequence_ui(&mut self, ui: &mut Ui) {
        if self.scene.structures().next().is_none() {
            ui.weak("No structure loaded.");
            return;
        }
        header::show(ui, self.sequence, self.scene);
        let (picked, focused) = self.strip(ui);
        if let Some(pick) = picked {
            self.select_residue(pick);
        }
        if let Some((id, residue)) = focused {
            self.focus_residue(id, residue);
        }
    }

    /// Draws the scrolling strip; returns what the pointer picked.
    fn strip(&mut self, ui: &mut Ui) -> (Option<Pick>, Option<(StructureId, u32)>) {
        let state = &mut *self.sequence;
        let scene: &Scene = self.scene;
        state.cache.retain_open(scene);
        let font = egui::TextStyle::Monospace.resolve(ui.style());
        let mut small = egui::TextStyle::Small.resolve(ui.style());
        small.size *= 0.85;
        let advance = ui.ctx().fonts_mut(|f| f.glyph_width(&font, 'W')) + 2.0;
        let row_h = ui.text_style_height(&egui::TextStyle::Monospace) + 6.0;
        state.cache.set_wake(ui.ctx());
        let env = state
            .cache
            .env(scene, state.antibody, state.track_on("conservation"));
        let blocks = layout::blocks(scene, &mut state.cache, &state.tracks, &env, row_h);
        let longest = blocks.iter().map(|b| b.row.columns()).max().unwrap_or(0);
        let height = blocks.last().map_or(0.0, |b| b.top + b.height);
        let label_width = label_width(ui, &blocks, &font, &small);
        let total = egui::vec2(label_width + longest as f32 * advance, height);
        let shortfall = height + FIT_SLACK - ui.available_height();
        let seen = (height as i32, ui.available_height() as i32);
        if shortfall > 0.0 && state.fitted != Some(seen) {
            state.fitted = Some(seen);
            *self.layout_request = Some(LayoutRequest::GrowSequence(shortfall));
        }
        let selected = rows::selected_residues(scene);
        let mut chips = std::collections::HashMap::new();
        for b in blocks.iter().filter(|b| b.row.ligands.is_empty()) {
            let id = b.row.structure;
            if let (None, Some(loaded)) = (chips.get(&id), scene.structure(id)) {
                chips.insert(id, state.cache.colors(id, loaded, state.color));
            }
        }

        let mut picked = None;
        let mut focused = None;
        egui::ScrollArea::both()
            .auto_shrink([false, false])
            .show_viewport(ui, |ui, viewport| {
                let (rect, response) = ui.allocate_exact_size(total, Sense::click());
                let painter = ui.painter_at(rect);
                let paint = draw::Paint {
                    painter: &painter,
                    origin: rect.min,
                    label_width,
                    label_x: rect.min.x + viewport.min.x,
                    panel: ui.visuals().panel_fill,
                    soft: state.color.categorical(),
                    advance,
                    row_h,
                    font,
                    small,
                    text: ui.visuals().text_color(),
                    dim: ui.visuals().weak_text_color(),
                    highlight: ui.visuals().selection.bg_fill,
                };
                let first_col =
                    ((viewport.min.x - label_width) / advance).floor().max(0.0) as usize;
                let last_col = ((viewport.max.x - label_width) / advance).ceil().max(0.0) as usize;
                let first = layout::block_at(&blocks, viewport.min.y).unwrap_or(0);
                for b in blocks
                    .iter()
                    .skip(first)
                    .take_while(|b| b.top < viewport.max.y)
                {
                    let Some(loaded) = scene.structure(b.row.structure) else {
                        continue;
                    };
                    let flags = selected
                        .as_ref()
                        .filter(|(id, _)| *id == b.row.structure)
                        .map(|(_, f)| f.as_slice());
                    let colors = chips.get(&b.row.structure).and_then(|c| c.as_deref());
                    let top = &loaded.structure.topology;
                    draw::block(
                        &paint,
                        b,
                        top,
                        colors.map(|c| c.as_slice()),
                        flags,
                        first_col..last_col,
                    );
                }

                let cols = layout::Columns {
                    label_width,
                    label_left: viewport.min.x,
                    advance,
                    row_h,
                };
                let hit = |pos: egui::Pos2| layout::hit(&blocks, pos - rect.min, &cols);
                match response.hover_pos().and_then(hit) {
                    Some(layout::Hit::Residue { block, residue }) => tooltip::residue(
                        &response,
                        scene,
                        &mut state.cache,
                        &blocks[block],
                        residue,
                    ),
                    Some(layout::Hit::Label { block }) => {
                        tooltip::chain(&response, scene, &blocks[block])
                    }
                    Some(layout::Hit::Ligands { block, group }) => {
                        tooltip::ligands(&response, scene, &blocks[block], group)
                    }
                    None => {}
                }
                if let Some(target) = response.interact_pointer_pos().and_then(hit) {
                    if let Some((structure, residues)) = clicked_residues(&blocks, target) {
                        if response.double_clicked() {
                            focused = Some((structure, residues[0]));
                        } else if response.clicked() {
                            let add = ui.input(|i| i.modifiers.command);
                            picked = Some(Pick {
                                structure,
                                residues,
                                add,
                            });
                        }
                    }
                }
            });
        (picked, focused)
    }

    fn select_residue(&mut self, pick: Pick) {
        let atoms: Option<Vec<u32>> = self.scene.structure(pick.structure).map(|s| {
            let top = &s.structure.topology;
            pick.residues
                .iter()
                .flat_map(|&r| top.residues[r as usize].atoms.clone())
                .collect()
        });
        if let Some(atoms) = atoms {
            self.select_atoms(pick.structure, &atoms, pick.add);
        }
    }

    fn focus_residue(&mut self, id: StructureId, residue: u32) {
        let first = self
            .scene
            .structure(id)
            .map(|s| s.structure.topology.residues[residue as usize].atoms.start);
        if let Some(atom) = first {
            self.zoom_to_residue(id, atom);
        }
    }
}
