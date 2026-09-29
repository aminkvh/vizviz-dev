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

mod cache;
mod color;
mod command;
mod draw;
mod header;
mod layout;
mod providers;
mod rows;
mod tooltip;
mod tracks;

use egui::{Sense, Ui};
use vv_scene::{Scene, StructureId};

use crate::ui::AppUi;
use color::SeqColor;

/// Room for "<structure> <chain>" at the left of every row.
const LABEL_WIDTH: f32 = 96.0;

/// What the strip shows, and what it has computed for it.
pub struct SequenceState {
    pub color: SeqColor,
    /// Ids of the enabled tracks (`tracks::PROVIDERS`).
    pub tracks: Vec<&'static str>,
    pub legend: bool,
    cache: cache::Cache,
}

impl Default for SequenceState {
    fn default() -> Self {
        Self {
            color: SeqColor::None,
            tracks: Vec::new(),
            legend: false,
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

/// A click on the strip: the residue, and whether it adds to the selection.
struct Pick {
    structure: StructureId,
    residue: u32,
    add: bool,
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
        let blocks = layout::blocks(scene, &mut state.cache, &state.tracks, row_h);
        let longest = blocks
            .iter()
            .map(|b| b.row.residues.len())
            .max()
            .unwrap_or(0);
        let height = blocks.last().map_or(0.0, |b| b.top + b.height);
        let total = egui::vec2(LABEL_WIDTH + longest as f32 * advance, height);
        let selected = rows::selected_residues(scene);
        let mut chips = std::collections::HashMap::new();
        for b in &blocks {
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
                    label_width: LABEL_WIDTH,
                    advance,
                    row_h,
                    font,
                    small,
                    text: ui.visuals().text_color(),
                    dim: ui.visuals().weak_text_color(),
                    highlight: ui.visuals().selection.bg_fill,
                };
                let first_col =
                    ((viewport.min.x - LABEL_WIDTH) / advance).floor().max(0.0) as usize;
                let last_col = ((viewport.max.x - LABEL_WIDTH) / advance).ceil().max(0.0) as usize;
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

                let hit = |pos: egui::Pos2| -> Option<(usize, usize)> {
                    let local = pos - rect.min;
                    if local.x < LABEL_WIDTH {
                        return None;
                    }
                    let i = layout::block_at(&blocks, local.y)?;
                    let column = ((local.x - LABEL_WIDTH) / advance).floor() as usize;
                    (column < blocks[i].row.residues.len()).then_some((i, column))
                };
                if let Some((i, column)) = response.hover_pos().and_then(hit) {
                    let b = &blocks[i];
                    let residue = b.row.residues.start + column as u32;
                    tooltip::show(&response, scene, &state.cache, b, residue);
                }
                if let Some((i, column)) = response.interact_pointer_pos().and_then(hit) {
                    let structure = blocks[i].row.structure;
                    let residue = blocks[i].row.residues.start + column as u32;
                    if response.double_clicked() {
                        focused = Some((structure, residue));
                    } else if response.clicked() {
                        let add = ui.input(|i| i.modifiers.command);
                        picked = Some(Pick {
                            structure,
                            residue,
                            add,
                        });
                    }
                }
            });
        (picked, focused)
    }

    fn select_residue(&mut self, pick: Pick) {
        let atoms: Option<Vec<u32>> = self.scene.structure(pick.structure).map(|s| {
            s.structure.topology.residues[pick.residue as usize]
                .atoms
                .clone()
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
