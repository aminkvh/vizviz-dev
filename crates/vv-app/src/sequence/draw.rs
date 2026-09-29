//! Painting the strip: letters, residue chips and the track glyphs, for the
//! visible columns only.

use std::ops::Range;

use egui::{pos2, vec2, Align2, Color32, FontId, Painter, Pos2, Rect, Shape, Stroke};
use vv_core::Topology;

use super::color::{from_packed, text_on};
use super::layout::{Block, TrackLine};
use super::rows::{one_letter, residue_number};
use super::tracks::Glyph;

/// Everything constant across one frame's painting.
pub struct Paint<'a> {
    pub painter: &'a Painter,
    pub origin: Pos2,
    pub label_width: f32,
    pub advance: f32,
    pub row_h: f32,
    pub font: FontId,
    pub small: FontId,
    pub text: Color32,
    pub dim: Color32,
    pub highlight: Color32,
}

impl Paint<'_> {
    fn cell(&self, block: &Block, column: usize, offset: f32, height: f32) -> Rect {
        let x = self.origin.x + self.label_width + column as f32 * self.advance;
        let y = self.origin.y + block.top + offset;
        Rect::from_min_size(pos2(x, y), vec2(self.advance, height))
    }
}

/// One block's sequence row and track rows over `columns`.
pub fn block(
    p: &Paint,
    b: &Block,
    top: &Topology,
    chips: Option<&[Color32]>,
    selected: Option<&[bool]>,
    columns: Range<usize>,
) {
    sequence_label(p, b);
    for column in columns.start..columns.end.min(b.row.residues.len()) {
        let residue = b.row.residues.start as usize + column;
        let cell = p.cell(b, column, 0.0, p.row_h);
        let chip = chips
            .map(|c| c[residue])
            .filter(|c| *c != Color32::TRANSPARENT);
        letter(
            p,
            cell,
            top,
            residue,
            chip,
            selected.is_some_and(|s| s[residue]),
        );
    }
    for line in &b.tracks {
        track_label(p, b, line);
        for column in columns.start..columns.end.min(b.row.residues.len()) {
            let residue = b.row.residues.start + column as u32;
            glyph(p, b, line, top, residue, column);
        }
    }
}

fn sequence_label(p: &Paint, b: &Block) {
    p.painter.text(
        pos2(p.origin.x + 4.0, p.origin.y + b.top + p.row_h * 0.5),
        Align2::LEFT_CENTER,
        &b.row.label,
        p.font.clone(),
        p.text,
    );
}

fn track_label(p: &Paint, b: &Block, line: &TrackLine) {
    p.painter.text(
        pos2(
            p.origin.x + 14.0,
            p.origin.y + b.top + line.offset + line.height * 0.5,
        ),
        Align2::LEFT_CENTER,
        line.provider.label(),
        p.small.clone(),
        p.dim,
    );
}

fn letter(
    p: &Paint,
    cell: Rect,
    top: &Topology,
    residue: usize,
    chip: Option<Color32>,
    selected: bool,
) {
    if selected {
        p.painter.rect_filled(cell.shrink(1.0), 2.0, p.highlight);
    }
    let (glyph, mut color) = match one_letter(top.residue_name(residue)) {
        Some(glyph) => (glyph, p.text),
        None => ('x', p.dim),
    };
    if let Some(fill) = chip {
        let inset = if selected { 2.5 } else { 1.0 };
        p.painter
            .rect_filled(cell.shrink2(vec2(inset, inset + 0.5)), 2.0, fill);
        color = text_on(fill);
    }
    p.painter.text(
        cell.center(),
        Align2::CENTER_CENTER,
        glyph,
        p.font.clone(),
        color,
    );
}

fn glyph(p: &Paint, b: &Block, line: &TrackLine, top: &Topology, residue: u32, column: usize) {
    let kind = line.data.kind(residue);
    if kind == 0 {
        return;
    }
    let cell = p.cell(b, column, line.offset, line.height);
    let color = from_packed(line.data.color(kind));
    match line.data.glyph {
        Glyph::Bar => {
            p.painter
                .rect_filled(cell.shrink2(vec2(0.0, 1.5)), 1.0, color);
        }
        Glyph::Structure => {
            structure_glyph(p, line, residue, kind, cell, color, b.row.residues.end)
        }
        Glyph::Ticks => tick(p, cell, &residue_number(top, residue)),
        Glyph::Edge => edge(p, cell, kind, color),
    }
}

fn structure_glyph(
    p: &Paint,
    line: &TrackLine,
    residue: u32,
    kind: u8,
    cell: Rect,
    color: Color32,
    end: u32,
) {
    match kind {
        1 => {
            p.painter
                .rect_filled(cell.shrink2(vec2(0.0, 0.5)), 2.0, color);
        }
        2 if residue + 1 >= end || line.data.kind(residue + 1) != 2 => arrow_head(p, cell, color),
        2 => {
            p.painter
                .rect_filled(cell.shrink2(vec2(0.0, 1.5)), 0.0, color);
        }
        _ => {
            p.painter
                .rect_filled(cell.shrink2(vec2(0.0, cell.height() * 0.38)), 0.0, color);
        }
    }
}

/// A strand's final residue: a triangle pointing at the C-terminus.
fn arrow_head(p: &Paint, cell: Rect, color: Color32) {
    let (l, r, t, b, m) = (
        cell.left(),
        cell.right(),
        cell.top(),
        cell.bottom(),
        cell.center().y,
    );
    p.painter.add(Shape::convex_polygon(
        vec![pos2(l, t + 0.5), pos2(r, m), pos2(l, b - 0.5)],
        color,
        Stroke::NONE,
    ));
}

fn tick(p: &Paint, cell: Rect, number: &str) {
    let x = cell.center().x;
    p.painter.line_segment(
        [pos2(x, cell.top()), pos2(x, cell.top() + 3.0)],
        Stroke::new(1.0, p.dim),
    );
    p.painter.text(
        pos2(x, cell.top() + 3.0),
        Align2::CENTER_TOP,
        number,
        p.small.clone(),
        p.dim,
    );
}

fn edge(p: &Paint, cell: Rect, kind: u8, color: Color32) {
    let x = if kind == 1 {
        cell.left()
    } else {
        cell.right() - 2.5
    };
    let bar = Rect::from_min_size(pos2(x, cell.top()), vec2(2.5, cell.height()));
    p.painter.rect_filled(bar, 1.0, color);
}
