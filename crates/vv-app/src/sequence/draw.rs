//! Painting the strip: letters, residue chips and the track glyphs, for the
//! visible columns only.

use std::ops::Range;

use egui::{pos2, vec2, Align2, Color32, FontId, Painter, Pos2, Rect, Shape, Stroke};
use vv_core::Topology;

use super::color::{from_packed, readable_on, tint};
use super::layout::{Block, TrackLine, TICK_HEIGHT};
use super::rows::{one_letter, residue_number};
use super::tracks::Glyph;

/// Everything constant across one frame's painting.
pub struct Paint<'a> {
    pub painter: &'a Painter,
    pub origin: Pos2,
    pub label_width: f32,
    /// Left edge of the label column: pinned to the visible area.
    pub label_x: f32,
    pub panel: Color32,
    /// Whether residue chips are drawn as soft tints.
    pub soft: bool,
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
        for column in columns.start..columns.end.min(b.row.residues.len()) {
            let residue = b.row.residues.start + column as u32;
            glyph(p, b, line, top, residue, column);
        }
        badges(p, b, line);
    }
    labels(p, b);
}

/// The label column over the letters, so it stays put while they scroll.
fn labels(p: &Paint, b: &Block) {
    let top = p.origin.y + b.top;
    let pane = Rect::from_min_size(pos2(p.label_x, top), vec2(p.label_width, b.height));
    p.painter.rect_filled(pane, 0.0, p.panel);
    p.painter.text(
        pos2(p.label_x + 4.0, top + p.row_h * 0.5),
        Align2::LEFT_CENTER,
        &b.row.label,
        p.font.clone(),
        p.text,
    );
    for line in &b.tracks {
        p.painter.text(
            pos2(
                p.label_x + TRACK_INDENT,
                top + line.offset + line.height * 0.5,
            ),
            Align2::LEFT_CENTER,
            line.provider.label(),
            p.small.clone(),
            p.dim,
        );
    }
}

/// Where a track's label starts, right of the chain label's edge.
pub const TRACK_INDENT: f32 = 14.0;

/// Tags such as a domain's chain type, on the bar row of their residue.
fn badges(p: &Paint, b: &Block, line: &TrackLine) {
    for (residue, text) in line.data.badges() {
        if !b.row.residues.contains(&residue) {
            continue;
        }
        let column = (residue - b.row.residues.start) as usize;
        let bar = p
            .cell(
                b,
                column,
                line.offset + TICK_HEIGHT,
                line.height - TICK_HEIGHT,
            )
            .shrink2(vec2(0.0, 1.5));
        let galley = p
            .painter
            .layout_no_wrap(text.to_string(), p.small.clone(), p.dim);
        let pill = Rect::from_min_size(
            bar.min + vec2(2.0, 0.0),
            vec2(galley.size().x + 8.0, bar.height()),
        );
        let fill = Color32::from_rgb(0x4A, 0x55, 0x66);
        p.painter.rect_filled(pill, 3.0, fill);
        let color = readable_on(fill, Color32::WHITE);
        p.painter.galley_with_override_text_color(
            pill.min + vec2(4.0, (pill.height() - galley.size().y) * 0.5),
            galley,
            color,
        );
    }
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
    if let Some(fill) = chip.map(|c| if p.soft { tint(c, p.panel) } else { c }) {
        let inset = if selected { 2.5 } else { 1.0 };
        p.painter
            .rect_filled(cell.shrink2(vec2(inset, inset + 0.5)), 2.0, fill);
        color = readable_on(fill, p.text);
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
        Glyph::LabeledBar => labeled_bar(p, line, cell, residue, color),
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

fn labeled_bar(p: &Paint, line: &TrackLine, cell: Rect, residue: u32, color: Color32) {
    let (ticks, bar) = cell.split_top_bottom_at_y(cell.top() + TICK_HEIGHT);
    if let Some(text) = line.data.tick(residue) {
        tick(p, ticks, text);
    }
    p.painter
        .rect_filled(bar.shrink2(vec2(0.0, 1.5)), 1.0, color);
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
