//! CPU text rasterization for baking labels and captions into exported
//! screenshots. The live viewport draws the same text through egui's own
//! vector renderer instead (`ui::labels_and_captions_ui`); this module
//! exists only because `state::export_screenshot` renders to a raw pixel
//! buffer with no egui pass over it at all, so nothing egui-only reaches
//! the exported file otherwise.
//!
//! Font: Exo 2 Light, SIL Open Font License 1.1, bundled at
//! `assets/fonts/` (license text alongside it) -- a clean, legible sans
//! that needs no external font file or system font lookup at runtime.

use ab_glyph::{Font, FontRef, GlyphId, PxScale, ScaleFont};

static FONT_BYTES: &[u8] = include_bytes!("../assets/fonts/Exo2-Light.otf");

fn font() -> FontRef<'static> {
    FontRef::try_from_slice(FONT_BYTES).expect("bundled font must parse")
}

/// Walks the glyph outlines of `text` (Exo 2 Light) at `size_px`,
/// left-aligned with its top-left corner at `(x, y)`, calling `plot` with
/// each covered pixel's unclipped position and coverage (0, 1]. Shared by
/// `draw_text` (RGBA8 buffers) and the splash's packed-u32 blending
/// (`splash.rs`).
pub(crate) fn walk_glyphs(
    x: f32,
    y: f32,
    text: &str,
    size_px: f32,
    mut plot: impl FnMut(i32, i32, f32),
) {
    let font = font();
    let scale = PxScale::from(size_px);
    let scaled = font.as_scaled(scale);
    let mut cursor_x = x;
    let baseline_y = y + scaled.ascent();
    let mut previous: Option<GlyphId> = None;
    for ch in text.chars() {
        let id = font.glyph_id(ch);
        if let Some(prev) = previous {
            cursor_x += scaled.kern(prev, id);
        }
        let glyph = id.with_scale_and_position(scale, ab_glyph::point(cursor_x, baseline_y));
        if let Some(outlined) = font.outline_glyph(glyph) {
            let bounds = outlined.px_bounds();
            outlined.draw(|gx, gy, coverage| {
                if coverage > 0.0 {
                    plot(
                        bounds.min.x as i32 + gx as i32,
                        bounds.min.y as i32 + gy as i32,
                        coverage,
                    );
                }
            });
        }
        cursor_x += scaled.h_advance(id);
        previous = Some(id);
    }
}

/// The rendered width in pixels of `text` at `size_px`, for centering a
/// line before drawing it.
pub(crate) fn text_width(text: &str, size_px: f32) -> f32 {
    let font = font();
    let scaled = font.as_scaled(PxScale::from(size_px));
    let mut width = 0.0;
    let mut previous: Option<GlyphId> = None;
    for ch in text.chars() {
        let id = font.glyph_id(ch);
        if let Some(prev) = previous {
            width += scaled.kern(prev, id);
        }
        width += scaled.h_advance(id);
        previous = Some(id);
    }
    width
}

/// The font's ascent + descent at `size_px`: the vertical room one line
/// of text needs, for stacking lines without measuring glyphs.
pub(crate) fn line_height(size_px: f32) -> f32 {
    let font = font();
    let scaled = font.as_scaled(PxScale::from(size_px));
    scaled.ascent() - scaled.descent()
}

/// Draws `text` into an RGBA8 `buffer` of `width x height`, left-aligned
/// with its top-left corner at `(x, y)` in pixels, alpha-blending each
/// glyph over the existing pixels so it composites onto the rendered
/// scene rather than punching a hole in it. Silently clips at the buffer
/// edges and off-buffer positions.
pub fn draw_text(
    buffer: &mut [u8],
    width: u32,
    height: u32,
    x: f32,
    y: f32,
    text: &str,
    size_px: f32,
    color: [u8; 3],
) {
    walk_glyphs(x, y, text, size_px, |px, py, coverage| {
        if px < 0 || py < 0 || px as u32 >= width || py as u32 >= height {
            return;
        }
        let idx = ((py as u32 * width + px as u32) * 4) as usize;
        for c in 0..3 {
            let bg = buffer[idx + c] as f32;
            let fg = color[c] as f32;
            buffer[idx + c] = (fg * coverage + bg * (1.0 - coverage)).round() as u8;
        }
    });
}

/// A filled circle, alpha-blended the same way `draw_text` is -- the
/// small anchor marker `ui::labels_and_captions_ui` draws at each labeled
/// atom in the live view, reproduced here so an exported figure matches.
pub fn draw_dot(
    buffer: &mut [u8],
    width: u32,
    height: u32,
    cx: f32,
    cy: f32,
    radius: f32,
    color: [u8; 3],
) {
    let r = radius.ceil() as i32;
    let (cxi, cyi) = (cx.round() as i32, cy.round() as i32);
    for dy in -r..=r {
        for dx in -r..=r {
            let dist = ((dx * dx + dy * dy) as f32).sqrt();
            // One pixel of soft edge instead of a hard-edged circle.
            let coverage = (radius + 0.5 - dist).clamp(0.0, 1.0);
            if coverage <= 0.0 {
                continue;
            }
            let (px, py) = (cxi + dx, cyi + dy);
            if px < 0 || py < 0 || px as u32 >= width || py as u32 >= height {
                continue;
            }
            let idx = ((py as u32 * width + px as u32) * 4) as usize;
            for c in 0..3 {
                let bg = buffer[idx + c] as f32;
                let fg = color[c] as f32;
                buffer[idx + c] = (fg * coverage + bg * (1.0 - coverage)).round() as u8;
            }
        }
    }
}

/// A dashed line from `(x0, y0)` to `(x1, y1)`, `dash` pixels on and
/// `gap` off, stamped as small dots -- a measurement's line, as
/// `ui::labels_and_captions_ui` draws it live.
#[allow(clippy::too_many_arguments)]
pub fn draw_dashed_line(
    buffer: &mut [u8],
    width: u32,
    height: u32,
    (x0, y0): (f32, f32),
    (x1, y1): (f32, f32),
    [dash, gap]: [f32; 2],
    color: [u8; 3],
) {
    let length = ((x1 - x0).powi(2) + (y1 - y0).powi(2)).sqrt();
    let mut s = 0.0;
    while s <= length {
        if s % (dash + gap) < dash {
            let t = s / length.max(1e-6);
            draw_dot(
                buffer,
                width,
                height,
                x0 + (x1 - x0) * t,
                y0 + (y1 - y0) * t,
                0.75,
                color,
            );
        }
        s += 0.5;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn draw_text_darkens_some_pixels_without_touching_alpha() {
        let mut buffer = vec![200u8; 64 * 16 * 4];
        for px in buffer.chunks_exact_mut(4) {
            px[3] = 255;
        }
        draw_text(&mut buffer, 64, 16, 2.0, 2.0, "Ag", 12.0, [0, 0, 0]);
        let changed = buffer.chunks_exact(4).filter(|px| px[0] < 200).count();
        assert!(changed > 0, "expected at least one darkened pixel");
        assert!(
            buffer.chunks_exact(4).all(|px| px[3] == 255),
            "alpha must stay opaque"
        );
    }

    #[test]
    fn draw_text_clips_at_buffer_edges_without_panicking() {
        let mut buffer = vec![0u8; 8 * 8 * 4];
        // Text starting near the right/bottom edge and running off it.
        draw_text(
            &mut buffer,
            8,
            8,
            6.0,
            6.0,
            "Wide text",
            20.0,
            [255, 255, 255],
        );
        // Negative position too.
        draw_text(
            &mut buffer,
            8,
            8,
            -50.0,
            -50.0,
            "Off screen",
            20.0,
            [255, 255, 255],
        );
    }

    #[test]
    fn draw_dot_covers_its_center_pixel() {
        let mut buffer = vec![0u8; 10 * 10 * 4];
        for px in buffer.chunks_exact_mut(4) {
            px[3] = 255;
        }
        draw_dot(&mut buffer, 10, 10, 5.0, 5.0, 3.0, [255, 200, 50]);
        let idx = ((5 * 10 + 5) * 4) as usize;
        assert!(buffer[idx] > 200, "center pixel should be near full color");
    }
}
