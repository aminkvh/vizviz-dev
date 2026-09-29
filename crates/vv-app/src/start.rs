//! The start card's logo texture (`ui.rs`'s `start_card_ui`): the mark
//! and wordmark PNG (`splash::full_logo`) decoded once per theme and
//! uploaded to the GPU, not re-decoded every frame the card is on
//! screen.

use std::sync::OnceLock;

use crate::theme::ThemeMode;

/// The mark and wordmark as a texture, matching `mode`. Decoded and
/// uploaded once per theme and reused for the life of the process.
pub fn logo_texture(ctx: &egui::Context, mode: ThemeMode) -> Option<egui::TextureHandle> {
    static DARK: OnceLock<Option<egui::TextureHandle>> = OnceLock::new();
    static LIGHT: OnceLock<Option<egui::TextureHandle>> = OnceLock::new();
    let (cell, name) = match mode {
        ThemeMode::Dark => (&DARK, "logo-full-dark"),
        ThemeMode::Light => (&LIGHT, "logo-full-light"),
    };
    cell.get_or_init(|| {
        let (rgba, w, h) = crate::splash::full_logo(mode)?;
        let image = egui::ColorImage::from_rgba_unmultiplied([w as usize, h as usize], &rgba);
        Some(ctx.load_texture(name, image, egui::TextureOptions::LINEAR))
    })
    .clone()
}
