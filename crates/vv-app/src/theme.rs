//! Design tokens (DESIGN.md): the one palette, spacing scale, type scale
//! and radii every panel uses; nothing else sets a colour, size or
//! spacing. Colours come from the vizviz logo: navy bonds, deep teal
//! centre, azure and violet atoms, periwinkle wings. Dark chrome is
//! navy-slate so the molecule is the brightest thing on screen; the light
//! theme is for figures and daylight. Text meets WCAG AA (4.5:1) on every
//! surface it is drawn on in both themes (`tests` checks the pairs).

use egui::{Color32, CornerRadius, FontFamily, FontId, Stroke, TextStyle};
use serde::{Deserialize, Serialize};

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub enum ThemeMode {
    #[default]
    Dark,
    Light,
}

impl ThemeMode {
    pub fn name(self) -> &'static str {
        match self {
            ThemeMode::Dark => "dark",
            ThemeMode::Light => "light",
        }
    }

    pub fn parse(word: &str) -> Option<ThemeMode> {
        match word {
            "dark" => Some(ThemeMode::Dark),
            "light" => Some(ThemeMode::Light),
            _ => None,
        }
    }
}

/// UI preferences kept with the saved layout, not in sessions: they are
/// about the person, not the figure.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct UiPrefs {
    #[serde(default)]
    pub theme: ThemeMode,
    /// The ribbon shows only its tab names.
    #[serde(default)]
    pub ribbon_collapsed: bool,
    /// Set by the quit dialog's "Don't remind me again" checkbox: quitting
    /// with unsaved changes no longer asks first. `confirmquit on` clears
    /// it.
    #[serde(default)]
    pub dont_confirm_quit: bool,
    /// Set by the start card's "Don't show again" checkbox. `welcome on`
    /// clears it.
    #[serde(default)]
    pub hide_start_card: bool,
    /// Structures and sessions opened from the start card's Recent list,
    /// most recent first; capped at `MAX_RECENTS`.
    #[serde(default)]
    pub recent_files: Vec<std::path::PathBuf>,
}

const MAX_RECENTS: usize = 8;

impl UiPrefs {
    /// Adds `path` to the front of `recent_files`, deduplicating and
    /// capping the list.
    pub fn add_recent(&mut self, path: std::path::PathBuf) {
        self.recent_files.retain(|p| p != &path);
        self.recent_files.insert(0, path);
        self.recent_files.truncate(MAX_RECENTS);
    }
}

/// The semantic colours of one theme (DESIGN.md, "Colour").
#[derive(Clone, Copy, Debug)]
pub struct Tokens {
    /// Behind docked panels.
    pub bg: Color32,
    /// Panels and dialogs.
    pub surface: Color32,
    /// Buttons, inputs and cards.
    pub surface_raised: Color32,
    /// The hover fill of anything clickable.
    pub hover: Color32,
    pub border: Color32,
    pub text: Color32,
    /// Secondary text: captions, hints, placeholders.
    pub text_muted: Color32,
    /// The primary action, selected state, links and focus ring: azure.
    pub primary: Color32,
    pub on_primary: Color32,
    pub danger: Color32,
    /// Done or confirmed: the logo's teal.
    pub success: Color32,
    /// Slow or partial results.
    pub warning: Color32,
}

impl Tokens {
    pub fn of(mode: ThemeMode) -> Tokens {
        match mode {
            ThemeMode::Dark => Tokens {
                bg: Color32::from_rgb(0x0B, 0x14, 0x24),
                surface: Color32::from_rgb(0x12, 0x1F, 0x36),
                surface_raised: Color32::from_rgb(0x1A, 0x2C, 0x4B),
                hover: Color32::from_rgb(0x24, 0x3C, 0x63),
                border: Color32::from_rgb(0x2C, 0x43, 0x69),
                text: Color32::from_rgb(0xE6, 0xED, 0xF7),
                text_muted: Color32::from_rgb(0x9F, 0xB0, 0xC8),
                primary: Color32::from_rgb(0x4D, 0xA3, 0xFF),
                on_primary: Color32::from_rgb(0x06, 0x11, 0x1F),
                danger: Color32::from_rgb(0xFF, 0x7A, 0x7A),
                success: Color32::from_rgb(0x2B, 0xB5, 0xA6),
                warning: Color32::from_rgb(0xF2, 0xB8, 0x4B),
            },
            ThemeMode::Light => Tokens {
                bg: Color32::from_rgb(0xDD, 0xE5, 0xF1),
                surface: Color32::from_rgb(0xF6, 0xF8, 0xFC),
                surface_raised: Color32::from_rgb(0xFF, 0xFF, 0xFF),
                hover: Color32::from_rgb(0xE3, 0xEC, 0xF9),
                border: Color32::from_rgb(0xBF, 0xCC, 0xE0),
                text: Color32::from_rgb(0x13, 0x29, 0x4B),
                text_muted: Color32::from_rgb(0x45, 0x56, 0x73),
                primary: Color32::from_rgb(0x1A, 0x66, 0xD4),
                on_primary: Color32::from_rgb(0xFF, 0xFF, 0xFF),
                danger: Color32::from_rgb(0xB4, 0x23, 0x18),
                success: Color32::from_rgb(0x15, 0x7F, 0x75),
                warning: Color32::from_rgb(0x8A, 0x55, 0x00),
            },
        }
    }

    /// The tokens egui is currently drawing with.
    pub fn current(ctx: &egui::Context) -> Tokens {
        Tokens::of(match ctx.theme() {
            egui::Theme::Dark => ThemeMode::Dark,
            egui::Theme::Light => ThemeMode::Light,
        })
    }
}

/// The spacing scale (logical px): 4, 8, 12, 16, 24, 32 and nothing else.
pub mod space {
    /// Icon to its label; label to its field.
    pub const TIGHT: f32 = 4.0;
    /// Between controls in a row or column.
    pub const GAP: f32 = 8.0;
    /// Inside cards and dialogs; a control to its help text.
    pub const INSET: f32 = 12.0;
    /// Panel padding; between a panel's sections.
    pub const PAD: f32 = 16.0;
    /// Dialog padding; around empty states.
    pub const WIDE: f32 = 24.0;
    pub const WIDEST: f32 = 32.0;
}

/// Corner radii: controls, and containers.
pub mod radius {
    pub const CONTROL: u8 = 4;
    pub const CONTAINER: u8 = 8;
}

/// Height of every button, input and select; also their hit target.
pub const CONTROL_HEIGHT: f32 = 28.0;
/// A ribbon group's main action: icon over label.
pub const LARGE_CONTROL_HEIGHT: f32 = 56.0;

/// The six type styles (DESIGN.md, "Type"). Label is Body, strong.
pub mod text {
    pub const CAPTION: f32 = 11.0;
    pub const BODY: f32 = 13.0;
    pub const SECTION: f32 = 15.0;
    pub const TITLE: f32 = 18.0;
    pub const MONO: f32 = 12.0;
    /// Icon glyphs: in icon-only buttons, and over a large button's label.
    pub const ICON: f32 = 18.0;
    pub const LARGE_ICON: f32 = 24.0;
}

/// egui's name for the Title style (Caption, Body, Section and Mono map
/// onto its Small, Body, Heading and Monospace).
pub fn title_style() -> TextStyle {
    TextStyle::Name("title".into())
}

/// Loads the icon font (Phosphor, MIT) next to egui's own fonts, so icon
/// glyphs mix into any label.
pub fn install_fonts(ctx: &egui::Context) {
    let mut fonts = egui::FontDefinitions::default();
    egui_phosphor::add_to_fonts(&mut fonts, egui_phosphor::Variant::Regular);
    ctx.set_fonts(fonts);
}

/// Applies `mode`'s tokens and the type and spacing scales to egui.
pub fn apply(ctx: &egui::Context, mode: ThemeMode) {
    let t = Tokens::of(mode);
    let mut visuals = match mode {
        ThemeMode::Dark => egui::Visuals::dark(),
        ThemeMode::Light => egui::Visuals::light(),
    };
    let control = CornerRadius::same(radius::CONTROL);
    visuals.panel_fill = t.surface;
    visuals.window_fill = t.surface;
    visuals.window_stroke = Stroke::new(1.0, t.border);
    visuals.window_corner_radius = CornerRadius::same(radius::CONTAINER);
    visuals.menu_corner_radius = CornerRadius::same(radius::CONTAINER);
    visuals.extreme_bg_color = t.bg;
    visuals.faint_bg_color = t.surface_raised;
    visuals.code_bg_color = t.surface_raised;
    visuals.hyperlink_color = t.primary;
    visuals.error_fg_color = t.danger;
    visuals.warn_fg_color = t.warning;
    visuals.weak_text_color = Some(t.text_muted);
    visuals.selection.bg_fill = t.primary;
    visuals.selection.stroke = Stroke::new(2.0, t.on_primary);
    let w = &mut visuals.widgets;
    w.noninteractive.bg_fill = t.surface;
    w.noninteractive.weak_bg_fill = t.surface;
    w.noninteractive.bg_stroke = Stroke::new(1.0, t.border);
    w.noninteractive.fg_stroke = Stroke::new(1.0, t.text);
    w.noninteractive.corner_radius = control;
    for (state, fill) in [
        (&mut w.inactive, t.surface_raised),
        (&mut w.hovered, t.hover),
        (&mut w.active, t.hover),
        (&mut w.open, t.hover),
    ] {
        state.bg_fill = fill;
        state.weak_bg_fill = fill;
        state.bg_stroke = Stroke::new(1.0, t.border);
        state.fg_stroke = Stroke::new(1.0, t.text);
        state.corner_radius = control;
        state.expansion = 0.0;
    }
    // Pressed: the hover fill with a primary border.
    w.active.bg_stroke = Stroke::new(1.0, t.primary);
    let theme = match mode {
        ThemeMode::Dark => egui::Theme::Dark,
        ThemeMode::Light => egui::Theme::Light,
    };
    ctx.set_theme(theme);
    ctx.set_visuals_of(theme, visuals);
    ctx.all_styles_mut(|style| {
        style.text_styles = [
            (TextStyle::Small, FontId::proportional(text::CAPTION)),
            (TextStyle::Body, FontId::proportional(text::BODY)),
            (TextStyle::Button, FontId::proportional(text::BODY)),
            (TextStyle::Heading, FontId::proportional(text::SECTION)),
            (title_style(), FontId::proportional(text::TITLE)),
            (
                TextStyle::Monospace,
                FontId::new(text::MONO, FontFamily::Monospace),
            ),
        ]
        .into();
        let s = &mut style.spacing;
        s.item_spacing = egui::vec2(space::GAP, space::TIGHT);
        s.button_padding = egui::vec2(space::GAP, space::TIGHT);
        s.interact_size = egui::vec2(CONTROL_HEIGHT, CONTROL_HEIGHT);
        s.icon_spacing = space::TIGHT;
        s.window_margin = egui::Margin::same(space::WIDE as i8);
        s.menu_margin = egui::Margin::same(space::TIGHT as i8);
        s.indent = space::PAD;
    });
}

/// The dock's panels, tab bars and splitters in `mode`'s tokens. Tab
/// bodies have no margin: the viewport fills its panel to the edge, and
/// every other panel pads itself (`space::PAD`).
pub fn dock_style(style: &egui::Style, mode: ThemeMode) -> egui_dock::Style {
    let t = Tokens::of(mode);
    let mut dock = egui_dock::Style::from_egui(style);
    dock.dock_area_padding = None;
    dock.main_surface_border_stroke = Stroke::NONE;
    dock.separator.width = 1.0;
    dock.separator.color_idle = t.border;
    dock.separator.color_hovered = t.primary;
    dock.separator.color_dragged = t.primary;
    dock.tab_bar.bg_fill = t.bg;
    dock.tab_bar.hline_color = t.border;
    dock.tab_bar.height = CONTROL_HEIGHT;
    dock.tab_bar.corner_radius = CornerRadius::ZERO;
    // Overflowing tabs still scroll (wheel or drag); no bar across them.
    dock.tab_bar.show_scroll_bar_on_overflow = false;
    let tab = |fill: Color32, text: Color32| egui_dock::TabInteractionStyle {
        outline_color: Color32::TRANSPARENT,
        corner_radius: CornerRadius::ZERO,
        bg_fill: fill,
        text_color: text,
    };
    dock.tab.active = tab(t.surface, t.text);
    dock.tab.focused = tab(t.surface, t.text);
    dock.tab.hovered = tab(t.hover, t.text);
    dock.tab.inactive = tab(t.bg, t.text_muted);
    dock.tab.active_with_kb_focus = tab(t.surface, t.primary);
    dock.tab.focused_with_kb_focus = tab(t.surface, t.primary);
    dock.tab.inactive_with_kb_focus = tab(t.bg, t.primary);
    dock.tab.hline_below_active_tab_name = true;
    dock.tab.tab_body.inner_margin = egui::Margin::ZERO;
    dock.tab.tab_body.bg_fill = t.surface;
    dock.tab.tab_body.stroke = Stroke::NONE;
    dock.tab.tab_body.corner_radius = CornerRadius::ZERO;
    // The viewport's own tab bar is force-hidden every frame
    // (`layout::hide_viewport_tab_bar`); zero these so hiding it leaves
    // no reveal triangle or drag strip in its place (no header, no ▼,
    // no +).
    dock.buttons.show_tab_bar_size = 0.0;
    dock.tab.tab_body.hidden_tab_bar_drag_height = Some(0.0);
    dock
}

#[cfg(test)]
mod tests {
    use super::*;

    /// WCAG 2 contrast ratio of two sRGB colours.
    fn contrast(a: Color32, b: Color32) -> f32 {
        let lum = |c: Color32| {
            let ch = |v: u8| {
                let s = v as f32 / 255.0;
                if s <= 0.04045 {
                    s / 12.92
                } else {
                    ((s + 0.055) / 1.055).powf(2.4)
                }
            };
            0.2126 * ch(c.r()) + 0.7152 * ch(c.g()) + 0.0722 * ch(c.b())
        };
        let (x, y) = (lum(a), lum(b));
        (x.max(y) + 0.05) / (x.min(y) + 0.05)
    }

    #[test]
    fn text_meets_wcag_aa_on_every_surface() {
        for mode in [ThemeMode::Dark, ThemeMode::Light] {
            let t = Tokens::of(mode);
            for (fg, bg, what) in [
                (t.text, t.bg, "text on bg"),
                (t.text, t.surface, "text on surface"),
                (t.text, t.surface_raised, "text on raised"),
                (t.text, t.hover, "text on hover"),
                (t.text_muted, t.surface, "muted text on surface"),
                (t.text_muted, t.surface_raised, "muted text on raised"),
                (t.on_primary, t.primary, "text on primary"),
                (t.primary, t.surface, "primary on surface"),
                (t.danger, t.surface, "danger on surface"),
                (t.success, t.surface, "success on surface"),
                (t.warning, t.surface, "warning on surface"),
            ] {
                let c = contrast(fg, bg);
                assert!(c >= 4.5, "{mode:?}: {what} is {c:.2}:1");
            }
        }
    }
}
