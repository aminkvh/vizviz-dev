//! Represent ▸ Coloring: one popover holding every way to color the
//! current rep. Schemes by what the structure is, properties with a
//! legend and an editable ramp and range, one solid color, and the
//! per-target overrides over any of them. Every click runs a `color`
//! command, so it is undoable, logged and scriptable like a typed one.

use egui::{Color32, Rect, RichText, Sense, Ui, Vec2};
use vv_render::color::ramp_color;
use vv_render::coloring::{legend, Legend};
use vv_scene::{ColorOverride, ColorScheme, Property, PropertyKind, Ramp};

use crate::ribbon::{base_command, with_current_override, COLORINGS, COLORING_GROUPS};
use crate::theme::{space, text, Tokens};
use crate::ui::AppUi;
use crate::widgets::{self, Variant};

const WIDTH: f32 = 360.0;
const MAX_HEIGHT: f32 = 680.0;

/// What the popover shows, copied out of the scene so its widgets can
/// queue commands while it is drawn.
struct Snapshot {
    coloring: ColorScheme,
    legend: Option<Legend>,
    channels: Vec<String>,
    overrides: Vec<ColorOverride>,
    glycan: bool,
}

fn snapshot(app: &AppUi<'_>) -> Option<Snapshot> {
    let loaded = app.scene.structure(app.current()?)?;
    let rep = loaded.rep();
    Some(Snapshot {
        coloring: rep.coloring.clone(),
        legend: legend(loaded, &rep.coloring),
        channels: loaded.values.keys().cloned().collect(),
        overrides: loaded.color_overrides.clone(),
        glycan: rep.representation == vv_scene::Representation::Glycan,
    })
}

pub(crate) fn coloring_popover(app: &mut AppUi<'_>, ui: &mut Ui) {
    ui.set_width(WIDTH);
    let Some(shot) = snapshot(app) else {
        widgets::caption(ui, "Open a structure to color it");
        return;
    };
    let mut lines = Vec::new();
    egui::ScrollArea::vertical()
        .max_height(MAX_HEIGHT)
        .show(ui, |ui| {
            if shot.glycan {
                widgets::caption(ui, "SNFG colours are fixed");
            } else {
                schemes_section(ui, &shot, &mut lines);
                properties_section(ui, &shot, &mut lines);
                uniform_section(ui, &shot, &mut lines);
            }
            overrides_section(ui, &shot, &mut lines);
        });
    for line in lines {
        let line = with_current_override(line, app);
        app.run_command_logged(&line);
    }
}

fn section(ui: &mut Ui, title: &str, open: bool, body: impl FnOnce(&mut Ui)) {
    egui::CollapsingHeader::new(RichText::new(title).strong())
        .id_salt(("coloring-section", title))
        .default_open(open)
        .show(ui, body);
}

/// `items` as toggle buttons in `columns`, the one whose command is `now`
/// selected; a click queues its command.
fn chips(ui: &mut Ui, items: &[(&str, &str)], now: &str, columns: usize, lines: &mut Vec<String>) {
    egui::Grid::new(ui.id().with(items.as_ptr() as usize))
        .spacing(Vec2::splat(space::TIGHT))
        .show(ui, |ui| {
            for (i, (label, command)) in items.iter().enumerate() {
                let r = ui
                    .selectable_label(*command == now, *label)
                    .on_hover_text(*command);
                if r.clicked() {
                    lines.push((*command).into());
                }
                if i % columns == columns - 1 {
                    ui.end_row();
                }
            }
        });
}

fn group_items(group: usize) -> &'static [(&'static str, &'static str)] {
    let start = COLORING_GROUPS[group].1;
    let end = COLORING_GROUPS
        .get(group + 1)
        .map_or(COLORINGS.len(), |(_, s)| *s);
    &COLORINGS[start..end]
}

fn schemes_section(ui: &mut Ui, shot: &Snapshot, lines: &mut Vec<String>) {
    section(
        ui,
        COLORING_GROUPS[0].0,
        !matches!(shot.coloring, ColorScheme::Property(_)),
        |ui| {
            chips(ui, group_items(0), &base_command(&shot.coloring), 3, lines);
        },
    );
}

fn properties_section(ui: &mut Ui, shot: &Snapshot, lines: &mut Vec<String>) {
    let now = base_command(&shot.coloring);
    section(
        ui,
        COLORING_GROUPS[1].0,
        matches!(shot.coloring, ColorScheme::Property(_)),
        |ui| {
            chips(ui, group_items(1), &now, 2, lines);
            channel_chips(ui, shot, lines);
            if let (ColorScheme::Property(p), Some(legend)) = (&shot.coloring, &shot.legend) {
                ui.add_space(space::GAP);
                property_controls(ui, p, legend, lines);
            }
        },
    );
}

/// Solvent-accessible area (computed on first use) and every attached
/// value channel, as colorings.
fn channel_chips(ui: &mut Ui, shot: &Snapshot, lines: &mut Vec<String>) {
    let current = match &shot.coloring {
        ColorScheme::Property(Property {
            kind: PropertyKind::Values(name),
            ..
        }) => Some(name.as_str()),
        _ => None,
    };
    ui.add_space(space::TIGHT);
    ui.horizontal_wrapped(|ui| {
        let sasa = ui
            .selectable_label(current == Some("sasa"), "Solvent-accessible area")
            .on_hover_text("sasa, then color values sasa");
        if sasa.clicked() {
            if !shot.channels.iter().any(|c| c == "sasa") {
                lines.push("sasa".into());
            }
            lines.push("color values sasa".into());
        }
        for name in shot.channels.iter().filter(|c| *c != "sasa") {
            let r = ui
                .selectable_label(current == Some(name), format!("Values: {name}"))
                .on_hover_text(format!("color values {name}"));
            if r.clicked() {
                lines.push(format!("color values {name}"));
            }
        }
    });
}

fn property_line(p: &Property) -> String {
    format!("color {}", p.name())
}

/// The legend, then the range's two ends and the ramp, each editable.
fn property_controls(ui: &mut Ui, p: &Property, legend: &Legend, lines: &mut Vec<String>) {
    ui.label(RichText::new(&legend.title).strong());
    ui.add_space(space::TIGHT);
    legend_bar(ui, legend, ui.available_width());
    ui.horizontal(|ui| {
        if let Some(range) = range_fields(ui, p, legend) {
            let mut edited = p.clone();
            edited.range = Some(range);
            lines.push(property_line(&edited));
        }
    });
    ui.horizontal(|ui| {
        ui.label("Ramp");
        if let Some(ramp) = ramp_select(ui, legend.ramp) {
            let mut edited = p.clone();
            edited.ramp = Some(ramp);
            lines.push(property_line(&edited));
        }
        let tuned = p.range.is_some() || p.ramp.is_some();
        if tuned && widgets::button(ui, "", "Reset", Variant::Secondary).clicked() {
            lines.push(format!("color {}", p.kind.name()));
        }
    });
    legend_overlay_switch(ui);
}

fn ramp_select(ui: &mut Ui, current: Ramp) -> Option<Ramp> {
    let labels: Vec<&str> = Ramp::ALL.iter().map(|r| r.label()).collect();
    let at = Ramp::ALL.iter().position(|r| *r == current);
    widgets::select(ui, "coloring-ramp", &labels, at, "Ramp").map(|i| Ramp::ALL[i])
}

/// The range's two drag fields; the range once an edit ends. The value
/// being edited is kept across frames under the property's own name, so a
/// drag does not snap back while the command waits for its end.
fn range_fields(ui: &mut Ui, p: &Property, legend: &Legend) -> Option<[f32; 2]> {
    let id = ui.id().with(("range", p.name()));
    let [mut lo, mut hi] = ui
        .data(|d| d.get_temp::<[f32; 2]>(id))
        .unwrap_or([legend.min, legend.max]);
    let speed = f64::from((hi - lo).abs() / 200.0).max(0.001);
    let edit = |ui: &mut Ui, value: &mut f32| {
        let r = ui.add(egui::DragValue::new(value).speed(speed).max_decimals(3));
        (r.changed(), r.drag_stopped() || r.lost_focus())
    };
    let (lo_changed, lo_done) = edit(ui, &mut lo);
    widgets::caption(ui, legend.low);
    let (hi_changed, hi_done) = ui
        .with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            let ends = edit(ui, &mut hi);
            widgets::caption(ui, legend.high);
            ends
        })
        .inner;
    if lo_done || hi_done {
        ui.data_mut(|d| d.remove::<[f32; 2]>(id));
        return ([lo, hi] != [legend.min, legend.max]).then_some([lo, hi]);
    }
    if lo_changed || hi_changed {
        ui.data_mut(|d| d.insert_temp(id, [lo, hi]));
    }
    None
}

fn unpack(packed: u32) -> Color32 {
    Color32::from_rgb(packed as u8, (packed >> 8) as u8, (packed >> 16) as u8)
}

/// `legend`'s ramp as a horizontal gradient in `rect`.
fn paint_gradient(painter: &egui::Painter, rect: Rect, legend: &Legend) {
    const STEPS: usize = 48;
    let mut mesh = egui::Mesh::default();
    for i in 0..=STEPS {
        let t = i as f32 / STEPS as f32;
        let color = unpack(ramp_color(legend.ramp, t));
        let x = rect.left() + t * rect.width();
        mesh.colored_vertex(egui::pos2(x, rect.top()), color);
        mesh.colored_vertex(egui::pos2(x, rect.bottom()), color);
        if i > 0 {
            let k = (2 * i) as u32;
            mesh.add_triangle(k - 2, k - 1, k);
            mesh.add_triangle(k - 1, k + 1, k);
        }
    }
    painter.add(mesh);
}

fn legend_bar(ui: &mut Ui, legend: &Legend, width: f32) {
    let (rect, _) = ui.allocate_exact_size(Vec2::new(width, space::PAD), Sense::hover());
    paint_gradient(ui.painter(), rect, legend);
    let border = Tokens::current(ui.ctx()).border;
    ui.painter().rect_stroke(
        rect,
        0.0,
        egui::Stroke::new(1.0, border),
        egui::StrokeKind::Outside,
    );
}

const LEGEND_OVERLAY: &str = "coloring-legend-overlay";

fn legend_overlay_on(ctx: &egui::Context) -> bool {
    ctx.data(|d| d.get_temp(egui::Id::new(LEGEND_OVERLAY)))
        .unwrap_or(true)
}

fn legend_overlay_switch(ui: &mut Ui) {
    let on = legend_overlay_on(ui.ctx());
    ui.add_space(space::TIGHT);
    if widgets::switch(ui, on, "Legend in the viewport").clicked() {
        ui.data_mut(|d| d.insert_temp(egui::Id::new(LEGEND_OVERLAY), !on));
    }
}

fn uniform_section(ui: &mut Ui, shot: &Snapshot, lines: &mut Vec<String>) {
    let solid = matches!(shot.coloring, ColorScheme::Constant(_));
    section(ui, "Uniform", solid, |ui| {
        let id = egui::Id::new("coloring-uniform");
        let mut color = match shot.coloring {
            ColorScheme::Constant([r, g, b]) => Color32::from_rgb(r, g, b),
            _ => ui
                .data(|d| d.get_temp::<Color32>(id))
                .unwrap_or(Color32::GRAY),
        };
        if widgets::color_editor(ui, id.with("editor"), &mut color) {
            ui.data_mut(|d| d.insert_temp(id, color));
            lines.push(format!("color {}", widgets::hex(color)));
        }
    });
}

/// (label, command word, hint) of each override target kind.
const TARGETS: [(&str, &str, &str); 5] = [
    ("Element", "element", "C"),
    ("Atom name", "name", "CA"),
    ("Residue name", "resname", "HEM"),
    ("Atom index", "atom", "123"),
    ("Selection", "sel", "chain A and helix"),
];

/// The add form's state, kept across frames.
#[derive(Clone)]
struct OverrideForm {
    kind: usize,
    value: String,
    color: Color32,
    picking: bool,
}

impl Default for OverrideForm {
    fn default() -> Self {
        Self {
            kind: 1,
            value: String::new(),
            color: Color32::from_rgb(0xFF, 0x00, 0xFF),
            picking: false,
        }
    }
}

fn overrides_section(ui: &mut Ui, shot: &Snapshot, lines: &mut Vec<String>) {
    section(ui, "Overrides", !shot.overrides.is_empty(), |ui| {
        widgets::caption(ui, "Colors laid over the coloring on every rep");
        for o in &shot.overrides {
            override_row(ui, o, lines);
        }
        if shot.overrides.len() > 1
            && widgets::button(ui, "", "Remove all", Variant::Ghost).clicked()
        {
            lines.push("color unset all".into());
        }
        ui.add_space(space::TIGHT);
        override_form(ui, lines);
    });
}

fn override_row(ui: &mut Ui, o: &ColorOverride, lines: &mut Vec<String>) {
    let [r, g, b] = o.color;
    ui.horizontal(|ui| {
        let (swatch, _) = ui.allocate_exact_size(Vec2::splat(space::PAD), Sense::hover());
        let border = Tokens::current(ui.ctx()).border;
        ui.painter().rect(
            swatch,
            egui::CornerRadius::same(crate::theme::radius::CONTROL),
            Color32::from_rgb(r, g, b),
            egui::Stroke::new(1.0, border),
            egui::StrokeKind::Inside,
        );
        ui.label(o.target.to_string());
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            let remove =
                widgets::icon_button(ui, egui_phosphor::regular::X, false).on_hover_text("Remove");
            if remove.clicked() {
                lines.push(format!("color unset {}", o.target));
            }
        });
    });
}

fn override_form(ui: &mut Ui, lines: &mut Vec<String>) {
    let id = egui::Id::new("coloring-override-form");
    let mut form: OverrideForm = ui.data(|d| d.get_temp(id)).unwrap_or_default();
    let labels: Vec<&str> = TARGETS.iter().map(|t| t.0).collect();
    ui.horizontal(|ui| {
        if let Some(i) = widgets::select(ui, "override-kind", &labels, Some(form.kind), "Target") {
            form.kind = i;
        }
        ui.add(
            egui::TextEdit::singleline(&mut form.value)
                .hint_text(TARGETS[form.kind].2)
                .desired_width(f32::INFINITY),
        );
    });
    ui.horizontal(|ui| {
        let (rect, swatch) =
            ui.allocate_exact_size(Vec2::splat(crate::theme::CONTROL_HEIGHT), Sense::click());
        let border = Tokens::current(ui.ctx()).border;
        ui.painter().rect(
            rect.shrink(space::TIGHT),
            egui::CornerRadius::same(crate::theme::radius::CONTROL),
            form.color,
            egui::Stroke::new(1.0, border),
            egui::StrokeKind::Outside,
        );
        if swatch.on_hover_text("Choose the color").clicked() {
            form.picking = !form.picking;
        }
        ui.label(RichText::new(widgets::hex(form.color)).size(text::MONO));
        let ready = !form.value.trim().is_empty();
        let add = ui.add_enabled_ui(ready, |ui| widgets::button(ui, "", "Add", Variant::Primary));
        if add.inner.clicked() {
            let word = TARGETS[form.kind].1;
            let color = widgets::hex(form.color);
            lines.push(format!("color set {word} {} {color}", form.value.trim()));
            form.value.clear();
        }
    });
    if form.picking {
        widgets::color_editor(ui, id.with("editor"), &mut form.color);
    }
    ui.data_mut(|d| d.insert_temp(id, form));
}

/// The legend of the current rep's property coloring in the viewport's
/// bottom-right corner, unless switched off in the popover.
pub(crate) fn legend_overlay(app: &AppUi<'_>, ui: &Ui, viewport: Rect, background: Color32) {
    if !legend_overlay_on(ui.ctx()) {
        return;
    }
    let light = background.r() as u32 + background.g() as u32 + background.b() as u32 > 3 * 128;
    let ink = Tokens::of(if light {
        crate::theme::ThemeMode::Light
    } else {
        crate::theme::ThemeMode::Dark
    })
    .text;
    let Some(legend) = app
        .current()
        .and_then(|id| app.scene.structure(id))
        .and_then(|loaded| legend(loaded, &loaded.rep().coloring))
    else {
        return;
    };
    let size = Vec2::new(160.0, 10.0);
    let bar = Rect::from_min_size(
        egui::pos2(viewport.right() - 16.0 - size.x, viewport.bottom() - 28.0),
        size,
    );
    let painter = ui.painter_at(viewport);
    paint_gradient(&painter, bar, &legend);
    painter.rect_stroke(
        bar,
        0.0,
        egui::Stroke::new(1.0, ink),
        egui::StrokeKind::Outside,
    );
    let font = egui::FontId::proportional(text::CAPTION);
    painter.text(
        bar.left_top() + Vec2::new(0.0, -4.0),
        egui::Align2::LEFT_BOTTOM,
        &legend.title,
        font.clone(),
        ink,
    );
    painter.text(
        bar.left_bottom() + Vec2::new(0.0, 3.0),
        egui::Align2::LEFT_TOP,
        format_end(legend.min),
        font.clone(),
        ink,
    );
    painter.text(
        bar.right_bottom() + Vec2::new(0.0, 3.0),
        egui::Align2::RIGHT_TOP,
        format_end(legend.max),
        font,
        ink,
    );
}

/// A legend end with as many decimals as it needs, at most two.
fn format_end(value: f32) -> String {
    let text = format!("{value:.2}");
    text.trim_end_matches('0').trim_end_matches('.').to_owned()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn legend_ends_drop_needless_zeros() {
        assert_eq!(format_end(4.5), "4.5");
        assert_eq!(format_end(-3.0), "-3");
        assert_eq!(format_end(0.126), "0.13");
        assert_eq!(format_end(0.0), "0");
    }

    #[test]
    fn every_coloring_group_is_reachable_from_the_popover() {
        let shown: usize = (0..COLORING_GROUPS.len())
            .map(|g| group_items(g).len())
            .sum();
        assert_eq!(shown, COLORINGS.len(), "each button belongs to one group");
        let commands: Vec<&str> = COLORINGS.iter().map(|(_, c)| *c).collect();
        let reached = |scheme: &ColorScheme| commands.contains(&base_command(scheme).as_str());
        for scheme in ColorScheme::plain() {
            assert!(reached(&scheme), "{scheme:?} has no button");
        }
        for kind in PropertyKind::builtin() {
            assert!(
                reached(&ColorScheme::from(kind.clone())),
                "{kind:?} has no button"
            );
        }
    }
}
