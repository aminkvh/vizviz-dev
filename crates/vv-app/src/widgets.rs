//! The components of DESIGN.md, each with all of its states (hover,
//! pressed, selected, disabled, focus) drawn from `theme`'s tokens, so
//! every panel looks and behaves the same. Anything that needs a new look
//! gets a component here first.

use egui::{
    Color32, CornerRadius, FontId, Response, RichText, Sense, Stroke, StrokeKind, Ui, Vec2,
    WidgetInfo, WidgetType,
};

use crate::theme::{radius, space, text, Tokens, CONTROL_HEIGHT, LARGE_CONTROL_HEIGHT};

/// A button's role (DESIGN.md, "Button").
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Variant {
    /// The one action the panel or dialog is for.
    Primary,
    Secondary,
    /// No fill until hovered: ribbon and toolbar actions.
    Ghost,
    /// Destructive; always confirmed by the caller.
    Danger,
}

/// Fill, border and text of a clickable control in its current state.
fn look(
    ui: &Ui,
    response: &Response,
    variant: Variant,
    selected: bool,
) -> (Color32, Stroke, Color32) {
    let t = Tokens::current(ui.ctx());
    let hovered = response.hovered() && ui.is_enabled();
    let pressed = response.is_pointer_button_down_on();
    if selected || variant == Variant::Primary {
        let fill = if hovered {
            t.primary.gamma_multiply(0.85)
        } else {
            t.primary
        };
        return (fill, Stroke::NONE, t.on_primary);
    }
    let fill = match (variant, hovered) {
        (Variant::Danger, true) => t.danger,
        (_, true) => t.hover,
        (Variant::Ghost | Variant::Danger, false) => Color32::TRANSPARENT,
        (_, false) => t.surface_raised,
    };
    let stroke = if pressed {
        Stroke::new(1.0, t.primary)
    } else if variant == Variant::Secondary {
        Stroke::new(1.0, t.border)
    } else {
        Stroke::NONE
    };
    let fg = match (variant, hovered) {
        (Variant::Danger, true) => t.on_primary,
        (Variant::Danger, false) => t.danger,
        _ => t.text,
    };
    (fill, stroke, fg)
}

/// The 2 px primary ring around a control with keyboard focus.
fn focus_ring(ui: &Ui, response: &Response, corner: u8) {
    if response.has_focus() {
        let t = Tokens::current(ui.ctx());
        ui.painter().rect_stroke(
            response.rect.expand(2.0),
            CornerRadius::same(corner + 2),
            Stroke::new(2.0, t.primary),
            StrokeKind::Outside,
        );
    }
}

/// A regular button: 32 tall, `icon` (may be empty) then `label`.
pub fn button(ui: &mut Ui, icon: &str, label: &str, variant: Variant) -> Response {
    button_selected(ui, icon, label, variant, false)
}

/// A button showing a selected state (a ribbon action that is current).
pub fn button_selected(
    ui: &mut Ui,
    icon: &str,
    label: &str,
    variant: Variant,
    selected: bool,
) -> Response {
    let painter = ui.painter();
    let font = FontId::proportional(text::BODY);
    let glyph = (!icon.is_empty())
        .then(|| painter.layout_no_wrap(icon.into(), font.clone(), Color32::PLACEHOLDER));
    let words = (!label.is_empty())
        .then(|| painter.layout_no_wrap(label.into(), font, Color32::PLACEHOLDER));
    let glyph_w = glyph.as_ref().map_or(0.0, |g| g.size().x);
    let words_w = words.as_ref().map_or(0.0, |w| w.size().x);
    let between = if glyph.is_some() && words.is_some() {
        space::GAP
    } else {
        0.0
    };
    let content = glyph_w + between + words_w;
    let size = Vec2::new(
        (content + 2.0 * space::GAP).max(CONTROL_HEIGHT),
        CONTROL_HEIGHT,
    );
    let (rect, response) = ui.allocate_exact_size(size, Sense::click());
    response.widget_info(|| WidgetInfo::labeled(WidgetType::Button, ui.is_enabled(), label));
    if ui.is_rect_visible(rect) {
        let (fill, stroke, fg) = look(ui, &response, variant, selected);
        let corner = CornerRadius::same(radius::CONTROL);
        ui.painter()
            .rect(rect, corner, fill, stroke, StrokeKind::Inside);
        let left = rect.center().x - content / 2.0;
        if let Some(glyph) = glyph {
            let center = egui::pos2(left + glyph_w / 2.0, rect.center().y);
            paint_centered(ui.painter(), glyph, center, fg);
        }
        if let Some(words) = words {
            let at = egui::pos2(
                left + glyph_w + between,
                rect.center().y - words.size().y / 2.0,
            );
            ui.painter().galley(at, words, fg);
        }
        focus_ring(ui, &response, radius::CONTROL);
    }
    response
}

/// A ribbon tab: selected when open; `contextual` (shown only while it
/// applies, like Trajectory) in the success colour.
pub fn tab(ui: &mut Ui, label: &str, selected: bool, contextual: bool) -> Response {
    let response = button_selected(ui, "", label, Variant::Ghost, selected);
    if contextual && !selected {
        let t = Tokens::current(ui.ctx());
        ui.painter().hline(
            response.rect.x_range().shrink(space::GAP),
            response.rect.bottom() - 2.0,
            Stroke::new(2.0, t.success),
        );
    }
    response
}

/// An icon-only button, 32 x 32; the caller adds the tooltip that says
/// what it does.
pub fn icon_button(ui: &mut Ui, glyph: &str, selected: bool) -> Response {
    icon_button_marked(ui, glyph, selected, false)
}

/// The side of the corner marker's square, which is also its hit area.
pub const CORNER_MARK: f32 = 10.0;

/// The square in a button's bottom-right corner that holds the marker
/// telling it opens a flyout of variants.
pub fn corner_zone(rect: egui::Rect) -> egui::Rect {
    egui::Rect::from_min_max(rect.max - Vec2::splat(CORNER_MARK), rect.max)
}

/// `icon_button`, with a small triangle in its bottom-right corner when
/// `flyout`: the button stands for a group of variants.
pub fn icon_button_marked(ui: &mut Ui, glyph: &str, selected: bool, flyout: bool) -> Response {
    let (rect, response) = ui.allocate_exact_size(Vec2::splat(CONTROL_HEIGHT), Sense::click());
    response
        .widget_info(|| WidgetInfo::selected(WidgetType::Button, ui.is_enabled(), selected, glyph));
    if ui.is_rect_visible(rect) {
        let (fill, stroke, fg) = look(ui, &response, Variant::Ghost, selected);
        ui.painter().rect(
            rect,
            CornerRadius::same(radius::CONTROL),
            fill,
            stroke,
            StrokeKind::Inside,
        );
        let galley =
            ui.painter()
                .layout_no_wrap(glyph.into(), FontId::proportional(text::ICON), fg);
        paint_centered(ui.painter(), galley, rect.center(), fg);
        if flyout {
            paint_corner_triangle(ui.painter(), rect, fg);
        }
        focus_ring(ui, &response, radius::CONTROL);
    }
    response
}

fn paint_corner_triangle(painter: &egui::Painter, rect: egui::Rect, color: Color32) {
    let corner = rect.max - Vec2::splat(3.0);
    painter.add(egui::Shape::convex_polygon(
        vec![
            corner,
            corner - Vec2::new(5.0, 0.0),
            corner - Vec2::new(0.0, 5.0),
        ],
        color,
        Stroke::NONE,
    ));
}

/// Joined icon-only buttons, one of them selected, each with its tooltip.
/// Returns the option clicked.
pub fn icon_segmented(ui: &mut Ui, options: &[(&str, &str)], selected: usize) -> Option<usize> {
    let mut picked = None;
    ui.scope(|ui| {
        ui.spacing_mut().item_spacing.x = 0.0;
        ui.horizontal(|ui| {
            for (i, (glyph, tip)) in options.iter().enumerate() {
                let r = icon_button(ui, glyph, selected == i).on_hover_text(*tip);
                if r.clicked() {
                    picked = Some(i);
                }
            }
        });
    });
    picked
}

/// Paints `galley` with its drawn shape, not its line box, centred on
/// `center`: icon glyphs sit off the text baseline, so line-box centring
/// leaves them visibly high or low.
fn paint_centered(
    painter: &egui::Painter,
    galley: std::sync::Arc<egui::Galley>,
    center: egui::Pos2,
    color: Color32,
) {
    let offset = center - galley.mesh_bounds.center();
    painter.galley(offset.to_pos2(), galley, color);
}

/// A ribbon group's main action: the icon over the label, 56 tall and as
/// wide as its label needs.
pub fn large_button(ui: &mut Ui, icon: &str, label: &str, selected: bool) -> Response {
    let painter = ui.painter();
    let glyph = painter.layout_no_wrap(
        icon.into(),
        FontId::proportional(text::LARGE_ICON),
        Color32::PLACEHOLDER,
    );
    let words = painter.layout_no_wrap(
        label.into(),
        FontId::proportional(text::BODY),
        Color32::PLACEHOLDER,
    );
    let width = glyph.size().x.max(words.size().x) + 2.0 * space::GAP;
    let (rect, response) = ui.allocate_exact_size(
        Vec2::new(width.max(LARGE_CONTROL_HEIGHT), LARGE_CONTROL_HEIGHT),
        Sense::click(),
    );
    response.widget_info(|| WidgetInfo::labeled(WidgetType::Button, ui.is_enabled(), label));
    if ui.is_rect_visible(rect) {
        let (fill, stroke, fg) = look(ui, &response, Variant::Ghost, selected);
        let corner = CornerRadius::same(radius::CONTROL);
        ui.painter()
            .rect(rect, corner, fill, stroke, StrokeKind::Inside);
        let gap = space::TIGHT;
        let top = rect.center().y - (glyph.size().y + gap + words.size().y) / 2.0;
        let x = |g: &egui::Galley| rect.center().x - g.size().x / 2.0;
        ui.painter()
            .galley(egui::pos2(x(&glyph), top), glyph.clone(), fg);
        ui.painter()
            .galley(egui::pos2(x(&words), top + glyph.size().y + gap), words, fg);
        focus_ring(ui, &response, radius::CONTROL);
    }
    response
}

/// A ghost button ending in "▾" (UX_PATTERNS' popover pattern: "Lights ▾",
/// "Clip ▾", "Layout ▾"...): its click opens a popover the caller anchors
/// on the returned response.
pub fn disclosure_button(ui: &mut Ui, icon: &str, label: &str) -> Response {
    button(
        ui,
        icon,
        &format!("{label} {}", egui_phosphor::regular::CARET_DOWN),
        Variant::Ghost,
    )
}

/// An on/off effect with its strength one click away.
/// Returns the switch's response (the caller flips the setting when
/// clicked) and, only while `on`, the small ▾ button's response to
/// anchor a strength popover on.
pub fn switch_disclosure(ui: &mut Ui, on: bool, label: &str) -> (Response, Option<Response>) {
    ui.horizontal(|ui| {
        let sw = switch(ui, on, label);
        let more = on.then(|| icon_button(ui, egui_phosphor::regular::CARET_DOWN, false));
        (sw, more)
    })
    .inner
}

/// An on/off setting: a switch and its label. Returns the response of the
/// whole row; the caller flips the setting when it is clicked.
pub fn switch(ui: &mut Ui, on: bool, label: &str) -> Response {
    let font = FontId::proportional(text::BODY);
    let galley = ui
        .painter()
        .layout_no_wrap(label.into(), font, Color32::PLACEHOLDER);
    let track = Vec2::new(28.0, 16.0);
    let size = Vec2::new(
        track.x + space::GAP + galley.size().x + space::GAP,
        CONTROL_HEIGHT,
    );
    let (rect, response) = ui.allocate_exact_size(size, Sense::click());
    response.widget_info(|| WidgetInfo::selected(WidgetType::Checkbox, ui.is_enabled(), on, label));
    if ui.is_rect_visible(rect) {
        let t = Tokens::current(ui.ctx());
        if response.hovered() && ui.is_enabled() {
            ui.painter()
                .rect_filled(rect, CornerRadius::same(radius::CONTROL), t.hover);
        }
        let track_rect = egui::Rect::from_min_size(
            egui::pos2(rect.left() + space::TIGHT, rect.center().y - track.y / 2.0),
            track,
        );
        let fill = if on { t.primary } else { t.surface_raised };
        ui.painter().rect(
            track_rect,
            CornerRadius::same(8),
            fill,
            Stroke::new(1.0, if on { t.primary } else { t.border }),
            StrokeKind::Inside,
        );
        let knob_x = if on {
            track_rect.right() - track.y / 2.0
        } else {
            track_rect.left() + track.y / 2.0
        };
        let knob = if on { t.on_primary } else { t.text_muted };
        ui.painter()
            .circle_filled(egui::pos2(knob_x, track_rect.center().y), 5.0, knob);
        let at = egui::pos2(
            track_rect.right() + space::GAP,
            rect.center().y - galley.size().y / 2.0,
        );
        ui.painter().galley(at, galley, t.text);
        focus_ring(ui, &response, radius::CONTROL);
    }
    response
}

/// One of a few options as joined buttons, the current one selected;
/// wraps onto more rows when it has to. Returns the option clicked.
pub fn segmented(ui: &mut Ui, options: &[&str], selected: Option<usize>) -> Option<usize> {
    let mut picked = None;
    ui.horizontal_wrapped(|ui| {
        ui.spacing_mut().item_spacing = Vec2::new(0.0, 0.0);
        for (i, name) in options.iter().enumerate() {
            let r = button_selected(ui, "", name, Variant::Ghost, selected == Some(i));
            if r.clicked() {
                picked = Some(i);
            }
        }
    });
    picked
}

/// One of many options as a drop-down with readable names; `placeholder`
/// shows when none is current. Returns the option picked.
pub fn select(
    ui: &mut Ui,
    id: impl std::hash::Hash + std::fmt::Debug,
    options: &[&str],
    selected: Option<usize>,
    placeholder: &str,
) -> Option<usize> {
    let mut picked = None;
    egui::ComboBox::from_id_salt(id)
        .selected_text(selected.map_or(placeholder, |i| options[i]))
        .height(400.0)
        .show_ui(ui, |ui| {
            for (i, name) in options.iter().enumerate() {
                if ui.selectable_label(selected == Some(i), *name).clicked() {
                    picked = Some(i);
                }
            }
        });
    picked
}

/// A "⋯" button whose menu `add_contents` fills; a menu-kind popup, so a
/// nested `ui.menu_button` inside it opens as a submenu.
pub fn menu(ui: &mut Ui, add_contents: impl FnOnce(&mut Ui)) -> Response {
    let r = icon_button(ui, egui_phosphor::regular::DOTS_THREE, false).on_hover_text("More");
    egui::Popup::menu(&r).show(|ui| {
        ui.set_min_width(160.0);
        add_contents(ui);
    });
    r
}

/// A row of a list (a structure, a rep): its text on one line, cut with
/// "…" to the width left after the row's own buttons, which `trailing`
/// adds from the right. Selected rows take the primary fill. Returns the
/// text's response (click to pick the row).
pub fn list_row(
    ui: &mut Ui,
    words: &str,
    selected: bool,
    trailing: impl FnOnce(&mut Ui),
) -> Response {
    ui.horizontal(|ui| {
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            trailing(ui);
            let width = ui.available_width().max(space::PAD);
            let (rect, response) =
                ui.allocate_exact_size(Vec2::new(width, CONTROL_HEIGHT), Sense::click());
            response.widget_info(|| {
                WidgetInfo::selected(
                    WidgetType::SelectableLabel,
                    ui.is_enabled(),
                    selected,
                    words,
                )
            });
            if ui.is_rect_visible(rect) {
                let (fill, stroke, fg) = look(ui, &response, Variant::Ghost, selected);
                ui.painter().rect(
                    rect,
                    CornerRadius::same(radius::CONTROL),
                    fill,
                    stroke,
                    StrokeKind::Inside,
                );
                let mut job = egui::text::LayoutJob::simple_singleline(
                    words.into(),
                    FontId::proportional(text::BODY),
                    fg,
                );
                job.wrap = egui::text::TextWrapping::truncate_at_width(width - 2.0 * space::GAP);
                let galley = ui.painter().layout_job(job);
                let at = rect.left_center() + Vec2::new(space::GAP, -galley.size().y / 2.0);
                ui.painter().galley(at, galley, fg);
                focus_ring(ui, &response, radius::CONTROL);
            }
            response.on_hover_text(words)
        })
        .inner
    })
    .inner
}

/// A list row's "⋯" menu (UX_PATTERNS "Contextual actions"): opens a
/// small popover listing each `(label, danger)` item as a button, top to
/// bottom; returns the index clicked and the "⋯" button's own response
/// (for a caller that anchors a second popover, e.g. a rep row's Options,
/// on the same spot -- `Popup::from_response(..).open_memory(None)`, not
/// `from_toggle_button_response`, or its own click would fight this
/// menu's for the same toggle). `danger` marks a destructive item (e.g.
/// "Close") in the danger colour, DESIGN.md's Danger button.
pub fn menu_button(ui: &mut Ui, items: &[(&str, bool)]) -> (Option<usize>, Response) {
    let r = icon_button(ui, egui_phosphor::regular::DOTS_THREE, false).on_hover_text("More");
    let mut picked = None;
    egui::Popup::from_toggle_button_response(&r)
        .close_behavior(egui::PopupCloseBehavior::CloseOnClick)
        .show(|ui| {
            ui.set_min_width(160.0);
            for (i, (label, danger)) in items.iter().enumerate() {
                let variant = if *danger {
                    Variant::Danger
                } else {
                    Variant::Ghost
                };
                if button(ui, "", label, variant).clicked() {
                    picked = Some(i);
                }
            }
        });
    (picked, r)
}

/// A labelled slider: the label on its own line (never cut off), the
/// slider across the panel with its value, two decimals and `unit`.
pub fn slider(
    ui: &mut Ui,
    label: &str,
    value: &mut f32,
    range: std::ops::RangeInclusive<f32>,
    unit: &str,
) -> Response {
    ui.label(label);
    // The handle is sized from the interact height; a whole control's
    // would make it a slab.
    ui.spacing_mut().interact_size.y = space::PAD + space::TIGHT;
    let value_width = 64.0;
    ui.spacing_mut().slider_width =
        (ui.available_width() - value_width - space::GAP).max(space::WIDEST);
    ui.add(
        egui::Slider::new(value, range)
            .fixed_decimals(2)
            .suffix(unit),
    )
}

/// Whether `slider`'s response reflects an actual drag or keystroke,
/// not just `.changed()`: a `fixed_decimals` slider fires one spurious
/// `changed` the first time an irrational value (e.g. an angle computed
/// via `atan2`) is displayed, from rounding to its shown decimals, with
/// no interaction at all. Callers that only care about deliberate moves
/// (Look ▸ Lights ▾'s "Custom" label) should check this instead.
pub fn slider_moved(response: &Response) -> bool {
    response.changed() && (response.dragged() || response.drag_stopped() || response.has_focus())
}

/// The viewport's top-left mode chip: `words` names the
/// current mode, then a "?" that lists every mouse binding. Returns the
/// label's response (the caller adds a frame-timings hover to it) and
/// the "?" button's (the caller anchors the bindings popover on it).
pub fn mode_chip(ui: &mut Ui, at: egui::Pos2, words: &str) -> (Response, Response) {
    let t = Tokens::current(ui.ctx());
    let area = egui::Rect::from_min_size(
        at + Vec2::splat(space::GAP),
        Vec2::new(420.0, CONTROL_HEIGHT),
    );
    let mut child = ui.new_child(egui::UiBuilder::new().max_rect(area));
    egui::Frame::NONE
        .fill(t.surface.gamma_multiply(0.85))
        .corner_radius(CornerRadius::same(radius::CONTROL))
        .inner_margin(egui::Margin::symmetric(
            space::GAP as i8,
            space::TIGHT as i8,
        ))
        .show(&mut child, |ui| {
            ui.horizontal(|ui| {
                let label = ui.label(RichText::new(words).small().color(t.text_muted));
                let help = ui.small_button(RichText::new("?").small().color(t.text_muted));
                (label, help)
            })
            .inner
        })
        .inner
}

/// What happened to a toast this frame.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum ToastOutcome {
    None,
    Dismissed,
    /// Its action button ran this command (e.g. "undo").
    Action(&'static str),
}

/// The notice bar along the bottom of `rect`: what happened, with an
/// icon in the primary (info) or danger (error) colour, an optional
/// action button (UX_PATTERNS "Feedback": "Undo, Show"), and a button
/// that dismisses it.
pub fn toast(
    ui: &mut Ui,
    rect: egui::Rect,
    message: &str,
    error: bool,
    action: Option<(&'static str, &'static str)>,
) -> ToastOutcome {
    use egui_phosphor::regular as icon;
    let t = Tokens::current(ui.ctx());
    let height = CONTROL_HEIGHT + 2.0 * space::TIGHT;
    let bar = egui::Rect::from_min_max(egui::pos2(rect.min.x, rect.max.y - height), rect.max);
    let mut child = ui.new_child(egui::UiBuilder::new().max_rect(bar));
    let (glyph, color) = if error {
        (icon::WARNING_CIRCLE, t.danger)
    } else {
        (icon::INFO, t.primary)
    };
    egui::Frame::NONE
        .fill(t.surface_raised)
        .stroke(Stroke::new(1.0, color))
        .inner_margin(egui::Margin::symmetric(
            space::GAP as i8,
            space::TIGHT as i8,
        ))
        .show(&mut child, |ui| {
            ui.set_width(bar.width() - 2.0 * space::GAP);
            ui.horizontal(|ui| {
                ui.label(RichText::new(glyph).color(color));
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    let mut outcome = ToastOutcome::None;
                    if button(ui, icon::X, "", Variant::Ghost)
                        .on_hover_text("Dismiss")
                        .clicked()
                    {
                        outcome = ToastOutcome::Dismissed;
                    }
                    if let Some((label, command)) = action {
                        if button(ui, "", label, Variant::Ghost).clicked() {
                            outcome = ToastOutcome::Action(command);
                        }
                    }
                    ui.with_layout(egui::Layout::left_to_right(egui::Align::Center), |ui| {
                        ui.add(egui::Label::new(RichText::new(message).color(t.text)).truncate())
                            .on_hover_text(message);
                    });
                    outcome
                })
                .inner
            })
            .inner
        })
        .inner
}

/// The whole viewport tinted as a drop target, with what dropping does.
pub fn drop_target(ui: &Ui, rect: egui::Rect, words: &str) {
    let t = Tokens::current(ui.ctx());
    let painter = ui.painter_at(rect);
    painter.rect_filled(rect, CornerRadius::ZERO, t.primary.gamma_multiply(0.25));
    painter.rect_stroke(
        rect.shrink(space::GAP),
        CornerRadius::same(radius::CONTAINER),
        Stroke::new(2.0, t.primary),
        StrokeKind::Inside,
    );
    painter.text(
        rect.center(),
        egui::Align2::CENTER_CENTER,
        words,
        FontId::proportional(text::TITLE),
        t.text,
    );
}

/// A form popover's bottom row (UX_PATTERNS "Create"): the primary
/// `verb` (or Enter) submits, "Cancel" (or the caller's own Escape
/// check) backs out. `Some(true)`/`Some(false)` once either fires; the
/// caller runs its command and closes the popover (`ui.close()`) on
/// either answer, same as a plain click would.
pub fn form_actions(ui: &mut Ui, verb: &str) -> Option<bool> {
    let enter = ui.input(|i| i.key_pressed(egui::Key::Enter));
    let mut answer = enter.then_some(true);
    ui.horizontal(|ui| {
        if button(ui, "", verb, Variant::Primary).clicked() {
            answer = Some(true);
        }
        if button(ui, "", "Cancel", Variant::Secondary).clicked() {
            answer = Some(false);
        }
    });
    answer
}

/// What the person did with a dialog this frame.
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Answer {
    Pending,
    Cancel,
    Confirm,
}

/// A modal dialog (DESIGN.md, "Dialog"): the title, `body`, then
/// "Cancel" and the `verb` button bottom-right. Escape cancels; Enter
/// confirms.
pub fn dialog(
    ctx: &egui::Context,
    id: &str,
    title: &str,
    verb: &str,
    variant: Variant,
    body: impl FnOnce(&mut Ui),
) -> Answer {
    match dialog_with(ctx, id, title, &[(verb, variant)], body) {
        Some(Some(_)) => Answer::Confirm,
        Some(None) => Answer::Cancel,
        None => Answer::Pending,
    }
}

/// A dialog with several actions, the last the primary one (Enter picks
/// it): `Some(Some(i))` for action `i`, `Some(None)` for Cancel or
/// Escape, `None` while it is open.
pub fn dialog_with(
    ctx: &egui::Context,
    id: &str,
    title: &str,
    actions: &[(&str, Variant)],
    body: impl FnOnce(&mut Ui),
) -> Option<Option<usize>> {
    let t = Tokens::current(ctx);
    let modal = egui::Modal::new(egui::Id::new(id))
        .frame(
            egui::Frame::NONE
                .fill(t.surface)
                .stroke(Stroke::new(1.0, t.border))
                .corner_radius(CornerRadius::same(radius::CONTAINER))
                .inner_margin(space::WIDE),
        )
        .show(ctx, |ui| {
            ui.set_max_width(420.0);
            ui.label(RichText::new(title).text_style(crate::theme::title_style()));
            ui.add_space(space::INSET);
            body(ui);
            ui.add_space(space::WIDE);
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                let mut picked = None;
                for (i, (label, variant)) in actions.iter().enumerate().rev() {
                    if button(ui, "", label, *variant).clicked() {
                        picked = Some(Some(i));
                    }
                }
                if button(ui, "", "Cancel", Variant::Secondary).clicked() {
                    picked = Some(None);
                }
                picked
            })
            .inner
        });
    let (enter, escape) = ctx.input(|i| {
        (
            i.key_pressed(egui::Key::Enter),
            i.key_pressed(egui::Key::Escape),
        )
    });
    match modal.inner {
        None if escape || modal.should_close() => Some(None),
        None if enter => Some(actions.len().checked_sub(1)),
        answer => answer,
    }
}

/// `slider` for a whole number (a frame index): no decimals.
pub fn slider_int(
    ui: &mut Ui,
    label: &str,
    value: &mut usize,
    range: std::ops::RangeInclusive<usize>,
) -> bool {
    ui.label(label);
    ui.spacing_mut().interact_size.y = space::PAD + space::TIGHT;
    let value_width = 64.0;
    ui.spacing_mut().slider_width =
        (ui.available_width() - value_width - space::GAP).max(space::WIDEST);
    ui.add(egui::Slider::new(value, range)).changed()
}

/// A compass-style direction pad (a Substance/Marmoset-style light
/// gizmo): the angle around a dot is its azimuth (0 = top/forward,
/// clockwise like a clock face), its distance from the centre is the
/// zenith angle (90 - elevation) -- centre is straight up, rim is
/// straight down -- so every direction lands on exactly one spot, unlike
/// an orthographic projection of the light vector (front and back would
/// overlap there). Draws every `lights` entry's dot (`(azimuth,
/// elevation, colour, enabled)`, dim while off), the `selected` one
/// ringed; dragging it returns its new (azimuth, elevation) in degrees.
pub fn direction_gizmo(
    ui: &mut Ui,
    size: f32,
    lights: &[(f32, f32, Color32, bool)],
    selected: usize,
) -> Option<(f32, f32)> {
    let t = Tokens::current(ui.ctx());
    let (rect, response) = ui.allocate_exact_size(Vec2::splat(size), Sense::click_and_drag());
    let center = rect.center();
    let radius = size * 0.5 - 6.0;
    if ui.is_rect_visible(rect) {
        let painter = ui.painter();
        painter.circle_filled(center, radius, t.surface_raised);
        painter.circle_stroke(center, radius, Stroke::new(1.0, t.border));
        painter.circle_stroke(
            center,
            radius * 0.5,
            Stroke::new(1.0, t.border.gamma_multiply(0.6)),
        );
        painter.line_segment(
            [center, center - Vec2::new(0.0, radius)],
            Stroke::new(1.0, t.border),
        );
        for (i, &(az, el, color, enabled)) in lights.iter().enumerate() {
            let p = center + gizmo_offset(az, el, radius);
            let r = if i == selected { 6.0 } else { 4.0 };
            let c = if enabled {
                color
            } else {
                color.gamma_multiply(0.35)
            };
            if i == selected {
                painter.circle_stroke(p, r + 2.0, Stroke::new(2.0, t.primary));
            }
            painter.circle_filled(p, r, c);
            painter.circle_stroke(p, r, Stroke::new(1.0, t.bg));
        }
    }
    response
        .interact_pointer_pos()
        .map(|pos| gizmo_angles(pos - center, radius))
}

/// The gizmo dot's offset from its centre for `(azimuth, elevation)`
/// degrees at `radius`; [`gizmo_angles`] inverts it.
fn gizmo_offset(azimuth: f32, elevation: f32, radius: f32) -> Vec2 {
    let theta = azimuth.to_radians();
    let r = (90.0 - elevation).clamp(0.0, 180.0) / 180.0 * radius;
    Vec2::new(theta.sin(), -theta.cos()) * r
}

/// Inverts [`gizmo_offset`]: the (azimuth, elevation) degrees a dot at
/// `offset` from the gizmo's centre represents, clamped to the disc.
fn gizmo_angles(offset: Vec2, radius: f32) -> (f32, f32) {
    let radius = radius.max(1.0);
    let r = offset.length().min(radius);
    let azimuth = offset.x.atan2(-offset.y).to_degrees();
    let elevation = 90.0 - r / radius * 180.0;
    (azimuth, elevation)
}

/// Key/value rows (DESIGN.md, "Table"): labels muted and right-aligned
/// in their own column, values wrapped and selectable.
pub fn table(ui: &mut Ui, id: &str, rows: &[(&str, String)]) {
    let t = Tokens::current(ui.ctx());
    let font = FontId::proportional(text::BODY);
    let label_width = rows
        .iter()
        .map(|(key, _)| {
            ui.painter()
                .layout_no_wrap((*key).into(), font.clone(), Color32::PLACEHOLDER)
                .size()
                .x
        })
        .fold(0.0, f32::max);
    egui::Grid::new(id)
        .num_columns(2)
        .spacing([space::INSET, space::TIGHT])
        .show(ui, |ui| {
            for (key, value) in rows {
                ui.allocate_ui_with_layout(
                    Vec2::new(label_width, text::BODY * 1.4),
                    egui::Layout::right_to_left(egui::Align::TOP),
                    |ui| ui.label(RichText::new(*key).color(t.text_muted)),
                );
                ui.add(egui::Label::new(value).wrap().selectable(true));
                ui.end_row();
            }
        });
}

/// Colors offered first in every color picker: backgrounds from black to
/// white, then a hue wheel for coloring (scene data, not chrome).
const SWATCHES: [Color32; 16] = [
    Color32::from_rgb(0x00, 0x00, 0x00),
    Color32::from_rgb(0x17, 0x17, 0x1C),
    Color32::from_rgb(0x0B, 0x14, 0x24),
    Color32::from_rgb(0x5A, 0x62, 0x70),
    Color32::from_rgb(0xD2, 0xD8, 0xE5),
    Color32::from_rgb(0xE6, 0xE6, 0xCC),
    Color32::from_rgb(0xF4, 0xF4, 0xF4),
    Color32::from_rgb(0xFF, 0xFF, 0xFF),
    Color32::from_rgb(0xE5, 0x48, 0x4D),
    Color32::from_rgb(0xF2, 0x9A, 0x3A),
    Color32::from_rgb(0xF2, 0xD0, 0x4B),
    Color32::from_rgb(0x46, 0xA7, 0x58),
    Color32::from_rgb(0x2B, 0xB5, 0xA6),
    Color32::from_rgb(0x3E, 0x8E, 0xF7),
    Color32::from_rgb(0x7B, 0x5C, 0xE6),
    Color32::from_rgb(0xD4, 0x5C, 0xB8),
];

/// A color, as its swatch; a click opens common colors, a full picker
/// and a hex field. True when the color changed.
pub fn color_picker(ui: &mut Ui, color: &mut Color32) -> bool {
    let t = Tokens::current(ui.ctx());
    let (rect, response) = ui.allocate_exact_size(Vec2::splat(CONTROL_HEIGHT), Sense::click());
    if ui.is_rect_visible(rect) {
        let swatch = rect.shrink(space::TIGHT);
        let stroke = if response.hovered() {
            t.primary
        } else {
            t.border
        };
        ui.painter().rect(
            swatch,
            CornerRadius::same(radius::CONTROL),
            *color,
            Stroke::new(1.0, stroke),
            StrokeKind::Outside,
        );
        focus_ring(ui, &response, radius::CONTROL);
    }
    let mut changed = false;
    egui::Popup::from_toggle_button_response(&response)
        .close_behavior(egui::PopupCloseBehavior::CloseOnClickOutside)
        .show(|ui| changed = color_editor(ui, response.id, color));
    changed
}

/// Common colors, a full picker and a hex field, drawn in place (inside a
/// popover, where `color_picker`'s own popup would close its parent).
/// True when the color changed.
pub fn color_editor(ui: &mut Ui, id: egui::Id, color: &mut Color32) -> bool {
    const COLUMNS: f32 = 8.0;
    let width = COLUMNS * space::WIDE + (COLUMNS - 1.0) * space::TIGHT;
    ui.set_width(width);
    let mut changed = swatch_grid(ui, color);
    ui.add_space(space::GAP);
    changed |= hsv_picker(ui, id, color, width);
    ui.add_space(space::GAP);
    changed |= hex_field(ui, id.with("hex"), color, width);
    changed
}

fn swatch_grid(ui: &mut Ui, color: &mut Color32) -> bool {
    let t = Tokens::current(ui.ctx());
    let mut changed = false;
    ui.horizontal_wrapped(|ui| {
        ui.spacing_mut().item_spacing = Vec2::splat(space::TIGHT);
        for swatch in SWATCHES {
            let (r, pick) = ui.allocate_exact_size(Vec2::splat(space::WIDE), Sense::click());
            let current = swatch == *color;
            ui.painter().rect(
                r,
                CornerRadius::same(radius::CONTROL),
                swatch,
                Stroke::new(
                    if current { 2.0 } else { 1.0 },
                    if current { t.primary } else { t.border },
                ),
                StrokeKind::Inside,
            );
            if pick.clicked() {
                *color = swatch;
                changed = true;
            }
        }
    });
    changed
}

/// A saturation/value field over a hue bar, both `width` wide. The hue is
/// kept across frames under `id`, since grey and black carry none.
fn hsv_picker(ui: &mut Ui, id: egui::Id, color: &mut Color32, width: f32) -> bool {
    use egui::ecolor::Hsva;
    let mut hsva = ui
        .data(|d| d.get_temp::<Hsva>(id))
        .filter(|h| h.to_srgb() == [color.r(), color.g(), color.b()])
        .unwrap_or_else(|| Hsva::from(*color));
    let before = hsva;
    let (field, drag) =
        ui.allocate_exact_size(Vec2::new(width, width * 0.6), Sense::click_and_drag());
    if let Some(p) = drag.interact_pointer_pos() {
        hsva.s = ((p.x - field.left()) / field.width()).clamp(0.0, 1.0);
        hsva.v = 1.0 - ((p.y - field.top()) / field.height()).clamp(0.0, 1.0);
    }
    ui.add_space(space::TIGHT);
    let (bar, hue_drag) =
        ui.allocate_exact_size(Vec2::new(width, space::PAD), Sense::click_and_drag());
    if let Some(p) = hue_drag.interact_pointer_pos() {
        hsva.h = ((p.x - bar.left()) / bar.width()).clamp(0.0, 1.0);
    }
    paint_sv_field(ui.painter(), field, hsva);
    paint_hue_bar(ui.painter(), bar, hsva.h);
    ui.data_mut(|d| d.insert_temp(id, hsva));
    if hsva == before {
        return false;
    }
    let [r, g, b] = hsva.to_srgb();
    *color = Color32::from_rgb(r, g, b);
    true
}

fn paint_sv_field(painter: &egui::Painter, rect: egui::Rect, hsva: egui::ecolor::Hsva) {
    let hue = egui::ecolor::Hsva::new(hsva.h, 1.0, 1.0, 1.0);
    let mut mesh = egui::Mesh::default();
    // White to the pure hue across, then transparent to black down.
    for (layer, [tl, tr, bl, br]) in [
        [
            Color32::WHITE,
            Color32::from(hue),
            Color32::WHITE,
            Color32::from(hue),
        ],
        [
            Color32::TRANSPARENT,
            Color32::TRANSPARENT,
            Color32::BLACK,
            Color32::BLACK,
        ],
    ]
    .into_iter()
    .enumerate()
    {
        let base = (layer * 4) as u32;
        for (pos, c) in [
            (rect.left_top(), tl),
            (rect.right_top(), tr),
            (rect.left_bottom(), bl),
            (rect.right_bottom(), br),
        ] {
            mesh.colored_vertex(pos, c);
        }
        mesh.add_triangle(base, base + 1, base + 2);
        mesh.add_triangle(base + 1, base + 3, base + 2);
    }
    painter.add(mesh);
    let at = egui::pos2(
        rect.left() + hsva.s * rect.width(),
        rect.top() + (1.0 - hsva.v) * rect.height(),
    );
    let ring = if hsva.v > 0.5 {
        Color32::BLACK
    } else {
        Color32::WHITE
    };
    painter.circle_stroke(at, space::TIGHT + 1.0, Stroke::new(1.5, ring));
}

fn paint_hue_bar(painter: &egui::Painter, rect: egui::Rect, hue: f32) {
    const STOPS: usize = 6;
    let mut mesh = egui::Mesh::default();
    for i in 0..=STOPS {
        let h = i as f32 / STOPS as f32;
        let c = Color32::from(egui::ecolor::Hsva::new(h, 1.0, 1.0, 1.0));
        let x = rect.left() + h * rect.width();
        mesh.colored_vertex(egui::pos2(x, rect.top()), c);
        mesh.colored_vertex(egui::pos2(x, rect.bottom()), c);
        if i > 0 {
            let k = (2 * i) as u32;
            mesh.add_triangle(k - 2, k - 1, k);
            mesh.add_triangle(k - 1, k + 1, k);
        }
    }
    painter.add(mesh);
    let x = rect.left() + hue * rect.width();
    painter.vline(x, rect.y_range(), Stroke::new(2.0, Color32::WHITE));
}

fn hex_field(ui: &mut Ui, id: egui::Id, color: &mut Color32, width: f32) -> bool {
    let mut hex = ui
        .data_mut(|d| d.get_temp::<String>(id))
        .unwrap_or_else(|| self::hex(*color));
    let edit = ui.add(egui::TextEdit::singleline(&mut hex).desired_width(width));
    let mut changed = false;
    if edit.changed() {
        if let Some(parsed) = parse_hex(&hex) {
            *color = parsed;
            changed = true;
        }
        ui.data_mut(|d| d.insert_temp(id, hex));
    } else if !edit.has_focus() {
        ui.data_mut(|d| d.remove::<String>(id));
    }
    changed
}

/// `#RRGGBB` (the `#` optional) as a color.
fn parse_hex(text: &str) -> Option<Color32> {
    let digits = text.trim().trim_start_matches('#');
    if digits.len() != 6 {
        return None;
    }
    let v = u32::from_str_radix(digits, 16).ok()?;
    Some(Color32::from_rgb((v >> 16) as u8, (v >> 8) as u8, v as u8))
}

/// A color as the command language writes it.
pub fn hex(color: Color32) -> String {
    format!("#{:02X}{:02X}{:02X}", color.r(), color.g(), color.b())
}

/// A panel section's heading.
pub fn section(ui: &mut Ui, title: &str) {
    ui.add_space(space::TIGHT);
    ui.label(
        RichText::new(title)
            .text_style(egui::TextStyle::Heading)
            .strong(),
    );
}

/// A caption under a field or group: small and muted.
pub fn caption(ui: &mut Ui, words: &str) -> Response {
    let t = Tokens::current(ui.ctx());
    ui.label(RichText::new(words).small().color(t.text_muted))
}

/// What goes here when nothing does yet: an icon, one line, and the one
/// action that fills it. True when the action is clicked.
pub fn empty_state(ui: &mut Ui, icon: &str, line: &str, action: &str) -> bool {
    let t = Tokens::current(ui.ctx());
    let mut clicked = false;
    ui.vertical_centered(|ui| {
        ui.add_space(space::WIDE);
        ui.label(
            RichText::new(icon)
                .size(text::TITLE * 2.0)
                .color(t.text_muted),
        );
        ui.add_space(space::GAP);
        ui.label(RichText::new(line).color(t.text_muted));
        ui.add_space(space::INSET);
        clicked = button(ui, "", action, Variant::Primary).clicked();
        ui.add_space(space::WIDE);
    });
    clicked
}

/// `n` things, in words: "1 atom", "4,779 atoms" (never "atom(s)").
pub fn count(n: usize, noun: &str) -> String {
    let digits = n.to_string();
    let mut grouped = String::new();
    for (i, c) in digits.chars().enumerate() {
        if i > 0 && (digits.len() - i) % 3 == 0 {
            grouped.push(',');
        }
        grouped.push(c);
    }
    let plural = if n == 1 { "" } else { "s" };
    format!("{grouped} {noun}{plural}")
}

#[cfg(test)]
mod tests {
    use super::{count, gizmo_angles, gizmo_offset};

    #[test]
    fn counts_read_as_words() {
        assert_eq!(count(1, "atom"), "1 atom");
        assert_eq!(count(0, "atom"), "0 atoms");
        assert_eq!(count(4779, "atom"), "4,779 atoms");
        assert_eq!(count(3_968_189, "atom"), "3,968,189 atoms");
    }

    #[test]
    fn gizmo_centre_is_straight_up_and_rim_is_straight_down() {
        assert!(gizmo_offset(0.0, 90.0, 100.0).length() < 1e-4);
        let rim = gizmo_offset(0.0, -90.0, 100.0);
        assert!((rim.length() - 100.0).abs() < 1e-3);
    }

    #[test]
    fn gizmo_azimuth_reads_clockwise_from_the_top() {
        // Straight right (az 90, elevation 0 -- the equator, half radius
        // out) should land on the gizmo's +x side, level with centre.
        let right = gizmo_offset(90.0, 0.0, 100.0);
        assert!(
            (right.x - 50.0).abs() < 1e-3 && right.y.abs() < 1e-3,
            "{right:?}"
        );
    }

    #[test]
    fn gizmo_angles_round_trip_through_its_offset() {
        for &(az, el) in &[
            (0.0, 0.0),
            (45.0, 30.0),
            (-90.0, -45.0),
            (179.0, 45.0),
            (-179.0, -45.0),
            (30.0, 89.0),
        ] {
            let offset = gizmo_offset(az, el, 100.0);
            let (az2, el2) = gizmo_angles(offset, 100.0);
            assert!((el - el2).abs() < 1e-2, "elevation {el} -> {el2}");
            let d = (az - az2).abs();
            assert!(d < 1e-1 || d > 359.9, "azimuth {az} -> {az2}");
        }
    }
}
