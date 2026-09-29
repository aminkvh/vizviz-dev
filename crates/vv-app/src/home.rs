//! The Home tab's selection controls, each a group body of captioned icon
//! cells: the select-tool button with its flyout, the pick-level choice
//! and the quick-select grid. Each returns the command line to run (so a
//! click here is also a typeable command) and fills the ribbon's body
//! height itself.

use egui::{Id, PopupCloseBehavior, Ui};
use egui_phosphor::regular as icon;

use crate::keys::shortcut_for;
use crate::select_tool::{classify_press, SelectLevel, SelectShape, ToolPress};
use crate::ui::AppUi;
use crate::widgets::{self, Variant};

struct ShapeItem {
    shape: SelectShape,
    icon: &'static str,
    label: &'static str,
    tip: &'static str,
    command: &'static str,
}

const SHAPES: [ShapeItem; 4] = [
    ShapeItem {
        shape: SelectShape::Click,
        icon: icon::CURSOR,
        label: "Click",
        tip: "Click an atom to select it; dragging orbits",
        command: "selectmode shape click",
    },
    ShapeItem {
        shape: SelectShape::Box,
        icon: icon::SELECTION,
        label: "Box",
        tip: "Drag a rectangle around what to select",
        command: "selectmode shape box",
    },
    ShapeItem {
        shape: SelectShape::Circle,
        icon: icon::CIRCLE_DASHED,
        label: "Circle",
        tip: "Drag a circle out from its center",
        command: "selectmode shape circle",
    },
    ShapeItem {
        shape: SelectShape::Lasso,
        icon: icon::LASSO,
        label: "Lasso",
        tip: "Draw a freehand loop around what to select",
        command: "selectmode shape lasso",
    },
];

struct LevelItem {
    level: SelectLevel,
    label: &'static str,
    icon: &'static str,
    tip: &'static str,
    command: &'static str,
}

const LEVELS: [LevelItem; 4] = [
    LevelItem {
        level: SelectLevel::Atom,
        label: "Atom",
        icon: icon::DOT_OUTLINE,
        tip: "Pick atoms",
        command: "selectmode level atom",
    },
    LevelItem {
        level: SelectLevel::Residue,
        label: "Residue",
        icon: icon::HEXAGON,
        tip: "Pick whole residues",
        command: "selectmode level residue",
    },
    LevelItem {
        level: SelectLevel::Chain,
        label: "Chain",
        icon: icon::LINK_SIMPLE,
        tip: "Pick whole chains",
        command: "selectmode level chain",
    },
    LevelItem {
        level: SelectLevel::Molecule,
        label: "Molecule",
        icon: icon::FLASK,
        tip: "Pick whole molecules (connected atoms)",
        command: "selectmode level molecule",
    },
];

struct QuickItem {
    icon: &'static str,
    label: &'static str,
    name: &'static str,
    replace: &'static str,
    add: &'static str,
}

const QUICK: [QuickItem; 8] = [
    QuickItem {
        icon: icon::SPIRAL,
        label: "Protein",
        name: "protein",
        replace: "select protein",
        add: "select add protein",
    },
    QuickItem {
        icon: icon::PILL,
        label: "Ligand",
        name: "ligands",
        replace: "select ligand",
        add: "select add ligand",
    },
    QuickItem {
        icon: icon::DROP,
        label: "Water",
        name: "water",
        replace: "select water",
        add: "select add water",
    },
    QuickItem {
        icon: icon::LIGHTNING,
        label: "Ion",
        name: "ions",
        replace: "select ion",
        add: "select add ion",
    },
    QuickItem {
        icon: icon::DROP_HALF,
        label: "Lipid",
        name: "lipids",
        replace: "select lipid",
        add: "select add lipid",
    },
    QuickItem {
        icon: icon::DNA,
        label: "Nucleic",
        name: "nucleic acids",
        replace: "select nucleic",
        add: "select add nucleic",
    },
    QuickItem {
        icon: icon::SHAPES,
        label: "Glycan",
        name: "glycans",
        replace: "select glycan",
        add: "select add glycan",
    },
    QuickItem {
        icon: icon::TEST_TUBE,
        label: "Additive",
        name: "additives",
        replace: "select additive",
        add: "select add additive",
    },
];

const INVERT: &str = "select invert";
const CLEAR: &str = "clear";
const CYCLE_SHAPE: &str = "selectmode shape next";

/// Every command the tool group can run; the first is what its key tip
/// runs.
pub fn tool_commands() -> Vec<&'static str> {
    ["selectmode shape toggle", CYCLE_SHAPE, "panel selection"]
        .into_iter()
        .chain(SHAPES.iter().map(|s| s.command))
        .collect()
}

/// Every command the pick-level group can run; the first is what its key
/// tip runs.
pub fn level_commands() -> Vec<&'static str> {
    LEVELS.iter().map(|l| l.command).collect()
}

/// Every command the quick-select group can run; the first is what its
/// key tip runs.
pub fn quick_commands() -> Vec<&'static str> {
    QUICK
        .iter()
        .map(|q| q.replace)
        .chain(QUICK.iter().map(|q| q.add))
        .chain([INVERT, CLEAR])
        .collect()
}

/// Runs `add_cells` vertically centered in the ribbon body, `height` tall,
/// with the cells 4 px apart.
fn in_body<R>(ui: &mut Ui, height: f32, add_cells: impl FnOnce(&mut Ui) -> R) -> R {
    ui.vertical(|ui| {
        ui.spacing_mut().item_spacing = egui::Vec2::splat(crate::theme::space::TIGHT);
        ui.add_space((crate::ribbon::BODY_HEIGHT - height) / 2.0);
        add_cells(ui)
    })
    .inner
}

/// The select-tool button and the expression shortcut.
pub fn tool_row(app: &mut AppUi<'_>, ui: &mut Ui) -> Option<String> {
    in_body(ui, widgets::CELL_TALL, |ui| {
        ui.horizontal(|ui| {
            let tool = tool_button(app, ui);
            let typed = expression_button(ui);
            tool.or(typed)
        })
        .inner
    })
}

/// Opens the Selections panel with its expression field focused, for
/// typing a selection instead of picking one.
fn expression_button(ui: &mut Ui) -> Option<String> {
    let response = widgets::captioned_button(
        ui,
        icon::CURSOR_TEXT,
        "Expression",
        widgets::CELL_TALL,
        false,
        false,
    )
    .on_hover_text("Type a selection expression: opens the Selections panel and focuses its field");
    if !response.clicked() {
        return None;
    }
    ui.ctx()
        .memory_mut(|m| m.request_focus(Id::new(crate::ui::SELECTION_FIELD_ID)));
    Some("panel selection".into())
}

/// The pick-level choice: one cell per level, the current one selected.
pub fn level_row(app: &mut AppUi<'_>, ui: &mut Ui) -> Option<String> {
    let now = app.view.select_tool.level;
    in_body(ui, widgets::CELL_TALL, |ui| {
        ui.horizontal(|ui| {
            let mut line = None;
            for l in &LEVELS {
                let cell = widgets::captioned_button(
                    ui,
                    l.icon,
                    l.label,
                    widgets::CELL_TALL,
                    l.level == now,
                    false,
                );
                if cell.on_hover_text(l.tip).clicked() {
                    line = Some(l.command.to_owned());
                }
            }
            line
        })
        .inner
    })
}

/// The quick-select classes in two rows, then invert and clear stacked.
pub fn quick_row(app: &mut AppUi<'_>, ui: &mut Ui) -> Option<String> {
    let has_structure = app.scene.structures().next().is_some();
    let add = ui.input(|i| i.modifiers.shift);
    let height = 2.0 * widgets::CELL_COMPACT + crate::theme::space::TIGHT;
    ui.add_enabled_ui(has_structure, |ui| {
        in_body(ui, height, |ui| {
            ui.horizontal_top(|ui| {
                let classes = class_grid(ui, add);
                classes.or(invert_and_clear(ui))
            })
            .inner
        })
    })
    .inner
}

fn class_grid(ui: &mut Ui, add: bool) -> Option<String> {
    ui.vertical(|ui| {
        let mut line = None;
        for row in QUICK.chunks(QUICK.len() / 2) {
            ui.horizontal(|ui| {
                for q in row {
                    let tip = format!("Select {}. Shift-click adds to the selection.", q.name);
                    let cell = widgets::captioned_button(
                        ui,
                        q.icon,
                        q.label,
                        widgets::CELL_COMPACT,
                        false,
                        false,
                    );
                    if cell.on_hover_text(tip).clicked() {
                        line = Some(if add { q.add } else { q.replace }.to_owned());
                    }
                }
            });
        }
        line
    })
    .inner
}

fn invert_and_clear(ui: &mut Ui) -> Option<String> {
    let cells = [
        (
            icon::SELECTION_INVERSE,
            "Invert",
            "Invert the selection",
            INVERT,
        ),
        (
            icon::SELECTION_SLASH,
            "Clear",
            "Clear the selection (Esc)",
            CLEAR,
        ),
    ];
    ui.vertical(|ui| {
        let mut line = None;
        for (glyph, label, tip, command) in cells {
            let cell =
                widgets::captioned_button(ui, glyph, label, widgets::CELL_COMPACT, false, false);
            if cell.on_hover_text(tip).clicked() {
                line = Some(command.to_owned());
            }
        }
        line
    })
    .inner
}

/// The grouped tool button: a click switches the tool on or off, a press
/// on its corner or a long press opens the flyout of shapes (right-click
/// too). The flyout outlives the press that opened it.
fn tool_button(app: &AppUi<'_>, ui: &mut Ui) -> Option<String> {
    let tool = app.view.select_tool;
    let item = SHAPES.iter().find(|s| s.shape == tool.shape)?;
    let response = widgets::captioned_button(
        ui,
        item.icon,
        item.label,
        widgets::CELL_TALL,
        tool.draws(),
        true,
    )
    .on_hover_text(tool_tip(item));
    let ctx = ui.ctx().clone();
    let by_press = flyout_requested(&ctx, &response);
    let opened_by_press = release_after_flyout_press(&ctx, &response, by_press);
    let toggle = (response.clicked() && !opened_by_press).then(|| "selectmode shape toggle".into());
    let open = by_press || response.secondary_clicked();
    let picked = show_flyout(&ctx, &response, open, opened_by_press, tool.shape);
    toggle.or(picked)
}

/// True on the release that ends a press which opened the flyout, so that
/// release does not also click the button. The flag lives in egui's temp
/// memory between frames.
fn release_after_flyout_press(
    ctx: &egui::Context,
    response: &egui::Response,
    by_press: bool,
) -> bool {
    let flag = Id::new("home-select-flyout-from-press");
    if by_press {
        ctx.data_mut(|d| d.insert_temp(flag, true));
    }
    if response.is_pointer_button_down_on() {
        return false;
    }
    ctx.data_mut(|d| d.remove_temp::<bool>(flag))
        .unwrap_or(false)
        && response.clicked()
}

fn show_flyout(
    ctx: &egui::Context,
    anchor: &egui::Response,
    open: bool,
    just_opened: bool,
    current: SelectShape,
) -> Option<String> {
    let shown = egui::Popup::from_response(anchor)
        .open_memory(open.then_some(egui::SetOpenCommand::Bool(true)))
        .close_behavior(PopupCloseBehavior::IgnoreClicks)
        .show(|ui| flyout(ui, current))?;
    let picked = shown.inner;
    if picked.is_some() || (!just_opened && shown.response.clicked_elsewhere()) {
        egui::Popup::close_id(ctx, egui::Popup::default_response_id(anchor));
    }
    picked
}

fn tool_tip(item: &ShapeItem) -> String {
    let cycle = shortcut_for(CYCLE_SHAPE).map_or(String::new(), |k| format!(" {k} cycles shapes."));
    format!(
        "Select tool ({}): click to switch it on or off.{cycle} Hold, or click the corner, for all shapes.",
        item.label.to_lowercase()
    )
}

/// Whether the press on the button is one that opens the flyout: on its
/// corner marker, or held long enough.
fn flyout_requested(ctx: &egui::Context, response: &egui::Response) -> bool {
    if !response.is_pointer_button_down_on() {
        return false;
    }
    ctx.request_repaint_after(std::time::Duration::from_millis(50));
    let (held, origin) = ctx.input(|i| {
        (
            i.time - i.pointer.press_start_time().unwrap_or(i.time),
            i.pointer.press_origin(),
        )
    });
    let on_corner = origin.is_some_and(|p| widgets::corner_zone(response.rect).contains(p));
    classify_press(held, on_corner) == ToolPress::OpenFlyout
}

fn flyout(ui: &mut Ui, current: SelectShape) -> Option<String> {
    ui.spacing_mut().item_spacing.y = 0.0;
    let mut picked = None;
    for item in &SHAPES {
        let r = widgets::button_selected(
            ui,
            item.icon,
            item.label,
            Variant::Ghost,
            item.shape == current,
        );
        if r.on_hover_text(item.tip).clicked() {
            picked = Some(item.command.to_owned());
        }
    }
    picked
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tables_cover_every_shape_and_level_once_in_order() {
        let shapes: Vec<_> = SHAPES.iter().map(|s| s.shape).collect();
        assert_eq!(shapes, SelectShape::ALL);
        let levels: Vec<_> = LEVELS.iter().map(|l| l.level).collect();
        assert_eq!(levels, SelectLevel::ALL);
        for s in &SHAPES {
            assert_eq!(s.command, format!("selectmode shape {}", s.shape.word()));
        }
        for l in &LEVELS {
            assert_eq!(l.command, format!("selectmode level {}", l.level.word()));
        }
    }

    #[test]
    fn quick_select_pairs_a_replace_and_an_add_form() {
        for q in &QUICK {
            let expr = q.replace.strip_prefix("select ").unwrap();
            assert_eq!(q.add, format!("select add {expr}"));
        }
    }

    #[test]
    fn icons_and_tooltips_are_unique() {
        let mut glyphs: Vec<&str> = SHAPES.iter().map(|s| s.icon).collect();
        glyphs.extend(LEVELS.iter().map(|l| l.icon));
        glyphs.extend(QUICK.iter().map(|q| q.icon));
        glyphs.extend([
            icon::SELECTION_INVERSE,
            icon::SELECTION_SLASH,
            icon::CURSOR_TEXT,
        ]);
        let mut tips: Vec<&str> = SHAPES.iter().map(|s| s.tip).collect();
        tips.extend(LEVELS.iter().map(|l| l.tip));
        for list in [&mut glyphs, &mut tips] {
            let n = list.len();
            list.sort_unstable();
            list.dedup();
            assert_eq!(list.len(), n, "a repeated icon or tooltip");
        }
    }

    #[test]
    fn every_home_command_is_known_and_captions_are_unique() {
        for c in tool_commands()
            .into_iter()
            .chain(level_commands())
            .chain(quick_commands())
        {
            crate::commands::validate(c).unwrap();
        }
        let mut labels: Vec<&str> = LEVELS.iter().map(|l| l.label).collect();
        labels.extend(QUICK.iter().map(|q| q.label));
        labels.extend(["Invert", "Clear", "Expression"]);
        let n = labels.len();
        labels.sort_unstable();
        labels.dedup();
        assert_eq!(labels.len(), n, "a repeated caption");
        assert_eq!(QUICK.len() % 2, 0, "the class grid is two rows");
    }
}
