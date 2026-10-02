//! The Home tab's selection controls, each a group body of captioned icon
//! cells: the select tool with its pick level, the by-type quick-select
//! grid, the modify grid and the view buttons (the interface group is in
//! `home_interface`). Each returns the command line to run (so a click here
//! is also a typeable command) and fills the ribbon's body height itself.

use egui::{Id, PopupCloseBehavior, Ui};
use egui_phosphor::regular as icon;

use crate::keys::shortcut_for;
use crate::select_tool::{
    classify_press, combine_for, Combine, SelectLevel, SelectShape, ToolPress,
};
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
    tip: &'static str,
    command: &'static str,
}

const LEVELS: [LevelItem; 4] = [
    LevelItem {
        level: SelectLevel::Atom,
        label: "Atom",
        tip: "Pick atoms",
        command: "selectmode level atom",
    },
    LevelItem {
        level: SelectLevel::Residue,
        label: "Residue",
        tip: "Pick whole residues",
        command: "selectmode level residue",
    },
    LevelItem {
        level: SelectLevel::Chain,
        label: "Chain",
        tip: "Pick whole chains",
        command: "selectmode level chain",
    },
    LevelItem {
        level: SelectLevel::Molecule,
        label: "Molecule",
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

/// The classes behind the quick-select grid's "More…" cell.
const MORE: [QuickItem; 4] = [
    QuickItem {
        icon: icon::ROWS,
        label: "Membrane",
        name: "membrane lipids",
        replace: "select membrane",
        add: "select add membrane",
    },
    QuickItem {
        icon: icon::PUZZLE_PIECE,
        label: "Cofactor",
        name: "cofactors",
        replace: "select cofactor",
        add: "select add cofactor",
    },
    QuickItem {
        icon: icon::WAVES,
        label: "Solvent",
        name: "water and ions",
        replace: "select solvent",
        add: "select add solvent",
    },
    QuickItem {
        icon: icon::TREE_STRUCTURE,
        label: "Polymer",
        name: "protein and nucleic acid",
        replace: "select polymer",
        add: "select add polymer",
    },
];

/// Keys `selectmode flyout` forces open through `RibbonState::pending_popover`.
pub const SHAPES_FLYOUT: &str = "home.shapes";
pub const MORE_FLYOUT: &str = "home.more";

const CLEAR: &str = "clear";
const CYCLE_SHAPE: &str = "selectmode shape next";

/// Every command the select group can run; the first is what its key tip
/// runs.
pub fn select_commands() -> Vec<&'static str> {
    ["selectmode shape toggle", CYCLE_SHAPE]
        .into_iter()
        .chain(SHAPES.iter().map(|s| s.command))
        .chain(LEVELS.iter().map(|l| l.command))
        .collect()
}

/// Every command the by-type group can run; the first is what its key tip
/// runs.
pub fn quick_commands() -> Vec<&'static str> {
    let classes = || QUICK.iter().chain(&MORE);
    classes()
        .map(|q| q.replace)
        .chain(classes().map(|q| q.add))
        .collect()
}

/// Runs `add_cells` vertically centered in the ribbon body, `height` tall,
/// with the cells 4 px apart.
pub(crate) fn in_body<R>(ui: &mut Ui, height: f32, add_cells: impl FnOnce(&mut Ui) -> R) -> R {
    ui.vertical(|ui| {
        ui.spacing_mut().item_spacing = egui::Vec2::splat(crate::theme::space::TIGHT);
        ui.add_space((crate::ribbon::BODY_HEIGHT - height) / 2.0);
        add_cells(ui)
    })
    .inner
}

/// The height of two compact cells and the gap between them.
pub(crate) fn two_rows() -> f32 {
    2.0 * widgets::CELL_COMPACT + crate::theme::space::TIGHT
}

pub(crate) struct CommandCell {
    pub icon: &'static str,
    pub label: &'static str,
    pub tip: &'static str,
    pub command: &'static str,
}

const VIEW: [CommandCell; 2] = [
    CommandCell {
        icon: icon::FRAME_CORNERS,
        label: "Reset",
        tip: "Reset the view to frame everything",
        command: "view reset",
    },
    CommandCell {
        icon: icon::CAMERA,
        label: "Screenshot",
        tip: "Save an image of the viewport",
        command: "screenshot",
    },
];

const MODIFY: [CommandCell; 8] = [
    CommandCell {
        icon: icon::SELECTION_INVERSE,
        label: "Invert",
        tip: "Select what is not selected",
        command: "select invert",
    },
    CommandCell {
        icon: icon::HEXAGON,
        label: "+Residue",
        tip: "Grow the selection to whole residues",
        command: "select expand residue",
    },
    CommandCell {
        icon: icon::LINK_SIMPLE,
        label: "+Chain",
        tip: "Grow the selection to whole chains",
        command: "select expand chain",
    },
    CommandCell {
        icon: icon::FLASK,
        label: "+Molecule",
        tip: "Grow the selection to whole molecules (connected atoms)",
        command: "select expand molecule",
    },
    CommandCell {
        icon: icon::ARROWS_OUT_LINE_HORIZONTAL,
        label: "Grow",
        tip: "Add the residue on each side of every selected run",
        command: "select grow",
    },
    CommandCell {
        icon: icon::ARROWS_IN_LINE_HORIZONTAL,
        label: "Shrink",
        tip: "Drop the end residue of every selected run",
        command: "select shrink",
    },
    CommandCell {
        icon: icon::ARROWS_OUT_CARDINAL,
        label: "+Shell",
        tip: "Add everything within 5 \u{c5} of the selection",
        command: "select expand within 5",
    },
    CommandCell {
        icon: icon::ARROWS_IN_CARDINAL,
        label: "\u{2212}Shell",
        tip: "Peel off atoms within 5 \u{c5} of what is not selected",
        command: "select shrink within 5",
    },
];

/// Every command the modify group can run; the first is what its key tip
/// runs.
pub fn modify_commands() -> Vec<&'static str> {
    MODIFY.iter().map(|c| c.command).chain([CLEAR]).collect()
}

/// Every command the view group can run; the first is what its key tip
/// runs.
pub fn view_commands() -> Vec<&'static str> {
    VIEW.iter().map(|c| c.command).collect()
}

pub fn view_row(_app: &mut AppUi<'_>, ui: &mut Ui) -> Option<String> {
    in_body(ui, widgets::CELL_TALL, |ui| {
        ui.horizontal(|ui| {
            let mut line = None;
            for c in &VIEW {
                let cell = widgets::captioned_button(
                    ui,
                    c.icon,
                    c.label,
                    widgets::CELL_TALL,
                    false,
                    false,
                );
                if cell.on_hover_text(c.tip).clicked() {
                    line = Some(c.command.to_owned());
                }
            }
            line
        })
        .inner
    })
}

/// The select tool button with the pick level beside it: both set what a
/// click or drag in the viewport selects.
pub fn select_row(app: &mut AppUi<'_>, ui: &mut Ui) -> Option<String> {
    in_body(ui, widgets::CELL_TALL, |ui| {
        ui.horizontal(|ui| {
            let tool = tool_button(app, ui);
            let level = pick_level(app, ui);
            tool.or(level)
        })
        .inner
    })
}

/// "Pick": one joined button per level, the current one selected.
fn pick_level(app: &AppUi<'_>, ui: &mut Ui) -> Option<String> {
    let now = app.view.select_tool.level;
    ui.vertical(|ui| {
        ui.spacing_mut().item_spacing.y = crate::theme::space::TIGHT;
        let used = crate::theme::CONTROL_HEIGHT + crate::theme::text::CAPTION + 6.0;
        ui.add_space((widgets::CELL_TALL - used) / 2.0);
        widgets::caption(ui, "Pick: what a click or drag selects");
        ui.horizontal(|ui| {
            ui.spacing_mut().item_spacing.x = 0.0;
            let mut line = None;
            for l in &LEVELS {
                let button =
                    widgets::button_selected(ui, "", l.label, Variant::Ghost, l.level == now);
                if button.on_hover_text(l.tip).clicked() {
                    line = Some(l.command.to_owned());
                }
            }
            line
        })
        .inner
    })
    .inner
}

/// By type: the quick-select classes in two rows beside the "More…" cell.
pub fn quick_row(app: &mut AppUi<'_>, ui: &mut Ui) -> Option<String> {
    let has_structure = app.scene.structures().next().is_some();
    let how = ui.input(|i| combine_for(i.modifiers, true));
    ui.add_enabled_ui(has_structure, |ui| {
        in_body(ui, two_rows(), |ui| {
            ui.horizontal_top(|ui| {
                let classes = class_grid(ui, how);
                classes.or(more_cell(app, ui, how))
            })
            .inner
        })
    })
    .inner
}

/// A grid of compact cells, `per_row` to a row; `cell` draws cell `i` and
/// returns the command a click chose.
pub(crate) fn compact_grid(
    ui: &mut Ui,
    count: usize,
    per_row: usize,
    mut cell: impl FnMut(&mut Ui, usize) -> Option<String>,
) -> Option<String> {
    ui.vertical(|ui| {
        let mut line = None;
        for row in 0..count.div_ceil(per_row) {
            ui.horizontal(|ui| {
                for i in row * per_row..((row + 1) * per_row).min(count) {
                    line = cell(ui, i).or(line.take());
                }
            });
        }
        line
    })
    .inner
}

/// The command a quick-select cell runs: replace, Shift adds, Ctrl removes.
fn class_line(item: &QuickItem, how: Combine) -> String {
    match how {
        Combine::Replace => item.replace.to_owned(),
        Combine::Add => item.add.to_owned(),
        Combine::Subtract => {
            let expr = item.replace.strip_prefix("select ").unwrap_or(item.replace);
            format!("select remove {expr}")
        }
    }
}

fn class_tip(item: &QuickItem) -> String {
    format!(
        "Select {}. Shift adds to the selection, Ctrl takes them out.",
        item.name
    )
}

/// Takes the deep link for the flyout `key`, if one is pending.
pub(crate) fn take_deep_link(app: &mut AppUi<'_>, key: &'static str) -> bool {
    let asked = app.ribbon.pending_popover == Some(key);
    if asked {
        app.ribbon.pending_popover = None;
    }
    asked
}

/// The "More…" cell, vertically centered beside the grid: a flyout of the
/// classes that have no cell of their own.
fn more_cell(app: &mut AppUi<'_>, ui: &mut Ui, how: Combine) -> Option<String> {
    ui.add_space(crate::theme::space::TIGHT);
    let cell = ui
        .vertical(|ui| {
            ui.add_space((two_rows() - widgets::CELL_TALL) / 2.0);
            widgets::captioned_button(
                ui,
                icon::DOTS_THREE,
                "More…",
                widgets::CELL_TALL,
                false,
                true,
            )
            .on_hover_text("More classes: membrane, cofactor, solvent, polymer")
        })
        .inner;
    let popup = Id::new(MORE_FLYOUT);
    if take_deep_link(app, MORE_FLYOUT) {
        egui::Popup::open_id(ui.ctx(), popup);
    }
    egui::Popup::from_toggle_button_response(&cell)
        .id(popup)
        .close_behavior(PopupCloseBehavior::CloseOnClick)
        .show(|ui| more_flyout(ui, how))?
        .inner
}

fn more_flyout(ui: &mut Ui, how: Combine) -> Option<String> {
    ui.spacing_mut().item_spacing.y = 0.0;
    let mut picked = None;
    for item in &MORE {
        let r = widgets::button_selected(ui, item.icon, item.label, Variant::Ghost, false);
        if r.on_hover_text(class_tip(item)).clicked() {
            picked = Some(class_line(item, how));
        }
    }
    picked
}

fn class_grid(ui: &mut Ui, how: Combine) -> Option<String> {
    compact_grid(ui, QUICK.len(), QUICK.len() / 2, |ui, i| {
        let q = &QUICK[i];
        let cell =
            widgets::captioned_button(ui, q.icon, q.label, widgets::CELL_COMPACT, false, false);
        cell.on_hover_text(class_tip(q))
            .clicked()
            .then(|| class_line(q, how))
    })
}

/// Modify: invert, grow to whole residues, chains or molecules, grow or
/// shrink by a residue or a shell, and clear.
pub fn modify_row(app: &mut AppUi<'_>, ui: &mut Ui) -> Option<String> {
    let has_selection = app.scene.active_selection().is_some();
    let has_structure = app.scene.structures().next().is_some();
    ui.add_enabled_ui(has_structure, |ui| {
        in_body(ui, two_rows(), |ui| {
            ui.horizontal_top(|ui| {
                let grid = compact_grid(ui, MODIFY.len(), MODIFY.len() / 2, |ui, i| {
                    let c = &MODIFY[i];
                    let needs_selection = c.command != "select invert";
                    ui.add_enabled_ui(has_selection || !needs_selection, |ui| {
                        let cell = widgets::captioned_button(
                            ui,
                            c.icon,
                            c.label,
                            widgets::CELL_COMPACT,
                            false,
                            false,
                        );
                        cell.on_hover_text(c.tip)
                            .clicked()
                            .then(|| c.command.to_owned())
                    })
                    .inner
                });
                grid.or(clear_cell(ui, has_selection))
            })
            .inner
        })
    })
    .inner
}

fn clear_cell(ui: &mut Ui, has_selection: bool) -> Option<String> {
    ui.add_space(crate::theme::space::TIGHT);
    ui.vertical(|ui| {
        ui.add_space((two_rows() - widgets::CELL_TALL) / 2.0);
        ui.add_enabled_ui(has_selection, |ui| {
            let cell = widgets::captioned_button(
                ui,
                icon::SELECTION_SLASH,
                "Clear",
                widgets::CELL_TALL,
                false,
                false,
            );
            cell.on_hover_text("Clear the selection (Esc)")
                .clicked()
                .then(|| CLEAR.to_owned())
        })
        .inner
    })
    .inner
}

/// The grouped tool button: a click switches the tool on or off, a press
/// on its corner or a long press opens the flyout of shapes (right-click
/// too). The flyout outlives the press that opened it.
fn tool_button(app: &mut AppUi<'_>, ui: &mut Ui) -> Option<String> {
    let deep_link = take_deep_link(app, SHAPES_FLYOUT);
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
    let open = by_press || response.secondary_clicked() || deep_link;
    let picked = show_flyout(
        &ctx,
        &response,
        open,
        opened_by_press || deep_link,
        tool.shape,
    );
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
        for q in QUICK.iter().chain(&MORE) {
            let expr = q.replace.strip_prefix("select ").unwrap();
            assert_eq!(q.add, format!("select add {expr}"));
        }
    }

    #[test]
    fn shift_adds_and_ctrl_removes_a_class() {
        let protein = &QUICK[0];
        assert_eq!(class_line(protein, Combine::Replace), "select protein");
        assert_eq!(class_line(protein, Combine::Add), "select add protein");
        assert_eq!(
            class_line(protein, Combine::Subtract),
            "select remove protein"
        );
    }

    #[test]
    fn icons_and_tooltips_are_unique() {
        let mut glyphs: Vec<&str> = SHAPES.iter().map(|s| s.icon).collect();
        glyphs.extend(QUICK.iter().chain(&MORE).map(|q| q.icon));
        glyphs.extend(VIEW.iter().chain(&MODIFY).map(|c| c.icon));
        glyphs.extend([
            icon::SELECTION_SLASH,
            icon::DOTS_THREE,
            crate::home_interface::CHAINS_ICON,
            crate::home_interface::POCKET_ICON,
        ]);
        let mut tips: Vec<&str> = SHAPES.iter().map(|s| s.tip).collect();
        tips.extend(LEVELS.iter().map(|l| l.tip));
        tips.extend(VIEW.iter().chain(&MODIFY).map(|c| c.tip));
        for list in [&mut glyphs, &mut tips] {
            let n = list.len();
            list.sort_unstable();
            list.dedup();
            assert_eq!(list.len(), n, "a repeated icon or tooltip");
        }
    }

    #[test]
    fn every_home_command_is_known_and_captions_are_unique() {
        for c in select_commands()
            .into_iter()
            .chain(quick_commands())
            .chain(modify_commands())
            .chain(crate::home_interface::commands())
            .chain(view_commands())
        {
            crate::commands::validate(c).unwrap();
        }
        let mut labels: Vec<&str> = LEVELS.iter().map(|l| l.label).collect();
        labels.extend(QUICK.iter().chain(&MORE).map(|q| q.label));
        labels.extend(VIEW.iter().chain(&MODIFY).map(|c| c.label));
        labels.extend(["Clear", "More…", "Chains…", "Pocket"]);
        let n = labels.len();
        labels.sort_unstable();
        labels.dedup();
        assert_eq!(labels.len(), n, "a repeated caption");
        assert_eq!(QUICK.len() % 2, 0, "the class grid is two rows");
        assert_eq!(MODIFY.len() % 2, 0, "the modify grid is two rows");
    }
}
