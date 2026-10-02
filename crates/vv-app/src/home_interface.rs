//! The Home tab's Interface group: pick the residues where two parts of a
//! structure touch. "Chains…" opens a small form (two sides, a distance),
//! "Pocket" selects the residues around every ligand in one click. Both
//! run `select interface ...`.

use egui::{Id, PopupCloseBehavior, Ui};
use egui_phosphor::regular as icon;

use crate::home::{in_body, take_deep_link};
use crate::ui::AppUi;
use crate::widgets;

pub const CHAINS_ICON: &str = icon::GIT_MERGE;
pub const POCKET_ICON: &str = icon::CROSSHAIR;

/// Keys `RibbonState::pending_popover` for the form.
pub const FORM_FLYOUT: &str = "home.interface";

const POCKET: &str = "select interface polymer with ligand within 5";
const EXAMPLE: &str = "select interface chain A to chain B within 5";

/// Every command the group can run; the first is what its key tip runs.
pub fn commands() -> Vec<&'static str> {
    vec![POCKET, EXAMPLE]
}

#[derive(Clone)]
struct Form {
    from: String,
    to: String,
    both: bool,
    reach: f32,
}

impl Default for Form {
    fn default() -> Self {
        Self {
            from: String::new(),
            to: String::new(),
            both: false,
            reach: 5.0,
        }
    }
}

/// The command for a filled form, or `None` while "From" is blank.
fn line(form: &Form) -> Option<String> {
    let from = form.from.trim();
    if from.is_empty() {
        return None;
    }
    let to = form.to.trim();
    let reach = form.reach;
    Some(match (to.is_empty(), form.both) {
        (true, _) => format!("select interface {from} within {reach}"),
        (false, false) => format!("select interface {from} to {to} within {reach}"),
        (false, true) => format!("select interface {from} with {to} within {reach}"),
    })
}

pub fn interface_row(app: &mut AppUi<'_>, ui: &mut Ui) -> Option<String> {
    let has_structure = app.scene.structures().next().is_some();
    ui.add_enabled_ui(has_structure, |ui| {
        in_body(ui, widgets::CELL_TALL, |ui| {
            ui.horizontal(|ui| {
                let form = chains_cell(app, ui);
                form.or(pocket_cell(ui))
            })
            .inner
        })
    })
    .inner
}

fn pocket_cell(ui: &mut Ui) -> Option<String> {
    let cell =
        widgets::captioned_button(ui, POCKET_ICON, "Pocket", widgets::CELL_TALL, false, false);
    cell.on_hover_text("Residues within 5 \u{c5} of any ligand, and the ligand itself")
        .clicked()
        .then(|| POCKET.to_owned())
}

fn chains_cell(app: &mut AppUi<'_>, ui: &mut Ui) -> Option<String> {
    let cell =
        widgets::captioned_button(ui, CHAINS_ICON, "Chains…", widgets::CELL_TALL, false, true)
            .on_hover_text("Residues where one part of the structure touches another");
    let popup = Id::new(FORM_FLYOUT);
    if take_deep_link(app, FORM_FLYOUT) {
        egui::Popup::open_id(ui.ctx(), popup);
    }
    let choices = side_choices(app);
    egui::Popup::from_toggle_button_response(&cell)
        .id(popup)
        .close_behavior(PopupCloseBehavior::CloseOnClickOutside)
        .show(|ui| form_ui(ui, &choices))?
        .inner
}

/// "chain A", "chain B", ... of the current structure, then the classes.
fn side_choices(app: &AppUi<'_>) -> Vec<String> {
    let mut choices: Vec<String> = app
        .current()
        .and_then(|id| app.scene.structure(id))
        .map(|loaded| {
            let top = &loaded.structure.topology;
            let mut names: Vec<&str> = (0..top.chains.len()).map(|c| top.chain_name(c)).collect();
            names.dedup();
            names.iter().map(|n| format!("chain {n}")).collect()
        })
        .unwrap_or_default();
    choices.sort();
    choices.dedup();
    choices.extend(["protein", "nucleic", "ligand", "lipid", "glycan"].map(String::from));
    choices
}

fn form_ui(ui: &mut Ui, choices: &[String]) -> Option<String> {
    let key = Id::new(FORM_FLYOUT).with("form");
    let mut form: Form = ui.data(|d| d.get_temp(key)).unwrap_or_default();
    ui.set_max_width(320.0);
    widgets::section(ui, "Interface");
    side_field(ui, "From", "e.g. chain A", &mut form.from, choices);
    side_field(ui, "To", "any other chain", &mut form.to, choices);
    ui.checkbox(&mut form.both, "Include the other side too");
    ui.horizontal(|ui| {
        ui.label("Within");
        ui.add(
            egui::DragValue::new(&mut form.reach)
                .range(1.0..=20.0)
                .suffix(" \u{c5}"),
        );
    });
    let picked = line(&form);
    ui.data_mut(|d| d.insert_temp(key, form));
    let answer = widgets::form_actions(ui, "Select interface")?;
    if answer {
        ui.close();
        return picked;
    }
    ui.close();
    None
}

/// A labelled expression field with a drop-down that fills it.
fn side_field(ui: &mut Ui, label: &str, hint: &str, text: &mut String, choices: &[String]) {
    ui.horizontal(|ui| {
        ui.label(label);
        ui.add(
            egui::TextEdit::singleline(text)
                .hint_text(hint)
                .desired_width(150.0),
        );
        let names: Vec<&str> = choices.iter().map(String::as_str).collect();
        if let Some(i) = widgets::select(ui, ("interface-side", label), &names, None, "Pick…") {
            *text = choices[i].clone();
        }
    });
}

#[cfg(test)]
mod tests {
    use super::*;

    fn form(from: &str, to: &str, both: bool) -> Form {
        Form {
            from: from.into(),
            to: to.into(),
            both,
            reach: 4.0,
        }
    }

    #[test]
    fn the_form_builds_each_interface_command() {
        assert_eq!(line(&form("", "chain B", false)), None);
        assert_eq!(
            line(&form("chain A", "", false)).unwrap(),
            "select interface chain A within 4"
        );
        assert_eq!(
            line(&form("chain A", "chain B", false)).unwrap(),
            "select interface chain A to chain B within 4"
        );
        assert_eq!(
            line(&form(" ligand ", " chain B ", true)).unwrap(),
            "select interface ligand with chain B within 4"
        );
    }

    #[test]
    fn the_defaults_are_the_documented_five_angstroms() {
        assert_eq!(Form::default().reach, 5.0);
        assert!(POCKET.ends_with("within 5") && EXAMPLE.ends_with("within 5"));
    }
}
