//! The top of the Selections panel: what is selected now (counts by chain
//! and class, with the save and clear actions) and the saved sets. The
//! same selection Home builds; every control here runs a command Home or
//! the console can run too.

use egui::{Id, Ui};
use egui_phosphor::regular as icon;
use vv_core::fixedbitset::FixedBitSet;
use vv_core::{Element, ResidueClass};
use vv_scene::{LoadedStructure, StructureId};

use crate::selection_outline::Highlight;
use crate::ui::AppUi;
use crate::widgets::{self, count, Variant};

/// What the panel and the Inspector show of the active selection: atoms
/// per element, residues per chain and per class. Built when the
/// selection changes, not every frame: a 4M-atom selection is a 4M-bit
/// mask.
#[derive(Default)]
pub struct SelectionSummary {
    /// The selection's mask, by address.
    pub(crate) key: Option<usize>,
    pub(crate) by_element: Vec<(&'static str, usize)>,
    pub(crate) residues: usize,
    /// Residues touched, by chain name.
    pub(crate) by_chain: Vec<(String, usize)>,
    /// Residues touched, by class label, in `ResidueClass::ALL` order.
    pub(crate) by_class: Vec<(&'static str, usize)>,
}

impl SelectionSummary {
    pub(crate) fn of(mask: &FixedBitSet, loaded: &LoadedStructure, key: usize) -> Self {
        let top = &loaded.structure.topology;
        let mut by_element = std::collections::BTreeMap::<&'static str, usize>::new();
        let mut by_chain = std::collections::BTreeMap::<&str, usize>::new();
        let mut classes = [0usize; ResidueClass::COUNT];
        let mut residues = 0;
        let mut last = usize::MAX;
        for atom in mask.ones() {
            *by_element
                .entry(Element::symbol(top.element[atom]))
                .or_default() += 1;
            let r = top.residue_index[atom] as usize;
            if r != last {
                last = r;
                residues += 1;
                *by_chain
                    .entry(top.chain_name(top.residues[r].chain as usize))
                    .or_default() += 1;
                classes[top.residue_class(r).index()] += 1;
            }
        }
        Self {
            key: Some(key),
            by_element: by_element.into_iter().collect(),
            residues,
            by_chain: by_chain
                .into_iter()
                .map(|(name, n)| (name.to_owned(), n))
                .collect(),
            by_class: ResidueClass::ALL
                .iter()
                .filter(|c| classes[c.index()] > 0)
                .map(|c| (c.label(), classes[c.index()]))
                .collect(),
        }
    }
}

/// `"A 120 · B 80"`.
fn joined<S: AsRef<str>>(items: &[(S, usize)]) -> String {
    items
        .iter()
        .map(|(name, n)| format!("{} {n}", name.as_ref()))
        .collect::<Vec<_>>()
        .join(" \u{b7} ")
}

fn save_name_key() -> Id {
    Id::new("selection-panel-save-name")
}

impl AppUi<'_> {
    /// The current selection's summary and the saved sets.
    pub(crate) fn selection_overview_ui(&mut self, ui: &mut Ui) {
        self.selection_card_ui(ui);
        self.saved_sets_ui(ui);
    }

    fn selection_card_ui(&mut self, ui: &mut Ui) {
        let active = self.scene.active_selection().map(|a| {
            let label = self.scene.structure(a.structure).map(|l| l.label.clone());
            (a.mask.count_ones(..), label.unwrap_or_default())
        });
        let Some((atoms, label)) = active.filter(|(n, _)| *n > 0) else {
            widgets::caption(ui, "Nothing selected. Click atoms, or use the Home tab.");
            return;
        };
        self.refresh_selection_summary();
        let s = &*self.selection_summary;
        ui.label(format!(
            "{} \u{b7} {} in {label}",
            count(atoms, "atom"),
            count(s.residues, "residue")
        ));
        widgets::caption(ui, &format!("Chains: {}", joined(&s.by_chain)));
        widgets::caption(ui, &format!("Classes: {}", joined(&s.by_class)));
        self.save_row_ui(ui);
    }

    /// A name field with Save, and Clear.
    fn save_row_ui(&mut self, ui: &mut Ui) {
        let key = save_name_key();
        let mut name: String = ui.data(|d| d.get_temp(key)).unwrap_or_default();
        let mut save = false;
        ui.horizontal(|ui| {
            let field = ui.add(
                egui::TextEdit::singleline(&mut name)
                    .hint_text("Set name")
                    .desired_width(110.0),
            );
            let enter = field.lost_focus() && ui.input(|i| i.key_pressed(egui::Key::Enter));
            let ready = !name.trim().is_empty();
            let button = ui.add_enabled_ui(ready, |ui| {
                widgets::button(ui, icon::FLOPPY_DISK, "Save set", Variant::Secondary)
            });
            save = ready && (button.inner.clicked() || enter);
            if widgets::button(ui, icon::SELECTION_SLASH, "Clear", Variant::Ghost).clicked() {
                self.run_command_logged("clear");
            }
        });
        if save {
            self.run_command_logged(&format!("saveset {}", name.trim()));
            name.clear();
        }
        ui.data_mut(|d| d.insert_temp(key, name));
    }

    fn saved_sets_ui(&mut self, ui: &mut Ui) {
        let sets: Vec<(String, StructureId, vv_scene::Mask)> = self
            .scene
            .selection_sets()
            .iter()
            .map(|s| (s.name.clone(), s.structure, s.mask.clone()))
            .collect();
        if sets.is_empty() {
            return;
        }
        widgets::section(ui, "Saved sets");
        let active = self.scene.active_selection().map(|a| a.mask.clone());
        for (name, id, mask) in sets {
            let current = active
                .as_ref()
                .is_some_and(|a| std::sync::Arc::ptr_eq(a, &mask));
            let words = format!("{name}  \u{b7}  {}", count(mask.count_ones(..), "atom"));
            let mut delete = false;
            let row = widgets::list_row(ui, &words, current, |ui| {
                delete = widgets::button(ui, icon::X, "", Variant::Ghost)
                    .on_hover_text("Delete this set")
                    .clicked();
            });
            let hovered = ui
                .input(|i| i.pointer.hover_pos())
                .is_some_and(|p| row.rect.contains(p));
            if hovered {
                *self.row_highlight = Some(Highlight::Atoms(id, mask));
            }
            if delete {
                self.run_command_logged(&format!("deleteset {name}"));
            } else if row.clicked() {
                self.run_command_logged(&format!("useset {name}"));
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use vv_scene::{Command, CommandHistory, Scene};

    #[test]
    fn the_summary_counts_residues_by_chain_and_class() {
        let mut scene = Scene::new();
        let mut history = CommandHistory::new(10);
        let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../fixtures/small/1AKE.pdb");
        history
            .dispatch(&mut scene, Command::LoadStructure { path })
            .unwrap();
        let loaded = scene.structures().next().unwrap().1;
        let mask = loaded
            .select("chain A and protein and resid 1-10", 0)
            .unwrap();
        let s = SelectionSummary::of(&mask, loaded, 7);
        assert_eq!(s.residues, 10);
        assert_eq!(s.by_chain, [("A".to_owned(), 10)]);
        assert_eq!(s.by_class, [("protein", 10)]);
        assert_eq!(
            s.by_element.iter().map(|e| e.1).sum::<usize>(),
            mask.count_ones(..)
        );
        assert_eq!(joined(&s.by_class), "protein 10");
    }
}
