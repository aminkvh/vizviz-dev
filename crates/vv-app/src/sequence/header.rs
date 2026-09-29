//! The Sequence panel's header: the coloring choice, the Tracks menu, and
//! the key of the enabled tracks.

use std::sync::Arc;

use egui::{vec2, Sense, Ui};
use vv_core::antibody::{CdrDefinition, Scheme};
use vv_scene::Scene;

use super::color::{from_packed, SCHEMES};
use super::tracks::{Glyph, TrackData, PROVIDERS};
use super::SequenceState;
use crate::widgets;

pub fn show(ui: &mut Ui, state: &mut SequenceState, scene: &Scene) {
    ui.horizontal_wrapped(|ui| {
        color_choice(ui, state);
        tracks_menu(ui, state);
        if state.track_on("antibody") {
            antibody_choices(ui, state);
        }
        if state.legend {
            key_inline(ui, state, scene);
        } else if !state.tracks.is_empty() {
            key_button(ui, state, scene);
        }
    });
}

fn color_choice(ui: &mut Ui, state: &mut SequenceState) {
    ui.label("Color");
    let labels: Vec<&str> = SCHEMES.iter().map(|s| s.2).collect();
    let current = SCHEMES.iter().position(|s| s.0 == state.color);
    if let Some(i) = widgets::select(ui, "sequence-color", &labels, current, "None") {
        state.color = SCHEMES[i].0;
    }
}

fn antibody_choices(ui: &mut Ui, state: &mut SequenceState) {
    ui.label("Numbering");
    let names: Vec<&str> = Scheme::ALL.iter().map(|s| s.name()).collect();
    let current = Scheme::ALL.iter().position(|s| *s == state.antibody.scheme);
    if let Some(i) = widgets::select(ui, "sequence-ab-scheme", &names, current, "Kabat") {
        state.antibody.scheme = Scheme::ALL[i];
    }
    ui.label("CDRs");
    let names: Vec<&str> = CdrDefinition::ALL.iter().map(|d| d.name()).collect();
    let current = CdrDefinition::ALL
        .iter()
        .position(|d| *d == state.antibody.cdr);
    if let Some(i) = widgets::select(ui, "sequence-ab-cdr", &names, current, "Kabat") {
        state.antibody.cdr = CdrDefinition::ALL[i];
    }
}

fn tracks_menu(ui: &mut Ui, state: &mut SequenceState) {
    let label = match state.tracks.len() {
        0 => "Tracks".to_string(),
        n => format!("Tracks ({n})"),
    };
    let button = ui.button(format!("{label} \u{25BE}"));
    egui::Popup::from_toggle_button_response(&button)
        .close_behavior(egui::PopupCloseBehavior::CloseOnClickOutside)
        .show(|ui| {
            for p in PROVIDERS {
                let mut on = state.track_on(p.id());
                if ui.checkbox(&mut on, p.label()).changed() {
                    state.set_track(p.id(), on);
                }
            }
            ui.separator();
            ui.checkbox(&mut state.legend, "Show key inline");
        });
}

/// The enabled tracks with something to show, and their legends, from the
/// first structure that has marks for each.
fn key_rows(state: &mut SequenceState, scene: &Scene) -> Vec<(&'static str, Arc<TrackData>)> {
    let enabled: Vec<_> = PROVIDERS
        .iter()
        .filter(|p| state.track_on(p.id()))
        .collect();
    let mut rows = Vec::new();
    for provider in enabled {
        let data = scene.structures().find_map(|(id, loaded)| {
            let data = state.cache.track(id, loaded, *provider, state.antibody);
            let all = 0..loaded.structure.topology.residue_count() as u32;
            data.any_in(&all).then_some(data)
        });
        if let Some(data) = data.filter(|d| d.glyph != Glyph::Ticks) {
            rows.push((provider.label(), data));
        }
    }
    rows
}

fn key_row(ui: &mut Ui, label: &str, data: &TrackData) {
    ui.weak(format!("{label}:"));
    for entry in &data.legend {
        swatch(ui, from_packed(entry.color));
        ui.weak(&entry.label);
    }
}

fn key_inline(ui: &mut Ui, state: &mut SequenceState, scene: &Scene) {
    for (label, data) in key_rows(state, scene) {
        ui.separator();
        key_row(ui, label, &data);
    }
}

fn key_button(ui: &mut Ui, state: &mut SequenceState, scene: &Scene) {
    let button = ui.button("Key \u{25BE}");
    egui::Popup::from_toggle_button_response(&button)
        .close_behavior(egui::PopupCloseBehavior::CloseOnClickOutside)
        .show(|ui| {
            ui.set_max_width(460.0);
            for (label, data) in key_rows(state, scene) {
                ui.horizontal_wrapped(|ui| key_row(ui, label, &data));
            }
        });
}

fn swatch(ui: &mut Ui, color: egui::Color32) {
    let (rect, _) = ui.allocate_exact_size(vec2(9.0, 9.0), Sense::hover());
    ui.painter().rect_filled(rect, 2.0, color);
}
