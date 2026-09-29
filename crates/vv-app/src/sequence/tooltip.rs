//! The hover card of a residue: what it is, and what each visible track
//! says about it.

use egui::{Response, RichText, Ui};
use vv_scene::Scene;

use super::cache::Cache;
use super::layout::Block;
use super::rows::describe_residue;

pub fn show(response: &Response, scene: &Scene, cache: &Cache, block: &Block, residue: u32) {
    let Some(loaded) = scene.structure(block.row.structure) else {
        return;
    };
    let top = &loaded.structure.topology;
    let summary = summary_line(top, cache, block, residue);
    response.clone().on_hover_ui(|ui| {
        ui.label(RichText::new(describe_residue(top, residue)).strong());
        ui.weak(summary);
        for line in &block.tracks {
            if let Some(text) = line.data.describe(residue) {
                track_entry(ui, line.provider.label(), &text);
            }
        }
    });
}

fn track_entry(ui: &mut Ui, label: &str, text: &str) {
    ui.separator();
    ui.label(RichText::new(label).small().weak());
    for line in text.lines() {
        ui.label(line);
    }
}

/// Atom count, mean B-factor and, once computed, relative SASA.
fn summary_line(top: &vv_core::Topology, cache: &Cache, block: &Block, residue: u32) -> String {
    let rec = &top.residues[residue as usize];
    let atoms = rec.atoms.len();
    let b: f32 = rec
        .atoms
        .clone()
        .map(|a| top.b_factor[a as usize])
        .sum::<f32>()
        / atoms.max(1) as f32;
    let mut line = format!("{atoms} atoms · B {b:.1} Å²");
    if let Some(rel) = cache.cached_sasa(block.row.structure, residue) {
        line.push_str(&format!(" · {:.0}% exposed", rel * 100.0));
    }
    line
}
