//! Hover cards of the strip: a residue (what it is, and what each visible
//! track says about it), a chain label (its bulk properties) and a ligand
//! chip.

use egui::{Response, RichText, Ui};
use vv_scene::Scene;

use super::cache::Cache;
use super::chain_props;
use super::layout::Block;
use super::rows::describe_residue;

/// Ligand residues listed by name in a chip's card.
const LISTED_LIGANDS: usize = 8;

pub fn residue(response: &Response, scene: &Scene, cache: &mut Cache, block: &Block, residue: u32) {
    let Some(loaded) = scene.structure(block.row.structure) else {
        return;
    };
    let top = &loaded.structure.topology;
    let exposed = cache.sasa_at(block.row.structure, loaded, residue);
    let summary = summary_line(top, exposed, residue);
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

pub fn chain(response: &Response, scene: &Scene, block: &Block) {
    let Some(loaded) = scene.structure(block.row.structure) else {
        return;
    };
    let top = &loaded.structure.topology;
    let props = chain_props::of_chain(top, block.row.residues.clone());
    response.clone().on_hover_ui(|ui| {
        ui.label(RichText::new(&block.row.label).strong());
        match &props {
            Some(p) => {
                for line in chain_props::lines(p) {
                    ui.label(line);
                }
            }
            None => {
                ui.weak("Not a protein chain: no protein properties.");
            }
        }
    });
}

pub fn ligands(response: &Response, scene: &Scene, block: &Block, group: usize) {
    let Some(loaded) = scene.structure(block.row.structure) else {
        return;
    };
    let top = &loaded.structure.topology;
    let group = &block.row.ligands[group];
    response.clone().on_hover_ui(|ui| {
        ui.label(RichText::new(group.text()).strong());
        for &r in group.residues.iter().take(LISTED_LIGANDS) {
            ui.label(describe_residue(top, r));
        }
        let more = group.residues.len().saturating_sub(LISTED_LIGANDS);
        if more > 0 {
            ui.weak(format!("and {more} more"));
        }
        ui.weak("Click to select all");
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
fn summary_line(top: &vv_core::Topology, exposed: Option<f32>, residue: u32) -> String {
    let rec = &top.residues[residue as usize];
    let atoms = rec.atoms.len();
    let b: f32 = rec
        .atoms
        .clone()
        .map(|a| top.b_factor[a as usize])
        .sum::<f32>()
        / atoms.max(1) as f32;
    let mut line = format!("{atoms} atoms · B {b:.1} Å²");
    if let Some(rel) = exposed {
        line.push_str(&format!(" · {:.0}% exposed", rel * 100.0));
    }
    line
}
