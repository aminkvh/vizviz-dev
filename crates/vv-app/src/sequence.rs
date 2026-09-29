//! The sequence strip: one row of one-letter residue codes per chain, for
//! every loaded structure, docked above the viewport by default. Click a
//! residue to select its atoms, Ctrl+click to add it to the selection,
//! double-click to zoom to it with its 5 Å neighbourhood (as the viewport
//! does), hover for its name and number. Drawn straight with the painter and
//! virtualized in both directions, so a 2.4M-atom structure costs only the
//! glyphs actually on screen.

use egui::{Align2, Rect, Sense, Ui};
use vv_scene::{Scene, StructureId};

use crate::ui::AppUi;

/// Room for "<structure> <chain>" at the left of every row.
const LABEL_WIDTH: f32 = 96.0;

/// One-letter code for a residue name: upper case for amino acids, lower
/// case for nucleotides, `None` for anything else (water, ligands, ions).
pub fn one_letter(name: &str) -> Option<char> {
    Some(match name {
        "ALA" => 'A',
        "ARG" => 'R',
        "ASN" => 'N',
        "ASP" => 'D',
        "CYS" => 'C',
        "GLN" => 'Q',
        "GLU" => 'E',
        "GLY" => 'G',
        "HIS" => 'H',
        "ILE" => 'I',
        "LEU" => 'L',
        "LYS" => 'K',
        "MET" | "MSE" => 'M',
        "PHE" => 'F',
        "PRO" => 'P',
        "SER" => 'S',
        "THR" => 'T',
        "TRP" => 'W',
        "TYR" => 'Y',
        "VAL" => 'V',
        "SEC" => 'U',
        "PYL" => 'O',
        "A" | "DA" => 'a',
        "C" | "DC" => 'c',
        "G" | "DG" => 'g',
        "U" => 'u',
        "DT" => 't',
        "I" | "DI" => 'i',
        _ => return None,
    })
}

struct Row {
    structure: StructureId,
    label: String,
    residues: std::ops::Range<u32>,
}

fn rows_of(scene: &Scene) -> Vec<Row> {
    let mut rows = Vec::new();
    for (id, loaded) in scene.structures() {
        for (name, residues) in chain_rows(&loaded.structure.topology) {
            rows.push(Row {
                structure: id,
                label: format!("{} {name}", loaded.label),
                residues,
            });
        }
    }
    rows
}

/// One `(chain name, residues)` row per chain name. A PDB `TER` splits a
/// chain letter into polymer and water/ligand records: adjacent records
/// with one name join into one row, and a later non-polymer record that
/// repeats a name is left out rather than shown as a second "A".
fn chain_rows(top: &vv_core::Topology) -> Vec<(String, std::ops::Range<u32>)> {
    let mut rows: Vec<(String, std::ops::Range<u32>)> = Vec::new();
    for (ci, chain) in top.chains.iter().enumerate() {
        let name = top.chain_name(ci);
        let same = rows.iter_mut().rev().find(|r| r.0 == name);
        match same {
            Some(row) if row.1.end == chain.residues.start => row.1.end = chain.residues.end,
            Some(_) if !is_polymer(top, chain.residues.start) => {}
            _ => rows.push((name.to_string(), chain.residues.clone())),
        }
    }
    rows
}

fn is_polymer(top: &vv_core::Topology, residue: u32) -> bool {
    use vv_core::residue_class::ResidueClass;
    matches!(
        top.residue_class(residue as usize),
        ResidueClass::Protein | ResidueClass::Nucleic
    )
}

/// Residues of the actively selected structure that contain at least one
/// selected atom. Recomputed every frame: it's one pass over the selected
/// atoms, not over the structure.
fn selected_residues(scene: &Scene) -> Option<(StructureId, Vec<bool>)> {
    let active = scene.active_selection()?;
    let top = &scene.structure(active.structure)?.structure.topology;
    let mut flags = vec![false; top.residue_count()];
    for atom in active.mask.ones() {
        flags[top.residue_index[atom] as usize] = true;
    }
    Some((active.structure, flags))
}

impl AppUi<'_> {
    pub(crate) fn sequence_ui(&mut self, ui: &mut Ui) {
        let rows = rows_of(self.scene);
        if rows.is_empty() {
            ui.weak("No structure loaded.");
            return;
        }
        let font = egui::TextStyle::Monospace.resolve(ui.style());
        let advance = ui.ctx().fonts_mut(|f| f.glyph_width(&font, 'W')) + 2.0;
        let row_h = ui.text_style_height(&egui::TextStyle::Monospace) + 6.0;
        let longest = rows.iter().map(|r| r.residues.len()).max().unwrap_or(0);
        let total = egui::vec2(
            LABEL_WIDTH + longest as f32 * advance,
            rows.len() as f32 * row_h,
        );
        let selected = selected_residues(self.scene);
        let scene: &Scene = self.scene;

        let mut clicked: Option<(StructureId, u32, bool)> = None;
        let mut focused: Option<(StructureId, u32)> = None;
        egui::ScrollArea::both()
            .auto_shrink([false, false])
            .show_viewport(ui, |ui, viewport| {
                let (rect, response) = ui.allocate_exact_size(total, Sense::click());
                let origin = rect.min;
                let text = ui.visuals().text_color();
                let dim = ui.visuals().weak_text_color();
                let highlight = ui.visuals().selection.bg_fill;
                let painter = ui.painter_at(rect);

                // Only the rows and columns inside `viewport` get drawn.
                let first_row = (viewport.min.y / row_h).floor().max(0.0) as usize;
                let last_row = ((viewport.max.y / row_h).ceil() as usize).min(rows.len());
                let first_col =
                    ((viewport.min.x - LABEL_WIDTH) / advance).floor().max(0.0) as usize;
                let last_col = ((viewport.max.x - LABEL_WIDTH) / advance).ceil().max(0.0) as usize;

                for (r, row) in rows.iter().enumerate().take(last_row).skip(first_row) {
                    let top_y = origin.y + r as f32 * row_h;
                    painter.text(
                        egui::pos2(origin.x + 4.0, top_y + row_h * 0.5),
                        Align2::LEFT_CENTER,
                        &row.label,
                        font.clone(),
                        text,
                    );
                    let Some(loaded) = scene.structure(row.structure) else {
                        continue;
                    };
                    let top = &loaded.structure.topology;
                    let flags = selected
                        .as_ref()
                        .filter(|(id, _)| *id == row.structure)
                        .map(|(_, flags)| flags);
                    for c in first_col..last_col.min(row.residues.len()) {
                        let residue = row.residues.start as usize + c;
                        let cell = Rect::from_min_size(
                            egui::pos2(origin.x + LABEL_WIDTH + c as f32 * advance, top_y),
                            egui::vec2(advance, row_h),
                        );
                        if flags.is_some_and(|f| f[residue]) {
                            painter.rect_filled(cell.shrink(1.0), 2.0, highlight);
                        }
                        let (glyph, color) = match one_letter(top.residue_name(residue)) {
                            Some(glyph) => (glyph, text),
                            None => ('x', dim),
                        };
                        painter.text(
                            cell.center(),
                            Align2::CENTER_CENTER,
                            glyph,
                            font.clone(),
                            color,
                        );
                    }
                }

                // Which (row, column) a pointer position lands on, if any.
                let hit = |pos: egui::Pos2| -> Option<(usize, usize)> {
                    let local = pos - origin;
                    if local.x < LABEL_WIDTH || local.y < 0.0 {
                        return None;
                    }
                    let r = (local.y / row_h).floor() as usize;
                    let c = ((local.x - LABEL_WIDTH) / advance).floor() as usize;
                    (r < rows.len() && c < rows[r].residues.len()).then_some((r, c))
                };
                if let Some((r, c)) = response.hover_pos().and_then(hit) {
                    let row = &rows[r];
                    let residue = row.residues.start as usize + c;
                    if let Some(loaded) = scene.structure(row.structure) {
                        let top = &loaded.structure.topology;
                        let rec = &top.residues[residue];
                        let ins = if rec.ins_code == 0 {
                            String::new()
                        } else {
                            (rec.ins_code as char).to_string()
                        };
                        response.clone().on_hover_text(format!(
                            "{} {}{ins}  ({} atoms)",
                            top.residue_name(residue),
                            rec.auth_seq_id,
                            rec.atoms.len()
                        ));
                    }
                }
                if let Some((r, c)) = response.interact_pointer_pos().and_then(hit) {
                    let residue = (rows[r].structure, rows[r].residues.start + c as u32);
                    if response.double_clicked() {
                        focused = Some(residue);
                    } else if response.clicked() {
                        let add = ui.input(|i| i.modifiers.command);
                        clicked = Some((residue.0, residue.1, add));
                    }
                }
            });

        if let Some((id, residue, add)) = clicked {
            let atoms: Option<Vec<u32>> = self.scene.structure(id).map(|s| {
                s.structure.topology.residues[residue as usize]
                    .atoms
                    .clone()
                    .collect()
            });
            if let Some(atoms) = atoms {
                self.select_atoms(id, &atoms, add);
            }
        }
        if let Some((id, residue)) = focused {
            let first = self
                .scene
                .structure(id)
                .map(|s| s.structure.topology.residues[residue as usize].atoms.start);
            if let Some(atom) = first {
                self.zoom_to_residue(id, atom);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn one_letter_codes_cover_the_standard_residues() {
        assert_eq!(one_letter("ALA"), Some('A'));
        assert_eq!(one_letter("TRP"), Some('W'));
        assert_eq!(one_letter("MSE"), Some('M'));
        assert_eq!(one_letter("DA"), Some('a'));
        assert_eq!(one_letter("U"), Some('u'));
        assert_eq!(one_letter("HOH"), None);
        assert_eq!(one_letter("HEM"), None);
    }

    #[test]
    fn a_ter_split_chain_letter_is_one_sequence_row() {
        let atom = |serial: u32, res: &str, chain: char, seq: u32| {
            format!(
                "{:<6}{serial:>5}  CA  {res:<4}{chain}{seq:>4}    {:8.3}{:8.3}{:8.3}  1.00  0.00           C
",
                if res == "HOH" { "HETATM" } else { "ATOM" },
                seq as f32 * 3.8,
                0.0,
                0.0
            )
        };
        let text = atom(1, "ALA", 'A', 1)
            + &atom(2, "GLY", 'A', 2)
            + "TER
" + &atom(3, "HOH", 'A', 3);
        let top = vv_io::pdb::parse(text.as_bytes()).unwrap().topology.clone();
        assert_eq!(top.chain_count(), 2);
        assert_eq!(chain_rows(&top), [("A".to_string(), 0..3)]);
    }
}
