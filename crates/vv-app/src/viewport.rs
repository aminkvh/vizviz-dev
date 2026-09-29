//! Pure helpers for viewport picking and the double-click residue zoom,
//! kept free of `AppUi` so they're unit-testable without a live scene.
//! `ui.rs`'s `viewport_ui` is the only caller.

use std::ops::Range;

use vv_core::glam::Vec3;
use vv_core::Topology;
use vv_scene::StructureId;

/// Tracks the pointer's last picked pixel so `viewport_ui` only re-picks
/// (a blocking GPU readback, `Renderer::pick`) when it actually moved,
/// not on every incidental repaint (trajectory playback, a fading toast,
/// ...) while the pointer sits still.
#[derive(Default)]
pub struct HoverPick {
    last_pixel: Option<(u32, u32)>,
    /// The residue-outline target and the measure chain's live rubber
    /// band both read this; `None` when nothing is hovered.
    pub atom: Option<(StructureId, u32)>,
}

impl HoverPick {
    /// `pixel` is `None` when nothing should be hovered right now (the
    /// pointer left the viewport, or a mouse button is held). Returns
    /// whether it differs from last time -- the caller re-picks only
    /// then, and should clear/set `atom` accordingly.
    pub fn moved(&mut self, pixel: Option<(u32, u32)>) -> bool {
        let changed = pixel != self.last_pixel;
        self.last_pixel = pixel;
        changed
    }
}

/// The render-target pixel a pointer position corresponds to, clamped
/// inside a `width x height` canvas -- shared by click, double-click and
/// hover picking so the three agree on the same mapping.
pub fn pointer_pixel(pos: egui::Pos2, rect: egui::Rect, width: u32, height: u32) -> (u32, u32) {
    let local = pos - rect.min;
    let x = (local.x / rect.width().max(1.0) * width as f32) as u32;
    let y = (local.y / rect.height().max(1.0) * height as f32) as u32;
    (
        x.min(width.saturating_sub(1)),
        y.min(height.saturating_sub(1)),
    )
}

/// The atom range of the residue containing `atom`.
pub fn residue_atoms(top: &Topology, atom: u32) -> Range<u32> {
    let r = top.residue_index[atom as usize] as usize;
    top.residues[r].atoms.clone()
}

/// A selection expression identifying exactly the residue containing
/// `atom` -- by chain and author residue number, the only handle
/// SELECTION.md's grammar gives a single residue (an insertion code
/// can't be expressed, so a residue sharing chain+resid with another
/// via one is picked ambiguously; accepted as a rare edge case).
pub fn residue_expr(top: &Topology, atom: u32) -> String {
    let res = &top.residues[top.residue_index[atom as usize] as usize];
    format!(
        "chain {} and resid {}",
        top.chain_name(res.chain as usize),
        res.auth_seq_id
    )
}

/// "Whole residues within `radius` Å" of the residue containing `atom`.
/// SELECTION.md lists `same residue as` under "Not yet"; `byres` is its
/// supported equivalent, and `within` is already inclusive of its own
/// expr, so this also includes the residue itself.
pub fn neighborhood_expr(top: &Topology, atom: u32, radius: f32) -> String {
    format!("byres within {radius} of ({})", residue_expr(top, atom))
}

/// Whether `selection` is one `neighborhood_expr` wrote, so a new
/// double-click can replace the layer the last one added.
pub fn is_neighborhood_expr(selection: &str) -> bool {
    selection.starts_with("byres within ") && selection.contains(" of (chain ")
}

/// Centroid and bounding radius of `atoms`' positions; `(ZERO, 0.0)` for
/// an empty range.
pub fn bounding_sphere(positions: &[Vec3], atoms: Range<u32>) -> (Vec3, f32) {
    let points: Vec<Vec3> = atoms.map(|a| positions[a as usize]).collect();
    if points.is_empty() {
        return (Vec3::ZERO, 0.0);
    }
    let center = points.iter().copied().sum::<Vec3>() / points.len() as f32;
    let radius = points
        .iter()
        .map(|p| p.distance(center))
        .fold(0.0f32, f32::max);
    (center, radius)
}

#[cfg(test)]
mod tests {
    use super::*;
    use vv_core::builder::{AtomRow, TopologyBuilder};
    use vv_core::{Element, Structure};

    #[test]
    fn a_generated_neighborhood_is_recognised_and_user_layers_are_not() {
        let s = fixture();
        assert!(is_neighborhood_expr(&neighborhood_expr(
            &s.topology,
            0,
            5.0
        )));
        assert!(!is_neighborhood_expr("all"));
        assert!(!is_neighborhood_expr("within 5 of resname HEM"));
    }

    /// Two residues, one two-atom ("A" 42) and one one-atom ("B" 7), for
    /// `residue_expr`/`neighborhood_expr`/`bounding_sphere`.
    fn fixture() -> Structure {
        let mut b = TopologyBuilder::new();
        let rows: [(&str, &str, &str, i32, Vec3); 3] = [
            ("CA", "ALA", "A", 42, Vec3::new(0.0, 0.0, 0.0)),
            ("CB", "ALA", "A", 42, Vec3::new(4.0, 0.0, 0.0)),
            ("O", "HOH", "B", 7, Vec3::new(10.0, 0.0, 0.0)),
        ];
        for (i, &(name, comp, asym, resid, position)) in rows.iter().enumerate() {
            let mut n = [b' '; 4];
            n[..name.len()].copy_from_slice(name.as_bytes());
            b.push(&AtomRow {
                element: Element::UNKNOWN,
                name: n,
                serial: i as u32 + 1,
                alt_loc: 0,
                comp,
                asym,
                auth_asym: asym,
                seq_id: resid,
                auth_seq_id: resid,
                ins_code: 0,
                entity: 1,
                position,
                occupancy: 1.0,
                b_factor: 0.0,
                charge: 0,
                hetero: false,
            });
        }
        b.finish().unwrap()
    }

    #[test]
    fn pointer_pixel_scales_and_clamps() {
        let rect = egui::Rect::from_min_size(egui::pos2(10.0, 10.0), egui::vec2(100.0, 50.0));
        assert_eq!(
            pointer_pixel(egui::pos2(60.0, 35.0), rect, 200, 100),
            (100, 50)
        );
        // At the rect's far edge, clamped inside the canvas (not `width`).
        assert_eq!(
            pointer_pixel(egui::pos2(110.0, 60.0), rect, 200, 100),
            (199, 99)
        );
    }

    #[test]
    fn hover_pick_reports_a_change_only_once() {
        let mut hover = HoverPick::default();
        assert!(hover.moved(Some((1, 1))), "first pixel is always a change");
        assert!(!hover.moved(Some((1, 1))), "same pixel is not a change");
        assert!(hover.moved(None), "leaving is a change");
        assert!(!hover.moved(None), "already gone is not a change");
    }

    #[test]
    fn residue_atoms_spans_just_that_residue() {
        let s = fixture();
        assert_eq!(residue_atoms(&s.topology, 0), 0..2);
        assert_eq!(residue_atoms(&s.topology, 1), 0..2);
        assert_eq!(residue_atoms(&s.topology, 2), 2..3);
    }

    #[test]
    fn residue_expr_names_chain_and_auth_resid() {
        let s = fixture();
        assert_eq!(residue_expr(&s.topology, 0), "chain A and resid 42");
        assert_eq!(residue_expr(&s.topology, 2), "chain B and resid 7");
    }

    #[test]
    fn neighborhood_expr_wraps_it_in_byres_within() {
        let s = fixture();
        assert_eq!(
            neighborhood_expr(&s.topology, 0, 5.0),
            "byres within 5 of (chain A and resid 42)"
        );
        // The grammar this builds against: `within` must parse its whole
        // parenthesised operand, and `byres` must accept `within`'s result.
        vv_core::select::parse(&neighborhood_expr(&s.topology, 0, 5.0)).expect("valid selection");
    }

    #[test]
    fn bounding_sphere_of_two_atoms() {
        let s = fixture();
        let coords = s.frame(0);
        let (center, radius) = bounding_sphere(coords.positions(), residue_atoms(&s.topology, 0));
        assert_eq!(center, Vec3::new(2.0, 0.0, 0.0));
        assert_eq!(radius, 2.0);
    }

    #[test]
    fn bounding_sphere_of_no_atoms_is_a_point() {
        assert_eq!(bounding_sphere(&[], 0..0), (Vec3::ZERO, 0.0));
    }
}
