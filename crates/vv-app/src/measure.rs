//! The viewport's Ctrl+click measurement chain: a run of atoms picked
//! while Ctrl stays held, growing a distance into an angle into a
//! dihedral (`vv_scene::Measurement`'s own 2-4 atom range). `ui.rs`'s
//! `viewport_ui` drives the state machine and owns the Ctrl/Escape
//! cancel policy; this module only tracks the chain and paints its
//! pending (not yet a `Measurement`) preview.

use egui::{Color32, Painter, Pos2, Stroke};
use vv_scene::{Measurement, StructureId};

/// What `MeasureChain::push` did, for the caller to turn into commands.
pub enum ChainStep {
    /// Fewer than 2 atoms so far: nothing to show.
    Pending,
    /// The clicked atom was already in the chain: nothing changed.
    Unchanged,
    /// Replace `hide` (the chain's previous measurement, if any -- an
    /// undo-able "hide the 2-atom distance, show the 3-atom angle") with
    /// `show`.
    Show {
        hide: Option<Measurement>,
        show: Measurement,
    },
}

/// The atoms picked so far, in order, on one structure. Starts over on a
/// different structure or once full (a dihedral is the last kind
/// `Measurement` supports); `cancel` drops it without touching whatever
/// measurement is already shown.
#[derive(Default)]
pub struct MeasureChain {
    structure: Option<StructureId>,
    atoms: Vec<u32>,
}

impl MeasureChain {
    /// Adds `atom` of `id` to the chain.
    pub fn push(&mut self, id: StructureId, atom: u32) -> ChainStep {
        if self.structure != Some(id) || self.atoms.len() >= 4 {
            self.structure = Some(id);
            self.atoms.clear();
        }
        if self.atoms.contains(&atom) {
            return ChainStep::Unchanged;
        }
        let hide = Measurement::new(self.atoms.clone());
        self.atoms.push(atom);
        match Measurement::new(self.atoms.clone()) {
            Some(show) => ChainStep::Show { hide, show },
            None => ChainStep::Pending,
        }
    }

    /// Drops the pending chain (Ctrl released, or Escape) without
    /// touching any measurement already dispatched.
    pub fn cancel(&mut self) {
        self.structure = None;
        self.atoms.clear();
    }

    pub fn is_pending(&self) -> bool {
        !self.atoms.is_empty()
    }

    /// The first atom picked (the anchor marker) and the chain's own
    /// structure, once there's at least one.
    pub fn anchor(&self) -> Option<(StructureId, u32)> {
        self.structure.zip(self.atoms.first().copied())
    }

    /// The most recently picked atom -- the rubber band's origin.
    pub fn tip(&self) -> Option<(StructureId, u32)> {
        self.structure.zip(self.atoms.last().copied())
    }
}

/// Ring radius and colour for the anchor marker (`draw_pending`).
const ANCHOR_RADIUS: f32 = 6.0;
const CHAIN_COLOR: Color32 = Color32::from_rgb(0x00, 0xE0, 0x40);

/// The pending chain's live preview: a ring on the anchor, and a rubber
/// band from the last picked atom to the pointer with the live distance
/// to `hover` (the hovered atom's position, if any). Screen-space points
/// only -- `viewport_ui` projects world positions first.
pub fn draw_pending(
    painter: &Painter,
    anchor: Pos2,
    tip: Pos2,
    pointer: Pos2,
    hover: Option<Pos2>,
) {
    painter.circle_stroke(anchor, ANCHOR_RADIUS, Stroke::new(2.0, CHAIN_COLOR));
    painter.line_segment([tip, pointer], Stroke::new(1.5, CHAIN_COLOR));
    if let Some(hover) = hover {
        painter.circle_stroke(hover, ANCHOR_RADIUS * 0.7, Stroke::new(1.5, CHAIN_COLOR));
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn id(n: u32) -> StructureId {
        StructureId::from_raw(n)
    }

    #[test]
    fn first_atom_is_pending_with_no_measurement() {
        let mut chain = MeasureChain::default();
        assert!(matches!(chain.push(id(1), 5), ChainStep::Pending));
        assert_eq!(chain.anchor(), Some((id(1), 5)));
        assert_eq!(chain.tip(), Some((id(1), 5)));
    }

    #[test]
    fn second_atom_shows_a_distance_with_nothing_to_hide() {
        let mut chain = MeasureChain::default();
        chain.push(id(1), 5);
        match chain.push(id(1), 6) {
            ChainStep::Show { hide, show } => {
                assert_eq!(hide, None);
                assert_eq!(show, Measurement::new(vec![5, 6]).unwrap());
            }
            _ => panic!("expected Show"),
        }
        assert_eq!(chain.anchor(), Some((id(1), 5)));
        assert_eq!(chain.tip(), Some((id(1), 6)));
    }

    #[test]
    fn third_and_fourth_atoms_replace_distance_then_angle() {
        let mut chain = MeasureChain::default();
        chain.push(id(1), 1);
        chain.push(id(1), 2);
        match chain.push(id(1), 3) {
            ChainStep::Show { hide, show } => {
                assert_eq!(hide, Measurement::new(vec![1, 2]));
                assert_eq!(show, Measurement::new(vec![1, 2, 3]).unwrap());
            }
            _ => panic!("expected Show"),
        }
        match chain.push(id(1), 4) {
            ChainStep::Show { hide, show } => {
                assert_eq!(hide, Measurement::new(vec![1, 2, 3]));
                assert_eq!(show, Measurement::new(vec![1, 2, 3, 4]).unwrap());
            }
            _ => panic!("expected Show"),
        }
    }

    #[test]
    fn a_fifth_click_starts_over_without_hiding_the_dihedral() {
        let mut chain = MeasureChain::default();
        for atom in 1..=4 {
            chain.push(id(1), atom);
        }
        assert!(matches!(chain.push(id(1), 9), ChainStep::Pending));
        assert_eq!(chain.anchor(), Some((id(1), 9)));
    }

    #[test]
    fn repicking_the_same_atom_is_unchanged() {
        let mut chain = MeasureChain::default();
        chain.push(id(1), 1);
        chain.push(id(1), 2);
        assert!(matches!(chain.push(id(1), 2), ChainStep::Unchanged));
    }

    #[test]
    fn a_different_structure_starts_a_fresh_chain() {
        let mut chain = MeasureChain::default();
        chain.push(id(1), 1);
        chain.push(id(1), 2);
        assert!(matches!(chain.push(id(2), 7), ChainStep::Pending));
        assert_eq!(chain.anchor(), Some((id(2), 7)));
    }

    #[test]
    fn cancel_clears_the_chain_but_is_silent_about_measurements() {
        let mut chain = MeasureChain::default();
        chain.push(id(1), 1);
        chain.push(id(1), 2);
        chain.cancel();
        assert!(!chain.is_pending());
        assert_eq!(chain.anchor(), None);
    }
}
