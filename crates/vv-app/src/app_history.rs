//! View settings and camera jumps join the scene's own undo
//! history from the user's perspective, without touching
//! `vv_scene::CommandHistory` (another agent's crate). A minimal
//! app-side wrapper was the cleanest of the options considered: `Scene`
//! is a marker recorded here in lockstep with every successful
//! `AppUi::dispatch` (see its doc), so undo/redo of a marker just
//! replays on `vv_scene::CommandHistory`, which still owns the real
//! scene-edit stack; `View`/`Camera` carry their own before/after
//! snapshot and are applied directly, since neither ever touches
//! `Scene`. Ordinary camera drag/pan/zoom/wheel are deliberately never
//! recorded here (only the three explicit "jump" call sites are) --
//! see `State::redraw`'s `camera_jump` and `AppUi::commit_view_edit`.

use std::collections::VecDeque;

use vv_render::Camera;

use crate::ui::ViewSettings;

pub enum Edit {
    /// A successful `AppUi::dispatch`; undo/redo replay on
    /// `vv_scene::CommandHistory`, which this is kept in lockstep with.
    Scene,
    View(Box<ViewSettings>, Box<ViewSettings>),
    Camera(Box<Camera>, Box<Camera>),
}

/// A capped double-ended history, the same shape as `vv_scene::
/// CommandHistory` (see its own doc for why capped) and the same
/// default cap, so the two stay roughly in step under ordinary use.
pub struct AppHistory {
    done: VecDeque<Edit>,
    undone: Vec<Edit>,
    cap: usize,
}

impl AppHistory {
    pub fn new(cap: usize) -> Self {
        Self {
            done: VecDeque::new(),
            undone: Vec::new(),
            cap: cap.max(1),
        }
    }

    /// Records `edit`, dropping the redo stack like any editor.
    pub fn push(&mut self, edit: Edit) {
        self.undone.clear();
        self.done.push_back(edit);
        if self.done.len() > self.cap {
            self.done.pop_front();
        }
    }

    pub fn pop_undo(&mut self) -> Option<Edit> {
        self.done.pop_back()
    }

    pub fn pop_redo(&mut self) -> Option<Edit> {
        self.undone.pop()
    }

    pub fn push_undone(&mut self, edit: Edit) {
        self.undone.push(edit);
    }

    pub fn push_done(&mut self, edit: Edit) {
        self.done.push_back(edit);
    }

    pub fn can_undo(&self) -> bool {
        !self.done.is_empty()
    }

    pub fn can_redo(&self) -> bool {
        !self.undone.is_empty()
    }
}

impl Default for AppHistory {
    fn default() -> Self {
        Self::new(100)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn is_scene(edit: &Edit) -> bool {
        matches!(edit, Edit::Scene)
    }

    /// Mirrors what `AppUi::app_undo`/`app_redo` actually do: every pop
    /// from one stack is pushed onto the other, keeping them symmetric.
    fn undo(h: &mut AppHistory) -> Option<Edit> {
        let edit = h.pop_undo()?;
        h.push_undone(match &edit {
            Edit::Scene => Edit::Scene,
            Edit::View(b, a) => Edit::View(b.clone(), a.clone()),
            Edit::Camera(b, a) => Edit::Camera(b.clone(), a.clone()),
        });
        Some(edit)
    }

    fn redo(h: &mut AppHistory) -> Option<Edit> {
        let edit = h.pop_redo()?;
        h.push_done(match &edit {
            Edit::Scene => Edit::Scene,
            Edit::View(b, a) => Edit::View(b.clone(), a.clone()),
            Edit::Camera(b, a) => Edit::Camera(b.clone(), a.clone()),
        });
        Some(edit)
    }

    #[test]
    fn undo_then_redo_round_trips_in_lifo_order() {
        let mut h = AppHistory::new(100);
        h.push(Edit::Scene); // 1st
        h.push(Edit::View(Box::default(), Box::default())); // 2nd
        h.push(Edit::Scene); // 3rd

        // Undo pops most-recent-first: 3rd (Scene), 2nd (View), 1st (Scene).
        assert!(is_scene(&undo(&mut h).unwrap()));
        assert!(!is_scene(&undo(&mut h).unwrap()));
        assert!(is_scene(&undo(&mut h).unwrap()));
        assert!(undo(&mut h).is_none(), "nothing left to undo");

        // Redo replays the same three, in the order they were undone.
        assert!(is_scene(&redo(&mut h).unwrap()));
        assert!(!is_scene(&redo(&mut h).unwrap()));
        assert!(is_scene(&redo(&mut h).unwrap()));
        assert!(redo(&mut h).is_none(), "nothing left to redo");
    }

    #[test]
    fn a_fresh_push_drops_the_redo_stack() {
        let mut h = AppHistory::new(100);
        h.push(Edit::Scene);
        let undone = h.pop_undo().unwrap();
        h.push_undone(undone);
        assert!(h.can_redo());
        h.push(Edit::Scene);
        assert!(!h.can_redo(), "a fresh push must drop the redo stack");
    }

    #[test]
    fn cap_drops_the_oldest_entry() {
        let mut h = AppHistory::new(2);
        h.push(Edit::Scene);
        h.push(Edit::Scene);
        h.push(Edit::Scene);
        assert!(h.pop_undo().is_some());
        assert!(h.pop_undo().is_some());
        assert!(
            h.pop_undo().is_none(),
            "the cap should have dropped the first push"
        );
    }

    #[test]
    fn camera_and_view_edits_carry_their_own_before_and_after() {
        let before = Camera::framing(glam::Vec3::ZERO, 1.0);
        let mut after = before.clone();
        after.distance *= 2.0;
        let mut h = AppHistory::new(100);
        h.push(Edit::Camera(
            Box::new(before.clone()),
            Box::new(after.clone()),
        ));
        match h.pop_undo().unwrap() {
            Edit::Camera(b, a) => {
                assert_eq!(*b, before);
                assert_eq!(*a, after);
            }
            _ => panic!("expected a Camera edit"),
        }
    }
}
