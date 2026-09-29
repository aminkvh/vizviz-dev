//! Undo/redo. `Undo` and `Redo` are deliberately not `Command` variants —
//! they are operations the dispatcher performs, not things that themselves
//! get recorded (undoing would otherwise enter its own history).

use std::collections::VecDeque;

use crate::command::{Command, SceneError};
use crate::scene::Scene;

struct Entry {
    command: Command,
    inverse: Command,
}

/// A capped double-ended history: `dispatch` grows `done`, dropping the
/// oldest entry past `cap` (a `LoadStructure` entry keeps its whole
/// `Structure` alive via its `RestoreStructure` inverse once undone past,
/// so an uncapped history on a long session would leak memory).
pub struct CommandHistory {
    done: VecDeque<Entry>,
    undone: Vec<Entry>,
    cap: usize,
    /// Changes made so far: every dispatch, undo and redo.
    edits: u64,
}

impl CommandHistory {
    pub fn new(cap: usize) -> Self {
        Self {
            done: VecDeque::new(),
            undone: Vec::new(),
            cap: cap.max(1),
            edits: 0,
        }
    }

    /// Applies `command`, and on success records it (dropping the redo
    /// stack, same as any editor). On error the scene and history are
    /// untouched — a failed command never enters history.
    pub fn dispatch(&mut self, scene: &mut Scene, command: Command) -> Result<(), SceneError> {
        let inverse = command.clone().apply(scene)?;
        self.edits += 1;
        self.undone.clear();
        self.done.push_back(Entry { command, inverse });
        if self.done.len() > self.cap {
            self.done.pop_front();
        }
        Ok(())
    }

    /// `Ok(false)` when there is nothing to undo (not an error).
    pub fn undo(&mut self, scene: &mut Scene) -> Result<bool, SceneError> {
        let Some(entry) = self.done.pop_back() else {
            return Ok(false);
        };
        // Redo replays what the inverse handed back, not the original
        // command: for a load that is `RestoreStructure` (same id, no
        // re-parse) rather than loading the file again under a fresh id,
        // which would strand any selection set referring to the old id.
        let command = entry.inverse.clone().apply(scene)?;
        self.edits += 1;
        self.undone.push(Entry {
            command,
            inverse: entry.inverse,
        });
        Ok(true)
    }

    /// `Ok(false)` when there is nothing to redo (not an error).
    pub fn redo(&mut self, scene: &mut Scene) -> Result<bool, SceneError> {
        let Some(entry) = self.undone.pop() else {
            return Ok(false);
        };
        let inverse = entry.command.clone().apply(scene)?;
        self.edits += 1;
        self.done.push_back(Entry {
            command: entry.command,
            inverse,
        });
        Ok(true)
    }

    pub fn can_undo(&self) -> bool {
        !self.done.is_empty()
    }

    pub fn can_redo(&self) -> bool {
        !self.undone.is_empty()
    }

    /// How many changes have been made; the scene changed since a moment
    /// when this differs from its value then (an undo back to it counts
    /// as a change too, as in most editors).
    pub fn edits(&self) -> u64 {
        self.edits
    }
}

impl Default for CommandHistory {
    fn default() -> Self {
        Self::new(100)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::selection::mask_of_one;
    use crate::Command;

    #[test]
    fn cap_drops_oldest_entries() {
        let mut scene = Scene::new();
        let mut history = CommandHistory::new(2);
        let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../fixtures/small/1CRN.cif");
        let mut ids = Vec::new();
        for _ in 0..3 {
            history
                .dispatch(&mut scene, Command::LoadStructure { path: path.clone() })
                .unwrap();
            ids.push(scene.structures().last().unwrap().0);
        }
        // Only the last 2 loads are undoable; the first one's slot can no
        // longer be closed via history (it was pushed out of the cap).
        assert!(history.undo(&mut scene).unwrap());
        assert!(history.undo(&mut scene).unwrap());
        assert!(!history.undo(&mut scene).unwrap());
        assert_eq!(scene.structures().count(), 1);
        assert_eq!(scene.structures().next().unwrap().0, ids[0]);
    }

    #[test]
    fn redo_of_an_undone_load_restores_the_same_id_without_reparsing() {
        let mut scene = Scene::new();
        let mut history = CommandHistory::new(10);
        let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../fixtures/small/1CRN.cif");
        history
            .dispatch(&mut scene, Command::LoadStructure { path })
            .unwrap();
        let id = scene.structures().next().unwrap().0;
        assert!(history.undo(&mut scene).unwrap());
        assert_eq!(scene.structures().count(), 0);
        assert!(history.redo(&mut scene).unwrap());
        let ids: Vec<_> = scene.structures().map(|(id, _)| id).collect();
        assert_eq!(ids, vec![id], "redo must bring back the same slot");
        // The cycle stays undoable afterwards.
        assert!(history.undo(&mut scene).unwrap());
        assert_eq!(scene.structures().count(), 0);
        assert!(history.redo(&mut scene).unwrap());
        assert_eq!(scene.structures().count(), 1);
    }

    #[test]
    fn every_dispatch_undo_and_redo_counts_as_an_edit() {
        let mut scene = Scene::new();
        let mut history = CommandHistory::new(10);
        let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../fixtures/small/1CRN.cif");
        assert_eq!(history.edits(), 0);
        history
            .dispatch(&mut scene, Command::LoadStructure { path: path.clone() })
            .unwrap();
        let saved = history.edits();
        assert!(history.undo(&mut scene).unwrap());
        assert!(history.redo(&mut scene).unwrap());
        assert_eq!(history.edits(), saved + 2);
        // A failed command changes nothing, so it is no edit.
        let missing = Command::LoadStructure {
            path: path.with_file_name("missing.cif"),
        };
        assert!(history.dispatch(&mut scene, missing).is_err());
        assert_eq!(history.edits(), saved + 2);
    }

    #[test]
    fn redo_of_select_is_a_no_op_error_free_when_stack_empty() {
        let mut scene = Scene::new();
        let mut history = CommandHistory::new(10);
        assert!(!history.redo(&mut scene).unwrap());
        let n = 10;
        let _ = mask_of_one(n, 0); // exercised in command.rs tests too
    }
}
