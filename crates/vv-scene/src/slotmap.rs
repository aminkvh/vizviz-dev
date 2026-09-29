//! A minimal slot map: stable ids that survive removal and undo/redo.
//!
//! Unlike a generational slot map, indices are never reused — closing a
//! structure leaves a `None` tombstone rather than freeing the slot for
//! reuse. That trades a few bytes per closed structure (irrelevant for a
//! desktop session) for a real simplification: an `Id` is valid forever
//! once allocated, so undo/redo never has to worry about a stale id
//! aliasing a newer value at the same slot.

use std::hash::{Hash, Hasher};
use std::marker::PhantomData;

pub struct Id<T> {
    index: u32,
    marker: PhantomData<fn() -> T>,
}

impl<T> Id<T> {
    /// For (de)serializing ids as plain indices (session files, tests).
    pub fn from_raw(index: u32) -> Self {
        Self {
            index,
            marker: PhantomData,
        }
    }

    pub fn to_raw(self) -> u32 {
        self.index
    }
}

// Manual impls: `derive` would require `T: Clone/Copy/...`, but `Id<T>` is
// just a `u32` regardless of `T`.
impl<T> Clone for Id<T> {
    fn clone(&self) -> Self {
        *self
    }
}
impl<T> Copy for Id<T> {}
impl<T> PartialEq for Id<T> {
    fn eq(&self, other: &Self) -> bool {
        self.index == other.index
    }
}
impl<T> Eq for Id<T> {}
impl<T> Hash for Id<T> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.index.hash(state);
    }
}
impl<T> std::fmt::Debug for Id<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Id({})", self.index)
    }
}

#[derive(Clone, Debug)]
pub struct SlotMap<T> {
    slots: Vec<Option<T>>,
}

// Manual, not derived: `#[derive(Default)]` would require `T: Default`,
// which an empty `Vec<Option<T>>` never needs.
impl<T> Default for SlotMap<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T> SlotMap<T> {
    pub fn new() -> Self {
        Self { slots: Vec::new() }
    }

    pub fn insert(&mut self, value: T) -> Id<T> {
        let index = self.slots.len() as u32;
        self.slots.push(Some(value));
        Id::from_raw(index)
    }

    /// Puts `value` back into a specific (previously emptied) slot. Used
    /// only to undo a `remove`; panics on a slot that was never allocated
    /// or is already occupied, since that indicates a logic error, not
    /// user-reachable state.
    pub fn restore(&mut self, id: Id<T>, value: T) {
        let slot = &mut self.slots[id.index as usize];
        assert!(slot.is_none(), "restore into an occupied slot");
        *slot = Some(value);
    }

    pub fn remove(&mut self, id: Id<T>) -> Option<T> {
        self.slots.get_mut(id.index as usize).and_then(Option::take)
    }

    pub fn get(&self, id: Id<T>) -> Option<&T> {
        self.slots.get(id.index as usize).and_then(Option::as_ref)
    }

    pub fn get_mut(&mut self, id: Id<T>) -> Option<&mut T> {
        self.slots
            .get_mut(id.index as usize)
            .and_then(Option::as_mut)
    }

    pub fn contains(&self, id: Id<T>) -> bool {
        self.get(id).is_some()
    }

    pub fn iter(&self) -> impl DoubleEndedIterator<Item = (Id<T>, &T)> {
        self.slots
            .iter()
            .enumerate()
            .filter_map(|(i, s)| s.as_ref().map(|v| (Id::from_raw(i as u32), v)))
    }

    pub fn len(&self) -> usize {
        self.slots.iter().filter(|s| s.is_some()).count()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn insert_get_remove_restore_round_trip() {
        let mut m: SlotMap<&str> = SlotMap::new();
        let a = m.insert("a");
        let b = m.insert("b");
        assert_eq!(m.get(a), Some(&"a"));
        assert_eq!(m.len(), 2);

        let removed = m.remove(a).unwrap();
        assert_eq!(removed, "a");
        assert_eq!(m.get(a), None);
        assert_eq!(m.len(), 1);
        // `b` keeps its id after `a`'s slot is emptied.
        assert_eq!(m.get(b), Some(&"b"));

        m.restore(a, "a again");
        assert_eq!(m.get(a), Some(&"a again"));
        assert_eq!(m.len(), 2);
    }

    #[test]
    fn ids_from_different_maps_do_not_collide_by_construction() {
        let mut m: SlotMap<i32> = SlotMap::new();
        let a = m.insert(1);
        let b = m.insert(2);
        assert_ne!(a, b);
        assert_eq!(a, Id::from_raw(a.to_raw()));
    }

    #[test]
    #[should_panic(expected = "occupied slot")]
    fn restore_into_occupied_slot_panics() {
        let mut m: SlotMap<i32> = SlotMap::new();
        let a = m.insert(1);
        m.restore(a, 2);
    }
}
