//! Selections: which atoms of which structure are "current" or named.

use std::sync::Arc;

use fixedbitset::FixedBitSet;

use crate::StructureId;

/// An immutable, cheaply-cloned atom mask. Cloning a selection (into undo
/// history, into a named set) clones the `Arc`, not the bits.
pub type Mask = Arc<FixedBitSet>;

pub fn mask_of_one(atom_count: usize, atom: u32) -> Mask {
    let mut bits = FixedBitSet::with_capacity(atom_count);
    bits.insert(atom as usize);
    Arc::new(bits)
}

pub fn empty_mask(atom_count: usize) -> Mask {
    Arc::new(FixedBitSet::with_capacity(atom_count))
}

/// A saved, named selection. Distinct from the transient "active" selection
/// (see `Scene::active`), which drives the inspector as the user clicks
/// around and is not itself named until explicitly saved.
#[derive(Clone, Debug, PartialEq)]
pub struct SelectionSet {
    pub name: String,
    pub structure: StructureId,
    pub mask: Mask,
    /// The expression the set was made from, if it came from one (see
    /// `ActiveSelection::expr`). A parsed expression is what a session
    /// file will want to store, not the bitset.
    pub expr: Option<String>,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn one_atom_mask_has_exactly_one_bit() {
        let m = mask_of_one(100, 42);
        assert_eq!(m.count_ones(..), 1);
        assert!(m[42]);
    }

    #[test]
    fn empty_mask_has_no_bits() {
        assert_eq!(empty_mask(50).count_ones(..), 0);
    }
}
