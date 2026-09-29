//! String interning for repeated names (residue types, chain ids).

use std::collections::HashMap;

/// Handle to an interned string. Only meaningful with the `Interner` that made it.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug, Default, PartialOrd, Ord)]
#[repr(transparent)]
pub struct InternId(pub u32);

#[derive(Clone, Debug, Default)]
pub struct Interner {
    strings: Vec<Box<str>>,
    index: HashMap<Box<str>, InternId>,
}

impl Interner {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn intern(&mut self, s: &str) -> InternId {
        if let Some(&id) = self.index.get(s) {
            return id;
        }
        let id = InternId(self.strings.len() as u32);
        let boxed: Box<str> = s.into();
        self.strings.push(boxed.clone());
        self.index.insert(boxed, id);
        id
    }

    pub fn lookup(&self, s: &str) -> Option<InternId> {
        self.index.get(s).copied()
    }

    pub fn get(&self, id: InternId) -> &str {
        &self.strings[id.0 as usize]
    }

    pub fn len(&self) -> usize {
        self.strings.len()
    }

    pub fn is_empty(&self) -> bool {
        self.strings.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn interning_is_idempotent_and_ordered() {
        let mut i = Interner::new();
        let ala = i.intern("ALA");
        let gly = i.intern("GLY");
        assert_eq!(i.intern("ALA"), ala);
        assert_ne!(ala, gly);
        assert_eq!(i.get(ala), "ALA");
        assert_eq!(i.get(gly), "GLY");
        assert_eq!(i.lookup("GLY"), Some(gly));
        assert_eq!(i.lookup("SER"), None);
        assert_eq!(i.len(), 2);
    }
}
