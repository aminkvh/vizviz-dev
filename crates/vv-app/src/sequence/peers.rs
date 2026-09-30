//! The protein chains of every loaded structure, for tracks that compare
//! a chain with the others (conservation).

use std::hash::{Hash, Hasher};

use vv_scene::{Scene, StructureId};

use super::rows::{chain_rows, is_protein, letters_and_breaks};

pub struct PeerChain {
    pub structure: StructureId,
    pub name: String,
    pub letters: Vec<u8>,
}

#[derive(Default)]
pub struct Peers {
    pub chains: Vec<PeerChain>,
}

impl Peers {
    pub fn of(scene: &Scene) -> Self {
        let mut chains = Vec::new();
        for (id, loaded) in scene.structures() {
            let top = &loaded.structure.topology;
            for (name, residues) in chain_rows(top) {
                if residues.is_empty() || !is_protein(top, residues.start) {
                    continue;
                }
                chains.push(PeerChain {
                    structure: id,
                    name,
                    letters: letters_and_breaks(top, residues).0,
                });
            }
        }
        Peers { chains }
    }

    /// The sequences of every chain but `(structure, name)`, in load order.
    pub fn others(&self, structure: StructureId, name: &str) -> Vec<&[u8]> {
        self.chains
            .iter()
            .filter(|c| !(c.structure == structure && c.name == name))
            .map(|c| c.letters.as_slice())
            .collect()
    }
}

/// Changes whenever a structure is opened, closed or replaced: what
/// [`Peers`] would be built from, without building it.
pub fn fingerprint(scene: &Scene) -> u64 {
    let mut hasher = std::collections::hash_map::DefaultHasher::new();
    for (id, loaded) in scene.structures() {
        let top = &loaded.structure.topology;
        (id, top.residue_count(), top.atom_count()).hash(&mut hasher);
    }
    hasher.finish()
}
