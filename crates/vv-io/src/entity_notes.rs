//! What the file said about each written entity (description, source
//! organism, branch descriptors), found through the entity id the file
//! gave its chains, so a read, write, read keeps it.
//!
//! An mmCIF chain carries its `label_entity_id`; a PDB chain has none, so
//! the `CHAIN:` lists of `COMPND` (kept as `entity_poly.pdbx_strand_id`)
//! stand in.

use vv_core::{AnnotationCategory, Topology};

use crate::polymer_layout::{Entity, Layout};

/// Source tables whose rows follow an entity; each names it in `entity_id`.
pub(crate) const ENTITY_TABLES: [&str; 4] = [
    "entity_src_gen",
    "entity_src_nat",
    "pdbx_entity_branch",
    "pdbx_entity_branch_descriptor",
];

pub(crate) struct EntityNotes<'a> {
    t: &'a Topology,
    /// Per written entity, the id it had in the file.
    original: Vec<Option<String>>,
}

impl<'a> EntityNotes<'a> {
    pub(crate) fn new(t: &'a Topology, layout: &Layout) -> Self {
        let original = layout
            .entities
            .iter()
            .map(|e| original_id(t, layout, e))
            .collect();
        Self { t, original }
    }

    /// The `entity` table's `item` for written entity `entity`.
    pub(crate) fn entity_value(&self, entity: usize, item: &str) -> Option<&'a str> {
        let cat = self.t.annotations.category("entity")?;
        let id = self.original[entity].as_deref()?;
        let row = (0..cat.rows.len()).find(|&r| cat.get("id", r) == Some(id))?;
        cat.get(item, row)
    }

    /// The file's `category` and its rows about written entity `entity`.
    pub(crate) fn rows(
        &self,
        category: &str,
        entity: usize,
    ) -> Option<(&'a AnnotationCategory, Vec<usize>)> {
        let cat = self.t.annotations.category(category)?;
        let id = self.original[entity].as_deref()?;
        let rows = (0..cat.rows.len())
            .filter(|&r| cat.get("entity_id", r) == Some(id))
            .collect();
        Some((cat, rows))
    }

    /// Source organism: expressed (`entity_src_gen`) or natural.
    pub(crate) fn organism(&self, entity: usize) -> Option<&'a str> {
        [
            ("entity_src_gen", "pdbx_gene_src_scientific_name"),
            ("entity_src_nat", "pdbx_organism_scientific"),
        ]
        .into_iter()
        .find_map(|(category, item)| {
            let (cat, rows) = self.rows(category, entity)?;
            rows.into_iter().find_map(|r| cat.get(item, r))
        })
    }
}

fn original_id(t: &Topology, layout: &Layout, entity: &Entity) -> Option<String> {
    let chain = &t.chains[layout.segments[*entity.segments.first()?].chain as usize];
    if chain.entity != 0 {
        return Some(chain.entity.to_string());
    }
    if !entity.kind.is_polymer() {
        return None;
    }
    let strands = t.annotations.category("entity_poly")?;
    let auth = t.names.get(chain.auth_asym);
    (0..strands.rows.len()).find_map(|r| {
        let listed = strands.get("pdbx_strand_id", r)?;
        listed
            .split(',')
            .any(|s| s.trim() == auth)
            .then(|| strands.get("entity_id", r).map(str::to_owned))?
    })
}
