//! The mmCIF entity tables as a statement of which residues are polymer.
//!
//! `_entity.type` says whether an entity is a polymer, non-polymer,
//! branched or water; `_entity_poly.type` says what a polymer is made of;
//! `_chem_comp.type` classifies a monomer when the polymer type does not
//! (a polysaccharide, `other`). Every `label_asym_id` belongs to exactly
//! one entity, so a chain record inherits its entity's kind, and a
//! modified residue inside a polymer entity is polymer whatever its name.

use std::collections::HashMap;

use vv_core::{PolymerHint, Topology};

use crate::cif::Category;
use crate::float::parse_i32;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum EntityKind {
    Protein,
    Nucleic,
    /// A polymer of neither (polysaccharide, `other`): decided per monomer.
    OtherPolymer,
    NonPolymer,
    Water,
}

#[derive(Default)]
pub(crate) struct EntityInfo {
    entities: HashMap<u16, EntityKind>,
    comp_hint: HashMap<String, PolymerHint>,
    /// Canonical one-letter sequence of each protein entity.
    sequences: HashMap<u16, String>,
}

/// `_entity_poly.type` -> what the polymer is made of.
fn polymer_kind(poly_type: &str) -> EntityKind {
    if poly_type.starts_with("polypeptide") || poly_type == "cyclic-pseudo-peptide" {
        EntityKind::Protein
    } else if poly_type.starts_with("polydeoxyribonucleotide")
        || poly_type.starts_with("polyribonucleotide")
    {
        EntityKind::Nucleic
    } else {
        EntityKind::OtherPolymer
    }
}

/// `_chem_comp.type` -> polymer family, or `Unknown` for a saccharide,
/// `non-polymer` or `other`.
fn comp_type_hint(comp_type: &str) -> PolymerHint {
    let t = comp_type.to_ascii_lowercase();
    if t.contains("peptide") && !t.contains("non-polymer") {
        PolymerHint::Protein
    } else if t.contains("dna") || t.contains("rna") {
        PolymerHint::Nucleic
    } else {
        PolymerHint::Unknown
    }
}

fn entity_id(cat: &Category<'_>, row: usize, item: &str) -> Option<u16> {
    cat.get(row, item).and_then(parse_i32).map(|v| v as u16)
}

impl EntityInfo {
    pub(crate) fn collect(cats: &HashMap<String, Category<'_>>) -> Self {
        let mut info = Self::default();
        if let Some(entity) = cats.get("entity") {
            for row in 0..entity.rows.len() {
                let kind = match entity.get_str(row, "type") {
                    Some("polymer") => EntityKind::OtherPolymer,
                    Some("water") => EntityKind::Water,
                    Some(_) => EntityKind::NonPolymer,
                    None => continue,
                };
                if let Some(id) = entity_id(entity, row, "id") {
                    info.entities.insert(id, kind);
                }
            }
        }
        if let Some(poly) = cats.get("entity_poly") {
            for row in 0..poly.rows.len() {
                let (Some(id), Some(t)) =
                    (entity_id(poly, row, "entity_id"), poly.get_str(row, "type"))
                else {
                    continue;
                };
                let kind = polymer_kind(t);
                info.entities.insert(id, kind);
                let code = poly.get_str(row, "pdbx_seq_one_letter_code_can");
                if let (EntityKind::Protein, Some(code)) = (kind, code) {
                    let letters = code.chars().filter(|c| c.is_ascii_alphabetic());
                    info.sequences.insert(id, letters.collect());
                }
            }
        }
        if let Some(comp) = cats.get("chem_comp") {
            for row in 0..comp.rows.len() {
                if let (Some(id), Some(t)) = (comp.get_str(row, "id"), comp.get_str(row, "type")) {
                    info.comp_hint.insert(id.to_string(), comp_type_hint(t));
                }
            }
        }
        info
    }

    fn residue_hint(&self, entity: u16, comp: &str) -> PolymerHint {
        match self.entities.get(&entity) {
            Some(EntityKind::Protein) => PolymerHint::Protein,
            Some(EntityKind::Nucleic) => PolymerHint::Nucleic,
            Some(EntityKind::OtherPolymer) => self.comp_hint.get(comp).copied().unwrap_or_default(),
            Some(EntityKind::NonPolymer) => PolymerHint::NonPolymer,
            Some(EntityKind::Water) => PolymerHint::Water,
            None => PolymerHint::Unknown,
        }
    }

    /// The deposited sequence of each chain record's entity, or empty when
    /// the file gives none.
    pub(crate) fn full_sequences(&self, topology: &Topology) -> Vec<String> {
        if self.sequences.is_empty() {
            return Vec::new();
        }
        topology
            .chains
            .iter()
            .map(|chain| self.sequences.get(&chain.entity).cloned().unwrap_or_default())
            .collect()
    }

    /// One hint per residue, or empty when the file has no entity table.
    pub(crate) fn hints(&self, topology: &Topology) -> Vec<PolymerHint> {
        if self.entities.is_empty() {
            return Vec::new();
        }
        topology
            .residues
            .iter()
            .map(|res| {
                let chain = &topology.chains[res.chain as usize];
                self.residue_hint(chain.entity, topology.names.get(res.comp))
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn polymer_types_map_to_a_family() {
        assert_eq!(polymer_kind("polypeptide(L)"), EntityKind::Protein);
        assert_eq!(polymer_kind("polydeoxyribonucleotide"), EntityKind::Nucleic);
        assert_eq!(
            polymer_kind("polydeoxyribonucleotide/polyribonucleotide hybrid"),
            EntityKind::Nucleic
        );
        assert_eq!(polymer_kind("polysaccharide(D)"), EntityKind::OtherPolymer);
    }

    #[test]
    fn component_types_map_to_a_family() {
        assert_eq!(comp_type_hint("L-peptide linking"), PolymerHint::Protein);
        assert_eq!(comp_type_hint("RNA linking"), PolymerHint::Nucleic);
        assert_eq!(
            comp_type_hint("DNA OH 5 prime terminus"),
            PolymerHint::Nucleic
        );
        assert_eq!(comp_type_hint("non-polymer"), PolymerHint::Unknown);
        assert_eq!(comp_type_hint("D-saccharide"), PolymerHint::Unknown);
    }
}
