//! Which written residues are polymer, and how they group into entities:
//! the writers' side of `Topology::polymer_hint` (`mmcif_entity` and
//! `pdb_seqres` are the readers' side).
//!
//! A chain record is cut into segments of one kind so a polymer and the
//! ligands or waters that share its chain record come back as separate
//! entities, and every kind survives a read.

use std::collections::HashMap;
use std::ops::Range;

use vv_core::fixedbitset::FixedBitSet;
use vv_core::{InternId, PolymerHint, ResidueClass, Topology};

use crate::write::residue_atoms;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub(crate) enum Kind {
    Protein,
    Nucleic,
    /// Two or more glycan residues in a row.
    Branched,
    NonPolymer,
    Water,
}

impl Kind {
    pub(crate) fn is_polymer(self) -> bool {
        matches!(self, Kind::Protein | Kind::Nucleic)
    }
}

/// The file's hint decides polymer and water; otherwise the residue class.
fn residue_kind(t: &Topology, residue: usize) -> Kind {
    match t.polymer_hint.get(residue).copied().unwrap_or_default() {
        PolymerHint::Protein => Kind::Protein,
        PolymerHint::Nucleic => Kind::Nucleic,
        PolymerHint::Water => Kind::Water,
        PolymerHint::NonPolymer | PolymerHint::Unknown => match t.residue_class(residue) {
            ResidueClass::Protein => Kind::Protein,
            ResidueClass::Nucleic => Kind::Nucleic,
            ResidueClass::Water => Kind::Water,
            ResidueClass::Glycan => Kind::Branched,
            _ => Kind::NonPolymer,
        },
    }
}

/// A run of consecutive written residues of one chain record and kind.
pub(crate) struct Segment {
    pub chain: u32,
    pub residues: Range<u32>,
    pub kind: Kind,
}

pub(crate) struct Entity {
    pub kind: Kind,
    /// The monomer sequence; a water entity holds its first residue's name.
    pub comps: Vec<InternId>,
    /// Chain records that carry it, by segment.
    pub segments: Vec<usize>,
}

pub(crate) struct Layout {
    pub segments: Vec<Segment>,
    pub entities: Vec<Entity>,
    /// Per segment, an index into `entities`.
    pub entity_of: Vec<usize>,
    /// Per residue, its segment, or `usize::MAX` when none of its atoms is written.
    segment_of: Vec<usize>,
}

impl Layout {
    pub(crate) fn new(t: &Topology, mask: Option<&FixedBitSet>) -> Self {
        let (segments, segment_of) = segments(t, mask);
        let (entities, entity_of) = group_entities(t, &segments);
        Self {
            segments,
            entities,
            entity_of,
            segment_of,
        }
    }

    pub(crate) fn segment(&self, residue: usize) -> Option<usize> {
        Some(self.segment_of[residue]).filter(|&s| s != usize::MAX)
    }

    /// One `label_asym_id` per segment: the chain record's own label, with a
    /// numeric suffix where a chain record was cut or a label repeats.
    pub(crate) fn asym_ids(&self, t: &Topology) -> Vec<String> {
        let mut used = std::collections::HashSet::new();
        self.segments
            .iter()
            .map(|s| {
                let name = t.names.get(t.chains[s.chain as usize].label_asym);
                let id = std::iter::once(name.to_string())
                    .chain((2..).map(|n| format!("{name}{n}")))
                    .find(|id| !used.contains(id))
                    .expect("unbounded candidates");
                used.insert(id.clone());
                id
            })
            .collect()
    }
}

fn segments(t: &Topology, mask: Option<&FixedBitSet>) -> (Vec<Segment>, Vec<usize>) {
    let mut out: Vec<Segment> = Vec::new();
    let mut segment_of = vec![usize::MAX; t.residues.len()];
    for (i, res) in t.residues.iter().enumerate() {
        if residue_atoms(res, mask).is_empty() {
            continue;
        }
        let kind = residue_kind(t, i);
        let extends = out.last().is_some_and(|s| {
            let prev = &t.residues[s.residues.end as usize - 1];
            s.chain == res.chain
                && s.kind == kind
                && (kind != Kind::NonPolymer || prev.comp == res.comp)
        });
        if extends {
            out.last_mut().expect("checked").residues.end = i as u32 + 1;
        } else {
            out.push(Segment {
                chain: res.chain,
                residues: i as u32..i as u32 + 1,
                kind,
            });
        }
        segment_of[i] = out.len() - 1;
    }
    for s in &mut out {
        if s.kind == Kind::Branched && s.residues.len() == 1 {
            s.kind = Kind::NonPolymer;
        }
    }
    (out, segment_of)
}

/// The monomer sequence of the written residues of `segment`. Of a
/// polymer's microheterogeneous residues (one number, several names) only
/// the first counts, as one position of the chain.
pub(crate) fn sequence(t: &Topology, segment: &Segment) -> Vec<InternId> {
    let r = segment.residues.clone();
    let residues = &t.residues[r.start as usize..r.end as usize];
    let mut out = Vec::with_capacity(residues.len());
    let mut prev: Option<(i32, u8)> = None;
    for res in residues {
        let key = (res.seq_id, res.ins_code);
        if segment.kind.is_polymer() && prev == Some(key) {
            continue;
        }
        prev = Some(key);
        out.push(res.comp);
    }
    out
}

/// Segments with the same kind and monomer sequence are one entity (all
/// waters are one, a non-polymer is by its residue name).
fn group_entities(t: &Topology, segments: &[Segment]) -> (Vec<Entity>, Vec<usize>) {
    let mut index: HashMap<(Kind, Vec<InternId>), usize> = HashMap::new();
    let mut entities: Vec<Entity> = Vec::new();
    let mut entity_of = Vec::with_capacity(segments.len());
    for (i, segment) in segments.iter().enumerate() {
        let mut comps = sequence(t, segment);
        if matches!(segment.kind, Kind::Water | Kind::NonPolymer) {
            comps.truncate(1);
        }
        let key = if segment.kind == Kind::Water {
            Vec::new()
        } else {
            comps.clone()
        };
        let e = *index.entry((segment.kind, key)).or_insert_with(|| {
            entities.push(Entity {
                kind: segment.kind,
                comps,
                segments: Vec::new(),
            });
            entities.len() - 1
        });
        entities[e].segments.push(i);
        entity_of.push(e);
    }
    (entities, entity_of)
}
