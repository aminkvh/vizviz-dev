//! The viewport's selection outline: which drawn pick ids belong to the
//! active selection, a panel row's highlight, or what a click on the
//! hovered atom would pick (`Renderer::set_selection`) -- all three share
//! the same outline, unioned per atom.

use std::sync::Arc;

use vv_core::fixedbitset::FixedBitSet;
use vv_render::{ItemSelection, Renderer};
use vv_scene::{ActiveSelection, LoadedStructure, Mask, RepId, Scene, StructureId};

use crate::gpu_cache::DrawSource;
use crate::select_tool::{picked_atoms, SelectLevel};

/// What a hovered panel row asks the viewport to outline.
#[derive(Clone, Debug)]
pub enum Highlight {
    /// One rep's atoms.
    Rep(StructureId, RepId),
    /// Every atom of a structure.
    Structure(StructureId),
    /// A saved set's atoms.
    Atoms(StructureId, Mask),
}

type HighlightKey = (StructureId, usize, usize);

impl Highlight {
    /// Identity for change detection: the mask by address, never by content.
    fn key(&self) -> HighlightKey {
        match self {
            Self::Rep(id, rep) => (*id, 0, rep.0 as usize),
            Self::Structure(id) => (*id, 1, 0),
            Self::Atoms(id, mask) => (*id, 2, Arc::as_ptr(mask) as usize),
        }
    }

    fn covers_source(&self, source: &DrawSource) -> bool {
        match self {
            Self::Rep(id, rep) => (*id, *rep) == (source.id, source.rep),
            Self::Structure(id) => *id == source.id,
            Self::Atoms(..) => false,
        }
    }

    fn atoms_of(&self, id: StructureId) -> Option<&Mask> {
        match self {
            Self::Atoms(set, mask) if *set == id => Some(mask),
            _ => None,
        }
    }
}

/// What the renderer's selection bits were last built from: the active
/// mask by address, the row highlight, the hovered atom and pick level,
/// and each draw source's structure and maps by address.
#[derive(Default, PartialEq)]
pub struct OutlineKey {
    mask: Option<usize>,
    highlight: Option<HighlightKey>,
    hover: Option<(StructureId, u32, SelectLevel)>,
    sources: Vec<(StructureId, usize, usize)>,
}

impl OutlineKey {
    pub fn of(
        scene: &Scene,
        sources: &[DrawSource],
        highlight: Option<&Highlight>,
        hover: Option<(StructureId, u32)>,
        level: SelectLevel,
    ) -> Self {
        Self {
            mask: scene
                .active_selection()
                .map(|a| Arc::as_ptr(&a.mask) as usize),
            highlight: highlight.map(Highlight::key),
            hover: hover.map(|(id, atom)| (id, atom, level)),
            sources: sources
                .iter()
                .map(|s| (s.id, address(&s.atom_map), address(&s.bond_atoms)))
                .collect(),
        }
    }
}

/// Rebuilds `renderer`'s selection bits when the selection, the row
/// highlight, the hovered atom, the pick level or the draw list changed
/// since `last`.
pub fn sync(
    renderer: &mut Renderer,
    last: &mut OutlineKey,
    scene: &Scene,
    sources: &[DrawSource],
    highlight: Option<&Highlight>,
    hover: Option<(StructureId, u32)>,
    level: SelectLevel,
) {
    let key = OutlineKey::of(scene, sources, highlight, hover, level);
    if key != *last {
        renderer.set_selection(item_selections(scene, sources, highlight, hover, level));
        *last = key;
    }
}

fn address<T>(map: &Option<Arc<T>>) -> usize {
    map.as_ref().map_or(0, |m| Arc::as_ptr(m) as usize)
}

/// One entry per draw source, in pick order; empty when nothing is
/// selected, highlighted or hovered.
fn item_selections(
    scene: &Scene,
    sources: &[DrawSource],
    highlight: Option<&Highlight>,
    hover: Option<(StructureId, u32)>,
    level: SelectLevel,
) -> Vec<ItemSelection> {
    let active = scene.active_selection();
    if active.is_none() && highlight.is_none() && hover.is_none() {
        return Vec::new();
    }
    // Once, not per draw source: a molecule-level pick walks the bond graph.
    let picked = hover.and_then(|(id, atom)| {
        let loaded = scene.structure(id)?;
        Some((id, picked_atoms(loaded, atom, level)))
    });
    sources
        .iter()
        .map(|source| {
            let Some(loaded) = scene.structure(source.id) else {
                return ItemSelection::default();
            };
            let marks = Marks {
                active: active.filter(|a| a.structure == source.id),
                whole: highlight.is_some_and(|h| h.covers_source(source)),
                set: highlight.and_then(|h| h.atoms_of(source.id)),
                hover: picked
                    .as_ref()
                    .filter(|(id, _)| *id == source.id)
                    .map(|(_, mask)| mask),
            };
            if marks.is_empty() {
                return ItemSelection::default();
            }
            item_selection(loaded, &marks, source)
        })
        .collect()
}

/// The reasons one structure's atoms are outlined.
struct Marks<'a> {
    active: Option<&'a ActiveSelection>,
    /// Every atom of the source (a hovered rep or structure row).
    whole: bool,
    set: Option<&'a Mask>,
    /// What a click on the hovered atom would pick.
    hover: Option<&'a FixedBitSet>,
}

impl Marks<'_> {
    fn is_empty(&self) -> bool {
        self.active.is_none() && !self.whole && self.set.is_none() && self.hover.is_none()
    }

    fn marks(&self, atom: u32) -> bool {
        let at = atom as usize;
        self.whole
            || self.active.is_some_and(|a| a.mask.contains(at))
            || self.set.is_some_and(|m| m.contains(at))
            || self.hover.is_some_and(|h| h.contains(at))
    }
}

/// A bond counts as selected when both of its atoms are.
fn item_selection(
    loaded: &LoadedStructure,
    marks: &Marks<'_>,
    source: &DrawSource,
) -> ItemSelection {
    let atom_count = source
        .atom_map
        .as_ref()
        .map_or(loaded.structure.atom_count(), |m| m.len());
    let pairs = source
        .bond_atoms
        .as_deref()
        .map_or(loaded.bonds.pairs.as_slice(), |p| p.as_slice());
    ItemSelection {
        atoms: ItemSelection::bits(atom_count, |i| marks.marks(source.atom(i as u32))),
        bonds: ItemSelection::bits(pairs.len(), |k| {
            let [a, b] = pairs[k];
            marks.marks(a) && marks.marks(b)
        }),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use vv_scene::{Command, CommandHistory};

    fn source(id: StructureId, rep: RepId) -> DrawSource {
        DrawSource {
            id,
            rep,
            atom_map: None,
            bond_atoms: None,
        }
    }

    fn loaded_1ake() -> (Scene, StructureId) {
        let mut scene = Scene::new();
        let mut history = CommandHistory::new(10);
        let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../fixtures/small/1AKE.pdb");
        history
            .dispatch(&mut scene, Command::LoadStructure { path })
            .unwrap();
        let id = scene.structures().next().unwrap().0;
        (scene, id)
    }

    fn outlined(picked: &[ItemSelection]) -> u32 {
        picked[0].atoms.iter().map(|w| w.count_ones()).sum()
    }

    #[test]
    fn a_highlighted_source_selects_every_atom_regardless_of_the_active_mask() {
        let (scene, id) = loaded_1ake();
        let rep = scene.structure(id).unwrap().rep().id;
        let other_rep = RepId(rep.0 + 1);
        let sources = vec![source(id, rep), source(id, other_rep)];
        let level = SelectLevel::Atom;

        // Nothing selected, nothing highlighted: no work, no sources.
        assert!(item_selections(&scene, &sources, None, None, level).is_empty());

        // Highlighting a rep selects its own source fully, and leaves the
        // other rep's source untouched -- distinguished by `DrawSource::rep`,
        // not just the structure id.
        let highlight = Highlight::Rep(id, rep);
        let picked = item_selections(&scene, &sources, Some(&highlight), None, level);
        assert_eq!(picked.len(), 2);
        let atom_count = scene.structure(id).unwrap().structure.atom_count();
        assert_eq!(picked[0].atoms, ItemSelection::bits(atom_count, |_| true));
        assert_eq!(
            picked[1],
            ItemSelection::default(),
            "not highlighted, not active: the cheap empty default, not a zero-filled bitset"
        );
    }

    #[test]
    fn hovering_outlines_what_a_click_would_pick_at_the_level() {
        let (scene, id) = loaded_1ake();
        let sources = vec![source(id, scene.structure(id).unwrap().rep().id)];
        let at = |level| {
            outlined(&item_selections(
                &scene,
                &sources,
                None,
                Some((id, 0)),
                level,
            ))
        };
        let total = scene.structure(id).unwrap().structure.atom_count() as u32;
        let (atom, residue, chain) = (
            at(SelectLevel::Atom),
            at(SelectLevel::Residue),
            at(SelectLevel::Chain),
        );
        assert_eq!(atom, 1);
        assert!(atom < residue && residue < chain && chain < total);
    }

    #[test]
    fn structure_and_set_rows_outline_their_atoms() {
        let (scene, id) = loaded_1ake();
        let sources = vec![source(id, scene.structure(id).unwrap().rep().id)];
        let total = scene.structure(id).unwrap().structure.atom_count();
        let level = SelectLevel::Atom;
        let row = Highlight::Structure(id);
        let all = item_selections(&scene, &sources, Some(&row), None, level);
        assert_eq!(outlined(&all) as usize, total);

        let mut few = FixedBitSet::with_capacity(total);
        few.insert_range(0..7);
        let set = Highlight::Atoms(id, Arc::new(few));
        let some = item_selections(&scene, &sources, Some(&set), None, level);
        assert_eq!(outlined(&some), 7);
    }
}
