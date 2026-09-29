//! The viewport's selection outline: which drawn pick ids belong to the
//! active selection, a Structures-panel row's hover/current highlight, or
//! the hovered atom's residue (`Renderer::set_selection`) -- all three
//! share the same outline, unioned per atom.

use std::sync::Arc;

use vv_render::{ItemSelection, Renderer};
use vv_scene::{ActiveSelection, LoadedStructure, RepId, Scene, StructureId};

use crate::gpu_cache::DrawSource;

/// What the renderer's selection bits were last built from: the active
/// mask by address, the row highlight, the hovered atom, and each draw
/// source's structure and maps by address.
#[derive(Default, PartialEq)]
pub struct OutlineKey {
    mask: Option<usize>,
    highlight: Option<(StructureId, RepId)>,
    hover: Option<(StructureId, u32)>,
    sources: Vec<(StructureId, usize, usize)>,
}

impl OutlineKey {
    pub fn of(
        scene: &Scene,
        sources: &[DrawSource],
        highlight: Option<(StructureId, RepId)>,
        hover: Option<(StructureId, u32)>,
    ) -> Self {
        Self {
            mask: scene
                .active_selection()
                .map(|a| Arc::as_ptr(&a.mask) as usize),
            highlight,
            hover,
            sources: sources
                .iter()
                .map(|s| (s.id, address(&s.atom_map), address(&s.bond_atoms)))
                .collect(),
        }
    }
}

/// Rebuilds `renderer`'s selection bits when the selection, the row
/// highlight, the hovered atom, or the draw list changed since `last`.
pub fn sync(
    renderer: &mut Renderer,
    last: &mut OutlineKey,
    scene: &Scene,
    sources: &[DrawSource],
    highlight: Option<(StructureId, RepId)>,
    hover: Option<(StructureId, u32)>,
) {
    let key = OutlineKey::of(scene, sources, highlight, hover);
    if key != *last {
        renderer.set_selection(item_selections(scene, sources, highlight, hover));
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
    highlight: Option<(StructureId, RepId)>,
    hover: Option<(StructureId, u32)>,
) -> Vec<ItemSelection> {
    let active = scene.active_selection();
    if active.is_none() && highlight.is_none() && hover.is_none() {
        return Vec::new();
    }
    sources
        .iter()
        .map(|source| {
            let highlighted = highlight == Some((source.id, source.rep));
            let active = active.filter(|a| a.structure == source.id);
            let hover_atom = hover.filter(|&(id, _)| id == source.id).map(|(_, a)| a);
            if active.is_none() && !highlighted && hover_atom.is_none() {
                return ItemSelection::default();
            }
            match scene.structure(source.id) {
                Some(loaded) => item_selection(loaded, active, highlighted, hover_atom, source),
                None => ItemSelection::default(),
            }
        })
        .collect()
}

/// A bond counts as selected when both of its atoms are.
fn item_selection(
    loaded: &LoadedStructure,
    active: Option<&ActiveSelection>,
    highlighted: bool,
    hover_atom: Option<u32>,
    source: &DrawSource,
) -> ItemSelection {
    // Hovering outlines the whole residue, not just the one atom under
    // the pointer, so it reads as "this is what a click would pick".
    let hover_range =
        hover_atom.map(|atom| crate::viewport::residue_atoms(&loaded.structure.topology, atom));
    let selected = |atom: u32| {
        highlighted
            || active.is_some_and(|a| a.mask.contains(atom as usize))
            || hover_range.as_ref().is_some_and(|r| r.contains(&atom))
    };
    let atom_count = source
        .atom_map
        .as_ref()
        .map_or(loaded.structure.atom_count(), |m| m.len());
    let pairs = source
        .bond_atoms
        .as_deref()
        .map_or(loaded.bonds.pairs.as_slice(), |p| p.as_slice());
    ItemSelection {
        atoms: ItemSelection::bits(atom_count, |i| selected(source.atom(i as u32))),
        bonds: ItemSelection::bits(pairs.len(), |k| {
            let [a, b] = pairs[k];
            selected(a) && selected(b)
        }),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn source(id: StructureId, rep: RepId) -> DrawSource {
        DrawSource {
            id,
            rep,
            atom_map: None,
            bond_atoms: None,
        }
    }

    #[test]
    fn a_highlighted_source_selects_every_atom_regardless_of_the_active_mask() {
        use vv_scene::{Command, CommandHistory};
        let mut scene = Scene::new();
        let mut history = CommandHistory::new(10);
        history
            .dispatch(
                &mut scene,
                Command::LoadStructure {
                    path: std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
                        .join("../../fixtures/small/1CRN.pdb"),
                },
            )
            .unwrap();
        let (id, loaded) = scene.structures().next().unwrap();
        let rep = loaded.rep().id;
        let other_rep = RepId(rep.0 + 1);
        let sources = vec![source(id, rep), source(id, other_rep)];

        // Nothing selected, nothing highlighted: no work, no sources.
        assert!(item_selections(&scene, &sources, None, None).is_empty());

        // Highlighting a rep selects its own source fully, and leaves the
        // other rep's source untouched -- distinguished by `DrawSource::rep`,
        // not just the structure id.
        let picked = item_selections(&scene, &sources, Some((id, rep)), None);
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
    fn hovering_an_atom_outlines_its_whole_residue() {
        use vv_scene::{Command, CommandHistory};
        let mut scene = Scene::new();
        let mut history = CommandHistory::new(10);
        history
            .dispatch(
                &mut scene,
                Command::LoadStructure {
                    path: std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
                        .join("../../fixtures/small/1CRN.pdb"),
                },
            )
            .unwrap();
        let (id, loaded) = scene.structures().next().unwrap();
        let rep = loaded.rep().id;
        let sources = vec![source(id, rep)];
        let residue = crate::viewport::residue_atoms(&loaded.structure.topology, 0);
        assert!(residue.len() > 1, "fixture needs a multi-atom residue");

        let picked = item_selections(&scene, &sources, None, Some((id, 0)));
        let expected = ItemSelection::bits(loaded.structure.atom_count(), |i| {
            residue.contains(&(i as u32))
        });
        assert_eq!(picked[0].atoms, expected);
    }
}
