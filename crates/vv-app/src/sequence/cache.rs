//! Computed annotations, kept until their structure or frame changes so
//! the strip only pays for drawing on ordinary frames.

use std::collections::HashMap;
use std::sync::Arc;

use egui::Color32;
use vv_scene::{ColorScheme, LoadedStructure, Scene, StructureId};

use super::color::{residue_colors, ColorInput, SeqColor};
use super::rows::chain_rows;
use super::tracks::{AntibodySettings, TrackContext, TrackData, TrackProvider};

/// Structures larger than this skip the solvent-accessibility scheme: its
/// cost grows with atom count and runs on the UI thread.
const SASA_ATOM_LIMIT: usize = 250_000;

/// The shape of a structure at one frame: what invalidates a result.
#[derive(Clone, Copy, PartialEq, Eq)]
struct Stamp {
    residues: usize,
    atoms: usize,
    frame: usize,
    antibody: Option<AntibodySettings>,
}

impl Stamp {
    /// `frame` counts only when the result depends on it.
    fn of(loaded: &LoadedStructure, frame: usize, per_frame: bool) -> Self {
        Self {
            residues: loaded.structure.topology.residue_count(),
            atoms: loaded.structure.topology.atom_count(),
            frame: if per_frame { frame } else { 0 },
            antibody: None,
        }
    }
}

struct Cached<T> {
    stamp: Stamp,
    value: T,
}

struct ColorEntry {
    scheme: SeqColor,
    view: ColorScheme,
    colors: Arc<Vec<Color32>>,
}

#[derive(Default)]
pub struct Cache {
    tracks: HashMap<(StructureId, &'static str), Cached<Arc<TrackData>>>,
    colors: HashMap<StructureId, Cached<ColorEntry>>,
    sasa: HashMap<StructureId, Cached<Arc<Vec<f32>>>>,
}

impl Cache {
    /// Forgets structures that are no longer open.
    pub fn retain_open(&mut self, scene: &Scene) {
        let open = |id: &StructureId| scene.structure(*id).is_some();
        self.tracks.retain(|(id, _), _| open(id));
        self.colors.retain(|id, _| open(id));
        self.sasa.retain(|id, _| open(id));
    }

    pub fn track(
        &mut self,
        id: StructureId,
        loaded: &LoadedStructure,
        provider: &dyn TrackProvider,
        antibody: AntibodySettings,
    ) -> Arc<TrackData> {
        let mut stamp = Stamp::of(loaded, loaded.frame, provider.per_frame());
        stamp.antibody = provider.uses_antibody_settings().then_some(antibody);
        let key = (id, provider.id());
        if let Some(hit) = self.tracks.get(&key).filter(|c| c.stamp == stamp) {
            return hit.value.clone();
        }
        let coords = loaded.structure.frame(loaded.frame);
        let rows = chain_rows(&loaded.structure.topology);
        let ctx = TrackContext {
            loaded,
            positions: coords.positions(),
            rows: &rows,
            antibody,
        };
        let value = Arc::new(provider.compute(&ctx));
        self.tracks.insert(
            key,
            Cached {
                stamp,
                value: value.clone(),
            },
        );
        value
    }

    /// Relative SASA per residue; `None` for a structure too large to
    /// compute on the UI thread.
    fn rel_sasa(&mut self, id: StructureId, loaded: &LoadedStructure) -> Option<Arc<Vec<f32>>> {
        let top = &loaded.structure.topology;
        if top.atom_count() > SASA_ATOM_LIMIT {
            return None;
        }
        let stamp = Stamp::of(loaded, loaded.frame, true);
        if let Some(hit) = self.sasa.get(&id).filter(|c| c.stamp == stamp) {
            return Some(hit.value.clone());
        }
        let coords = loaded.structure.frame(loaded.frame);
        let value = Arc::new(vv_core::seqfeat::relative_sasa(
            top,
            coords.positions(),
            vv_core::ses::WATER_PROBE,
        ));
        self.sasa.insert(
            id,
            Cached {
                stamp,
                value: value.clone(),
            },
        );
        Some(value)
    }

    /// One color per residue under `scheme`, or `None` when it draws none.
    pub fn colors(
        &mut self,
        id: StructureId,
        loaded: &LoadedStructure,
        scheme: SeqColor,
    ) -> Option<Arc<Vec<Color32>>> {
        if scheme == SeqColor::None {
            return None;
        }
        let stamp = Stamp::of(loaded, loaded.frame, scheme.per_frame());
        let view = loaded.reps[loaded.current_rep.min(loaded.reps.len() - 1)]
            .coloring
            .clone();
        let fresh = |c: &&Cached<ColorEntry>| {
            c.stamp == stamp
                && c.value.scheme == scheme
                && (scheme != SeqColor::View || c.value.view == view)
        };
        if let Some(hit) = self.colors.get(&id).filter(fresh) {
            return Some(hit.value.colors.clone());
        }
        let rel = (scheme == SeqColor::Sasa)
            .then(|| self.rel_sasa(id, loaded))
            .flatten();
        let coords = loaded.structure.frame(loaded.frame);
        let input = ColorInput {
            frame: loaded.frame,
            loaded,
            positions: coords.positions(),
            rel_sasa: rel.as_deref().map(|v| v.as_slice()),
        };
        let colors = Arc::new(residue_colors(scheme, &input));
        let entry = ColorEntry {
            scheme,
            view,
            colors: colors.clone(),
        };
        self.colors.insert(
            id,
            Cached {
                stamp,
                value: entry,
            },
        );
        Some(colors)
    }

    /// The cached relative SASA of residue `r`, if already computed.
    pub fn cached_sasa(&self, id: StructureId, r: u32) -> Option<f32> {
        let v = *self.sasa.get(&id)?.value.get(r as usize)?;
        (!v.is_nan()).then_some(v)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sequence::color::SCHEMES;
    use vv_scene::{Command, CommandHistory};

    fn open(fixture: &str) -> (Scene, StructureId) {
        let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../fixtures")
            .join(fixture);
        let mut scene = Scene::new();
        let mut history = CommandHistory::new(10);
        history
            .dispatch(&mut scene, Command::LoadStructure { path })
            .unwrap();
        let id = scene.structures().next().unwrap().0;
        (scene, id)
    }

    #[test]
    fn every_scheme_colors_some_residues_and_is_cached() {
        let (scene, id) = open("small/1CRN.pdb");
        let loaded = scene.structure(id).unwrap();
        let n = loaded.structure.topology.residue_count();
        let mut cache = Cache::default();
        for (scheme, word, _) in SCHEMES {
            let colors = cache.colors(id, loaded, scheme);
            if scheme == SeqColor::None {
                assert!(colors.is_none());
                continue;
            }
            let colors = colors.unwrap();
            assert_eq!(colors.len(), n, "{word}");
            assert!(colors.iter().any(|c| *c != Color32::TRANSPARENT), "{word}");
            let again = cache.colors(id, loaded, scheme).unwrap();
            assert!(Arc::ptr_eq(&colors, &again), "{word} recomputed");
        }
    }

    #[test]
    fn burial_records_relative_sasa_for_the_tooltip() {
        let (scene, id) = open("small/1CRN.pdb");
        let loaded = scene.structure(id).unwrap();
        let mut cache = Cache::default();
        assert_eq!(cache.cached_sasa(id, 0), None);
        cache.colors(id, loaded, SeqColor::Sasa);
        assert!(cache.cached_sasa(id, 0).is_some());
    }

    #[test]
    fn tracks_are_computed_once_and_dropped_with_their_structure() {
        let (mut scene, id) = open("small/1CRN.pdb");
        let mut cache = Cache::default();
        let provider = crate::sequence::tracks::provider("disulfide").unwrap();
        let first = cache.track(
            id,
            scene.structure(id).unwrap(),
            provider,
            Default::default(),
        );
        let second = cache.track(
            id,
            scene.structure(id).unwrap(),
            provider,
            Default::default(),
        );
        assert!(Arc::ptr_eq(&first, &second));
        let mut history = CommandHistory::new(10);
        history
            .dispatch(&mut scene, Command::CloseStructure { id })
            .unwrap();
        cache.retain_open(&scene);
        assert!(cache.tracks.is_empty());
    }
}
