//! Computed annotations, kept until their structure or frame changes so
//! the strip only pays for drawing on ordinary frames. Results that take
//! long (solvent accessibility, reading a file's sequences, UniProt) are
//! computed on worker threads (`background.rs`) and appear when done.

use std::collections::HashMap;
use std::sync::Arc;

use egui::Color32;
use vv_core::seqfeat::EntityChain;
use vv_scene::{ColorScheme, LoadedStructure, Scene, StructureId};

use super::background::{Jobs, Progress};
use super::color::{residue_colors, ColorInput, SeqColor};
use super::peers::{self, Peers};
use super::rows::{chain_rows, ligand_groups, LigandGroup};
use super::tracks::{AntibodySettings, Extras, TrackContext, TrackData, TrackProvider};
use super::uniprot;

/// Structures larger than this skip solvent accessibility: its cost grows
/// with atom count.
const SASA_ATOM_LIMIT: usize = 1_000_000;

/// The shape of a structure at one frame: what invalidates a result.
#[derive(Clone, Copy, PartialEq, Eq)]
struct Stamp {
    residues: usize,
    atoms: usize,
    frame: usize,
    antibody: Option<AntibodySettings>,
    /// Which background inputs had arrived, one bit each.
    ready: u8,
    peers: u64,
}

impl Stamp {
    /// `frame` is the frame number the result depends on, `0` if none.
    fn of(loaded: &LoadedStructure, frame: usize) -> Self {
        Self {
            residues: loaded.structure.topology.residue_count(),
            atoms: loaded.structure.topology.atom_count(),
            frame,
            antibody: None,
            ready: 0,
            peers: 0,
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

/// What the strip needs beyond the structure to compute a track.
pub struct Env {
    pub antibody: AntibodySettings,
    peers: Arc<Peers>,
    peers_stamp: u64,
}

/// A result from a worker thread, as far as it has got.
pub enum Fetch<T> {
    Pending,
    /// It cannot be had: too large, no source, offline, or the worker died.
    Unavailable,
    Ready(Arc<Option<T>>),
}

impl<T> Fetch<T> {
    pub fn value(&self) -> Option<&T> {
        match self {
            Fetch::Ready(v) => v.as_ref().as_ref(),
            _ => None,
        }
    }

    pub fn settled(&self) -> bool {
        !matches!(self, Fetch::Pending)
    }
}

type Job<T> = Jobs<StructureId, Stamp, Option<T>>;

#[derive(Default)]
pub struct Cache {
    tracks: HashMap<(StructureId, &'static str), Cached<Arc<TrackData>>>,
    colors: HashMap<StructureId, Cached<ColorEntry>>,
    ligands: HashMap<StructureId, Cached<Arc<Vec<LigandGroup>>>>,
    sasa: Job<Vec<f32>>,
    entity: Job<HashMap<String, EntityChain>>,
    uniprot: Job<uniprot::Data>,
    peers: Arc<Peers>,
    peers_stamp: u64,
    wake: Option<egui::Context>,
}

fn settle<T: Send + 'static, W: FnOnce() -> Option<T> + Send + 'static>(
    jobs: &mut Job<T>,
    wake: Option<&egui::Context>,
    id: StructureId,
    stamp: Stamp,
    prepare: impl FnOnce() -> W,
) -> Fetch<T> {
    match jobs.get(id, stamp, wake, prepare) {
        Some(value) => Fetch::Ready(value),
        None if jobs.progress(id, stamp) == Progress::Failed => Fetch::Unavailable,
        None => Fetch::Pending,
    }
}

impl Cache {
    /// The context worker threads wake when they finish.
    pub fn set_wake(&mut self, ctx: &egui::Context) {
        if self.wake.is_none() {
            self.wake = Some(ctx.clone());
        }
    }

    /// Forgets structures that are no longer open.
    pub fn retain_open(&mut self, scene: &Scene) {
        let open = |id: &StructureId| scene.structure(*id).is_some();
        self.tracks.retain(|(id, _), _| open(id));
        self.colors.retain(|id, _| open(id));
        self.ligands.retain(|id, _| open(id));
        self.sasa.retain(open);
        self.entity.retain(open);
        self.uniprot.retain(open);
    }

    /// The structure's non-polymer molecules by name.
    pub fn ligands(&mut self, id: StructureId, loaded: &LoadedStructure) -> Arc<Vec<LigandGroup>> {
        let stamp = Stamp::of(loaded, 0);
        if let Some(hit) = self.ligands.get(&id).filter(|c| c.stamp == stamp) {
            return hit.value.clone();
        }
        let value = Arc::new(ligand_groups(&loaded.structure.topology));
        self.ligands.insert(
            id,
            Cached {
                stamp,
                value: value.clone(),
            },
        );
        value
    }

    /// The environment tracks are computed in; the other chains are only
    /// gathered when `with_peers` (a track compares against them).
    pub fn env(&mut self, scene: &Scene, antibody: AntibodySettings, with_peers: bool) -> Env {
        if with_peers {
            let now = peers::fingerprint(scene);
            if now != self.peers_stamp {
                self.peers = Arc::new(Peers::of(scene));
                self.peers_stamp = now;
            }
        }
        Env {
            antibody,
            peers: self.peers.clone(),
            peers_stamp: if with_peers { self.peers_stamp } else { 0 },
        }
    }

    pub fn track(
        &mut self,
        id: StructureId,
        loaded: &LoadedStructure,
        provider: &dyn TrackProvider,
        env: &Env,
    ) -> Arc<TrackData> {
        let inputs = provider.inputs();
        let sasa = if inputs.sasa {
            self.rel_sasa(id, loaded)
        } else {
            Fetch::Unavailable
        };
        let entity = if inputs.entity {
            self.entity(id, loaded)
        } else {
            Fetch::Unavailable
        };
        let uniprot = if inputs.uniprot {
            self.uniprot(id, loaded)
        } else {
            Fetch::Unavailable
        };
        let mut stamp = Stamp::of(loaded, provider.frame_key(loaded, loaded.frame));
        stamp.antibody = provider.uses_antibody_settings().then_some(env.antibody);
        stamp.ready = u8::from(sasa.settled() && inputs.sasa)
            | u8::from(entity.settled() && inputs.entity) << 1
            | u8::from(uniprot.settled() && inputs.uniprot) << 2;
        stamp.peers = if inputs.peers { env.peers_stamp } else { 0 };
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
            antibody: env.antibody,
            extras: Extras {
                rel_sasa: sasa.value().map(|v| v.as_slice()),
                entity: entity.value(),
                uniprot: uniprot.value(),
                peers: inputs.peers.then_some(&*env.peers),
                structure: Some(id),
                model: loaded.frame as u32 + 1,
            },
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

    /// Relative SASA per residue; computed on a worker thread, and never
    /// for a structure too large to.
    fn rel_sasa(&mut self, id: StructureId, loaded: &LoadedStructure) -> Fetch<Vec<f32>> {
        let top = loaded.structure.topology.clone();
        if top.atom_count() > SASA_ATOM_LIMIT {
            return Fetch::Unavailable;
        }
        let stamp = Stamp::of(loaded, loaded.frame);
        settle(&mut self.sasa, self.wake.as_ref(), id, stamp, || {
            let positions = loaded.structure.frame(loaded.frame).positions().to_vec();
            move || {
                Some(vv_core::seqfeat::relative_sasa(
                    &top,
                    &positions,
                    vv_core::ses::WATER_PROBE,
                ))
            }
        })
    }

    /// The chains' full sequences, read from the structure's file on a
    /// worker thread.
    fn entity(
        &mut self,
        id: StructureId,
        loaded: &LoadedStructure,
    ) -> Fetch<HashMap<String, EntityChain>> {
        let Some(path) = loaded.path.as_deref() else {
            return Fetch::Unavailable;
        };
        let stamp = Stamp::of(loaded, 0);
        settle(&mut self.entity, self.wake.as_ref(), id, stamp, || {
            let path = path.to_path_buf();
            move || vv_io::seqdata::entity::read(&path)
        })
    }

    /// UniProt features for the structure's entry, fetched on a worker
    /// thread; unavailable offline unless cached.
    fn uniprot(&mut self, id: StructureId, loaded: &LoadedStructure) -> Fetch<uniprot::Data> {
        let path = loaded.path.as_deref();
        let entry = uniprot::pdb_id(&loaded.structure.topology, path);
        let (Some(path), Some(entry)) = (path, entry) else {
            return Fetch::Unavailable;
        };
        let stamp = Stamp::of(loaded, 0);
        settle(&mut self.uniprot, self.wake.as_ref(), id, stamp, || {
            let path = path.to_path_buf();
            move || uniprot::load(&entry, &path)
        })
    }

    /// Where the UniProt fetch for a structure stands, starting it if it
    /// has not begun.
    pub fn uniprot_status(&mut self, id: StructureId, loaded: &LoadedStructure) -> Fetch<()> {
        match self.uniprot(id, loaded) {
            Fetch::Pending => Fetch::Pending,
            Fetch::Ready(data) if data.is_some() => Fetch::Ready(Arc::new(Some(()))),
            _ => Fetch::Unavailable,
        }
    }

    /// One color per residue under `scheme`, or `None` when it draws none
    /// (or, for burial, until the accessible areas are computed).
    pub fn colors(
        &mut self,
        id: StructureId,
        loaded: &LoadedStructure,
        scheme: SeqColor,
    ) -> Option<Arc<Vec<Color32>>> {
        if scheme == SeqColor::None {
            return None;
        }
        let frame = if scheme.per_frame() { loaded.frame } else { 0 };
        let stamp = Stamp::of(loaded, frame);
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
        let rel = match scheme {
            SeqColor::Sasa => match self.rel_sasa(id, loaded) {
                Fetch::Pending => return None,
                other => Some(other),
            },
            _ => None,
        };
        let coords = loaded.structure.frame(loaded.frame);
        let input = ColorInput {
            frame: loaded.frame,
            loaded,
            positions: coords.positions(),
            rel_sasa: rel.as_ref().and_then(|f| f.value()).map(|v| v.as_slice()),
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

    /// Relative SASA of residue `r`: starts the computation if needed and
    /// answers once it is done.
    pub fn sasa_at(&mut self, id: StructureId, loaded: &LoadedStructure, r: u32) -> Option<f32> {
        let fetch = self.rel_sasa(id, loaded);
        let v = *fetch.value()?.get(r as usize)?;
        (!v.is_nan()).then_some(v)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sequence::color::SCHEMES;
    use std::time::{Duration, Instant};
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

    /// Calls `step` until it yields, for results that arrive from a worker.
    fn until<T>(mut step: impl FnMut() -> Option<T>) -> T {
        let start = Instant::now();
        loop {
            if let Some(v) = step() {
                return v;
            }
            assert!(start.elapsed() < Duration::from_secs(60), "never arrived");
            std::thread::sleep(Duration::from_millis(5));
        }
    }

    #[test]
    fn every_scheme_colors_some_residues_and_is_cached() {
        let (scene, id) = open("small/1CRN.pdb");
        let loaded = scene.structure(id).unwrap();
        let n = loaded.structure.topology.residue_count();
        let mut cache = Cache::default();
        for (scheme, word, _) in SCHEMES {
            let colors = until(|| {
                let c = cache.colors(id, loaded, scheme);
                (c.is_some() || scheme == SeqColor::None).then_some(c)
            });
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
    fn burial_is_computed_off_the_calling_thread_then_answers_the_tooltip() {
        let (scene, id) = open("small/1CRN.pdb");
        let loaded = scene.structure(id).unwrap();
        let mut cache = Cache::default();
        assert_eq!(cache.sasa_at(id, loaded, 0), None);
        let exposed = until(|| cache.sasa_at(id, loaded, 0));
        assert!((0.0..=1.5).contains(&exposed));
    }

    #[test]
    fn the_burial_track_fills_in_when_the_areas_arrive() {
        let (scene, id) = open("small/1CRN.pdb");
        let loaded = scene.structure(id).unwrap();
        let mut cache = Cache::default();
        let env = cache.env(&scene, Default::default(), false);
        let provider = crate::sequence::tracks::provider("burial").unwrap();
        let first = cache.track(id, loaded, provider, &env);
        let all = 0..loaded.structure.topology.residue_count() as u32;
        assert!(!first.any_in(&all));
        let ready = until(|| {
            let t = cache.track(id, loaded, provider, &env);
            t.any_in(&all).then_some(t)
        });
        assert!((0..46).all(|r| ready.kind(r) != 0));
    }

    #[test]
    fn tracks_are_computed_once_and_dropped_with_their_structure() {
        let (mut scene, id) = open("small/1CRN.pdb");
        let mut cache = Cache::default();
        let env = cache.env(&scene, Default::default(), false);
        let provider = crate::sequence::tracks::provider("disulfide").unwrap();
        let first = cache.track(id, scene.structure(id).unwrap(), provider, &env);
        let second = cache.track(id, scene.structure(id).unwrap(), provider, &env);
        assert!(Arc::ptr_eq(&first, &second));
        let mut history = CommandHistory::new(10);
        history
            .dispatch(&mut scene, Command::CloseStructure { id })
            .unwrap();
        cache.retain_open(&scene);
        assert!(cache.tracks.is_empty());
    }
}
