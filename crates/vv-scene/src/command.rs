//! Commands: the only way `Scene` changes. Every command's `apply` returns
//! its own inverse, computed from the state it just overwrote — undo never
//! snapshots the whole scene, it replays one command.
//!
//! A few variants (`RestoreStructure`, `RestoreSelection`,
//! `RestoreSelectionSet`) are never constructed by callers; they only ever
//! appear as an inverse. They are still ordinary `Command`s so undo/redo
//! needs no special case.

use std::path::{Path, PathBuf};
use std::sync::Arc;

use crate::scene::{
    ActiveSelection, Caption, ColorScheme, LoadedStructure, Material, Measurement, Rep, RepId,
    Representation, Scene, StructureId,
};
use crate::selection::{Mask, SelectionSet};
use crate::values::ValueChannel;
use vv_core::Structure;

#[derive(thiserror::Error, Debug)]
pub enum SceneError {
    #[error("no structure loaded at {0:?}")]
    NoSuchStructure(StructureId),
    #[error("structure {0:?} has no rep {1:?}")]
    NoSuchRep(StructureId, RepId),
    #[error("a structure keeps at least one rep; hide it instead")]
    LastRep(StructureId),
    #[error("no active selection to save")]
    NoActiveSelection,
    #[error("no selection set named `{0}`")]
    NoSuchSelectionSet(String),
    #[error("a selection set named `{0}` already exists")]
    SelectionSetNameTaken(String),
    #[error("failed to load {path}: {source}")]
    Load {
        path: PathBuf,
        source: vv_io::ParseError,
    },
    #[error("failed to read trajectory {path}: {source}")]
    LoadTrajectory {
        path: PathBuf,
        source: vv_io::trajectory::TrajectoryError,
    },
    #[error(
        "trajectory has {trajectory} atoms per frame, topology {topology} has {topology_atoms}"
    )]
    TrajectoryAtomMismatch {
        topology: PathBuf,
        topology_atoms: usize,
        trajectory: usize,
    },
    #[error("{0}")]
    Structure(#[from] vv_core::StructureError),
    #[error("bad selection `{expr}`: {source}")]
    Selection {
        expr: String,
        source: vv_core::SelectError,
    },
    #[error("value channel `{name}` has {actual} atoms, structure has {expected}")]
    ValuesLength {
        name: String,
        expected: usize,
        actual: usize,
    },
    #[error("no value channel named `{0}`")]
    NoSuchValues(String),
    #[error("atom {atom} is out of range: structure {id:?} has {atoms} atoms")]
    NoSuchAtom {
        id: StructureId,
        atom: u32,
        atoms: usize,
    },
    #[error("no caption named `{0}`")]
    NoSuchCaption(String),
    #[error("frame {frame} out of range: {id:?} has {count} frame(s)")]
    FrameOutOfRange {
        id: StructureId,
        frame: usize,
        count: usize,
    },
}

#[derive(Clone, Debug)]
pub enum Command {
    LoadStructure {
        path: PathBuf,
    },
    /// Loads a topology file plus a trajectory as one multi-frame
    /// structure, its frames streamed from the file: every frame lands in
    /// the same `Scene` slot `LoadStructure` would use for a single-frame
    /// file, so the Timeline tab and `frame()`/`frame_count()` need no
    /// special case. The topology is PDB, mmCIF, PSF or PRMTOP
    /// (`vv_io::load_topology`); the trajectory is DCD, XTC, TRR, AMBER
    /// NetCDF, or a PDB/mmCIF coordinate file read as one frame per model
    /// (`vv_io::trajectory`). A PSF/PRMTOP topology's own bond list is
    /// used verbatim (`vv_core::bonds` module doc).
    LoadTrajectory {
        topology: PathBuf,
        trajectory: PathBuf,
    },
    /// Replaces structure `id`'s coordinates with the frames of
    /// `trajectory`, streamed like `LoadTrajectory`. The topology, bonds,
    /// reps, selections, labels and id stay; the current frame goes back to
    /// 0. The trajectory must have the structure's atom count.
    AttachTrajectory {
        id: StructureId,
        trajectory: PathBuf,
    },
    CloseStructure {
        id: StructureId,
    },
    /// The structure's display label (Structures panel row). Doesn't
    /// touch its file `path`.
    SetStructureLabel {
        id: StructureId,
        label: String,
    },
    /// Whole-structure show/hide (Structures panel row's eye), independent
    /// of any rep's own `ShowRep`. See `LoadedStructure::visible`'s doc for
    /// what hiding skips.
    ShowStructure {
        id: StructureId,
        visible: bool,
    },
    SetRepresentation {
        id: StructureId,
        rep: RepId,
        representation: Representation,
    },
    SetColoring {
        id: StructureId,
        rep: RepId,
        coloring: ColorScheme,
    },
    SetMaterial {
        id: StructureId,
        rep: RepId,
        material: Material,
    },
    /// Inserts `rep` at `index` (clamped) and makes it current.
    AddRep {
        id: StructureId,
        index: usize,
        rep: Rep,
    },
    /// Refuses to remove a structure's last rep (hide it instead).
    RemoveRep {
        id: StructureId,
        rep: RepId,
    },
    /// The rep's atoms, as a selection expression (docs/SELECTION.md);
    /// checked for syntax when applied.
    SetRepSelection {
        id: StructureId,
        rep: RepId,
        selection: String,
    },
    ShowRep {
        id: StructureId,
        rep: RepId,
        visible: bool,
    },
    /// Sets option `name` of a rep (`Representation::options`), or back
    /// to its default (`None`).
    SetRepOption {
        id: StructureId,
        rep: RepId,
        name: String,
        value: Option<f32>,
    },
    SetCurrentRep {
        id: StructureId,
        rep: RepId,
    },
    /// Several commands as one undo step, applied in order (e.g. a style
    /// preset setting every structure's material). If one fails, the
    /// ones before it are undone and the error returned.
    Batch(Vec<Command>),
    /// Sets the current coordinate-set frame, undoably — for a
    /// deliberate, one-shot set from a script, the console, or Python;
    /// see [`LoadedStructure::frame`]'s doc for why playback/scrubbing
    /// use `Scene::set_frame_live` instead of this.
    SetFrame {
        id: StructureId,
        frame: usize,
    },
    /// Attach (`Some`) or remove (`None`) the per-atom value channel
    /// `name` on a structure. Undo restores what was there before.
    SetValues {
        id: StructureId,
        name: String,
        channel: Option<ValueChannel>,
    },
    /// Select exactly `mask` (a pick, or a mask computed elsewhere).
    Select {
        id: StructureId,
        mask: Mask,
    },
    /// Select the atoms matching a selection expression such as
    /// `chain A and name CA` (grammar: `vv_core::select`, docs/SELECTION.md).
    /// The expression is remembered on the selection.
    SelectExpr {
        id: StructureId,
        expr: String,
    },
    ClearSelection,
    SaveSelectionSet {
        name: String,
    },
    DeleteSelectionSet {
        name: String,
    },
    /// Renames a selection set (Selection panel). Refuses a collision
    /// with another set's name.
    RenameSelectionSet {
        old_name: String,
        new_name: String,
    },
    /// Sets (`Some`) or clears (`None`) the structure-anchored text label
    /// on one atom (docs/UI_DESIGN.md's "structure-anchored" annotation
    /// type). Billboarded in the viewport, following that atom's current
    /// frame position. Undo restores what was there before.
    SetLabel {
        id: StructureId,
        atom: u32,
        text: Option<String>,
    },
    /// Shows (`shown`) or removes a measurement on structure `id`. Undo
    /// puts back what was there.
    SetMeasurement {
        id: StructureId,
        measurement: Measurement,
        shown: bool,
    },
    /// Adds a free-floating screen-space text caption, not tied to any
    /// structure (docs/UI_DESIGN.md's "screen-space" annotation type).
    /// Overwrites any existing caption of the same name.
    SetCaption {
        caption: Caption,
    },
    DeleteCaption {
        name: String,
    },

    /// Inverse of `AttachTrajectory` (and of itself): puts back a
    /// coordinate source, its trajectory path and the current frame. Not
    /// for callers to construct.
    RestoreCoordinates {
        id: StructureId,
        structure: Box<Structure>,
        trajectory: Option<PathBuf>,
        frame: usize,
    },
    /// Inverse of `CloseStructure`. Not for callers to construct.
    RestoreStructure {
        id: StructureId,
        structure: Box<LoadedStructure>,
        prior_active: Option<ActiveSelection>,
    },
    /// Inverse of every selection change (`Select`, `SelectExpr`,
    /// `ClearSelection`): puts back the whole previous active selection,
    /// expression included. Not for callers to construct.
    RestoreSelection {
        active: Option<ActiveSelection>,
    },
    /// Inverse of `SaveSelectionSet`/`DeleteSelectionSet`. Not for callers
    /// to construct; overwrites (or creates) the set named `set.name`.
    RestoreSelectionSet {
        set: SelectionSet,
    },
    /// Inverse of `SetCaption`/`DeleteCaption`. Not for callers to
    /// construct; overwrites (or creates) the caption named `caption.name`.
    RestoreCaption {
        caption: Caption,
    },
}

/// `loaded`'s topology with the frames of `trajectory`, streamed.
fn streamed_onto(loaded: &LoadedStructure, trajectory: &Path) -> Result<Structure, SceneError> {
    let reader = vv_io::trajectory::Trajectory::open(trajectory).map_err(|source| {
        SceneError::LoadTrajectory {
            path: trajectory.to_path_buf(),
            source,
        }
    })?;
    let topology = loaded.structure.topology.clone();
    if reader.atom_count() != topology.atom_count() {
        return Err(SceneError::TrajectoryAtomMismatch {
            topology: loaded
                .path
                .clone()
                .unwrap_or_else(|| loaded.label.clone().into()),
            topology_atoms: topology.atom_count(),
            trajectory: reader.atom_count(),
        });
    }
    Ok(Structure::streamed(
        topology,
        Box::new(reader),
        vv_core::DEFAULT_FRAME_BUDGET,
    )?)
}

fn rep_mut(scene: &mut Scene, id: StructureId, rep: RepId) -> Result<&mut Rep, SceneError> {
    let loaded = scene
        .structures_mut()
        .get_mut(id)
        .ok_or(SceneError::NoSuchStructure(id))?;
    let index = loaded
        .rep_index(rep)
        .ok_or(SceneError::NoSuchRep(id, rep))?;
    Ok(&mut loaded.reps[index])
}

impl Command {
    /// Applies this command to `scene`, returning the command that undoes it.
    pub fn apply(self, scene: &mut Scene) -> Result<Command, SceneError> {
        let inverse = match self {
            Command::LoadStructure { path } => {
                let structure = vv_io::load(&path).map_err(|source| SceneError::Load {
                    path: path.clone(),
                    source,
                })?;
                let bonds =
                    vv_core::bonds::perceive(&structure.topology, structure.frame(0).positions());
                let label = structure_label(&structure.topology.id, &path);
                let loaded = LoadedStructure::new(structure, Some(path), label, bonds);
                let id = scene.structures_mut().insert(loaded);
                Command::CloseStructure { id }
            }

            Command::LoadTrajectory {
                topology,
                trajectory,
            } => {
                let base_topology =
                    vv_io::load_topology(&topology).map_err(|source| SceneError::Load {
                        path: topology.clone(),
                        source,
                    })?;
                let reader =
                    vv_io::trajectory::Trajectory::open(&trajectory).map_err(|source| {
                        SceneError::LoadTrajectory {
                            path: trajectory.clone(),
                            source,
                        }
                    })?;
                if reader.atom_count() != base_topology.atom_count() {
                    return Err(SceneError::TrajectoryAtomMismatch {
                        topology: topology.clone(),
                        topology_atoms: base_topology.atom_count(),
                        trajectory: reader.atom_count(),
                    });
                }
                // Frames are read as they are shown, so a trajectory larger
                // than memory plays.
                let structure = vv_core::Structure::streamed(
                    base_topology,
                    Box::new(reader),
                    vv_core::DEFAULT_FRAME_BUDGET,
                )?;
                let bonds =
                    vv_core::bonds::perceive(&structure.topology, structure.frame(0).positions());
                let label = structure_label(&structure.topology.id, &topology);
                let mut loaded = LoadedStructure::new(structure, Some(topology), label, bonds);
                loaded.trajectory = Some(trajectory);
                let id = scene.structures_mut().insert(loaded);
                Command::CloseStructure { id }
            }

            Command::AttachTrajectory { id, trajectory } => {
                let loaded = scene
                    .structures_mut()
                    .get_mut(id)
                    .ok_or(SceneError::NoSuchStructure(id))?;
                let structure = streamed_onto(loaded, &trajectory)?;
                Command::RestoreCoordinates {
                    id,
                    structure: Box::new(std::mem::replace(&mut loaded.structure, structure)),
                    trajectory: loaded.trajectory.replace(trajectory),
                    frame: std::mem::take(&mut loaded.frame),
                }
            }

            Command::RestoreCoordinates {
                id,
                structure,
                trajectory,
                frame,
            } => {
                let loaded = scene
                    .structures_mut()
                    .get_mut(id)
                    .ok_or(SceneError::NoSuchStructure(id))?;
                Command::RestoreCoordinates {
                    id,
                    structure: Box::new(std::mem::replace(&mut loaded.structure, *structure)),
                    trajectory: std::mem::replace(&mut loaded.trajectory, trajectory),
                    frame: std::mem::replace(&mut loaded.frame, frame),
                }
            }

            Command::CloseStructure { id } => {
                let removed = scene
                    .structures_mut()
                    .remove(id)
                    .ok_or(SceneError::NoSuchStructure(id))?;
                let prior_active = if scene.active_selection().map(|a| a.structure) == Some(id) {
                    scene.set_active(None)
                } else {
                    None
                };
                Command::RestoreStructure {
                    id,
                    structure: Box::new(removed),
                    prior_active,
                }
            }

            Command::RestoreStructure {
                id,
                structure,
                prior_active,
            } => {
                scene.structures_mut().restore(id, *structure);
                if prior_active.is_some() {
                    scene.set_active(prior_active);
                }
                Command::CloseStructure { id }
            }

            Command::SetStructureLabel { id, label } => {
                let loaded = scene
                    .structures_mut()
                    .get_mut(id)
                    .ok_or(SceneError::NoSuchStructure(id))?;
                let old = std::mem::replace(&mut loaded.label, label);
                Command::SetStructureLabel { id, label: old }
            }

            Command::ShowStructure { id, visible } => {
                let loaded = scene
                    .structures_mut()
                    .get_mut(id)
                    .ok_or(SceneError::NoSuchStructure(id))?;
                let old = std::mem::replace(&mut loaded.visible, visible);
                Command::ShowStructure { id, visible: old }
            }

            Command::SetRepresentation {
                id,
                rep,
                representation,
            } => {
                let r = rep_mut(scene, id, rep)?;
                let old = std::mem::replace(&mut r.representation, representation);
                Command::SetRepresentation {
                    id,
                    rep,
                    representation: old,
                }
            }

            Command::SetColoring { id, rep, coloring } => {
                let r = rep_mut(scene, id, rep)?;
                let old = std::mem::replace(&mut r.coloring, coloring);
                Command::SetColoring {
                    id,
                    rep,
                    coloring: old,
                }
            }

            Command::SetMaterial { id, rep, material } => {
                let r = rep_mut(scene, id, rep)?;
                let old = std::mem::replace(&mut r.material, material);
                Command::SetMaterial {
                    id,
                    rep,
                    material: old,
                }
            }

            Command::AddRep { id, index, rep } => {
                let loaded = scene
                    .structures_mut()
                    .get_mut(id)
                    .ok_or(SceneError::NoSuchStructure(id))?;
                let index = index.min(loaded.reps.len());
                let rep_id = rep.id;
                loaded.next_rep_id = loaded.next_rep_id.max(rep.id.0 + 1);
                loaded.reps.insert(index, rep);
                loaded.current_rep = index;
                Command::RemoveRep { id, rep: rep_id }
            }

            Command::RemoveRep { id, rep } => {
                let loaded = scene
                    .structures_mut()
                    .get_mut(id)
                    .ok_or(SceneError::NoSuchStructure(id))?;
                let index = loaded
                    .rep_index(rep)
                    .ok_or(SceneError::NoSuchRep(id, rep))?;
                if loaded.reps.len() == 1 {
                    return Err(SceneError::LastRep(id));
                }
                let removed = loaded.reps.remove(index);
                loaded.current_rep = loaded.current_rep.min(loaded.reps.len() - 1);
                Command::AddRep {
                    id,
                    index,
                    rep: removed,
                }
            }

            Command::SetRepSelection { id, rep, selection } => {
                vv_core::select::parse(&selection).map_err(|source| SceneError::Selection {
                    expr: selection.clone(),
                    source,
                })?;
                let r = rep_mut(scene, id, rep)?;
                let old = std::mem::replace(&mut r.selection, selection);
                Command::SetRepSelection {
                    id,
                    rep,
                    selection: old,
                }
            }

            Command::SetRepOption {
                id,
                rep,
                name,
                value,
            } => {
                let r = rep_mut(scene, id, rep)?;
                let old = match value {
                    Some(v) => r.options.insert(name.clone(), v),
                    None => r.options.remove(&name),
                };
                Command::SetRepOption {
                    id,
                    rep,
                    name,
                    value: old,
                }
            }
            Command::ShowRep { id, rep, visible } => {
                let r = rep_mut(scene, id, rep)?;
                let old = std::mem::replace(&mut r.visible, visible);
                Command::ShowRep {
                    id,
                    rep,
                    visible: old,
                }
            }

            Command::SetCurrentRep { id, rep } => {
                let loaded = scene
                    .structures_mut()
                    .get_mut(id)
                    .ok_or(SceneError::NoSuchStructure(id))?;
                let index = loaded
                    .rep_index(rep)
                    .ok_or(SceneError::NoSuchRep(id, rep))?;
                let old = loaded.reps[loaded.current_rep].id;
                loaded.current_rep = index;
                Command::SetCurrentRep { id, rep: old }
            }

            Command::Batch(commands) => {
                let mut inverses = Vec::with_capacity(commands.len());
                for command in commands {
                    match command.apply(scene) {
                        Ok(inverse) => inverses.push(inverse),
                        Err(e) => {
                            for inverse in inverses.into_iter().rev() {
                                let _ = inverse.apply(scene);
                            }
                            return Err(e);
                        }
                    }
                }
                inverses.reverse();
                Command::Batch(inverses)
            }

            Command::SetFrame { id, frame } => {
                let loaded = scene
                    .structures_mut()
                    .get_mut(id)
                    .ok_or(SceneError::NoSuchStructure(id))?;
                let count = loaded.structure.frame_count();
                if frame >= count {
                    return Err(SceneError::FrameOutOfRange { id, frame, count });
                }
                let old = std::mem::replace(&mut loaded.frame, frame);
                Command::SetFrame { id, frame: old }
            }

            Command::SetValues { id, name, channel } => {
                let loaded = scene
                    .structures_mut()
                    .get_mut(id)
                    .ok_or(SceneError::NoSuchStructure(id))?;
                let atoms = loaded.structure.atom_count();
                if let Some(c) = &channel {
                    if c.atoms() != atoms {
                        return Err(SceneError::ValuesLength {
                            name,
                            expected: atoms,
                            actual: c.atoms(),
                        });
                    }
                }
                let old = match channel {
                    Some(c) => loaded.values.insert(name.clone(), c),
                    None => Some(
                        loaded
                            .values
                            .remove(&name)
                            .ok_or_else(|| SceneError::NoSuchValues(name.clone()))?,
                    ),
                };
                Command::SetValues {
                    id,
                    name,
                    channel: old,
                }
            }

            Command::Select { id, mask } => {
                if !scene.structures_mut().contains(id) {
                    return Err(SceneError::NoSuchStructure(id));
                }
                let old = scene.set_active(Some(ActiveSelection {
                    structure: id,
                    mask,
                    expr: None,
                }));
                Command::RestoreSelection { active: old }
            }

            Command::SelectExpr { id, expr } => {
                let loaded = scene.structure(id).ok_or(SceneError::NoSuchStructure(id))?;
                let structure = &loaded.structure;
                let bits =
                    vv_core::select(&structure.topology, structure.frame(0).positions(), &expr)
                        .map_err(|source| SceneError::Selection {
                            expr: expr.clone(),
                            source,
                        })?;
                let old = scene.set_active(Some(ActiveSelection {
                    structure: id,
                    mask: Arc::new(bits),
                    expr: Some(expr),
                }));
                Command::RestoreSelection { active: old }
            }

            Command::ClearSelection => {
                let old = scene.set_active(None);
                Command::RestoreSelection { active: old }
            }

            Command::RestoreSelection { active } => {
                let old = scene.set_active(active);
                Command::RestoreSelection { active: old }
            }

            Command::SaveSelectionSet { name } => {
                let active = scene
                    .active_selection()
                    .cloned()
                    .ok_or(SceneError::NoActiveSelection)?;
                let set = SelectionSet {
                    name: name.clone(),
                    structure: active.structure,
                    mask: active.mask,
                    expr: active.expr,
                };
                replace_selection_set(scene, set)
            }

            Command::DeleteSelectionSet { name } => {
                let sets = scene.selection_sets_mut();
                let index = sets
                    .iter()
                    .position(|s| s.name == name)
                    .ok_or_else(|| SceneError::NoSuchSelectionSet(name.clone()))?;
                let removed = sets.remove(index);
                Command::RestoreSelectionSet { set: removed }
            }

            Command::RestoreSelectionSet { set } => replace_selection_set(scene, set),

            Command::RenameSelectionSet { old_name, new_name } => {
                if scene
                    .selection_sets()
                    .iter()
                    .any(|s| s.name == new_name && s.name != old_name)
                {
                    return Err(SceneError::SelectionSetNameTaken(new_name));
                }
                let sets = scene.selection_sets_mut();
                let index = sets
                    .iter()
                    .position(|s| s.name == old_name)
                    .ok_or_else(|| SceneError::NoSuchSelectionSet(old_name.clone()))?;
                sets[index].name = new_name.clone();
                Command::RenameSelectionSet {
                    old_name: new_name,
                    new_name: old_name,
                }
            }

            Command::SetLabel { id, atom, text } => {
                let loaded = scene
                    .structures_mut()
                    .get_mut(id)
                    .ok_or(SceneError::NoSuchStructure(id))?;
                let atoms = loaded.structure.atom_count();
                if atom as usize >= atoms {
                    return Err(SceneError::NoSuchAtom { id, atom, atoms });
                }
                let old = match text {
                    Some(t) => loaded.labels.insert(atom, t),
                    None => loaded.labels.remove(&atom),
                };
                Command::SetLabel {
                    id,
                    atom,
                    text: old,
                }
            }

            Command::SetMeasurement {
                id,
                measurement,
                shown,
            } => {
                let loaded = scene
                    .structures_mut()
                    .get_mut(id)
                    .ok_or(SceneError::NoSuchStructure(id))?;
                let atoms = loaded.structure.atom_count();
                if let Some(&atom) = measurement.atoms().iter().find(|&&a| a as usize >= atoms) {
                    return Err(SceneError::NoSuchAtom { id, atom, atoms });
                }
                let list = &mut loaded.measurements;
                let was = list.iter().position(|m| *m == measurement);
                match (was, shown) {
                    (None, true) => list.push(measurement.clone()),
                    (Some(index), false) => {
                        list.remove(index);
                    }
                    _ => {}
                }
                Command::SetMeasurement {
                    id,
                    measurement,
                    shown: was.is_some(),
                }
            }

            Command::SetCaption { caption } => replace_caption(scene, caption),

            Command::DeleteCaption { name } => {
                let list = scene.captions_mut();
                let index = list
                    .iter()
                    .position(|c| c.name == name)
                    .ok_or_else(|| SceneError::NoSuchCaption(name.clone()))?;
                let removed = list.remove(index);
                Command::RestoreCaption { caption: removed }
            }

            Command::RestoreCaption { caption } => replace_caption(scene, caption),
        };
        scene.touch();
        Ok(inverse)
    }
}

/// A loaded structure's display label: the file's own id, or the file
/// name (gz and one extension stripped) when the id is empty or the RCSB
/// assembly-file placeholder `XXXX`.
fn structure_label(topology_id: &str, path: &std::path::Path) -> String {
    if !topology_id.is_empty() && !topology_id.eq_ignore_ascii_case("XXXX") {
        return topology_id.to_string();
    }
    let name = path
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_default();
    name.trim_end_matches(".gz")
        .rsplit_once('.')
        .map_or(name.clone(), |(stem, _)| stem.to_string())
}

/// Inserts `set`, overwriting any existing set of the same name, and
/// returns the command that undoes exactly that (restore the old one, or
/// delete this one if there wasn't an old one).
fn replace_selection_set(scene: &mut Scene, set: SelectionSet) -> Command {
    let sets = scene.selection_sets_mut();
    match sets.iter().position(|s| s.name == set.name) {
        Some(index) => {
            let old = std::mem::replace(&mut sets[index], set);
            Command::RestoreSelectionSet { set: old }
        }
        None => {
            let name = set.name.clone();
            sets.push(set);
            Command::DeleteSelectionSet { name }
        }
    }
}

/// Inserts `caption`, overwriting any existing caption of the same name,
/// and returns the command that undoes exactly that.
fn replace_caption(scene: &mut Scene, caption: Caption) -> Command {
    let list = scene.captions_mut();
    match list.iter().position(|c| c.name == caption.name) {
        Some(index) => {
            let old = std::mem::replace(&mut list[index], caption);
            Command::RestoreCaption { caption: old }
        }
        None => {
            let name = caption.name.clone();
            list.push(caption);
            Command::DeleteCaption { name }
        }
    }
}

#[cfg(test)]
mod tests {

    #[test]
    fn a_batch_is_one_undo_step() {
        let mut scene = Scene::new();
        let mut history = CommandHistory::default();
        load(&mut scene, &mut history);
        load(&mut scene, &mut history);
        let ids: Vec<_> = scene.structures().map(|(id, _)| id).collect();
        let (a, b) = (ids[0], ids[1]);
        history
            .dispatch(
                &mut scene,
                Command::Batch(vec![
                    Command::SetMaterial {
                        id: a,
                        rep: RepId(0),
                        material: Material::Glossy,
                    },
                    Command::SetMaterial {
                        id: b,
                        rep: RepId(0),
                        material: Material::Textbook,
                    },
                ]),
            )
            .unwrap();
        assert_eq!(
            scene.structure(b).unwrap().rep().material,
            Material::Textbook
        );
        history.undo(&mut scene).unwrap();
        assert_eq!(scene.structure(a).unwrap().rep().material, Material::Opaque);
        assert_eq!(scene.structure(b).unwrap().rep().material, Material::Opaque);
        history.redo(&mut scene).unwrap();
        assert_eq!(scene.structure(a).unwrap().rep().material, Material::Glossy);
    }

    use super::*;
    use crate::history::CommandHistory;
    use crate::selection::mask_of_one;

    fn small_cif_path() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/small/1CRN.cif")
    }

    fn small_fixture(name: &str) -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../fixtures/small")
            .join(name)
    }

    fn load(scene: &mut Scene, history: &mut CommandHistory) -> StructureId {
        history
            .dispatch(
                scene,
                Command::LoadStructure {
                    path: small_cif_path(),
                },
            )
            .unwrap();
        scene.structures().next().unwrap().0
    }

    #[test]
    fn load_trajectory_produces_a_multi_frame_structure() {
        let mut scene = Scene::new();
        let mut history = CommandHistory::new(100);
        history
            .dispatch(
                &mut scene,
                Command::LoadTrajectory {
                    topology: small_fixture("1CRN.pdb"),
                    trajectory: small_fixture("1CRN_traj.dcd"),
                },
            )
            .unwrap();
        let (id, loaded) = scene.structures().next().unwrap();
        assert_eq!(loaded.label, "1CRN");
        assert_eq!(loaded.structure.atom_count(), 327);
        assert_eq!(loaded.structure.frame_count(), 3);
        let base = loaded.structure.frame(0).positions()[0];
        let last = loaded.structure.frame(2).positions()[0];
        assert!((last.x - base.x - 2.0).abs() < 1e-4);

        // Undo removes it exactly like closing a single-frame structure.
        assert!(history.undo(&mut scene).unwrap());
        assert_eq!(scene.structures().count(), 0);
        assert!(history.redo(&mut scene).unwrap());
        assert_eq!(scene.structure(id).unwrap().structure.frame_count(), 3);
    }

    #[test]
    fn xtc_trr_and_netcdf_trajectories_load_like_dcd() {
        let traj = |name: &str| {
            std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
                .join("../../fixtures/traj")
                .join(name)
        };
        for name in ["1crn.xtc", "1crn.trr", "1crn.nc"] {
            let mut scene = Scene::new();
            let mut history = CommandHistory::new(100);
            history
                .dispatch(
                    &mut scene,
                    Command::LoadTrajectory {
                        topology: small_fixture("1CRN.pdb"),
                        trajectory: traj(name),
                    },
                )
                .unwrap();
            let (_, loaded) = scene.structures().next().unwrap();
            assert_eq!(loaded.structure.frame_count(), 10, "{name}");
        }
    }

    #[test]
    fn load_trajectory_rejects_an_atom_count_mismatch() {
        let mut scene = Scene::new();
        let mut history = CommandHistory::new(100);
        let err = history
            .dispatch(
                &mut scene,
                Command::LoadTrajectory {
                    topology: small_cif_path(),
                    trajectory: small_fixture("sample.dcd"),
                },
            )
            .unwrap_err();
        assert!(matches!(
            err,
            SceneError::TrajectoryAtomMismatch {
                topology_atoms: 327,
                trajectory: 5,
                ..
            }
        ));
        assert_eq!(scene.structures().count(), 0);
    }

    #[test]
    fn attach_trajectory_keeps_the_structure_and_undoes() {
        let mut scene = Scene::new();
        let mut history = CommandHistory::new(100);
        let id = load(&mut scene, &mut history);
        history
            .dispatch(
                &mut scene,
                Command::SetLabel {
                    id,
                    atom: 3,
                    text: Some("here".into()),
                },
            )
            .unwrap();
        history
            .dispatch(
                &mut scene,
                Command::AttachTrajectory {
                    id,
                    trajectory: small_fixture("1CRN_traj.dcd"),
                },
            )
            .unwrap();
        assert_eq!(scene.structures().count(), 1);
        let loaded = scene.structure(id).unwrap();
        assert_eq!(loaded.structure.frame_count(), 3);
        assert_eq!(loaded.frame, 0);
        assert_eq!(loaded.labels.get(&3).map(String::as_str), Some("here"));
        assert_eq!(loaded.trajectory, Some(small_fixture("1CRN_traj.dcd")));

        assert!(history.undo(&mut scene).unwrap());
        let loaded = scene.structure(id).unwrap();
        assert_eq!(loaded.structure.frame_count(), 1);
        assert_eq!(loaded.trajectory, None);
        assert!(history.redo(&mut scene).unwrap());
        assert_eq!(scene.structure(id).unwrap().structure.frame_count(), 3);
    }

    #[test]
    fn attach_trajectory_rejects_an_atom_count_mismatch() {
        let mut scene = Scene::new();
        let mut history = CommandHistory::new(100);
        let id = load(&mut scene, &mut history);
        let err = history
            .dispatch(
                &mut scene,
                Command::AttachTrajectory {
                    id,
                    trajectory: small_fixture("sample.dcd"),
                },
            )
            .unwrap_err();
        assert!(matches!(
            err,
            SceneError::TrajectoryAtomMismatch {
                topology_atoms: 327,
                trajectory: 5,
                ..
            }
        ));
        assert_eq!(scene.structure(id).unwrap().structure.frame_count(), 1);
    }

    fn write_temp(name: &str, content: &str) -> PathBuf {
        let dir = std::env::temp_dir().join("vizviz-tests");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join(name);
        std::fs::write(&path, content).unwrap();
        path
    }

    /// 5 atoms, matching `fixtures/small/sample.dcd`'s atom count exactly.
    /// Its bonds have nothing to do with that fixture's own coordinates
    /// (atoms ~17 A apart) -- see the next test.
    const FIVE_ATOM_PSF: &str = "PSF\n\
\n\
       0 !NTITLE\n\
\n\
       5 !NATOM\n\
       1 A    1    LIG  C1   C      0.000000       12.0110           0\n\
       2 A    1    LIG  C2   C      0.000000       12.0110           0\n\
       3 A    1    LIG  C3   C      0.000000       12.0110           0\n\
       4 A    1    LIG  C4   C      0.000000       12.0110           0\n\
       5 A    1    LIG  C5   C      0.000000       12.0110           0\n\
\n\
       4 !NBOND: bonds\n\
       1       2       2       3       3       4       4       5\n";

    #[test]
    fn loadtrajectory_accepts_a_psf_topology_with_a_real_dcd_and_keeps_its_bonds_verbatim() {
        let psf = write_temp("md_bonds.psf", FIVE_ATOM_PSF);
        let mut scene = Scene::new();
        let mut history = CommandHistory::new(100);
        history
            .dispatch(
                &mut scene,
                Command::LoadTrajectory {
                    topology: psf,
                    trajectory: small_fixture("sample.dcd"),
                },
            )
            .unwrap();
        let (_, loaded) = scene.structures().next().unwrap();
        assert_eq!(loaded.structure.atom_count(), 5);
        assert_eq!(loaded.structure.frame_count(), 4);
        // sample.dcd's atoms sit ~17 A apart: geometric perception alone
        // would find nothing. The PSF's own chain of bonds must survive.
        assert_eq!(loaded.bonds.pairs, vec![[0, 1], [1, 2], [2, 3], [3, 4]]);
    }

    /// One fixed-column wwPDB `ATOM` record (columns as `vv_io::pdb`
    /// reads them); only the coordinates matter here, since a PDB used as
    /// `loadtrajectory`'s trajectory argument contributes positions only,
    /// never the topology.
    fn pdb_atom_line(serial: u32, x: f32) -> String {
        let mut line = vec![b' '; 80];
        line[0..6].copy_from_slice(b"ATOM  ");
        line[6..11].copy_from_slice(format!("{serial:>5}").as_bytes());
        line[12..16].copy_from_slice(b"C1  ");
        line[17..20].copy_from_slice(b"LIG");
        line[21] = b'A';
        line[22..26].copy_from_slice(b"   1");
        line[30..38].copy_from_slice(format!("{x:>8.3}").as_bytes());
        line[38..46].copy_from_slice(b"   0.000");
        line[46..54].copy_from_slice(b"   0.000");
        line[54..60].copy_from_slice(b"  1.00");
        line[60..66].copy_from_slice(b"  0.00");
        String::from_utf8(line).unwrap()
    }

    #[test]
    fn loadtrajectory_accepts_a_pdb_as_the_coordinate_source() {
        let psf = write_temp("with_pdb_coords.psf", FIVE_ATOM_PSF);
        let mut text = String::new();
        for (i, x) in [0.0, 1.5, 3.0, 4.5, 6.0].into_iter().enumerate() {
            text.push_str(&pdb_atom_line(i as u32 + 1, x));
            text.push('\n');
        }
        text.push_str("END\n");
        let pdb = write_temp("with_pdb_coords.pdb", &text);
        let mut scene = Scene::new();
        let mut history = CommandHistory::new(100);
        history
            .dispatch(
                &mut scene,
                Command::LoadTrajectory {
                    topology: psf,
                    trajectory: pdb,
                },
            )
            .unwrap();
        let (_, loaded) = scene.structures().next().unwrap();
        assert_eq!(loaded.structure.atom_count(), 5);
        assert_eq!(loaded.structure.frame_count(), 1);
        assert_eq!(loaded.structure.frame(0).positions()[1].x, 1.5);
        assert_eq!(loaded.bonds.pairs, vec![[0, 1], [1, 2], [2, 3], [3, 4]]);
    }

    #[test]
    fn loadtrajectory_reports_an_atom_count_mismatch_for_a_psf_topology() {
        let psf = write_temp(
            "mismatch.psf",
            "PSF\n\n       0 !NTITLE\n\n       3 !NATOM\n\
       1 A    1    LIG  C1   C      0.000000       12.0110           0\n\
       2 A    1    LIG  C2   C      0.000000       12.0110           0\n\
       3 A    1    LIG  C3   C      0.000000       12.0110           0\n\
\n       2 !NBOND: bonds\n       1       2       2       3\n",
        );
        let mut scene = Scene::new();
        let mut history = CommandHistory::new(100);
        let err = history
            .dispatch(
                &mut scene,
                Command::LoadTrajectory {
                    topology: psf,
                    trajectory: small_fixture("sample.dcd"),
                },
            )
            .unwrap_err();
        assert!(matches!(
            err,
            SceneError::TrajectoryAtomMismatch {
                topology_atoms: 3,
                trajectory: 5,
                ..
            }
        ));
    }

    #[test]
    fn load_close_undo_redo_round_trips() {
        let mut scene = Scene::new();
        let mut history = CommandHistory::new(100);
        let id = load(&mut scene, &mut history);
        assert_eq!(scene.structures().count(), 1);
        assert_eq!(scene.structure(id).unwrap().label, "1CRN");

        history
            .dispatch(&mut scene, Command::CloseStructure { id })
            .unwrap();
        assert_eq!(scene.structures().count(), 0);

        assert!(history.undo(&mut scene).unwrap());
        assert_eq!(scene.structures().count(), 1);
        assert_eq!(scene.structure(id).unwrap().label, "1CRN");

        assert!(history.redo(&mut scene).unwrap());
        assert_eq!(scene.structures().count(), 0);

        // `done` now holds both the load and the close; two more undos
        // (close, then load) exhaust it.
        assert!(history.undo(&mut scene).unwrap());
        assert_eq!(scene.structures().count(), 1);
        assert!(history.undo(&mut scene).unwrap());
        assert_eq!(scene.structures().count(), 0);
        assert!(!history.undo(&mut scene).unwrap(), "nothing left to undo");
    }

    #[test]
    fn renaming_a_structure_undoes_back_to_its_old_label() {
        let mut scene = Scene::new();
        let mut history = CommandHistory::new(100);
        let id = load(&mut scene, &mut history);
        history
            .dispatch(
                &mut scene,
                Command::SetStructureLabel {
                    id,
                    label: "My protein".into(),
                },
            )
            .unwrap();
        assert_eq!(scene.structure(id).unwrap().label, "My protein");
        assert!(history.undo(&mut scene).unwrap());
        assert_eq!(scene.structure(id).unwrap().label, "1CRN");
        assert!(history.redo(&mut scene).unwrap());
        assert_eq!(scene.structure(id).unwrap().label, "My protein");
    }

    #[test]
    fn selecting_something_new_clears_the_previous_redo_stack() {
        let mut scene = Scene::new();
        let mut history = CommandHistory::new(100);
        let id = load(&mut scene, &mut history);
        history.undo(&mut scene).unwrap();
        assert!(history.can_redo());
        let id2 = load(&mut scene, &mut history);
        assert!(
            !history.can_redo(),
            "a fresh command should drop the redo stack"
        );
        assert_eq!(scene.structures().count(), 1);
        // The reloaded structure gets a new id (the old slot stays a
        // tombstone), which is the point of never reusing indices.
        assert_ne!(id, id2);
    }

    #[test]
    fn select_and_clear_undo_restores_prior_active_selection() {
        let mut scene = Scene::new();
        let mut history = CommandHistory::new(100);
        let id = load(&mut scene, &mut history);
        let n = scene.structure(id).unwrap().structure.atom_count();

        assert!(scene.active_selection().is_none());
        history
            .dispatch(
                &mut scene,
                Command::Select {
                    id,
                    mask: mask_of_one(n, 0),
                },
            )
            .unwrap();
        assert_eq!(scene.active_selection().unwrap().mask.count_ones(..), 1);

        history
            .dispatch(
                &mut scene,
                Command::Select {
                    id,
                    mask: mask_of_one(n, 5),
                },
            )
            .unwrap();
        assert!(scene.active_selection().unwrap().mask[5]);

        history.undo(&mut scene).unwrap();
        assert!(
            scene.active_selection().unwrap().mask[0],
            "undo should restore the first pick"
        );

        history.undo(&mut scene).unwrap();
        assert!(
            scene.active_selection().is_none(),
            "undo should restore no-selection"
        );
    }

    #[test]
    fn select_by_expression_remembers_the_expression_and_undoes() {
        let mut scene = Scene::new();
        let mut history = CommandHistory::new(100);
        let id = load(&mut scene, &mut history);
        let residues = scene
            .structure(id)
            .unwrap()
            .structure
            .topology
            .residue_count();

        history
            .dispatch(
                &mut scene,
                Command::SelectExpr {
                    id,
                    expr: "name CA".into(),
                },
            )
            .unwrap();
        let active = scene.active_selection().unwrap();
        assert_eq!(active.mask.count_ones(..), residues, "one CA per residue");
        assert_eq!(active.expr.as_deref(), Some("name CA"));

        // The expression travels into a saved set.
        history
            .dispatch(&mut scene, Command::SaveSelectionSet { name: "ca".into() })
            .unwrap();
        assert_eq!(
            scene.selection_set("ca").unwrap().expr.as_deref(),
            Some("name CA")
        );

        // A pick has no expression; undoing it brings the expression back.
        history
            .dispatch(
                &mut scene,
                Command::Select {
                    id,
                    mask: mask_of_one(residues, 0),
                },
            )
            .unwrap();
        assert_eq!(scene.active_selection().unwrap().expr, None);
        history.undo(&mut scene).unwrap();
        assert_eq!(
            scene.active_selection().unwrap().expr.as_deref(),
            Some("name CA")
        );

        // A bad expression is an error and leaves the selection alone.
        let err = history
            .dispatch(
                &mut scene,
                Command::SelectExpr {
                    id,
                    expr: "chain A and".into(),
                },
            )
            .unwrap_err();
        assert!(matches!(err, SceneError::Selection { .. }), "{err}");
        assert_eq!(
            scene.active_selection().unwrap().mask.count_ones(..),
            residues
        );
    }

    #[test]
    fn save_selection_set_then_delete_then_undo_round_trips() {
        let mut scene = Scene::new();
        let mut history = CommandHistory::new(100);
        let id = load(&mut scene, &mut history);
        let n = scene.structure(id).unwrap().structure.atom_count();
        history
            .dispatch(
                &mut scene,
                Command::Select {
                    id,
                    mask: mask_of_one(n, 3),
                },
            )
            .unwrap();
        history
            .dispatch(
                &mut scene,
                Command::SaveSelectionSet {
                    name: "pocket".into(),
                },
            )
            .unwrap();
        assert_eq!(
            scene.selection_set("pocket").unwrap().mask.count_ones(..),
            1
        );

        // Overwriting an existing set must be undoable back to the old mask.
        history
            .dispatch(
                &mut scene,
                Command::Select {
                    id,
                    mask: mask_of_one(n, 7),
                },
            )
            .unwrap();
        history
            .dispatch(
                &mut scene,
                Command::SaveSelectionSet {
                    name: "pocket".into(),
                },
            )
            .unwrap();
        assert!(scene.selection_set("pocket").unwrap().mask[7]);
        history.undo(&mut scene).unwrap();
        assert!(
            scene.selection_set("pocket").unwrap().mask[3],
            "undo should restore old mask"
        );

        history
            .dispatch(
                &mut scene,
                Command::DeleteSelectionSet {
                    name: "pocket".into(),
                },
            )
            .unwrap();
        assert!(scene.selection_set("pocket").is_none());
        history.undo(&mut scene).unwrap();
        assert!(
            scene.selection_set("pocket").is_some(),
            "undo should restore the deleted set"
        );
    }

    #[test]
    fn renaming_a_selection_set_undoes_and_refuses_a_collision() {
        let mut scene = Scene::new();
        let mut history = CommandHistory::new(100);
        let id = load(&mut scene, &mut history);
        let n = scene.structure(id).unwrap().structure.atom_count();
        for name in ["pocket", "loop"] {
            history
                .dispatch(
                    &mut scene,
                    Command::Select {
                        id,
                        mask: mask_of_one(n, 3),
                    },
                )
                .unwrap();
            history
                .dispatch(&mut scene, Command::SaveSelectionSet { name: name.into() })
                .unwrap();
        }
        history
            .dispatch(
                &mut scene,
                Command::RenameSelectionSet {
                    old_name: "pocket".into(),
                    new_name: "active site".into(),
                },
            )
            .unwrap();
        assert!(scene.selection_set("pocket").is_none());
        assert!(scene.selection_set("active site").is_some());
        assert!(history.undo(&mut scene).unwrap());
        assert!(scene.selection_set("pocket").is_some());
        assert!(scene.selection_set("active site").is_none());

        let err = history
            .dispatch(
                &mut scene,
                Command::RenameSelectionSet {
                    old_name: "pocket".into(),
                    new_name: "loop".into(),
                },
            )
            .unwrap_err();
        assert!(matches!(err, SceneError::SelectionSetNameTaken(n) if n == "loop"));

        // Renaming a set to its own name is not a collision.
        history
            .dispatch(
                &mut scene,
                Command::RenameSelectionSet {
                    old_name: "pocket".into(),
                    new_name: "pocket".into(),
                },
            )
            .unwrap();
    }

    #[test]
    fn errors_do_not_corrupt_history() {
        let mut scene = Scene::new();
        let mut history = CommandHistory::new(100);
        let bogus = crate::slotmap::Id::from_raw(999);
        let err = history
            .dispatch(&mut scene, Command::CloseStructure { id: bogus })
            .unwrap_err();
        assert!(matches!(err, SceneError::NoSuchStructure(_)));
        assert!(
            !history.can_undo(),
            "a failed command must not enter history"
        );

        let err = history
            .dispatch(&mut scene, Command::SaveSelectionSet { name: "x".into() })
            .unwrap_err();
        assert!(matches!(err, SceneError::NoActiveSelection));
    }

    #[test]
    fn set_representation_round_trips() {
        let mut scene = Scene::new();
        let mut history = CommandHistory::new(100);
        let id = load(&mut scene, &mut history);
        assert_eq!(
            scene.structure(id).unwrap().rep().representation,
            Representation::Lines
        );
        history
            .dispatch(
                &mut scene,
                Command::SetRepresentation {
                    id,
                    rep: RepId(0),
                    representation: Representation::BallAndStick,
                },
            )
            .unwrap();
        assert_eq!(
            scene.structure(id).unwrap().rep().representation,
            Representation::BallAndStick
        );
        history.undo(&mut scene).unwrap();
        assert_eq!(
            scene.structure(id).unwrap().rep().representation,
            Representation::Lines
        );
    }

    #[test]
    fn show_structure_round_trips_independent_of_show_rep() {
        let mut scene = Scene::new();
        let mut history = CommandHistory::new(100);
        let id = load(&mut scene, &mut history);
        assert!(scene.structure(id).unwrap().visible);
        history
            .dispatch(&mut scene, Command::ShowStructure { id, visible: false })
            .unwrap();
        assert!(!scene.structure(id).unwrap().visible);
        // A rep's own visibility is untouched by the structure-level flag.
        assert!(scene.structure(id).unwrap().rep().visible);
        history.undo(&mut scene).unwrap();
        assert!(scene.structure(id).unwrap().visible);
    }

    #[test]
    fn set_frame_round_trips_validates_and_is_undoable() {
        let mut scene = Scene::new();
        let mut history = CommandHistory::new(100);
        history
            .dispatch(
                &mut scene,
                Command::LoadTrajectory {
                    topology: small_fixture("1CRN.pdb"),
                    trajectory: small_fixture("1CRN_traj.dcd"),
                },
            )
            .unwrap();
        let id = scene.structures().next().unwrap().0;
        assert_eq!(scene.structure(id).unwrap().frame, 0);

        history
            .dispatch(&mut scene, Command::SetFrame { id, frame: 2 })
            .unwrap();
        assert_eq!(scene.structure(id).unwrap().frame, 2);
        history.undo(&mut scene).unwrap();
        assert_eq!(scene.structure(id).unwrap().frame, 0);
        history.redo(&mut scene).unwrap();
        assert_eq!(scene.structure(id).unwrap().frame, 2);

        // Out of range (the fixture has 3 frames, 0..=2): an error, and
        // the frame is left unchanged -- a rejected command must not
        // enter history either.
        let err = history
            .dispatch(&mut scene, Command::SetFrame { id, frame: 3 })
            .unwrap_err();
        assert!(matches!(
            err,
            SceneError::FrameOutOfRange {
                frame: 3,
                count: 3,
                ..
            }
        ));
        assert_eq!(scene.structure(id).unwrap().frame, 2);
    }

    #[test]
    fn set_frame_live_bypasses_history_and_clamps() {
        let mut scene = Scene::new();
        let mut history = CommandHistory::new(100);
        history
            .dispatch(
                &mut scene,
                Command::LoadTrajectory {
                    topology: small_fixture("1CRN.pdb"),
                    trajectory: small_fixture("1CRN_traj.dcd"),
                },
            )
            .unwrap();
        let id = scene.structures().next().unwrap().0;
        assert!(history.can_undo(), "the load itself is undoable");

        scene.set_frame_live(id, 1);
        assert_eq!(scene.structure(id).unwrap().frame, 1);
        // No new history entry was created: undoing still reverts the
        // load (removing the structure), not just a frame change.
        assert!(history.undo(&mut scene).unwrap());
        assert!(scene.structure(id).is_none());
        history.redo(&mut scene).unwrap();
        let id = scene.structures().next().unwrap().0;

        // Clamped to the last valid frame, not an error or a panic --
        // playback advancing past the end should never crash.
        scene.set_frame_live(id, 999);
        assert_eq!(scene.structure(id).unwrap().frame, 2);

        // A structure that doesn't exist is a silent no-op.
        scene.set_frame_live(StructureId::from_raw(9999), 1);
    }

    #[test]
    fn set_coloring_round_trips() {
        let mut scene = Scene::new();
        let mut history = CommandHistory::new(100);
        let id = load(&mut scene, &mut history);
        assert_eq!(
            scene.structure(id).unwrap().rep().coloring,
            ColorScheme::Element
        );
        history
            .dispatch(
                &mut scene,
                Command::SetColoring {
                    id,
                    rep: RepId(0),
                    coloring: ColorScheme::BFactor,
                },
            )
            .unwrap();
        assert_eq!(
            scene.structure(id).unwrap().rep().coloring,
            ColorScheme::BFactor
        );
        history.undo(&mut scene).unwrap();
        assert_eq!(
            scene.structure(id).unwrap().rep().coloring,
            ColorScheme::Element
        );
    }

    #[test]
    fn set_label_then_clear_then_undo_round_trips() {
        let mut scene = Scene::new();
        let mut history = CommandHistory::new(100);
        let id = load(&mut scene, &mut history);
        history
            .dispatch(
                &mut scene,
                Command::SetLabel {
                    id,
                    atom: 3,
                    text: Some("Thr1 N".into()),
                },
            )
            .unwrap();
        assert_eq!(
            scene
                .structure(id)
                .unwrap()
                .labels
                .get(&3)
                .map(String::as_str),
            Some("Thr1 N")
        );

        // Overwriting an existing label must be undoable back to the old text.
        history
            .dispatch(
                &mut scene,
                Command::SetLabel {
                    id,
                    atom: 3,
                    text: Some("renamed".into()),
                },
            )
            .unwrap();
        history.undo(&mut scene).unwrap();
        assert_eq!(
            scene
                .structure(id)
                .unwrap()
                .labels
                .get(&3)
                .map(String::as_str),
            Some("Thr1 N")
        );

        // Clearing (text: None) removes it; undo restores it.
        history
            .dispatch(
                &mut scene,
                Command::SetLabel {
                    id,
                    atom: 3,
                    text: None,
                },
            )
            .unwrap();
        assert!(scene.structure(id).unwrap().labels.is_empty());
        history.undo(&mut scene).unwrap();
        assert_eq!(
            scene
                .structure(id)
                .unwrap()
                .labels
                .get(&3)
                .map(String::as_str),
            Some("Thr1 N")
        );
    }

    #[test]
    fn set_label_rejects_an_out_of_range_atom() {
        let mut scene = Scene::new();
        let mut history = CommandHistory::new(100);
        let id = load(&mut scene, &mut history);
        let atoms = scene.structure(id).unwrap().structure.atom_count() as u32;
        let err = history
            .dispatch(
                &mut scene,
                Command::SetLabel {
                    id,
                    atom: atoms,
                    text: Some("out of range".into()),
                },
            )
            .unwrap_err();
        assert!(matches!(err, SceneError::NoSuchAtom { atom, .. } if atom == atoms));
        // A rejected command must not have entered the history: the load
        // that set up this test is still the only undoable step.
        assert!(history.can_undo());
        assert!(history.undo(&mut scene).unwrap());
        assert!(!history.can_undo());
    }

    #[test]
    fn caption_save_delete_undo_round_trips() {
        let mut scene = Scene::new();
        let mut history = CommandHistory::new(100);
        history
            .dispatch(
                &mut scene,
                Command::SetCaption {
                    caption: Caption {
                        name: "title".into(),
                        x: 0.5,
                        y: 0.05,
                        text: "Crambin, 1CRN".into(),
                    },
                },
            )
            .unwrap();
        assert_eq!(scene.captions().len(), 1);
        assert_eq!(scene.captions()[0].text, "Crambin, 1CRN");

        // Overwriting an existing caption (same name) must be undoable.
        history
            .dispatch(
                &mut scene,
                Command::SetCaption {
                    caption: Caption {
                        name: "title".into(),
                        x: 0.5,
                        y: 0.05,
                        text: "renamed".into(),
                    },
                },
            )
            .unwrap();
        history.undo(&mut scene).unwrap();
        assert_eq!(scene.captions()[0].text, "Crambin, 1CRN");

        history
            .dispatch(
                &mut scene,
                Command::DeleteCaption {
                    name: "title".into(),
                },
            )
            .unwrap();
        assert!(scene.captions().is_empty());
        history.undo(&mut scene).unwrap();
        assert_eq!(scene.captions()[0].text, "Crambin, 1CRN");
    }
}
