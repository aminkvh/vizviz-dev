//! Session files: everything needed to get back to where you were,
//! as JSON. A session stores *references* (file paths, selection
//! expressions, names), never the structure's own atom data: a 10M-atom
//! structure is a one-line path, and a saved selection is the expression
//! that made it, re-evaluated on load. A picked selection with no
//! expression falls back to its atom indices. The one exception is value
//! channels (`vv_scene::values`): numbers computed outside vizviz have no
//! path to re-run, so they are written as `.npy` sidecars next to the
//! session file and referenced by path, same as a structure.
//!
//! The document half (structures, sets, active selection, each
//! structure's current frame) lives here. The app adds its own `view`
//! block (camera, style) as opaque JSON; `vv-scene` carries it without
//! looking inside, so a session saved from the app round-trips through a
//! headless Python `Session` intact.
//!
//! A structure's current frame moved from the app's own opaque `view`
//! block to here (`SavedStructure::frame`) once `Command::SetFrame` made
//! it real document state — see `vv_scene::scene::LoadedStructure::frame`'s
//! doc. A session saved by the version of this file that kept frame in
//! `view` still loads (the key is simply absent from `structures[i]`, so
//! `#[serde(default)]` gives `0`), but its own recorded frame is not
//! restored: carrying that old key's restoration forward would mean two
//! places a loaded session's frame could come from, exactly the
//! confusion moving it here was meant to end. A session saved by *this*
//! version always round-trips its frame correctly.

use std::path::{Path, PathBuf};
use std::sync::Arc;

use fixedbitset::FixedBitSet;
use serde::{Deserialize, Serialize};

use crate::command::Command;
use crate::history::CommandHistory;
use crate::scene::{
    Caption, ColorScheme, Material, Measurement, Rep, RepId, Representation, Scene, StructureId,
};
use crate::values::ValueChannel;

/// 2: structures carry `reps` (several representations each; a schema-1
/// file's single representation/coloring/material still loads, as one
/// rep). 3: rep options, atom labels and captions. 4: the app's `view`
/// block also carries the dock layout at save time and playback state --
/// still opaque here, so a schema-3 file simply lacks them, and loading
/// one leaves both as they already are.
pub const SCHEMA_VERSION: u32 = 4;

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct SessionFile {
    pub schema_version: u32,
    /// Structures in scene order; selections refer to them by index.
    pub structures: Vec<SavedStructure>,
    #[serde(default)]
    pub selection_sets: Vec<SavedSet>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub active: Option<SavedSelection>,
    /// The app's view state (camera, style, frames). Opaque here.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub view: Option<serde_json::Value>,
    /// Screen captions; absent before schema 3.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub captions: Vec<SavedCaption>,
}

/// `scene::Caption` as saved: position is a fraction of the viewport.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct SavedCaption {
    pub name: String,
    pub x: f32,
    pub y: f32,
    pub text: String,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct SavedStructure {
    pub path: PathBuf,
    /// The trajectory its frames came from, loaded with it again.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub trajectory: Option<PathBuf>,
    pub label: String,
    /// Whole-structure show/hide (`LoadedStructure::visible`); absent
    /// before this field existed, so `#[serde(default = "shown")]`.
    #[serde(default = "shown")]
    pub visible: bool,
    /// How the structure is drawn; empty in a schema-1 file, which has
    /// the single `representation`/`coloring`/`material` below instead.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub reps: Vec<SavedRep>,
    #[serde(default)]
    pub current_rep: usize,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub representation: String,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub coloring: String,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub material: String,
    /// Value channels attached to this structure, each as a `.npy`
    /// sidecar path (see the module docs).
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub values: Vec<SavedValues>,
    /// Measurements on screen, as their atoms.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub measurements: Vec<Vec<u32>>,
    /// `LoadedStructure::frame` — absent (so `0`) in a file saved before
    /// frame became document state; see the module doc.
    #[serde(default)]
    pub frame: usize,
    /// Atom labels, by atom; absent before schema 3.
    #[serde(default, skip_serializing_if = "std::collections::BTreeMap::is_empty")]
    pub labels: std::collections::BTreeMap<u32, String>,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct SavedRep {
    pub selection: String,
    pub representation: String,
    pub coloring: String,
    pub material: String,
    #[serde(default = "shown")]
    pub visible: bool,
    /// `Rep::options`; absent before schema 3.
    #[serde(default, skip_serializing_if = "std::collections::BTreeMap::is_empty")]
    pub options: std::collections::BTreeMap<String, f32>,
}

fn shown() -> bool {
    true
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct SavedValues {
    pub name: String,
    pub path: PathBuf,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct SavedSelection {
    /// Index into `SessionFile::structures`.
    pub structure: usize,
    /// Re-evaluated on load when present (docs/SELECTION.md).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub expr: Option<String>,
    /// Explicit atom indices, only when there is no expression.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub atoms: Option<Vec<u32>>,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct SavedSet {
    pub name: String,
    #[serde(flatten)]
    pub selection: SavedSelection,
}

#[derive(Debug, thiserror::Error)]
pub enum SessionError {
    #[error("cannot {action} {path}: {source}")]
    Io {
        action: &'static str,
        path: PathBuf,
        source: std::io::Error,
    },
    #[error("{path} is not a vizviz session file: {source}")]
    Format {
        path: PathBuf,
        source: serde_json::Error,
    },
    #[error("{path} uses session schema {found}; this build reads schema {supported}")]
    Schema {
        path: PathBuf,
        found: u32,
        supported: u32,
    },
    #[error("structure #{index} in the session was never loaded (no path)")]
    Unsaved { index: usize },
}

pub fn representation_name(r: Representation) -> &'static str {
    match r {
        Representation::Spacefill => "spacefill",
        Representation::BallAndStick => "ball_and_stick",
        Representation::Tube => "tube",
        Representation::Cartoon => "cartoon",
        Representation::GaussianSurface => "gaussian_surface",
        Representation::SkinSurface => "skin_surface",
        Representation::Ses => "ses",
        Representation::Sas => "sas",
        Representation::Sticks => "sticks",
        Representation::Lines => "lines",
        Representation::Glycan => "glycan",
    }
}

pub fn coloring_name(c: &ColorScheme) -> String {
    c.name()
}

fn parse_representation(name: &str) -> Option<Representation> {
    match name {
        "spacefill" => Some(Representation::Spacefill),
        "ball_and_stick" => Some(Representation::BallAndStick),
        "tube" => Some(Representation::Tube),
        "cartoon" => Some(Representation::Cartoon),
        "gaussian_surface" => Some(Representation::GaussianSurface),
        "skin_surface" => Some(Representation::SkinSurface),
        "ses" => Some(Representation::Ses),
        "sas" => Some(Representation::Sas),
        // "licorice": a session saved before the Sticks rename.
        "sticks" | "licorice" => Some(Representation::Sticks),
        "lines" => Some(Representation::Lines),
        "glycan" => Some(Representation::Glycan),
        _ => None,
    }
}

/// A saved rep as a `Rep` for structure `id`, or `None` (with a warning)
/// when a name in it is unknown.
fn restore_rep(
    scene: &Scene,
    id: StructureId,
    index: usize,
    saved: &SavedRep,
    warnings: &mut Vec<String>,
) -> Option<Rep> {
    let Some(representation) = parse_representation(&saved.representation) else {
        warnings.push(format!(
            "rep {index}: unknown representation `{}`",
            saved.representation
        ));
        return None;
    };
    let Some(coloring) = parse_coloring(&saved.coloring) else {
        warnings.push(format!(
            "rep {index}: unknown coloring `{}`",
            saved.coloring
        ));
        return None;
    };
    let Some(material) = Material::parse(&saved.material) else {
        warnings.push(format!(
            "rep {index}: unknown material `{}`",
            saved.material
        ));
        return None;
    };
    let next = scene.structure(id)?.next_rep_id;
    Some(Rep {
        id: RepId(next),
        selection: saved.selection.clone(),
        representation,
        coloring,
        material,
        visible: saved.visible,
        options: saved.options.clone(),
    })
}

fn parse_coloring(name: &str) -> Option<ColorScheme> {
    ColorScheme::parse(name)
}

/// Describes `scene` as a session, with `view` attached verbatim. Paths
/// are made absolute so the file works from any working directory. Value
/// channels are written as `.npy` files under `values_dir` (created if
/// there are any).
pub fn capture(
    scene: &Scene,
    view: Option<serde_json::Value>,
    values_dir: &Path,
) -> Result<SessionFile, SessionError> {
    let values_dir = std::path::absolute(values_dir).unwrap_or_else(|_| values_dir.to_owned());
    if scene.structures().any(|(_, l)| !l.values.is_empty()) {
        std::fs::create_dir_all(&values_dir).map_err(|source| SessionError::Io {
            action: "create",
            path: values_dir.clone(),
            source,
        })?;
    }
    let ids: Vec<StructureId> = scene.structures().map(|(id, _)| id).collect();
    let index_of = |id: StructureId| ids.iter().position(|i| *i == id);
    let mut structures = Vec::new();
    for (index, (_, loaded)) in scene.structures().enumerate() {
        let path = loaded.path.clone().ok_or(SessionError::Unsaved { index })?;
        let path = std::path::absolute(&path).unwrap_or(path);
        let mut values = Vec::new();
        for (chan_index, (name, channel)) in loaded.values.iter().enumerate() {
            let chan_path = values_dir.join(format!("{index}_{chan_index}.npy"));
            crate::values::write_values(&chan_path, channel).map_err(|source| {
                SessionError::Io {
                    action: "write",
                    path: chan_path.clone(),
                    source,
                }
            })?;
            values.push(SavedValues {
                name: name.clone(),
                path: chan_path,
            });
        }
        structures.push(SavedStructure {
            path,
            trajectory: loaded
                .trajectory
                .as_ref()
                .map(|t| std::path::absolute(t).unwrap_or_else(|_| t.clone())),
            label: loaded.label.clone(),
            visible: loaded.visible,
            reps: loaded
                .reps
                .iter()
                .map(|r| SavedRep {
                    selection: r.selection.clone(),
                    representation: representation_name(r.representation).into(),
                    coloring: coloring_name(&r.coloring),
                    material: r.material.name().into(),
                    visible: r.visible,
                    options: r.options.clone(),
                })
                .collect(),
            current_rep: loaded.current_rep,
            representation: String::new(),
            coloring: String::new(),
            material: String::new(),
            values,
            measurements: loaded
                .measurements
                .iter()
                .map(|m| m.atoms().to_vec())
                .collect(),
            frame: loaded.frame,
            labels: loaded.labels.clone(),
        });
    }
    let saved = |structure: StructureId, expr: &Option<String>, mask: &FixedBitSet| {
        index_of(structure).map(|structure| SavedSelection {
            structure,
            expr: expr.clone(),
            atoms: expr
                .is_none()
                .then(|| mask.ones().map(|i| i as u32).collect()),
        })
    };
    let selection_sets = scene
        .selection_sets()
        .iter()
        .filter_map(|s| {
            saved(s.structure, &s.expr, &s.mask).map(|selection| SavedSet {
                name: s.name.clone(),
                selection,
            })
        })
        .collect();
    let active = scene
        .active_selection()
        .and_then(|a| saved(a.structure, &a.expr, &a.mask));
    let captions = scene
        .captions()
        .iter()
        .map(|c| SavedCaption {
            name: c.name.clone(),
            x: c.x,
            y: c.y,
            text: c.text.clone(),
        })
        .collect();
    Ok(SessionFile {
        schema_version: SCHEMA_VERSION,
        structures,
        selection_sets,
        active,
        view,
        captions,
    })
}

pub fn save(
    scene: &Scene,
    path: &Path,
    view: Option<serde_json::Value>,
) -> Result<(), SessionError> {
    let file = capture(scene, view, &values_sidecar_dir(path))?;
    let json = serde_json::to_string_pretty(&file).expect("session serializes");
    std::fs::write(path, json).map_err(|source| SessionError::Io {
        action: "write",
        path: path.to_owned(),
        source,
    })
}

/// Where `save` writes value-channel sidecars for a session at `path`:
/// `work.vviz` -> `work.values/`, next to it.
fn values_sidecar_dir(path: &Path) -> PathBuf {
    let mut name = path.file_stem().unwrap_or_default().to_os_string();
    name.push(".values");
    path.with_file_name(name)
}

pub fn read(path: &Path) -> Result<SessionFile, SessionError> {
    let text = std::fs::read_to_string(path).map_err(|source| SessionError::Io {
        action: "read",
        path: path.to_owned(),
        source,
    })?;
    let file: SessionFile = serde_json::from_str(&text).map_err(|source| SessionError::Format {
        path: path.to_owned(),
        source,
    })?;
    if file.schema_version > SCHEMA_VERSION {
        return Err(SessionError::Schema {
            path: path.to_owned(),
            found: file.schema_version,
            supported: SCHEMA_VERSION,
        });
    }
    Ok(file)
}

/// What `apply` produced: the id each saved structure got (`None` if its
/// file could not be loaded) and human-readable warnings for anything
/// skipped. A session with a missing file still loads the rest.
#[derive(Debug, Default)]
pub struct Applied {
    pub ids: Vec<Option<StructureId>>,
    pub warnings: Vec<String>,
}

/// Reads a value-channel sidecar back for the just-loaded structure `id`,
/// shaped to its current atom and frame count.
fn load_saved_values(path: &Path, scene: &Scene, id: StructureId) -> Result<ValueChannel, String> {
    let loaded = scene
        .structure(id)
        .ok_or_else(|| "structure is not loaded".to_string())?;
    let (atoms, frames) = (
        loaded.structure.atom_count(),
        loaded.structure.frame_count(),
    );
    let (data, shape) = crate::values::read_values(path).map_err(|e| e.to_string())?;
    ValueChannel::new(data, &shape, atoms, frames).map_err(|e| e.to_string())
}

/// Replays `file` into `scene` through the command bus: closes what is
/// open, loads each structure, restores representations, colorings,
/// selection sets, and the active selection. Every step is an ordinary
/// command, so a loaded session is undoable step by step.
pub fn apply(file: &SessionFile, scene: &mut Scene, history: &mut CommandHistory) -> Applied {
    let mut out = Applied::default();
    let open: Vec<StructureId> = scene.structures().map(|(id, _)| id).collect();
    for id in open {
        if let Err(e) = history.dispatch(scene, Command::CloseStructure { id }) {
            out.warnings.push(e.to_string());
        }
    }
    for (index, s) in file.structures.iter().enumerate() {
        let load = match &s.trajectory {
            Some(trajectory) => Command::LoadTrajectory {
                topology: s.path.clone(),
                trajectory: trajectory.clone(),
            },
            None => Command::LoadStructure {
                path: s.path.clone(),
            },
        };
        let id = match history.dispatch(scene, load) {
            Ok(()) => scene.structures().next_back().map(|(id, _)| id),
            Err(e) => {
                out.warnings.push(format!("structure #{index}: {e}"));
                None
            }
        };
        if let Some(id) = id {
            if !s.visible {
                if let Err(e) =
                    history.dispatch(scene, Command::ShowStructure { id, visible: false })
                {
                    out.warnings.push(format!("{}: {e}", s.label));
                }
            }
            let saved_reps = if s.reps.is_empty() {
                // Schema 1: one rep, from the structure's own fields.
                vec![SavedRep {
                    selection: "all".into(),
                    representation: s.representation.clone(),
                    coloring: s.coloring.clone(),
                    material: if s.material.is_empty() {
                        Material::default().name().into()
                    } else {
                        s.material.clone()
                    },
                    visible: true,
                    options: Default::default(),
                }]
            } else {
                s.reps.clone()
            };
            let mut attempted = std::collections::HashSet::new();
            for sv in &s.values {
                attempted.insert(sv.name.as_str());
                match load_saved_values(&sv.path, scene, id) {
                    Ok(channel) => {
                        let _ = history.dispatch(
                            scene,
                            Command::SetValues {
                                id,
                                name: sv.name.clone(),
                                channel: Some(channel),
                            },
                        );
                    }
                    Err(e) => out
                        .warnings
                        .push(format!("{}: value channel `{}`: {e}", s.label, sv.name)),
                }
            }
            for (i, saved) in saved_reps.iter().enumerate() {
                let rep = match restore_rep(scene, id, i, saved, &mut out.warnings) {
                    Some(rep) => rep,
                    None => continue,
                };
                if let ColorScheme::Values(name) = &rep.coloring {
                    let restored = scene
                        .structure(id)
                        .is_some_and(|l| l.values.contains_key(name));
                    if !restored && !attempted.contains(name.as_str()) {
                        out.warnings.push(format!(
                            "{}: value channel `{name}` is not available; run `values {name} PATH` to see it again",
                            s.label
                        ));
                    }
                }
                if i == 0 {
                    // The load made rep 0 already; bring it to the saved one.
                    let first = scene.structure(id).expect("loaded").reps[0].id;
                    let mut edits = vec![
                        Command::SetRepresentation {
                            id,
                            rep: first,
                            representation: rep.representation,
                        },
                        Command::SetColoring {
                            id,
                            rep: first,
                            coloring: rep.coloring.clone(),
                        },
                        Command::SetMaterial {
                            id,
                            rep: first,
                            material: rep.material,
                        },
                        Command::SetRepSelection {
                            id,
                            rep: first,
                            selection: rep.selection.clone(),
                        },
                        Command::ShowRep {
                            id,
                            rep: first,
                            visible: rep.visible,
                        },
                    ];
                    edits.extend(
                        rep.options
                            .iter()
                            .map(|(name, &value)| Command::SetRepOption {
                                id,
                                rep: first,
                                name: name.clone(),
                                value: Some(value),
                            }),
                    );
                    if let Err(e) = history.dispatch(scene, Command::Batch(edits)) {
                        out.warnings.push(format!("{}: rep 0: {e}", s.label));
                    }
                } else if let Err(e) =
                    history.dispatch(scene, Command::AddRep { id, index: i, rep })
                {
                    out.warnings.push(format!("{}: rep {i}: {e}", s.label));
                }
            }
            let current = scene
                .structure(id)
                .and_then(|l| l.reps.get(s.current_rep))
                .map(|r| r.id);
            if let Some(rep) = current {
                let _ = history.dispatch(scene, Command::SetCurrentRep { id, rep });
            }
            for atoms in &s.measurements {
                let shown = Measurement::new(atoms.clone())
                    .ok_or_else(|| format!("{atoms:?} is not 2 to 4 atoms"))
                    .and_then(|measurement| {
                        history
                            .dispatch(
                                scene,
                                Command::SetMeasurement {
                                    id,
                                    measurement,
                                    shown: true,
                                },
                            )
                            .map_err(|e| e.to_string())
                    });
                if let Err(e) = shown {
                    out.warnings.push(format!("{}: measurement: {e}", s.label));
                }
            }
            for (&atom, text) in &s.labels {
                let label = Command::SetLabel {
                    id,
                    atom,
                    text: Some(text.clone()),
                };
                if let Err(e) = history.dispatch(scene, label) {
                    out.warnings.push(format!("{}: label: {e}", s.label));
                }
            }
            if s.frame != 0 {
                // Out of range when the trajectory changed on disk or did
                // not load: warn rather than quietly show frame 0.
                if let Err(e) = history.dispatch(scene, Command::SetFrame { id, frame: s.frame }) {
                    out.warnings.push(format!("{}: {e}", s.label));
                }
            }
        }
        out.ids.push(id);
    }

    for c in &file.captions {
        let caption = Caption {
            name: c.name.clone(),
            x: c.x,
            y: c.y,
            text: c.text.clone(),
        };
        if let Err(e) = history.dispatch(scene, Command::SetCaption { caption }) {
            out.warnings.push(format!("caption {}: {e}", c.name));
        }
    }

    let select = |scene: &mut Scene, history: &mut CommandHistory, sel: &SavedSelection| {
        let id = out
            .ids
            .get(sel.structure)
            .copied()
            .flatten()
            .ok_or_else(|| format!("structure #{} is not loaded", sel.structure))?;
        let command = match (&sel.expr, &sel.atoms) {
            (Some(expr), _) => Command::SelectExpr {
                id,
                expr: expr.clone(),
            },
            (None, Some(atoms)) => {
                let n = scene.structure(id).map_or(0, |s| s.structure.atom_count());
                let mut bits = FixedBitSet::with_capacity(n);
                for &a in atoms {
                    if (a as usize) < n {
                        bits.insert(a as usize);
                    }
                }
                Command::Select {
                    id,
                    mask: Arc::new(bits),
                }
            }
            (None, None) => return Err("selection has neither expression nor atoms".into()),
        };
        history.dispatch(scene, command).map_err(|e| e.to_string())
    };
    for set in &file.selection_sets {
        match select(scene, history, &set.selection) {
            Ok(()) => {
                let _ = history.dispatch(
                    scene,
                    Command::SaveSelectionSet {
                        name: set.name.clone(),
                    },
                );
            }
            Err(e) => out
                .warnings
                .push(format!("selection set `{}`: {e}", set.name)),
        }
    }
    match &file.active {
        Some(active) => {
            if let Err(e) = select(scene, history, active) {
                out.warnings.push(format!("active selection: {e}"));
                let _ = history.dispatch(scene, Command::ClearSelection);
            }
        }
        None => {
            let _ = history.dispatch(scene, Command::ClearSelection);
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::selection::mask_of_one;

    fn fixture(name: &str) -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../fixtures/small")
            .join(name)
    }

    fn temp_path(name: &str) -> PathBuf {
        std::env::temp_dir().join(format!("vizviz_session_test_{}_{name}", std::process::id()))
    }

    #[test]
    fn representation_names_round_trip_and_licorice_is_a_hidden_alias() {
        for r in [
            Representation::Spacefill,
            Representation::BallAndStick,
            Representation::Tube,
            Representation::Cartoon,
            Representation::GaussianSurface,
            Representation::SkinSurface,
            Representation::Ses,
            Representation::Sas,
            Representation::Sticks,
            Representation::Lines,
            Representation::Glycan,
        ] {
            assert_eq!(parse_representation(representation_name(r)), Some(r));
        }
        assert_eq!(representation_name(Representation::Sticks), "sticks");
        assert_eq!(
            parse_representation("licorice"),
            Some(Representation::Sticks),
            "a session saved before the rename"
        );
    }

    #[test]
    fn a_session_round_trips_structures_sets_and_view() {
        let mut scene = Scene::new();
        let mut history = CommandHistory::default();
        crate::script::run_script(
            &mut scene,
            &mut history,
            &format!(
                "load {}\nload {}\nselect chain A and name CA\nsaveset ca\ncolor chain\nrep ballstick 0\nmaterial glossy 0\naddrep cartoon protein",
                fixture("4HHB.cif").display(),
                fixture("1CRN.cif").display()
            ),
        );
        // A picked selection (no expression) on the first structure.
        let first = scene.structures().next().unwrap().0;
        history
            .dispatch(
                &mut scene,
                Command::Select {
                    id: first,
                    mask: mask_of_one(4779, 7),
                },
            )
            .unwrap();
        // A value channel (the FastSASA round trip), colored by.
        let sasa = ValueChannel::new(vec![1.0; 4779], &[4779], 4779, 1).unwrap();
        history
            .dispatch(
                &mut scene,
                Command::SetValues {
                    id: first,
                    name: "sasa".into(),
                    channel: Some(sasa),
                },
            )
            .unwrap();
        history
            .dispatch(
                &mut scene,
                Command::SetColoring {
                    id: first,
                    rep: RepId(0),
                    coloring: ColorScheme::Values("sasa".into()),
                },
            )
            .unwrap();

        let view = serde_json::json!({"style": "white", "zoom": 1.5});
        let path = temp_path("roundtrip.json");
        save(&scene, &path, Some(view.clone())).unwrap();

        let file = read(&path).unwrap();
        assert_eq!(file.schema_version, SCHEMA_VERSION);
        assert_eq!(file.structures.len(), 2);
        assert!(file.structures[0].path.is_absolute());
        assert_eq!(file.structures[0].reps[0].representation, "ball_and_stick");
        assert_eq!(file.structures[0].reps[0].material, "glossy");
        assert_eq!(file.structures[1].reps[0].material, "opaque");
        assert_eq!(file.structures[0].reps[0].coloring, "values:sasa");
        assert_eq!(file.structures[0].values.len(), 1);
        assert_eq!(file.structures[0].values[0].name, "sasa");
        let values_path = &file.structures[0].values[0].path;
        assert!(values_path.is_absolute());
        assert!(values_path.exists());
        assert_eq!(
            values_path.parent().and_then(|p| p.file_name()),
            temp_path("roundtrip.values").file_name()
        );
        assert_eq!(file.structures[1].reps[0].coloring, "chain");
        assert_eq!(file.structures[1].reps.len(), 2);
        assert_eq!(file.structures[1].reps[1].selection, "protein");
        assert_eq!(file.structures[1].current_rep, 1);
        assert!(file.structures[1].values.is_empty());
        assert_eq!(file.selection_sets[0].name, "ca");
        assert_eq!(file.selection_sets[0].selection.structure, 1);
        assert_eq!(
            file.selection_sets[0].selection.expr.as_deref(),
            Some("chain A and name CA")
        );
        assert_eq!(file.selection_sets[0].selection.atoms, None);
        let active = file.active.as_ref().unwrap();
        assert_eq!((active.structure, active.expr.as_deref()), (0, None));
        assert_eq!(active.atoms.as_deref(), Some(&[7][..]));
        assert_eq!(file.view, Some(view));

        // Into a fresh scene that already has something open: it is
        // replaced, and everything comes back.
        let mut scene2 = Scene::new();
        let mut history2 = CommandHistory::default();
        crate::script::run_line(
            &mut scene2,
            &mut history2,
            &format!("load {}", fixture("1UBQ.cif").display()),
        )
        .unwrap();
        let applied = apply(&file, &mut scene2, &mut history2);
        assert!(applied.warnings.is_empty(), "{:?}", applied.warnings);
        let loaded: Vec<_> = scene2.structures().map(|(_, s)| s.label.clone()).collect();
        assert_eq!(loaded, ["4HHB", "1CRN"]);
        let (id0, id1) = (applied.ids[0].unwrap(), applied.ids[1].unwrap());
        assert_eq!(
            scene2.structure(id0).unwrap().rep().representation,
            Representation::BallAndStick
        );
        assert_eq!(
            scene2.structure(id0).unwrap().rep().material,
            Material::Glossy
        );
        assert_eq!(
            scene2.structure(id0).unwrap().rep().coloring,
            ColorScheme::Values("sasa".into())
        );
        let first = scene2.structure(id0).unwrap();
        let (name, channel) = first.values_for(&first.rep().coloring).unwrap();
        assert_eq!(name, "sasa");
        assert_eq!(channel.frame(0), vec![1.0; 4779].as_slice());
        let second = scene2.structure(id1).unwrap();
        assert_eq!(second.reps.len(), 2);
        assert_eq!(second.reps[0].coloring, ColorScheme::Chain);
        assert_eq!(second.rep().representation, Representation::Cartoon);
        assert_eq!(second.rep().selection, "protein");
        let set = scene2.selection_set("ca").unwrap();
        assert_eq!((set.structure, set.mask.count_ones(..)), (id1, 46));
        assert_eq!(set.expr.as_deref(), Some("chain A and name CA"));
        let active = scene2.active_selection().unwrap();
        assert_eq!(active.structure, id0);
        assert!(active.mask[7] && active.mask.count_ones(..) == 1);
        // Undo walks the whole load back, one command at a time.
        while history2.undo(&mut scene2).unwrap() {}
        assert_eq!(scene2.structures().count(), 0);
        let _ = std::fs::remove_file(path);
        let _ = std::fs::remove_dir_all(temp_path("roundtrip.values"));
    }

    /// A 2-model PDB (the exact fixture `vv-io`'s own
    /// `parses_records_models_and_connectivity` test uses), so
    /// `LoadStructure` alone produces a real multi-frame structure --
    /// what a session file actually restores through. A
    /// `LoadTrajectory`-sourced structure's frame does *not* round-trip
    /// (see the `apply()` comment this test also checks): its saved
    /// `path` is the topology file alone, so reloading always comes back
    /// single-frame.
    const TWO_MODEL_PDB: &str = "HEADER    PLANT PROTEIN                           30-APR-81   1CRN
MODEL        1
ATOM      1  N   THR A   1      17.047  14.099   3.625  1.00 13.79           N
ATOM      2  CA  THR A   1      16.967  12.784   4.338  1.00 10.80           C
ENDMDL
MODEL        2
ATOM      1  N   THR A   1      17.147  14.099   3.625  1.00 13.79           N
ATOM      2  CA  THR A   1      17.067  12.784   4.338  1.00 10.80           C
ENDMDL
END
";

    #[test]
    fn a_structures_frame_round_trips_through_a_session() {
        let pdb_path = temp_path("two_model.pdb");
        std::fs::write(&pdb_path, TWO_MODEL_PDB).unwrap();

        let mut scene = Scene::new();
        let mut history = CommandHistory::default();
        history
            .dispatch(
                &mut scene,
                Command::LoadStructure {
                    path: pdb_path.clone(),
                },
            )
            .unwrap();
        let id = scene.structures().next().unwrap().0;
        assert_eq!(scene.structure(id).unwrap().structure.frame_count(), 2);
        history
            .dispatch(&mut scene, Command::SetFrame { id, frame: 1 })
            .unwrap();

        let path = temp_path("frame_roundtrip.json");
        save(&scene, &path, None).unwrap();
        let file = read(&path).unwrap();
        assert_eq!(file.structures[0].frame, 1);

        let mut scene2 = Scene::new();
        let mut history2 = CommandHistory::default();
        let applied = apply(&file, &mut scene2, &mut history2);
        assert!(applied.warnings.is_empty(), "{:?}", applied.warnings);
        let id2 = applied.ids[0].unwrap();
        assert_eq!(scene2.structure(id2).unwrap().frame, 1);

        // A file saved before `frame` existed (the key absent from JSON,
        // not just `0`) still loads, at frame 0 -- `#[serde(default)]`,
        // not a parse failure.
        let mut old_json: serde_json::Value =
            serde_json::from_str(&std::fs::read_to_string(&path).unwrap()).unwrap();
        old_json["structures"][0]
            .as_object_mut()
            .unwrap()
            .remove("frame");
        std::fs::write(&path, serde_json::to_string(&old_json).unwrap()).unwrap();
        let old_file = read(&path).unwrap();
        assert_eq!(old_file.structures[0].frame, 0);

        let _ = std::fs::remove_file(path);
        let _ = std::fs::remove_file(pdb_path);
    }

    #[test]
    fn a_hidden_structure_stays_hidden_through_a_session_round_trip() {
        let mut scene = Scene::new();
        let mut history = CommandHistory::default();
        crate::script::run_line(
            &mut scene,
            &mut history,
            &format!("load {}", fixture("1CRN.cif").display()),
        )
        .unwrap();
        let id = scene.structures().next().unwrap().0;
        history
            .dispatch(&mut scene, Command::ShowStructure { id, visible: false })
            .unwrap();

        let path = temp_path("hidden_roundtrip.json");
        save(&scene, &path, None).unwrap();
        let file = read(&path).unwrap();
        assert!(!file.structures[0].visible);

        let mut scene2 = Scene::new();
        let mut history2 = CommandHistory::default();
        let applied = apply(&file, &mut scene2, &mut history2);
        assert!(applied.warnings.is_empty(), "{:?}", applied.warnings);
        let id2 = applied.ids[0].unwrap();
        assert!(!scene2.structure(id2).unwrap().visible);

        // A file saved before `visible` existed loads as shown.
        let mut old_json: serde_json::Value =
            serde_json::from_str(&std::fs::read_to_string(&path).unwrap()).unwrap();
        old_json["structures"][0]
            .as_object_mut()
            .unwrap()
            .remove("visible");
        std::fs::write(&path, serde_json::to_string(&old_json).unwrap()).unwrap();
        let old_file = read(&path).unwrap();
        assert!(old_file.structures[0].visible);

        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn a_trajectory_its_frame_and_measurements_survive_a_session_round_trip() {
        let mut scene = Scene::new();
        let mut history = CommandHistory::default();
        let run = |scene: &mut Scene, history: &mut CommandHistory, line: String| {
            crate::script::run_line(scene, history, &line).unwrap()
        };
        run(
            &mut scene,
            &mut history,
            format!("load {}", fixture("1CRN.pdb").display()),
        );
        run(
            &mut scene,
            &mut history,
            format!("attachtrajectory {}", fixture("1CRN_traj.dcd").display()),
        );
        let id = scene.structures().next().unwrap().0;
        history
            .dispatch(&mut scene, Command::SetFrame { id, frame: 2 })
            .unwrap();
        crate::script::run_line(&mut scene, &mut history, "measure 0 1 2").unwrap();
        crate::script::run_line(&mut scene, &mut history, "label 5 the fifth").unwrap();
        crate::script::run_line(&mut scene, &mut history, "caption title 0.1 0.2 Crambin").unwrap();

        let path = temp_path("trajectory_frame.json");
        save(&scene, &path, None).unwrap();
        let file = read(&path).unwrap();
        assert!(file.structures[0].trajectory.is_some());

        let mut scene2 = Scene::new();
        let mut history2 = CommandHistory::default();
        let applied = apply(&file, &mut scene2, &mut history2);
        assert!(applied.warnings.is_empty(), "{:?}", applied.warnings);
        let loaded = scene2.structure(applied.ids[0].unwrap()).unwrap();
        assert_eq!(loaded.structure.frame_count(), 3);
        assert_eq!(loaded.frame, 2);
        assert_eq!(loaded.measurements.len(), 1);
        assert_eq!(loaded.measurements[0].atoms(), [0, 1, 2]);
        assert_eq!(loaded.labels.get(&5).map(String::as_str), Some("the fifth"));
        assert_eq!(scene2.captions(), scene.captions());

        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn coloring_by_a_never_attached_channel_warns_on_load() {
        let mut scene = Scene::new();
        let mut history = CommandHistory::default();
        crate::script::run_line(
            &mut scene,
            &mut history,
            &format!("load {}", fixture("1CRN.cif").display()),
        )
        .unwrap();
        let id = scene.structures().next().unwrap().0;
        history
            .dispatch(
                &mut scene,
                Command::SetColoring {
                    id,
                    rep: RepId(0),
                    coloring: ColorScheme::Values("sasa".into()),
                },
            )
            .unwrap();
        let file = capture(&scene, None, &temp_path("no_such_values")).unwrap();
        assert!(file.structures[0].values.is_empty());

        let mut scene2 = Scene::new();
        let mut history2 = CommandHistory::default();
        let applied = apply(&file, &mut scene2, &mut history2);
        assert_eq!(applied.warnings.len(), 1, "{:?}", applied.warnings);
        assert!(
            applied.warnings[0].contains("`sasa`"),
            "{:?}",
            applied.warnings
        );
        assert_eq!(
            scene2
                .structure(applied.ids[0].unwrap())
                .unwrap()
                .rep()
                .coloring,
            ColorScheme::Values("sasa".into())
        );
    }

    #[test]
    fn a_missing_sidecar_file_is_a_warning_not_a_crash() {
        let mut scene = Scene::new();
        let mut history = CommandHistory::default();
        crate::script::run_line(
            &mut scene,
            &mut history,
            &format!("load {}", fixture("1CRN.cif").display()),
        )
        .unwrap();
        let id = scene.structures().next().unwrap().0;
        let atoms = scene.structure(id).unwrap().structure.atom_count();
        history
            .dispatch(
                &mut scene,
                Command::SetValues {
                    id,
                    name: "sasa".into(),
                    channel: Some(ValueChannel::new(vec![0.0; atoms], &[atoms], atoms, 1).unwrap()),
                },
            )
            .unwrap();
        let values_dir = temp_path("deleted_values");
        let file = capture(&scene, None, &values_dir).unwrap();
        assert_eq!(file.structures[0].values.len(), 1);
        std::fs::remove_dir_all(&values_dir).unwrap();

        let mut scene2 = Scene::new();
        let mut history2 = CommandHistory::default();
        let applied = apply(&file, &mut scene2, &mut history2);
        assert_eq!(applied.warnings.len(), 1, "{:?}", applied.warnings);
        assert!(
            applied.warnings[0].contains("`sasa`"),
            "{:?}",
            applied.warnings
        );
        assert!(scene2
            .structure(applied.ids[0].unwrap())
            .unwrap()
            .values
            .is_empty());
    }

    #[test]
    fn missing_files_are_warnings_not_errors() {
        let file = SessionFile {
            schema_version: SCHEMA_VERSION,
            structures: vec![
                SavedStructure {
                    path: PathBuf::from("/no/such/file.cif"),
                    trajectory: None,
                    label: "gone".into(),
                    visible: true,
                    reps: Vec::new(),
                    current_rep: 0,
                    representation: "spacefill".into(),
                    coloring: "element".into(),
                    material: String::new(),
                    values: Vec::new(),
                    measurements: Vec::new(),
                    frame: 0,
                    labels: Default::default(),
                },
                SavedStructure {
                    path: fixture("1CRN.cif"),
                    trajectory: None,
                    label: "1CRN".into(),
                    visible: true,
                    reps: Vec::new(),
                    current_rep: 0,
                    representation: "spacefill".into(),
                    coloring: "element".into(),
                    material: String::new(),
                    values: Vec::new(),
                    measurements: Vec::new(),
                    frame: 0,
                    labels: Default::default(),
                },
            ],
            selection_sets: vec![SavedSet {
                name: "lost".into(),
                selection: SavedSelection {
                    structure: 0,
                    expr: Some("all".into()),
                    atoms: None,
                },
            }],
            active: Some(SavedSelection {
                structure: 1,
                expr: Some("name CA".into()),
                atoms: None,
            }),
            view: None,
            captions: Vec::new(),
        };
        let mut scene = Scene::new();
        let mut history = CommandHistory::default();
        let applied = apply(&file, &mut scene, &mut history);
        assert_eq!(applied.ids[0], None);
        assert!(applied.ids[1].is_some());
        assert_eq!(applied.warnings.len(), 2, "{:?}", applied.warnings);
        assert!(applied.warnings[0].contains("structure #0"));
        assert!(applied.warnings[1].contains("`lost`"));
        assert_eq!(scene.selection_sets().len(), 0);
        assert_eq!(scene.active_selection().unwrap().mask.count_ones(..), 46);
    }

    #[test]
    fn read_rejects_garbage_and_future_schemas() {
        let path = temp_path("garbage.json");
        std::fs::write(&path, "{\"schema_version\": 1, \"structures\": 5}").unwrap();
        assert!(matches!(read(&path), Err(SessionError::Format { .. })));
        std::fs::write(&path, "{\"schema_version\": 99, \"structures\": []}").unwrap();
        assert!(matches!(
            read(&path),
            Err(SessionError::Schema { found: 99, .. })
        ));
        let _ = std::fs::remove_file(&path);
        assert!(matches!(
            read(Path::new("/nope.json")),
            Err(SessionError::Io { .. })
        ));
    }
}
