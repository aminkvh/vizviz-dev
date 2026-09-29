//! `vizviz.Session`: the `vv-scene` command bus from Python.
//!
//! Every mutation goes through `CommandHistory::dispatch`, exactly as a
//! click in the desktop app does, so `undo()`/`redo()` work and a script
//! can't put the scene into a state the GUI couldn't. Structure ids are
//! the plain `u32` slot indices of `vv_scene::StructureId` (never reused
//! within a session, so a stale id fails loudly instead of aliasing).

use std::path::PathBuf;
use std::sync::Arc;

use fixedbitset::FixedBitSet;
use numpy::{PyArray1, PyArrayMethods, PyReadonlyArray1, PyReadonlyArray2, PyUntypedArrayMethods};
use pyo3::exceptions::{PyIndexError, PyKeyError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use vv_scene::{
    ColorScheme, Command, CommandHistory, LoadedStructure, Mask, Representation, Scene, SceneError,
    StructureId, ValueChannel,
};

use crate::structure::Structure;

fn sid(id: u32) -> StructureId {
    StructureId::from_raw(id)
}

fn scene_err(e: SceneError) -> PyErr {
    match e {
        SceneError::NoSuchStructure(id) => {
            PyKeyError::new_err(format!("no structure with id {}", id.to_raw()))
        }
        e @ SceneError::NoSuchRep(..) => PyKeyError::new_err(e.to_string()),
        e @ SceneError::LastRep(_) => PyValueError::new_err(e.to_string()),
        SceneError::NoSuchSelectionSet(name) => {
            PyKeyError::new_err(format!("no selection set named {name:?}"))
        }
        SceneError::NoActiveSelection => PyValueError::new_err("nothing is selected"),
        e @ SceneError::SelectionSetNameTaken(_) => PyValueError::new_err(e.to_string()),
        SceneError::Load { path, source } => crate::parse_err(&path, source),
        e @ SceneError::LoadTrajectory { .. } => PyValueError::new_err(e.to_string()),
        e @ SceneError::TrajectoryAtomMismatch { .. } => PyValueError::new_err(e.to_string()),
        e @ SceneError::Structure(_) => PyValueError::new_err(e.to_string()),
        SceneError::Selection { source, .. } => crate::select_err(source),
        SceneError::NoSuchValues(name) => {
            PyKeyError::new_err(format!("no value channel named {name:?}"))
        }
        e @ SceneError::ValuesLength { .. } => PyValueError::new_err(e.to_string()),
        e @ SceneError::NoSuchAtom { .. } => PyIndexError::new_err(e.to_string()),
        SceneError::NoSuchCaption(name) => {
            PyKeyError::new_err(format!("no caption named {name:?}"))
        }
        e @ SceneError::FrameOutOfRange { .. } => PyIndexError::new_err(e.to_string()),
    }
}

/// Delegates to the same alias-tolerant parser the console/`exec()` use
/// (`vv_scene::script::parse_representation`), so `set_representation`
/// accepts everything `rep`/`Session.exec("rep ...")` does.
fn parse_representation(name: &str) -> PyResult<Representation> {
    vv_scene::parse_representation(name).map_err(|_| {
        PyValueError::new_err(format!(
            "unknown representation {name:?}; expected \"spacefill\", \"ball_and_stick\", \"sticks\", \"lines\", \"sas\", \"tube\", \"cartoon\", \"gaussian_surface\", or \"skin_surface\""
        ))
    })
}

fn parse_scene_coloring(name: &str) -> PyResult<ColorScheme> {
    ColorScheme::parse(name).ok_or_else(|| {
        PyValueError::new_err(format!(
            "unknown coloring {name:?}; expected \"element\", \"chain\", \"b_factor\", or \"values:NAME\""
        ))
    })
}

/// A boolean mask of length `atom_count` (NumPy array or list of bools),
/// or any sequence of atom indices (a list, range, or integer array).
fn mask_from_python(atom_count: usize, atoms: &Bound<'_, PyAny>) -> PyResult<Mask> {
    let mut bits = FixedBitSet::with_capacity(atom_count);
    let bools: Option<Vec<bool>> = if let Ok(mask) = atoms.extract::<PyReadonlyArray1<bool>>() {
        Some(mask.as_array().iter().copied().collect())
    } else {
        // `bool` extraction rejects plain ints, so `[0, 1, 2]` falls through
        // to the index path while `[True, False, ...]` is treated as a mask.
        atoms.extract::<Vec<bool>>().ok()
    };
    if let Some(mask) = bools {
        if mask.len() != atom_count {
            return Err(PyValueError::new_err(format!(
                "mask has {} entries but the structure has {atom_count} atoms",
                mask.len()
            )));
        }
        for (i, on) in mask.into_iter().enumerate() {
            if on {
                bits.insert(i);
            }
        }
    } else {
        let indices: Vec<i64> = atoms.extract().map_err(|_| {
            PyTypeError::new_err("atoms must be a boolean mask or a sequence of atom indices")
        })?;
        for i in indices {
            if i < 0 || i as usize >= atom_count {
                return Err(PyIndexError::new_err(format!(
                    "atom index {i} out of range for {atom_count} atoms"
                )));
            }
            bits.insert(i as usize);
        }
    }
    Ok(Arc::new(bits))
}

fn indices_of<'py>(py: Python<'py>, mask: &Mask) -> Bound<'py, PyArray1<u32>> {
    PyArray1::from_iter(py, mask.ones().map(|i| i as u32))
}

#[pyclass(name = "Session", module = "vizviz")]
pub struct Session {
    scene: Scene,
    history: CommandHistory,
}

impl Session {
    fn dispatch(&mut self, command: Command) -> PyResult<()> {
        self.history
            .dispatch(&mut self.scene, command)
            .map_err(scene_err)
    }

    fn loaded(&self, id: u32) -> PyResult<&LoadedStructure> {
        self.scene
            .structure(sid(id))
            .ok_or_else(|| scene_err(SceneError::NoSuchStructure(sid(id))))
    }

    fn newest(&self) -> Option<StructureId> {
        self.scene.structures().last().map(|(id, _)| id)
    }
}

#[pymethods]
impl Session {
    /// `undo_limit` caps the history (a closed structure stays in memory
    /// until its close is pushed out of the undo stack).
    #[new]
    #[pyo3(signature = (undo_limit = 100))]
    fn new(undo_limit: usize) -> Self {
        Self {
            scene: Scene::new(),
            history: CommandHistory::new(undo_limit),
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "Session(structures={}, selection_sets={}, version={})",
            self.scene.structures().count(),
            self.scene.selection_sets().len(),
            self.scene.version()
        )
    }

    // ---- structures ------------------------------------------------------

    /// Loads a file and returns its structure id. Releases the GIL while
    /// parsing.
    fn load(&mut self, py: Python<'_>, path: PathBuf) -> PyResult<u32> {
        let (scene, history) = (&mut self.scene, &mut self.history);
        py.detach(|| history.dispatch(scene, Command::LoadStructure { path }))
            .map_err(scene_err)?;
        Ok(self.newest().expect("just loaded").to_raw())
    }

    /// Loads a topology file plus a trajectory (DCD, XTC, TRR, NetCDF) as one multi-frame
    /// structure. Releases the GIL while parsing.
    fn load_trajectory(
        &mut self,
        py: Python<'_>,
        topology: PathBuf,
        trajectory: PathBuf,
    ) -> PyResult<u32> {
        let (scene, history) = (&mut self.scene, &mut self.history);
        py.detach(|| {
            history.dispatch(
                scene,
                Command::LoadTrajectory {
                    topology,
                    trajectory,
                },
            )
        })
        .map_err(scene_err)?;
        Ok(self.newest().expect("just loaded").to_raw())
    }

    fn close(&mut self, id: u32) -> PyResult<()> {
        self.dispatch(Command::CloseStructure { id: sid(id) })
    }

    /// Ids of the currently loaded structures, oldest first.
    #[getter]
    fn structures(&self) -> Vec<u32> {
        self.scene.structures().map(|(id, _)| id.to_raw()).collect()
    }

    /// The structure behind `id`. Shares memory with the session (no
    /// copy) and stays valid after the session closes or drops it.
    fn structure(&self, id: u32) -> PyResult<Structure> {
        Ok(Structure::new(self.loaded(id)?.structure.clone()))
    }

    fn label(&self, id: u32) -> PyResult<String> {
        Ok(self.loaded(id)?.label.clone())
    }

    fn path(&self, id: u32) -> PyResult<Option<PathBuf>> {
        Ok(self.loaded(id)?.path.clone())
    }

    /// The current rep's style: `"spacefill"`, `"ball_and_stick"`,
    /// `"tube"`, `"cartoon"`, `"gaussian_surface"` or `"skin_surface"`.
    fn representation(&self, id: u32) -> PyResult<&'static str> {
        Ok(vv_scene::session::representation_name(
            self.loaded(id)?.rep().representation,
        ))
    }

    fn set_representation(&mut self, id: u32, representation: &str) -> PyResult<()> {
        let representation = parse_representation(representation)?;
        let rep = self.loaded(id)?.rep().id;
        self.dispatch(Command::SetRepresentation {
            id: sid(id),
            rep,
            representation,
        })
    }

    /// The current coordinate-set index (0-based) of structure `id`.
    /// Always `0` for a single-frame structure.
    fn frame(&self, id: u32) -> PyResult<usize> {
        Ok(self.loaded(id)?.frame)
    }

    /// Sets structure `id`'s current frame. Undoable; this is the value
    /// `save_session` persists. Errors if `frame` is out of range.
    fn set_frame(&mut self, id: u32, frame: usize) -> PyResult<()> {
        self.dispatch(Command::SetFrame { id: sid(id), frame })
    }

    /// The current rep's colouring: `"element"`, `"chain"`, `"b_factor"`,
    /// or `"values:NAME"`, among others.
    fn coloring(&self, id: u32) -> PyResult<String> {
        Ok(self.loaded(id)?.rep().coloring.name())
    }

    /// Same names as `coloring()`; `"values:NAME"` colors by an attached
    /// value channel (see `set_values`).
    fn set_coloring(&mut self, id: u32, coloring: &str) -> PyResult<()> {
        let coloring = parse_scene_coloring(coloring)?;
        let rep = self.loaded(id)?.rep().id;
        self.dispatch(Command::SetColoring {
            id: sid(id),
            rep,
            coloring,
        })
    }

    // ---- per-atom value channels -----------------------------------------

    /// Attaches `values` to structure `id` under `name`: one number per
    /// atom or per residue, or `(frames, n)` of those for a trajectory.
    /// Then `set_coloring(id, "values:NAME")` (or `exec("color values
    /// NAME")`) colors by it, with one color scale over all frames.
    /// Undoable. `save_session` persists channels as `.npy` sidecars next
    /// to the session file.
    fn set_values(&mut self, id: u32, name: &str, values: &Bound<'_, PyAny>) -> PyResult<()> {
        let loaded = self.loaded(id)?;
        let (data, shape) = values_and_shape(values)?;
        let channel = ValueChannel::for_structure(data, &shape, &loaded.structure)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        self.dispatch(Command::SetValues {
            id: sid(id),
            name: name.to_owned(),
            channel: Some(channel),
        })
    }

    /// Names of the value channels attached to structure `id`.
    fn values(&self, id: u32) -> PyResult<Vec<String>> {
        Ok(self.loaded(id)?.values.keys().cloned().collect())
    }

    /// A copy of channel `name` of structure `id`: `(atoms,)` float32, or
    /// `(frames, atoms)` for a per-frame channel.
    fn get_values<'py>(&self, py: Python<'py>, id: u32, name: &str) -> PyResult<Bound<'py, PyAny>> {
        let loaded = self.loaded(id)?;
        let channel = loaded
            .values
            .get(name)
            .ok_or_else(|| PyKeyError::new_err(format!("no value channel named {name:?}")))?;
        let flat = PyArray1::from_slice(py, channel.data());
        if channel.frames() > 1 {
            Ok(flat
                .reshape([channel.frames(), channel.atoms()])?
                .into_any())
        } else {
            Ok(flat.into_any())
        }
    }

    /// Detaches channel `name` from structure `id`. Undoable.
    fn remove_values(&mut self, id: u32, name: &str) -> PyResult<()> {
        self.dispatch(Command::SetValues {
            id: sid(id),
            name: name.to_owned(),
            channel: None,
        })
    }

    // ---- selection -------------------------------------------------------

    /// Makes `atoms` of structure `id` the active selection. `atoms` is a
    /// selection expression string (`"chain A and name CA"`, see
    /// docs/SELECTION.md), a boolean mask of length `atom_count`, or a
    /// sequence of atom indices.
    fn select(&mut self, py: Python<'_>, id: u32, atoms: &Bound<'_, PyAny>) -> PyResult<()> {
        if let Ok(expr) = atoms.extract::<String>() {
            let (scene, history) = (&mut self.scene, &mut self.history);
            let id = sid(id);
            return py
                .detach(|| history.dispatch(scene, Command::SelectExpr { id, expr }))
                .map_err(scene_err);
        }
        let n = self.loaded(id)?.structure.atom_count();
        let mask = mask_from_python(n, atoms)?;
        self.dispatch(Command::Select { id: sid(id), mask })
    }

    fn clear_selection(&mut self) -> PyResult<()> {
        self.dispatch(Command::ClearSelection)
    }

    /// `(structure_id, indices)` of the active selection, or `None`.
    #[getter]
    fn selection<'py>(&self, py: Python<'py>) -> Option<(u32, Bound<'py, PyArray1<u32>>)> {
        self.scene
            .active_selection()
            .map(|a| (a.structure.to_raw(), indices_of(py, &a.mask)))
    }

    /// The expression the active selection was made with, or `None` if it
    /// came from a mask, indices, or a pick.
    #[getter]
    fn selection_expr(&self) -> Option<String> {
        self.scene.active_selection().and_then(|a| a.expr.clone())
    }

    /// The expression the saved set `name` was made with, if any.
    fn selection_set_expr(&self, name: &str) -> PyResult<Option<String>> {
        let set = self
            .scene
            .selection_set(name)
            .ok_or_else(|| scene_err(SceneError::NoSuchSelectionSet(name.into())))?;
        Ok(set.expr.clone())
    }

    /// Saves the active selection under `name` (overwriting any set with
    /// that name).
    fn save_selection_set(&mut self, name: &str) -> PyResult<()> {
        self.dispatch(Command::SaveSelectionSet { name: name.into() })
    }

    fn delete_selection_set(&mut self, name: &str) -> PyResult<()> {
        self.dispatch(Command::DeleteSelectionSet { name: name.into() })
    }

    /// Names of the saved selection sets.
    #[getter]
    fn selection_sets(&self) -> Vec<String> {
        self.scene
            .selection_sets()
            .iter()
            .map(|s| s.name.clone())
            .collect()
    }

    /// `(structure_id, indices)` of the saved set `name`.
    fn selection_set<'py>(
        &self,
        py: Python<'py>,
        name: &str,
    ) -> PyResult<(u32, Bound<'py, PyArray1<u32>>)> {
        let set = self
            .scene
            .selection_set(name)
            .ok_or_else(|| scene_err(SceneError::NoSuchSelectionSet(name.into())))?;
        Ok((set.structure.to_raw(), indices_of(py, &set.mask)))
    }

    // ---- scripting -------------------------------------------------------

    /// Runs text commands (`"load x.cif; select chain A; color chain"`,
    /// one per line or `;`-separated; see docs/COMMANDS.md) and returns
    /// each command's result line. Stops at the first failure and raises
    /// `ValueError` naming the command that failed. Only document verbs
    /// exist here; view verbs (style, screenshot) belong to the app.
    fn exec(&mut self, py: Python<'_>, script: &str) -> PyResult<Vec<String>> {
        let (scene, history) = (&mut self.scene, &mut self.history);
        let (outputs, failure) = py.detach(|| vv_scene::run_script(scene, history, script));
        match failure {
            None => Ok(outputs),
            Some((line, e)) => Err(PyValueError::new_err(format!("`{line}`: {e}"))),
        }
    }

    /// Writes a session file (docs/COMMANDS.md): structure paths,
    /// representations, colorings, selection sets, and the active
    /// selection. The desktop app reads the same files (and adds its
    /// camera and style, which a headless session has no use for).
    fn save_session(&self, path: PathBuf) -> PyResult<()> {
        vv_scene::session::save(&self.scene, &path, None)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    /// Replaces what is loaded with a session file's contents. Returns
    /// warnings for anything that could not be restored (a moved file, a
    /// selection on it); a missing file is a warning, not an error.
    fn load_session(&mut self, py: Python<'_>, path: PathBuf) -> PyResult<Vec<String>> {
        let file =
            vv_scene::session::read(&path).map_err(|e| PyValueError::new_err(e.to_string()))?;
        let (scene, history) = (&mut self.scene, &mut self.history);
        let applied = py.detach(|| vv_scene::session::apply(&file, scene, history));
        Ok(applied.warnings)
    }

    /// The commands `exec` understands, as `(id, usage, help)` tuples: an
    /// agent-readable manifest of what a headless session can do.
    #[staticmethod]
    fn commands() -> Vec<(&'static str, &'static str, &'static str)> {
        vv_scene::COMMAND_SPECS
            .iter()
            .map(|s| (s.id, s.usage, s.help))
            .collect()
    }

    // ---- history ---------------------------------------------------------

    /// Undoes the last command; `False` if there was nothing to undo.
    fn undo(&mut self) -> PyResult<bool> {
        self.history.undo(&mut self.scene).map_err(scene_err)
    }

    /// Redoes the last undone command; `False` if there was nothing to redo.
    fn redo(&mut self) -> PyResult<bool> {
        self.history.redo(&mut self.scene).map_err(scene_err)
    }

    #[getter]
    fn can_undo(&self) -> bool {
        self.history.can_undo()
    }

    #[getter]
    fn can_redo(&self) -> bool {
        self.history.can_redo()
    }

    /// Bumped by every applied command (including undo/redo); compare
    /// against a remembered value to know whether anything changed.
    #[getter]
    fn version(&self) -> u64 {
        self.scene.version()
    }

    // ---- rendering -------------------------------------------------------

    /// Renders structure `id` (default: the most recently loaded one) with
    /// its current coloring, using the CPU backend (spacefill only), and
    /// returns an `(height, width, 4)` uint8 RGBA image.
    #[cfg(feature = "render")]
    #[pyo3(signature = (width, height, *, id = None, style = "dark_presentation",
                        yaw = 0.0, pitch = 0.0, zoom = 1.0, background = None))]
    #[allow(clippy::too_many_arguments)]
    fn render<'py>(
        &self,
        py: Python<'py>,
        width: u32,
        height: u32,
        id: Option<u32>,
        style: &str,
        yaw: f32,
        pitch: f32,
        zoom: f32,
        background: Option<(u8, u8, u8, u8)>,
    ) -> PyResult<Bound<'py, numpy::PyArray3<u8>>> {
        let id = match id {
            Some(id) => id,
            None => self
                .newest()
                .ok_or_else(|| PyValueError::new_err("no structure is loaded"))?
                .to_raw(),
        };
        let loaded = self.loaded(id)?;
        if matches!(
            loaded.rep().representation,
            Representation::BallAndStick | Representation::Tube | Representation::Cartoon
        ) {
            PyErr::warn(
                py,
                py.get_type::<pyo3::exceptions::PyUserWarning>().as_any(),
                c"the CPU render backend draws bonds/tubes/cartoon ribbons as spacefill (spacefill, gaussian_surface, and skin_surface render as themselves)",
                1,
            )?;
        }
        let options = crate::render::RenderOptions {
            width,
            height,
            style,
            yaw,
            pitch,
            zoom,
            background,
        };
        let colors = crate::render::colors_of(loaded, &loaded.rep().coloring, 0);
        crate::render::render_to_array(
            py,
            &loaded.structure,
            colors,
            loaded.rep().representation,
            &options,
        )
    }
}

/// A 1-D or 2-D numeric array (or nested lists) as f32 with its shape.
fn values_and_shape(values: &Bound<'_, PyAny>) -> PyResult<(Vec<f32>, Vec<usize>)> {
    macro_rules! try_2d {
        ($t:ty) => {
            if let Ok(a) = values.extract::<PyReadonlyArray2<$t>>() {
                let shape = a.shape().to_vec();
                return Ok((a.as_array().iter().map(|&v| v as f32).collect(), shape));
            }
        };
    }
    try_2d!(f32);
    try_2d!(f64);
    try_2d!(i64);
    try_2d!(i32);
    if let Ok(rows) = values.extract::<Vec<Vec<f32>>>() {
        let cols = rows.first().map_or(0, Vec::len);
        if rows.iter().any(|r| r.len() != cols) {
            return Err(PyValueError::new_err("values rows have different lengths"));
        }
        let shape = vec![rows.len(), cols];
        return Ok((rows.into_iter().flatten().collect(), shape));
    }
    let column = crate::render::scalar_column(values)?;
    let n = column.len();
    Ok((column, vec![n]))
}
