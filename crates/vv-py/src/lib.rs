//! `vizviz._core`: the compiled half of the `vizviz` Python package.
//!
//! Three things live here. `load()` + [`Structure`] expose a parsed file's
//! per-atom columns as read-only, zero-copy NumPy views. [`Session`] wraps
//! the `vv-scene` command bus (load/close/select/undo/redo) so a script
//! does exactly what a click in the desktop app does. Behind the `render`
//! feature, `Structure.render()`/`Session.render()` draw a frame with the
//! CPU backend — no GPU, no window — and `write_png()` saves it.

mod arrays;
mod session;
mod structure;

#[cfg(feature = "render")]
mod render;

use std::path::{Path, PathBuf};

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

pub use session::Session;
pub use structure::Structure;

/// Element symbol for an atomic number, or `None` if out of range.
#[pyfunction]
fn element_symbol(atomic_number: u8) -> Option<&'static str> {
    vv_core::Element::from_atomic_number(atomic_number).map(|e| e.symbol())
}

/// Reads an mmCIF or PDB file (optionally gzip-compressed) into a
/// [`Structure`]. The GIL is released while parsing, so other Python
/// threads keep running during a multi-second load.
#[pyfunction]
fn load(py: Python<'_>, path: PathBuf) -> PyResult<Structure> {
    let inner = py
        .detach(|| vv_io::load(&path))
        .map_err(|e| parse_err(&path, e))?;
    Ok(Structure::new(inner))
}

/// A selection-expression error, with the byte span of the offending
/// token so a caller can point at it.
pub(crate) fn select_err(e: vv_core::SelectError) -> PyErr {
    PyValueError::new_err(format!(
        "{} (at characters {}..{})",
        e.message, e.span.start, e.span.end
    ))
}

/// Downloads a PDB entry into the per-user cache (or reuses the cached
/// copy) and returns its path; `assembly=N` fetches biological assembly
/// N instead of the asymmetric unit. Pair with `load()`.
#[pyfunction]
#[pyo3(signature = (id, assembly = None))]
fn fetch(py: Python<'_>, id: &str, assembly: Option<u32>) -> PyResult<PathBuf> {
    let which = match assembly {
        Some(n) => vv_io::fetch::Assembly::Biological(n.max(1)),
        None => vv_io::fetch::Assembly::AsymmetricUnit,
    };
    let cache = vv_io::fetch::cache_dir()
        .ok_or_else(|| PyValueError::new_err("no cache directory (no home directory?)"))?;
    py.detach(|| vv_io::fetch::fetch(id, which, &cache))
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Reads many files in parallel (rayon, GIL released) and returns the
/// structures in the same order. Fails on the first unreadable file,
/// naming it. For millions of files, call this in batches and drop each
/// batch's structures before the next: memory is the limit, not time.
#[pyfunction]
fn load_many(py: Python<'_>, paths: Vec<PathBuf>) -> PyResult<Vec<Structure>> {
    use rayon::prelude::*;
    let loaded: Vec<Result<vv_core::Structure, (PathBuf, vv_io::ParseError)>> = py.detach(|| {
        paths
            .par_iter()
            .map(|path| vv_io::load(path).map_err(|e| (path.clone(), e)))
            .collect()
    });
    loaded
        .into_iter()
        .map(|r| match r {
            Ok(s) => Ok(Structure::new(s)),
            Err((path, e)) => Err(parse_err(&path, e)),
        })
        .collect()
}

/// I/O failures become the matching built-in `OSError` subclass
/// (`FileNotFoundError`, `PermissionError`, ...); everything else about a
/// file we could open but not understand is a `ValueError`.
pub(crate) fn parse_err(path: &Path, e: vv_io::ParseError) -> PyErr {
    match e {
        vv_io::ParseError::Io { source, .. } => {
            std::io::Error::new(source.kind(), format!("{}: {source}", path.display())).into()
        }
        other => PyValueError::new_err(format!("{}: {other}", path.display())),
    }
}

/// Same split as `parse_err`: an I/O failure becomes the matching
/// `OSError` subclass, an unrecognized extension a `ValueError`.
pub(crate) fn save_err(e: vv_io::SaveError) -> PyErr {
    match e {
        vv_io::SaveError::Io { path, source } => {
            std::io::Error::new(source.kind(), format!("{path}: {source}")).into()
        }
        other => PyValueError::new_err(other.to_string()),
    }
}

#[pymodule]
fn _core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add("HAS_RENDER", cfg!(feature = "render"))?;
    m.add_function(wrap_pyfunction!(element_symbol, m)?)?;
    m.add_function(wrap_pyfunction!(load, m)?)?;
    m.add_function(wrap_pyfunction!(load_many, m)?)?;
    m.add_function(wrap_pyfunction!(fetch, m)?)?;
    m.add_function(wrap_pyfunction!(arrays::from_columns, m)?)?;
    m.add_function(wrap_pyfunction!(arrays::guess_elements, m)?)?;
    m.add_class::<Structure>()?;
    m.add_class::<Session>()?;
    #[cfg(feature = "render")]
    m.add_function(wrap_pyfunction!(render::write_png, m)?)?;
    Ok(())
}
