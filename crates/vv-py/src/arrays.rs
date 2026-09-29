//! `vizviz.from_arrays`: build a `Structure` from plain NumPy columns, so
//! coordinates from any other package (MDAnalysis, mdtraj, RDKit, a
//! simulation engine, your own code) can be rendered and analysed without
//! going through a file.
//!
//! The Python wrapper in `python/vizviz/interop.py` does the friendly
//! normalization (strings to interned indices, elements from symbols or
//! atom names, defaults); this is the strict, typed entry point it calls.

use numpy::{PyArray1, PyFixedString, PyReadonlyArray1, PyReadonlyArray3, PyUntypedArrayMethods};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use vv_core::glam::Vec3;
use vv_core::{AtomRow, Element, TopologyBuilder};

use crate::structure::Structure;

pub struct Columns<'py> {
    pub positions: PyReadonlyArray3<'py, f32>,
    pub element: PyReadonlyArray1<'py, u8>,
    pub name: PyReadonlyArray1<'py, PyFixedString<4>>,
    pub residue_names: Vec<String>,
    pub residue_name_index: PyReadonlyArray1<'py, u32>,
    pub residue_ids: PyReadonlyArray1<'py, i32>,
    pub chain_names: Vec<String>,
    pub chain_index: PyReadonlyArray1<'py, u32>,
    pub auth_chain_names: Vec<String>,
    pub auth_chain_index: PyReadonlyArray1<'py, u32>,
    pub b_factor: PyReadonlyArray1<'py, f32>,
    pub occupancy: PyReadonlyArray1<'py, f32>,
    pub hetero: PyReadonlyArray1<'py, bool>,
}

/// Builds the structure. Residues are runs of atoms with the same
/// (chain, residue id, residue name); chains are runs of the same chain
/// name. Atoms are taken in the order given, so a residue split by other
/// atoms becomes two residues, which is valid and visible rather than
/// silently reordered.
#[allow(clippy::too_many_arguments)]
#[pyfunction]
#[pyo3(name = "_from_columns")]
pub fn from_columns(
    positions: PyReadonlyArray3<'_, f32>,
    element: PyReadonlyArray1<'_, u8>,
    name: PyReadonlyArray1<'_, PyFixedString<4>>,
    residue_names: Vec<String>,
    residue_name_index: PyReadonlyArray1<'_, u32>,
    residue_ids: PyReadonlyArray1<'_, i32>,
    chain_names: Vec<String>,
    chain_index: PyReadonlyArray1<'_, u32>,
    auth_chain_names: Vec<String>,
    auth_chain_index: PyReadonlyArray1<'_, u32>,
    b_factor: PyReadonlyArray1<'_, f32>,
    occupancy: PyReadonlyArray1<'_, f32>,
    hetero: PyReadonlyArray1<'_, bool>,
    id: &str,
    title: &str,
) -> PyResult<Structure> {
    let c = Columns {
        positions,
        element,
        name,
        residue_names,
        residue_name_index,
        residue_ids,
        chain_names,
        chain_index,
        auth_chain_names,
        auth_chain_index,
        b_factor,
        occupancy,
        hetero,
    };
    build(&c, id, title)
}

fn build(c: &Columns<'_>, id: &str, title: &str) -> PyResult<Structure> {
    let shape = c.positions.shape();
    let (frames, atoms) = (shape[0], shape[1]);
    if shape[2] != 3 {
        return Err(PyValueError::new_err(format!(
            "positions must be (frames, atoms, 3), got {shape:?}"
        )));
    }
    if frames == 0 {
        return Err(PyValueError::new_err("positions has no frames"));
    }
    let check = |what: &str, len: usize| {
        if len == atoms {
            Ok(())
        } else {
            Err(PyValueError::new_err(format!(
                "{what} has {len} entries, positions has {atoms} atoms"
            )))
        }
    };
    check("elements", c.element.len())?;
    check("names", c.name.len())?;
    check("residue_names", c.residue_name_index.len())?;
    check("residue_ids", c.residue_ids.len())?;
    check("chain_ids", c.chain_index.len())?;
    check("auth_chain_ids", c.auth_chain_index.len())?;
    check("b_factors", c.b_factor.len())?;
    check("occupancies", c.occupancy.len())?;
    check("hetero", c.hetero.len())?;

    let positions = c.positions.as_array();
    let element = c.element.as_array();
    let name = c.name.as_array();
    let res_idx = c.residue_name_index.as_array();
    let res_id = c.residue_ids.as_array();
    let chain_idx = c.chain_index.as_array();
    let auth_idx = c.auth_chain_index.as_array();
    let b_factor = c.b_factor.as_array();
    let occupancy = c.occupancy.as_array();
    let hetero = c.hetero.as_array();

    let mut builder = TopologyBuilder::with_capacity(atoms);
    for a in 0..atoms {
        let comp = c
            .residue_names
            .get(res_idx[a] as usize)
            .ok_or_else(|| PyValueError::new_err("residue name index out of range"))?;
        let asym = c
            .chain_names
            .get(chain_idx[a] as usize)
            .ok_or_else(|| PyValueError::new_err("chain index out of range"))?;
        let auth_asym = c
            .auth_chain_names
            .get(auth_idx[a] as usize)
            .ok_or_else(|| PyValueError::new_err("author chain index out of range"))?;
        let mut padded = name[a].0;
        for b in padded.iter_mut() {
            if *b == 0 {
                *b = b' ';
            }
        }
        let z = element[a];
        let element = Element::from_atomic_number(z).unwrap_or(Element::UNKNOWN);
        let row = AtomRow {
            element,
            name: padded,
            serial: a as u32 + 1,
            alt_loc: 0,
            comp,
            asym,
            auth_asym,
            seq_id: res_id[a],
            auth_seq_id: res_id[a],
            ins_code: 0,
            entity: 0,
            position: Vec3::new(
                positions[[0, a, 0]],
                positions[[0, a, 1]],
                positions[[0, a, 2]],
            ),
            occupancy: occupancy[a],
            b_factor: b_factor[a],
            charge: 0,
            hetero: hetero[a],
        };
        builder.push(&row);
    }
    builder.topology.id = id.to_owned();
    builder.topology.title = title.to_owned();
    let extra: Vec<Vec<Vec3>> = (1..frames)
        .map(|f| {
            (0..atoms)
                .map(|a| {
                    Vec3::new(
                        positions[[f, a, 0]],
                        positions[[f, a, 1]],
                        positions[[f, a, 2]],
                    )
                })
                .collect()
        })
        .collect();
    let inner = builder
        .finish_with_frames(extra)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    Ok(Structure::new(inner))
}

/// Atomic numbers guessed from `|S4` atom names (`CA` -> 6, `FE` -> 26),
/// the rule PDB files without an element column are read with.
#[pyfunction]
#[pyo3(name = "_guess_elements")]
pub fn guess_elements<'py>(
    py: Python<'py>,
    names: PyReadonlyArray1<'py, PyFixedString<4>>,
) -> Bound<'py, PyArray1<u8>> {
    PyArray1::from_iter(
        py,
        names
            .as_array()
            .iter()
            .map(|n| Element::from_atom_name(&n.0).atomic_number()),
    )
}
