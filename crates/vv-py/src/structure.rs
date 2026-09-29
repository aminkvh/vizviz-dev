//! `vizviz.Structure`: a parsed structure whose per-atom columns are
//! exposed to Python as read-only, zero-copy NumPy views.
//!
//! Why this is sound: `vv_core::Structure` owns its columns through `Arc`s
//! and never mutates them after construction (the parsers call
//! `Arc::get_mut` only on a structure nothing else has seen yet). Each
//! view is made with `PyArray::borrow_from_array`, which records the Python
//! `Structure` object as the array's `base`, so the array keeps the Rust
//! data alive for as long as it exists — even after the user drops their
//! last reference to the `Structure`. Views are marked non-writeable so
//! Python code can't reach into memory the renderer and other views read.
//!
//! Per-residue and per-chain tables are *copies*: they are records, not
//! columns, and 1000x smaller than the atom arrays, so a copy costs
//! milliseconds even on a 10M-atom structure.

use std::path::PathBuf;
use std::sync::Arc;

use numpy::ndarray::{Array2, ArrayView1, ArrayView2};
use numpy::{
    Element as NpElement, PyArray1, PyArray2, PyArrayMethods, PyFixedString, PyReadonlyArray2,
};
use pyo3::exceptions::{PyIOError, PyIndexError, PyValueError};
use pyo3::prelude::*;
use vv_core::analysis;
use vv_core::fixedbitset::FixedBitSet;
use vv_core::glam::Vec3;
use vv_core::{CoordSet, ExplicitBondKind, InternId, Structure as CoreStructure, StructureError};

#[pyclass(name = "Structure", module = "vizviz", frozen)]
pub struct Structure {
    pub(crate) inner: CoreStructure,
}

impl Structure {
    pub(crate) fn new(inner: CoreStructure) -> Self {
        Self { inner }
    }
}

/// A read-only NumPy view of `data`, kept alive by `owner`.
fn view1<'py, T: NpElement>(owner: &Bound<'py, Structure>, data: &[T]) -> Bound<'py, PyArray1<T>> {
    let view = ArrayView1::from(data);
    // SAFETY: `data` lives inside `owner`'s Arc-owned, never-mutated
    // columns (module docs); `owner` becomes the array's base object, so
    // the memory outlives every view.
    let array = unsafe { PyArray1::borrow_from_array(&view, owner.clone().into_any()) };
    array.readwrite().make_nonwriteable();
    array
}

/// Same as [`view1`] for a flat column viewed as `(len / cols, cols)`.
fn view2<'py, T: NpElement>(
    owner: &Bound<'py, Structure>,
    data: &[T],
    cols: usize,
) -> Bound<'py, PyArray2<T>> {
    let view = ArrayView2::from_shape((data.len() / cols, cols), data)
        .expect("flat column length is a multiple of cols");
    // SAFETY: as in `view1`.
    let array = unsafe { PyArray2::borrow_from_array(&view, owner.clone().into_any()) };
    array.readwrite().make_nonwriteable();
    array
}

/// An owned `(len / cols, cols)` array (for the small record tables).
fn owned2<T: NpElement>(py: Python<'_>, data: Vec<T>, cols: usize) -> Bound<'_, PyArray2<T>> {
    let rows = data.len() / cols;
    PyArray2::from_owned_array(
        py,
        Array2::from_shape_vec((rows, cols), data).expect("rows * cols entries"),
    )
}

/// Index tuples of width `N` from a `(k, N)` integer array or a sequence
/// of sequences, bounds-checked against `atom_count`.
fn index_tuples<const N: usize>(
    obj: &Bound<'_, PyAny>,
    atom_count: usize,
    what: &str,
) -> PyResult<Vec<[u32; N]>> {
    let rows: Vec<Vec<i64>> = if let Ok(a) = obj.extract::<PyReadonlyArray2<i64>>() {
        a.as_array()
            .rows()
            .into_iter()
            .map(|r| r.to_vec())
            .collect()
    } else if let Ok(a) = obj.extract::<PyReadonlyArray2<u32>>() {
        a.as_array()
            .rows()
            .into_iter()
            .map(|r| r.iter().map(|&v| v as i64).collect())
            .collect()
    } else if let Ok(a) = obj.extract::<PyReadonlyArray2<i32>>() {
        a.as_array()
            .rows()
            .into_iter()
            .map(|r| r.iter().map(|&v| v as i64).collect())
            .collect()
    } else {
        obj.extract::<Vec<Vec<i64>>>().map_err(|_| {
            PyValueError::new_err(format!(
                "{what} must be a (k, {N}) integer array or a list of {N}-tuples"
            ))
        })?
    };
    rows.into_iter()
        .map(|row| {
            if row.len() != N {
                return Err(PyValueError::new_err(format!(
                    "{what} rows must have {N} atom indices, got {}",
                    row.len()
                )));
            }
            let mut out = [0u32; N];
            for (o, &i) in out.iter_mut().zip(&row) {
                if i < 0 || i as usize >= atom_count {
                    return Err(PyIndexError::new_err(format!(
                        "atom index {i} out of range for {atom_count} atoms"
                    )));
                }
                *o = i as u32;
            }
            Ok(out)
        })
        .collect()
}

/// A group of atoms: a selection expression, or a sequence / array of
/// atom indices. Deduplicated and sorted.
fn atom_group(structure: &CoreStructure, obj: &Bound<'_, PyAny>, what: &str) -> PyResult<Vec<u32>> {
    let mut atoms: Vec<u32> = if let Ok(expr) = obj.extract::<String>() {
        vv_core::select(&structure.topology, structure.frame(0).positions(), &expr)
            .map_err(crate::select_err)?
            .ones()
            .map(|i| i as u32)
            .collect()
    } else {
        let n = structure.atom_count();
        let indices: Vec<i64> = obj.extract().map_err(|_| {
            PyValueError::new_err(format!(
                "{what} must be a selection expression or a sequence of atom indices"
            ))
        })?;
        indices
            .into_iter()
            .map(|i| {
                if i < 0 || i as usize >= n {
                    Err(PyIndexError::new_err(format!(
                        "atom index {i} out of range for {n} atoms"
                    )))
                } else {
                    Ok(i as u32)
                }
            })
            .collect::<PyResult<_>>()?
    };
    atoms.sort_unstable();
    atoms.dedup();
    Ok(atoms)
}

/// Frame indices to work on: all of them when `frames` is `None`.
fn frame_list(structure: &CoreStructure, frames: Option<Vec<usize>>) -> PyResult<Vec<usize>> {
    let n = structure.frame_count();
    match frames {
        None => Ok((0..n).collect()),
        Some(list) => {
            for &f in &list {
                if f >= n {
                    return Err(PyIndexError::new_err(format!(
                        "frame {f} out of range: structure has {n} frame(s)"
                    )));
                }
            }
            Ok(list)
        }
    }
}

/// Frames `frames` of `structure`, read if they are streamed; an
/// unreadable frame is an error, never a stand-in.
fn read_frames(
    structure: &CoreStructure,
    frames: &[usize],
) -> Result<Vec<Arc<CoordSet>>, StructureError> {
    frames.iter().map(|&f| structure.try_frame(f)).collect()
}

fn read_error(e: StructureError) -> PyErr {
    PyIOError::new_err(e.to_string())
}

fn contact_arrays<'py>(
    py: Python<'py>,
    contacts: &[analysis::Contact],
) -> (Bound<'py, PyArray2<u32>>, Bound<'py, PyArray1<f32>>) {
    let flat: Vec<u32> = contacts.iter().flat_map(|c| [c.a, c.b]).collect();
    let pairs = owned2(py, flat, 2);
    let distances = PyArray1::from_iter(py, contacts.iter().map(|c| c.distance));
    (pairs, distances)
}

fn check_cutoff(cutoff: f32) -> PyResult<()> {
    if !(cutoff.is_finite() && cutoff > 0.0) {
        return Err(PyValueError::new_err("cutoff must be a positive number"));
    }
    Ok(())
}

/// `(rows, cols)` from a flat row-major vector; `cols == 0` gives an
/// empty `(rows, 0)` array.
fn owned2_shape<T: NpElement>(
    py: Python<'_>,
    data: Vec<T>,
    rows: usize,
    cols: usize,
) -> Bound<'_, PyArray2<T>> {
    PyArray2::from_owned_array(
        py,
        Array2::from_shape_vec((rows, cols), data).expect("rows * cols entries"),
    )
}

fn bond_kind_code(kind: ExplicitBondKind) -> u8 {
    match kind {
        ExplicitBondKind::Covalent => 0,
        ExplicitBondKind::Disulfide => 1,
        ExplicitBondKind::Metal => 2,
        ExplicitBondKind::Other => 3,
    }
}

#[pymethods]
impl Structure {
    // ---- scalars -------------------------------------------------------

    #[getter]
    fn atom_count(&self) -> usize {
        self.inner.atom_count()
    }

    #[getter]
    fn residue_count(&self) -> usize {
        self.inner.topology.residue_count()
    }

    #[getter]
    fn chain_count(&self) -> usize {
        self.inner.topology.chain_count()
    }

    /// Number of coordinate sets (models / trajectory frames).
    #[getter]
    fn frame_count(&self) -> usize {
        self.inner.frame_count()
    }

    /// The entry id from the file (`_entry.id` / PDB HEADER), or `""`.
    #[getter]
    fn id(&self) -> String {
        self.inner.topology.id.clone()
    }

    #[getter]
    fn title(&self) -> String {
        self.inner.topology.title.clone()
    }

    fn __len__(&self) -> usize {
        self.inner.atom_count()
    }

    fn __repr__(&self) -> String {
        let t = &self.inner.topology;
        format!(
            "Structure(id={:?}, atoms={}, residues={}, chains={}, frames={})",
            t.id,
            self.inner.atom_count(),
            t.residue_count(),
            t.chain_count(),
            self.inner.frame_count()
        )
    }

    // ---- per-atom columns: zero-copy, read-only ------------------------

    /// `(atom_count, 3)` float32 coordinates of frame 0, in the file's
    /// units (Angstroms for mmCIF/PDB).
    #[getter]
    fn positions<'py>(this: Bound<'py, Self>) -> Bound<'py, PyArray2<f32>> {
        let s = this.get();
        let frame = s.inner.resident_frame(0).expect("frame 0 stays in memory");
        let flat: &[f32] = bytemuck::cast_slice(frame.positions());
        view2(&this, flat, 3)
    }

    /// `(atom_count, 3)` float32 coordinates of frame `index`.
    fn frame<'py>(this: Bound<'py, Self>, index: usize) -> PyResult<Bound<'py, PyArray2<f32>>> {
        let s = this.get();
        if index >= s.inner.frame_count() {
            return Err(PyIndexError::new_err(format!(
                "frame {index} out of range: structure has {} frame(s)",
                s.inner.frame_count()
            )));
        }
        // A view where the frame stays in memory, a copy of a streamed one.
        if let Some(frame) = s.inner.resident_frame(index) {
            let flat: &[f32] = bytemuck::cast_slice(frame.positions());
            return Ok(view2(&this, flat, 3));
        }
        let frame = s.inner.try_frame(index).map_err(read_error)?;
        let flat: &[f32] = bytemuck::cast_slice(frame.positions());
        Ok(owned2(this.py(), flat.to_vec(), 3))
    }

    /// `(atom_count,)` uint8 atomic numbers; 0 means unknown.
    #[getter]
    fn element<'py>(this: Bound<'py, Self>) -> Bound<'py, PyArray1<u8>> {
        let data: &[u8] = bytemuck::cast_slice(&this.get().inner.topology.element);
        view1(&this, data)
    }

    /// `(atom_count,)` `S4` atom names as stored in the file, space padded
    /// (`b" CA "`). `numpy.char.strip(s.name)` gives clean bytes.
    #[getter]
    fn name<'py>(this: Bound<'py, Self>) -> Bound<'py, PyArray1<PyFixedString<4>>> {
        let names: &[[u8; 4]] = &this.get().inner.topology.name;
        // SAFETY: `PyFixedString<4>` is `#[repr(transparent)]` over `[u8; 4]`.
        let names: &[PyFixedString<4>] =
            unsafe { std::slice::from_raw_parts(names.as_ptr().cast(), names.len()) };
        view1(&this, names)
    }

    /// `(atom_count,)` uint32 atom serial numbers from the file.
    #[getter]
    fn serial<'py>(this: Bound<'py, Self>) -> Bound<'py, PyArray1<u32>> {
        view1(&this, &this.get().inner.topology.serial)
    }

    /// `(atom_count,)` float32 B-factors (0 when the file has none).
    #[getter]
    fn b_factor<'py>(this: Bound<'py, Self>) -> Bound<'py, PyArray1<f32>> {
        view1(&this, &this.get().inner.topology.b_factor)
    }

    /// `(atom_count,)` float32 occupancies.
    #[getter]
    fn occupancy<'py>(this: Bound<'py, Self>) -> Bound<'py, PyArray1<f32>> {
        view1(&this, &this.get().inner.topology.occupancy)
    }

    /// `(atom_count,)` uint8 alternate-location ids as ASCII; 0 = none.
    #[getter]
    fn alt_loc<'py>(this: Bound<'py, Self>) -> Bound<'py, PyArray1<u8>> {
        view1(&this, &this.get().inner.topology.alt_loc)
    }

    /// `(atom_count,)` int8 formal charges.
    #[getter]
    fn charge<'py>(this: Bound<'py, Self>) -> Bound<'py, PyArray1<i8>> {
        view1(&this, &this.get().inner.topology.charge)
    }

    /// `(atom_count,)` uint8 bit flags; bit 0 set = HETATM record.
    #[getter]
    fn flags<'py>(this: Bound<'py, Self>) -> Bound<'py, PyArray1<u8>> {
        view1(&this, &this.get().inner.topology.flags)
    }

    /// `(atom_count,)` uint32 index of each atom's residue.
    #[getter]
    fn residue_index<'py>(this: Bound<'py, Self>) -> Bound<'py, PyArray1<u32>> {
        view1(&this, &this.get().inner.topology.residue_index)
    }

    // ---- per-residue and per-chain tables: small copies ----------------

    /// `(residue_count, 2)` uint32 `[start, end)` atom ranges.
    #[getter]
    fn residue_atoms<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<u32>> {
        let flat = self
            .inner
            .topology
            .residues
            .iter()
            .flat_map(|r| [r.atoms.start, r.atoms.end])
            .collect();
        owned2(py, flat, 2)
    }

    /// `(residue_count,)` uint32 chain index of each residue.
    #[getter]
    fn residue_chain<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<u32>> {
        PyArray1::from_iter(py, self.inner.topology.residues.iter().map(|r| r.chain))
    }

    /// `(residue_count,)` int32 `label_seq_id` (mmCIF) / residue number (PDB).
    #[getter]
    fn residue_seq_id<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<i32>> {
        PyArray1::from_iter(py, self.inner.topology.residues.iter().map(|r| r.seq_id))
    }

    /// `(residue_count,)` int32 author residue numbers.
    #[getter]
    fn residue_auth_seq_id<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<i32>> {
        PyArray1::from_iter(
            py,
            self.inner.topology.residues.iter().map(|r| r.auth_seq_id),
        )
    }

    /// `(residue_count,)` uint8 insertion codes as ASCII; 0 = none.
    #[getter]
    fn residue_ins_code<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<u8>> {
        PyArray1::from_iter(py, self.inner.topology.residues.iter().map(|r| r.ins_code))
    }

    /// `(residue_count,)` uint8 secondary structure from the file's
    /// records: 0 unknown, 1 coil, 2 helix, 3 strand.
    #[getter]
    fn residue_ss<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<u8>> {
        PyArray1::from_iter(py, self.inner.topology.residues.iter().map(|r| r.ss as u8))
    }

    /// `(residue_count,)` uint32 index into `names` of each residue's
    /// type (e.g. `"ALA"`). `residue_names` is the string form.
    #[getter]
    fn residue_comp<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<u32>> {
        PyArray1::from_iter(py, self.inner.topology.residues.iter().map(|r| r.comp.0))
    }

    /// Residue type names, one per residue (`["VAL", "LEU", ...]`).
    #[getter]
    fn residue_names(&self) -> Vec<String> {
        let t = &self.inner.topology;
        (0..t.residue_count())
            .map(|i| t.residue_name(i).to_owned())
            .collect()
    }

    /// `(chain_count, 2)` uint32 `[start, end)` residue ranges.
    #[getter]
    fn chain_residues<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<u32>> {
        let flat = self
            .inner
            .topology
            .chains
            .iter()
            .flat_map(|c| [c.residues.start, c.residues.end])
            .collect();
        owned2(py, flat, 2)
    }

    /// `(chain_count,)` uint16 entity ids (mmCIF `label_entity_id`).
    #[getter]
    fn chain_entity<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<u16>> {
        PyArray1::from_iter(py, self.inner.topology.chains.iter().map(|c| c.entity))
    }

    /// Chain names (`label_asym_id`), one per chain.
    #[getter]
    fn chain_names(&self) -> Vec<String> {
        let t = &self.inner.topology;
        (0..t.chain_count())
            .map(|i| t.chain_name(i).to_owned())
            .collect()
    }

    /// Author chain names (`auth_asym_id`), one per chain.
    #[getter]
    fn chain_auth_names(&self) -> Vec<String> {
        let t = &self.inner.topology;
        t.chains
            .iter()
            .map(|c| t.names.get(c.auth_asym).to_owned())
            .collect()
    }

    /// The interned string table that `residue_comp` indexes into.
    #[getter]
    fn names(&self) -> Vec<String> {
        let names = &self.inner.topology.names;
        (0..names.len())
            .map(|i| names.get(InternId(i as u32)).to_owned())
            .collect()
    }

    // ---- annotations ---------------------------------------------------

    /// The file's header annotations as `{category: [row, ...]}`, each row
    /// a `{item: value}` dict, in the file's own mmCIF vocabulary (PDB
    /// headers are mapped onto the same names). Empty values are omitted.
    #[getter]
    fn annotations<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, pyo3::types::PyDict>> {
        use pyo3::types::{PyDict, PyList};
        let out = PyDict::new(py);
        for category in &self.inner.topology.annotations.categories {
            let rows = PyList::empty(py);
            for row in &category.rows {
                let d = PyDict::new(py);
                for (item, value) in category.items.iter().zip(row) {
                    if !value.is_empty() {
                        d.set_item(item, value)?;
                    }
                }
                rows.append(d)?;
            }
            out.set_item(&category.name, rows)?;
        }
        Ok(out)
    }

    /// The usual summary fields, ready to print: title, method,
    /// resolution (float, Angstroms), deposition_date, organism, keywords,
    /// citation_title, doi, uniprot (list), entities (list of (id,
    /// description)). Missing fields are `None`.
    #[getter]
    fn info<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, pyo3::types::PyDict>> {
        use pyo3::types::PyDict;
        let a = &self.inner.topology.annotations;
        let d = PyDict::new(py);
        d.set_item("id", &self.inner.topology.id)?;
        d.set_item("title", a.title())?;
        d.set_item("method", a.method())?;
        d.set_item("resolution", a.resolution())?;
        d.set_item("deposition_date", a.deposition_date())?;
        d.set_item("organism", a.organism())?;
        d.set_item("keywords", a.keywords())?;
        d.set_item("citation_title", a.citation_title())?;
        d.set_item("doi", a.doi())?;
        d.set_item("uniprot", a.uniprot_accessions())?;
        d.set_item("entities", a.entities())?;
        Ok(d)
    }

    // ---- bonds ---------------------------------------------------------

    /// `(bond_count, 2)` uint32 atom pairs (`a < b`) from geometric bond
    /// perception on frame 0 (see docs/FORMATS.md). Computed on each call.
    fn bonds<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<u32>> {
        let table = py.detach(|| {
            vv_core::bonds::perceive(&self.inner.topology, self.inner.frame(0).positions())
        });
        let flat = table.pairs.iter().flat_map(|p| [p[0], p[1]]).collect();
        owned2(py, flat, 2)
    }

    /// `(k, 2)` uint32 atom pairs the file listed explicitly (mmCIF
    /// `_struct_conn`, PDB CONECT).
    #[getter]
    fn explicit_bonds<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<u32>> {
        let flat = self
            .inner
            .topology
            .explicit_bonds
            .iter()
            .flat_map(|b| b.atoms)
            .collect();
        owned2(py, flat, 2)
    }

    /// `(k,)` uint8 kind of each explicit bond: 0 covalent, 1 disulfide,
    /// 2 metal coordination, 3 other.
    #[getter]
    fn explicit_bond_kinds<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<u8>> {
        PyArray1::from_iter(
            py,
            self.inner
                .topology
                .explicit_bonds
                .iter()
                .map(|b| bond_kind_code(b.kind)),
        )
    }

    // ---- scalar lookups and analysis ----------------------------------

    /// Trimmed name of one atom (`"CA"`).
    fn atom_name(&self, atom: usize) -> PyResult<String> {
        self.check_atom(atom)?;
        Ok(self.inner.topology.atom_name(atom).to_owned())
    }

    fn residue_name(&self, residue: usize) -> PyResult<String> {
        if residue >= self.inner.topology.residue_count() {
            return Err(PyIndexError::new_err(format!(
                "residue {residue} out of range for {} residues",
                self.inner.topology.residue_count()
            )));
        }
        Ok(self.inner.topology.residue_name(residue).to_owned())
    }

    fn chain_name(&self, chain: usize) -> PyResult<String> {
        if chain >= self.inner.topology.chain_count() {
            return Err(PyIndexError::new_err(format!(
                "chain {chain} out of range for {} chains",
                self.inner.topology.chain_count()
            )));
        }
        Ok(self.inner.topology.chain_name(chain).to_owned())
    }

    /// Indices of the atoms matching a selection expression such as
    /// `"chain A and name CA"` or `"within 5 of resname HEM"`; see
    /// docs/SELECTION.md for the grammar. A bad expression raises
    /// `ValueError` naming the offending characters.
    fn select<'py>(&self, py: Python<'py>, expr: &str) -> PyResult<Bound<'py, PyArray1<u32>>> {
        let structure = &self.inner;
        let bits = py
            .detach(|| vv_core::select(&structure.topology, structure.frame(0).positions(), expr))
            .map_err(crate::select_err)?;
        Ok(PyArray1::from_iter(py, bits.ones().map(|i| i as u32)))
    }

    /// Angle at `b` between `a-b` and `c-b`, in degrees.
    #[pyo3(signature = (a, b, c, frame = 0))]
    fn angle(&self, a: usize, b: usize, c: usize, frame: usize) -> PyResult<f32> {
        for i in [a, b, c] {
            self.check_atom(i)?;
        }
        let frame = frame_list(&self.inner, Some(vec![frame]))?[0];
        Ok(analysis::angle(
            self.inner.frame(frame).positions(),
            a,
            b,
            c,
        ))
    }

    /// Dihedral (torsion) of `a-b-c-d` about `b-c`, in degrees, IUPAC
    /// sign (the same convention as MDAnalysis).
    #[pyo3(signature = (a, b, c, d, frame = 0))]
    fn dihedral(&self, a: usize, b: usize, c: usize, d: usize, frame: usize) -> PyResult<f32> {
        for i in [a, b, c, d] {
            self.check_atom(i)?;
        }
        let frame = frame_list(&self.inner, Some(vec![frame]))?[0];
        Ok(analysis::dihedral(
            self.inner.frame(frame).positions(),
            a,
            b,
            c,
            d,
        ))
    }

    /// Distances for many atom pairs at once: `pairs` is `(k, 2)`; the
    /// result is `(frames, k)` float32, one row per frame (all frames by
    /// default). Frames run in parallel with the GIL released.
    #[pyo3(signature = (pairs, frames = None))]
    fn distances<'py>(
        &self,
        py: Python<'py>,
        pairs: &Bound<'py, PyAny>,
        frames: Option<Vec<usize>>,
    ) -> PyResult<Bound<'py, PyArray2<f32>>> {
        let pairs = index_tuples::<2>(pairs, self.inner.atom_count(), "pairs")?;
        let frames = frame_list(&self.inner, frames)?;
        let structure = &self.inner;
        let out = py
            .detach(|| -> Result<_, StructureError> {
                let held = read_frames(structure, &frames)?;
                let slices: Vec<&[Vec3]> = held.iter().map(|f| f.positions()).collect();
                let mut out = vec![0.0f32; slices.len() * pairs.len()];
                analysis::distances_frames(&slices, &pairs, &mut out);
                Ok(out)
            })
            .map_err(read_error)?;
        Ok(owned2_shape(py, out, frames.len(), pairs.len()))
    }

    /// Angles for many atom triples: `triples` is `(k, 3)`; result is
    /// `(frames, k)` float32 degrees.
    #[pyo3(signature = (triples, frames = None))]
    fn angles<'py>(
        &self,
        py: Python<'py>,
        triples: &Bound<'py, PyAny>,
        frames: Option<Vec<usize>>,
    ) -> PyResult<Bound<'py, PyArray2<f32>>> {
        let triples = index_tuples::<3>(triples, self.inner.atom_count(), "triples")?;
        let frames = frame_list(&self.inner, frames)?;
        let structure = &self.inner;
        let out = py
            .detach(|| -> Result<_, StructureError> {
                let held = read_frames(structure, &frames)?;
                let slices: Vec<&[Vec3]> = held.iter().map(|f| f.positions()).collect();
                let mut out = vec![0.0f32; slices.len() * triples.len()];
                analysis::angles_frames(&slices, &triples, &mut out);
                Ok(out)
            })
            .map_err(read_error)?;
        Ok(owned2_shape(py, out, frames.len(), triples.len()))
    }

    /// Dihedrals for many atom quadruples: `quads` is `(k, 4)`; result is
    /// `(frames, k)` float32 degrees, IUPAC sign.
    #[pyo3(signature = (quads, frames = None))]
    fn dihedrals<'py>(
        &self,
        py: Python<'py>,
        quads: &Bound<'py, PyAny>,
        frames: Option<Vec<usize>>,
    ) -> PyResult<Bound<'py, PyArray2<f32>>> {
        let quads = index_tuples::<4>(quads, self.inner.atom_count(), "quads")?;
        let frames = frame_list(&self.inner, frames)?;
        let structure = &self.inner;
        let out = py
            .detach(|| -> Result<_, StructureError> {
                let held = read_frames(structure, &frames)?;
                let slices: Vec<&[Vec3]> = held.iter().map(|f| f.positions()).collect();
                let mut out = vec![0.0f32; slices.len() * quads.len()];
                analysis::dihedrals_frames(&slices, &quads, &mut out);
                Ok(out)
            })
            .map_err(read_error)?;
        Ok(owned2_shape(py, out, frames.len(), quads.len()))
    }

    /// Atom pairs from group `a` to group `b` within `cutoff` Angstroms in
    /// one frame: `(pairs (k, 2) uint32, distances (k,) float32)`, sorted.
    /// Groups are selection expressions or index sequences.
    #[pyo3(signature = (a, b, cutoff, frame = 0))]
    fn contacts<'py>(
        &self,
        py: Python<'py>,
        a: &Bound<'py, PyAny>,
        b: &Bound<'py, PyAny>,
        cutoff: f32,
        frame: usize,
    ) -> PyResult<(Bound<'py, PyArray2<u32>>, Bound<'py, PyArray1<f32>>)> {
        check_cutoff(cutoff)?;
        let (ga, gb) = (
            atom_group(&self.inner, a, "a")?,
            atom_group(&self.inner, b, "b")?,
        );
        let frame = frame_list(&self.inner, Some(vec![frame]))?[0];
        let structure = &self.inner;
        let found = py.detach(|| {
            let mut out = Vec::new();
            analysis::contacts_into(
                structure.frame(frame).positions(),
                &ga,
                &gb,
                cutoff,
                &mut out,
            );
            out
        });
        Ok(contact_arrays(py, &found))
    }

    /// Number of `a`/`b` contacts within `cutoff` in each frame (all
    /// frames by default): `(frames,)` uint32. The time-series form; it
    /// never materializes the pairs.
    #[pyo3(signature = (a, b, cutoff, frames = None))]
    fn contact_counts<'py>(
        &self,
        py: Python<'py>,
        a: &Bound<'py, PyAny>,
        b: &Bound<'py, PyAny>,
        cutoff: f32,
        frames: Option<Vec<usize>>,
    ) -> PyResult<Bound<'py, PyArray1<u32>>> {
        check_cutoff(cutoff)?;
        let (ga, gb) = (
            atom_group(&self.inner, a, "a")?,
            atom_group(&self.inner, b, "b")?,
        );
        let frames = frame_list(&self.inner, frames)?;
        let structure = &self.inner;
        let counts = py
            .detach(|| -> Result<_, StructureError> {
                let held = read_frames(structure, &frames)?;
                let slices: Vec<&[Vec3]> = held.iter().map(|f| f.positions()).collect();
                Ok(analysis::contact_counts_frames(&slices, &ga, &gb, cutoff))
            })
            .map_err(read_error)?;
        Ok(PyArray1::from_vec(py, counts))
    }

    /// Every pair within one group closer than `cutoff`, each once with
    /// `a < b`: `(pairs (k, 2) uint32, distances (k,) float32)`. `atoms`
    /// defaults to the whole structure.
    #[pyo3(signature = (cutoff, atoms = None, frame = 0))]
    fn neighbors<'py>(
        &self,
        py: Python<'py>,
        cutoff: f32,
        atoms: Option<&Bound<'py, PyAny>>,
        frame: usize,
    ) -> PyResult<(Bound<'py, PyArray2<u32>>, Bound<'py, PyArray1<f32>>)> {
        check_cutoff(cutoff)?;
        let group = match atoms {
            Some(obj) => atom_group(&self.inner, obj, "atoms")?,
            None => (0..self.inner.atom_count() as u32).collect(),
        };
        let frame = frame_list(&self.inner, Some(vec![frame]))?[0];
        let structure = &self.inner;
        let found = py.detach(|| {
            let mut out = Vec::new();
            analysis::neighbor_pairs_into(
                structure.frame(frame).positions(),
                &group,
                cutoff,
                &mut out,
            );
            out
        });
        Ok(contact_arrays(py, &found))
    }

    /// The distinct residue pairs behind atom `pairs` (`(k, 2)`), as
    /// `(m, 2)` uint32 residue indices, sorted; same-residue pairs dropped.
    fn residue_pairs<'py>(
        &self,
        py: Python<'py>,
        pairs: &Bound<'py, PyAny>,
    ) -> PyResult<Bound<'py, PyArray2<u32>>> {
        let pairs = index_tuples::<2>(pairs, self.inner.atom_count(), "pairs")?;
        let contacts: Vec<analysis::Contact> = pairs
            .iter()
            .map(|[a, b]| analysis::Contact {
                a: *a,
                b: *b,
                distance: 0.0,
            })
            .collect();
        let residues = analysis::residue_pairs(&self.inner.topology, &contacts);
        let flat: Vec<u32> = residues.iter().flat_map(|p| *p).collect();
        Ok(owned2_shape(py, flat, residues.len(), 2))
    }

    /// Distance between atoms `a` and `b` in `frame`, in the file's units.
    #[pyo3(signature = (a, b, frame = 0))]
    fn distance(&self, a: usize, b: usize, frame: usize) -> PyResult<f32> {
        self.check_atom(a)?;
        self.check_atom(b)?;
        if frame >= self.inner.frame_count() {
            return Err(PyIndexError::new_err(format!(
                "frame {frame} out of range: structure has {} frame(s)",
                self.inner.frame_count()
            )));
        }
        Ok(vv_core::distance(self.inner.frame(frame).positions(), a, b))
    }

    /// Writes the structure to `path`: format from its extension
    /// (`.pdb`/`.ent`, `.cif`/`.mmcif`/`.pdbx`, `.xyz`, `.pqr`, `.gro`,
    /// any of them optionally `.gz`). `selection` is a selection
    /// expression or a sequence of atom indices (default: every atom);
    /// `frames` are coordinate-set indices to write, `MODEL`/`ENDMDL`- or
    /// `pdbx_PDB_model_num`-wrapped when there is more than one (default:
    /// every frame). Returns one warning per thing the format couldn't
    /// represent exactly (e.g. a PDB chain id longer than one character).
    #[pyo3(signature = (path, selection = None, frames = None))]
    fn save(
        &self,
        py: Python<'_>,
        path: PathBuf,
        selection: Option<&Bound<'_, PyAny>>,
        frames: Option<Vec<usize>>,
    ) -> PyResult<Vec<String>> {
        let mask = match selection {
            None => None,
            Some(obj) => {
                let atoms = atom_group(&self.inner, obj, "selection")?;
                if atoms.is_empty() {
                    return Err(PyValueError::new_err("selection matches no atoms"));
                }
                let mut bits = FixedBitSet::with_capacity(self.inner.atom_count());
                bits.extend(atoms.into_iter().map(|a| a as usize));
                Some(bits)
            }
        };
        let frame_indices = frame_list(&self.inner, frames)?;
        let structure = &self.inner;
        py.detach(|| {
            let opts = vv_io::SaveOptions {
                atoms: mask.as_ref(),
                frames: &frame_indices,
            };
            vv_io::save(structure, &path, &opts)
        })
        .map_err(crate::save_err)
    }

    /// Renders frame 0 with the CPU backend (spacefill only) and returns
    /// an `(height, width, 4)` uint8 RGBA image. Color by a named scheme,
    /// by `values` (one number per atom, e.g. a SASA from another package,
    /// on the blue-white-red ramp over `range` or its own min..max), or by
    /// explicit `colors` (an `(atoms, 3)` or `(atoms, 4)` uint8 array).
    /// See docs/PYTHON.md.
    #[cfg(feature = "render")]
    #[pyo3(signature = (width, height, *, coloring = "element", values = None, colors = None,
                        range = None, style = "dark_presentation", yaw = 0.0, pitch = 0.0,
                        zoom = 1.0, background = None))]
    #[allow(clippy::too_many_arguments)]
    fn render<'py>(
        &self,
        py: Python<'py>,
        width: u32,
        height: u32,
        coloring: &str,
        values: Option<&Bound<'py, PyAny>>,
        colors: Option<&Bound<'py, PyAny>>,
        range: Option<(f32, f32)>,
        style: &str,
        yaw: f32,
        pitch: f32,
        zoom: f32,
        background: Option<(u8, u8, u8, u8)>,
    ) -> PyResult<Bound<'py, numpy::PyArray3<u8>>> {
        let options = crate::render::RenderOptions {
            width,
            height,
            style,
            yaw,
            pitch,
            zoom,
            background,
        };
        let colors =
            crate::render::colors_from_python(&self.inner, coloring, values, colors, range)?;
        crate::render::render_spacefill_to_array(py, &self.inner, colors, &options)
    }
}

impl Structure {
    fn check_atom(&self, atom: usize) -> PyResult<()> {
        let n = self.inner.atom_count();
        if atom >= n {
            return Err(PyIndexError::new_err(format!(
                "atom {atom} out of range for {n} atoms"
            )));
        }
        Ok(())
    }
}
