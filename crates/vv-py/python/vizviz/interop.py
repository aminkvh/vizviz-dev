"""Bring your own arrays: build a `Structure` from NumPy columns, and move
structures to and from MDAnalysis without touching a file.

vizviz never imports MDAnalysis itself; `from_mdanalysis` duck-types the
AtomGroup attributes it needs and `to_mdanalysis` imports MDAnalysis only
when called. Missing optional columns (elements, B-factors, chain ids)
get sensible defaults instead of errors.
"""

from __future__ import annotations

import numpy as np

from . import _core

__all__ = ["from_arrays", "from_mdanalysis", "to_mdanalysis"]

_SYMBOLS = {_core.element_symbol(z).upper(): z for z in range(1, 119)}


def _atomic_numbers(elements, names):
    """`(atoms,)` uint8 atomic numbers from symbols, numbers, or (when
    `elements` is None) a guess from atom names, the way PDB files
    without an element column are read."""
    if elements is None:
        return _core._guess_elements(np.asarray(names).astype("S4"))
    arr = np.asarray(elements)
    if arr.dtype.kind in "iu":
        return arr.astype(np.uint8)
    if arr.dtype.kind == "O":  # MDAnalysis hands out object arrays of str
        arr = arr.astype("U")
    if arr.dtype.kind in "SU":
        symbols = np.char.upper(np.char.strip(arr.astype("U2")))
        return np.fromiter(
            (_SYMBOLS.get(s, 0) for s in symbols), dtype=np.uint8, count=len(symbols)
        )
    raise TypeError("elements must be symbols ('C', 'Fe'), atomic numbers, or None")


def _column(value, n, default, dtype, what):
    if value is None:
        return np.full(n, default, dtype=dtype)
    arr = np.asarray(value)
    if arr.shape != (n,):
        raise ValueError(f"{what} has shape {arr.shape}, expected ({n},)")
    return np.ascontiguousarray(arr, dtype=dtype)


def _strings(value, n, default, what):
    if value is None:
        return np.full(n, default, dtype=f"U{max(1, len(default))}")
    arr = np.asarray(value)
    if arr.shape != (n,):
        raise ValueError(f"{what} has shape {arr.shape}, expected ({n},)")
    return np.char.strip(arr.astype("U"))


def from_arrays(
    positions,
    *,
    names=None,
    elements=None,
    residue_names=None,
    residue_ids=None,
    chain_ids=None,
    auth_chain_ids=None,
    b_factors=None,
    occupancies=None,
    hetero=None,
    id="",
    title="",
):
    """A `Structure` from per-atom columns.

    `positions` is `(atoms, 3)` or `(frames, atoms, 3)` in Angstrom. Every
    other argument is optional and per atom: `names` (`"CA"`), `elements`
    (symbols or atomic numbers; guessed from names when omitted),
    `residue_names` (`"ALA"`), `residue_ids` (ints), `chain_ids` (`"A"`;
    the mmCIF `label_asym_id`), `auth_chain_ids` (the author / PDB chain
    letter the selection language's `chain` matches; defaults to
    `chain_ids`), `b_factors`, `occupancies`, and `hetero` (bools, True
    for ligands, ions and water). Residues are consecutive runs of the
    same chain, residue id and residue name, so atoms should be grouped
    the way a PDB file lists them.
    """
    pos = np.asarray(positions, dtype=np.float32)
    if pos.ndim == 2:
        pos = pos[None]
    if pos.ndim != 3 or pos.shape[-1] != 3:
        raise ValueError(f"positions must be (atoms, 3) or (frames, atoms, 3), got {pos.shape}")
    pos = np.ascontiguousarray(pos)
    n = pos.shape[1]

    names_u = _strings(names, n, "X", "names")
    res_names = _strings(residue_names, n, "UNK", "residue_names")
    chains = _strings(chain_ids, n, "A", "chain_ids")
    auth = chains if auth_chain_ids is None else _strings(auth_chain_ids, n, "A", "auth_chain_ids")
    res_uniq, res_inv = np.unique(res_names, return_inverse=True)
    chain_uniq, chain_inv = np.unique(chains, return_inverse=True)
    auth_uniq, auth_inv = np.unique(auth, return_inverse=True)

    return _core._from_columns(
        pos,
        _atomic_numbers(elements, names_u),
        np.ascontiguousarray(names_u.astype("S4")),
        [str(s) for s in res_uniq],
        np.ascontiguousarray(res_inv.reshape(-1), dtype=np.uint32),
        _column(residue_ids, n, 1, np.int32, "residue_ids")
        if residue_ids is not None
        else np.ones(n, dtype=np.int32),
        [str(s) for s in chain_uniq],
        np.ascontiguousarray(chain_inv.reshape(-1), dtype=np.uint32),
        [str(s) for s in auth_uniq],
        np.ascontiguousarray(auth_inv.reshape(-1), dtype=np.uint32),
        _column(b_factors, n, 0.0, np.float32, "b_factors"),
        _column(occupancies, n, 1.0, np.float32, "occupancies"),
        _column(hetero, n, False, np.bool_, "hetero"),
        str(id),
        str(title),
    )


def _optional(group, attr):
    """A topology attribute of an MDAnalysis group, or None when the
    topology does not carry it (MDAnalysis raises NoDataError, an
    AttributeError subclass)."""
    try:
        return getattr(group, attr)
    except AttributeError:
        return None


def from_mdanalysis(atoms, *, frames=None, id=None):
    """A `Structure` from an MDAnalysis `Universe` or `AtomGroup`.

    Takes the group's current frame, or `frames` (a slice, a list of
    frame indices, or `"all"`) from its trajectory as the structure's
    frames. Names, residue names and ids come from the topology; elements
    from the `elements` attribute when present, otherwise guessed from
    names; chains from `chainIDs` when present, otherwise `segids`;
    B-factors from `tempfactors` when present. MDAnalysis is not imported
    by vizviz, only duck-typed.
    """
    group = atoms.atoms  # a Universe and an AtomGroup both have .atoms
    if frames is None:
        positions = np.asarray(group.positions, dtype=np.float32)[None]
    else:
        trajectory = group.universe.trajectory
        if isinstance(frames, str) and frames == "all":
            frames = slice(None)
        if isinstance(frames, slice):
            indices = range(*frames.indices(len(trajectory)))
        else:
            indices = [int(i) for i in frames]
        positions = np.empty((len(indices), len(group), 3), dtype=np.float32)
        for k, i in enumerate(indices):
            trajectory[i]
            positions[k] = group.positions
    segids = _optional(group, "segids")
    chain_ids = _optional(group, "chainIDs")
    chains = segids if segids is not None else chain_ids
    record_types = _optional(group, "record_types")
    hetero = None if record_types is None else np.asarray(record_types) == "HETATM"
    label = id if id is not None else ""
    return from_arrays(
        positions,
        names=_optional(group, "names"),
        elements=_optional(group, "elements"),
        residue_names=_optional(group, "resnames"),
        residue_ids=_optional(group, "resids"),
        chain_ids=chains,
        auth_chain_ids=chain_ids,
        b_factors=_optional(group, "tempfactors"),
        occupancies=_optional(group, "occupancies"),
        hetero=hetero,
        id=label,
    )


def to_mdanalysis(structure, *, bonds=False):
    """An MDAnalysis `Universe` with this structure's topology and all
    its frames (an in-memory trajectory). `bonds=True` also runs vizviz's
    bond perception and attaches the pairs. Requires MDAnalysis."""
    import MDAnalysis as mda
    from MDAnalysis.coordinates.memory import MemoryReader

    s = structure
    n = s.atom_count
    u = mda.Universe.empty(
        n,
        n_residues=s.residue_count,
        n_segments=s.chain_count,
        atom_resindex=np.asarray(s.residue_index, dtype=np.int64),
        residue_segindex=np.asarray(s.residue_chain, dtype=np.int64),
        trajectory=False,
    )
    names = np.char.strip(np.asarray(s.name).astype("U4"))
    symbols = np.array([_core.element_symbol(z) or "" for z in range(256)], dtype="U2")
    u.add_TopologyAttr("names", names)
    u.add_TopologyAttr("types", names)
    u.add_TopologyAttr("elements", symbols[np.asarray(s.element)])
    u.add_TopologyAttr("tempfactors", np.asarray(s.b_factor, dtype=np.float64))
    u.add_TopologyAttr("occupancies", np.asarray(s.occupancy, dtype=np.float64))
    hetero = (np.asarray(s.flags) & 1) != 0
    u.add_TopologyAttr("record_types", np.where(hetero, "HETATM", "ATOM"))
    u.add_TopologyAttr("resnames", list(s.residue_names))
    u.add_TopologyAttr("resids", np.asarray(s.residue_auth_seq_id, dtype=np.int64))
    u.add_TopologyAttr("segids", list(s.chain_auth_names))
    chain_of_atom = np.asarray(s.residue_chain)[np.asarray(s.residue_index)]
    u.add_TopologyAttr("chainIDs", np.asarray(s.chain_auth_names, dtype="U")[chain_of_atom])
    if bonds:
        u.add_TopologyAttr("bonds", [tuple(map(int, p)) for p in s.bonds()])
    coords = np.stack([np.asarray(s.frame(f), dtype=np.float32) for f in range(s.frame_count)])
    u.load_new(coords, format=MemoryReader)
    return u
