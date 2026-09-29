"""`vizviz.load` and `Structure`: zero-copy NumPy views over a parsed file."""

import gc

import numpy as np
import pytest
import vizviz
from conftest import FIXTURES

PER_ATOM_COLUMNS = {
    "element": np.uint8,
    "serial": np.uint32,
    "b_factor": np.float32,
    "occupancy": np.float32,
    "alt_loc": np.uint8,
    "charge": np.int8,
    "flags": np.uint8,
    "residue_index": np.uint32,
}


def test_counts_and_repr(hemoglobin):
    s = hemoglobin
    assert s.atom_count == 4779  # human deoxyhemoglobin, 4 chains + 4 hemes + waters
    assert len(s) == s.atom_count
    assert s.frame_count == 1
    assert s.chain_count >= 4
    assert s.residue_count > 500
    assert s.id == "4HHB"
    assert "4HHB" in repr(s)
    assert isinstance(s, vizviz.Structure)


def test_positions_is_a_read_only_zero_copy_view(hemoglobin):
    s = hemoglobin
    p1 = s.positions
    p2 = s.positions
    assert p1.shape == (s.atom_count, 3)
    assert p1.dtype == np.float32
    assert p1.flags.c_contiguous
    # Two accesses give two array objects over the *same* memory, whose
    # owner is the Structure itself (no intermediate copy).
    assert p1 is not p2
    assert np.shares_memory(p1, p2)
    assert p1.base is s
    assert not p1.flags.writeable
    with pytest.raises(ValueError):
        p1[0, 0] = 1.0
    # A real coordinate from the file (N of VAL A 1).
    assert np.allclose(p1[0], [19.323, 29.727, 42.781], atol=1e-3)


@pytest.mark.parametrize("column", sorted(PER_ATOM_COLUMNS))
def test_every_per_atom_column_is_a_zero_copy_view(hemoglobin, column):
    s = hemoglobin
    a = getattr(s, column)
    assert a.shape == (s.atom_count,)
    assert a.dtype == PER_ATOM_COLUMNS[column]
    assert a.base is s
    assert not a.flags.writeable
    assert np.shares_memory(a, getattr(s, column))


def test_view_keeps_the_structure_alive():
    positions = vizviz.load(FIXTURES / "1CRN.cif").positions
    gc.collect()
    # The Structure was dropped by Python, but the array's `base` still
    # holds it, so the memory is intact.
    assert positions.shape[0] > 0
    assert np.isfinite(positions).all()
    assert positions.base is not None


def test_frames(hemoglobin):
    s = hemoglobin
    assert np.shares_memory(s.frame(0), s.positions)
    with pytest.raises(IndexError):
        s.frame(1)


def test_names_and_elements(hemoglobin):
    s = hemoglobin
    assert s.name.dtype == np.dtype("S4")
    assert s.name.base is s
    assert np.char.strip(s.name)[0] == b"N"
    assert s.atom_name(0) == "N"
    assert s.atom_name(1) == "CA"
    assert s.element[0] == 7
    assert vizviz.element_symbol(int(s.element[1])) == "C"
    assert (s.element > 0).all()
    # Fe in the hemes.
    assert (s.element == 26).sum() == 4


def test_b_factors_and_flags_are_real_data(hemoglobin):
    s = hemoglobin
    assert s.b_factor.std() > 1.0
    assert np.isfinite(s.b_factor).all()
    hetero = (s.flags & 1) != 0
    assert hetero.any() and not hetero.all()
    assert (s.occupancy > 0).all()


def test_residue_and_chain_tables_are_consistent(hemoglobin):
    s = hemoglobin
    ra = s.residue_atoms
    assert ra.shape == (s.residue_count, 2) and ra.dtype == np.uint32
    assert ra[0, 0] == 0 and ra[-1, 1] == s.atom_count
    assert (ra[:, 1] > ra[:, 0]).all()
    assert (ra[1:, 0] == ra[:-1, 1]).all(), "residues tile the atom array"
    assert (s.residue_index[ra[:, 0]] == np.arange(s.residue_count)).all()

    cr = s.chain_residues
    assert cr.shape == (s.chain_count, 2)
    assert cr[0, 0] == 0 and cr[-1, 1] == s.residue_count
    assert (cr[1:, 0] == cr[:-1, 1]).all(), "chains tile the residue array"
    assert (s.residue_chain[cr[:, 0]] == np.arange(s.chain_count)).all()

    assert s.residue_names[0] == "VAL"
    assert s.residue_name(0) == "VAL"
    assert s.names[s.residue_comp[0]] == "VAL"
    assert s.residue_seq_id[0] == 1
    assert s.residue_auth_seq_id[0] == 1
    assert s.residue_ins_code.dtype == np.uint8
    assert set(np.unique(s.residue_ss)) <= {0, 1, 2, 3}
    assert s.chain_names[0] == "A"
    assert s.chain_name(0) == "A"
    assert len(s.chain_auth_names) == s.chain_count
    assert s.chain_entity.shape == (s.chain_count,)
    with pytest.raises(IndexError):
        s.residue_name(s.residue_count)
    with pytest.raises(IndexError):
        s.chain_name(s.chain_count)


def test_annotations_and_info(hemoglobin, hemoglobin_pdb):
    info = hemoglobin.info
    assert "haemoglobin" in info["title"].lower()
    assert info["method"] == "X-RAY DIFFRACTION"
    assert info["resolution"] == pytest.approx(1.74, abs=0.01)
    assert info["deposition_date"]
    assert "P69905" in info["uniprot"] and "P68871" in info["uniprot"]
    assert len(info["entities"]) >= 2
    ann = hemoglobin.annotations
    assert ann["struct"][0]["title"] == info["title"]
    assert ann["exptl"][0]["method"] == "X-RAY DIFFRACTION"
    assert isinstance(ann["entity"], list) and "pdbx_description" in ann["entity"][0]
    # The PDB header lands in the same vocabulary.
    pdb = hemoglobin_pdb.info
    assert pdb["method"] == info["method"]
    assert pdb["resolution"] == pytest.approx(info["resolution"], abs=0.01)
    assert hemoglobin_pdb.annotations["entry"][0]["id"] == "4HHB"


def test_bonds(hemoglobin):
    s = hemoglobin
    bonds = s.bonds()
    assert bonds.ndim == 2 and bonds.shape[1] == 2 and bonds.dtype == np.uint32
    assert len(bonds) > s.atom_count * 0.8
    assert (bonds[:, 0] < bonds[:, 1]).all()
    assert bonds.max() < s.atom_count
    # N-CA of the first residue is a bond, ~1.46 A.
    assert (bonds == [0, 1]).all(axis=1).any()
    p = s.positions
    lengths = np.linalg.norm(p[bonds[:, 0]] - p[bonds[:, 1]], axis=1)
    assert lengths.max() < 2.5
    assert s.explicit_bonds.shape[1] == 2
    assert s.explicit_bond_kinds.shape == (len(s.explicit_bonds),)


def test_distance(hemoglobin):
    s = hemoglobin
    p = s.positions
    expected = float(np.linalg.norm(p[0] - p[1]))
    assert s.distance(0, 1) == pytest.approx(expected, rel=1e-5)
    assert s.distance(1, 0) == s.distance(0, 1)
    assert s.distance(0, 0) == 0.0
    with pytest.raises(IndexError):
        s.distance(0, s.atom_count)
    with pytest.raises(IndexError):
        s.distance(0, 1, frame=1)


def test_select_expressions(hemoglobin):
    s = hemoglobin
    ca = s.select("chain A and name CA")
    assert ca.dtype == np.uint32 and len(ca) == 141
    assert (np.char.strip(s.name)[ca] == b"CA").all()
    assert len(s.select("resname HEM")) == 172
    fe = s.select("element FE")
    assert (s.element[fe] == 26).all() and len(fe) == 4
    near = s.select("within 3 of element FE")
    assert set(fe.tolist()) <= set(near.tolist()) and len(near) > len(fe)
    assert len(s.select("all")) == s.atom_count
    assert len(s.select("none")) == 0
    assert len(s.select("protein")) + len(s.select("not protein")) == s.atom_count
    with pytest.raises(ValueError, match="characters"):
        s.select("chain A and")
    with pytest.raises(ValueError):
        s.select("bogus")


def test_pdb_and_mmcif_agree(hemoglobin, hemoglobin_pdb):
    cif, pdb = hemoglobin, hemoglobin_pdb
    assert cif.atom_count == pdb.atom_count
    assert np.allclose(cif.positions, pdb.positions, atol=1e-3)
    assert (cif.element == pdb.element).all()
    assert np.allclose(cif.b_factor, pdb.b_factor, atol=1e-2)


def test_save_round_trips_pdb_and_mmcif(hemoglobin, tmp_path):
    s = hemoglobin
    pdb_path = tmp_path / "4hhb.pdb"
    warnings = s.save(pdb_path)
    assert isinstance(warnings, list)
    reloaded = vizviz.load(pdb_path)
    assert reloaded.atom_count == s.atom_count
    assert np.allclose(reloaded.positions, s.positions, atol=1e-3)

    cif_path = tmp_path / "4hhb_out.cif"
    s.save(cif_path)
    reloaded_cif = vizviz.load(cif_path)
    assert reloaded_cif.atom_count == s.atom_count


def test_save_selection_writes_only_those_atoms(hemoglobin, tmp_path):
    s = hemoglobin
    out = tmp_path / "ca_only.pdb"
    s.save(out, selection="chain A and name CA")
    reloaded = vizviz.load(out)
    assert reloaded.atom_count == 141
    assert (np.char.strip(reloaded.name) == b"CA").all()

    out_idx = tmp_path / "by_index.pdb"
    s.save(out_idx, selection=[0, 1, 2])
    assert vizviz.load(out_idx).atom_count == 3


def test_save_frames_selects_coordinate_sets(hemoglobin, tmp_path):
    s = hemoglobin
    out = tmp_path / "frame0.pdb"
    s.save(out, frames=[0])
    assert vizviz.load(out).frame_count == 1


def test_save_errors(hemoglobin, tmp_path):
    s = hemoglobin
    with pytest.raises(ValueError, match="no atoms"):
        s.save(tmp_path / "empty.pdb", selection="resname XYZ")
    with pytest.raises(ValueError):
        s.save(tmp_path / "out.unknownext")
    with pytest.raises(IndexError):
        s.save(tmp_path / "out.pdb", frames=[999])


def test_load_errors(tmp_path):
    with pytest.raises(FileNotFoundError):
        vizviz.load(FIXTURES / "does-not-exist.cif")
    junk = tmp_path / "notes.txt"
    junk.write_text("hello")
    with pytest.raises(ValueError):
        vizviz.load(junk)
    empty = tmp_path / "empty.pdb"
    empty.write_text("HEADER    EMPTY\nEND\n")
    with pytest.raises(ValueError):
        vizviz.load(empty)
