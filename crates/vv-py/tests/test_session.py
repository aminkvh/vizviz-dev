"""`vizviz.Session`: the command bus (load/select/undo/redo) from Python."""

import numpy as np
import pytest
import vizviz
from conftest import FIXTURES

HEMOGLOBIN = FIXTURES / "4HHB.cif"


@pytest.fixture
def session():
    return vizviz.Session()


def test_load_close_and_ids(session):
    assert session.structures == []
    sid = session.load(HEMOGLOBIN)
    assert isinstance(sid, int)
    assert session.structures == [sid]
    assert session.label(sid) == "4HHB"
    assert session.path(sid) == HEMOGLOBIN
    assert "structures=1" in repr(session)

    session.close(sid)
    assert session.structures == []
    with pytest.raises(KeyError):
        session.structure(sid)
    with pytest.raises(KeyError):
        session.close(sid)

    # Ids are never reused, so a reload gets a fresh one.
    sid2 = session.load(HEMOGLOBIN)
    assert sid2 != sid


def test_load_trajectory_produces_a_multi_frame_structure(session):
    sid = session.load_trajectory(FIXTURES / "1CRN.pdb", FIXTURES / "1CRN_traj.dcd")
    assert isinstance(sid, int)
    s = session.structure(sid)
    assert s.atom_count == 327
    assert s.frame_count == 3


def test_load_trajectory_rejects_an_atom_count_mismatch(session):
    with pytest.raises(ValueError):
        session.load_trajectory(HEMOGLOBIN, FIXTURES / "1CRN_traj.dcd")
    assert session.structures == []


def test_set_frame_is_undoable_and_validates(session):
    sid = session.load_trajectory(FIXTURES / "1CRN.pdb", FIXTURES / "1CRN_traj.dcd")
    assert session.frame(sid) == 0
    session.set_frame(sid, 2)
    assert session.frame(sid) == 2
    assert session.undo()
    assert session.frame(sid) == 0
    with pytest.raises(IndexError):
        session.set_frame(sid, 99)
    assert session.frame(sid) == 0
    # The `frame N` console verb dispatches the same command.
    out = session.exec("frame 1")[0]
    assert out.endswith("frame 1")
    assert session.frame(sid) == 1


def test_structure_shares_memory_with_the_session(session, hemoglobin):
    sid = session.load(HEMOGLOBIN)
    s = session.structure(sid)
    assert s.atom_count == hemoglobin.atom_count
    assert np.allclose(s.positions, hemoglobin.positions)
    # Two handles to the same loaded structure view the same Arc-owned memory.
    assert np.shares_memory(s.positions, session.structure(sid).positions)
    # ... and survive the session letting go of it.
    positions = s.positions
    session.close(sid)
    del s, session
    assert np.isfinite(positions).all()


def test_load_error_is_not_recorded(session, tmp_path):
    with pytest.raises(FileNotFoundError):
        session.load(tmp_path / "missing.cif")
    assert not session.can_undo
    assert session.structures == []


def test_representation_and_coloring(session):
    sid = session.load(HEMOGLOBIN)
    assert session.representation(sid) == "spacefill"
    assert session.coloring(sid) == "element"
    session.set_representation(sid, "ball_and_stick")
    session.set_coloring(sid, "chain")
    assert session.representation(sid) == "ball_and_stick"
    assert session.coloring(sid) == "chain"
    session.set_representation(sid, "tube")
    assert session.representation(sid) == "tube"
    assert session.exec("rep ballstick")[0].endswith("drawn as ballstick")
    # set_representation delegates to the same alias-tolerant parser as
    # exec()/the console, so it accepts an alias like "ribbon" too.
    session.set_representation(sid, "ribbon")
    assert session.representation(sid) == "cartoon"
    with pytest.raises(ValueError):
        session.set_representation(sid, "not-a-representation")
    with pytest.raises(ValueError):
        session.set_coloring(sid, "plaid")
    assert all(session.undo() for _ in range(5))
    assert session.representation(sid) == "spacefill"
    assert session.coloring(sid) == "element"


def test_set_coloring_accepts_the_consoles_space_syntax_too(session):
    # `color values NAME` (space-separated) works at the console; direct
    # Python calls used to require `"values:NAME"` only. Both now parse to
    # the same ColorScheme.
    sid = session.load(HEMOGLOBIN)
    atoms = session.structure(sid).atom_count
    session.set_values(sid, "sasa", np.zeros(atoms, dtype=np.float32))
    session.set_coloring(sid, "values sasa")
    assert session.coloring(sid) == "values:sasa"


def test_select_by_indices_and_by_mask(session):
    sid = session.load(HEMOGLOBIN)
    n = session.structure(sid).atom_count
    assert session.selection is None

    session.select(sid, [2, 0, 1, 1])
    got_id, indices = session.selection
    assert got_id == sid
    assert indices.dtype == np.uint32
    assert indices.tolist() == [0, 1, 2]

    session.select(sid, np.array([5, 7], dtype=np.int64))
    assert session.selection[1].tolist() == [5, 7]

    session.select(sid, range(10, 13))
    assert session.selection[1].tolist() == [10, 11, 12]

    mask = np.zeros(n, dtype=bool)
    mask[[3, 4]] = True
    session.select(sid, mask)
    assert session.selection[1].tolist() == [3, 4]

    # A strided view of a mask is fine too.
    wide = np.zeros(2 * n, dtype=bool)
    wide[2 * 9] = True
    session.select(sid, wide[::2])
    assert session.selection[1].tolist() == [9]

    session.select(sid, [True] + [False] * (n - 1))
    assert session.selection[1].tolist() == [0]

    session.clear_selection()
    assert session.selection is None


def test_select_errors(session):
    sid = session.load(HEMOGLOBIN)
    n = session.structure(sid).atom_count
    with pytest.raises(ValueError):
        session.select(sid, np.zeros(n - 1, dtype=bool))
    with pytest.raises(IndexError):
        session.select(sid, [n])
    with pytest.raises(IndexError):
        session.select(sid, [-1])
    with pytest.raises(TypeError):
        session.select(sid, 3.5)
    with pytest.raises(KeyError):
        session.select(sid + 100, [0])
    with pytest.raises(ValueError, match="characters"):
        session.select(sid, "chain A and")
    assert session.selection is None


def test_select_by_expression(session, hemoglobin):
    sid = session.load(HEMOGLOBIN)
    session.select(sid, "chain A and name CA")
    got_id, indices = session.selection
    assert got_id == sid and len(indices) == 141
    assert session.selection_expr == "chain A and name CA"
    assert np.array_equal(indices, hemoglobin.select("chain A and name CA"))

    session.save_selection_set("ca")
    assert session.selection_set_expr("ca") == "chain A and name CA"
    session.select(sid, [0])
    assert session.selection_expr is None
    session.undo()
    assert session.selection_expr == "chain A and name CA"
    with pytest.raises(KeyError):
        session.selection_set_expr("nope")


def test_selection_sets(session):
    sid = session.load(HEMOGLOBIN)
    with pytest.raises(ValueError):
        session.save_selection_set("pocket")  # nothing selected
    session.select(sid, [1, 2, 3])
    session.save_selection_set("pocket")
    assert session.selection_sets == ["pocket"]
    got_id, indices = session.selection_set("pocket")
    assert got_id == sid and indices.tolist() == [1, 2, 3]

    session.select(sid, [8])
    session.save_selection_set("pocket")  # overwrite
    assert session.selection_set("pocket")[1].tolist() == [8]
    session.undo()
    assert session.selection_set("pocket")[1].tolist() == [1, 2, 3]

    session.delete_selection_set("pocket")
    assert session.selection_sets == []
    with pytest.raises(KeyError):
        session.selection_set("pocket")
    with pytest.raises(KeyError):
        session.delete_selection_set("pocket")
    session.undo()
    assert session.selection_sets == ["pocket"]


def test_undo_redo_and_version(session):
    v0 = session.version
    assert not session.can_undo and not session.can_redo
    assert session.undo() is False and session.redo() is False

    sid = session.load(HEMOGLOBIN)
    v1 = session.version
    assert v1 > v0
    session.select(sid, [0])
    session.clear_selection()
    assert session.selection is None

    assert session.undo() is True
    assert session.selection[1].tolist() == [0]
    assert session.can_redo
    assert session.redo() is True
    assert session.selection is None

    assert session.undo() and session.undo() and session.undo()
    assert session.structures == []
    assert not session.can_undo
    assert session.redo()
    assert session.structures == [sid]
    assert session.version > v1


def test_exec_runs_text_commands(session):
    out = session.exec(f"load {HEMOGLOBIN}\nselect chain A and name CA; saveset ca\ncolor chain")
    assert len(out) == 4
    assert "4HHB" in out[0] and "4779 atoms" in out[0]
    assert out[1] == "selected 141 atom(s)"
    sid = session.structures[0]
    assert session.coloring(sid) == "chain"
    assert session.selection_sets == ["ca"]
    assert session.selection_set_expr("ca") == "chain A and name CA"
    assert session.exec("undo")[0] == "undone"
    assert session.coloring(sid) == "element"
    assert "#0 4HHB" in session.exec("structures")[0]

    with pytest.raises(ValueError, match="`bogus`"):
        session.exec("clear; bogus; clear")
    with pytest.raises(ValueError, match="no structure with id 9"):
        session.exec("color chain 9")
    assert session.exec("") == []
    assert session.exec("# just a comment") == []


def test_session_file_round_trip(session, tmp_path):
    sid = session.load(HEMOGLOBIN)
    session.select(sid, "chain A and name CA")
    session.save_selection_set("ca")
    session.set_coloring(sid, "chain")
    session.select(sid, [3, 5])  # a picked selection: stored as indices
    path = tmp_path / "work.vizviz.json"
    session.save_session(path)
    text = path.read_text()
    assert '"schema_version": 4' in text and "chain A and name CA" in text

    other = vizviz.Session()
    other.load(FIXTURES / "1CRN.cif")  # gets replaced
    assert other.load_session(path) == []
    (new_id,) = other.structures
    assert other.label(new_id) == "4HHB"
    assert other.coloring(new_id) == "chain"
    assert other.selection_set_expr("ca") == "chain A and name CA"
    assert other.selection_set("ca")[1].tolist() == session.selection_set("ca")[1].tolist()
    assert other.selection[1].tolist() == [3, 5]

    # Through exec too, and a moved file is a warning, not a crash.
    assert other.exec(f"savesession {tmp_path / 'again.json'}")[0].startswith("saved session")
    broken = tmp_path / "broken.json"
    broken.write_text(text.replace("4HHB.cif", "moved.cif"))
    warnings = other.load_session(broken)
    assert warnings and "structure #0" in warnings[0]
    assert other.structures == []
    with pytest.raises(ValueError):
        other.load_session(tmp_path / "missing.json")


def test_session_file_persists_value_channels(session, tmp_path):
    sid = session.load(HEMOGLOBIN)
    atoms = session.structure(sid).atom_count
    sasa = np.linspace(0.0, 1.0, atoms, dtype=np.float32)
    session.set_values(sid, "sasa", sasa)
    session.set_coloring(sid, "values:sasa")
    path = tmp_path / "with_values.vizviz.json"
    session.save_session(path)
    values_dir = tmp_path / "with_values.vizviz.values"
    assert values_dir.is_dir()
    assert list(values_dir.glob("*.npy"))

    other = vizviz.Session()
    assert other.load_session(path) == []
    (new_id,) = other.structures
    assert other.coloring(new_id) == "values:sasa"
    assert other.values(new_id) == ["sasa"]
    assert np.allclose(other.get_values(new_id, "sasa"), sasa)

    # A deleted sidecar is a warning, not a crash; the coloring stays put.
    import shutil

    shutil.rmtree(values_dir)
    other2 = vizviz.Session()
    warnings = other2.load_session(path)
    assert warnings and "sasa" in warnings[0]
    (new_id2,) = other2.structures
    assert other2.values(new_id2) == []
    assert other2.coloring(new_id2) == "values:sasa"


def test_commands_manifest():
    commands = vizviz.Session.commands()
    ids = [c[0] for c in commands]
    assert "load" in ids and "select" in ids and "undo" in ids
    assert "style" not in ids, "view verbs belong to the app, not a headless session"
    for cid, usage, help_text in commands:
        assert usage.startswith(cid) and help_text


def test_undo_limit_caps_history():
    session = vizviz.Session(undo_limit=1)
    sid = session.load(HEMOGLOBIN)
    session.select(sid, [0])
    assert session.undo() is True
    assert session.undo() is False, "the load fell off the capped history"
    assert session.structures == [sid]
