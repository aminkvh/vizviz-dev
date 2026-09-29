"""Piping through other packages: structures from arrays, per-atom value
channels back in, and the MDAnalysis bridge (skipped when MDAnalysis is
not installed)."""

import numpy as np
import pytest
import vizviz
from conftest import FIXTURES

HEMOGLOBIN = FIXTURES / "4HHB.cif"


def columns_of(s):
    """The per-atom columns another package would hand us."""
    return dict(
        names=s.name,
        elements=s.element,
        residue_names=np.asarray(s.residue_names)[s.residue_index],
        residue_ids=s.residue_auth_seq_id[s.residue_index],
        chain_ids=np.asarray(s.chain_names)[s.residue_chain[s.residue_index]],
        auth_chain_ids=np.asarray(s.chain_auth_names)[s.residue_chain[s.residue_index]],
        b_factors=s.b_factor,
        hetero=(s.flags & 1) != 0,
    )


def test_from_arrays_rebuilds_a_parsed_structure(hemoglobin):
    s = hemoglobin
    t = vizviz.from_arrays(s.positions, id="copy", **columns_of(s))
    assert t.id == "copy"
    assert (t.atom_count, t.residue_count, t.chain_count) == (
        s.atom_count,
        s.residue_count,
        s.chain_count,
    )
    assert np.array_equal(t.positions, s.positions)
    assert np.array_equal(t.element, s.element)
    assert np.array_equal(t.name, s.name)
    assert np.array_equal(t.residue_index, s.residue_index)
    assert np.array_equal(t.b_factor, s.b_factor)
    assert t.residue_names == s.residue_names
    assert t.chain_names == s.chain_names
    assert t.chain_auth_names == s.chain_auth_names
    # The selection language and the analysis kernels see the same atoms.
    expr = "chain A and resname HIS and name NE2"
    assert np.array_equal(t.select(expr), s.select(expr))
    assert t.contacts(t.select("chain A"), t.select("chain B"), 4.0)[0].shape == (
        s.contacts(s.select("chain A"), s.select("chain B"), 4.0)[0].shape
    )


def test_from_arrays_defaults_and_element_guessing():
    pos = np.array([[0.0, 0, 0], [1.5, 0, 0], [3.0, 0, 0]])
    m = vizviz.from_arrays(pos, names=["N", "CA", "FE"])
    assert m.element.tolist() == [7, 6, 26]
    assert m.residue_names == ["UNK"] and m.chain_names == ["A"]
    assert m.atom_name(1) == "CA"
    by_symbol = vizviz.from_arrays(pos, elements=["c", "Fe", "Zn"])
    assert by_symbol.element.tolist() == [6, 26, 30]
    bare = vizviz.from_arrays(np.zeros((4, 3)))
    assert bare.atom_count == 4 and bare.element.tolist() == [0, 0, 0, 0]
    frames = vizviz.from_arrays(np.random.default_rng(0).random((6, 5, 3)))
    assert frames.frame_count == 6 and frames.atom_count == 5


def test_from_arrays_groups_residues_and_chains_by_runs():
    pos = np.zeros((6, 3))
    m = vizviz.from_arrays(
        pos,
        names=["N", "CA", "N", "CA", "N", "CA"],
        residue_names=["GLY", "GLY", "ALA", "ALA", "GLY", "GLY"],
        residue_ids=[1, 1, 2, 2, 1, 1],
        chain_ids=["A", "A", "A", "A", "B", "B"],
    )
    assert m.residue_count == 3 and m.chain_count == 2
    assert m.residue_index.tolist() == [0, 0, 1, 1, 2, 2]
    assert m.residue_names == ["GLY", "ALA", "GLY"]
    assert m.chain_names == ["A", "B"]
    assert m.select("chain B").tolist() == [4, 5]


def test_from_arrays_rejects_mismatched_columns():
    with pytest.raises(ValueError, match="names"):
        vizviz.from_arrays(np.zeros((3, 3)), names=["N", "CA"])
    with pytest.raises(ValueError, match="positions"):
        vizviz.from_arrays(np.zeros((3, 2)))
    with pytest.raises(ValueError, match="b_factors"):
        vizviz.from_arrays(np.zeros((3, 3)), b_factors=[1.0])


@pytest.mark.skipif(not vizviz.HAS_RENDER, reason="built without the render feature")
def test_render_by_values_and_explicit_colors(hemoglobin):
    s = hemoglobin
    score = np.linspace(0.0, 1.0, s.atom_count)
    by_values = s.render(48, 32, values=score)
    by_range = s.render(48, 32, values=score, range=(-10.0, 10.0))
    assert by_values.shape == (32, 48, 4) and not np.array_equal(by_values, by_range)
    per_residue = np.linspace(0.0, 1.0, s.residue_count)
    assert np.array_equal(
        s.render(48, 32, values=per_residue), s.render(48, 32, values=per_residue[s.residue_index])
    )
    # A constant column and an explicit uniform color draw the same shape.
    flat = s.render(48, 32, values=np.zeros(s.atom_count))
    red = np.tile(np.array([[200, 30, 30]], dtype=np.uint8), (s.atom_count, 1))
    by_colors = s.render(48, 32, colors=red)
    rgba = np.concatenate([red, np.full((s.atom_count, 1), 255, np.uint8)], axis=1)
    assert np.array_equal(by_colors, s.render(48, 32, colors=rgba))
    assert (flat[..., 3] > 0).sum() == (by_colors[..., 3] > 0).sum()
    with pytest.raises(ValueError, match="entries"):
        s.render(16, 16, values=score[:-1])
    with pytest.raises(ValueError, match="colors"):
        s.render(16, 16, colors=red[:, :2])


def test_session_value_channels_round_trip():
    session = vizviz.Session()
    sid = session.load(HEMOGLOBIN)
    n = session.structure(sid).atom_count
    score = np.arange(n, dtype=np.float64)
    session.set_values(sid, "score", score)
    assert session.values(sid) == ["score"]
    assert np.array_equal(session.get_values(sid, "score"), score.astype(np.float32))

    session.set_coloring(sid, "values:score")
    assert session.coloring(sid) == "values:score"
    assert session.exec("values")[0].startswith("score: 0.000 ..")
    if vizviz.HAS_RENDER:
        img = session.render(48, 32)
        assert img.shape == (32, 48, 4)

    # Per-frame channels take (frames, atoms); 4HHB has one frame.
    session.set_values(sid, "score", score[None, :] * 2)
    assert session.get_values(sid, "score").shape == (n,)
    # One value per residue is expanded to that residue's atoms.
    s = session.structure(sid)
    session.set_values(sid, "per_res", np.arange(s.residue_count, dtype=np.float32))
    expanded = session.get_values(sid, "per_res")
    assert expanded.shape == (n,) and expanded[-1] == s.residue_count - 1
    assert np.array_equal(expanded, s.residue_index.astype(np.float32))
    with pytest.raises(ValueError, match="per atom"):
        session.set_values(sid, "bad", score[:10])
    with pytest.raises(KeyError):
        session.get_values(sid, "nope")
    with pytest.raises(KeyError):
        session.remove_values(sid, "nope")

    # Undo walks back through set_values / set_coloring like anything else.
    session.remove_values(sid, "score")
    assert session.values(sid) == ["per_res"]
    assert session.undo() and session.values(sid) == ["per_res", "score"]
    assert session.exec("color values other")[0].endswith("attaches it)")
    assert session.coloring(sid) == "values:other"


def test_value_coloring_survives_a_session_file(tmp_path):
    """The FastSASA-style pipe (see examples/fastsasa_pipe.py), through a
    saved and reloaded session: the channel is written as a `.npy`
    sidecar and comes back attached, not just the coloring."""
    session = vizviz.Session()
    sid = session.load(HEMOGLOBIN)
    sasa = np.ones(session.structure(sid).atom_count)
    session.set_values(sid, "sasa", sasa)
    session.set_coloring(sid, "values:sasa")
    path = tmp_path / "s.vizviz.json"
    session.save_session(path)
    assert '"coloring": "values:sasa"' in path.read_text()
    other = vizviz.Session()
    assert other.load_session(path) == []
    (new_id,) = other.structures
    assert other.coloring(new_id) == "values:sasa"
    assert other.values(new_id) == ["sasa"]
    assert np.allclose(other.get_values(new_id, "sasa"), sasa)


def test_from_arrays_is_fast_enough_for_big_inputs():
    """2M atoms through from_arrays should be a second-ish, not a minute."""
    import time

    n = 2_000_000
    rng = np.random.default_rng(1)
    pos = rng.random((n, 3), dtype=np.float32) * 300
    names = np.tile(np.array(["N", "CA", "C", "O"]), n // 4)
    resid = np.repeat(np.arange(n // 4), 4)
    t0 = time.perf_counter()
    m = vizviz.from_arrays(pos, names=names, residue_ids=resid)
    dt = time.perf_counter() - t0
    assert m.atom_count == n and m.residue_count == n // 4
    assert dt < 20, f"from_arrays took {dt:.1f}s for 2M atoms"


def test_mdanalysis_round_trip(hemoglobin):
    pytest.importorskip("MDAnalysis", reason="MDAnalysis not installed")
    s = hemoglobin
    u = vizviz.to_mdanalysis(s)
    assert len(u.atoms) == s.atom_count
    assert len(u.residues) == s.residue_count
    assert len(u.segments) == s.chain_count
    assert u.atoms.names[1] == "CA"
    assert u.atoms.elements[0] == "N"
    assert np.allclose(u.atoms.positions, s.positions)
    ca = u.select_atoms("protein and name CA")
    assert len(ca) == len(s.select("protein and name CA"))

    back = vizviz.from_mdanalysis(u)
    assert back.atom_count == s.atom_count and back.residue_count == s.residue_count
    assert np.array_equal(back.element, s.element)
    assert back.residue_names == s.residue_names
    assert np.allclose(back.positions, s.positions)

    # An AtomGroup gives a sub-structure; a slice of frames a trajectory.
    chain_a = vizviz.from_mdanalysis(u.select_atoms("segid A"), id="A")
    assert chain_a.atom_count == len(s.select("chain A")) and chain_a.id == "A"
    multi = vizviz.from_arrays(np.random.default_rng(2).random((5, 20, 3)))
    u2 = vizviz.to_mdanalysis(multi)
    assert len(u2.trajectory) == 5
    assert vizviz.from_mdanalysis(u2, frames="all").frame_count == 5
    assert vizviz.from_mdanalysis(u2, frames=slice(0, 5, 2)).frame_count == 3
    assert vizviz.from_mdanalysis(u2, frames=[4]).frame_count == 1


def test_mdanalysis_reads_a_trajectory_as_vizviz_does():
    mda = pytest.importorskip("MDAnalysis", reason="MDAnalysis not installed")
    topology, trajectory = FIXTURES / "1CRN.pdb", FIXTURES / "1CRN_traj.dcd"
    session = vizviz.Session()
    ours = session.structure(session.load_trajectory(topology, trajectory))
    theirs = vizviz.from_mdanalysis(mda.Universe(topology, trajectory), frames="all")
    assert theirs.frame_count == ours.frame_count > 1
    for f in (0, ours.frame_count - 1):
        assert np.allclose(theirs.frame(f), ours.frame(f), atol=1e-3)


@pytest.mark.parametrize(
    "expr",
    ["protein", "backbone", "name CA", "resname CYS", "resid 1-10", "chain A and not backbone"],
)
def test_selections_match_mdanalysis(hemoglobin, expr):
    pytest.importorskip("MDAnalysis", reason="MDAnalysis not installed")
    u = vizviz.to_mdanalysis(hemoglobin)
    # vizviz's backbone includes the terminal OXT; MDAnalysis's does not.
    theirs = expr.replace("chain A", "segid A").replace(
        "backbone", "(backbone or (protein and name OXT))"
    )
    assert list(hemoglobin.select(expr)) == list(u.select_atoms(theirs).indices)
