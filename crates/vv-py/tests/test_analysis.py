"""Analysis kernels: checked against plain NumPy, not against themselves."""

import numpy as np
import pytest
import vizviz
from conftest import FIXTURES


def numpy_dihedral(p0, p1, p2, p3):
    """The standard formula (as used by MDAnalysis); IUPAC sign."""
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1 = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    return np.degrees(np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w)))


def backbone(s, residue, name):
    start, end = s.residue_atoms[residue]
    names = np.char.strip(s.name[start:end])
    return int(start + np.flatnonzero(names == name.encode())[0])


def test_angle_and_dihedral_match_numpy(hemoglobin):
    s = hemoglobin
    p = s.positions.astype(np.float64)
    n, ca, c = (backbone(s, 5, x) for x in ("N", "CA", "C"))
    u, v = p[n] - p[ca], p[c] - p[ca]
    expected = np.degrees(np.arccos(np.dot(u, v) / np.linalg.norm(u) / np.linalg.norm(v)))
    assert s.angle(n, ca, c) == pytest.approx(expected, abs=1e-3)
    assert 100 < s.angle(n, ca, c) < 122  # tau is chemistry

    c_prev = backbone(s, 4, "C")
    phi = s.dihedral(c_prev, n, ca, c)
    assert phi == pytest.approx(numpy_dihedral(p[c_prev], p[n], p[ca], p[c]), abs=1e-3)
    assert -180 <= phi <= 180
    with pytest.raises(IndexError):
        s.angle(0, 1, s.atom_count)


def test_batch_measurements_have_one_row_per_frame(hemoglobin):
    s = hemoglobin
    quads = np.array(
        [
            [backbone(s, r - 1, "C"), *(backbone(s, r, x) for x in ("N", "CA", "C"))]
            for r in range(1, 30)
        ]
    )
    d = s.dihedrals(quads)
    assert d.shape == (1, len(quads)) and d.dtype == np.float32
    p = s.positions.astype(np.float64)
    ref = [numpy_dihedral(*(p[i] for i in q)) for q in quads]
    assert np.allclose(d[0], ref, atol=1e-3)
    # Same numbers through the scalar entry point and through a list.
    assert s.dihedrals(quads.tolist())[0, 0] == pytest.approx(s.dihedral(*quads[0]), abs=1e-5)

    pairs = np.array([[0, 1], [0, 2], [10, 20]])
    dist = s.distances(pairs)
    assert dist.shape == (1, 3)
    assert np.allclose(dist[0], np.linalg.norm(p[pairs[:, 0]] - p[pairs[:, 1]], axis=1), atol=1e-4)
    assert s.distances(pairs, frames=[0]).shape == (1, 3)
    with pytest.raises(IndexError):
        s.distances(pairs, frames=[1])
    with pytest.raises(ValueError):
        s.distances([[0, 1, 2]])
    assert s.angles(np.zeros((0, 3), dtype=np.int64)).shape == (1, 0)


def test_contacts_and_neighbors_against_brute_force(hemoglobin):
    s = hemoglobin
    p = s.positions.astype(np.float64)
    a = s.select("chain A and protein")
    b = s.select("chain B and protein")
    pairs, dist = s.contacts(a, b, 4.0)
    d = np.linalg.norm(p[a][:, None, :] - p[b][None, :, :], axis=2)
    ia, ib = np.nonzero(d <= 4.0)
    expected = sorted(zip(a[ia].tolist(), b[ib].tolist(), strict=True))
    assert pairs.dtype == np.uint32 and dist.dtype == np.float32
    assert list(map(tuple, pairs.tolist())) == expected
    assert np.allclose(dist, d[ia, ib][np.lexsort((b[ib], a[ia]))], atol=1e-4)

    # Expressions work as groups, and the counts agree with the pairs.
    pairs2, _ = s.contacts("chain A and protein", "chain B and protein", 4.0)
    assert np.array_equal(pairs, pairs2)
    counts = s.contact_counts("chain A and protein", "chain B and protein", 4.0)
    assert counts.shape == (1,) and counts[0] == len(pairs)

    residues = s.residue_pairs(pairs)
    ri = s.residue_index
    expected_res = sorted(
        {(min(ri[i], ri[j]), max(ri[i], ri[j])) for i, j in expected if ri[i] != ri[j]}
    )
    assert list(map(tuple, residues.tolist())) == expected_res

    # Fe coordinates the proximal histidine at ~2.1 A: exactly four.
    fe_pairs, fe_dist = s.contacts("element FE", "resname HIS and name NE2", 2.6)
    assert len(fe_pairs) == 4 and fe_dist.max() < 2.4

    # Neighbor pairs within 1.9 A of heavy protein atoms are the bonds.
    heavy = s.select("protein and not hydrogen")
    nb, _ = s.neighbors(1.9, heavy)
    assert (nb[:, 0] < nb[:, 1]).all()
    bonds = {tuple(x) for x in s.bonds().tolist()}
    covered = sum(tuple(x) in bonds for x in nb.tolist())
    assert covered > 0.95 * len(nb)
    with pytest.raises(ValueError):
        s.neighbors(0.0)
    with pytest.raises(ValueError):
        s.contacts("chain A and", b, 4.0)


def test_load_many_is_parallel_and_ordered():
    paths = [FIXTURES / n for n in ("1CRN.cif", "4HHB.cif", "1UBQ.pdb", "1CRN.pdb")]
    loaded = vizviz.load_many(paths)
    assert [s.id for s in loaded] == ["1CRN", "4HHB", "1UBQ", "1CRN"]
    assert loaded[1].atom_count == 4779
    with pytest.raises(FileNotFoundError):
        vizviz.load_many([FIXTURES / "1CRN.cif", FIXTURES / "nope.cif"])
    assert vizviz.load_many([]) == []
