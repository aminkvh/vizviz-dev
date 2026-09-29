"""Writes the trajectory fixtures: crambin (1CRN) with 10 frames of a
seeded random walk, as XTC, TRR and AMBER NetCDF, and for each the
coordinates MDAnalysis reads back from that file (the reference the vv-io
readers are tested against), as raw little-endian f32 (frames x atoms x 3).

    py -3 fixtures/traj/make.py
"""

from pathlib import Path

import MDAnalysis as mda
import numpy as np

here = Path(__file__).parent
u = mda.Universe(str(here.parent / "small" / "1CRN.pdb"))
rng = np.random.default_rng(7)
start = u.atoms.positions.copy()
frames = [
    start + np.cumsum(rng.normal(0, 0.3, (k + 1, *start.shape)), axis=0)[-1] for k in range(10)
]

for ext in ["xtc", "trr", "nc"]:
    path = here / f"1crn.{ext}"
    with mda.Writer(str(path), n_atoms=u.atoms.n_atoms) as w:
        for k, xyz in enumerate(frames):
            u.atoms.positions = xyz
            u.trajectory.ts.time = 2.0 * k
            w.write(u.atoms)
    back = mda.Universe(str(here.parent / "small" / "1CRN.pdb"), str(path))
    ref = np.stack([ts.positions.copy() for ts in back.trajectory]).astype("<f4")
    ref.tofile(here / f"1crn.{ext}.f32")
    print(path.name, path.stat().st_size, "bytes,", ref.shape)
