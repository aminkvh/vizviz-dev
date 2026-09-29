# Analysis

Native Rust kernels in `vv_core::analysis`, reached from Python
(`Structure` methods), the app (inspector, `measure`, `contacts`), and any
future trajectory reader. Rule from docs/VALIDATION.md: nothing ships
without a check against an independent reference.

## Design rule: assume 40,000 frames or 2,000,000 structures

Every kernel is written for batch use first:

- Inputs are a coordinate slice plus index arrays, never a `Structure`.
  The same call serves one crystal structure, one NMR model, or a frame
  streamed from a trajectory.
- The `*_into` form writes into a caller-provided buffer and allocates
  nothing per call; the scalar form is a convenience over it.
- The `*_frames` form runs frames in parallel with rayon.
- Contacts and neighbor lists use a uniform grid (the same one behind
  `within` in selections), built over the smaller group, walked in
  parallel over the larger.

Measured in release mode on a 16-thread laptop (`cargo test -p vv-io
--test analysis --release -- --ignored --nocapture`):

| Workload | Time |
|---|---|
| 40,000 frames x 573 phi dihedrals (22.9 M values) | 78 ms |
| 40,000 frames, contact count chain A vs B of 4HHB at 4 Å | 276 ms (7 µs/frame) |
| Neighbor list, 4HHB (4,779 atoms) at 4 Å, 29 K pairs | 4 ms |
| Neighbor list, 3J3Q (2.44 M atoms) at 4 Å, 14.3 M pairs | 0.97 s |

For millions of separate files, `vizviz.load_many(paths)` parses in
parallel with the GIL released; do it in batches, since memory, not
time, is the limit.

## What exists

| Kernel | Python | Reference it was checked against |
|---|---|---|
| distance | `s.distance(a, b)`, `s.distances(pairs, frames=None)` | NumPy |
| angle (degrees) | `s.angle(a, b, c)`, `s.angles(triples, frames=None)` | NumPy; backbone N-CA-C is 100-122° on 1CRN |
| dihedral (degrees, IUPAC sign) | `s.dihedral(a, b, c, d)`, `s.dihedrals(quads, frames=None)` | NumPy; crambin's helix has phi -100..-40, psi -70..-10 |
| contacts between two groups | `s.contacts(a, b, cutoff, frame=0)` → pairs, distances | brute-force O(n·m) on 4HHB; Fe–His NE2 coordination = exactly 4 |
| contact counts per frame | `s.contact_counts(a, b, cutoff, frames=None)` | equals `len(contacts)` per frame |
| neighbor list within one group | `s.neighbors(cutoff, atoms=None, frame=0)` | pairs under 1.9 Å are the perceived bonds |
| residue pairs from atom pairs | `s.residue_pairs(pairs)` | brute force |

Groups (`a`, `b`, `atoms`) are selection expressions or index sequences.
Multi-frame results are `(frames, k)` float32, one row per frame.

In the app: Ctrl+click 2, 3 or 4 atoms and the inspector shows the
distance, angle or dihedral in click order; `measure A B [C [D]]` and
`contacts 4 chain A | chain B` do the same from the command line.

## Not yet, and how each will arrive

- **RMSD / RMSF / alignment**: needs a second frame or structure. Multi-
  model files load today, so RMSD across an ensemble is the next kernel;
  trajectories follow the trajectory reader.
- **Trajectory readers (XTC, DCD, TRR)**: out-of-core streaming, frame
  windows, never "load the whole file". The kernels above already take
  frames as slices for this reason.
- **Secondary structure, SASA, hydrogen bonds, salt bridges**: each ships
  only with a validation against a reference tool.
- **Heavier analyses** (docking scores, pocket detection, MD setup) come
  as external packages behind the plugin/container layer, not as native
  kernels: that is the modularity rule in README.md. Lahuta is an
  inspiration for the contact/interaction vocabulary, not a dependency.

## Piping through other packages (what exists today)

The in-process form of that rule works now and needs no plugin layer:
`Structure` columns out as NumPy arrays, a per-atom result back in as a
value channel (`Session.set_values`, the `values` command,
`render(values=...)`). FastSASA runs on a vizviz structure unchanged
(`examples/fastsasa_pipe.py`, docs/PYTHON.md); MDAnalysis objects go in
and out through `from_mdanalysis` / `to_mdanalysis`. What this does not
cover is the out-of-process case (a container, a Slurm job, a tool with
no Python API): that is the manifest layer in docs/CONTAINERS.md, still
unbuilt.
