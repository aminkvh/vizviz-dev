# Roadmap

Open work, roughly in priority order. Not a schedule.

## Rendering at scale

- Hierarchical level of detail: distant clusters of atoms drawn as one
  proxy, so frame cost stops scaling with atom count.
- Sparse brick volumes for the Gaussian surface (finer voxels up close,
  cheaper rebuilds during trajectory playback).
- GPU-generated cartoon with level of detail.
- Skin and solvent-excluded surfaces at multi-million-atom scale: GPU
  build and per-patch culling.
- Half-resolution ambient occlusion for integrated GPUs and 4K.
- Lit shading for the distant point level of detail.
- Surface colouring as a distance-weighted blend across overlapping atoms
  instead of nearest atom.
- Hardware ray tracing as an optional offline render backend.

## Analysis and interop

- Colormap and range controls per value channel in the UI.
- A `to_mdtraj` counterpart to `to_mdanalysis`.
- Writing edited annotations back to mmCIF/PDB.
- Secondary structure (native or wrapped) recomputed per frame as a live
  value channel.
- Cryo-EM density and volume rendering, with map fetch (EMDB, 2Fo-Fc)
  next to a fetched entry.

## Extensibility

- Out-of-process plugin manifest for tools without a Python API
  ([Containers](docs/CONTAINERS.md)), starting with permissively licensed
  structure-prediction and docking tools.
- WASM plugin host for sandboxed, native-speed third-party plugins.
- In-app Python console mirroring the command bus.
- Additional HPC schedulers (PBS/Torque, LSF) beyond Slurm.

## Ideas, unscoped

- Fetch by ID into an editable, styleable scene (PubChem, RCSB).
- A notebook-style panel for reproducible figure and analysis cells.
- Several AI agents driving the app concurrently.
- Plugin registry and installer on top of the manifest format.
