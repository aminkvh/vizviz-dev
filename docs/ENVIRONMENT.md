# Environment & dependency strategy

Goal: clone the repo on Windows, Linux, or macOS and get a working,
GPU-capable build without hand-debugging toolchain mismatches — and the
same lockfiles work on an HPC login node.

## Rust

- `rust-toolchain.toml` pins an exact compiler; `rustup` installs it on
  first `cargo` call. Windows also needs the VS Build Tools C++ workload.
- `cargo deny check` enforces the dependency license rule from
  CONTRIBUTING.md (config in `deny.toml`).
- GPU/UI crates (`wgpu`, `winit`, `egui*`, `egui_dock`) are upgraded together,
  one PR per egui release — they only work in matching versions.

## Python: `uv`

`uv` manages the Python version (`.python-version`) and the environment;
the Python package lives in `crates/vv-py` and is built by `maturin`:

```bash
uv sync --project crates/vv-py     # builds the extension, installs dev deps
uv run --project crates/vv-py pytest crates/vv-py/tests
```

`uv sync` rebuilds the extension whenever a Rust source in the workspace
changes (`[tool.uv] cache-keys` in `crates/vv-py/pyproject.toml`). To run
`cargo check/clippy -p vv-py` directly, point pyo3 at a Python 3.10+:
`PYO3_PYTHON=crates/vv-py/.venv/Scripts/python.exe` (or `bin/python`).

`pixi` (conda-forge) is deferred until we need conda-only scientific
packages; v1's Python dependencies are all on PyPI.

## GPU / driver handling

- Baseline renderer targets wgpu's native backends (Vulkan on
  Windows/Linux, Metal on macOS, DX12 as Windows fallback) — must always
  work, no vendor lock-in.
- Optional acceleration (hardware RT, CUDA for future compute plugins) is
  probed at startup; the app reports what's active/unavailable, never a
  silent slowdown.
- Minimum driver versions documented here once the renderer's real
  requirements are known.

## Containers

- **Dev/CI:** Docker (not needed for v1).
- **HPC:** Apptainer/Singularity — most clusters don't grant Docker's
  required privileges. See [HPC.md](HPC.md) and [CONTAINERS.md](CONTAINERS.md).

## License boundary

Core is Apache-2.0, so GPL-licensed tools (MDAnalysis, GPL-2.0+; GNINA and
LAMMPS default builds) are never linked in — they run as separate-process
plugins over IPC or in containers, in their own environment. LGPL
libraries (GROMACS, OpenMM's GPU platforms, MDTraj) may only be
dynamically linked.
