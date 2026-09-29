# Contributing

Locked in so far:

- **License:** Apache-2.0, plain.
- **Core:** Rust, `wgpu` + `egui`, PyO3 Python bindings as first-class.

## Provenance: DCO, not a CLA

No CLA is needed. Contributions need a **DCO sign-off**
(`git commit -s`): a statement that you wrote the contribution or have the
right to submit it.

## Dependency license rule

Core dependencies must be permissive; `cargo deny check` enforces it
([deny.toml](deny.toml)).

- **Permissive** (MIT/BSD/Apache-2.0, or dual-licensed with a permissive
  option): fine to link into the core.
- **LGPL** (e.g. GROMACS, OpenMM's GPU platforms): dynamic linking only,
  never static or vendored. A container still counts as distributing the
  binary, so ship the license text and source offer.
- **GPL or non-commercial** (e.g. NAMD, AMBER's `pmemd`, GNINA in its
  default build): never linked, in any form. Integrate as a
  separate-process plugin over IPC or container, or leave as
  bring-your-own for tools whose license bars redistribution (see
  [docs/CONTAINERS.md](docs/CONTAINERS.md)).
- **Check what actually gets installed, not just the tool's own LICENSE.**
  A clean top-level license can still pull in a restrictively licensed
  dependency.

## Code quality

- **Rust:** `rustfmt` and `clippy` with warnings denied; `cargo deny check`
  for licenses and advisories.
- **Python:** `ruff` for lint and format.
- Pre-commit hooks run both locally before CI does
  (`uv tool install pre-commit && pre-commit install`).
- Comments state only what the code cannot: a non-obvious reason, a
  constraint, a units contract, a citation. Keep functions small and
  single-purpose.

## Testing

Three tiers:

1. **Unit tests**: per-module correctness.
2. **Integration tests**: real file formats and small structures, no
   mocked internals.
3. **Scientific regression tests**: see [docs/VALIDATION.md](docs/VALIDATION.md).
   A passing unit test proves the code runs, not that the number is right.

CI runs the full suite on Windows, Linux and macOS.

## Read first

- [docs/VERSIONING.md](docs/VERSIONING.md)
- [docs/ENVIRONMENT.md](docs/ENVIRONMENT.md)
- [docs/HPC.md](docs/HPC.md)
- [docs/VALIDATION.md](docs/VALIDATION.md)
- [docs/UI_DESIGN.md](docs/UI_DESIGN.md)
