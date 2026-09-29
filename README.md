# vizviz

A fast molecular viewer for proteins, nucleic acids, glycans and their
trajectories, for Windows, Linux and macOS.

![Hemoglobin, path-traced by vizviz](docs/gallery/v02_hemoglobin_traced.jpg)

> **A personal project. Experimental and under active development.**
> Things change without notice, features are unfinished, and rough edges
> are expected. Use it, break it, and tell me what you find, but do not
> depend on it yet.

## Why

I study biomolecules for a living and wanted one viewer that stays smooth
on multi-million-atom structures and long trajectories, produces
publication-quality figures without a second program, and can be driven
from a script or Python as easily as from a mouse. vizviz is my attempt at
that.

The name is the sound a fly makes in Farsi.

## What it does

- Loads PDB, mmCIF and common trajectory formats (DCD, XTC, TRR, NetCDF);
  trajectories stream from disk, so they can be larger than memory.
- Draws spheres, sticks, cartoons, glycan symbols and several molecular
  surfaces, in as many layers per structure as you like.
- Path-traced still images, GIF and frame-sequence movies from a
  multi-track timeline.
- Selections by shape (box, circle, lasso) at atom, residue or chain
  level, or by a small expression language.
- Measurements, contacts, secondary structure and surface area.
- One command language behind every button, so anything in the interface
  can be scripted, run with `--exec`, or driven from Python.

## Install

From a checkout (needs the [Rust toolchain](https://rustup.rs)):

```bash
cargo install --path crates/vv-app --locked
```

Then:

```bash
vizviz path/to/structure.cif
vizviz --exec "load fixtures/small/4HHB.cif; color chain; screenshot fig.png; quit"
vizviz --exec help
```

You can also drop files onto the window or use File ▸ Fetch from RCSB.

## Documentation

- [Commands](docs/COMMANDS.md)
- [Selections](docs/SELECTION.md)
- [Movies](docs/MOVIE.md)
- [Python](docs/PYTHON.md)
- [File formats](docs/FORMATS.md)
- [Driving a live window](docs/LIVE.md)
- [Design system](DESIGN.md), [Contributing](CONTRIBUTING.md)

## License

Apache-2.0. See [LICENSE](LICENSE) and [NOTICE](NOTICE).
