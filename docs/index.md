# vizviz

Fast, open molecular visualization and analysis for Windows, Linux and
macOS. Millions of atoms at interactive frame rates, on a laptop's
integrated GPU too; publication figures path-traced in seconds; every
action a command you can script, from Python or an AI agent.

![Hemoglobin, path-traced](gallery/v02_hemoglobin_traced.jpg)

## Install

Download the installer for your system from the
[releases page](https://github.com/aminkvh/vizviz-dev/releases), or with
Rust installed:

```bash
cargo install --git https://github.com/aminkvh/vizviz-dev vizviz
```

The Python package (NumPy in and out, headless rendering, driving the
app): build the wheel from a checkout. See [Python](PYTHON.md).

## A first figure

```bash
vizviz --exec "fetch 4HHB; rep cartoon; color chain; addrep sticks resname HEM; lighting full; render hemoglobin.png 3840x2160 64; quit"
```

Or open a file, pick a style from the ribbon, and press Render. Every
button runs a command from the same language ([Commands](COMMANDS.md)).

## What it does

- **Big structures.** 8GLV, four million atoms, at 30 frames per second
  or better at 1080p in the atom styles, cartoon and Gaussian surface on
  a laptop RTX GPU, and on an integrated GPU for most of them.
- **Representations.** Spacefill, ball-and-stick, sticks, lines,
  solvent-accessible spheres, tube, cartoon, Gaussian, skin and
  solvent-excluded surfaces; several per structure, each with its own
  selection, colouring and material.
- **Trajectories.** DCD, XTC, TRR and AMBER NetCDF, streamed from disk.
- **Analysis.** Distances, angles and dihedrals that follow the frame,
  contacts, and solvent-accessible surface area computed as
  [FastSASA](https://github.com/aminkvh/FastSASA) does.
- **Figures.** Screenshots, SVG, and a GPU path tracer with real shadows
  and ambient occlusion.
- **Scripting.** `--exec` scripts, a Python package, and a running window
  that Python or an MCP-speaking agent can drive ([Driving the
  app](LIVE.md)).
