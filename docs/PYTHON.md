
`import vizviz` gives you the same parsers, data model, and command bus
the desktop app uses, with NumPy arrays that point straight at the Rust
memory. Aimed at anyone who can write a NumPy script.

## Install

No PyPI release yet; build a wheel from a checkout (needs `uv` and a
Rust toolchain) and install it into whatever environment you like:

```bash
uv run --project crates/vv-py maturin build --release --manifest-path crates/vv-py/Cargo.toml --out dist
pip install dist/vizviz-*.whl        # one abi3 wheel covers Python 3.10+
```

For hacking on the bindings themselves, `uv sync --project crates/vv-py`
builds them in place and `uv run --project crates/vv-py pytest` tests them.

## Quickstart

```python
import numpy as np
import vizviz

s = vizviz.load("fixtures/small/4HHB.cif")
print(s)  # Structure(id="4HHB", atoms=4779, residues=801, chains=14, frames=1)

ca = np.char.strip(s.name) == b"CA"  # per-atom boolean mask
chain_a = s.residue_chain[s.residue_index] == 0  # residue table -> per atom
centroid = s.positions[ca & chain_a].mean(axis=0)

img = s.render(1200, 900, coloring="chain", style="publication_white")
vizviz.write_png("4hhb.png", img)
```

## `Structure`

`vizviz.load(path)` reads mmCIF or PDB (optionally `.gz`) and returns a
`Structure`. `vizviz.fetch("4HHB")` downloads an entry from RCSB into the
per-user cache and returns its path (`assembly=1` for the biological
assembly); `vizviz.load_many(paths)` parses many files in parallel. Its per-atom columns are **zero-copy, read-only NumPy
views** of the Rust data:

| Attribute | dtype, shape | Meaning |
|---|---|---|
| `positions` | float32 `(n, 3)` | coordinates of frame 0, in Å |
| `frame(i)` | float32 `(n, 3)` | coordinates of frame `i` |
| `element` | uint8 `(n,)` | atomic number; 0 = unknown |
| `name` | `S4` `(n,)` | atom name as in the file, space padded |
| `serial` | uint32 `(n,)` | atom serial number |
| `b_factor`, `occupancy` | float32 `(n,)` | |
| `alt_loc` | uint8 `(n,)` | ASCII alt-loc id; 0 = none |
| `charge` | int8 `(n,)` | formal charge |
| `flags` | uint8 `(n,)` | bit 0 = HETATM |
| `residue_index` | uint32 `(n,)` | residue of each atom |

Rules for views:

- They are read-only (`arr.flags.writeable` is `False`); `np.array(arr)`
  makes an editable copy.
- They keep the `Structure` alive (`arr.base is s`), so dropping `s` never
  invalidates an array you still hold.
- Each access returns a new array object over the same memory.

Per-residue and per-chain tables are small copies:

| Attribute | dtype, shape |
|---|---|
| `residue_atoms` | uint32 `(r, 2)`, `[start, end)` atom range |
| `residue_chain` | uint32 `(r,)` |
| `residue_seq_id`, `residue_auth_seq_id` | int32 `(r,)` |
| `residue_ins_code` | uint8 `(r,)`, ASCII; 0 = none |
| `residue_ss` | uint8 `(r,)`: 0 unknown, 1 coil, 2 helix, 3 strand |
| `residue_comp` | uint32 `(r,)`, index into `names` |
| `residue_names` | `list[str]` |
| `chain_residues` | uint32 `(c, 2)`, `[start, end)` residue range |
| `chain_names`, `chain_auth_names` | `list[str]` |
| `chain_entity` | uint16 `(c,)` |
| `names` | `list[str]`, the interned string table |

Also: `atom_count`, `residue_count`, `chain_count`, `frame_count`, `id`,
`title`, `atom_name(i)`, `residue_name(i)`, `chain_name(i)`,
`distance(a, b, frame=0)`, `bonds()` (perceived, uint32 `(k, 2)`),
`explicit_bonds`, `explicit_bond_kinds`.

`info` is a dict of the usual header fields (title, method, resolution,
deposition_date, organism, keywords, citation_title, doi, uniprot,
entities); `annotations` is everything the parser kept, as
`{category: [row, ...]}` in mmCIF vocabulary. PDB headers are mapped onto
the same names, so `s.info["resolution"]` works for both formats.

`select(expr)` returns the uint32 indices matching a selection expression
(grammar in [SELECTION.md](SELECTION.md)):

```python
ca = s.select("chain A and name CA")
pocket = s.select("byres within 4 of resname HEM and protein")
```

`save(path, selection=None, frames=None)` writes a structure file, format
from `path`'s extension (`.pdb`/`.ent`, `.cif`/`.mmcif`/`.pdbx`, `.xyz`,
`.pqr`, `.gro`; any of them optionally `.gz`). `selection` is a selection
expression or a sequence of atom indices (default: every atom); `frames`
are coordinate-set indices to write (default: every frame). Returns a
list of warnings, e.g. a PDB chain id that had to be remapped to one
character:

```python
s.save("pocket.pdb", selection="byres within 4 of resname HEM")
s.save("model0.cif", frames=[0])
```

## Analysis

Batch-first kernels (see [ANALYSIS.md](ANALYSIS.md) for the design and
measured speed). Multi-frame results have one row per frame; groups are
selection expressions or index sequences.

```python
s.angle(n, ca, c)  # degrees
s.dihedral(c_prev, n, ca, c)  # degrees, IUPAC sign
s.dihedrals(quads)  # (frames, k) float32
pairs, dist = s.contacts("chain A", "chain B", 4.0)
counts = s.contact_counts("chain A", "chain B", 4.0)  # (frames,) uint32
pairs, dist = s.neighbors(4.0, "protein and not hydrogen")
s.residue_pairs(pairs)  # (m, 2) residue indices
structures = vizviz.load_many(paths)  # parallel parse, GIL released
```

## Bring your own arrays

`vizviz.from_arrays(positions, *, names, elements, residue_names,
residue_ids, chain_ids, auth_chain_ids, b_factors, occupancies, hetero,
id, title)` builds a `Structure` from per-atom NumPy columns; only
`positions` (`(atoms, 3)` or `(frames, atoms, 3)`, Å) is required.
Elements are guessed from atom names when omitted. Residues are runs of
the same chain + residue id + residue name, in the order given. 2M atoms
take about a second.

`vizviz.from_mdanalysis(universe_or_atomgroup, frames=None)` and
`vizviz.to_mdanalysis(structure)` move structures across without a file
(`frames="all"` or a slice pulls a trajectory). vizviz never imports
MDAnalysis: the first reads plain attributes off whatever you pass in,
the second imports the MDAnalysis you installed yourself, when called.
That keeps the license boundary in docs/ENVIRONMENT.md intact.

## Piping through another package

The contract is "NumPy in, NumPy out": vizviz hands over positions and
names, the other package computes, and one number per atom comes back as
a **value channel** that both renderers color by. Nothing in the other
package has to know vizviz exists. With [FastSASA](https://github.com/aminkvh/FastSASA),
unchanged (the full script is `examples/fastsasa_pipe.py`):

```python
import fastsasa, fastsasa_adapters as fa, numpy as np, vizviz

s = vizviz.load(vizviz.fetch("4HHB"))
names = np.char.strip(s.name.astype("U4"))
resn = np.asarray(s.residue_names)[s.residue_index]
radii = fa.load_radius_config().radii(resn, names, default=1.7)  # FastSASA's own table
total, atom = fastsasa.sasa(s.positions.astype(np.float64), radii, atom_sasa=True)

img = s.render(1200, 900, values=atom[0])  # color by any per-atom number
np.save("sasa.npy", atom[0])  # or hand it to the app:
#   vizviz --exec "fetch 4HHB; values sasa sasa.npy; color values sasa"
```

In a `Session`, `set_values(id, "sasa", atom)` attaches the channel and
`set_coloring(id, "values:sasa")` (or `exec("color values sasa")`)
colors by it. One value per residue is expanded to the residue's atoms
(what a residue SASA or conservation score needs, and what reads well on
a tube); a `(frames, n)` array gives a per-frame channel with one color
scale over the whole trajectory. `save_session` writes each channel as a
`.npy` sidecar next to the session file, so `color values sasa` survives
a save/load round trip.

## `Session`: the command bus

Every click in the desktop app is a command; a `Session` dispatches the
same commands, with undo/redo.

```python
sess = vizviz.Session()
sid = sess.load("fixtures/small/4HHB.cif")
sess.select(sid, np.flatnonzero(sess.structure(sid).element == 26))  # the 4 Fe
sess.save_selection_set("iron")
sess.set_coloring(sid, "b_factor")
sess.undo()
sess.redo()
img = sess.render(800, 600)
```

| Method or property | Notes |
|---|---|
| `load(path) -> id`, `close(id)`, `structures` | ids are ints, never reused within a session |
| `load_trajectory(topology, dcd) -> id` | a topology file plus a DCD trajectory as one multi-frame structure |
| `structure(id)` | a `Structure` sharing memory with the session |
| `label(id)`, `path(id)` | |
| `representation(id)`, `set_representation(id, r)` | `"spacefill"`, `"ball_and_stick"`, `"tube"`, `"cartoon"`, `"gaussian_surface"`, or `"skin_surface"`; `Session.render()` (the CPU backend) draws `spacefill`/`gaussian_surface`/`skin_surface` as themselves, and `ball_and_stick`/`tube`/`cartoon` as spacefill (a warning) — no bond cylinders or ribbon triangles in this backend |
| `frame(id)`, `set_frame(id, n)` | the current coordinate-set index (0-based) of a multi-frame structure; `set_frame` is undoable and raises `IndexError` if `n` is out of range. This is the frame `save_session` persists |
| `coloring(id)`, `set_coloring(id, c)` | `"element"`, `"chain"`, `"b_factor"`, or `"values:NAME"` |
| `set_values(id, name, array)`, `values(id)`, `get_values(id, name)`, `remove_values(id, name)` | value channels: `(atoms,)` or `(residues,)`, or `(frames, n)` of either; undoable; persisted as `.npy` sidecars by `save_session` |
| `select(id, atoms)`, `clear_selection()`, `selection`, `selection_expr` | `atoms`: an expression string, a bool mask, or an index sequence; `selection` is `(id, indices)` or `None`; `selection_expr` is the expression it came from, if any |
| `save_selection_set(name)`, `delete_selection_set(name)`, `selection_sets`, `selection_set(name)`, `selection_set_expr(name)` | a set remembers its expression |
| `undo()`, `redo()`, `can_undo`, `can_redo`, `version` | `version` bumps on every applied command |
| `exec(script) -> list[str]`, `Session.commands()` | text commands, one per line or `;`-separated ([COMMANDS.md](COMMANDS.md)); `commands()` lists them as `(id, usage, help)` |
| `save_session(path)`, `load_session(path) -> list[str]` | session files shared with the desktop app (structure and trajectory paths, selections as expressions, value channels as `.npy` sidecars, each structure's frame — no structure atom data); `load_session` returns warnings for anything it could not restore |

Errors: unknown id or set name → `KeyError`; bad option strings →
`ValueError`; file problems → `FileNotFoundError` / `OSError` /
`ValueError`. A failed command never enters the history.

## Headless rendering

`Structure.render(width, height, *, coloring, values, colors, range,
style, yaw, pitch, zoom, background)` and `Session.render(width, height,
*, id, style, ...)` return an `(height, width, 4)` uint8 RGBA array;
`vizviz.write_png(path, img)` saves it. `values` (one number per atom
or per residue) colors on the blue-white-red ramp over `range` or the
column's own min..max, exactly as the B-factor coloring does; `colors`
is an `(atoms, 3)` or `(atoms, 4)` uint8 array used as is.

- Uses the CPU backend (`vv-cpu`), so it works with no GPU at all (HPC
  login nodes, CI). Spacefill only for now; a ball-and-stick structure
  renders as spacefill with a `UserWarning`.
- The camera auto-frames the structure; `yaw`/`pitch` are radians on top
  of that, `zoom > 1` moves the camera away.
- Styles: `dark_presentation` (`dark`), `publication_white` (`white`),
  `glossy`, `flat_cel` (`flat`). `background=(r, g, b, a)` overrides the
  style's own.
- `vizviz.HAS_RENDER` is `False` if the package was built without the
  `render` Cargo feature.

## Not yet

Trajectory file readers, GPU rendering in a headless `Session`, a
`to_mdtraj` twin. See ROADMAP.md at the repository root. (Driving the live
window from Python, GPU rendering included, is `vizviz.connect()`:
[LIVE.md](LIVE.md).)
