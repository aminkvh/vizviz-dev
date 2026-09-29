# Changelog

## 0.2.0

**Big structures at interactive rates.** A 4-million-atom cryo-EM
structure (8GLV) draws at 30 fps or better at 1080p on a laptop RTX GPU
in the atom styles, cartoon and Gaussian surface, and on an integrated
GPU for most of them. Two-phase occlusion culling, per-residue cartoon
detail and conservative depth keep the cost off the atom count
([benchmarks](benchmarks/README.md)).

**Representations.** Several per structure, each with its own selection,
colouring and material: spacefill, ball-and-stick, sticks, lines,
solvent-accessible spheres, tube, cartoon (DSSP, nucleic-acid ladders),
Gaussian, skin and analytic solvent-excluded surfaces, and 3D-SNFG
glycan symbols (Thieker et al. 2016) detected from the structure's own
bonds. Colour by element, chain, secondary structure, residue type or
name, B-factor, occupancy, hydrophobicity, several residue-type
palettes, or any attached per-atom value channel.

**Look.** A set of materials from opaque to glass to toon, lighting
presets, up to four directional lights with colour and direction,
gradient backgrounds, screen-space and world-space ambient occlusion,
shadows, depth cue, outlines, depth of field, perspective or
orthographic projection, and a clip plane that caps every solid and
surface. Transparent materials blend order-independently
(McGuire and Bavoil 2013) on every representation.

**Bonds and selections.** Bonds come from the wwPDB Chemical Component
Dictionary's residue templates first and from distances, with valence
caps and a clash floor, only for the rest; metals bond through cofactor
templates or explicit records. Selection keywords recognise CHARMM,
Amber, GLYCAM and GROMACS residue and water naming, plus `ion` and
`glycan`.

**Formats.** PDB and mmCIF (also gzipped, or fetched by ID), CHARMM/NAMD
PSF and AMBER PRMTOP topologies, DCD, XTC, TRR and AMBER NetCDF
trajectories streamed from disk. Structures export to PDB, mmCIF, XYZ,
PQR and GRO.

**Trajectories and analysis.** Playback and scrubbing, a trajectory
attached to an open structure, distances, angles and dihedrals that
follow the frame, contacts, and solvent-accessible surface area
(Shrake-Rupley, ProtOr radii) as a per-atom value channel.

**Figures and movies.** `render` path-traces the view (real shadows and
ambient occlusion, 4K in seconds) alongside screenshots and SVG. The
Movie panel lays rotations, zooms, camera moves and trajectory playback on
a multi-track timeline and exports a PNG sequence, plus an MP4 when
`ffmpeg` is available ([Movies](docs/MOVIE.md)).

**Interface.** A ribbon of five tabs (File, Represent, Look, Analyze,
View) with key tips, a quick-access bar, dockable panels and saved
workspaces, a command palette, a Preferences dialog, a view cube and
scale bar, light and dark themes, undo for scene and view edits, and
toasts that report every change. The design rules are in
[DESIGN.md](DESIGN.md).

**Scripting.** One versioned command language ([Commands](docs/COMMANDS.md),
[Versioning](docs/VERSIONING.md)) behind the console, `--exec` scripts,
session files (`.vviz`), the Python package, and `vizviz --listen`, which
lets Python (`vizviz.connect()`) and AI agents (an MCP server,
`python -m vizviz.mcp`) drive a running window
([Driving the app](docs/LIVE.md)).

## 0.1.0

First release: mmCIF/PDB loading, spacefill and ball-and-stick on the GPU,
selections, sessions, the command language, and the Python package.
