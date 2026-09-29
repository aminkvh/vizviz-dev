# File formats

All readers are written from the format specifications, not from other
tools' code (see CONTRIBUTING.md). `vv_io::load` picks the format from the
extension (`.cif`/`.mmcif`/`.pdbx`, `.pdb`/`.ent`, optional `.gz`) or from
the content, memory-maps the file, and decompresses gzip on the fly.

`vv_io::load_topology` additionally reads `.psf` and `.prmtop`/`.parm7`
(text, optional `.gz`) as a bare `Topology` with no coordinates -- for
`loadtrajectory`'s topology argument (docs/COMMANDS.md), which pairs it
with a coordinate file or trajectory. Loading one of these directly
(`loadstructure`, `vv_io::load`) still reports "unknown format": a PSF or
PRMTOP alone is not a valid structure.

## mmCIF (PDBx)

- The `_atom_site` table is parsed in parallel: the file is split at line
  boundaries into chunks, each chunk fills its own columnar builder, and
  the builders are stitched back together (a residue split across a chunk
  seam is merged). No per-atom allocation; a hand-written decimal parser
  handles coordinates. About 8 M atoms/s on 16 threads.
- Rows are read one per line (the fast path). A row with the wrong field
  count triggers a token-based re-read that accepts rows wrapped over
  several lines; if the token count is still not a whole number of rows,
  the error names the line.
- Column order is taken from the loop header; only `Cartn_x/y/z`,
  `label_comp_id`, `label_asym_id` and an atom-name column are required.
- `label_seq_id` is the residue key; `auth_seq_id` is kept alongside.
  Chains come from `label_asym_id` (author chain kept as `auth_asym`).
- Models: the first `pdbx_PDB_model_num` defines the topology; further
  models with the same atom count become extra coordinate frames.
- Alternate locations are kept and tagged, never dropped.
- Also read: `_entry.id`, `_struct.title`, `_struct_conf` /
  `_struct_sheet_range` (secondary structure onto residues), and
  `_struct_conn` (disulfide, covalent, metal bonds resolved to atom
  indices). Assemblies are not expanded yet.

## PDB (legacy)

- Fixed columns per wwPDB format v3.3; hybrid-36 serial and residue
  numbers; `MODEL`/`ENDMDL` frames; `CONECT` explicit bonds; `HELIX` and
  `SHEET`; `HEADER` id and `TITLE`.
- Element from columns 77-78, falling back to the atom name (a two-letter
  symbol left-justified at column 13, else the first letter).
- Serials and residue numbers beyond 99999 / 9999 are read as hybrid-36
  (or plain decimal), never as garbage; a tool that wraps them modulo
  100000 / 10000 (GROMACS-style) is read as written.

## PSF (CHARMM/NAMD/X-PLOR)

- A full topology (`vv_io::psf::parse`), not just bonds: fields are
  whitespace-separated, not column-fixed, so the standard, EXT (wider
  serial/name columns) and CHEQ (two extra charge-equilibration columns)
  variants all parse the same way.
- Per atom: segid -> chain, resid -> residue (a trailing insertion
  letter, e.g. `82A`, is split off), resname, atom name, charge, mass.
  The atom "type" field (a force-field type code, e.g. `CT1`) is read
  past, not kept -- nothing in this data model consumes it.
- Element from mass (nearest standard atomic weight within 0.5 Da,
  [`Element::from_mass`]), the atom name as fallback when no tabulated
  mass is close enough (a virtual site/lone pair, or hydrogen mass
  repartitioning, which shifts mass between a heavy atom and its
  hydrogens and can defeat the mass match for both).
- Charge is a fractional force-field partial charge; `Topology::charge`
  is a small formal-charge column (as PDB uses it), so only a value
  within 0.1 e of an integer -- a monatomic ion, not a polarized
  covalent atom -- is kept, everything else reads as neutral.
- The `!NBOND` section becomes `Topology::md_bonds`, used verbatim (see
  "How bonds are decided" below).
- No coordinates: a PSF alone has no positions, so `parse` returns a
  `Topology`, not a `Structure` -- pair it with a coordinate file or
  trajectory of the same atom count via `loadtrajectory`
  (`vv_io::load_topology`).

## PRMTOP (AMBER)

- A full topology (`vv_io::prmtop::parse`) from its `%FLAG`/`%FORMAT`
  sections: a section's FORTRAN format (`20a4`, `5E16.8`, `10I8`, ...)
  gives its field width; every data line up to the next `%FLAG` is
  concatenated first, since the line-wrapping itself carries no meaning,
  then sliced into fields (a `%COMMENT` line some AmberTools versions add
  between `%FLAG` and `%FORMAT` is skipped).
- Per atom: `ATOM_NAME`, `CHARGE` (AMBER's internal units, ÷18.2223 to
  elementary charge, same 0.1 e formal-charge rule as PSF above), `MASS`,
  residue from `RESIDUE_LABEL`/`RESIDUE_POINTER`. Element from
  `ATOMIC_NUMBER` when that section is present (every AmberTools >= 12
  file), else from mass, the atom name as further fallback -- the same
  rule PSF uses.
- Chains from `ATOMS_PER_MOLECULE`: each solute molecule (before
  `SOLVENT_POINTERS`' first solvent molecule) gets its own chain, and
  every solvent molecule is lumped into one trailing chain, since a
  solvated system's thousands of individually numbered waters are not
  thousands of meaningful chains. `ATOMS_PER_MOLECULE` without
  `SOLVENT_POINTERS`: every molecule is its own chain. Neither section:
  one chain.
- Both bond sections (`BONDS_WITHOUT_HYDROGEN`, `BONDS_INC_HYDROGEN`)
  become `Topology::md_bonds`, the bond-type column discarded (vizviz
  does not model bond order).
- No coordinates, same as PSF: pair with a coordinate file or trajectory
  of the same atom count (typically the matching `.inpcrd`/`.rst7` or a
  NetCDF/DCD trajectory from the same run).

## Specification checklist (PDB v3.3 and PDBx/mmCIF)

Every row is pinned by a test in `crates/vv-io/tests/formats_audit.rs`
(hand-written snippets, so the rule is visible next to its test).
"Yes" means handled as specified.

### PDB

| Rule | Status | How |
|---|---|---|
| ATOM/HETATM columns (serial 7-11, name 13-16, altLoc 17, resName 18-20, chain 22, resSeq 23-26, iCode 27, xyz 31-54, occupancy, B, element 77-78, charge 79-80) | yes | fixed slices; a record ending before column 54 is a line-numbered error, not a panic |
| 4-character residue names (column 21 borrowed) | yes | columns 18-21 read as the name; the writer emits the same |
| Segment id (73-76) | partial | used as the chain name only when column 22 is blank |
| Alternate locations | yes | all conformers kept, tagged in `alt_loc`; bonds, SASA and `altloc` selection are alt-aware (same policy as mmCIF) |
| Insertion codes | yes | `(chain, seq, iCode, name)` is the residue key; `52` and `52A` are distinct residues |
| Microheterogeneity (`ASER`/`BTHR` at one number) | yes | residue name is part of the key, so each variant is its own residue |
| Element: column 77-78, else name | yes | `D` is hydrogen; blank falls back to name (`CA` at column 13 is calcium, ` CA ` carbon) |
| Charge `2+` / `1-` | yes | also accepts the sign-first `+2` some writers emit |
| Hybrid-36 serials and residue numbers | yes | read and written (`pdb::hybrid36`, `pdb::encode_hybrid36`) |
| TER | yes | starts a new chain record even when the letter repeats (MD systems with one letter); chain selection, chain colours and the sequence strip treat records sharing a name as one chain |
| END | yes | blocks after an `END` are extra frames when atom count and names match the first block (MD output without MODEL); a mismatching block is a second entry, so reading stops and the first is kept |
| MODEL/ENDMDL | yes | first model is the topology, later models with the same atom count are frames; others are dropped |
| CONECT | yes | repeated partner = bond order, as the old convention |
| LINK, SSBOND | yes | resolved by (chain, number, iCode, atom name); SSBOND is a disulfide between the two `SG`; LINK to a metal is a metal bond; duplicates of a CONECT bond are merged |
| HELIX, SHEET | yes | range walked in file order from first to last residue, so insertion-code residues are inside; HETATM ligands reusing a number are not |
| CRYST1 | yes | kept as `cell` / `symmetry` annotations (and written back) |
| HEADER, TITLE, COMPND, SOURCE, EXPDTA, REMARK 2, JRNL, KEYWDS, DBREF | yes | mapped onto mmCIF category names (see Annotations) |
| ANISOU, SEQRES, MASTER, SITE, REMARK 350 | ignored | no field in the data model; skipped without error |
| Assemblies (REMARK 350 BIOMT) | no | use the `fetch` assembly files |

### mmCIF

| Rule | Status | How |
|---|---|---|
| Loop columns in any order | yes | roles come from the header; optional columns default sensibly |
| Quoted values with spaces and quotes, `;` text fields | yes | quote closes only before whitespace (`O5'` and `'A B'`); text fields read in headers |
| `.` and `?` | yes | both are null everywhere |
| Rows wrapped over lines | yes | token-based re-read on a field-count mismatch |
| `_atom_site` as key/value pairs (single atom) | yes | converted to one row |
| Several data blocks | yes | the first block that has an `_atom_site` table is used |
| Gzip | yes | detected by magic bytes |
| BinaryCIF, MMTF | no | reported as an unknown format |
| label_* vs auth_* | yes | one rule: chain record = `label_asym_id` (`auth_asym_id` kept), `seq_id` = `label_seq_id` else `auth_seq_id` (non-polymers have none), `auth_seq_id` always kept; the PDB reader sets both from the same column, so cross-format comparisons key on the auth values |
| `label_alt_id`, `pdbx_PDB_ins_code`, `pdbx_formal_charge`, `type_symbol` | yes | as above; missing `type_symbol` falls back to the atom name |
| `pdbx_PDB_model_num` | yes | first model's number is the topology; other models with the same atom count are frames |
| `struct_conf`, `struct_sheet_range` | yes | ranges by label chain and number (insertion codes keep residues that share a number apart) |
| `struct_conn` (`disulf`, `metalc`, `covale*`) | yes | resolved to atom indices; `hydrog`, `saltbr`, `modres` are not bonds |
| `chem_comp_bond` orders | yes | non-single orders only, by atom name within each residue |
| `entity`, `entity_poly`, `chem_comp.type` | annotations only | polymer vs non-polymer comes from residue names (`vv_core::select` classes), not from these tables |
| Atom names longer than 4 characters | no | truncated to the 4-byte name column |

### Policies

- Alternate locations are kept, never dropped or merged: every conformer
  is an atom with its `alt_loc` letter and occupancy, in both readers.
  Choosing a conformer is a selection (`altloc A`), not a load-time filter,
  so atom counts equal the file's row counts.
- Extra models must match the first model's atom count; a model that does
  not (a truncated file) is dropped silently rather than misaligned.

### Validation

`cargo test -p vv-io --test formats_audit` runs the snippet, round-trip
(read, write, read: names, elements, altlocs, charges, coordinates,
secondary structure) and cross-format tests (1CRN, 1AKE, 4HHB as `.pdb`
and `.cif`: atoms, residues, coordinates to 1e-3, helix/strand residues).
With `-- --ignored` it also checks entries downloaded into the gitignored
`fixtures/real/` (`<ID>.pdb` and `<ID>.cif` from files.rcsb.org for 1IGY
antibody insertion codes, 3NIR alternate locations, 2K39 116 NMR models,
1OKC and 2RH1 membrane proteins with lipids, 1BNA DNA, 1ZNI; `4V6X.cif`
for the >99999-atom hybrid-36 round trip). Parsed atom and model counts
are asserted against the raw `ATOM`/`HETATM`/`ENDMDL` counts in the files.

## Writers

`vv_io::save(structure, path, &SaveOptions { atoms, frames })` picks the
format from `path`'s extension (gzip-compressed when it ends in `.gz`)
and writes it; `atoms` restricts to a subset (`None` is every atom),
`frames` are the coordinate sets to write. The `savestructure` command
(docs/COMMANDS.md) and Python's `Structure.save()` (docs/PYTHON.md) are
built on it; each format's own module (`pdb_write`, `mmcif_write`, ...)
can also be called directly, as `cargo xtask synth` does for benchmark
files.

- **PDB** (`pdb_write`): fixed columns per wwPDB v3.3. Atom name
  placement matches the reader's own heuristic (the element symbol fills
  columns 13-14, right-justified for one letter, left for two, with any
  remaining name characters in 15-16). Residue names take columns 18-21
  (four characters, as the reader accepts). `HELIX`/`SHEET` are written
  from the residues' secondary structure (a strand is a one-strand sheet).
  `TER` closes each polymer run
  (chain change or the first `HETATM` after a run of `ATOM`s). `CONECT`
  covers every explicit bond (deduplicated, written both directions);
  a disulfide additionally gets an `SSBOND` record, which the reader
  turns back into a disulfide bond. `CRYST1` is written when `_cell`'s six
  geometry values all parse. Atom serials and residue numbers wider than
  5/4 decimal digits switch to the hybrid-36 extension
  (`pdb::encode_hybrid36`, inverse of `pdb::hybrid36`) rather than
  wrapping. Chains are written by author chain name, so a chain's
  ligands and waters share its letter as in the PDB archive; a name
  longer than one character takes a free `A-Za-z0-9` letter (a shared
  `?` once that 62-character pool runs out), with a warning.
- **mmCIF** (`mmcif_write`): one `_atom_site` loop, `_struct_conf` /
  `_struct_sheet_range` for helices and strands, one `label_asym_id` per
  chain record (a repeated chain name gets a numeric suffix; `auth_asym_id`
  stays the letter), each requested frame
  its own `pdbx_PDB_model_num` (the original frame index + 1); explicit
  bonds as a `_struct_conn` loop. No column-width limit, so atom counts
  and chain names never need remapping.
- **XYZ** (`xyz_write`), **PQR** (`pqr_write`), **GRO** (`gro_write`):
  minimal, no readers in this crate (round-trip tests parse the text
  directly). PQR's charge is the file's formal integer charge, not a
  force-field partial charge; its radius is `Element::vdw_radius`, not a
  per-atom-type radius. GRO is nm rather than the rest of this crate's
  Å, with GROMACS's own residue/atom-number wraparound (modulo
  100000, not hybrid-36) and a box line from the exported atoms'
  bounding box (no unit cell is part of the data model).

## How bonds are decided

Not a file feature but runs right after loading (`vv_core::bonds::perceive`).
The problem with deciding bonds by distance alone: two atoms that are just
close together (a tight crystal packing, a compressed MD frame) look
exactly like a real bond, and a real bond in a slightly distorted or
low-resolution structure can look longer than the geometric cutoff allows.
So instead of one rule, there are four, tried in order:

1. **Templates.** A residue named ALA, HOH, DA, NAG, HEM, ... (see
   `vv_core::bonds::templates`, generated from the PDB Chemical Component
   Dictionary) bonds exactly the atom-name pairs that residue's chemistry
   says it has, not whichever atoms happen to be close. This is what
   keeps two unrelated atoms of the same residue from bonding just
   because a bad model, or thermal noise, put them near each other.
2. **Named links between residues.** A peptide bond (`C` of residue *i* to
   `N` of residue *i*+1), a nucleic backbone bond (`O3'` to `P`), a
   disulfide (`SG` to `SG`, ≤ 2.3 Å), and a glycosidic bond (a glycan's
   anomeric carbon to an acceptor oxygen, or to an N-glycosylation `ND2`)
   are the *only* ways two different residues of a polymer end up
   bonded. A distance check still applies, but only within the window
   that kind of bond actually has.
3. **Distance, for everything a template doesn't name** -- ligands,
   unknown or modified residues. Same idea as before (covalent radii sum
   plus 0.45 Å), plus two guards the old version didn't have: a per-element
   cap on simultaneous bonds (carbon 4, nitrogen 4, oxygen 2, sulfur 6,
   ...), and a floor that rejects a pair much closer than any real bond
   between those two elements (a clash, not a bond). Never applied to a
   metal atom or a lone ion: those only bond through a cofactor template
   (heme's Fe-N) or an explicit file bond, never by proximity -- this is
   what stops a metal ion from bonding a nearby water.
4. **Explicit file bonds** (`CONECT`, `_struct_conn`, `LINK`) are always
   merged in, on top of the above.

If a topology came from an MD file that already lists its own bonds
(PSF, PRMTOP; see `vv_io::psf`, `vv_io::prmtop`), perception is skipped
and those bonds are used as they are.

See docs/VALIDATION.md for how this is checked and current timings.

## Annotations (header metadata)

Both readers keep the file's annotations on `Topology::annotations`
(`vv_core::Annotations`): a list of categories, each with its item names
and rows, stored verbatim so the data model never changes when one more
field matters. mmCIF keeps `entry`, `struct`, `struct_keywords`, `exptl`,
`refine` (resolution and R factors), `cell`, `symmetry`, `citation`,
`citation_author`, `pdbx_database_status`, `em_3d_reconstruction`,
`entity`, `entity_src_gen`, `entity_src_nat`, `struct_ref`,
`pdbx_struct_assembly`, `pdbx_struct_assembly_gen`, `database_2`,
`audit_author`, `pdbx_audit_revision_history` (author lists capped at 50
rows). PDB HEADER, TITLE, COMPND, SOURCE, EXPDTA, REMARK 2, JRNL, KEYWDS
and DBREF are mapped onto the same category and item names, so
`annotations.resolution()` reads either format. A malformed header block
is ignored, never a parse error. Shown in the app's Info tab and the
`info` command; in Python as `Structure.info` / `Structure.annotations`.
Editing and writing annotations back is not implemented.
