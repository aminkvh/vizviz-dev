# File formats

All readers are written from the format specifications, not from other
tools' code (see CONTRIBUTING.md). `vv_io::load` picks the format from the
extension (`.cif`/`.mmcif`/`.pdbx`, `.bcif`, `.pdb`/`.ent`, optional `.gz`)
or from the content, memory-maps the file, and decompresses gzip on the fly.

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
- Atom names longer than four characters are kept whole
  (`Topology::atom_name`); `Topology::name` holds the first four bytes.
- `type_symbol` `D` is deuterium: the element is hydrogen, the isotope is
  the `flags::DEUTERIUM` bit, so every hydrogen rule still applies.
- The PDBx dictionary has no segment id. The local item
  `_atom_site.vizviz_segid` carries one (written only when a structure
  has segment ids, read back as the chain record's segid).
- `_entity.type`, `_entity_poly.type` and `_chem_comp.type` say which
  residues are polymer (see "Polymer status" below).
- Also read: `_entry.id`, `_struct.title`, `_struct_conf` /
  `_struct_sheet_range` (secondary structure onto residues), and
  `_struct_conn` (disulfide, covalent, metal bonds resolved to atom
  indices). Assemblies are not expanded yet.

## BinaryCIF

- Read only (`vv_io::bcif`), from the BinaryCIF specification:
  MessagePack (`vv_io::msgpack`, a small reader of the msgpack.org
  format), then each column's encodings undone last to first
  (`ByteArray`, `FixedPoint`, `IntervalQuantization`, `RunLength`,
  `Delta`, `IntegerPacking`, `StringArray`; masks give `.` and `?`).
- In the first data block with an `_atom_site` table, only the columns
  the mmCIF reader has a role for are decoded (in parallel, one task per
  column) and fed row by row, in parallel chunks, through the same
  column-to-atom mapping as the mmCIF text reader (`atom_site` module), so
  atoms, chains, segids, deuterium and models agree by construction. The
  other categories the reader consumes (entity tables, secondary
  structure, `struct_conn`, `chem_comp_bond`, the annotation categories)
  are decoded to a small CIF text and go through the mmCIF reader's own
  category code; every other category is never decoded. No text copy of
  the atoms is made: 4V6X (237,685 atoms, 12 MB `.bcif`) loads in about
  75 ms against 140 ms for its `.cif.gz` (release build, 7 runs, minimum;
  the earlier text re-emit took 400 ms).
- `fetch ID bcif` (`vv_io::fetch::fetch_bcif`) downloads the asymmetric
  unit from `models.rcsb.org/ID.bcif`.
- MMTF is not read: RCSB stopped serving it in July 2024 and recommends
  BinaryCIF instead.

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
| Segment id (73-76) | yes | kept on the chain record (`Topology::segid`, `segname` selection) and written back (a segid over four characters is cut, with one writer warning); a change of segid within one chain letter starts a new chain record; the chain name is the segid only when column 22 is blank |
| Alternate locations | yes | all conformers kept, tagged in `alt_loc`; bonds, SASA and `altloc` selection are alt-aware (same policy as mmCIF); which one draws is the `altloc` display policy (see Policies) |
| Insertion codes | yes | `(chain, seq, iCode, name)` is the residue key; `52` and `52A` are distinct residues |
| Microheterogeneity (`ASER`/`BTHR` at one number) | yes | residue name is part of the key, so each variant is its own residue |
| Element: column 77-78, else name | yes | `D` is deuterium (hydrogen with the `DEUTERIUM` flag, written back as `D`); blank falls back to name (`CA` at column 13 is calcium, ` CA ` carbon); with the column blank, deuterium is read from the name only for a lone `D` right-justified at column 14 (` D  `, ` DA `) or a filled four-character locant (`DD21`, `DG12`); `DY  `, `DUM `, ` CD ` and every name of a filled element column are unaffected (no heavy atom is aligned like that) |
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
| SEQRES | yes | polymer status only (see "Polymer status"); the sequence itself is not kept, but the writer regenerates `SEQRES` from the polymer residues (13 names per line, columns per v3.3) |
| ANISOU, MASTER, SITE, REMARK 350 | ignored | no field in the data model; skipped without error |
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
| BinaryCIF | yes | `.bcif`, decoded straight into the mmCIF reader's atom rows (see BinaryCIF) |
| MMTF | no | retired by RCSB in 2024; reported as an unknown format |
| label_* vs auth_* | yes | one rule: chain record = `label_asym_id` (`auth_asym_id` kept), `seq_id` = `label_seq_id` else `auth_seq_id` (non-polymers have none), `auth_seq_id` always kept; the PDB reader sets both from the same column, so cross-format comparisons key on the auth values |
| `label_alt_id`, `pdbx_PDB_ins_code`, `pdbx_formal_charge`, `type_symbol` | yes | as above; missing `type_symbol` falls back to the atom name |
| `pdbx_PDB_model_num` | yes | first model's number is the topology; other models with the same atom count are frames |
| `struct_conf`, `struct_sheet_range` | yes | ranges by label chain and number (insertion codes keep residues that share a number apart) |
| `struct_conn` (`disulf`, `metalc`, `covale*`) | yes | resolved to atom indices; `hydrog`, `saltbr`, `modres` are not bonds |
| `chem_comp_bond` orders | yes | non-single orders only, by atom name within each residue |
| `entity`, `entity_poly`, `chem_comp.type` | yes | decide polymer vs non-polymer per residue (see "Polymer status"); `entity` is also kept as an annotation |
| Atom names longer than 4 characters | yes | kept whole; the PDB writer truncates with a warning |
| `type_symbol` `D` | yes | deuterium flag on a hydrogen |
| Segment id | local item | `_atom_site.vizviz_segid` (the dictionary defines none); a change of segid within one `label_asym_id` starts a new chain record, as in PDB |
| Entity tables written | yes | `_entity`, `_entity_poly` (type, canonical one-letter sequence, strand ids), `_pdbx_entity_nonpoly`, `_chem_comp`; `label_entity_id` matches; see Writers |

### Policies

- Alternate locations are kept, never dropped or merged: every conformer
  is an atom with its `alt_loc` letter and occupancy, in both readers.
  Choosing a conformer to reach is a selection (`altloc A`), not a
  load-time filter, so atom counts equal the file's row counts. Which
  conformer *draws* is the per-structure `altloc` display policy
  (`altloc first|all|LABEL`, `vv_core::altloc`): `first` (the default) shows
  per residue the conformer with the largest summed occupancy, ties to the
  lowest label; untagged atoms always show; bonds and surfaces are built
  from the drawn atoms. A label a residue lacks falls back to `first`.
- Polymer status ("Polymer status" below) is a per-residue hint from the
  file; without one, or where it says nothing, the name tables decide.
- Extra models must match the first model's atom count; a model that does
  not (a truncated file) is dropped silently rather than misaligned.

### Polymer status

`Topology::polymer_hint` holds one hint per residue, set by the reader
before classes are assigned (`vv_core::residue_class`):

- **mmCIF and BinaryCIF.** A chain record inherits its entity's kind
  (`_entity.type`, `_entity_poly.type`). A polypeptide or polynucleotide
  entity makes every residue in it protein or nucleic, modified residues
  included; a non-polymer, branched or water entity makes its residues
  never protein or nucleic (a free amino acid ligand is a small molecule,
  water is water). A polymer of another type (polysaccharide, `other`)
  falls back to each residue's `_chem_comp.type`.
- **PDB.** A residue named in its chain's `SEQRES` is polymer, protein or
  nucleic by the majority of that chain's names. A `HETATM` residue is
  not polymer unless it is a modified residue the `SEQRES` names (`MSE`);
  a standard residue is written `ATOM` inside a polymer, so a `HETATM` one
  is a free ligand. Chains without `SEQRES` (MD output) give no hint.
- **Writers.** Both writers state the same thing back (see "Writers"), so
  a read, write, read cycle keeps every residue's class and, up to the
  distinction between water, non-polymer and unstated, its hint. mmCIF
  writes cut a chain record into one `label_asym_id` per run of one kind,
  so a polymer and the ligands sharing its chain record come back as
  separate chain records. A standard-residue `HETATM` in a chain absent
  from `SEQRES` is read as a free ligand.
- A `non-polymer` hint only demotes a protein or nucleic name reading;
  glycan, lipid and ion names keep their class. Residues with no hint
  keep the name-table class.

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

`cargo test -p vv-io --test polymer_roundtrip` writes 4HHB, 1AKE and 1CRN
(`.pdb` and `.cif` sources) through both writers and compares classes and
hints, plus a free amino acid ligand snippet and the written entity and
`SEQRES` text.

`cargo test -p vv-io --test bcif` reads 1CRN, 1AKE and 4HHB as `.bcif` and
checks atoms, names, coordinates to 1e-3, chains, elements, classes, hints,
residues, bonds and secondary structure against the `.cif` and `.pdb`
(`-- --ignored` adds 4V6X `.bcif` against `.cif.gz`). `--test atom_identity`
covers deuterium, long names, segids (`fixtures/small/md_segments.pdb`),
entity and `SEQRES` polymer status and the altloc policy; with `--ignored`
it uses `fixtures/real/` downloads of the neutron entries 3KCJ (H and D
mixed) and 5A93 (`.pdb`, `.cif`, `.bcif`) and 3NIR.

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
  wrapping. Deuterium is written as element `D`. The segment id fills
  columns 73-76 (cut to four characters, with one warning counting the
  chain records affected). `SEQRES` is written for each polymer chain from
  its residue sequence (13 three-character names per line; a residue
  number shared by several names counts once), before `HELIX`. An atom name over four
  characters is cut to four; when the cut collides with another name of
  the same residue and alternate location, its last character becomes a
  counter (`1`-`9`, `A`-`Z`) so names stay unique, with one warning that
  counts them. Chains are written by author chain name, so a chain's
  ligands and waters share its letter as in the PDB archive; a name
  longer than one character takes a free `A-Za-z0-9` letter (a shared
  `?` once that 62-character pool runs out), with a warning.
- **mmCIF** (`mmcif_write`): one `_atom_site` loop, `_struct_conf` /
  `_struct_sheet_range` for helices and strands, one `label_asym_id` per
  chain record (a repeated chain name gets a numeric suffix; `auth_asym_id`
  stays the letter), each requested frame
  its own `pdbx_PDB_model_num` (the original frame index + 1); explicit
  bonds as a `_struct_conn` loop. No column-width limit, so atom counts,
  chain names and atom names are lossless. Deuterium is `type_symbol` `D`;
  segment ids go in the local `_atom_site.vizviz_segid` column when any
  chain has one. Entities: consecutive residues of one kind (protein or
  nucleic polymer from the hint, else the residue class; branched for two
  or more glycans; water; other non-polymer, one run per residue name)
  form one `label_asym_id`; asyms with the same kind and monomer sequence
  share an entity id. `_entity` (`polymer`, `non-polymer`, `branched`,
  `water`), `_entity_poly` (`polypeptide(L)`, `polyribonucleotide`,
  `polydeoxyribonucleotide` or the hybrid; `pdbx_seq_one_letter_code_can`,
  parent letter for common variants, else `X`/`N`; `pdbx_strand_id` the
  author chains), `_pdbx_entity_nonpoly` and `_chem_comp` (`id`, `type`)
  precede `_atom_site`, whose `label_entity_id` matches. There is no
  BinaryCIF writer.
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
