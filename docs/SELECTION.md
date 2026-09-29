# Selection expressions

One small language, used everywhere a selection is typed: the Selection
panel in the app, `Session.select(id, "...")` and `Structure.select("...")`
in Python, and (later) the CLI. Implemented in `vv_core::select`.

## Examples

| Expression | Selects |
|---|---|
| `chain A` | every atom of chain A |
| `chain A and name CA` | the alpha carbons of chain A |
| `resname HEM` | the heme groups |
| `resid 1-10 and chain A` | residues 1 to 10 of chain A (author numbering) |
| `within 5 of resname HEM` | every atom within 5 Å of any heme atom (heme included) |
| `byres within 4 of resname HEM and protein` | whole protein residues touching a heme |
| `protein and not backbone` | side chains |
| `bfactor > 50` | poorly ordered atoms |
| `not polymer` | everything that is not protein or nucleic |
| `membrane and within 4 of protein` | the membrane lipids touching the protein |
| `ligand or ion` | cofactors, drugs, bound lipids and metal ions |
| `lipid and ligand` | lipids bound in pockets, not part of a membrane |

## Grammar

Keywords are case-insensitive. `not` binds tighter than `and`, which
binds tighter than `or`; use parentheses to override.

| Form | Meaning |
|---|---|
| `all`, `none` | everything, nothing |
| `protein`, `nucleic`, `lipid`, `glycan`, `water`, `ion` | by molecule class, one class per residue: see "Molecule classes" below |
| `ligand`, `membrane`, `additive`, `cofactor` | by role, which can overlap the classes: see "Class vs role" below |
| `solvent`, `polymer` | `water or ion`, `protein or nucleic` |
| `hetero` | atoms from `HETATM` records (a file flag, unrelated to the classes and roles: a modified residue or an attached glycan is `hetero` but not `ligand`) |
| `backbone`, `sidechain`, `hydrogen` | by atom role |
| `helix`, `strand`, `coil` | secondary structure from the file's records |
| `chain A B` | chain ids (label or author; case-sensitive) |
| `resname ALA GLY` | residue names |
| `name CA CB` | atom names |
| `element C N` | element symbols |
| `resid 1-10 25` | author residue numbers, single or inclusive ranges (`-` or `:`) |
| `seqid 1-10` | label sequence numbers (mmCIF `label_seq_id`) |
| `serial 100-200` | atom serial numbers from the file |
| `index 0-9` | 0-based atom indices |
| `altloc A` | alternate-location ids |
| `bfactor < 20`, `occupancy >= 0.5` | numeric compare: `< <= > >= == !=` |
| `within 5 of <expr>` | atoms within 5 Å of any atom of `<expr>` (inclusive; `<expr>` is included) |
| `byres <expr>` | expand `<expr>` to whole residues |
| `not e`, `a and b`, `a or b`, `( e )` | boolean logic |

A bad expression is an error that names the offending characters, for
example `unknown keyword `foo` (at characters 8..11)`.

## Performance

`within` builds a uniform grid over the inner set and tests every atom in
parallel. On a 10M-atom structure, `within 5 of index 0-99999` takes
about 170 ms on a 16-thread laptop; the per-residue and per-chain
keywords fill whole atom ranges at once and are faster.

## Molecule classes

Every residue gets exactly one class when a structure loads (one pass,
parallel, cached as a byte per residue; nothing is recomputed per frame or
per selection): `protein`, `nucleic`, `lipid`, `glycan`, `water`, `ion`,
`small molecule` (a non-polymer molecule in none of the others: a
cofactor, a drug, a phosphate; no keyword of its own, select it with
`ligand or additive` minus lipids and glycans) or `other` (dummy atoms and
residues with no known element; no keyword). The Info tab lists the
counts, and `color class` colours by class. Two layers decide:

1. **Names.** Tables of residue names, each group citing the force field
   or dictionary it comes from (`vv_core::residue_class`): standard and
   modified amino acids with CHARMM/Amber/GROMACS protonation variants and
   `N`/`C`-terminal templates; nucleotides (DNA/RNA, terminal and
   force-field variants); water models; monatomic ions; glycans (the
   `glycan` table below); lipids from CHARMM36, Amber Lipid21/Lipid17,
   Slipids, GROMOS, Martini and the PDB Chemical Component Dictionary
   (phospholipids, sterols, fatty acids, detergents). Amber Lipid21 splits
   one lipid into a head-group residue plus one residue per tail (`PC`,
   `PA`, `OL`, ...); every piece is `lipid`. A few of those codes are also
   element symbols (`PA`, `LA`, `AR`): they are lipid only when the residue
   has more than one atom.
2. **Structure**, for a name no table knows, from the bonds of just those
   residues and their chain neighbours (perceived from the first frame, or
   an MD topology's own bond list): a lone metal or halide atom is an ion;
   one O with up to two H is water; a residue joined by a peptide bond
   (`C(i-1)-N(i)`) to protein is protein, by a phosphodiester bond
   (`O3'(i-1)-P(i)`) to nucleic is nucleic, so modified residues follow
   their chain (any other bond, such as a heme iron to a histidine, does
   not count); a residue with a run of at least 8 unbranched carbons bonded
   only to C and H **and** an acyl (ester, acid, amide) or phosphate/sulfate
   head, or a carbon/oxygen skeleton of at least 4 rings and 27 carbons
   (a sterol), is lipid; everything else is a small molecule.

A monatomic ion code (`CO`, `NI`, `CA`) is an ion only for a one-atom
residue. Limits of the structural layer: without hydrogens or bond orders
a chain cannot be checked for double bonds, and with hydrogens a
polyunsaturated tail (linoleoyl and beyond) has no 8-carbon saturated run
(oleoyl passes on its distal segment); a drug with a long alkyl chain and an
ester passes as lipid; alkyl glycosides and polyethylene glycol
detergents (ether or alcohol heads only) come out small molecules unless
their code is in the table; bile-acid-sized steroids (under 27 carbons)
are small molecules; without first-frame positions or MD bonds unknown
residues stay small molecules. `DGN` is deliberately in no table (D-glutamine in the
PDB, a DNA guanine terminus in Amber) and neither is `SPM` (spermine in
the PDB, possibly a sphingomyelin head in Amber): structure decides them.

## Class vs role

A class is chemistry (what the residue is); a role is what it does in this
structure. Every residue has one class and any number of roles, so a
palmitate in a protein pocket is a `lipid` (class) and a `ligand` (role),
while the same palmitate in a bilayer is a `lipid` and `membrane`. Roles
are decided once at load, after the classes (`vv_core::residue_class::roles`).

| Keyword | A residue has it when |
|---|---|
| `ligand` | it is a non-polymer, non-solvent molecule that is not an additive: a small molecule (drug, cofactor, unknown name), a lipid that is not `membrane`, or a glycan not covalently attached to protein or nucleic acid. A ligand bonded to the polymer through a non-peptide bond (a covalent inhibitor, a palmitoylated cysteine's lipid) is still a ligand. |
| `membrane` | it is a lipid in a contact network of at least 12 lipid residues (heavy atoms within 4 Å of each other): a bilayer patch or a micelle. Bound lipids form networks of a few. |
| `additive` | it is a small molecule whose code is a common crystallization reagent: sulfate, phosphate, glycerol, ethylene glycol, PEG fragments, MPD, acetate, formate, citrate, Tris, HEPES, MES, DMSO, ... (`role_names.rs` lists them with sources). |
| `cofactor` | it is a small molecule named as an enzyme cofactor or prosthetic group (heme variants, FAD/FMN, NAD(P), SAM/SAH, PLP, TPP, CoA, B12, iron-sulfur clusters, chlorophylls, quinones, ...). Every cofactor is also a `ligand`. |

`protein`, `nucleic`, `water`, `ion` and `other` residues have no role.
Glycan attachment uses the MD topology's bonds when the file has them, else
any two heavy atoms of different residues closer than 1.9 Å; a whole glycan
tree is attached when any of its sugars is, and a free oligosaccharide (a
bound sugar, a lectin ligand) is a ligand.

Examples: `ligand and lipid` (lipids sitting in pockets or grooves),
`ligand or additive` (every non-polymer non-solvent molecule),
`cofactor`, `membrane and within 5 of protein`.

Limits. The lists are by name only, so a glycerol that really binds in an
active site is an `additive`, and a novel cofactor under an unlisted code is
a plain `ligand`. The membrane test counts at most six contacts per lipid,
so two dense clusters of under 12 lipids each that touch through a single
atom pair can both read as bound. A nanodisc or crystal with fewer than 12
modelled lipids has none marked `membrane`. Without coordinates (a bare
topology) roles fall back to the class: every lipid is `membrane`, no
glycan is a `ligand`.

## Naming conventions

Different MD packages name residues and atoms differently; the classes
above, `backbone` and element inference (used when a file has no element
column) all recognize every convention below, so a structure prepared for
any of these tools selects the same way. The PDB reader takes residue
names of up to four characters (columns 18-21), so `POPC` and `TIP3` load
intact; the PDB writer still keeps three, so a saved membrane reloads with
`POP` and falls to the structural layer.

| Convention | `protein` | `nucleic` | `water` | `ion` |
|---|---|---|---|---|
| PDB/mmCIF CCD | the 20 standard names, `MSE`, `HYP`, `SEC`, `PYL`, `SEP`, `TPO`, `PTR`, `CSO`, `ASX`, `GLX`, `UNK` | `A C G U I` and `D`-prefixed DNA, plus common modified bases (`PSU`, `5MC`, `7MG`, ...) | `HOH`, `WAT`, `H2O`, `DOD` | element-symbol residue names (`NA`, `CL`, `ZN`, `FE`, ...) |
| CHARMM | + `HSD`/`HSE`/`HSP` (neutral delta/epsilon, doubly protonated histidine) | (PDB names) | + `TIP3` (atoms `OH2`/`H1`/`H2`) | + `SOD`, `CLA`, `POT`, `CAL` |
| Amber | + `HID`/`HIE`/`HIP` (histidine), `CYX` (disulfide cystine), `CYM` (thiolate), `ASH`/`GLH` (protonated Asp/Glu), `LYN` (deprotonated Lys); a 4-letter `N`/`C`-prefixed name whose last 3 letters are a protein code (`NALA`, `CHIS`, ...) is that residue's N-/C-terminal template | + `A5`/`A3`/`DA5`/`DA3`/... 5'/3'-terminal templates | `WAT` | + charge-suffixed `atomic_ions.lib` names (`Na+`, `Cl-`, `Mg2+`, ...) |
| GROMACS | mostly Amber's or plain PDB names | mostly Amber's or plain PDB names | `SOL` | plain element-symbol names (`NA`, `CL`, ...) |
| GLYCAM (glycoproteins) | + `NLN`, `OLS`, `OLT`, `ZOLS`, `ZOLT` (the amino acid a glycan attaches to) | - | - | - |

OPLS-AA isn't its own row: it is almost always run through GROMACS or a
GROMACS-compatible tool, and keeps GROMACS's residue/atom names rather
than its own; the GROMACS row already covers it.

`T3P`, `T4P`, `T5P` (TIP3P/TIP4P/TIP5P under the names some MD packages
give them), `SPC`, `SPCE`, `OPC` and the Martini bead `W` are water too.

`glycan` is separate: it matches a residue name from `vv_core::glycan`'s
own table (GLYCAM, CHARMM and PDB CCD monosaccharide codes, alpha and
beta forms), documented alongside the `Glycan` representation
(docs/COMMANDS.md).

`backbone` additionally recognizes CHARMM's own C-terminal oxygen names
`OT1`/`OT2` alongside PDB's `OXT`.

Element inference (`vv_core::Element::from_atom_name`, used by the
legacy PDB reader and as mmCIF's own fallback) also handles conventions
that don't fit the table above: a two-letter symbol match is trusted only
when it isn't obviously a hydrogen's locant (`HG11`, `HE21`: Val/Ile's
and Gln's own branched hydrogens read as hydrogen, not mercury or
helium), and a single-atom residue's element is read from its *residue*
name when that residue is a recognized ion (`CA`, `CD`, `SOD`, ...),
correctly telling e.g. a calcium ion from an alpha carbon even when a
file's own atom-name column alignment doesn't follow the PDB convention
precisely.

## Not yet

`same residue as`, `around` (exclusive `within`), regular expressions in
names, and selections spanning several structures.
