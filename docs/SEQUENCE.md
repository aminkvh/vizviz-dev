# Sequence coloring and tracks

The Sequence panel shows one row of one-letter codes per chain. Its header
adds two optional layers: a **color** for the letters and **tracks**, thin
rows of annotation under each chain. Both are off until you choose them.
Click a residue to select it, Ctrl+click to add, double-click to zoom;
hover for its name, number and what every visible track says about it.

## Color

`sequence color SCHEME`, or the Color menu.

| Scheme | Colors residues by |
|---|---|
| `view` | the 3D view's color of the structure's current representation |
| `ss` | secondary structure (helix, strand, turn), as the cartoon draws it |
| `chemistry` | acidic, basic, polar, nonpolar |
| `hydrophobicity` | Kyte-Doolittle hydropathy |
| `bfactor` | mean B-factor of the residue, blue to red |
| `sasa` | burial: solvent-accessible area relative to the residue's maximum (Tien et al. 2013), white exposed to blue buried. Computed in the background; skipped above 1,000,000 atoms |
| `charge` | Asp/Glu red, Lys/Arg blue, His pale blue |
| `clustal`, `zappo`, `taylor` | the usual alignment palettes |

## Tracks

`sequence track NAME on|off`, or the Tracks menu. A track only appears
under a chain that has something to show. The Key button lists the colors
of the enabled tracks.

| Track | Shows |
|---|---|
| `ss` | helix (thick bar), strand (arrow), turn (line) |
| `numbering` | author numbers every tenth residue, insertion codes kept |
| `missing` | residues of the chain's full sequence (mmCIF `_pdbx_poly_seq_scheme`, PDB `SEQRES`) with no coordinates, as an edge marker where the gap sits; hover lists them. The modeled residues are aligned to the full sequence (Needleman and Wunsch 1970, affine gaps), so insertion codes and numbering that is not ascending place correctly. Zero-occupancy residues have coordinates and are not missing. For an ensemble the gap list is the shown model's. Without the full sequence (BinaryCIF, no file) the file's unobserved list is placed by author number |
| `burial` | relative accessible area in four classes (under 10%, 10-25%, 25-50%, 50% or more exposed), computed in the background; appears when ready. Hover any residue for its exact "% exposed" |
| `conservation` | against every other protein chain of the loaded structures that is the same sequence or at least 30% identical after pairwise alignment. Shading is an entropy score, 1 minus the Shannon entropy over log2 20, with Henikoff and Henikoff (1994) sequence weights and a gap penalty, the properties Valdar (2002, Proteins 48:227) asks of a conservation score; the darkest class is a residue identical in every chain. Hover: how many chains share the residue, the score, the residues found. No network |
| `uniprot` | features of the chain's UniProt entry: domains and regions, signal and transmembrane segments, modifications, other annotated sites (cleavage, glycation; muted), and active or binding sites in red (later kinds are drawn over earlier; the Key names each color). The mapping of chain residues to UniProt positions is PDBe SIFTS; features come from the UniProt REST API. Fetched in the background and cached on disk (`uniprot` under the fetch cache), so a repeat or offline use needs no network; when nothing is cached and the network is down the header says so. Off by default and left out of `sequence tracks all`: `sequence uniprot on\|off` or the header checkbox |
| `variants` | UniProt's natural variants as their own row, off by default because a well-studied protein has hundreds (hemoglobin beta over 250). `sequence uniprot variants on\|off` or the header's variants checkbox; hover gives the substitution and UniProt's note |
| `disulfide` | bonded cysteines; partners share a color |
| `glycan` | N-X-S/T sequons (N-X-C, rare) with no glycan modeled, and Asn, Ser or Thr that carry one |
| `liability` | deamidation (NG, NS, NT), isomerization (DG, DS, DT), fragmentation (DP), integrin binding (RGD), unpaired Cys, N-terminal Gln/Glu, Met and Trp oxidation |
| `ligand` | residues within 4 Å of a bound ligand |
| `interface` | residues within 4.5 Å of another chain |
| `altloc` | residues with alternate locations |
| `modified` | non-standard amino acids |
| `antibody` | variable domains (heavy, kappa, lambda) with the numbering scheme's labels above framework and CDR bars, and a Heavy/Kappa/Lambda badge where each domain starts. Hover for `H52A · CDR-H2 (Kabat)`. See [Antibodies](ANTIBODY.md) |

Motifs never span a numbering gap. Motifs are risks to check, not
predictions: most NG sites do not deamidate, and a sequon is necessary but
not sufficient for glycosylation.

### Antibody settings

With the `antibody` track on, the header adds a "Numbers from" menu (native, or the optional external ANARCI,
`sequence antibody backend anarci`; see [Antibodies](ANTIBODY.md#external-backend); or the scheme authors' public web
service Abnum, `sequence antibody backend abnum`, Kabat/Chothia/Martin only and sent over plain HTTP, see
[Antibodies](ANTIBODY.md#scheme-authors-program-as-a-backend)) and two menus: the numbering
(`sequence antibody scheme kabat|chothia|imgt|martin|aho`; `enhancedchothia` is another name for `martin`) and the CDR
definition (`sequence antibody cdr kabat|chothia|imgt|contact|north`). They
are independent, so Kabat numbers can sit under Chothia loops. Numbers are
labeled at every tenth position and at insertion codes (`52A`), thinned
where labels would overlap; hover any residue for its own. `select cdr h3`
selects loops by the same definitions ([Selections](SELECTION.md)).

## Rows and ligands

Only polymer chains get a sequence row. Each structure's ligands, sugars
and ions are summarized in one `<structure> ligands` row of chips such as
`NAG x3`; click a chip to select every residue of that name (Ctrl+click to
add). Water is not listed.

## Chain properties

Hover a chain's label, or run `sequence props`, for its modeled residues,
mass (average residue masses), net charge at pH 7 and isoelectric point
(Henderson-Hasselbalch with the pK set of Bjellqvist et al. 1993, 1994,
including the residue-specific terminal pKs), and the extinction
coefficient at 280 nm from Trp, Tyr and cystine (Pace et al. 1995). When
Cys pairing changes it, both values are shown.

## Layout and saving

Chain and track labels stay pinned at the left while the sequence scrolls
sideways, and the label column widens to fit the longest one. Turning on
tracks grows the panel (up to 70% of its split) so they fit; dragging the
split afterwards is respected. Residue colors from categorical schemes are
soft tints of the panel color, with text chosen for WCAG AA contrast in
both themes; continuous ramps keep their full colors. Color, enabled
tracks and antibody settings are saved in sessions.

## Adding a track

A track is a `TrackProvider` (`crates/vv-app/src/sequence/tracks.rs`)
listed in `PROVIDERS`. It returns one mark per residue plus a legend, is
cached per structure (and per frame if it depends on geometry), and gets
its menu entry and `sequence track` name from its `id`. Kernels shared
with scripts live in `vv_core::seqfeat`.
