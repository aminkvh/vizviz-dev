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
| `sasa` | burial: solvent-accessible area relative to the residue's maximum (Tien et al. 2013), white exposed to blue buried. Skipped above 250,000 atoms |
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
| `missing` | residues in the file's sequence with no coordinates (`REMARK 465`, mmCIF `_pdbx_unobs_or_zero_occ_residues`), as an edge marker where the gap sits; hover lists them |
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

With the `antibody` track on, the header adds two menus: the numbering
(`sequence antibody scheme kabat|chothia|imgt|martin`) and the CDR
definition (`sequence antibody cdr kabat|chothia|imgt|contact|north`). They
are independent, so Kabat numbers can sit under Chothia loops. Numbers are
labeled at every tenth position and at insertion codes (`52A`), thinned
where labels would overlap; hover any residue for its own. `select cdr h3`
selects loops by the same definitions ([Selections](SELECTION.md)).

## Layout and saving

Chain and track labels stay pinned at the left while the sequence scrolls
sideways, and the label column widens to fit the longest one. Turning on
tracks grows the panel (up to half of its split) so they fit; dragging the
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
