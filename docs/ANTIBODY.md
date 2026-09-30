# Antibody variable domains

`vv_core::antibody` finds antibody variable domains (heavy, kappa, lambda,
including nanobodies) and T-cell receptor alpha and beta variable domains in
a one-letter protein sequence, numbers them, and labels the framework and
CDR regions. Pure Rust, no I/O, no runtime data files.

```rust
use vv_core::antibody::{find_domains, CdrDefinition, Scheme};

for domain in find_domains(seq) {                      // empty for non-antibodies
    // domain.chain: Heavy | Kappa | Lambda; domain.start..domain.end in `seq`
    for a in domain.annotate(Scheme::Kabat, CdrDefinition::Chothia) {
        // a.index (into seq), a.label ("52A"), a.region (Fr1..Cdr3..Fr4)
    }
}
```

| Item | Meaning |
|---|---|
| `find_domains(seq) -> Vec<Domain>` | Every antibody domain, in sequence order (an scFv gives two). Receptors are not returned. |
| `find_variable_domains(seq)` | The same, plus `TcrAlpha` and `TcrBeta` domains. Each domain goes to the profile it fits best. |
| `Domain::{chain, start, end, score, confidence}` | `confidence` is the alignment score over a perfect framework match (0 to 1). Antibodies are reported from 0.30 (real ones sit at 0.42 to 0.97), receptors from 0.15. |
| `Domain::numbering(Scheme) -> Vec<(usize, Label)>` | `Scheme` is `Imgt`, `Kabat`, `Chothia`, `Martin` (enhanced Chothia, Abhinandan & Martin 2008) or `Aho`. IMGT insertion letters stand for `.1`, `.2`, ... (`112A` is 112.1). Receptors are always numbered in IMGT. |
| `Domain::annotate(Scheme, CdrDefinition) -> Vec<Annotation>` | Label plus `Region`, with CDR edges from `Kabat`, `Chothia`, `Imgt`, `Contact` or `North`. The two choices are independent; receptors use the IMGT loops. |
| `find_domains_with(seq, min_confidence)` | Same, with your own cutoff. |
| `find_in_residues(&Topology, range)` | Domains of a chain's protein residues; offsets are from the range start. |
| `chain_domains(&Topology)` | Domains of every chain, computed once per topology and cached on it (`Topology::antibody`). |
| `cdr_residues(&Topology, CdrDefinition)` | Every CDR residue of a structure, chains in parallel. |

In the app, the same detection drives the Sequence panel's `antibody` track
([Sequence](SEQUENCE.md)) and the `cdr` selection keyword (`select cdr h3`,
Kabat unless a definition is named; [Selections](SELECTION.md)). Detection
runs once per structure and again only when the numbering or definition
changes, never per frame.

`find_domains` is single threaded and allocation-light; run it from rayon
for bulk work.

## How it works

1. **Candidate window.** A V domain has two Cys 40 to 110 residues apart
   followed by an FR4 `[WF]G.G`-like motif. Sequences without that are
   rejected without alignment.
2. **Framework profiles.** One position-specific scoring profile per chain
   family covers only the framework columns in the IMGT frame (FR1 1-26,
   FR2 39-55, FR3 66-104, FR4 118-128; 93 columns). Scores are log-odds
   from clustered wwPDB sequences with BLOSUM62 pseudocounts. There are
   five profiles for antibodies and receptors: heavy, kappa, lambda, and
   receptor alpha and beta.
3. **Alignment.** Local alignment with affine gaps inside the framework.
   Columns a germline leaves empty (heavy 10 and 73; kappa 73, 81, 82;
   lambda 10, 73, 81, 82) are free to skip. CDRs are free-length loops
   priced by a soft length prior; their residues are never scored, so any
   CDR length works. The best profile wins and must place Cys at 23 and 104
   and Trp at 41. Leftover stretches are searched again, which finds the
   second domain of an scFv. Receptor profiles are tried only when no
   antibody profile reaches confidence 0.5.
4. **Numbering.** IMGT is exact from the loop lengths. Kabat, Chothia,
   Martin and AHo are label rules over the IMGT frame.

### IMGT loops

Loops shorter than the maximum leave gaps at the top of the loop, with the
extra residue of an odd length on the N-terminal side; longer loops grow
symmetrically around the axis (Lefranc et al. 2003). CDR1 27-38 (axis
32|33), CDR2 56-65 (60|61), CDR3 105-117 (111|112; gaps 111, 112, 110,
113, ...; insertions 112.1, 111.1, 112.2, ...).

### Kabat, Chothia and Martin from IMGT

| Stretch | IMGT positions | Labels | Insertions after | Deletions, in order |
|---|---|---|---|---|
| Heavy FR1 | 1-9, 11-26 | 1-9, 10-25 | (Martin: 8) | |
| H1 | 27-40 | 26-35 | 35 (Kabat), 31 (Chothia, Martin) | Kabat 35 34 33 32; Chothia, Martin 31 30 29 28 |
| Heavy FR2 | 41-54 | 36-49 | | Martin 44 43 42 |
| H2 | 55-74 | 50-65 | 52 | 53 54 55 56 57 |
| Heavy FR3 | 75-91, 92-94, 95-104 | 66-82, 82A-C, 83-92 (Martin: 66-72, 72A-C, 73-92) | | |
| H3 | 105-117 | 93-102 | 100 | 100 99 98 97 96 95 |
| Heavy FR4 | 118-128 | 103-113 | | |
| Light FR1 | 1-23 | 1-23 (Kabat, Chothia lambda: no 10; Martin lambda: no 7) | | |
| L1 | 24-40 | 24-34 | 27 (Kabat), 30 (Chothia, Martin) | Kabat 28 29 30 31; Chothia 31 32 33 34; Martin 30 29 28 27 |
| Light FR2 | 41-55 | 35-49 | (Martin: 40) | Martin 41 |
| L2 | 56-69 | 50-56 | 54 (Kabat, Chothia), 52 (Martin) | Kabat, Chothia 54 53; Martin 52 51 |
| Light FR3 | 70-72, 74-80, 83-104 | 57-59, 60-66, 67-88 (Martin: 57-88) | (Martin: 68) | Martin 68 |
| L3 | 105-117 | 89-97 | 95 | 95 94 93 92 |
| Light FR4 | 118-128 | kappa 98-108; lambda 98-106, 106A, 107 (also in Martin) | | |

H1, H2 and L1, L2 loops include the framework columns next to them
(IMGT 39-40 are Kabat H34-H35 and L33-L34).

**Martin** (Abhinandan & Martin 2008) is Chothia plus indel sites in the
framework. The published label list puts the heavy FR3 insertion at H72
(H72A-C) instead of Kabat's H82, and adds sites at H8, H42, L40A/L41 and
L68. Its effect is on almost every chain: on the held-out set 392 of 392
heavy, 61 of 61 lambda and 34 of 325 kappa domains are numbered
differently from Chothia. The order in which further residues are deleted
is not in the paper; the orders above are the ones the scheme authors'
own numbering program produced (next section), except heavy FR2, where the
program's few cases (residues 43 and 44 deleted) disagree with the
published site H42, so it deletes from 44.

**Deletion order** (which labels a loop gives up when shorter than its
base run) was measured, not derived: in structures renumbered by the
scheme authors' program, loops shorter than the base run drop labels in
one fixed order per scheme and loop, at every length seen (H3 lengths 4 to
9, L1 7 to 10, L3 5 to 8 and so on, at least three domains each). The
earlier tables were the modal choices of depositors, who mix conventions;
the measured orders raise agreement with the reference structures from
94.2% to 96.3% (Kabat) and from 91.8% to 96.2% (Chothia).

### CDR definitions

Boundaries as stated by their authors, in the numbering they were
published in:

| | H1 | H2 | H3 | L1 | L2 | L3 |
|---|---|---|---|---|---|---|
| Kabat (Kabat et al. 1991) | H31-H35B | H50-H65 | H95-H102 | L24-L34 | L50-L56 | L89-L97 |
| Chothia (Al-Lazikani et al. 1997), Chothia numbering | H26-H32 | H52-H56 | H95-H102 | L24-L34 | L50-L56 | L89-L97 |
| IMGT (Lefranc et al. 2003), IMGT numbering | 27-38 | 56-65 | 105-117 | 27-38 | 56-65 | 105-117 |
| Contact (MacCallum et al. 1996), Kabat numbering | H30-H35B | H47-H58 | H93-H101 | L30-L36 | L46-L55 | L89-L96 |
| North (North et al. 2011), Chothia numbering | H23-H35 | H50-H58 | H95-H102 | L24-L34 | L50-L56 | L89-L97 |

Insertion sites differ by scheme (Chothia & Lesk 1987; Kabat: L27A-F,
H35A-B; Chothia: L30A-F, H31A-B), which is why the same CDR spans
different residue numbers.

## AHo

AHo (Honegger & Plückthun 2001) gives every V domain positions 1 to 149,
with fixed gap sites instead of insertion letters. `Scheme::Aho` covers
antibodies (the paper's receptor columns are not modelled).

Rules, from the paper's text and Figure 1:

| Stretch | AHo positions | Rule |
|---|---|---|
| FR1 | 1-7, 9-26 | Heavy and lambda have no residue at 8; kappa fills 8. Cys at 23. |
| CDR1 outer loop | 27-30 | One gap at 28 (heavy, lambda), two at 27-28 (kappa). |
| CDR1 | 31 | One residue always sits at 31 (the hydrophobic one that packs between the sheets). |
| CDR1 inner loop | 32-40 | Gap centred on 36; W at 43. |
| FR2 | 41-57 | IMGT + 2. |
| CDR2 | 58-68 | Gap centred on 63. |
| FR3 | 69-106 | IMGT + 3 for IMGT 66-72, then + 2. Cys at 106. Light chains have no residues at 85-86. |
| CDR3 | 107-138 | The first two and the last residue sit at 107, 108, 138; the rest fill 109-137 with the gap centred on 123. |
| FR4 | 139-149 | IMGT + 21; first Gly of the FR4 bulge at 140. |

The paper does not say everything. Read from Figure 1(b), the alignment of
14 known structures, and pinned by `tests/antibody_aho.rs` (all 14 reproduce
column for column): the number of outer-loop residues; which side of a gap
takes the odd residue (right of the gap in CDR1 inner and CDR2, left in
CDR3); the split of CDR3 around 123; and the gap order in the 83-88
hairpin (86, then 85). Not in the paper or the figure, and not guessed
silently: the lambda outer-loop rule (three residues when the IMGT CDR1
has seven or more, one for six or fewer, fitted to the three lambda rows
of the figure; germlines of other lengths are untested), and every
position beyond a stretch's capacity (36A, 63A, 123A, 75A for an IMGT 73
residue), which take insertion letters.

Validation, 778 held-out domains: numbered 1 to 149 and strictly rising
(778, no insertion letters), Cys 23 (778), Trp 43 (778), Cys 106 (778),
`[WF]` at 139 (777) and Gly at 140 (777). The exceptions are the FR4
variants `WA` (1GHF heavy) and `RG` (7KN7 light). Light CDR2 leaves the
paper's eight-column gap 59-66 in 385 of 385 chains with a seven-residue
L2 (one more has nine); heavy CDR2 leaves the paper's "one to four" in 391
of 392 (the last has eight, a nanobody-like loop).

## T-cell receptors

Alpha and beta variable domains get their own profiles, built the same way
as the antibody ones: framework columns from wwPDB entries (CC0), CDRs
discarded. The seeds are the chains of 454 entries (`tcr_seed_ids.txt`,
X-ray entries whose polymer description names a T cell receptor) whose
depositors numbered the anchors as IMGT does (Cys 23, Trp 41, Cys 104):
104 alpha and 106 beta chains, 35 and 40 clusters at 90% identity
(`gen_tcr_seeds.py`). Every fourth entry in sorted-ID order, and 1AO7,
1MI5, 2BNR, 3HG1 and 5HHO, are held out: 154 entries.

Chains of one type differ in which FR3 columns they fill, so alpha has two
profiles (IMGT 69-73 empty, or 71-77 empty with 84A-84C filled) and beta
two (73 and 82 empty, or 82 only). Gamma and delta chains are not modelled:
no gamma or delta chain in the entries carries IMGT numbering. Their
variable domains come out as alpha or are missed.

Chain type comes from the constant domain (`PAVYQL` for alpha, `VFPPEVAV`
for beta) where the chain has one, since deposited names are sometimes
swapped; otherwise from the name.

| Held-out chains | Found as alpha | Found as beta | Antibody | Not found |
|---|---|---|---|---|
| Alpha (254) | 250 | 2 (1NFD) | 0 | 2 |
| Beta (253) | 0 | 249 | 0 | 4 |
| Gamma or delta (6) | 1 | 1 | 0 | 4 |
| Antibody chains of receptor entries (7) | 0 | 0 | 6 | 1 |
| Other chains of receptor entries (590) | 0 | 1 (3O9W) | 0 | 589 |
| Antibody entries (974 chains) | 2 (1NFD) | 2 (1NFD) | 778 domains | 215 chains |

Of the alpha row, the two "found as beta" are the mouse N15 beta chains of
1NFD, which the entry names "alpha-beta"; they are also the four receptor
domains among the antibody entries (1NFD holds an H57 Fab). The chain
counted under 3O9W is a chimeric receptor. The misses are two
constant-domain-only chains, one chain with Cys 104 mutated to Trp (3W0W),
and gamma or delta chimeras. No antibody is reported as a receptor, no
receptor as an antibody, and no alpha as beta. Eight beta chains of V-gene
families the seeds lack (Mel5, 5CC7) score 0.19 to 0.21, below the antibody
cutoff; the receptor cutoff of 0.15 keeps them, and nothing else scores
above 0 against a receptor profile.

Agreement with depositors is measured where they numbered all five anchors
the IMGT way (Cys 23, Trp 41, Cys 104, and Phe and Gly at 118 and 119): 46
domains in 17 entries.

| | Domains | Whole domain identical | Framework residues | Loop residues |
|---|---|---|---|---|
| Alpha | 19 | 10 | 1,536 of 1,636 (93.9%) | 416 of 455 (91.4%) |
| Beta | 27 | 13 | 2,408 of 2,410 (99.9%) | 603 of 639 (94.4%) |

Every alpha framework miss lies in FR3 (positions 65-84): the two chains
that number 84A-84C are not reproduced (the smaller profile loses to the
larger), and seven more are shifted by up to four columns. Beta misses one
residue in two chains (position 73). In the seed entries about a third of
the chains with only the three inner anchors put their CDR residues
left-justified instead of around the loop axis; those follow neither the
IMGT rule nor the majority and are left out of the comparison. Independent check: the A6 receptor
(1AO7) gets the CDRs published for it (alpha DRGSQS, IYSNGD, AVTTDSWGKLQ;
beta MNHEY, SVGAGI, ASRPGLAGGRPEQY), asserted in `tests/antibody.rs`.

## Data and licences

- **Antibody profiles** are built from framework segments (CDRs discarded)
  of 582 wwPDB entries listed in `seed_ids.txt`, clustered at 90% identity
  into 227 heavy, 73 kappa and 36 lambda examples in `seeds.rs`. wwPDB data
  is CC0. `gen_seeds.py` regenerates the table; it is not part of the build.
- **Nothing is copied from a germline database or another tool.** The
  profile-HMM route (Dunbar & Deane 2016) builds its HMMs from IMGT/GENE-DB
  germline alignments. IMGT data was licensed CC BY-NC-ND 4.0 with a
  financial arrangement for the private sector; reports dated July 2026
  (ELIXIR Core Data Resource announcement, a public blog summary) say it is
  now CC BY 4.0. We could not load imgt.org to confirm the wording, and we
  do not depend on it. The Dunbar & Deane paper states GPLv3 while its
  repository lists BSD-3-Clause; we read neither the code nor the HMMs.
- **Scheme rules** were re-derived from the publications below, the
  Honegger and Martin numbering charts, and the numbering of real
  structures, then checked against real depositions and against structures
  renumbered by the scheme authors' own program.

## Validation

`fixtures/real/holdout/` holds 260 entries that did not seed the profiles
(RCSB queries for Fab, nanobody and Fv entries, X-ray, 3.0 Å or better,
oldest first, minus the seed list; downloaded again in 2026-09, so the
counts differ slightly from earlier runs). `cargo test --release -p vv-core
--test antibody_pdb -- --ignored --nocapture` reprints these numbers.

Ground truth for the first table is the depositors' own residue numbers in
entries whose author numbering places the four conserved anchors (Cys, Trp,
Cys, FR4 Trp/Phe) at their Kabat/Chothia numbers. About half of the domains
qualify; the others are numbered 1..n or by a private convention.

| | Held-out |
|---|---|
| Domains found (chains scanned) | 778 (974) |
| Author numbering Kabat/Chothia-framed | 391 |
| Whole domain identical, Kabat | 321 (82.1%) |
| Whole domain identical, best of Kabat/Chothia | 330 (84.4%) |
| Residues agreeing, best scheme per domain | 99.53% |
| Heavy / kappa / lambda whole-domain | 127/148, 189/212, 14/31 |

No held-out heavy chain was deposited with H72A: all 159 that use an
insertion there number it H82A-C, so Martin cannot be checked against
depositors. It is checked against the scheme authors' program instead
(`tests/antibody_reference.rs`, needs `fixtures/real/reference/`, one PDB
file per Fv of the 2019 non-redundant set, seed entries excluded). Whole
domains identical:

| | Domains | Identical | Residues |
|---|---|---|---|
| Kabat | 3,625 | 3,490 (96.3%) | 99.67% |
| Chothia | 3,617 | 3,481 (96.2%) | 99.67% |
| Martin | 3,605 | 3,539 (98.2%) | 99.70% |

Martin: heavy 1,777 of 1,801, kappa 1,414 of 1,434, lambda 348 of 370.
What the misses are, from a read of the held-out depositor set:

- Lambda FR4: depositors split, 13 with `106, 106A, 107` (used here,
  Kabat's published rule) against 17 with plain `106, 107, 108`, of 31
  lambda domains with Kabat-framed numbering (one more ends at 106). Kabat's
  rule stays the default. The scheme authors' program also uses 106A.
- CDR-H3 insertions lettered from the far end (`100J, 100K` for the last
  two inserted residues); the common convention, `100A...`, is used here.
- Depositors who numbered by Kabat up to some point and sequentially after
  (trastuzumab heavy chain 1N8Z is Kabat through H52, then 1..n).
- A few framework insertions and deletions numbered at a different site.
- In the reference set, long L2 loops (lambda), a few deleted framework
  residues and domains typed differently (a kappa-like FR4 on a chain typed
  lambda) are what remains.

Trastuzumab reproduces the published Kabat CDRs exactly (H1 DTYIH, H2
RIYPTNGYTRYADSVKG, H3 WGGDGFYAMDY, L1 RASQDVNTAVA, L2 SASFLYS, L3
QQHYTTPPT) under `Kabat`, and the CDRs of every other definition are
asserted in `tests/antibody.rs`. In 1N8Z the light chain matches the
depositor's Kabat numbers residue for residue.

Detection over the 974 chains of the antibody set: 778 domains found, 215
chains without one, no antibody chain reported as a receptor (the four
receptor domains are 1NFD). Other immunoglobulin-fold proteins (MHC heavy
chains, beta-2 microglobulin) produce no domain: 589 of the 590 other
chains of the receptor entries; the exception is a chimeric receptor.

Timing, release build, one thread: 10,000 chains (1.9 M residues, 7,993
domains found) in 4.2 s, 420 µs per chain (330 µs before the receptor
profiles).

### `cdr` selection

`select cdr ...` runs on the domains cached on the topology
(`chain_domains`), not on a fresh detection. Release build, one `select cdr
h3` query on a synthetic topology (one CA per residue) built from the
deposited chains of two Fab-decorated capsids:

| Entry | Chains (residues) | No cache, 16 threads | First, 16 threads | Repeat | No cache, 1 thread | Repeat, 1 thread |
|---|---|---|---|---|---|---|
| 12US (AAV2 capsid, 60 Fab) | 180 (46,500) | 14.4 ms | 16.4 ms | 0.74 ms | 62.2 ms | 0.49 ms |
| 4UDF (parechovirus, 60 Fab) | 240 (38,580) | 7.0 ms | 7.6 ms | 0.51 ms | 30.3 ms | 0.40 ms |

The cache is per `Topology` and starts empty on a clone; edit `chains` or
`residues` in place only before the first query, or call
`AntibodyCache::clear`.

## Not done

- The AHo columns for receptors are not implemented.
- **Gamma and delta receptor chains** have no profile (see above).
- Profiles come from PDB Fabs and nanobodies (mostly human, mouse, camelid,
  some rabbit) and from human and mouse receptors; shark VNAR and unusual
  germlines are untested.
- The IMGT framework gap columns for lambda FR3 are assumed equal to
  kappa's, and IMGT output itself was not compared with IMGT's own
  software; receptor IMGT numbering is compared with depositors only.
- The two-profile split of receptor alpha is not learned well enough to
  reproduce the 84A-84C numbering.

## References

Kabat EA, Wu TT, Perry HM, Gottesman KS, Foeller C. Sequences of Proteins
of Immunological Interest, 5th ed. NIH 1991. Chothia C, Lesk AM. J Mol
Biol 196:901 (1987). Al-Lazikani B, Lesk AM, Chothia C. J Mol Biol
273:927 (1997). Abhinandan KR, Martin ACR. Mol Immunol 45:3832 (2008).
Lefranc MP et al. Dev Comp Immunol 27:55 (2003). Honegger A, Plückthun A.
J Mol Biol 309:657 (2001). MacCallum RM, Martin ACR, Thornton JM. J Mol
Biol 262:732 (1996). North B, Lehmann A, Dunbrack RL. J Mol Biol 406:228
(2011). Dunbar J, Deane CM. Bioinformatics 32:298 (2016). Garboczi DN et
al. Nature 384:134 (1996) (the A6 receptor).
