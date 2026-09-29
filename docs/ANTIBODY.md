# Antibody variable domains

`vv_core::antibody` finds antibody variable domains (heavy, kappa, lambda,
including nanobodies) in a one-letter protein sequence, numbers them, and
labels the framework and CDR regions. Pure Rust, no I/O, no runtime data
files.

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
| `find_domains(seq) -> Vec<Domain>` | Every domain, in sequence order (an scFv gives two). |
| `Domain::{chain, start, end, score, confidence}` | `confidence` is the alignment score over a perfect framework match (0 to 1). Reported domains are at least 0.30; real ones sit at 0.42 to 0.97. |
| `Domain::numbering(Scheme) -> Vec<(usize, Label)>` | `Scheme` is `Imgt`, `Kabat`, `Chothia` or `Martin`. IMGT insertion letters stand for `.1`, `.2`, ... (`112A` is 112.1). |
| `Domain::annotate(Scheme, CdrDefinition) -> Vec<Annotation>` | Label plus `Region`, with CDR edges from `Kabat`, `Chothia`, `Imgt`, `Contact` or `North`. The two choices are independent. |
| `find_domains_with(seq, min_confidence)` | Same, with your own cutoff. |
| `find_in_residues(&Topology, range)` | Domains of a chain's protein residues; offsets are from the range start. |
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
   type covers only the framework columns in the IMGT frame (FR1 1-26,
   FR2 39-55, FR3 66-104, FR4 118-128; 93 columns). Scores are log-odds
   from clustered wwPDB sequences with BLOSUM62 pseudocounts.
3. **Alignment.** Local alignment with affine gaps inside the framework.
   Columns a germline leaves empty (heavy 10 and 73; kappa 73, 81, 82;
   lambda 10, 73, 81, 82) are free to skip. CDRs are free-length loops
   priced by a soft length prior; their residues are never scored, so any
   CDR length works. The best of the three profiles wins and must place
   Cys at 23 and 104 and Trp at 41. Leftover stretches are searched again,
   which finds the second domain of an scFv.
4. **Numbering.** IMGT is exact from the loop lengths. Kabat, Chothia and
   Martin are label rules over the IMGT frame: fixed offsets in the
   framework and, for each loop, a base run of labels, an insertion site
   and a deletion order.

### IMGT loops

Loops shorter than the maximum leave gaps at the top of the loop, with the
extra residue of an odd length on the N-terminal side; longer loops grow
symmetrically around the axis (Lefranc et al. 2003). CDR1 27-38 (axis
32|33), CDR2 56-65 (60|61), CDR3 105-117 (111|112; gaps 111, 112, 110,
113, ...; insertions 112.1, 111.1, 112.2, ...).

### Kabat, Chothia and Martin from IMGT

| Stretch | IMGT positions | Labels | Insertions after | Deletions, in order |
|---|---|---|---|---|
| Heavy FR1 | 1-9, 11-26 | 1-9, 10-25 | | |
| H1 | 27-40 | 26-35 | 35 (Kabat), 31 (Chothia, Martin) | 32 33 31 34 30 35 |
| Heavy FR2 | 41-54 | 36-49 | | |
| H2 | 55-74 | 50-65 | 52 | 54 53 55 52 56 51 |
| Heavy FR3 | 75-91, 92-94, 95-104 | 66-82, 82A-C, 83-92 | | |
| H3 | 105-117 | 93-102 | 100 | 100 99 98 97 96 95 |
| Heavy FR4 | 118-128 | 103-113 | | |
| Light FR1 | 1-23 | 1-23 (lambda: no 10) | | |
| L1 | 24-40 | 24-34 | 27 (Kabat), 30 (Chothia, Martin) | 28 29 30 27 31 |
| Light FR2 | 41-55 | 35-49 | | |
| L2 | 56-69 | 50-56 | 54 (Kabat, Chothia), 52 (Martin) | 54 53 55 52 |
| Light FR3 | 70-72, 74-80, 83-104 | 57-59, 60-66, 67-88 | | |
| L3 | 105-117 | 89-97 | 95 | 95 94 96 93 |
| Light FR4 | 118-128 | kappa 98-108; lambda 98-106, 106A, 107 | | |

H1, H2 and L1, L2 loops include the framework columns next to them
(IMGT 39-40 are Kabat H34-H35 and L33-L34). The deletion orders are the
modal choices in real depositions (see validation), not derived from a
paper; they only matter for loops shorter than the base run.

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

## Data and licences

- **Profiles** are built from framework segments (CDRs discarded) of 582
  wwPDB entries listed in `seed_ids.txt`, clustered at 90% identity into
  227 heavy, 73 kappa and 36 lambda examples in `seeds.rs`. wwPDB data is
  CC0. `gen_seeds.py` regenerates the table; it is not part of the build.
- **Nothing is copied from a germline database or another tool.** The
  profile-HMM route (Dunbar & Deane 2016) builds its HMMs from IMGT/GENE-DB
  germline alignments. IMGT data was licensed CC BY-NC-ND 4.0 with a
  financial arrangement for the private sector; reports dated July 2026
  (ELIXIR Core Data Resource announcement, a public blog summary) say it is
  now CC BY 4.0. We could not load imgt.org to confirm the wording, and we
  do not depend on it. The Dunbar & Deane paper states GPLv3 while its
  repository lists BSD-3-Clause; we read neither the code nor the HMMs.
- **Scheme rules** were re-derived from the publications below and from
  the IMGT and Honegger numbering charts, then checked against real
  depositions.

## Validation

Ground truth is the depositors' own residue numbers in wwPDB entries whose
author numbering places the four conserved anchors (Cys, Trp, Cys, FR4
Trp/Phe) at their Kabat/Chothia numbers, so the rest is comparable.
About half of the held-out domains qualify; the others are numbered 1..n or
by a private convention.

`fixtures/real/holdout/` holds 260 entries that did not seed the profiles
(RCSB queries for Fab, nanobody and Fv entries, X-ray, 3.0 Å or better,
oldest first, minus the seed list). `cargo test --release -p vv-core
--test antibody_pdb -- --ignored --nocapture` reprints these numbers.

| | Held-out | Seed entries |
|---|---|---|
| Domains found (chains scanned) | 747 (931) | 757 (965) |
| Author numbering Kabat/Chothia-framed | 390 | 228 |
| Whole domain identical, Kabat | 307 (78.7%) | 202 |
| Whole domain identical, best of Kabat/Chothia | 316 (81.0%) | 211 (92.5%) |
| Residues agreeing, best scheme per domain | 99.34% | 99.60% |
| Heavy / kappa / lambda whole-domain | 125/154, 180/209, 11/27 | 82/90, 120/122, 9/16 |

Chothia and Martin agree with each other on every entry because they only
differ at unusual loop lengths that the entries do not contain, so Martin
is untested. What the misses are, from a read of every one in the held-out
set:

- Lambda FR4: depositors split evenly between `106, 106A, 107` (used
  here, Kabat and Honegger's table) and plain `106, 107, 108`.
- CDR-H3 insertions lettered from the far end (`100J, 100K` for the last
  two inserted residues); the common convention, `100A...`, is used here
  (162 of 181 heavy chains with inserted H3 residues, both sets pooled).
- Depositors who numbered by Kabat up to some point and sequentially after
  (trastuzumab heavy chain 1N8Z is Kabat through H52, then 1..n).
- A few framework insertions and deletions numbered at a different site.

Trastuzumab reproduces the published Kabat CDRs exactly (H1 DTYIH, H2
RIYPTNGYTRYADSVKG, H3 WGGDGFYAMDY, L1 RASQDVNTAVA, L2 SASFLYS, L3
QQHYTTPPT) under `Kabat`, and the CDRs of every other definition are
asserted in `tests/antibody.rs`. In 1N8Z the light chain matches the
depositor's Kabat numbers residue for residue.

Detection over the 1,896 chains of both sets: every chain that got a
domain is an antibody fragment (checked by name; nanobodies, diabodies and
V-domain chimeras included), and the 16 chains that carry an antibody word
in their name but got none are constant domains, Fc receptors, V-set
immune receptors (TIM-3, VSIG4) and PDZ domains. Other immunoglobulin-fold
proteins (TCR alpha and beta from 16 complexes, MHC, beta-2 microglobulin,
PD-1, PD-L1, CD4, CD8, CD2, CTLA-4, Fc, titin, fibronectin) produce no
domain; the best sub-threshold hit is a TCR beta chain at 0.18, against
0.42 for the weakest antibody.

Timing, release build, one thread: 10,000 chains (1.88 M residues, 7,945
domains found) in 3.3 s, 330 µs per chain.

## Not done

- **AHo** (Honegger & Plückthun 2001) is not implemented. Its framework
  conserved C23, W43 (IMGT 41), C106 (IMGT 104) and G140 (IMGT 119) line
  up with the IMGT frame, but the loop occupancy rules (positions 28, 36, 63, 74-75,
  85-86 and 123 are the gap sites) are only in the paper's figure and would
  be guesswork. Wolfguy is not implemented either.
- **T-cell receptors** are rejected, not numbered; they need their own
  profiles (their FR3 length varies).
- **Framework insertions** beyond the germline columns (Martin's H72 and
  H8 sites) are not modeled; such chains still align but get approximate
  labels there.
- Profiles come from PDB Fabs and nanobodies (mostly human, mouse, camelid,
  some rabbit); shark VNAR and unusual germlines are untested.
- The IMGT framework gap columns for lambda FR3 are assumed equal to
  kappa's, and IMGT output itself was not compared with IMGT's own
  software.

## References

Kabat EA, Wu TT, Perry HM, Gottesman KS, Foeller C. Sequences of Proteins
of Immunological Interest, 5th ed. NIH 1991. Chothia C, Lesk AM. J Mol
Biol 196:901 (1987). Al-Lazikani B, Lesk AM, Chothia C. J Mol Biol
273:927 (1997). Abhinandan KR, Martin ACR. Mol Immunol 45:3832 (2008).
Lefranc MP et al. Dev Comp Immunol 27:55 (2003). Honegger A, Plückthun A.
J Mol Biol 309:657 (2001). MacCallum RM, Martin ACR, Thornton JM. J Mol
Biol 262:732 (1996). North B, Lehmann A, Dunbrack RL. J Mol Biol 406:228
(2011). Dunbar J, Deane CM. Bioinformatics 32:298 (2016).
