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
| Heavy FR1 | 1-9, 11-26 | 1-9, 10-25 | 6 (Martin 7), only for a chain too long for the stretch | Residues missing inside the chain: 10 9 11 8 (all three schemes) |
| H1 | 27-40 | 26-35 | 35 (Kabat), 31 (Chothia, Martin) | Kabat 35 34 33 32; Chothia, Martin 31 30 29 28 |
| Heavy FR2 | 41-54 | 36-49 | | One missing residue: where the aligner left a column empty; two or more: 44 43 42 |
| H2 | 55-74 | 50-65 | 52 | 53 54 55 56 57 |
| Heavy FR3 | 75-104 | 66-92 and insertions, counted | 82 (Martin 72) | Where the aligner left columns empty |
| H3 | 105-117 | 93-102 | 100 | 100 99 98 97 96 95 |
| Heavy FR4 | 118-128 | 103-113 | | From the end (truncated FR4) |
| Light FR1 | 1-23 | 1-23 | | At most 1-4 missing at the start, then Kabat, Chothia 10 9 8 7 6 5; Martin 7 6 5 8 9 10 (lambda always lacks one: 10, Martin 7) |
| L1 | 24-40 | 24-34 | 27 (Kabat), 30 (Chothia, Martin) | Kabat 28 29 30 31; Chothia 31 32 33 34; Martin 30 29 28 27 |
| Light FR2 | 41-55 | 35-49 | (Martin: 40) | Martin 41 |
| L2 | 56-69 | 50-56 | 54 (Kabat, Chothia), 52 (Martin) | Kabat, Chothia 54 53; Martin 52 51 |
| Light FR3 | 70-104 | 57-88 and insertions, counted | 66 (Martin 68) | Where the aligner left columns empty; Martin, failing that, 68 |
| L3 | 105-117 | 89-97 | 95 | 95 94 93 92 |
| Light FR4 | 118-128 | 98-108; with a lambda-type J segment 98-106, 106A, 107 (also in Martin) | | From the end |

H1, H2 and L1, L2 loops include the framework columns next to them
(IMGT 39-40 are Kabat H34-H35 and L33-L34).

Frameworks are counted, not indexed by IMGT column, the way the scheme
authors' program numbers them: a heavy FR3 of 27 residues has no 82A, one
of 28 has 82A, and so on to 82E; the light FR3 gets 66A-C likewise. Where
residues are missing, the aligner says which columns are empty (FR2, FR3).
In FR1 it says only where the chain starts: one or two residues missing
inside FR1 go to the scheme's own site whatever the sequence says. A light
chain loses at most positions 1-4 at the start (a shorter stub still starts
at label 5 in the reference set, `LIGHT_MAX_LEADING`) and the rest inside
FR1; a heavy chain has no such limit. A residue in front of the aligned
start belongs to the domain when fewer than 23 residues precede the first
conserved Cys: the chain then holds nothing but framework there, and the
reference numbers every one of them (`claim_leading`). A heavy chain with
more framework residues than the stretch has positions gets its extra one
at 6A (Martin 7A), through a spare column that costs 3 bits. FR4 is
cut from its end. A lambda V domain followed by a kappa-type J segment (Lys
or Arg at IMGT 127, as in `...TKLEIKR`) gets no 106A, as in the reference
set. In the aligner, a real gap right after a free column pays the full
opening price, and kappa and lambda (after IMGT 82) and heavy (after 94)
have three spare insertion columns that cost little; without them a long
framework spills into the neighbouring CDR.

**Martin** (Abhinandan & Martin 2008) is Chothia plus indel sites in the
framework. The published label list puts the heavy FR3 insertion at H72
(H72A-C) instead of Kabat's H82, and adds sites at H8, H42, L40A/L41 and
L68. The reference program's output shows no H8 deletion in FR1 (it
deletes at 10 like the others), and its heavy FR1 insertions sit at 7A
(Kabat and Chothia put them at 6A). Its effect is on almost every chain: on the held-out set 392 of 392
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

`fixtures/real/holdout/` holds entries that did not seed the profiles (RCSB
queries for Fab, nanobody and Fv entries, X-ray, 3.0 Å or better, oldest
first, minus the seed list). The directory was rebuilt on 2026-09-30 with
the same queries, since the earlier copy and its id list were gone: 257
entries now, against 260 then. The first table is the earlier download (778
domains, which the AHo, receptor and detection figures below also use); the
second is the current one. `cargo test --release -p vv-core --test
antibody_pdb -- --ignored --nocapture` reprints the second.

Ground truth for both tables is the depositors' own residue numbers in
entries whose author numbering places the four conserved anchors (Cys, Trp,
Cys, FR4 Trp/Phe) at their Kabat/Chothia numbers. About a third to a half
of the domains qualify; the others are numbered 1..n or by a private
convention.

| | Earlier download |
|---|---|
| Domains found (chains scanned) | 778 (974) |
| Author numbering Kabat/Chothia-framed | 391 |
| Whole domain identical, Kabat | 321 (82.1%) |
| Whole domain identical, best of Kabat/Chothia | 330 (84.4%) |
| Residues agreeing, best scheme per domain | 99.53% |
| Heavy / kappa / lambda whole-domain | 127/148, 189/212, 14/31 |

| | Current download, before the reference work | Now |
|---|---|---|
| Domains found (chains scanned) | 580 (849) | 580 (849) |
| Author numbering Kabat/Chothia-framed | 198 | 198 |
| Whole domain identical, Kabat | 182 | 179 |
| Whole domain identical, best of Kabat/Chothia | 182 | 179 |
| Residues agreeing, best scheme per domain | 99.76% | 99.58% |
| Heavy / kappa whole-domain | 50/61, 132/137 | 47/61, 132/137 |

No framed lambda domain is in the current download. The three domains that
no longer agree are explained under the reference comparison below.

No held-out heavy chain was deposited with H72A: all 159 that use an
insertion there number it H82A-C, so Martin cannot be checked against
depositors. It is checked against the scheme authors' program instead.

### Agreement with the scheme authors' program

`tests/antibody_reference.rs` and `tests/antibody_reference_exact.rs` read
`fixtures/real/reference/{kabat,chothia,martin}/*.pdb`, one PDB file per Fv
of the non-redundant set. The files come from the antibody structure
database of the Martin group (Ferdous & Martin 2018, Database 2018:bay040),
snapshot of 2019-07-26. Re-download:

```
base=http://www.abybank.org/abdb/Data
for s in Kabat Chothia Martin; do
  curl -O $base/NR_LH_Combined_$s.tar.bz2     # about 78 MB each
  tar xjf NR_LH_Combined_$s.tar.bz2           # NR_LH_Combined_$s/<ID>_<n>.pdb
done
```

Put the three directories under `fixtures/real/reference/` as `kabat`,
`chothia` and `martin`, or point `VIZVIZ_FIXTURES_REAL` at another
`fixtures/real`. The bundles hold 1,934, 1,930 and 1,924 files; 3,625,
3,617 and 3,605 domains remain after the seed entries are dropped. The
`NR_LH_Protein_*` bundles are a subset and were not used. The site now
also offers newer snapshots as ZIP files; these numbers are for the 2019
bundles. The scheme authors' program itself is not needed: the files are
its output.

Whole domains identical (`cargo test --release -p vv-core --test
antibody_reference -- --ignored --nocapture` prints the same counts):

| | Domains | Identical before | Identical now | Residues now |
|---|---|---|---|---|
| Kabat | 3,625 | 3,490 (96.3%) | 3,605 (99.4%) | 99.88% |
| Chothia | 3,617 | 3,481 (96.2%) | 3,596 (99.4%) | 99.88% |
| Martin | 3,605 | 3,539 (98.2%) | 3,585 (99.4%) | 99.88% |

Three sets are reported (`antibody_reference_exact.rs`). The **raw** set is
the table above. The **comparison set** drops what the reference itself did
not scheme-number or numbered inconsistently: the chains its `REMARK 950`
lines type `A` (antigen; deposited numbers kept) and every copy of a domain
that it numbers differently from an identical one (same chain sequence and
domain sequence, different labels; the majority numbering is kept, and with
no majority all copies go). The reference files carry only renumbered
labels, so "identical" means identical sequence, not identical deposited
numbers. 14 chains are typed antigen; no two identical chains are numbered
differently, so nothing else is dropped. (An earlier reading of the files
took one residue per CA atom and so skipped residues built without one:
4XCF_1:H has a Leu 99 with only its C and O, 1MFE_1:H a Gly 42 with only
its N, 1QFW_1:L a residue 27 without a CA. Dropping them shifted every
label after and made those chains look like reference contradictions or
loop anomalies; read whole, they follow the rules.) Across the roughly
1,800 non-antigen heavy domains of each scheme the labels of H3 depend on
its length alone, one pattern for each length from 4 to 30 in every scheme
(`reference_loops_follow_one_pattern_per_length`), with no exception, and
the same holds for H1, L1 and L3 (kappa and lambda). Our loops reproduce
the pattern at every length. Both halves split the sorted non-seed
PDB ids of the three bundles together, even index and odd index.

| | Raw | Comparison | Comparison, even half | Comparison, odd half |
|---|---|---|---|---|
| Kabat, before | 3,602 / 3,625 | 3,595 / 3,609 | 1,809 / 1,818 | 1,786 / 1,791 |
| Kabat, now | 3,608 / 3,625 | 3,602 / 3,611 | 1,811 / 1,818 | 1,791 / 1,793 |
| Chothia, before | 3,593 / 3,617 | 3,587 / 3,601 | 1,807 / 1,816 | 1,780 / 1,785 |
| Chothia, now | 3,599 / 3,617 | 3,594 / 3,603 | 1,809 / 1,816 | 1,785 / 1,787 |
| Martin, before | 3,582 / 3,605 | 3,577 / 3,589 | 1,801 / 1,808 | 1,776 / 1,781 |
| Martin, now | 3,588 / 3,605 | 3,584 / 3,591 | 1,803 / 1,808 | 1,781 / 1,783 |

The comparison set is not at 100%: 11 domains (Kabat, Chothia) or 9
(Martin) remain, listed in `RESIDUAL` of the exact test, which fails on any
other difference and on a listed one that starts to agree. The held-out
depositor report did not move (99.58% of residues, Kabat/Chothia-framed
domains, heavy 47 of 61 whole-domain identical).

Mismatching domains before and how many agree now, by the region of the
first difference (Kabat; Chothia within one domain of it):

| Region | Kabat | Martin | What fixed it |
|---|---|---|---|
| Heavy FR3 | 43, 41 agree | 3, 0 | Count the residues; insert after 82 (Martin 72) |
| Light FR3 | 35, 31 | 8, 4 | Count the residues, insert at 66; no cheap gap after a free column |
| Lambda FR4 | 9, 9 | 9, 9 | 106A only with a lambda-type J |
| Lambda L2 (FR3 of 35 residues) | 11, 11 | 9, 9 | Spare columns after IMGT 82 |
| Light FR1 | 19, 9 | 13, 9 | Canonical gap site |
| Heavy FR1 | 6, 0 | 8, 2 | Canonical gap site |
| Heavy FR2 | 3, 2 | 5, 3 | One gap by the aligner, several from 44 |
| Heavy H2 | 6, 1 | 6, 1 | Spare columns after IMGT 94 |
| Heavy FR4 | 0, 0 | 1, 1 | Cut from the end |
| L1, H3 | 3, 0 | 4, 0 | |
| Total | 135, 104 | 66, 39 | |

### Placement by consensus

The rules above fix one labeling per domain. The domains still wrong were
all framework indel *placements* among near-equivalent alternatives, so a
placement step (`placement.rs`) lets alternatives compete on similarity to a
consensus of the reference numbering itself.

**Consensus.** Per scheme and label, residue counts over the heavy chains of
the even half of the comparison set (`consensus.rs`, written by
`tests/antibody_consensus_gen.rs`; labels 36-49 in every scheme, and 57-92
with their letters in Martin's only). A residue at a label scores log2 of its frequency there, with
background-weighted pseudocounts (weight 5), over its background frequency;
labels without data score zero. Light chains and other labels have no table.

**Candidates and prices.** Candidates are labelings around the rule's. Each
label used or skipped unlike the rule costs bits, and so does each insertion
letter more or fewer; a candidate must beat the rule by more than float noise,
so the rule is the tie-break. Two places are open, everywhere else the price
is prohibitive: FR2 (36-49), where the block of skipped labels may slide
(4 bits per label changed, one contiguous block so the residues keep their
order), and the CDR-H2 tail to FR3 (50-92), a dynamic program over labels 50-65
skippable at 30 bits and letters at the H2 and FR3 sites movable at 8 bits per
letter (at most 6 above the rule's count). That split is decided once, in
Martin's frame (FR3 letters after 72), and the result relabeled to Kabat's or
Chothia's 82: the split is a property of the residues, and the schemes
disagreed on it for 4YDL_1:H (H2 of 19 residues and FR3 of 30 in Kabat and
Chothia, 17 and 32 in Martin and in the reference), each running the program
on its own letter rows. It is the only one of 1,815 heavy domains whose H2/FR3
residue counts differ between schemes; no other domain moved.

**How it was chosen.** A first version opened every label and letter under one
price. At every price either it lost whole domains (the reference keeps the
rules' gap site in FR3, FR1 and FR4 whatever the sequence says: tens to
hundreds more errors) or fixed nothing; where it fixed a domain the free optimum often was
not the reference's. Opening only FR2 and the H2/FR3 boundary and shutting the
rest gave no loss at 4 bits per label in FR2 (2 lost three domains), 20 to 50 in
the tail, and the same result with the letter price from 4 to 8; the pseudocount weight was the one setting
that mattered (1 to 5 fixed the Martin boundary case, 10 did not; 30 fixed
nothing), taken from the middle of the plateau. Lambda FR1 start placements were
tried as a third family (leading run plus one block); it either lost kappa and
lambda stubs or left the three lambda cases unchanged and was dropped.

**Result.** Fixed: 4LLV_3:H (FR2 block, G at 44) and 1QFW_1:H (H2 tail: 65
skipped, Lys at 66) in all three schemes, and 4YDL_1:H in all three (Kabat and
Chothia through the shared split above). The even half gained 4LLV_3:H (all
schemes) and 4YDL_1:H, the odd half 1QFW_1:H; the consensus comes from the even
half only, so the odd half's gain is validation. No domain that agreed before
disagrees now.

**Residual domains.** `VIZVIZ_EVIDENCE=1` prints, for each one, the rule's,
our and the reference's labels over the differing stretch and every placement
of its residues onto the base labels between the stretch's ends, scored by the
even-half consensus (bits). In every enumerated case the reference's labeling
is among the placements; the question is where it ranks.

| Domain | Site | Reference rank among placements | Consensus |
|---|---|---|---|
| 4LLV_3:L, 3UTZ_1:L, 5EOC_2:L, 5VTA_2:L, 6BPC_1:L | kappa FR3 gap site | 2 of 10 (11.0 vs 11.1), 2 of 2, 2 of 2, 2 of 3, 3 of 3 | **Rule not found.** Not a tie: no other chain holds the window (`FSGSGTDFTL`, `RFSGSGGTDF`, `TLNIPVEEEDAA`, `LTRVEAEDAA`, `DLAYFC`). Of the 12 comparison-set kappa domains with a missing FR3 residue the reference's gap sits at 66, 66-67, 68 (three), 74-75, 77 (three), 83 and 82-85; we agree in 7. The consensus prefers the rule by 4 to 14 bits in four of the five, and ties 4LLV_3:L |
| 5CEY_1:L, 6NNJ_1:L | lambda start | 1 of 5 (-4.8 vs -12.5; Martin 5CEY by 0.1 bit) | **Rule not found.** The reference is consistent (the same `YVRPLSVALG` is Y4 in both), so ours is the miss. `SYVRPLSVALG` (5CEY_2:L) is S5, `VRPLSVALG` (4FQ2_1:L) is V5, `YVSPLSVALG` (4R26_1:L) is Y5: no fixed start or leading count reproduces the four. A consensus over starts reproduces this pair and loses stubs (see above) |
| 1OAY_2:L | lambda start | 9 of 36 (6.4 vs 14.2 best) | **Rule not found, fault in our start.** `AVVTQESALTT` siblings (1MFE_1:L, 1OAX_1:L, 1IND_1:L) start at 2 and leave 1 and 10 empty; with Gln missing the reference still starts at 2 and drops 9 and 10. The aligner starts at 4 and fills 10; the reference leaves 10 empty in all 375 comparison-set lambda domains, and we do in all but this one. Forbidding 10 alone gives a third labeling, not the reference's |
| 3GK8_1:H | heavy start (`AVHLQG` numbered 3-6, 6A, 6B) | letters, not enumerated; -10.3 vs 0.5 | **Rule not found.** One chain, two letters in FR1 after a start at 3; nothing like it elsewhere, the consensus prefers the rule |

Fixed this round: 4YDL_1:H in Kabat and Chothia (shared H2/FR3 split above).
Resolved by reading every residue (above): 4XCF_1:H, 1MFE_1:H and
1QFW_1:L agree in all three schemes, and 4XAW_1:H is back in the set. No
remaining miss is a tie: none has a second chain holding its window.
5CEY_2:L is not a twin of 5CEY_1:L (it has an extra leading Ser).

Conclusion for the comparison set: 99.75% (Kabat, 3,602 of 3,611), 99.75%
(Chothia, 3,594 of 3,603), 99.81% (Martin, 3,584 of 3,591); every remaining
miss is listed in `RESIDUAL` with its placement evidence.

**What the stubs show.** Residues in front of our domain were the largest
single cause in the previous round, and they show that the reference does
not number a domain alone: it numbers the chain. Of the 27 domains that
were listed before, 12 were a chain start. 4G6A_1:H and 4G6A_2:H, for
instance, differ in one leading `V` (ours started after it); 4N0Y_1:H, 4G6A_1:H
and 2HH0_1:H begin `VQLLEQSG` and carry 6A where 1R70_1:H, `VKLLEQSG`, does
not: one residue decides between a leading gap with an insertion and a
straight run, with the insertion column at 3 bits below a match. Seven
light chains start with a stub of one to seven residues that the reference
numbers from 5 (6AOD_1:L, 5WB9_1:L, 4JY6_1:L, 4UOM_1:L, 5FYL_1:L, 4R26_1:L,
6MCO_1:L). The reference chains hold only the variable domain plus a few
constant residues (no tags); the library numbers a domain after a tag as
before.

The numbering does not vary. `real_chains_number_identically_every_run_and_thread_count`
numbers every reference and held-out chain twice on one thread and on
pools of 1 and 4 threads, and compares domain ranges, score bits and the
labels of every scheme; `numbering_is_identical_every_run_and_thread_count`
does the same without fixtures. Profiles are built with ordered maps,
equal profile scores go to the later profile of a fixed table, and the one
parallel loop (`chain_domains`) collects in chain order.

On the rebuilt held-out set, whole domains identical under Kabat went from
182 to 179 of 198. The three are heavy chains of 1MRD, 1MRE and 1MRF, which
their depositors number with a three-residue framework deletion at 73-75
and 82A-C, where the reference program counts residues (66-92, no
insertion). The reference is the ground truth for this check.

What the misses are, from a read of the held-out depositor set:

- Lambda FR4: depositors split between `106, 106A, 107` (used here,
  Kabat's published rule) and plain `106, 107, 108`. The scheme authors'
  program uses 106A with a lambda-type J segment, and so do we.
- CDR-H3 insertions lettered from the far end (`100J, 100K` for the last
  two inserted residues); the common convention, `100A...`, is used here.
- Depositors who numbered by Kabat up to some point and sequentially after
  (trastuzumab heavy chain 1N8Z is Kabat through H52, then 1..n).
- A few framework insertions and deletions numbered at a different site.

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

Timing, release build, one thread: 10,000 chains (1.7 M residues, 6,862
domains found) in 4.4 s, 437 µs per chain.

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

## External backend

Native numbering is the default and needs nothing installed. For a second
opinion, `sequence antibody backend anarci` (header menu "Numbers from")
takes the numbers from ANARCI (Dunbar & Deane 2016, Bioinformatics 32:298;
BSD-3-Clause), which the user installs; nothing of it is bundled or linked.
Regions still come from the CDR definitions chosen under "CDRs".

Install (needs HMMER 3): `conda create -n anarci -c conda-forge -c bioconda
python=3.10 hmmer anarci`. On Windows, ANARCI runs inside WSL; install it in
the Linux distribution the same way.

How vizviz finds it, first match wins: `sequence antibody exe PATH` (saved in
sessions), the `VIZVIZ_ANARCI` environment variable, `ANARCI` on the PATH,
then `ANARCI` inside WSL (Windows only, on the distribution's login PATH). A
conda install is usually not on that PATH: point at it with a `wsl:` prefix,
for example `VIZVIZ_ANARCI=wsl:/home/me/miniconda3/envs/anarci/bin/ANARCI`
(paths without spaces; its directory is added to the PATH for the run).

All protein chains of a structure (80 residues or more) go in one FASTA, one
ANARCI run per scheme, on a background thread; switching scheme or CDR
definition afterwards does not rerun it. The track is empty while the header
says "ANARCI running...". If the program is missing or fails, the header says
why and the track shows native numbering.

What differs from native: gamma and delta receptor chains are dropped (no
native equivalent), and receptors are numbered in IMGT only, as natively.
ANARCI does not number shark VNAR. Its AHo range can end one residue after
its IMGT range; the extra residue is left unlabelled. The `select cdr`
expression still uses native domains whatever the backend
(`antibody::external::cdr_residues_external` is the ready bridge).

Checked against ANARCI 2024.05.21 with HMMER 3.4 (conda, under WSL). Hand-made
trastuzumab (heavy + light) agrees with native in all 227 residues in every
scheme; 1N8Z in 227 of 228 (native's domain is one residue longer than
ANARCI's). On the reference set above (seed entries dropped), same domains
found by both, ANARCI against the scheme authors' program, as the native
columns above:

| | Domains | Native identical | ANARCI identical | Residues native | Residues ANARCI | Native = ANARCI residues |
|---|---|---|---|---|---|---|
| Kabat | 3,625 | 3,602 (99.4%) | 3,180 (87.7%) | 99.87% | 99.63% | 99.73% |
| Chothia | 3,617 | 3,593 (99.3%) | 3,171 (87.7%) | 99.87% | 99.59% | 99.68% |
| Martin | 3,605 | 3,582 (99.4%) | 3,054 (84.7%) | 99.87% | 99.40% | 99.50% |

ANARCI also finds 17 domains native does not, and native finds none that
ANARCI misses. Residues are those inside both domain ranges; "identical"
needs the whole ANARCI domain to match the reference. Reprint with `cargo
test --release -p vv-core --test external_numbering -- --ignored --nocapture`
with `VIZVIZ_ANARCI` set (the run takes a while: one ANARCI process per
scheme over about 3,600 chains).

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
