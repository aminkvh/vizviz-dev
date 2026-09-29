# Scientific validation policy

Non-negotiable. Any analysis method this project ships (RMSD/RMSF,
secondary structure, contact maps, SASA, docking scores, trajectory stats,
etc.) must be demonstrably correct before it ships, and stay correct as
code changes. A beautiful renderer showing a wrong number is worse than no
number.

## The rule

**No analysis feature merges without a recorded validation.** Each method
needs, committed alongside it:

- A **reference source** — a published reference value, or an established
  tool's output on the same fixed input (e.g. check our RMSD against
  MDAnalysis's or MDTraj's on the same trajectory).
- A **fixed input fixture** committed to the repo, not "trust me, I ran it
  once locally."
- A **regression test** that fails CI if output drifts past a stated,
  justified tolerance.
- A short **validation note**: method, reference, tolerance, date, who
  validated it.

## Prefer wrapping over reimplementing

Default to wrapping validated existing libraries (MDAnalysis, MDTraj,
Biotite, RDKit, OpenMM) via the plugin/IPC boundary rather than
reimplementing analysis algorithms from scratch. Write our own only when
performance requires it or licensing forces it (see ENVIRONMENT.md).
Reference tools that are GPL or non-commercial can still
be used offline in the dev/CI pipeline to generate validation fixtures —
that doesn't make them a runtime dependency of the shipped product.

## What this isn't

- Not a ban on looking at existing implementations of standard, published
  algorithms (Kabsch RMSD, DSSP-style secondary structure) — those are
  established method, not proprietary IP.
- Not required for rendering/UI features — only for anything producing a
  number or classification a researcher might cite or act on.

## Tracking

| Method | Reference | Check | Where | Date |
|---|---|---|---|---|
| Bond perception (`vv_core::bonds::perceive`) | Standard-residue templates from the wwPDB Chemical Component Dictionary (`_chem_comp_bond`, fetched 2026-09-25); named links (peptide, nucleic, disulfide ≤ 2.3 Å, glycosidic) from standard biochemistry; remaining distance rule uses Cordero et al. 2008 covalent radii + 0.45 Å tolerance, a per-element valence cap, and a 0.6×-ideal clash floor | On 1CRN: every residue has N–CA, CA–C, C–O; all 45 peptide bonds with no extras; exactly 3 disulfides; no atom > 4 bonds; heavy-atom bond count matches an independent spanning-tree-plus-rings count exactly. Same independent count on 1UBQ (exact) and 4HHB (within 2%, the gap being real crystallographic disorder in a few exposed side-chain tips, verified by distance). 4 heme irons each with exactly 5 bonds (4 porphyrin N + proximal His). All 21 `LINK`-recorded glycosidic bonds in 6X3Z present, and recovered by the geometric rule alone with `explicit_bonds` cleared. Unit tests for the standard-residue-vs-clash distinction (two different residues 1.7 Å apart never bond; the same geometry does once neither has a template), metal/single-ion exclusion, valence-cap priority-to-shortest, and verbatim MD topology bonds. | `crates/vv-io/tests/fixtures.rs`, `crates/vv-io/tests/glycan.rs`, `crates/vv-core/src/bonds/mod.rs`, `crates/vv-core/src/bonds/templates.rs` | 2026-09-25 |
| mmCIF / PDB readers | Same entries in both formats from RCSB | Atom-for-atom agreement (positions, elements, names, serials, secondary structure) on 1CRN, 1UBQ, 4HHB, 1AKE; writer round-trip of a 20k-atom synthetic structure | `crates/vv-io/tests/fixtures.rs` | 2026-09-18 |
| PDB / mmCIF writers (`pdb_write`, `mmcif_write`) | Round-trip through this project's own readers, wwPDB v3.3 and PDBx/mmCIF column/category layout | Full topology (names, residues, chains, elements, coordinates within 0.001 Å, B-factors, occupancies, `CONECT`/`_struct_conn` bonds) preserved on 1CRN, 4HHB, 1UBQ; a selection writes only those atoms; multi-model round trip; a >99,999-atom synthetic structure round-trips through both formats, PDB via hybrid-36 serials; every fixed-column PDB record line is exactly 80 columns | `crates/vv-io/tests/save.rs` | 2026-09-25 |

Still to do: compare bond tables against an independent tool (gemmi or
Biotite, offline in CI) on a larger fixture set, cover bond orders, and
extend the metal-cofactor template list past HEM/HEC (e.g. iron-sulfur
clusters, chlorophyll).

### Bond perception timings

Release mode, min of 5 runs
(`cargo test -p vv-io --test fixtures --release -- --ignored --nocapture
bond_perception_timing`), on an 8-core laptop CPU. Bond counts are exact.
Sub-millisecond runs carry thread-pool wake-up jitter.

| Structure | Atoms | Bonds | Time |
|---|---|---|---|
| 4HHB | 4,779 | 4,647 | about 1 ms |
| 6X3Z (glycoprotein) | 17,365 | 17,829 | about 2 ms |
| Synthetic protein-like | 1,000,000 | 1,442,618 | about 180 ms |

A distance-only rule is about 1.7x faster on the synthetic system but
makes 2,191,743 bonds there: roughly 750,000 are false contacts from
unbounded degree and no clash floor. The valence cap's
"priority to shortest" fill is the extra cost, and it is near zero on real
structures (6X3Z: 14 of 17,365 atoms need it). The result is pinned by
`crates/vv-io/tests/bonds_pin.rs` (exact count, first and last pair, and a
hash of the sorted bond list) on 1CRN, 1UBQ, 4HHB, 6X3Z and a 20k-atom
synthetic system.
