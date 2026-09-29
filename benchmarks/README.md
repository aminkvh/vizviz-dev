# Benchmarks

`vv-bench` renders a synthetic protein-density structure headlessly through
scripted orbits and reports frame-time percentiles. Results are committed
here per machine so regressions are visible in review; they are never gated
in hosted CI (no GPUs, too noisy).

```bash
cargo run --release -p vv-bench -- --atoms 10000000 --frames 60 \
    --json benchmarks/results/<date>_<gpu>_<atoms>_<setting>.json
cargo run --release -p vv-bench -- --help      # presets, --adapter, --adaptive, --png
```

Presets (camera distance in bounding radii): `overview` 2.74 (structure
framed), `mid` 1.6, `close` 1.15 (just outside the surface), `inside` 0.5
(camera inside the structure — the overdraw worst case).

The v1 bar: 10M atoms, 1920x1080, p95 GPU time ≤ 33 ms on a laptop GPU,
with every atom drawn as a ray-cast sphere (no point fallback).

## Results (p95 GPU ms per preset)

| Date | GPU | Atoms | Setting | overview | mid | close | inside |
|---|---|---|---|---|---|---|---|
| 2026-09-18 | RTX 3060 Laptop (Vulkan) | 10M | no occlusion culling | 26.8 | 49.6 | 55.8 | 44.4 |
| 2026-09-18 | RTX 3060 Laptop (Vulkan) | 10M | occlusion culling | **9.9** | **14.1** | **12.1** | 35.2 (p50 3.0) |
| 2026-09-18 | AMD Radeon iGPU (Vulkan) | 10M | occlusion culling | 28.9 | 37.4 | 34.3 | 67.5 (p50 13.1) |
| 2026-09-18 | AMD Radeon iGPU (Vulkan) | 10M | occlusion + adaptive LOD | 26.6 | 34.5 | 31.2 | 89.0 (p50 11.3) |
| 2026-09-18 | AMD Radeon iGPU (Vulkan) | 2M | occlusion culling | 10.6 | 19.3 | 19.5 | 51.6 (p50 9.3) |
| 2026-09-18 | RTX 3060 Laptop (Vulkan) | 3J3Q 2.44M (real) | spacefill | 2.8 | 2.9 | 3.2 | 2.8 |
| 2026-09-18 | RTX 3060 Laptop (Vulkan) | 3J3Q 2.44M (real) | ball-and-stick, 2.5M bonds | 1.0 | 1.4 | 1.7 | 1.2 |
| 2026-09-18 | RTX 3060 Laptop (Vulkan) | 10M | ball-and-stick, 21.8M bonds | 4.3 | 8.6 | 8.3 | 9.0 |

Real structures: `cargo xtask fetch` downloads 4V6X and 3J3Q (HIV-1 capsid,
2,440,800 atoms) into `fixtures/large/`; `--file` benchmarks them. Loading
3J3Q from gzip takes 0.8 s, bond perception 0.24 s, GPU upload 43 ms.

Files: `results/<date>_<gpu>_<atoms>_<setting>.json`.

**Molecule classification (2026-09-29).** One class per residue, computed
at structure construction (`cargo test --release -p vv-io --test
residue_class_timing -- --ignored --nocapture`): a synthetic membrane of
100k lipids and 2M TIP3 waters (9.0M atoms, 2.1M residues) classifies in
4 ms when every name is known, and in 33 ms with 5,000 unknown-name
residues (their bonds perceived in small per-residue windows).

**8GLV, cartoon and Gaussian added (2026-09-24).** Same setup, 120 frames per orbit (3°
per frame), `--look`, perspective. p95 ms, overview / mid / close / inside:

| GPU | spacefill | ball-and-stick | cartoon | Gaussian surface |
|---|---|---|---|---|
| RTX 3060 Laptop | 5.1 / 6.8 / 7.2 / 5.5 | 2.9 / 3.1 / 4.6 / 4.2 | 5.3 / 3.8 / 4.4 / 3.5 | 1.8 / 2.5 / 3.0 / 3.2 |
| AMD iGPU | 15.7 / 25.0 / 26.9 / 21.1 | 17.5 / 25.0 / **35.9** / 32.9 | 10.7 / 18.3 / 21.9 / 18.8 | 8.5 / 14.6 / 17.5 / 18.6 |

What did it: conservative depth for sphere and cylinder impostors, 64-atom
cluster culling, bond lines only past 3 px, half-resolution AO, and a
cartoon built from stored cross-sections with per-residue LOD (it now fits
the iGPU). Files: `results/2026-09-24_*_m6.json`.

**8GLV (2026-09-23), the canonical real structure from now on**
(3,968,189 atoms, cryo-EM; `--file 8GLV.cif`). `--look` turns on what the
app draws by default (AO, depth cue), and from this date the GPU number
spans the whole frame: both occlusion phases and the post passes. p95 ms:

| GPU | Size | spacefill | ball-and-stick | cartoon | Gaussian surface |
|---|---|---|---|---|---|
| RTX 3060 Laptop | 1080p | 5.4 / 9.1 / 12.6 / 14.0 | 2.9 / 3.0 / 4.2 / 4.6 | 4.5 / 4.5 / 4.6 / 4.4 | 1.4 / 2.6 / 3.1 / 3.3 |
| RTX 3060 Laptop | 4K | 10.5 / 23.9 / **38.5 / 49.8** | 4.6 / 8.6 / 15.7 / 16.3 | 12.5 / 6.5 / 8.4 / 9.2 | 3.9 / 8.1 / 10.7 / 13.1 |
| AMD iGPU | 1080p | 17.9 / 30.6 / **41.8 / 45.2** | 19.8 / 26.0 / 33.0 / 32.1 | out of memory | 8.8 / 15.8 / 19.0 / 19.7 |

(overview / mid / close / inside; bold misses 30 fps.) Files:
`results/2026-09-23_*_8GLV_*_look.json`.

The old `inside` p95 spikes (previous-frame depth pyramid failing to cull
under fast rotation) are gone with two-phase occlusion culling, which
also stopped it from culling *visible* atoms while the camera moved.

**Transparency (2026-09-25).** Weighted-blended OIT extended to every
representation (see docs/RENDERING.md's "Transparency"), 1M-atom
synthetic spacefill, `--material opaque` vs `--material transparent`
(opacity 0.30), perspective, no `--look`. p95 GPU ms,
overview / mid / close / inside:

| GPU | opaque | transparent |
|---|---|---|
| RTX 3060 Laptop (Vulkan) | 2.8 / 3.0 / 2.5 / 3.8 | 16.0 / 26.8 / 30.6 / 25.4 |
| AMD iGPU (Vulkan) | 13.8 / 13.7 / 11.7 / 17.3 | 50.8 / 103.9 / 130.7 / 115.2 |

Transparent is 5-11x slower at this size: glass writes no depth, so the
occlusion pass has nothing of its own to cull against and every atom's
own layer survives and gets fragment-shaded every frame (`visible_quads`
hits the full 1M at `overview`, vs ~460k opaque) -- inherent to any OIT
scheme, not specific to this implementation, and why an opaque-only scene
must never pay for it (`opaque_only_render_is_unaffected_by_an_unseen_
glass_atom` in `crates/vv-render/tests/headless.rs`). 4HHB (4,779 atoms)
spacefill transparent + an opaque cartoon, 1400x900, wall-clock (`cargo
test -p vv-render --release --test headless -- --ignored
frame_time_4hhb_transparent_spacefill_and_opaque_cartoon --nocapture`):
RTX 3060 Laptop 1.53 -> 1.93 ms/frame, AMD iGPU 4.94 -> 7.02 ms/frame
(opaque -> transparent) -- a realistic small-molecule scene barely
notices it. Files: `results/2026-09-25_*_1M_*_spacefill.json`. (Rechecked
after glass self-occlusion culling below: this test uses
`RenderSettings::default()`, i.e. `fast_glass: false`, so it never
engages the cull and these numbers stay within run-to-run noise of the
original -- expected, and correct: `fast_glass` targets the large, dense
case below, which a still small-molecule render never was.)

**Glass self-occlusion culling (2026-09-25, after).** Same setup and
`--material transparent`, after culling a glass atom against the
previous frame's nearest-glass depth once it estimates that atom
contributes less than 1/255 once composited (docs/RENDERING.md's
"Transparency", "Glass self-occlusion culling"; commit `e8c5727` is the
baseline this built on). **First shipped unconditionally; after merge,
a before/after comparison of `close.png` showed a real, visible look
change on a densely overlapping scene (spheres leave gaps a ray passes
through that the isotropic-density estimate doesn't model, so some
genuine contribution gets dropped, reading darker and more saturated)
-- see "Interaction-time only" in docs/RENDERING.md.** Now gated behind
`RenderSettings::fast_glass` (`--fast-glass` here), off by default: the
`transparent` column below is what a still view, a screenshot, an
export, or `--material transparent` *without* `--fast-glass` renders --
identical to *before this feature existed at all* (`visible_quads`
matches exactly). `--fast-glass` is what `vv-app` turns on while the
camera is actively moving or a trajectory is playing, never for a still
frame. p95 GPU ms, overview / mid / close / inside:

| GPU | opaque (unchanged) | transparent (default, exact) | transparent `--fast-glass` | speedup while interacting |
|---|---|---|---|---|
| RTX 3060 Laptop (Vulkan) | 2.9 / 2.9 / 2.5 / 3.8 | 16.0 / 26.8 / 30.6 / 25.4 (unchanged) | 4.6 / 7.4 / 7.0 / 8.0 | 3.5x / 3.6x / 4.4x / 3.2x |
| AMD iGPU (Vulkan) | 14.1 / 14.0 / 12.4 / 16.9 | 50.8 / 103.9 / 130.7 / 115.2 (unchanged) | 15.7 / 29.7 / 34.4 / 38.6 | 3.3x / 3.6x / 3.8x / 3.0x |

The "transparent (default, exact)" column is the original `--material
transparent` numbers a few sections up, unchanged: `fast_glass` off
computes a `0` margin for every item regardless of material, so this
path is identical by construction, not a separate measurement worth its
own noise -- confirmed by `visible_quads` matching exactly at every
preset (1,000,000/847,221/604,009/235,823 both GPUs, checked against a
fresh `--material transparent` run without `--fast-glass`). `visible_
quads` at `overview`, RTX 3060 Laptop: 1,000,000 (default) -> 134,217
(`--fast-glass`) -- the atom count dropped, not just its cost per pixel,
which is what the sub-linear-with-resolution scaling measured before
picking this fix over half-resolution accumulation predicted (docs/
RENDERING.md). The opaque columns are the same `--material opaque` run
rechecked after the change; `visible_quads` for every opaque preset
matches the `before` run exactly (459180/293380/138913/891 RTX,
461301/292011/138925/891 AMD), and `opaque_only_render_is_unaffected_
by_an_unseen_glass_atom` (pixel-exact), `glass_self_occlusion_culling_
is_exact_when_fast_glass_is_off`, and the full `cargo test -p vv-render
--release --test headless` suite pass on both adapters (`VIZVIZ_TEST_
ADAPTER=AMD` for the iGPU); `--fast-glass` numbers are within
measurement noise of the first version of this fix (RTX `overview` 4.59
-> 4.6 here), confirming the phase-2 `occluded_bits` change above is
pure waste removal, not a behaviour change. Files: `results/2026-09-25_
*_1M_opaque_spacefill_after.json`, `results/2026-09-25_*_1M_transparent_
spacefill.json` (the unchanged default path), `results/2026-09-25_*_1M_
transparent_spacefill_fastglass.json`. Images: `images/{still,
interacting}/{overview,close}.png` (not committed:
regenerate with
`vv-bench --atoms 1000000 --material transparent --preset
overview,close --png <dir>`, adding `--fast-glass` for the interacting
set).
