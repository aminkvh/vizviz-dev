# Rendering pipeline

How `vv-render` draws millions of atoms per frame. Read this before touching
the shaders in `crates/vv-render/src/shaders/`.

## Data on the GPU

Atoms live in pages of up to 2^22 atoms (64 MiB), each a storage buffer of
`vec4<f32>` = position + display radius, plus an RGBA8 color buffer. Radii
are derived on upload (`vv_core::Element::vdw_radius`) and never change.
Colors are derived under a `ColorScheme` (element, chain, or b-factor;
`crates/vv-render/src/color.rs`) and rewritten in place — `Page.colors`
carries `COPY_DST` for exactly this — by `GpuStructure::recolor` whenever
the scheme changes, with no re-upload of positions or radii. Per page there
are also two index buffers and one indirect-arguments buffer the cull pass
fills in; the CPU never touches per-atom data after upload otherwise.

There are two full-frame color textures, not one: `color` (the opaque
pass's raw output; `Renderer::color_view`/`read_color`, what screenshots
read) and `fxaa` (the post-process output; `Renderer::display_view`, what's
actually shown). See "Style presets" and "Anti-aliasing" below for why
they're kept separate.

## One frame

1. **Cull (compute, 256 threads/workgroup, one thread per atom).**
   Frustum test against six planes → occlusion test against the depth
   pyramid of the previous frame → size classification: projected radius
   above a per-atom jittered threshold → *quad*, else *point*. Survivors
   take a slot from a workgroup-local atomic counter; one thread claims the
   workgroup's range with a single global atomic add. Output: compacted
   index lists and `DrawIndirectArgs`.
2. **Opaque (render pass, two indirect draws).** Quads are ray-cast sphere
   impostors: a screen-aligned billboard per atom, the fragment shader
   intersects the eye ray with the sphere, writes the true depth via
   `frag_depth`, and shades. Points are one vertex per atom with depth from
   the vertex. Targets: sRGB color, view-space normal (for future SSAO),
   reversed-Z `Depth32Float` (compare `Greater`, clear 0).
3. **Post-process (render pass, one full-screen triangle).** Reads `color`,
   writes `fxaa` — FXAA if `RenderSettings.fxaa` is set (the default),
   otherwise a plain copy, so `display_view` is always the same texture
   regardless of the setting (see "Anti-aliasing" below). Always runs; not
   bracketed by a GPU timestamp (see "Why it is fast" below).
4. **Depth pyramid (compute).** Level 0 copies the depth buffer; each
   further level is the min (= farthest, in reversed-Z) of a 2x2 block,
   folding odd leftover rows/columns into the last texel so the pyramid
   never claims a region is nearer than it is. Used by step 1 next frame.

GPU timestamps bracket steps 1, 2, and 4 (`Renderer::last_times()` reports
them as `cull_ms`/`draw_ms`/`hiz_ms`); step 3 is not timestamped.

When anything drawn is transparent, two more passes run, both skipped
entirely otherwise (see "Transparency" below for why and how): a "glass
id" pass before step 2 gives every transparent atom, bond, cartoon and
glycan vertex its own pick id wherever nothing opaque ends up covering
it; a "glass" pass after step 3 draws every transparent fragment's
weighted-blended contribution and composites it onto `outlined`. Neither
is timestamped (see "Transparency"'s own note on this).

There is no swapchain-blit step in the live app: `vv-app` registers
`display_view()` directly as an egui texture and lets egui's own pass
composite it into the window, so `display_view()`/`read_color`/
`read_display_color` are how anything (screenshots, tests, the viewport)
ever sees a frame. `Renderer::blit()` exists for presenting straight to a
swapchain without egui in the loop, but nothing in this repo calls it yet.

## Why it is fast

- Impostors cost fragments, not triangles: 10M small spheres are
  fragment-bound (~25 ms on an RTX 3060 laptop at 1080p), and instanced
  vs. non-instanced quads made no difference.
- `frag_depth` disables early-Z, so hidden atoms would run the full
  fragment shader before being rejected. **Occlusion culling in the
  compute pass is what makes 10M atoms fit in 10–14 ms**: with a dense
  structure only the front layer survives (see benchmarks/README.md).
- The point tier and the adaptive threshold (`AdaptiveLod`) are the safety
  valve for weaker GPUs; the threshold is jittered per atom so a structure
  whose atoms all project to the same size degrades smoothly.
- The post-process pass (FXAA or its passthrough) is one full-screen
  triangle sampling a handful of texels per pixel — cheap next to a cull
  or draw pass over millions of atoms, which is why it isn't worth a
  dedicated GPU timestamp slot (`AdaptiveLod` steers off `cull_ms +
  draw_ms + hiz_ms`, not total frame time, so it doesn't see this pass
  either; it hasn't needed to).

## Correctness rules (tested in `crates/vv-render/tests/headless.rs`)

- A framed structure must survive frustum culling entirely.
- With a static camera, occlusion culling must not change a single pixel.
- An opaque-only scene must render identically whether or not the glass
  pipelines exist and whether or not an (unseen) glass item is also in
  the scene: adding transparency support must cost an opaque scene
  nothing (see "Transparency" below, and
  `opaque_only_render_is_unaffected_by_an_unseen_glass_atom`).
- Weighted-blended OIT must not depend on submission order or on which
  of two overlapping glass atoms happens to be nearer the camera
  (`overlapping_transparent_atoms_blend_the_same_regardless_of_draw_
  order`/`..._after_a_180_degree_rotation`).
- Run the GPU tests on a specific adapter with `VIZVIZ_TEST_ADAPTER=<name>`.
  Do this on at least one NVIDIA and one AMD device before merging shader
  changes: a shared-memory prefix scan that was correct on NVIDIA silently
  dropped wave tails on AMD, which is why compaction uses atomics.

## Style presets

`StylePreset` (`crates/vv-render/src/style.rs`) is four named points in a
small material space — ambient/diffuse/specular/shininess plus an optional
toon-band count — read every frame from a dedicated uniform at
`@group(0) @binding(2)` in `shaders/draw.wgsl` (declared only there, not in
the shared `atoms.wgsl`, since the cull and Hi-Z shaders never shade
anything). `shade()` (spheres, cylinders) and `shade_point()` (the point
tier) both read it, so switching presets never splits the two LOD tiers.
`toon_bands > 0` quantizes the diffuse term (`floor(diffuse * bands) /
bands`) for the flat/cel look; every other preset leaves it at `0` for
continuous shading. A preset also suggests a background
(`StylePreset::background()`), but `RenderSettings.background` stays an
independent field — vv-app seeds it from the preset when the user switches,
but a screenshot export will want to override it without touching the
shading model. This is a session view preference (`vv-app`'s
`ViewSettings`), not document state: no `vv_scene` mirror, no `Command`
variant, unlike `ColorScheme`.

**Tonemapping.** `shade()` ends with `tonemap()`, an ACES filmic fit
(`x(2.51x+0.03) / (x(2.43x+0.59)+0.14)`) selected per preset through
`StyleUniform.tonemap`. It runs *inside* the fragment shader, before the
8-bit sRGB target quantizes, which is the only place values above 1.0
still exist: the Glossy preset's specular term reaches ~1.6 and would
otherwise clip to a flat white disc. Dark Presentation and Glossy use it;
Publication White and Flat/Cel stay linear so a picked color prints as
picked. The CPU backend calls the same Rust function
(`style::tonemap_channel`), so both backends agree.

**Outlines.** A post pass (`shaders/outline.wgsl`) between the opaque
pass and FXAA reads `color`, the normal buffer and the depth buffer and
darkens a pixel when any 4-neighbour is a depth jump (reversed-Z
linearized to view depth; background counts as infinitely far, so the
silhouette lands on the object's own rim, never on the background) or a
normal crease (`dot < 0.5`, which catches where two overlapping atoms
meet even though depth is continuous there). Thresholds are in
`renderer.rs` (`OUTLINE_*`); the depth one grows with view depth so a
distant structure does not dissolve into noise. The pass always runs,
uniform-gated (a copy when off), so the texture chain
`color -> outlined -> display` never changes identity; `read_color()`
returns the outlined image, so screenshots and SSAA export include the
lines, and export doubles the line width to survive the 2x downsample.
Flat/Cel turns outlines on by default (`StylePreset::outline`); the app
exposes the toggle as `outline on|off`.

**Selection outline.** The active selection is a green halo
(`selection.rs`, `shaders/selection_outline.wgsl`). The pass
reads this frame's pick-id target and two bitsets over the frame-wide
atom and bond ids (`Renderer::set_selection`, uploaded only when the
selection or the draw list changes), and colours each unselected pixel
with a selected pixel within `RenderSettings::selection_width` pixels
(the app uses 2 logical px). So it follows any representation that
writes ids, is hidden wherever the selected atoms are, and moves with
them in the same frame. It is drawn after FXAA onto the display only:
screenshots and exports never show it. A bond counts as selected when
both its atoms are. Not outlined: the Gaussian surface (writes no ids)
and SES patches drawn as points when tiny.

## Anti-aliasing

No MSAA: `frag_depth` impostors don't resolve well under it (the
rasterizer's coverage is for the billboard quad, not the ray-cast sphere
inside it). Two different answers instead, for two different jobs:

- **FXAA**, for the live viewport. A single full-screen post-process pass
  (`fs_fxaa` in `shaders/blit.wgsl`, following the structure of Timothy
  Lottes' public NVIDIA FXAA whitepaper — a published technique, not code
  taken from any viewer) reads `color` and writes `fxaa`
  (`Renderer::display_view`). Luma edge-detection finds a local contrast
  step, estimates its direction from the four diagonal neighbors, and
  blends along it; flat regions (`range < max(0.0312, lMax * 0.125)`) are
  left untouched to avoid softening real detail. `RenderSettings.fxaa`
  (default on) picks between this and a plain copy — never whether the
  pass runs at all, since `display_view` is registered with egui once and
  switching which *texture* it points at would need the same
  re-registration a resize does (`ensure_viewport_texture`/
  `maybe_resize_viewport`).
- **Live supersampling** (`supersample 1.5|2`, Scene Settings): the app
  asks the renderer for a target that many times larger than the
  viewport tab and lets egui's linear texture filter average it down.
  Real subpixel coverage for thin tubes and sticks, at 2.25x or 4x the
  pixel cost; off by default so the 10M-atom fps bar still holds.
- **2x SSAA**, for screenshot export (see "Style presets" for why exports
  and the live viewport are allowed to diverge here). Real subpixel
  information from rendering at double resolution beats FXAA's single-
  sample edge guess, so exports render at 2x and downsample in linear
  light (`crates/vv-render/src/postprocess.rs`) rather than running FXAA —
  `read_color()` reads the pre-FXAA (but post-outline) texture, and
  `vv-bench` explicitly sets `fxaa: false` and `outline: false` so
  committed benchmark numbers don't silently include untimestamped passes
  (see "Why it is fast" below).

## Bonds and representations

Ball-and-stick scales atom radii by 0.25 and draws bonds as ray-cast
capped cylinders (radius 0.15 Å): a second compute entry point culls
bonds by their bounding sphere (frustum, occlusion, and "thinner than
half a pixel" — the atoms cover it) into a third indirect draw; the
vertex shader builds an oriented billboard along the bond's on-screen
direction, and the fragment shader intersects the eye ray with the finite
cylinder and its end caps, colouring each half by its atom. Bonds are
stored per page as page-relative index pairs; the rare bond crossing a
page boundary is dropped. `Representation` in `RenderSettings` switches
between spacefill and ball-and-stick without re-uploading anything.

**Tube (cartoon-lite).** No new shader: `vv_core::backbone` picks one
trace atom per polymer residue (`CA`, or `P`/`C4'` for nucleotides) and
splits the chain wherever consecutive trace atoms are too far apart to be
bonded (4.5 Å for protein, 8.5 Å for nucleic acid). The trace then feeds
`vv_core::cartoon::tube_plan`, the same ribbon-cartoon mesh path the
Cartoon representation uses (`vv_core::cartoon::plan`/`CartoonPlan`,
`vv_render::CartoonGpu`, `shaders/cartoon.wgsl`/`cartoon_frame.wgsl`/
`cartoon_lod.wgsl`): every residue is a round cross-section (roundness 1,
half-width = half-thickness), so it draws as a real swept tube instead of
spheres on a stick, with the radius interpolated smoothly residue to
residue -- same GPU code, same picking, clipping, LOD and path-traced
mesh as a ribbon, just a circle instead of a superellipse. The radius is
`TUBE_RADIUS` (0.3 Å) or the `radius` rep option; with `radius_by` set to
putty, each residue's radius instead comes from
`vv_core::backbone::putty_radius`: linear from `radius_min` to `radius`
over the B-factor range of the tube's own atoms
(`gpu_cache::tube_radius_fn`), baked into the plan once per rebuild since
B-factor doesn't change across frames -- only the spline moves
(`CartoonGpu::set_frame`), same as a playing cartoon. The mesh itself has
no end caps (`vv_core::cartoon`'s scope, shared with the ribbon); a tube
closes each chain end, break or selection cut with a plain sphere at that
end's own radius (`CartoonMesh::loose_ends`), cheap because a Catmull-Rom
spline passes exactly through its control points, so the end section's
center already equals the end atom's own position. The non-polymer
content goes up beside it as described under "Cartoon and tube scene
composition" below. Picking maps
back to real atoms the same way a cartoon does (each mesh section
remembers its nearest trace atom), so picking anywhere on the tube
selects that residue's CA and picking a heme stick selects the heme
atoms. The app builds this lazily the first time a structure is shown as
a tube and keeps it (`crates/vv-app/src/gpu_cache.rs`); a moving frame
only re-evaluates the spline, never rebuilds the mesh, so the section
count under the GPU buffers never changes. Helices and strands are not
drawn as ribbons or arrows for Tube: that is the Cartoon representation,
which needs validated secondary-structure assignment (docs/VALIDATION.md).

## Cartoon and tube scene composition

A cartoon or tube also draws what a reader expects beside the polymer,
decided per residue from its class and roles (`vv_core::companion`,
`vv_core::residue_class`) and limited to the rep's selection. Defaults
(each a `repopt` on/off choice):

| Content | Default | Drawn as |
| --- | --- | --- |
| ligands, cofactors (`ligands`) | on | ball-and-stick |
| ions (`ions`) | on | spheres at half van der Waals radius |
| glycans (`glycans`) | on | thin sticks |
| lipids (`lipids`) | on | bound: thin sticks; membrane: lines |
| water (`water`) | off | ball-and-stick |
| crystallization additives (`additives`) | off | ball-and-stick |

Each kind is its own piece of geometry (`gpu_cache/companions.rs`), so
picking, selection outlines, transparency, the occlusion volume and the
path tracer treat it like any other atoms.

A cartoon draws each nucleotide's base (`bases`, `gpu_cache/bases.rs`,
`vv_core::bases`). `plate` (default) is a slab on the base's real ring
atoms (nine for a purine, six for a pyrimidine, fan-triangulated from
the ring centroid, 0.4 A thick), with sticks from the backbone trace atom
to C1' and on to the glycosidic atom; it goes through the glycan mesh
path, so a pick on a plate names a ring atom. `ladder` draws the same
stems and one rung per base pair, C1' to C1' in the two bases' colours.
Pairs are found geometrically: a purine and a pyrimidine whose N1 and N3
are within 3.4 A (or with the wobble contacts N1-O2 and O6-N3), rings
within 40 degrees of parallel and 2 A of co-planar; each base pairs at
most once. `stick` is the older single stick to the pairing atom. Slabs
take their colour from the rep's coloring at the glycosidic atom, so
`color nucleotide` gives the usual per-base colours. Under `element`
coloring that atom is always nitrogen, so each plate is instead coloured
by the base's identity (the nucleotide scheme). Works under cartoon and tube.

## Interaction overlay

`interactions hbond|metal|saltbridge on|off` (`vv_core::interactions`,
`gpu_cache/interactions.rs`) draws each contact as a run of short
cylinders through the ball-and-stick pipeline, so the dashes are
depth-tested and reach path-traced renders and SVG export (one `<line>` per dash, depth-sorted with the atoms). A metal contact is a metal to
an N, O or S of another residue within mean + 3 SD of that metal-ligand
pair's distance in high-resolution PDB entries (Zheng et al. 2008, J Inorg
Biochem 102:1765-1776, Table 3, PDB-HR; e.g. Zn-N 2.40, Zn-S 2.50, Zn-O
2.68 A). The paper gives distributions, so the 3 SD tolerance is this
project's choice. Pairings the table lacks (Mn-S, Ca-N, ...) take that
metal's widest tabulated cutoff; a metal the table lacks uses the paper's
3 A contact radius.

## Transparency

A material with opacity below 1 (`Material::is_glass`: opacity < 1 or
`transmode` set) draws through a second path instead of the opaque one
above, on every representation: spheres/impostors, bonds and sticks,
cartoon, tube and glycan meshes, and the Gaussian/skin/SES surfaces.

**Technique: weighted-blended OIT (McGuire & Bavoil, JCGT 2013),
already the SES/skin/Gaussian glass path's own technique, now shared by
everything.** Every glass fragment adds its own premultiplied colour and
alpha, weighted toward the nearer layers, into one `Rgba16Float` `accum`
target (additive blend) and multiplies one `R8Unorm` `reveal` target by
`1 - alpha` (see `glass_lit`/`glass_point`/`glass_accum` in
`shaders/shading.wgsl`); `shaders/glass_resolve.wgsl` composites
`accum.rgb / accum.a` over the opaque image with coverage `1 - reveal`.
Considered against the alternatives named in the task this shipped
under:

- **Depth peeling** (draw the whole scene once per transparent layer,
  peeling off the nearest not-yet-drawn surface each pass) needs one
  full extra impostor pass *per layer*: a spacefill structure with
  dozens of overlapping transparent atoms along one ray would cost tens
  of extra opaque-sized passes. Weighted-blended OIT costs exactly one
  extra pass regardless of how many layers overlap.
- **The pre-existing nearest-layer-only path** (an owner buffer picks
  the one nearest glass fragment and discards the rest) is *why*
  overlapping transparent atoms used to pop when rotating: whichever
  atom happens to be nearest changes discontinuously as the camera
  moves, and every other overlapping layer's colour and alpha are simply
  thrown away. Weighted-blended OIT sums every layer instead, so the
  result is a continuous function of each layer's own depth and colour,
  order-independent by construction (`overlapping_transparent_atoms_
  blend_the_same_regardless_of_draw_order`/`..._after_a_180_degree_
  rotation` in `crates/vv-render/tests/headless.rs`).

It is an approximation (the resolved colour is a depth-weighted average,
not an exact per-pixel sort), which is exactly the tradeoff McGuire &
Bavoil describe. Surfaces already used it, so extending it to every
representation avoids a second scheme.

**Frame changes (`Renderer::render_all`), all skipped when nothing
drawn is glass so an opaque-only scene pays nothing new:**

1. **"Glass id" pass**, before the opaque pass, only for representations
   that carry a frame-wide pick id (spheres, bonds, cartoons, glycans --
   not the Gaussian/skin/SES surfaces, see "Picking" below): draws every
   glass item's id into `id_view` (the very first write to it that
   frame) depth-tested against a scratch buffer cleared to the far
   plane, so overlapping glass atoms still resolve to the nearest one.
2. **Opaque pass**: `id_view`'s load op becomes `Load` instead of
   `Clear` when the id pass ran, so an opaque hit always overwrites a
   glass one wherever it draws — regardless of relative depth, which is
   the picking policy below, not a depth test.
3. Every glass item's own `_glass` fragment entry point (`fs_sphere_
   glass`, `fs_cartoon_glass`, ...) draws into `accum`/`reveal` after the
   outline pass, depth-tested against the opaque depth (`Greater`,
   write disabled) so opaque geometry always occludes glass and glass
   never occludes anything, then `glass_resolve.wgsl` composites onto
   `outlined`.
4. Skin/SES surfaces keep their own owner-buffer prepass first (a patch
   surface tiles many overlapping billboards, so it needs a "which
   patch of *this* surface is nearest" step the other representations
   don't); Gaussian surfaces need no prepass at all (one full-screen
   ray-march is already exactly one layer). This is unchanged and still
   restricted to one glass surface total: if two different transparent
   Skin/SES surfaces overlap, only the nearer one's patches contribute
   — a real, documented gap the task that shipped this left as-is
   (rare: it needs two separate glass surfaces on screen at once, not
   two overlapping transparent atoms, which is the common case and is
   unrestricted).

Neither new pass carries a GPU timestamp: the timer's 8 query slots
(`timing.rs`) are already all spoken for by the existing passes, and
adding a ninth would mean re-plumbing `PassTimes`' public shape for a
path that costs an opaque scene nothing to begin with. Time it by
wall-clock instead (`cargo test -p vv-render --release --test headless
-- --ignored frame_time_4hhb_transparent_spacefill_and_opaque_cartoon
--nocapture`) or with `vv-bench --material transparent`
(benchmarks/README.md's "Transparency" entry).

**Glass self-occlusion culling.** A fully transparent scene is 5-11x
slower than the same scene opaque (RTX 3060 laptop 2.8/3.0/2.5/3.8 ->
15.97/26.75/30.62/25.38 ms p95, AMD iGPU 13.8/13.7/11.7/17.3 ->
50.8/103.9/130.7/115.2 ms p95, 1M-atom synthetic spacefill,
`overview`/`mid`/`close`/`inside`; `benchmarks/results/2026-09-25_*_1M_*`).
Cause: cull's occlusion test (`shaders/cull.wgsl`'s `occluded`) checks
survivors against the *opaque* Hi-Z pyramid, which glass never
contributes to (it never writes depth, on purpose -- see "Frame changes"
above) — so it already culls a glass atom sitting behind opaque geometry
for free, but an all-glass region has nothing in that pyramid to test
against, and every atom in it survives cull and gets fragment-shaded
every frame. Measured before picking a fix: at 1920x1080 the untimed
"glass id" pass is ~20% of the gap (a full impostor pass over every
glass atom, `frag_depth` disabling early-Z, same as the opaque pass
without occlusion culling to shrink it); the rest scales *sub-linearly*
with render target resolution (960x540 -> 1920x1080 -> 3840x2160 gave
~1.8x -> ~2.5x, not the ~4x a purely fragment-bound cost would show),
meaning the dominant cost is bound to the *number of atoms drawn*, not
pixels shaded — ruling out a purely resolution-based fix (half-resolution
accumulation with bilateral upsample) as the main lever, since roughly
half that cost does not shrink with the render target.

So the fix cuts the number of glass atoms drawn, the same lever
occlusion culling already uses for opaque: **cull a glass atom against
*last frame's* nearest-glass depth**, the "glass id" pass's own
`glass_depth` target, copied out (`hiz_copy_pipeline`, reused verbatim,
one texel, no mip chain) right after that pass and before it is
overwritten as Skin/SES scratch (`Renderer::render_all`, the "glass hiz"
compute pass). `shaders/cull.wgsl`'s `glass_hidden` culls an atom whose
own depth lies more than a *margin* behind that reference, sampling only
the centre texel of its own screen footprint (cheaper than a
conservative multi-texel query like `occluded`'s, and exact for a single
atom's footprint, but `cull_clusters` also calls it on a whole cluster's
bounding sphere, whose footprint spans more than one pixel — a cluster
straddling a glass silhouette reads only its centre's reference, part of
what the quality-bound test below measures). It is gated by the same
`RenderSettings::occlusion_culling` flag as the opaque test
(`culled_by_depth`) and fed into the *same* two-phase cluster/atom cull
pass opaque items already go through — no new pass, no new indirect
buffers, no separate disocclusion handling to reason about.

The margin (`vv_render::scene::glass_cull_margin`, one `PageParams`
field per item) approximates "beyond this many layers, weighted-blended
OIT's contribution rounds to zero": `K(opacity) = ceil(ln(1/255) / ln(1
- opacity))` layers of transmittance `(1 - opacity)` cross the 1/255
threshold, times a mean atom spacing estimated from packing the item's
atom count into its bounding sphere (`(volume / atom_count) ^ (1/3)`, an
isotropic-density approximation), plus one extra layer of slack. **Both
halves are approximations, not exact bounds**: real structures are not
isotropic
(a flat sheet of atoms and a packed globule with the same count and
bounding sphere have very different true spacing), and weighted-blended
OIT's weight function is not a clean per-layer transmittance product to
begin with (see "Technique" above) — so `K * spacing` is a plausible
estimate, not a derived guarantee. `tests/headless.rs`'s
`glass_self_occlusion_culling_stays_within_a_measured_pixel_error_bound`
measures the actual pixel error against the uncontrolled render on a
dense synthetic structure instead of trusting the derivation (observed:
mean 2.1/255, max 108/255 at 0.35 opacity, RTX 3060 laptop) and is the
thing to rerun if the formula ever changes.

**Known limitation**: unlike the opaque Hi-Z, this reference has no
phase-2 revalidation — a glass atom wrongly culled because last frame's
`glass_depth` suggested deeper coverage than is actually there this
frame stays wrong for the whole frame, not just until phase 2 catches
it, so a moving camera can show one frame of drift before the next
frame's copy corrects it (the extra layer of slack in the margin exists
to absorb exactly this). `glass_self_occlusion_culling_has_no_popping_
across_a_180_degree_sweep` bounds how much *extra* frame-to-frame
jumpiness culling adds on top of the baseline camera motion already
causes (observed: +1.6/255 mean over a 36-step half turn) rather than
asserting an absolute bound, since ordinary parallax already moves
plenty of pixels frame to frame on its own. A scene that gains its first
glass item this frame, or one where the previous frame drew no glass at
that pixel, reads a `0` reference and never culls (safe default, not a
special case in the shader).

Only `occluded` firing sets an atom's `occluded_bits` for phase 2 to
re-test (`cull`/`cull_bonds`'s `opaque_hidden`, `glass_only_hidden` in
`shaders/cull.wgsl`): a glass-culled atom would just read the same
unchanged `glass_hiz` again for the same answer, so re-testing it wastes
a phase-2 dispatch entry for nothing.

**Interaction-time only.** The isotropic-density approximation above
does not model real packing: spheres leave gaps between them a ray can
pass through, so "cull everything past `K` layers of mean spacing" drops
some real, visible contribution from atoms the isotropic estimate placed
too deep, and the composited pixel reads darker and more saturated than
the exact render wherever that happens (visible on a 1M-atom scene at
`close`: mean 4.3/255, 32,125 of 2,073,600 pixels off by more than
32/255, comparing `benchmarks/images/still/close.png` and `interacting/
close.png` -- `vv-bench --atoms 1000000 --material transparent --preset
close --png`, `--frames 10`, the 10th frame after a cold start, so
`glass_hiz` is warmed up; a different frame count or warm-up length
will read a different number, which is why this and the coordinator's
independently observed 8,342 need not match exactly). That is too
visible for a still image, so it is gated behind
`RenderSettings::fast_glass` (default `false`, independent of
`occlusion_culling`, which still runs its exact opaque test either way):
`false` computes a `0` margin for every item regardless of material and
skips the "glass hiz" copy entirely, so a still view, a screenshot or
export, and any direct `Renderer::render`/`render_all` caller that
doesn't opt in all render exactly, at zero extra cost over plain
occlusion culling. `vv-app` (`State::render_scene` in `crates/vv-app/
src/state.rs`) sets it from `self.camera != self.previous_camera ||
self.timeline.playing` — the camera's own equality against last frame's,
not a separate drag/animate flag ui.rs would have to maintain, since
`Camera` already carries every knob (pan, zoom, dolly, rotation,
projection) a still view needs to hold steady — and always `false` for
`export_screenshot`'s one-shot render. `glass_self_occlusion_culling_
is_exact_when_fast_glass_is_off` covers the one sharp edge this needs to
get right: turning `fast_glass` off must stop approximating immediately,
not once a stale "glass hiz" reference from an earlier `fast_glass: true`
frame happens to get overwritten — the gate is on the margin computation
itself, checked by priming that reference and then asserting the next
render matches the uncontrolled one (within OIT's own draw-order float
noise, ~0.05/255 mean; see the test's own comment for why weighted-
blended OIT is never quite bit-exact between two renders of a busy
overlap regardless of this feature).

Effect at the same 1M-atom benchmark (RTX 3060 laptop, p95 ms,
`overview`/`mid`/`close`/`inside`), `fast_glass: true`: 15.97/26.75/
30.62/25.38 -> 4.59/7.36/6.97/7.96 (`visible_quads` 1,000,000 ->
134,217/72,407/51,074/10,497, confirming the atom count is what dropped,
not just its cost per pixel). `fast_glass: false` (the default) matches
the *before* numbers and `visible_quads` exactly, by construction. See
benchmarks/README.md's "Transparency" entry for the full before/after
table on both GPUs.

**Picking and the selection outline.** A click resolves to the nearest
*opaque* thing under the cursor if there is one, or the nearest
*transparent* one only if nothing opaque is there at all (`id_view`
stays exactly what the opaque pass wrote, everywhere it drew) — not
"whichever is nearer": an opaque atom behind a glass one is still what
gets picked, matching the pre-existing Skin/SES convention ("clicks pick
what is behind"), now extended so a glass-only region (nothing opaque
behind it) is finally pickable too, instead of always missing. The
selection outline (`selection.rs`) needs no change for this: it already
just reads whatever is in `id_view`. A glass Gaussian/Skin/SES surface
stays click-through, as before.

**The path tracer (`render`, `shaders/path_trace.wgsl`) already handled
opacity correctly and needed no change**: every camera ray tests
`random() >= material_alpha(...)` at each hit and, on failure, continues
straight through with no refraction (unbiased stochastic transmission,
not an approximation), for every primitive it hits — spheres, cylinders,
triangles, SES and skin patches alike. The rasterized
`glass()` uses `material_alpha` directly, with no Fresnel brightening at
grazing angles, so both renderers agree.

## CPU backend (`vv-cpu`)

Not part of the v1 GPU pipeline above — a separate, much simpler renderer
for machines with no usable GPU at all (an HPC login node, a headless box
without llvmpipe/WARP, a laptop with a broken driver), not for everyday
use. Scalar ray-sphere intersection, one thread per image row via `rayon`,
a bounding-volume hierarchy (`vv-cpu/src/bvh.rs`, median-split top-down
over atom bounding spheres, no SAH) for nearest-hit queries. Spacefill
only so far — no bonds, no occlusion culling, no FXAA/SSAA.

Deliberately reuses `vv_render::{Camera, colors_for, StylePreset}` rather
than re-deriving camera or color math, so a CPU and GPU spacefill render
of the same file and camera agree (checked directly: `vv-cpu/tests/
agreement.rs` renders both and asserts >=95% pixel agreement on
background vs. drawn). Shading (`vv-cpu/src/renderer.rs`'s `shade()`) is a
hand-written twin of `shade()` in `draw.wgsl`, not a shared function —
WGSL and Rust can't share code, so if you tune one, check the other or
the two backends will silently drift apart.

Brute force (no BVH, every atom tested against every pixel) measured 2.0
fps for 4HHB (4,779 atoms) at 512x512 on a 16-thread laptop CPU — 5x short
of the target below, which is what motivated the BVH. Measured results
(`vv-cpu/tests/bench.rs`, run manually: `cargo test -p vv-cpu --release
--test bench -- --ignored --nocapture`):

| Structure | Atoms | Size | fps (16 threads) |
|---|---|---|---|
| 4HHB | 4,779 | 512x512 | ~57 (target: >=10) |
| 3J3Q | 2,440,800 | 512x512 | ~18 (informational, no target) |

The BVH's O(log n)-per-pixel scaling is why 3J3Q (510x more atoms than
4HHB) only costs ~3x more time, not 510x — brute force would have been
unusable at that size. Not a substitute for the GPU pipeline (10M atoms
at 10-14 ms there, not 512x512 at 512x512), and not SIMD-widened; revisit
the scalar inner loop if a real use case needs more than this gives.

## Known limits

- Occlusion uses last frame's depth with a half-radius safety margin; a
  camera jump gives one conservative (slower) frame, and fast rotation
  inside a dense structure shows p95 spikes (`inside` benchmark preset).
- A structure with no b-factors recorded (or all equal, e.g. every
  synthetic/benchmark structure) colors solid white under `ColorScheme::
  BFactor` rather than erroring — `by_b_factor`'s `min == max` case. Not a
  bug, but worth knowing before wondering why "by b-factor" looks blank.
- `picking_agrees_with_the_color_buffer_under_occlusion_culling` and
  `picking_hits_the_framed_center_and_misses_the_empty_corner` have each
  shown the same kind of intermittent single-pixel Hi-Z boundary flake as
  `occlusion_culling_does_not_change_the_image` above (observed on NVIDIA
  only, roughly 1 run in 5-10 when the full `--test headless` suite runs
  back to back, always a viewport-corner background check, never
  reproduces in isolation). Neither is given a tolerance, since both are
  exactly the tests that would catch picking or frustum-cull coverage
  actually disagreeing with what is drawn — loosen one only with real
  evidence the two concerns can't be told apart, not just to silence the
  flake. If this keeps recurring as more headless tests are added, it may
  be worth an adapter-level fix (e.g. serializing GPU-heavy tests) rather
  than chasing it test by test.
