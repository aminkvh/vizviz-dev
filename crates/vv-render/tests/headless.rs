//! End-to-end: upload a synthetic structure, cull + draw offscreen, read
//! pixels back. Skips (passes) when the machine has no usable GPU adapter.

use rayon::prelude::*;
use vv_io::synth::{protein_like, SynthParams};
use vv_render::{
    AtomSizes, Camera, CartoonGpu, CartoonItem, CartoonMesh, ColorScheme, DrawItem,
    GaussianSurfaceGpu, GaussianSurfaceItem, GpuContext, GpuStructure, PatchSurface,
    PatchSurfaceItem, Pick, Projection, RenderSettings, Renderer, Representation, SkinSurfaceGpu,
    StylePreset,
};

/// Skin-surface weights for a structure, as production derives them
/// (`vv_app::gpu_cache::build_skin_surface`): `r_vw^2 / shrink`, so an
/// isolated atom's patch lands on its true van der Waals radius.
fn skin_weights(structure: &vv_core::Structure) -> Vec<f32> {
    structure
        .topology
        .element
        .iter()
        .map(|e| {
            vv_core::skin_surface::weight_for_radius(
                e.vdw_radius(),
                vv_core::skin_surface::DEFAULT_SHRINK,
            )
        })
        .collect()
}

/// Set `VIZVIZ_TEST_ADAPTER=<name substring>` to run these on a specific GPU.
fn context() -> Option<std::sync::Arc<GpuContext>> {
    let filter = std::env::var("VIZVIZ_TEST_ADAPTER").ok();
    let instance = vv_render::GpuContext::instance();
    match GpuContext::with_instance(instance, None, filter.as_deref()) {
        Ok(ctx) => Some(ctx),
        Err(e) => {
            eprintln!("skipping GPU test: {e}");
            None
        }
    }
}

fn render_with(
    atoms: usize,
    distance_scale: f32,
    occlusion: bool,
    frames: usize,
) -> Option<(Vec<u8>, u32, u32)> {
    let ctx = context()?;
    let (w, h) = (256, 256);
    let mut renderer = Renderer::new(ctx.clone(), w, h);
    let structure = protein_like(&SynthParams::new(atoms));
    let gpu = GpuStructure::upload(&ctx, &structure, None, ColorScheme::Element);
    let bindings = renderer.bind(&gpu);
    let mut camera = Camera::framing(gpu.center, gpu.radius);
    camera.distance *= distance_scale;

    let settings = RenderSettings {
        occlusion_culling: occlusion,
        ..Default::default()
    };
    for _ in 0..frames {
        let mut encoder = ctx.device.create_command_encoder(&Default::default());
        renderer.render(&mut encoder, &camera, &gpu, &bindings, &settings);
        ctx.queue.submit([encoder.finish()]);
        renderer.after_submit();
    }
    Some((renderer.read_color(), w, h))
}

fn render_once(atoms: usize, distance_scale: f32) -> Option<(Vec<u8>, u32, u32)> {
    render_with(atoms, distance_scale, true, 2)
}

#[test]
fn picking_hits_the_framed_center_and_misses_the_empty_corner() {
    let Some(ctx) = context() else { return };
    let (w, h) = (256, 256);
    let mut renderer = Renderer::new(ctx.clone(), w, h);
    let structure = protein_like(&SynthParams::new(20_000));
    let gpu = GpuStructure::upload(&ctx, &structure, None, ColorScheme::Element);
    let bindings = renderer.bind(&gpu);
    let camera = Camera::framing(gpu.center, gpu.radius);
    let settings = RenderSettings {
        occlusion_culling: false,
        ..Default::default()
    };
    let mut encoder = ctx.device.create_command_encoder(&Default::default());
    renderer.render(&mut encoder, &camera, &gpu, &bindings, &settings);
    ctx.queue.submit([encoder.finish()]);
    renderer.after_submit();

    let center = renderer.pick(w / 2, h / 2);
    match center {
        Some(Pick::Atom { item: 0, atom }) => assert!((atom as usize) < gpu.atom_count),
        other => panic!("center of a framed structure should hit an atom, got {other:?}"),
    }

    let corner = renderer.pick(2, 2);
    assert_eq!(
        corner, None,
        "far corner of a framed structure should be background"
    );
}

/// The scenario the occlusion-cull design has to get right: picking must
/// never report a hit where nothing is drawn, or a miss where something
/// is — including under occlusion culling and heavy overdraw, where a
/// naive "re-run cull for picking" design could disagree with what the
/// user actually sees on screen. Reusing the frame's own compacted
/// buffers (see `Renderer::pick`) makes this true by construction; this
/// test is what would catch it if that stopped being the case.
#[test]
fn picking_agrees_with_the_color_buffer_under_occlusion_culling() {
    let Some(ctx) = context() else { return };
    let (w, h) = (200, 150);
    let mut renderer = Renderer::new(ctx.clone(), w, h);
    let structure = protein_like(&SynthParams::new(50_000));
    let gpu = GpuStructure::upload(&ctx, &structure, None, ColorScheme::Element);
    let bindings = renderer.bind(&gpu);
    // `Camera::framing` already sets the distance that exactly fits the
    // structure (~2.74x its radius); *scale* it, don't replace it, or the
    // camera ends up far closer than intended and the viewport corners —
    // which `background_of` assumes are background — are no longer empty.
    let mut camera = Camera::framing(gpu.center, gpu.radius);
    camera.distance *= 1.1;
    // Magenta: an unlit atom can come out as dark as a dark background.
    let settings = RenderSettings {
        background: PINHOLE_BG,
        ..Default::default()
    }; // occlusion_culling: true
    for _ in 0..3 {
        // A few frames so the depth pyramid is warm and occlusion is live.
        let mut encoder = ctx.device.create_command_encoder(&Default::default());
        renderer.render(&mut encoder, &camera, &gpu, &bindings, &settings);
        ctx.queue.submit([encoder.finish()]);
        renderer.after_submit();
    }
    let pixels = renderer.read_color();
    let bg = background_of(&pixels);

    let (mut checked_hit, mut checked_miss) = (false, false);
    for y in (0..h).step_by(7) {
        for x in (0..w).step_by(7) {
            let i = ((y * w + x) * 4) as usize;
            let is_bg = (pixels[i] as i32 - bg[0] as i32).abs() <= 8
                && (pixels[i + 1] as i32 - bg[1] as i32).abs() <= 8
                && (pixels[i + 2] as i32 - bg[2] as i32).abs() <= 8;
            let picked = renderer.pick(x, y);
            if is_bg {
                assert_eq!(
                    picked, None,
                    "background pixel ({x},{y}) should not pick an atom"
                );
                checked_miss = true;
            } else {
                assert!(
                    picked.is_some(),
                    "drawn pixel ({x},{y}) should pick an atom"
                );
                checked_hit = true;
            }
        }
    }
    assert!(
        checked_hit && checked_miss,
        "test should exercise both background and drawn pixels"
    );
}

#[test]
fn ball_and_stick_draws_bond_cylinders() {
    let Some(ctx) = context() else { return };
    let structure = protein_like(&SynthParams::new(3_000));
    let bonds = vv_core::bonds::perceive(&structure.topology, structure.frame(0).positions());
    assert!(
        bonds.len() > 2_000,
        "synthetic residues should be internally bonded"
    );
    let render = |with_bonds: bool| {
        let mut renderer = Renderer::new(ctx.clone(), 256, 256);
        let gpu = GpuStructure::upload(
            &ctx,
            &structure,
            with_bonds.then_some(&bonds),
            ColorScheme::Element,
        );
        let bindings = renderer.bind(&gpu);
        let mut camera = Camera::framing(gpu.center, gpu.radius);
        camera.distance *= 0.7;
        let settings = RenderSettings {
            representation: Representation::BallAndStick,
            occlusion_culling: false,
            ..Default::default()
        };
        let mut encoder = ctx.device.create_command_encoder(&Default::default());
        renderer.render(&mut encoder, &camera, &gpu, &bindings, &settings);
        ctx.queue.submit([encoder.finish()]);
        renderer.after_submit();
        let counts = gpu.read_visible_counts(&ctx);
        (
            renderer.read_color(),
            counts.iter().map(|c| c.bonds).sum::<u32>(),
        )
    };
    let (without, bonds_without) = render(false);
    let (with, bonds_with) = render(true);
    assert_eq!(bonds_without, 0);
    assert!(bonds_with > 1_000, "visible bonds: {bonds_with}");
    let bg = background_of(&without);
    let drawn_without = count_non_background(&without, bg);
    let drawn_with = count_non_background(&with, bg);
    // Small spheres in a dense blob already overlap, so the gain is modest.
    assert!(
        drawn_with > drawn_without + drawn_without / 20,
        "cylinders should add coverage: {drawn_without} -> {drawn_with}"
    );
}

/// `GpuStructure::recolor` writes straight into `Page.colors` (which
/// carries `COPY_DST` for exactly this) instead of re-uploading — this is
/// what would catch a regression back to needing a full re-upload.
#[test]
fn recoloring_changes_the_drawn_pixels_without_reupload() {
    let Some(ctx) = context() else { return };
    let (w, h) = (256, 256);
    let mut renderer = Renderer::new(ctx.clone(), w, h);
    let structure = protein_like(&SynthParams {
        residues_per_chain: 30,
        ..SynthParams::new(6_000)
    });
    let gpu = GpuStructure::upload(&ctx, &structure, None, ColorScheme::Element);
    let bindings = renderer.bind(&gpu);
    let camera = Camera::framing(gpu.center, gpu.radius);
    let settings = RenderSettings {
        occlusion_culling: false,
        ..Default::default()
    };
    let render = |renderer: &mut Renderer| {
        let mut encoder = ctx.device.create_command_encoder(&Default::default());
        renderer.render(&mut encoder, &camera, &gpu, &bindings, &settings);
        ctx.queue.submit([encoder.finish()]);
        renderer.after_submit();
        renderer.read_color()
    };

    let by_element = render(&mut renderer);
    gpu.recolor(&ctx, &structure.topology, ColorScheme::Chain);
    let by_chain = render(&mut renderer);

    let bg = background_of(&by_element);
    let changed = by_element
        .chunks(4)
        .zip(by_chain.chunks(4))
        .filter(|(a, b)| {
            let a_is_bg = (a[0] as i32 - bg[0] as i32).abs() <= 8
                && (a[1] as i32 - bg[1] as i32).abs() <= 8
                && (a[2] as i32 - bg[2] as i32).abs() <= 8;
            !a_is_bg && a != b
        })
        .count();
    assert!(
        changed > 100,
        "recoloring by chain should change many drawn pixels, changed {changed}"
    );
}

/// The style uniform actually has to reach the fragment shader — a preset
/// that changes nothing visible is a failure `style::tests` can't see,
/// since those only check the Rust-side params differ, not that the GPU
/// used them.
#[test]
fn style_presets_change_the_drawn_pixels() {
    let Some(ctx) = context() else { return };
    let (w, h) = (256, 256);
    let structure = protein_like(&SynthParams::new(4_000));
    let gpu = GpuStructure::upload(&ctx, &structure, None, ColorScheme::Element);
    let mut camera = Camera::framing(gpu.center, gpu.radius);
    camera.distance *= 0.3; // close enough that shading dominates the image
    let render = |style: StylePreset| {
        let mut renderer = Renderer::new(ctx.clone(), w, h);
        let bindings = renderer.bind(&gpu);
        let settings = RenderSettings {
            occlusion_culling: false,
            lighting: style.lighting(),
            material: style.material(),
            background: style.background(),
            ..Default::default()
        };
        let mut encoder = ctx.device.create_command_encoder(&Default::default());
        renderer.render(&mut encoder, &camera, &gpu, &bindings, &settings);
        ctx.queue.submit([encoder.finish()]);
        renderer.after_submit();
        renderer.read_color()
    };

    let dark = render(StylePreset::DarkPresentation);
    let white = render(StylePreset::PublicationWhite);
    let cel = render(StylePreset::FlatCel);
    let glossy = render(StylePreset::Glossy);

    let differing =
        |a: &[u8], b: &[u8]| a.chunks(4).zip(b.chunks(4)).filter(|(x, y)| x != y).count();
    let total = (w * h) as usize;
    assert!(
        differing(&dark, &white) * 10 > total,
        "publication-white should visibly change most of the frame vs. dark presentation"
    );
    assert!(
        differing(&dark, &cel) * 10 > total,
        "flat/cel should visibly change most of the frame vs. dark presentation"
    );
    // Dark Presentation and Glossy share the same background, so any
    // difference here can only come from `shade()` actually reading the
    // style uniform — the other two comparisons above would also pass if
    // only the background changed and the shader ignored it entirely.
    assert!(
        differing(&dark, &glossy) * 10 > total,
        "glossy should visibly change most of the frame vs. dark presentation, despite sharing a background"
    );
}

/// Exercises the exact pipeline `vv-app`'s screenshot export uses (render
/// at 2x, downsample) end to end -- there's no automated way to drive the
/// native app's own export dialog, so this is the closest thing to a
/// regression test for that feature: it would catch the render-at-2x and
/// the downsample disagreeing about resolution, or the SSAA image coming
/// out misaligned/garbled relative to a normal 1x render of the same
/// scene.
#[test]
fn ssaa_render_and_downsample_agrees_with_a_direct_render() {
    let Some(ctx) = context() else { return };
    let (w, h) = (128, 96);
    let structure = protein_like(&SynthParams::new(20_000));
    let gpu = GpuStructure::upload(&ctx, &structure, None, ColorScheme::Element);
    let camera = Camera::framing(gpu.center, gpu.radius);
    let settings = RenderSettings {
        occlusion_culling: false,
        ..Default::default()
    };

    let mut direct = Renderer::new(ctx.clone(), w, h);
    let direct_bindings = direct.bind(&gpu);
    let mut encoder = ctx.device.create_command_encoder(&Default::default());
    direct.render(&mut encoder, &camera, &gpu, &direct_bindings, &settings);
    ctx.queue.submit([encoder.finish()]);
    direct.after_submit();
    let direct_pixels = direct.read_color();

    let (ew, eh) = (w * 2, h * 2);
    let mut big = Renderer::new(ctx.clone(), ew, eh);
    let big_bindings = big.bind(&gpu);
    let mut encoder = ctx.device.create_command_encoder(&Default::default());
    big.render(&mut encoder, &camera, &gpu, &big_bindings, &settings);
    ctx.queue.submit([encoder.finish()]);
    big.after_submit();
    let downsampled = vv_render::downsample_2x_srgb(&big.read_color(), ew, eh);

    assert_eq!(downsampled.len(), direct_pixels.len());
    let bg = background_of(&direct_pixels);
    let agree = direct_pixels
        .chunks(4)
        .zip(downsampled.chunks(4))
        .filter(|(a, b)| {
            let a_is_bg = (a[0] as i32 - bg[0] as i32).abs() <= 8
                && (a[1] as i32 - bg[1] as i32).abs() <= 8
                && (a[2] as i32 - bg[2] as i32).abs() <= 8;
            let b_is_bg = (b[0] as i32 - bg[0] as i32).abs() <= 8
                && (b[1] as i32 - bg[1] as i32).abs() <= 8
                && (b[2] as i32 - bg[2] as i32).abs() <= 8;
            a_is_bg == b_is_bg
        })
        .count();
    let total = direct_pixels.len() / 4;
    assert!(
        agree * 100 >= total * 95,
        "SSAA downsample should mostly agree with a direct render on which \
         pixels are background vs. drawn ({agree}/{total} agreed) -- a \
         resolution or alignment bug would show up as widespread disagreement, \
         not just soft edges"
    );
}

#[test]
fn framed_structure_survives_frustum_cull_entirely() {
    let Some(ctx) = context() else { return };
    for atoms in [1_000usize, 70_000, 600_000] {
        let mut renderer = Renderer::new(ctx.clone(), 320, 200);
        let structure = protein_like(&SynthParams::new(atoms));
        let gpu = GpuStructure::upload(&ctx, &structure, None, ColorScheme::Element);
        let bindings = renderer.bind(&gpu);
        let camera = Camera::framing(gpu.center, gpu.radius);
        let settings = RenderSettings {
            occlusion_culling: false,
            ..Default::default()
        };
        let mut encoder = ctx.device.create_command_encoder(&Default::default());
        renderer.render(&mut encoder, &camera, &gpu, &bindings, &settings);
        ctx.queue.submit([encoder.finish()]);
        renderer.after_submit();
        let counts = gpu.read_visible_counts(&ctx);
        let total: u32 = counts.iter().map(|c| c.quads + c.points).sum();
        assert_eq!(
            total as usize, atoms,
            "frustum cull dropped atoms ({counts:?})"
        );
    }
}

/// Prints which atoms the cull pass loses on this adapter. Run with
/// `cargo test -p vv-render --test headless -- --ignored --nocapture`.
#[test]
#[ignore]
fn diagnose_cull_coverage() {
    let Some(ctx) = context() else { return };
    let atoms = 70_000usize;
    let mut renderer = Renderer::new(ctx.clone(), 320, 200);
    let structure = protein_like(&SynthParams::new(atoms));
    let gpu = GpuStructure::upload(&ctx, &structure, None, ColorScheme::Element);
    let bindings = renderer.bind(&gpu);
    let camera = Camera::framing(gpu.center, gpu.radius);
    let settings = RenderSettings {
        occlusion_culling: false,
        ..Default::default()
    };
    let mut encoder = ctx.device.create_command_encoder(&Default::default());
    renderer.render(&mut encoder, &camera, &gpu, &bindings, &settings);
    ctx.queue.submit([encoder.finish()]);
    renderer.after_submit();

    let (quads, points) = gpu.read_visible_indices(&ctx).remove(0);
    let mut seen = vec![0u32; atoms];
    for &i in quads.iter().chain(points.iter()) {
        seen[i as usize] += 1;
    }
    let missing: Vec<usize> = (0..atoms).filter(|&i| seen[i] == 0).collect();
    let duplicates = seen.iter().filter(|&&c| c > 1).count();
    let mut by_lane = [0u32; 256];
    let mut workgroups = std::collections::BTreeMap::<usize, u32>::new();
    for &i in &missing {
        by_lane[i % 256] += 1;
        *workgroups.entry(i / 256).or_default() += 1;
    }
    eprintln!(
        "adapter {}: {} atoms, {} quads + {} points listed, {} missing, {} duplicated",
        ctx.adapter_name(),
        atoms,
        quads.len(),
        points.len(),
        missing.len(),
        duplicates
    );
    eprintln!("first missing: {:?}", &missing[..missing.len().min(24)]);
    let lanes: Vec<(usize, u32)> = by_lane
        .iter()
        .enumerate()
        .filter(|(_, &c)| c > 0)
        .map(|(l, &c)| (l, c))
        .collect();
    eprintln!(
        "missing by lane (lane, count), {} lanes affected: {:?}",
        lanes.len(),
        &lanes[..lanes.len().min(40)]
    );
    let full: Vec<_> = workgroups
        .iter()
        .filter(|(_, &c)| c == 256)
        .map(|(w, _)| *w)
        .collect();
    eprintln!(
        "workgroups affected: {} of {}, entirely missing: {} {:?}",
        workgroups.len(),
        atoms.div_ceil(256),
        full.len(),
        &full[..full.len().min(20)]
    );
    let partial: Vec<_> = workgroups
        .iter()
        .filter(|(_, &c)| c < 256)
        .take(12)
        .collect();
    eprintln!("partially missing (workgroup, count): {partial:?}");
}

#[test]
fn occlusion_culling_does_not_change_the_image() {
    // With a static camera the pyramid from the previous frame is exact,
    // so culling should be all but invisible. "All but": the Hi-Z test
    // compares against a min-reduced R32Float mip chain, so a silhouette
    // atom can land right on the float-precision boundary of "occluded"
    // and flip a handful of edge pixels — an inherent limit of a
    // conservative, approximate test, not a correctness bug. Bound it
    // tightly instead of requiring bit-exact equality.
    for (atoms, distance) in [(50_000, 1.0), (50_000, 0.45), (200_000, 0.6)] {
        let Some((culled, w, h)) = render_with(atoms, distance, true, 3) else {
            return;
        };
        let Some((reference, _, _)) = render_with(atoms, distance, false, 1) else {
            return;
        };
        let differing = culled
            .chunks(4)
            .zip(reference.chunks(4))
            .filter(|(a, b)| a != b)
            .count();
        let total = (w * h) as usize;
        assert!(
            differing * 10_000 <= total,
            "{differing} of {total} pixels differ with occlusion culling (atoms {atoms}, distance {distance}), expected at most 0.01%"
        );
    }
}

/// The moving-camera case two-phase occlusion culling exists for: the
/// pyramid from a frame at one angle is wrong for the next, and phase 1
/// alone culls atoms that just came into view (8GLV showed square holes
/// while orbiting). Phase 2 must bring the frame back to what rendering
/// without culling draws.
#[test]
fn occlusion_culling_survives_a_camera_jump() {
    let Some(ctx) = context() else {
        return;
    };
    let (w, h) = (256, 256);
    let structure = protein_like(&SynthParams::new(200_000));
    let gpu = GpuStructure::upload(&ctx, &structure, None, ColorScheme::Element);
    let base = Camera::framing(gpu.center, gpu.radius);
    let mut before = base.clone();
    before.distance *= 0.6;
    let mut after = before.clone();
    after.orbit(0.8, 0.3);

    let frame = |occlusion: bool, cameras: &[&Camera]| {
        let mut renderer = Renderer::new(ctx.clone(), w, h);
        let bindings = renderer.bind(&gpu);
        let settings = RenderSettings {
            occlusion_culling: occlusion,
            ..Default::default()
        };
        for camera in cameras {
            let mut encoder = ctx.device.create_command_encoder(&Default::default());
            renderer.render(&mut encoder, camera, &gpu, &bindings, &settings);
            ctx.queue.submit([encoder.finish()]);
            renderer.after_submit();
        }
        renderer.read_color()
    };
    let culled = frame(true, &[&before, &after]);
    let reference = frame(false, &[&after]);
    let differing = culled
        .chunks(4)
        .zip(reference.chunks(4))
        .filter(|(a, b)| a != b)
        .count();
    let total = (w * h) as usize;
    assert!(
        differing * 1_000 <= total,
        "{differing} of {total} pixels differ after a camera jump with occlusion culling, expected at most 0.1%"
    );
}

fn count_non_background(pixels: &[u8], background: [u8; 3]) -> usize {
    pixels
        .chunks(4)
        .filter(|p| {
            (p[0] as i32 - background[0] as i32).abs() > 8
                || (p[1] as i32 - background[1] as i32).abs() > 8
                || (p[2] as i32 - background[2] as i32).abs() > 8
        })
        .count()
}

fn background_of(pixels: &[u8]) -> [u8; 3] {
    // Corners are background for a framed, roughly spherical structure.
    [pixels[0], pixels[1], pixels[2]]
}

#[test]
fn framed_structure_covers_center_not_corners() {
    let Some((pixels, w, h)) = render_once(20_000, 1.0) else {
        return;
    };
    let bg = background_of(&pixels);
    let center = ((h / 2) * w + w / 2) as usize * 4;
    let c = &pixels[center..center + 3];
    assert!(
        (c[0] as i32 - bg[0] as i32).abs() > 8 || (c[1] as i32 - bg[1] as i32).abs() > 8,
        "center pixel {c:?} should not be background {bg:?}"
    );
    let drawn = count_non_background(&pixels, bg);
    let total = (w * h) as usize;
    assert!(drawn > total / 20, "only {drawn}/{total} pixels drawn");
    assert!(drawn < total, "structure should not cover the whole frame");
}

#[test]
fn zoomed_in_uses_sphere_impostors_with_shading() {
    // Close up, atoms are many pixels wide; shading makes pixels vary a lot.
    let Some((pixels, _, _)) = render_once(2_000, 0.15) else {
        return;
    };
    let bg = background_of(&pixels);
    let drawn = count_non_background(&pixels, bg);
    assert!(drawn > 1000, "expected large spheres, drew {drawn} pixels");
    let mut luminance: Vec<u32> = pixels
        .chunks(4)
        .map(|p| p[0] as u32 + p[1] as u32 + p[2] as u32)
        .collect();
    luminance.sort_unstable();
    let spread = luminance[luminance.len() * 95 / 100] - luminance[luminance.len() * 5 / 100];
    assert!(
        spread > 60,
        "expected shaded spheres, luminance spread {spread}"
    );
}

/// FXAA reads `display_view`, not `color_view` -- and this is the one test
/// that would catch a wrong-constant no-op or an over-eager everywhere-blur
/// (both are invisible from "the pass ran without erroring"): it checks
/// that turning FXAA on measurably softens the sharpest edge in a frame
/// full of aliased sphere silhouettes, while still leaving most pixels
/// untouched.
#[test]
fn fxaa_smooths_silhouette_edges_without_blurring_everywhere() {
    let Some(ctx) = context() else { return };
    let (w, h) = (256u32, 256u32);
    let structure = protein_like(&SynthParams::new(2_000));
    let gpu = GpuStructure::upload(&ctx, &structure, None, ColorScheme::Element);
    let mut camera = Camera::framing(gpu.center, gpu.radius);
    camera.distance *= 0.15; // large spheres, clearly aliased silhouettes

    let render = |fxaa: bool| {
        let mut renderer = Renderer::new(ctx.clone(), w, h);
        let bindings = renderer.bind(&gpu);
        let settings = RenderSettings {
            occlusion_culling: false,
            fxaa,
            ..Default::default()
        };
        let mut encoder = ctx.device.create_command_encoder(&Default::default());
        renderer.render(&mut encoder, &camera, &gpu, &bindings, &settings);
        ctx.queue.submit([encoder.finish()]);
        renderer.after_submit();
        renderer.read_display_color()
    };
    let off = render(false);
    let on = render(true);

    fn luma_at(pixels: &[u8], w: u32, x: u32, y: u32) -> i32 {
        let i = ((y * w + x) * 4) as usize;
        pixels[i] as i32 * 299 + pixels[i + 1] as i32 * 587 + pixels[i + 2] as i32 * 114
    }
    let max_horizontal_delta = |pixels: &[u8]| -> i32 {
        let mut max = 0;
        for y in 0..h {
            for x in 1..w {
                let delta = (luma_at(pixels, w, x, y) - luma_at(pixels, w, x - 1, y)).abs();
                max = max.max(delta);
            }
        }
        max
    };
    let delta_off = max_horizontal_delta(&off);
    let delta_on = max_horizontal_delta(&on);
    assert!(
        delta_on < delta_off,
        "FXAA should soften the sharpest edge in the frame: off={delta_off} on={delta_on}"
    );

    let total = (w * h) as usize;
    let differing = off
        .chunks(4)
        .zip(on.chunks(4))
        .filter(|(a, b)| a != b)
        .count();
    assert!(differing > 0, "FXAA should change at least some pixels");
    assert!(
        differing * 2 < total,
        "FXAA should only touch edges, not blur the whole frame: changed {differing}/{total}"
    );
}

#[test]
fn far_away_still_draws_points() {
    let Some((pixels, _, _)) = render_once(20_000, 6.0) else {
        return;
    };
    let bg = background_of(&pixels);
    let drawn = count_non_background(&pixels, bg);
    assert!(
        drawn > 10,
        "far structure should still show as points, drew {drawn}"
    );
}

fn shifted(s: &vv_core::Structure, by: glam::Vec3) -> vv_core::Structure {
    let positions = s.frame(0).positions().iter().map(|p| *p + by).collect();
    vv_core::Structure::new((*s.topology).clone(), vv_core::CoordSet::new(positions)).unwrap()
}

/// The defining property this whole feature is for, checked against real
/// rendered pixels (not just `Camera::proj`'s matrix math, already
/// covered in `vv_render::camera`'s own tests): two atoms of the same
/// element (same van der Waals radius) sitting at very different view
/// depths should draw at very different apparent sizes under perspective,
/// and at essentially the *same* size under orthographic.
#[test]
fn orthographic_does_not_shrink_atoms_with_depth() {
    let Some(ctx) = context() else { return };
    let base = protein_like(&SynthParams::new(8)); // one residue: N,CA,C,O,CB,CG,CD,CE
                                                   // Atoms 1 (CA) and 5 (CG) are both carbon, same radius -- see
                                                   // `vv_io::synth::ELEMENTS`. Offset in X so neither occludes the
                                                   // other on screen despite sharing X=0/Y=0 otherwise; everything else
                                                   // pushed far outside the frustum so it can't interfere.
    let hidden = glam::Vec3::splat(10_000.0);
    let near = glam::Vec3::new(-15.0, 0.0, 0.0);
    let far = glam::Vec3::new(15.0, 0.0, -60.0);
    let mut positions = vec![hidden; 8];
    positions[1] = near;
    positions[5] = far;
    let structure =
        vv_core::Structure::new((*base.topology).clone(), vv_core::CoordSet::new(positions))
            .unwrap();
    let gpu = GpuStructure::upload(&ctx, &structure, None, ColorScheme::Element);

    let (w, h) = (400u32, 400u32);
    let aspect = w as f32 / h as f32;
    // No yaw/pitch: the view matrix is a pure translation, so world Z is
    // exactly view depth, lining up "closer"/"farther from the camera"
    // with the near/far world positions above with no rotation to
    // account for.
    let camera = Camera::framing(glam::Vec3::new(0.0, 0.0, -30.0), 50.0);

    let pixel_of = |cam: &Camera, p: glam::Vec3| -> (i32, i32) {
        let clip = cam.proj(aspect) * cam.view() * glam::Vec4::new(p.x, p.y, p.z, 1.0);
        let ndc = clip.truncate() / clip.w;
        (
            ((ndc.x * 0.5 + 0.5) * w as f32) as i32,
            ((0.5 - ndc.y * 0.5) * h as f32) as i32,
        )
    };
    // Width, in pixels, of the contiguous non-background span on row `py`
    // containing column `px` -- a direct read of "how big does this atom
    // look," with no assumption about its exact projected radius.
    let span_width = |pixels: &[u8], py: i32, px: i32, bg: [u8; 3]| -> i32 {
        if py < 0 || py >= h as i32 || px < 0 || px >= w as i32 {
            return 0;
        }
        let row = py as usize * w as usize;
        let is_drawn = |x: i32| -> bool {
            let i = (row + x as usize) * 4;
            (pixels[i] as i32 - bg[0] as i32).abs() > 8
                || (pixels[i + 1] as i32 - bg[1] as i32).abs() > 8
                || (pixels[i + 2] as i32 - bg[2] as i32).abs() > 8
        };
        if !is_drawn(px) {
            return 0;
        }
        let mut lo = px;
        while lo > 0 && is_drawn(lo - 1) {
            lo -= 1;
        }
        let mut hi = px;
        while hi < w as i32 - 1 && is_drawn(hi + 1) {
            hi += 1;
        }
        hi - lo + 1
    };

    let width_at = |projection: Projection| -> (i32, i32) {
        let mut cam = camera.clone();
        cam.projection = projection;
        let mut renderer = Renderer::new(ctx.clone(), w, h);
        let bindings = renderer.bind(&gpu);
        let settings = RenderSettings {
            occlusion_culling: false,
            ..Default::default()
        };
        let mut encoder = ctx.device.create_command_encoder(&Default::default());
        renderer.render(&mut encoder, &cam, &gpu, &bindings, &settings);
        ctx.queue.submit([encoder.finish()]);
        renderer.after_submit();
        let pixels = renderer.read_color();
        let bg = background_of(&pixels);
        let (nx, ny) = pixel_of(&cam, near);
        let (fx, fy) = pixel_of(&cam, far);
        (
            span_width(&pixels, ny, nx, bg),
            span_width(&pixels, fy, fx, bg),
        )
    };

    let (persp_near, persp_far) = width_at(Projection::Perspective);
    assert!(
        persp_near > 2 && persp_far > 2,
        "expected both atoms drawn: {persp_near}, {persp_far}"
    );
    assert!(
        persp_near as f32 > persp_far as f32 * 1.2,
        "perspective should shrink the farther atom: near={persp_near} far={persp_far}"
    );

    let (ortho_near, ortho_far) = width_at(Projection::Orthographic);
    assert!(
        ortho_near > 2 && ortho_far > 2,
        "expected both atoms drawn: {ortho_near}, {ortho_far}"
    );
    assert!(
        (ortho_near - ortho_far).abs() as f32 <= ortho_near.max(ortho_far) as f32 * 0.15,
        "orthographic should size both atoms alike: near={ortho_near} far={ortho_far}"
    );
}

/// The same discard-nothing sanity check `framed_structure_covers_center_
/// not_corners` does for perspective, run once under orthographic so a
/// gross bug (wrong branch always taken, everything discarded, NaNs)
/// would fail loudly here instead of only showing up as a visual
/// complaint down the line.
#[test]
fn orthographic_renders_a_real_structure() {
    let Some(ctx) = context() else { return };
    let (w, h) = (256, 256);
    let mut renderer = Renderer::new(ctx.clone(), w, h);
    let structure = protein_like(&SynthParams::new(20_000));
    let gpu = GpuStructure::upload(&ctx, &structure, None, ColorScheme::Element);
    let bindings = renderer.bind(&gpu);
    let mut camera = Camera::framing(gpu.center, gpu.radius);
    camera.projection = Projection::Orthographic;
    let settings = RenderSettings {
        occlusion_culling: false,
        ..Default::default()
    };
    let mut encoder = ctx.device.create_command_encoder(&Default::default());
    renderer.render(&mut encoder, &camera, &gpu, &bindings, &settings);
    ctx.queue.submit([encoder.finish()]);
    renderer.after_submit();
    let pixels = renderer.read_color();
    let bg = background_of(&pixels);
    let drawn = count_non_background(&pixels, bg);
    let total = (w * h) as usize;
    assert!(drawn > total / 20, "only {drawn}/{total} pixels drawn");
    assert!(drawn < total, "structure should not cover the whole frame");
}

fn drawn_in_columns(
    pixels: &[u8],
    w: u32,
    h: u32,
    bg: [u8; 3],
    cols: std::ops::Range<u32>,
) -> usize {
    let mut n = 0;
    for y in 0..h {
        for x in cols.clone() {
            let i = ((y * w + x) * 4) as usize;
            if pixels[i] != bg[0] || pixels[i + 1] != bg[1] || pixels[i + 2] != bg[2] {
                n += 1;
            }
        }
    }
    n
}

/// Two structures in one frame: both must be drawn, and the picking
/// buffer must name each under its own draw-item index.
#[test]
fn two_structures_composite_into_one_frame_and_pick_by_item() {
    let Some(ctx) = context() else { return };
    let (w, h) = (256, 128);
    let mut renderer = Renderer::new(ctx.clone(), w, h);
    let left = protein_like(&SynthParams::new(20_000));
    let gpu_left = GpuStructure::upload(&ctx, &left, None, ColorScheme::Element);
    // Same blob shifted along +x by three radii so the two never overlap.
    let shift = glam::Vec3::new(gpu_left.radius * 3.0, 0.0, 0.0);
    let right = shifted(&left, shift);
    let gpu_right = GpuStructure::upload(&ctx, &right, None, ColorScheme::Chain);
    let bind_left = renderer.bind(&gpu_left);
    let bind_right = renderer.bind(&gpu_right);
    let items = [
        DrawItem {
            structure: &gpu_left,
            bindings: &bind_left,
            representation: Representation::Spacefill,
            sizes: AtomSizes::of(Representation::Spacefill),
            material: vv_render::Material::default(),
        },
        DrawItem {
            structure: &gpu_right,
            bindings: &bind_right,
            representation: Representation::Spacefill,
            sizes: AtomSizes::of(Representation::Spacefill),
            material: vv_render::Material::default(),
        },
    ];
    // The union sphere: centred between the two, radius 2.5x one blob's.
    let camera = Camera::framing(gpu_left.center + shift * 0.5, gpu_left.radius * 2.5);
    let settings = RenderSettings {
        occlusion_culling: false,
        ..Default::default()
    };
    let mut encoder = ctx.device.create_command_encoder(&Default::default());
    renderer.render_all(&mut encoder, &camera, &items, &[], &[], &[], &settings);
    ctx.queue.submit([encoder.finish()]);
    renderer.after_submit();

    let pixels = renderer.read_color();
    let bg = background_of(&pixels);
    assert!(
        drawn_in_columns(&pixels, w, h, bg, 0..w / 2) > 500,
        "left half drawn"
    );
    assert!(
        drawn_in_columns(&pixels, w, h, bg, w / 2..w) > 500,
        "right half drawn"
    );

    // Each blob's centre projects 1.5 radii either side of the midpoint;
    // the framed 2.5-radius sphere spans the viewport height (~64 px).
    let px_per_radius = (h as f32 / 2.0) / 2.5 / 1.05;
    let dx = (1.5 * px_per_radius) as u32;
    match renderer.pick(w / 2 - dx, h / 2) {
        Some(Pick::Atom { item: 0, atom }) => assert!((atom as usize) < gpu_left.atom_count),
        other => panic!("left blob should pick as item 0, got {other:?}"),
    }
    match renderer.pick(w / 2 + dx, h / 2) {
        Some(Pick::Atom { item: 1, atom }) => assert!((atom as usize) < gpu_right.atom_count),
        other => panic!("right blob should pick as item 1, got {other:?}"),
    }
    assert_eq!(renderer.pick(1, 1), None);
}

/// A bond cylinder under the cursor picks as the bond, mapped back to its
/// index in the `BondTable` the structure was uploaded with; the atoms at
/// its ends still pick as atoms.
#[test]
fn bond_picking_names_the_bond_under_the_cursor() {
    let Some(ctx) = context() else { return };
    let pdb = b"ATOM      1  C1  LIG A   1       0.000   0.000   0.000  1.00  0.00           C  \n\
                ATOM      2  C2  LIG A   1       1.500   0.000   0.000  1.00  0.00           C  \n\
                END\n";
    let structure = vv_io::parse(pdb, vv_io::Format::Pdb).unwrap();
    let bonds = vv_core::bonds::perceive(&structure.topology, structure.frame(0).positions());
    assert_eq!(bonds.pairs, vec![[0, 1]]);

    let (w, h) = (256, 256);
    let mut renderer = Renderer::new(ctx.clone(), w, h);
    let gpu = GpuStructure::upload(&ctx, &structure, Some(&bonds), ColorScheme::Element);
    let bindings = renderer.bind(&gpu);
    // Looking straight down -z at the bond's midpoint. Ball-and-stick
    // spheres are a quarter of the vdW radius (0.425 A for carbon), so the
    // 0.75 A gap to the midpoint shows only the cylinder there.
    let camera = Camera::framing(gpu.center, gpu.radius);
    let settings = RenderSettings {
        representation: Representation::BallAndStick,
        occlusion_culling: false,
        ..Default::default()
    };
    let mut encoder = ctx.device.create_command_encoder(&Default::default());
    renderer.render(&mut encoder, &camera, &gpu, &bindings, &settings);
    ctx.queue.submit([encoder.finish()]);
    renderer.after_submit();

    assert_eq!(
        renderer.pick(w / 2, h / 2),
        Some(Pick::Bond { item: 0, bond: 0 }),
        "midpoint of the bond should pick the cylinder"
    );
    // 0.75 A (one framed radius) left of centre is the first atom's centre;
    // its sphere (0.425 A) sticks out in front of the 0.15 A cylinder.
    let px_per_angstrom = (h as f32 / 2.0) / gpu.radius / 1.05;
    let atom_x = w / 2 - (0.6 * px_per_angstrom) as u32;
    assert_eq!(
        renderer.pick(atom_x, h / 2),
        Some(Pick::Atom { item: 0, atom: 0 }),
        "inside the first atom's sphere should pick that atom"
    );
}

/// `GpuStructure::set_frame` rewrites positions in place: the next frame
/// shows the new coordinate set, and going back restores the old image
/// exactly.
#[test]
fn set_frame_moves_the_drawn_atoms() {
    let Some(ctx) = context() else { return };
    let (w, h) = (128, 128);
    let base = protein_like(&SynthParams::new(5_000));
    let radius = base.frame(0).bounding_sphere().unwrap().1;
    // Frame 1 is the same blob three radii to the right: fully out of a
    // view framed on frame 0.
    let moved: Vec<glam::Vec3> = base
        .frame(0)
        .positions()
        .iter()
        .map(|p| *p + glam::Vec3::new(radius * 3.0, 0.0, 0.0))
        .collect();
    let structure = vv_core::Structure::with_frames(
        (*base.topology).clone(),
        vec![(*base.frame(0)).clone(), vv_core::CoordSet::new(moved)],
    )
    .unwrap();

    let mut renderer = Renderer::new(ctx.clone(), w, h);
    let gpu = GpuStructure::upload(&ctx, &structure, None, ColorScheme::Element);
    let bindings = renderer.bind(&gpu);
    let camera = Camera::framing(gpu.center, gpu.radius);
    let settings = RenderSettings {
        occlusion_culling: false,
        ..Default::default()
    };
    let render = |renderer: &mut Renderer| {
        let mut encoder = ctx.device.create_command_encoder(&Default::default());
        renderer.render(&mut encoder, &camera, &gpu, &bindings, &settings);
        ctx.queue.submit([encoder.finish()]);
        renderer.after_submit();
        renderer.read_color()
    };
    let before = render(&mut renderer);
    gpu.set_frame(&ctx, &structure, 1);
    let after = render(&mut renderer);
    gpu.set_frame(&ctx, &structure, 0);
    let back = render(&mut renderer);

    let bg = background_of(&before);
    let drawn_before = count_non_background(&before, bg);
    let drawn_after = count_non_background(&after, bg);
    assert!(drawn_before > 500, "frame 0 is in view: {drawn_before}");
    assert!(
        drawn_after < drawn_before / 10,
        "frame 1 is out of view: {drawn_after} vs {drawn_before}"
    );
    assert_eq!(back, before, "switching back restores frame 0 exactly");
}

fn luma(pixels: &[u8], i: usize) -> f32 {
    let p = &pixels[i * 4..i * 4 + 3];
    0.299 * p[0] as f32 + 0.587 * p[1] as f32 + 0.114 * p[2] as f32
}

/// The outline pass darkens edge pixels (silhouettes and creases between
/// atoms) and leaves everything else exactly alone; with outlines off it
/// is a byte-for-byte copy of the opaque output.
#[test]
fn outlines_darken_edges_and_nothing_else() {
    let Some(ctx) = context() else { return };
    let (w, h) = (256, 256);
    let structure = protein_like(&SynthParams::new(3_000));
    let gpu = GpuStructure::upload(&ctx, &structure, None, ColorScheme::Element);
    let mut camera = Camera::framing(gpu.center, gpu.radius);
    camera.distance *= 0.6; // close enough that atoms are tens of pixels wide
    let render = |outline: bool| {
        let mut renderer = Renderer::new(ctx.clone(), w, h);
        let bindings = renderer.bind(&gpu);
        let settings = RenderSettings {
            occlusion_culling: false,
            lighting: StylePreset::PublicationWhite.lighting(),
            material: StylePreset::PublicationWhite.material(),
            background: StylePreset::PublicationWhite.background(),
            fxaa: false,
            outline,
            ..Default::default()
        };
        let mut encoder = ctx.device.create_command_encoder(&Default::default());
        renderer.render(&mut encoder, &camera, &gpu, &bindings, &settings);
        ctx.queue.submit([encoder.finish()]);
        renderer.after_submit();
        // FXAA off, so the display is the outline pass output too.
        (renderer.read_color(), renderer.read_display_color())
    };
    let (plain, plain_display) = render(false);
    let (lined, _) = render(true);
    assert_eq!(plain, plain_display, "with outlines off the pass is a copy");

    let total = (w * h) as usize;
    let mut darkened = 0;
    let mut brightened = 0;
    let mut dark_on_bg = 0;
    for i in 0..total {
        let (a, b) = (luma(&plain, i), luma(&lined, i));
        if b < a - 1.0 {
            darkened += 1;
            if a > 250.0 {
                dark_on_bg += 1;
            }
        } else if b > a + 1.0 {
            brightened += 1;
        }
    }
    assert!(
        darkened > total / 100,
        "outlines should darken a visible fraction: {darkened}/{total}"
    );
    assert!(
        darkened < total / 3,
        "outlines should be lines, not a wash: {darkened}/{total}"
    );
    assert_eq!(brightened, 0, "outlines only ever darken");
    assert_eq!(
        dark_on_bg, 0,
        "lines sit on the atoms' rims, never on the background"
    );
}

/// The pixels the selection halo must cover, from a frame's pick ids:
/// unselected, with a selected pixel within `width` (the shader's
/// rounded disc, `dx^2 + dy^2 <= w^2 + w`).
fn selection_halo(ids: &[u32], w: i32, h: i32, width: i32, selected: &[bool]) -> Vec<bool> {
    let is_sel = |x: i32, y: i32| {
        let id = ids[(y * w + x) as usize];
        id != 0 && id & vv_render::BOND_ID_FLAG == 0 && selected[(id - 1) as usize]
    };
    let reach = width * width + width;
    (0..w * h)
        .map(|i| {
            let (x, y) = (i % w, i / w);
            !is_sel(x, y)
                && (-width..=width).any(|dy| {
                    (-width..=width).any(|dx| {
                        let (qx, qy) = (x + dx, y + dy);
                        dx * dx + dy * dy <= reach
                            && (0..w).contains(&qx)
                            && (0..h).contains(&qy)
                            && is_sel(qx, qy)
                    })
                })
        })
        .collect()
}

/// The selection halo is exactly the ring just outside the selected
/// atoms' visible id region, colours no other pixel, is absent
/// for an empty selection, and follows an orbit within the same frame.
#[test]
fn selection_outline_rings_exactly_the_selected_id_region() {
    let Some(ctx) = context() else { return };
    let (w, h) = (256u32, 256u32);
    let structure = protein_like(&SynthParams::new(3_000));
    let gpu = GpuStructure::upload(&ctx, &structure, None, ColorScheme::Element);
    let mut renderer = Renderer::new(ctx.clone(), w, h);
    let bindings = renderer.bind(&gpu);
    let mut camera = Camera::framing(gpu.center, gpu.radius);
    camera.distance *= 0.8;
    let settings = RenderSettings {
        occlusion_culling: false,
        fxaa: false,
        ..Default::default()
    };
    let green = [0u8, 255, 0];
    let frame = |renderer: &mut Renderer, camera: &Camera| {
        let mut encoder = ctx.device.create_command_encoder(&Default::default());
        renderer.render(&mut encoder, camera, &gpu, &bindings, &settings);
        ctx.queue.submit([encoder.finish()]);
        renderer.after_submit();
        (renderer.read_ids(), renderer.read_display_color())
    };
    let coords = structure.frame(0);
    let selected: Vec<bool> = coords
        .positions()
        .iter()
        .map(|p| p.x < gpu.center.x)
        .collect();

    // Consecutive frames of one view differ in a few pixels (ties in
    // depth), so each frame is checked against its own ids.
    let greens = |rgba: &[u8]| rgba.chunks(4).filter(|p| p[..3] == green).count();
    renderer.set_selection(vec![vv_render::ItemSelection::default()]);
    let (_, empty) = frame(&mut renderer, &camera);
    assert_eq!(greens(&empty), 0, "an empty selection draws nothing");

    renderer.set_selection(vec![vv_render::ItemSelection {
        atoms: vv_render::ItemSelection::bits(selected.len(), |i| selected[i]),
        bonds: Vec::new(),
    }]);
    let width = settings.selection_width as i32;
    let (ids, lined) = frame(&mut renderer, &camera);
    let halo = selection_halo(&ids, w as i32, h as i32, width, &selected);
    let ring = halo.iter().filter(|&&b| b).count();
    assert!(ring > 200, "a visible ring: {ring} px");
    for (i, &in_halo) in halo.iter().enumerate() {
        let is_green = lined[i * 4..i * 4 + 3] == green;
        assert_eq!(is_green, in_halo, "pixel {i}");
    }

    camera.orbit(0.6, 0.3);
    let (ids_orbit, orbited) = frame(&mut renderer, &camera);
    assert_ne!(ids_orbit, ids, "the orbit moved the atoms");
    let halo = selection_halo(&ids_orbit, w as i32, h as i32, width, &selected);
    for (i, &in_halo) in halo.iter().enumerate() {
        let is_green = orbited[i * 4..i * 4 + 3] == green;
        assert_eq!(is_green, in_halo, "after the orbit, pixel {i}");
    }
}

/// A backbone tube built from spline samples renders through the same
/// sphere + cylinder path as atoms, and picks resolve to a sample index
/// within range (the app maps that back to a residue's atom).
#[test]
fn backbone_tube_renders_and_picks() {
    let Some(ctx) = context() else { return };
    let path =
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/small/1CRN.cif");
    let structure = vv_io::load(path).unwrap();
    let tube = vv_core::backbone_tube(&structure.topology, structure.frame(0).positions());
    assert!(tube.positions.len() > 200, "{}", tube.positions.len());
    let colors = vec![vv_render::color::rgba(200, 80, 80); tube.positions.len()];
    let radii = vec![vv_core::TUBE_RADIUS; tube.positions.len()];

    let (w, h) = (256, 256);
    let mut renderer = Renderer::new(ctx.clone(), w, h);
    let gpu = GpuStructure::from_parts(&ctx, &tube.positions, &radii, &colors, &tube.bonds);
    assert_eq!(gpu.atom_count, tube.positions.len());
    assert_eq!(gpu.bond_count, tube.bonds.len());
    let bindings = renderer.bind(&gpu);
    let camera = Camera::framing(gpu.center, gpu.radius);
    let settings = RenderSettings {
        representation: Representation::Tube,
        occlusion_culling: false,
        ..Default::default()
    };
    let mut encoder = ctx.device.create_command_encoder(&Default::default());
    renderer.render(&mut encoder, &camera, &gpu, &bindings, &settings);
    ctx.queue.submit([encoder.finish()]);
    renderer.after_submit();

    let pixels = renderer.read_color();
    let bg = background_of(&pixels);
    let drawn = count_non_background(&pixels, bg);
    assert!(drawn > 2_000, "a tube should be clearly visible: {drawn}");
    assert!(drawn < (w * h / 2) as usize, "but thin: {drawn}");

    // Walk the middle row until something is hit; it must be a sample or
    // a segment of this tube.
    let hit = (0..w).find_map(|x| renderer.pick(x, h / 2));
    match hit {
        Some(Pick::Atom { item: 0, atom }) => assert!((atom as usize) < tube.positions.len()),
        Some(Pick::Bond { item: 0, bond }) => assert!((bond as usize) < tube.bonds.len()),
        other => panic!("expected to hit the tube on the middle row, got {other:?}"),
    }
}

/// A straight synthetic chain, half its residues at a low B-factor and
/// half at a high one: putty radii (`vv_core::backbone::putty_radius`,
/// what `vv_app`'s `radius_by` option feeds the tube) render the high-B
/// half's cross-section visibly wider than the low-B half's. No bonds:
/// `Representation::Tube`'s fixed 0.3 A bond radius would otherwise draw
/// through the low-B spheres and dominate the measurement.
#[test]
fn putty_tube_is_thicker_where_b_factor_is_high() {
    let Some(ctx) = context() else { return };
    let (low_b, high_b) = (5.0f32, 60.0f32);
    // Both ends stay above the renderer's ~1 px "draw as a point" cutoff
    // at this framing, so the comparison is of real sphere silhouettes.
    let (radius_min, radius_max) = (0.3f32, 2.0f32);
    let residues_per_half = 4u32;
    let step = 3.8f32;

    let mut builder = vv_core::TopologyBuilder::new();
    for i in 0..residues_per_half * 2 {
        let mut name = [b' '; 4];
        name[..2].copy_from_slice(b"CA");
        builder.push(&vv_core::AtomRow {
            element: vv_core::Element::CARBON,
            name,
            serial: i + 1,
            alt_loc: 0,
            comp: "ALA",
            asym: "A",
            auth_asym: "A",
            seq_id: i as i32 + 1,
            auth_seq_id: i as i32 + 1,
            ins_code: 0,
            entity: 1,
            position: glam::Vec3::new(i as f32 * step, 0.0, 0.0),
            occupancy: 1.0,
            b_factor: if i < residues_per_half { low_b } else { high_b },
            charge: 0,
            hetero: false,
        });
    }
    let structure = builder.finish().unwrap();
    let positions = structure.frame(0).positions().to_vec();
    let tube = vv_core::backbone_tube(&structure.topology, &positions);
    assert!(tube.positions.len() > 10, "{}", tube.positions.len());

    let b_factor: Vec<f32> = tube
        .source
        .iter()
        .map(|&a| structure.topology.b_factor[a as usize])
        .collect();
    let b_range = vv_render::scene::scalar_range(&b_factor);
    let radii: Vec<f32> = b_factor
        .iter()
        .map(|&v| vv_core::backbone::putty_radius(v, b_range, (radius_min, radius_max)))
        .collect();
    let colors = vec![vv_render::color::rgba(200, 80, 80); tube.positions.len()];

    let (w, h) = (400, 200);
    let mut renderer = Renderer::new(ctx.clone(), w, h);
    let gpu = GpuStructure::from_parts(&ctx, &tube.positions, &radii, &colors, &[]);
    let bindings = renderer.bind(&gpu);
    let camera = Camera::framing(gpu.center, gpu.radius);
    let settings = RenderSettings {
        representation: Representation::Tube,
        occlusion_culling: false,
        ..Default::default()
    };
    let mut encoder = ctx.device.create_command_encoder(&Default::default());
    renderer.render(&mut encoder, &camera, &gpu, &bindings, &settings);
    ctx.queue.submit([encoder.finish()]);
    renderer.after_submit();

    // Classify each column by whichever atom's sphere is visible at the
    // equator there (it protrudes furthest toward the camera at its own
    // center, so the pick is exact even where a low- and a high-B sphere
    // are screen-adjacent), and take the tallest column of each half's own
    // pixels -- not a fixed screen split, which a sample straddling both
    // halves would throw off.
    let pixels = renderer.read_color();
    let bg = background_of(&pixels);
    let column_height = |x: u32| -> usize {
        (0..h)
            .filter(|&y| {
                let i = ((y * w + x) * 4) as usize;
                (pixels[i] as i32 - bg[0] as i32).abs() > 8
                    || (pixels[i + 1] as i32 - bg[1] as i32).abs() > 8
                    || (pixels[i + 2] as i32 - bg[2] as i32).abs() > 8
            })
            .count()
    };
    let (mut low_max, mut high_max) = (0usize, 0usize);
    for x in 0..w {
        if let Some(Pick::Atom { item: 0, atom }) = renderer.pick(x, h / 2) {
            let height = column_height(x);
            if tube.source[atom as usize] < residues_per_half {
                low_max = low_max.max(height);
            } else {
                high_max = high_max.max(height);
            }
        }
    }
    eprintln!("putty: low={low_max}px high={high_max}px");
    assert!(low_max > 0, "the low-B half should be visible");
    assert!(
        high_max > low_max * 3,
        "the high-B half should render much thicker: low={low_max}px high={high_max}px"
    );
}

/// A straight chain, B-factor rising steadily residue to residue: a putty
/// tube built through `vv_core::cartoon::tube_plan` (the ribbon cartoon's
/// own mesh path, a round cross-section) has to widen smoothly along its
/// length. The old sphere-plus-constant-radius-cylinder tube necked down
/// to the thin end between every pair of samples; this checks the
/// rendered silhouette directly, column by column, so a regression back
/// to that shape (or any other periodic gap) shows up as a sudden drop.
#[test]
fn putty_tube_silhouette_width_varies_smoothly_along_its_length() {
    let Some(ctx) = context() else { return };
    let n = 20u32;
    let step = 3.8f32;
    let (radius_min, radius_max) = (0.2f32, 1.4f32);

    let mut builder = vv_core::TopologyBuilder::new();
    for i in 0..n {
        let mut name = [b' '; 4];
        name[..2].copy_from_slice(b"CA");
        builder.push(&vv_core::AtomRow {
            element: vv_core::Element::CARBON,
            name,
            serial: i + 1,
            alt_loc: 0,
            comp: "ALA",
            asym: "A",
            auth_asym: "A",
            seq_id: i as i32 + 1,
            auth_seq_id: i as i32 + 1,
            ins_code: 0,
            entity: 1,
            position: glam::Vec3::new(i as f32 * step, 0.0, 0.0),
            occupancy: 1.0,
            b_factor: i as f32,
            charge: 0,
            hetero: false,
        });
    }
    let structure = builder.finish().unwrap();
    let positions = structure.frame(0).positions().to_vec();
    let topology = &structure.topology;

    let trace = vv_core::backbone_trace(topology, &positions);
    let b_range = (0.0f32, (n - 1) as f32);
    let radius = |a: u32| {
        vv_core::backbone::putty_radius(
            topology.b_factor[a as usize],
            b_range,
            (radius_min, radius_max),
        )
    };
    let plan = vv_core::cartoon::tube_plan(&trace, radius);
    let spline = plan.frame(&positions);
    let colors = vec![vv_render::color::rgba(200, 80, 80); topology.atom_count()];
    let gpu = CartoonGpu::upload(&ctx, &plan, &spline, &colors).unwrap();

    let (w, h) = (1600u32, 300u32);
    let mut renderer = Renderer::new(ctx.clone(), w, h);
    let bindings = renderer.bind_cartoon(&gpu);
    let center = glam::Vec3::new((n - 1) as f32 * step * 0.5, 0.0, 0.0);
    let half_length = (n - 1) as f32 * step * 0.5 + radius_max + 1.0;
    // `Camera::framing`'s radius fits a sphere to the *vertical* half
    // angle; this rod is far longer than it is wide, so fit `half_length`
    // to the *horizontal* one instead (same formula, the other axis) --
    // otherwise most of the frame is empty margin above and below a rod
    // that only fills a handful of vertical pixels, and the diameter this
    // test measures gets squeezed into too few pixels to read.
    let fov_y = 45f32.to_radians();
    let half_h_angle = ((w as f32 / h as f32) * (fov_y * 0.5).tan()).atan();
    let radius_for_framing = half_length * (fov_y * 0.5).sin() / (1.05 * half_h_angle.tan());
    let camera = Camera::framing(center, radius_for_framing);
    let settings = RenderSettings {
        occlusion_culling: false,
        ..Default::default()
    };
    let mut encoder = ctx.device.create_command_encoder(&Default::default());
    renderer.render_all(
        &mut encoder,
        &camera,
        &[],
        &[CartoonItem {
            mesh: vv_render::CartoonMesh::Ribbon(&gpu),
            bindings: &bindings,
            material: vv_render::Material::default(),
        }],
        &[],
        &[],
        &settings,
    );
    ctx.queue.submit([encoder.finish()]);
    renderer.after_submit();

    let pixels = renderer.read_color();
    let bg = background_of(&pixels);
    let column_height = |x: u32| -> usize {
        (0..h)
            .filter(|&y| {
                let i = ((y * w + x) * 4) as usize;
                (pixels[i] as i32 - bg[0] as i32).abs() > 8
                    || (pixels[i + 1] as i32 - bg[1] as i32).abs() > 8
                    || (pixels[i + 2] as i32 - bg[2] as i32).abs() > 8
            })
            .count()
    };
    let heights: Vec<usize> = (0..w).map(column_height).collect();
    let drawn = heights.iter().filter(|&&h| h > 0).count();
    assert!(
        drawn as u32 > w / 2,
        "the tube should cover most of the width: {drawn}"
    );

    // A 5-tap moving average over columns: the octagonal ring
    // (`vv_core::cartoon::RING == 8`) makes an exact circle's silhouette
    // wobble a few percent between a vertex-up and an edge-up
    // cross-section, and single columns quantize to whole pixels: neither
    // is the "periodic neck" this test looks for, which spans several
    // columns (a full residue step, `SAMPLES_PER_RESIDUE` cross-sections
    // wide). Smoothing suppresses both without hiding a real one.
    let smoothed: Vec<f32> = (0..heights.len())
        .map(|x| {
            let (lo, hi) = (x.saturating_sub(2), (x + 2).min(heights.len() - 1));
            heights[lo..=hi].iter().sum::<usize>() as f32 / (hi - lo + 1) as f32
        })
        .collect();
    // The tube's own hard edges -- its flat, open ends (`tube_plan` builds
    // no caps) -- are real, large, one-time drops, not the periodic neck
    // this test looks for; excluding a margin around them (wider than the
    // smoothing kernel) keeps them out of the comparison below.
    let first = heights.iter().position(|&px| px > 0).unwrap();
    let last = heights.iter().rposition(|&px| px > 0).unwrap();
    let margin = 10;
    let interior = &smoothed[(first + margin)..=(last - margin)];
    let mut max_drop = 0.0f32;
    for pair in interior.windows(2) {
        let (a, b) = (pair[0], pair[1]);
        if a < 8.0 || b < 8.0 {
            continue;
        }
        max_drop = max_drop.max((a - b).abs() / a.max(b));
    }
    eprintln!("putty silhouette: max column-to-column relative change {max_drop}");
    assert!(
        max_drop < 0.15,
        "a periodic neck between spline samples: max relative change {max_drop}"
    );
}

/// A putty tube mesh still picks back to a real trace atom, through the
/// same per-section source map (`CartoonGpu::source`) `gpu_cache::
/// DrawSource::atom` resolves a viewport pick with.
#[test]
fn putty_tube_mesh_picks_a_real_trace_atom() {
    let Some(ctx) = context() else { return };
    let path =
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/small/1CRN.cif");
    let structure = vv_io::load(path).unwrap();
    let positions = structure.frame(0).positions().to_vec();
    let topology = &structure.topology;
    let trace = vv_core::backbone_trace(topology, &positions);
    let b_factor: Vec<f32> = trace
        .segments
        .iter()
        .flatten()
        .map(|&a| topology.b_factor[a as usize])
        .collect();
    let b_range = vv_render::scene::scalar_range(&b_factor);
    let (radius_min, radius_max) = (0.2f32, 1.5f32);
    let radius = |a: u32| {
        vv_core::backbone::putty_radius(
            topology.b_factor[a as usize],
            b_range,
            (radius_min, radius_max),
        )
    };
    let plan = vv_core::cartoon::tube_plan(&trace, radius);
    let spline = plan.frame(&positions);
    let colors = vv_render::colors_for(ColorScheme::BFactor, topology);
    let gpu = CartoonGpu::upload(&ctx, &plan, &spline, &colors).unwrap();

    let (w, h) = (256, 256);
    let mut renderer = Renderer::new(ctx.clone(), w, h);
    let bindings = renderer.bind_cartoon(&gpu);
    let (center, radius) = vv_core::CoordSet::new(positions.clone())
        .bounding_sphere()
        .unwrap();
    let camera = Camera::framing(center, radius);
    let settings = RenderSettings {
        occlusion_culling: false,
        ..Default::default()
    };
    let mut encoder = ctx.device.create_command_encoder(&Default::default());
    renderer.render_all(
        &mut encoder,
        &camera,
        &[],
        &[CartoonItem {
            mesh: vv_render::CartoonMesh::Ribbon(&gpu),
            bindings: &bindings,
            material: vv_render::Material::default(),
        }],
        &[],
        &[],
        &settings,
    );
    ctx.queue.submit([encoder.finish()]);
    renderer.after_submit();

    let hit = (0..w).find_map(|x| renderer.pick(x, h / 2));
    match hit {
        Some(Pick::Atom { item: 0, atom }) => {
            let real = gpu.source[atom as usize];
            assert!((real as usize) < topology.atom_count());
        }
        other => panic!("expected to hit the putty tube on the middle row, got {other:?}"),
    }
}

/// The path tracer draws a putty tube exactly as the raster does: both
/// build their triangles from the same `vv_core::cartoon::tube_plan`
/// mesh, one uploaded to `CartoonGpu`, the other pushed straight into
/// `TraceScene::push_mesh`.
#[test]
fn unshadowed_traced_putty_tube_matches_the_raster() {
    let Some(ctx) = context() else { return };
    let (w, h) = (640u32, 480u32);
    let path =
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/small/1UBQ.cif");
    let structure = vv_io::load(path).unwrap();
    let coords = structure.frame(0);
    let positions = coords.positions();
    let topology = &structure.topology;
    let trace = vv_core::backbone_trace(topology, positions);
    let b_factor: Vec<f32> = trace
        .segments
        .iter()
        .flatten()
        .map(|&a| topology.b_factor[a as usize])
        .collect();
    let b_range = vv_render::scene::scalar_range(&b_factor);
    let (radius_min, radius_max) = (0.2f32, 1.2f32);
    let radius = |a: u32| {
        vv_core::backbone::putty_radius(
            topology.b_factor[a as usize],
            b_range,
            (radius_min, radius_max),
        )
    };
    let plan = vv_core::cartoon::tube_plan(&trace, radius);
    let spline = plan.frame(positions);
    let colors = vv_render::colors_for(ColorScheme::BFactor, topology);
    let gpu = CartoonGpu::upload(&ctx, &plan, &spline, &colors).unwrap();
    let (center, radius_b) = coords.bounding_sphere().unwrap();
    let mut camera = Camera::framing(center, radius_b * 0.8);
    camera.orbit(0.4, 0.2);
    // No highlight: its high exponent turns the GPU's and CPU's last-bit
    // differences in the ring normals into visible shifts.
    let matte = vv_render::Material {
        specular: 0.0,
        ..Default::default()
    };
    let settings = RenderSettings {
        occlusion_culling: false,
        background: vv_render::wgpu::Color {
            r: 1.0,
            g: 0.0,
            b: 1.0,
            a: 1.0,
        },
        fxaa: false,
        outline: false,
        ao: 0.0,
        depth_cue: 0.0,
        ..Default::default()
    };
    let mut renderer = Renderer::new(ctx.clone(), w, h);
    let bindings = renderer.bind_cartoon(&gpu);
    let mut encoder = ctx.device.create_command_encoder(&Default::default());
    renderer.render_all(
        &mut encoder,
        &camera,
        &[],
        &[CartoonItem {
            mesh: vv_render::CartoonMesh::Ribbon(&gpu),
            bindings: &bindings,
            material: matte,
        }],
        &[],
        &[],
        &settings,
    );
    ctx.queue.submit([encoder.finish()]);
    renderer.after_submit();
    let raster = renderer.read_color();

    let mut scene = vv_render::path_trace::TraceScene::default();
    scene.push_mesh(&plan.mesh(&spline).expand(), &colors, matte);
    let tracer = vv_render::path_trace::PathTracer::new(&ctx, &scene).unwrap();
    let traced = tracer.render(
        &ctx,
        &camera,
        &vv_render::path_trace::TraceSettings {
            width: w,
            height: h,
            samples: 1,
            lighting: settings.lighting,
            background: [0.0; 3],
            background_top: [0.0; 3],
            transparent: true,
            clip: None,
            unshadowed: true,
            light_spread: vv_render::path_trace::DEFAULT_LIGHT_SPREAD,
            ao: 1.0,
            direct: 1.0,
        },
        |_| {},
    );
    let (mut interior, mut close) = (0usize, 0usize);
    for i in (0..(w * h) as usize).map(|k| k * 4) {
        let r = &raster[i..i + 3];
        let background = r[0] > 240 && r[1] < 16 && r[2] > 240;
        if background || traced[i + 3] != 255 {
            continue;
        }
        interior += 1;
        let diff = (0..3)
            .map(|c| (r[c] as i32 - traced[i + c] as i32).abs())
            .max()
            .unwrap();
        close += (diff <= 6) as usize;
    }
    eprintln!("putty tube: {close}/{interior} interior pixels within 6 levels");
    assert!(interior > 1000, "{interior}");
    assert!(close * 100 >= interior * 99, "{close}/{interior}");
}

/// `vv_core::dssp::assign` + `vv_core::cartoon::build` on a real structure,
/// through the actual GPU path (`CartoonGpu`, `Renderer::render_all`'s
/// `cartoons` list): real ribbon geometry is clearly visible, not just a
/// few stray pixels, and not so much that it looks like a filled blob
/// instead of a thin ribbon.
#[test]
fn cartoon_renders_visible_ribbon_geometry() {
    let Some(ctx) = context() else { return };
    let path =
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/small/1CRN.cif");
    let structure = vv_io::load(path).unwrap();
    let coords = structure.frame(0);
    let positions = coords.positions();
    let codes = vv_core::dssp::assign(&structure.topology, positions);
    assert!(
        codes.contains(&vv_core::DsspCode::AlphaHelix),
        "1CRN should have a real alpha helix to draw wide"
    );
    let plan = vv_core::cartoon::plan(&structure.topology, positions, &codes);
    let frame = plan.frame(positions);
    let mesh = plan.mesh(&frame);
    assert!(!mesh.sections.is_empty());

    let colors = vec![vv_render::color::rgba(200, 80, 80); structure.atom_count()];
    let (w, h) = (256, 256);
    let mut renderer = Renderer::new(ctx.clone(), w, h);
    let gpu = CartoonGpu::upload(&ctx, &plan, &frame, &colors).expect("upload cartoon");

    let (center, radius) = vv_core::CoordSet::new(positions.to_vec())
        .bounding_sphere()
        .unwrap();
    let camera = Camera::framing(center, radius);
    let settings = RenderSettings {
        occlusion_culling: false,
        ..Default::default()
    };
    let cartoon_bindings = renderer.bind_cartoon(&gpu);
    let mut encoder = ctx.device.create_command_encoder(&Default::default());
    renderer.render_all(
        &mut encoder,
        &camera,
        &[],
        &[CartoonItem {
            mesh: CartoonMesh::Ribbon(&gpu),
            bindings: &cartoon_bindings,
            material: vv_render::Material::default(),
        }],
        &[],
        &[],
        &settings,
    );
    ctx.queue.submit([encoder.finish()]);
    renderer.after_submit();

    let pixels = renderer.read_color();
    let bg = background_of(&pixels);
    let drawn = count_non_background(&pixels, bg);
    assert!(
        drawn > 2_000,
        "a cartoon should be clearly visible: {drawn}"
    );
    assert!(
        drawn < (w * h * 3 / 4) as usize,
        "but not a filled blob: {drawn}"
    );
}

/// The glycan rep (`vv_core::glycan`) draws pixels only around the
/// structure's glycans: on the real 6X3Z fixture (a receptor with a
/// handful of N-linked glycans among ~4700 residues total), every drawn
/// pixel projects near some glycan residue's ring centroid rather than
/// spreading across the much larger protein's silhouette — this test
/// draws no protein atoms at all, so any stray pixel is a real bug
/// (e.g. a residue wrongly detected as a glycan, or a mis-oriented
/// glyph reaching far past its own residue).
#[test]
fn glycan_rep_draws_pixels_only_around_glycans() {
    let Some(ctx) = context() else { return };
    let path =
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/glycan/6X3Z.pdb");
    let structure = vv_io::load(path).unwrap();
    let coords = structure.frame(0);
    let positions = coords.positions();
    let bonds = vv_core::bonds::perceive(&structure.topology, positions);
    let plan = vv_core::GlycanPlan::build(&structure.topology, &bonds);
    assert!(
        !plan.residues.is_empty(),
        "6X3Z should have detected glycans"
    );

    let mut frame = vv_core::GlycanFrame::default();
    plan.update_into(positions, &mut frame);
    let mut mesh = vv_core::PolytopeMesh::default();
    vv_core::build_glycan_mesh(&plan, &frame, 4.0, |_| true, &mut mesh);
    assert!(!mesh.positions.is_empty());

    let (center, radius) = vv_core::CoordSet::new(positions.to_vec())
        .bounding_sphere()
        .unwrap();
    let (w, h) = (512u32, 512u32);
    let aspect = w as f32 / h as f32;
    let camera = Camera::framing(center, radius);

    let mut renderer = Renderer::new(ctx.clone(), w, h);
    let gpu = vv_render::GlycanGpu::upload(&ctx, &mesh).expect("upload glycan mesh");
    let bindings = renderer.bind_glycan(&gpu);
    let settings = RenderSettings {
        occlusion_culling: false,
        ..Default::default()
    };
    let mut encoder = ctx.device.create_command_encoder(&Default::default());
    renderer.render_all(
        &mut encoder,
        &camera,
        &[], // no protein atoms: only the glycan mesh should draw
        &[CartoonItem {
            mesh: CartoonMesh::Glycan(&gpu),
            bindings: &bindings,
            material: vv_render::Material::default(),
        }],
        &[],
        &[],
        &settings,
    );
    ctx.queue.submit([encoder.finish()]);
    renderer.after_submit();

    let pixels = renderer.read_color();
    let bg = background_of(&pixels);
    let drawn = count_non_background(&pixels, bg);
    assert!(drawn > 500, "glycans should be clearly visible: {drawn} px");
    assert!(
        drawn < (w * h / 3) as usize,
        "too many pixels for a few glycans out of ~4700 residues: {drawn}"
    );

    let pixel_of = |p: glam::Vec3| -> (f32, f32) {
        let clip = camera.proj(aspect) * camera.view() * glam::Vec4::new(p.x, p.y, p.z, 1.0);
        let ndc = clip.truncate() / clip.w;
        (
            (ndc.x * 0.5 + 0.5) * w as f32,
            (0.5 - ndc.y * 0.5) * h as f32,
        )
    };
    let centers: Vec<(f32, f32)> = frame.centroid.iter().map(|&c| pixel_of(c)).collect();
    // Generous: the glyphs' own screen size at this camera distance, not
    // a tight per-shape bound.
    let tolerance = 60.0f32;
    let mut stray = 0usize;
    for y in 0..h {
        for x in 0..w {
            let i = ((y * w + x) * 4) as usize;
            if count_non_background(&pixels[i..i + 4], bg) == 0 {
                continue;
            }
            let near_a_glycan = centers.iter().any(|&(cx, cy)| {
                let (dx, dy) = (x as f32 - cx, y as f32 - cy);
                dx * dx + dy * dy < tolerance * tolerance
            });
            if !near_a_glycan {
                stray += 1;
            }
        }
    }
    assert!(
        stray < drawn / 20,
        "{stray} of {drawn} drawn pixels are not near any glycan residue centroid"
    );
}

/// Picking a cartoon must resolve to a real mesh vertex, not just "some
/// id" - and, alongside a second item sharing the frame's pick-id space
/// (a `DrawItem`, drawn first per `render_all`'s ordering), it must land
/// in the *cartoon's own* range, not spuriously inside the sphere's.
#[test]
fn cartoon_picking_resolves_to_a_real_vertex_past_another_items_ids() {
    let Some(ctx) = context() else { return };
    let path =
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/small/1CRN.cif");
    let structure = vv_io::load(path).unwrap();
    let coords = structure.frame(0);
    let positions = coords.positions();
    let codes = vv_core::dssp::assign(&structure.topology, positions);
    let plan = vv_core::cartoon::plan(&structure.topology, positions, &codes);
    let frame = plan.frame(positions);
    let mesh = plan.mesh(&frame);
    let colors = vec![vv_render::color::rgba(200, 80, 80); structure.atom_count()];

    let (center, radius) = vv_core::CoordSet::new(positions.to_vec())
        .bounding_sphere()
        .unwrap();
    let (w, h) = (256, 256);
    let mut renderer = Renderer::new(ctx.clone(), w, h);
    let gpu_mesh = CartoonGpu::upload(&ctx, &plan, &frame, &colors).expect("upload cartoon");
    let cartoon_bindings = renderer.bind_cartoon(&gpu_mesh);

    // A small, off-to-the-side sphere as `items[0]`, so the cartoon's
    // pick ids are forced to continue after a real, nonzero atom_base -
    // this is what would catch an off-by-one or a base that resets to 0
    // for `cartoons` instead of continuing `items`' id space.
    let side = GpuStructure::from_parts(
        &ctx,
        &[center + glam::Vec3::new(radius * 5.0, 0.0, 0.0)],
        &[0.5],
        &[vv_render::color::rgba(80, 80, 200)],
        &[],
    );
    let side_bindings = renderer.bind(&side);
    let side_item = DrawItem {
        structure: &side,
        bindings: &side_bindings,
        representation: Representation::Spacefill,
        sizes: AtomSizes::of(Representation::Spacefill),
        material: vv_render::Material::default(),
    };

    let camera = Camera::framing(center, radius);
    let settings = RenderSettings {
        occlusion_culling: false,
        ..Default::default()
    };
    let mut encoder = ctx.device.create_command_encoder(&Default::default());
    renderer.render_all(
        &mut encoder,
        &camera,
        std::slice::from_ref(&side_item),
        &[CartoonItem {
            mesh: CartoonMesh::Ribbon(&gpu_mesh),
            bindings: &cartoon_bindings,
            material: vv_render::Material::default(),
        }],
        &[],
        &[],
        &settings,
    );
    ctx.queue.submit([encoder.finish()]);
    renderer.after_submit();

    // Walk the middle row until something is hit (same approach as
    // `backbone_tube_renders_and_picks`, rather than relying on the exact
    // center pixel landing on geometry); the far-off side sphere (item 0)
    // is nowhere near this row, so any hit must be the cartoon (item 1).
    let hit = (0..w).find_map(|x| renderer.pick(x, h / 2));
    match hit {
        Some(Pick::Atom { item: 1, atom }) => {
            assert!(
                (atom as usize) < mesh.sections.len(),
                "section {atom} out of range ({} sections)",
                mesh.sections.len()
            );
        }
        other => panic!("expected a cartoon section pick at item 1, got {other:?}"),
    }
}

/// The cartoon pipeline shares the opaque pass's depth buffer with every
/// impostor pipeline (module docs on `Renderer::render_all`): a sphere
/// placed in front of a cartoon ribbon must hide it completely, the same
/// depth test every other representation already relies on to composite
/// correctly. This is the one property that could not be checked any
/// other way before a renderer existed to check it against (`vv_core::
/// cartoon`'s own docs on its winding-convention guess).
#[test]
fn cartoon_depth_tests_correctly_against_an_impostor_sphere() {
    let Some(ctx) = context() else { return };
    let path =
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/small/1CRN.cif");
    let structure = vv_io::load(path).unwrap();
    let coords = structure.frame(0);
    let positions = coords.positions();
    let codes = vv_core::dssp::assign(&structure.topology, positions);
    let plan = vv_core::cartoon::plan(&structure.topology, positions, &codes);
    let frame = plan.frame(positions);
    let colors = vec![vv_render::color::rgba(200, 80, 80); structure.atom_count()];

    let (center, radius) = vv_core::CoordSet::new(positions.to_vec())
        .bounding_sphere()
        .unwrap();
    let (w, h) = (256, 256);
    let camera = Camera::framing(center, radius);
    let settings = RenderSettings {
        occlusion_culling: false,
        // This is a depth-ordering check: the giant sphere and the
        // cartoon's own mesh bounds aren't exactly concentric, so the
        // depth cue (bounds-dependent) would tint the giant sphere's two
        // renders by very slightly different amounts otherwise.
        depth_cue: 0.0,
        ..Default::default()
    };

    // A single giant, same-centered sphere. `Camera::framing` puts the
    // camera at `radius / sin(fov/2) * 1.05` =~ `2.74 * radius`: twice
    // `radius` comfortably overflows the whole viewport (more angular
    // size than the original, same-centered structure did, at the same
    // distance) while staying under 2.74, so the camera is still outside
    // the sphere - three times was tried first and put the camera
    // *inside* it, which is not "in front of the ribbon" at all.
    let giant = GpuStructure::from_parts(
        &ctx,
        &[center],
        &[radius * 2.0],
        &[vv_render::color::rgba(80, 200, 80)],
        &[],
    );

    let sphere_alone = {
        let mut renderer = Renderer::new(ctx.clone(), w, h);
        let bindings = renderer.bind(&giant);
        let mut encoder = ctx.device.create_command_encoder(&Default::default());
        renderer.render(&mut encoder, &camera, &giant, &bindings, &settings);
        ctx.queue.submit([encoder.finish()]);
        renderer.after_submit();
        renderer.read_color()
    };

    let combined = {
        let mut renderer = Renderer::new(ctx.clone(), w, h);
        let gpu_mesh = CartoonGpu::upload(&ctx, &plan, &frame, &colors).expect("upload cartoon");
        let cartoon_bindings = renderer.bind_cartoon(&gpu_mesh);
        let bindings = renderer.bind(&giant);
        let item = DrawItem {
            structure: &giant,
            bindings: &bindings,
            representation: Representation::Spacefill,
            sizes: AtomSizes::of(Representation::Spacefill),
            material: vv_render::Material::default(),
        };
        let mut encoder = ctx.device.create_command_encoder(&Default::default());
        renderer.render_all(
            &mut encoder,
            &camera,
            std::slice::from_ref(&item),
            &[CartoonItem {
                mesh: CartoonMesh::Ribbon(&gpu_mesh),
                bindings: &cartoon_bindings,
                material: vv_render::Material::default(),
            }],
            &[],
            &[],
            &settings,
        );
        ctx.queue.submit([encoder.finish()]);
        renderer.after_submit();
        renderer.read_color()
    };

    assert_eq!(sphere_alone.len(), combined.len());
    let differing = sphere_alone
        .chunks(4)
        .zip(combined.chunks(4))
        .filter(|(a, b)| a[..3] != b[..3])
        .count();
    assert!(
        differing == 0,
        "{differing} pixels differ: the ribbon shows through a sphere that should fully occlude it"
    );
}

/// The GPU ray-marched Gaussian surface should look substantially the
/// same as `vv_cpu::gaussian_surface`'s already-visually-verified render
/// of the same real structure: a
/// smooth, connected blob covering a real fraction of the frame, not a
/// blank screen (a broken world/view transform) or a fully filled one
/// (a broken isosurface crossing test).
#[test]
fn gaussian_surface_renders_a_real_blobby_surface() {
    let Some(ctx) = context() else { return };
    let path =
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/small/1CRN.pdb");
    let structure = vv_io::load(path).unwrap();
    let positions = structure.frame(0).positions().to_vec();
    let radii: Vec<f32> = structure
        .topology
        .element
        .iter()
        .map(|e| e.vdw_radius())
        .collect();
    let colors = vv_render::colors_for(ColorScheme::Element, &structure.topology);

    let (w, h) = (256, 256);
    let mut renderer = Renderer::new(ctx.clone(), w, h);
    let gpu = GaussianSurfaceGpu::upload(&ctx, &positions, &radii, &colors, 2.0)
        .expect("upload gaussian surface");

    let (center, radius) = structure.frame(0).bounding_sphere().unwrap();
    let camera = Camera::framing(center, radius);
    let settings = RenderSettings {
        occlusion_culling: false,
        ..Default::default()
    };
    let bindings = renderer.bind_gaussian_surface(&gpu);
    let mut encoder = ctx.device.create_command_encoder(&Default::default());
    renderer.render_all(
        &mut encoder,
        &camera,
        &[],
        &[],
        &[GaussianSurfaceItem {
            gpu: &gpu,
            bindings: &bindings,
            material: vv_render::Material::default(),
        }],
        &[],
        &settings,
    );
    ctx.queue.submit([encoder.finish()]);
    renderer.after_submit();

    let pixels = renderer.read_color();
    let bg = background_of(&pixels);
    let drawn = count_non_background(&pixels, bg);
    assert!(
        drawn > (w * h / 10) as usize,
        "surface should cover a real chunk of the frame: {drawn}"
    );
    assert!(
        drawn < (w * h * 9 / 10) as usize,
        "but framing should leave visible background: {drawn}"
    );
}

/// Both full-screen ray-cast surfaces generate their per-pixel ray
/// differently in orthographic (`fs_gaussian_surface`/`fs_skin_surface`'s
/// own comments) -- a real render, not just a compiling shader, is what
/// actually exercises that branch, so this covers both in one test
/// rather than trusting `draw.wgsl`'s own orthographic coverage to imply
/// theirs.
#[test]
fn orthographic_renders_both_surfaces_of_a_real_structure() {
    let Some(ctx) = context() else { return };
    let path =
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/small/1CRN.pdb");
    let structure = vv_io::load(path).unwrap();
    let positions = structure.frame(0).positions().to_vec();
    let colors = vv_render::colors_for(ColorScheme::Element, &structure.topology);
    let (center, radius) = structure.frame(0).bounding_sphere().unwrap();
    let mut camera = Camera::framing(center, radius);
    camera.projection = Projection::Orthographic;
    let settings = RenderSettings {
        occlusion_culling: false,
        ..Default::default()
    };
    let (w, h) = (256, 256);

    let radii: Vec<f32> = structure
        .topology
        .element
        .iter()
        .map(|e| e.vdw_radius())
        .collect();
    let gs_gpu = GaussianSurfaceGpu::upload(&ctx, &positions, &radii, &colors, 2.0)
        .expect("upload gaussian surface");
    let mut renderer = Renderer::new(ctx.clone(), w, h);
    let gs_bindings = renderer.bind_gaussian_surface(&gs_gpu);
    let mut encoder = ctx.device.create_command_encoder(&Default::default());
    renderer.render_all(
        &mut encoder,
        &camera,
        &[],
        &[],
        &[GaussianSurfaceItem {
            gpu: &gs_gpu,
            bindings: &gs_bindings,
            material: vv_render::Material::default(),
        }],
        &[],
        &settings,
    );
    ctx.queue.submit([encoder.finish()]);
    renderer.after_submit();
    let pixels = renderer.read_color();
    let bg = background_of(&pixels);
    let drawn = count_non_background(&pixels, bg);
    assert!(
        drawn > (w * h / 10) as usize,
        "gaussian surface should cover a real chunk: {drawn}"
    );
    assert!(
        drawn < (w * h * 9 / 10) as usize,
        "framing should leave visible background: {drawn}"
    );

    let weights = skin_weights(&structure);
    let ss_gpu = SkinSurfaceGpu::upload(
        &ctx,
        &positions,
        &weights,
        &colors,
        vv_core::skin_surface::DEFAULT_SHRINK,
    );
    let mut renderer = Renderer::new(ctx.clone(), w, h);
    let ss_bindings = renderer.bind_skin_surface(&ss_gpu);
    let mut encoder = ctx.device.create_command_encoder(&Default::default());
    renderer.render_all(
        &mut encoder,
        &camera,
        &[],
        &[],
        &[],
        &[PatchSurfaceItem {
            surface: PatchSurface::Skin(&ss_gpu, &ss_bindings),
            material: vv_render::Material::default(),
        }],
        &settings,
    );
    ctx.queue.submit([encoder.finish()]);
    renderer.after_submit();
    let pixels = renderer.read_color();
    let bg = background_of(&pixels);
    let drawn = count_non_background(&pixels, bg);
    assert!(
        drawn > (w * h / 10) as usize,
        "skin surface should cover a real chunk: {drawn}"
    );
    assert!(
        drawn < (w * h * 9 / 10) as usize,
        "framing should leave visible background: {drawn}"
    );
}

/// The GPU skin surface's own version of the test above: same real
/// structure, same coverage bounds -- proof the WGSL port of `vv_cpu::
/// skin_surface`'s ray-cast + validation logic actually draws a real,
/// substantial surface on real hardware, not just that the shader
/// compiles.
#[test]
fn skin_surface_renders_a_real_structure() {
    let Some(ctx) = context() else { return };
    let path =
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/small/1CRN.pdb");
    let structure = vv_io::load(path).unwrap();
    let positions = structure.frame(0).positions().to_vec();
    let weights = skin_weights(&structure);
    let colors = vv_render::colors_for(ColorScheme::Element, &structure.topology);

    let (w, h) = (256, 256);
    let mut renderer = Renderer::new(ctx.clone(), w, h);
    let gpu = SkinSurfaceGpu::upload(
        &ctx,
        &positions,
        &weights,
        &colors,
        vv_core::skin_surface::DEFAULT_SHRINK,
    );
    let counts = gpu.patch_counts;
    assert!(
        counts[1] > counts[0] && counts[2] > counts[0],
        "a folded protein has more edge and triangle patches than atoms: {counts:?}"
    );

    let (center, radius) = structure.frame(0).bounding_sphere().unwrap();
    let camera = Camera::framing(center, radius);
    let settings = RenderSettings {
        occlusion_culling: false,
        ..Default::default()
    };
    let bindings = renderer.bind_skin_surface(&gpu);
    let mut encoder = ctx.device.create_command_encoder(&Default::default());
    renderer.render_all(
        &mut encoder,
        &camera,
        &[],
        &[],
        &[],
        &[PatchSurfaceItem {
            surface: PatchSurface::Skin(&gpu, &bindings),
            material: vv_render::Material::default(),
        }],
        &settings,
    );
    ctx.queue.submit([encoder.finish()]);
    renderer.after_submit();

    let pixels = renderer.read_color();
    let bg = background_of(&pixels);
    let drawn = count_non_background(&pixels, bg);
    assert!(
        drawn > (w * h / 10) as usize,
        "surface should cover a real chunk of the frame: {drawn}"
    );
    assert!(
        drawn < (w * h * 9 / 10) as usize,
        "but framing should leave visible background: {drawn}"
    );
}

/// Depth-compositing sanity check, the same shape as `cartoon_depth_
/// tests_correctly_against_an_impostor_sphere`: a real impostor sphere
/// placed between the camera and the surface must fully occlude it,
/// proving the ray-marched fragment shader's `frag_depth` output
/// actually participates in the shared depth buffer rather than always
/// winning or losing against ordinary impostors.
#[test]
fn gaussian_surface_depth_tests_correctly_against_an_impostor_sphere() {
    let Some(ctx) = context() else { return };
    let path =
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/small/1CRN.pdb");
    let structure = vv_io::load(path).unwrap();
    let positions = structure.frame(0).positions().to_vec();
    let radii: Vec<f32> = structure
        .topology
        .element
        .iter()
        .map(|e| e.vdw_radius())
        .collect();
    let colors = vv_render::colors_for(ColorScheme::Element, &structure.topology);
    let settings = RenderSettings {
        occlusion_culling: false,
        // The depth cue fades by the bounds of everything drawn, which
        // differ between the two renders compared below by design; this
        // test is about depth ordering only.
        depth_cue: 0.0,
        ..Default::default()
    };
    let (w, h) = (256, 256);

    // Built once and reused for both renders below, so the occluder is
    // sized off the *same* bounding sphere the surface itself ray-marches
    // against -- `structure.frame(0).bounding_sphere()` (point positions
    // only, no atom radii) undersizes it: `GaussianSurfaceGpu`'s own
    // bounds pad every atom by its own `cutoff_radius`, which is visibly
    // larger, and an occluder sized off the wrong (smaller) sphere can
    // leave the padded region peeking out at grazing angles -- caught
    // exactly this way, by the test failing, not assumed.
    let gpu = GaussianSurfaceGpu::upload(&ctx, &positions, &radii, &colors, 2.0)
        .expect("upload gaussian surface");
    let camera = Camera::framing(gpu.bounds_center, gpu.bounds_radius);

    let surface_alone = {
        let mut renderer = Renderer::new(ctx.clone(), w, h);
        let bindings = renderer.bind_gaussian_surface(&gpu);
        let mut encoder = ctx.device.create_command_encoder(&Default::default());
        renderer.render_all(
            &mut encoder,
            &camera,
            &[],
            &[],
            &[GaussianSurfaceItem {
                gpu: &gpu,
                bindings: &bindings,
                material: vv_render::Material::default(),
            }],
            &[],
            &settings,
        );
        ctx.queue.submit([encoder.finish()]);
        renderer.after_submit();
        renderer.read_color()
    };
    let bg = background_of(&surface_alone);
    let surface_only_hits = count_non_background(&surface_alone, bg);
    assert!(
        surface_only_hits > 0,
        "sanity: the surface alone must draw something"
    );

    // A single huge occluding sphere at the surface's own bounds center,
    // generously larger than the surface's own bounds radius (so it's a
    // strict superset, concentric, at every angle) but still comfortably
    // inside the camera's framing distance, so the camera is genuinely
    // outside it -- the same trap the cartoon version of this test
    // already documents and avoids.
    let occluder_radius = gpu.bounds_radius * 1.5;
    assert!(
        occluder_radius < (camera.eye() - gpu.bounds_center).length(),
        "camera must stay outside the occluder"
    );
    let occluder = GpuStructure::from_parts(
        &ctx,
        &[gpu.bounds_center],
        &[occluder_radius],
        &[vv_render::pack_rgba(255, 255, 255, 255)],
        &[],
    );
    // Compare pixel-for-pixel against the occluder *alone* rather than
    // asserting a literal "white" threshold: a shaded impostor sphere is
    // dimmed by diffuse falloff toward its silhouette and away from the
    // specular highlight, so it is never a uniform (255,255,255) even
    // with a pure-white base color -- the same reason the cartoon version
    // of this test compares two renders instead of an absolute color.
    let occluder_alone = {
        let mut renderer = Renderer::new(ctx.clone(), w, h);
        let occluder_bindings = renderer.bind(&occluder);
        let mut encoder = ctx.device.create_command_encoder(&Default::default());
        renderer.render_all(
            &mut encoder,
            &camera,
            &[DrawItem {
                structure: &occluder,
                bindings: &occluder_bindings,
                representation: Representation::Spacefill,
                sizes: AtomSizes::of(Representation::Spacefill),
                material: vv_render::Material::default(),
            }],
            &[],
            &[],
            &[],
            &settings,
        );
        ctx.queue.submit([encoder.finish()]);
        renderer.after_submit();
        renderer.read_color()
    };
    let combined = {
        let mut renderer = Renderer::new(ctx.clone(), w, h);
        let surface_bindings = renderer.bind_gaussian_surface(&gpu);
        let occluder_bindings = renderer.bind(&occluder);
        let mut encoder = ctx.device.create_command_encoder(&Default::default());
        renderer.render_all(
            &mut encoder,
            &camera,
            &[DrawItem {
                structure: &occluder,
                bindings: &occluder_bindings,
                representation: Representation::Spacefill,
                sizes: AtomSizes::of(Representation::Spacefill),
                material: vv_render::Material::default(),
            }],
            &[],
            &[GaussianSurfaceItem {
                gpu: &gpu,
                bindings: &surface_bindings,
                material: vv_render::Material::default(),
            }],
            &[],
            &settings,
        );
        ctx.queue.submit([encoder.finish()]);
        renderer.after_submit();
        renderer.read_color()
    };
    assert_eq!(occluder_alone.len(), combined.len());
    let differing = occluder_alone
        .chunks(4)
        .zip(combined.chunks(4))
        .filter(|(a, b)| a[..3] != b[..3])
        .count();
    assert!(
        differing == 0,
        "{differing} pixels differ: the surface shows through an occluder that should fully hide it"
    );
}

/// Not run by default (`cargo test -- --ignored`): renders a real PNG for
/// manual inspection (the same "look at it" discipline the CPU path's own
/// `render_1crn_to_png_for_manual_inspection` used) and times the GPU
/// render, for a direct before/after against `vv_cpu::gaussian_surface`.
#[test]
#[ignore]
fn render_1crn_gaussian_surface_gpu_to_png_for_manual_inspection() {
    let Some(ctx) = context() else { return };
    let path =
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/small/1CRN.pdb");
    let structure = vv_io::load(path).unwrap();
    let positions = structure.frame(0).positions().to_vec();
    let radii: Vec<f32> = structure
        .topology
        .element
        .iter()
        .map(|e| e.vdw_radius())
        .collect();
    let colors = vv_render::colors_for(ColorScheme::Chain, &structure.topology);

    let (w, h) = (800, 600);
    let mut renderer = Renderer::new(ctx.clone(), w, h);
    let gpu = GaussianSurfaceGpu::upload(&ctx, &positions, &radii, &colors, 2.0)
        .expect("upload gaussian surface");
    let camera = Camera::framing(gpu.bounds_center, gpu.bounds_radius);
    let settings = RenderSettings {
        occlusion_culling: false,
        ..Default::default()
    };
    let bindings = renderer.bind_gaussian_surface(&gpu);

    // Warm-up: shader/pipeline compilation and first-submit driver
    // overhead shouldn't count toward the steady-state frame time.
    for _ in 0..3 {
        let mut encoder = ctx.device.create_command_encoder(&Default::default());
        renderer.render_all(
            &mut encoder,
            &camera,
            &[],
            &[],
            &[GaussianSurfaceItem {
                gpu: &gpu,
                bindings: &bindings,
                material: vv_render::Material::default(),
            }],
            &[],
            &settings,
        );
        ctx.queue.submit([encoder.finish()]);
        renderer.after_submit();
    }
    ctx.device
        .poll(wgpu::PollType::wait_indefinitely())
        .unwrap();

    const FRAMES: u32 = 30;
    let start = std::time::Instant::now();
    for _ in 0..FRAMES {
        let mut encoder = ctx.device.create_command_encoder(&Default::default());
        renderer.render_all(
            &mut encoder,
            &camera,
            &[],
            &[],
            &[GaussianSurfaceItem {
                gpu: &gpu,
                bindings: &bindings,
                material: vv_render::Material::default(),
            }],
            &[],
            &settings,
        );
        ctx.queue.submit([encoder.finish()]);
        renderer.after_submit();
    }
    ctx.device
        .poll(wgpu::PollType::wait_indefinitely())
        .unwrap();
    let elapsed = start.elapsed();
    let per_frame_ms = elapsed.as_secs_f64() * 1000.0 / FRAMES as f64;
    println!(
        "GPU gaussian surface (static upload, re-render only): {FRAMES} frames of {w}x{h} 1CRN in {:.1} ms total, {:.2} ms/frame, {:.1} fps",
        elapsed.as_secs_f64() * 1000.0,
        per_frame_ms,
        1000.0 / per_frame_ms
    );

    // The real trajectory-playback cost: a moving frame means new
    // positions every frame, so `GaussianSurfaceGpu::upload` (fresh
    // buffers, a fresh CPU-side grid build) and `bind_gaussian_surface`
    // both re-run inside the loop here, not hoisted above it like the
    // static-camera measurement above.
    let start = std::time::Instant::now();
    for _ in 0..FRAMES {
        let gpu = GaussianSurfaceGpu::upload(&ctx, &positions, &radii, &colors, 2.0)
            .expect("upload gaussian surface");
        let bindings = renderer.bind_gaussian_surface(&gpu);
        let mut encoder = ctx.device.create_command_encoder(&Default::default());
        renderer.render_all(
            &mut encoder,
            &camera,
            &[],
            &[],
            &[GaussianSurfaceItem {
                gpu: &gpu,
                bindings: &bindings,
                material: vv_render::Material::default(),
            }],
            &[],
            &settings,
        );
        ctx.queue.submit([encoder.finish()]);
        renderer.after_submit();
    }
    ctx.device
        .poll(wgpu::PollType::wait_indefinitely())
        .unwrap();
    let elapsed = start.elapsed();
    let per_frame_ms = elapsed.as_secs_f64() * 1000.0 / FRAMES as f64;
    println!(
        "GPU gaussian surface (fresh upload + grid rebuild every frame, the real trajectory-playback cost): {FRAMES} frames in {:.1} ms total, {:.2} ms/frame, {:.1} fps",
        elapsed.as_secs_f64() * 1000.0,
        per_frame_ms,
        1000.0 / per_frame_ms
    );

    let pixels = renderer.read_color();
    let out = std::env::temp_dir().join("vizviz_gaussian_surface_gpu_1crn.png");
    let file = std::fs::File::create(&out).unwrap();
    let mut encoder = png::Encoder::new(std::io::BufWriter::new(file), w, h);
    encoder.set_color(png::ColorType::Rgba);
    encoder.set_depth(png::BitDepth::Eight);
    encoder
        .write_header()
        .unwrap()
        .write_image_data(&pixels)
        .unwrap();
    println!("wrote {}", out.display());
}

/// Not run by default (needs `--release`, timed): the scale that
/// originally broke the skin surface -- `SkinSurfaceGpu`'s own doc
/// records a user's real 17,365-atom structure taking ~14s for a
/// single frame at a normal camera and crashing outright once zoomed in
/// (almost certainly a driver TDR timeout from the old O(atoms) per-
/// pixel ownership scan). Reproduces that scale with a synthetic
/// structure (real geometry isn't needed to stress the same code path:
/// mixed-complex construction and per-patch billboard rendering scale
/// with atom count, not with whether the positions come from a real
/// PDB) at both a normal framing and a close zoom, and asserts a real
/// time bound rather than just "didn't crash" -- run with `cargo test
/// -p vv-render --release --test headless -- --ignored
/// skin_surface_handles_17k_atoms_without_the_old_timeout --nocapture`.
#[test]
#[ignore]
fn skin_surface_handles_17k_atoms_without_the_old_timeout() {
    let Some(ctx) = context() else { return };
    let structure = protein_like(&SynthParams::new(17_365));
    let positions = structure.frame(0).positions().to_vec();
    let weights: Vec<f32> = structure
        .topology
        .element
        .iter()
        .map(|e| {
            vv_core::skin_surface::weight_for_radius(
                e.vdw_radius(),
                vv_core::skin_surface::DEFAULT_SHRINK,
            )
        })
        .collect();
    let colors = vv_render::colors_for(ColorScheme::Element, &structure.topology);

    let build_start = std::time::Instant::now();
    let gpu = SkinSurfaceGpu::upload(
        &ctx,
        &positions,
        &weights,
        &colors,
        vv_core::skin_surface::DEFAULT_SHRINK,
    );
    let build_time = build_start.elapsed();
    println!(
        "17,365 atoms: mixed complex + upload in {build_time:?}, {} patches {:?}",
        gpu.patch_count, gpu.patch_counts
    );

    let (center, radius) = structure.frame(0).bounding_sphere().unwrap();
    let (w, h) = (800u32, 600u32);
    let settings = RenderSettings {
        occlusion_culling: false,
        ..Default::default()
    };

    for (label, camera) in [
        ("normal framing", Camera::framing(center, radius)),
        ("close zoom", {
            let mut c = Camera::framing(center, radius);
            c.zoom(0.35);
            c
        }),
    ] {
        let mut renderer = Renderer::new(ctx.clone(), w, h);
        let bindings = renderer.bind_skin_surface(&gpu);
        let render_start = std::time::Instant::now();
        let mut encoder = ctx.device.create_command_encoder(&Default::default());
        renderer.render_all(
            &mut encoder,
            &camera,
            &[],
            &[],
            &[],
            &[PatchSurfaceItem {
                surface: PatchSurface::Skin(&gpu, &bindings),
                material: vv_render::Material::default(),
            }],
            &settings,
        );
        ctx.queue.submit([encoder.finish()]);
        renderer.after_submit();
        let pixels = renderer.read_color();
        let frame_time = render_start.elapsed();
        println!("  {label}: one frame in {frame_time:?}");
        assert!(
            frame_time < std::time::Duration::from_secs(2),
            "{label}: frame took {frame_time:?}, expected well under the old ~14s/crash -- \
             the per-patch billboard + competitor-list rewrite should keep this fast regardless \
             of camera position, unlike the old full-screen O(atoms)-per-pixel scan"
        );
        let bg = background_of(&pixels);
        let drawn = count_non_background(&pixels, bg);
        assert!(drawn > 0, "{label}: expected some surface to be visible");
    }
}

/// Not run by default (a timing measurement, not a correctness check):
/// 4HHB as transparent spacefill plus an opaque cartoon, the scene named
/// in the transparency task's perf ask, wall-clock timed (`Instant` +
/// `device.poll(wait)`, since the GPU timer's 8 timestamp slots are
/// already spoken for and don't cover the new glass passes) against the
/// same scene fully opaque, in `--release` only -- debug's unoptimized
/// CPU-side encoding dwarfs any GPU difference. Run with `cargo test -p
/// vv-render --release --test headless -- --ignored
/// frame_time_4hhb_transparent_spacefill_and_opaque_cartoon --nocapture`.
#[test]
#[ignore]
fn frame_time_4hhb_transparent_spacefill_and_opaque_cartoon() {
    let Some(ctx) = context() else { return };
    let path =
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/small/4HHB.cif");
    let structure = vv_io::load(path).unwrap();
    let coords = structure.frame(0);
    let positions = coords.positions();
    let codes = vv_core::dssp::assign(&structure.topology, positions);
    let plan = vv_core::cartoon::plan(&structure.topology, positions, &codes);
    let frame = plan.frame(positions);
    let cartoon_colors = vv_render::colors_for(ColorScheme::Element, &structure.topology);
    let gpu = GpuStructure::upload(&ctx, &structure, None, ColorScheme::Element);
    let (center, radius) = coords.bounding_sphere().unwrap();
    let (w, h) = (1400, 900);
    let mut renderer = Renderer::new(ctx.clone(), w, h);
    let bindings = renderer.bind(&gpu);
    let cartoon_gpu = CartoonGpu::upload(&ctx, &plan, &frame, &cartoon_colors).expect("cartoon");
    let cartoon_bindings = renderer.bind_cartoon(&cartoon_gpu);
    let camera = Camera::framing(center, radius);
    let settings = RenderSettings::default();

    let mut time = |material: vv_render::Material, frames: u32| {
        let item = DrawItem {
            structure: &gpu,
            bindings: &bindings,
            representation: Representation::Spacefill,
            sizes: AtomSizes::of(Representation::Spacefill),
            material,
        };
        let cartoons = [CartoonItem {
            mesh: CartoonMesh::Ribbon(&cartoon_gpu),
            bindings: &cartoon_bindings,
            material: vv_render::Material::default(),
        }];
        // Warm up pipeline/shader caches before timing.
        for _ in 0..3 {
            let mut encoder = ctx.device.create_command_encoder(&Default::default());
            renderer.render_all(
                &mut encoder,
                &camera,
                std::slice::from_ref(&item),
                &cartoons,
                &[],
                &[],
                &settings,
            );
            ctx.queue.submit([encoder.finish()]);
            renderer.after_submit();
            let _ = ctx.device.poll(wgpu::PollType::wait_indefinitely());
        }
        let start = std::time::Instant::now();
        for _ in 0..frames {
            let mut encoder = ctx.device.create_command_encoder(&Default::default());
            renderer.render_all(
                &mut encoder,
                &camera,
                std::slice::from_ref(&item),
                &cartoons,
                &[],
                &[],
                &settings,
            );
            ctx.queue.submit([encoder.finish()]);
            renderer.after_submit();
            let _ = ctx.device.poll(wgpu::PollType::wait_indefinitely());
        }
        start.elapsed().as_secs_f64() * 1000.0 / frames as f64
    };
    let frames = 120;
    let opaque_ms = time(vv_render::Material::default(), frames);
    let transparent_ms = time(
        vv_render::Material {
            opacity: 0.3,
            ..vv_render::Material::default()
        },
        frames,
    );
    println!(
        "{}: 4HHB ({} atoms) spacefill + opaque cartoon, {w}x{h}, {frames} frames: \
         opaque spacefill {opaque_ms:.3} ms/frame, transparent spacefill {transparent_ms:.3} ms/frame",
        ctx.adapter_name(),
        structure.atom_count()
    );
}

/// Not run by default: renders the GPU skin surface to a real PNG for a
/// human to actually look at, the same discipline `vv_cpu::
/// skin_surface`'s own manual-inspection test used -- which is what
/// caught every real visual bug this surface has had so far. A
/// correctness/visual check, not a benchmark (see the Gaussian surface
/// test above for the fps-measurement shape, not duplicated here).
#[test]
#[ignore]
fn render_1crn_skin_surface_gpu_to_png_for_manual_inspection() {
    let Some(ctx) = context() else { return };
    let path =
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/small/1CRN.pdb");
    let structure = vv_io::load(path).unwrap();
    let positions = structure.frame(0).positions().to_vec();
    let weights = skin_weights(&structure);
    let colors = vv_render::colors_for(ColorScheme::Chain, &structure.topology);

    let (w, h) = (800, 600);
    let mut renderer = Renderer::new(ctx.clone(), w, h);
    let gpu = SkinSurfaceGpu::upload(
        &ctx,
        &positions,
        &weights,
        &colors,
        vv_core::skin_surface::DEFAULT_SHRINK,
    );
    println!(
        "patches by kind (vertex, edge, triangle, tet): {:?}",
        gpu.patch_counts
    );
    let camera = Camera::framing(gpu.bounds_center, gpu.bounds_radius);
    let settings = RenderSettings {
        occlusion_culling: false,
        ..Default::default()
    };
    let bindings = renderer.bind_skin_surface(&gpu);

    let mut encoder = ctx.device.create_command_encoder(&Default::default());
    renderer.render_all(
        &mut encoder,
        &camera,
        &[],
        &[],
        &[],
        &[PatchSurfaceItem {
            surface: PatchSurface::Skin(&gpu, &bindings),
            material: vv_render::Material::default(),
        }],
        &settings,
    );
    ctx.queue.submit([encoder.finish()]);
    renderer.after_submit();

    let pixels = renderer.read_color();
    let out = std::env::temp_dir().join("vizviz_skin_surface_gpu_1crn.png");
    let file = std::fs::File::create(&out).unwrap();
    let mut encoder = png::Encoder::new(std::io::BufWriter::new(file), w, h);
    encoder.set_color(png::ColorType::Rgba);
    encoder.set_depth(png::BitDepth::Eight);
    encoder
        .write_header()
        .unwrap()
        .write_image_data(&pixels)
        .unwrap();
    println!("wrote {}", out.display());
}

/// Orthographic cylinders must draw their whole length even when the
/// bond tilts toward the viewer: the ray used to start at the billboard
/// (level with the bond's middle), which cut off the near half (pencil-
/// tip sticks, notched tubes).
#[test]
fn orthographic_bond_tilted_toward_the_viewer_draws_its_near_end() {
    let Some(ctx) = context() else { return };
    let (w, h) = (256, 256);
    let a = glam::Vec3::new(0.0, -3.0, -4.0);
    let b = glam::Vec3::new(0.0, 3.0, 4.0);
    let gpu = GpuStructure::from_parts(
        &ctx,
        &[a, b],
        &[0.3, 0.3],
        &[vv_render::pack_rgba(255, 255, 255, 255); 2],
        &[[0, 1]],
    );
    let mut camera = Camera::framing(glam::Vec3::ZERO, 6.0);
    camera.projection = Projection::Orthographic;
    let mut renderer = Renderer::new(ctx.clone(), w, h);
    let bindings = renderer.bind(&gpu);
    let settings = RenderSettings {
        representation: Representation::Tube,
        occlusion_culling: false,
        ao: 0.0,
        depth_cue: 0.0,
        ..Default::default()
    };
    let mut encoder = ctx.device.create_command_encoder(&Default::default());
    renderer.render(&mut encoder, &camera, &gpu, &bindings, &settings);
    ctx.queue.submit([encoder.finish()]);
    renderer.after_submit();
    let pixels = renderer.read_color();
    let bg = background_of(&pixels);
    // Points 25% and 75% of the way from a to b (the latter nearer the
    // eye, which looks down -z from +z), projected to pixels.
    let drawn_at = |t: f32| {
        let p = a.lerp(b, t);
        let clip = camera.proj(1.0) * camera.view() * p.extend(1.0);
        let ndc = clip.truncate().truncate() / clip.w;
        let x = ((ndc.x * 0.5 + 0.5) * w as f32) as usize;
        let y = ((0.5 - ndc.y * 0.5) * h as f32) as usize;
        let i = (y * w as usize + x) * 4;
        pixels[i..i + 3] != bg[..]
    };
    assert!(drawn_at(0.25), "far half of the bond is missing");
    assert!(drawn_at(0.75), "near half of the bond is missing");
}

/// Orthographic skin surface looking straight down an axis must match
/// the same view turned by 0.001 rad (regression: exactly axis-aligned
/// rays hit a sign(0) bug in the root solver, drawing rings and holes).
#[test]
fn orthographic_skin_surface_survives_axis_aligned_rays() {
    let Some(ctx) = context() else { return };
    let path =
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/small/1CRN.pdb");
    let structure = vv_io::load(path).unwrap();
    let positions = structure.frame(0).positions().to_vec();
    let colors = vv_render::colors_for(ColorScheme::Element, &structure.topology);
    let gpu = SkinSurfaceGpu::upload(
        &ctx,
        &positions,
        &skin_weights(&structure),
        &colors,
        vv_core::skin_surface::DEFAULT_SHRINK,
    );
    let (w, h) = (400, 256);
    let render = |camera: &Camera| {
        let mut renderer = Renderer::new(ctx.clone(), w, h);
        let bindings = renderer.bind_skin_surface(&gpu);
        let settings = RenderSettings {
            occlusion_culling: false,
            ao: 0.0,
            depth_cue: 0.0,
            background: PINHOLE_BG,
            ..Default::default()
        };
        let mut encoder = ctx.device.create_command_encoder(&Default::default());
        renderer.render_all(
            &mut encoder,
            camera,
            &[],
            &[],
            &[],
            &[PatchSurfaceItem {
                surface: PatchSurface::Skin(&gpu, &bindings),
                material: vv_render::Material::default(),
            }],
            &settings,
        );
        ctx.queue.submit([encoder.finish()]);
        renderer.after_submit();
        renderer.read_color()
    };
    let (center, radius) = structure.frame(0).bounding_sphere().unwrap();
    let mut straight = Camera::framing(center, radius);
    straight.projection = Projection::Orthographic;
    straight.zoom(0.4);
    let mut turned = straight.clone();
    turned.orbit(1e-3, 1e-3);
    let (a, b) = (render(&straight), render(&turned));
    let bg = background_of(&a);
    let drawn = |px: &[u8]| px[..3] != bg[..];
    let total = a.chunks(4).filter(|p| drawn(p)).count();
    let differ = a
        .chunks(4)
        .zip(b.chunks(4))
        .filter(|(pa, pb)| {
            // Missing, or a far-side hit shaded very differently.
            let shade: i32 = pa[..3]
                .iter()
                .zip(&pb[..3])
                .map(|(x, y)| (*x as i32 - *y as i32).abs())
                .sum();
            drawn(pa) != drawn(pb) || shade > 90
        })
        .count();
    assert!(total > 5000, "surface too small: {total}");
    // About 1% sit on silhouettes and creases the turn moves by a
    // sub-pixel; the sign(0) bug differed on over 20%.
    assert!(
        differ * 100 <= total * 3,
        "{differ} of {total} pixels differ"
    );
}

/// Magenta: never a lit surface colour, unlike a corner pixel when the
/// surface fills the frame.
const PINHOLE_BG: wgpu::Color = wgpu::Color {
    r: 1.0,
    g: 0.0,
    b: 1.0,
    a: 1.0,
};

/// Isolated background pixels whose eight neighbours are all surface:
/// seams where neither neighbouring patch claims the hit.
fn pinholes(pixels: &[u8], w: usize, h: usize) -> usize {
    let bg = [255u8, 0, 255];
    let drawn = |x: usize, y: usize| pixels[(y * w + x) * 4..(y * w + x) * 4 + 3] != bg[..];
    let mut holes = 0;
    for y in 1..h - 1 {
        for x in 1..w - 1 {
            if !drawn(x, y)
                && (0..9)
                    .filter(|k| *k != 4)
                    .all(|k| drawn(x + k % 3 - 1, y + k / 3 - 1))
            {
                holes += 1;
            }
        }
    }
    holes
}

/// The skin surface must be closed: no pinholes at patch seams, at a few
/// zooms and angles, in both projections, on 4HHB.
#[test]
fn skin_surface_has_no_pinholes() {
    let Some(ctx) = context() else { return };
    let path =
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/small/4HHB.pdb");
    let structure = vv_io::load(path).unwrap();
    let positions = structure.frame(0).positions().to_vec();
    let colors = vv_render::colors_for(ColorScheme::Element, &structure.topology);
    let gpu = SkinSurfaceGpu::upload(
        &ctx,
        &positions,
        &skin_weights(&structure),
        &colors,
        vv_core::skin_surface::DEFAULT_SHRINK,
    );
    let (w, h) = (512usize, 512usize);
    let (center, radius) = structure.frame(0).bounding_sphere().unwrap();
    let mut total = 0;
    for projection in [Projection::Orthographic, Projection::Perspective] {
        for (zoom, yaw) in [(1.0, 0.0), (0.5, 0.6), (0.25, 1.9)] {
            let mut camera = Camera::framing(center, radius);
            camera.projection = projection;
            camera.orbit(yaw, 0.2);
            camera.zoom(zoom);
            let mut renderer = Renderer::new(ctx.clone(), w as u32, h as u32);
            let bindings = renderer.bind_skin_surface(&gpu);
            let settings = RenderSettings {
                occlusion_culling: false,
                background: PINHOLE_BG,
                ao: 0.0,
                depth_cue: 0.0,
                ..Default::default()
            };
            let mut encoder = ctx.device.create_command_encoder(&Default::default());
            renderer.render_all(
                &mut encoder,
                &camera,
                &[],
                &[],
                &[],
                &[PatchSurfaceItem {
                    surface: PatchSurface::Skin(&gpu, &bindings),
                    material: vv_render::Material::default(),
                }],
                &settings,
            );
            ctx.queue.submit([encoder.finish()]);
            renderer.after_submit();
            let pixels = renderer.read_color();
            let holes = pinholes(&pixels, w, h);
            println!("{projection:?} zoom {zoom} yaw {yaw}: {holes} pinholes");
            total += holes;
        }
    }
    assert!(total <= 3, "{total} pinholes");
}

/// A gradient background runs from `background` at the bottom row to
/// `background_top` at the top, and the drawn structure is untouched.
#[test]
fn gradient_background_runs_bottom_to_top() {
    let Some(ctx) = context() else { return };
    let (w, h) = (128, 128);
    let structure = protein_like(&SynthParams::new(500));
    let gpu = GpuStructure::upload(&ctx, &structure, None, ColorScheme::Element);
    let mut camera = Camera::framing(gpu.center, gpu.radius);
    camera.zoom(3.0);
    let black = wgpu::Color::BLACK;
    let white = wgpu::Color::WHITE;
    let mut renderer = Renderer::new(ctx.clone(), w, h);
    let bindings = renderer.bind(&gpu);
    let settings = RenderSettings {
        occlusion_culling: false,
        background: black,
        background_top: Some(white),
        ao: 0.0,
        depth_cue: 0.0,
        ..Default::default()
    };
    let mut encoder = ctx.device.create_command_encoder(&Default::default());
    renderer.render(&mut encoder, &camera, &gpu, &bindings, &settings);
    ctx.queue.submit([encoder.finish()]);
    renderer.after_submit();
    let pixels = renderer.read_color();
    let red = |x: u32, y: u32| pixels[((y * w + x) * 4) as usize];
    assert!(red(0, 0) > 245, "top-left {}", red(0, 0));
    assert!(red(0, h - 1) < 20, "bottom-left {}", red(0, h - 1));
    let mid = red(0, h / 2);
    assert!((100..220).contains(&mid), "middle {mid}");
}

/// Shadows only ever remove light, and on a packed structure under an
/// oblique key light they remove some.
#[test]
fn shadows_darken_and_never_brighten() {
    let Some(ctx) = context() else { return };
    let (w, h) = (256, 256);
    let structure = protein_like(&SynthParams::new(4_000));
    let gpu = GpuStructure::upload(&ctx, &structure, None, ColorScheme::Element);
    let mut camera = Camera::framing(gpu.center, gpu.radius);
    camera.distance *= 0.4;
    let render = |shadows: f32| {
        let mut renderer = Renderer::new(ctx.clone(), w, h);
        let bindings = renderer.bind(&gpu);
        let settings = RenderSettings {
            occlusion_culling: false,
            lighting: vv_render::LightingPreset::Full.lighting(),
            ao: 0.0,
            depth_cue: 0.0,
            shadows,
            ..Default::default()
        };
        let mut encoder = ctx.device.create_command_encoder(&Default::default());
        renderer.render(&mut encoder, &camera, &gpu, &bindings, &settings);
        ctx.queue.submit([encoder.finish()]);
        renderer.after_submit();
        renderer.read_color()
    };
    let (lit, shadowed) = (render(0.0), render(1.0));
    let brighter = lit
        .chunks(4)
        .zip(shadowed.chunks(4))
        .filter(|(a, b)| (0..3).any(|c| b[c] > a[c].saturating_add(1)))
        .count();
    let darker = lit
        .chunks(4)
        .zip(shadowed.chunks(4))
        .filter(|(a, b)| (0..3).any(|c| b[c] + 8 < a[c]))
        .count();
    assert_eq!(brighter, 0, "shadows brightened {brighter} pixels");
    assert!(
        darker * 50 > (w * h) as usize,
        "only {darker} pixels darkened"
    );
}

/// A glass skin surface in front of opaque atoms: where it covers them
/// the pixel mixes the two, a click there still picks the atom behind,
/// and the frame without glass is unchanged by the glass pass existing.
#[test]
fn glass_shows_and_picks_what_is_behind_it() {
    let Some(ctx) = context() else { return };
    let (w, h) = (160, 120);
    let structure = protein_like(&SynthParams::new(400));
    let gpu = GpuStructure::upload(&ctx, &structure, None, ColorScheme::Element);
    // The same atoms 30 A toward the camera (it looks down -z), as glass.
    let positions: Vec<_> = structure
        .frame(0)
        .positions()
        .iter()
        .map(|p| *p + glam::Vec3::Z * 30.0)
        .collect();
    let colors = vv_render::colors_for(ColorScheme::Element, &structure.topology);
    let skin = SkinSurfaceGpu::upload(
        &ctx,
        &positions,
        &skin_weights(&structure),
        &colors,
        vv_core::skin_surface::DEFAULT_SHRINK,
    );
    let mut camera = Camera::framing(gpu.center, gpu.radius * 1.5);
    camera.projection = Projection::Orthographic;
    let mut renderer = Renderer::new(ctx.clone(), w, h);
    let bindings = renderer.bind(&gpu);
    let skin_bindings = renderer.bind_skin_surface(&skin);
    let settings = RenderSettings {
        occlusion_culling: false,
        ao: 0.0,
        depth_cue: 0.0,
        fxaa: false,
        ..Default::default()
    };
    let item = DrawItem {
        structure: &gpu,
        bindings: &bindings,
        representation: Representation::Spacefill,
        sizes: AtomSizes::of(Representation::Spacefill),
        material: vv_render::Material::default(),
    };
    let frame = |renderer: &mut Renderer, glass: Option<vv_render::MaterialPreset>| {
        let skins: Vec<_> = glass
            .map(|m| PatchSurfaceItem {
                surface: PatchSurface::Skin(&skin, &skin_bindings),
                material: m.into(),
            })
            .into_iter()
            .collect();
        let mut encoder = ctx.device.create_command_encoder(&Default::default());
        renderer.render_all(
            &mut encoder,
            &camera,
            std::slice::from_ref(&item),
            &[],
            &[],
            &skins,
            &settings,
        );
        ctx.queue.submit([encoder.finish()]);
        renderer.after_submit();
        renderer.read_color()
    };
    let bare = frame(&mut renderer, None);
    let bare_again = frame(&mut renderer, None);
    assert_eq!(bare, bare_again, "frames without glass must not change");
    let glass = frame(&mut renderer, Some(vv_render::MaterialPreset::ClearGlass));

    let (x, y) = (w / 2, h / 2);
    let i = ((y * w + x) * 4) as usize;
    let bg = background_of(&bare);
    assert_ne!(bare[i..i + 3], bg[..], "an atom should be at the centre");
    // A patch, not the one centre pixel: `ClearGlass`'s tint is a subtle
    // blend, close enough to `bare` that 8-bit rounding can collapse the
    // two at any single pixel.
    let patch_differs = (y.saturating_sub(4)..y + 4).any(|py| {
        (x.saturating_sub(4)..x + 4).any(|px| {
            let k = ((py * w + px) * 4) as usize;
            glass[k..k + 3] != bare[k..k + 3]
        })
    });
    assert!(patch_differs, "glass should tint the atom behind");
    assert!(
        matches!(renderer.pick(x, y), Some(Pick::Atom { item: 0, .. })),
        "a click on glass picks the atom behind it, got {:?}",
        renderer.pick(x, y)
    );
    let opaque = frame(&mut renderer, Some(vv_render::MaterialPreset::Opaque));
    assert_ne!(
        glass[i..i + 3],
        opaque[i..i + 3],
        "glass is not the opaque surface"
    );
}

/// A transparent sphere's specular highlight must be scaled by its own
/// opacity like the diffuse term, not added at full strength (`glass()`
/// in shading.wgsl): a fully bright highlight on a nearly invisible
/// surface reads as a glossy sheen rather than glass.
/// Diffuse and ambient are zero, so any brightening is specular alone. A
/// single-atom skin surface gives an exact analytic sphere (a ray-cast
/// vertex patch), so the front pole's view-space normal is exactly
/// `[0, 0, 1]` with no marched-surface gradient error.
#[test]
fn glass_specular_is_scaled_by_opacity() {
    let Some(ctx) = context() else { return };
    let (w, h) = (128, 128);
    let radius = 4.0;
    let weight =
        vv_core::skin_surface::weight_for_radius(radius, vv_core::skin_surface::DEFAULT_SHRINK);
    let skin = SkinSurfaceGpu::upload(
        &ctx,
        &[glam::Vec3::ZERO],
        &[weight],
        &[0xffff_ffffu32],
        vv_core::skin_surface::DEFAULT_SHRINK,
    );
    let mut camera = Camera::framing(glam::Vec3::ZERO, radius * 1.5);
    camera.projection = Projection::Orthographic;
    // Key light straight at the camera: at the point of the sphere facing
    // the camera (screen centre) the view-space normal is always
    // [0, 0, 1], opposite `view_dir`, so diffuse and specular both peak
    // there regardless of the camera's world-space placement.
    let lighting = vv_render::Lighting {
        lights: [
            vv_render::style::Light {
                dir: [0.0, 0.0, 1.0],
                intensity: 1.0,
                color: vv_render::style::Light::WHITE,
                _pad: 0.0,
            },
            vv_render::style::Light::default(),
            vv_render::style::Light::default(),
            vv_render::style::Light::default(),
        ],
        ambient_color: vv_render::style::Light::WHITE,
        ambient: 0.0,
        tonemap: vv_render::style::TONEMAP_NONE,
        count: 1.0,
        _pad: [0.0, 0.0],
    };
    let background = wgpu::Color {
        r: 0.5,
        g: 0.5,
        b: 0.5,
        a: 1.0,
    };
    let settings = RenderSettings {
        occlusion_culling: false,
        lighting,
        background,
        ao: 0.0,
        depth_cue: 0.0,
        fxaa: false,
        ..Default::default()
    };
    let mut renderer = Renderer::new(ctx.clone(), w, h);
    let bindings = renderer.bind_skin_surface(&skin);
    // Opacity picked so the fixed-point render target's colour, once
    // divided by a small alpha and clamped back to 1.0 before blending,
    // is where the *old*, unscaled-specular formula would blow the
    // highlight out to flat white: specular (0.55) clears alpha (~0.48),
    // triggering that clamp, while diffuse stays zero so any brightening
    // is specular alone.
    let dark = vv_render::Material {
        ambient: 0.0,
        diffuse: 0.0,
        specular: 0.0,
        shininess: 1.0,
        toon_bands: 0.0,
        opacity: 0.46,
        outline: 0.0,
        outline_width: 0.0,
        transmode: 0.0,
        rim: 0.0,
        _pad: [0.0; 2],
    };
    let glossy = vv_render::Material {
        specular: 0.55,
        ..dark
    };
    let mut render = |material: vv_render::Material| {
        let mut encoder = ctx.device.create_command_encoder(&Default::default());
        renderer.render_all(
            &mut encoder,
            &camera,
            &[],
            &[],
            &[],
            &[PatchSurfaceItem {
                surface: PatchSurface::Skin(&skin, &bindings),
                material,
            }],
            &settings,
        );
        ctx.queue.submit([encoder.finish()]);
        renderer.after_submit();
        renderer.read_color()
    };
    let base = render(dark);
    let glass = render(glossy);
    let pixel = |p: &[u8], x: u32, y: u32| p[((y * w + x) * 4) as usize] as i32;
    let (cx, cy) = (w / 2, h / 2);
    let bg = pixel(&base, 0, 0);
    let base_c = pixel(&base, cx, cy);
    let glass_c = pixel(&glass, cx, cy);
    assert!(
        base_c < bg,
        "opacity should darken the background some: bg {bg} base {base_c}"
    );
    assert!(
        glass_c > base_c,
        "a highlight should still lighten the pixel some: base {base_c} glass {glass_c}"
    );
    assert!(
        glass_c < 210,
        "specular must be scaled by opacity, not clamped to full brightness first: {glass_c}"
    );
}

/// `count` atoms, one per residue so nothing bonds them, at exact
/// `position(i)` world coordinates -- unlike `protein_like`, whose atoms
/// scatter through a packed sphere, the order-independence and blend-math
/// tests below need atoms placed by hand.
fn exact_atoms(count: usize, position: impl Fn(usize) -> glam::Vec3) -> vv_core::Structure {
    let mut b = vv_core::TopologyBuilder::with_capacity(count);
    for i in 0..count {
        b.push(&vv_core::AtomRow {
            element: vv_core::Element::CARBON,
            name: *b"C   ",
            serial: i as u32 + 1,
            alt_loc: 0,
            comp: "LIG",
            asym: "A",
            auth_asym: "A",
            seq_id: i as i32 + 1,
            auth_seq_id: i as i32 + 1,
            ins_code: 0,
            entity: 0,
            position: position(i),
            occupancy: 1.0,
            b_factor: 0.0,
            charge: 0,
            hetero: true,
        });
    }
    b.finish().expect("exact_atoms")
}

/// Packed RGBA8 (`unpack_color`'s layout: R low byte .. A high byte), opaque.
const RED: u32 = 0xFF0000FF;
const GREEN: u32 = 0xFF00FF00;
const WHITE: u32 = 0xFFFFFFFF;

/// A spacefill `DrawItem` for one of the order-independence/blend tests
/// below. A free function, not a closure, so it stays generic over the
/// borrow's lifetime at each call site instead of being pinned to one.
fn spacefill_item<'a>(
    structure: &'a GpuStructure,
    bindings: &'a vv_render::renderer::PageBindings,
    material: vv_render::Material,
) -> DrawItem<'a> {
    DrawItem {
        structure,
        bindings,
        representation: Representation::Spacefill,
        sizes: AtomSizes::of(Representation::Spacefill),
        material,
    }
}

/// A transparent atom in front of an opaque one shows both, blended more
/// toward the back atom's colour as the front atom's opacity drops --
/// the material row of docs/COMMANDS.md's `material` command, extended
/// from Gaussian/skin/SES surfaces to spacefill impostors.
#[test]
fn transparent_sphere_in_front_shows_both_atoms_blended_by_opacity() {
    let Some(ctx) = context() else { return };
    let (w, h) = (128, 128);
    let front = exact_atoms(1, |_| glam::Vec3::new(0.0, 0.0, 5.0));
    let back = exact_atoms(1, |_| glam::Vec3::new(0.0, 0.0, -5.0));
    let gpu_front = GpuStructure::upload_colored(&ctx, &front, None, &[RED]);
    let gpu_back = GpuStructure::upload_colored(&ctx, &back, None, &[GREEN]);
    let mut renderer = Renderer::new(ctx.clone(), w, h);
    let bindings_front = renderer.bind(&gpu_front);
    let bindings_back = renderer.bind(&gpu_back);
    let camera = Camera::framing(glam::Vec3::ZERO, 20.0);
    let settings = RenderSettings {
        occlusion_culling: false,
        ao: 0.0,
        depth_cue: 0.0,
        fxaa: false,
        outline: false,
        ..Default::default()
    };
    // No specular: a highlight anywhere near this on-axis scene's dead
    // centre would broadcast equally into every channel and swamp the
    // green-channel opacity signal this test reads.
    let flat = vv_render::Material {
        specular: 0.0,
        ..vv_render::Material::default()
    };
    let mut render = |front_material: vv_render::Material| {
        let items = [
            DrawItem {
                structure: &gpu_front,
                bindings: &bindings_front,
                representation: Representation::Spacefill,
                sizes: AtomSizes::of(Representation::Spacefill),
                material: front_material,
            },
            DrawItem {
                structure: &gpu_back,
                bindings: &bindings_back,
                representation: Representation::Spacefill,
                sizes: AtomSizes::of(Representation::Spacefill),
                material: flat,
            },
        ];
        let mut encoder = ctx.device.create_command_encoder(&Default::default());
        renderer.render_all(&mut encoder, &camera, &items, &[], &[], &[], &settings);
        ctx.queue.submit([encoder.finish()]);
        renderer.after_submit();
        renderer.read_color()
    };
    let (cx, cy) = (w / 2, h / 2);
    let i = ((cy * w + cx) * 4) as usize;
    let mut green_at = Vec::new();
    for opacity in [1.0, 0.8, 0.4, 0.1] {
        let pixels = render(vv_render::Material { opacity, ..flat });
        green_at.push((opacity, pixels[i], pixels[i + 1]));
    }
    let (_, opaque_r, opaque_g) = green_at[0];
    assert!(
        opaque_r > 150 && opaque_g < 40,
        "a fully opaque front atom should read as plain red, got r={opaque_r} g={opaque_g}"
    );
    for pair in green_at.windows(2) {
        let (a_op, _, a_g) = pair[0];
        let (b_op, _, b_g) = pair[1];
        assert!(
            b_g > a_g,
            "green should rise monotonically as the front atom's opacity drops: {a_op} -> {a_g}, {b_op} -> {b_g}"
        );
    }
}

/// An opaque sphere in front of a transparent one hides it completely,
/// exactly as if the transparent atom were not drawn at all.
#[test]
fn opaque_sphere_in_front_hides_a_transparent_one_behind() {
    let Some(ctx) = context() else { return };
    let (w, h) = (128, 128);
    let front = exact_atoms(1, |_| glam::Vec3::new(0.0, 0.0, 5.0));
    let back = exact_atoms(1, |_| glam::Vec3::new(0.0, 0.0, -5.0));
    let gpu_front = GpuStructure::upload_colored(&ctx, &front, None, &[GREEN]);
    let gpu_back = GpuStructure::upload_colored(&ctx, &back, None, &[RED]);
    let mut renderer = Renderer::new(ctx.clone(), w, h);
    let bindings_front = renderer.bind(&gpu_front);
    let bindings_back = renderer.bind(&gpu_back);
    let camera = Camera::framing(glam::Vec3::ZERO, 20.0);
    let settings = RenderSettings {
        occlusion_culling: false,
        ao: 0.0,
        depth_cue: 0.0,
        fxaa: false,
        outline: false,
        ..Default::default()
    };
    let front_item = DrawItem {
        structure: &gpu_front,
        bindings: &bindings_front,
        representation: Representation::Spacefill,
        sizes: AtomSizes::of(Representation::Spacefill),
        material: vv_render::Material::default(),
    };
    let back_item = DrawItem {
        structure: &gpu_back,
        bindings: &bindings_back,
        representation: Representation::Spacefill,
        sizes: AtomSizes::of(Representation::Spacefill),
        material: vv_render::Material {
            opacity: 0.4,
            ..vv_render::Material::default()
        },
    };
    let mut render = |items: &[DrawItem]| {
        let mut encoder = ctx.device.create_command_encoder(&Default::default());
        renderer.render_all(&mut encoder, &camera, items, &[], &[], &[], &settings);
        ctx.queue.submit([encoder.finish()]);
        renderer.after_submit();
        renderer.read_color()
    };
    let with_hidden_atom = render(&[front_item, back_item]);
    let front_alone = render(std::slice::from_ref(&front_item));
    assert_eq!(
        with_hidden_atom, front_alone,
        "a transparent atom fully hidden behind an opaque one must not change a single pixel"
    );
}

/// Weighted-blended OIT sums every layer's contribution, so it cannot
/// depend on which order the caller submits overlapping glass items in
/// -- unlike the nearest-layer-only Skin/SES glass path, whose owner
/// buffer picks a single winner (docs/RENDERING.md).
#[test]
fn overlapping_transparent_atoms_blend_the_same_regardless_of_draw_order() {
    let Some(ctx) = context() else { return };
    let (w, h) = (128, 128);
    let a = exact_atoms(1, |_| glam::Vec3::new(0.0, 0.0, 4.0));
    let b = exact_atoms(1, |_| glam::Vec3::new(0.0, 0.0, -4.0));
    let gpu_a = GpuStructure::upload_colored(&ctx, &a, None, &[RED]);
    let gpu_b = GpuStructure::upload_colored(&ctx, &b, None, &[GREEN]);
    let mut renderer = Renderer::new(ctx.clone(), w, h);
    let ba = renderer.bind(&gpu_a);
    let bb = renderer.bind(&gpu_b);
    let camera = Camera::framing(glam::Vec3::ZERO, 20.0);
    let settings = RenderSettings {
        occlusion_culling: false,
        ao: 0.0,
        depth_cue: 0.0,
        fxaa: false,
        outline: false,
        ..Default::default()
    };
    let material = vv_render::Material {
        opacity: 0.5,
        ..vv_render::Material::default()
    };
    let mut render = |items: &[DrawItem]| {
        let mut encoder = ctx.device.create_command_encoder(&Default::default());
        renderer.render_all(&mut encoder, &camera, items, &[], &[], &[], &settings);
        ctx.queue.submit([encoder.finish()]);
        renderer.after_submit();
        renderer.read_color()
    };
    let forward = render(&[
        spacefill_item(&gpu_a, &ba, material),
        spacefill_item(&gpu_b, &bb, material),
    ]);
    let reversed = render(&[
        spacefill_item(&gpu_b, &bb, material),
        spacefill_item(&gpu_a, &ba, material),
    ]);
    assert_eq!(
        forward, reversed,
        "the same two overlapping glass atoms must blend identically regardless of submission order"
    );
}

/// The same two overlapping glass atoms, seen from the opposite side (a
/// literal 180-degree camera rotation about the pair's midpoint) blend to
/// the same pixel: which one happens to be nearer must not change the
/// weighted-OIT result, so rotating the view never makes an overlap
/// "pop". Also confirms it is really two layers merging (alpha
/// `1 - (1-a)^2`), not the nearer one winning alone.
#[test]
fn overlapping_transparent_atoms_blend_the_same_after_a_180_degree_rotation() {
    let Some(ctx) = context() else { return };
    let (w, h) = (128, 128);
    let d = 5.0;
    let a = exact_atoms(1, |_| glam::Vec3::new(0.0, 0.0, d));
    let b = exact_atoms(1, |_| glam::Vec3::new(0.0, 0.0, -d));
    let gpu_a = GpuStructure::upload_colored(&ctx, &a, None, &[WHITE]);
    let gpu_b = GpuStructure::upload_colored(&ctx, &b, None, &[WHITE]);
    let mut renderer = Renderer::new(ctx.clone(), w, h);
    let ba = renderer.bind(&gpu_a);
    let bb = renderer.bind(&gpu_b);
    let settings = RenderSettings {
        occlusion_culling: false,
        ao: 0.0,
        depth_cue: 0.0,
        fxaa: false,
        outline: false,
        background: wgpu::Color::BLACK,
        ..Default::default()
    };
    let material = vv_render::Material {
        opacity: 0.5,
        ..vv_render::Material::default()
    };
    let mut render = |camera: &Camera, items: &[DrawItem]| {
        let mut encoder = ctx.device.create_command_encoder(&Default::default());
        renderer.render_all(&mut encoder, camera, items, &[], &[], &[], &settings);
        ctx.queue.submit([encoder.finish()]);
        renderer.after_submit();
        renderer.read_color()
    };
    let mut camera = Camera::framing(glam::Vec3::ZERO, 20.0);
    let both = [
        spacefill_item(&gpu_a, &ba, material),
        spacefill_item(&gpu_b, &bb, material),
    ];
    let front = render(&camera, &both);
    camera.orientation = glam::Quat::from_rotation_y(std::f32::consts::PI);
    let rotated = render(&camera, &both);
    let (cx, cy) = (w / 2, h / 2);
    let i = ((cy * w + cx) * 4) as usize;
    for c in 0..3 {
        assert!(
            (front[i + c] as i32 - rotated[i + c] as i32).abs() <= 3,
            "the same pair seen from the opposite side must blend the same: {:?} vs {:?}",
            &front[i..i + 3],
            &rotated[i..i + 3]
        );
    }
    camera.orientation = glam::Quat::IDENTITY;
    let single = render(
        &camera,
        std::slice::from_ref(&spacefill_item(&gpu_a, &ba, material)),
    );
    let coverage = |p: &[u8]| (0..3).map(|c| p[i + c] as i32).sum::<i32>();
    assert!(
        coverage(&front) > coverage(&single),
        "two overlapping glass layers should cover more of the black background than one: two {:?} one {:?}",
        &front[i..i + 3],
        &single[i..i + 3]
    );
}

/// A transparent atom hidden behind opaque geometry never changes an
/// opaque-only render: the same scene rendered with and without an
/// additional, fully-occluded glass atom must be pixel-identical, so
/// adding the glass pipelines never costs an opaque-only scene anything
/// it did not already pay (see docs/RENDERING.md's "Why it is fast").
#[test]
fn opaque_only_render_is_unaffected_by_an_unseen_glass_atom() {
    let Some(ctx) = context() else { return };
    let (w, h) = (192, 160);
    let structure = protein_like(&SynthParams::new(5_000));
    let gpu = GpuStructure::upload(&ctx, &structure, None, ColorScheme::Element);
    // Dead centre of the structure and far enough behind it (not an
    // arbitrary world position) so it lands on the same screen pixels as
    // the opaque atoms already covering the centre, deep enough behind
    // them that occlusion is unambiguous.
    let behind = gpu.center - glam::Vec3::Z * gpu.radius * 20.0;
    let hidden = exact_atoms(1, |_| behind);
    let gpu_hidden = GpuStructure::upload_colored(&ctx, &hidden, None, &[RED]);
    let mut renderer = Renderer::new(ctx.clone(), w, h);
    let bindings = renderer.bind(&gpu);
    let bindings_hidden = renderer.bind(&gpu_hidden);
    let camera = Camera::framing(gpu.center, gpu.radius);
    let settings = RenderSettings {
        occlusion_culling: false,
        // Fixed, not the default auto-bounds: `render_all` folds every
        // item's bounds (the hidden atom's included) into the depth
        // cue's range when this is `None`, which would shift every
        // pixel's fog term by *adding* the item, independent of
        // anything to do with transparency -- pin it so the only
        // difference between the two renders is the glass pass itself.
        scene_bounds: Some((gpu.center, gpu.radius)),
        ..Default::default()
    };
    // Spacefill, not ball-and-stick: full-radius spheres pack solidly
    // enough that the structure's centre has no background gap for the
    // hidden atom to peek through, unlike ball-and-stick's thin sticks.
    let item = DrawItem {
        structure: &gpu,
        bindings: &bindings,
        representation: Representation::Spacefill,
        sizes: AtomSizes::of(Representation::Spacefill),
        material: vv_render::Material::default(),
    };
    let glass_item = DrawItem {
        structure: &gpu_hidden,
        bindings: &bindings_hidden,
        representation: Representation::Spacefill,
        sizes: AtomSizes::of(Representation::Spacefill),
        material: vv_render::Material {
            opacity: 0.3,
            ..vv_render::Material::default()
        },
    };
    let mut render = |items: &[DrawItem]| {
        let mut encoder = ctx.device.create_command_encoder(&Default::default());
        renderer.render_all(&mut encoder, &camera, items, &[], &[], &[], &settings);
        ctx.queue.submit([encoder.finish()]);
        renderer.after_submit();
        renderer.read_color()
    };
    let without_glass = render(std::slice::from_ref(&item));
    let with_unseen_glass = render(&[item, glass_item]);
    assert_eq!(
        without_glass, with_unseen_glass,
        "an off-screen glass atom must not change one opaque pixel"
    );
}

/// Picking policy for a glass atom with nothing opaque behind it: the
/// "glass id" pass in `renderer.rs` gives it its own id, so a click on
/// empty space that only a transparent atom covers still picks that atom
/// (docs/RENDERING.md). An opaque atom drawn afterward always overwrites
/// it (`glass_shows_and_picks_what_is_behind_it` covers that half).
#[test]
fn a_transparent_atom_is_picked_when_nothing_opaque_is_behind_it() {
    let Some(ctx) = context() else { return };
    let (w, h) = (64, 64);
    let atom = exact_atoms(1, |_| glam::Vec3::ZERO);
    let gpu = GpuStructure::upload_colored(&ctx, &atom, None, &[RED]);
    let mut renderer = Renderer::new(ctx.clone(), w, h);
    let bindings = renderer.bind(&gpu);
    let camera = Camera::framing(glam::Vec3::ZERO, 5.0);
    let settings = RenderSettings {
        occlusion_culling: false,
        ..Default::default()
    };
    let item = DrawItem {
        structure: &gpu,
        bindings: &bindings,
        representation: Representation::Spacefill,
        sizes: AtomSizes::of(Representation::Spacefill),
        material: vv_render::Material {
            opacity: 0.3,
            ..vv_render::Material::default()
        },
    };
    let mut encoder = ctx.device.create_command_encoder(&Default::default());
    renderer.render_all(
        &mut encoder,
        &camera,
        std::slice::from_ref(&item),
        &[],
        &[],
        &[],
        &settings,
    );
    ctx.queue.submit([encoder.finish()]);
    renderer.after_submit();
    let (cx, cy) = (w / 2, h / 2);
    assert!(
        matches!(renderer.pick(cx, cy), Some(Pick::Atom { item: 0, atom: 0 })),
        "clicking a glass-only atom should pick it, got {:?}",
        renderer.pick(cx, cy)
    );
    assert_eq!(
        renderer.pick(0, 0),
        None,
        "background past the glass atom should still miss"
    );
}

/// Glass self-occlusion culling (`vv_render::scene::glass_cull_margin`,
/// `shaders/cull.wgsl`'s `glass_hidden`) drops a glass atom once it
/// estimates that atom's own weighted-blended OIT contribution is below
/// 1/255 once composited -- an approximation of the exact (uncontrolled)
/// render, not a match to it, visible enough on a dense scene that it is
/// only ever on with `fast_glass: true` (`RenderSettings`'s doc): the
/// interactive path, not a still view. Bound the actual pixel error it
/// introduces instead of trusting the derivation (docs/RENDERING.md's
/// "Transparency", "Glass self-occlusion culling").
#[test]
fn glass_self_occlusion_culling_stays_within_a_measured_pixel_error_bound() {
    let Some(ctx) = context() else { return };
    let (w, h) = (256, 256);
    let structure = protein_like(&SynthParams::new(20_000));
    let gpu = GpuStructure::upload(&ctx, &structure, None, ColorScheme::Element);
    let mut renderer = Renderer::new(ctx.clone(), w, h);
    let bindings = renderer.bind(&gpu);
    let mut camera = Camera::framing(gpu.center, gpu.radius);
    camera.distance *= 0.6; // real depth complexity to cull against
    let material = vv_render::Material {
        opacity: 0.35,
        ..vv_render::Material::default()
    };
    let item = spacefill_item(&gpu, &bindings, material);
    let mut render = |occlusion_culling: bool, frames: usize| {
        let settings = RenderSettings {
            occlusion_culling,
            // Only ever matters while `occlusion_culling` is also on
            // (`glass_only_hidden` checks `cam.occlusion` first), so
            // leaving it on here doesn't affect the exact reference below.
            fast_glass: true,
            ao: 0.0,
            depth_cue: 0.0,
            fxaa: false,
            outline: false,
            ..Default::default()
        };
        let mut pixels = Vec::new();
        for _ in 0..frames {
            let mut encoder = ctx.device.create_command_encoder(&Default::default());
            renderer.render_all(
                &mut encoder,
                &camera,
                std::slice::from_ref(&item),
                &[],
                &[],
                &[],
                &settings,
            );
            ctx.queue.submit([encoder.finish()]);
            renderer.after_submit();
            pixels = renderer.read_color();
        }
        pixels
    };
    // The glass depth reference is last frame's (`glass_hidden`'s doc),
    // so culling only starts engaging from the second frame on.
    let culled = render(true, 5);
    let reference = render(false, 1);
    let (mut sum_err, mut max_err, mut n) = (0u64, 0u8, 0u64);
    for (a, b) in culled.chunks(4).zip(reference.chunks(4)) {
        for c in 0..3 {
            let d = a[c].abs_diff(b[c]);
            sum_err += d as u64;
            max_err = max_err.max(d);
            n += 1;
        }
    }
    let mean_err = sum_err as f64 / n as f64;
    eprintln!(
        "glass self-occlusion culling: mean per-channel error {mean_err:.3}/255, max {max_err}/255"
    );
    assert!(
        mean_err <= 6.0,
        "mean per-channel error {mean_err:.3}/255 exceeds the bound"
    );
    assert!(
        max_err <= 235,
        "a single pixel differs by {max_err}/255 from the uncontrolled render"
    );
}

/// The glass self-occlusion cull reference is last frame's, with no
/// phase-2 revalidation the way the opaque Hi-Z gets
/// (`vv_render::scene::glass_cull_margin`'s doc), so a wrongly-culled
/// atom stays wrong for a whole frame instead of being caught before it
/// is drawn. Sweep the camera through a dense glass scene and compare
/// how much consecutive frames differ against the same sweep with
/// culling off: ordinary camera motion (parallax, silhouette shift)
/// already changes pixels frame to frame, so a real regression is extra
/// jumpiness *culling itself* adds on top of that baseline (an atom's
/// worth of colour popping in or out), not the baseline's own size.
#[test]
fn glass_self_occlusion_culling_has_no_popping_across_a_180_degree_sweep() {
    let Some(ctx) = context() else { return };
    let (w, h) = (192, 192);
    let structure = protein_like(&SynthParams::new(20_000));
    let gpu = GpuStructure::upload(&ctx, &structure, None, ColorScheme::Element);
    let mut renderer = Renderer::new(ctx.clone(), w, h);
    let bindings = renderer.bind(&gpu);
    let material = vv_render::Material {
        opacity: 0.35,
        ..vv_render::Material::default()
    };
    let item = spacefill_item(&gpu, &bindings, material);
    let mut base = Camera::framing(gpu.center, gpu.radius);
    base.distance *= 0.6;
    const STEPS: usize = 36;
    let mut sweep = |occlusion_culling: bool| {
        let settings = RenderSettings {
            occlusion_culling,
            fast_glass: true,
            ao: 0.0,
            depth_cue: 0.0,
            fxaa: false,
            outline: false,
            ..Default::default()
        };
        let mut previous: Option<Vec<u8>> = None;
        let mut worst_step_mean = 0.0f64;
        for step in 0..=STEPS {
            let angle = std::f32::consts::PI * step as f32 / STEPS as f32;
            let mut camera = base.clone();
            camera.orientation = glam::Quat::from_rotation_y(angle);
            let mut encoder = ctx.device.create_command_encoder(&Default::default());
            renderer.render_all(
                &mut encoder,
                &camera,
                std::slice::from_ref(&item),
                &[],
                &[],
                &[],
                &settings,
            );
            ctx.queue.submit([encoder.finish()]);
            renderer.after_submit();
            let pixels = renderer.read_color();
            if let Some(prev) = &previous {
                let mut sum = 0u64;
                for (a, b) in pixels.chunks(4).zip(prev.chunks(4)) {
                    for c in 0..3 {
                        sum += a[c].abs_diff(b[c]) as u64;
                    }
                }
                let mean = sum as f64 / (pixels.len() / 4 * 3) as f64;
                worst_step_mean = worst_step_mean.max(mean);
            }
            previous = Some(pixels);
        }
        worst_step_mean
    };
    let worst_culled = sweep(true);
    let worst_reference = sweep(false);
    eprintln!(
        "glass self-occlusion culling: worst consecutive-frame mean error \
         {worst_culled:.3}/255 culled vs {worst_reference:.3}/255 reference over {STEPS} steps"
    );
    assert!(
        worst_culled <= worst_reference + 5.0,
        "culling adds {:.3}/255 of consecutive-frame jumpiness on top of the \
         {worst_reference:.3}/255 baseline from camera motion alone -- looks like popping",
        worst_culled - worst_reference
    );
}

/// `fast_glass: false` (`RenderSettings::default`, what a still view,
/// screenshot, export or `Renderer::render` caller gets unless it opts
/// in) must never approximate glass, even if a stale "glass hiz"
/// reference from an earlier `fast_glass: true` frame is still sitting
/// in the texture: the margin computation itself is gated, not just
/// whether that texture gets refreshed, so turning `fast_glass` off
/// mid-session falls back to exact immediately, not once the old
/// reference happens to clear.
#[test]
fn glass_self_occlusion_culling_is_exact_when_fast_glass_is_off() {
    let Some(ctx) = context() else { return };
    let (w, h) = (256, 256);
    let structure = protein_like(&SynthParams::new(20_000));
    let gpu = GpuStructure::upload(&ctx, &structure, None, ColorScheme::Element);
    let mut renderer = Renderer::new(ctx.clone(), w, h);
    let bindings = renderer.bind(&gpu);
    let mut camera = Camera::framing(gpu.center, gpu.radius);
    camera.distance *= 0.6;
    let material = vv_render::Material {
        opacity: 0.35,
        ..vv_render::Material::default()
    };
    let item = spacefill_item(&gpu, &bindings, material);
    let mut render = |occlusion_culling: bool, fast_glass: bool, frames: usize| {
        let settings = RenderSettings {
            occlusion_culling,
            fast_glass,
            ao: 0.0,
            depth_cue: 0.0,
            fxaa: false,
            outline: false,
            ..Default::default()
        };
        let mut pixels = Vec::new();
        for _ in 0..frames {
            let mut encoder = ctx.device.create_command_encoder(&Default::default());
            renderer.render_all(
                &mut encoder,
                &camera,
                std::slice::from_ref(&item),
                &[],
                &[],
                &[],
                &settings,
            );
            ctx.queue.submit([encoder.finish()]);
            renderer.after_submit();
            pixels = renderer.read_color();
        }
        pixels
    };
    render(true, true, 3); // prime a "glass hiz" reference, then abandon it
    let still = render(true, false, 4);
    let exact_reference = render(false, false, 1);
    // Not bit-exact: weighted-blended OIT sums overlapping glass
    // fragments in whatever order the cull pass's atomic compaction
    // happens to produce, which the GPU does not guarantee is the same
    // list order every dispatch, so floating-point addition alone gives
    // a few pixels a difference of 1-2/255 between *any* two renders of
    // this scene -- present even for two `occlusion_culling: false`
    // renders back to back. Bound tightly around that noise floor
    // instead: an order of magnitude below the mean-2/max-108 error the
    // quality-bound test measures with `fast_glass: true` actually
    // culling something.
    let (mut sum_err, mut max_err) = (0u64, 0u8);
    for (a, b) in still.iter().zip(&exact_reference) {
        let d = a.abs_diff(*b);
        sum_err += d as u64;
        max_err = max_err.max(d);
    }
    let mean_err = sum_err as f64 / still.len() as f64;
    eprintln!("fast_glass off: mean error {mean_err:.4}/255, max {max_err}/255");
    assert!(
        mean_err <= 0.5 && max_err <= 8,
        "fast_glass: false differs from the exact reference by mean {mean_err:.4}/255 max \
         {max_err}/255 -- more than OIT's own draw-order noise floor, looks like it is still \
         approximating glass"
    );
}

/// A transparent cartoon ribbon in front of an opaque backdrop shows the
/// backdrop through it, exactly like a transparent sphere -- the
/// weighted-blended OIT path extends to real triangle meshes
/// (`fs_cartoon_glass`), not just impostors.
#[test]
fn transparent_cartoon_shows_the_opaque_atom_behind_it() {
    let Some(ctx) = context() else { return };
    let path =
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/small/1CRN.cif");
    let structure = vv_io::load(path).unwrap();
    let coords = structure.frame(0);
    let positions = coords.positions();
    let codes = vv_core::dssp::assign(&structure.topology, positions);
    let plan = vv_core::cartoon::plan(&structure.topology, positions, &codes);
    let frame = plan.frame(positions);
    let colors = vec![vv_render::color::rgba(220, 30, 30); structure.atom_count()];
    let (center, radius) = vv_core::CoordSet::new(positions.to_vec())
        .bounding_sphere()
        .unwrap();
    let (w, h) = (256, 256);
    let camera = Camera::framing(center, radius);
    let settings = RenderSettings {
        occlusion_culling: false,
        depth_cue: 0.0,
        ao: 0.0,
        fxaa: false,
        outline: false,
        ..Default::default()
    };
    // A backdrop sphere well behind the ribbon along the view axis, not
    // concentric with it: same-centred and bigger (as
    // `cartoon_depth_tests_correctly_against_an_impostor_sphere` uses)
    // would *enclose* the ribbon and hide it, the opposite of a backdrop.
    // Orthographic has no perspective falloff, so pushing it back does
    // not need a bigger radius to keep covering the same screen area.
    let giant = GpuStructure::from_parts(
        &ctx,
        &[center + camera.forward() * radius * 5.0],
        &[radius * 1.5],
        &[vv_render::color::rgba(30, 200, 30)],
        &[],
    );
    let mut renderer = Renderer::new(ctx.clone(), w, h);
    let gpu_mesh = CartoonGpu::upload(&ctx, &plan, &frame, &colors).expect("upload cartoon");
    let cartoon_bindings = renderer.bind_cartoon(&gpu_mesh);
    let bindings = renderer.bind(&giant);
    let backdrop_item = DrawItem {
        structure: &giant,
        bindings: &bindings,
        representation: Representation::Spacefill,
        sizes: AtomSizes::of(Representation::Spacefill),
        material: vv_render::Material::default(),
    };
    let mut render = |cartoons: &[CartoonItem]| {
        let mut encoder = ctx.device.create_command_encoder(&Default::default());
        renderer.render_all(
            &mut encoder,
            &camera,
            std::slice::from_ref(&backdrop_item),
            cartoons,
            &[],
            &[],
            &settings,
        );
        ctx.queue.submit([encoder.finish()]);
        renderer.after_submit();
        renderer.read_color()
    };
    let backdrop_only = render(&[]);
    let opaque = render(&[CartoonItem {
        mesh: CartoonMesh::Ribbon(&gpu_mesh),
        bindings: &cartoon_bindings,
        material: vv_render::Material::default(),
    }]);
    let glass = render(&[CartoonItem {
        mesh: CartoonMesh::Ribbon(&gpu_mesh),
        bindings: &cartoon_bindings,
        material: vv_render::Material {
            opacity: 0.3,
            ..vv_render::Material::default()
        },
    }]);
    let ribbon_pixel = (0..(w * h) as usize)
        .find(|&idx| {
            let p = idx * 4;
            opaque[p..p + 3] != backdrop_only[p..p + 3]
        })
        .expect("the opaque ribbon should draw over the backdrop somewhere");
    let p = ribbon_pixel * 4;
    assert_ne!(
        glass[p..p + 3],
        opaque[p..p + 3],
        "a transparent ribbon should let the backdrop show through"
    );
    assert!(
        glass[p + 1] > opaque[p + 1],
        "green (the backdrop's colour) should rise once the ribbon is transparent: opaque {:?} glass {:?}",
        &opaque[p..p + 3],
        &glass[p..p + 3]
    );
}

/// Screen-space AO barely darkens a cartoon (a ribbon's occluders sit
/// behind it, out of the depth buffer); the occlusion volume must.
#[test]
fn occlusion_volume_darkens_a_cartoon() {
    let Some(ctx) = context() else { return };
    let path =
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/small/4HHB.pdb");
    let structure = vv_io::load(path).unwrap();
    let coords = structure.frame(0);
    let positions = coords.positions();
    let codes = vv_core::dssp::assign(&structure.topology, positions);
    let plan = vv_core::cartoon::plan(&structure.topology, positions, &codes);
    let frame = plan.frame(positions);
    let mesh = plan.mesh(&frame);
    let colors = vec![vv_render::color::rgba(200, 200, 200); structure.atom_count()];
    let gpu = CartoonGpu::upload(&ctx, &plan, &frame, &colors).expect("upload cartoon");
    let spheres: Vec<_> = mesh
        .sections
        .iter()
        .map(|s| (s.center, 0.5 * (s.half_width + s.half_thickness)))
        .collect();
    let volume = std::sync::Arc::new(vv_render::OcclusionVolume::build(&ctx, &spheres).unwrap());

    let (center, radius) = vv_core::CoordSet::new(positions.to_vec())
        .bounding_sphere()
        .unwrap();
    let camera = Camera::framing(center, radius);
    let (w, h) = (256, 256);
    let mut renderer = Renderer::new(ctx.clone(), w, h);
    let bindings = renderer.bind_cartoon(&gpu);
    let frame = |renderer: &mut Renderer, ao: f32| {
        let settings = RenderSettings {
            occlusion_culling: false,
            background: PINHOLE_BG,
            depth_cue: 0.0,
            ao,
            ..Default::default()
        };
        let mut encoder = ctx.device.create_command_encoder(&Default::default());
        renderer.render_all(
            &mut encoder,
            &camera,
            &[],
            &[CartoonItem {
                mesh: CartoonMesh::Ribbon(&gpu),
                bindings: &bindings,
                material: vv_render::Material::default(),
            }],
            &[],
            &[],
            &settings,
        );
        ctx.queue.submit([encoder.finish()]);
        renderer.after_submit();
        renderer.read_color()
    };
    let darkening = |plain: &[u8], shaded: &[u8]| {
        let (mut sum, mut n) = (0i64, 0i64);
        for (a, b) in plain.chunks(4).zip(shaded.chunks(4)) {
            if a[..3] != [255, 0, 255] {
                sum += (0..3).map(|c| a[c] as i64 - b[c] as i64).sum::<i64>();
                n += 1;
            }
        }
        sum as f64 / n.max(1) as f64
    };
    let plain = frame(&mut renderer, 0.0);
    let screen_only = frame(&mut renderer, 2.0);
    renderer.set_occlusion_volume(Some(volume));
    let with_volume = frame(&mut renderer, 2.0);
    let (ss, vol) = (
        darkening(&plain, &screen_only),
        darkening(&plain, &with_volume),
    );
    // Measured 65 -> 95 on 4HHB.
    assert!(
        vol > 1.3 * ss,
        "volume adds too little: {ss:.1} -> {vol:.1}"
    );
}

/// Half of a synthetic structure plus its copy turned 180 degrees about
/// the vertical axis through `pivot`, all atoms one colour: it looks the
/// same from the front and the back.
fn twofold(half: usize) -> (vv_core::Structure, glam::Vec3) {
    let source = protein_like(&SynthParams::new(2 * half));
    let mut topology = (*source.topology).clone();
    topology.element.copy_within(..half, half);
    let first = source.frame(0).positions()[..half].to_vec();
    let pivot = first.iter().copied().sum::<glam::Vec3>() / half as f32;
    let turned = first.iter().map(|&p| {
        let d = p - pivot;
        pivot + glam::Vec3::new(-d.x, d.y, -d.z)
    });
    let positions = first.iter().copied().chain(turned).collect();
    let structure = vv_core::Structure::new(topology, vv_core::CoordSet::new(positions)).unwrap();
    (structure, pivot)
}

/// The occlusion volume is world space, so it must darken a structure the
/// same whichever side faces the camera: a twofold-symmetric structure
/// seen from the front and from the back shades alike. Mips whose sizes
/// round down misplace coarse occupancy toward one corner of the grid, a
/// dark side that turns with the structure (26% of pixels differed).
#[test]
fn occlusion_volume_shades_both_sides_of_a_symmetric_structure_alike() {
    let Some(ctx) = context() else { return };
    let (w, h) = (256, 256);
    let (structure, pivot) = twofold(3_000);
    let colors = vec![vv_render::color::rgba(200, 200, 200); structure.atom_count()];
    let gpu = GpuStructure::upload_colored(&ctx, &structure, None, &colors);
    let spheres: Vec<_> = structure
        .frame(0)
        .positions()
        .iter()
        .zip(&structure.topology.element)
        .map(|(&p, e)| (p, e.vdw_radius()))
        .collect();
    let volume = vv_render::OcclusionVolume::build(&ctx, &spheres).unwrap();
    let mut renderer = Renderer::new(ctx.clone(), w, h);
    renderer.set_occlusion_volume(Some(std::sync::Arc::new(volume)));
    let bindings = renderer.bind(&gpu);
    let style = StylePreset::default();
    let settings = RenderSettings {
        occlusion_culling: false,
        ao: style.ao(),
        depth_cue: style.depth_cue(),
        ..Default::default()
    };
    let item = DrawItem {
        structure: &gpu,
        bindings: &bindings,
        representation: Representation::default(),
        sizes: AtomSizes::of(Representation::default()),
        material: style.material(),
    };
    let mut front = Camera::framing(pivot, gpu.radius);
    front.target = pivot;
    let mut back = front.clone();
    back.orbit(std::f32::consts::PI, 0.0);
    let items = std::slice::from_ref(&item);
    let a = render_items(&ctx, &mut renderer, &front, items, &[], &settings);
    let b = render_items(&ctx, &mut renderer, &back, items, &[], &settings);

    let (mut differing, mut signed) = (0usize, 0i64);
    for (p, q) in a.chunks(4).zip(b.chunks(4)) {
        if (0..3).any(|c| p[c].abs_diff(q[c]) > 2) {
            differing += 1;
        }
        signed += (0..3).map(|c| p[c] as i64 - q[c] as i64).sum::<i64>();
    }
    let total = (w * h) as usize;
    let mean = signed as f64 / total as f64;
    // Measured 0.6% and 0.06; the rest is the cones' tangent
    // frame, which the turn does not carry along.
    assert!(
        differing * 20 < total && mean.abs() < 1.0,
        "front and back differ: {differing} of {total} pixels, mean {mean:.2}"
    );
}

/// A clip plane through a packed structure's centre, facing the camera:
/// the cut shows capped cross-sections (the centre stays covered and
/// pickable, not a hole to the background), and the frame changes.
#[test]
fn clip_plane_caps_what_it_cuts() {
    let Some(ctx) = context() else { return };
    let (w, h) = (128, 128);
    let structure = protein_like(&SynthParams::new(3_000));
    let gpu = GpuStructure::upload(&ctx, &structure, None, ColorScheme::Element);
    let mut camera = Camera::framing(gpu.center, gpu.radius);
    camera.projection = Projection::Orthographic;
    let mut renderer = Renderer::new(ctx.clone(), w, h);
    let bindings = renderer.bind(&gpu);
    let forward = camera.orientation * glam::Vec3::NEG_Z;
    let frame = |renderer: &mut Renderer, clip: Option<[f32; 4]>| {
        let settings = RenderSettings {
            occlusion_culling: false,
            background: PINHOLE_BG,
            clip,
            ..Default::default()
        };
        let mut encoder = ctx.device.create_command_encoder(&Default::default());
        renderer.render(&mut encoder, &camera, &gpu, &bindings, &settings);
        ctx.queue.submit([encoder.finish()]);
        renderer.after_submit();
        renderer.read_color()
    };
    let whole = frame(&mut renderer, None);
    let plane = forward.extend(-forward.dot(gpu.center)).to_array();
    let cut = frame(&mut renderer, Some(plane));
    let (x, y) = (w / 2, h / 2);
    let i = ((y * w + x) * 4) as usize;
    assert_ne!(
        cut[i..i + 3],
        [255, 0, 255],
        "the cut centre is capped, not a hole"
    );
    assert!(renderer.pick(x, y).is_some(), "a cap is pickable");
    let differing = whole
        .chunks(4)
        .zip(cut.chunks(4))
        .filter(|(a, b)| a != b)
        .count();
    assert!(
        differing * 4 > (w * h) as usize,
        "only {differing} pixels changed"
    );
}

/// Cavities and the gap between atoms stay open: a clip plane through two
/// separated atoms caps each one where it straddles the plane, but does
/// not paint the empty span between them.
#[test]
fn clip_plane_does_not_fill_the_gap_between_atoms() {
    let Some(ctx) = context() else { return };
    let base = protein_like(&SynthParams::new(8));
    let hidden = glam::Vec3::splat(10_000.0);
    let mut positions = vec![hidden; 8];
    // Atoms 1 and 5 are both carbon (same van der Waals radius): far
    // enough apart in X that their discs, and the gap between them, never
    // overlap on screen.
    positions[1] = glam::Vec3::new(-15.0, 0.0, 0.0);
    positions[5] = glam::Vec3::new(15.0, 0.0, 0.0);
    let structure =
        vv_core::Structure::new((*base.topology).clone(), vv_core::CoordSet::new(positions))
            .unwrap();
    let gpu = GpuStructure::upload(&ctx, &structure, None, ColorScheme::Element);

    let (w, h) = (400u32, 400u32);
    let aspect = w as f32 / h as f32;
    let mut camera = Camera::framing(glam::Vec3::ZERO, 30.0);
    camera.projection = Projection::Orthographic;
    let pixel_of = |p: glam::Vec3| -> (u32, u32) {
        let clip = camera.proj(aspect) * camera.view() * glam::Vec4::new(p.x, p.y, p.z, 1.0);
        let ndc = clip.truncate() / clip.w;
        (
            ((ndc.x * 0.5 + 0.5) * w as f32) as u32,
            ((0.5 - ndc.y * 0.5) * h as f32) as u32,
        )
    };

    let mut renderer = Renderer::new(ctx.clone(), w, h);
    let bindings = renderer.bind(&gpu);
    // Through both atoms' centres (Z=0), facing the camera: the near half
    // of each sphere is cut away and shows its equatorial cap.
    let forward = camera.orientation * glam::Vec3::NEG_Z;
    let plane = forward.extend(-forward.dot(glam::Vec3::ZERO)).to_array();
    let settings = RenderSettings {
        occlusion_culling: false,
        background: PINHOLE_BG,
        clip: Some(plane),
        ..Default::default()
    };
    let mut encoder = ctx.device.create_command_encoder(&Default::default());
    renderer.render(&mut encoder, &camera, &gpu, &bindings, &settings);
    ctx.queue.submit([encoder.finish()]);
    renderer.after_submit();
    let pixels = renderer.read_color();

    let is_background = |p: glam::Vec3| -> bool {
        let (x, y) = pixel_of(p);
        let i = ((y * w + x) * 4) as usize;
        pixels[i..i + 3] == [255, 0, 255]
    };
    assert!(
        !is_background(glam::Vec3::new(-15.0, 0.0, 0.0)),
        "the left atom's cap should show, not a hole"
    );
    assert!(
        !is_background(glam::Vec3::new(15.0, 0.0, 0.0)),
        "the right atom's cap should show, not a hole"
    );
    assert!(
        is_background(glam::Vec3::ZERO),
        "the empty gap between the atoms must not be painted over"
    );
}

/// The viewport's atom-sphere caps (`draw.wgsl`'s `fs_sphere`) and
/// `render`'s (`path_trace.wgsl`'s `hit_sphere`) agree on where the clip
/// plane leaves a cap versus background: same technique
/// (`shading.wgsl`'s `kept_span`), so coverage should match almost
/// everywhere; shadows and AO are off (`unshadowed`) so only the cut
/// itself, not sampling noise, can disagree.
#[test]
fn traced_spacefill_clip_matches_the_raster_coverage() {
    let Some(ctx) = context() else { return };
    let (w, h) = (240u32, 240u32);
    let structure = protein_like(&SynthParams::new(400));
    let positions = structure.frame(0).positions().to_vec();
    let radii: Vec<f32> = structure
        .topology
        .element
        .iter()
        .map(|e| e.vdw_radius())
        .collect();
    let colors = vv_render::colors_for(ColorScheme::Element, &structure.topology);
    let (center, radius) = structure.frame(0).bounding_sphere().unwrap();
    let mut camera = Camera::framing(center, radius);
    camera.orbit(0.4, 0.2);
    let forward = camera.orientation * glam::Vec3::NEG_Z;
    let plane = forward.extend(-forward.dot(camera.target)).to_array();

    let settings = RenderSettings {
        occlusion_culling: false,
        background: PINHOLE_BG,
        fxaa: false,
        outline: false,
        ao: 0.0,
        depth_cue: 0.0,
        clip: Some(plane),
        ..Default::default()
    };
    let gpu = GpuStructure::upload(&ctx, &structure, None, ColorScheme::Element);
    let mut renderer = Renderer::new(ctx.clone(), w, h);
    let bindings = renderer.bind(&gpu);
    let mut encoder = ctx.device.create_command_encoder(&Default::default());
    renderer.render(&mut encoder, &camera, &gpu, &bindings, &settings);
    ctx.queue.submit([encoder.finish()]);
    renderer.after_submit();
    let raster = renderer.read_color();

    let mut scene = vv_render::path_trace::TraceScene::default();
    scene.push_atoms(
        &positions,
        &radii,
        &colors,
        0..positions.len(),
        std::iter::empty(),
        0.0,
        vv_render::Material::default(),
    );
    let tracer = vv_render::path_trace::PathTracer::new(&ctx, &scene).unwrap();
    let trace_settings = vv_render::path_trace::TraceSettings {
        width: w,
        height: h,
        samples: 32,
        lighting: vv_render::Lighting::default(),
        background: [1.0, 0.0, 1.0],
        background_top: [1.0, 0.0, 1.0],
        transparent: false,
        clip: Some(plane),
        unshadowed: true,
        light_spread: vv_render::path_trace::DEFAULT_LIGHT_SPREAD,
        ao: 1.0,
        direct: 1.0,
    };
    let traced = tracer.render(&ctx, &camera, &trace_settings, |_| {});

    let mut coverage = 0;
    let mut agree = 0;
    for i in (0..(w * h) as usize).map(|k| k * 4) {
        let raster_bg = raster[i..i + 3] == [255, 0, 255];
        let traced_bg = traced[i] > 240 && traced[i + 1] < 16 && traced[i + 2] > 240;
        if raster_bg && traced_bg {
            continue;
        }
        coverage += 1;
        agree += (raster_bg == traced_bg) as usize;
    }
    assert!(
        coverage > 500,
        "too little non-background coverage: {coverage}"
    );
    assert!(
        agree * 100 >= coverage * 95,
        "{agree}/{coverage} pixels agree on cut coverage"
    );
}

/// Depth of field blurs what is away from the focus and leaves a frame
/// with it off untouched.
#[test]
fn depth_of_field_blurs_only_when_on() {
    let Some(ctx) = context() else { return };
    let (w, h) = (128, 128);
    let structure = protein_like(&SynthParams::new(3_000));
    let gpu = GpuStructure::upload(&ctx, &structure, None, ColorScheme::Element);
    let camera = Camera::framing(gpu.center, gpu.radius);
    let mut renderer = Renderer::new(ctx.clone(), w, h);
    let bindings = renderer.bind(&gpu);
    let frame = |renderer: &mut Renderer, dof: f32| {
        let settings = RenderSettings {
            occlusion_culling: false,
            dof,
            ..Default::default()
        };
        let mut encoder = ctx.device.create_command_encoder(&Default::default());
        renderer.render(&mut encoder, &camera, &gpu, &bindings, &settings);
        ctx.queue.submit([encoder.finish()]);
        renderer.after_submit();
        renderer.read_color()
    };
    let sharp = frame(&mut renderer, 0.0);
    let blurred = frame(&mut renderer, 1.0);
    assert_eq!(
        sharp,
        frame(&mut renderer, 0.0),
        "dof 0 must not change the frame"
    );
    // Blur lowers local contrast: sum of neighbour differences drops.
    let contrast = |px: &[u8]| -> i64 {
        (0..(w * h - 1) as usize)
            .map(|i| (px[i * 4] as i64 - px[(i + 1) * 4] as i64).abs())
            .sum()
    };
    let (a, b) = (contrast(&sharp), contrast(&blurred));
    assert!(b * 10 < a * 9, "contrast {a} -> {b}");
}

/// Two reps of the same atoms in one frame, each with its own draw state:
/// spacefill by element and ball-and-stick by chain (inside the spheres,
/// so the pair looks exactly like spacefill alone). Sharing per-draw
/// state between them drew both with whichever rep wrote last.
#[test]
fn two_reps_of_one_structure_keep_their_own_state() {
    let Some(ctx) = context() else { return };
    let (w, h) = (160, 120);
    let structure = protein_like(&SynthParams::new(2_000));
    let bonds = vv_core::bonds::perceive(&structure.topology, structure.frame(0).positions());
    let gpu = GpuStructure::upload(&ctx, &structure, Some(&bonds), ColorScheme::Element);
    let chain_colors = vv_render::colors_for(ColorScheme::Chain, &structure.topology);
    let chain_state = vv_render::DrawState::new(&ctx, &gpu, &chain_colors);
    let camera = Camera::framing(gpu.center, gpu.radius);
    let mut renderer = Renderer::new(ctx.clone(), w, h);
    let spacefill_bindings = renderer.bind(&gpu);
    let sticks_bindings = renderer.bind_state(&gpu, &chain_state);
    let spacefill = DrawItem {
        structure: &gpu,
        bindings: &spacefill_bindings,
        representation: Representation::Spacefill,
        sizes: AtomSizes::of(Representation::Spacefill),
        material: vv_render::Material::default(),
    };
    let sticks = DrawItem {
        structure: &gpu,
        bindings: &sticks_bindings,
        representation: Representation::BallAndStick,
        sizes: AtomSizes::of(Representation::BallAndStick),
        material: vv_render::MaterialPreset::Glossy.into(),
    };
    let settings = RenderSettings {
        occlusion_culling: false,
        ..Default::default()
    };
    let frame = |renderer: &mut Renderer, items: &[DrawItem]| {
        for _ in 0..2 {
            let mut encoder = ctx.device.create_command_encoder(&Default::default());
            renderer.render_all(&mut encoder, &camera, items, &[], &[], &[], &settings);
            ctx.queue.submit([encoder.finish()]);
            renderer.after_submit();
        }
        renderer.read_color()
    };
    let alone = frame(&mut renderer, std::slice::from_ref(&spacefill));
    let sticks_alone = frame(&mut renderer, std::slice::from_ref(&sticks));
    assert_ne!(alone, sticks_alone, "the two reps should look different");
    let both = frame(&mut renderer, &[spacefill, sticks]);
    let differing = alone
        .chunks(4)
        .zip(both.chunks(4))
        .filter(|(a, b)| a != b)
        .count();
    assert!(
        differing * 1000 < (w * h) as usize,
        "{differing} pixels changed when the sticks joined"
    );
    let (x, y) = (w / 2, h / 2);
    assert!(
        matches!(renderer.pick(x, y), Some(Pick::Atom { item: 0, .. })),
        "the spheres are in front: {:?}",
        renderer.pick(x, y)
    );
}

/// The SES shader against its CPU twin (`vv_core::ses::Ses::cast`, itself
/// checked against an independent oracle): per pixel, the same ray must
/// hit or miss on both, in both projections.
#[test]
fn ses_surface_matches_the_cpu_cast() {
    let Some(ctx) = context() else { return };
    let (w, h) = (120u32, 90u32);
    let path =
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/small/1CRN.pdb");
    let structure = vv_io::load(&path).unwrap();
    let positions = structure.frame(0).positions().to_vec();
    let radii: Vec<f32> = structure
        .topology
        .element
        .iter()
        .map(|e| e.vdw_radius())
        .collect();
    let colors = vv_render::colors_for(ColorScheme::Element, &structure.topology);
    let ses = vv_core::ses::build(&positions, &radii, 1.4);
    let layout = vv_render::SesLayout::new(&ses, &positions, &radii);
    let gpu = vv_render::SesGpu::upload(&ctx, &layout, &colors).unwrap();
    let mut renderer = Renderer::new(ctx.clone(), w, h);
    let bindings = renderer.bind_ses_surface(&gpu);
    let settings = RenderSettings {
        occlusion_culling: true,
        background: vv_render::wgpu::Color {
            r: 1.0,
            g: 0.0,
            b: 1.0,
            a: 1.0,
        },
        fxaa: false,
        outline: false,
        ao: 0.0,
        depth_cue: 0.0,
        ..Default::default()
    };
    for projection in [Projection::Perspective, Projection::Orthographic] {
        let mut camera = Camera::framing(gpu.bounds_center, gpu.bounds_radius * 0.6);
        camera.orbit(0.6, 0.3);
        camera.projection = projection;
        let mut encoder = ctx.device.create_command_encoder(&Default::default());
        renderer.render_all(
            &mut encoder,
            &camera,
            &[],
            &[],
            &[],
            &[PatchSurfaceItem {
                surface: PatchSurface::Ses(&gpu, &bindings),
                material: vv_render::Material::default(),
            }],
            &settings,
        );
        ctx.queue.submit([encoder.finish()]);
        renderer.after_submit();
        let pixels = renderer.read_color();

        let view_proj_inv = (camera.proj(w as f32 / h as f32) * camera.view()).inverse();
        let eye = camera.eye();
        let forward = (camera.orientation * glam::Vec3::NEG_Z).normalize();
        let (mut agree, mut drawn) = (0usize, 0usize);
        for y in 0..h {
            for x in 0..w {
                let ndc = glam::Vec2::new(
                    (x as f32 + 0.5) / w as f32 * 2.0 - 1.0,
                    1.0 - (y as f32 + 0.5) / h as f32 * 2.0,
                );
                let on = view_proj_inv.project_point3(ndc.extend(0.5));
                let (origin, dir) = match projection {
                    Projection::Perspective => (eye, (on - eye).normalize()),
                    Projection::Orthographic => (on - forward * 1000.0, forward),
                };
                let cpu = ses.cast(&positions, &radii, origin, dir, 0.0).is_some();
                let i = ((y * w + x) * 4) as usize;
                let p = &pixels[i..i + 3];
                let gpu_hit = !(p[0] > 240 && p[1] < 16 && p[2] > 240);
                agree += (cpu == gpu_hit) as usize;
                drawn += gpu_hit as usize;
            }
        }
        let total = (w * h) as usize;
        eprintln!("{projection:?}: {agree}/{total} agree, {drawn} drawn");
        assert!(
            drawn > total / 5,
            "{projection:?}: only {drawn} pixels drawn"
        );
        assert!(
            agree * 1000 >= total * 995,
            "{projection:?}: {agree}/{total} agree"
        );
    }
}

/// Every patch's bounding sphere, in `Ses::cast`'s order: convex, tori,
/// concave.
fn ses_patch_bounds(
    ses: &vv_core::ses::Ses,
    positions: &[glam::Vec3],
    radii: &[f32],
) -> Vec<(glam::Vec3, f32)> {
    let rp = ses.probe_radius;
    let convex = ses.convex.iter().map(|c| {
        let a = c.atom as usize;
        (positions[a], radii[a])
    });
    let tori = ses.tori.iter().map(|t| {
        let [i, j] = t.atoms.map(|a| a as usize);
        let c =
            vv_core::ses::circle(positions[i], radii[i] + rp, positions[j], radii[j] + rp).unwrap();
        vv_core::ses::torus_bound(&c, radii[i], radii[j], rp, ses.torus_arc(t, &c))
    });
    let concave = ses.probes.iter().map(|p| (p.center, rp));
    convex.chain(tori).chain(concave).collect()
}

/// `Ses::cast`, testing only the patches whose bounds the ray passes: the
/// nearest hit as (t, whether it faces the ray), fast enough per pixel.
fn ses_cast_culled(
    ses: &vv_core::ses::Ses,
    bounds: &[(glam::Vec3, f32)],
    positions: &[glam::Vec3],
    radii: &[f32],
    o: glam::Vec3,
    d: glam::Vec3,
    t_min: f32,
) -> Option<(f32, bool)> {
    let (nc, nt) = (ses.convex.len(), ses.tori.len());
    bounds
        .iter()
        .enumerate()
        .filter(|(_, &(c, r))| (c - o).reject_from_normalized(d).length_squared() <= r * r)
        .filter_map(|(k, _)| {
            if k < nc {
                ses.hit_convex(&ses.convex[k], positions, radii, o, d, t_min)
            } else if k < nc + nt {
                ses.hit_torus(&ses.tori[k - nc], positions, radii, o, d, t_min)
            } else {
                ses.hit_concave(k - nc - nt, positions, o, d, t_min)
            }
        })
        .min_by(|a, b| a.0.total_cmp(&b.0))
        .map(|(t, n)| (t, n.dot(d) < 0.0))
}

/// Holes in the SES viewport, per pixel against the CPU cast's depth, at
/// views that lost patches: fillets skipped by their own small bounds
/// (or drawn as points) while their atoms cut caps open, and patch bounds
/// holding the eye. Pixels on the CPU surface's creases and silhouettes,
/// and where the eye sits inside the excluded volume, are left out.
#[test]
fn ses_surface_has_no_holes() {
    let Some(ctx) = context() else { return };
    let software = ctx.adapter.get_info().device_type == vv_render::wgpu::DeviceType::Cpu;
    let (w, h) = (320u32, 240u32);
    let cases = [
        ("4HHB.cif", Projection::Orthographic, 1.0, 1),
        ("4HHB.cif", Projection::Orthographic, 0.6, 2),
        ("4HHB.cif", Projection::Perspective, 0.6, 1),
        ("1CRN.pdb", Projection::Perspective, 0.12, 4),
    ];
    let settings = RenderSettings {
        background: vv_render::wgpu::Color {
            r: 1.0,
            g: 0.0,
            b: 1.0,
            a: 1.0,
        },
        fxaa: false,
        outline: false,
        ao: 0.0,
        depth_cue: 0.0,
        ..Default::default()
    };
    for (file, projection, zoom, view) in cases {
        let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../fixtures/small")
            .join(file);
        let structure = vv_io::load(&path).unwrap();
        let positions = structure.frame(0).positions().to_vec();
        let radii: Vec<f32> = structure
            .topology
            .element
            .iter()
            .map(|e| e.vdw_radius())
            .collect();
        let colors = vv_render::colors_for(ColorScheme::Element, &structure.topology);
        let ses = vv_core::ses::build(&positions, &radii, vv_core::ses::WATER_PROBE);
        let bounds = ses_patch_bounds(&ses, &positions, &radii);
        let layout = vv_render::SesLayout::new(&ses, &positions, &radii);
        let gpu = vv_render::SesGpu::upload(&ctx, &layout, &colors).unwrap();
        let mut renderer = Renderer::new(ctx.clone(), w, h);
        let bindings = renderer.bind_ses_surface(&gpu);

        let (center, radius) = structure.frame(0).bounding_sphere().unwrap();
        let mut camera = Camera::framing(center, radius);
        camera.projection = projection;
        camera.orbit(view as f32 * 0.83, (view as f32 * 1.37).sin() * 0.9);
        camera.zoom(zoom);
        // The second frame's occlusion cull reads the first's depth.
        for _ in 0..2 {
            let mut encoder = ctx.device.create_command_encoder(&Default::default());
            renderer.render_all(
                &mut encoder,
                &camera,
                &[],
                &[],
                &[],
                &[PatchSurfaceItem {
                    surface: PatchSurface::Ses(&gpu, &bindings),
                    material: vv_render::Material::default(),
                }],
                &settings,
            );
            ctx.queue.submit([encoder.finish()]);
            renderer.after_submit();
        }
        let depth = renderer.read_depth();

        let proj = camera.proj(w as f32 / h as f32);
        let eye = camera.eye();
        let (right, up) = (
            camera.orientation * glam::Vec3::X,
            camera.orientation * glam::Vec3::Y,
        );
        let forward = camera.orientation * glam::Vec3::NEG_Z;
        let tan = (camera.fov_y * 0.5).tan();
        let cut = camera.near_cut().unwrap_or(0.0);
        // View depth of the CPU's front-facing hit through each pixel centre,
        // past the camera's near cut.
        let cpu: Vec<Option<f32>> = (0..w * h)
            .into_par_iter()
            .map(|k| {
                let x = ((k % w) as f32 + 0.5) / w as f32 * 2.0 - 1.0;
                let y = 1.0 - ((k / w) as f32 + 0.5) / h as f32 * 2.0;
                let screen = right * (x * tan * w as f32 / h as f32) + up * (y * tan);
                let (o, d) = match projection {
                    Projection::Perspective => (eye, (screen + forward).normalize()),
                    Projection::Orthographic => (eye + screen * camera.distance, forward),
                };
                ses_cast_culled(
                    &ses,
                    &bounds,
                    &positions,
                    &radii,
                    o,
                    d,
                    cut / d.dot(forward),
                )
                .and_then(|(t, front)| front.then_some(t * d.dot(forward)))
            })
            .collect();
        let (p22, p32, p23, p33) = (proj.z_axis.z, proj.w_axis.z, proj.z_axis.w, proj.w_axis.w);
        let gpu_depth = |z: f32| (z > 0.0).then(|| (z * p33 - p32) / (z * p23 - p22));

        let (mut smooth, mut holes, mut off) = (0usize, 0usize, 0usize);
        for y in 1..h - 1 {
            for x in 1..w - 1 {
                let k = (y * w + x) as usize;
                let Some(want) = cpu[k] else { continue };
                let flat = (-1i32..=1).all(|dy| {
                    (-1i32..=1).all(|dx| {
                        let n = ((y as i32 + dy) * w as i32 + x as i32 + dx) as usize;
                        cpu[n].is_some_and(|t| (t - want).abs() < 1.0)
                    })
                });
                if !flat {
                    continue;
                }
                smooth += 1;
                match gpu_depth(depth[k]) {
                    None => holes += 1,
                    Some(got) => off += ((got - want).abs() > 0.3) as usize,
                }
            }
        }
        let case = format!("{file} {projection:?} zoom {zoom} view {view}");
        eprintln!("{case}: {holes} holes, {off}/{smooth} off");
        assert!(
            smooth > (w * h / 5) as usize,
            "{case}: {smooth} smooth pixels"
        );
        if software {
            // Lavapipe leaves a few dozen scattered pixels open at seams
            // that hardware rasterizers close; the regressions this guards
            // against opened thousands.
            assert!(holes * 1000 <= smooth, "{case}: {holes} holes");
        } else {
            assert_eq!(holes, 0, "{case}");
        }
        assert!(off * 1000 <= smooth, "{case}: {off}/{smooth} off");
    }
}

/// The SES baked into a volume (`GaussianSurfaceGpu::ses`) against the
/// analytic cast: per pixel the same ray hits or misses, except along the
/// silhouette, where the volume is off by up to about a voxel.
#[test]
fn ses_volume_matches_the_cpu_cast() {
    let Some(ctx) = context() else { return };
    let (w, h) = (120u32, 90u32);
    let path =
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/small/1CRN.pdb");
    let structure = vv_io::load(&path).unwrap();
    let positions = structure.frame(0).positions().to_vec();
    let radii: Vec<f32> = structure
        .topology
        .element
        .iter()
        .map(|e| e.vdw_radius())
        .collect();
    let colors = vv_render::colors_for(ColorScheme::Element, &structure.topology);
    let ses = vv_core::ses::build(&positions, &radii, 1.4);
    let gpu = GaussianSurfaceGpu::ses(&ctx, &ses, &positions, &radii, &colors).unwrap();
    let mut renderer = Renderer::new(ctx.clone(), w, h);
    let bindings = renderer.bind_gaussian_surface(&gpu);
    let settings = RenderSettings {
        occlusion_culling: false,
        background: vv_render::wgpu::Color {
            r: 1.0,
            g: 0.0,
            b: 1.0,
            a: 1.0,
        },
        fxaa: false,
        outline: false,
        ao: 0.0,
        depth_cue: 0.0,
        ..Default::default()
    };
    for projection in [Projection::Perspective, Projection::Orthographic] {
        let mut camera = Camera::framing(gpu.bounds_center, gpu.bounds_radius * 0.6);
        camera.orbit(0.6, 0.3);
        camera.projection = projection;
        let mut encoder = ctx.device.create_command_encoder(&Default::default());
        renderer.render_all(
            &mut encoder,
            &camera,
            &[],
            &[],
            &[GaussianSurfaceItem {
                gpu: &gpu,
                bindings: &bindings,
                material: vv_render::Material::default(),
            }],
            &[],
            &settings,
        );
        ctx.queue.submit([encoder.finish()]);
        renderer.after_submit();
        let pixels = renderer.read_color();

        let view_proj_inv = (camera.proj(w as f32 / h as f32) * camera.view()).inverse();
        let eye = camera.eye();
        let forward = (camera.orientation * glam::Vec3::NEG_Z).normalize();
        let (mut agree, mut drawn) = (0usize, 0usize);
        for y in 0..h {
            for x in 0..w {
                let ndc = glam::Vec2::new(
                    (x as f32 + 0.5) / w as f32 * 2.0 - 1.0,
                    1.0 - (y as f32 + 0.5) / h as f32 * 2.0,
                );
                let on = view_proj_inv.project_point3(ndc.extend(0.5));
                let (origin, dir) = match projection {
                    Projection::Perspective => (eye, (on - eye).normalize()),
                    Projection::Orthographic => (on - forward * 1000.0, forward),
                };
                let cpu = ses.cast(&positions, &radii, origin, dir, 0.0).is_some();
                let i = ((y * w + x) * 4) as usize;
                let p = &pixels[i..i + 3];
                let gpu_hit = !(p[0] > 240 && p[1] < 16 && p[2] > 240);
                agree += (cpu == gpu_hit) as usize;
                drawn += gpu_hit as usize;
            }
        }
        let total = (w * h) as usize;
        eprintln!("{projection:?}: {agree}/{total} agree, {drawn} drawn");
        assert!(
            drawn > total / 5,
            "{projection:?}: only {drawn} pixels drawn"
        );
        assert!(
            agree * 100 >= total * 98,
            "{projection:?}: {agree}/{total} agree"
        );
    }
}

/// The path tracer with every shadow and occlusion ray let through shades
/// exactly as the viewport does: same rays from the same camera (both
/// projections), normals, colours, lights and tonemap. Compared away from
/// silhouettes, which the tracer antialiases and the raster does not.
#[test]
fn unshadowed_path_trace_matches_the_raster() {
    let Some(ctx) = context() else { return };
    let (w, h) = (160u32, 120u32);
    let path =
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/small/1CRN.pdb");
    let structure = vv_io::load(&path).unwrap();
    let positions = structure.frame(0).positions().to_vec();
    let colors = vv_render::colors_for(ColorScheme::Element, &structure.topology);
    let gpu = GpuStructure::upload(&ctx, &structure, None, ColorScheme::Element);
    let mut renderer = Renderer::new(ctx.clone(), w, h);
    let bindings = renderer.bind(&gpu);
    let settings = RenderSettings {
        occlusion_culling: false,
        background: vv_render::wgpu::Color {
            r: 1.0,
            g: 0.0,
            b: 1.0,
            a: 1.0,
        },
        fxaa: false,
        outline: false,
        ao: 0.0,
        depth_cue: 0.0,
        ..Default::default()
    };
    let scene = vv_render::path_trace::TraceScene {
        spheres: positions
            .iter()
            .zip(&structure.topology.element)
            .zip(&colors)
            .map(
                |((&center, e), &color)| vv_render::path_trace::TraceSphere {
                    center,
                    radius: e.vdw_radius(),
                    color,
                    material: 0,
                },
            )
            .collect(),
        cylinders: Vec::new(),
        materials: vec![settings.material],
        ..Default::default()
    };
    let tracer = vv_render::path_trace::PathTracer::new(&ctx, &scene).unwrap();
    for projection in [Projection::Perspective, Projection::Orthographic] {
        let mut camera = Camera::framing(gpu.center, gpu.radius);
        camera.orbit(0.6, 0.3);
        camera.projection = projection;
        let mut encoder = ctx.device.create_command_encoder(&Default::default());
        renderer.render(&mut encoder, &camera, &gpu, &bindings, &settings);
        ctx.queue.submit([encoder.finish()]);
        renderer.after_submit();
        let raster = renderer.read_color();
        let traced = tracer.render(
            &ctx,
            &camera,
            &vv_render::path_trace::TraceSettings {
                width: w,
                height: h,
                samples: 1,
                lighting: settings.lighting,
                background: [0.0; 3],
                background_top: [0.0; 3],
                transparent: true,
                clip: None,
                unshadowed: true,
                light_spread: vv_render::path_trace::DEFAULT_LIGHT_SPREAD,
                ao: 1.0,
                direct: 1.0,
            },
            |_| {},
        );
        let (mut interior, mut close) = (0usize, 0usize);
        for i in (0..(w * h) as usize).map(|k| k * 4) {
            let r = &raster[i..i + 3];
            let background = r[0] > 240 && r[1] < 16 && r[2] > 240;
            if background || traced[i + 3] != 255 {
                continue;
            }
            interior += 1;
            let diff = (0..3)
                .map(|c| (r[c] as i32 - traced[i + c] as i32).abs())
                .max()
                .unwrap();
            close += (diff <= 6) as usize;
        }
        eprintln!("{projection:?}: {close}/{interior} interior pixels within 6 levels");
        assert!(
            interior > (w * h / 5) as usize,
            "{projection:?}: {interior}"
        );
        assert!(
            close * 100 >= interior * 99,
            "{projection:?}: {close}/{interior}"
        );
    }
}

/// Shadow and occlusion rays only take light away: traced with them,
/// crambin is darker overall than traced without, and no pixel gets
/// brighter beyond sampling noise.
#[test]
fn shadows_and_occlusion_only_darken() {
    let Some(ctx) = context() else { return };
    let (w, h) = (128u32, 96u32);
    let path =
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/small/1CRN.pdb");
    let structure = vv_io::load(&path).unwrap();
    let positions = structure.frame(0).positions().to_vec();
    let colors = vv_render::colors_for(ColorScheme::Element, &structure.topology);
    let radii: Vec<f32> = structure
        .topology
        .element
        .iter()
        .map(|e| e.vdw_radius())
        .collect();
    let mut scene = vv_render::path_trace::TraceScene::default();
    scene.push_atoms(
        &positions,
        &radii,
        &colors,
        0..positions.len(),
        std::iter::empty(),
        0.0,
        vv_render::Material::default(),
    );
    let tracer = vv_render::path_trace::PathTracer::new(&ctx, &scene).unwrap();
    let (lo, hi) = positions.iter().fold(
        (glam::Vec3::splat(f32::MAX), glam::Vec3::splat(f32::MIN)),
        |(lo, hi), &p| (lo.min(p), hi.max(p)),
    );
    let mut camera = Camera::framing((lo + hi) * 0.5, (hi - lo).length() * 0.5);
    camera.orbit(0.6, 0.3);
    let render = |unshadowed| {
        tracer.render(
            &ctx,
            &camera,
            &vv_render::path_trace::TraceSettings {
                width: w,
                height: h,
                samples: 64,
                lighting: vv_render::LightingPreset::Full.lighting(),
                background: [0.0; 3],
                background_top: [0.0; 3],
                transparent: true,
                clip: None,
                unshadowed,
                light_spread: vv_render::path_trace::DEFAULT_LIGHT_SPREAD,
                ao: 1.0,
                direct: 1.0,
            },
            |_| {},
        )
    };
    let (open, shaded) = (render(true), render(false));
    let (mut sum_open, mut sum_shaded, mut brighter, mut covered) = (0u64, 0u64, 0usize, 0usize);
    for i in (0..(w * h) as usize).map(|k| k * 4) {
        if open[i + 3] != 255 {
            continue;
        }
        covered += 1;
        let (a, b) = (
            open[i..i + 3].iter().map(|&c| c as u64).sum::<u64>(),
            shaded[i..i + 3].iter().map(|&c| c as u64).sum::<u64>(),
        );
        sum_open += a;
        sum_shaded += b;
        brighter += (b > a + 3 * 12) as usize;
    }
    eprintln!("mean {sum_open} -> {sum_shaded} over {covered} pixels, {brighter} brighter");
    assert!(covered > (w * h / 10) as usize);
    assert!(
        sum_shaded * 100 < sum_open * 90,
        "{sum_open} -> {sum_shaded}"
    );
    assert!(brighter * 100 < covered, "{brighter}/{covered} brighter");
}

/// The traced cartoon (the raster's own mesh, `CartoonMesh::expand`)
/// with shadow rays off shades as the raster cartoon does.
#[test]
fn unshadowed_traced_cartoon_matches_the_raster() {
    let Some(ctx) = context() else { return };
    let (w, h) = (640u32, 480u32);
    let path =
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/small/1UBQ.cif");
    let structure = vv_io::load(path).unwrap();
    let coords = structure.frame(0);
    let positions = coords.positions();
    let codes = vv_core::dssp::assign(&structure.topology, positions);
    let plan = vv_core::cartoon::plan(&structure.topology, positions, &codes);
    let spline = plan.frame(positions);
    let colors = vv_render::colors_for(ColorScheme::Rainbow, &structure.topology);
    let gpu = CartoonGpu::upload(&ctx, &plan, &spline, &colors).unwrap();
    let (center, radius) = coords.bounding_sphere().unwrap();
    let mut camera = Camera::framing(center, radius * 0.8);
    camera.orbit(0.4, 0.2);
    // No highlight: its high exponent turns the GPU's and CPU's last-bit
    // differences in the ring normals into visible shifts.
    let matte = vv_render::Material {
        specular: 0.0,
        ..Default::default()
    };
    let settings = RenderSettings {
        occlusion_culling: false,
        background: vv_render::wgpu::Color {
            r: 1.0,
            g: 0.0,
            b: 1.0,
            a: 1.0,
        },
        fxaa: false,
        outline: false,
        ao: 0.0,
        depth_cue: 0.0,
        ..Default::default()
    };
    let mut renderer = Renderer::new(ctx.clone(), w, h);
    let bindings = renderer.bind_cartoon(&gpu);
    let mut encoder = ctx.device.create_command_encoder(&Default::default());
    renderer.render_all(
        &mut encoder,
        &camera,
        &[],
        &[CartoonItem {
            mesh: CartoonMesh::Ribbon(&gpu),
            bindings: &bindings,
            material: matte,
        }],
        &[],
        &[],
        &settings,
    );
    ctx.queue.submit([encoder.finish()]);
    renderer.after_submit();
    let raster = renderer.read_color();

    let mut scene = vv_render::path_trace::TraceScene::default();
    scene.push_mesh(&plan.mesh(&spline).expand(), &colors, matte);
    let tracer = vv_render::path_trace::PathTracer::new(&ctx, &scene).unwrap();
    let traced = tracer.render(
        &ctx,
        &camera,
        &vv_render::path_trace::TraceSettings {
            width: w,
            height: h,
            samples: 1,
            lighting: settings.lighting,
            background: [0.0; 3],
            background_top: [0.0; 3],
            transparent: true,
            clip: None,
            unshadowed: true,
            light_spread: vv_render::path_trace::DEFAULT_LIGHT_SPREAD,
            ao: 1.0,
            direct: 1.0,
        },
        |_| {},
    );
    let (mut interior, mut close) = (0usize, 0usize);
    for i in (0..(w * h) as usize).map(|k| k * 4) {
        let r = &raster[i..i + 3];
        let background = r[0] > 240 && r[1] < 16 && r[2] > 240;
        if background || traced[i + 3] != 255 {
            continue;
        }
        interior += 1;
        let diff = (0..3)
            .map(|c| (r[c] as i32 - traced[i + c] as i32).abs())
            .max()
            .unwrap();
        close += (diff <= 6) as usize;
    }
    eprintln!("{close}/{interior} interior pixels within 6 levels");
    assert!(interior > 1500, "{interior}");
    assert!(close * 100 >= interior * 98, "{close}/{interior}");
}

/// The traced SES and skin surface (their patches, hit by the viewport's
/// own WGSL) with shadow rays off shade as the viewport's; with shadows
/// on, no pixel gets brighter (what self-intersection on the concave and
/// saddle patches would show as speckle).
#[test]
fn traced_ses_and_skin_match_the_raster() {
    let Some(ctx) = context() else { return };
    let (w, h) = (320u32, 240u32);
    let path =
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/small/1CRN.pdb");
    let structure = vv_io::load(&path).unwrap();
    let coords = structure.frame(0);
    let positions = coords.positions().to_vec();
    let radii: Vec<f32> = structure
        .topology
        .element
        .iter()
        .map(|e| e.vdw_radius())
        .collect();
    let colors = vv_render::colors_for(ColorScheme::Element, &structure.topology);
    let (center, radius) = coords.bounding_sphere().unwrap();
    let mut camera = Camera::framing(center, radius * 0.9);
    camera.orbit(0.5, 0.3);
    let matte = vv_render::Material {
        specular: 0.0,
        ..Default::default()
    };
    let settings = RenderSettings {
        occlusion_culling: false,
        background: vv_render::wgpu::Color {
            r: 1.0,
            g: 0.0,
            b: 1.0,
            a: 1.0,
        },
        fxaa: false,
        outline: false,
        ao: 0.0,
        depth_cue: 0.0,
        ..Default::default()
    };
    let ses = vv_core::ses::build(&positions, &radii, vv_core::ses::WATER_PROBE);
    let layout = vv_render::SesLayout::new(&ses, &positions, &radii);
    let ses_gpu = vv_render::SesGpu::upload(&ctx, &layout, &colors).unwrap();
    let weights = skin_weights(&structure);
    let complex = vv_core::skin_surface::build_complex(
        &positions,
        &weights,
        vv_core::skin_surface::DEFAULT_SHRINK,
    );
    let skin_gpu = SkinSurfaceGpu::from_complex(&ctx, &complex, &positions, &weights, &colors);
    for kind in ["ses", "skin"] {
        let mut renderer = Renderer::new(ctx.clone(), w, h);
        let mut scene = vv_render::path_trace::TraceScene::default();
        let raster = if kind == "ses" {
            scene.push_ses(&ses, &positions, &radii, &colors, matte);
            let bindings = renderer.bind_ses_surface(&ses_gpu);
            raster_patches(
                &ctx,
                &mut renderer,
                &camera,
                PatchSurface::Ses(&ses_gpu, &bindings),
                matte,
                &settings,
            )
        } else {
            scene.push_skin(&complex, &positions, &weights, &colors, matte);
            let bindings = renderer.bind_skin_surface(&skin_gpu);
            raster_patches(
                &ctx,
                &mut renderer,
                &camera,
                PatchSurface::Skin(&skin_gpu, &bindings),
                matte,
                &settings,
            )
        };
        let tracer = vv_render::path_trace::PathTracer::new(&ctx, &scene).unwrap();
        let trace = |unshadowed, samples| {
            tracer.render(
                &ctx,
                &camera,
                &vv_render::path_trace::TraceSettings {
                    width: w,
                    height: h,
                    samples,
                    lighting: settings.lighting,
                    background: [0.0; 3],
                    background_top: [0.0; 3],
                    transparent: true,
                    clip: None,
                    unshadowed,
                    light_spread: vv_render::path_trace::DEFAULT_LIGHT_SPREAD,
                    ao: 1.0,
                    direct: 1.0,
                },
                |_| {},
            )
        };
        let open = trace(true, 1);
        let (mut interior, mut close) = (0usize, 0usize);
        for i in (0..(w * h) as usize).map(|k| k * 4) {
            let r = &raster[i..i + 3];
            let background = r[0] > 240 && r[1] < 16 && r[2] > 240;
            if background || open[i + 3] != 255 {
                continue;
            }
            interior += 1;
            let diff = (0..3)
                .map(|c| (r[c] as i32 - open[i + c] as i32).abs())
                .max()
                .unwrap();
            close += (diff <= 6) as usize;
        }
        eprintln!("{kind}: {close}/{interior} interior pixels within 6 levels");
        assert!(interior > (w * h / 5) as usize, "{kind}: {interior}");
        assert!(close * 100 >= interior * 99, "{kind}: {close}/{interior}");

        let (open, shaded) = (trace(true, 32), trace(false, 32));
        let brighter = (0..(w * h) as usize)
            .map(|k| k * 4)
            .filter(|&i| open[i + 3] == 255)
            .filter(|&i| (0..3).any(|c| shaded[i + c] as i32 > open[i + c] as i32 + 12))
            .count();
        eprintln!("{kind}: {brighter} pixels brighter with shadows");
        assert!(
            brighter * 1000 < interior,
            "{kind}: {brighter}/{interior} brighter"
        );
    }
}

/// The traced Gaussian surface (`vv_core::gaussian_mesh`) covers the
/// pixels the viewport's ray-marched volume does: per pixel, both hit or
/// both miss, except along silhouettes, where the viewport's voxels show.
#[test]
fn traced_gaussian_mesh_covers_the_raster_surface() {
    let Some(ctx) = context() else { return };
    let (w, h) = (320u32, 240u32);
    let path =
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/small/1CRN.pdb");
    let structure = vv_io::load(&path).unwrap();
    let positions = structure.frame(0).positions().to_vec();
    let radii: Vec<f32> = structure
        .topology
        .element
        .iter()
        .map(|e| e.vdw_radius())
        .collect();
    let colors = vv_render::colors_for(ColorScheme::Element, &structure.topology);
    let blob = vv_core::gaussian_surface::DEFAULT_BLOB_FACTOR;
    let gpu = GaussianSurfaceGpu::upload(&ctx, &positions, &radii, &colors, blob).unwrap();
    let mut renderer = Renderer::new(ctx.clone(), w, h);
    let bindings = renderer.bind_gaussian_surface(&gpu);
    let settings = RenderSettings {
        occlusion_culling: false,
        background: vv_render::wgpu::Color {
            r: 1.0,
            g: 0.0,
            b: 1.0,
            a: 1.0,
        },
        fxaa: false,
        outline: false,
        ao: 0.0,
        depth_cue: 0.0,
        ..Default::default()
    };
    let mut camera = Camera::framing(gpu.bounds_center, gpu.bounds_radius * 0.4);
    camera.orbit(0.5, 0.2);
    let mut encoder = ctx.device.create_command_encoder(&Default::default());
    renderer.render_all(
        &mut encoder,
        &camera,
        &[],
        &[],
        &[GaussianSurfaceItem {
            gpu: &gpu,
            bindings: &bindings,
            material: vv_render::Material::default(),
        }],
        &[],
        &settings,
    );
    ctx.queue.submit([encoder.finish()]);
    renderer.after_submit();
    let raster = renderer.read_color();

    let mesh = vv_core::gaussian_mesh::mesh(
        &positions,
        &radii,
        blob,
        vv_render::scene::GAUSSIAN_EPSILON,
        0.5,
    );
    let mut scene = vv_render::path_trace::TraceScene::default();
    scene.push_mesh(&mesh, &colors, vv_render::Material::default());
    let traced = vv_render::path_trace::PathTracer::new(&ctx, &scene)
        .unwrap()
        .render(
            &ctx,
            &camera,
            &vv_render::path_trace::TraceSettings {
                width: w,
                height: h,
                samples: 1,
                lighting: settings.lighting,
                background: [0.0; 3],
                background_top: [0.0; 3],
                transparent: true,
                clip: None,
                unshadowed: true,
                light_spread: vv_render::path_trace::DEFAULT_LIGHT_SPREAD,
                ao: 1.0,
                direct: 1.0,
            },
            |_| {},
        );
    let (mut agree, mut drawn) = (0usize, 0usize);
    for i in (0..(w * h) as usize).map(|k| k * 4) {
        let r = &raster[i..i + 3];
        let raster_hit = !(r[0] > 240 && r[1] < 16 && r[2] > 240);
        let traced_hit = traced[i + 3] > 127;
        agree += (raster_hit == traced_hit) as usize;
        drawn += traced_hit as usize;
    }
    let total = (w * h) as usize;
    eprintln!("{agree}/{total} agree, {drawn} drawn");
    assert!(drawn > total / 5, "{drawn}");
    assert!(agree * 100 >= total * 97, "{agree}/{total}");
}

/// One frame of a patch surface alone, read back.
fn raster_patches(
    ctx: &GpuContext,
    renderer: &mut Renderer,
    camera: &Camera,
    surface: PatchSurface,
    material: vv_render::Material,
    settings: &RenderSettings,
) -> Vec<u8> {
    let mut encoder = ctx.device.create_command_encoder(&Default::default());
    renderer.render_all(
        &mut encoder,
        camera,
        &[],
        &[],
        &[],
        &[PatchSurfaceItem { surface, material }],
        settings,
    );
    ctx.queue.submit([encoder.finish()]);
    renderer.after_submit();
    renderer.read_color()
}

/// The baked, ray-marched Gaussian surface against the exact field
/// (`vv_core::gaussian_surface::density`, every atom, no cutoff) marched
/// on the CPU along the same rays: per pixel, both hit or both miss,
/// except along silhouettes, where a voxel's worth of difference shows.
#[test]
fn gaussian_surface_matches_the_exact_field() {
    let Some(ctx) = context() else { return };
    let (w, h) = (64u32, 48u32);
    let path =
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/small/1CRN.pdb");
    let structure = vv_io::load(&path).unwrap();
    let positions = structure.frame(0).positions().to_vec();
    let radii: Vec<f32> = structure
        .topology
        .element
        .iter()
        .map(|e| e.vdw_radius())
        .collect();
    let colors = vv_render::colors_for(ColorScheme::Element, &structure.topology);
    let blob = vv_core::gaussian_surface::DEFAULT_BLOB_FACTOR;
    let gpu = GaussianSurfaceGpu::upload(&ctx, &positions, &radii, &colors, blob).unwrap();
    let mut renderer = Renderer::new(ctx.clone(), w, h);
    let bindings = renderer.bind_gaussian_surface(&gpu);
    let settings = RenderSettings {
        occlusion_culling: false,
        background: vv_render::wgpu::Color {
            r: 1.0,
            g: 0.0,
            b: 1.0,
            a: 1.0,
        },
        fxaa: false,
        outline: false,
        ao: 0.0,
        depth_cue: 0.0,
        ..Default::default()
    };
    let mut camera = Camera::framing(gpu.bounds_center, gpu.bounds_radius * 0.3);
    camera.orbit(0.5, 0.2);
    let mut encoder = ctx.device.create_command_encoder(&Default::default());
    renderer.render_all(
        &mut encoder,
        &camera,
        &[],
        &[],
        &[GaussianSurfaceItem {
            gpu: &gpu,
            bindings: &bindings,
            material: vv_render::Material::default(),
        }],
        &[],
        &settings,
    );
    ctx.queue.submit([encoder.finish()]);
    renderer.after_submit();
    let pixels = renderer.read_color();

    let blobs: Vec<_> = positions
        .iter()
        .zip(&radii)
        .map(|(&center, &radius)| vv_core::gaussian_surface::Blob { center, radius })
        .collect();
    let view_proj_inv = (camera.proj(w as f32 / h as f32) * camera.view()).inverse();
    let eye = camera.eye();
    let forward = (camera.orientation * glam::Vec3::NEG_Z).normalize();
    let (mut agree, mut drawn) = (0usize, 0usize);
    for y in 0..h {
        for x in 0..w {
            let ndc = glam::Vec2::new(
                (x as f32 + 0.5) / w as f32 * 2.0 - 1.0,
                1.0 - (y as f32 + 0.5) / h as f32 * 2.0,
            );
            let on = view_proj_inv.project_point3(ndc.extend(0.5));
            let (origin, dir) = match camera.projection {
                Projection::Perspective => (eye, (on - eye).normalize()),
                Projection::Orthographic => (on - forward * 1000.0, forward),
            };
            // March the exact field through the surface's bounding sphere.
            let oc = origin - gpu.bounds_center;
            let b = oc.dot(dir);
            let disc = b * b - (oc.length_squared() - gpu.bounds_radius * gpu.bounds_radius);
            let cpu = disc > 0.0 && {
                let (t0, t1) = (-b - disc.sqrt(), -b + disc.sqrt());
                let mut t = t0.max(0.0);
                let mut hit = false;
                while t < t1 && !hit {
                    hit = vv_core::gaussian_surface::density(origin + dir * t, &blobs, blob) >= 1.0;
                    t += 0.15;
                }
                hit
            };
            let i = ((y * w + x) * 4) as usize;
            let p = &pixels[i..i + 3];
            let gpu_hit = !(p[0] > 240 && p[1] < 16 && p[2] > 240);
            agree += (cpu == gpu_hit) as usize;
            drawn += gpu_hit as usize;
        }
    }
    let total = (w * h) as usize;
    eprintln!("{agree}/{total} agree, {drawn} drawn");
    assert!(drawn > total / 5, "only {drawn} drawn");
    assert!(agree * 100 >= total * 98, "{agree}/{total} agree");
}

/// A Gaussian surface baked at a new frame in place (`set_frame`) and
/// back draws exactly as the original upload; atoms leaving the volume's
/// box are refused, for a fresh upload.
#[test]
fn a_gaussian_surface_bakes_new_frames_in_place() {
    let Some(ctx) = context() else { return };
    let path =
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/small/1CRN.pdb");
    let structure = vv_io::load(&path).unwrap();
    let start = structure.frame(0).positions().to_vec();
    let radii: Vec<f32> = structure
        .topology
        .element
        .iter()
        .map(|e| e.vdw_radius())
        .collect();
    let colors = vv_render::colors_for(ColorScheme::Element, &structure.topology);
    let mut gpu = GaussianSurfaceGpu::upload(&ctx, &start, &radii, &colors, 2.0).unwrap();
    let (w, h) = (160u32, 120u32);
    let mut renderer = Renderer::new(ctx.clone(), w, h);
    let camera = Camera::framing(gpu.bounds_center, gpu.bounds_radius * 0.7);
    let settings = RenderSettings {
        occlusion_culling: false,
        ..Default::default()
    };
    let mut draw = |gpu: &GaussianSurfaceGpu| {
        let bindings = renderer.bind_gaussian_surface(gpu);
        let mut encoder = ctx.device.create_command_encoder(&Default::default());
        renderer.render_all(
            &mut encoder,
            &camera,
            &[],
            &[],
            &[GaussianSurfaceItem {
                gpu,
                bindings: &bindings,
                material: vv_render::Material::default(),
            }],
            &[],
            &settings,
        );
        ctx.queue.submit([encoder.finish()]);
        renderer.after_submit();
        renderer.read_color()
    };
    let before = draw(&gpu);
    let wiggled: Vec<glam::Vec3> = start
        .iter()
        .map(|p| *p + glam::Vec3::new((p.y * 0.3).sin(), 0.0, (p.x * 0.2).cos()) * 1.5)
        .collect();
    assert!(gpu.set_frame(&ctx, &wiggled));
    let moved = draw(&gpu);
    assert!(gpu.set_frame(&ctx, &start));
    let back = draw(&gpu);
    assert_ne!(before, moved, "the new frame is drawn");
    assert_eq!(before, back, "and the old one again, exactly");
    let away: Vec<glam::Vec3> = start.iter().map(|p| *p + glam::Vec3::X * 50.0).collect();
    assert!(!gpu.set_frame(&ctx, &away), "outside the box");
}

/// A playing trajectory moves a cartoon with `CartoonGpu::set_frame`
/// (the spline uploaded, `shaders/cartoon_frame.wgsl` rebuilding the
/// sections): it must draw what uploading the new frame does.
#[test]
fn a_moved_cartoon_draws_like_a_fresh_one() {
    let Some(ctx) = context() else { return };
    let path =
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/small/1UBQ.cif");
    let structure = vv_io::load(path).unwrap();
    let coords = structure.frame(0);
    let before = coords.positions();
    let codes = vv_core::dssp::assign(&structure.topology, before);
    let plan = vv_core::cartoon::plan(&structure.topology, before, &codes);
    // A bend and a shift, as a trajectory frame might bring.
    let after: Vec<glam::Vec3> = before
        .iter()
        .map(|p| *p + glam::Vec3::new((p.y * 0.2).sin(), 0.5, (p.x * 0.15).cos()))
        .collect();
    let (start, moved) = (plan.frame(before), plan.frame(&after));
    let colors = vec![vv_render::color::rgba(200, 80, 80); structure.atom_count()];
    let (center, radius) = vv_core::CoordSet::new(after.clone())
        .bounding_sphere()
        .unwrap();
    let camera = Camera::framing(center, radius);
    let settings = RenderSettings {
        occlusion_culling: false,
        ..Default::default()
    };
    let (w, h) = (192, 192);
    let draw = |gpu: &CartoonGpu| {
        let mut renderer = Renderer::new(ctx.clone(), w, h);
        let bindings = renderer.bind_cartoon(gpu);
        let mut encoder = ctx.device.create_command_encoder(&Default::default());
        renderer.render_all(
            &mut encoder,
            &camera,
            &[],
            &[CartoonItem {
                mesh: CartoonMesh::Ribbon(gpu),
                bindings: &bindings,
                material: vv_render::Material::default(),
            }],
            &[],
            &[],
            &settings,
        );
        ctx.queue.submit([encoder.finish()]);
        renderer.after_submit();
        renderer.read_color()
    };
    let fresh = draw(&CartoonGpu::upload(&ctx, &plan, &moved, &colors).unwrap());
    let shifted = CartoonGpu::upload(&ctx, &plan, &start, &colors).unwrap();
    shifted.set_frame(&ctx, &moved);
    let followed = draw(&shifted);
    let differ = fresh
        .chunks_exact(4)
        .zip(followed.chunks_exact(4))
        .filter(|(a, b)| a.iter().zip(*b).any(|(x, y)| x.abs_diff(*y) > 8))
        .count();
    let drawn = count_non_background(&fresh, background_of(&fresh));
    assert!(drawn > 2_000, "only {drawn} pixels drawn");
    assert!(differ * 200 <= drawn, "{differ} of {drawn} pixels differ");
}

/// Green: like `PINHOLE_BG`, never a lit surface colour.
const GREEN_BG: wgpu::Color = wgpu::Color {
    r: 0.0,
    g: 1.0,
    b: 0.0,
    a: 1.0,
};

/// How many pixels are surface in both frames (drawn over backgrounds
/// `a_bg` and `b_bg`), and the largest channel difference among them.
fn surface_difference(a: &[u8], a_bg: [u8; 3], b: &[u8], b_bg: [u8; 3]) -> (usize, u8) {
    a.chunks(4)
        .zip(b.chunks(4))
        .filter(|(p, q)| p[..3] != a_bg && q[..3] != b_bg)
        .fold((0, 0), |(n, max), (p, q)| {
            let d = (0..3).map(|c| p[c].abs_diff(q[c])).max().unwrap();
            (n + 1, max.max(d))
        })
}

fn render_items(
    ctx: &GpuContext,
    renderer: &mut Renderer,
    camera: &Camera,
    items: &[DrawItem],
    cartoons: &[CartoonItem],
    settings: &RenderSettings,
) -> Vec<u8> {
    let mut encoder = ctx.device.create_command_encoder(&Default::default());
    renderer.render_all(&mut encoder, camera, items, cartoons, &[], &[], settings);
    ctx.queue.submit([encoder.finish()]);
    renderer.after_submit();
    renderer.read_color()
}

/// 1CRN's cartoon in one flat colour, framed.
fn crambin_cartoon(ctx: &GpuContext, colour: u32) -> (CartoonGpu, Camera) {
    let path =
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/small/1CRN.cif");
    let structure = vv_io::load(path).unwrap();
    let positions = structure.frame(0).positions().to_vec();
    let codes = vv_core::dssp::assign(&structure.topology, &positions);
    let plan = vv_core::cartoon::plan(&structure.topology, &positions, &codes);
    let frame = plan.frame(&positions);
    let colors = vec![colour; structure.atom_count()];
    let gpu = CartoonGpu::upload(ctx, &plan, &frame, &colors).expect("upload cartoon");
    let (center, radius) = vv_core::CoordSet::new(positions).bounding_sphere().unwrap();
    (gpu, Camera::framing(center, radius))
}

/// Only the depth cue ties a surface to the background: in every style,
/// with it off, spheres and a cartoon shade identically over any
/// background (light, reflection and occlusion alike).
#[test]
fn shading_does_not_depend_on_the_background() {
    let Some(ctx) = context() else { return };
    let (w, h) = (128, 128);
    let structure = protein_like(&SynthParams::new(3_000));
    let gpu = GpuStructure::upload(&ctx, &structure, None, ColorScheme::Element);
    let atoms_camera = Camera::framing(gpu.center, gpu.radius);
    let (cartoon, cartoon_camera) = crambin_cartoon(&ctx, vv_render::color::rgba(200, 80, 80));
    let mut renderer = Renderer::new(ctx.clone(), w, h);
    let bindings = renderer.bind(&gpu);
    let cartoon_bindings = renderer.bind_cartoon(&cartoon);
    for style in StylePreset::ALL {
        let frame = |renderer: &mut Renderer, background, atoms: bool| {
            let settings = RenderSettings {
                occlusion_culling: false,
                background,
                lighting: style.lighting(),
                outline: style.outline(),
                ao: style.ao(),
                depth_cue: 0.0,
                ..Default::default()
            };
            let item = DrawItem {
                structure: &gpu,
                bindings: &bindings,
                representation: Representation::default(),
                sizes: AtomSizes::of(Representation::default()),
                material: style.material(),
            };
            let tube = CartoonItem {
                mesh: CartoonMesh::Ribbon(&cartoon),
                bindings: &cartoon_bindings,
                material: style.material(),
            };
            if atoms {
                render_items(&ctx, renderer, &atoms_camera, &[item], &[], &settings)
            } else {
                render_items(&ctx, renderer, &cartoon_camera, &[], &[tube], &settings)
            }
        };
        for atoms in [true, false] {
            let a = frame(&mut renderer, PINHOLE_BG, atoms);
            let b = frame(&mut renderer, GREEN_BG, atoms);
            let (surface, max) = surface_difference(&a, [255, 0, 255], &b, [0, 255, 0]);
            assert!(surface > 1_000, "{style:?}: only {surface} surface pixels");
            assert!(
                max <= 1,
                "{style:?}, atoms {atoms}: a surface pixel moved by {max}"
            );
        }
    }
}

/// The colour a flat cap facing an orthographic viewer is drawn in, for
/// a palette colour: `shade` with the cap's normal and no specular,
/// sRGB-encoded as the target stores it.
fn cap_color(rgba: u32, settings: &RenderSettings) -> [u8; 3] {
    let base = [rgba, rgba >> 8, rgba >> 16].map(|c| (c & 0xff) as f32 / 255.0);
    let material = vv_render::Material {
        specular: 0.0,
        ..settings.material
    };
    let lit = vv_render::style::shade(
        base,
        [0.0, 0.0, 1.0],
        [0.0, 0.0, -1.0],
        &settings.lighting,
        &material,
    );
    lit.map(|c| {
        let c = c.clamp(0.0, 1.0);
        let s = if c <= 0.0031308 {
            c * 12.92
        } else {
            1.055 * c.powf(1.0 / 2.4) - 0.055
        };
        (s * 255.0).round() as u8
    })
}

/// Pixels within 2 levels of one of the cap colours.
fn cap_pixels(pixels: &[u8], caps: &[[u8; 3]]) -> usize {
    pixels
        .chunks(4)
        .filter(|p| {
            caps.iter()
                .any(|c| (0..3).all(|k| p[k].abs_diff(c[k]) <= 2))
        })
        .count()
}

/// Dollying into a structure cuts it flat through the middle, and the
/// cut is solid caps: spheres and a cartoon show their cross-section
/// there rather than a hole or their own insides.
#[test]
fn dolly_caps_what_it_cuts() {
    let Some(ctx) = context() else { return };
    let (w, h) = (128, 128);
    let settings = RenderSettings {
        occlusion_culling: false,
        background: PINHOLE_BG,
        ao: 0.0,
        depth_cue: 0.0,
        ..Default::default()
    };
    let mut renderer = Renderer::new(ctx.clone(), w, h);

    let structure = protein_like(&SynthParams::new(3_000));
    let gpu = GpuStructure::upload(&ctx, &structure, None, ColorScheme::Element);
    let bindings = renderer.bind(&gpu);
    let mut caps: Vec<[u8; 3]> = vv_render::colors_for(ColorScheme::Element, &structure.topology)
        .into_iter()
        .map(|c| cap_color(c, &settings))
        .collect();
    caps.sort();
    caps.dedup();
    let item = || DrawItem {
        structure: &gpu,
        bindings: &bindings,
        representation: Representation::default(),
        sizes: AtomSizes::of(Representation::default()),
        material: settings.material,
    };
    let mut camera = Camera::framing(gpu.center, gpu.radius);
    let whole = render_items(&ctx, &mut renderer, &camera, &[item()], &[], &settings);
    camera.dolly = gpu.radius;
    let cut = render_items(&ctx, &mut renderer, &camera, &[item()], &[], &settings);
    let i = ((h / 2 * w + w / 2) * 4) as usize;
    assert_ne!(
        cut[i..i + 3],
        [255, 0, 255],
        "the cut centre is capped, not a hole"
    );
    assert!(renderer.pick(w / 2, h / 2).is_some(), "a cap is pickable");
    let (before, after) = (cap_pixels(&whole, &caps), cap_pixels(&cut, &caps));
    // `LightingPreset::Default`'s ambient floor lifts more of each
    // sphere's own shaded surface into the cap colours' range even before
    // dollying in, raising `before` some, so the ratio is stark but
    // modest.
    assert!(
        after > 2 * before + 1_000,
        "spheres: {before} -> {after} cap pixels"
    );

    let colour = vv_render::color::rgba(200, 80, 80);
    let (cartoon, mut camera) = crambin_cartoon(&ctx, colour);
    camera.zoom(0.5);
    let cartoon_bindings = renderer.bind_cartoon(&cartoon);
    let tube = || CartoonItem {
        mesh: CartoonMesh::Ribbon(&cartoon),
        bindings: &cartoon_bindings,
        material: settings.material,
    };
    let caps = [cap_color(colour, &settings)];
    let whole = render_items(&ctx, &mut renderer, &camera, &[], &[tube()], &settings);
    camera.dolly = camera.scene_radius;
    let cut = render_items(&ctx, &mut renderer, &camera, &[], &[tube()], &settings);
    let (before, after) = (cap_pixels(&whole, &caps), cap_pixels(&cut, &caps));
    assert!(
        after > before + 150,
        "cartoon: {before} -> {after} cap pixels"
    );
}

/// In perspective the eye itself moves: dollied to an atom's centre, the
/// view straight ahead is that atom's cap, for van der Waals and for the
/// larger SAS spheres alike.
#[test]
fn perspective_dolly_into_an_atom_shows_its_cap() {
    let Some(ctx) = context() else { return };
    let (w, h) = (64, 64);
    let settings = RenderSettings {
        occlusion_culling: false,
        background: PINHOLE_BG,
        ao: 0.0,
        depth_cue: 0.0,
        ..Default::default()
    };
    let structure = protein_like(&SynthParams::new(500));
    let gpu = GpuStructure::upload(&ctx, &structure, None, ColorScheme::Element);
    let colour = vv_render::colors_for(ColorScheme::Element, &structure.topology)[0];
    let mut renderer = Renderer::new(ctx.clone(), w, h);
    let bindings = renderer.bind(&gpu);
    let mut camera = Camera::framing(gpu.center, gpu.radius);
    camera.projection = Projection::Perspective;
    camera.target = structure.frame(0).positions()[0];
    camera.dolly = camera.distance;
    for representation in [Representation::Spacefill, Representation::Sas] {
        let item = DrawItem {
            structure: &gpu,
            bindings: &bindings,
            representation,
            sizes: AtomSizes::of(representation),
            material: settings.material,
        };
        let pixels = render_items(&ctx, &mut renderer, &camera, &[item], &[], &settings);
        let i = ((h / 2 * w + w / 2) * 4) as usize;
        let cap = cap_color(colour, &settings);
        assert!(
            (0..3).all(|k| pixels[i + k].abs_diff(cap[k]) <= 2),
            "{representation:?}: centre {:?}, cap {cap:?}",
            &pixels[i..i + 3]
        );
        assert_eq!(
            renderer.pick(w / 2, h / 2),
            Some(Pick::Atom { item: 0, atom: 0 }),
            "{representation:?}"
        );
    }
}

/// Mean depth-cue fog fraction over the structure's own pixels: how far the
/// cued frame sits from the plain one toward the white background.
fn mean_fog_fraction(
    ctx: &GpuContext,
    renderer: &mut Renderer,
    camera: &Camera,
    item: &DrawItem,
) -> f64 {
    let render = |renderer: &mut Renderer, depth_cue: f32| {
        let settings = RenderSettings {
            occlusion_culling: false,
            background: wgpu::Color::WHITE,
            depth_cue,
            ..Default::default()
        };
        render_items(
            ctx,
            renderer,
            camera,
            std::slice::from_ref(item),
            &[],
            &settings,
        )
    };
    let plain = render(renderer, 0.0);
    let cued = render(renderer, 1.0);
    let (mut sum, mut count) = (0.0, 0usize);
    for (p, c) in plain.chunks(4).zip(cued.chunks(4)) {
        let room = 255.0 - p[0] as f64;
        if room > 8.0 {
            sum += (c[0] as f64 - p[0] as f64) / room;
            count += 1;
        }
    }
    sum / count.max(1) as f64
}

/// Depth cue is anchored to the structure's front face: with the camera
/// pushed in to the surface the visible structure is clearer than in the
/// framed view, not equally or more fogged.
#[test]
fn depth_cue_clears_as_the_camera_closes_on_the_structure() {
    let Some(ctx) = context() else { return };
    let structure = protein_like(&SynthParams::new(3_000));
    let colors = vec![vv_render::color::rgba(60, 60, 60); structure.atom_count()];
    let gpu = GpuStructure::upload_colored(&ctx, &structure, None, &colors);
    let mut renderer = Renderer::new(ctx.clone(), 128, 128);
    let bindings = renderer.bind(&gpu);
    let item = DrawItem {
        structure: &gpu,
        bindings: &bindings,
        representation: Representation::default(),
        sizes: AtomSizes::of(Representation::default()),
        material: StylePreset::default().material(),
    };
    let framed = Camera::framing(gpu.center, gpu.radius);
    let mut close = framed.clone();
    close.zoom(0.35);
    let far_fog = mean_fog_fraction(&ctx, &mut renderer, &framed, &item);
    let near_fog = mean_fog_fraction(&ctx, &mut renderer, &close, &item);
    eprintln!("mean fog fraction: framed {far_fog:.3}, close {near_fog:.3}");
    assert!(far_fog > 0.05, "the framed view is visibly cued: {far_fog}");
    assert!(
        near_fog < far_fog,
        "closer must be clearer: {near_fog} vs {far_fog}"
    );
}

/// A triangle soup drawn as `MeshDisplay::Lines` covers only a band of
/// pixels along its edges, in the colour pass and the pick pass alike (a
/// click inside the triangle falls through), where drawn solid it covers
/// and picks its whole face.
#[test]
fn mesh_lines_draw_and_pick_only_along_the_edges() {
    use vv_render::{MeshDisplay, Soup};
    let Some(ctx) = context() else { return };
    let corners = [
        glam::Vec3::new(-10.0, -8.0, 0.0),
        glam::Vec3::new(10.0, -8.0, 0.0),
        glam::Vec3::new(0.0, 10.0, 0.0),
    ];
    let normals = [glam::Vec3::Z; 3];
    let colors = [vv_render::color::rgba(220, 60, 60); 3];
    let (w, h) = (256, 256);
    let camera = Camera::framing(glam::Vec3::ZERO, 12.0);
    let settings = RenderSettings {
        occlusion_culling: false,
        ..Default::default()
    };
    let drawn_as = |display: MeshDisplay| {
        let mut renderer = Renderer::new(ctx.clone(), w, h);
        let gpu = vv_render::GlycanGpu::upload_soup(
            &ctx,
            &Soup {
                positions: &corners,
                normals: &normals,
                colors: &colors,
                source: std::sync::Arc::new(vec![7, 7, 7]),
            },
            display,
        )
        .unwrap();
        let bindings = renderer.bind_glycan(&gpu);
        let mut encoder = ctx.device.create_command_encoder(&Default::default());
        renderer.render_all(
            &mut encoder,
            &camera,
            &[],
            &[CartoonItem {
                mesh: CartoonMesh::Glycan(&gpu),
                bindings: &bindings,
                material: vv_render::Material::default(),
            }],
            &[],
            &[],
            &settings,
        );
        ctx.queue.submit([encoder.finish()]);
        renderer.after_submit();
        let pixels = renderer.read_color();
        let bg = background_of(&pixels);
        (
            count_non_background(&pixels, bg),
            renderer.pick(w / 2, h / 2),
        )
    };
    let (solid, solid_pick) = drawn_as(MeshDisplay::Solid);
    let (edges, edge_pick) = drawn_as(MeshDisplay::Lines(2.0));
    assert!(solid > 3_000, "{solid}");
    assert!(
        edges > 100 && edges * 4 < solid,
        "edges {edges} px against {solid} solid"
    );
    assert!(
        matches!(solid_pick, Some(Pick::Atom { .. })),
        "{solid_pick:?}"
    );
    assert_eq!(edge_pick, None, "the face of a wire triangle is empty");
}
