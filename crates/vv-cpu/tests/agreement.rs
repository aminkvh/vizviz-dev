//! The one test that matters most for a second render backend: does it
//! agree with the first? Renders the same synthetic structure and camera
//! through both `vv-cpu` and `vv-render`, and checks they agree on which
//! pixels are background vs. drawn. This is the test that would catch a
//! ray-generation sign error, a flipped V coordinate, a wrong FOV
//! interpretation, or a radius mismatch -- the class of bug a CPU-only
//! "something is drawn" test cannot see (see the vv-render FXAA/SSAA
//! commits this session for the same reasoning applied to those passes).
//!
//! Skips (passes) when the machine has no usable GPU adapter, same as
//! vv-render's own headless tests.

use vv_io::synth::{protein_like, SynthParams};
use vv_render::{
    Camera, ColorScheme, GpuContext, GpuStructure, PatchSurface, PatchSurfaceItem, RenderSettings,
    Renderer, SkinSurfaceGpu, StylePreset,
};

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

const MAGENTA: vv_render::wgpu::Color = vv_render::wgpu::Color {
    r: 1.0,
    g: 0.0,
    b: 1.0,
    a: 1.0,
};

fn is_background(pixels: &[u8], bg: [u8; 3], i: usize) -> bool {
    let p = &pixels[i * 4..i * 4 + 3];
    (p[0] as i32 - bg[0] as i32).abs() <= 8
        && (p[1] as i32 - bg[1] as i32).abs() <= 8
        && (p[2] as i32 - bg[2] as i32).abs() <= 8
}

/// `vv-render`'s color targets (`color`, `outlined`, `fxaa_color` --
/// `COLOR_FORMAT` in `renderer.rs`) are `Rgba8UnormSrgb`: wgpu treats a
/// fragment shader's `shade()` output as linear light and gamma-encodes
/// it on write, and `read_color`'s raw byte copy returns those encoded
/// bytes as-is. `vv_cpu::renderer::shade()` computes the identical
/// formula but writes its result straight to bytes with no such encode
/// -- neither backend linearizes the 8-bit element/chain color table
/// first, so `shade()`'s own output was never meant to be "linear light"
/// in the physical sense; the GPU's sRGB target just applies an extra
/// step CPU's software rasterizer has no equivalent of. Confirmed by
/// direct measurement (a single face-on-lit atom): GPU's byte matches
/// sRGB-encoding CPU's byte far more closely than it matches CPU's raw
/// byte, and the mismatch is uniform across the whole surface (mean
/// delta ~51/255 on every backend-agreed hit pixel, not just at
/// silhouettes), never in the reverse direction (CPU never calls a pixel
/// background where GPU calls it a hit -- a one-sided brightening does
/// exactly that).
///
/// This makes CPU and GPU bytes not directly comparable for a pixel a
/// backend actually shaded, only for background pixels (copied verbatim
/// by both, never gamma-mapped, and already observed to agree byte for
/// byte). Applying this same encode to every CPU pixel that isn't
/// exactly the background color it was asked to fill (a shaded pixel
/// that happens to land on those exact bytes is skipped too, aliased
/// with real background -- harmless at the 99%+ agreement this reaches,
/// but worth knowing if this helper is ever reused somewhere tighter)
/// makes the comparison apples to apples without changing what either
/// renderer actually draws -- this is a test-fairness fix, not a rendering fix:
/// changing the app's actual on-screen colors or changing vv-cpu's
/// headless/Python renders to match is a separate, user-visible decision
/// this function does not make.
fn as_gpu_would_show(cpu_pixels: &[u8], cpu_background: [u8; 4]) -> Vec<u8> {
    fn srgb_encode_u8(b: u8) -> u8 {
        let c = b as f32 / 255.0;
        let encoded = if c <= 0.003_130_8 {
            12.92 * c
        } else {
            1.055 * c.powf(1.0 / 2.4) - 0.055
        };
        (encoded.clamp(0.0, 1.0) * 255.0).round() as u8
    }
    let mut out = cpu_pixels.to_vec();
    for px in out.chunks_mut(4) {
        if px == cpu_background {
            continue;
        }
        px[0] = srgb_encode_u8(px[0]);
        px[1] = srgb_encode_u8(px[1]);
        px[2] = srgb_encode_u8(px[2]);
    }
    out
}

#[test]
fn cpu_and_gpu_agree_on_background_vs_drawn() {
    let Some(ctx) = context() else { return };
    let (w, h) = (128u32, 128u32);
    let structure = protein_like(&SynthParams::new(5_000));
    let style = StylePreset::default();

    let gpu = GpuStructure::upload(&ctx, &structure, None, ColorScheme::Element);
    let mut camera = Camera::framing(gpu.center, gpu.radius);
    // vv-cpu casts perspective rays only; compare like with like.
    camera.projection = vv_render::Projection::Perspective;
    camera.distance *= 1.1; // real background at the corners, not just tight framing

    let mut renderer = Renderer::new(ctx.clone(), w, h);
    let bindings = renderer.bind(&gpu);
    let settings = RenderSettings {
        occlusion_culling: false, // vv-cpu has no occlusion culling to compare against
        lighting: style.lighting(),
        material: style.material(),
        background: style.background(),
        // Nor ambient occlusion or depth cue: this compares geometry.
        ao: 0.0,
        depth_cue: 0.0,
        ..Default::default()
    };
    let mut encoder = ctx.device.create_command_encoder(&Default::default());
    renderer.render(&mut encoder, &camera, &gpu, &bindings, &settings);
    ctx.queue.submit([encoder.finish()]);
    renderer.after_submit();
    let gpu_pixels = renderer.read_color();

    let cpu_background = {
        let c = style.background();
        [
            (c.r * 255.0).round() as u8,
            (c.g * 255.0).round() as u8,
            (c.b * 255.0).round() as u8,
            255,
        ]
    };
    let cpu_pixels = vv_cpu::render_spacefill(
        &structure,
        ColorScheme::Element,
        &camera,
        w,
        h,
        style,
        cpu_background,
    );

    // Each backend's own top-left corner as its background reference,
    // rather than assuming the two clear colors round-trip identically
    // through wgpu's sRGB target -- classification agreement is what
    // matters here, not byte-identical background color.
    let gpu_bg = [gpu_pixels[0], gpu_pixels[1], gpu_pixels[2]];
    let cpu_bg = [cpu_pixels[0], cpu_pixels[1], cpu_pixels[2]];

    let total = (w * h) as usize;
    let mut agree = 0usize;
    let (mut checked_hit, mut checked_miss) = (false, false);
    for i in 0..total {
        let gpu_is_bg = is_background(&gpu_pixels, gpu_bg, i);
        let cpu_is_bg = is_background(&cpu_pixels, cpu_bg, i);
        if gpu_is_bg == cpu_is_bg {
            agree += 1;
        }
        if gpu_is_bg {
            checked_miss = true;
        } else {
            checked_hit = true;
        }
    }
    assert!(
        checked_hit && checked_miss,
        "test should exercise both background and drawn pixels"
    );
    assert!(
        agree * 100 >= total * 95,
        "CPU and GPU should mostly agree on background vs. drawn ({agree}/{total} agreed) -- \
         a ray-generation or radius bug would show up as widespread disagreement, not just \
         soft edges"
    );
}

/// Focused repro for a real disagreement found via the test below at
/// zoom 0.35: patch 3091 (a Triangle patch, atoms [30,310,311] of
/// 1CRN.pdb) is hit by a CPU ray solidly inside its mixed cell
/// (membership margins of order 1, nowhere near `MEMBERSHIP_EPS`) at
/// several consecutive pixels, but GPU shows background there even
/// though the billboard's own angular coverage comfortably contains
/// those ray directions. Aims a camera straight at that patch's bound
/// sphere from a modest distance and dumps both renders for direct
/// visual comparison, isolated from the rest of the structure's
/// depth-compositing.
#[test]
#[ignore]
fn diagnose_triangle_patch_gap() {
    let Some(ctx) = context() else { return };
    let path =
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/small/1CRN.pdb");
    let structure = vv_io::load(&path).unwrap();
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
    let scene = vv_cpu::skin_surface::SkinSurfaceScene::new(
        positions.clone(),
        weights.clone(),
        colors.clone(),
        vv_core::skin_surface::DEFAULT_SHRINK,
    );
    let idx = scene
        .complex()
        .patches
        .iter()
        .position(|p| {
            p.kind == vv_core::skin_surface::PatchKind::Triangle && {
                let mut m = p.members().to_vec();
                m.sort_unstable();
                m == [30, 310, 311]
            }
        })
        .expect("patch 3091 should still exist");
    let p = &scene.complex().patches[idx];
    println!(
        "patch {idx}: bound_center={:?} bound_radius={} weight={} axis={:?}",
        p.bound_center, p.bound_radius, p.weight, p.axis
    );

    let (w, h) = (128u32, 128u32);
    let style = StylePreset::default();
    let mut camera = Camera::framing(p.bound_center, p.bound_radius * 2.5);
    // vv-cpu casts perspective rays only; compare like with like.
    camera.projection = vv_render::Projection::Perspective;

    let gpu = SkinSurfaceGpu::upload(
        &ctx,
        &positions,
        &weights,
        &colors,
        vv_core::skin_surface::DEFAULT_SHRINK,
    );
    let mut renderer = Renderer::new(ctx.clone(), w, h);
    let bindings = renderer.bind_skin_surface(&gpu);
    let settings = RenderSettings {
        occlusion_culling: false,
        lighting: style.lighting(),
        material: style.material(),
        background: style.background(),
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
            material: style.material(),
        }],
        &settings,
    );
    ctx.queue.submit([encoder.finish()]);
    renderer.after_submit();
    let gpu_pixels = renderer.read_color();

    let cpu_background = {
        let c = style.background();
        [
            (c.r * 255.0).round() as u8,
            (c.g * 255.0).round() as u8,
            (c.b * 255.0).round() as u8,
            255,
        ]
    };
    let cpu_pixels = scene.render(&camera, w, h, style, cpu_background);

    let dump = |name: &str, pixels: &[u8]| {
        let out = std::env::temp_dir().join(name);
        let file = std::fs::File::create(&out).unwrap();
        let mut encoder = png::Encoder::new(std::io::BufWriter::new(file), w, h);
        encoder.set_color(png::ColorType::Rgba);
        encoder.set_depth(png::BitDepth::Eight);
        encoder
            .write_header()
            .unwrap()
            .write_image_data(pixels)
            .unwrap();
        println!("wrote {}", out.display());
    };
    dump("triangle_gap_gpu.png", &gpu_pixels);
    dump("triangle_gap_cpu.png", &cpu_pixels);
}

/// The skin surface's own version of the test above, on a real structure
/// zoomed in close: the CPU ray cast (`vv_cpu::skin_surface`) and the
/// shader's port of the same math must agree on where the surface is.
/// Perspective only (the CPU backend casts no orthographic rays). A
/// magenta background keeps dark shading from passing for background;
/// the CPU image is sRGB-encoded first (`as_gpu_would_show`).
#[test]
fn cpu_and_gpu_agree_on_skin_surface_background_vs_drawn() {
    let Some(ctx) = context() else { return };
    let (w, h) = (200u32, 150u32);
    let path =
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/small/1CRN.pdb");
    let structure = vv_io::load(&path).unwrap();
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
    let style = StylePreset::default();
    let (center, radius) = structure.frame(0).bounding_sphere().unwrap();

    let mut camera = Camera::framing(center, radius);
    // vv-cpu casts perspective rays only; compare like with like.
    camera.projection = vv_render::Projection::Perspective;
    // Closer, and the surface fills the frame: no background left to compare.
    camera.zoom(0.5);

    let gpu = SkinSurfaceGpu::upload(
        &ctx,
        &positions,
        &weights,
        &colors,
        vv_core::skin_surface::DEFAULT_SHRINK,
    );
    let mut renderer = Renderer::new(ctx.clone(), w, h);
    let bindings = renderer.bind_skin_surface(&gpu);
    let settings = RenderSettings {
        occlusion_culling: false,
        lighting: style.lighting(),
        material: style.material(),
        // Magenta: never a lit surface colour, so dark shading can't pass
        // for background.
        background: MAGENTA,
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
            material: style.material(),
        }],
        &settings,
    );
    ctx.queue.submit([encoder.finish()]);
    renderer.after_submit();
    let gpu_pixels = renderer.read_color();

    let cpu_background = {
        let c = MAGENTA;
        [
            (c.r * 255.0).round() as u8,
            (c.g * 255.0).round() as u8,
            (c.b * 255.0).round() as u8,
            255,
        ]
    };
    let cpu_pixels_raw = vv_cpu::skin_surface::SkinSurfaceScene::new(
        positions,
        weights,
        colors,
        vv_core::skin_surface::DEFAULT_SHRINK,
    )
    .render(&camera, w, h, style, cpu_background);
    let cpu_pixels = as_gpu_would_show(&cpu_pixels_raw, cpu_background);

    let (gpu_bg, cpu_bg) = ([255, 0, 255], [255, 0, 255]);

    let total = (w * h) as usize;
    let mut agree = 0usize;
    let (mut checked_hit, mut checked_miss) = (false, false);
    for i in 0..total {
        let gpu_is_bg = is_background(&gpu_pixels, gpu_bg, i);
        let cpu_is_bg = is_background(&cpu_pixels, cpu_bg, i);
        if gpu_is_bg == cpu_is_bg {
            agree += 1;
        }
        if gpu_is_bg {
            checked_miss = true;
        } else {
            checked_hit = true;
        }
    }
    assert!(
        checked_hit && checked_miss,
        "test should exercise both background and drawn pixels"
    );
    assert!(
        agree * 1000 >= total * 998,
        "CPU and GPU skin surfaces disagree on coverage: {agree}/{total} pixels agree"
    );
}
