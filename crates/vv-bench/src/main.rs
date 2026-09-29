//! Headless frame-time harness.
//!
//! Renders a synthetic structure through scripted camera orbits at several
//! zoom levels and reports frame-time percentiles as JSON. No window, so it
//! runs anywhere a GPU adapter exists, and every frame is waited on, so the
//! CPU number is an upper bound on real frame time.

use std::path::PathBuf;
use std::time::Instant;

use serde::Serialize;
use vv_io::synth::{protein_like, SynthParams};
use vv_render::{
    colors_for, AdaptiveLod, AtomSizes, Camera, CartoonGpu, CartoonItem, ColorScheme, DrawItem,
    GaussianSurfaceGpu, GaussianSurfaceItem, GpuContext, GpuStructure, PatchSurface,
    PatchSurfaceItem, RenderSettings, Renderer, Representation, SesGpu, SesLayout, SkinSurfaceGpu,
};

const USAGE: &str = "\
usage: vv-bench [options]
  --atoms N          synthetic atom count (default 10000000)
  --file PATH        benchmark a structure file instead of a synthetic one
  --rep NAME         spacefill, ballstick, cartoon, gaussian, skin, ses, sesvolume (default spacefill)
  --limit N          keep only the first N atoms of the structure (scaling tests)
  --frames N         frames per preset after warm-up (default 120)
  --size WxH         render target size (default 1920x1080)
  --preset a,b,c     subset of overview,mid,close,inside (default all)
  --threshold PX     fixed quad/point threshold in pixels (default 1.0)
  --adaptive         let the LOD controller move the threshold
  --no-occlusion     disable depth-pyramid occlusion culling
  --fast-glass       cull glass atoms against last frame's own nearest-
                     glass depth (RenderSettings::fast_glass); off by
                     default like the app's still view, on while it's
                     interacting -- see docs/RENDERING.md, Transparency
  --clip             cut the structure with a plane through its centre
  --projection P     perspective (default, as in every committed result) or orthographic
  --look             include the post passes the app runs (AO, depth cue)
  --shadows          also screen-space shadows (with Full lighting)
  --material NAME    material for the drawn structure (default opaque; glass1 etc. are transparent)
  --adapter NAME     use the first GPU whose name contains NAME
  --png DIR          write the last frame of each preset as PNG
  --trace SAMPLES    path-trace one image per preset at --size (atom reps,
                     ses, skin or gaussian) with SAMPLES per pixel, timed;
                     writes PNGs to --png
  --json FILE        write results to FILE instead of stdout";

struct Args {
    atoms: usize,
    file: Option<PathBuf>,
    rep: BenchRep,
    limit: Option<usize>,
    frames: usize,
    width: u32,
    height: u32,
    presets: Vec<Preset>,
    threshold: f32,
    adaptive: bool,
    occlusion: bool,
    fast_glass: bool,
    look: bool,
    shadows: bool,
    material: vv_render::MaterialPreset,
    orthographic: bool,
    adapter: Option<String>,
    trace: Option<u32>,
    png: Option<PathBuf>,
    json: Option<PathBuf>,
    /// A clip plane through the structure's centre, view-facing at each
    /// preset's fixed orbit start (not re-facing every frame, unlike the
    /// app's `clip view`): isolates the cut's own cost from camera motion.
    clip: bool,
}

/// What gets drawn: an impostor representation, or one of the derived
/// geometries the app builds in `vv-app/src/gpu_cache.rs` (built here the
/// same way, so build time is measured too).
#[derive(Clone, Copy, Debug, PartialEq)]
enum BenchRep {
    Atoms(Representation),
    Cartoon,
    Gaussian,
    Skin,
    Ses,
    SesVolume,
}

#[derive(Clone, Copy, Debug, PartialEq, Serialize)]
#[serde(rename_all = "lowercase")]
enum Preset {
    /// Whole structure framed in view.
    Overview,
    /// Structure overflows the frame; camera well outside the surface.
    Mid,
    /// Camera just outside the surface, atoms tens of pixels wide.
    Close,
    /// Camera inside the structure: the overdraw worst case.
    Inside,
}

impl Preset {
    fn all() -> Vec<Preset> {
        vec![Preset::Overview, Preset::Mid, Preset::Close, Preset::Inside]
    }

    fn parse(s: &str) -> Option<Preset> {
        match s {
            "overview" => Some(Preset::Overview),
            "mid" => Some(Preset::Mid),
            "close" => Some(Preset::Close),
            "inside" => Some(Preset::Inside),
            _ => None,
        }
    }

    /// Camera distance from the centre as a multiple of the bounding radius.
    /// The framing distance is ~2.74 radii, so anything below 1.0 is inside.
    fn distance_radii(self) -> f32 {
        match self {
            Preset::Overview => 2.74,
            Preset::Mid => 1.6,
            Preset::Close => 1.15,
            Preset::Inside => 0.5,
        }
    }
}

#[derive(Serialize)]
struct Percentiles {
    p50: f32,
    p95: f32,
    p99: f32,
}

#[derive(Serialize)]
struct PresetResult {
    preset: Preset,
    frames: usize,
    threshold_px: f32,
    cpu_ms: Percentiles,
    gpu_ms: Option<Percentiles>,
    cull_ms_p50: Option<f32>,
    draw_ms_p50: Option<f32>,
    hiz_ms_p50: Option<f32>,
    post_ms_p50: Option<f32>,
    visible_quads: u32,
    visible_points: u32,
    visible_bonds: u32,
}

#[derive(Serialize)]
struct Report {
    adapter: String,
    source: String,
    atoms: usize,
    bonds: usize,
    representation: String,
    /// CPU + upload time for derived geometry (cartoon, surfaces).
    build_ms: Option<f32>,
    width: u32,
    height: u32,
    adaptive: bool,
    occlusion: bool,
    results: Vec<PresetResult>,
}

fn parse_args() -> Result<Args, String> {
    let mut args = Args {
        atoms: 10_000_000,
        file: None,
        rep: BenchRep::Atoms(Representation::Spacefill),
        limit: None,
        frames: 120,
        width: 1920,
        height: 1080,
        presets: Preset::all(),
        threshold: 1.0,
        adaptive: false,
        occlusion: true,
        fast_glass: false,
        look: false,
        shadows: false,
        material: vv_render::MaterialPreset::Opaque,
        orthographic: false,
        adapter: None,
        trace: None,
        png: None,
        json: None,
        clip: false,
    };
    let mut it = std::env::args().skip(1);
    while let Some(flag) = it.next() {
        let mut value = || it.next().ok_or_else(|| format!("{flag} needs a value"));
        match flag.as_str() {
            "--atoms" => args.atoms = value()?.parse().map_err(|e| format!("--atoms: {e}"))?,
            "--file" => args.file = Some(PathBuf::from(value()?)),
            "--rep" => {
                args.rep = match value()?.as_str() {
                    "spacefill" => BenchRep::Atoms(Representation::Spacefill),
                    "ballstick" => BenchRep::Atoms(Representation::BallAndStick),
                    "cartoon" => BenchRep::Cartoon,
                    "gaussian" | "gaussiansurface" => BenchRep::Gaussian,
                    "skin" | "skinsurface" => BenchRep::Skin,
                    "ses" => BenchRep::Ses,
                    "sesvolume" => BenchRep::SesVolume,
                    other => return Err(format!("--rep: unknown representation `{other}`")),
                }
            }
            "--limit" => args.limit = Some(value()?.parse().map_err(|e| format!("--limit: {e}"))?),
            "--trace" => args.trace = Some(value()?.parse().map_err(|e| format!("--trace: {e}"))?),
            "--frames" => args.frames = value()?.parse().map_err(|e| format!("--frames: {e}"))?,
            "--size" => {
                let v = value()?;
                let (w, h) = v.split_once('x').ok_or("--size expects WxH")?;
                args.width = w.parse().map_err(|e| format!("--size: {e}"))?;
                args.height = h.parse().map_err(|e| format!("--size: {e}"))?;
            }
            "--preset" => {
                args.presets = value()?
                    .split(',')
                    .map(|p| Preset::parse(p).ok_or_else(|| format!("unknown preset `{p}`")))
                    .collect::<Result<_, _>>()?;
            }
            "--threshold" => {
                args.threshold = value()?.parse().map_err(|e| format!("--threshold: {e}"))?
            }
            "--adaptive" => args.adaptive = true,
            "--no-occlusion" => args.occlusion = false,
            "--fast-glass" => args.fast_glass = true,
            "--clip" => args.clip = true,
            "--look" => args.look = true,
            "--shadows" => args.shadows = true,
            "--material" => {
                let name = value()?;
                args.material = vv_render::MaterialPreset::parse(&name)
                    .ok_or(format!("--material: unknown `{name}`"))?;
            }
            "--projection" => {
                args.orthographic = match value()?.as_str() {
                    "perspective" => false,
                    "orthographic" | "ortho" => true,
                    other => return Err(format!("--projection: unknown `{other}`")),
                }
            }
            "--adapter" => args.adapter = Some(value()?),
            "--png" => args.png = Some(PathBuf::from(value()?)),
            "--json" => args.json = Some(PathBuf::from(value()?)),
            "-h" | "--help" => {
                println!("{USAGE}");
                std::process::exit(0);
            }
            other => return Err(format!("unknown option `{other}`\n\n{USAGE}")),
        }
    }
    Ok(args)
}

fn percentiles(samples: &mut [f32]) -> Percentiles {
    samples.sort_by(|a, b| a.total_cmp(b));
    let at = |q: f32| samples[((samples.len() - 1) as f32 * q).round() as usize];
    Percentiles {
        p50: at(0.50),
        p95: at(0.95),
        p99: at(0.99),
    }
}

fn main() {
    let args = match parse_args() {
        Ok(a) => a,
        Err(e) => {
            eprintln!("{e}");
            std::process::exit(2);
        }
    };

    let instance = vv_render::GpuContext::instance();
    let ctx = match GpuContext::with_instance(instance, None, args.adapter.as_deref()) {
        Ok(ctx) => ctx,
        Err(e) => {
            eprintln!("{e}");
            std::process::exit(1);
        }
    };
    eprintln!("adapter: {}", ctx.adapter_name());

    let t = Instant::now();
    let (structure, source) = match &args.file {
        Some(path) => {
            let s = match vv_io::load(path) {
                Ok(s) => s,
                Err(e) => {
                    eprintln!("{e}");
                    std::process::exit(1);
                }
            };
            eprintln!(
                "loaded {} atoms from {} in {:.2?}",
                s.atom_count(),
                path.display(),
                t.elapsed()
            );
            (s, path.display().to_string())
        }
        None => {
            let s = protein_like(&SynthParams::new(args.atoms));
            eprintln!("generated {} atoms in {:.2?}", s.atom_count(), t.elapsed());
            let label = format!("synthetic {} atoms", s.atom_count());
            (s, label)
        }
    };
    let bonds = match args.rep {
        BenchRep::Atoms(Representation::BallAndStick) => {
            let t = Instant::now();
            let b = vv_core::bonds::perceive(&structure.topology, structure.frame(0).positions());
            eprintln!("{} bonds perceived in {:.2?}", b.len(), t.elapsed());
            Some(b)
        }
        BenchRep::Atoms(Representation::Tube) => unreachable!("not parsed"),
        _ => None,
    };
    let t = Instant::now();
    let gpu = GpuStructure::upload(&ctx, &structure, bonds.as_ref(), ColorScheme::Element);
    eprintln!(
        "uploaded in {:.2?} ({} pages, {} bonds)",
        t.elapsed(),
        gpu.pages.len(),
        gpu.bond_count
    );

    let mut renderer = Renderer::new(ctx.clone(), args.width, args.height);
    let bindings = renderer.bind(&gpu);
    if args.look {
        // The app's long-range AO and shadows; atoms stand in for every
        // representation here (the app uses cartoon sections for cartoons).
        let t = Instant::now();
        let radii = vv_render::scene::vdw_radii(&structure.topology.element);
        let spheres: Vec<_> = structure
            .frame(0)
            .positions()
            .iter()
            .zip(&radii)
            .map(|(p, r)| (*p, *r))
            .collect();
        let volume = vv_render::OcclusionVolume::build(&ctx, &spheres).map(std::sync::Arc::new);
        eprintln!("occlusion volume in {:.2?}", t.elapsed());
        renderer.set_occlusion_volume(volume);
    }

    // Derived geometry, built the way `vv-app/src/gpu_cache.rs` builds it.
    let all_colors = colors_for(ColorScheme::Element, &structure.topology);
    let n = args.limit.unwrap_or(usize::MAX).min(structure.atom_count());
    let coords = structure.frame(0);
    let positions = &coords.positions()[..n];
    let elements = &structure.topology.element[..n];
    let colors = &all_colors[..n];
    let t = Instant::now();
    let mut cartoon = None;
    let mut gaussian = None;
    let mut skin = None;
    let mut ses = None;
    match args.rep {
        BenchRep::Atoms(_) => {}
        BenchRep::Cartoon => {
            let coords = structure.frame(0);
            let all_positions = coords.positions();
            let codes = vv_core::cartoon::secondary_structure(
                &structure.topology,
                all_positions,
                structure.frame_count() == 1,
            );
            eprintln!("dssp in {:.2?}", t.elapsed());
            let plan = vv_core::cartoon::plan(&structure.topology, all_positions, &codes);
            let frame = plan.frame(all_positions);
            eprintln!(
                "cartoon planned at {:.2?}: {} sections, {} joins",
                t.elapsed(),
                plan.recipes.len(),
                plan.joins.len()
            );
            let gpu = CartoonGpu::upload(&ctx, &plan, &frame, &all_colors).unwrap_or_else(|e| {
                eprintln!("cartoon: {e}");
                std::process::exit(1);
            });
            let bind = renderer.bind_cartoon(&gpu);
            cartoon = Some((gpu, bind));
        }
        BenchRep::Gaussian => {
            let radii: Vec<f32> = elements.iter().map(|e| e.vdw_radius()).collect();
            let gpu = GaussianSurfaceGpu::upload(
                &ctx,
                positions,
                &radii,
                colors,
                vv_core::gaussian_surface::DEFAULT_BLOB_FACTOR,
            )
            .unwrap_or_else(|e| {
                eprintln!("gaussian surface: {e}");
                std::process::exit(1);
            });
            // Include the GPU bake: this is the per-frame cost when a
            // trajectory frame changes.
            let _ = ctx
                .device
                .poll(vv_render::wgpu::PollType::wait_indefinitely());
            eprintln!(
                "gaussian volume: {:?} voxels at {:.2} A",
                gpu.base_params.dims, gpu.base_params.voxel
            );
            let bind = renderer.bind_gaussian_surface(&gpu);
            gaussian = Some((gpu, bind));
        }
        BenchRep::Skin => {
            let shrink = vv_core::skin_surface::DEFAULT_SHRINK;
            let weights: Vec<f32> = elements
                .iter()
                .map(|e| vv_core::skin_surface::weight_for_radius(e.vdw_radius(), shrink))
                .collect();
            let gpu = SkinSurfaceGpu::upload(&ctx, positions, &weights, colors, shrink);
            eprintln!(
                "skin patches by kind (v, e, t, tet): {:?}",
                gpu.patch_counts
            );
            let bind = renderer.bind_skin_surface(&gpu);
            skin = Some((gpu, bind));
        }
        BenchRep::Ses => {
            let radii: Vec<f32> = elements.iter().map(|e| e.vdw_radius()).collect();
            let built = vv_core::ses::build(positions, &radii, vv_core::ses::WATER_PROBE);
            eprintln!(
                "ses built in {:.2?}: {} convex, {} tori, {} concave",
                t.elapsed(),
                built.convex.len(),
                built.tori.len(),
                built.probes.len()
            );
            let layout = SesLayout::new(&built, positions, &radii);
            let gpu = SesGpu::upload(&ctx, &layout, colors).unwrap_or_else(|_| {
                eprintln!("ses: out of GPU memory");
                std::process::exit(1);
            });
            let bind = renderer.bind_ses_surface(&gpu);
            ses = Some((gpu, bind));
        }
        BenchRep::SesVolume => {
            let radii: Vec<f32> = elements.iter().map(|e| e.vdw_radius()).collect();
            let built = vv_core::ses::build(positions, &radii, vv_core::ses::WATER_PROBE);
            eprintln!("ses built in {:.2?}", t.elapsed());
            let baking = std::time::Instant::now();
            let gpu = GaussianSurfaceGpu::ses(&ctx, &built, positions, &radii, colors)
                .unwrap_or_else(|_| {
                    eprintln!("ses volume: out of GPU memory");
                    std::process::exit(1);
                });
            let _ = ctx
                .device
                .poll(vv_render::wgpu::PollType::wait_indefinitely());
            eprintln!(
                "ses volume: {:?} voxels at {:.2} A, baked in {:.2?}",
                gpu.base_params.dims,
                gpu.base_params.voxel,
                baking.elapsed()
            );
            let bind = renderer.bind_gaussian_surface(&gpu);
            gaussian = Some((gpu, bind));
        }
    }
    let build_ms = (!matches!(args.rep, BenchRep::Atoms(_))).then(|| {
        let ms = t.elapsed().as_secs_f32() * 1e3;
        eprintln!("built {:?} over {n} atoms in {ms:.1} ms", args.rep);
        ms
    });
    let representation = match args.rep {
        BenchRep::Atoms(r) => r,
        _ => Representation::Spacefill,
    };
    if let Some(samples) = args.trace {
        let t = Instant::now();
        let vdw: Vec<f32> = elements.iter().map(|e| e.vdw_radius()).collect();
        let material = vv_render::Material::default();
        let mut scene = vv_render::path_trace::TraceScene::default();
        match args.rep {
            BenchRep::Atoms(rep) => {
                let radii: Vec<f32> = vdw.iter().map(|&r| rep.atom_radius(r)).collect();
                scene.push_atoms(
                    positions,
                    &radii,
                    colors,
                    0..n,
                    bonds
                        .iter()
                        .flat_map(|b| b.pairs.iter().copied())
                        .filter(|p| (p[1] as usize) < n),
                    rep.bond_radius(),
                    material,
                );
            }
            BenchRep::Ses => {
                let built = vv_core::ses::build(positions, &vdw, vv_core::ses::WATER_PROBE);
                scene.push_ses(&built, positions, &vdw, colors, material);
            }
            BenchRep::Skin => {
                let shrink = vv_core::skin_surface::DEFAULT_SHRINK;
                let weights: Vec<f32> = vdw
                    .iter()
                    .map(|&r| vv_core::skin_surface::weight_for_radius(r, shrink))
                    .collect();
                let complex = vv_core::skin_surface::build_complex(positions, &weights, shrink);
                scene.push_skin(&complex, positions, &weights, colors, material);
            }
            BenchRep::Gaussian => {
                let mesh = vv_core::gaussian_mesh::mesh(
                    positions,
                    &vdw,
                    vv_core::gaussian_surface::DEFAULT_BLOB_FACTOR,
                    vv_render::scene::GAUSSIAN_EPSILON,
                    0.5,
                );
                scene.push_mesh(&mesh, colors, material);
            }
            BenchRep::Cartoon | BenchRep::SesVolume => {
                eprintln!("--trace draws atoms, ses, skin and gaussian");
                std::process::exit(2);
            }
        }
        let tracer = vv_render::path_trace::PathTracer::new(&ctx, &scene).unwrap_or_else(|_| {
            eprintln!("trace: out of GPU memory");
            std::process::exit(1);
        });
        eprintln!(
            "trace scene: {} spheres, {} cylinders, {} triangles, built in {:.2?}",
            scene.spheres.len(),
            scene.cylinders.len(),
            scene.triangles.len(),
            t.elapsed()
        );
        // Framed on what is traced (`--limit` may cut the structure).
        let (lo, hi) = positions.iter().fold(
            (
                vv_core::glam::Vec3::splat(f32::MAX),
                vv_core::glam::Vec3::splat(f32::MIN),
            ),
            |(lo, hi), &p| (lo.min(p), hi.max(p)),
        );
        let (center, radius) = ((lo + hi) * 0.5, ((hi - lo).length() * 0.5).max(1.0));
        for preset in &args.presets {
            let mut camera = Camera::framing(center, radius);
            camera.distance = radius * preset.distance_radii();
            camera.projection = if args.orthographic {
                vv_render::Projection::Orthographic
            } else {
                vv_render::Projection::Perspective
            };
            let settings = vv_render::path_trace::TraceSettings {
                width: args.width,
                height: args.height,
                samples,
                lighting: vv_render::Lighting::default(),
                background: [0.02; 3],
                background_top: [0.08; 3],
                transparent: true,
                clip: args.clip.then(|| [0.0, 0.0, 1.0, -center.z]),
                unshadowed: false,
                light_spread: vv_render::path_trace::DEFAULT_LIGHT_SPREAD,
                ao: 1.0,
                direct: 1.0,
            };
            let t = Instant::now();
            let image = tracer.render(&ctx, &camera, &settings, |_| {});
            eprintln!(
                "{preset:?}: traced {}x{} at {samples} spp in {:.2?}",
                args.width,
                args.height,
                t.elapsed()
            );
            if let Some(dir) = &args.png {
                let path = dir.join(format!("trace_{preset:?}.png").to_lowercase());
                let file = std::fs::File::create(&path).expect("create png");
                let mut encoder =
                    png::Encoder::new(std::io::BufWriter::new(file), args.width, args.height);
                encoder.set_color(png::ColorType::Rgba);
                encoder.set_depth(png::BitDepth::Eight);
                encoder
                    .write_header()
                    .and_then(|mut w| w.write_image_data(&image))
                    .expect("write png");
                eprintln!("wrote {}", path.display());
            }
        }
        return;
    }

    let warmup = 10;
    let mut results = Vec::new();

    for preset in &args.presets {
        let mut settings = RenderSettings {
            representation,
            quad_px_threshold: args.threshold,
            occlusion_culling: args.occlusion,
            fast_glass: args.fast_glass,
            // Benchmarks measure the rendering pipeline, not the
            // post-process filter on top of it; FXAA isn't GPU-timestamped
            // (see docs/RENDERING.md), so leaving it on would silently
            // fold into these numbers as unaccounted wall-clock time.
            fxaa: false,
            // Same reason as fxaa: keep committed numbers free of post passes.
            outline: false,
            outline_width: 1,
            ao: if args.look { 1.0 } else { 0.0 },
            depth_cue: if args.look { 0.5 } else { 0.0 },
            shadows: if args.shadows { 1.0 } else { 0.0 },
            lighting: if args.shadows {
                vv_render::LightingPreset::Full.lighting()
            } else {
                vv_render::Lighting::default()
            },
            clip: args.clip.then(|| [0.0, 0.0, 1.0, -gpu.center.z]),
            ..Default::default()
        };
        let mut lod = AdaptiveLod {
            threshold_px: args.threshold,
            ..Default::default()
        };
        let mut base = Camera::framing(gpu.center, gpu.radius);
        // The app defaults to orthographic; committed results are perspective.
        base.projection = if args.orthographic {
            vv_render::Projection::Orthographic
        } else {
            vv_render::Projection::Perspective
        };
        let mut cpu = Vec::with_capacity(args.frames);
        let mut gpu_total = Vec::new();
        let mut cull = Vec::new();
        let mut draw = Vec::new();
        let mut hiz = Vec::new();
        let mut post = Vec::new();
        let total = warmup + args.frames;

        for i in 0..total {
            let phase = i as f32 / args.frames as f32 * std::f32::consts::TAU;
            let mut camera = base.clone();
            camera.distance = gpu.radius * preset.distance_radii();
            camera.orientation = Camera::angles(phase, 0.3 * phase.sin());

            let start = Instant::now();
            let mut encoder = ctx.device.create_command_encoder(&Default::default());
            let items: Vec<DrawItem> = match args.rep {
                BenchRep::Atoms(r) => vec![DrawItem {
                    structure: &gpu,
                    bindings: &bindings,
                    representation: r,
                    sizes: AtomSizes::of(r),
                    material: args.material.into(),
                }],
                _ => Vec::new(),
            };
            let cartoons: Vec<CartoonItem> = cartoon
                .iter()
                .map(|(mesh, bindings)| CartoonItem {
                    mesh: vv_render::CartoonMesh::Ribbon(mesh),
                    bindings,
                    material: args.material.into(),
                })
                .collect();
            let gaussians: Vec<GaussianSurfaceItem> = gaussian
                .iter()
                .map(|(gpu, bindings)| GaussianSurfaceItem {
                    gpu,
                    bindings,
                    material: args.material.into(),
                })
                .collect();
            let skins: Vec<PatchSurfaceItem> = skin
                .iter()
                .map(|(gpu, bindings)| PatchSurface::Skin(gpu, bindings))
                .chain(
                    ses.iter()
                        .map(|(gpu, bindings)| PatchSurface::Ses(gpu, bindings)),
                )
                .map(|surface| PatchSurfaceItem {
                    surface,
                    material: args.material.into(),
                })
                .collect();
            renderer.render_all(
                &mut encoder,
                &camera,
                &items,
                &cartoons,
                &gaussians,
                &skins,
                &settings,
            );
            ctx.queue.submit([encoder.finish()]);
            renderer.after_submit();
            let _ = ctx
                .device
                .poll(vv_render::wgpu::PollType::wait_indefinitely());
            let cpu_ms = start.elapsed().as_secs_f32() * 1e3;

            if let Some(times) = renderer.last_times() {
                let frame_ms = times.total_ms();
                if args.adaptive {
                    settings.quad_px_threshold = lod.update(frame_ms);
                }
                if i >= warmup {
                    gpu_total.push(frame_ms);
                    cull.push(times.cull_ms);
                    draw.push(times.draw_ms);
                    hiz.push(times.hiz_ms);
                    post.push(times.post_ms);
                }
            }
            if i >= warmup {
                cpu.push(cpu_ms);
            }
        }

        let counts = gpu.read_visible_counts(&ctx);
        let visible_quads = counts.iter().map(|c| c.quads).sum();
        let visible_points = counts.iter().map(|c| c.points).sum();
        let visible_bonds = counts.iter().map(|c| c.bonds).sum();

        if let Some(dir) = &args.png {
            std::fs::create_dir_all(dir).expect("create png dir");
            let path = dir.join(format!("{:?}.png", preset).to_lowercase());
            write_png(&path, &renderer.read_color(), args.width, args.height);
            eprintln!("wrote {}", path.display());
        }

        let result = PresetResult {
            preset: *preset,
            frames: cpu.len(),
            threshold_px: settings.quad_px_threshold,
            cpu_ms: percentiles(&mut cpu),
            gpu_ms: (!gpu_total.is_empty()).then(|| percentiles(&mut gpu_total)),
            cull_ms_p50: (!cull.is_empty()).then(|| percentiles(&mut cull).p50),
            draw_ms_p50: (!draw.is_empty()).then(|| percentiles(&mut draw).p50),
            hiz_ms_p50: (!hiz.is_empty()).then(|| percentiles(&mut hiz).p50),
            post_ms_p50: (!post.is_empty()).then(|| percentiles(&mut post).p50),
            visible_quads,
            visible_points,
            visible_bonds,
        };
        eprintln!(
            "{:?}: cpu p50 {:.2} ms p95 {:.2} ms | gpu p50 {} p95 {} | quads {} points {} bonds {} | threshold {:.1}px",
            preset,
            result.cpu_ms.p50,
            result.cpu_ms.p95,
            result
                .gpu_ms
                .as_ref()
                .map(|g| format!("{:.2} ms", g.p50))
                .unwrap_or_else(|| "n/a".into()),
            result
                .gpu_ms
                .as_ref()
                .map(|g| format!("{:.2} ms", g.p95))
                .unwrap_or_else(|| "n/a".into()),
            visible_quads,
            visible_points,
            visible_bonds,
            result.threshold_px
        );
        results.push(result);
    }

    let report = Report {
        adapter: ctx.adapter_name(),
        source,
        atoms: gpu.atom_count,
        bonds: gpu.bond_count,
        representation: format!("{:?}", args.rep),
        build_ms,
        width: args.width,
        height: args.height,
        adaptive: args.adaptive,
        occlusion: args.occlusion,
        results,
    };
    let json = serde_json::to_string_pretty(&report).expect("serialize report");
    match &args.json {
        Some(path) => std::fs::write(path, json).expect("write json"),
        None => println!("{json}"),
    }
}

fn write_png(path: &std::path::Path, rgba: &[u8], width: u32, height: u32) {
    let file = std::fs::File::create(path).expect("create png");
    let mut encoder = png::Encoder::new(std::io::BufWriter::new(file), width, height);
    encoder.set_color(png::ColorType::Rgba);
    encoder.set_depth(png::BitDepth::Eight);
    let mut writer = encoder.write_header().expect("png header");
    writer.write_image_data(rgba).expect("png data");
}
