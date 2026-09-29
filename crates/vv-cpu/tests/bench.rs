//! Not a correctness test: memorialized performance numbers, run manually
//! (`cargo test -p vv-cpu --release --test bench -- --ignored --nocapture`)
//! rather than on every `cargo test`, since they're meant to be read, not
//! just pass/fail. See docs/RENDERING.md for what these numbers mean (and
//! don't mean) for real use.
//!
//! Each builds the `CpuScene` (and so its BVH) once, outside the timed
//! loop, matching how a real caller would use it across many frames of the
//! same structure as the camera moves.

use std::time::Instant;

use vv_cpu::CpuScene;
use vv_render::{Camera, ColorScheme, StylePreset};

fn fixture(name: &str) -> std::path::PathBuf {
    std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(format!("../../fixtures/{name}"))
}

/// Renders `frames` times after one warm-up render, returns (ms/frame, fps).
fn measure(scene: &CpuScene, camera: &Camera, w: u32, h: u32, style: StylePreset) -> (f64, f64) {
    let background = [23, 23, 28, 255]; // matches StylePreset::DarkPresentation, roughly
    let _ = scene.render(camera, w, h, style, background); // warm up: page faults, first-touch
    let frames = 30;
    let start = Instant::now();
    for _ in 0..frames {
        let _ = scene.render(camera, w, h, style, background);
    }
    let per_frame_ms = start.elapsed().as_secs_f64() * 1000.0 / frames as f64;
    (per_frame_ms, 1000.0 / per_frame_ms)
}

/// The concrete target this crate exists to hit: interactive framerate on
/// a small, common structure size, on a laptop CPU, with no GPU at all.
/// Brute force (no BVH) measured 2.0 fps here on this machine -- 5x short
/// -- which is what motivated the BVH in src/bvh.rs; that got to ~58 fps,
/// a ~28x speedup on exactly this case.
#[test]
#[ignore]
fn bench_4hhb_512x512() {
    let structure = vv_io::load(fixture("small/4HHB.cif")).expect("load 4HHB");
    let (center, radius) = structure
        .frame(0)
        .bounding_sphere()
        .expect("4HHB has atoms");
    let camera = Camera::framing(center, radius);
    let scene = CpuScene::from_structure(&structure, ColorScheme::Element);
    let atoms = scene.atom_count();

    let (per_frame_ms, fps) = measure(&scene, &camera, 512, 512, StylePreset::default());
    eprintln!(
        "4HHB ({atoms} atoms) at 512x512, BVH-accelerated CPU, {} threads: \
         {per_frame_ms:.2} ms/frame, {fps:.1} fps",
        rayon::current_num_threads()
    );
    assert!(
        fps >= 10.0,
        "target: 4HHB at 512x512 should render at >=10 fps on CPU (BVH build excluded, \
         this machine's {} rayon threads); got {fps:.1} fps",
        rayon::current_num_threads()
    );
}

/// Informational only, no target: where does this fall over? 3J3Q (2.44M
/// atoms) is the large real-world fixture used throughout this project's
/// GPU benchmarks (see benchmarks/README.md) -- this is not remotely the
/// use case vv-cpu targets (that's the GPU renderer's job, at 10-14 ms for
/// 10M atoms), but the number is worth having on record rather than
/// guessing at where brute-force-per-pixel BVH traversal (still O(log n)
/// per pixel, but with a much taller tree and much larger images expected
/// at this scale) stops being usable.
#[test]
#[ignore]
fn bench_3j3q_512x512() {
    let structure = vv_io::load(fixture("large/3J3Q.cif.gz")).expect(
        "load 3J3Q -- run `cargo xtask fetch` first if this fixture hasn't been downloaded",
    );
    let (center, radius) = structure
        .frame(0)
        .bounding_sphere()
        .expect("3J3Q has atoms");
    let camera = Camera::framing(center, radius);
    let scene = CpuScene::from_structure(&structure, ColorScheme::Element);
    let atoms = scene.atom_count();

    let (per_frame_ms, fps) = measure(&scene, &camera, 512, 512, StylePreset::default());
    eprintln!(
        "3J3Q ({atoms} atoms) at 512x512, BVH-accelerated CPU, {} threads: \
         {per_frame_ms:.2} ms/frame, {fps:.1} fps (informational; no target at this scale)",
        rayon::current_num_threads()
    );
}
