//! CPU ray-marcher for `vv_core::gaussian_surface`'s density field -- the
//! twin of `renderer.rs`'s sphere-impostor path, but for a smooth,
//! ray-castable isosurface instead of many separate spheres, matching
//! this project's "no meshing" choice already made for spheres,
//! cylinders and tubes.
//!
//! Per-step atom lookups go through a `vv_core::spatial::Grid` whose
//! cell edge is the largest atom's `gaussian_surface::cutoff_radius`,
//! rather than every atom in the structure -- the same "measure the
//! brute-force cost, then grid it" progression this crate's own
//! sphere-impostor path took with its BVH (see `renderer.rs`'s module
//! doc and `CpuScene`). That cell size is what makes one query radius
//! (the cell edge) provably sufficient: any atom farther away than that
//! is farther than its own cutoff, so it contributes less than the
//! grid's `epsilon` and is safe to skip.

use glam::Vec3;
use rayon::prelude::*;
use vv_core::gaussian_surface::{self, Blob};
use vv_core::spatial::Grid;
use vv_core::Structure;
use vv_render::{Camera, ColorScheme, StylePreset};

use crate::renderer::{camera_basis, shade, unpack_color};

/// The isovalue every atom's own van der Waals surface sits exactly on
/// (`gaussian_surface`'s per-atom contribution is constructed so this is
/// always `1.0`, independent of `blob_factor`) -- the natural default
/// for a "hugs the atoms" Gaussian surface.
pub const DEFAULT_ISOVALUE: f32 = 1.0;

/// How finely a ray through the scene's bounding sphere is sampled
/// before bisection refines the crossing. Fixed step count (not a fixed
/// world-space step size) so cost per ray is bounded regardless of
/// structure size.
const MARCH_STEPS: usize = 256;
/// Bisection refinement steps once a crossing is bracketed: each halves
/// the bracket, so 10 steps resolve to roughly 1/1024 of one march step.
const REFINE_STEPS: usize = 10;

/// The "negligible past this" threshold the scene's grid cell size is
/// derived from (see the module doc).
const GRID_EPSILON: f32 = 0.01;

/// Atoms plus a spatial grid over their derived blobs and a bounding
/// sphere, ready to ray-march from any camera. Build once per
/// structure/frame and reuse across camera moves -- there's no
/// incremental update, matching `CpuScene`'s own "build once, render
/// many" contract. Rebuilding is cheap relative to a render (both are
/// `O(atoms)`, single-pass), so a moving trajectory frame is a fresh
/// `GaussianSurfaceScene::new` per frame, not an in-place update.
pub struct GaussianSurfaceScene {
    positions: Vec<Vec3>,
    blobs: Vec<Blob>,
    grid: Grid,
    colors: Vec<u32>,
    bounds_center: Vec3,
    bounds_radius: f32,
    blob_factor: f32,
    isovalue: f32,
}

impl GaussianSurfaceScene {
    /// `positions`, `radii`, and `colors` (packed as `vv_render::color::
    /// rgba` does) must have the same length, one entry per atom.
    pub fn new(positions: Vec<Vec3>, radii: Vec<f32>, colors: Vec<u32>, blob_factor: f32) -> Self {
        assert_eq!(positions.len(), radii.len());
        assert_eq!(positions.len(), colors.len());
        assert!(!positions.is_empty(), "need at least one atom");
        let blobs: Vec<Blob> = positions
            .iter()
            .zip(&radii)
            .map(|(&center, &radius)| Blob { center, radius })
            .collect();
        let (bounds_center, bounds_radius) = bounding_sphere(&positions, &radii, blob_factor);
        let indices: Vec<u32> = (0..positions.len() as u32).collect();
        let grid = Grid::build(
            &positions,
            &indices,
            gaussian_surface::grid_cell_size(&radii, blob_factor, GRID_EPSILON),
        );
        Self {
            positions,
            blobs,
            grid,
            colors,
            bounds_center,
            bounds_radius,
            blob_factor,
            isovalue: DEFAULT_ISOVALUE,
        }
    }

    /// Every blob that can contribute more than the grid's epsilon at
    /// `point`, into `out` (cleared first): the atoms within one cell
    /// edge, which is at least every atom's own cutoff radius.
    fn nearby(&self, point: Vec3, out: &mut Vec<Blob>) {
        out.clear();
        let radius = self.grid.cell_size();
        self.grid
            .for_each_within(&self.positions, point, radius, |i, _| {
                out.push(self.blobs[i as usize]);
            });
    }

    /// Van der Waals radii and colors under `scheme`, the same derivation
    /// `CpuScene::from_structure` uses for the sphere path, so a spacefill
    /// and a Gaussian-surface render of the same file agree on radii and
    /// colors (just not on how the surface between atoms is drawn).
    pub fn from_structure(structure: &Structure, scheme: ColorScheme, blob_factor: f32) -> Self {
        let colors = vv_render::colors_for(scheme, &structure.topology);
        Self::from_structure_colored(structure, colors, blob_factor)
    }

    /// `from_structure` with one packed color per atom already decided by
    /// the caller (a value column, an explicit color array from Python).
    pub fn from_structure_colored(
        structure: &Structure,
        colors: Vec<u32>,
        blob_factor: f32,
    ) -> Self {
        let positions = structure.frame(0).positions().to_vec();
        let radii = structure
            .topology
            .element
            .iter()
            .map(|e| e.vdw_radius())
            .collect();
        Self::new(positions, radii, colors, blob_factor)
    }

    pub fn atom_count(&self) -> usize {
        self.positions.len()
    }

    /// Renders this scene by ray-marching the density field, into a
    /// tightly packed RGBA8 buffer the same shape `CpuScene::render`
    /// produces. Parallelized by image row via `rayon`, same as the
    /// sphere path.
    pub fn render(
        &self,
        camera: &Camera,
        width: u32,
        height: u32,
        style: StylePreset,
        background: [u8; 4],
    ) -> Vec<u8> {
        let (lighting, material) = (style.lighting(), style.material());
        let (eye, forward, right, up) = camera_basis(camera);
        let aspect = width as f32 / height.max(1) as f32;
        let tan_half_fov = (camera.fov_y * 0.5).tan();

        let mut buffer = vec![0u8; (width as usize) * (height as usize) * 4];
        buffer
            .par_chunks_mut(width as usize * 4)
            .enumerate()
            .for_each(|(y, row)| {
                // One scratch buffer per row, reused for every pixel and
                // every march/refine step in that row -- the grid query
                // runs in the hottest of hot loops, so this keeps it
                // allocation-free rather than a fresh `Vec` per query.
                let mut nearby = Vec::new();
                let ndc_y = (1.0 - 2.0 * (y as f32 + 0.5) / height as f32) * tan_half_fov;
                for x in 0..width as usize {
                    let ndc_x =
                        (2.0 * (x as f32 + 0.5) / width as f32 - 1.0) * aspect * tan_half_fov;
                    let dir = (forward + right * ndc_x + up * ndc_y).normalize();
                    let px = x * 4;
                    match self.march(eye, dir, &mut nearby) {
                        None => row[px..px + 4].copy_from_slice(&background),
                        Some((normal, nearest_atom)) => {
                            let base = unpack_color(self.colors[nearest_atom]);
                            let shaded = shade(base, normal, dir, &lighting, &material);
                            row[px] = (shaded.x.clamp(0.0, 1.0) * 255.0).round() as u8;
                            row[px + 1] = (shaded.y.clamp(0.0, 1.0) * 255.0).round() as u8;
                            row[px + 2] = (shaded.z.clamp(0.0, 1.0) * 255.0).round() as u8;
                            row[px + 3] = 255;
                        }
                    }
                }
            });
        buffer
    }

    /// Marches one ray through the scene's bounding sphere; `None` if
    /// the isosurface is never crossed. On a hit, returns the outward
    /// normal (`-gradient`, since `gradient` points toward increasing
    /// density, i.e. into the surface) and the index of the atom nearest
    /// the hit point, used to color it. Coloring by nearest atom rather
    /// than a distance-weighted blend of overlapping atoms is a
    /// deliberate phase-1 simplification -- see the module doc. `nearby`
    /// is caller-owned scratch space for the grid query, reused across
    /// every step of this ray (and every other ray in the caller's row)
    /// rather than allocated fresh each time.
    fn march(&self, eye: Vec3, dir: Vec3, nearby: &mut Vec<Blob>) -> Option<(Vec3, usize)> {
        let (t_enter, t_exit) = self.bounds_interval(eye, dir)?;
        let step = (t_exit - t_enter) / MARCH_STEPS as f32;
        let mut prev_t = t_enter;
        self.nearby(eye + dir * prev_t, nearby);
        let mut prev_d = gaussian_surface::density(eye + dir * prev_t, &*nearby, self.blob_factor);
        for i in 1..=MARCH_STEPS {
            let t = t_enter + step * i as f32;
            self.nearby(eye + dir * t, nearby);
            let d = gaussian_surface::density(eye + dir * t, &*nearby, self.blob_factor);
            if prev_d < self.isovalue && d >= self.isovalue {
                let mut lo = prev_t;
                let mut hi = t;
                for _ in 0..REFINE_STEPS {
                    let mid = 0.5 * (lo + hi);
                    self.nearby(eye + dir * mid, nearby);
                    let dm = gaussian_surface::density(eye + dir * mid, &*nearby, self.blob_factor);
                    if dm < self.isovalue {
                        lo = mid;
                    } else {
                        hi = mid;
                    }
                }
                let hit = eye + dir * hi;
                self.nearby(hit, nearby);
                let grad = gaussian_surface::gradient(hit, &*nearby, self.blob_factor);
                let normal = (-grad).normalize_or_zero();
                return Some((normal, self.nearest_atom(hit)));
            }
            prev_t = t;
            prev_d = d;
        }
        None
    }

    /// The ray's entry/exit parametric distances through the scene's
    /// bounding sphere, entry clamped to `0` (a camera already inside
    /// the sphere still marches from itself, not from behind it).
    fn bounds_interval(&self, eye: Vec3, dir: Vec3) -> Option<(f32, f32)> {
        let oc = eye - self.bounds_center;
        let b = oc.dot(dir);
        let c = oc.length_squared() - self.bounds_radius * self.bounds_radius;
        let disc = b * b - c;
        if disc < 0.0 {
            return None;
        }
        let sqrt_disc = disc.sqrt();
        let t_exit = -b + sqrt_disc;
        if t_exit <= 0.0 {
            return None;
        }
        let t_enter = (-b - sqrt_disc).max(0.0);
        Some((t_enter, t_exit))
    }

    fn nearest_atom(&self, point: Vec3) -> usize {
        self.positions
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| {
                (point - **a)
                    .length_squared()
                    .total_cmp(&(point - **b).length_squared())
            })
            .map(|(i, _)| i)
            .expect("scene has at least one atom (checked in new())")
    }
}

/// A sphere guaranteed to contain the whole field within `epsilon` of
/// zero at its boundary (each atom's own [`gaussian_surface::
/// cutoff_radius`] past its position), so marching stops at a boundary
/// where the surface truly has ended, not an arbitrarily chosen box.
fn bounding_sphere(positions: &[Vec3], radii: &[f32], blob_factor: f32) -> (Vec3, f32) {
    let (min, max) = positions.iter().fold(
        (Vec3::splat(f32::INFINITY), Vec3::splat(f32::NEG_INFINITY)),
        |(min, max), &p| (min.min(p), max.max(p)),
    );
    let center = (min + max) * 0.5;
    let radius = positions
        .iter()
        .zip(radii)
        .map(|(&p, &r)| {
            (p - center).length() + gaussian_surface::cutoff_radius(r, blob_factor, 0.05)
        })
        .fold(0.0f32, f32::max);
    (center, radius.max(1.0))
}

/// One-shot convenience over `GaussianSurfaceScene` for callers that
/// don't need to render more than once -- see `CpuScene`'s own
/// `render_spacefill` for the same caveat about repeated calls.
pub fn render_gaussian_surface(
    structure: &Structure,
    scheme: ColorScheme,
    blob_factor: f32,
    camera: &Camera,
    width: u32,
    height: u32,
    style: StylePreset,
    background: [u8; 4],
) -> Vec<u8> {
    GaussianSurfaceScene::from_structure(structure, scheme, blob_factor)
        .render(camera, width, height, style, background)
}

#[cfg(test)]
mod tests {
    use super::*;

    const WHITE: u32 = 0xFFFFFFFF; // little-endian RGBA8: r=g=b=a=0xFF

    fn white_scene(positions: Vec<Vec3>, radii: Vec<f32>) -> GaussianSurfaceScene {
        let n = positions.len();
        GaussianSurfaceScene::new(
            positions,
            radii,
            vec![WHITE; n],
            gaussian_surface::DEFAULT_BLOB_FACTOR,
        )
    }

    #[test]
    fn a_single_atom_is_hit_dead_center_and_missed_at_the_corners() {
        let scene = white_scene(vec![Vec3::ZERO], vec![1.5]);
        let camera = Camera::framing(Vec3::ZERO, 1.5);
        let (w, h) = (64, 64);
        let buf = scene.render(&camera, w, h, StylePreset::ALL[0], [10, 10, 10, 255]);
        let px = |x: u32, y: u32| {
            let i = ((y * w + x) * 4) as usize;
            (buf[i], buf[i + 1], buf[i + 2], buf[i + 3])
        };
        let center = px(w / 2, h / 2);
        assert_ne!(center, (10, 10, 10, 255), "center should hit the surface");
        let corner = px(0, 0);
        assert_eq!(
            corner,
            (10, 10, 10, 255),
            "corner should miss, showing background"
        );
    }

    #[test]
    fn two_close_atoms_render_a_single_blended_blob_wider_than_either_atom_alone() {
        // The visual proof of `gaussian_surface`'s own blending test:
        // count hit pixels along the horizontal midline for one atom vs.
        // two overlapping atoms side by side -- the pair should cover
        // more than either atom would need alone, because the field
        // between them lifts above the isovalue too.
        let radius = 1.0;
        let one = white_scene(vec![Vec3::ZERO], vec![radius]);
        let two = white_scene(
            vec![Vec3::new(-0.9, 0.0, 0.0), Vec3::new(0.9, 0.0, 0.0)],
            vec![radius, radius],
        );
        let camera = Camera::framing(Vec3::ZERO, 3.0);
        let (w, h) = (128, 32);
        let bg = [0, 0, 0, 255];
        let count_hits = |buf: &[u8]| {
            (0..w)
                .filter(|&x| {
                    let i = ((h / 2 * w + x) * 4) as usize;
                    (buf[i], buf[i + 1], buf[i + 2], buf[i + 3]) != (bg[0], bg[1], bg[2], bg[3])
                })
                .count()
        };
        let one_hits = count_hits(&one.render(&camera, w, h, StylePreset::ALL[0], bg));
        let two_hits = count_hits(&two.render(&camera, w, h, StylePreset::ALL[0], bg));
        assert!(
            two_hits > one_hits,
            "two atoms ({two_hits}px) should cover more than one ({one_hits}px)"
        );
    }

    #[test]
    fn grid_accelerated_density_matches_brute_force_within_epsilon() {
        // The grid trades an `epsilon`-sized truncation (one skipped
        // atom's worst-case contribution) for never scanning every atom;
        // several near-boundary atoms can each contribute up to epsilon,
        // so allow a few multiples of it, still tight against the
        // isovalue of 1.0.
        let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../fixtures/small/1CRN.pdb");
        let structure = vv_io::load(&path).unwrap();
        let scene = GaussianSurfaceScene::from_structure(
            &structure,
            ColorScheme::Element,
            gaussian_surface::DEFAULT_BLOB_FACTOR,
        );
        let mut buf = Vec::new();
        let mut max_abs_diff = 0.0f32;
        for (i, &p) in scene.positions.iter().enumerate() {
            if i % 7 != 0 {
                continue;
            }
            let probe = p + Vec3::new(0.3, -0.2, 0.4);
            scene.nearby(probe, &mut buf);
            let via_grid = gaussian_surface::density(probe, &buf, scene.blob_factor);
            let brute = gaussian_surface::density(probe, &scene.blobs, scene.blob_factor);
            max_abs_diff = max_abs_diff.max((via_grid - brute).abs());
        }
        assert!(
            max_abs_diff < GRID_EPSILON * 3.0,
            "max diff on real structure: {max_abs_diff}"
        );
    }

    #[test]
    fn a_real_structure_renders_without_panicking_and_covers_the_center() {
        let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../fixtures/small/1CRN.pdb");
        let structure = vv_io::load(&path).unwrap();
        let (center, radius) = structure.frame(0).bounding_sphere().unwrap();
        let camera = Camera::framing(center, radius);
        let scene = GaussianSurfaceScene::from_structure(
            &structure,
            ColorScheme::Element,
            gaussian_surface::DEFAULT_BLOB_FACTOR,
        );
        assert_eq!(scene.atom_count(), structure.atom_count());
        let (w, h) = (200, 150);
        let bg = [30, 30, 30, 255];
        let buf = scene.render(&camera, w, h, StylePreset::ALL[0], bg);
        let i = ((h / 2 * w + w / 2) * 4) as usize;
        assert_ne!(
            (buf[i], buf[i + 1], buf[i + 2]),
            (bg[0], bg[1], bg[2]),
            "a structure framed to fill the view should hit something dead center"
        );
    }

    /// Not run by default (`cargo test -- --ignored`): renders a real PNG
    /// to disk for a human to actually look at, the same "verify by
    /// opening the file" discipline the rest of this session's work used
    /// (screenshots, SVGs) rather than trusting pixel-count assertions
    /// alone as proof the surface looks right.
    #[test]
    #[ignore]
    fn render_1crn_to_png_for_manual_inspection() {
        let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../fixtures/small/1CRN.pdb");
        let structure = vv_io::load(&path).unwrap();
        let (center, radius) = structure.frame(0).bounding_sphere().unwrap();
        let camera = Camera::framing(center, radius);
        let scene = GaussianSurfaceScene::from_structure(
            &structure,
            ColorScheme::Chain,
            gaussian_surface::DEFAULT_BLOB_FACTOR,
        );
        let (w, h) = (800, 600);
        let buf = scene.render(&camera, w, h, StylePreset::ALL[0], [85, 85, 93, 255]);

        let out = std::env::temp_dir().join("vizviz_gaussian_surface_1crn.png");
        let file = std::fs::File::create(&out).unwrap();
        let mut encoder = png::Encoder::new(std::io::BufWriter::new(file), w, h);
        encoder.set_color(png::ColorType::Rgba);
        encoder.set_depth(png::BitDepth::Eight);
        encoder
            .write_header()
            .unwrap()
            .write_image_data(&buf)
            .unwrap();
        println!("wrote {}", out.display());
    }
}
