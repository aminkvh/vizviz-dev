//! Ray-sphere impostor rendering: the CPU twin of `vv_render`'s
//! `fs_sphere`/`shade()` (see `crates/vv-render/src/shaders/draw.wgsl`).
//! Kept structurally identical to that shader on purpose — if you tune one,
//! check the other, or the two backends will silently drift apart.

use glam::Vec3;
use rayon::prelude::*;
use vv_core::Structure;
use vv_render::style::{Lighting, Material};
use vv_render::{Camera, ColorScheme, StylePreset};

use crate::bvh::Bvh;

/// Unpacks a little-endian `0xAABBGGRR` color (see `vv_render::color::rgba`)
/// into linear-order `[r, g, b]` in `0.0..=1.0`.
pub(crate) fn unpack_color(c: u32) -> Vec3 {
    Vec3::new(
        (c & 0xff) as f32,
        ((c >> 8) & 0xff) as f32,
        ((c >> 16) & 0xff) as f32,
    ) / 255.0
}

/// `vv_render::style::shade` on `Vec3`s.
pub(crate) fn shade(
    base: Vec3,
    n: Vec3,
    view_dir: Vec3,
    lighting: &Lighting,
    material: &Material,
) -> Vec3 {
    Vec3::from(vv_render::style::shade(
        base.into(),
        n.into(),
        view_dir.into(),
        lighting,
        material,
    ))
}

/// Nearest positive `t` where the ray `origin + t*dir` (`dir` must be unit
/// length) meets the sphere at `center` with the given `radius`, or `None`.
/// Same quadratic as `fs_sphere` in `draw.wgsl`, just solved on the CPU:
/// `a == 1` since `dir` is normalized, so this is the two-root formula with
/// the near root preferred (camera rays always start outside every atom).
pub(crate) fn intersect_sphere(origin: Vec3, dir: Vec3, center: Vec3, radius: f32) -> Option<f32> {
    let oc = origin - center;
    let b = oc.dot(dir);
    let c = oc.length_squared() - radius * radius;
    let disc = b * b - c;
    if disc < 0.0 {
        return None;
    }
    let sqrt_disc = disc.sqrt();
    let near = -b - sqrt_disc;
    if near > 0.0 {
        return Some(near);
    }
    let far = -b + sqrt_disc;
    (far > 0.0).then_some(far)
}

/// The camera's forward/right/up basis, derived the same way `Camera::view`
/// builds its view matrix (`f = Camera::forward`, `right =
/// cross(f, world_up)`, `up = cross(right, f)`) so ray directions land in
/// the same orientation `view()`/`proj()` would put the scene in.
pub(crate) fn camera_basis(camera: &Camera) -> (Vec3, Vec3, Vec3, Vec3) {
    // The CPU backend casts perspective rays only, so it takes the
    // perspective eye even when the camera is orthographic (whose eye is
    // pushed out beyond the scene; see `Camera::eye_distance`).
    let perspective = Camera {
        projection: vv_render::Projection::Perspective,
        ..camera.clone()
    };
    let eye = perspective.eye();
    let forward = camera.forward();
    let right = camera.orientation * Vec3::X;
    let up = camera.orientation * Vec3::Y;
    (eye, forward, right, up)
}

/// Atoms plus a prebuilt spatial index, ready to render from any camera.
/// Build once per structure and reuse across every frame the camera moves
/// through — the BVH build is the expensive one-time cost a real caller
/// shouldn't pay every frame (mirrors `vv_render::GpuStructure::upload`
/// once, `Renderer::render` many times).
pub struct CpuScene {
    positions: Vec<Vec3>,
    radii: Vec<f32>,
    colors: Vec<u32>,
    bvh: Bvh,
}

impl CpuScene {
    /// `positions`, `radii`, and `colors` (packed as `vv_render::color::
    /// rgba` does) must have the same length, one entry per atom.
    pub fn new(positions: Vec<Vec3>, radii: Vec<f32>, colors: Vec<u32>) -> Self {
        assert_eq!(positions.len(), radii.len());
        assert_eq!(positions.len(), colors.len());
        let bvh = Bvh::build(&positions, &radii);
        Self {
            positions,
            radii,
            colors,
            bvh,
        }
    }

    /// Van der Waals radii (spacefill, no representation scaling yet) and
    /// colors under `scheme`, both derived the same way `vv_render::
    /// GpuStructure::upload` does, so a CPU and GPU spacefill render of the
    /// same file should look the same.
    pub fn from_structure(structure: &Structure, scheme: ColorScheme) -> Self {
        let colors = vv_render::colors_for(scheme, &structure.topology);
        Self::from_structure_colored(structure, colors)
    }

    /// `from_structure` with one packed RGBA8 color per atom decided by the
    /// caller (a value column, an explicit color array from Python).
    pub fn from_structure_colored(structure: &Structure, colors: Vec<u32>) -> Self {
        assert_eq!(colors.len(), structure.topology.atom_count());
        let positions = structure.frame(0).positions().to_vec();
        let radii = structure
            .topology
            .element
            .iter()
            .map(|e| e.vdw_radius())
            .collect();
        Self::new(positions, radii, colors)
    }

    pub fn atom_count(&self) -> usize {
        self.positions.len()
    }

    /// Renders this scene as ray-cast sphere impostors, spacefill only (no
    /// bonds yet — see docs/RENDERING.md), into a tightly packed RGBA8
    /// buffer the same shape `vv_render::Renderer::read_color` returns.
    /// Parallelized by image row via `rayon`.
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
                let ndc_y = (1.0 - 2.0 * (y as f32 + 0.5) / height as f32) * tan_half_fov;
                for x in 0..width as usize {
                    let ndc_x =
                        (2.0 * (x as f32 + 0.5) / width as f32 - 1.0) * aspect * tan_half_fov;
                    let dir = (forward + right * ndc_x + up * ndc_y).normalize();

                    let px = x * 4;
                    match self.bvh.nearest_hit(eye, dir, &self.positions, &self.radii) {
                        None => row[px..px + 4].copy_from_slice(&background),
                        Some((t, i)) => {
                            let hit = eye + dir * t;
                            let n = (hit - self.positions[i]) / self.radii[i];
                            let base = unpack_color(self.colors[i]);
                            let shaded = shade(base, n, dir, &lighting, &material);
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
}

/// One-shot convenience over `CpuScene` for callers that don't need to
/// render more than once (tests, a single screenshot). Building a fresh
/// `CpuScene` (and so a fresh BVH) on every call is wasteful across
/// multiple frames of the same structure — use `CpuScene` directly then.
pub fn render_spheres(
    positions: &[Vec3],
    radii: &[f32],
    colors: &[u32],
    camera: &Camera,
    width: u32,
    height: u32,
    style: StylePreset,
    background: [u8; 4],
) -> Vec<u8> {
    CpuScene::new(positions.to_vec(), radii.to_vec(), colors.to_vec())
        .render(camera, width, height, style, background)
}

/// One-shot convenience over `CpuScene::from_structure` — see the same
/// caveat as `render_spheres` about repeated calls.
pub fn render_spacefill(
    structure: &Structure,
    scheme: ColorScheme,
    camera: &Camera,
    width: u32,
    height: u32,
    style: StylePreset,
    background: [u8; 4],
) -> Vec<u8> {
    CpuScene::from_structure(structure, scheme).render(camera, width, height, style, background)
}

/// `render_spacefill` with one packed RGBA8 color per atom instead of a
/// scheme (`vv_render::colors_from_scalar` turns a value column into one).
pub fn render_spacefill_colored(
    structure: &Structure,
    colors: Vec<u32>,
    camera: &Camera,
    width: u32,
    height: u32,
    style: StylePreset,
    background: [u8; 4],
) -> Vec<u8> {
    CpuScene::from_structure_colored(structure, colors)
        .render(camera, width, height, style, background)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ray_through_center_hits_at_distance_minus_radius() {
        let origin = Vec3::new(0.0, 0.0, 10.0);
        let dir = Vec3::new(0.0, 0.0, -1.0);
        let t = intersect_sphere(origin, dir, Vec3::ZERO, 2.0).expect("should hit");
        assert!((t - 8.0).abs() < 1e-5, "expected t=8.0, got {t}");
    }

    #[test]
    fn ray_missing_the_sphere_returns_none() {
        let origin = Vec3::new(5.0, 0.0, 10.0);
        let dir = Vec3::new(0.0, 0.0, -1.0);
        assert_eq!(intersect_sphere(origin, dir, Vec3::ZERO, 2.0), None);
    }

    #[test]
    fn ray_starting_inside_the_sphere_hits_the_far_side() {
        let origin = Vec3::ZERO;
        let dir = Vec3::new(0.0, 0.0, -1.0);
        let t = intersect_sphere(origin, dir, Vec3::ZERO, 2.0).expect("should hit exiting");
        assert!((t - 2.0).abs() < 1e-5, "expected t=2.0, got {t}");
    }

    #[test]
    fn sphere_entirely_behind_the_ray_is_not_hit() {
        let origin = Vec3::new(0.0, 0.0, 10.0);
        let dir = Vec3::new(0.0, 0.0, -1.0);
        // A sphere behind the camera, along the same line.
        assert_eq!(
            intersect_sphere(origin, dir, Vec3::new(0.0, 0.0, 20.0), 2.0),
            None
        );
    }

    #[test]
    fn a_single_centered_sphere_covers_the_middle_not_the_corners() {
        let (w, h) = (64u32, 64u32);
        let camera = Camera::framing(Vec3::ZERO, 5.0);
        let pixels = render_spheres(
            &[Vec3::ZERO],
            &[5.0],
            &[vv_render::color::rgba(255, 0, 0)],
            &camera,
            w,
            h,
            StylePreset::default(),
            [0, 0, 0, 255],
        );
        let center = ((h / 2 * w + w / 2) * 4) as usize;
        assert_ne!(
            &pixels[center..center + 3],
            &[0, 0, 0],
            "center pixel should be the sphere, not background"
        );
        let corner = 0usize;
        assert_eq!(
            &pixels[corner..corner + 3],
            &[0, 0, 0],
            "far corner of a framed sphere should be background"
        );
    }

    #[test]
    fn reusing_a_scene_across_frames_gives_the_same_image_as_building_fresh() {
        let camera = Camera::framing(Vec3::ZERO, 5.0);
        let scene = CpuScene::new(
            vec![Vec3::ZERO],
            vec![5.0],
            vec![vv_render::color::rgba(0, 255, 0)],
        );
        let a = scene.render(&camera, 32, 32, StylePreset::default(), [1, 2, 3, 255]);
        let b = render_spheres(
            &[Vec3::ZERO],
            &[5.0],
            &[vv_render::color::rgba(0, 255, 0)],
            &camera,
            32,
            32,
            StylePreset::default(),
            [1, 2, 3, 255],
        );
        assert_eq!(a, b);
    }
}
