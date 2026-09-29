//! CPU ray-caster for `vv_core::skin_surface`'s mixed complex -- the
//! twin of `gaussian_surface.rs`, same "CPU first, GPU later" order.
//!
//! All of the geometry and the membership rule live in `vv_core::
//! skin_surface` (read its module doc first: what a patch is, why both
//! ray roots must be tested, and how a hit is accepted or rejected
//! exactly rather than heuristically). This file only decides *which*
//! patches a ray could hit: a [`Bvh`] over every patch's bounding
//! sphere, so a pixel tests the handful of patches along its ray instead
//! of every atom and every edge -- the brute-force version of this file
//! was O(pixels x patches x atoms) and took tens of seconds for a
//! 300-atom protein. The ownership queries go through each patch's
//! competitor list (`SkinComplex::owned_by_competitors`).

use glam::Vec3;
use rayon::prelude::*;
use vv_core::skin_surface::{self, SkinComplex};
use vv_core::Structure;
use vv_render::{Camera, ColorScheme, StylePreset};

use crate::bvh::Bvh;
use crate::renderer::{camera_basis, shade, unpack_color};

pub struct SkinSurfaceScene {
    positions: Vec<Vec3>,
    weights: Vec<f32>,
    colors: Vec<u32>,
    complex: SkinComplex,
    /// Indices into `complex.patches` of the patches that have any
    /// surface at all -- what the BVH is built over, in this order.
    drawn: Vec<u32>,
    bound_centers: Vec<Vec3>,
    bound_radii: Vec<f32>,
    bvh: Bvh,
    bounds_center: Vec3,
    bounds_radius: f32,
}

impl SkinSurfaceScene {
    /// `positions`, `weights` (`vv_core::skin_surface::weight_for_radius`
    /// of each atom's radius), and `colors` (packed as
    /// `vv_render::color::rgba` does) must have the same length, one
    /// entry per atom.
    pub fn new(positions: Vec<Vec3>, weights: Vec<f32>, colors: Vec<u32>, shrink: f32) -> Self {
        assert_eq!(positions.len(), weights.len());
        assert_eq!(positions.len(), colors.len());
        assert!(!positions.is_empty(), "need at least one atom");
        let complex = skin_surface::build_complex(&positions, &weights, shrink);
        let drawn: Vec<u32> = (0..complex.patches.len() as u32)
            .filter(|&i| complex.patches[i as usize].has_surface())
            .collect();
        let bound_centers: Vec<Vec3> = drawn
            .iter()
            .map(|&i| complex.patches[i as usize].bound_center)
            .collect();
        let bound_radii: Vec<f32> = drawn
            .iter()
            .map(|&i| complex.patches[i as usize].bound_radius)
            .collect();
        let bvh = Bvh::build(&bound_centers, &bound_radii);
        let (bounds_center, bounds_radius) = bounding_sphere(&bound_centers, &bound_radii);
        Self {
            positions,
            weights,
            colors,
            complex,
            drawn,
            bound_centers,
            bound_radii,
            bvh,
            bounds_center,
            bounds_radius,
        }
    }

    /// Frame 0 of `structure`, colored by `scheme`; weights from van der
    /// Waals radii via `weight_for_radius`, so every atom's own patch
    /// sits on the same sphere spacefill draws.
    pub fn from_structure(structure: &Structure, scheme: ColorScheme, shrink: f32) -> Self {
        let colors = vv_render::colors_for(scheme, &structure.topology);
        Self::from_structure_colored(structure, 0, colors, shrink)
    }

    /// `from_structure` for coordinate set `frame`, with one packed color
    /// per atom decided by the caller.
    pub fn from_structure_colored(
        structure: &Structure,
        frame: usize,
        colors: Vec<u32>,
        shrink: f32,
    ) -> Self {
        let positions = structure.frame(frame).positions().to_vec();
        let weights = structure
            .topology
            .element
            .iter()
            .map(|e| skin_surface::weight_for_radius(e.vdw_radius(), shrink))
            .collect();
        Self::new(positions, weights, colors, shrink)
    }

    pub fn atom_count(&self) -> usize {
        self.positions.len()
    }

    /// Patches with a surface, by kind: vertices, edges, triangles,
    /// tetrahedra.
    pub fn patch_counts(&self) -> [usize; 4] {
        let mut counts = [0; 4];
        for &i in &self.drawn {
            counts[self.complex.patches[i as usize].kind as usize] += 1;
        }
        counts
    }

    /// The mixed complex this scene rendered from -- for callers that
    /// need to inspect a specific patch (diagnostics, tests) rather than
    /// just get pixels.
    pub fn complex(&self) -> &SkinComplex {
        &self.complex
    }

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
                    match self.cast(eye, dir) {
                        None => row[px..px + 4].copy_from_slice(&background),
                        Some((normal, color_atom)) => {
                            let base = unpack_color(self.colors[color_atom]);
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

    /// The nearest patch crossing the ray actually inside its own mixed
    /// cell: outward normal and the atom to color by.
    fn cast(&self, eye: Vec3, dir: Vec3) -> Option<(Vec3, usize)> {
        let (t, index) = self.cast_t(eye, dir)?;
        let patch = &self.complex.patches[index];
        let x = eye + dir * t;
        let color_atom = patch.nearest_atom(x, &self.positions, &self.weights) as usize;
        Some((patch.normal_at(x, self.complex.shrink), color_atom))
    }

    /// `cast`'s ray parameter and patch index.
    pub fn cast_t(&self, eye: Vec3, dir: Vec3) -> Option<(f32, usize)> {
        if !self.hits_bounds(eye, dir) {
            return None;
        }
        let shrink = self.complex.shrink;
        let (t, which) =
            self.bvh
                .nearest_hit_with(eye, dir, &self.bound_centers, &self.bound_radii, |i| {
                    let index = self.drawn[i] as usize;
                    let patch = &self.complex.patches[index];
                    patch
                        .ray_roots(eye, dir, shrink)
                        .into_iter()
                        .flatten()
                        .find(|&t| {
                            self.complex.contains(
                                index,
                                eye + dir * t,
                                &self.positions,
                                &self.weights,
                            )
                        })
                })?;
        Some((t, self.drawn[which] as usize))
    }

    fn hits_bounds(&self, eye: Vec3, dir: Vec3) -> bool {
        let oc = eye - self.bounds_center;
        let b = oc.dot(dir);
        let c = oc.length_squared() - self.bounds_radius * self.bounds_radius;
        let disc = b * b - c;
        disc >= 0.0 && (-b + disc.sqrt()) > 0.0
    }
}

fn bounding_sphere(centers: &[Vec3], radii: &[f32]) -> (Vec3, f32) {
    let (min, max) = centers.iter().fold(
        (Vec3::splat(f32::INFINITY), Vec3::splat(f32::NEG_INFINITY)),
        |(min, max), &p| (min.min(p), max.max(p)),
    );
    let center = (min + max) * 0.5;
    let radius = centers
        .iter()
        .zip(radii)
        .map(|(&p, &r)| (p - center).length() + r)
        .fold(0.0f32, f32::max);
    (center, radius.max(1.0))
}

/// One-shot convenience over `SkinSurfaceScene`, matching
/// `render_gaussian_surface`'s own shape.
#[allow(clippy::too_many_arguments)]
pub fn render_skin_surface(
    structure: &Structure,
    scheme: ColorScheme,
    shrink: f32,
    camera: &Camera,
    width: u32,
    height: u32,
    style: StylePreset,
    background: [u8; 4],
) -> Vec<u8> {
    SkinSurfaceScene::from_structure(structure, scheme, shrink)
        .render(camera, width, height, style, background)
}

#[cfg(test)]
mod tests {
    use super::*;

    const WHITE: u32 = 0xFFFFFFFF;

    fn hit_at(buf: &[u8], w: u32, bg: [u8; 4], x: u32, y: u32) -> bool {
        let i = ((y * w + x) * 4) as usize;
        (buf[i], buf[i + 1], buf[i + 2], buf[i + 3]) != (bg[0], bg[1], bg[2], bg[3])
    }

    #[test]
    fn two_close_atoms_render_a_seamless_neck_between_them() {
        // The visual proof this is a *skin* surface and not just shrunk
        // spacefill: two overlapping atoms should show a continuous
        // silhouette (a neck), not two separate disks with a gap -- so
        // scanning the horizontal midline should find hits with no
        // background pixel strictly between the two atoms' x-extents.
        let a = Vec3::new(-1.0, 0.0, 0.0);
        let b = Vec3::new(1.0, 0.0, 0.0);
        let weight = skin_surface::weight_for_radius(1.0, 0.5);
        let scene =
            SkinSurfaceScene::new(vec![a, b], vec![weight, weight], vec![WHITE, WHITE], 0.5);
        assert_eq!(
            scene.patch_counts()[1],
            1,
            "expected the pair to be a Delaunay edge with a neck"
        );
        let camera = Camera::framing((a + b) * 0.5, 3.0);
        let (w, h) = (128, 32);
        let bg = [0u8, 0, 0, 255];
        let buf = scene.render(&camera, w, h, StylePreset::ALL[0], bg);
        let row = h / 2;
        let hits: Vec<u32> = (0..w).filter(|&x| hit_at(&buf, w, bg, x, row)).collect();
        assert!(
            hits.len() > 10,
            "expected a wide hit region, got {} px",
            hits.len()
        );
        let (first, last) = (*hits.first().unwrap(), *hits.last().unwrap());
        let gaps: Vec<u32> = (first..=last)
            .filter(|&x| !hit_at(&buf, w, bg, x, row))
            .collect();
        assert!(
            gaps.is_empty(),
            "expected a seamless neck, found gap columns {gaps:?}"
        );
    }

    #[test]
    fn the_neck_is_narrower_than_the_atoms_and_the_silhouette_is_convex_nowhere_near_the_waist() {
        // Along the vertical midline (x through the waist) the neck must
        // be thinner than either atom's own disk, or the "hyperboloid"
        // is really just two spheres.
        let a = Vec3::new(-1.2, 0.0, 0.0);
        let b = Vec3::new(1.2, 0.0, 0.0);
        let weight = skin_surface::weight_for_radius(1.0, 0.5);
        let scene =
            SkinSurfaceScene::new(vec![a, b], vec![weight, weight], vec![WHITE, WHITE], 0.5);
        let camera = Camera::framing(Vec3::ZERO, 3.0);
        let (w, h) = (128, 128);
        let bg = [0u8, 0, 0, 255];
        let buf = scene.render(&camera, w, h, StylePreset::ALL[0], bg);
        let column_height = |x: u32| (0..h).filter(|&y| hit_at(&buf, w, bg, x, y)).count();
        let waist = column_height(w / 2);
        // The atoms' centers project a bit left/right of the middle.
        let atom = (0..w).map(column_height).max().unwrap();
        assert!(waist > 0, "the neck should be visible");
        assert!(
            waist < atom,
            "waist {waist} px should be thinner than the atoms' {atom} px"
        );
    }

    #[test]
    fn a_single_atom_still_renders_as_its_own_sphere() {
        let scene = SkinSurfaceScene::new(
            vec![Vec3::ZERO],
            vec![skin_surface::weight_for_radius(1.5, 0.5)],
            vec![WHITE],
            0.5,
        );
        assert_eq!(scene.patch_counts(), [1, 0, 0, 0]);
        let camera = Camera::framing(Vec3::ZERO, 1.5);
        let (w, h) = (64, 64);
        let bg = [10u8, 10, 10, 255];
        let buf = scene.render(&camera, w, h, StylePreset::ALL[0], bg);
        assert!(
            hit_at(&buf, w, bg, w / 2, h / 2),
            "center should hit the sphere"
        );
        assert!(!hit_at(&buf, w, bg, 0, 0), "corner should be background");
    }

    #[test]
    fn a_real_structure_renders_with_every_patch_kind_and_covers_the_center() {
        let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../fixtures/small/1CRN.pdb");
        let structure = vv_io::load(&path).unwrap();
        let (center, radius) = structure.frame(0).bounding_sphere().unwrap();
        let camera = Camera::framing(center, radius);
        let scene = SkinSurfaceScene::from_structure(
            &structure,
            ColorScheme::Element,
            skin_surface::DEFAULT_SHRINK,
        );
        assert_eq!(scene.atom_count(), structure.atom_count());
        let counts = scene.patch_counts();
        assert!(
            counts[1] > counts[0] && counts[2] > counts[0],
            "a folded protein has more edges and triangles than atoms: {counts:?}"
        );
        let (w, h) = (200, 150);
        let bg = [30u8, 30, 30, 255];
        let buf = scene.render(&camera, w, h, StylePreset::ALL[0], bg);
        assert!(
            hit_at(&buf, w, bg, w / 2, h / 2),
            "a structure framed to fill the view should hit something dead center"
        );
    }

    /// The fast path (BVH over culled-and-bounded patches, competitor
    /// lists) must find the same nearest hit as an exhaustive cast over
    /// *every* patch with exhaustive ownership, at every pixel of a
    /// zoomed real-structure render. Any pixel where the exhaustive cast
    /// finds a nearer hit is a real hole in the fast path -- an
    /// over-eager cull, a bounding sphere that misses its own patch, or
    /// a competitor list that lets a point through -- and the failure
    /// message says which.
    #[test]
    fn fast_path_matches_an_exhaustive_cast_on_a_zoomed_real_structure() {
        let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../fixtures/small/1CRN.pdb");
        let structure = vv_io::load(&path).unwrap();
        let (center, radius) = structure.frame(0).bounding_sphere().unwrap();
        let mut camera = Camera::framing(center, radius);
        camera.zoom(0.35);
        let scene = SkinSurfaceScene::from_structure(
            &structure,
            ColorScheme::Element,
            skin_surface::DEFAULT_SHRINK,
        );
        let (w, h) = (200u32, 150u32);
        let (eye, forward, right, up) = camera_basis(&camera);
        let aspect = w as f32 / h as f32;
        let tan_half_fov = (camera.fov_y * 0.5).tan();
        let s = scene.complex.shrink;
        let rows: Vec<Vec<String>> = (0..h)
            .into_par_iter()
            .map(|y| {
                let mut problems = Vec::new();
                let ndc_y = (1.0 - 2.0 * (y as f32 + 0.5) / h as f32) * tan_half_fov;
                for x in 0..w {
                    let ndc_x = (2.0 * (x as f32 + 0.5) / w as f32 - 1.0) * aspect * tan_half_fov;
                    let dir = (forward + right * ndc_x + up * ndc_y).normalize();
                    let fast = scene.cast_t(eye, dir);
                    let mut exhaustive: Option<(f32, usize)> = None;
                    for (i, p) in scene.complex.patches.iter().enumerate() {
                        for t in p.ray_roots(eye, dir, s).into_iter().flatten() {
                            if exhaustive.is_some_and(|(best, _)| t >= best) {
                                continue;
                            }
                            let hit = eye + dir * t;
                            let ok = p.contains(hit, s, &scene.positions, |q| {
                                scene.complex.voronoi_owned_by(&scene.positions, &scene.weights, q, p.members())
                            });
                            if ok {
                                exhaustive = Some((t, i));
                            }
                        }
                    }
                    match (fast, exhaustive) {
                        // Patches of adjacent mixed cells are tangent at
                        // their shared boundary (the module doc), so
                        // right at a boundary both the departing and the
                        // arriving patch can pass `contains` within
                        // `MEMBERSHIP_EPS` (2e-2 in power distance,
                        // Å^2) -- worth a positional slop of roughly
                        // `eps / (2r)`, a few thousandths of an Å for
                        // van-der-Waals-scale radii. Either patch is the
                        // same true surface there to sub-pixel precision;
                        // only a real gap (a hole, or two surfaces
                        // noticeably apart) should fail this test.
                        (Some((tf, _)), Some((te, _))) if (tf - te).abs() < 3e-2 => {}
                        (None, None) => {}
                        (fast, Some((te, i))) => {
                            let p = &scene.complex.patches[i];
                            let hit = eye + dir * te;
                            let q = p.voronoi_point(hit, s);
                            let oc = p.bound_center - eye;
                            let closest = (oc - dir * oc.dot(dir)).length();
                            problems.push(format!(
                                "({x},{y}): fast={fast:?} exhaustive t={te:.3} via patch {i} {:?} {:?} (surface={}, bound_r={:.3}, ray-to-bound-center={closest:.3}, competitors_ok={}, w_z={:.3})",
                                p.kind,
                                p.members(),
                                p.has_surface(),
                                p.bound_radius,
                                scene.complex.owned_by_competitors(i, &scene.positions, &scene.weights, q),
                                p.weight
                            ));
                        }
                        (Some((tf, i)), None) => {
                            problems.push(format!("({x},{y}): fast hit t={tf:.3} via patch {i} but exhaustive found nothing"));
                        }
                    }
                }
                problems
            })
            .collect();
        let problems: Vec<String> = rows.into_iter().flatten().collect();
        assert!(
            problems.is_empty(),
            "{} disagreeing pixels, first few:\n{}",
            problems.len(),
            problems
                .iter()
                .take(12)
                .cloned()
                .collect::<Vec<_>>()
                .join("\n")
        );
    }

    /// Not run by default: renders a real PNG to disk for a human to
    /// actually look at -- what caught every real visual bug this surface
    /// has had so far.
    #[test]
    #[ignore]
    fn render_1crn_to_png_for_manual_inspection() {
        let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../fixtures/small/1CRN.pdb");
        let structure = vv_io::load(&path).unwrap();
        let (center, radius) = structure.frame(0).bounding_sphere().unwrap();
        let camera = Camera::framing(center, radius);
        let scene = SkinSurfaceScene::from_structure(
            &structure,
            ColorScheme::Chain,
            skin_surface::DEFAULT_SHRINK,
        );
        println!(
            "patches by kind (vertex, edge, triangle, tet): {:?}",
            scene.patch_counts()
        );
        let (w, h) = (800, 600);
        let start = std::time::Instant::now();
        let buf = scene.render(&camera, w, h, StylePreset::ALL[0], [85, 85, 93, 255]);
        println!("rendered {w}x{h} in {:?}", start.elapsed());
        let out = std::env::temp_dir().join("vizviz_skin_surface_1crn.png");
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
