//! `render`: a path-traced image of what the viewport shows
//! (`vv_render::path_trace`), rendered on a worker thread so the window
//! stays responsive, with its progress in the notice bar.

use std::path::PathBuf;
use std::sync::atomic::{AtomicU32, Ordering};
use std::sync::{mpsc, Arc};
use std::time::Instant;

use glam::Vec3;
use vv_render::path_trace::{PathTracer, TraceScene, TraceSettings};
use vv_render::{Camera, GpuContext, Representation as GpuRepresentation};
use vv_scene::{LoadedStructure, Representation, Scene};

use crate::gpu_cache::{
    atom_arrays, colors_of, keep_mask, rep_sizes, select_atoms, skin_weights, tube_ligands,
    tube_mesh_with_density, without_water, Atoms, GpuCache, Waker, SKIN_MAX_ATOMS,
};

/// Grid spacing a Gaussian surface is meshed at, in angstroms: the
/// vertices then move onto the exact surface, so this only sets how fine
/// its silhouettes are.
const GAUSSIAN_MESH_VOXEL: f32 = 0.5;

/// A render's geometry density and sample count, traded off against
/// time: `High` is the default, chosen to stay close to a `draft`
/// render's time while visibly smoothing cartoon and tube silhouettes
/// over the viewport's own (interactive) tessellation and reducing
/// noise (`light_spread` is exposed for soft shadows but not varied by
/// quality -- see its doc).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum Quality {
    Draft,
    #[default]
    High,
    Ultra,
}

impl Quality {
    pub const ALL: [Quality; 3] = [Quality::Draft, Quality::High, Quality::Ultra];

    pub fn name(self) -> &'static str {
        match self {
            Quality::Draft => "draft",
            Quality::High => "high",
            Quality::Ultra => "ultra",
        }
    }

    pub fn parse(word: &str) -> Option<Quality> {
        Self::ALL.into_iter().find(|q| q.name() == word)
    }

    /// Pixel samples used when the `render` command is given none.
    pub fn default_samples(self) -> u32 {
        match self {
            Quality::Draft => 16,
            Quality::High => 64,
            Quality::Ultra => 160,
        }
    }

    /// Cross-section polygon vertex count for cartoon and tube meshes
    /// (`vv_core::cartoon::RING` is the viewport's fixed 8).
    pub fn ring(self) -> usize {
        match self {
            Quality::Draft => 8,
            Quality::High => 16,
            Quality::Ultra => 24,
        }
    }

    /// Spline samples per residue for cartoon and tube meshes
    /// (`vv_core::cartoon::SAMPLES_PER_RESIDUE` is the viewport's fixed 6).
    pub fn samples_per_residue(self) -> usize {
        match self {
            Quality::Draft => 6,
            Quality::High => 12,
            Quality::Ultra => 18,
        }
    }

    /// Key and fill light cone half-angle, radians. Tested up to 0.15 (7x
    /// [`vv_render::path_trace::DEFAULT_LIGHT_SPREAD`]) on a cartoon
    /// crevice and a ball-and-stick contact shadow at up to 160 spp: no
    /// visible penumbra either way, since a molecule's blockers sit too
    /// close to what they shade for the cone width to matter. Every tier
    /// uses the default until a scene turns up where widening it visibly
    /// helps; the field stays real (`TraceSettings::light_spread`) for
    /// that case or a future UI slider.
    pub fn light_spread(self) -> f32 {
        vv_render::path_trace::DEFAULT_LIGHT_SPREAD
    }
}

/// What `render` asked for.
pub struct RenderRequest {
    pub path: PathBuf,
    pub width: u32,
    pub height: u32,
    pub samples: u32,
    pub transparent: bool,
    pub quality: Quality,
}

/// What a render draws: every visible rep at its structure's current
/// frame, coloured and made as the viewport draws it.
pub struct TraceInput {
    pub scene: TraceScene,
    /// Surfaces, built on the render's worker: they take seconds.
    surfaces: Vec<SurfaceBuild>,
    /// A line for each visible rep left out (a cartoon still building).
    pub skipped: Vec<String>,
}

impl TraceInput {
    pub fn is_empty(&self) -> bool {
        let s = &self.scene;
        s.spheres.is_empty()
            && s.cylinders.is_empty()
            && s.triangles.is_empty()
            && self.surfaces.is_empty()
    }
}

/// A surface rep, from the atoms, radii and colours the viewport's is
/// built from. The SES and skin builds are deterministic, so building
/// them again gives the surface on screen; the Gaussian surface, which
/// the viewport ray-marches in a volume, is meshed.
struct SurfaceBuild {
    kind: Representation,
    positions: Vec<Vec3>,
    radii: Vec<f32>,
    colors: Vec<u32>,
    material: vv_render::Material,
    /// The rep's probe, shrink or blobbiness, whichever its kind has.
    parameter: f32,
}

impl SurfaceBuild {
    fn push_into(self, scene: &mut TraceScene) {
        let (p, r, c, m) = (&self.positions, &self.radii, &self.colors, self.material);
        match self.kind {
            Representation::Ses => {
                let ses = vv_core::ses::build(p, r, self.parameter);
                scene.push_ses(&ses, p, r, c, m);
            }
            Representation::SkinSurface => {
                let weights = skin_weights(r, self.parameter);
                let complex = vv_core::skin_surface::build_complex(p, &weights, self.parameter);
                scene.push_skin(&complex, p, &weights, c, m);
            }
            _ => {
                let mesh = vv_core::gaussian_mesh::mesh(
                    p,
                    r,
                    self.parameter,
                    vv_render::scene::GAUSSIAN_EPSILON,
                    GAUSSIAN_MESH_VOXEL,
                );
                scene.push_mesh(&mesh, c, m);
            }
        }
    }
}

/// A surface rep's probe radius (SES), shrink (skin) or blobbiness
/// (Gaussian), as set or by default.
fn surface_parameter(rep: &vv_scene::Rep) -> f32 {
    let name = match rep.representation {
        Representation::Ses => "probe",
        Representation::SkinSurface => "shrink",
        _ => "blob",
    };
    rep.option(name).expect("every surface has its parameter")
}

/// Everything `scene` shows at `quality`'s geometry density; cartoons and
/// tubes rebuilt denser than the mesh `cache` draws (`CartoonPlan::
/// redensify`, `tube_mesh_with_density`), everything else unaffected by
/// `quality` (already analytic spheres and cylinders, or its own fixed
/// mesh resolution).
pub fn trace_scene(scene: &Scene, cache: &GpuCache, quality: Quality) -> TraceInput {
    let mut traced = TraceScene::default();
    let mut surfaces = Vec::new();
    let mut skipped = Vec::new();
    for (id, loaded) in scene.structures() {
        if !loaded.visible {
            continue;
        }
        let frame = loaded
            .frame
            .min(loaded.structure.frame_count().saturating_sub(1));
        let coords = loaded.structure.frame(frame);
        let positions = coords.positions();
        for (index, rep) in loaded.reps.iter().enumerate() {
            if !rep.visible {
                continue;
            }
            let atoms = match select_atoms(loaded, rep, frame) {
                Ok(atoms) => atoms,
                Err(e) => {
                    skipped.push(e);
                    continue;
                }
            };
            let n = loaded.structure.atom_count();
            let colors = colors_of(loaded, &rep.coloring, frame);
            let material = vv_render::Material::from(rep.material);
            if rep.representation == Representation::Cartoon {
                let Some((plan, spline)) = cache.cartoon(id) else {
                    skipped.push(format!(
                        "{} rep {index} (cartoon still building)",
                        loaded.label
                    ));
                    continue;
                };
                let dense = plan.redensify(quality.samples_per_residue());
                let mesh = match keep_mask(&atoms, n) {
                    None => dense.mesh(&spline),
                    Some(keep) => dense.filter(|a| keep[a as usize]).mesh(&spline),
                };
                traced.push_mesh(&mesh.expand_with_ring(quality.ring()), &colors, material);
                continue;
            }
            if matches!(
                rep.representation,
                Representation::Ses | Representation::SkinSurface | Representation::GaussianSurface
            ) {
                // The SES is of the molecule, not its waters (as drawn).
                let of = match rep.representation {
                    Representation::Ses => without_water(&loaded.structure, &atoms),
                    _ => atoms,
                };
                let (positions, radii, colors) =
                    atom_arrays(&loaded.structure, &of, frame, &colors);
                if rep.representation == Representation::SkinSurface
                    && positions.len() > SKIN_MAX_ATOMS
                {
                    skipped.push(format!(
                        "{} rep {index} (skin surface too large)",
                        loaded.label
                    ));
                    continue;
                }
                let parameter = surface_parameter(rep);
                // The tracer takes one probe radius (and one shrink) per
                // render: the first surface's.
                let first = surfaces
                    .iter()
                    .find(|s: &&SurfaceBuild| s.kind == rep.representation)
                    .map(|s| s.parameter);
                if rep.representation != Representation::GaussianSurface
                    && first.is_some_and(|f| f != parameter)
                {
                    skipped.push(format!(
                        "{} rep {index} (a second probe or shrink value in one render)",
                        loaded.label
                    ));
                    continue;
                }
                surfaces.push(SurfaceBuild {
                    kind: rep.representation,
                    positions,
                    radii,
                    colors,
                    material,
                    parameter,
                });
                continue;
            }
            if rep.representation == Representation::Tube {
                push_tube(
                    &mut traced,
                    loaded,
                    rep,
                    &atoms,
                    frame,
                    &colors,
                    material,
                    quality,
                );
                continue;
            }
            if rep.representation == Representation::Glycan {
                let plan = vv_core::GlycanPlan::build(&loaded.structure.topology, &loaded.bonds);
                let mut glycan_frame = vv_core::GlycanFrame::default();
                plan.update_into(positions, &mut glycan_frame);
                let keep = keep_mask(&atoms, n);
                let keep_residue = |residue: u32| match &keep {
                    None => true,
                    Some(mask) => plan
                        .residues
                        .iter()
                        .find(|res| res.residue == residue)
                        .is_some_and(|res| mask[res.anomeric() as usize]),
                };
                let size = rep.option("size").unwrap_or(4.0);
                let link_radius = rep.option("radius").unwrap_or(0.5);
                let mut mesh = vv_core::PolytopeMesh::default();
                vv_core::build_glycan_mesh(&plan, &glycan_frame, size, keep_residue, &mut mesh);
                vv_core::build_glycan_linkage_mesh(
                    &plan,
                    &glycan_frame,
                    link_radius,
                    keep_residue,
                    &mut mesh,
                );
                traced.push_glycan_mesh(&mesh, material);
                continue;
            }
            let sizes = rep_sizes(rep);
            let radii: Vec<f32> = loaded
                .structure
                .topology
                .element
                .iter()
                .map(|e| sizes.atom_radius(e.vdw_radius()))
                .collect();
            let keep = keep_mask(&atoms, n).unwrap_or_else(|| vec![true; n]);
            // Non-`Single` bonds draw as strands (below) instead of the
            // ordinary single cylinder push_atoms would give them.
            let single_bonds = loaded
                .bonds
                .pairs
                .iter()
                .enumerate()
                .filter(|&(i, &[a, b])| {
                    keep[a as usize]
                        && keep[b as usize]
                        && loaded.bonds.order_of(i as u32) == vv_core::BondOrder::Single
                })
                .map(|(_, &pair)| pair);
            let material_index = traced.materials.len() as u32;
            traced.push_atoms(
                positions,
                &radii,
                &colors,
                (0..n).filter(|&a| keep[a]),
                single_bonds,
                sizes.bond_radius,
                material,
            );
            push_bond_strands(
                &mut traced,
                &loaded.bonds,
                positions,
                &colors,
                &keep,
                sizes.bond_radius,
                material_index,
            );
        }
    }
    TraceInput {
        scene: traced,
        surfaces,
        skipped,
    }
}

/// The extra parallel-strand cylinders of `bonds`'s non-`Single` entries
/// (see `gpu_cache::build_strand_geometry`, which draws the same thing on
/// the raster path from the same `vv_core::bond_strands` data), for a rep
/// whose bonds are drawn at all (`bond_radius > 0`) and only where both
/// atoms pass `keep`. `material_index` is the index `push_atoms` (called
/// just before this) registered the rep's material at.
fn push_bond_strands(
    traced: &mut TraceScene,
    bonds: &vv_core::BondTable,
    positions: &[Vec3],
    colors: &[u32],
    keep: &[bool],
    bond_radius: f32,
    material_index: u32,
) {
    if bond_radius <= 0.0 || bonds.orders.is_empty() {
        return;
    }
    let adjacency = vv_core::adjacency(bonds, positions.len());
    let separation = bond_radius * crate::gpu_cache::STRAND_SEPARATION_SCALE;
    let radius = bond_radius * crate::gpu_cache::STRAND_RADIUS_SCALE;
    for strand in vv_core::bond_strands(bonds, &adjacency, separation) {
        let [a, b] = strand.atoms;
        if !keep[a as usize] || !keep[b as usize] {
            continue;
        }
        let (pa, pb) = vv_core::bond_strand_endpoints(&strand, positions);
        traced.cylinders.push(vv_render::path_trace::TraceCylinder {
            a: pa,
            b: pb,
            radius,
            color_a: colors[a as usize],
            color_b: colors[b as usize],
            material: material_index,
        });
    }
}

/// A tube rep as the viewport draws it, but at `quality`'s density
/// (`gpu_cache::tube_mesh_with_density`), with rounded end caps and
/// ligands as ball-and-stick.
fn push_tube(
    traced: &mut TraceScene,
    loaded: &LoadedStructure,
    rep: &vv_scene::Rep,
    atoms: &Atoms,
    frame: usize,
    colors: &[u32],
    material: vv_render::Material,
    quality: Quality,
) {
    let keep = keep_mask(atoms, loaded.structure.atom_count());
    let (plan, spline) =
        tube_mesh_with_density(loaded, rep, &keep, frame, quality.samples_per_residue());
    let mesh = plan.mesh(&spline);
    traced.push_mesh(&mesh.expand_with_ring(quality.ring()), colors, material);

    let ends = mesh.loose_ends();
    if !ends.is_empty() {
        let coords = loaded.structure.frame(frame);
        let positions = coords.positions();
        let cap_pos: Vec<Vec3> = ends.iter().map(|&(a, _)| positions[a as usize]).collect();
        let cap_radii: Vec<f32> = ends.iter().map(|&(_, r)| r).collect();
        let cap_colors: Vec<u32> = ends.iter().map(|&(a, _)| colors[a as usize]).collect();
        traced.push_atoms(
            &cap_pos,
            &cap_radii,
            &cap_colors,
            0..cap_pos.len(),
            std::iter::empty(),
            0.0,
            material,
        );
    }

    let ligands = tube_ligands(&loaded.structure, &loaded.bonds.pairs, &keep);
    let coords = loaded.structure.frame(frame);
    let sticks = GpuRepresentation::BallAndStick;
    let radii: Vec<f32> = loaded
        .structure
        .topology
        .element
        .iter()
        .map(|e| sticks.atom_radius(e.vdw_radius()))
        .collect();
    traced.push_atoms(
        coords.positions(),
        &radii,
        colors,
        ligands.atoms.iter().map(|&a| a as usize),
        ligands.real_pairs.iter().copied(),
        sticks.bond_radius(),
        material,
    );
}

/// A render running on a worker thread.
pub struct RenderJob {
    pub path: PathBuf,
    /// Fraction done, as f32 bits.
    progress: Arc<AtomicU32>,
    result: mpsc::Receiver<Result<Vec<u8>, String>>,
    pub width: u32,
    pub height: u32,
    pub samples: u32,
    pub started: Instant,
}

impl RenderJob {
    /// Builds the surfaces and the hierarchy and renders `input` through
    /// `camera` on a worker; `waker` asks the window for a frame whenever
    /// a tile is done.
    pub fn start(
        ctx: Arc<GpuContext>,
        input: TraceInput,
        camera: Camera,
        settings: TraceSettings,
        path: PathBuf,
        waker: Waker,
    ) -> Self {
        let progress = Arc::new(AtomicU32::new(0));
        let (tx, result) = mpsc::channel();
        let seen = progress.clone();
        std::thread::spawn(move || {
            let mut scene = input.scene;
            for surface in input.surfaces {
                surface.push_into(&mut scene);
            }
            let image = PathTracer::new(&ctx, &scene)
                .map_err(|e| format!("render: {e}"))
                .map(|tracer| {
                    tracer.render(&ctx, &camera, &settings, |f| {
                        seen.store(f.to_bits(), Ordering::Relaxed);
                        waker();
                    })
                });
            let _ = tx.send(image);
            waker();
        });
        Self {
            path,
            progress,
            result,
            width: settings.width,
            height: settings.height,
            samples: settings.samples,
            started: Instant::now(),
        }
    }

    pub fn progress(&self) -> f32 {
        f32::from_bits(self.progress.load(Ordering::Relaxed))
    }

    /// The finished image (straight-alpha sRGB RGBA8), once there is one.
    pub fn poll(&self) -> Option<Result<Vec<u8>, String>> {
        match self.result.try_recv() {
            Ok(done) => Some(done),
            Err(mpsc::TryRecvError::Empty) => None,
            Err(mpsc::TryRecvError::Disconnected) => Some(Err("render: the worker died".into())),
        }
    }
}
