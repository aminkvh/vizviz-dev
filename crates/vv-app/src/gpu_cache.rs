//! Bridges `vv_scene::Scene` (CPU document model) to the GPU. Every open
//! structure's visible reps are drawn every frame, composited into one
//! image by `Renderer::render_all`; this cache keeps the GPU copies in
//! step with the scene (uploads, recolors, frame changes, drops) without
//! re-uploading anything that hasn't changed.
//!
//! Per structure: its atom geometry (uploaded once, only if some rep
//! draws all of its atoms), and its cartoon mesh (DSSP always runs on the
//! whole structure, never on a selection). Per rep: which atoms its
//! selection picks, its colours, and its style's geometry -- a draw state
//! over the structure's geometry when it draws every atom, else a subset
//! of its own; the tube or cartoon filtered to its atoms; its own
//! surface. Colours are always computed on the whole structure and then
//! picked out, so a selection never changes how an atom is coloured.
//!
//! Derived geometry keeps a map back to real atoms, so a pick on a tube
//! sample, a subset sphere or a cartoon section still names an atom of
//! the structure. Cartoon and skin-surface builds take seconds on
//! multi-million-atom structures, so their CPU halves run on a worker
//! thread (`Job`) and the window is woken to upload them.

pub(crate) mod bases;
pub(crate) mod companions;
pub(crate) mod interactions;
mod surface_mesh;

use std::collections::HashMap;
use std::sync::mpsc;
use std::sync::Arc;

use bases::Bases;
use companions::Companions;
use interactions::Overlay;

use glam::Vec3;
use rayon::prelude::*;
use vv_core::cartoon::{CartoonFrame, CartoonPlan};
use vv_core::skin_surface::SkinComplex;
use vv_core::ResidueClass;
use vv_core::Structure;
use vv_render::renderer::PageBindings;
use vv_render::scene::vdw_radii;
use vv_render::{
    AtomSizes, Camera, CartoonGpu, CartoonItem, DrawItem, DrawState, GaussianSurfaceGpu,
    GaussianSurfaceItem, GpuContext, GpuStructure, OcclusionVolume, OutOfGpuMemory, PatchSurface,
    PatchSurfaceItem, Renderer, Representation as GpuRepresentation, SesBindings, SesGpu,
    SesLayout, SkinSurfaceGpu,
};
use vv_scene::{
    ColorOverride, ColorScheme as SceneColorScheme, LoadedStructure, Rep, RepId,
    Representation as SceneRepresentation, Scene, StructureId, ValueChannel,
};

/// A rep's atom and bond sizes: its style's, with its options applied
/// (`vv_scene::Representation::options`).
pub(crate) fn rep_sizes(rep: &Rep) -> AtomSizes {
    sizes_of(rep.representation, &rep.options)
}

fn sizes_of(
    representation: SceneRepresentation,
    set: &std::collections::BTreeMap<String, f32>,
) -> AtomSizes {
    let mut sizes = AtomSizes::of(map_representation(representation));
    let option = |name| representation.option(set, name);
    match representation {
        SceneRepresentation::Spacefill | SceneRepresentation::BallAndStick => {
            sizes.radius_scale = option("scale").unwrap_or(sizes.radius_scale);
            sizes.bond_radius = option("bond").unwrap_or(sizes.bond_radius);
        }
        SceneRepresentation::Sas => {
            sizes.radius_offset = option("probe").unwrap_or(sizes.radius_offset);
        }
        SceneRepresentation::Sticks => {
            let stick = option("bond").unwrap_or(sizes.bond_radius);
            sizes.radius_offset = stick;
            sizes.bond_radius = stick;
        }
        // Tube's own radius (constant or by B-factor) is baked into its
        // mesh (`gpu_cache::tube_radius_fn`), not read from `AtomSizes`;
        // its end caps and ligand sticks size themselves.
        _ => {}
    }
    sizes
}

pub fn map_representation(r: SceneRepresentation) -> GpuRepresentation {
    match r {
        SceneRepresentation::Spacefill => GpuRepresentation::Spacefill,
        SceneRepresentation::BallAndStick => GpuRepresentation::BallAndStick,
        SceneRepresentation::Tube => GpuRepresentation::Tube,
        SceneRepresentation::Sas => GpuRepresentation::Sas,
        SceneRepresentation::Sticks => GpuRepresentation::Sticks,
        SceneRepresentation::Lines => GpuRepresentation::Lines,
        // Drawn through the cartoon and surface lists, never as a
        // `DrawItem`; `Spacefill` only keeps the match exhaustive.
        SceneRepresentation::Cartoon => GpuRepresentation::Spacefill,
        SceneRepresentation::GaussianSurface => GpuRepresentation::Spacefill,
        SceneRepresentation::SkinSurface => GpuRepresentation::Spacefill,
        SceneRepresentation::Ses => GpuRepresentation::Spacefill,
        SceneRepresentation::Glycan => GpuRepresentation::Spacefill,
    }
}

pub use vv_render::coloring::colors_of;

/// What a rep's colours were last computed from, so `sync` only rewrites
/// them when that changed: the scheme, the channel it names (by identity,
/// so re-attaching a channel under the same name recolors), and the frame
/// when the colours are per-frame.
#[derive(Clone)]
struct ColorKey {
    coloring: SceneColorScheme,
    overrides: Vec<ColorOverride>,
    channel: Option<ValueChannel>,
    frame: usize,
}

impl ColorKey {
    fn of(loaded: &LoadedStructure, coloring: &SceneColorScheme, frame: usize) -> Self {
        let channel = loaded.values_for(coloring).map(|(_, c)| c.clone());
        // Per-frame colors: a per-frame channel, secondary structure on
        // a trajectory (DSSP on every frame), or overrides whose
        // selections can move with the atoms.
        let multi_frame = loaded.structure.frame_count() > 1;
        let frame = match &channel {
            Some(c) if c.frames() > 1 => frame,
            None if multi_frame
                && (*coloring == SceneColorScheme::SecondaryStructure
                    || !loaded.color_overrides.is_empty()) =>
            {
                frame
            }
            _ => 0,
        };
        Self {
            coloring: coloring.clone(),
            overrides: loaded.color_overrides.clone(),
            channel,
            frame,
        }
    }

    fn same(&self, other: &Self) -> bool {
        self.coloring == other.coloring
            && self.overrides == other.overrides
            && self.frame == other.frame
            && match (&self.channel, &other.channel) {
                (None, None) => true,
                (Some(a), Some(b)) => std::sync::Arc::ptr_eq(a.data(), b.data()),
                _ => false,
            }
    }
}

/// The structure ribbon actions apply to when they need one: the
/// most recently loaded structure still open.
pub fn current_structure(scene: &Scene) -> Option<StructureId> {
    scene.structures().next_back().map(|(id, _)| id)
}

/// A camera pivoting on `id`'s own centre: the mean of its current
/// frame's atom positions, not `bounding_sphere`'s
/// AABB midpoint, which can sit well clear of an unevenly shaped
/// molecule (a branched glycan, a compact core with a floppy tail) --
/// sized to its farthest atom from that centre. `None` if `id` isn't
/// loaded or has no atoms.
pub fn frame_structure(scene: &Scene, id: StructureId) -> Option<Camera> {
    let loaded = scene.structure(id)?;
    let frame = loaded.structure.frame(loaded.frame);
    let center = frame.centroid()?;
    let radius = frame
        .positions()
        .iter()
        .map(|p| p.distance(center))
        .fold(0.0f32, f32::max);
    Some(Camera::framing(center, radius))
}

/// Where one drawn item's picks map back to. `atom_map[i]` is the real
/// atom behind local sphere `i` (identity when `None`); `bond_atoms[k]`
/// is the real atom pair behind local cylinder `k` (the structure's own
/// bond table when `None`).
#[derive(Clone, Debug)]
pub struct DrawSource {
    pub id: StructureId,
    /// The rep this source draws, so a Structures-panel row can highlight
    /// exactly its own atoms in the viewport (`selection_outline`).
    pub rep: RepId,
    pub atom_map: Option<Arc<Vec<u32>>>,
    pub bond_atoms: Option<Arc<Vec<[u32; 2]>>>,
}

impl DrawSource {
    pub fn atom(&self, local: u32) -> u32 {
        self.atom_map
            .as_ref()
            .and_then(|m| m.get(local as usize).copied())
            .unwrap_or(local)
    }
}

/// The atoms `keep` marks, in order, and the bonds between two of them:
/// re-indexed into that list, and as real atom pairs.
fn subset_of(keep: &[bool], bonds: &[[u32; 2]]) -> LigandSubset {
    let atoms: Vec<u32> = (0..keep.len() as u32)
        .filter(|&a| keep[a as usize])
        .collect();
    let mut local = vec![u32::MAX; keep.len()];
    for (i, &a) in atoms.iter().enumerate() {
        local[a as usize] = i as u32;
    }
    let (local_bonds, real_pairs): (Vec<[u32; 2]>, Vec<[u32; 2]>) = bonds
        .iter()
        .filter(|[a, b]| keep[*a as usize] && keep[*b as usize])
        .map(|&[a, b]| ([local[a as usize], local[b as usize]], [a, b]))
        .unzip();
    LigandSubset {
        atoms,
        local_bonds,
        real_pairs,
    }
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct LigandSubset {
    pub atoms: Vec<u32>,
    pub local_bonds: Vec<[u32; 2]>,
    pub real_pairs: Vec<[u32; 2]>,
}

/// Perpendicular gap between a multi-order bond's strands, and each
/// strand's own cylinder radius, both proportional to the rep's ordinary
/// bond radius so a thin (ball-and-stick) and a thick (licorice) style
/// each draw proportioned strands.
pub(crate) const STRAND_SEPARATION_SCALE: f32 = 2.4;
pub(crate) const STRAND_RADIUS_SCALE: f32 = 0.5;

/// Extra parallel-strand cylinders for a rep's non-`Single` bonds
/// (`vv_core::bond_strands`): its own small `GpuStructure` of virtual,
/// offset endpoint spheres joined pairwise into strands, rebuilt from
/// live atom positions every frame rather than cached, so a played
/// trajectory's strands track their bond and ring plane exactly
/// (`vv_core::bond_strand_endpoints`'s own doc). `bond_atoms[k]` is the
/// real atom pair behind local cylinder `k`, for `DrawSource` picking --
/// the same contract `Derived`/`LigandSubset::real_pairs` use.
///
/// The endpoints' own radius is `0`: a screenshot check across ball-and-
/// stick, licorice and lines found no stray dot at any of them (a
/// zero-radius atom's impostor billboard collapses to a degenerate
/// quad, `shaders/draw.wgsl`'s `vs_sphere`; whether that fully explains
/// every case, versus the point sitting behind the coincident strand
/// cylinder it ends, wasn't traced further -- flag it if a future
/// offset/geometry combination ever puts an endpoint somewhere no other
/// opaque geometry covers it).
struct StrandGeometry {
    gpu: GpuStructure,
    bindings: PageBindings,
    strands: Vec<vv_core::BondStrand>,
    bond_atoms: Arc<Vec<[u32; 2]>>,
}

/// Interleaved `(a, b)` endpoint positions of every strand, at `positions`.
fn strand_positions(strands: &[vv_core::BondStrand], positions: &[Vec3]) -> Vec<Vec3> {
    strands
        .iter()
        .flat_map(|s| {
            let (a, b) = vv_core::bond_strand_endpoints(s, positions);
            [a, b]
        })
        .collect()
}

/// Interleaved `(colors[a], colors[b])` per strand: split coloring at the
/// midpoint, matching how the shader colors an ordinary bond cylinder.
fn strand_colors(strands: &[vv_core::BondStrand], all: &[u32]) -> Vec<u32> {
    strands
        .iter()
        .flat_map(|s| [all[s.atoms[0] as usize], all[s.atoms[1] as usize]])
        .collect()
}

/// `bonds.pairs` with every non-`Single` bond removed -- what the
/// ordinary center-cylinder bond list draws instead of the full table,
/// so a double/triple/aromatic bond's strands (`build_strand_geometry`)
/// are the only cylinders it draws, not one too many.
fn single_bond_pairs(bonds: &vv_core::BondTable) -> Vec<[u32; 2]> {
    bonds
        .pairs
        .iter()
        .enumerate()
        .filter(|&(i, _)| bonds.order_of(i as u32) == vv_core::BondOrder::Single)
        .map(|(_, &pair)| pair)
        .collect()
}

/// `strands`'s endpoints as their own small `GpuStructure` -- the core
/// both `build_strand_geometry` (every non-`Single` bond) and `Derived::
/// atoms` (only those within one subset) upload through.
fn strand_geometry_from(
    ctx: &GpuContext,
    renderer: &Renderer,
    strands: Vec<vv_core::BondStrand>,
    positions: &[Vec3],
    colors: &[u32],
) -> StrandGeometry {
    let pos = strand_positions(&strands, positions);
    let col = strand_colors(&strands, colors);
    let radii = vec![0.0f32; pos.len()];
    let local_bonds: Vec<[u32; 2]> = (0..strands.len() as u32)
        .map(|k| [2 * k, 2 * k + 1])
        .collect();
    let gpu = GpuStructure::from_parts(ctx, &pos, &radii, &col, &local_bonds);
    let bindings = renderer.bind(&gpu);
    let bond_atoms = Arc::new(strands.iter().map(|s| s.atoms).collect());
    StrandGeometry {
        gpu,
        bindings,
        strands,
        bond_atoms,
    }
}

/// `bonds`'s strands (`separation` apart), or `None` when there are none
/// to draw (the common all-`Single` case costs nothing beyond the empty
/// check).
fn build_strand_geometry(
    ctx: &GpuContext,
    renderer: &Renderer,
    bonds: &vv_core::BondTable,
    adjacency: &vv_core::Adjacency,
    positions: &[Vec3],
    separation: f32,
    colors: &[u32],
) -> Option<StrandGeometry> {
    if bonds.orders.is_empty() {
        return None;
    }
    let strands = vv_core::bond_strands(bonds, adjacency, separation);
    Some(strand_geometry_from(
        ctx, renderer, strands, positions, colors,
    ))
}

/// What a `Derived` piece needs to also draw its non-`Single` bonds as
/// strands: the structure-wide bond table and adjacency (both already
/// built once and cached per structure, `Entry::adjacency`) and this
/// rep's strand separation. A piece with no real chemistry bonds of its
/// own (a tube's end caps, a nucleic ladder's synthetic rungs) passes
/// `None` instead rather than paying for a lookup that can only ever
/// come up empty.
struct StrandInput<'a> {
    bonds: &'a vv_core::BondTable,
    adjacency: &'a vv_core::Adjacency,
    separation: f32,
}

/// Whether each of `real_pairs` (e.g. a `LigandSubset::real_pairs`) is a
/// `Single` bond in `bonds`. A pair `bonds` doesn't know at all (never
/// true for a real subset's own bonds, since they come from `bonds`
/// itself) counts as `Single` too, so it still draws its ordinary
/// cylinder rather than silently vanishing.
fn is_single(bonds: &vv_core::BondTable, real_pairs: &[[u32; 2]]) -> Vec<bool> {
    real_pairs
        .iter()
        .map(|&pair| {
            bonds.pairs.binary_search(&pair).map_or(true, |i| {
                bonds.order_of(i as u32) == vv_core::BondOrder::Single
            })
        })
        .collect()
}

/// `strands` restricted to those whose real atom pair is one of `sorted_
/// real_pairs` (a subset's own bonds, e.g. `LigandSubset::real_pairs`,
/// which stays sorted -- a filtered subsequence of the sorted `BondTable`
/// it was built from).
fn strands_within(
    strands: Vec<vv_core::BondStrand>,
    sorted_real_pairs: &[[u32; 2]],
) -> Vec<vv_core::BondStrand> {
    strands
        .into_iter()
        .filter(|s| sorted_real_pairs.binary_search(&s.atoms).is_ok())
        .collect()
}

impl StrandGeometry {
    fn set_frame(&self, ctx: &GpuContext, positions: &[Vec3]) {
        let pos = strand_positions(&self.strands, positions);
        let radii = vec![0.0f32; pos.len()];
        self.gpu.set_instances(ctx, &pos, &radii);
    }

    fn recolor(&self, ctx: &GpuContext, all_colors: &[u32]) {
        self.gpu
            .set_colors(ctx, &strand_colors(&self.strands, all_colors));
    }
}

/// Spheres and cylinders of their own for some of a structure's atoms (a
/// selection, the tube's samples, the ligand sticks), with the map back.
struct Derived {
    gpu: GpuStructure,
    bindings: PageBindings,
    atom_map: Arc<Vec<u32>>,
    bond_atoms: Arc<Vec<[u32; 2]>>,
    /// This piece's own non-`Single` bonds, drawn as strands instead of
    /// (not alongside) their ordinary cylinder; `None` when `strand_input`
    /// was `None` or none of `subset`'s bonds are non-`Single`.
    strands: Option<StrandGeometry>,
}

impl Derived {
    /// This piece's colors picked out of the structure-wide `all`.
    fn colors(&self, all: &[u32]) -> Vec<u32> {
        self.atom_map.iter().map(|&a| all[a as usize]).collect()
    }

    /// `subset`'s atoms of `structure` at `frame`, van der Waals radii,
    /// and -- given `strand_input` -- its own non-`Single` bonds split
    /// into strands the same way `AllAtoms`'s shared geometry is.
    fn atoms(
        ctx: &GpuContext,
        renderer: &Renderer,
        structure: &Structure,
        frame: usize,
        subset: LigandSubset,
        all_colors: &[u32],
        strand_input: Option<StrandInput>,
    ) -> Self {
        let coords = structure.frame(frame);
        let positions = coords.positions();
        let LigandSubset {
            atoms,
            local_bonds,
            real_pairs,
        } = subset;

        let (local_bonds, real_pairs, strands) = match strand_input {
            Some(si) => {
                let single = is_single(si.bonds, &real_pairs);
                let (kept_local, kept_real): (Vec<_>, Vec<_>) = local_bonds
                    .into_iter()
                    .zip(&real_pairs)
                    .zip(single)
                    .filter(|&(_, keep)| keep)
                    .map(|((l, &r), _)| (l, r))
                    .unzip();
                let subset_strands = strands_within(
                    vv_core::bond_strands(si.bonds, si.adjacency, si.separation),
                    &real_pairs,
                );
                let strands = (!subset_strands.is_empty()).then(|| {
                    strand_geometry_from(ctx, renderer, subset_strands, positions, all_colors)
                });
                (kept_local, kept_real, strands)
            }
            None => (local_bonds, real_pairs, None),
        };

        let pos: Vec<Vec3> = atoms.iter().map(|&a| positions[a as usize]).collect();
        let radii: Vec<f32> = atoms
            .iter()
            .map(|&a| structure.topology.element[a as usize].vdw_radius())
            .collect();
        let colors: Vec<u32> = atoms.iter().map(|&a| all_colors[a as usize]).collect();
        let gpu = GpuStructure::from_parts(ctx, &pos, &radii, &colors, &local_bonds);
        Derived {
            bindings: renderer.bind(&gpu),
            gpu,
            atom_map: Arc::new(atoms),
            bond_atoms: Arc::new(real_pairs),
            strands,
        }
    }

    /// Moves this subset's atoms, and its strands, to `frame`.
    fn set_frame(&self, ctx: &GpuContext, structure: &Structure, frame: usize) {
        let coords = structure.frame(frame);
        let positions = coords.positions();
        let pos: Vec<Vec3> = self
            .atom_map
            .iter()
            .map(|&a| positions[a as usize])
            .collect();
        let radii: Vec<f32> = self
            .atom_map
            .iter()
            .map(|&a| structure.topology.element[a as usize].vdw_radius())
            .collect();
        self.gpu.set_instances(ctx, &pos, &radii);
        if let Some(s) = &self.strands {
            s.set_frame(ctx, positions);
        }
    }

    /// Recolors this subset and its strands from the structure-wide `all`.
    fn recolor(&self, ctx: &GpuContext, all: &[u32]) {
        self.gpu.set_colors(ctx, &self.colors(all));
        if let Some(s) = &self.strands {
            s.recolor(ctx, all);
        }
    }
}

/// A tube rep as the viewport draws it: a round-cross-section mesh
/// through the same GPU path as the ribbon cartoon (`gpu_cache::
/// tube_mesh`), rounded end caps, and ligand sticks.
struct TubeGeometry {
    /// The rep's selection, already filtered in: fixed section/join/span
    /// counts for the geometry's lifetime (B-factor doesn't change across
    /// frames), only the spline moves per frame.
    plan: CartoonPlan,
    gpu: CartoonGpu,
    bindings: vv_render::CartoonBindings,
    /// Rounded end caps (`vv_core::cartoon::CartoonMesh::loose_ends`), and
    /// their radii: fixed, like `plan`'s.
    caps: Option<Derived>,
    cap_radii: Vec<f32>,
    bases: Bases,
    companions: Companions,
}

/// A tube rep's per-atom radius (`vv_core::backbone::putty_radius`):
/// `radius` for every atom, or (`radius_by` >= 0.5) linear from
/// `radius_min` to the `putty`-scaled thick end over the B-factor range of `source` -- the
/// tube's own drawn atoms, not the whole structure, so a selection reads
/// by its own contrast.
fn tube_radius_fn<'a>(
    topology: &'a vv_core::Topology,
    source: &[u32],
    rep: &Rep,
) -> impl Fn(u32) -> f32 + 'a {
    let radius = rep.option("radius").unwrap_or(vv_core::TUBE_RADIUS);
    let putty = rep.option("radius_by").unwrap_or(0.0) >= 0.5;
    let radius_min = rep.option("radius_min").unwrap_or(radius);
    let strength = rep.option("putty").unwrap_or(1.0);
    let thick = radius_min + strength * (radius - radius_min);
    let b_range = if putty {
        let b_factor: Vec<f32> = source
            .iter()
            .map(|&a| topology.b_factor.get(a as usize).copied().unwrap_or(0.0))
            .collect();
        vv_render::scene::scalar_range(&b_factor)
    } else {
        (0.0, 0.0)
    };
    move |a: u32| {
        if !putty {
            return radius;
        }
        let b = topology.b_factor.get(a as usize).copied().unwrap_or(0.0);
        vv_core::backbone::putty_radius(b, b_range, (radius_min, thick))
    }
}

/// The tube's plan and spline at `frame`: the whole backbone trace,
/// radius by `tube_radius_fn`, filtered to `keep`. Shared by the viewport
/// (`build_tube`, below) and the path tracer (`render::push_tube`), so
/// both draw the same tube, at the viewport's section density.
pub(crate) fn tube_mesh(
    loaded: &LoadedStructure,
    rep: &Rep,
    keep: &Option<Vec<bool>>,
    frame: usize,
) -> (CartoonPlan, CartoonFrame) {
    tube_mesh_with_density(
        loaded,
        rep,
        keep,
        frame,
        vv_core::cartoon::SAMPLES_PER_RESIDUE,
    )
}

/// As [`tube_mesh`], with the spline samples per residue given explicitly
/// -- the path tracer's own, denser tube, rebuilt fresh (no DSSP, so this
/// is cheap unlike a ribbon cartoon's [`vv_core::cartoon::CartoonPlan::
/// redensify`]).
pub(crate) fn tube_mesh_with_density(
    loaded: &LoadedStructure,
    rep: &Rep,
    keep: &Option<Vec<bool>>,
    frame: usize,
    samples_per_residue: usize,
) -> (CartoonPlan, CartoonFrame) {
    let structure = &loaded.structure;
    let trace = vv_core::backbone_trace(&structure.topology, loaded.drawn_positions(0).positions());
    let keep = &loaded.trace_keep(keep);
    let all_trace_atoms = trace.segments.iter().flatten().copied();
    let source: Vec<u32> = match keep {
        None => all_trace_atoms.collect(),
        Some(keep) => all_trace_atoms.filter(|&a| keep[a as usize]).collect(),
    };
    let radius = tube_radius_fn(&structure.topology, &source, rep);
    let mut plan = vv_core::cartoon::tube_plan_with_density(&trace, radius, samples_per_residue);
    if let Some(keep) = keep {
        plan = plan.filter(|a| keep[a as usize]);
    }
    let spline = plan.frame(loaded.drawn_positions(frame).positions());
    (plan, spline)
}

/// Rounded-cap spheres for a tube's mesh (`CartoonMesh::loose_ends`): one
/// per open end, at its own radius. `None` when the mesh has no sections
/// (an empty tube).
fn build_tube_caps(
    ctx: &GpuContext,
    renderer: &Renderer,
    structure: &Structure,
    frame: usize,
    ends: &[(u32, f32)],
    all_colors: &[u32],
) -> (Option<Derived>, Vec<f32>) {
    if ends.is_empty() {
        return (None, Vec::new());
    }
    let atoms: Vec<u32> = ends.iter().map(|&(a, _)| a).collect();
    let radii: Vec<f32> = ends.iter().map(|&(_, r)| r).collect();
    let coords = structure.frame(frame);
    let positions = coords.positions();
    let pos: Vec<Vec3> = atoms.iter().map(|&a| positions[a as usize]).collect();
    let colors: Vec<u32> = atoms.iter().map(|&a| all_colors[a as usize]).collect();
    let gpu = GpuStructure::from_parts(ctx, &pos, &radii, &colors, &[]);
    let derived = Derived {
        bindings: renderer.bind(&gpu),
        gpu,
        atom_map: Arc::new(atoms),
        bond_atoms: Arc::new(Vec::new()),
        strands: None,
    };
    (Some(derived), radii)
}

#[allow(clippy::too_many_arguments)]
fn build_tube(
    ctx: &GpuContext,
    renderer: &Renderer,
    loaded: &LoadedStructure,
    rep: &Rep,
    keep: Option<Vec<bool>>,
    all_colors: &[u32],
    frame: usize,
    adjacency: Option<&vv_core::Adjacency>,
) -> Result<TubeGeometry, OutOfGpuMemory> {
    let structure = &loaded.structure;
    let (plan, spline) = tube_mesh(loaded, rep, &keep, frame);
    let gpu = CartoonGpu::upload(ctx, &plan, &spline, all_colors)?;
    let bindings = renderer.bind_cartoon(&gpu);
    let ends = plan.mesh(&spline).loose_ends();
    let (caps, cap_radii) = build_tube_caps(ctx, renderer, structure, frame, &ends, all_colors);
    let style = bases::BaseStyle::of(rep);
    let bases = Bases::build(ctx, renderer, loaded, &keep, frame, all_colors, style)?;
    let traced = companions::traced_atoms(&plan);
    let drawn = companions::Drawn {
        keep: &keep,
        traced: &traced,
    };
    let companions = Companions::build(
        ctx, renderer, loaded, rep, drawn, frame, all_colors, adjacency,
    );
    Ok(TubeGeometry {
        plan,
        gpu,
        bindings,
        caps,
        cap_radii,
        bases,
        companions,
    })
}

impl TubeGeometry {
    fn set_frame(
        &mut self,
        ctx: &GpuContext,
        renderer: &Renderer,
        loaded: &LoadedStructure,
        bases_of: (&Option<Vec<bool>>, bases::BaseStyle, &[u32]),
        frame: usize,
    ) {
        let structure = &loaded.structure;
        let spline = self.plan.frame(loaded.drawn_positions(frame).positions());
        self.gpu.set_frame(ctx, &spline);
        if let Some(caps) = &self.caps {
            let coords = structure.frame(frame);
            let positions = coords.positions();
            let pos: Vec<Vec3> = caps
                .atom_map
                .iter()
                .map(|&a| positions[a as usize])
                .collect();
            caps.gpu.set_instances(ctx, &pos, &self.cap_radii);
        }
        let (keep, style, colors) = bases_of;
        self.bases =
            Bases::build(ctx, renderer, loaded, keep, frame, colors, style).unwrap_or_default();
        self.companions.set_frame(ctx, structure, frame);
    }
}

fn oom_notice(what: &str, label: &str) -> String {
    format!("{label}: not enough GPU memory for the {what}; try `rep lines` or a smaller selection")
}

/// A CPU build running on a worker thread for coordinate set `frame`.
struct Job<T> {
    frame: usize,
    result: mpsc::Receiver<T>,
}

/// What finishing a job means for the window: it may be idle
/// (`ControlFlow::Wait`), so the worker asks it for a frame.
pub type Waker = Arc<dyn Fn() + Send + Sync>;

fn spawn_job<T: Send + 'static>(
    frame: usize,
    waker: Option<Waker>,
    work: impl FnOnce() -> T + Send + 'static,
) -> Job<T> {
    let (tx, result) = mpsc::channel();
    std::thread::spawn(move || {
        // The receiver is gone if the structure closed meanwhile.
        let _ = tx.send(work());
        if let Some(wake) = waker {
            wake();
        }
    });
    Job { frame, result }
}

/// A finished job's result; `Err(())` if the worker died (panicked).
fn poll_job<T>(job: &mut Option<Job<T>>) -> Option<Result<(usize, T), ()>> {
    let j = job.as_ref()?;
    match j.result.try_recv() {
        Ok(value) => {
            let frame = j.frame;
            *job = None;
            Some(Ok((frame, value)))
        }
        Err(mpsc::TryRecvError::Empty) => None,
        Err(mpsc::TryRecvError::Disconnected) => {
            *job = None;
            Some(Err(()))
        }
    }
}

/// Above this, a skin surface is refused with a notice: the CPU build is
/// ~50 us/atom (minutes at 4M atoms) and drawing its patches drops well
/// below 30 fps up close.
/// The Gaussian surface handles these sizes.
pub const SKIN_MAX_ATOMS: usize = 200_000;

struct GaussianSurfaceEntry {
    gpu: GaussianSurfaceGpu,
    bindings: wgpu::BindGroup,
    /// The coarse volume a playing trajectory re-bakes every frame
    /// (`GAUSSIAN_PLAYBACK_VOXELS`); the full one when paused.
    playback: bool,
}

struct SesEntry {
    gpu: SesGpu,
    bindings: SesBindings,
    /// The rep's atoms without water, which the surface is of.
    atoms: Atoms,
    frame: usize,
}

/// `atoms` (every atom when `None`) without water: a molecular surface is
/// of the molecule, not the crystal waters on it.
pub(crate) fn without_water(structure: &Structure, atoms: &Atoms) -> Atoms {
    let t = &structure.topology;
    let mut water = vec![false; t.atom_count()];
    for (r, rec) in t.residues.iter().enumerate() {
        if t.residue_class(r) == ResidueClass::Water {
            for a in rec.atoms.clone() {
                water[a as usize] = true;
            }
        }
    }
    if !water.contains(&true) {
        return atoms.clone();
    }
    let kept: Vec<u32> = match atoms {
        None => (0..water.len() as u32)
            .filter(|&a| !water[a as usize])
            .collect(),
        Some(list) => list
            .iter()
            .copied()
            .filter(|&a| !water[a as usize])
            .collect(),
    };
    Some(Arc::new(kept))
}

/// Runs `f`, printing its wall time under `VIZVIZ_TIMING` (`state.rs`'s
/// startup-timing convention) for profiling a background surface build.
fn timed<T>(what: &str, atoms: usize, f: impl FnOnce() -> T) -> T {
    if std::env::var_os("VIZVIZ_TIMING").is_none() {
        return f();
    }
    let t0 = std::time::Instant::now();
    let result = f();
    eprintln!(
        "{what}: {atoms} atoms in {:.0} ms",
        t0.elapsed().as_secs_f64() * 1e3
    );
    result
}

/// What a finished SES job hands back: the upload (a second or more of
/// copying at millions of patches, so it happens on the worker too) and
/// the colours it was made with.
type SesBuilt = (Result<SesGpu, OutOfGpuMemory>, Vec<u32>);

/// Builds the SES of `positions`, lays it out and uploads it, on a
/// worker.
fn spawn_ses(
    ctx: Arc<GpuContext>,
    positions: Vec<Vec3>,
    radii: Vec<f32>,
    colors: Vec<u32>,
    probe: f32,
    frame: usize,
    waker: Option<Waker>,
) -> Job<SesBuilt> {
    let count = positions.len();
    spawn_job(frame, waker, move || {
        let ses = timed("SES build", count, || {
            vv_core::ses::build(&positions, &radii, probe)
        });
        let layout = timed("SES layout", count, || {
            SesLayout::new(&ses, &positions, &radii)
        });
        let gpu = timed("SES gpu upload", count, || {
            SesGpu::upload(&ctx, &layout, &colors)
        });
        (gpu, colors)
    })
}

struct SkinSurfaceEntry {
    gpu: SkinSurfaceGpu,
    bindings: wgpu::BindGroup,
    /// The rep's atoms this was built from: a recolor or a pick must use
    /// these, not `RepEntry::atoms`, which a newer selection may already
    /// have replaced while this (carried-over, still showing) surface
    /// rebuilds.
    atoms: Atoms,
    frame: usize,
}

struct CartoonEntry {
    gpu: CartoonGpu,
    bindings: vv_render::CartoonBindings,
    /// The structure's plan this one was filtered from, and the spline
    /// (frame) it shows.
    from: Arc<CartoonPlan>,
    spline: Arc<CartoonFrame>,
    /// With cylinder helices, `from` restyled: the plan whose frame moves
    /// the axes with the atoms.
    cylinders: Option<CartoonPlan>,
    bases: Bases,
    companions: Companions,
}

/// A glycan rep's SNFG shapes and linkage cylinders: unlike a cartoon,
/// there is no worker-thread build to wait on (detection and linkage
/// resolution run once, synchronously -- a few hundred residues against
/// an already-perceived bond table) and no GPU-side rebuild -- the whole
/// mesh is regenerated on the CPU and re-uploaded whenever the frame
/// moves or `size`/`radius` changes; see [`vv_core::glycan`].
struct GlycanEntry {
    gpu: vv_render::GlycanGpu,
    bindings: vv_render::CartoonBindings,
    plan: Arc<vv_core::GlycanPlan>,
    /// Reused across frames (`GlycanPlan::update_into`): after the first
    /// frame sizes its buffers, replaying a trajectory allocates nothing
    /// here.
    frame_buf: vv_core::GlycanFrame,
    size: f32,
    radius: f32,
}

/// A structure's cartoon: planned (secondary structure and all) at one
/// frame, its spline at another. A playing trajectory moves the spline
/// every frame and keeps the plan; paused, the plan catches up on a
/// worker (DSSP for that frame).
struct CartoonModel {
    plan: Arc<CartoonPlan>,
    planned_at: usize,
    spline: Arc<CartoonFrame>,
    spline_at: usize,
}

/// A rep's atoms: all of them (`None`), or the listed ones.
pub(crate) type Atoms = Option<Arc<Vec<u32>>>;

/// Positions, van der Waals radii and colours of `atoms` at `frame`.
pub(crate) fn atom_arrays(
    structure: &Structure,
    atoms: &Atoms,
    frame: usize,
    all_colors: &[u32],
) -> (Vec<Vec3>, Vec<f32>, Vec<u32>) {
    let coords = structure.frame(frame);
    let positions = coords.positions();
    let pick = |a: usize| {
        (
            positions[a],
            structure.topology.element[a].vdw_radius(),
            all_colors[a],
        )
    };
    let rows: Vec<(Vec3, f32, u32)> = match atoms {
        None => (0..positions.len()).map(pick).collect(),
        Some(list) => list.iter().map(|&a| pick(a as usize)).collect(),
    };
    (
        rows.iter().map(|r| r.0).collect(),
        rows.iter().map(|r| r.1).collect(),
        rows.iter().map(|r| r.2).collect(),
    )
}

/// The positions of `atoms` (all when `None`) at `frame`.
fn atom_positions(structure: &Structure, atoms: &Atoms, frame: usize) -> Vec<Vec3> {
    use rayon::prelude::*;
    let coords = structure.frame(frame);
    let positions = coords.positions();
    match atoms {
        None => positions.to_vec(),
        Some(list) => list.par_iter().map(|&a| positions[a as usize]).collect(),
    }
}

/// Skin-surface weights for van der Waals `radii` (`vv_core::skin_surface::
/// weight_for_radius`, Chavent, Levy & Maigret 2008 eq. 1): an isolated
/// atom's patch stays at its true radius for any shrink.
pub(crate) fn skin_weights(radii: &[f32], shrink: f32) -> Vec<f32> {
    radii
        .iter()
        .map(|&r| vv_core::skin_surface::weight_for_radius(r, shrink))
        .collect()
}

/// Starts building the mixed complex of `positions` on a worker.
fn spawn_skin(
    positions: Vec<Vec3>,
    radii: &[f32],
    shrink: f32,
    frame: usize,
    waker: Option<Waker>,
) -> Job<SkinComplex> {
    let weights = skin_weights(radii, shrink);
    let count = positions.len();
    spawn_job(frame, waker, move || {
        timed("skin surface build", count, || {
            vv_core::skin_surface::build_complex(&positions, &weights, shrink)
        })
    })
}

/// Starts DSSP + the cartoon plan for the whole `structure` at `frame`
/// on a worker. Secondary structure: the file's records for a single
/// structure, DSSP at that frame for a trajectory.
fn spawn_cartoon(
    structure: &Structure,
    positions: Vec<Vec3>,
    frame: usize,
    waker: Option<Waker>,
) -> Job<(CartoonPlan, CartoonFrame)> {
    let structure = structure.clone();
    spawn_job(frame, waker, move || plan_cartoon(&structure, &positions))
}

fn plan_cartoon(structure: &Structure, positions: &[Vec3]) -> (CartoonPlan, CartoonFrame) {
    let single_frame = structure.frame_count() == 1;
    let codes = vv_core::cartoon::secondary_structure(&structure.topology, positions, single_frame);
    let plan = vv_core::cartoon::plan(&structure.topology, positions, &codes);
    let spline = plan.frame(positions);
    (plan, spline)
}

/// The cartoon `rep` draws when it asks for cylinder helices (`repopt
/// helix cylinder`): the structure's `plan` and `spline` with each helix
/// straightened, sections `samples` to a residue step. `None` for ribbon
/// helices, which draw the structure's own.
pub(crate) fn cylinder_cartoon(
    rep: &Rep,
    plan: &CartoonPlan,
    spline: &CartoonFrame,
    positions: &[Vec3],
    samples: usize,
) -> Option<(CartoonPlan, CartoonFrame)> {
    if rep.option("helix").unwrap_or(0.0) < 0.5 {
        return None;
    }
    let styled = plan.with_cylinder_helices(spline, samples);
    let frame = styled.frame(positions);
    Some((styled, frame))
}

/// What a rep draws, by style.
enum RepGeometry {
    /// Not built yet (or waiting on a worker).
    Pending,
    /// Built: the selection gives this style nothing to draw.
    Empty,
    /// Every atom: a draw state over the structure's own geometry.
    AllAtoms {
        state: DrawState,
        bindings: PageBindings,
        /// Extra parallel-strand cylinders for double/triple/aromatic
        /// bonds, drawn instead of (not alongside) their ordinary single
        /// cylinder; `None` when the rep draws no bonds (Spacefill, Sas)
        /// or the structure has none. See `build_strand_geometry`.
        strands: Option<StrandGeometry>,
    },
    /// Some atoms, as their own spheres and cylinders.
    SomeAtoms(Derived),
    Tube(Box<TubeGeometry>),
    Cartoon(Box<CartoonEntry>),
    Gaussian(GaussianSurfaceEntry),
    Skin(SkinSurfaceEntry),
    Ses(SesEntry),
    Glycan(GlycanEntry),
    /// A surface as its mesh edges or vertices (`repopt surface mesh|dots`).
    Mesh(surface_mesh::SurfaceMesh),
}

struct RepEntry {
    representation: SceneRepresentation,
    /// `Rep::options` the geometry was built with.
    options: std::collections::BTreeMap<String, f32>,
    /// The selection text and, for a geometric selection (`within`), the
    /// frame it was evaluated at.
    selection: (String, usize),
    atoms: Atoms,
    colors: ColorKey,
    /// The coordinate set the geometry holds.
    frame: usize,
    geometry: RepGeometry,
    skin_job: Option<Job<SkinComplex>>,
    mesh_job: Option<Job<vv_core::cartoon::ExpandedMesh>>,
    ses_job: Option<Job<SesBuilt>>,
    /// Out of GPU memory, or refused (skin surface too large): draws
    /// nothing, and nothing retries until the rep changes.
    failed: bool,
}

struct Entry {
    /// The structure's atoms, uploaded when a rep first draws all of them.
    /// Its bonds are `loaded.bonds.pairs` filtered to `Single` (see
    /// `single_pairs`) whenever any bond isn't -- a non-`Single` bond
    /// draws only as its `StrandGeometry` strands, never also as this
    /// ordinary center cylinder, or double/triple bonds would draw one
    /// cylinder too many, overlapping the strands instead of leaving the
    /// gap between them the task asks for.
    gpu: Option<GpuStructure>,
    /// The exact bonds `gpu` was built from when that differs from
    /// `loaded.bonds.pairs` (some bond is non-`Single`): `None` keeps the
    /// all-`Single` case's `DrawSource::bond_atoms` at `None` too (a pick
    /// indexes `loaded.bonds.pairs` directly, byte-identical to before
    /// this feature). `Some` when filtered, so picking a single bond
    /// after a removed one still resolves to the right atom pair.
    single_pairs: Option<Arc<Vec<[u32; 2]>>>,
    /// `loaded.bonds`'s adjacency, for strand plane lookups
    /// (`build_strand_geometry`); built once and shared by every rep,
    /// since it depends only on connectivity, never on a rep's style.
    adjacency: Option<Arc<vv_core::Adjacency>>,
    /// The coordinate set `gpu` holds.
    frame: usize,
    /// `Structure::frames_id` of the coordinates everything here was built
    /// from; a different one (an attached trajectory) discards the entry.
    coords: usize,
    /// The altloc policy `reps` were selected under; a change discards the entry.
    altloc: vv_core::altloc::AltlocPolicy,
    /// The whole structure's cartoon, and its plan being built.
    cartoon: Option<CartoonModel>,
    cartoon_job: Option<Job<(CartoonPlan, CartoonFrame)>>,
    /// The contact dashes, when any overlay is on.
    overlay: Option<Overlay>,
    reps: HashMap<RepId, RepEntry>,
}

/// The structure's bond adjacency, built on first use and shared by every
/// rep; `None` when every bond is `Single` (no strand needs it).
fn adjacency_of(
    cached: &mut Option<Arc<vv_core::Adjacency>>,
    loaded: &LoadedStructure,
) -> Option<Arc<vv_core::Adjacency>> {
    (!loaded.bonds.orders.is_empty()).then(|| {
        cached
            .get_or_insert_with(|| {
                Arc::new(vv_core::adjacency(
                    &loaded.bonds,
                    loaded.structure.atom_count(),
                ))
            })
            .clone()
    })
}

/// Whether a selection depends on coordinates, so it is re-evaluated when
/// the frame changes.
fn is_geometric(selection: &str) -> bool {
    selection.contains("within") || selection.contains("around")
}

/// The atoms `rep` selects in `loaded` at `frame`.
pub(crate) fn select_atoms(
    loaded: &LoadedStructure,
    rep: &Rep,
    frame: usize,
) -> Result<Atoms, String> {
    let shown = loaded.shown_atoms();
    let mut bits = if rep.selects_all() {
        match &shown {
            Some(shown) => shown.clone(),
            None => return Ok(None),
        }
    } else {
        loaded
            .select(&rep.selection, frame)
            .map_err(|e| format!("{}: rep selection `{}`: {e}", loaded.label, rep.selection))?
    };
    if let Some(shown) = &shown {
        bits.intersect_with(shown);
    }
    Ok(Some(Arc::new(bits.ones().map(|a| a as u32).collect())))
}

pub(crate) fn keep_mask(atoms: &Atoms, atom_count: usize) -> Option<Vec<bool>> {
    atoms.as_ref().map(|list| {
        let mut keep = vec![false; atom_count];
        for &a in list.iter() {
            keep[a as usize] = true;
        }
        keep
    })
}

#[derive(Default)]
pub struct GpuCache {
    entries: HashMap<StructureId, Entry>,
    waker: Option<Waker>,
    /// Notices for the user (e.g. a refused skin surface), drained by
    /// `take_messages`.
    messages: Vec<String>,
    /// Progress notes (e.g. a background build started), shown as info
    /// rather than errors; drained by `take_infos`.
    infos: Vec<String>,
    /// What the occlusion volume was built from, and the volume.
    occlusion_key: Vec<OcclusionKey>,
    occlusion: Option<Arc<OcclusionVolume>>,
    /// A volume being baked on a worker, and what it is built from.
    occlusion_job: Option<Job<Option<Arc<OcclusionVolume>>>>,
    occlusion_pending: Vec<OcclusionKey>,
}

/// Everything about one rep that changes what it puts in the occlusion
/// volume.
#[derive(Clone, PartialEq)]
struct OcclusionKey {
    id: StructureId,
    rep: RepId,
    representation: SceneRepresentation,
    selection: (String, usize),
    frame: usize,
    /// Which built geometry the proxies come from; 0 while pending.
    geometry: usize,
}

impl RepEntry {
    fn geometry_id(&self) -> usize {
        match &self.geometry {
            RepGeometry::Pending => 0,
            RepGeometry::Cartoon(c) => Arc::as_ptr(&c.from) as usize,
            RepGeometry::Tube(t) => Arc::as_ptr(&t.gpu.source) as usize,
            _ => 1,
        }
    }
}

/// Spheres standing in for what a rep draws: atoms at their drawn
/// radius, the tube's trace, the cartoon's cross-sections.
fn occlusion_proxies(
    loaded: &LoadedStructure,
    rep: &RepEntry,
    frame: usize,
    out: &mut Vec<(Vec3, f32)>,
) {
    let coords = loaded.structure.frame(frame);
    let positions = coords.positions();
    let elements = &loaded.structure.topology.element;
    let mut atoms = |sizes: AtomSizes, list: &Atoms| match list {
        None => out.par_extend(
            positions
                .par_iter()
                .zip(elements)
                .map(|(p, e)| (*p, sizes.atom_radius(e.vdw_radius()))),
        ),
        Some(list) => out.par_extend(list.par_iter().map(|&a| {
            (
                positions[a as usize],
                sizes.atom_radius(elements[a as usize].vdw_radius()),
            )
        })),
    };
    match &rep.geometry {
        // A wire mesh hides nothing behind it.
        RepGeometry::Pending | RepGeometry::Empty | RepGeometry::Mesh(_) => {}
        RepGeometry::AllAtoms { .. } | RepGeometry::SomeAtoms(_) => {
            atoms(sizes_of(rep.representation, &rep.options), &rep.atoms)
        }
        RepGeometry::Gaussian(_) | RepGeometry::Skin(_) | RepGeometry::Ses(_) => {
            atoms(AtomSizes::of(GpuRepresentation::Spacefill), &rep.atoms)
        }
        RepGeometry::Tube(t) => {
            out.extend(
                t.gpu
                    .source
                    .iter()
                    .map(|&a| (positions[a as usize], vv_core::TUBE_RADIUS)),
            );
            if let Some(caps) = &t.caps {
                out.extend(
                    caps.atom_map
                        .iter()
                        .zip(&t.cap_radii)
                        .map(|(&a, &r)| (positions[a as usize], r)),
                );
            }
            bases::push_stick_proxies(&t.bases, positions, out);
            t.companions.proxies(positions, elements, out);
        }
        RepGeometry::Cartoon(c) => {
            // One sphere per residue on the spline, as wide as a ribbon.
            let keep = loaded.trace_keep(&keep_mask(&rep.atoms, positions.len()));
            let radius = 0.5 * vv_core::cartoon::RIBBON_WIDTH;
            out.extend(
                c.from
                    .slot_atoms()
                    .iter()
                    .zip(&c.spline.controls)
                    .filter(|(&a, _)| keep.as_ref().is_none_or(|k| k[a as usize]))
                    .map(|(_, &p)| (p, radius)),
            );
            bases::push_stick_proxies(&c.bases, positions, out);
            c.companions.proxies(positions, elements, out);
        }
        RepGeometry::Glycan(g) => {
            // One sphere per shape, at its ring centroid.
            out.extend(g.frame_buf.centroid.iter().map(|&p| (p, 0.5 * g.size)));
        }
    }
}

/// The frame `sync`'s per-structure checks are keyed on: `live` while
/// `visible`, else `cached` (the entry's last synced frame) held still, so
/// every frame/spline/colour/upload comparison below reads "unchanged"
/// and does no work. Its own function so the freeze rule is checked
/// without a GPU context (`sync` itself needs one to run at all).
fn sync_frame(visible: bool, cached: usize, live: usize) -> usize {
    if visible {
        live
    } else {
        cached
    }
}

impl GpuCache {
    /// Called by a worker when a background build finishes.
    pub fn set_waker(&mut self, waker: Waker) {
        self.waker = Some(waker);
    }

    /// Rebuilds the occlusion volume (`vv_render::OcclusionVolume`) when
    /// what is drawn opaque changed since the last call: reps, styles,
    /// selections, frames, glass, or geometry finishing a background
    /// build. The bake runs on a worker (tens of milliseconds at millions
    /// of atoms, every frame of a playing trajectory), the old volume
    /// staying in use meanwhile. Returns a finished volume, for
    /// `Renderer::set_occlusion_volume`. Camera moves never rebuild it.
    pub fn update_occlusion(
        &mut self,
        scene: &Scene,
        ctx: &Arc<GpuContext>,
    ) -> Option<Option<Arc<OcclusionVolume>>> {
        let mut finished = None;
        match poll_job(&mut self.occlusion_job) {
            Some(Ok((_, volume))) => {
                self.occlusion = volume;
                self.occlusion_key = std::mem::take(&mut self.occlusion_pending);
                finished = Some(self.occlusion.clone());
            }
            Some(Err(())) => self.occlusion_pending.clear(),
            None => {}
        }
        let mut key = Vec::new();
        for (id, loaded) in scene.structures() {
            if !loaded.visible {
                continue;
            }
            let Some(entry) = self.entries.get(&id) else {
                continue;
            };
            for rep in loaded.reps.iter().filter(|r| r.visible) {
                let glass = vv_render::Material::from(rep.material).is_glass()
                    && matches!(
                        rep.representation,
                        SceneRepresentation::GaussianSurface
                            | SceneRepresentation::SkinSurface
                            | SceneRepresentation::Ses
                    );
                let Some(r) = entry.reps.get(&rep.id).filter(|_| !glass) else {
                    continue;
                };
                key.push(OcclusionKey {
                    id,
                    rep: rep.id,
                    representation: r.representation,
                    selection: r.selection.clone(),
                    frame: r.frame,
                    geometry: r.geometry_id(),
                });
            }
        }
        // Up to date, or a bake is running (the next call starts one for
        // whatever is drawn by then).
        if key == self.occlusion_key || self.occlusion_job.is_some() {
            return finished;
        }
        let mut spheres = Vec::new();
        for k in &key {
            let (Some(loaded), Some(entry)) = (scene.structure(k.id), self.entries.get(&k.id))
            else {
                continue;
            };
            if let Some(rep) = entry.reps.get(&k.rep) {
                occlusion_proxies(loaded, rep, rep.frame, &mut spheres);
            }
        }
        let ctx = ctx.clone();
        self.occlusion_pending = key;
        self.occlusion_job = Some(spawn_job(0, self.waker.clone(), move || {
            OcclusionVolume::build(&ctx, &spheres).map(Arc::new)
        }));
        finished
    }

    /// The current occlusion volume (an export's own renderer needs it).
    pub fn occlusion(&self) -> Option<Arc<OcclusionVolume>> {
        self.occlusion.clone()
    }

    /// Notices produced since the last call.
    pub fn take_infos(&mut self) -> Vec<String> {
        std::mem::take(&mut self.infos)
    }

    pub fn take_messages(&mut self) -> Vec<String> {
        std::mem::take(&mut self.messages)
    }

    /// Whether any cartoon or skin surface is still building.
    /// Structure `id`'s cartoon as the viewport draws it: its plan and
    /// the spline at the current frame, once built.
    pub fn cartoon(&self, id: StructureId) -> Option<(Arc<CartoonPlan>, Arc<CartoonFrame>)> {
        let model = self.entries.get(&id)?.cartoon.as_ref()?;
        Some((model.plan.clone(), model.spline.clone()))
    }

    pub fn building(&self) -> bool {
        self.occlusion_job.is_some()
            || self.entries.values().any(|e| {
                e.cartoon_job.is_some()
                    || e.reps.values().any(|r| {
                        r.skin_job.is_some() || r.ses_job.is_some() || r.mesh_job.is_some()
                    })
            })
    }

    /// Whether `rep`'s own geometry is still building in the background
    /// (a Structures panel row's loading spinner): skin/SES jobs are
    /// tracked per rep; a cartoon plan is tracked per structure (one
    /// cartoon shape can back several cartoon reps), so `is_cartoon`
    /// tells this whether that structure-wide job belongs to `rep`.
    pub fn rep_building(&self, id: StructureId, rep: RepId, is_cartoon: bool) -> bool {
        let Some(entry) = self.entries.get(&id) else {
            return false;
        };
        (is_cartoon && entry.cartoon_job.is_some())
            || entry.reps.get(&rep).is_some_and(|r| {
                r.skin_job.is_some() || r.ses_job.is_some() || r.mesh_job.is_some()
            })
    }

    /// Brings the GPU copies in step with `scene`: drops what closed,
    /// uploads what is new, and rebuilds, recolours or moves what changed.
    /// Call once per frame before rendering; cheap when nothing changed.
    /// `playing`: a trajectory is playing, so cartoons keep their
    /// secondary structure and only follow the atoms.
    ///
    /// A hidden structure (`LoadedStructure::visible`) costs nothing here:
    /// `sync_frame` pins the frame every check below is keyed on to
    /// whatever it already was, so no rebuild, recolour or upload is
    /// judged necessary while hidden -- `sync_cartoon`/`sync_rep` still
    /// run (so a job already in flight is still polled and applied, never
    /// left to finish into a void), they just find nothing changed. Un-
    /// hiding lets the real frame back in on the next call, so every
    /// check fires once and the structure catches up in a single sync,
    /// not one step per frame missed while hidden.
    pub fn sync(
        &mut self,
        scene: &Scene,
        ctx: &Arc<GpuContext>,
        renderer: &Renderer,
        playing: bool,
    ) {
        // Out of `self` while the helpers below post notices to it.
        let mut entries = std::mem::take(&mut self.entries);
        entries.retain(|id, _| scene.structure(*id).is_some());
        let waker = self.waker.clone();
        for (id, loaded) in scene.structures() {
            let live_frame = loaded
                .frame
                .min(loaded.structure.frame_count().saturating_sub(1));
            let coords = loaded.structure.frames_id();
            if entries
                .get(&id)
                .is_some_and(|e| e.coords != coords || e.altloc != loaded.altloc)
            {
                entries.remove(&id);
            }
            let entry = entries.entry(id).or_insert_with(|| Entry {
                gpu: None,
                single_pairs: None,
                adjacency: None,
                frame: live_frame,
                coords,
                altloc: loaded.altloc,
                cartoon: None,
                cartoon_job: None,
                overlay: None,
                reps: HashMap::new(),
            });
            let frame = sync_frame(loaded.visible, entry.frame, live_frame);
            entry.reps.retain(|r, _| loaded.rep_index(*r).is_some());
            self.sync_cartoon(entry, loaded, frame, playing, &waker);
            for rep in loaded.reps.iter().filter(|r| r.visible) {
                self.sync_rep(entry, loaded, rep, frame, playing, ctx, renderer, &waker);
            }
            Overlay::sync(&mut entry.overlay, ctx, renderer, loaded, frame);
            // The shared atoms move only while a visible rep draws them (a
            // playing cartoon would otherwise re-upload every atom too).
            let drawn = loaded.reps.iter().filter(|r| r.visible).any(|rep| {
                entry
                    .reps
                    .get(&rep.id)
                    .is_some_and(|r| matches!(r.geometry, RepGeometry::AllAtoms { .. }))
            });
            if entry.frame != frame && drawn {
                if let Some(gpu) = &entry.gpu {
                    gpu.set_frame(ctx, &loaded.structure, frame);
                }
                let coords = loaded.structure.frame(frame);
                let positions = coords.positions();
                for r in entry.reps.values() {
                    if let RepGeometry::AllAtoms {
                        strands: Some(s), ..
                    } = &r.geometry
                    {
                        s.set_frame(ctx, positions);
                    }
                }
                entry.frame = frame;
            }
        }
        self.entries = entries;
    }

    /// The whole structure's cartoon while any visible rep draws one:
    /// planned on a worker, its spline moved with the frame at once, and
    /// re-planned (secondary structure) for a new frame once playback
    /// stops.
    fn sync_cartoon(
        &mut self,
        entry: &mut Entry,
        loaded: &LoadedStructure,
        frame: usize,
        playing: bool,
        waker: &Option<Waker>,
    ) {
        let wanted = loaded
            .reps
            .iter()
            .any(|r| r.visible && r.representation == SceneRepresentation::Cartoon);
        if !wanted {
            // Dropped rather than kept stale: switching back must never
            // show a cartoon from an older frame.
            entry.cartoon = None;
            entry.cartoon_job = None;
            return;
        }
        match poll_job(&mut entry.cartoon_job) {
            Some(Ok((built, (plan, spline)))) => {
                entry.cartoon = Some(CartoonModel {
                    plan: Arc::new(plan),
                    planned_at: built,
                    spline: Arc::new(spline),
                    spline_at: built,
                })
            }
            Some(Err(())) => self.messages.push("cartoon build failed".into()),
            None => {}
        }
        if let Some(model) = &mut entry.cartoon {
            if model.spline_at != frame {
                model.spline =
                    Arc::new(model.plan.frame(loaded.drawn_positions(frame).positions()));
                model.spline_at = frame;
            }
        }
        // One job at a time, and none while playing: a new plan waits for
        // the frame to stay put.
        let planned_at = entry.cartoon.as_ref().map(|c| c.planned_at);
        let replan = planned_at.is_none() || (planned_at != Some(frame) && !playing);
        if replan && entry.cartoon_job.is_none() {
            if entry.cartoon.is_none() {
                self.infos.push(format!(
                    "building cartoon for {}; it appears when ready",
                    loaded.label
                ));
            }
            let positions = loaded.drawn_positions(frame).positions().to_vec();
            entry.cartoon_job = Some(spawn_cartoon(
                &loaded.structure,
                positions,
                frame,
                waker.clone(),
            ));
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn sync_rep(
        &mut self,
        entry: &mut Entry,
        loaded: &LoadedStructure,
        rep: &Rep,
        frame: usize,
        playing: bool,
        ctx: &Arc<GpuContext>,
        renderer: &Renderer,
        waker: &Option<Waker>,
    ) {
        let selection = (
            rep.selection.clone(),
            if is_geometric(&rep.selection) {
                frame
            } else {
                0
            },
        );
        let colors_key = ColorKey::of(loaded, &rep.coloring, frame);
        let stale = entry.reps.get(&rep.id).is_none_or(|r| {
            r.representation != rep.representation
                || r.selection != selection
                || r.options != rep.options
        });
        if stale {
            let atoms = match select_atoms(loaded, rep, frame) {
                Ok(atoms) => atoms,
                Err(message) => {
                    self.messages.push(message);
                    Some(Arc::new(Vec::new()))
                }
            };
            let previous = entry.reps.remove(&rep.id);
            // Skin and SES build in the background (seconds at 100K+
            // atoms): a selection or option change on the same
            // representation keeps the old surface on screen -- stale but
            // still valid -- until the match arms below finish a rebuild,
            // so it never flashes to nothing while it rebuilds.
            // `skin_job`/`ses_job` are dropped either way: a job in flight
            // was building for the old selection or options. Every other
            // representation rebuilds synchronously right here (cheap) or
            // keys its own freshness off geometry already being `Pending`
            // (Tube, Cartoon), so it always resets on `stale`.
            let carries_over = matches!(
                rep.representation,
                SceneRepresentation::SkinSurface | SceneRepresentation::Ses
            ) || surface_mesh::style(rep).is_some();
            let geometry = match previous {
                Some(p) if carries_over && p.representation == rep.representation => p.geometry,
                _ => RepGeometry::Pending,
            };
            entry.reps.insert(
                rep.id,
                RepEntry {
                    representation: rep.representation,
                    selection,
                    options: rep.options.clone(),
                    atoms,
                    colors: colors_key.clone(),
                    frame,
                    geometry,
                    skin_job: None,
                    ses_job: None,
                    mesh_job: None,
                    failed: false,
                },
            );
        }
        let r = entry.reps.get_mut(&rep.id).expect("inserted above");
        if r.failed {
            return;
        }
        let recolor = !r.colors.same(&colors_key);
        let moved = r.frame != frame;
        let mut all_colors: Option<Vec<u32>> = None;
        let mut colors = || -> Vec<u32> {
            all_colors
                .get_or_insert_with(|| colors_of(loaded, &rep.coloring, frame))
                .clone()
        };
        let structure = &loaded.structure;
        match rep.representation {
            SceneRepresentation::GaussianSurface
            | SceneRepresentation::SkinSurface
            | SceneRepresentation::Ses
                if surface_mesh::style(rep).is_some() =>
            {
                let style = surface_mesh::style(rep).expect("guarded");
                self.sync_surface_mesh(
                    r,
                    loaded,
                    rep,
                    style,
                    (frame, recolor),
                    ctx,
                    renderer,
                    waker,
                    &mut colors,
                );
            }
            SceneRepresentation::Spacefill
            | SceneRepresentation::BallAndStick
            | SceneRepresentation::Sas
            | SceneRepresentation::Sticks
            | SceneRepresentation::Lines => match &r.geometry {
                RepGeometry::Pending => {
                    r.geometry = match &r.atoms {
                        None => {
                            // A non-`Single` bond draws only as its
                            // `StrandGeometry` strands (below), never
                            // also as this ordinary center cylinder --
                            // else a double bond would draw 3 overlapping
                            // cylinders instead of 2 separated ones.
                            // Computed once (`entry.gpu.is_none()`, the
                            // same "build once" gate `get_or_insert_with`
                            // uses below) and cached on `Entry`, not
                            // recomputed for every rep drawing all atoms.
                            if entry.gpu.is_none() && !loaded.bonds.orders.is_empty() {
                                entry.single_pairs =
                                    Some(Arc::new(single_bond_pairs(&loaded.bonds)));
                            }
                            let bond_pairs: &[[u32; 2]] = entry
                                .single_pairs
                                .as_deref()
                                .map_or(&loaded.bonds.pairs[..], |v| v.as_slice());
                            let gpu = entry.gpu.get_or_insert_with(|| {
                                let coords = structure.frame(frame);
                                let positions = coords.positions();
                                GpuStructure::geometry(
                                    ctx,
                                    positions,
                                    &vdw_radii(&structure.topology.element),
                                    bond_pairs,
                                )
                            });
                            let state = DrawState::new(ctx, gpu, &colors());
                            let bindings = renderer.bind_state(gpu, &state);
                            let bond_radius = rep_sizes(rep).bond_radius;
                            // The `orders.is_empty()` check keeps the
                            // common all-`Single` case from ever building
                            // an `Adjacency` (O(bonds)) it wouldn't use.
                            let strands = if bond_radius > 0.0 && !loaded.bonds.orders.is_empty() {
                                let coords = structure.frame(frame);
                                let adjacency = entry.adjacency.get_or_insert_with(|| {
                                    Arc::new(vv_core::adjacency(
                                        &loaded.bonds,
                                        structure.atom_count(),
                                    ))
                                });
                                build_strand_geometry(
                                    ctx,
                                    renderer,
                                    &loaded.bonds,
                                    adjacency,
                                    coords.positions(),
                                    bond_radius * STRAND_SEPARATION_SCALE,
                                    &colors(),
                                )
                            } else {
                                None
                            };
                            RepGeometry::AllAtoms {
                                state,
                                bindings,
                                strands,
                            }
                        }
                        Some(list) => {
                            let keep = keep_mask(&Some(list.clone()), structure.atom_count())
                                .expect("a list");
                            let bond_radius = rep_sizes(rep).bond_radius;
                            let strand_input =
                                (bond_radius > 0.0 && !loaded.bonds.orders.is_empty()).then(|| {
                                    let adjacency = entry.adjacency.get_or_insert_with(|| {
                                        Arc::new(vv_core::adjacency(
                                            &loaded.bonds,
                                            structure.atom_count(),
                                        ))
                                    });
                                    StrandInput {
                                        bonds: &loaded.bonds,
                                        adjacency,
                                        separation: bond_radius * STRAND_SEPARATION_SCALE,
                                    }
                                });
                            RepGeometry::SomeAtoms(Derived::atoms(
                                ctx,
                                renderer,
                                structure,
                                frame,
                                subset_of(&keep, &loaded.bonds.pairs),
                                &colors(),
                                strand_input,
                            ))
                        }
                    };
                }
                RepGeometry::AllAtoms { state, strands, .. } if recolor => {
                    state.set_colors(ctx, &colors());
                    if let Some(s) = strands {
                        s.recolor(ctx, &colors());
                    }
                }
                RepGeometry::SomeAtoms(d) => {
                    if recolor {
                        d.recolor(ctx, &colors());
                    }
                    if moved {
                        d.set_frame(ctx, structure, frame);
                    }
                }
                _ => {}
            },
            // A recolor rebuilds the mesh from scratch, like Cartoon:
            // `CartoonGpu` has no incremental recolor (its sections'
            // colors are baked in at upload).
            SceneRepresentation::Tube => match &mut r.geometry {
                RepGeometry::Tube(t) if !recolor => {
                    if moved {
                        let keep = keep_mask(&r.atoms, structure.atom_count());
                        let style = bases::BaseStyle::of(rep);
                        t.set_frame(ctx, renderer, loaded, (&keep, style, &colors()), frame);
                    }
                }
                _ => {
                    let keep = keep_mask(&r.atoms, structure.atom_count());
                    let adjacency = adjacency_of(&mut entry.adjacency, loaded);
                    match build_tube(
                        ctx,
                        renderer,
                        loaded,
                        rep,
                        keep,
                        &colors(),
                        frame,
                        adjacency.as_deref(),
                    ) {
                        Ok(geometry) => r.geometry = RepGeometry::Tube(Box::new(geometry)),
                        Err(OutOfGpuMemory) => {
                            r.geometry = RepGeometry::Pending;
                            r.failed = true;
                            self.messages.push(oom_notice("tube", &loaded.label));
                        }
                    }
                }
            },
            SceneRepresentation::Cartoon => {
                let Some(model) = &entry.cartoon else {
                    return;
                };
                if let RepGeometry::Cartoon(c) = &mut r.geometry {
                    if Arc::ptr_eq(&c.from, &model.plan) && !recolor {
                        // Same plan: move the spline, bases and companions.
                        if !Arc::ptr_eq(&c.spline, &model.spline) {
                            match &c.cylinders {
                                None => c.gpu.set_frame(ctx, &model.spline),
                                Some(styled) => {
                                    let positions = loaded.drawn_positions(frame);
                                    c.gpu.set_frame(ctx, &styled.frame(positions.positions()));
                                }
                            }
                            c.spline = model.spline.clone();
                            c.companions.set_frame(ctx, structure, frame);
                            let keep = keep_mask(&r.atoms, structure.atom_count());
                            let style = bases::BaseStyle::of(rep);
                            match Bases::build(
                                ctx,
                                renderer,
                                loaded,
                                &keep,
                                frame,
                                &colors(),
                                style,
                            ) {
                                Ok(bases) => c.bases = bases,
                                Err(OutOfGpuMemory) => c.bases = Bases::default(),
                            }
                        }
                        r.colors = colors_key;
                        r.frame = frame;
                        return;
                    }
                }
                let keep = keep_mask(&r.atoms, structure.atom_count());
                let cylinders = cylinder_cartoon(
                    rep,
                    &model.plan,
                    &model.spline,
                    loaded.drawn_positions(frame).positions(),
                    vv_core::cartoon::SAMPLES_PER_RESIDUE,
                );
                let (plan, spline) = match &cylinders {
                    Some((plan, spline)) => (plan, spline),
                    None => (&*model.plan, &*model.spline),
                };
                let filtered = match &loaded.trace_keep(&keep) {
                    None => plan.clone(),
                    Some(keep) => plan.filter(|a| keep[a as usize]),
                };
                let all = colors();
                let adjacency = adjacency_of(&mut entry.adjacency, loaded);
                let built = CartoonGpu::upload(ctx, &filtered, spline, &all).and_then(|gpu| {
                    let style = bases::BaseStyle::of(rep);
                    let bases = Bases::build(ctx, renderer, loaded, &keep, frame, &all, style)?;
                    Ok((gpu, bases))
                });
                match built {
                    Ok((gpu, bases)) => {
                        let bindings = renderer.bind_cartoon(&gpu);
                        let traced = companions::traced_atoms(&filtered);
                        let drawn = companions::Drawn {
                            keep: &keep,
                            traced: &traced,
                        };
                        let companions = Companions::build(
                            ctx,
                            renderer,
                            loaded,
                            rep,
                            drawn,
                            frame,
                            &all,
                            adjacency.as_deref(),
                        );
                        r.geometry = RepGeometry::Cartoon(Box::new(CartoonEntry {
                            gpu,
                            bindings,
                            from: model.plan.clone(),
                            spline: model.spline.clone(),
                            cylinders: cylinders.map(|(plan, _)| plan),
                            bases,
                            companions,
                        }));
                    }
                    Err(OutOfGpuMemory) => {
                        r.geometry = RepGeometry::Pending;
                        r.failed = true;
                        self.messages.push(oom_notice("cartoon", &loaded.label));
                    }
                }
            }
            SceneRepresentation::GaussianSurface => {
                // A playing trajectory re-bakes a coarser volume every
                // frame; the full one comes back when it stops.
                let detail_changed =
                    matches!(&r.geometry, RepGeometry::Gaussian(g) if g.playback != playing);
                if matches!(r.geometry, RepGeometry::Pending) || recolor || moved || detail_changed
                {
                    // A new frame of the same atoms: only their positions,
                    // baked into the same volume.
                    if let RepGeometry::Gaussian(g) = &mut r.geometry {
                        if moved
                            && !recolor
                            && !detail_changed
                            && g.gpu
                                .set_frame(ctx, &atom_positions(structure, &r.atoms, frame))
                        {
                            r.colors = colors_key;
                            r.frame = frame;
                            return;
                        }
                    }
                    let (positions, radii, surface_colors) =
                        atom_arrays(structure, &r.atoms, frame, &colors());
                    if positions.is_empty() {
                        r.geometry = RepGeometry::Empty;
                        r.colors = colors_key;
                        r.frame = frame;
                        return;
                    }
                    let voxels = if playing {
                        vv_render::scene::GAUSSIAN_PLAYBACK_VOXELS
                    } else {
                        vv_render::scene::GAUSSIAN_MAX_VOXELS
                    };
                    match GaussianSurfaceGpu::upload_within(
                        ctx,
                        &positions,
                        &radii,
                        &surface_colors,
                        rep.option("blob")
                            .unwrap_or(vv_core::gaussian_surface::DEFAULT_BLOB_FACTOR),
                        voxels,
                    ) {
                        Ok(gpu) => {
                            let bindings = renderer.bind_gaussian_surface(&gpu);
                            r.geometry = RepGeometry::Gaussian(GaussianSurfaceEntry {
                                gpu,
                                bindings,
                                playback: playing,
                            });
                        }
                        Err(OutOfGpuMemory) => {
                            r.geometry = RepGeometry::Pending;
                            r.failed = true;
                            self.messages
                                .push(oom_notice("Gaussian surface", &loaded.label));
                        }
                    }
                }
            }
            SceneRepresentation::SkinSurface => {
                let count = r.atoms.as_ref().map_or(structure.atom_count(), |a| a.len());
                if count == 0 {
                    r.geometry = RepGeometry::Empty;
                    r.skin_job = None;
                    return;
                }
                if count > SKIN_MAX_ATOMS {
                    // Whatever was showing no longer matches `r.atoms`: drop
                    // it rather than leave a surface of the wrong selection.
                    r.geometry = RepGeometry::Pending;
                    r.failed = true;
                    self.messages.push(format!(
                        "skin surface is limited to {SKIN_MAX_ATOMS} atoms for now ({} has {count} \
                         in this rep); try `rep gaussiansurface` or a smaller selection.",
                        loaded.label
                    ));
                    return;
                }
                match poll_job(&mut r.skin_job) {
                    Some(Ok((built, complex))) => {
                        let (positions, radii, surface_colors) =
                            atom_arrays(structure, &r.atoms, built, &colors());
                        // Runs on the main thread (it needs `ctx`'s device):
                        // a stall here, unlike the build above, is a frame
                        // hitch, not a background wait.
                        let gpu = timed("skin surface gpu upload", positions.len(), || {
                            SkinSurfaceGpu::from_complex(
                                ctx,
                                &complex,
                                &positions,
                                &skin_weights(&radii, complex.shrink),
                                &surface_colors,
                            )
                        });
                        let bindings = renderer.bind_skin_surface(&gpu);
                        r.geometry = RepGeometry::Skin(SkinSurfaceEntry {
                            gpu,
                            bindings,
                            atoms: r.atoms.clone(),
                            frame: built,
                        });
                    }
                    Some(Err(())) => self.messages.push("skin surface build failed".into()),
                    None => {}
                }
                // `stale` (the selection or an option just changed) makes
                // whatever is showing outdated even when its frame still
                // matches; keep drawing it (no flash to nothing) until the
                // rebuild below replaces it.
                let outdated =
                    stale || !matches!(&r.geometry, RepGeometry::Skin(s) if s.frame == frame);
                if outdated {
                    if r.skin_job.is_none() {
                        if matches!(r.geometry, RepGeometry::Pending) {
                            self.infos.push(format!(
                                "building skin surface for {} ({count} atoms); it appears when ready",
                                loaded.label
                            ));
                        }
                        let (positions, radii, _) =
                            atom_arrays(structure, &r.atoms, frame, &colors());
                        let shrink = rep
                            .option("shrink")
                            .unwrap_or(vv_core::skin_surface::DEFAULT_SHRINK);
                        r.skin_job =
                            Some(spawn_skin(positions, &radii, shrink, frame, waker.clone()));
                    }
                    // Else a rebuild is already in flight: wait for it.
                } else if recolor {
                    // Colours are per-atom, not baked into the patches
                    // (unlike SES): a rewrite of this one buffer is enough,
                    // no need to rebuild the complex's GPU upload. `s.atoms`
                    // (not `r.atoms`, which a newer selection may already
                    // have replaced) is what this surface was built from.
                    if let RepGeometry::Skin(s) = &r.geometry {
                        let (_, _, surface_colors) =
                            atom_arrays(structure, &s.atoms, s.frame, &colors());
                        timed("skin surface recolor", s.gpu.atom_count as usize, || {
                            s.gpu.set_colors(ctx, &surface_colors)
                        });
                    }
                }
            }
            SceneRepresentation::Ses => {
                let atoms = without_water(structure, &r.atoms);
                match poll_job(&mut r.ses_job) {
                    Some(Ok((built, (upload, used)))) => {
                        match upload {
                            Ok(gpu) => {
                                // Recoloured while it was building.
                                let (_, _, now) = atom_arrays(structure, &atoms, built, &colors());
                                if now != used {
                                    gpu.set_colors(ctx, &now);
                                }
                                let bindings = renderer.bind_ses_surface(&gpu);
                                r.geometry = RepGeometry::Ses(SesEntry {
                                    gpu,
                                    bindings,
                                    atoms: atoms.clone(),
                                    frame: built,
                                });
                            }
                            Err(OutOfGpuMemory) => {
                                // Whatever was showing no longer matches
                                // `atoms`: drop it rather than leave a
                                // surface of the wrong selection.
                                r.geometry = RepGeometry::Pending;
                                r.failed = true;
                                self.messages.push(format!(
                                    "{}: not enough GPU memory for the molecular surface; try `rep gaussiansurface` or a smaller selection",
                                    loaded.label
                                ));
                                return;
                            }
                        }
                    }
                    Some(Err(())) => self.messages.push("SES build failed".into()),
                    None => {}
                }
                // `stale` (the selection or an option just changed) makes
                // whatever is showing outdated even when its frame still
                // matches; keep drawing it (no flash to nothing) until the
                // rebuild below replaces it.
                let outdated =
                    stale || !matches!(&r.geometry, RepGeometry::Ses(s) if s.frame == frame);
                if outdated {
                    if r.ses_job.is_none() {
                        if matches!(r.geometry, RepGeometry::Pending) {
                            let count = atoms.as_ref().map_or(structure.atom_count(), |a| a.len());
                            self.infos.push(format!(
                                "building the molecular surface for {} ({count} atoms); it appears \
                                 when ready",
                                loaded.label
                            ));
                        }
                        let (positions, radii, surface_colors) =
                            atom_arrays(structure, &atoms, frame, &colors());
                        r.ses_job = Some(spawn_ses(
                            ctx.clone(),
                            positions,
                            radii,
                            surface_colors,
                            rep.option("probe").unwrap_or(vv_core::ses::WATER_PROBE),
                            frame,
                            waker.clone(),
                        ));
                    }
                    // Else a rebuild is already in flight; wait for it.
                } else if recolor {
                    if let RepGeometry::Ses(s) = &r.geometry {
                        let (_, _, surface_colors) =
                            atom_arrays(structure, &s.atoms, s.frame, &colors());
                        s.gpu.set_colors(ctx, &surface_colors);
                    }
                }
            }
            SceneRepresentation::Glycan => {
                let size = rep.option("size").unwrap_or(4.0);
                let link_radius = rep.option("radius").unwrap_or(0.5);
                if !moved {
                    if let RepGeometry::Glycan(g) = &r.geometry {
                        if g.size == size && g.radius == link_radius {
                            r.colors = colors_key;
                            r.frame = frame;
                            return;
                        }
                    }
                }
                let plan = match &r.geometry {
                    RepGeometry::Glycan(g) => g.plan.clone(),
                    _ => Arc::new(vv_core::GlycanPlan::build(
                        &structure.topology,
                        &loaded.bonds,
                    )),
                };
                let mut frame_buf = match &mut r.geometry {
                    RepGeometry::Glycan(g) if Arc::ptr_eq(&g.plan, &plan) => {
                        std::mem::take(&mut g.frame_buf)
                    }
                    _ => vv_core::GlycanFrame::default(),
                };
                plan.update_into(loaded.drawn_positions(frame).positions(), &mut frame_buf);

                let keep_atoms = loaded.trace_keep(&keep_mask(&r.atoms, structure.atom_count()));
                let keep_residue = |residue: u32| match &keep_atoms {
                    None => true,
                    Some(mask) => plan
                        .residues
                        .iter()
                        .find(|res| res.residue == residue)
                        .is_some_and(|res| mask[res.anomeric() as usize]),
                };

                let mut mesh = vv_core::PolytopeMesh::default();
                vv_core::build_glycan_mesh(&plan, &frame_buf, size, keep_residue, &mut mesh);
                vv_core::build_glycan_linkage_mesh(
                    &plan,
                    &frame_buf,
                    link_radius,
                    keep_residue,
                    &mut mesh,
                );
                match vv_render::GlycanGpu::upload(ctx, &mesh) {
                    Ok(gpu) => {
                        let bindings = renderer.bind_glycan(&gpu);
                        r.geometry = RepGeometry::Glycan(GlycanEntry {
                            gpu,
                            bindings,
                            plan,
                            frame_buf,
                            size,
                            radius: link_radius,
                        });
                    }
                    Err(OutOfGpuMemory) => {
                        r.failed = true;
                        self.messages.push(format!(
                            "{}: not enough GPU memory for the glycan rep",
                            loaded.label
                        ));
                    }
                }
            }
        }
        r.colors = colors_key;
        r.frame = frame;
    }

    /// The draw list for this frame: each structure's visible reps in
    /// order, split by kind into atom `items`, `cartoons`,
    /// `gaussian_surfaces` and `skin_surfaces`. `sources` parallels
    /// `items`, *then* `cartoons`, *then* `skin_surfaces` concatenated in
    /// that order -- exactly how `Renderer::render_all` assigns pick ids --
    /// so a `Pick::Atom { item, atom }` indexes straight into it. Gaussian
    /// surfaces are not pickable, so they have no `sources` entry.
    pub fn draw_items<'a>(
        &'a self,
        scene: &Scene,
    ) -> (
        Vec<DrawSource>,
        Vec<DrawItem<'a>>,
        Vec<CartoonItem<'a>>,
        Vec<GaussianSurfaceItem<'a>>,
        Vec<PatchSurfaceItem<'a>>,
    ) {
        let mut sources = Vec::new();
        let mut items = Vec::new();
        let mut cartoon_sources = Vec::new();
        let mut cartoons = Vec::new();
        let mut gaussian_surfaces = Vec::new();
        let mut skin_sources = Vec::new();
        let mut skin_surfaces = Vec::new();
        for (id, loaded) in scene.structures() {
            if !loaded.visible {
                continue;
            }
            let Some(entry) = self.entries.get(&id) else {
                continue;
            };
            for rep in loaded.reps.iter().filter(|r| r.visible) {
                let Some(r) = entry.reps.get(&rep.id) else {
                    continue;
                };
                let material = vv_render::Material::from(rep.material);
                let mut derived = |d: &'a Derived, representation, sizes: AtomSizes| {
                    sources.push(DrawSource {
                        id,
                        rep: rep.id,
                        atom_map: Some(d.atom_map.clone()),
                        bond_atoms: Some(d.bond_atoms.clone()),
                    });
                    items.push(DrawItem {
                        structure: &d.gpu,
                        bindings: &d.bindings,
                        representation,
                        sizes,
                        material,
                    });
                    if let Some(s) = &d.strands {
                        sources.push(DrawSource {
                            id,
                            rep: rep.id,
                            atom_map: None,
                            bond_atoms: Some(s.bond_atoms.clone()),
                        });
                        items.push(DrawItem {
                            structure: &s.gpu,
                            bindings: &s.bindings,
                            representation,
                            sizes: AtomSizes {
                                radius_scale: 0.0,
                                radius_offset: 0.0,
                                bond_radius: sizes.bond_radius * STRAND_RADIUS_SCALE,
                            },
                            material,
                        });
                    }
                };
                let own = rep_sizes(rep);
                match &r.geometry {
                    RepGeometry::Pending => {}
                    RepGeometry::AllAtoms {
                        bindings, strands, ..
                    } => {
                        let Some(gpu) = &entry.gpu else { continue };
                        sources.push(DrawSource {
                            id,
                            rep: rep.id,
                            atom_map: None,
                            // `Some` exactly when `gpu` was built from a
                            // filtered bond list (see `Entry::
                            // single_pairs`), so a pick's `bond` index
                            // still resolves to the right atom pair.
                            bond_atoms: entry.single_pairs.clone(),
                        });
                        items.push(DrawItem {
                            structure: gpu,
                            bindings,
                            representation: map_representation(rep.representation),
                            sizes: own,
                            material,
                        });
                        if let Some(s) = strands {
                            sources.push(DrawSource {
                                id,
                                rep: rep.id,
                                atom_map: None,
                                bond_atoms: Some(s.bond_atoms.clone()),
                            });
                            items.push(DrawItem {
                                structure: &s.gpu,
                                bindings: &s.bindings,
                                representation: map_representation(rep.representation),
                                sizes: AtomSizes {
                                    radius_scale: 0.0,
                                    radius_offset: 0.0,
                                    bond_radius: own.bond_radius * STRAND_RADIUS_SCALE,
                                },
                                material,
                            });
                        }
                    }
                    RepGeometry::SomeAtoms(d) => {
                        derived(d, map_representation(rep.representation), own)
                    }
                    RepGeometry::Tube(t) => {
                        // An empty mesh's buffers are placeholders its
                        // shaders reject; a ligand-only tube draws just its
                        // sticks.
                        if t.gpu.section_count > 0 {
                            cartoon_sources.push(DrawSource {
                                id,
                                rep: rep.id,
                                atom_map: Some(t.gpu.source.clone()),
                                bond_atoms: None,
                            });
                            cartoons.push(CartoonItem {
                                mesh: vv_render::CartoonMesh::Ribbon(&t.gpu),
                                bindings: &t.bindings,
                                material,
                            });
                        }
                        if let Some(caps) = &t.caps {
                            derived(
                                caps,
                                GpuRepresentation::Spacefill,
                                AtomSizes::of(GpuRepresentation::Spacefill),
                            );
                        }
                        if let Some((source, item)) = t.bases.plates_draw(id, rep.id, material) {
                            cartoon_sources.push(source);
                            cartoons.push(item);
                        }
                        if let Some(sticks) = t.bases.sticks() {
                            derived(sticks, GpuRepresentation::Sticks, bases::stick_sizes());
                        }
                        for (d, representation, sizes) in t.companions.draw() {
                            derived(d, representation, sizes);
                        }
                    }
                    RepGeometry::Cartoon(c) => {
                        if c.gpu.section_count > 0 {
                            cartoon_sources.push(DrawSource {
                                id,
                                rep: rep.id,
                                atom_map: Some(c.gpu.source.clone()),
                                bond_atoms: None,
                            });
                            cartoons.push(CartoonItem {
                                mesh: vv_render::CartoonMesh::Ribbon(&c.gpu),
                                bindings: &c.bindings,
                                material,
                            });
                        }
                        if let Some((source, item)) = c.bases.plates_draw(id, rep.id, material) {
                            cartoon_sources.push(source);
                            cartoons.push(item);
                        }
                        if let Some(sticks) = c.bases.sticks() {
                            derived(sticks, GpuRepresentation::Sticks, bases::stick_sizes());
                        }
                        for (d, representation, sizes) in c.companions.draw() {
                            derived(d, representation, sizes);
                        }
                    }
                    RepGeometry::Glycan(g) if g.gpu.vertex_count > 0 => {
                        cartoon_sources.push(DrawSource {
                            id,
                            rep: rep.id,
                            atom_map: Some(g.gpu.source.clone()),
                            bond_atoms: None,
                        });
                        cartoons.push(CartoonItem {
                            mesh: vv_render::CartoonMesh::Glycan(&g.gpu),
                            bindings: &g.bindings,
                            material,
                        });
                    }
                    RepGeometry::Mesh(m) => match &m.drawn {
                        surface_mesh::Drawn::Edges { gpu, bindings } => {
                            cartoon_sources.push(DrawSource {
                                id,
                                rep: rep.id,
                                atom_map: Some(gpu.source.clone()),
                                bond_atoms: None,
                            });
                            cartoons.push(CartoonItem {
                                mesh: vv_render::CartoonMesh::Glycan(gpu),
                                bindings,
                                material,
                            });
                        }
                        surface_mesh::Drawn::Dots { points, radius } => derived(
                            points,
                            GpuRepresentation::Spacefill,
                            surface_mesh::SurfaceMesh::dot_sizes(*radius),
                        ),
                    },
                    RepGeometry::Glycan(_) | RepGeometry::Empty => {}
                    RepGeometry::Gaussian(g) => gaussian_surfaces.push(GaussianSurfaceItem {
                        gpu: &g.gpu,
                        bindings: &g.bindings,
                        material,
                    }),
                    RepGeometry::Skin(s) => {
                        // A skin pick names the patch's atom within the
                        // atom list this surface was built from (`s.atoms`,
                        // not `r.atoms`: a newer selection may already have
                        // replaced it while this one still shows).
                        skin_sources.push(DrawSource {
                            id,
                            rep: rep.id,
                            atom_map: s.atoms.clone(),
                            bond_atoms: None,
                        });
                        skin_surfaces.push(PatchSurfaceItem {
                            surface: PatchSurface::Skin(&s.gpu, &s.bindings),
                            material,
                        });
                    }
                    RepGeometry::Ses(s) => {
                        skin_sources.push(DrawSource {
                            id,
                            rep: rep.id,
                            atom_map: s.atoms.clone(),
                            bond_atoms: None,
                        });
                        skin_surfaces.push(PatchSurfaceItem {
                            surface: PatchSurface::Ses(&s.gpu, &s.bindings),
                            material,
                        });
                    }
                }
            }
            if let Some(overlay) = &entry.overlay {
                let (source, item) = overlay.draw(id);
                sources.push(source);
                items.push(item);
            }
        }
        sources.extend(cartoon_sources);
        sources.extend(skin_sources);
        (sources, items, cartoons, gaussian_surfaces, skin_surfaces)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sync_frame_pins_a_hidden_structures_frame_and_catches_up_at_once() {
        assert_eq!(sync_frame(true, 0, 5), 5, "visible tracks the live frame");
        assert_eq!(
            sync_frame(false, 2, 9),
            2,
            "hidden stays pinned to its cached frame"
        );
        // Un-hiding jumps straight to wherever playback moved on to: one
        // step, not one per frame missed while hidden.
        assert_eq!(sync_frame(true, 2, 40), 40);
    }

    #[test]
    fn a_rep_draws_one_conformer_unless_the_altloc_policy_says_all() {
        use vv_core::altloc::AltlocPolicy;
        use vv_scene::{Command, CommandHistory};
        let mut scene = Scene::new();
        let mut history = CommandHistory::new(10);
        let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../fixtures/small/1AKE.pdb");
        history
            .dispatch(&mut scene, Command::LoadStructure { path })
            .unwrap();
        let id = scene.structures().next().unwrap().0;
        let count = |scene: &Scene, selection: &str| {
            let loaded = scene.structure(id).unwrap();
            let mut rep = loaded.rep().clone();
            rep.selection = selection.to_string();
            select_atoms(loaded, &rep, 0)
                .unwrap()
                .map(|a| a.len())
                .unwrap_or(loaded.structure.atom_count())
        };
        let total = scene.structure(id).unwrap().structure.atom_count();
        assert!(count(&scene, "all") < total);
        history
            .dispatch(
                &mut scene,
                Command::SetAltloc {
                    id,
                    policy: AltlocPolicy::All,
                },
            )
            .unwrap();
        assert_eq!(count(&scene, "all"), total);
    }

    /// Four alanines whose second CA has conformers A (0.3, listed first,
    /// off the chain) and B (0.7, on it): `First` shows B.
    fn backbone_with_a_split_ca() -> (Scene, vv_scene::StructureId) {
        let atom = |serial: u32, alt: char, seq: u32, x: f32, y: f32, occ: f32| {
            format!(
                "ATOM  {serial:>5}  CA {alt}ALA A{seq:>4}    {x:>8.3}{y:>8.3}{:>8.3}{occ:>6.2}{:>6.2}           C\n",
                0.0, 0.0
            )
        };
        let text = [
            atom(1, ' ', 1, 0.0, 0.0, 1.0),
            atom(2, 'A', 2, 3.8, 2.0, 0.3),
            atom(3, 'B', 2, 3.8, 0.0, 0.7),
            atom(4, ' ', 3, 7.6, 0.0, 1.0),
            atom(5, ' ', 4, 11.4, 0.0, 1.0),
        ]
        .concat();
        static NEXT: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
        let n = NEXT.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        let name = format!("vv_split_ca_{}_{n}.pdb", std::process::id());
        let path = std::env::temp_dir().join(name);
        std::fs::write(&path, text).unwrap();
        let mut scene = Scene::new();
        vv_scene::CommandHistory::new(10)
            .dispatch(&mut scene, vv_scene::Command::LoadStructure { path })
            .unwrap();
        let id = scene.structures().next().unwrap().0;
        (scene, id)
    }

    fn all_keep(loaded: &LoadedStructure) -> Option<Vec<bool>> {
        let atoms = select_atoms(loaded, loaded.rep(), 0).unwrap();
        keep_mask(&atoms, loaded.structure.atom_count())
    }

    #[test]
    fn a_tube_traces_the_shown_conformer_of_a_split_backbone_atom() {
        let (scene, id) = backbone_with_a_split_ca();
        let loaded = scene.structure(id).unwrap();
        let (plan, spline) = tube_mesh(loaded, loaded.rep(), &all_keep(loaded), 0);
        assert_eq!(plan.spans.len(), 3, "no residue drops out of the tube");
        assert_eq!(spline.controls[1], Vec3::new(3.8, 0.0, 0.0));
    }

    #[test]
    fn a_cartoon_traces_the_shown_conformer_of_a_split_backbone_atom() {
        let (scene, id) = backbone_with_a_split_ca();
        let loaded = scene.structure(id).unwrap();
        let drawn = loaded.drawn_positions(0);
        let (plan, spline) = plan_cartoon(&loaded.structure, drawn.positions());
        let keep = loaded.trace_keep(&all_keep(loaded)).unwrap();
        let drawn_plan = plan.filter(|a| keep[a as usize]);
        assert_eq!(drawn_plan.spans.len(), plan.spans.len());
        assert_eq!(spline.controls[1], Vec3::new(3.8, 0.0, 0.0));
    }

    /// Two guanines whose first P and N1 are conformer A (0.3, listed
    /// first, hidden under `First`) with a B twin (0.7) elsewhere.
    fn dna_with_split_atoms() -> (Scene, vv_scene::StructureId) {
        let atom = |serial: u32, name: &str, alt: char, seq: u32, x: f32, occ: f32| {
            format!(
                "ATOM  {serial:>5} {name:<4}{alt} DG A{seq:>4}    {x:>8.3}{:>8.3}{:>8.3}{occ:>6.2}{:>6.2}           C
",
                0.0, 0.0, 0.0
            )
        };
        let text = [
            atom(1, "P", 'A', 1, 90.0, 0.3),
            atom(2, "P", 'B', 1, 1.0, 0.7),
            atom(3, "N9", ' ', 1, 2.0, 1.0),
            atom(4, "N1", 'A', 1, 91.0, 0.3),
            atom(5, "N1", 'B', 1, 3.0, 0.7),
            atom(6, "P", ' ', 2, 6.0, 1.0),
            atom(7, "N9", ' ', 2, 7.0, 1.0),
            atom(8, "N1", ' ', 2, 8.0, 1.0),
        ]
        .concat();
        let path = std::env::temp_dir().join(format!("vv_split_dna_{}.pdb", std::process::id()));
        std::fs::write(&path, text).unwrap();
        let mut scene = Scene::new();
        vv_scene::CommandHistory::new(10)
            .dispatch(&mut scene, vv_scene::Command::LoadStructure { path })
            .unwrap();
        let id = scene.structures().next().unwrap().0;
        (scene, id)
    }

    #[test]
    fn a_ladder_rung_uses_the_shown_conformer_of_a_split_base() {
        let (scene, id) = dna_with_split_atoms();
        let loaded = scene.structure(id).unwrap();
        let rungs = bases::stick_rungs(loaded, &all_keep(loaded));
        assert_eq!(rungs, vec![[1, 4], [5, 7]]);
        let x = loaded.structure.frame(0).positions()[rungs[0][1] as usize].x;
        assert_eq!(x, 3.0, "the hidden N1 at 91 is not used");
    }

    #[test]
    fn a_hidden_conformer_is_not_in_any_drawn_atom_list() {
        let (scene, id) = backbone_with_a_split_ca();
        let loaded = scene.structure(id).unwrap();
        let drawn = select_atoms(loaded, loaded.rep(), 0).unwrap().unwrap();
        assert!(!drawn.contains(&1), "conformer A of residue 2 is hidden");
        assert!(drawn.contains(&2));
    }

    #[test]
    fn a_hidden_conformer_is_not_drawn_or_pickable_in_any_rep() {
        use vv_scene::{Command, CommandHistory, Representation as R};
        let Some(ctx) = gpu_context() else { return };
        let renderer = vv_render::Renderer::new(ctx.clone(), 64, 64);
        let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../fixtures/small/1AKE.pdb");
        let mut scene = Scene::new();
        let mut history = CommandHistory::new(10);
        history
            .dispatch(&mut scene, Command::LoadStructure { path })
            .unwrap();
        let id = scene.structures().next().unwrap().0;
        let rep = scene.structure(id).unwrap().reps[0].id;
        let hidden: Vec<u32> = {
            let loaded = scene.structure(id).unwrap();
            let shown = loaded.shown_atoms().expect("1AKE has conformers");
            shown.zeroes().map(|a| a as u32).collect()
        };
        for representation in [R::Spacefill, R::Sticks, R::Lines, R::BallAndStick] {
            let set = Command::SetRepresentation {
                id,
                rep,
                representation,
            };
            history.dispatch(&mut scene, set).unwrap();
            let mut cache = GpuCache::default();
            cache.sync(&scene, &ctx, &renderer, false);
            let (sources, ..) = cache.draw_items(&scene);
            for source in &sources {
                let atoms = source.atom_map.iter().flat_map(|m| m.iter().copied());
                let bonded = source
                    .bond_atoms
                    .iter()
                    .flat_map(|p| p.iter().flatten().copied());
                assert!(
                    atoms.chain(bonded).all(|a| !hidden.contains(&a)),
                    "{representation:?} draws a hidden conformer"
                );
            }
        }
    }

    /// Playback (`Scene::set_frame_live`) keeps advancing a hidden
    /// structure's own `frame` -- `sync`'s per-structure checks, driven by
    /// `sync_frame`, never see any of that motion until it is shown again,
    /// at which point they see it all at once.
    #[test]
    fn a_hidden_structures_document_frame_advances_while_its_gpu_sync_stays_pinned() {
        use vv_scene::{Command, CommandHistory};
        let path = |name: &str| {
            std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
                .join("../../fixtures/small")
                .join(name)
        };
        let mut scene = Scene::new();
        let mut history = CommandHistory::new(10);
        history
            .dispatch(
                &mut scene,
                Command::LoadTrajectory {
                    topology: path("1CRN.pdb"),
                    trajectory: path("1CRN_traj.dcd"),
                },
            )
            .unwrap();
        let id = scene.structures().next().unwrap().0;
        history
            .dispatch(&mut scene, Command::ShowStructure { id, visible: false })
            .unwrap();

        // The frame `sync` last recorded for this structure before it was
        // hidden (0, since it was never shown at any other frame here).
        let synced = 0;
        let last_frame = scene.structure(id).unwrap().structure.frame_count() - 1;
        assert!(last_frame >= 2, "fixture needs at least 3 frames");
        for live in 1..=last_frame {
            scene.set_frame_live(id, live);
            let loaded = scene.structure(id).unwrap();
            assert_eq!(
                loaded.frame, live,
                "playback still advances the document frame"
            );
            assert_eq!(
                sync_frame(loaded.visible, synced, loaded.frame),
                synced,
                "but a hidden structure's gpu_cache sync stays pinned"
            );
        }

        history
            .dispatch(&mut scene, Command::ShowStructure { id, visible: true })
            .unwrap();
        let loaded = scene.structure(id).unwrap();
        assert_eq!(
            sync_frame(loaded.visible, synced, loaded.frame),
            last_frame,
            "showing it again catches up to the current frame in one step"
        );
    }

    /// A structure loaded twice under different ids frames identically
    /// either time: `frame_structure` only looks at the one it's given,
    /// unlike the old union-of-everything `frame_all` it replaced.
    #[test]
    fn frame_structure_pivots_on_its_own_centroid_not_the_other_structures() {
        use vv_scene::{Command, CommandHistory};
        let mut scene = Scene::new();
        let mut history = CommandHistory::new(10);
        let path = |name: &str| {
            std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
                .join("../../fixtures/small")
                .join(name)
        };
        history
            .dispatch(
                &mut scene,
                Command::LoadStructure {
                    path: path("1CRN.pdb"),
                },
            )
            .unwrap();
        history
            .dispatch(
                &mut scene,
                Command::LoadStructure {
                    path: path("4HHB.cif"),
                },
            )
            .unwrap();
        let ids: Vec<StructureId> = scene.structures().map(|(id, _)| id).collect();
        let [crn, hhb] = ids[..] else {
            panic!("2 structures")
        };

        let camera = frame_structure(&scene, hhb).unwrap();
        let centroid = scene
            .structure(hhb)
            .unwrap()
            .structure
            .frame(0)
            .centroid()
            .unwrap();
        assert_eq!(camera.target, centroid);
        // Framing 4HHB alone must not be pulled toward 1CRN's atoms.
        let crn_centroid = scene
            .structure(crn)
            .unwrap()
            .structure
            .frame(0)
            .centroid()
            .unwrap();
        assert!(camera.target.distance(crn_centroid) > 10.0);
    }

    #[test]
    fn companions_of_a_cartoon_are_hemes_only_until_additives_and_water_are_asked_for() {
        let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../fixtures/small/4HHB.cif");
        let structure = vv_io::load(path).unwrap();
        let bonds = vv_core::bonds::perceive(&structure.topology, structure.frame(0).positions());
        let t = &structure.topology;
        let mut rep = Rep::new(RepId(0), SceneRepresentation::Cartoon);
        let trace: Vec<u32> = vv_core::backbone_trace(t, structure.frame(0).positions())
            .segments
            .concat();
        let drawn = companions::Drawn {
            keep: &None,
            traced: &trace,
        };
        let names_of = |rep: &Rep| -> Vec<(vv_core::Companion, Vec<String>)> {
            companions::subsets(&structure, &bonds.pairs, rep, drawn)
                .into_iter()
                .map(|(kind, s)| {
                    let mut names: Vec<String> = s
                        .atoms
                        .iter()
                        .map(|&a| t.residue_name(t.residue_index[a as usize] as usize).into())
                        .collect();
                    names.dedup();
                    (kind, names)
                })
                .collect()
        };
        let found = names_of(&rep);
        assert_eq!(found.len(), 1, "{found:?}");
        assert_eq!(
            found[0],
            (vv_core::Companion::Ligand, vec!["HEM".to_string()])
        );
        rep.options.insert("additives".into(), 1.0);
        rep.options.insert("water".into(), 1.0);
        let kinds: Vec<_> = names_of(&rep).into_iter().map(|(k, _)| k).collect();
        assert!(kinds.contains(&vv_core::Companion::Additive));
        assert!(kinds.contains(&vv_core::Companion::Water));
        rep.options.insert("ligands".into(), 0.0);
        assert!(!names_of(&rep)
            .iter()
            .any(|(k, _)| *k == vv_core::Companion::Ligand));
        rep.options.insert("ligands".into(), 1.0);
        let subset = companions::subsets(&structure, &bonds.pairs, &rep, drawn)
            .into_iter()
            .find(|(k, _)| *k == vv_core::Companion::Ligand)
            .unwrap()
            .1;
        // 4 hemes of 43 atoms each.
        assert_eq!(subset.atoms.len(), 4 * 43);
        assert!(
            subset.local_bonds.len() > 4 * 40,
            "{}",
            subset.local_bonds.len()
        );
        assert_eq!(subset.local_bonds.len(), subset.real_pairs.len());
        for ([la, lb], [ra, rb]) in subset.local_bonds.iter().zip(&subset.real_pairs) {
            assert_eq!(subset.atoms[*la as usize], *ra);
            assert_eq!(subset.atoms[*lb as usize], *rb);
        }
    }

    /// A selection of a run of residues, one lone residue, a heme and the
    /// waters, drawn as a cartoon: the run is ribbon, the lone residue (no
    /// neighbour to span to) and the heme are licorice, and the selected
    /// waters are licorice too because the selection names them.
    #[test]
    fn a_mixed_selection_draws_what_the_ribbon_cannot_as_licorice() {
        let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../fixtures/small/4HHB.cif");
        let structure = vv_io::load(path).unwrap();
        let bonds = vv_core::bonds::perceive(&structure.topology, structure.frame(0).positions());
        let loaded = LoadedStructure::new(structure, None, "4HHB".into(), bonds);
        let t = &loaded.structure.topology;
        let mut rep = Rep::new(RepId(0), SceneRepresentation::Cartoon);
        rep.selection = "chain A and (resid 10-14 50 or resname HEM or water)".into();
        let selected = loaded.select(&rep.selection, 0).unwrap();
        let keep: Option<Vec<bool>> =
            Some((0..t.atom_count()).map(|a| selected.contains(a)).collect());

        let coords = loaded.structure.frame(0);
        let codes = vv_core::cartoon::secondary_structure(t, coords.positions(), true);
        let plan = vv_core::cartoon::plan(t, coords.positions(), &codes);
        let mask = keep.clone().unwrap();
        let plan = plan.filter(|a| mask[a as usize]);
        assert!(!plan.recipes.is_empty(), "the run of residues is a ribbon");

        let traced = companions::traced_atoms(&plan);
        let drawn = companions::Drawn {
            keep: &keep,
            traced: &traced,
        };
        let pieces = companions::subsets(&loaded.structure, &loaded.bonds.pairs, &rep, drawn);
        let names_of = |kind| -> Vec<String> {
            let found = pieces.iter().find(|(k, _)| *k == kind);
            let mut names: Vec<String> = found
                .map(|(_, s)| &s.atoms)
                .into_iter()
                .flatten()
                .map(|&a| t.residue_name(t.residue_index[a as usize] as usize).into())
                .collect();
            names.dedup();
            names
        };
        assert_eq!(names_of(vv_core::Companion::Ligand), ["HEM"]);
        assert_eq!(names_of(vv_core::Companion::Untraced).len(), 1);
        assert_eq!(names_of(vv_core::Companion::Water), ["HOH"]);
        let ribbon_residues = (10..=14).count();
        let traced_residues: std::collections::HashSet<u32> = traced
            .iter()
            .map(|&a| t.residue_index[a as usize])
            .collect();
        assert_eq!(traced_residues.len(), ribbon_residues);
    }

    /// Regression for the tube's `radius_by` (putty) option: the plan's
    /// sections must actually carry a varying radius when it's on, and a
    /// flat one when it's off (`tube_radius_fn`, `vv_core::cartoon::
    /// tube_plan`).
    #[test]
    fn tube_radius_by_varies_with_b_factor_and_constant_does_not() {
        let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../fixtures/small/1UBQ.cif");
        let structure = vv_io::load(path).unwrap();
        let bonds = vv_core::bonds::perceive(&structure.topology, structure.frame(0).positions());
        let loaded = LoadedStructure::new(structure, None, "1UBQ".into(), bonds);
        let mut rep = Rep::new(RepId(0), SceneRepresentation::Tube);

        rep.options.insert("radius_by".into(), 1.0);
        let (putty, _) = tube_mesh(&loaded, &rep, &None, 0);
        let (min, max) = putty
            .recipes
            .iter()
            .map(|r| r.half_width)
            .fold((f32::MAX, f32::MIN), |(lo, hi), r| (lo.min(r), hi.max(r)));
        assert!(
            max - min > 0.05,
            "putty radius should vary with B-factor, got {min}..{max}"
        );

        rep.options.insert("radius_by".into(), 0.0);
        let (constant, _) = tube_mesh(&loaded, &rep, &None, 0);
        let first = constant.recipes[0].half_width;
        assert!(
            constant
                .recipes
                .iter()
                .all(|r| (r.half_width - first).abs() < 1e-6),
            "constant radius should not vary"
        );
    }

    /// `Set VIZVIZ_TEST_ADAPTER=<name substring>` to run this on a
    /// specific GPU; skips (passes) with no usable adapter, like
    /// `vv-render`'s own headless tests.
    fn gpu_context() -> Option<Arc<GpuContext>> {
        let filter = std::env::var("VIZVIZ_TEST_ADAPTER").ok();
        let instance = GpuContext::instance();
        GpuContext::with_instance(instance, None, filter.as_deref()).ok()
    }

    /// `single_bond_pairs` drops exactly the bonds `orders` names and
    /// nothing else -- the invariant that keeps the ordinary center
    /// cylinder from ever drawing a non-`Single` bond a `StrandGeometry`
    /// is about to draw as strands (else a double bond would show 3
    /// overlapping cylinders instead of 2 separated ones).
    #[test]
    fn single_bond_pairs_drops_exactly_the_non_single_bonds() {
        let bonds = vv_core::BondTable {
            pairs: vec![[0, 1], [1, 2], [2, 3], [3, 4]],
            orders: vec![
                (1, vv_core::BondOrder::Double),
                (3, vv_core::BondOrder::Triple),
            ],
        };
        let singles = single_bond_pairs(&bonds);
        assert_eq!(singles.len(), bonds.pairs.len() - bonds.orders.len());
        assert_eq!(singles, vec![[0, 1], [2, 3]]);
        for &[a, b] in &singles {
            let i = bonds.pairs.binary_search(&[a, b]).unwrap();
            assert_eq!(bonds.order_of(i as u32), vv_core::BondOrder::Single);
        }
    }

    /// A double bond's strand cylinder picks back to its real atom pair,
    /// through the exact same `Pick::Bond { bond }` -> `DrawSource::
    /// bond_atoms` path a selection-derived ligand stick uses (`ui.rs`,
    /// `selection_outline.rs`) -- `build_strand_geometry` is what feeds it.
    #[test]
    fn a_double_bonds_strand_picks_back_to_its_real_atom_pair() {
        let Some(ctx) = gpu_context() else { return };
        let mut renderer = vv_render::Renderer::new(ctx.clone(), 64, 64);
        let bonds = vv_core::BondTable {
            pairs: vec![[0, 1]],
            orders: vec![(0, vv_core::BondOrder::Double)],
        };
        let adjacency = vv_core::adjacency(&bonds, 2);
        let positions = [Vec3::new(-0.75, 0.0, 0.0), Vec3::new(0.75, 0.0, 0.0)];
        let colors = [0xFFFFFFFFu32, 0xFFFFFFFFu32];
        // A tiny separation next to the strand radius: the two strands
        // overlap on screen, so the framed center pixel is guaranteed to
        // land on one of them regardless of which side of the (arbitrary,
        // no-neighbour) offset axis the camera ends up looking along.
        let strands = build_strand_geometry(
            &ctx, &renderer, &bonds, &adjacency, &positions, 0.02, &colors,
        )
        .expect("a double bond has strands");
        assert_eq!(strands.gpu.bond_count, 2, "double bond: 2 strands");

        let camera = vv_render::Camera::framing(strands.gpu.center, strands.gpu.radius.max(1.0));
        let settings = vv_render::RenderSettings {
            representation: vv_render::Representation::BallAndStick,
            occlusion_culling: false,
            ..Default::default()
        };
        let mut encoder = ctx.device.create_command_encoder(&Default::default());
        renderer.render(
            &mut encoder,
            &camera,
            &strands.gpu,
            &strands.bindings,
            &settings,
        );
        ctx.queue.submit([encoder.finish()]);
        renderer.after_submit();

        let pick = renderer.pick(32, 32);
        let vv_render::Pick::Bond { item: 0, bond } = pick.expect("center pixel hits a strand")
        else {
            panic!("expected a bond pick, got {pick:?}");
        };
        assert_eq!(
            strands.bond_atoms[bond as usize],
            [0, 1],
            "the strand's local bond maps back to the real C=C pair"
        );
    }

    /// The wiring `sync` builds, not just the helpers it's built from:
    /// loading a real structure with double bonds (1CRN has PHE/TYR rings
    /// and every backbone C=O) and drawing it ball-and-stick, the shared
    /// `entry.gpu`'s own bond count must be *smaller* than the full bond
    /// table -- non-`Single` bonds excluded, drawn only by their strands
    /// -- and `entry.single_pairs` must account for exactly the
    /// difference. This is the exact check that would have caught
    /// `entry.gpu` being built from the unfiltered `loaded.bonds.pairs`.
    #[test]
    fn sync_filters_non_single_bonds_out_of_the_shared_geometry() {
        use vv_scene::{Command, CommandHistory, Representation as SceneRepresentation, Scene};
        let Some(ctx) = gpu_context() else { return };
        let renderer = vv_render::Renderer::new(ctx.clone(), 64, 64);
        let mut scene = Scene::new();
        let mut history = CommandHistory::new(10);
        let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../fixtures/small/1CRN.cif");
        history
            .dispatch(&mut scene, Command::LoadStructure { path })
            .unwrap();
        let id = scene.structures().next().unwrap().0;
        let rep_id = scene.structure(id).unwrap().reps[0].id;
        history
            .dispatch(
                &mut scene,
                Command::SetRepresentation {
                    id,
                    rep: rep_id,
                    representation: SceneRepresentation::BallAndStick,
                },
            )
            .unwrap();

        let loaded = scene.structure(id).unwrap();
        assert!(
            !loaded.bonds.orders.is_empty(),
            "1CRN's PHE/TYR rings and backbone carbonyls must give it some"
        );

        let mut cache = GpuCache::default();
        cache.sync(&scene, &ctx, &renderer, false);

        let entry = cache.entries.get(&id).expect("synced");
        let gpu = entry
            .gpu
            .as_ref()
            .expect("ball-and-stick uploads all atoms");
        let single_pairs = entry
            .single_pairs
            .as_ref()
            .expect("some bond is non-Single, so this must be Some");
        assert_eq!(gpu.bond_table_len, single_pairs.len());
        assert!(
            single_pairs.len() < loaded.bonds.pairs.len(),
            "the shared geometry must have fewer bonds than the full table"
        );
        assert_eq!(
            single_pairs.len(),
            loaded.bonds.pairs.len() - loaded.bonds.orders.len(),
            "exactly the non-Single bonds are excluded"
        );
    }

    /// The same check as `sync_filters_non_single_bonds_out_of_the_
    /// shared_geometry`, but for a selection layer (`RepGeometry::
    /// SomeAtoms`): PHE's ring and backbone carbonyl give it non-Single
    /// bonds, so a `resname PHE` ball-and-stick layer must draw fewer
    /// bonds on its own `Derived::gpu` than the residue actually has, and
    /// carry a `Derived::strands` for the rest, still picking back to the
    /// real atom pairs via `Derived::bond_atoms`.
    #[test]
    fn sync_filters_non_single_bonds_out_of_a_selection_layer() {
        use vv_scene::{Command, CommandHistory, Representation as SceneRepresentation, Scene};
        let Some(ctx) = gpu_context() else { return };
        let renderer = vv_render::Renderer::new(ctx.clone(), 64, 64);
        let mut scene = Scene::new();
        let mut history = CommandHistory::new(10);
        let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../fixtures/small/1CRN.cif");
        history
            .dispatch(&mut scene, Command::LoadStructure { path })
            .unwrap();
        let id = scene.structures().next().unwrap().0;
        let rep_id = scene.structure(id).unwrap().reps[0].id;
        history
            .dispatch(
                &mut scene,
                Command::SetRepresentation {
                    id,
                    rep: rep_id,
                    representation: SceneRepresentation::BallAndStick,
                },
            )
            .unwrap();
        history
            .dispatch(
                &mut scene,
                Command::SetRepSelection {
                    id,
                    rep: rep_id,
                    selection: "resname PHE".into(),
                },
            )
            .unwrap();

        let mut cache = GpuCache::default();
        cache.sync(&scene, &ctx, &renderer, false);

        let entry = cache.entries.get(&id).expect("synced");
        let r = entry.reps.get(&rep_id).expect("the rep synced");
        let RepGeometry::SomeAtoms(d) = &r.geometry else {
            panic!("expected a selection layer, got {:?}", r.geometry_id());
        };
        let strands = d
            .strands
            .as_ref()
            .expect("PHE's ring and C=O are non-Single");
        assert!(!strands.strands.is_empty());
        // Every strand's real atom pair is one of PHE's bonds, confirming
        // `strands_within` kept only bonds inside this selection.
        for pair in strands.bond_atoms.iter() {
            assert!(
                d.atom_map.contains(&pair[0]) && d.atom_map.contains(&pair[1]),
                "strand {pair:?} must be between two atoms of the selection"
            );
        }
        // The main geometry's own bonds are all Single: none of them is
        // also one of the strands' real pairs.
        for &pair in d.bond_atoms.iter() {
            assert!(
                !strands.bond_atoms.contains(&pair),
                "{pair:?} drawn both as an ordinary cylinder and as strands"
            );
        }
    }

    /// Loads 1CRN, makes its rep `representation` with `surface` set to
    /// `display` (1 mesh, 2 dots), and syncs until its background builds
    /// are done.
    fn synced_surface(
        ctx: &Arc<GpuContext>,
        representation: SceneRepresentation,
        display: f32,
    ) -> (Scene, StructureId, GpuCache) {
        use vv_scene::{Command, CommandHistory};
        let renderer = vv_render::Renderer::new(ctx.clone(), 64, 64);
        let mut scene = Scene::new();
        let mut history = CommandHistory::new(10);
        let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../fixtures/small/1CRN.cif");
        history
            .dispatch(&mut scene, Command::LoadStructure { path })
            .unwrap();
        let id = scene.structures().next().unwrap().0;
        let rep = scene.structure(id).unwrap().reps[0].id;
        let set = |name: &str, value: f32| Command::SetRepOption {
            id,
            rep,
            name: name.into(),
            value: Some(value),
        };
        history
            .dispatch(
                &mut scene,
                Command::SetRepresentation {
                    id,
                    rep,
                    representation,
                },
            )
            .unwrap();
        history
            .dispatch(&mut scene, set("surface", display))
            .unwrap();
        let mut cache = GpuCache::default();
        for _ in 0..600 {
            cache.sync(&scene, ctx, &renderer, false);
            if !cache.building() {
                break;
            }
            std::thread::sleep(std::time::Duration::from_millis(20));
        }
        (scene, id, cache)
    }

    /// Each surface drawn as `mesh` is a triangle soup of edges whose every
    /// vertex names a real atom of the structure, and as `dots` a point
    /// per mesh vertex doing the same; the solid surface lists none.
    #[test]
    fn every_surface_draws_as_a_mesh_or_dots_that_picks_real_atoms() {
        let Some(ctx) = gpu_context() else { return };
        for representation in [
            SceneRepresentation::GaussianSurface,
            SceneRepresentation::SkinSurface,
            SceneRepresentation::Ses,
        ] {
            let (scene, id, cache) = synced_surface(&ctx, representation, 1.0);
            let atoms = scene.structure(id).unwrap().structure.atom_count() as u32;
            let (sources, items, cartoons, ..) = cache.draw_items(&scene);
            assert!(
                items.is_empty(),
                "{representation:?}: mesh edges are not atoms"
            );
            assert_eq!(cartoons.len(), 1, "{representation:?}");
            let vv_render::CartoonMesh::Glycan(soup) = cartoons[0].mesh else {
                panic!("{representation:?}: edges draw as a triangle soup");
            };
            assert!(soup.vertex_count > 3_000, "{representation:?}");
            assert_eq!(soup.vertex_count % 3, 0);
            assert_eq!(soup.display, vv_render::MeshDisplay::Lines(1.5));
            let source = sources.last().unwrap().atom_map.as_ref().unwrap();
            assert_eq!(source.len(), soup.vertex_count as usize);
            assert!(source.iter().all(|&a| a < atoms), "{representation:?}");

            let (scene, _, cache) = synced_surface(&ctx, representation, 2.0);
            let (sources, items, cartoons, ..) = cache.draw_items(&scene);
            assert!(cartoons.is_empty(), "{representation:?}: dots are atoms");
            assert_eq!(items.len(), 1);
            let map = sources[0].atom_map.as_ref().unwrap();
            assert!(map.len() > 1_000 && map.iter().all(|&a| a < atoms));

            let (scene, _, cache) = synced_surface(&ctx, representation, 0.0);
            let (_, items, cartoons, ..) = cache.draw_items(&scene);
            assert!(items.is_empty() && cartoons.is_empty());
        }
    }

    #[test]
    fn cylinder_helices_are_a_cartoon_option_that_straightens_helices() {
        let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../fixtures/small/1CRN.cif");
        let structure = vv_io::load(path).unwrap();
        let positions = structure.frame(0).positions().to_vec();
        let (plan, spline) = plan_cartoon(&structure, &positions);
        let mut rep = Rep::new(RepId(0), SceneRepresentation::Cartoon);
        let samples = vv_core::cartoon::SAMPLES_PER_RESIDUE;
        assert!(cylinder_cartoon(&rep, &plan, &spline, &positions, samples).is_none());
        rep.options.insert("helix".into(), 1.0);
        let (styled, frame) = cylinder_cartoon(&rep, &plan, &spline, &positions, samples).unwrap();
        let widest = |plan: &CartoonPlan| {
            plan.recipes
                .iter()
                .map(|r| r.half_width)
                .fold(0.0, f32::max)
        };
        assert!(widest(&styled) > widest(&plan), "helices became cylinders");
        assert_eq!(frame.controls.len(), spline.controls.len());
        assert_ne!(
            frame.controls, spline.controls,
            "controls moved to the axes"
        );
    }
}
