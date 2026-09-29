//! Structures uploaded to the GPU as pages of storage buffers.

use std::sync::atomic::{AtomicBool, Ordering};

use bytemuck::{Pod, Zeroable};
use glam::Vec3;
use rayon::prelude::*;
use vv_core::{BondTable, CoordSet, Structure, Topology};
use wgpu::util::DeviceExt;

use crate::color::{self, ColorScheme};
use crate::context::GpuContext;
use crate::style::Material;

/// Packed RGBA8 colors for every atom in `topology`, in atom order, under
/// `scheme`. Shared by `GpuStructure::upload` (first paint) and `recolor`
/// (switching schemes without re-uploading positions/radii) — and, being
/// plain CPU math with no wgpu types, also by `vv-cpu`, so the two render
/// backends can never disagree on what a given scheme's colors are.
pub fn colors_for(scheme: ColorScheme, topology: &Topology) -> Vec<u32> {
    match scheme {
        ColorScheme::Element => topology
            .element
            .iter()
            .map(|e| color::by_element(*e))
            .collect(),
        ColorScheme::Chain => {
            let by_name = topology.chain_name_index();
            (0..topology.atom_count())
                .map(|a| color::by_chain(by_name[topology.chain_of_atom(a) as usize]))
                .collect()
        }
        ColorScheme::BFactor => colors_from_scalar(&topology.b_factor),
        ColorScheme::SecondaryStructure => {
            let codes: Vec<vv_core::dssp::DsspCode> = topology
                .residues
                .iter()
                .map(|r| match r.ss {
                    vv_core::SecondaryStructure::Helix => vv_core::dssp::DsspCode::AlphaHelix,
                    vv_core::SecondaryStructure::Strand => vv_core::dssp::DsspCode::Strand,
                    _ => vv_core::dssp::DsspCode::Coil,
                })
                .collect();
            colors_for_ss(topology, &codes)
        }
        ColorScheme::ResidueType => (0..topology.atom_count())
            .map(|a| {
                let r = topology.residue_index[a] as usize;
                color::by_residue_type(topology.residue_name(r))
                    .unwrap_or_else(|| color::by_element(topology.element[a]))
            })
            .collect(),
        ColorScheme::Rainbow => {
            // Position of each polymer residue along its chain, 0 at the
            // first; waters and ligands (often in the same chain) keep
            // element colors instead of stretching the ramp.
            let polymer = |r: usize| {
                matches!(
                    topology.residue_class(r),
                    vv_core::ResidueClass::Protein | vv_core::ResidueClass::Nucleic
                )
            };
            let mut t = vec![None; topology.residues.len()];
            for chain in &topology.chains {
                let members: Vec<u32> = chain
                    .residues
                    .clone()
                    .filter(|&r| polymer(r as usize))
                    .collect();
                let n = members.len().max(2) as f32 - 1.0;
                for (k, r) in members.into_iter().enumerate() {
                    t[r as usize] = Some(k as f32 / n);
                }
            }
            (0..topology.atom_count())
                .map(|a| match t[topology.residue_index[a] as usize] {
                    Some(t) => color::rainbow(t),
                    None => color::by_element(topology.element[a]),
                })
                .collect()
        }
        ColorScheme::ResidueName => (0..topology.atom_count())
            .map(|a| {
                color::by_residue_name(topology.residue_name(topology.residue_index[a] as usize))
            })
            .collect(),
        ColorScheme::Occupancy => colors_from_scalar(&topology.occupancy),
        ColorScheme::Hydrophobicity => by_residue_scale(topology, color::kyte_doolittle, |v| {
            color::by_hydropathy(v, color::KD_SCALE)
        }),
        ColorScheme::WimleyWhite => by_residue_scale(topology, color::wimley_white, |v| {
            color::by_hydropathy(-v, color::WW_SCALE)
        }),
        ColorScheme::HelixPropensity => by_residue_scale(
            topology,
            color::helix_propensity,
            color::by_helix_propensity,
        ),
        ColorScheme::StrandPropensity => by_residue_scale(
            topology,
            color::strand_propensity,
            color::by_strand_propensity,
        ),
        ColorScheme::TurnPropensity => {
            by_residue_scale(topology, color::turn_propensity, color::by_turn_propensity)
        }
        ColorScheme::BuriedIndex => {
            by_residue_scale(topology, color::buried_index, color::by_buried_index)
        }
        ColorScheme::Zappo => by_residue_lookup(topology, color::by_zappo),
        ColorScheme::Taylor => by_residue_lookup(topology, color::by_taylor),
        ColorScheme::Clustal => by_residue_lookup(topology, color::by_clustal),
        ColorScheme::Nucleotide => by_residue_lookup(topology, color::by_nucleotide),
        ColorScheme::PurinePyrimidine => by_residue_lookup(topology, color::by_purine_pyrimidine),
        ColorScheme::Class => (0..topology.atom_count())
            .map(|a| color::by_class(topology.residue_class(topology.residue_index[a] as usize)))
            .collect(),
        ColorScheme::SegmentName => {
            let by_segid = !topology.segids.is_empty()
                && (0..topology.chain_count()).any(|c| !topology.segid(c).is_empty());
            let mut seen: Vec<u32> = Vec::new();
            let chain_segment: Vec<u32> = (0..topology.chain_count())
                .map(|c| {
                    let id = if by_segid {
                        topology.segids[c].0
                    } else {
                        topology.chains[c].auth_asym.0
                    };
                    seen.iter().position(|&s| s == id).unwrap_or_else(|| {
                        seen.push(id);
                        seen.len() - 1
                    }) as u32
                })
                .collect();
            (0..topology.atom_count())
                .map(|a| color::by_chain(chain_segment[topology.chain_of_atom(a) as usize]))
                .collect()
        }
        ColorScheme::Constant(c) => vec![c; topology.atom_count()],
        ColorScheme::Hetero => {
            let by_name = topology.chain_name_index();
            (0..topology.atom_count())
                .map(|a| {
                    if topology.element[a] == vv_core::Element::CARBON {
                        color::by_chain(by_name[topology.chain_of_atom(a) as usize])
                    } else {
                        color::by_element(topology.element[a])
                    }
                })
                .collect()
        }
    }
}

/// Every atom's residue name through `scale` (a Chou-Fasman-style table,
/// `None` for a non-standard residue), then `ramp`; atoms `scale` doesn't
/// cover keep their element colour.
fn by_residue_scale(
    topology: &Topology,
    scale: impl Fn(&str) -> Option<f32>,
    ramp: impl Fn(f32) -> u32,
) -> Vec<u32> {
    (0..topology.atom_count())
        .map(|a| {
            let r = topology.residue_index[a] as usize;
            scale(topology.residue_name(r))
                .map(&ramp)
                .unwrap_or_else(|| color::by_element(topology.element[a]))
        })
        .collect()
}

/// Every atom's residue name through a direct name -> colour lookup
/// (`None` for a residue it doesn't cover, which keeps its element colour).
fn by_residue_lookup(topology: &Topology, lookup: impl Fn(&str) -> Option<u32>) -> Vec<u32> {
    (0..topology.atom_count())
        .map(|a| {
            let r = topology.residue_index[a] as usize;
            lookup(topology.residue_name(r))
                .unwrap_or_else(|| color::by_element(topology.element[a]))
        })
        .collect()
}

/// Packed RGBA8 colors for `vv_scene::ColorScheme::Fragment`: one hue per connected component of `bonds`. Kept separate
/// from `colors_for` because it needs the structure's bond graph, which
/// `Topology` alone doesn't carry -- callers that have a `BondTable` at
/// hand (e.g. `vv-app`'s `gpu_cache`) call this directly instead, the same
/// way `colors_for_ss` takes a DSSP assignment `colors_for` alone can't
/// produce.
pub fn colors_for_fragments(topology: &Topology, bonds: &BondTable) -> Vec<u32> {
    bonds
        .fragments(topology.atom_count())
        .into_iter()
        .map(color::by_chain)
        .collect()
}

/// Secondary-structure colors from `codes` (one per residue, e.g. DSSP
/// on the current frame). Atoms outside protein residues (ligands,
/// water, ions) keep their element color, so they stand out.
pub fn colors_for_ss(topology: &Topology, codes: &[vv_core::dssp::DsspCode]) -> Vec<u32> {
    (0..topology.atom_count())
        .map(|a| {
            let r = topology.residue_index[a] as usize;
            if color::by_residue_type(topology.residue_name(r)).is_some() {
                color::by_secondary_structure(
                    codes
                        .get(r)
                        .copied()
                        .unwrap_or(vv_core::dssp::DsspCode::Coil),
                )
            } else {
                color::by_element(topology.element[a])
            }
        })
        .collect()
}

/// Packed RGBA8 colors for one per-atom scalar (B-factor, a SASA or
/// conservation score from another package, ...) on the blue-white-red
/// ramp, normalized to `[min, max]`. NaNs are ignored for the range and
/// drawn as the midpoint. The one function every "color by a number" path
/// goes through, so a value column from Python looks exactly like the
/// built-in B-factor coloring of the same numbers.
pub fn colors_from_scalar(values: &[f32]) -> Vec<u32> {
    let (min, max) = scalar_range(values);
    values
        .iter()
        .map(|&v| color::by_b_factor(if v.is_nan() { 0.5 * (min + max) } else { v }, min, max))
        .collect()
}

/// `colors_from_scalar` with an explicit range, for values that should
/// keep the same colors across frames or structures.
pub fn colors_from_scalar_in(values: &[f32], min: f32, max: f32) -> Vec<u32> {
    values
        .iter()
        .map(|&v| color::by_b_factor(if v.is_nan() { 0.5 * (min + max) } else { v }, min, max))
        .collect()
}

/// `(min, max)` over the finite entries of `values`; `(0, 0)` when there
/// are none, which `by_b_factor` draws as the uniform midpoint color.
pub fn scalar_range(values: &[f32]) -> (f32, f32) {
    let (min, max) = values
        .iter()
        .filter(|v| v.is_finite())
        .fold((f32::MAX, f32::MIN), |(lo, hi), &b| (lo.min(b), hi.max(b)));
    if min > max {
        (0.0, 0.0)
    } else {
        (min, max)
    }
}

/// Packs `(r, g, b, a)` bytes the way the color buffers store them.
pub fn pack_rgba(r: u8, g: u8, b: u8, a: u8) -> u32 {
    u32::from_le_bytes([r, g, b, a])
}

/// Inverse of `pack_rgba`, as normalized floats - what a vertex color
/// attribute needs (`shaders/atoms.wgsl`'s `unpack_color` does the same
/// thing in WGSL, for shaders that read packed colors from a buffer
/// instead of a vertex attribute).
pub fn unpack_rgba(c: u32) -> [f32; 4] {
    c.to_le_bytes().map(|b| b as f32 / 255.0)
}

/// A GPU allocation failed: the device ran out of memory (8GLV's cartoon
/// or Gaussian volume on an integrated GPU). Callers fall back to a
/// cheaper representation instead of the process panicking, which is
/// what an uncaptured wgpu error does.
#[derive(Debug, Clone, Copy, PartialEq, Eq, thiserror::Error)]
#[error("out of GPU memory")]
pub struct OutOfGpuMemory;

/// Runs `f` (resource creation) inside an out-of-memory error scope.
fn catch_oom<T>(device: &wgpu::Device, f: impl FnOnce() -> T) -> Result<T, OutOfGpuMemory> {
    let scope = device.push_error_scope(wgpu::ErrorFilter::OutOfMemory);
    let value = f();
    match pollster::block_on(scope.pop()) {
        None => Ok(value),
        Some(_) => Err(OutOfGpuMemory),
    }
}

/// `create_buffer_init`, but reporting out-of-memory. Not mapped at
/// creation: mapping a buffer that failed to allocate would itself be an
/// error outside the scope. `contents.len()` must be a multiple of 4.
pub(crate) fn try_buffer_init(
    ctx: &GpuContext,
    label: &str,
    contents: &[u8],
    usage: wgpu::BufferUsages,
) -> Result<wgpu::Buffer, OutOfGpuMemory> {
    debug_assert_eq!(contents.len() % 4, 0);
    let buffer = catch_oom(&ctx.device, || {
        ctx.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some(label),
            size: (contents.len() as u64).max(4),
            usage: usage | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        })
    })?;
    if !contents.is_empty() {
        ctx.queue.write_buffer(&buffer, 0, contents);
    }
    Ok(buffer)
}

/// Atoms per page. 2^22 atoms × 16 B = 64 MiB, under wgpu's default
/// 128 MiB storage-binding limit with room for the index buffers.
pub const PAGE_ATOMS: usize = 1 << 22;

/// `xyz` = position, `w` = van der Waals radius. Matches `array<vec4<f32>>`.
#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
pub struct AtomInstance {
    pub position: [f32; 3],
    pub radius: f32,
}

/// Three `DrawIndirectArgs` back to back: quads at 0, points at 16, bonds
/// at 32. All draws are non-instanced; the cull shaders accumulate vertex
/// counts (6 per quad or bond, 1 per point) and everything else is constant.
/// Draw args for occlusion phase 1 (quads, points, bond cylinders, bond
/// lines), then the same four for phase 2 at [`LATE_INDIRECT_OFFSET`]
/// (see shaders/cull.wgsl).
pub const INDIRECT_RESET: [u32; 32] = [
    0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, //
    0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0,
];
pub const LATE_INDIRECT_OFFSET: u64 = 64;
/// Consecutive atoms culled as one bounding sphere before per-atom tests.
/// File order keeps them spatially compact (a few residues).
pub const CLUSTER_ATOMS: u32 = 64;
pub const LINES_INDIRECT_OFFSET: u64 = 48;
pub const POINTS_INDIRECT_OFFSET: u64 = 16;
pub const BONDS_INDIRECT_OFFSET: u64 = 32;

/// Per-page shader parameters, rewritten by the renderer every frame.
/// Layout must match `PageParams` in `shaders/atoms.wgsl`.
///
/// These used to live in the camera uniform, which only works when one
/// structure is drawn per frame. With several structures composited into
/// one frame each needs its own representation, and every page needs its
/// own id offsets so the picking buffer can name any atom or bond in the
/// frame unambiguously (see `Renderer::pick`).
#[repr(C)]
#[derive(Clone, Copy, Debug, Default, PartialEq, Pod, Zeroable)]
pub struct PageParams {
    /// Multiplier on stored (van der Waals) atom radii.
    pub radius_scale: f32,
    /// Cylinder radius for bonds; `0` disables bonds.
    pub bond_radius: f32,
    /// Added to a page-relative atom index to get its frame-wide pick id.
    pub atom_id_base: u32,
    /// Added to a bond's structure-wide index to get its frame-wide pick id.
    pub bond_id_base: u32,
    /// Added to every drawn radius, after `radius_scale`.
    pub radius_offset: f32,
    /// Non-zero: every bond as a 1-px line.
    pub lines: u32,
    /// Size classification sees at least this radius: an SES fillet with
    /// a tiny bound stays a ray-cast quad, not a one-pixel point, while
    /// the atoms beside it still cut their caps open.
    pub min_quad_radius: f32,
    /// World-space view-depth margin for glass self-occlusion culling
    /// (`glass_cull_margin` below); `0` disables it (every opaque page,
    /// and any glass page too sparse to estimate a spacing from).
    pub glass_cull_margin: f32,
    pub material: Material,
}

/// World-space depth margin beyond which a glass atom's own weighted-
/// blended OIT contribution is below 1/255 once composited, so
/// `shaders/cull.wgsl`'s `glass_hidden` can cull it against the previous
/// frame's nearest-glass depth (docs/RENDERING.md's "Transparency": glass
/// self-occlusion culling). `K(opacity)` layers of transmittance `(1 -
/// opacity)` cross that threshold at `K = ceil(ln(1/255) / ln(1 -
/// opacity))`; `spacing` is the mean nearest-neighbour distance implied by
/// packing `atom_count` atoms into a sphere of `radius` (`(volume /
/// atom_count) ^ (1/3)`), an isotropic-density approximation of its own.
/// Both are approximations: real structures are not isotropic, and OIT's
/// weight function is not a clean
/// per-layer transmittance product -- the quality-bound test in
/// `tests/headless.rs` measures the resulting pixel error directly rather
/// than trusting this derivation. One extra layer of slack absorbs a
/// frame's worth of disocclusion drift while the camera moves -- this
/// reference has no phase-2 revalidation pass the way the opaque Hi-Z
/// does (`shaders/cull.wgsl`'s two-phase `PHASE`), so a margin this frame
/// finds too tight would stay wrong for the whole frame rather than being
/// caught before it is drawn. `0` when there is nothing to estimate a
/// spacing from.
pub fn glass_cull_margin(material: &Material, atom_count: usize, radius: f32) -> f32 {
    if !material.is_glass() || atom_count < 2 || radius <= 0.0 {
        return 0.0;
    }
    const LN_1_255: f32 = -5.541_263; // ln(1.0 / 255.0)
    let opacity = material.opacity.clamp(0.02, 0.98);
    let layers = (LN_1_255 / (1.0 - opacity).ln()).ceil().max(1.0) + 1.0;
    let volume = std::f32::consts::PI * (4.0 / 3.0) * radius.powi(3);
    let spacing = (volume / atom_count as f32).cbrt();
    layers * spacing
}

/// Up to [`PAGE_ATOMS`] atoms' geometry, shared by every rep drawing
/// them: positions and radii, bonds, cluster bounds. Frame changes
/// rewrite these; nothing per rep lives here.
pub struct Page {
    pub atom_count: u32,
    pub bond_count: u32,
    pub atoms: wgpu::Buffer,
    /// Page-relative atom index pairs; bonds crossing pages are dropped.
    pub bonds: wgpu::Buffer,
    /// For each entry of `bonds`, its index in the source `BondTable`, so
    /// a picked cylinder maps back to `pairs[i]` without a CPU-side copy
    /// of the (page-filtered) bond list.
    pub bond_ids: wgpu::Buffer,
    /// Per-cluster bounding spheres, rewritten with the atoms.
    pub cluster_bounds: wgpu::Buffer,
}

/// One rep's state for one page: its colours, per-item parameters, and
/// what the cull pass writes each frame. Apart from [`Page`] so several
/// reps can draw the same geometry in one frame, each with its own
/// radius, colours, material and survivors.
pub struct PageState {
    pub colors: wgpu::Buffer,
    pub visible_quads: wgpu::Buffer,
    pub visible_points: wgpu::Buffer,
    pub visible_bonds: wgpu::Buffer,
    pub indirect: wgpu::Buffer,
    /// Bits for the atoms, then the bonds, that occlusion phase 1
    /// rejected; cleared every frame.
    pub occluded: wgpu::Buffer,
    /// Per-cluster cull state, one `u32` per [`CLUSTER_ATOMS`] atoms
    /// (shaders/cull.wgsl), written every frame.
    pub clusters: wgpu::Buffer,
    /// A [`PageParams`] uniform.
    pub params: wgpu::Buffer,
}

/// A rep's draw state over all of a structure's pages
/// (`Renderer::bind_state`). A [`GpuStructure`] carries one of its own
/// (`state`); add more to draw the same atoms several ways at once.
pub struct DrawState {
    pub pages: Vec<PageState>,
}

impl DrawState {
    /// Fresh state for `structure`'s pages, one packed colour per atom.
    pub fn new(ctx: &GpuContext, structure: &GpuStructure, colors: &[u32]) -> Self {
        Self::try_new(ctx, structure, colors).expect("out of GPU memory")
    }

    pub fn try_new(
        ctx: &GpuContext,
        structure: &GpuStructure,
        colors: &[u32],
    ) -> Result<Self, OutOfGpuMemory> {
        assert_eq!(colors.len(), structure.atom_count);
        let pages = structure
            .pages
            .iter()
            .zip(colors.chunks(PAGE_ATOMS))
            .enumerate()
            .map(|(i, (page, chunk))| PageState::new(ctx, i, page, chunk))
            .collect::<Result<_, _>>()?;
        Ok(Self { pages })
    }

    /// Rewrites every page's colours in place, one packed colour per atom.
    pub fn set_colors(&self, ctx: &GpuContext, colors: &[u32]) {
        for (page, chunk) in self.pages.iter().zip(colors.chunks(PAGE_ATOMS)) {
            ctx.queue
                .write_buffer(&page.colors, 0, bytemuck::cast_slice(chunk));
        }
    }
}

pub struct GpuStructure {
    pub pages: Vec<Page>,
    /// The draw state `Renderer::bind` uses; other reps of the same atoms
    /// bring their own (`DrawState::new`).
    pub state: DrawState,
    pub atom_count: usize,
    /// Bonds actually uploaded (those crossing a page boundary are not).
    pub bond_count: usize,
    /// Length of the `BondTable` the bonds came from, which is the space
    /// `Page::bond_ids` indexes into.
    pub bond_table_len: usize,
    pub center: Vec3,
    pub radius: f32,
}

/// Survivor counts of the last cull pass for one page.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct VisibleCounts {
    pub quads: u32,
    pub points: u32,
    pub bonds: u32,
}

fn read_buffer(ctx: &GpuContext, src: &wgpu::Buffer, size: u64) -> Vec<u32> {
    let staging = ctx.device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("readback"),
        size,
        usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
        mapped_at_creation: false,
    });
    let mut encoder = ctx.device.create_command_encoder(&Default::default());
    encoder.copy_buffer_to_buffer(src, 0, &staging, 0, size);
    ctx.queue.submit([encoder.finish()]);
    let slice = staging.slice(..);
    let (tx, rx) = std::sync::mpsc::channel();
    slice.map_async(wgpu::MapMode::Read, move |r| {
        let _ = tx.send(r);
    });
    let _ = ctx.device.poll(wgpu::PollType::wait_indefinitely());
    rx.recv().expect("map callback").expect("map readback");
    let data = slice.get_mapped_range().expect("mapped range");
    let words: Vec<u32> = bytemuck::cast_slice(&data).to_vec();
    drop(data);
    staging.unmap();
    words
}

impl GpuStructure {
    /// Blocking readback of each page's indirect arguments. Diagnostics only.
    pub fn read_visible_counts(&self, ctx: &GpuContext) -> Vec<VisibleCounts> {
        self.state
            .pages
            .iter()
            .map(|page| {
                let w = read_buffer(ctx, &page.indirect, INDIRECT_RESET.len() as u64 * 4);
                VisibleCounts {
                    quads: (w[0] + w[16]) / 6,
                    points: w[4] + w[20],
                    // Cylinders and 1-px lines alike.
                    bonds: (w[8] + w[24]) / 6 + (w[12] + w[28]) / 2,
                }
            })
            .collect()
    }

    /// Blocking readback of every page's compacted index lists (quads,
    /// points), truncated to the survivor counts. Diagnostics only.
    pub fn read_visible_indices(&self, ctx: &GpuContext) -> Vec<(Vec<u32>, Vec<u32>)> {
        let counts = self.read_visible_counts(ctx);
        self.state
            .pages
            .iter()
            .zip(counts)
            .map(|(page, count)| {
                let read = |src: &wgpu::Buffer, n: u32| {
                    if n == 0 {
                        Vec::new()
                    } else {
                        read_buffer(ctx, src, n as u64 * 4)
                    }
                };
                (
                    read(&page.visible_quads, count.quads),
                    read(&page.visible_points, count.points),
                )
            })
            .collect()
    }

    /// Uploads frame 0 with van der Waals radii and `scheme` colors, plus
    /// `bonds` if given (needed for ball-and-stick).
    pub fn upload(
        ctx: &GpuContext,
        structure: &Structure,
        bonds: Option<&BondTable>,
        scheme: ColorScheme,
    ) -> Self {
        let colors = colors_for(scheme, &structure.topology);
        Self::upload_colored(ctx, structure, bonds, &colors)
    }

    /// `upload` with one packed color per atom already decided by the
    /// caller (a scheme, or a value column from outside).
    pub fn upload_colored(
        ctx: &GpuContext,
        structure: &Structure,
        bonds: Option<&BondTable>,
        colors: &[u32],
    ) -> Self {
        let coords = structure.frame(0);
        let positions = coords.positions();
        let radii = vdw_radii(&structure.topology.element);
        let pairs = bonds.map_or(&[][..], |b| b.pairs.as_slice());
        Self::from_parts(ctx, positions, &radii, colors, pairs)
    }

    /// Uploads arbitrary spheres (one position, radius, and packed color
    /// each) plus cylinders between `bonds` pairs. `upload` is this with
    /// van der Waals radii; the tube representation feeds it spline
    /// samples instead of atoms.
    pub fn from_parts(
        ctx: &GpuContext,
        positions: &[Vec3],
        radii: &[f32],
        colors: &[u32],
        bonds: &[[u32; 2]],
    ) -> Self {
        Self::try_from_parts(ctx, positions, radii, colors, bonds).expect("out of GPU memory")
    }

    pub fn try_from_parts(
        ctx: &GpuContext,
        positions: &[Vec3],
        radii: &[f32],
        colors: &[u32],
        bonds: &[[u32; 2]],
    ) -> Result<Self, OutOfGpuMemory> {
        assert_eq!(positions.len(), colors.len());
        let mut gpu = Self::try_geometry(ctx, positions, radii, bonds)?;
        gpu.state = DrawState::try_new(ctx, &gpu, colors)?;
        Ok(gpu)
    }

    /// Geometry only, with an empty draw state: for atoms drawn only
    /// through other states (`DrawState::new`, `Renderer::bind_state`),
    /// so no colour or cull buffers are allocated that nothing uses.
    pub fn geometry(
        ctx: &GpuContext,
        positions: &[Vec3],
        radii: &[f32],
        bonds: &[[u32; 2]],
    ) -> Self {
        Self::try_geometry(ctx, positions, radii, bonds).expect("out of GPU memory")
    }

    pub fn try_geometry(
        ctx: &GpuContext,
        positions: &[Vec3],
        radii: &[f32],
        bonds: &[[u32; 2]],
    ) -> Result<Self, OutOfGpuMemory> {
        assert_eq!(positions.len(), radii.len());
        let (center, radius) = CoordSet::new(positions.to_vec())
            .bounding_sphere()
            .unwrap_or((Vec3::ZERO, 1.0));

        let mut bond_count = 0;
        let pages: Vec<Page> = positions
            .chunks(PAGE_ATOMS)
            .zip(radii.chunks(PAGE_ATOMS))
            .enumerate()
            .map(|(i, (pos, rad))| {
                let base = (i * PAGE_ATOMS) as u32;
                let end = base + pos.len() as u32;
                let (page_bonds, page_bond_ids): (Vec<[u32; 2]>, Vec<u32>) = bonds
                    .iter()
                    .enumerate()
                    .filter(|(_, [a, c])| *a >= base && *a < end && *c >= base && *c < end)
                    .map(|(index, [a, c])| ([a - base, c - base], index as u32))
                    .unzip();
                bond_count += page_bonds.len();
                Page::new(ctx, i, pos, rad, &page_bonds, &page_bond_ids)
            })
            .collect::<Result<_, _>>()?;
        let state = DrawState { pages: Vec::new() };

        Ok(Self {
            pages,
            state,
            atom_count: positions.len(),
            bond_count,
            bond_table_len: bonds.len(),
            center,
            radius,
        })
    }

    /// Rewrites every page's positions and radii in place (colors and
    /// bonds are untouched). `Page.atoms` carries `COPY_DST` for this.
    /// Bounds (`center`/`radius`) stay as uploaded so the camera does not
    /// jump while scrubbing frames.
    pub fn set_instances(&self, ctx: &GpuContext, positions: &[Vec3], radii: &[f32]) {
        assert_eq!(positions.len(), self.atom_count);
        for (page, (pos, rad)) in self
            .pages
            .iter()
            .zip(positions.chunks(PAGE_ATOMS).zip(radii.chunks(PAGE_ATOMS)))
        {
            ctx.queue.write_buffer(
                &page.atoms,
                0,
                bytemuck::cast_slice(&instances_for(pos, rad)),
            );
            ctx.queue.write_buffer(
                &page.cluster_bounds,
                0,
                bytemuck::cast_slice(&cluster_bounds(pos, rad)),
            );
        }
    }

    /// `set_instances` from coordinate set `frame` of `structure` (the
    /// topology, and so the radii, is the same for every frame).
    pub fn set_frame(&self, ctx: &GpuContext, structure: &Structure, frame: usize) {
        let radii = vdw_radii(&structure.topology.element);
        self.set_instances(ctx, structure.frame(frame).positions(), &radii);
    }

    /// Rewrites the default state's colours in place, one packed colour
    /// per sphere.
    pub fn set_colors(&self, ctx: &GpuContext, colors: &[u32]) {
        assert_eq!(colors.len(), self.atom_count);
        self.state.set_colors(ctx, colors);
    }

    /// `set_colors` under `scheme` for a structure uploaded with `upload`.
    pub fn recolor(&self, ctx: &GpuContext, topology: &Topology, scheme: ColorScheme) {
        self.set_colors(ctx, &colors_for(scheme, topology));
    }
}

/// Van der Waals radius per atom, what `upload`/`set_frame` store.
pub fn vdw_radii(elements: &[vv_core::Element]) -> Vec<f32> {
    elements.par_iter().map(|e| e.vdw_radius()).collect()
}

/// A bounding sphere (xyz centre, w radius) per [`CLUSTER_ATOMS`] atoms,
/// at full van der Waals radius so it bounds every representation.
fn cluster_bounds(positions: &[Vec3], radii: &[f32]) -> Vec<[f32; 4]> {
    positions
        .par_chunks(CLUSTER_ATOMS as usize)
        .zip(radii.par_chunks(CLUSTER_ATOMS as usize))
        .map(|(chunk, radii)| {
            let (lo, hi) = chunk.iter().fold(
                (Vec3::splat(f32::MAX), Vec3::splat(f32::MIN)),
                |(lo, hi), &p| (lo.min(p), hi.max(p)),
            );
            let r = radii.iter().copied().fold(0.0f32, f32::max);
            let c = (lo + hi) * 0.5;
            [c.x, c.y, c.z, (hi - lo).length() * 0.5 + r]
        })
        .collect()
}

fn instances_for(positions: &[Vec3], radii: &[f32]) -> Vec<AtomInstance> {
    positions
        .par_iter()
        .zip(radii)
        .map(|(p, r)| AtomInstance {
            position: p.to_array(),
            radius: *r,
        })
        .collect()
}

/// One `vv_core::cartoon::CartoonMesh` vertex, 20 bytes: position,
/// normal, color. The vertex fetch unpacks the normal and color, so
/// `shaders/cartoon.wgsl` sees plain floats.
/// A cartoon's draw args before its LOD pass: 48 indices per instance
/// (`cartoon_join_indices`), no instances yet.
pub const CARTOON_INDIRECT_RESET: [u32; 5] = [vv_core::cartoon::RING as u32 * 6, 0, 0, 0, 0];

/// Index pattern for one cartoon join (instance): two rings of
/// `vv_core::cartoon::RING` vertices, two triangles per side, wound
/// counter-clockwise from outside (as `CartoonMesh::expand`).
pub fn cartoon_join_indices() -> Vec<u16> {
    let ring = vv_core::cartoon::RING as u16;
    (0..ring)
        .flat_map(|i| {
            let i2 = (i + 1) % ring;
            [i, i2, ring + i2, i, ring + i2, ring + i]
        })
        .collect()
}

/// One cartoon cross-section on the GPU (`shaders/cartoon.wgsl`, which
/// builds the ring vertices from it): 64 bytes, std430.
#[repr(C)]
#[derive(Clone, Copy, Debug, Pod, Zeroable)]
pub struct CartoonSectionGpu {
    pub center: [f32; 3],
    pub half_width: f32,
    pub across: [f32; 3],
    pub half_thickness: f32,
    pub up: [f32; 3],
    pub roundness: f32,
    /// Packed RGBA8, the same `u32` `colors_for` produces.
    pub color: u32,
    pub _pad: [u32; 3],
}

/// Per-cartoon uniform: the frame-wide pick id its first vertex gets
/// (`shaders/cartoon.wgsl` adds `@builtin(vertex_index)` to this, the
/// same "index + base, then +1 for the zero-is-nothing sentinel"
/// convention `atoms.wgsl`'s `PageParams` uses). Written fresh every
/// frame by `Renderer::render_all`, since a cartoon's place in the
/// frame-wide id space depends on how many atoms every `DrawItem` before
/// it claimed that frame.
#[repr(C)]
#[derive(Clone, Copy, Debug, Pod, Zeroable)]
pub struct CartoonParams {
    pub id_base: u32,
    pub _pad: [u32; 3],
    pub material: Material,
}

/// A cartoon on the GPU: its cross-sections and joins
/// (`vv_core::cartoon::CartoonPlan`); the vertex shader builds the ring
/// triangles, so there are no vertex or index buffers (8GLV: ~1 GB of
/// mesh became ~200 MB). `source` (one trace atom per section) maps a
/// pick back to an atom, as `DrawSource::atom` does for tubes.
pub struct CartoonGpu {
    pub sections: wgpu::Buffer,
    /// `[first, last]` section per residue step (`CartoonPlan::spans`).
    pub spans: wgpu::Buffer,
    /// Section pairs to draw this frame, written by the level-of-detail
    /// pass (`shaders/cartoon_lod.wgsl`); room for every join.
    pub pairs: wgpu::Buffer,
    /// `DrawIndexedIndirectArgs`; the LOD pass counts instances into it.
    pub indirect: wgpu::Buffer,
    /// Each section's place on the spline (`vv_core::cartoon::
    /// SectionRecipe`), and each residue slot's control point and guide:
    /// what `shaders/cartoon_frame.wgsl` rebuilds the sections from when
    /// the atoms move ([`CartoonGpu::set_frame`]).
    pub recipes: wgpu::Buffer,
    pub spline: wgpu::Buffer,
    /// A new frame's spline is uploaded but the sections not yet rebuilt.
    stale: AtomicBool,
    pub section_count: u32,
    pub span_count: u32,
    pub source: std::sync::Arc<Vec<u32>>,
    /// A [`CartoonParams`] uniform; rewritten every frame before the
    /// draw pass (see [`Renderer::render_all`]).
    pub params: wgpu::Buffer,
    /// A sphere around the mesh (framing, depth cue).
    pub bounds_center: Vec3,
    pub bounds_radius: f32,
}

/// Slot `s`'s control point at `2 s`, its guide at `2 s + 1`.
fn spline_data(frame: &vv_core::cartoon::CartoonFrame) -> Vec<[f32; 4]> {
    frame
        .controls
        .iter()
        .zip(&frame.guides)
        .flat_map(|(c, g)| [c.extend(0.0).to_array(), g.extend(0.0).to_array()])
        .collect()
}

impl CartoonGpu {
    /// `plan`'s sections at `frame` (its [`vv_core::cartoon::
    /// CartoonPlan::frame`]), each colored by `colors[section.source]`.
    pub fn upload(
        ctx: &GpuContext,
        plan: &vv_core::cartoon::CartoonPlan,
        frame: &vv_core::cartoon::CartoonFrame,
        colors: &[u32],
    ) -> Result<Self, OutOfGpuMemory> {
        let mesh = plan.mesh(frame);
        let (min, max) = mesh.sections.iter().fold(
            (Vec3::splat(f32::INFINITY), Vec3::splat(f32::NEG_INFINITY)),
            |(min, max), s| (min.min(s.center), max.max(s.center)),
        );
        let (bounds_center, bounds_radius) = if mesh.sections.is_empty() {
            (Vec3::ZERO, 1.0)
        } else {
            (
                (min + max) * 0.5,
                ((max - min).length() * 0.5 + 2.0).max(1.0),
            )
        };
        let sections: Vec<CartoonSectionGpu> = mesh
            .sections
            .iter()
            .map(|s| CartoonSectionGpu {
                center: s.center.to_array(),
                half_width: s.half_width,
                across: s.across.to_array(),
                half_thickness: s.half_thickness,
                up: s.up.to_array(),
                roundness: s.roundness,
                color: colors[s.source as usize],
                _pad: [0; 3],
            })
            .collect();
        let recipes: Vec<[u32; 4]> = plan
            .recipes
            .iter()
            .map(|r| [r.slot, r.t.to_bits(), r.first, r.end])
            .collect();
        let sections_buffer = try_buffer_init(
            ctx,
            "cartoon sections",
            bytemuck::cast_slice(&sections),
            wgpu::BufferUsages::STORAGE,
        )?;
        let spans = try_buffer_init(
            ctx,
            "cartoon spans",
            bytemuck::cast_slice(&mesh.spans),
            wgpu::BufferUsages::STORAGE,
        )?;
        let recipes = try_buffer_init(
            ctx,
            "cartoon recipes",
            bytemuck::cast_slice(&recipes),
            wgpu::BufferUsages::STORAGE,
        )?;
        let spline = try_buffer_init(
            ctx,
            "cartoon spline",
            bytemuck::cast_slice(&spline_data(frame)),
            wgpu::BufferUsages::STORAGE,
        )?;
        let pairs = catch_oom(&ctx.device, || {
            ctx.device.create_buffer(&wgpu::BufferDescriptor {
                label: Some("cartoon pairs"),
                size: (mesh.joins.len().max(1) * 8) as u64,
                usage: wgpu::BufferUsages::STORAGE,
                mapped_at_creation: false,
            })
        })?;
        let indirect = ctx
            .device
            .create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("cartoon indirect"),
                contents: bytemuck::cast_slice(&CARTOON_INDIRECT_RESET),
                usage: wgpu::BufferUsages::STORAGE
                    | wgpu::BufferUsages::INDIRECT
                    | wgpu::BufferUsages::COPY_DST,
            });
        let params = ctx
            .device
            .create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("cartoon params"),
                contents: bytemuck::bytes_of(&CartoonParams {
                    id_base: 0,
                    _pad: [0; 3],
                    material: Material::default(),
                }),
                usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            });
        Ok(Self {
            sections: sections_buffer,
            spans,
            pairs,
            indirect,
            recipes,
            spline,
            stale: AtomicBool::new(false),
            section_count: mesh.sections.len() as u32,
            span_count: mesh.spans.len() as u32,
            source: std::sync::Arc::new(mesh.sections.iter().map(|s| s.source).collect()),
            params,
            bounds_center,
            bounds_radius,
        })
    }

    /// Moves the cartoon to `frame` (of the plan it was uploaded with):
    /// uploads the spline; the renderer rebuilds the sections from it
    /// before the next draw.
    pub fn set_frame(&self, ctx: &GpuContext, frame: &vv_core::cartoon::CartoonFrame) {
        ctx.queue
            .write_buffer(&self.spline, 0, bytemuck::cast_slice(&spline_data(frame)));
        self.stale.store(true, Ordering::Relaxed);
    }

    /// Whether the sections need rebuilding from a new spline; clears it.
    pub(crate) fn take_stale(&self) -> bool {
        self.stale.swap(false, Ordering::Relaxed)
    }
}

/// One glycan-glyph mesh vertex (`shaders/glycan_mesh.wgsl`'s `Vertex`):
/// 32 bytes, std430. No atom index: the vertex shader stamps
/// `vertex_index` as the pick id, and [`GlycanGpu::source`] maps that
/// back to a real atom on the CPU side, the same way [`CartoonGpu::
/// source`] maps a section index.
#[repr(C)]
#[derive(Clone, Copy, Debug, Pod, Zeroable)]
pub struct GlycanVertexGpu {
    pub position: [f32; 3],
    pub _pad0: f32,
    pub normal: [f32; 3],
    pub color: u32,
}

/// A glycan rep's SNFG shapes on the GPU: a triangle soup rebuilt on the
/// CPU every frame (`vv_core::glycan::mesh::build_mesh`, a residue's
/// shape needs only its ring centroid and two reference points, not a
/// spline) and uploaded whole. Unlike [`CartoonGpu`], there is no
/// per-frame GPU rebuild to avoid: a few hundred residues' worth of
/// small glyphs is a fraction of a cartoon's mesh, so recreating the
/// buffer each frame is simpler and still well under budget.
pub struct GlycanGpu {
    pub vertices: wgpu::Buffer,
    pub vertex_count: u32,
    pub params: wgpu::Buffer,
    /// The real atom behind vertex `i` (its residue's anomeric carbon).
    pub source: std::sync::Arc<Vec<u32>>,
    pub bounds_center: Vec3,
    pub bounds_radius: f32,
}

impl GlycanGpu {
    /// `mesh`'s triangles (`vv_core::PolytopeMesh`), freshly built for
    /// the current frame.
    pub fn upload(ctx: &GpuContext, mesh: &vv_core::PolytopeMesh) -> Result<Self, OutOfGpuMemory> {
        let (min, max) = mesh.positions.iter().fold(
            (Vec3::splat(f32::INFINITY), Vec3::splat(f32::NEG_INFINITY)),
            |(min, max), &p| (min.min(p), max.max(p)),
        );
        let (bounds_center, bounds_radius) = if mesh.positions.is_empty() {
            (Vec3::ZERO, 1.0)
        } else {
            (
                (min + max) * 0.5,
                ((max - min).length() * 0.5 + 2.0).max(1.0),
            )
        };
        let verts: Vec<GlycanVertexGpu> = mesh
            .positions
            .iter()
            .zip(&mesh.normals)
            .zip(&mesh.colors)
            .map(|((&p, &n), &c)| GlycanVertexGpu {
                position: p.to_array(),
                _pad0: 0.0,
                normal: n.to_array(),
                color: color::rgba(c[0], c[1], c[2]),
            })
            .collect();
        let vertices = try_buffer_init(
            ctx,
            "glycan vertices",
            bytemuck::cast_slice(&verts),
            wgpu::BufferUsages::STORAGE,
        )?;
        let params = try_buffer_init(
            ctx,
            "glycan params",
            bytemuck::bytes_of(&CartoonParams {
                id_base: 0,
                _pad: [0; 3],
                material: Material::default(),
            }),
            wgpu::BufferUsages::UNIFORM,
        )?;
        Ok(Self {
            vertices,
            vertex_count: verts.len() as u32,
            params,
            source: std::sync::Arc::new(mesh.source_atom.clone()),
            bounds_center,
            bounds_radius,
        })
    }
}

/// The march's uniform (`shaders/gaussian_surface.wgsl`); layout matches
/// byte-for-byte (every `[T; 3]` is followed by a 4-byte scalar, which is
/// exactly how WGSL packs `vec3<T>` then a scalar). Only `view_inv`
/// changes per frame (`Renderer::render_all` rewrites it).
#[repr(C)]
#[derive(Clone, Copy, Debug, Pod, Zeroable)]
pub struct GaussianSurfaceParams {
    pub view_inv: [[f32; 4]; 4],
    /// World position of voxel (0, 0, 0).
    pub vol_origin: [f32; 3],
    /// Voxel spacing in Angstrom.
    pub voxel: f32,
    pub dims: [u32; 3],
    pub isovalue: f32,
    /// 8^3-voxel bricks per axis, for empty-space skipping.
    pub brick_dims: [u32; 3],
    pub blob_factor: f32,
    pub material: Material,
}

/// The bake's uniform (`Params` in `shaders/gaussian_density.wgsl`).
#[repr(C)]
#[derive(Clone, Copy, Debug, Pod, Zeroable)]
struct BakeParams {
    vol_origin: [f32; 3],
    voxel: f32,
    dims: [u32; 3],
    isovalue: f32,
    brick_dims: [u32; 3],
    blob_factor: f32,
    cutoff_ratio: f32,
    atom_count: u32,
    atom_bits: u32,
    _pad: u32,
}

/// The bake's compute pipelines, built once per device
/// (`GpuContext::gaussian_bake`).
pub(crate) struct GaussianBake {
    scatter: wgpu::ComputePipeline,
    resolve: wgpu::ComputePipeline,
    layout: wgpu::BindGroupLayout,
}

impl GaussianBake {
    fn new(device: &wgpu::Device) -> Self {
        let buffer = |binding, ty| wgpu::BindGroupLayoutEntry {
            binding,
            visibility: wgpu::ShaderStages::COMPUTE,
            ty: wgpu::BindingType::Buffer {
                ty,
                has_dynamic_offset: false,
                min_binding_size: None,
            },
            count: None,
        };
        let read = wgpu::BufferBindingType::Storage { read_only: true };
        let write = wgpu::BufferBindingType::Storage { read_only: false };
        let layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("gaussian bake"),
            entries: &[
                buffer(0, read),
                buffer(1, read),
                buffer(2, wgpu::BufferBindingType::Uniform),
                buffer(3, write),
                buffer(4, write),
                wgpu::BindGroupLayoutEntry {
                    binding: 5,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::StorageTexture {
                        access: wgpu::StorageTextureAccess::WriteOnly,
                        format: wgpu::TextureFormat::Rgba16Float,
                        view_dimension: wgpu::TextureViewDimension::D3,
                    },
                    count: None,
                },
                buffer(6, write),
            ],
        });
        let module = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("gaussian bake"),
            source: wgpu::ShaderSource::Wgsl(include_str!("shaders/gaussian_density.wgsl").into()),
        });
        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("gaussian bake"),
            bind_group_layouts: &[Some(&layout)],
            immediate_size: 0,
        });
        let pipeline = |entry: &str| {
            device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
                label: Some("gaussian bake"),
                layout: Some(&pipeline_layout),
                module: &module,
                entry_point: Some(entry),
                compilation_options: Default::default(),
                cache: None,
            })
        };
        Self {
            scatter: pipeline("scatter"),
            resolve: pipeline("resolve"),
            layout,
        }
    }
}

/// Voxel budget for one Gaussian surface's density volume (8 bytes per
/// voxel: RGBA16F). 32M voxels = 256 MB; 8GLV (4M atoms) fits at ~1 A.
pub const GAUSSIAN_MAX_VOXELS: f64 = 32.0e6;
/// Finest voxel spacing, used whenever the budget allows (small
/// structures). Well under an atom radius, so the trilinear surface
/// stays close to the exact one `vv_cpu::gaussian_surface` renders.
pub const GAUSSIAN_MIN_VOXEL: f32 = 0.35;

/// The Gaussian bake's bindings: atoms, colours, params, density, colour
/// bids (`shaders/gaussian_density.wgsl`, bindings 0-4), the volume and
/// the brick flags.
fn bake_bind_group(
    device: &wgpu::Device,
    bake: &GaussianBake,
    buffers: [&wgpu::Buffer; 5],
    volume: &wgpu::TextureView,
    bricks: &wgpu::Buffer,
) -> wgpu::BindGroup {
    let mut entries: Vec<wgpu::BindGroupEntry> = buffers
        .iter()
        .enumerate()
        .map(|(k, b)| wgpu::BindGroupEntry {
            binding: k as u32,
            resource: b.as_entire_binding(),
        })
        .collect();
    entries.push(wgpu::BindGroupEntry {
        binding: 5,
        resource: wgpu::BindingResource::TextureView(volume),
    });
    entries.push(wgpu::BindGroupEntry {
        binding: 6,
        resource: bricks.as_entire_binding(),
    });
    device.create_bind_group(&wgpu::BindGroupDescriptor {
        label: Some("gaussian bake"),
        layout: &bake.layout,
        entries: &entries,
    })
}

/// Voxels of a Gaussian surface re-baked every frame of a playing
/// trajectory: an eighth of [`GAUSSIAN_MAX_VOXELS`], voxels twice as wide,
/// an eighth of the bake (the full volume comes back when playback stops).
pub const GAUSSIAN_PLAYBACK_VOXELS: f64 = GAUSSIAN_MAX_VOXELS / 8.0;

/// Retries `upload` with half the voxel budget while the GPU is out of
/// memory, from `start` down to 1/8 of [`GAUSSIAN_MAX_VOXELS`].
pub(crate) fn with_voxel_budget<T>(
    start: f64,
    mut upload: impl FnMut(f64) -> Result<T, OutOfGpuMemory>,
) -> Result<T, OutOfGpuMemory> {
    let mut budget = start;
    loop {
        match upload(budget) {
            Err(OutOfGpuMemory) if budget > GAUSSIAN_MAX_VOXELS / 8.0 => budget /= 2.0,
            result => return result,
        }
    }
}

/// Where a baked volume sits: voxel (0, 0, 0) at `min`, `dims` voxels of
/// `voxel` Angstrom, and its 8^3 bricks.
#[derive(Clone, Copy, Debug)]
pub(crate) struct VolumeBox {
    pub min: glam::Vec3,
    pub max: glam::Vec3,
    pub voxel: f32,
    pub dims: glam::UVec3,
    pub brick_dims: glam::UVec3,
}

impl VolumeBox {
    /// The box around `positions` grown by `margin`, in the finest voxels
    /// within `max_voxels` and the device's 3D texture size.
    pub fn around(
        ctx: &GpuContext,
        positions: &[glam::Vec3],
        margin: f32,
        max_voxels: f64,
    ) -> Self {
        let (min, max) = positions.iter().fold(
            (
                glam::Vec3::splat(f32::INFINITY),
                glam::Vec3::splat(f32::NEG_INFINITY),
            ),
            |(min, max), &p| (min.min(p), max.max(p)),
        );
        let (min, max) = (
            min - glam::Vec3::splat(margin),
            max + glam::Vec3::splat(margin),
        );
        let extent = max - min;
        let max_dim = ctx.device.limits().max_texture_dimension_3d.min(2048) as f32;
        let budget_voxel =
            ((extent.x as f64 * extent.y as f64 * extent.z as f64) / max_voxels).cbrt() as f32;
        let voxel = budget_voxel
            .max(GAUSSIAN_MIN_VOXEL)
            .max(extent.max_element() / (max_dim - 1.0));
        let dims = (extent / voxel).ceil().as_uvec3() + 1;
        Self {
            min,
            max,
            voxel,
            dims,
            brick_dims: ((dims - 1).max(glam::UVec3::ONE) + 7) / 8,
        }
    }

    pub fn center(&self) -> glam::Vec3 {
        (self.min + self.max) * 0.5
    }

    pub fn radius(&self) -> f32 {
        ((self.max - self.min).length() * 0.5).max(1.0)
    }

    /// The RGBA16F volume a bake writes and the march samples, its view,
    /// the brick flags (zeroed: every brick starts empty), and the sampler.
    pub fn allocate(
        &self,
        ctx: &GpuContext,
    ) -> Result<
        (
            wgpu::Texture,
            wgpu::TextureView,
            wgpu::Buffer,
            wgpu::Sampler,
        ),
        OutOfGpuMemory,
    > {
        let device = &ctx.device;
        let volume = catch_oom(device, || {
            device.create_texture(&wgpu::TextureDescriptor {
                label: Some("surface volume"),
                size: wgpu::Extent3d {
                    width: self.dims.x,
                    height: self.dims.y,
                    depth_or_array_layers: self.dims.z,
                },
                mip_level_count: 1,
                sample_count: 1,
                dimension: wgpu::TextureDimension::D3,
                format: wgpu::TextureFormat::Rgba16Float,
                usage: wgpu::TextureUsages::STORAGE_BINDING | wgpu::TextureUsages::TEXTURE_BINDING,
                view_formats: &[],
            })
        })?;
        let view = volume.create_view(&Default::default());
        let b = self.brick_dims;
        let bricks = catch_oom(device, || {
            device.create_buffer(&wgpu::BufferDescriptor {
                label: Some("surface volume bricks"),
                size: (b.x as u64 * b.y as u64 * b.z as u64 * 4).max(4),
                usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
                mapped_at_creation: false,
            })
        })?;
        let sampler = device.create_sampler(&wgpu::SamplerDescriptor {
            label: Some("surface volume"),
            address_mode_u: wgpu::AddressMode::ClampToEdge,
            address_mode_v: wgpu::AddressMode::ClampToEdge,
            address_mode_w: wgpu::AddressMode::ClampToEdge,
            mag_filter: wgpu::FilterMode::Linear,
            min_filter: wgpu::FilterMode::Linear,
            ..Default::default()
        });
        Ok((volume, view, bricks, sampler))
    }
}

/// One Gaussian surface on the GPU: its density field baked into a 3D
/// texture (rgb = color of the dominant atom, a = density) by a compute
/// pass at upload, plus an occupancy flag per 8^3 brick. Drawn by a
/// full-screen ray march (`shaders/gaussian_surface.wgsl`) that costs one
/// texture fetch per step and skips empty bricks, so frame time follows
/// what's on screen, not atom count. Same field as `vv_core::
/// gaussian_surface` / `vv_cpu::gaussian_surface`, sampled at `voxel`.
pub struct GaussianSurfaceGpu {
    pub volume: wgpu::Texture,
    pub volume_view: wgpu::TextureView,
    pub sampler: wgpu::Sampler,
    pub bricks: wgpu::Buffer,
    pub params: wgpu::Buffer,
    /// Everything but `view_inv`, which `Renderer::render_all` fills in.
    pub base_params: GaussianSurfaceParams,
    pub atom_count: u32,
    /// A sphere around the surface, for framing.
    pub bounds_center: glam::Vec3,
    pub bounds_radius: f32,
    /// What the bake reads, kept to bake a new frame in place
    /// ([`GaussianSurfaceGpu::set_frame`]); `None` for a baked SES.
    pub(crate) inputs: Option<BakeInputs>,
}

/// A Gaussian surface's bake inputs on the GPU, and what they were built
/// from.
pub(crate) struct BakeInputs {
    atoms: wgpu::Buffer,
    /// Per voxel, filled by the scatter pass: density (fixed point) and
    /// the colour bid.
    density: wgpu::Buffer,
    owner: wgpu::Buffer,
    bind: wgpu::BindGroup,
    radii: Vec<f32>,
    /// The largest atom's cutoff.
    reach: f32,
    /// The volume's box: a frame whose atoms' cutoff spheres leave it
    /// needs a new volume.
    min: glam::Vec3,
    max: glam::Vec3,
}

/// Room around the atoms' cutoff spheres, so a trajectory's frames bake
/// into the same volume as the first.
const TRAJECTORY_ROOM: f32 = 4.0;

/// Contributions below this are left out (`vv_core::gaussian_surface::
/// cutoff_radius`).
pub const GAUSSIAN_EPSILON: f32 = 0.01;

fn atom_instances(positions: &[glam::Vec3], radii: &[f32]) -> Vec<AtomInstance> {
    use rayon::prelude::*;
    positions
        .par_iter()
        .zip(radii)
        .map(|(&p, &r)| AtomInstance {
            position: p.to_array(),
            radius: r,
        })
        .collect()
}

impl GaussianSurfaceGpu {
    /// `positions`, `radii`, and `colors` (packed as `colors_for` does)
    /// must have the same length, one entry per atom. `blob_factor` is
    /// `vv_core::gaussian_surface::DEFAULT_BLOB_FACTOR` unless tuning it
    /// is worth exposing. Records and submits the bake; it runs on the
    /// GPU before the first frame that draws this surface. Out of GPU
    /// memory, it retries with coarser voxels (down to 1/8 of the budget)
    /// before giving up.
    pub fn upload(
        ctx: &GpuContext,
        positions: &[glam::Vec3],
        radii: &[f32],
        colors: &[u32],
        blob_factor: f32,
    ) -> Result<Self, OutOfGpuMemory> {
        Self::upload_within(
            ctx,
            positions,
            radii,
            colors,
            blob_factor,
            GAUSSIAN_MAX_VOXELS,
        )
    }

    /// [`GaussianSurfaceGpu::upload`] into at most `max_voxels` voxels
    /// (e.g. [`GAUSSIAN_PLAYBACK_VOXELS`]).
    pub fn upload_within(
        ctx: &GpuContext,
        positions: &[glam::Vec3],
        radii: &[f32],
        colors: &[u32],
        blob_factor: f32,
        max_voxels: f64,
    ) -> Result<Self, OutOfGpuMemory> {
        with_voxel_budget(max_voxels, |budget| {
            Self::upload_with_budget(ctx, positions, radii, colors, blob_factor, budget)
        })
    }

    fn upload_with_budget(
        ctx: &GpuContext,
        positions: &[glam::Vec3],
        radii: &[f32],
        colors: &[u32],
        blob_factor: f32,
        max_voxels: f64,
    ) -> Result<Self, OutOfGpuMemory> {
        assert_eq!(positions.len(), radii.len());
        assert_eq!(positions.len(), colors.len());
        assert!(!positions.is_empty(), "need at least one atom");
        let device = &ctx.device;

        // Volume box: every atom's cutoff sphere.
        let reach = vv_core::gaussian_surface::grid_cell_size(radii, blob_factor, GAUSSIAN_EPSILON);
        let grid_box = VolumeBox::around(ctx, positions, reach + TRAJECTORY_ROOM, max_voxels);
        let VolumeBox {
            min,
            voxel,
            dims,
            brick_dims,
            ..
        } = grid_box;
        let base_params = GaussianSurfaceParams {
            view_inv: glam::Mat4::IDENTITY.to_cols_array_2d(),
            vol_origin: min.to_array(),
            voxel,
            dims: dims.to_array(),
            // `vv_core::gaussian_surface`'s contribution is exactly 1.0
            // at a lone atom's own surface, for any `blob_factor`.
            isovalue: 1.0,
            brick_dims: brick_dims.to_array(),
            blob_factor,
            material: Material::default(),
        };
        let atom_count = positions.len() as u32;
        let bake_params = BakeParams {
            vol_origin: min.to_array(),
            voxel,
            dims: dims.to_array(),
            isovalue: 1.0,
            brick_dims: brick_dims.to_array(),
            blob_factor,
            cutoff_ratio: vv_core::gaussian_surface::cutoff_radius(
                1.0,
                blob_factor,
                GAUSSIAN_EPSILON,
            ),
            atom_count,
            // Enough bits for every atom index; the rest rank colour bids.
            atom_bits: (32 - (atom_count.max(2) - 1).leading_zeros()).max(20),
            _pad: 0,
        };

        let atoms = try_buffer_init(
            ctx,
            "gaussian surface atoms",
            bytemuck::cast_slice(&atom_instances(positions, radii)),
            wgpu::BufferUsages::STORAGE,
        )?;
        let colors_buffer = try_buffer_init(
            ctx,
            "gaussian surface colors",
            bytemuck::cast_slice(colors),
            wgpu::BufferUsages::STORAGE,
        )?;
        let voxels = dims.x as u64 * dims.y as u64 * dims.z as u64;
        let per_voxel = |label| {
            try_buffer(
                ctx,
                label,
                voxels * 4,
                wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
            )
        };
        let density = per_voxel("gaussian surface density")?;
        let owner = per_voxel("gaussian surface colour bids")?;
        let bake_uniform = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("gaussian bake params"),
            contents: bytemuck::bytes_of(&bake_params),
            usage: wgpu::BufferUsages::UNIFORM,
        });
        let params = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("gaussian surface params"),
            contents: bytemuck::bytes_of(&base_params),
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
        });
        let (volume, volume_view, bricks, sampler) = grid_box.allocate(ctx)?;

        let bake = ctx.gaussian_bake.get_or_init(|| GaussianBake::new(device));
        let bind = bake_bind_group(
            device,
            bake,
            [&atoms, &colors_buffer, &bake_uniform, &density, &owner],
            &volume_view,
            &bricks,
        );
        let surface = Self {
            volume,
            volume_view,
            sampler,
            bricks,
            params,
            base_params,
            atom_count,
            bounds_center: grid_box.center(),
            bounds_radius: grid_box.radius(),
            inputs: Some(BakeInputs {
                atoms,
                density,
                owner,
                bind,
                radii: radii.to_vec(),
                reach,
                min: grid_box.min,
                max: grid_box.max,
            }),
        };
        surface.bake(ctx);
        Ok(surface)
    }

    /// Records and submits the bake: clear, scatter, resolve.
    fn bake(&self, ctx: &GpuContext) {
        let Some(inputs) = &self.inputs else {
            return;
        };
        let bake = ctx
            .gaussian_bake
            .get_or_init(|| GaussianBake::new(&ctx.device));
        let mut encoder = ctx
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                label: Some("gaussian bake"),
            });
        for buffer in [&inputs.density, &inputs.owner, &self.bricks] {
            encoder.clear_buffer(buffer, 0, None);
        }
        {
            let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: Some("gaussian bake"),
                timestamp_writes: None,
            });
            pass.set_bind_group(0, &inputs.bind, &[]);
            pass.set_pipeline(&bake.scatter);
            pass.dispatch_workgroups(self.atom_count.div_ceil(64), 1, 1);
            pass.set_pipeline(&bake.resolve);
            let [x, y, z] = self.base_params.dims;
            pass.dispatch_workgroups(x.div_ceil(4), y.div_ceil(4), z.div_ceil(4));
        }
        ctx.queue.submit([encoder.finish()]);
    }

    /// Bakes the same atoms, same colours, at new `positions` (a
    /// trajectory frame) into this volume, reusing every buffer and the
    /// texture. `false` when it cannot -- a different atom count, or atoms
    /// outside the volume's box -- and the caller uploads afresh.
    pub fn set_frame(&mut self, ctx: &GpuContext, positions: &[glam::Vec3]) -> bool {
        let Some(inputs) = &self.inputs else {
            return false;
        };
        if positions.len() != self.atom_count as usize {
            return false;
        }
        use rayon::prelude::*;
        let reach = glam::Vec3::splat(inputs.reach);
        let inside = positions
            .par_iter()
            .all(|&p| (p - reach).cmpge(inputs.min).all() && (p + reach).cmple(inputs.max).all());
        if !inside {
            return false;
        }
        ctx.queue.write_buffer(
            &inputs.atoms,
            0,
            bytemuck::cast_slice(&atom_instances(positions, &inputs.radii)),
        );
        self.bake(ctx);
        true
    }

    /// Voxels in the baked volume (diagnostics).
    pub fn voxel_count(&self) -> u64 {
        let [x, y, z] = self.base_params.dims;
        x as u64 * y as u64 * z as u64
    }
}

/// Uniform for `shaders/skin_surface.wgsl`: `view_inv` and `id_base` are
/// rewritten every frame by `Renderer::render_all`, the rest stays fixed
/// once uploaded.
#[repr(C)]
#[derive(Clone, Copy, Debug, Pod, Zeroable)]
pub struct SkinSurfaceParams {
    pub view_inv: [[f32; 4]; 4],
    pub shrink: f32,
    pub atom_count: u32,
    pub patch_count: u32,
    /// Frame-wide pick id of this structure's atom 0 (a skin-surface hit
    /// is picked as the atom whose patch it lies on, `id_base + atom +
    /// 1`, the same convention `PageParams::atom_id_base` follows).
    pub id_base: u32,
    pub material: Material,
}

/// One `vv_core::skin_surface::Patch` as `shaders/skin_patch.wgsl`
/// reads it: flat scalars (no `vec3<f32>`, same reason as
/// [`SkinSurfaceParams`]), 64 bytes, one per mixed cell that has any
/// surface in it. `s_weight` is `shrink * w_z`, the constant term of the
/// quadric, folded in once here instead of per fragment.
#[repr(C)]
#[derive(Clone, Copy, Debug, Pod, Zeroable)]
pub struct SkinPatchGpu {
    pub center: [f32; 3],
    pub s_weight: f32,
    pub axis: [f32; 3],
    /// `vv_core::skin_surface::PatchKind` as `0..=3`.
    pub kind: u32,
    pub bound_center: [f32; 3],
    pub bound_radius: f32,
    /// Unused slots are `u32::MAX`.
    pub atoms: [u32; 4],
}

impl SkinPatchGpu {
    pub(crate) fn from_patch(p: &vv_core::skin_surface::Patch, shrink: f32) -> Self {
        Self {
            center: p.center.to_array(),
            s_weight: shrink * p.weight,
            axis: p.axis.to_array(),
            kind: p.kind as u32,
            bound_center: p.bound_center.to_array(),
            bound_radius: p.bound_radius,
            atoms: p.atoms,
        }
    }
}

/// Atoms, colors, the mixed complex's patches, and the ownership grid
/// for one skin surface on the GPU (`shaders/skin_surface.wgsl`) -- the
/// WGSL twin of `vv_cpu::skin_surface::SkinSurfaceScene`, same geometry
/// and membership rule (`vv_core::skin_surface`, read its module doc
/// first).
///
/// Drawn as one screen-aligned billboard per patch (six vertices, like
/// `draw.wgsl`'s sphere impostors, sized from the patch's bounding
/// sphere), each fragment ray-casting only its own patch and letting the
/// depth buffer pick the nearest accepted hit across patches -- the same
/// structure as every other impostor representation here. The previous
/// version was a full-screen pass that looped over every atom and every
/// edge per pixel, O(pixels x patches), which is what made the skin
/// surface unusably slow at any real zoom or atom count. Only patches
/// with a surface are uploaded (`Patch::has_surface`, the build-time
/// visibility cull), each with its competitor list
/// (`SkinComplex::competitors`) for the membership test.
pub struct SkinSurfaceGpu {
    pub atoms: wgpu::Buffer,
    pub colors: wgpu::Buffer,
    pub patches: wgpu::Buffer,
    /// CSR over the uploaded patches, in `patches` order.
    pub competitor_starts: wgpu::Buffer,
    pub competitors: wgpu::Buffer,
    pub params: wgpu::Buffer,
    pub atom_count: u32,
    pub patch_count: u32,
    /// Patches by kind (vertex, edge, triangle, tetrahedron), diagnostics.
    pub patch_counts: [usize; 4],
    /// A sphere around every patch, for framing.
    pub bounds_center: glam::Vec3,
    pub bounds_radius: f32,
    /// Cached so `Renderer::render_all` can rewrite the `params` uniform
    /// every frame (fresh `view_inv`/`id_base`, same shrink).
    pub shrink: f32,
}

impl SkinSurfaceGpu {
    /// `positions`, `weights` (`vv_core::skin_surface::weight_for_radius`
    /// of each atom's radius), and `colors` (packed as `colors_for` does)
    /// must have the same length, one entry per atom. Builds the mixed
    /// complex itself (`vv_core::skin_surface::build_complex`).
    pub fn upload(
        ctx: &GpuContext,
        positions: &[glam::Vec3],
        weights: &[f32],
        colors: &[u32],
        shrink: f32,
    ) -> Self {
        assert_eq!(positions.len(), weights.len());
        assert_eq!(positions.len(), colors.len());
        assert!(!positions.is_empty(), "need at least one atom");
        let complex = vv_core::skin_surface::build_complex(positions, weights, shrink);
        Self::from_complex(ctx, &complex, positions, weights, colors)
    }

    /// Uploads an already-built complex (the app builds it on a worker
    /// thread: it takes seconds at 100K+ atoms). Same argument rules as
    /// [`Self::upload`].
    pub fn from_complex(
        ctx: &GpuContext,
        complex: &vv_core::skin_surface::SkinComplex,
        positions: &[glam::Vec3],
        weights: &[f32],
        colors: &[u32],
    ) -> Self {
        assert_eq!(positions.len(), weights.len());
        assert_eq!(positions.len(), colors.len());
        let shrink = complex.shrink;
        let atom_data: Vec<[f32; 4]> = positions
            .iter()
            .zip(weights)
            .map(|(&p, &w)| [p.x, p.y, p.z, w])
            .collect();
        let mut patch_counts = [0usize; 4];
        let mut patch_data: Vec<SkinPatchGpu> = Vec::new();
        let mut competitor_starts: Vec<u32> = vec![0];
        let mut competitor_data: Vec<u32> = Vec::new();
        for (i, p) in complex.patches.iter().enumerate() {
            if !p.has_surface() {
                continue;
            }
            patch_counts[p.kind as usize] += 1;
            patch_data.push(SkinPatchGpu::from_patch(p, shrink));
            competitor_data.extend_from_slice(complex.competitors(i));
            competitor_starts.push(competitor_data.len() as u32);
        }
        let (min, max) = patch_data.iter().fold(
            (
                glam::Vec3::splat(f32::INFINITY),
                glam::Vec3::splat(f32::NEG_INFINITY),
            ),
            |(min, max), p| {
                let c = glam::Vec3::from(p.bound_center);
                (min.min(c), max.max(c))
            },
        );
        let bounds_center = (min + max) * 0.5;
        let bounds_radius = patch_data
            .iter()
            .map(|p| (glam::Vec3::from(p.bound_center) - bounds_center).length() + p.bound_radius)
            .fold(0.0f32, f32::max)
            .max(1.0);

        let device = &ctx.device;
        let atoms = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("skin surface atoms"),
            contents: bytemuck::cast_slice(&atom_data),
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
        });
        let colors_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("skin surface colors"),
            contents: bytemuck::cast_slice(colors),
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
        });
        let competitor_starts_buffer =
            device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("skin surface competitor starts"),
                contents: bytemuck::cast_slice(&competitor_starts),
                usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
            });
        let competitors_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("skin surface competitors"),
            contents: bytemuck::cast_slice(if competitor_data.is_empty() {
                &[0u32][..]
            } else {
                &competitor_data
            }),
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
        });
        // Storage bindings can't be empty; a zeroed placeholder keeps a
        // patchless structure on the same layout (never drawn: see
        // `Renderer::render_all`'s `patch_count == 0` skip).
        let placeholder = SkinPatchGpu::zeroed();
        let patches = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("skin surface patches"),
            contents: if patch_data.is_empty() {
                bytemuck::bytes_of(&placeholder)
            } else {
                bytemuck::cast_slice(&patch_data)
            },
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
        });
        let params = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("skin surface params"),
            contents: bytemuck::bytes_of(&SkinSurfaceParams {
                view_inv: glam::Mat4::IDENTITY.to_cols_array_2d(),
                shrink,
                atom_count: positions.len() as u32,
                patch_count: patch_data.len() as u32,
                id_base: 0,
                material: Material::default(),
            }),
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
        });
        Self {
            atoms,
            colors: colors_buffer,
            patches,
            competitor_starts: competitor_starts_buffer,
            competitors: competitors_buffer,
            params,
            atom_count: positions.len() as u32,
            patch_count: patch_data.len() as u32,
            patch_counts,
            bounds_center,
            bounds_radius,
            shrink,
        }
    }

    /// Rewrites the colours in place, one packed colour per atom (the
    /// shader looks a patch's atom up in this buffer, unlike SES's
    /// per-patch colours): a recolor needs no other buffer to change.
    pub fn set_colors(&self, ctx: &GpuContext, colors: &[u32]) {
        assert_eq!(colors.len(), self.atom_count as usize);
        ctx.queue
            .write_buffer(&self.colors, 0, bytemuck::cast_slice(colors));
    }
}

impl Page {
    fn new(
        ctx: &GpuContext,
        index: usize,
        positions: &[Vec3],
        radii: &[f32],
        bonds: &[[u32; 2]],
        bond_ids: &[u32],
    ) -> Result<Self, OutOfGpuMemory> {
        let instances = instances_for(positions, radii);
        let label = |what: &str| format!("page {index} {what}");
        // Storage bindings can't be empty; one-entry placeholders keep
        // bondless pages on the same layout as bonded ones.
        let bond_ids: &[u32] = if bond_ids.is_empty() { &[0] } else { bond_ids };
        let bonds_or_placeholder: &[[u32; 2]] = if bonds.is_empty() { &[[0, 0]] } else { bonds };
        Ok(Self {
            atom_count: instances.len() as u32,
            bond_count: bonds.len() as u32,
            atoms: try_buffer_init(
                ctx,
                &label("atoms"),
                bytemuck::cast_slice(&instances),
                wgpu::BufferUsages::STORAGE,
            )?,
            bonds: try_buffer_init(
                ctx,
                &label("bonds"),
                bytemuck::cast_slice(bonds_or_placeholder),
                wgpu::BufferUsages::STORAGE,
            )?,
            bond_ids: try_buffer_init(
                ctx,
                &label("bond ids"),
                bytemuck::cast_slice(bond_ids),
                wgpu::BufferUsages::STORAGE,
            )?,
            cluster_bounds: try_buffer_init(
                ctx,
                &label("cluster bounds"),
                bytemuck::cast_slice(&cluster_bounds(positions, radii)),
                wgpu::BufferUsages::STORAGE,
            )?,
        })
    }
}

/// An uninitialised buffer of `size` bytes, reporting out-of-memory.
fn try_buffer(
    ctx: &GpuContext,
    label: &str,
    size: u64,
    usage: wgpu::BufferUsages,
) -> Result<wgpu::Buffer, OutOfGpuMemory> {
    catch_oom(&ctx.device, || {
        ctx.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some(label),
            size,
            usage,
            mapped_at_creation: false,
        })
    })
}

impl PageState {
    fn new(
        ctx: &GpuContext,
        index: usize,
        page: &Page,
        colors: &[u32],
    ) -> Result<Self, OutOfGpuMemory> {
        let (n, bonds) = (page.atom_count, page.bond_count);
        let label = |what: &str| format!("page {index} {what}");
        let index_buffer = |what: &str, count: u32| {
            try_buffer(
                ctx,
                &label(what),
                (count.max(1) as u64) * 4,
                wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
            )
        };
        Ok(Self {
            params: try_buffer_init(
                ctx,
                &label("params"),
                bytemuck::bytes_of(&PageParams::default()),
                wgpu::BufferUsages::UNIFORM,
            )?,
            colors: try_buffer_init(
                ctx,
                &label("colors"),
                bytemuck::cast_slice(colors),
                wgpu::BufferUsages::STORAGE,
            )?,
            visible_quads: index_buffer("visible quads", n)?,
            visible_points: index_buffer("visible points", n)?,
            visible_bonds: index_buffer("visible bonds", bonds)?,
            indirect: try_buffer_init(
                ctx,
                &label("indirect"),
                bytemuck::cast_slice(&INDIRECT_RESET),
                wgpu::BufferUsages::STORAGE
                    | wgpu::BufferUsages::INDIRECT
                    | wgpu::BufferUsages::COPY_SRC,
            )?,
            occluded: try_buffer(
                ctx,
                &label("occluded"),
                ((n as usize).div_ceil(32) + (bonds as usize).div_ceil(32)).max(1) as u64 * 4,
                wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
            )?,
            clusters: try_buffer(
                ctx,
                &label("clusters"),
                (n.div_ceil(CLUSTER_ATOMS).max(1) * 4) as u64,
                wgpu::BufferUsages::STORAGE,
            )?,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn small_cif_path(name: &str) -> std::path::PathBuf {
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../fixtures/small")
            .join(name)
    }

    #[test]
    fn segname_tells_apart_segments_that_share_a_chain_letter() {
        let s = vv_io::load(small_cif_path("md_segments.pdb")).unwrap();
        let t = &s.topology;
        let colors = colors_for(ColorScheme::SegmentName, t);
        let of_segid = |segid: &str| {
            let chain = (0..t.chain_count()).find(|&c| t.segid(c) == segid).unwrap();
            colors[t.residues[t.chains[chain].residues.start as usize]
                .atoms
                .start as usize]
        };
        assert_ne!(of_segid("PROA"), of_segid("PROB"));
        let plain = vv_io::load(small_cif_path("1CRN.pdb")).unwrap();
        let by_chain = colors_for(ColorScheme::SegmentName, &plain.topology);
        assert_eq!(by_chain[0], color::by_chain(0), "no segids: by auth chain");
    }

    /// 1UBQ's chain A holds 76 residues and then 58 waters: the ramp
    /// must span the protein alone (blue first residue, red last) and
    /// leave the waters in element colors.
    #[test]
    fn rainbow_spans_the_polymer_and_skips_water() {
        let s = vv_io::load(small_cif_path("1UBQ.cif")).unwrap();
        let t = &s.topology;
        let colors = colors_for(ColorScheme::Rainbow, t);
        let ca = |r: usize| {
            t.residues[r]
                .atoms
                .clone()
                .find(|&a| t.atom_name(a as usize) == "CA")
                .unwrap() as usize
        };
        let rgb = |c: u32| c.to_le_bytes();
        let first = rgb(colors[ca(0)]);
        let last = rgb(colors[ca(75)]);
        assert!(
            first[2] > first[0],
            "first residue should be blue: {first:?}"
        );
        assert!(last[0] > last[2], "last residue should be red: {last:?}");
        let water = (0..t.atom_count())
            .find(|&a| t.residue_name(t.residue_index[a] as usize) == "HOH")
            .unwrap();
        assert_eq!(colors[water], color::by_element(t.element[water]));
    }

    #[test]
    fn scalar_colors_match_b_factor_colors_and_skip_nans() {
        let values = [1.0, f32::NAN, 3.0, 2.0];
        assert_eq!(scalar_range(&values), (1.0, 3.0));
        let colors = colors_from_scalar(&values);
        assert_eq!(colors[1], colors[3], "NaN draws as the midpoint");
        assert_ne!(colors[0], colors[2]);
        assert_eq!(colors, colors_from_scalar_in(&values, 1.0, 3.0));
        assert_eq!(scalar_range(&[f32::NAN]), (0.0, 0.0));
        assert_eq!(pack_rgba(1, 2, 3, 4).to_le_bytes(), [1, 2, 3, 4]);
    }

    #[test]
    fn b_factor_coloring_varies_on_a_real_structure() {
        // 4HHB has real, non-uniform recorded b-factors; this is the one
        // hop in the coloring path (`topology.b_factor` -> min/max fold ->
        // per-atom color) that the headless GPU tests don't exercise, since
        // synthetic benchmark structures always report b_factor == 0.0
        // everywhere (see vv-io's synth module) and so always render white.
        let structure = vv_io::load(small_cif_path("4HHB.cif")).unwrap();
        let colors = colors_for(ColorScheme::BFactor, &structure.topology);
        assert_eq!(colors.len(), structure.topology.atom_count());
        assert!(
            colors.iter().any(|&c| c != colors[0]),
            "4HHB has varying b-factors; coloring should not be uniform"
        );
    }
}
