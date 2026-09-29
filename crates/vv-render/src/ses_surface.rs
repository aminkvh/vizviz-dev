//! The solvent-excluded surface on the GPU: `vv_core::ses::Ses` uploaded
//! as one 16-byte record per patch, drawn as one ray-cast billboard each
//! by `shaders/ses_surface.wgsl` (`shaders/ses_patch.wgsl` is the WGSL
//! twin of `Ses::cast`).
//!
//! Millions of patches are culled the way atoms are: each patch's
//! bounding sphere is an atom of [`SesGpu::bounds`], convex patches first
//! and in spatial order within a kind so runs of patches make compact
//! clusters, and the atom cull passes decide
//! which patches draw -- as billboards, or as points below a pixel.

use bytemuck::{Pod, Zeroable};
use glam::Vec3;
use rayon::prelude::*;
use vv_core::ses::{circle, torus_bound, Ses};

use crate::context::GpuContext;
use crate::scene::{try_buffer_init, GpuStructure, OutOfGpuMemory, PAGE_ATOMS};
use crate::style::Material;

/// Top two bits of a patch record's first word.
pub(crate) const KIND_CONVEX: u32 = 0;
pub(crate) const KIND_TORUS: u32 = 1 << 30;
pub(crate) const KIND_CONCAVE: u32 = 2 << 30;

/// Layout must match `SesParams` in `shaders/ses_surface.wgsl`.
#[repr(C)]
#[derive(Clone, Copy, Debug, Pod, Zeroable)]
pub struct SesParams {
    pub view_inv: [[f32; 4]; 4],
    pub probe_radius: f32,
    pub atom_count: u32,
    pub patch_count: u32,
    /// Frame-wide pick id of atom 0; a hit picks the atom it lies nearest.
    pub id_base: u32,
    pub material: Material,
}

/// A probe position as the shader reads it: centre, its three atoms and
/// its overlapping probes' range in `probe_neighbors`.
#[repr(C)]
#[derive(Clone, Copy, Debug, Pod, Zeroable)]
pub(crate) struct ProbeGpu {
    pub center: [f32; 3],
    pub atoms: [u32; 3],
    pub neighbors: [u32; 2],
}

/// Every patch of `ses` as (record, bounding sphere centre and radius,
/// the atom it shows as): convex patches, tori, then concave patches.
/// `positions` and `radii` (van der Waals) are the atoms `ses` was built
/// from.
pub(crate) fn patch_records(
    ses: &Ses,
    positions: &[Vec3],
    radii: &[f32],
) -> Vec<([u32; 4], Vec3, f32, u32)> {
    assert_eq!(positions.len(), radii.len());
    let rp = ses.probe_radius;
    let mut patches = Vec::with_capacity(ses.convex.len() + ses.tori.len() + ses.probes.len());
    patches.par_extend(ses.convex.par_iter().map(|c| {
        let a = c.atom as usize;
        (
            [KIND_CONVEX | c.atom, c.caps[0], c.caps[1], 0],
            positions[a],
            radii[a],
            c.atom,
        )
    }));
    patches.par_extend(ses.tori.par_iter().map(|t| {
        let [i, j] = t.atoms.map(|a| a as usize);
        let (center, radius) = circle(positions[i], radii[i] + rp, positions[j], radii[j] + rp)
            .map_or((positions[i], radii[i] + rp), |c| {
                torus_bound(&c, radii[i], radii[j], rp, ses.torus_arc(t, &c))
            });
        let [a, b] = t.probes;
        (
            [KIND_TORUS | t.atoms[0], t.atoms[1], a, b],
            center,
            radius,
            t.atoms[0],
        )
    }));
    patches.par_extend(
        ses.probes
            .par_iter()
            .enumerate()
            .map(|(k, p)| ([KIND_CONCAVE | k as u32, 0, 0, 0], p.center, rp, p.atoms[0])),
    );
    patches
}

pub(crate) fn probe_records(ses: &Ses) -> Vec<ProbeGpu> {
    ses.probes
        .iter()
        .enumerate()
        .map(|(k, p)| ProbeGpu {
            center: p.center.to_array(),
            atoms: p.atoms,
            neighbors: [
                ses.probe_neighbor_starts[k],
                ses.probe_neighbor_starts[k + 1],
            ],
        })
        .collect()
}

/// An SES laid out for the GPU, on the CPU: seconds of sorting and
/// packing at millions of patches, so it is built off the UI thread
/// (next to the SES itself) and only copied up by [`SesGpu::upload`].
pub struct SesLayout {
    atoms: Vec<[f32; 4]>,
    /// In drawing order: `[kind | a, b, c, d]`.
    records: Vec<[u32; 4]>,
    centers: Vec<Vec3>,
    bound_radii: Vec<f32>,
    /// The atom whose colour each patch shows as a point.
    patch_atoms: Vec<u32>,
    caps: Vec<u32>,
    probes: Vec<ProbeGpu>,
    probe_neighbors: Vec<u32>,
    probe_radius: f32,
    bounds_center: Vec3,
    bounds_radius: f32,
}

pub struct SesGpu {
    /// xyz = position, w = van der Waals radius.
    pub atoms: wgpu::Buffer,
    pub colors: wgpu::Buffer,
    /// Each page's patch records, as [`SesGpu::bounds`] pages them: no
    /// buffer outgrows a small GPU's local heap (256 MB on an AMD iGPU,
    /// where 8GLV's 18M patches in one buffer failed to allocate).
    pub patches: Vec<wgpu::Buffer>,
    pub caps: wgpu::Buffer,
    pub probes: wgpu::Buffer,
    pub probe_neighbors: wgpu::Buffer,
    pub params: wgpu::Buffer,
    /// Patch `k`'s bounding sphere is atom `k` here, coloured as the
    /// patch's first atom (what a patch below a pixel draws as).
    pub bounds: GpuStructure,
    /// Each page's first patch index (`PageView` in the shader).
    pub page_bases: Vec<wgpu::Buffer>,
    patch_atoms: Vec<u32>,
    pub atom_count: u32,
    pub patch_count: u32,
    pub probe_radius: f32,
    pub bounds_center: Vec3,
    pub bounds_radius: f32,
}

/// Interleaves the low 21 bits of `v` with two zero bits between each.
fn spread(v: u64) -> u64 {
    let mut x = v & 0x1f_ffff;
    x = (x | x << 32) & 0x1f00000000ffff;
    x = (x | x << 16) & 0x1f0000ff0000ff;
    x = (x | x << 8) & 0x100f00f00f00f00f;
    x = (x | x << 4) & 0x10c30c30c30c30c3;
    (x | x << 2) & 0x1249249249249249
}

/// Z-order key of `p` in the box `lo..lo + extent`.
fn morton(p: Vec3, lo: Vec3, extent: f32) -> u64 {
    let q = ((p - lo) / extent * 2_097_151.0).clamp(Vec3::ZERO, Vec3::splat(2_097_151.0));
    spread(q.x as u64) | spread(q.y as u64) << 1 | spread(q.z as u64) << 2
}

impl SesLayout {
    /// `positions` and `radii` (van der Waals) are per atom, the same atoms
    /// `ses` was built from.
    pub fn new(ses: &Ses, positions: &[Vec3], radii: &[f32]) -> Self {
        let rp = ses.probe_radius;
        let mut patches = patch_records(ses, positions, radii);
        let (lo, hi) = patches.iter().fold(
            (Vec3::splat(f32::INFINITY), Vec3::splat(f32::NEG_INFINITY)),
            |(lo, hi), p| (lo.min(p.1), hi.max(p.1)),
        );
        let extent = (hi - lo).max_element().max(1e-3);
        // Convex patches first: drawn first, they hide most fillet
        // fragments from the early depth test (close-up p95 halved).
        // Spatial order within a kind keeps clusters compact.
        patches.par_sort_by_cached_key(|p| (p.0[0] >> 30, morton(p.1, lo, extent)));
        let probes = probe_records(ses);
        Self {
            atoms: positions
                .iter()
                .zip(radii)
                .map(|(p, &r)| [p.x, p.y, p.z, r])
                .collect(),
            records: patches.iter().map(|p| p.0).collect(),
            centers: patches.iter().map(|p| p.1).collect(),
            bound_radii: patches.iter().map(|p| p.2).collect(),
            patch_atoms: patches.iter().map(|p| p.3).collect(),
            caps: ses.caps.clone(),
            probes,
            probe_neighbors: ses.probe_neighbors.clone(),
            probe_radius: rp,
            bounds_center: 0.5 * (lo + hi),
            bounds_radius: (0.5 * (hi - lo).length() + 2.0 * rp).max(1.0),
        }
    }
}

impl SesGpu {
    /// `colors` has one packed colour per atom of `layout`.
    pub fn upload(
        ctx: &GpuContext,
        layout: &SesLayout,
        colors: &[u32],
    ) -> Result<Self, OutOfGpuMemory> {
        assert_eq!(colors.len(), layout.atoms.len());
        // Storage bindings can't be empty.
        let storage = |label: &str, bytes: &[u8]| {
            let bytes = if bytes.is_empty() {
                &[0u8; 16][..]
            } else {
                bytes
            };
            try_buffer_init(ctx, label, bytes, wgpu::BufferUsages::STORAGE)
        };
        let patch_colors: Vec<u32> = layout
            .patch_atoms
            .iter()
            .map(|&a| colors[a as usize])
            .collect();
        let atom_count = layout.atoms.len() as u32;
        let patch_count = layout.records.len() as u32;
        let bounds = GpuStructure::try_from_parts(
            ctx,
            &layout.centers,
            &layout.bound_radii,
            &patch_colors,
            &[],
        )?;
        let page_bases = (0..bounds.pages.len())
            .map(|p| {
                try_buffer_init(
                    ctx,
                    "ses page base",
                    bytemuck::cast_slice(&[(p * PAGE_ATOMS) as u32, 0, 0, 0]),
                    wgpu::BufferUsages::UNIFORM,
                )
            })
            .collect::<Result<_, _>>()?;
        Ok(Self {
            atoms: storage("ses atoms", bytemuck::cast_slice(&layout.atoms))?,
            colors: storage("ses colors", bytemuck::cast_slice(colors))?,
            patches: (0..bounds.pages.len())
                .map(|p| {
                    let n = layout.records.len();
                    let page =
                        &layout.records[(p * PAGE_ATOMS).min(n)..((p + 1) * PAGE_ATOMS).min(n)];
                    storage("ses patches", bytemuck::cast_slice(page))
                })
                .collect::<Result<_, _>>()?,
            caps: storage("ses caps", bytemuck::cast_slice(&layout.caps))?,
            probes: storage("ses probes", bytemuck::cast_slice(&layout.probes))?,
            probe_neighbors: storage(
                "ses probe neighbors",
                bytemuck::cast_slice(&layout.probe_neighbors),
            )?,
            params: try_buffer_init(
                ctx,
                "ses params",
                bytemuck::bytes_of(&SesParams {
                    view_inv: glam::Mat4::IDENTITY.to_cols_array_2d(),
                    probe_radius: layout.probe_radius,
                    atom_count,
                    patch_count,
                    id_base: 0,
                    material: Material::default(),
                }),
                wgpu::BufferUsages::UNIFORM,
            )?,
            bounds,
            page_bases,
            patch_atoms: layout.patch_atoms.clone(),
            atom_count,
            patch_count,
            probe_radius: layout.probe_radius,
            bounds_center: layout.bounds_center,
            bounds_radius: layout.bounds_radius,
        })
    }

    /// New per-atom colours, same atoms.
    pub fn set_colors(&self, ctx: &GpuContext, colors: &[u32]) {
        assert_eq!(colors.len(), self.atom_count as usize);
        if colors.is_empty() {
            return;
        }
        ctx.queue
            .write_buffer(&self.colors, 0, bytemuck::cast_slice(colors));
        let patch_colors: Vec<u32> = self
            .patch_atoms
            .iter()
            .map(|&a| colors[a as usize])
            .collect();
        self.bounds.state.set_colors(ctx, &patch_colors);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn morton_keys_interleave_axes_and_keep_near_points_near() {
        assert_eq!(spread(0b111), 0b1001001);
        let lo = Vec3::ZERO;
        let a = morton(Vec3::new(1.0, 1.0, 1.0), lo, 100.0);
        let b = morton(Vec3::new(1.1, 1.0, 1.0), lo, 100.0);
        let far = morton(Vec3::new(90.0, 90.0, 90.0), lo, 100.0);
        assert!(a.abs_diff(b) < a.abs_diff(far));
    }
}
