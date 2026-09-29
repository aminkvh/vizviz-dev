//! A coarse world-space occupancy volume of everything drawn opaque, for
//! long-range ambient occlusion and large-scale shadows
//! (`shaders/ao.wgsl`). Screen-space AO only sees what is on screen near
//! a pixel, so a ribbon in front of another ribbon never darkens it; the
//! volume holds the whole scene, seen from any side.
//!
//! Built on the CPU from proxy spheres (the caller decides what stands in
//! for each representation) into about [`MAX_DIM`] voxels per axis,
//! with a mip chain of 2x2x2 averages for cone tracing. World space, so
//! turning the camera never rebuilds it; only a change in what is drawn
//! does.

use std::sync::atomic::{AtomicU32, Ordering};

use glam::Vec3;
use rayon::prelude::*;

use crate::context::GpuContext;

/// Voxels along the longest axis, before padding to [`ALIGN`].
pub const MAX_DIM: u32 = 160;
/// Voxels are never finer than this (Angstrom).
const MIN_VOXEL: f32 = 1.0;
/// Fixed-point scale for the atomic accumulation.
const FIXED: f32 = 65536.0;
/// Mips built. Enough for the widest footprint `shaders/ao.wgsl` samples
/// (12 voxels, lod 3.6); the sampler clamps anything coarser to the last.
const MIP_LEVELS: u32 = 5;
/// Every dimension is a multiple of this, so each built mip is exactly
/// half the one below: `ao.wgsl` maps a point to one `uvw` for every mip,
/// which lands on the same world point only when no mip rounds a size
/// down (else coarse mips stretch toward the grid's far corner).
const ALIGN: u32 = 1 << (MIP_LEVELS - 1);

/// The voxel grid around a box (corner `min`, size `extent`): its corner
/// and dimensions, at least one empty voxel on every side so cones leaving
/// the scene sample 0, padded to [`ALIGN`] evenly on both sides.
fn grid(min: Vec3, extent: Vec3, voxel: f32) -> (Vec3, [u32; 3]) {
    let inner = (extent / voxel).ceil().as_uvec3() + 2;
    let dims = inner.map(|d| d.next_multiple_of(ALIGN));
    let pad = (dims - inner).as_vec3() * 0.5;
    let origin = min - (pad.floor() + 1.0) * voxel;
    (origin, dims.to_array())
}

pub struct OcclusionVolume {
    pub(crate) view: wgpu::TextureView,
    /// World position of voxel (0, 0, 0)'s corner.
    pub origin: Vec3,
    pub voxel: f32,
    pub dims: [u32; 3],
    pub mip_count: u32,
}

impl OcclusionVolume {
    /// One empty voxel: what the renderer binds when there is no volume.
    pub(crate) fn empty(ctx: &GpuContext) -> Self {
        Self::upload(ctx, Vec3::ZERO, 1.0, vec![([1, 1, 1], vec![0.0])])
    }

    /// Occupancy (filled fraction of each voxel, 0..1) of `spheres`
    /// (centre, radius). `None` when there is nothing to hold.
    pub fn build(ctx: &GpuContext, spheres: &[(Vec3, f32)]) -> Option<Self> {
        let (min, max) = spheres.iter().fold(
            (Vec3::splat(f32::INFINITY), Vec3::splat(f32::NEG_INFINITY)),
            |(lo, hi), &(c, r)| (lo.min(c - r), hi.max(c + r)),
        );
        if !min.is_finite() || !max.is_finite() {
            return None;
        }
        let extent = max - min;
        let voxel = (extent.max_element() / MAX_DIM as f32).max(MIN_VOXEL);
        let (origin, dims) = grid(min, extent, voxel);
        let level0 = splat(spheres, origin, voxel, dims);

        let mut levels = vec![(dims, level0)];
        while (levels.len() as u32) < MIP_LEVELS {
            let (d, data) = levels.last().unwrap();
            levels.push(downsample(*d, data));
        }
        Some(Self::upload(ctx, origin, voxel, levels))
    }

    fn upload(
        ctx: &GpuContext,
        origin: Vec3,
        voxel: f32,
        levels: Vec<([u32; 3], Vec<f32>)>,
    ) -> Self {
        let dims = levels[0].0;
        let texture = ctx.device.create_texture(&wgpu::TextureDescriptor {
            label: Some("occlusion volume"),
            size: wgpu::Extent3d {
                width: dims[0],
                height: dims[1],
                depth_or_array_layers: dims[2],
            },
            mip_level_count: levels.len() as u32,
            sample_count: 1,
            dimension: wgpu::TextureDimension::D3,
            format: wgpu::TextureFormat::R8Unorm,
            usage: wgpu::TextureUsages::TEXTURE_BINDING | wgpu::TextureUsages::COPY_DST,
            view_formats: &[],
        });
        for (mip, (d, data)) in levels.iter().enumerate() {
            let bytes: Vec<u8> = data
                .iter()
                .map(|v| (v.clamp(0.0, 1.0) * 255.0).round() as u8)
                .collect();
            ctx.queue.write_texture(
                wgpu::TexelCopyTextureInfo {
                    texture: &texture,
                    mip_level: mip as u32,
                    origin: wgpu::Origin3d::ZERO,
                    aspect: wgpu::TextureAspect::All,
                },
                &bytes,
                wgpu::TexelCopyBufferLayout {
                    offset: 0,
                    bytes_per_row: Some(d[0]),
                    rows_per_image: Some(d[1]),
                },
                wgpu::Extent3d {
                    width: d[0],
                    height: d[1],
                    depth_or_array_layers: d[2],
                },
            );
        }
        Self {
            view: texture.create_view(&Default::default()),
            origin,
            voxel,
            dims,
            mip_count: levels.len() as u32,
        }
    }
}

/// Each sphere's volume spread over the voxels it covers: a sphere
/// smaller than a voxel adds its volume fraction trilinearly around its
/// centre; a larger one adds each covered voxel's approximate filled
/// fraction.
fn splat(spheres: &[(Vec3, f32)], origin: Vec3, voxel: f32, dims: [u32; 3]) -> Vec<f32> {
    let [nx, ny, nz] = dims.map(|d| d as usize);
    let cells: Vec<AtomicU32> = (0..nx * ny * nz).map(|_| AtomicU32::new(0)).collect();
    let add = |x: usize, y: usize, z: usize, v: f32| {
        if x < nx && y < ny && z < nz && v > 0.0 {
            cells[(z * ny + y) * nx + x].fetch_add((v * FIXED).round() as u32, Ordering::Relaxed);
        }
    };
    let voxel_volume = voxel * voxel * voxel;
    spheres.par_iter().for_each(|&(c, r)| {
        let g = (c - origin) / voxel - 0.5;
        if r < voxel {
            let fraction = 4.0 / 3.0 * std::f32::consts::PI * r * r * r / voxel_volume;
            let base = g.floor();
            let t = g - base;
            for corner in 0..8 {
                let (dx, dy, dz) = (corner & 1, (corner >> 1) & 1, (corner >> 2) & 1);
                let w = (if dx == 1 { t.x } else { 1.0 - t.x })
                    * (if dy == 1 { t.y } else { 1.0 - t.y })
                    * (if dz == 1 { t.z } else { 1.0 - t.z });
                let p = base + Vec3::new(dx as f32, dy as f32, dz as f32);
                if p.min_element() >= 0.0 {
                    add(p.x as usize, p.y as usize, p.z as usize, fraction * w);
                }
            }
        } else {
            let reach = r / voxel + 1.0;
            let lo = (g - reach).floor().max(Vec3::ZERO);
            let hi = (g + reach).ceil();
            for z in lo.z as usize..=hi.z as usize {
                for y in lo.y as usize..=hi.y as usize {
                    for x in lo.x as usize..=hi.x as usize {
                        let d = (Vec3::new(x as f32, y as f32, z as f32) - g).length() * voxel;
                        add(x, y, z, ((r - d) / voxel + 0.5).clamp(0.0, 1.0));
                    }
                }
            }
        }
    });
    cells
        .into_iter()
        .map(|c| (c.into_inner() as f32 / FIXED).min(1.0))
        .collect()
}

/// The next mip, sized as wgpu sizes mips (halved, rounded down): each
/// voxel the mean of the 2x2x2 below it.
fn downsample(dims: [u32; 3], data: &[f32]) -> ([u32; 3], Vec<f32>) {
    let [nx, ny, nz] = dims.map(|d| d as usize);
    let next = dims.map(|d| (d / 2).max(1));
    let [mx, my, mz] = next.map(|d| d as usize);
    let out = (0..mx * my * mz)
        .into_par_iter()
        .map(|i| {
            let (x, y, z) = (i % mx, (i / mx) % my, i / (mx * my));
            let mut sum = 0.0;
            for (dz, dy, dx) in (0..8).map(|k| (k >> 2, (k >> 1) & 1, k & 1)) {
                let (sx, sy, sz) = (2 * x + dx, 2 * y + dy, 2 * z + dz);
                if sx < nx && sy < ny && sz < nz {
                    sum += data[(sz * ny + sy) * nx + sx];
                }
            }
            sum / 8.0
        })
        .collect();
    (next, out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_big_sphere_fills_its_centre_and_leaves_the_corners_empty() {
        let dims = [16, 16, 16];
        let v = splat(&[(Vec3::splat(8.0), 5.0)], Vec3::ZERO, 1.0, dims);
        let at = |x: usize, y: usize, z: usize| v[(z * 16 + y) * 16 + x];
        assert_eq!(at(7, 7, 7), 1.0);
        assert_eq!(at(0, 0, 0), 0.0);
        let total: f32 = v.iter().sum();
        let exact = 4.0 / 3.0 * std::f32::consts::PI * 125.0;
        assert!((total - exact).abs() / exact < 0.1, "{total} vs {exact}");
    }

    #[test]
    fn every_mip_halves_exactly_and_the_box_keeps_a_margin() {
        let (min, extent, voxel) = (Vec3::new(-3.2, 7.0, 0.4), Vec3::new(64.3, 9.9, 171.0), 1.1);
        let (origin, dims) = grid(min, extent, voxel);
        for (d, (lo, hi)) in dims.iter().zip([
            (
                min.x - origin.x,
                origin.x + dims[0] as f32 * voxel - (min.x + extent.x),
            ),
            (
                min.y - origin.y,
                origin.y + dims[1] as f32 * voxel - (min.y + extent.y),
            ),
            (
                min.z - origin.z,
                origin.z + dims[2] as f32 * voxel - (min.z + extent.z),
            ),
        ]) {
            assert_eq!(d % ALIGN, 0, "{dims:?}");
            assert!(lo >= voxel * 0.999 && hi >= voxel * 0.999, "{lo} {hi}");
        }
    }

    #[test]
    fn small_spheres_keep_their_volume_and_mips_keep_the_mean() {
        let dims = [8, 8, 8];
        let spheres: Vec<_> = (0..20)
            .map(|i| {
                (
                    Vec3::new(1.3 + (i % 5) as f32 * 1.1, 2.7, 3.1 + (i / 5) as f32),
                    0.3,
                )
            })
            .collect();
        let v = splat(&spheres, Vec3::ZERO, 1.0, dims);
        let total: f32 = v.iter().sum();
        let exact = 20.0 * 4.0 / 3.0 * std::f32::consts::PI * 0.027;
        assert!((total - exact).abs() / exact < 0.02, "{total} vs {exact}");
        let (next, half) = downsample(dims, &v);
        assert_eq!(next, [4, 4, 4]);
        let mean0 = total / 512.0;
        let mean1 = half.iter().sum::<f32>() / 64.0;
        assert!((mean0 - mean1).abs() < 1e-5);
    }
}
