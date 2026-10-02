//! The solvent-excluded surface as triangles, for the mesh display. A point
//! is inside the excluded volume when it is farther than the probe radius
//! from every place the probe's centre can be (outside every atom grown by
//! the probe) -- the field [`crate::ses::Ses::probe_distance`] evaluates
//! patch by patch. Here it is found on a grid instead: the Euclidean
//! distance transform (Felzenszwalb & Huttenlocher 2012, Theory of
//! Computing 8:415) of the free grid points gives the distance to the
//! nearest probe position at every point, and surface nets
//! ([`crate::surface_net`]) mesh where it equals the probe radius. Accurate
//! to about a voxel; the viewport's solid surface is exact.

use glam::Vec3;
use rayon::prelude::*;

use crate::cartoon::ExpandedMesh;
use crate::spatial::Grid;
use crate::surface_net::{self, Isosurface, NetGrid};

/// Voxels a mesh may sample, at most (two 4-byte arrays at once).
pub const MAX_VOXELS: usize = 1 << 25;

/// Stands for "no probe position on this line yet"; far beyond any squared
/// distance in voxels, small enough that `FAR + q * q` stays exact in f32.
const FAR: f32 = 1e8;

/// The SES of atoms at `positions` with van der Waals `radii` for a probe
/// of radius `probe`, on a grid of `voxel` angstroms (raised as needed to
/// stay within [`MAX_VOXELS`]). Each vertex's `source` is the atom whose
/// surface is nearest it.
pub fn mesh(positions: &[Vec3], radii: &[f32], probe: f32, voxel: f32) -> ExpandedMesh {
    assert_eq!(positions.len(), radii.len());
    if positions.is_empty() {
        return ExpandedMesh::default();
    }
    let reach = radii.iter().copied().fold(0.0, f32::max) + probe;
    let (lo, hi) = positions.iter().zip(radii).fold(
        (Vec3::splat(f32::INFINITY), Vec3::splat(f32::NEG_INFINITY)),
        |(lo, hi), (&p, &r)| (lo.min(p - (r + probe)), hi.max(p + (r + probe))),
    );
    let grid = NetGrid::around(lo, hi, voxel, MAX_VOXELS);
    let mut field = probe_positions(&grid, positions, radii, probe);
    distance_transform(&grid, &mut field);
    let scale = grid.voxel;
    field.par_iter_mut().for_each(|d2| *d2 = d2.sqrt() * scale);
    let net = surface_net::build(
        &grid,
        &Isosurface {
            grid: &grid,
            values: &field,
            level: probe,
        },
    );
    let all: Vec<u32> = (0..positions.len() as u32).collect();
    let neighbors = Grid::build(positions, &all, reach);
    surface_net::mesh_of(net, |p| {
        surface_net::nearest_surface_atom(&neighbors, positions, radii, reach, p)
    })
}

/// `0` at every grid point outside all atoms grown by the probe (a place the
/// probe's centre can be), [`FAR`] elsewhere.
fn probe_positions(grid: &NetGrid, positions: &[Vec3], radii: &[f32], probe: f32) -> Vec<f32> {
    let [nx, ny, nz] = grid.dims;
    let grown: Vec<f32> = radii.iter().map(|r| r + probe).collect();
    let reach = grown.iter().copied().fold(0.0, f32::max);
    let mut by_z: Vec<u32> = (0..positions.len() as u32).collect();
    by_z.sort_by(|&a, &b| positions[a as usize].z.total_cmp(&positions[b as usize].z));
    let zs: Vec<f32> = by_z.iter().map(|&a| positions[a as usize].z).collect();
    let mut values = vec![0.0f32; nx * ny * nz];
    values
        .par_chunks_mut(nx * ny)
        .enumerate()
        .for_each(|(k, plane)| {
            let z = grid.origin.z + k as f32 * grid.voxel;
            let from = zs.partition_point(|&az| az < z - reach);
            let to = zs.partition_point(|&az| az <= z + reach);
            for &a in &by_z[from..to] {
                let (c, big) = (positions[a as usize], grown[a as usize]);
                let disc = big * big - (z - c.z) * (z - c.z);
                if disc <= 0.0 {
                    continue;
                }
                let r = disc.sqrt();
                let span = |center: f32, o: f32, n: usize| {
                    let first = ((center - r - o) / grid.voxel).ceil().max(0.0) as usize;
                    let last = (((center + r - o) / grid.voxel).floor().max(-1.0) + 1.0) as usize;
                    first..last.min(n)
                };
                for j in span(c.y, grid.origin.y, ny) {
                    let y = grid.origin.y + j as f32 * grid.voxel;
                    for i in span(c.x, grid.origin.x, nx) {
                        let x = grid.origin.x + i as f32 * grid.voxel;
                        if (x - c.x).powi(2) + (y - c.y).powi(2) + (z - c.z).powi(2) < big * big {
                            plane[i + nx * j] = FAR;
                        }
                    }
                }
            }
        });
    values
}

/// Turns `field` (`0` at sources, [`FAR`] elsewhere) into each point's
/// squared distance, in voxels, to the nearest source: the 1-D transform
/// along x, then y, then z.
fn distance_transform(grid: &NetGrid, field: &mut [f32]) {
    let [nx, ny, nz] = grid.dims;
    field
        .par_chunks_mut(nx)
        .for_each_init(|| Envelope::new(nx), |env, line| env.transform(line));
    field.par_chunks_mut(nx * ny).for_each_init(
        || (Envelope::new(ny), vec![0.0; ny]),
        |(env, column), plane| {
            for i in 0..nx {
                for (j, v) in column.iter_mut().enumerate() {
                    *v = plane[i + nx * j];
                }
                env.transform(column);
                for (j, v) in column.iter().enumerate() {
                    plane[i + nx * j] = *v;
                }
            }
        },
    );
    let columns: Vec<Vec<f32>> = (0..nx * ny)
        .into_par_iter()
        .map_init(
            || Envelope::new(nz),
            |env, ij| {
                let mut column: Vec<f32> = (0..nz).map(|k| field[ij + nx * ny * k]).collect();
                env.transform(&mut column);
                column
            },
        )
        .collect();
    for (ij, column) in columns.iter().enumerate() {
        for (k, &v) in column.iter().enumerate() {
            field[ij + nx * ny * k] = v;
        }
    }
}

/// Scratch for the lower-envelope pass over lines of one length.
struct Envelope {
    hull: Vec<usize>,
    cuts: Vec<f32>,
    out: Vec<f32>,
}

impl Envelope {
    fn new(n: usize) -> Self {
        Self {
            hull: vec![0; n],
            cuts: vec![0.0; n + 1],
            out: vec![0.0; n],
        }
    }

    /// `line[q]` becomes `min over p of (q - p)^2 + line[p]`.
    fn transform(&mut self, line: &mut [f32]) {
        let n = line.len();
        let (hull, cuts) = (&mut self.hull, &mut self.cuts);
        let meet = |q: usize, p: usize| {
            ((line[q] + (q * q) as f32) - (line[p] + (p * p) as f32))
                / (2.0 * (q as f32 - p as f32))
        };
        let mut k = 0;
        hull[0] = 0;
        cuts[0] = f32::NEG_INFINITY;
        cuts[1] = f32::INFINITY;
        for q in 1..n {
            let mut s = meet(q, hull[k]);
            while s <= cuts[k] {
                k -= 1;
                s = meet(q, hull[k]);
            }
            k += 1;
            hull[k] = q;
            cuts[k] = s;
            cuts[k + 1] = f32::INFINITY;
        }
        k = 0;
        for q in 0..n {
            while cuts[k + 1] < q as f32 {
                k += 1;
            }
            let p = hull[k];
            self.out[q] = ((q as f32 - p as f32).powi(2)) + line[p];
        }
        line.copy_from_slice(&self.out[..n]);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_lone_atoms_surface_is_its_van_der_waals_sphere() {
        let center = Vec3::new(1.0, -2.0, 3.0);
        let mesh = mesh(&[center], &[1.7], 1.4, 0.25);
        assert!(mesh.indices.len() > 500, "{}", mesh.indices.len());
        // The distance field is only known to a voxel, so the surface is
        // too: bumpy by under half a voxel, but outward everywhere.
        for (p, n) in mesh.positions.iter().zip(&mesh.normals) {
            assert!(((*p - center).length() - 1.7).abs() < 0.15, "{p}");
            assert!(n.dot((*p - center).normalize()) > 0.5, "{n}");
        }
    }

    #[test]
    fn two_overlapping_atoms_give_a_closed_outward_mesh_with_a_filled_crevice() {
        let positions = [Vec3::ZERO, Vec3::new(3.0, 0.0, 0.0)];
        let mesh = mesh(&positions, &[1.7, 1.7], 1.4, 0.25);
        let mut edges = std::collections::HashMap::new();
        let mut volume = 0.0;
        for t in mesh.indices.chunks_exact(3) {
            for k in 0..3 {
                let (a, b) = (t[k], t[(k + 1) % 3]);
                *edges.entry((a.min(b), a.max(b))).or_insert(0) += 1;
            }
            let [a, b, c] = [t[0], t[1], t[2]].map(|v| mesh.positions[v as usize]);
            volume += a.dot(b.cross(c)) / 6.0;
        }
        assert!(edges.values().all(|&n| n == 2), "open mesh");
        // The two spheres overlap by a 0.4 A lens: their union is 40.8 A^3,
        // and the filled crevice adds to it.
        assert!(volume > 40.8 && volume < 50.0, "{volume}");
        // The probe cannot reach into the crease between the atoms, so the
        // surface there is the probe's own torus, 3.1 A (atom plus probe)
        // from both centres: sqrt(3.1^2 - 1.5^2) - 1.4 = 1.31 A from the
        // axis, where the spheres alone would cut in to 0.8.
        let waist = mesh
            .positions
            .iter()
            .filter(|p| (p.x - 1.5).abs() < 0.3)
            .map(|p| (p.y * p.y + p.z * p.z).sqrt())
            .fold(f32::INFINITY, f32::min);
        assert!((waist - 1.31).abs() < 0.2, "{waist}");
    }

    #[test]
    fn the_distance_transform_of_two_sources_is_the_nearer_ones_squared_distance() {
        let grid = NetGrid {
            origin: Vec3::ZERO,
            voxel: 1.0,
            dims: [9, 4, 3],
        };
        let mut field = vec![FAR; grid.len()];
        field[grid.index(1, 1, 1)] = 0.0;
        field[grid.index(7, 2, 0)] = 0.0;
        distance_transform(&grid, &mut field);
        for k in 0..3usize {
            for j in 0..4usize {
                for i in 0..9usize {
                    let near = |a: [usize; 3]| {
                        (i as f32 - a[0] as f32).powi(2)
                            + (j as f32 - a[1] as f32).powi(2)
                            + (k as f32 - a[2] as f32).powi(2)
                    };
                    let want = near([1, 1, 1]).min(near([7, 2, 0]));
                    assert!((field[grid.index(i, j, k)] - want).abs() < 1e-3);
                }
            }
        }
    }
}
