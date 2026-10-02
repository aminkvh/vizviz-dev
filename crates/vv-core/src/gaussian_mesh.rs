//! The Gaussian surface as triangles, for renderers that take meshes (the
//! path tracer): naive surface nets (Gibson 1998, "Constrained elastic
//! surface nets") over the density sampled on a grid, then each vertex
//! moved onto the exact isosurface by Newton steps along the exact
//! gradient ([`crate::gaussian_surface::gradient`]), which also gives its
//! normal. So the grid only decides the triangles' size, not where the
//! surface is or how it is shaded.

use glam::Vec3;
use rayon::prelude::*;

use crate::cartoon::ExpandedMesh;
use crate::gaussian_surface::cutoff_radius;
use crate::spatial::Grid;
use crate::surface_net::{self, Isosurface, NetGrid};

/// An atom's field is 1 at its van der Waals surface
/// ([`crate::gaussian_surface`]); the surface is where the sum is 1.
pub const ISOVALUE: f32 = 1.0;
const NEWTON_STEPS: usize = 4;

/// Voxels a mesh may sample, at most (4 bytes each).
pub const MAX_VOXELS: usize = 1 << 26;

/// The Gaussian surface of atoms at `positions` with van der Waals
/// `radii`, each atom's field cut off where it drops below `epsilon`, on
/// a grid of `voxel` angstroms (raised as needed to stay within
/// [`MAX_VOXELS`]). Each vertex's `source` is the atom contributing most
/// there, the atom the viewport colours it by.
pub fn mesh(
    positions: &[Vec3],
    radii: &[f32],
    blob_factor: f32,
    epsilon: f32,
    voxel: f32,
) -> ExpandedMesh {
    assert_eq!(positions.len(), radii.len());
    let field = Field::new(positions, radii, blob_factor, epsilon);
    if positions.is_empty() {
        return ExpandedMesh::default();
    }
    let grid = SampleGrid::around(&field, voxel);
    let values = grid.sample(&field);
    let net = surface_net::build(
        &grid.net,
        &Isosurface {
            grid: &grid.net,
            values: &values,
            level: ISOVALUE,
        },
    );
    let indices = net.indices;
    let refined: Vec<(Vec3, Vec3, u32)> = net
        .cells
        .par_iter()
        .map(|&(_, p)| field.onto_surface(p, grid.net.voxel))
        .collect();
    let mut out = ExpandedMesh {
        positions: refined.iter().map(|v| v.0).collect(),
        normals: refined.iter().map(|v| v.1).collect(),
        indices,
        source: refined.iter().map(|v| v.2).collect(),
    };
    orient(&mut out);
    out
}

/// The exact field, with a neighbour grid for its sums.
struct Field<'a> {
    positions: &'a [Vec3],
    radii: &'a [f32],
    cutoffs: Vec<f32>,
    reach: f32,
    blob_factor: f32,
    neighbors: Grid,
}

impl<'a> Field<'a> {
    fn new(positions: &'a [Vec3], radii: &'a [f32], blob_factor: f32, epsilon: f32) -> Self {
        let cutoffs: Vec<f32> = radii
            .iter()
            .map(|&r| cutoff_radius(r, blob_factor, epsilon))
            .collect();
        let reach = cutoffs.iter().copied().fold(1e-3, f32::max);
        let all: Vec<u32> = (0..positions.len() as u32).collect();
        Self {
            positions,
            radii,
            cutoffs,
            reach,
            blob_factor,
            neighbors: Grid::build(positions, &all, reach),
        }
    }

    fn contribution(&self, a: usize, d2: f32) -> f32 {
        let r = self.radii[a];
        (-self.blob_factor * (d2 / (r * r) - 1.0)).exp()
    }

    /// Density, gradient and the most contributing atom at `p`, over the
    /// atoms within their cutoff (the same atoms the viewport sums).
    fn at(&self, p: Vec3) -> (f32, Vec3, u32) {
        let (mut f, mut g, mut best) = (0.0, Vec3::ZERO, (0.0, 0u32));
        self.neighbors
            .for_each_within(self.positions, p, self.reach, |a, d2| {
                let k = a as usize;
                if d2 > self.cutoffs[k] * self.cutoffs[k] {
                    return;
                }
                let c = self.contribution(k, d2);
                let r = self.radii[k];
                f += c;
                g += (p - self.positions[k]) * (-2.0 * self.blob_factor / (r * r) * c);
                if c > best.0 {
                    best = (c, a);
                }
            });
        (f, g, best.1)
    }

    /// `p` moved along the gradient onto the isosurface (by at most one
    /// voxel), its outward normal and its atom.
    fn onto_surface(&self, start: Vec3, voxel: f32) -> (Vec3, Vec3, u32) {
        let mut p = start;
        for _ in 0..NEWTON_STEPS {
            let (f, g, _) = self.at(p);
            let g2 = g.length_squared();
            if g2 < 1e-12 {
                break;
            }
            p -= g * ((f - ISOVALUE) / g2);
            let moved = p - start;
            if moved.length() > voxel {
                p = start + moved.normalize() * voxel;
            }
        }
        let (_, g, atom) = self.at(p);
        (p, (-g).normalize_or(Vec3::Z), atom)
    }
}

/// A grid a voxel beyond every atom's cutoff on each side, so the surface
/// closes inside it.
struct SampleGrid {
    net: NetGrid,
}

impl SampleGrid {
    fn around(field: &Field, voxel: f32) -> Self {
        let (lo, hi) = field.positions.iter().zip(&field.cutoffs).fold(
            (Vec3::splat(f32::INFINITY), Vec3::splat(f32::NEG_INFINITY)),
            |(lo, hi), (&p, &c)| (lo.min(p - c), hi.max(p + c)),
        );
        Self {
            net: NetGrid::around(lo, hi, voxel, MAX_VOXELS),
        }
    }

    /// The density at every grid point: each plane adds up the atoms
    /// whose cutoff reaches it, over the disc they cover there.
    fn sample(&self, field: &Field) -> Vec<f32> {
        let NetGrid {
            origin,
            voxel,
            dims: [nx, ny, nz],
        } = self.net;
        let mut by_z: Vec<u32> = (0..field.positions.len() as u32).collect();
        by_z.sort_by(|&a, &b| {
            field.positions[a as usize]
                .z
                .total_cmp(&field.positions[b as usize].z)
        });
        let zs: Vec<f32> = by_z
            .iter()
            .map(|&a| field.positions[a as usize].z)
            .collect();
        let mut values = vec![0.0f32; nx * ny * nz];
        values
            .par_chunks_mut(nx * ny)
            .enumerate()
            .for_each(|(k, plane)| {
                let z = origin.z + k as f32 * voxel;
                let from = zs.partition_point(|&az| az < z - field.reach);
                let to = zs.partition_point(|&az| az <= z + field.reach);
                for &a in &by_z[from..to] {
                    let a = a as usize;
                    let c = field.positions[a];
                    let cut = field.cutoffs[a];
                    let dz = z - c.z;
                    let disc = cut * cut - dz * dz;
                    if disc < 0.0 {
                        continue;
                    }
                    let r = disc.sqrt();
                    let range = |center: f32, o: f32, n: usize| {
                        let first = ((center - r - o) / voxel).ceil().max(0.0) as usize;
                        let last = (((center + r - o) / voxel).floor().max(-1.0) + 1.0) as usize;
                        first..last.min(n)
                    };
                    for j in range(c.y, origin.y, ny) {
                        let y = origin.y + j as f32 * voxel;
                        for i in range(c.x, origin.x, nx) {
                            let x = origin.x + i as f32 * voxel;
                            let d2 = (x - c.x).powi(2) + (y - c.y).powi(2) + dz * dz;
                            if d2 <= cut * cut {
                                plane[i + nx * j] += field.contribution(a, d2);
                            }
                        }
                    }
                }
            });
        values
    }
}

/// Winds every triangle counter-clockwise seen from outside (along its
/// vertices' normals), as the tracer reads which side a ray hit.
fn orient(mesh: &mut ExpandedMesh) {
    for t in mesh.indices.chunks_exact_mut(3) {
        let [a, b, c] = [t[0], t[1], t[2]].map(|v| v as usize);
        let p = &mesh.positions;
        let face = (p[b] - p[a]).cross(p[c] - p[a]);
        let n = mesh.normals[a] + mesh.normals[b] + mesh.normals[c];
        if face.dot(n) < 0.0 {
            t.swap(1, 2);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gaussian_surface::{density, Blob};

    #[test]
    fn a_lone_atoms_mesh_is_its_sphere() {
        let center = Vec3::new(1.0, -2.0, 3.0);
        let mesh = mesh(&[center], &[1.7], 2.0, 0.01, 0.4);
        assert!(mesh.indices.len() > 300, "{}", mesh.indices.len());
        for (p, n) in mesh.positions.iter().zip(&mesh.normals) {
            assert!(((*p - center).length() - 1.7).abs() < 1e-3, "{p}");
            assert!((*n - (*p - center).normalize()).length() < 1e-3);
        }
    }

    /// Vertices lie on the isosurface (the exact field less what the
    /// cutoff drops), the mesh is closed (every edge in exactly two
    /// triangles), and it faces outward (positive enclosed volume).
    #[test]
    fn three_blobs_mesh_is_closed_outward_and_on_the_surface() {
        let positions = [
            Vec3::ZERO,
            Vec3::new(2.4, 0.3, 0.0),
            Vec3::new(1.0, 2.0, 0.5),
        ];
        let radii = [1.7, 1.5, 1.2];
        let mesh = mesh(&positions, &radii, 2.0, 0.01, 0.3);
        let blobs: Vec<Blob> = positions
            .iter()
            .zip(radii)
            .map(|(&center, radius)| Blob { center, radius })
            .collect();
        for p in &mesh.positions {
            let f = density(*p, &blobs, 2.0);
            assert!((f - ISOVALUE).abs() < 0.03, "{f}");
        }
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
        let spheres: f32 = radii
            .iter()
            .map(|r| 4.0 / 3.0 * std::f32::consts::PI * r * r * r)
            .sum();
        assert!(
            volume > 0.5 * spheres && volume < spheres,
            "{volume} vs {spheres}"
        );
    }
}
