//! Surface nets (Gibson 1998, "Constrained elastic surface nets") on a
//! regular grid, for any surface that says which grid points are inside it
//! and where it crosses a grid edge: each cell the surface passes through
//! gets one vertex at the mean of its edges' crossings, and each crossed
//! edge becomes a quad joining the four cells around it. The Gaussian, solvent-excluded
//! and skin surface meshes ([`crate::gaussian_mesh`], [`crate::ses_mesh`],
//! [`crate::skin_mesh`]) differ only in how they answer those two questions.

use glam::Vec3;
use rayon::prelude::*;

use crate::cartoon::ExpandedMesh;
use crate::spatial::Grid;

/// Grid points `origin + voxel * (i, j, k)`.
#[derive(Clone, Copy, Debug)]
pub struct NetGrid {
    pub origin: Vec3,
    pub voxel: f32,
    pub dims: [usize; 3],
}

impl NetGrid {
    /// A grid covering `lo..hi` with a voxel of margin on every side, so a
    /// surface inside the box closes within it; `voxel` is raised as needed
    /// to stay within `max_voxels` points.
    pub fn around(lo: Vec3, hi: Vec3, voxel: f32, max_voxels: usize) -> Self {
        let extent = hi - lo;
        let fits = (extent.x * extent.y * extent.z / max_voxels as f32).cbrt();
        let voxel = voxel.max(fits * 1.01);
        let dims = ((extent / voxel).ceil() + 3.0)
            .to_array()
            .map(|d| d as usize);
        Self {
            origin: lo - Vec3::splat(voxel),
            voxel,
            dims,
        }
    }

    pub fn len(&self) -> usize {
        self.dims.iter().product()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn index(&self, i: usize, j: usize, k: usize) -> usize {
        i + self.dims[0] * (j + self.dims[1] * k)
    }

    pub fn point(&self, i: usize, j: usize, k: usize) -> Vec3 {
        self.origin + Vec3::new(i as f32, j as f32, k as f32) * self.voxel
    }
}

/// A closed surface on a [`NetGrid`].
pub trait Surface: Sync {
    /// Whether grid point `(i, j, k)` is inside the surface.
    fn inside(&self, i: usize, j: usize, k: usize) -> bool;

    /// Where the surface crosses the grid edge from `(i, j, k)` toward
    /// `+axis` (0, 1, 2 for x, y, z), as a fraction of the voxel. Asked
    /// only of edges whose ends differ in [`Surface::inside`].
    fn crossing(&self, axis: usize, i: usize, j: usize, k: usize) -> f32;
}

/// The level set `values >= level` of a scalar field sampled at every grid
/// point, crossing by linear interpolation between grid points.
pub struct Isosurface<'a> {
    pub grid: &'a NetGrid,
    pub values: &'a [f32],
    pub level: f32,
}

impl Surface for Isosurface<'_> {
    fn inside(&self, i: usize, j: usize, k: usize) -> bool {
        self.values[self.grid.index(i, j, k)] >= self.level
    }

    fn crossing(&self, axis: usize, i: usize, j: usize, k: usize) -> f32 {
        let to = [[i + 1, j, k], [i, j + 1, k], [i, j, k + 1]][axis];
        let a = self.values[self.grid.index(i, j, k)];
        let b = self.values[self.grid.index(to[0], to[1], to[2])];
        (self.level - a) / (b - a)
    }
}

/// The cell vertices, each keyed by its cell's grid index and ascending
/// by it, and the triangles joining them.
pub struct Net {
    pub cells: Vec<(usize, Vec3)>,
    pub indices: Vec<u32>,
}

const CORNERS: [[usize; 3]; 8] = [
    [0, 0, 0],
    [1, 0, 0],
    [0, 1, 0],
    [1, 1, 0],
    [0, 0, 1],
    [1, 0, 1],
    [0, 1, 1],
    [1, 1, 1],
];

/// A cell's edges as corner pairs, with the axis each runs along.
const EDGES: [([usize; 2], usize); 12] = [
    ([0, 1], 0),
    ([2, 3], 0),
    ([4, 5], 0),
    ([6, 7], 0),
    ([0, 2], 1),
    ([1, 3], 1),
    ([4, 6], 1),
    ([5, 7], 1),
    ([0, 4], 2),
    ([1, 5], 2),
    ([2, 6], 2),
    ([3, 7], 2),
];

/// The net of `surface` on `grid`, triangles wound counter-clockwise seen
/// from outside.
pub fn build<S: Surface>(grid: &NetGrid, surface: &S) -> Net {
    let inside = inside_flags(grid, surface);
    let cells = surface_cells(grid, surface, &inside);
    let indices = faces(grid, &inside, &cells);
    Net { cells, indices }
}

fn inside_flags<S: Surface>(grid: &NetGrid, surface: &S) -> Vec<bool> {
    let [nx, ny, _] = grid.dims;
    let mut flags = vec![false; grid.len()];
    flags
        .par_chunks_mut(nx * ny)
        .enumerate()
        .for_each(|(k, plane)| {
            for j in 0..ny {
                for i in 0..nx {
                    plane[i + nx * j] = surface.inside(i, j, k);
                }
            }
        });
    flags
}

/// Every cell the surface crosses, in index order, with the mean of its
/// edges' crossings.
fn surface_cells<S: Surface>(grid: &NetGrid, surface: &S, inside: &[bool]) -> Vec<(usize, Vec3)> {
    let [nx, ny, nz] = grid.dims;
    (0..nz - 1)
        .into_par_iter()
        .flat_map_iter(|k| {
            let mut found = Vec::new();
            for j in 0..ny - 1 {
                for i in 0..nx - 1 {
                    let flags = CORNERS.map(|[a, b, c]| inside[grid.index(i + a, j + b, k + c)]);
                    if flags.iter().all(|&x| x) || flags.iter().all(|&x| !x) {
                        continue;
                    }
                    let (mut sum, mut n) = (Vec3::ZERO, 0.0);
                    for ([a, b], axis) in EDGES {
                        if flags[a] != flags[b] {
                            let [ca, cb, cc] = CORNERS[a];
                            let t = surface.crossing(axis, i + ca, j + cb, k + cc);
                            sum += Vec3::from(CORNERS[a].map(|c| c as f32)) + Vec3::AXES[axis] * t;
                            n += 1.0;
                        }
                    }
                    found.push((
                        grid.index(i, j, k),
                        grid.point(i, j, k) + sum / n * grid.voxel,
                    ));
                }
            }
            found
        })
        .collect()
}

/// Two triangles per grid edge the surface crosses, joining the four cells
/// around it, wound so their normal points from inside to outside.
fn faces(grid: &NetGrid, inside: &[bool], cells: &[(usize, Vec3)]) -> Vec<u32> {
    let [nx, ny, nz] = grid.dims;
    let vertex = |i: usize, j: usize, k: usize| {
        let key = grid.index(i, j, k);
        cells
            .binary_search_by_key(&key, |c| c.0)
            .ok()
            .map(|v| v as u32)
    };
    (1..nz - 1)
        .into_par_iter()
        .flat_map_iter(|k| {
            let mut tris = Vec::new();
            for j in 1..ny - 1 {
                for i in 1..nx - 1 {
                    let here = inside[grid.index(i, j, k)];
                    // The edge to the next point along each axis, and the
                    // other two axes' cells around it.
                    let around: [([usize; 3], [[usize; 3]; 4]); 3] = [
                        (
                            [i + 1, j, k],
                            [[i, j - 1, k - 1], [i, j, k - 1], [i, j, k], [i, j - 1, k]],
                        ),
                        (
                            [i, j + 1, k],
                            [[i - 1, j, k - 1], [i - 1, j, k], [i, j, k], [i, j, k - 1]],
                        ),
                        (
                            [i, j, k + 1],
                            [[i - 1, j - 1, k], [i, j - 1, k], [i, j, k], [i - 1, j, k]],
                        ),
                    ];
                    for (axis, ([a, b, c], quad)) in around.into_iter().enumerate() {
                        if inside[grid.index(a, b, c)] == here {
                            continue;
                        }
                        let v = quad.map(|[x, y, z]| vertex(x, y, z));
                        if let [Some(v0), Some(v1), Some(v2), Some(v3)] = v {
                            let outward = if here { 1.0 } else { -1.0 } * Vec3::AXES[axis];
                            let p = |v: u32| cells[v as usize].1;
                            let facing = (p(v1) - p(v0)).cross(p(v2) - p(v0)).dot(outward);
                            if facing >= 0.0 {
                                tris.extend([v0, v1, v2, v0, v2, v3]);
                            } else {
                                tris.extend([v0, v2, v1, v0, v3, v2]);
                            }
                        }
                    }
                }
            }
            tris
        })
        .collect()
}

/// `net` as a mesh, normals from its triangles (area-weighted), each vertex
/// attributed to `source(position)`.
pub fn mesh_of(net: Net, source: impl Fn(Vec3) -> u32 + Sync) -> ExpandedMesh {
    let positions: Vec<Vec3> = net.cells.iter().map(|c| c.1).collect();
    let normals = vertex_normals(&positions, &net.indices);
    let source = positions.par_iter().map(|&p| source(p)).collect();
    ExpandedMesh {
        positions,
        normals,
        indices: net.indices,
        source,
    }
}

fn vertex_normals(positions: &[Vec3], indices: &[u32]) -> Vec<Vec3> {
    let mut normals = vec![Vec3::ZERO; positions.len()];
    for t in indices.chunks_exact(3) {
        let [a, b, c] = [t[0], t[1], t[2]].map(|v| v as usize);
        let area_normal = (positions[b] - positions[a]).cross(positions[c] - positions[a]);
        for v in [a, b, c] {
            normals[v] += area_normal;
        }
    }
    normals
        .into_iter()
        .map(|n| n.normalize_or(Vec3::Z))
        .collect()
}

/// The atom whose surface (`|p - center| - radius`) is nearest `p`, found
/// among those within `reach` of it.
pub fn nearest_surface_atom(
    neighbors: &Grid,
    positions: &[Vec3],
    radii: &[f32],
    reach: f32,
    p: Vec3,
) -> u32 {
    let mut best = (f32::INFINITY, 0u32);
    neighbors.for_each_within(positions, p, reach, |a, d2| {
        let gap = d2.sqrt() - radii[a as usize];
        if gap < best.0 {
            best = (gap, a);
        }
    });
    best.1
}
