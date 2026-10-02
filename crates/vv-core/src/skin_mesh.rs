//! The skin surface as triangles, for the mesh display. The surface is the
//! union of quadric patches, each valid only inside its own mixed cell
//! ([`crate::skin_surface`]), so there is no single field to sample.
//! Instead every grid line along x, y and z is cast through the patches near
//! it, using the exact ray roots and membership test the ray caster uses:
//! the sorted hits of a line say which of its grid points are inside (by
//! parity) and exactly where each grid edge is crossed, which is what
//! surface nets ([`crate::surface_net`]) needs. The vertices therefore lie
//! on the surface, not a voxel off it.

use glam::Vec3;
use rayon::prelude::*;

use crate::cartoon::ExpandedMesh;
use crate::skin_surface::{build_complex, weight_for_radius, SkinComplex};
use crate::spatial::Grid;
use crate::surface_net::{self, NetGrid, Surface};

/// Voxels a mesh may sample, at most.
pub const MAX_VOXELS: usize = 1 << 25;

/// Hits closer than this along a line are one hit seen by two patches
/// that share a boundary.
const SAME_HIT: f32 = 1e-3;

/// The skin surface of atoms at `positions` with van der Waals `radii` and
/// `shrink`, on a grid of `voxel` angstroms (raised as needed to stay within
/// [`MAX_VOXELS`]). Each vertex's `source` is the atom whose surface is
/// nearest it.
pub fn mesh(positions: &[Vec3], radii: &[f32], shrink: f32, voxel: f32) -> ExpandedMesh {
    assert_eq!(positions.len(), radii.len());
    let weights: Vec<f32> = radii
        .iter()
        .map(|&r| weight_for_radius(r, shrink))
        .collect();
    let complex = build_complex(positions, &weights, shrink);
    let Some((lo, hi)) = surface_bounds(&complex) else {
        return ExpandedMesh::default();
    };
    let grid = NetGrid::around(lo, hi, voxel, MAX_VOXELS);
    let lines = Lines::cast(&grid, &complex, positions, &weights);
    let net = surface_net::build(&grid, &lines.surface(&grid));
    let reach = 2.0 * radii.iter().copied().fold(0.0, f32::max) + 1.0;
    let all: Vec<u32> = (0..positions.len() as u32).collect();
    let neighbors = Grid::build(positions, &all, reach);
    surface_net::mesh_of(net, |p| {
        surface_net::nearest_surface_atom(&neighbors, positions, radii, reach, p)
    })
}

/// A box around every patch that has a surface.
fn surface_bounds(complex: &SkinComplex) -> Option<(Vec3, Vec3)> {
    complex
        .patches
        .iter()
        .filter(|p| p.has_surface())
        .map(|p| {
            (
                p.bound_center - Vec3::splat(p.bound_radius),
                p.bound_center + Vec3::splat(p.bound_radius),
            )
        })
        .reduce(|(lo, hi), (a, b)| (lo.min(a), hi.max(b)))
}

/// The other two axes, in index order, of lines along `axis`.
fn across(axis: usize) -> [usize; 2] {
    [(axis + 1) % 3, (axis + 2) % 3]
}

/// Where every grid line along each axis meets the surface: for axis `w`,
/// line `iu + dims[u] * iv` (`[u, v] = across(w)`) holds its hits' distances
/// from the grid's origin along `w`, ascending.
struct Lines {
    along: [Vec<Vec<f32>>; 3],
}

impl Lines {
    fn cast(grid: &NetGrid, complex: &SkinComplex, positions: &[Vec3], weights: &[f32]) -> Self {
        let along = [0, 1, 2].map(|axis| {
            let [u, v] = across(axis);
            let (du, dv) = (grid.dims[u], grid.dims[v]);
            let hits: Vec<(usize, f32)> = (0..complex.patches.len())
                .into_par_iter()
                .filter(|&p| complex.patches[p].has_surface())
                .flat_map_iter(|p| patch_hits(grid, complex, positions, weights, axis, p))
                .collect();
            let mut lines = vec![Vec::new(); du * dv];
            for (line, t) in hits {
                lines[line].push(t);
            }
            lines.par_iter_mut().for_each(settle);
            lines
        });
        Self { along }
    }

    fn surface(&self, grid: &NetGrid) -> HitSurface<'_> {
        let [nx, ny, nz] = grid.dims;
        let mut inside = vec![false; grid.len()];
        for j in 0..ny {
            for i in 0..nx {
                let hits = &self.along[2][i + nx * j];
                for k in 0..nz {
                    let below = hits.partition_point(|&t| t < k as f32 * grid.voxel);
                    inside[grid.index(i, j, k)] = below % 2 == 1;
                }
            }
        }
        HitSurface {
            grid: *grid,
            lines: &self.along,
            inside,
        }
    }
}

/// Every hit of the grid lines along `axis` with patch `p`, as (line, distance).
fn patch_hits(
    grid: &NetGrid,
    complex: &SkinComplex,
    positions: &[Vec3],
    weights: &[f32],
    axis: usize,
    p: usize,
) -> Vec<(usize, f32)> {
    let patch = &complex.patches[p];
    let [u, v] = across(axis);
    let (c, r) = (patch.bound_center, patch.bound_radius);
    let span = |a: usize| {
        let lo = ((c[a] - r - grid.origin[a]) / grid.voxel).ceil().max(0.0) as usize;
        let hi = (((c[a] + r - grid.origin[a]) / grid.voxel).floor().max(-1.0) + 1.0) as usize;
        lo..hi.min(grid.dims[a])
    };
    let dir = Vec3::AXES[axis];
    let mut hits = Vec::new();
    for iv in span(v) {
        for iu in span(u) {
            let mut start = grid.origin;
            start[u] += iu as f32 * grid.voxel;
            start[v] += iv as f32 * grid.voxel;
            let off = Vec3::new(start[u] - c[u], start[v] - c[v], 0.0).length();
            if off > r {
                continue;
            }
            for t in patch
                .ray_roots(start, dir, complex.shrink)
                .into_iter()
                .flatten()
            {
                if complex.contains(p, start + dir * t, positions, weights) {
                    hits.push((iu + grid.dims[u] * iv, t));
                }
            }
        }
    }
    hits
}

/// Sorts a line's hits, merges the doubled ones, and drops the line if it
/// still crosses an odd number of times (a hit lost at a seam): a surface
/// is closed, so a line enters it as often as it leaves.
fn settle(hits: &mut Vec<f32>) {
    hits.sort_by(f32::total_cmp);
    hits.dedup_by(|b, a| (*b - *a).abs() < SAME_HIT);
    if hits.len() % 2 == 1 {
        hits.clear();
    }
}

/// The surface as the hits of its grid lines.
struct HitSurface<'a> {
    grid: NetGrid,
    lines: &'a [Vec<Vec<f32>>; 3],
    inside: Vec<bool>,
}

impl Surface for HitSurface<'_> {
    fn inside(&self, i: usize, j: usize, k: usize) -> bool {
        self.inside[self.grid.index(i, j, k)]
    }

    fn crossing(&self, axis: usize, i: usize, j: usize, k: usize) -> f32 {
        let at = [i, j, k];
        let [u, v] = across(axis);
        let line = &self.lines[axis][at[u] + self.grid.dims[u] * at[v]];
        let from = at[axis] as f32 * self.grid.voxel;
        let hit = line.partition_point(|&t| t < from);
        match line.get(hit) {
            Some(&t) if t < from + self.grid.voxel => (t - from) / self.grid.voxel,
            _ => 0.5,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_lone_atoms_skin_is_its_van_der_waals_sphere_for_any_shrink() {
        let center = Vec3::new(0.5, -1.0, 2.0);
        for shrink in [0.3, 0.5, 0.8] {
            let mesh = mesh(&[center], &[1.7], shrink, 0.3);
            assert!(mesh.indices.len() > 300, "{}", mesh.indices.len());
            for (p, n) in mesh.positions.iter().zip(&mesh.normals) {
                assert!(((*p - center).length() - 1.7).abs() < 0.05, "{shrink} {p}");
                assert!(n.dot((*p - center).normalize()) > 0.95);
            }
        }
    }

    #[test]
    fn two_bonded_atoms_give_one_closed_outward_mesh() {
        let positions = [Vec3::ZERO, Vec3::new(1.5, 0.0, 0.0)];
        let mesh = mesh(&positions, &[1.7, 1.5], 0.5, 0.25);
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
        assert!(volume > 20.0, "{volume}");
    }
}
