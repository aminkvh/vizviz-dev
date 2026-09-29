//! Uniform grid for fixed-radius neighbour queries ("is any stored atom
//! within `r` of this point?").
//!
//! Built over a *subset* of atoms (the inner set of `within`), queried with
//! arbitrary points. Cells tile the subset's bounding box and are filled
//! with a counting sort, so each cell is a contiguous slice of atom indices.
//! The cell edge is at least the query radius, so a query only ever touches
//! the 3x3x3 block of cells around its own.

use glam::{IVec3, Vec3};

/// Upper bound on cells per axis. A tiny cutoff over a huge set (`within
/// 0.5 of all` on a 10M-atom box) would otherwise allocate a grid of
/// hundreds of millions of cells; past this bound cells grow instead.
const MAX_CELLS_PER_AXIS: f32 = 128.0;

#[derive(Clone, Debug)]
pub struct Grid {
    origin: Vec3,
    cell: f32,
    dims: IVec3,
    /// `starts[c]..starts[c + 1]` is the slice of `atoms` in cell `c`.
    starts: Vec<u32>,
    atoms: Vec<u32>,
}

/// Atoms from which [`Grid::build`] sorts in parallel.
const PARALLEL_BUILD: usize = 1 << 16;

impl Grid {
    /// Bins the atoms listed in `indices` into cells of edge `cell`
    /// (larger if needed to respect [`MAX_CELLS_PER_AXIS`]). `cell` must be
    /// positive and finite.
    pub fn build(positions: &[Vec3], indices: &[u32], cell: f32) -> Grid {
        assert!(cell > 0.0 && cell.is_finite(), "grid cell must be positive");
        let Some(&first) = indices.first() else {
            return Grid {
                origin: Vec3::ZERO,
                cell,
                dims: IVec3::ZERO,
                starts: vec![0],
                atoms: Vec::new(),
            };
        };
        let seed = positions[first as usize];
        let (min, max) = indices.iter().fold((seed, seed), |(lo, hi), &i| {
            let p = positions[i as usize];
            (lo.min(p), hi.max(p))
        });
        let cell = cell.max((max - min).max_element() / MAX_CELLS_PER_AXIS);
        let dims = ((max - min) / cell).floor().as_ivec3() + IVec3::ONE;
        let cell_count = (dims.x * dims.y * dims.z) as usize;

        let mut grid = Grid {
            origin: min,
            cell,
            dims,
            starts: vec![0; cell_count + 1],
            atoms: vec![0; indices.len()],
        };
        if indices.len() >= PARALLEL_BUILD {
            grid.fill_parallel(positions, indices);
            return grid;
        }
        // Counting sort: histogram, prefix sum, scatter.
        let cells: Vec<usize> = indices
            .iter()
            .map(|&i| grid.cell_index(grid.cell_of(positions[i as usize])))
            .collect();
        for &c in &cells {
            grid.starts[c + 1] += 1;
        }
        for c in 0..cell_count {
            grid.starts[c + 1] += grid.starts[c];
        }
        let mut cursor = grid.starts.clone();
        for (&atom, &c) in indices.iter().zip(&cells) {
            grid.atoms[cursor[c] as usize] = atom;
            cursor[c] += 1;
        }
        grid
    }

    /// The counting sort's result, by a parallel sort of (cell, input
    /// position) keys: the same cells, each in input order. Its scatter
    /// is random access, which dominates at millions of atoms.
    fn fill_parallel(&mut self, positions: &[Vec3], indices: &[u32]) {
        use rayon::prelude::*;
        let mut keys: Vec<u64> = indices
            .par_iter()
            .enumerate()
            .map(|(k, &i)| {
                let cell = self.cell_index(self.cell_of(positions[i as usize])) as u64;
                (cell << 32) | k as u64
            })
            .collect();
        keys.par_sort_unstable();
        self.atoms
            .par_iter_mut()
            .zip(&keys)
            .for_each(|(atom, &key)| *atom = indices[(key & 0xffff_ffff) as usize]);
        // Cell c starts at the first key in a cell >= c (empty cells
        // start where the next non-empty one does).
        let mut k = keys.len();
        for c in (0..self.starts.len()).rev() {
            while k > 0 && (keys[k - 1] >> 32) as usize >= c {
                k -= 1;
            }
            self.starts[c] = k as u32;
        }
    }

    /// Number of stored atoms.
    pub fn len(&self) -> usize {
        self.atoms.len()
    }

    pub fn is_empty(&self) -> bool {
        self.atoms.is_empty()
    }

    /// The box corner cell `(0, 0, 0)` starts at.
    pub fn origin(&self) -> Vec3 {
        self.origin
    }

    /// The cell edge actually used (at least what [`Grid::build`] was
    /// asked for; larger when [`MAX_CELLS_PER_AXIS`] forced it).
    pub fn cell_size(&self) -> f32 {
        self.cell
    }

    pub fn dims(&self) -> IVec3 {
        self.dims
    }

    /// The CSR layout, for uploading the grid as-is to a GPU: cell `c`
    /// (linear index `(z * dims.y + y) * dims.x + x`) holds
    /// `cell_atoms()[cell_starts()[c]..cell_starts()[c + 1]]`. A query
    /// point's cell is `floor((p - origin) / cell_size)`, which may fall
    /// outside `0..dims` for points outside the box; a 3x3x3 block
    /// around it, skipping out-of-range cells, is what
    /// [`Grid::for_each_within`] scans.
    pub fn cell_starts(&self) -> &[u32] {
        &self.starts
    }

    pub fn cell_atoms(&self) -> &[u32] {
        &self.atoms
    }

    /// The `(2k+1)^3` block of cells around `p`, clipped to the grid, as
    /// an inclusive `(lo, hi)` pair; `None` when `p` is more than `k`
    /// cells away from the box (so nothing stored can be within `k` cell
    /// edges of it).
    fn block_around(&self, p: Vec3, k: f32) -> Option<(IVec3, IVec3)> {
        let c = ((p - self.origin) / self.cell).floor();
        let lo = (c - k).max(Vec3::ZERO);
        let hi = (c + k).min((self.dims - IVec3::ONE).as_vec3());
        if lo.cmpgt(hi).any() {
            return None;
        }
        Some((lo.as_ivec3(), hi.as_ivec3()))
    }

    /// How many rings of cells a query of `radius` needs on each side:
    /// one when the radius fits in a cell edge (the common case), more
    /// for a larger radius -- so any radius is answered exactly, just at
    /// the cost of scanning a bigger block.
    fn rings(&self, radius: f32) -> f32 {
        (radius / self.cell).ceil().max(1.0)
    }

    /// True if any stored atom lies within `radius` (inclusive) of `p`.
    pub fn any_within(&self, positions: &[Vec3], p: Vec3, radius: f32) -> bool {
        let r2 = radius * radius;
        let Some((lo, hi)) = self.block_around(p, self.rings(radius)) else {
            return false;
        };
        for z in lo.z..=hi.z {
            for y in lo.y..=hi.y {
                for x in lo.x..=hi.x {
                    let cell = self.cell_index(IVec3::new(x, y, z));
                    let hit = self
                        .atoms_in(cell)
                        .iter()
                        .any(|&a| positions[a as usize].distance_squared(p) <= r2);
                    if hit {
                        return true;
                    }
                }
            }
        }
        false
    }

    /// Calls `f(atom, distance_squared)` for every stored atom within
    /// `radius` (inclusive) of `p`, in no particular order. The workhorse
    /// behind contacts and neighbour lists. Cheapest when `radius` is at
    /// most the cell edge given to [`Grid::build`] (a 3x3x3 block);
    /// larger radii scan proportionally bigger blocks.
    pub fn for_each_within(
        &self,
        positions: &[Vec3],
        p: Vec3,
        radius: f32,
        mut f: impl FnMut(u32, f32),
    ) {
        let r2 = radius * radius;
        let Some((lo, hi)) = self.block_around(p, self.rings(radius)) else {
            return;
        };
        for z in lo.z..=hi.z {
            for y in lo.y..=hi.y {
                for x in lo.x..=hi.x {
                    let cell = self.cell_index(IVec3::new(x, y, z));
                    for &a in self.atoms_in(cell) {
                        let d2 = positions[a as usize].distance_squared(p);
                        if d2 <= r2 {
                            f(a, d2);
                        }
                    }
                }
            }
        }
    }

    /// Cell coordinates of a point inside the bounding box. Clamped so that
    /// float rounding on the far face can never index past the last cell.
    fn cell_of(&self, p: Vec3) -> IVec3 {
        ((p - self.origin) / self.cell)
            .floor()
            .as_ivec3()
            .clamp(IVec3::ZERO, self.dims - IVec3::ONE)
    }

    fn cell_index(&self, c: IVec3) -> usize {
        ((c.z * self.dims.y + c.y) * self.dims.x + c.x) as usize
    }

    fn atoms_in(&self, cell: usize) -> &[u32] {
        &self.atoms[self.starts[cell] as usize..self.starts[cell + 1] as usize]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_parallel_build_bins_exactly_as_the_counting_sort() {
        let points = cloud(3 * PARALLEL_BUILD, 120.0, 11);
        // Every other atom, reversed: the input order must survive.
        let indices: Vec<u32> = (0..points.len() as u32).rev().step_by(2).collect();
        let big = Grid::build(&points, &indices, 3.0);
        assert!(indices.len() >= PARALLEL_BUILD);
        let mut cells: Vec<Vec<u32>> = vec![Vec::new(); big.starts.len() - 1];
        for &i in &indices {
            cells[big.cell_index(big.cell_of(points[i as usize]))].push(i);
        }
        let mut starts = vec![0u32];
        for c in &cells {
            starts.push(starts.last().unwrap() + c.len() as u32);
        }
        assert_eq!(big.starts, starts);
        assert_eq!(big.atoms, cells.concat());
    }

    /// Deterministic pseudo-random points in a `side`-Å cube (LCG, no
    /// external crate).
    fn cloud(n: usize, side: f32, seed: u64) -> Vec<Vec3> {
        let mut state = seed;
        let mut next = move || {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (state >> 40) as f32 / (1u64 << 24) as f32 * side
        };
        (0..n).map(|_| Vec3::new(next(), next(), next())).collect()
    }

    fn brute(positions: &[Vec3], indices: &[u32], p: Vec3, r: f32) -> bool {
        indices
            .iter()
            .any(|&i| positions[i as usize].distance_squared(p) <= r * r)
    }

    #[test]
    fn agrees_with_brute_force_on_a_random_cloud() {
        let positions = cloud(3000, 40.0, 7);
        // Store every third atom; query with all of them plus points
        // outside the box.
        let indices: Vec<u32> = (0..positions.len() as u32).step_by(3).collect();
        let outside = cloud(200, 60.0, 99)
            .into_iter()
            .map(|p| p - Vec3::splat(10.0))
            .collect::<Vec<_>>();
        let grid = Grid::build(&positions, &indices, 3.0);
        assert_eq!(grid.len(), indices.len());
        for &r in &[0.5f32, 1.7, 3.0] {
            let mut hits = 0;
            for &p in positions.iter().chain(&outside) {
                let expected = brute(&positions, &indices, p, r);
                assert_eq!(grid.any_within(&positions, p, r), expected, "r={r} p={p}");
                hits += usize::from(expected);
            }
            assert!(hits > 0, "test should exercise positive cases at r={r}");
        }
    }

    #[test]
    fn a_radius_larger_than_the_cell_still_answers_exactly() {
        let positions = cloud(2000, 30.0, 11);
        let indices: Vec<u32> = (0..positions.len() as u32).collect();
        let grid = Grid::build(&positions, &indices, 2.0);
        for &r in &[2.0f32, 3.5, 7.0] {
            for &p in positions.iter().step_by(37) {
                let mut via_grid: Vec<u32> = Vec::new();
                grid.for_each_within(&positions, p, r, |a, _| via_grid.push(a));
                via_grid.sort_unstable();
                let mut brute: Vec<u32> = indices
                    .iter()
                    .copied()
                    .filter(|&i| positions[i as usize].distance_squared(p) <= r * r)
                    .collect();
                brute.sort_unstable();
                assert_eq!(via_grid, brute, "r={r} p={p}");
            }
        }
    }

    #[test]
    fn radius_is_inclusive_and_crosses_cell_borders() {
        let positions = [Vec3::ZERO, Vec3::new(3.0, 4.0, 0.0)];
        let grid = Grid::build(&positions, &[1], 5.0);
        assert!(grid.any_within(&positions, Vec3::ZERO, 5.0));
        assert!(!grid.any_within(&positions, Vec3::ZERO, 4.99));
        // The query point sits well outside the one-cell box.
        assert!(grid.any_within(&positions, Vec3::new(3.0, 4.0, -4.9), 5.0));
        assert!(!grid.any_within(&positions, Vec3::new(3.0, 4.0, -5.1), 5.0));
        assert!(!grid.any_within(&positions, Vec3::new(1e6, 0.0, 0.0), 5.0));
    }

    #[test]
    fn empty_grid_matches_nothing() {
        let positions = [Vec3::ZERO];
        let grid = Grid::build(&positions, &[], 2.0);
        assert!(grid.is_empty());
        assert!(!grid.any_within(&positions, Vec3::ZERO, 2.0));
    }

    #[test]
    fn tiny_cell_over_a_large_box_stays_bounded() {
        // 0.01 Å cells over 1000 Å would be 10^15 cells; the cap makes this
        // finish instantly and still answer correctly.
        let positions = [Vec3::ZERO, Vec3::splat(1000.0), Vec3::new(500.0, 0.0, 0.0)];
        let grid = Grid::build(&positions, &[0, 1, 2], 0.01);
        assert!(grid.any_within(&positions, Vec3::new(500.005, 0.0, 0.0), 0.01));
        assert!(!grid.any_within(&positions, Vec3::new(500.05, 0.0, 0.0), 0.01));
    }
}
