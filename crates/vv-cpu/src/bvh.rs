//! A bounding-volume hierarchy over bounding spheres, for nearest-hit
//! ray queries. Built once per structure and reused across every frame the
//! camera moves through (see `CpuScene`) — brute force over every atom per
//! pixel measured at 2 fps for a 4,779-atom structure at 512x512 (see
//! `tests/bench.rs`), 5x short of the >=10 fps target, so this exists.
//!
//! Median-split top-down build (no SAH — atoms are roughly uniform in
//! density within one structure, so a cheaper heuristic is a fine trade),
//! recursive nearest-hit traversal with AABB-distance pruning. Children
//! aren't visited in near-to-far order, which leaves some pruning on the
//! table; revisit if `tests/bench.rs`'s target ever regresses.
//!
//! The primitives need not be spheres: [`Bvh::nearest_hit_with`] takes
//! the real intersection test as a closure and only uses the sphere as a
//! bound (the skin surface's quadric patches go through it that way).

use glam::Vec3;

const LEAF_SIZE: usize = 4;

enum NodeKind {
    Leaf { start: u32, end: u32 },
    Interior { left: u32, right: u32 },
}

struct Node {
    min: Vec3,
    max: Vec3,
    kind: NodeKind,
}

pub struct Bvh {
    nodes: Vec<Node>,
    /// Primitive indices reordered by the build; leaves reference
    /// contiguous ranges into this rather than the original order.
    indices: Vec<u32>,
    root: u32,
}

fn bounds(positions: &[Vec3], radii: &[f32], indices: &[u32]) -> (Vec3, Vec3) {
    indices.iter().fold(
        (Vec3::splat(f32::INFINITY), Vec3::splat(f32::NEG_INFINITY)),
        |(lo, hi), &i| {
            let i = i as usize;
            let r = Vec3::splat(radii[i]);
            (lo.min(positions[i] - r), hi.max(positions[i] + r))
        },
    )
}

impl Bvh {
    pub fn build(positions: &[Vec3], radii: &[f32]) -> Self {
        let mut indices: Vec<u32> = (0..positions.len() as u32).collect();
        let mut nodes = Vec::new();
        let root = if indices.is_empty() {
            nodes.push(Node {
                min: Vec3::ZERO,
                max: Vec3::ZERO,
                kind: NodeKind::Leaf { start: 0, end: 0 },
            });
            0
        } else {
            build_range(&mut indices, 0, positions, radii, &mut nodes)
        };
        Bvh {
            nodes,
            indices,
            root,
        }
    }

    /// Nearest sphere the ray `origin + t*dir` (`dir` unit length) hits, if
    /// any, as `(t, index)`.
    pub fn nearest_hit(
        &self,
        origin: Vec3,
        dir: Vec3,
        positions: &[Vec3],
        radii: &[f32],
    ) -> Option<(f32, usize)> {
        self.nearest_hit_with(origin, dir, positions, radii, |i| {
            crate::renderer::intersect_sphere(origin, dir, positions[i], radii[i])
        })
    }

    /// Nearest primitive along the ray, where `hit(i)` is the real
    /// intersection test for primitive `i` (returning its ray parameter,
    /// or `None`) and `positions[i]`/`radii[i]` only bound it: `hit` is
    /// called for each primitive whose bounding sphere the ray enters
    /// before the best hit so far, in no particular order.
    pub fn nearest_hit_with(
        &self,
        origin: Vec3,
        dir: Vec3,
        positions: &[Vec3],
        radii: &[f32],
        mut hit: impl FnMut(usize) -> Option<f32>,
    ) -> Option<(f32, usize)> {
        let inv_dir = Vec3::new(1.0 / dir.x, 1.0 / dir.y, 1.0 / dir.z);
        let mut best_t = f32::INFINITY;
        let mut best_i = usize::MAX;
        self.traverse(
            self.root,
            origin,
            dir,
            inv_dir,
            positions,
            radii,
            &mut hit,
            &mut best_t,
            &mut best_i,
        );
        (best_i != usize::MAX).then_some((best_t, best_i))
    }

    #[allow(clippy::too_many_arguments)]
    fn traverse(
        &self,
        node: u32,
        origin: Vec3,
        dir: Vec3,
        inv_dir: Vec3,
        positions: &[Vec3],
        radii: &[f32],
        hit: &mut impl FnMut(usize) -> Option<f32>,
        best_t: &mut f32,
        best_i: &mut usize,
    ) {
        let n = &self.nodes[node as usize];
        match ray_aabb(origin, inv_dir, n.min, n.max) {
            Some(t) if t < *best_t => {}
            _ => return,
        }
        match n.kind {
            NodeKind::Leaf { start, end } => {
                for &i in &self.indices[start as usize..end as usize] {
                    let i = i as usize;
                    match sphere_entry(origin, dir, positions[i], radii[i]) {
                        Some(entry) if entry < *best_t => {}
                        _ => continue,
                    }
                    if let Some(t) = hit(i) {
                        if t < *best_t {
                            *best_t = t;
                            *best_i = i;
                        }
                    }
                }
            }
            NodeKind::Interior { left, right } => {
                self.traverse(
                    left, origin, dir, inv_dir, positions, radii, hit, best_t, best_i,
                );
                self.traverse(
                    right, origin, dir, inv_dir, positions, radii, hit, best_t, best_i,
                );
            }
        }
    }
}

/// Where the ray enters the sphere (clamped to 0 for a ray starting
/// inside it), or `None` if it never does -- a conservative lower bound
/// on any hit inside that sphere.
fn sphere_entry(origin: Vec3, dir: Vec3, center: Vec3, radius: f32) -> Option<f32> {
    let oc = origin - center;
    let b = oc.dot(dir);
    let c = oc.length_squared() - radius * radius;
    let disc = b * b - c;
    if disc < 0.0 {
        return None;
    }
    let sqrt_disc = disc.sqrt();
    let far = -b + sqrt_disc;
    (far > 0.0).then(|| (-b - sqrt_disc).max(0.0))
}

/// Slab test. `None` if the ray misses the box entirely or the box is
/// entirely behind the ray; otherwise the near distance (clamped to 0, so
/// a ray starting inside the box reports "already there").
fn ray_aabb(origin: Vec3, inv_dir: Vec3, min: Vec3, max: Vec3) -> Option<f32> {
    let t1 = (min - origin) * inv_dir;
    let t2 = (max - origin) * inv_dir;
    let tmin = t1.min(t2).max_element().max(0.0);
    let tmax = t1.max(t2).min_element();
    (tmax >= tmin).then_some(tmin)
}

/// `base` is `indices`' absolute offset within the top-level array being
/// partitioned in place (like quicksort): `indices` here is always a
/// sub-slice of that same backing array, so a leaf's `[start, end)` range
/// is only correct once expressed in absolute terms, not relative to
/// whatever local slice this call happened to receive.
fn build_range(
    indices: &mut [u32],
    base: u32,
    positions: &[Vec3],
    radii: &[f32],
    nodes: &mut Vec<Node>,
) -> u32 {
    let (min, max) = bounds(positions, radii, indices);
    if indices.len() <= LEAF_SIZE {
        nodes.push(Node {
            min,
            max,
            kind: NodeKind::Leaf {
                start: base,
                end: base + indices.len() as u32,
            },
        });
        return (nodes.len() - 1) as u32;
    }

    let extent = max - min;
    let axis = if extent.x >= extent.y && extent.x >= extent.z {
        0
    } else if extent.y >= extent.z {
        1
    } else {
        2
    };
    indices.sort_unstable_by(|&a, &b| {
        positions[a as usize][axis].total_cmp(&positions[b as usize][axis])
    });
    let mid = indices.len() / 2;
    let (left_indices, right_indices) = indices.split_at_mut(mid);
    let left = build_range(left_indices, base, positions, radii, nodes);
    let right = build_range(right_indices, base + mid as u32, positions, radii, nodes);
    nodes.push(Node {
        min,
        max,
        kind: NodeKind::Interior { left, right },
    });
    (nodes.len() - 1) as u32
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn finds_the_nearest_of_several_colinear_spheres() {
        let positions = vec![
            Vec3::new(0.0, 0.0, -5.0),
            Vec3::new(0.0, 0.0, -10.0),
            Vec3::new(0.0, 0.0, -15.0),
        ];
        let radii = vec![1.0, 1.0, 1.0];
        let bvh = Bvh::build(&positions, &radii);
        let (t, i) = bvh
            .nearest_hit(Vec3::ZERO, Vec3::new(0.0, 0.0, -1.0), &positions, &radii)
            .expect("should hit the nearest sphere");
        assert_eq!(i, 0, "should hit the closest sphere, not a farther one");
        assert!((t - 4.0).abs() < 1e-4);
    }

    #[test]
    fn misses_when_nothing_is_along_the_ray() {
        let positions = vec![Vec3::new(10.0, 10.0, 10.0)];
        let radii = vec![1.0];
        let bvh = Bvh::build(&positions, &radii);
        assert!(bvh
            .nearest_hit(Vec3::ZERO, Vec3::new(0.0, 0.0, -1.0), &positions, &radii)
            .is_none());
    }

    #[test]
    fn empty_scene_never_hits() {
        let bvh = Bvh::build(&[], &[]);
        assert!(bvh
            .nearest_hit(Vec3::ZERO, Vec3::new(0.0, 0.0, -1.0), &[], &[])
            .is_none());
    }

    #[test]
    fn a_custom_hit_test_can_reject_the_nearest_bound_and_take_the_next() {
        // Two spheres along the ray; the closer one's real test says
        // "miss", so the farther one must win even though its bound is
        // entered later.
        let positions = vec![Vec3::new(0.0, 0.0, -5.0), Vec3::new(0.0, 0.0, -10.0)];
        let radii = vec![1.0, 1.0];
        let bvh = Bvh::build(&positions, &radii);
        let dir = Vec3::new(0.0, 0.0, -1.0);
        let (t, i) = bvh
            .nearest_hit_with(Vec3::ZERO, dir, &positions, &radii, |i| {
                (i == 1)
                    .then(|| {
                        crate::renderer::intersect_sphere(Vec3::ZERO, dir, positions[i], radii[i])
                    })
                    .flatten()
            })
            .unwrap();
        assert_eq!(i, 1);
        assert!((t - 9.0).abs() < 1e-4);
    }

    #[test]
    fn agrees_with_brute_force_on_a_random_scene() {
        // Deterministic pseudo-random cloud of overlapping spheres,
        // checked against an O(n) linear scan for many ray directions.
        let mut state = 0x1234_5678u32;
        let mut next = move || {
            state ^= state << 13;
            state ^= state >> 17;
            state ^= state << 5;
            state
        };
        let n = 300;
        let positions: Vec<Vec3> = (0..n)
            .map(|_| {
                let f = |bits: u32| (bits % 2000) as f32 / 100.0 - 10.0;
                Vec3::new(f(next()), f(next()), f(next()))
            })
            .collect();
        let radii: Vec<f32> = (0..n)
            .map(|_| (next() % 100) as f32 / 100.0 + 0.3)
            .collect();
        let bvh = Bvh::build(&positions, &radii);

        for _ in 0..200 {
            let dir = Vec3::new(
                (next() % 200) as f32 / 100.0 - 1.0,
                (next() % 200) as f32 / 100.0 - 1.0,
                (next() % 200) as f32 / 100.0 - 1.0,
            )
            .normalize_or(Vec3::Z);
            let origin = Vec3::new(0.0, 0.0, 30.0);

            let bvh_hit = bvh.nearest_hit(origin, dir, &positions, &radii);
            let brute_force = (0..n)
                .filter_map(|i| {
                    crate::renderer::intersect_sphere(origin, dir, positions[i], radii[i])
                        .map(|t| (t, i))
                })
                .min_by(|a, b| a.0.total_cmp(&b.0));

            match (bvh_hit, brute_force) {
                (None, None) => {}
                (Some((t_bvh, _)), Some((t_brute, _))) => {
                    assert!(
                        (t_bvh - t_brute).abs() < 1e-3,
                        "BVH and brute force disagree on hit distance: {t_bvh} vs {t_brute}"
                    );
                }
                other => panic!("BVH and brute force disagree on hit/miss: {other:?}"),
            }
        }
    }
}
