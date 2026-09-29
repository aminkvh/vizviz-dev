//! Gaussian ("blobby") molecular surface density field: Blinn's original
//! blobby model (J.F. Blinn, "A Generalization of Algebraic Surface
//! Drawing," ACM TOG 1(3):235-256, 1982, doi:10.1145/357306.357310),
//! the same per-atom falloff Krone, Stone, Ertl & Schulten's fast
//! Gaussian density surfaces use (EuroVis 2012 Short Papers, doi:10.2312/PE/EuroVisShort/EuroVisShort2012/067-071).
//!
//! Chosen as the first native surface kernel specifically because it needs
//! **no triangulation at all**, weighted or otherwise -- a smooth,
//! everywhere-differentiable scalar field is evaluated directly, so
//! there is no weighted-Delaunay/Voronoi blocker, no CGAL/geogram
//! license question, and no combinatorial-geometry correctness risk to
//! get wrong. It sidesteps the entire class of problem SES/SAS and the
//! skin surface both have.
//!
//! Both [`density`] and [`gradient`] are closed-form (exact derivatives,
//! not finite-differenced), which is what makes this kernel ray-castable
//! without meshing: a renderer can sphere-trace the isosurface directly
//! (`density(p) == isovalue`) and get an exact analytic normal from
//! `gradient` at the hit point, matching this project's "no meshing,
//! ray-cast the real surface" choice already made for spheres, cylinders
//! and tubes. `vv_cpu::gaussian_surface` does exactly that.

use glam::Vec3;

/// One atom's contribution to the field: center and radius (van der
/// Waals radius for a classic "soft VDW" surface; any positive radius
/// works, e.g. a probe-inflated radius for a smoother blend).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Blob {
    pub center: Vec3,
    pub radius: f32,
}

/// Default steepness (after Krone et al. 2012): steep enough that each atom's
/// field is negligible a couple of radii out (see [`cutoff_radius`]), so
/// a renderer can safely ignore distant atoms.
pub const DEFAULT_BLOB_FACTOR: f32 = 2.0;

/// One atom's Gaussian contribution at squared distance `d2` from its
/// center, with its own `radius`. Chosen so the value is exactly `1.0`
/// at `d2 == radius^2` (the atom's own van der Waals surface) regardless
/// of `blob_factor` -- the isovalue for a "hugs the VDW surface" render
/// is therefore always `1.0`, independent of how steep the blend is.
fn contribution(d2: f32, radius: f32, blob_factor: f32) -> f32 {
    (-blob_factor * (d2 / (radius * radius) - 1.0)).exp()
}

/// The field value at `point`: the sum of every blob's Gaussian
/// contribution. Isosurfaces of this field are the Gaussian ("blobby")
/// surface; nearby atoms blend smoothly where their fields overlap
/// (Blinn's whole point -- no explicit union/blend logic needed, it
/// falls out of summing continuous fields).
///
/// Generic over any `&Blob` iterator (a plain slice still works
/// unchanged) rather than requiring a `&[Blob]`, so a spatial index like
/// `vv_cpu`'s blob grid can pass a filtered, zero-allocation iterator of
/// just the nearby blobs instead of collecting them into a new slice
/// first.
pub fn density<'a>(
    point: Vec3,
    blobs: impl IntoIterator<Item = &'a Blob>,
    blob_factor: f32,
) -> f32 {
    blobs
        .into_iter()
        .map(|b| contribution((point - b.center).length_squared(), b.radius, blob_factor))
        .sum()
}

/// The exact gradient of [`density`] at `point` (sum of each blob's own
/// closed-form derivative -- linearity of differentiation, not a finite
/// difference). Points toward *increasing* density, i.e. toward the
/// nearest mass concentration; a renderer wants the outward surface
/// normal, which is `-gradient(...).normalize()`. Generic over any
/// `&Blob` iterator for the same reason [`density`] is.
pub fn gradient<'a>(
    point: Vec3,
    blobs: impl IntoIterator<Item = &'a Blob>,
    blob_factor: f32,
) -> Vec3 {
    blobs.into_iter().fold(Vec3::ZERO, |acc, b| {
        let delta = point - b.center;
        let inv_r2 = 1.0 / (b.radius * b.radius);
        let g = contribution(delta.length_squared(), b.radius, blob_factor);
        acc + delta * (-2.0 * blob_factor * inv_r2 * g)
    })
}

/// A conservative radius beyond which one isolated blob's own
/// contribution drops below `epsilon` -- solves `epsilon =
/// exp(-blob_factor * (d^2/r^2 - 1))` for `d`. Useful for a renderer or
/// spatial index deciding which atoms can be ignored near a given point.
/// Cell edge for a spatial grid over atoms of these `radii`: the largest
/// [`cutoff_radius`], so that one cell-edge-radius query around any point
/// reaches every atom able to contribute more than `epsilon` there (an
/// atom farther away than the cell edge is farther than its own cutoff).
/// Shared by the CPU (`vv_cpu::gaussian_surface`) and GPU
/// (`vv_render::scene::GaussianSurfaceGpu`) grids so the two agree.
pub fn grid_cell_size(radii: &[f32], blob_factor: f32, epsilon: f32) -> f32 {
    radii
        .iter()
        .map(|&r| cutoff_radius(r, blob_factor, epsilon))
        .fold(0.0f32, f32::max)
        .max(1e-3)
}

pub fn cutoff_radius(radius: f32, blob_factor: f32, epsilon: f32) -> f32 {
    debug_assert!(epsilon > 0.0 && epsilon < 1.0);
    let ratio_sq = 1.0 - epsilon.ln() / blob_factor;
    radius * ratio_sq.max(0.0).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn blob(center: Vec3, radius: f32) -> Blob {
        Blob { center, radius }
    }

    #[test]
    fn an_isolated_atoms_own_surface_is_exactly_the_isovalue_one() {
        // contribution() is constructed so d2 == radius^2 gives exactly
        // exp(0) == 1.0, for ANY blob_factor -- this is the load-bearing
        // property that lets a renderer use a single fixed isovalue.
        for blob_factor in [0.5, DEFAULT_BLOB_FACTOR, 5.0] {
            let b = blob(Vec3::new(1.0, 2.0, 3.0), 1.7);
            let on_surface = b.center + Vec3::X * b.radius;
            let d = density(on_surface, &[b], blob_factor);
            assert!((d - 1.0).abs() < 1e-5, "blob_factor {blob_factor}: {d}");
        }
    }

    #[test]
    fn density_decays_monotonically_away_from_a_single_atom() {
        let b = blob(Vec3::ZERO, 1.5);
        let mut last = f32::INFINITY;
        for i in 1..20 {
            let d = density(Vec3::X * (i as f32 * 0.3), &[b], DEFAULT_BLOB_FACTOR);
            assert!(d < last, "density should strictly decrease with distance");
            last = d;
        }
        // And it grows without bound approaching the center (not clamped).
        let center_density = density(Vec3::ZERO, &[b], DEFAULT_BLOB_FACTOR);
        assert!(center_density > 1.0);
    }

    #[test]
    fn two_overlapping_atoms_blend_above_either_ones_own_contribution() {
        // The whole point of a blobby/Gaussian surface: fields superpose,
        // so the midpoint between two overlapping atoms reads higher than
        // either atom's field would alone -- a real "merge," not a union
        // of two separate spheres.
        let a = blob(Vec3::new(-1.0, 0.0, 0.0), 1.5);
        let b = blob(Vec3::new(1.0, 0.0, 0.0), 1.5);
        let midpoint = Vec3::ZERO;
        let combined = density(midpoint, &[a, b], DEFAULT_BLOB_FACTOR);
        let alone = density(midpoint, &[a], DEFAULT_BLOB_FACTOR);
        assert!(
            combined > alone * 1.9,
            "combined {combined} vs alone {alone}"
        );
    }

    #[test]
    fn gradient_points_toward_the_only_atom() {
        let b = blob(Vec3::new(2.0, 0.0, 0.0), 1.0);
        let probe = Vec3::new(4.0, 0.0, 0.0);
        let g = gradient(probe, &[b], DEFAULT_BLOB_FACTOR);
        // Vector from probe to the atom center.
        let toward_atom = b.center - probe;
        assert!(
            g.dot(toward_atom) > 0.0,
            "gradient {g:?} should point toward the atom, not away"
        );
    }

    #[test]
    fn gradient_matches_a_finite_difference() {
        // Cross-check the closed-form derivative against a numeric one,
        // independent evidence the formula's algebra is right.
        let blobs = [
            blob(Vec3::new(0.3, -0.2, 0.1), 1.4),
            blob(Vec3::new(-0.6, 0.5, 0.2), 1.1),
        ];
        let p = Vec3::new(0.4, 0.1, -0.3);
        let h = 1e-3;
        let numeric = Vec3::new(
            (density(p + Vec3::X * h, &blobs, DEFAULT_BLOB_FACTOR)
                - density(p - Vec3::X * h, &blobs, DEFAULT_BLOB_FACTOR))
                / (2.0 * h),
            (density(p + Vec3::Y * h, &blobs, DEFAULT_BLOB_FACTOR)
                - density(p - Vec3::Y * h, &blobs, DEFAULT_BLOB_FACTOR))
                / (2.0 * h),
            (density(p + Vec3::Z * h, &blobs, DEFAULT_BLOB_FACTOR)
                - density(p - Vec3::Z * h, &blobs, DEFAULT_BLOB_FACTOR))
                / (2.0 * h),
        );
        let analytic = gradient(p, &blobs, DEFAULT_BLOB_FACTOR);
        assert!(
            (numeric - analytic).length() < 1e-2,
            "numeric {numeric:?} vs analytic {analytic:?}"
        );
    }

    #[test]
    fn cutoff_radius_is_conservative() {
        let r = cutoff_radius(1.5, DEFAULT_BLOB_FACTOR, 0.01);
        assert!(r > 1.5, "cutoff must be past the atom's own radius");
        let b = blob(Vec3::ZERO, 1.5);
        let d = density(Vec3::X * r, &[b], DEFAULT_BLOB_FACTOR);
        assert!((d - 0.01).abs() < 1e-4, "density at cutoff: {d}");
    }
}
