//! CPU renderers: ray-sphere impostors (`renderer`), a Gaussian-surface
//! ray-marcher (`gaussian_surface`), and a skin-surface ray-caster
//! (`skin_surface`).
//!
//! For machines with no usable GPU (an HPC login node, a headless CI box
//! without llvmpipe/WARP, a laptop with a broken driver) rather than for
//! everyday use — see docs/RENDERING.md for the fps envelope the sphere
//! path is actually good for. Scalar ray intersection, a BVH accelerating
//! nearest-hit queries (spheres and skin-surface patches alike), a
//! `vv_core::spatial::Grid` for the Gaussian march's per-step atom
//! lookups, parallelized by image row via `rayon`.
//!
//! Deliberately reuses `vv_render::{Camera, color::*, StylePreset}` rather
//! than re-deriving camera or shading math, so the two backends render the
//! same scene identically wherever they overlap.

mod bvh;
pub mod gaussian_surface;
pub mod renderer;
pub mod skin_surface;

pub use gaussian_surface::{render_gaussian_surface, GaussianSurfaceScene};
pub use renderer::{render_spacefill, render_spacefill_colored, render_spheres, CpuScene};
pub use skin_surface::{render_skin_surface, SkinSurfaceScene};
