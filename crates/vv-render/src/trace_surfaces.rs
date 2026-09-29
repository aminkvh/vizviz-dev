//! Analytic surfaces for the path tracer: the SES and skin patches the
//! viewport ray-casts, packed as the viewport packs them and hit by the
//! same WGSL (`ses_patch.wgsl`, `skin_patch.wgsl`), so a render shows
//! exactly the surface on screen.
//!
//! Every surface of a scene shares one set of arrays: each push appends
//! its atoms, probes and patches and offsets their cross-references.

use glam::Vec3;
use vv_core::ses::{Ses, FULL_CIRCLE};
use vv_core::skin_surface::SkinComplex;

use crate::path_trace::TraceScene;
use crate::scene::SkinPatchGpu;
use crate::ses_surface::{
    patch_records, probe_records, ProbeGpu, KIND_CONCAVE, KIND_CONVEX, KIND_TORUS,
};
use crate::style::Material;

/// A skin patch's record: this kind, and its index in `skin`.
pub(crate) const KIND_SKIN: u32 = 3 << 30;
const INDEX: u32 = (1 << 30) - 1;

/// The patches of every surface, and what their records point into.
#[derive(Clone, Debug)]
pub struct TraceSurfaces {
    /// SES records as `SesLayout` packs them, or `KIND_SKIN | index`.
    pub(crate) records: Vec<[u32; 4]>,
    /// Each record's bounding sphere.
    pub(crate) bounds: Vec<(Vec3, f32)>,
    /// xyz = position; w = van der Waals radius (SES) or weight (skin).
    pub(crate) atoms: Vec<[f32; 4]>,
    /// (colour, material, -, -) per atom.
    pub(crate) looks: Vec<[u32; 4]>,
    pub(crate) caps: Vec<u32>,
    pub(crate) probes: Vec<ProbeGpu>,
    pub(crate) probe_neighbors: Vec<u32>,
    pub(crate) skin: Vec<SkinPatchGpu>,
    /// CSR over `skin`, as `SkinSurfaceGpu` lays it out.
    pub(crate) competitor_starts: Vec<u32>,
    pub(crate) competitors: Vec<u32>,
    /// One probe radius for every SES and one shrink factor for every
    /// skin surface: the first pushed sets it.
    pub(crate) probe_radius: Option<f32>,
    pub(crate) shrink: Option<f32>,
}

impl Default for TraceSurfaces {
    fn default() -> Self {
        Self {
            records: Vec::new(),
            bounds: Vec::new(),
            atoms: Vec::new(),
            looks: Vec::new(),
            caps: Vec::new(),
            probes: Vec::new(),
            probe_neighbors: Vec::new(),
            skin: Vec::new(),
            competitor_starts: vec![0],
            competitors: Vec::new(),
            probe_radius: None,
            shrink: None,
        }
    }
}

impl TraceSurfaces {
    /// Appends atoms with their looks; returns the first one's index.
    fn push_atoms(&mut self, positions: &[Vec3], w: &[f32], colors: &[u32], material: u32) -> u32 {
        assert!(positions.len() == w.len() && w.len() == colors.len());
        let base = self.atoms.len() as u32;
        self.atoms.extend(
            positions
                .iter()
                .zip(w)
                .map(|(p, &w)| p.extend(w).to_array()),
        );
        self.looks
            .extend(colors.iter().map(|&c| [c, material, 0, 0]));
        base
    }
}

impl TraceScene {
    /// Adds the SES `ses` of atoms at `positions` with van der Waals
    /// `radii`, coloured per atom by `colors`.
    pub fn push_ses(
        &mut self,
        ses: &Ses,
        positions: &[Vec3],
        radii: &[f32],
        colors: &[u32],
        material: Material,
    ) {
        let m = self.materials.len() as u32;
        self.materials.push(material);
        let s = &mut self.surfaces;
        let probe = *s.probe_radius.get_or_insert(ses.probe_radius);
        assert_eq!(ses.probe_radius, probe, "one probe radius per scene");
        let atom = s.push_atoms(positions, radii, colors, m);
        let cap = s.caps.len() as u32;
        let probe = s.probes.len() as u32;
        let neighbor = s.probe_neighbors.len() as u32;
        s.caps.extend(ses.caps.iter().map(|&a| a + atom));
        s.probes
            .extend(probe_records(ses).into_iter().map(|p| ProbeGpu {
                atoms: p.atoms.map(|a| a + atom),
                neighbors: p.neighbors.map(|n| n + neighbor),
                ..p
            }));
        s.probe_neighbors
            .extend(ses.probe_neighbors.iter().map(|&p| p + probe));
        let arc_end = |p: u32| if p == FULL_CIRCLE { p } else { p + probe };
        for (r, center, radius, _) in patch_records(ses, positions, radii) {
            let kind = r[0] & !INDEX;
            let first = r[0] & INDEX;
            s.records.push(match kind {
                KIND_CONVEX => [kind | (first + atom), r[1] + cap, r[2] + cap, 0],
                KIND_TORUS => [
                    kind | (first + atom),
                    r[1] + atom,
                    arc_end(r[2]),
                    arc_end(r[3]),
                ],
                _ => {
                    debug_assert_eq!(kind, KIND_CONCAVE);
                    [kind | (first + probe), 0, 0, 0]
                }
            });
            s.bounds.push((center, radius));
        }
    }

    /// Adds the skin surface `complex` of atoms at `positions` with skin
    /// `weights` (`vv_core::skin_surface::weight_for_radius`), coloured per
    /// atom by `colors`.
    pub fn push_skin(
        &mut self,
        complex: &SkinComplex,
        positions: &[Vec3],
        weights: &[f32],
        colors: &[u32],
        material: Material,
    ) {
        let m = self.materials.len() as u32;
        self.materials.push(material);
        let s = &mut self.surfaces;
        let shrink = *s.shrink.get_or_insert(complex.shrink);
        assert_eq!(complex.shrink, shrink, "one shrink factor per scene");
        let atom = s.push_atoms(positions, weights, colors, m);
        let shift = |a: u32| if a == u32::MAX { a } else { a + atom };
        for (i, p) in complex.patches.iter().enumerate() {
            if !p.has_surface() {
                continue;
            }
            let mut gpu = SkinPatchGpu::from_patch(p, complex.shrink);
            gpu.atoms = gpu.atoms.map(shift);
            s.records.push([KIND_SKIN | s.skin.len() as u32, 0, 0, 0]);
            s.bounds
                .push((Vec3::from(gpu.bound_center), gpu.bound_radius));
            s.skin.push(gpu);
            s.competitors
                .extend(complex.competitors(i).iter().map(|&a| a + atom));
            s.competitor_starts.push(s.competitors.len() as u32);
        }
    }
}
