//! Ribbon/arrow cartoon geometry from a backbone trace and a secondary-
//! structure assignment, after M. Carson & C. Bugg, J. Mol. Graphics 4:121,
//! 1986 and Carson, J. Appl. Cryst. 24:958, 1991:
//!
//! - A Catmull-Rom spline through the CA trace (`backbone::catmull_rom`).
//!   Inside a strand the CAs are first averaged with their neighbours,
//!   which removes the pleat, so strands run straight and flat.
//! - **Orientation from the peptide plane.** The ribbon's width follows
//!   each residue's carbonyl (C to O) direction, projected perpendicular
//!   to the trace. In a sheet the carbonyls lie in the sheet plane, so
//!   strands lie flat in it; in a helix they run roughly along the axis,
//!   so the ribbon face points outward. Carbonyls alternate sign along a
//!   strand, so consecutive guides are sign-corrected *before* they are
//!   smoothed (averaging uncorrected ones cancels them). Residues without
//!   C and O (CA-only models, nucleic acids) fall back to the local
//!   curvature normal (used everywhere it would stand strands on edge:
//!   the pleat's curvature points out of the sheet).
//! - Cross-sections are [`RING`]-point superellipses: flat with rounded
//!   edges for helices and strands, round for coil, blended across a
//!   boundary. A strand ends in an arrowhead: a sharp step out to
//!   [`ARROW_FLARE`] times the ribbon width, then a taper to a point.
//! - Isolated bridges and runs too short to read (strand < 2, helix < 3
//!   residues) draw as coil, as the reference programs do.
//!
//! Scope: geometry only (`vv_render::CartoonGpu` uploads it). No caps at
//! chain breaks; nucleic acids get the coil tube.

use glam::Vec3;
use rayon::prelude::*;

use crate::backbone::{catmull_rom, trace, Trace};
use crate::dssp::DsspCode;
use crate::topology::{SecondaryStructure, Topology};

/// Ribbon width (the wide, flat face) for a helix or strand residue, in
/// Angstroms.
pub const RIBBON_WIDTH: f32 = 2.2;
/// Ribbon thickness (the narrow face), Angstroms.
pub const RIBBON_THICKNESS: f32 = 0.5;
/// Coil tube diameter (round: width = thickness).
pub const COIL_WIDTH: f32 = 0.6;
pub const COIL_THICKNESS: f32 = 0.6;
/// Spline samples per residue.
pub const SAMPLES_PER_RESIDUE: usize = 6;
/// Vertices per cross-section. Vertex 0 is the `+width` edge and vertex
/// `RING / 2` the `-width` edge, so the width is the distance between
/// them.
pub const RING: usize = 8;
/// The arrowhead's base, as a multiple of the ribbon width.
pub const ARROW_FLARE: f32 = 1.6;
/// Superellipse exponents: small = box-like, 1 = ellipse.
const RIBBON_ROUNDNESS: f32 = 0.35;
const COIL_ROUNDNESS: f32 = 1.0;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Shape {
    Helix,
    Strand,
    Coil,
}

/// What each residue draws as, after dropping fragments too short to
/// read: isolated bridges, strands under 2 and helices under 3 residues.
fn shapes(codes: &[DsspCode]) -> Vec<Shape> {
    let raw: Vec<Shape> = codes
        .iter()
        .map(|&c| match c {
            DsspCode::AlphaHelix | DsspCode::Helix3_10 | DsspCode::HelixPi => Shape::Helix,
            DsspCode::Strand => Shape::Strand,
            _ => Shape::Coil,
        })
        .collect();
    let mut out = raw.clone();
    let mut i = 0;
    while i < raw.len() {
        let mut j = i;
        while j + 1 < raw.len() && raw[j + 1] == raw[i] {
            j += 1;
        }
        let len = j - i + 1;
        let too_short = match raw[i] {
            Shape::Strand => len < 2,
            Shape::Helix => len < 3,
            Shape::Coil => false,
        };
        if too_short {
            out[i..=j].fill(Shape::Coil);
        }
        i = j + 1;
    }
    out
}

/// Width, thickness and roundness of a residue's cross-section.
fn profile(shape: Shape) -> (f32, f32, f32) {
    match shape {
        Shape::Helix | Shape::Strand => (RIBBON_WIDTH, RIBBON_THICKNESS, RIBBON_ROUNDNESS),
        Shape::Coil => (COIL_WIDTH, COIL_THICKNESS, COIL_ROUNDNESS),
    }
}

/// One cross-section: a superellipse of `half_width` along `across` and
/// `half_thickness` along `up` around `center`. The [`RING`] vertices are
/// derived from it ([`CartoonMesh::vertex`]; the GPU does the same in its
/// vertex shader), so only sections are stored.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CartoonSection {
    pub center: Vec3,
    pub across: Vec3,
    pub up: Vec3,
    pub half_width: f32,
    pub half_thickness: f32,
    /// Superellipse exponent: 1 is an ellipse, smaller flattens the faces.
    pub roundness: f32,
    /// Nearest trace atom, for picking and coloring.
    pub source: u32,
}

/// Cartoon geometry: cross-sections in trace order; `joins`, the
/// sections `k` that connect to `k + 1` with a tube of triangles (not
/// across a chain break); and `spans`, one `[first, last]` section range
/// per residue step. Spans share their end sections with their
/// neighbours, so a renderer may drop sections inside a span (level of
/// detail) without opening cracks.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct CartoonMesh {
    pub sections: Vec<CartoonSection>,
    pub joins: Vec<u32>,
    pub spans: Vec<[u32; 2]>,
}

impl CartoonMesh {
    /// Position and outward normal of vertex `i` (`0..RING`) of section
    /// `k`: vertex 0 is the `+across` edge, `RING / 2` the `-across` edge.
    pub fn vertex(&self, k: usize, i: usize) -> (Vec3, Vec3) {
        self.vertex_ring(k, i, RING)
    }

    /// As [`Self::vertex`], with the cross-section polygon's vertex count
    /// given explicitly instead of fixed at [`RING`] (`expand_with_ring`).
    fn vertex_ring(&self, k: usize, i: usize, ring: usize) -> (Vec3, Vec3) {
        let s = &self.sections[k];
        let theta = std::f32::consts::TAU * i as f32 / ring as f32;
        let (sin, cos) = theta.sin_cos();
        let pow = |x: f32, e: f32| x.signum() * x.abs().powf(e);
        let (a, b) = (s.half_width.max(1e-4), s.half_thickness.max(1e-4));
        let x = s.half_width * pow(cos, s.roundness);
        let y = s.half_thickness * pow(sin, s.roundness);
        // Gradient of the implicit superellipse at this parameter.
        let nx = pow(cos, 2.0 - s.roundness) / a;
        let ny = pow(sin, 2.0 - s.roundness) / b;
        (
            s.center + s.across * x + s.up * y,
            (s.across * nx + s.up * ny).normalize_or_zero(),
        )
    }

    /// The full triangle mesh (for export, and tests): vertex
    /// `k * RING + i` is `vertex(k, i)`, triangles wound counter-clockwise
    /// seen from outside.
    pub fn expand(&self) -> ExpandedMesh {
        self.expand_with_ring(RING)
    }

    /// As [`Self::expand`], tessellating each cross-section as an `ring`-gon
    /// instead of the viewport's fixed [`RING`] -- a denser polygon for a
    /// path-traced mesh, built once per render rather than every frame.
    pub fn expand_with_ring(&self, ring: usize) -> ExpandedMesh {
        let n = self.sections.len() * ring;
        let (mut positions, mut normals) = (Vec::with_capacity(n), Vec::with_capacity(n));
        for k in 0..self.sections.len() {
            for i in 0..ring {
                let (p, nrm) = self.vertex_ring(k, i, ring);
                positions.push(p);
                normals.push(nrm);
            }
        }
        let ring = ring as u32;
        let indices = self
            .joins
            .iter()
            .flat_map(|&k| {
                let (a, b) = (k * ring, (k + 1) * ring);
                (0..ring).flat_map(move |i| {
                    let i2 = (i + 1) % ring;
                    [a + i, a + i2, b + i2, a + i, b + i2, b + i]
                })
            })
            .collect();
        let source = self
            .sections
            .iter()
            .flat_map(|s| std::iter::repeat_n(s.source, ring as usize))
            .collect();
        ExpandedMesh {
            positions,
            normals,
            indices,
            source,
        }
    }

    /// Both ends of every run of connected spans (chain end, break or
    /// selection cut): the end's trace atom and half-width, for capping.
    pub fn loose_ends(&self) -> Vec<(u32, f32)> {
        let mut ends = Vec::new();
        let mut i = 0;
        while i < self.spans.len() {
            let mut j = i;
            while j + 1 < self.spans.len() && self.spans[j][1] == self.spans[j + 1][0] {
                j += 1;
            }
            let (start, end) = (self.spans[i][0], self.spans[j][1]);
            let (s0, s1) = (&self.sections[start as usize], &self.sections[end as usize]);
            ends.push((s0.source, s0.half_width));
            ends.push((s1.source, s1.half_width));
            i = j + 1;
        }
        ends
    }
}

/// [`CartoonMesh::expand`]'s output: [`RING`] vertices per section.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ExpandedMesh {
    pub positions: Vec<Vec3>,
    pub normals: Vec<Vec3>,
    pub indices: Vec<u32>,
    pub source: Vec<u32>,
}

/// The secondary structure to draw: the file's own HELIX/SHEET records
/// when it has any, else DSSP on
/// `positions`. File records only say helix/strand/coil, so DSSP still
/// supplies the finer kinds inside them (3-10 and pi helices within a
/// file helix, turns within file coil). Trajectories (`single_frame == false`) always use
/// DSSP: file records describe one conformation, not every frame.
pub fn secondary_structure(
    topology: &Topology,
    positions: &[Vec3],
    single_frame: bool,
) -> Vec<DsspCode> {
    let dssp = crate::dssp::assign(topology, positions);
    let has_records = topology
        .residues
        .iter()
        .any(|r| matches!(r.ss, SecondaryStructure::Helix | SecondaryStructure::Strand));
    if !single_frame || !has_records {
        return dssp;
    }
    topology
        .residues
        .iter()
        .zip(&dssp)
        .map(|(r, &d)| match (r.ss, d) {
            (SecondaryStructure::Helix, DsspCode::Helix3_10 | DsspCode::HelixPi) => d,
            (SecondaryStructure::Helix, _) => DsspCode::AlphaHelix,
            (SecondaryStructure::Strand, _) => DsspCode::Strand,
            (SecondaryStructure::Coil, DsspCode::Turn | DsspCode::Bend) => d,
            (SecondaryStructure::Coil, _) => DsspCode::Coil,
            (SecondaryStructure::Unknown, _) => DsspCode::None,
        })
        .collect()
}

/// How a cartoon is put together, fixed while its atoms move: the trace's
/// segments (chain breaks where they are when planned), each residue's
/// shape, and every cross-section's place along the spline and profile.
/// Only the spline's control points and guides follow the atoms
/// ([`CartoonPlan::frame`]), so a playing trajectory re-evaluates its
/// sections ([`section`]; on the GPU, `vv_render`'s cartoon frame pass)
/// instead of rebuilding the mesh.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct CartoonPlan {
    /// Trace atom of each residue slot, every segment's in turn.
    atoms: Vec<u32>,
    /// Each slot's residue's C and O atoms, when it has both.
    carbonyls: Vec<Option<[u32; 2]>>,
    shapes: Vec<Shape>,
    /// `[first, end)` slots of each segment.
    segments: Vec<[u32; 2]>,
    pub recipes: Vec<SectionRecipe>,
    /// As [`CartoonMesh::joins`] and [`CartoonMesh::spans`].
    pub joins: Vec<u32>,
    pub spans: Vec<[u32; 2]>,
}

/// Where one cross-section sits: `t` of the way from residue slot `slot`
/// to the next along the spline of the segment of slots `first..end`,
/// with its profile.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SectionRecipe {
    pub slot: u32,
    pub t: f32,
    pub first: u32,
    pub end: u32,
    pub half_width: f32,
    pub half_thickness: f32,
    pub roundness: f32,
    /// Nearest trace atom, for picking and colouring.
    pub source: u32,
}

/// One frame's spline: each residue slot's control point and smoothed
/// guide (the ribbon's width direction before it is made perpendicular to
/// the spline).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct CartoonFrame {
    pub controls: Vec<Vec3>,
    pub guides: Vec<Vec3>,
}

/// Plans the cartoon of every polymer segment of `topology`, segmented as
/// in `positions` and shaped by `codes` (one per residue, e.g. from
/// [`secondary_structure`]).
pub fn plan(topology: &Topology, positions: &[Vec3], codes: &[DsspCode]) -> CartoonPlan {
    let trace = trace(topology, positions);
    let mut plan = CartoonPlan::default();
    for segment in &trace.segments {
        let first = plan.atoms.len();
        let residues: Vec<usize> = segment
            .iter()
            .map(|&a| topology.residue_index[a as usize] as usize)
            .collect();
        let seg_codes: Vec<DsspCode> = residues.iter().map(|&r| codes[r]).collect();
        plan.atoms.extend(segment);
        plan.carbonyls
            .extend(residues.iter().map(|&r| carbonyl_atoms(topology, r)));
        plan.shapes.extend(shapes(&seg_codes));
        plan.segments.push([first as u32, plan.atoms.len() as u32]);
        plan.plan_segment(first, plan.atoms.len(), SAMPLES_PER_RESIDUE);
    }
    plan
}

/// Builds cartoon geometry for every polymer segment of `topology`, using
/// `codes` (one per residue, e.g. from [`secondary_structure`]) to shape
/// each residue's cross-section.
pub fn build(topology: &Topology, positions: &[Vec3], codes: &[DsspCode]) -> CartoonMesh {
    let plan = plan(topology, positions, codes);
    plan.mesh(&plan.frame(positions))
}

/// A round tube through `trace`'s segments with each residue's own
/// `radius` (a putty when it follows B-factor), built as a cartoon
/// plan so it shares the ribbon's mesh path: picking, clipping, level of
/// detail and path tracing.
pub fn tube_plan(trace: &Trace, radius: impl Fn(u32) -> f32) -> CartoonPlan {
    tube_plan_with_density(trace, radius, SAMPLES_PER_RESIDUE)
}

/// As [`tube_plan`], with the spline samples per residue given explicitly
/// instead of fixed at [`SAMPLES_PER_RESIDUE`] -- a denser tube for path
/// tracing, rebuilt fresh each render (cheap: no DSSP, unlike a ribbon).
pub fn tube_plan_with_density(
    trace: &Trace,
    radius: impl Fn(u32) -> f32,
    samples_per_residue: usize,
) -> CartoonPlan {
    let mut plan = CartoonPlan::default();
    for segment in &trace.segments {
        let first = plan.atoms.len();
        plan.atoms.extend(segment.iter().copied());
        plan.carbonyls.extend(segment.iter().map(|_| None));
        plan.shapes.extend(segment.iter().map(|_| Shape::Coil));
        plan.segments.push([first as u32, plan.atoms.len() as u32]);
        let radii: Vec<f32> = segment.iter().map(|&a| radius(a)).collect();
        plan.plan_tube_segment(first, plan.atoms.len(), &radii, samples_per_residue);
    }
    plan
}

/// Residue `r`'s C and O atoms, when it has both.
fn carbonyl_atoms(topology: &Topology, r: usize) -> Option<[u32; 2]> {
    let atoms = topology.residues[r].atoms.clone();
    let find = |name: &str| {
        atoms
            .clone()
            .find(|&a| topology.atom_name(a as usize) == name)
    };
    Some([find("C")?, find("O")?])
}

/// `v` with its component along unit `axis` removed, normalized; zero
/// when `v` is (nearly) parallel to it.
fn perpendicular(v: Vec3, axis: Vec3) -> Vec3 {
    (v - axis * v.dot(axis)).normalize_or_zero()
}

/// Any unit vector perpendicular to unit `axis`.
fn any_perpendicular(axis: Vec3) -> Vec3 {
    let helper = if axis.x.abs() < 0.9 { Vec3::X } else { Vec3::Y };
    perpendicular(helper, axis)
}

impl CartoonPlan {
    /// Residue slots (trace atoms).
    pub fn slot_count(&self) -> usize {
        self.atoms.len()
    }

    /// Each slot's trace atom (the spline's control points follow them).
    pub fn slot_atoms(&self) -> &[u32] {
        &self.atoms
    }

    /// The part drawn for the atoms `keep` accepts: every residue step
    /// (span) whose end sections' trace atoms are both kept, with its
    /// sections and joins. Spans still tile the joins, and neighbouring
    /// kept spans still share their end section. The spline (every slot)
    /// is unchanged, so frames made for the whole plan fit.
    pub fn filter(&self, keep: impl Fn(u32) -> bool) -> CartoonPlan {
        let mut new_index = vec![u32::MAX; self.recipes.len()];
        let mut out = CartoonPlan {
            atoms: self.atoms.clone(),
            carbonyls: self.carbonyls.clone(),
            shapes: self.shapes.clone(),
            segments: self.segments.clone(),
            ..CartoonPlan::default()
        };
        for &[a, b] in &self.spans {
            let (first, last) = (&self.recipes[a as usize], &self.recipes[b as usize]);
            if !(keep(first.source) && keep(last.source)) {
                continue;
            }
            for k in a..=b {
                if new_index[k as usize] == u32::MAX {
                    new_index[k as usize] = out.recipes.len() as u32;
                    out.recipes.push(self.recipes[k as usize]);
                }
            }
            out.joins.extend((a..b).map(|k| new_index[k as usize]));
            out.spans
                .push([new_index[a as usize], new_index[b as usize]]);
        }
        out
    }

    /// This plan's sections rebuilt at a different `samples_per_residue`
    /// (more for a denser path-traced mesh, without the viewport's per-
    /// frame budget): reuses the stored shapes and segments, so unlike
    /// [`plan`] this needs neither DSSP nor the topology again, only
    /// [`Self::plan_segment`]'s cheap arithmetic. The spline ([`Self::
    /// frame`]) does not depend on section density, so it is still valid
    /// for the result.
    pub fn redensify(&self, samples_per_residue: usize) -> CartoonPlan {
        let mut out = CartoonPlan {
            atoms: self.atoms.clone(),
            carbonyls: self.carbonyls.clone(),
            shapes: self.shapes.clone(),
            segments: self.segments.clone(),
            ..CartoonPlan::default()
        };
        for &[first, end] in &self.segments {
            out.plan_segment(first as usize, end as usize, samples_per_residue);
        }
        out
    }

    /// Sections along slots `first..end`, one segment, `samples_per_residue`
    /// cross-sections per residue step.
    fn plan_segment(&mut self, first: usize, end: usize, samples_per_residue: usize) {
        let n = end - first;
        let shape = &self.shapes[first..end];
        let strand_at = |i: usize| i < n && shape[i] == Shape::Strand;
        // Arrowheads: the segment from a strand's second-to-last residue
        // to its last flares out and tapers to a point.
        let arrow_start = |i: usize| strand_at(i) && strand_at(i + 1) && !strand_at(i + 2);
        let profiles: Vec<(f32, f32, f32)> = (0..n)
            .map(|i| {
                let (w, t, r) = profile(shape[i]);
                // The tip.
                if i > 0 && strand_at(i) && !strand_at(i + 1) && strand_at(i - 1) {
                    (0.0, t, r)
                } else {
                    (w, t, r)
                }
            })
            .collect();
        let atoms = &self.atoms[first..end];
        let recipe = |i: usize, t: f32, w: f32, th: f32, r: f32| SectionRecipe {
            slot: (first + i) as u32,
            t,
            first: first as u32,
            end: end as u32,
            half_width: w * 0.5,
            half_thickness: th * 0.5,
            roundness: r,
            source: if t < 0.5 { atoms[i] } else { atoms[i + 1] },
        };
        let base = self.recipes.len() as u32;
        let mut span_starts = Vec::with_capacity(n);
        for i in 0..n - 1 {
            span_starts.push(self.recipes.len() as u32);
            let (w0, t0, r0) = profiles[i];
            let (w1, t1, r1) = profiles[i + 1];
            for s in 0..samples_per_residue {
                let t = s as f32 / samples_per_residue as f32;
                let (th, r) = (t0 + (t1 - t0) * t, r0 + (r1 - r0) * t);
                if arrow_start(i) {
                    if s == 0 {
                        // The step out to the arrowhead's base.
                        self.recipes.push(recipe(i, 0.0, w0, th, r));
                    }
                    let flare = w0 * ARROW_FLARE;
                    self.recipes.push(recipe(i, t, flare * (1.0 - t), th, r));
                } else {
                    self.recipes.push(recipe(i, t, w0 + (w1 - w0) * t, th, r));
                }
            }
        }
        let (w, th, r) = profiles[n - 1];
        self.recipes.push(recipe(n - 2, 1.0, w, th, r));
        span_starts.push(self.recipes.len() as u32 - 1);
        self.joins.extend(base..self.recipes.len() as u32 - 1);
        self.spans
            .extend(span_starts.windows(2).map(|w| [w[0], w[1]]));
    }

    /// Round sections for tube segment `first..end`, the radius linear
    /// across each residue step as ribbon width is in `plan_segment`.
    fn plan_tube_segment(
        &mut self,
        first: usize,
        end: usize,
        radius: &[f32],
        samples_per_residue: usize,
    ) {
        let n = end - first;
        let atoms = &self.atoms[first..end];
        let recipe = |i: usize, t: f32, r: f32| SectionRecipe {
            slot: (first + i) as u32,
            t,
            first: first as u32,
            end: end as u32,
            half_width: r,
            half_thickness: r,
            roundness: 1.0,
            source: if t < 0.5 { atoms[i] } else { atoms[i + 1] },
        };
        let base = self.recipes.len() as u32;
        let mut span_starts = Vec::with_capacity(n);
        for i in 0..n - 1 {
            span_starts.push(self.recipes.len() as u32);
            let (r0, r1) = (radius[i], radius[i + 1]);
            for s in 0..samples_per_residue {
                let t = s as f32 / samples_per_residue as f32;
                self.recipes.push(recipe(i, t, r0 + (r1 - r0) * t));
            }
        }
        self.recipes.push(recipe(n - 2, 1.0, radius[n - 1]));
        span_starts.push(self.recipes.len() as u32 - 1);
        self.joins.extend(base..self.recipes.len() as u32 - 1);
        self.spans
            .extend(span_starts.windows(2).map(|w| [w[0], w[1]]));
    }

    /// The spline at `positions`: control points (CAs, with strand pleats
    /// averaged out) and guides (the carbonyl, or curvature where there is
    /// none, perpendicular to the trace, sign-consistent, then smoothed
    /// within helices and strands). Segments are independent, so they are
    /// done in parallel.
    pub fn frame(&self, positions: &[Vec3]) -> CartoonFrame {
        let per_segment: Vec<(Vec<Vec3>, Vec<Vec3>)> = self
            .segments
            .par_iter()
            .map(|&[first, end]| self.segment_frame(first as usize, end as usize, positions))
            .collect();
        let mut frame = CartoonFrame {
            controls: Vec::with_capacity(self.atoms.len()),
            guides: Vec::with_capacity(self.atoms.len()),
        };
        for (controls, guides) in per_segment {
            frame.controls.extend(controls);
            frame.guides.extend(guides);
        }
        frame
    }

    fn segment_frame(
        &self,
        first: usize,
        end: usize,
        positions: &[Vec3],
    ) -> (Vec<Vec3>, Vec<Vec3>) {
        let n = end - first;
        let shape = &self.shapes[first..end];
        let points: Vec<Vec3> = self.atoms[first..end]
            .iter()
            .map(|&a| positions[a as usize])
            .collect();
        let strand_at = |i: usize| i < n && shape[i] == Shape::Strand;
        let controls: Vec<Vec3> = (0..n)
            .map(|i| {
                if strand_at(i) && i > 0 && strand_at(i - 1) && strand_at(i + 1) {
                    (points[i - 1] + points[i] * 2.0 + points[i + 1]) * 0.25
                } else {
                    points[i]
                }
            })
            .collect();
        let at = |i: isize| controls[i.clamp(0, n as isize - 1) as usize];
        let tangent_at = |i: usize| (at(i as isize + 1) - at(i as isize - 1)).normalize_or_zero();

        let mut guides: Vec<Vec3> = (0..n)
            .map(|i| {
                let t = tangent_at(i);
                let raw = self.carbonyls[first + i]
                    .map(|[c, o]| positions[o as usize] - positions[c as usize])
                    .unwrap_or_else(|| {
                        let (prev, cur, next) =
                            (at(i as isize - 1), controls[i], at(i as isize + 1));
                        (next - cur).cross(cur - prev)
                    });
                perpendicular(raw, t)
            })
            .collect();
        for i in 0..n {
            if guides[i] == Vec3::ZERO {
                guides[i] = if i > 0 {
                    perpendicular(guides[i - 1], tangent_at(i))
                } else {
                    Vec3::ZERO
                };
            }
        }
        if let Some(found) = guides.iter().position(|g| *g != Vec3::ZERO) {
            for i in (0..found).rev() {
                guides[i] = perpendicular(guides[i + 1], tangent_at(i));
            }
        }
        for (i, guide) in guides.iter_mut().enumerate().take(n) {
            if *guide == Vec3::ZERO {
                *guide = any_perpendicular(if tangent_at(i) == Vec3::ZERO {
                    Vec3::X
                } else {
                    tangent_at(i)
                });
            }
        }
        for i in 1..n {
            if guides[i].dot(guides[i - 1]) < 0.0 {
                guides[i] = -guides[i];
            }
        }
        let smoothed = (0..n)
            .map(|i| {
                if shape[i] == Shape::Coil || i == 0 || i + 1 == n {
                    return guides[i];
                }
                let sum = guides[i - 1] + guides[i] * 2.0 + guides[i + 1];
                let g = perpendicular(sum, tangent_at(i));
                if g == Vec3::ZERO {
                    guides[i]
                } else {
                    g
                }
            })
            .collect();
        (controls, smoothed)
    }

    /// Every section at `frame` (from [`CartoonPlan::frame`]).
    pub fn mesh(&self, frame: &CartoonFrame) -> CartoonMesh {
        CartoonMesh {
            sections: self.recipes.par_iter().map(|r| section(r, frame)).collect(),
            joins: self.joins.clone(),
            spans: self.spans.clone(),
        }
    }
}

/// Cross-section `recipe` at `frame`: the Catmull-Rom point, its tangent,
/// and the guide made perpendicular to it. Mirrored by `vv_render`'s
/// `shaders/cartoon_frame.wgsl`.
pub fn section(recipe: &SectionRecipe, frame: &CartoonFrame) -> CartoonSection {
    let (i, first, last) = (
        recipe.slot as isize,
        recipe.first as isize,
        recipe.end as isize - 1,
    );
    let at = |k: isize| frame.controls[k.clamp(first, last) as usize];
    let (p0, p1, p2, p3) = (at(i - 1), at(i), at(i + 1), at(i + 2));
    let t = recipe.t;
    let center = catmull_rom(p0, p1, p2, p3, t);
    let d = 1e-3;
    let mut tangent = (catmull_rom(p0, p1, p2, p3, (t + d).min(1.0))
        - catmull_rom(p0, p1, p2, p3, (t - d).max(0.0)))
    .normalize_or_zero();
    if tangent == Vec3::ZERO {
        tangent = (p2 - p1).normalize_or_zero();
    }
    let guide = frame.guides[i as usize].lerp(frame.guides[(i + 1).min(last) as usize], t);
    let mut across = perpendicular(guide, tangent);
    if across == Vec3::ZERO {
        across = any_perpendicular(tangent);
    }
    CartoonSection {
        center,
        across,
        up: tangent.cross(across),
        half_width: recipe.half_width,
        half_thickness: recipe.half_thickness,
        roundness: recipe.roundness,
        source: recipe.source,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::builder::{AtomRow, TopologyBuilder};

    fn straight_chain(n: usize, codes: &[DsspCode]) -> (Topology, Vec<Vec3>) {
        assert_eq!(n, codes.len());
        let mut b = TopologyBuilder::new();
        for i in 0..n {
            let mut name = [b' '; 4];
            name[..2].copy_from_slice(b"CA");
            b.push(&AtomRow {
                element: crate::Element::CARBON,
                name,
                serial: i as u32 + 1,
                alt_loc: 0,
                comp: "ALA",
                asym: "A",
                auth_asym: "A",
                seq_id: i as i32 + 1,
                auth_seq_id: i as i32 + 1,
                ins_code: 0,
                entity: 1,
                position: Vec3::new(i as f32 * 3.8, 0.0, 0.0),
                occupancy: 1.0,
                b_factor: 0.0,
                charge: 0,
                hetero: false,
            });
        }
        let structure = b.finish().unwrap();
        let positions = structure.frame(0).positions().to_vec();
        (
            std::sync::Arc::try_unwrap(structure.topology).unwrap(),
            positions,
        )
    }

    /// A dead-straight chain has no curvature, so the naive per-residue
    /// normal is always the zero vector; this only checks the mesh comes
    /// out well-formed (right vertex/index counts, everything finite),
    /// not any particular orientation - a curved case does that below.
    #[test]
    fn a_coil_segment_has_the_expected_vertex_and_index_counts() {
        let codes = vec![DsspCode::Coil; 6];
        let (t, p) = straight_chain(6, &codes);
        let mesh = build(&t, &p, &codes).expand();
        let cross_sections = (6 - 1) * SAMPLES_PER_RESIDUE + 1;
        assert_eq!(mesh.positions.len(), cross_sections * RING);
        assert_eq!(mesh.normals.len(), mesh.positions.len());
        assert_eq!(mesh.source.len(), mesh.positions.len());
        assert_eq!(mesh.indices.len(), (cross_sections - 1) * RING * 6);
        assert!(mesh
            .indices
            .iter()
            .all(|&i| (i as usize) < mesh.positions.len()));
        assert!(mesh.positions.iter().all(|p| p.is_finite()));
    }

    /// `redensify` reuses the plan's stored shapes rather than re-running
    /// DSSP; this pins that it produces the same section count and the
    /// same boundary sections (every residue slot's `t == 0`, and the
    /// segment's final `t == 1`) as the original, whatever the density --
    /// only the samples in between are new.
    #[test]
    fn redensify_keeps_the_same_residue_boundaries_at_a_new_density() {
        let n = 9;
        let codes: Vec<DsspCode> = (0..n)
            .map(|i| {
                if (3..7).contains(&i) {
                    DsspCode::Strand
                } else {
                    DsspCode::Coil
                }
            })
            .collect();
        let (t, p) = straight_chain(n, &codes);
        let base = plan(&t, &p, &codes);
        let dense = base.redensify(11);

        // Every step contributes `samples_per_residue` recipes (plus one
        // extra at an arrowhead's step-out, the same at any density, plus
        // one final endpoint): the density difference alone explains the
        // recipe-count difference.
        assert_eq!(
            dense.recipes.len() - base.recipes.len(),
            (n - 1) * (11 - SAMPLES_PER_RESIDUE)
        );

        let at_t0 = |p: &CartoonPlan, slot: usize| {
            p.recipes
                .iter()
                .find(|r| r.slot as usize == slot && r.t == 0.0)
                .copied()
                .unwrap_or_else(|| panic!("no t=0 recipe at slot {slot}"))
        };
        for slot in 0..n - 1 {
            assert_eq!(at_t0(&base, slot), at_t0(&dense, slot), "slot {slot}");
        }
        assert_eq!(base.recipes.last(), dense.recipes.last());
    }

    #[test]
    fn helix_cross_sections_are_wider_than_coil_ones() {
        // A gentle curve (not a straight line) so the ribbon frame has a
        // real, non-degenerate normal to measure width along.
        let n = 10;
        let codes: Vec<DsspCode> = (0..n)
            .map(|i| {
                if (3..7).contains(&i) {
                    DsspCode::AlphaHelix
                } else {
                    DsspCode::Coil
                }
            })
            .collect();
        let mut b = TopologyBuilder::new();
        for i in 0..n {
            let mut name = [b' '; 4];
            name[..2].copy_from_slice(b"CA");
            let a = i as f32 * 0.35;
            b.push(&AtomRow {
                element: crate::Element::CARBON,
                name,
                serial: i as u32 + 1,
                alt_loc: 0,
                comp: "ALA",
                asym: "A",
                auth_asym: "A",
                seq_id: i + 1,
                auth_seq_id: i + 1,
                ins_code: 0,
                entity: 1,
                position: Vec3::new(10.0 * a.cos(), 10.0 * a.sin(), 0.0),
                occupancy: 1.0,
                b_factor: 0.0,
                charge: 0,
                hetero: false,
            });
        }
        let structure = b.finish().unwrap();
        let positions = structure.frame(0).positions().to_vec();
        let topology = std::sync::Arc::try_unwrap(structure.topology).unwrap();
        let mesh = build(&topology, &positions, &codes).expand();

        // Ring k's width: vertex 0 is the +width edge, RING/2 the -width.
        let width_at = |k: usize| {
            let base = k * RING;
            mesh.positions[base].distance(mesh.positions[base + RING / 2])
        };
        let per_residue = SAMPLES_PER_RESIDUE;
        let helix_width = width_at(5 * per_residue);
        let coil_width = width_at(per_residue);
        assert!(
            (helix_width - RIBBON_WIDTH).abs() < 1e-3,
            "helix width {helix_width}"
        );
        assert!(
            (coil_width - COIL_WIDTH).abs() < 1e-3,
            "coil width {coil_width}"
        );
        assert!(helix_width > coil_width * 2.0);
    }

    #[test]
    fn a_strand_ends_in_a_flared_arrowhead_that_tapers_to_a_point() {
        let n = 8;
        let codes: Vec<DsspCode> = (0..n)
            .map(|i| {
                if i >= 2 {
                    DsspCode::Strand
                } else {
                    DsspCode::Coil
                }
            })
            .collect();
        let (t, p) = straight_chain(n, &codes);
        let mesh = build(&t, &p, &codes).expand();
        let rings = mesh.positions.len() / RING;
        let width =
            |k: usize| mesh.positions[k * RING].distance(mesh.positions[k * RING + RING / 2]);
        let widest = (0..rings).map(width).fold(0.0f32, f32::max);
        assert!(
            (widest - RIBBON_WIDTH * ARROW_FLARE).abs() < 1e-3,
            "the arrowhead's base should flare to {} A, widest ring is {widest}",
            RIBBON_WIDTH * ARROW_FLARE
        );
        // The final ring collapses onto the tip along the width axis.
        assert!(width(rings - 1) < 1e-3, "tip width {}", width(rings - 1));
    }

    /// Spans tile the joins exactly: consecutive, end-to-end, covering
    /// every join once.
    #[test]
    fn spans_tile_the_joins() {
        let codes: Vec<DsspCode> = (0..10)
            .map(|i| {
                if (3..7).contains(&i) {
                    DsspCode::Strand
                } else {
                    DsspCode::Coil
                }
            })
            .collect();
        let (t, p) = straight_chain(10, &codes);
        let mesh = build(&t, &p, &codes);
        assert_eq!(mesh.spans.len(), 9);
        assert_eq!(mesh.spans[0][0], 0);
        assert_eq!(
            mesh.spans.last().unwrap()[1] as usize,
            mesh.sections.len() - 1
        );
        for w in mesh.spans.windows(2) {
            assert_eq!(w[0][1], w[1][0], "spans must share end sections");
        }
        let covered: u32 = mesh.spans.iter().map(|[a, b]| b - a).sum();
        assert_eq!(covered as usize, mesh.joins.len());
    }

    #[test]
    fn filtering_keeps_whole_residue_steps_and_their_tiling() {
        let codes: Vec<DsspCode> = (0..10)
            .map(|i| {
                if (3..7).contains(&i) {
                    DsspCode::Strand
                } else {
                    DsspCode::Coil
                }
            })
            .collect();
        let (t, p) = straight_chain(10, &codes);
        let plan = plan(&t, &p, &codes);
        let frame = plan.frame(&p);
        let mesh = plan.mesh(&frame);
        assert_eq!(plan.filter(|_| true), plan);
        assert!(plan.filter(|_| false).recipes.is_empty());

        // Keep the trace atoms of residues 2..=5: steps 2-3, 3-4, 4-5.
        let sources: Vec<u32> = mesh
            .spans
            .iter()
            .map(|&[a, _]| mesh.sections[a as usize].source)
            .collect();
        let kept: std::collections::HashSet<u32> = sources[2..6].iter().copied().collect();
        let part = plan.filter(|a| kept.contains(&a)).mesh(&frame);
        assert_eq!(part.spans.len(), 3);
        for w in part.spans.windows(2) {
            assert_eq!(w[0][1], w[1][0], "kept spans still share end sections");
        }
        let covered: u32 = part.spans.iter().map(|[a, b]| b - a).sum();
        assert_eq!(covered as usize, part.joins.len());
        assert_eq!(
            part.spans.last().unwrap()[1] as usize,
            part.sections.len() - 1
        );
    }

    #[test]
    fn short_fragments_and_isolated_bridges_draw_as_coil() {
        use DsspCode::*;
        let codes = [
            Coil, Bridge, Coil, Strand, Coil, AlphaHelix, AlphaHelix, Coil, Strand, Strand,
        ];
        let s = shapes(&codes);
        assert_eq!(s[1], Shape::Coil, "an isolated bridge");
        assert_eq!(s[3], Shape::Coil, "a one-residue strand");
        assert_eq!(s[5], Shape::Coil, "a two-residue helix");
        assert_eq!(s[8], Shape::Strand, "a two-residue strand is kept");
    }

    #[test]
    fn the_flip_guard_keeps_consecutive_cross_sections_from_reversing() {
        // A backbone that curves one way, then the other (an S-shape):
        // the naive per-residue normal flips sign at the inflection
        // without the flip guard, which would show up as corner 0 of one
        // cross-section suddenly sitting where corner 2 of the previous
        // one was. Checked by tracking that corner 0 moves smoothly
        // (small steps) all the way along, never jumping to the opposite
        // side of the ribbon.
        let n = 12;
        let codes = vec![DsspCode::AlphaHelix; n];
        let mut b = TopologyBuilder::new();
        for i in 0..n {
            let mut name = [b' '; 4];
            name[..2].copy_from_slice(b"CA");
            let x = i as f32 * 3.0;
            // Sine wave: curvature direction reverses every half period.
            let y = 8.0 * (x / 12.0).sin();
            b.push(&AtomRow {
                element: crate::Element::CARBON,
                name,
                serial: i as u32 + 1,
                alt_loc: 0,
                comp: "ALA",
                asym: "A",
                auth_asym: "A",
                seq_id: i as i32 + 1,
                auth_seq_id: i as i32 + 1,
                ins_code: 0,
                entity: 1,
                position: Vec3::new(x, y, 0.0),
                occupancy: 1.0,
                b_factor: 0.0,
                charge: 0,
                hetero: false,
            });
        }
        let structure = b.finish().unwrap();
        let positions = structure.frame(0).positions().to_vec();
        let topology = std::sync::Arc::try_unwrap(structure.topology).unwrap();
        let mesh = build(&topology, &positions, &codes).expand();

        let cross_sections = mesh.positions.len() / RING;
        let max_step = RIBBON_WIDTH; // a flip would jump by roughly 2x this
        for k in 1..cross_sections {
            let step = mesh.positions[k * RING].distance(mesh.positions[(k - 1) * RING]);
            assert!(
                step < max_step,
                "corner 0 jumped {step} A at cross-section {k}"
            );
        }
    }

    /// A necessary (not sufficient - see `push_cross_section`'s doc
    /// comment) precondition for a GPU pipeline to enable backface
    /// culling: every triangle's vertex order must itself agree with the
    /// outward direction its own vertex normals already point in.
    /// Concretely, `(pb-pa) x (pc-pa)` for triangle `(a, b, c)` should
    /// point the same way as `a`'s own vertex normal - the standard
    /// "counter-clockwise as seen from outside" convention almost every
    /// renderer assumes. If this ever fails on a triangle whose
    /// cross-section is *not* one of the excluded degenerate ones below,
    /// the bug is in `stitch`'s index order or `push_cross_section`'s
    /// corner order, not in `vv-render`'s pipeline setup.
    ///
    /// Excluded: triangles touching a cross-section with a zero-width (a
    /// strand's exact arrowhead tip, `control_width` returning `0.0`) or
    /// zero-thickness corner. There, two of the four corners are
    /// literally the same point, so `push_cross_section`'s vertex normal
    /// (`corner - center`) collapses onto the other axis alone and stops
    /// meaning "outward" for that vertex - a limitation of the normal
    /// proxy exactly at that singularity, confirmed on a straight,
    /// uncurved chain: every triangle away
    /// from a degenerate cross-section agreed perfectly (hundreds
    /// checked, zero disagreements); only ones touching the exact tip
    /// disagreed, all three with a near-zero dot product (magnitude
    /// under 3e-4) rather than a clean, confidently-wrong negative one -
    /// the signature of a degenerate proxy, not a real winding flip.
    #[test]
    fn triangle_winding_agrees_with_outward_vertex_normals() {
        let n = 16;
        let codes: Vec<DsspCode> = (0..n)
            .map(|i| {
                if (2..6).contains(&i) {
                    DsspCode::Strand
                } else {
                    DsspCode::AlphaHelix
                }
            })
            .collect();
        let mut b = TopologyBuilder::new();
        for i in 0..n {
            let mut name = [b' '; 4];
            name[..2].copy_from_slice(b"CA");
            let t = i as f32 * 0.4;
            // A curve with both in-plane bending and out-of-plane rise, so
            // the check exercises every axis of the local frame, not just
            // a flat circle.
            let position = Vec3::new(10.0 * t.cos(), 10.0 * t.sin(), t * 1.5);
            b.push(&AtomRow {
                element: crate::Element::CARBON,
                name,
                serial: i as u32 + 1,
                alt_loc: 0,
                comp: "ALA",
                asym: "A",
                auth_asym: "A",
                seq_id: i + 1,
                auth_seq_id: i + 1,
                ins_code: 0,
                entity: 1,
                position,
                occupancy: 1.0,
                b_factor: 0.0,
                charge: 0,
                hetero: false,
            });
        }
        let structure = b.finish().unwrap();
        let positions = structure.frame(0).positions().to_vec();
        let topology = std::sync::Arc::try_unwrap(structure.topology).unwrap();
        let mesh = build(&topology, &positions, &codes).expand();

        // A ring is degenerate when it has (nearly) zero width - the
        // arrowhead's tip - where the superellipse normal stops meaning
        // "outward". Triangles joining two rings at the same spline point
        // (the arrowhead's step out) are flat walls facing along the
        // trace, which radial vertex normals don't describe either.
        let ring_of = |v: usize| v / RING;
        let center = |r: usize| {
            mesh.positions[r * RING..(r + 1) * RING]
                .iter()
                .copied()
                .sum::<Vec3>()
                / RING as f32
        };
        let ring_width =
            |r: usize| mesh.positions[r * RING].distance(mesh.positions[r * RING + RING / 2]);
        let cross_section_degenerate = |v: usize| ring_width(ring_of(v)) < 1e-2;
        let is_step = |a: usize, b: usize| {
            ring_of(a) != ring_of(b) && center(ring_of(a)).distance(center(ring_of(b))) < 1e-4
        };

        let mut checked = 0;
        let mut skipped = 0;
        for tri in mesh.indices.chunks_exact(3) {
            let (a, bb, c) = (tri[0] as usize, tri[1] as usize, tri[2] as usize);
            if cross_section_degenerate(a)
                || cross_section_degenerate(bb)
                || cross_section_degenerate(c)
                || is_step(a, bb)
                || is_step(a, c)
            {
                skipped += 1;
                continue;
            }
            let (pa, pb, pc) = (mesh.positions[a], mesh.positions[bb], mesh.positions[c]);
            let face_normal = (pb - pa).cross(pc - pa);
            assert!(
                face_normal.length_squared() > 1e-8,
                "triangle ({a},{bb},{c}) has near-zero area despite a non-degenerate cross-section"
            );
            checked += 1;
            assert!(
                face_normal.normalize().dot(mesh.normals[a]) > 0.0,
                "triangle ({a},{bb},{c}) winds against its own vertex normal"
            );
        }
        assert!(checked > 100, "too few checked triangles: {checked}");
        // The strand arrowhead's tip touches a handful of triangles on
        // both sides (the exact tip cross-section, plus the ones either
        // side of it); if this balloons, something broader than one tip
        // is degenerating.
        assert!(
            skipped > 0 && skipped <= 4 * RING * 2,
            "{skipped} skipped triangles"
        );
    }

    /// A straight chain, radius rising from one end to the other: every
    /// ring is a perfect circle (`tube_plan` sets `roundness = 1`,
    /// `half_width == half_thickness`) whose radius matches the given
    /// function exactly at each residue and interpolates smoothly (no
    /// jump bigger than one residue step's share of the total change)
    /// between them.
    #[test]
    fn tube_plan_is_a_smoothly_varying_round_tube() {
        let n = 6;
        let codes = vec![DsspCode::Coil; n];
        let (t, p) = straight_chain(n, &codes);
        let trace = crate::backbone::trace(&t, &p);
        let radius = |a: u32| if a == 0 { 0.2 } else { 1.0 };
        let plan = tube_plan(&trace, radius);
        let frame = plan.frame(&p);
        let mesh = plan.mesh(&frame);

        for s in &mesh.sections {
            assert!((s.half_width - s.half_thickness).abs() < 1e-6);
            assert_eq!(s.roundness, 1.0);
        }
        assert!((mesh.sections[0].half_width - 0.2).abs() < 1e-6);
        assert!((mesh.sections.last().unwrap().half_width - 1.0).abs() < 1e-6);
        let max_step = 0.8 / SAMPLES_PER_RESIDUE as f32 + 1e-6;
        for w in mesh.sections.windows(2) {
            assert!(
                (w[1].half_width - w[0].half_width).abs() <= max_step,
                "radius jumped from {} to {}",
                w[0].half_width,
                w[1].half_width
            );
        }
    }

    /// Two disjoint runs (a chain break, or a selection cut after
    /// `CartoonPlan::filter`): `loose_ends` reports each run's own two
    /// termini, not the whole plan's first and last section.
    #[test]
    fn loose_ends_are_reported_per_disjoint_run() {
        let trace = Trace {
            segments: vec![vec![0, 1, 2], vec![5, 6, 7]],
        };
        let plan = tube_plan(&trace, |_| 0.5);
        let positions: Vec<Vec3> = (0..8)
            .map(|i| Vec3::new(i as f32 * 3.8, 0.0, 0.0))
            .collect();
        let frame = plan.frame(&positions);
        let mesh = plan.mesh(&frame);
        assert_eq!(
            mesh.loose_ends(),
            vec![(0, 0.5), (2, 0.5), (5, 0.5), (7, 0.5)]
        );
    }
}
