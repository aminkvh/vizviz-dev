//! Offline path-traced images: the scene's spheres, bond cylinders,
//! triangle meshes (cartoons, Gaussian surfaces) and SES and skin patches
//! (`crate::trace_surfaces`) ray-traced on the GPU
//! (`shaders/path_trace.wgsl`) with real shadows
//! and ambient occlusion -- the
//! viewport's own shading (`shaders/shading.wgsl`'s `light_surface`),
//! with each light seen only where a shadow ray reaches it and the
//! ambient light only where an occlusion ray escapes to the sky.
//!
//! The image is rendered a tile at a time, every sample of a tile before
//! the next, so no submit runs long enough for the OS to reset the GPU
//! and nothing is image-sized on the GPU at any resolution. Primitives
//! sit in a bounding-volume hierarchy built on the CPU (binned surface
//! area heuristic).

use bytemuck::{Pod, Zeroable};
use glam::{Vec3, Vec4};
use wgpu::util::DeviceExt;

use crate::camera::{Camera, Projection, KEEP_ALL};
use crate::context::GpuContext;
use crate::scene::{try_buffer_init, OutOfGpuMemory};
use crate::style::{Lighting, Material};
use crate::trace_surfaces::TraceSurfaces;

/// What a path-traced image is made of.
#[derive(Clone, Debug, Default)]
pub struct TraceScene {
    pub spheres: Vec<TraceSphere>,
    pub cylinders: Vec<TraceCylinder>,
    /// Triangles, as indices into `vertices`.
    pub vertices: Vec<TraceVertex>,
    pub triangles: Vec<[u32; 3]>,
    /// Indexed by the primitives' `material`.
    pub materials: Vec<Material>,
    /// SES and skin patches (`push_ses`, `push_skin`).
    pub surfaces: TraceSurfaces,
}

impl TraceScene {
    /// Adds atoms drawn with `material`: a sphere of `radii[a]` for each
    /// atom `a` of `atoms` with a positive radius, and a cylinder of
    /// `bond_radius` for each of `bonds` (when positive), each half
    /// coloured as its atom.
    pub fn push_atoms(
        &mut self,
        positions: &[Vec3],
        radii: &[f32],
        colors: &[u32],
        atoms: impl IntoIterator<Item = usize>,
        bonds: impl IntoIterator<Item = [u32; 2]>,
        bond_radius: f32,
        material: Material,
    ) {
        let m = self.materials.len() as u32;
        self.materials.push(material);
        self.spheres.extend(
            atoms
                .into_iter()
                .filter(|&a| radii[a] > 0.0)
                .map(|a| TraceSphere {
                    center: positions[a],
                    radius: radii[a],
                    color: colors[a],
                    material: m,
                }),
        );
        if bond_radius > 0.0 {
            self.cylinders.extend(bonds.into_iter().map(|[i, j]| {
                let (i, j) = (i as usize, j as usize);
                TraceCylinder {
                    a: positions[i],
                    b: positions[j],
                    radius: bond_radius,
                    color_a: colors[i],
                    color_b: colors[j],
                    material: m,
                }
            }));
        }
    }
}

impl TraceScene {
    /// Adds a cartoon's triangles (`vv_core::cartoon::CartoonMesh::
    /// expand`) drawn with `material`, each vertex coloured as
    /// `colors[source]` of the atom its section belongs to.
    pub fn push_mesh(
        &mut self,
        mesh: &vv_core::cartoon::ExpandedMesh,
        colors: &[u32],
        material: Material,
    ) {
        let m = self.materials.len() as u32;
        self.materials.push(material);
        let base = self.vertices.len() as u32;
        self.vertices.extend(
            mesh.positions
                .iter()
                .zip(&mesh.normals)
                .zip(&mesh.source)
                .map(|((&position, &normal), &atom)| TraceVertex {
                    position,
                    normal,
                    color: colors[atom as usize],
                    material: m,
                }),
        );
        self.triangles.extend(
            mesh.indices
                .chunks_exact(3)
                .map(|t| [base + t[0], base + t[1], base + t[2]]),
        );
    }

    /// Adds a glycan rep's SNFG shapes and linkage cylinders
    /// (`vv_core::glycan::mesh::build_mesh`/`build_linkage_mesh`), drawn
    /// with `material`. Unlike [`Self::push_mesh`], colors are already
    /// baked per vertex (the source script's fixed SNFG palette, not a
    /// per-atom color scheme), so there is no `colors` lookup.
    pub fn push_glycan_mesh(&mut self, mesh: &vv_core::PolytopeMesh, material: Material) {
        let m = self.materials.len() as u32;
        self.materials.push(material);
        let base = self.vertices.len() as u32;
        self.vertices.extend(
            mesh.positions
                .iter()
                .zip(&mesh.normals)
                .zip(&mesh.colors)
                .map(|((&position, &normal), &c)| TraceVertex {
                    position,
                    normal,
                    color: crate::color::rgba(c[0], c[1], c[2]),
                    material: m,
                }),
        );
        self.triangles.extend(
            (0..mesh.positions.len() as u32 / 3)
                .map(|t| [base + t * 3, base + t * 3 + 1, base + t * 3 + 2]),
        );
    }
}

/// A mesh vertex: shading normal, colour and material are interpolated
/// (normal, colour) or taken from the first vertex (material).
#[derive(Clone, Copy, Debug)]
pub struct TraceVertex {
    pub position: Vec3,
    pub normal: Vec3,
    pub color: u32,
    pub material: u32,
}

#[derive(Clone, Copy, Debug)]
pub struct TraceSphere {
    pub center: Vec3,
    pub radius: f32,
    /// Packed sRGB, as `colors_for` makes them.
    pub color: u32,
    pub material: u32,
}

/// A cylinder from `a` to `b`, each half coloured as its end (a bond
/// between two atoms), open-ended: the atoms' spheres or the neighbouring
/// bonds close it.
#[derive(Clone, Copy, Debug)]
pub struct TraceCylinder {
    pub a: Vec3,
    pub b: Vec3,
    pub radius: f32,
    pub color_a: u32,
    pub color_b: u32,
    pub material: u32,
}

/// How to render.
#[derive(Clone, Copy, Debug)]
pub struct TraceSettings {
    pub width: u32,
    pub height: u32,
    /// Samples per pixel: antialiasing, soft shadows and ambient
    /// occlusion all converge with it.
    pub samples: u32,
    pub lighting: Lighting,
    /// Bottom and top of the background gradient, linear RGB; ignored
    /// when `transparent`.
    pub background: [f32; 3],
    pub background_top: [f32; 3],
    /// Background pixels get alpha 0, edges partial alpha.
    pub transparent: bool,
    /// World-space clip plane, as `RenderSettings::clip`.
    pub clip: Option<[f32; 4]>,
    /// Every shadow and occlusion ray counted as reaching the light: the
    /// viewport's shading exactly, for checking the tracer against it.
    pub unshadowed: bool,
    /// Every live light's shadow cone half-angle (soft shadows), in
    /// radians; [`DEFAULT_LIGHT_SPREAD`] matches the viewport's near-hard
    /// shadows.
    pub light_spread: f32,
    /// Scales the sky term before `light_surface_seen` sees it, on top of
    /// `lighting.ambient`; 1.0 is parity with the viewport. [`RENDER_AO`]
    /// is the render default.
    pub ao: f32,
    /// Scales the key+fill terms, specular included; 1.0 is parity with
    /// the viewport. [`RENDER_DIRECT`] is the render default.
    pub direct: f32,
}

/// Default for [`TraceSettings::ao`]: doubles the current lighting
/// preset's ambient so crevices read with strong occlusion depth.
pub const RENDER_AO: f32 = 2.0;
/// Default for [`TraceSettings::direct`]: cuts the key+fill (and their
/// specular) to about a third, so they no longer wash out the
/// ambient-occlusion shading above.
pub const RENDER_DIRECT: f32 = 0.35;

/// Layout must match `TraceParams` in `shaders/path_trace.wgsl`.
#[repr(C)]
#[derive(Clone, Copy, Debug, Pod, Zeroable)]
struct TraceParams {
    view_inv: [[f32; 4]; 4],
    view: [[f32; 4]; 4],
    background: [f32; 4],
    background_top: [f32; 4],
    /// proj[0][0], proj[1][1], orthographic (0/1), transparent (0/1).
    lens: [f32; 4],
    size: [u32; 2],
    tile_origin: [u32; 2],
    /// This dispatch adds samples `first_sample..first_sample +
    /// sample_count` (the first overwrites); `resolve` divides by
    /// `first_sample + sample_count`.
    first_sample: u32,
    sample_count: u32,
    unshadowed: u32,
    /// Light cone half-angle (soft shadows), in radians.
    light_spread: f32,
    scene_radius: f32,
    /// Of every SES and skin surface (`TraceSurfaces`).
    probe_radius: f32,
    shrink: f32,
    ao_scale: f32,
    direct_scale: f32,
    _pad: [f32; 3],
}

/// Layout must match `Camera` in `shaders/path_trace.wgsl` (only what
/// `shading.wgsl` reads).
#[repr(C)]
#[derive(Clone, Copy, Debug, Pod, Zeroable)]
struct ClipUniform {
    clip: [f32; 4],
    cut: [f32; 4],
}

/// Pixels per tile side; one tile, all its samples, per submit.
const TILE: u32 = 256;
/// Samples per submit within a tile.
const SAMPLES_PER_SUBMIT: u32 = 16;
/// [`TraceSettings::light_spread`]'s default: about the sun's apparent
/// size, doubled -- the viewport's own near-hard shadows. Widening it
/// softens shadow edges but measured ~2x render time at 3x this angle
/// (4HHB cartoon, `Default` lighting); A/B'd up to 7x wider with no
/// visible penumbra on molecular scenes (blockers sit too close to what
/// they shade), so quality presets don't vary it -- see `vv-app`'s
/// `Quality::light_spread`.
pub const DEFAULT_LIGHT_SPREAD: f32 = 0.02;
/// Primitives per BVH leaf, at most.
const LEAF_MAX: usize = 4;
/// SAH bins per axis.
const BINS: usize = 12;

/// One BVH node as the shader reads it: a leaf holds primitives
/// `refs[first..first + count]`, an inner node its children at `first`
/// and `first + 1` (`count == 0`).
#[repr(C)]
#[derive(Clone, Copy, Debug, Pod, Zeroable)]
struct Node {
    min: [f32; 3],
    first: u32,
    max: [f32; 3],
    count: u32,
}

/// Kind of a primitive reference, in its top two bits (spheres 0).
const CYLINDER: u32 = 1 << 30;
const TRIANGLE: u32 = 2 << 30;
/// An SES or skin patch: a record of `TraceSurfaces`.
const PATCH: u32 = 3 << 30;

fn triangle_bounds(scene: &TraceScene, t: [u32; 3]) -> (Vec3, Vec3) {
    let p = t.map(|v| scene.vertices[v as usize].position);
    (p[0].min(p[1]).min(p[2]), p[0].max(p[1]).max(p[2]))
}

fn cylinder_bounds(c: &TraceCylinder) -> (Vec3, Vec3) {
    let r = Vec3::splat(c.radius);
    (c.a.min(c.b) - r, c.a.max(c.b) + r)
}

struct Prim {
    reference: u32,
    min: Vec3,
    max: Vec3,
    centroid: Vec3,
}

fn area(min: Vec3, max: Vec3) -> f32 {
    let e = (max - min).max(Vec3::ZERO);
    e.x * e.y + e.y * e.z + e.z * e.x
}

/// Top-down binned-SAH build (Wald 2007). Returns the flattened nodes
/// and the primitive references in leaf order.
fn build_bvh(scene: &TraceScene) -> (Vec<Node>, Vec<u32>) {
    let mut prims: Vec<Prim> = scene
        .spheres
        .iter()
        .enumerate()
        .map(|(k, s)| Prim {
            reference: k as u32,
            min: s.center - Vec3::splat(s.radius),
            max: s.center + Vec3::splat(s.radius),
            centroid: s.center,
        })
        .chain(scene.cylinders.iter().enumerate().map(|(k, c)| {
            let (min, max) = cylinder_bounds(c);
            Prim {
                reference: k as u32 | CYLINDER,
                min,
                max,
                centroid: (c.a + c.b) * 0.5,
            }
        }))
        .chain(scene.triangles.iter().enumerate().map(|(k, &t)| {
            let (min, max) = triangle_bounds(scene, t);
            Prim {
                reference: k as u32 | TRIANGLE,
                min,
                max,
                centroid: (min + max) * 0.5,
            }
        }))
        .chain(
            scene
                .surfaces
                .bounds
                .iter()
                .enumerate()
                .map(|(k, &(center, radius))| Prim {
                    reference: k as u32 | PATCH,
                    min: center - Vec3::splat(radius),
                    max: center + Vec3::splat(radius),
                    centroid: center,
                }),
        )
        .collect();
    let mut nodes = vec![Node::zeroed()];
    if prims.is_empty() {
        return (nodes, vec![0]);
    }
    // (node, start, end) still to split.
    let mut todo = vec![(0usize, 0usize, prims.len())];
    while let Some((node, start, end)) = todo.pop() {
        let span = &mut prims[start..end];
        let (min, max) = span.iter().fold(
            (Vec3::splat(f32::INFINITY), Vec3::splat(f32::NEG_INFINITY)),
            |(lo, hi), p| (lo.min(p.min), hi.max(p.max)),
        );
        nodes[node].min = min.to_array();
        nodes[node].max = max.to_array();
        let split = if span.len() > LEAF_MAX {
            best_split(span)
        } else {
            None
        };
        let Some((axis, cut)) = split else {
            nodes[node].first = start as u32;
            nodes[node].count = (end - start) as u32;
            continue;
        };
        let mut mid = start + partition(span, |p| p.centroid[axis] < cut);
        if mid == start || mid == end {
            // Rounding put every centroid on one side: split at the median.
            let half = span.len() / 2;
            span.select_nth_unstable_by(half, |a, b| a.centroid[axis].total_cmp(&b.centroid[axis]));
            mid = start + half;
        }
        let children = nodes.len();
        nodes.push(Node::zeroed());
        nodes.push(Node::zeroed());
        nodes[node].first = children as u32;
        nodes[node].count = 0;
        todo.push((children, start, mid));
        todo.push((children + 1, mid, end));
    }
    (nodes, prims.iter().map(|p| p.reference).collect())
}

/// Moves the primitives `keep` accepts to the front; returns how many.
fn partition(span: &mut [Prim], keep: impl Fn(&Prim) -> bool) -> usize {
    let mut front = 0;
    for k in 0..span.len() {
        if keep(&span[k]) {
            span.swap(front, k);
            front += 1;
        }
    }
    front
}

/// The cheapest (axis, centroid cut) by the surface area heuristic, or
/// `None` when no split beats a leaf (or the centroids coincide).
fn best_split(span: &[Prim]) -> Option<(usize, f32)> {
    let (lo, hi) = span.iter().fold(
        (Vec3::splat(f32::INFINITY), Vec3::splat(f32::NEG_INFINITY)),
        |(lo, hi), p| (lo.min(p.centroid), hi.max(p.centroid)),
    );
    let mut best: Option<(f32, usize, f32)> = None;
    for axis in 0..3 {
        let extent = hi[axis] - lo[axis];
        if extent <= 1e-6 {
            continue;
        }
        let mut count = [0usize; BINS];
        let mut bounds = [(Vec3::splat(f32::INFINITY), Vec3::splat(f32::NEG_INFINITY)); BINS];
        for p in span {
            let b = (((p.centroid[axis] - lo[axis]) / extent * BINS as f32) as usize).min(BINS - 1);
            count[b] += 1;
            bounds[b] = (bounds[b].0.min(p.min), bounds[b].1.max(p.max));
        }
        // Cost of cutting after bin k: left and right areas times counts.
        let mut right = [(0usize, 0.0f32); BINS];
        let (mut n, mut b) = (
            0,
            (Vec3::splat(f32::INFINITY), Vec3::splat(f32::NEG_INFINITY)),
        );
        for k in (1..BINS).rev() {
            n += count[k];
            b = (b.0.min(bounds[k].0), b.1.max(bounds[k].1));
            right[k] = (n, area(b.0, b.1));
        }
        let (mut n, mut b) = (
            0,
            (Vec3::splat(f32::INFINITY), Vec3::splat(f32::NEG_INFINITY)),
        );
        for k in 0..BINS - 1 {
            n += count[k];
            b = (b.0.min(bounds[k].0), b.1.max(bounds[k].1));
            let (rn, ra) = right[k + 1];
            if n == 0 || rn == 0 {
                continue;
            }
            let cost = n as f32 * area(b.0, b.1) + rn as f32 * ra;
            if best.is_none_or(|(c, _, _)| cost < c) {
                best = Some((cost, axis, lo[axis] + extent * (k + 1) as f32 / BINS as f32));
            }
        }
    }
    let (min, max) = span.iter().fold(
        (Vec3::splat(f32::INFINITY), Vec3::splat(f32::NEG_INFINITY)),
        |(lo, hi), p| (lo.min(p.min), hi.max(p.max)),
    );
    // A leaf costs one test per primitive; a split one traversal step
    // (about one primitive test) plus its children's expected tests.
    let leaf = span.len() as f32;
    best.filter(|(cost, _, _)| span.len() > 16 || 1.0 + cost / area(min, max) < leaf)
        .map(|(_, axis, cut)| (axis, cut))
}

/// The traced scene on the GPU, and the pipeline that renders it.
pub struct PathTracer {
    pipeline: wgpu::ComputePipeline,
    resolve: wgpu::ComputePipeline,
    bind: wgpu::BindGroup,
    params: wgpu::Buffer,
    lighting: wgpu::Buffer,
    clip: wgpu::Buffer,
    pixels: wgpu::Buffer,
    readback: wgpu::Buffer,
    /// Radius of a sphere around everything (occlusion rays stop there).
    scene_radius: f32,
    probe_radius: f32,
    shrink: f32,
}

impl PathTracer {
    /// Uploads `scene` and builds its hierarchy.
    pub fn new(ctx: &GpuContext, scene: &TraceScene) -> Result<Self, OutOfGpuMemory> {
        let device = &ctx.device;
        let (nodes, refs) = build_bvh(scene);
        let root = nodes[0];
        let scene_radius = (Vec3::from(root.max) - Vec3::from(root.min))
            .length()
            .max(1.0);
        let spheres: Vec<[f32; 4]> = scene
            .spheres
            .iter()
            .map(|s| s.center.extend(s.radius).to_array())
            .collect();
        let sphere_looks: Vec<[u32; 2]> = scene
            .spheres
            .iter()
            .map(|s| [s.color, s.material])
            .collect();
        let cylinders: Vec<[f32; 8]> = scene
            .cylinders
            .iter()
            .map(|c| {
                let [ax, ay, az] = c.a.to_array();
                let [bx, by, bz] = c.b.to_array();
                [ax, ay, az, c.radius, bx, by, bz, 0.0]
            })
            .collect();
        let cylinder_looks: Vec<[u32; 4]> = scene
            .cylinders
            .iter()
            .map(|c| [c.color_a, c.color_b, c.material, 0])
            .collect();
        let vertices: Vec<[f32; 8]> = scene
            .vertices
            .iter()
            .map(|v| {
                let [px, py, pz] = v.position.to_array();
                let [nx, ny, nz] = v.normal.to_array();
                [
                    px,
                    py,
                    pz,
                    f32::from_bits(v.color),
                    nx,
                    ny,
                    nz,
                    f32::from_bits(v.material),
                ]
            })
            .collect();
        let triangles: Vec<[u32; 4]> = scene
            .triangles
            .iter()
            .map(|t| [t[0], t[1], t[2], 0])
            .collect();
        let materials: Vec<Material> = if scene.materials.is_empty() {
            vec![Material::default()]
        } else {
            scene.materials.clone()
        };
        // An empty binding still needs one element of the largest kind (a
        // skin patch, 64 bytes). One bigger than a binding may be would
        // fail validation, not allocation: refused as out of memory too.
        let max_binding = device.limits().max_storage_buffer_binding_size;
        let storage = |label: &str, bytes: &[u8]| {
            let bytes = if bytes.is_empty() {
                &[0u8; 64][..]
            } else {
                bytes
            };
            if bytes.len() as u64 > max_binding {
                return Err(OutOfGpuMemory);
            }
            try_buffer_init(ctx, label, bytes, wgpu::BufferUsages::STORAGE)
        };
        let nodes = storage("trace nodes", bytemuck::cast_slice(&nodes))?;
        let refs = storage("trace refs", bytemuck::cast_slice(&refs))?;
        let spheres = storage("trace spheres", bytemuck::cast_slice(&spheres))?;
        let sphere_looks = storage("trace sphere looks", bytemuck::cast_slice(&sphere_looks))?;
        let cylinders = storage("trace cylinders", bytemuck::cast_slice(&cylinders))?;
        let cylinder_looks = storage(
            "trace cylinder looks",
            bytemuck::cast_slice(&cylinder_looks),
        )?;
        let materials = storage("trace materials", bytemuck::cast_slice(&materials))?;
        let vertices = storage("trace vertices", bytemuck::cast_slice(&vertices))?;
        let triangles = storage("trace triangles", bytemuck::cast_slice(&triangles))?;
        let surf = &scene.surfaces;
        let patches = storage("trace patches", bytemuck::cast_slice(&surf.records))?;
        let surface_atoms = storage("trace surface atoms", bytemuck::cast_slice(&surf.atoms))?;
        let surface_looks = storage("trace surface looks", bytemuck::cast_slice(&surf.looks))?;
        let caps = storage("trace caps", bytemuck::cast_slice(&surf.caps))?;
        let probes = storage("trace probes", bytemuck::cast_slice(&surf.probes))?;
        let probe_neighbors = storage(
            "trace probe neighbors",
            bytemuck::cast_slice(&surf.probe_neighbors),
        )?;
        let skin = storage("trace skin patches", bytemuck::cast_slice(&surf.skin))?;
        let competitor_starts = storage(
            "trace competitor starts",
            bytemuck::cast_slice(&surf.competitor_starts),
        )?;
        let competitors = storage("trace competitors", bytemuck::cast_slice(&surf.competitors))?;
        let uniform = |label: &str, bytes: &[u8]| {
            device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some(label),
                contents: bytes,
                usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            })
        };
        let params = uniform("trace params", bytemuck::bytes_of(&TraceParams::zeroed()));
        let lighting = uniform("trace lighting", bytemuck::bytes_of(&Lighting::default()));
        let clip = uniform("trace clip", bytemuck::bytes_of(&ClipUniform::zeroed()));
        let pixel_count = (TILE * TILE) as u64;
        let accum = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("trace accum"),
            size: pixel_count * 16,
            usage: wgpu::BufferUsages::STORAGE,
            mapped_at_creation: false,
        });
        let pixels = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("trace pixels"),
            size: pixel_count * 4,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
            mapped_at_creation: false,
        });
        let readback = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("trace readback"),
            size: pixel_count * 4,
            usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        let entry = |binding, ty| wgpu::BindGroupLayoutEntry {
            binding,
            visibility: wgpu::ShaderStages::COMPUTE,
            ty: wgpu::BindingType::Buffer {
                ty,
                has_dynamic_offset: false,
                min_binding_size: None,
            },
            count: None,
        };
        let read = wgpu::BufferBindingType::Storage { read_only: true };
        let write = wgpu::BufferBindingType::Storage { read_only: false };
        let layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("path trace"),
            entries: &[
                entry(0, wgpu::BufferBindingType::Uniform),
                entry(1, wgpu::BufferBindingType::Uniform),
                entry(2, wgpu::BufferBindingType::Uniform),
                entry(3, read),
                entry(4, read),
                entry(5, read),
                entry(6, read),
                entry(7, read),
                entry(8, read),
                entry(9, read),
                entry(10, write),
                entry(11, write),
                entry(12, read),
                entry(13, read),
                entry(14, read),
                entry(15, read),
                entry(16, read),
                entry(17, read),
                entry(18, read),
                entry(19, read),
                entry(20, read),
                entry(21, read),
                entry(22, read),
            ],
        });
        let buffers = [
            &params,
            &clip,
            &lighting,
            &nodes,
            &refs,
            &spheres,
            &sphere_looks,
            &cylinders,
            &cylinder_looks,
            &materials,
            &accum,
            &pixels,
            &vertices,
            &triangles,
            &patches,
            &surface_atoms,
            &surface_looks,
            &caps,
            &probes,
            &probe_neighbors,
            &skin,
            &competitor_starts,
            &competitors,
        ];
        let entries: Vec<wgpu::BindGroupEntry> = buffers
            .iter()
            .enumerate()
            .map(|(k, b)| wgpu::BindGroupEntry {
                binding: k as u32,
                resource: b.as_entire_binding(),
            })
            .collect();
        let bind = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("path trace"),
            layout: &layout,
            entries: &entries,
        });
        let module = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("path trace"),
            source: wgpu::ShaderSource::Wgsl(
                format!(
                    "{}\n{}\n{}\n{}",
                    include_str!("shaders/path_trace.wgsl"),
                    include_str!("shaders/shading.wgsl"),
                    include_str!("shaders/ses_patch.wgsl"),
                    include_str!("shaders/skin_patch.wgsl")
                )
                .into(),
            ),
        });
        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("path trace"),
            bind_group_layouts: &[Some(&layout)],
            immediate_size: 0,
        });
        let pipeline = |entry: &str| {
            device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
                label: Some("path trace"),
                layout: Some(&pipeline_layout),
                module: &module,
                entry_point: Some(entry),
                compilation_options: Default::default(),
                cache: None,
            })
        };
        Ok(Self {
            pipeline: pipeline("trace_tile"),
            resolve: pipeline("resolve"),
            bind,
            params,
            lighting,
            clip,

            pixels,
            readback,
            scene_radius,
            probe_radius: scene
                .surfaces
                .probe_radius
                .unwrap_or(vv_core::ses::WATER_PROBE),
            shrink: scene
                .surfaces
                .shrink
                .unwrap_or(vv_core::skin_surface::DEFAULT_SHRINK),
        })
    }

    /// Renders through `camera` into straight-alpha sRGB RGBA8, rows top
    /// to bottom. `progress` hears the fraction done after each tile.
    pub fn render(
        &self,
        ctx: &GpuContext,
        camera: &Camera,
        settings: &TraceSettings,
        mut progress: impl FnMut(f32),
    ) -> Vec<u8> {
        let (w, h) = (settings.width.max(1), settings.height.max(1));
        let proj = camera.proj(w as f32 / h as f32);
        let view = camera.view();
        ctx.queue
            .write_buffer(&self.lighting, 0, bytemuck::bytes_of(&settings.lighting));
        ctx.queue.write_buffer(
            &self.clip,
            0,
            bytemuck::bytes_of(&ClipUniform {
                clip: settings
                    .clip
                    .map_or(KEEP_ALL, |plane| camera.view_clip(plane)),
                cut: camera
                    .near_cut_plane()
                    .map_or(KEEP_ALL, |plane| camera.view_clip(plane)),
            }),
        );
        let mut image = vec![0u8; (w * h * 4) as usize];
        let samples = settings.samples.max(1);
        let tiles_x = w.div_ceil(TILE);
        let tiles_y = h.div_ceil(TILE);
        let total = (tiles_x * tiles_y) as f32;
        for ty in 0..tiles_y {
            for tx in 0..tiles_x {
                let origin = [tx * TILE, ty * TILE];
                let mut first = 0;
                while first < samples {
                    let count = SAMPLES_PER_SUBMIT.min(samples - first);
                    let params = TraceParams {
                        view_inv: view.inverse().to_cols_array_2d(),
                        view: view.to_cols_array_2d(),
                        background: Vec3::from(settings.background).extend(1.0).to_array(),
                        background_top: Vec3::from(settings.background_top).extend(1.0).to_array(),
                        lens: Vec4::new(
                            proj.x_axis.x,
                            proj.y_axis.y,
                            (camera.projection == Projection::Orthographic) as u32 as f32,
                            settings.transparent as u32 as f32,
                        )
                        .to_array(),
                        size: [w, h],
                        tile_origin: origin,
                        first_sample: first,
                        sample_count: count,
                        unshadowed: settings.unshadowed as u32,
                        light_spread: settings.light_spread,
                        scene_radius: self.scene_radius,
                        probe_radius: self.probe_radius,
                        shrink: self.shrink,
                        ao_scale: settings.ao,
                        direct_scale: settings.direct,
                        _pad: [0.0; 3],
                    };
                    ctx.queue
                        .write_buffer(&self.params, 0, bytemuck::bytes_of(&params));
                    let mut encoder = ctx.device.create_command_encoder(&Default::default());
                    {
                        let mut pass = encoder.begin_compute_pass(&Default::default());
                        pass.set_bind_group(0, &self.bind, &[]);
                        pass.set_pipeline(&self.pipeline);
                        pass.dispatch_workgroups(TILE / 8, TILE / 8, 1);
                        if first + count == samples {
                            pass.set_pipeline(&self.resolve);
                            pass.dispatch_workgroups(TILE / 8, TILE / 8, 1);
                        }
                    }
                    first += count;
                    if first == samples {
                        encoder.copy_buffer_to_buffer(
                            &self.pixels,
                            0,
                            &self.readback,
                            0,
                            (TILE * TILE * 4) as u64,
                        );
                    }
                    ctx.queue.submit([encoder.finish()]);
                }
                self.read_tile(ctx, origin, w, h, &mut image);
                progress((ty * tiles_x + tx + 1) as f32 / total);
            }
        }
        image
    }

    /// Copies the resolved tile at `origin` into `image` (`w` x `h`).
    fn read_tile(&self, ctx: &GpuContext, origin: [u32; 2], w: u32, h: u32, image: &mut [u8]) {
        let slice = self.readback.slice(..);
        slice.map_async(wgpu::MapMode::Read, |_| {});
        let _ = ctx.device.poll(wgpu::PollType::wait_indefinitely());
        {
            let data = slice.get_mapped_range().expect("mapped trace tile");
            for row in 0..TILE.min(h - origin[1]) {
                let cols = TILE.min(w - origin[0]) as usize;
                let src = (row * TILE * 4) as usize;
                let dst = (((origin[1] + row) * w + origin[0]) * 4) as usize;
                image[dst..dst + cols * 4].copy_from_slice(&data[src..src + cols * 4]);
            }
        }
        self.readback.unmap();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    const INDEX: u32 = (1 << 30) - 1;

    fn sphere(center: Vec3, radius: f32) -> TraceSphere {
        TraceSphere {
            center,
            radius,
            color: 0,
            material: 0,
        }
    }

    /// Every primitive is in exactly one leaf, inside its node's box.
    #[test]
    fn the_hierarchy_holds_every_primitive_once_inside_its_boxes() {
        let mut scene = TraceScene::default();
        for k in 0..1000 {
            let f = k as f32;
            scene.spheres.push(sphere(
                Vec3::new((f * 0.37).sin() * 30.0, (f * 0.11).cos() * 30.0, f * 0.05),
                1.0 + (k % 3) as f32 * 0.3,
            ));
        }
        scene.cylinders.push(TraceCylinder {
            a: Vec3::ZERO,
            b: Vec3::X * 5.0,
            radius: 0.2,
            color_a: 0,
            color_b: 0,
            material: 0,
        });
        for p in [Vec3::ZERO, Vec3::Y * 3.0, Vec3::Z * 2.0] {
            scene.vertices.push(TraceVertex {
                position: p,
                normal: Vec3::X,
                color: 0,
                material: 0,
            });
        }
        scene.triangles.push([0, 1, 2]);
        let (nodes, refs) = build_bvh(&scene);
        let mut seen = vec![0; refs.len()];
        let mut stack = vec![0usize];
        while let Some(n) = stack.pop() {
            let node = nodes[n];
            let (lo, hi) = (Vec3::from(node.min), Vec3::from(node.max));
            if node.count == 0 {
                for c in [node.first as usize, node.first as usize + 1] {
                    let child = nodes[c];
                    assert!(Vec3::from(child.min).cmpge(lo - 1e-4).all());
                    assert!(Vec3::from(child.max).cmple(hi + 1e-4).all());
                    stack.push(c);
                }
                continue;
            }
            assert!(node.count as usize <= 16);
            for k in node.first..node.first + node.count {
                seen[k as usize] += 1;
                let r = refs[k as usize];
                let (min, max) = match r & !INDEX {
                    CYLINDER => cylinder_bounds(&scene.cylinders[(r & INDEX) as usize]),
                    TRIANGLE => triangle_bounds(&scene, scene.triangles[(r & INDEX) as usize]),
                    _ => {
                        let s = scene.spheres[r as usize];
                        (s.center - s.radius, s.center + s.radius)
                    }
                };
                assert!(min.cmpge(lo - 1e-4).all() && max.cmple(hi + 1e-4).all());
            }
        }
        assert!(seen.iter().all(|&n| n == 1));
        let mut all: Vec<u32> = refs.clone();
        all.sort_unstable();
        all.dedup();
        assert_eq!(all.len(), 1002);
    }
}
