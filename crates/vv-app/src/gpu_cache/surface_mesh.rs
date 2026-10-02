//! A surface drawn as its triangle mesh in place of the solid: the mesh's
//! edges as lines of a set pixel width, or its vertices as dots. The mesh
//! is built on a worker thread (`vv_core::{gaussian_mesh, ses_mesh,
//! skin_mesh}`); picking and colouring go by each vertex's nearest atom.

use super::*;
use vv_core::cartoon::ExpandedMesh;
use vv_render::{MeshDisplay, Soup};

/// More triangles than this are refused with a notice rather than drawn.
const MAX_TRIANGLES: usize = 4_000_000;
/// A dot's radius, in angstroms, per pixel of `line` width.
const DOT_RADIUS_PER_WIDTH: f32 = 0.1;
/// The mesh's grid spacing for a small structure, in angstroms; it grows
/// with the cube root of the atom count so a larger structure meshes in
/// about the same time.
const VOXEL: f32 = 0.6;
const VOXEL_FULL_ATOMS: f32 = 20_000.0;

/// How a surface rep draws its mesh; the width is the `line` option.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) enum Style {
    Edges(f32),
    Dots(f32),
}

/// How `rep` displays its surface, `None` for the solid (and for any rep
/// that is not a surface).
pub(super) fn style(rep: &Rep) -> Option<Style> {
    let width = rep.option("line")?;
    match rep.option("surface")? as u32 {
        1 => Some(Style::Edges(width)),
        2 => Some(Style::Dots(width)),
        _ => None,
    }
}

/// What a mesh was built from: it stays valid while this does not change,
/// whatever the style drawing it.
#[derive(Clone, PartialEq)]
struct Key {
    kind: SceneRepresentation,
    /// Probe radius, blobbiness or shrink.
    parameter: f32,
    atoms: Atoms,
    frame: usize,
}

pub(super) enum Drawn {
    Edges {
        gpu: vv_render::GlycanGpu,
        bindings: vv_render::CartoonBindings,
    },
    Dots {
        points: Derived,
        /// Angstroms.
        radius: f32,
    },
}

pub(super) struct SurfaceMesh {
    key: Key,
    style: Style,
    mesh: Arc<ExpandedMesh>,
    pub(super) drawn: Drawn,
}

fn surface_parameter(rep: &Rep) -> f32 {
    let name = match rep.representation {
        SceneRepresentation::Ses => "probe",
        SceneRepresentation::SkinSurface => "shrink",
        _ => "blob",
    };
    rep.option(name).expect("every surface has its parameter")
}

/// The mesh of the surface `kind` of `positions`, on a worker.
fn spawn_mesh(
    kind: SceneRepresentation,
    positions: Vec<Vec3>,
    radii: Vec<f32>,
    parameter: f32,
    frame: usize,
    waker: Option<Waker>,
) -> Job<ExpandedMesh> {
    let count = positions.len();
    let voxel = VOXEL * (count as f32 / VOXEL_FULL_ATOMS).cbrt().max(1.0);
    spawn_job(frame, waker, move || {
        timed("surface mesh", count, || match kind {
            SceneRepresentation::Ses => {
                vv_core::ses_mesh::mesh(&positions, &radii, parameter, voxel)
            }
            SceneRepresentation::SkinSurface => {
                vv_core::skin_mesh::mesh(&positions, &radii, parameter, voxel)
            }
            _ => vv_core::gaussian_mesh::mesh(
                &positions,
                &radii,
                parameter,
                vv_render::scene::GAUSSIAN_EPSILON,
                voxel,
            ),
        })
    })
}

impl SurfaceMesh {
    fn upload(
        ctx: &GpuContext,
        renderer: &Renderer,
        key: Key,
        style: Style,
        mesh: Arc<ExpandedMesh>,
        all_colors: &[u32],
    ) -> Result<Self, OutOfGpuMemory> {
        let real = |local: u32| key.atoms.as_ref().map_or(local, |a| a[local as usize]);
        let drawn = match style {
            Style::Edges(width) => {
                let (positions, normals, colors, source) = soup(&mesh, real, all_colors);
                let gpu = vv_render::GlycanGpu::upload_soup(
                    ctx,
                    &Soup {
                        positions: &positions,
                        normals: &normals,
                        colors: &colors,
                        source: Arc::new(source),
                    },
                    MeshDisplay::Lines(width),
                )?;
                let bindings = renderer.bind_glycan(&gpu);
                Drawn::Edges { gpu, bindings }
            }
            Style::Dots(width) => {
                let atoms: Vec<u32> = mesh.source.iter().map(|&s| real(s)).collect();
                let colors: Vec<u32> = atoms.iter().map(|&a| all_colors[a as usize]).collect();
                let radii = vec![0.0; mesh.positions.len()];
                let gpu = GpuStructure::from_parts(ctx, &mesh.positions, &radii, &colors, &[]);
                Drawn::Dots {
                    points: Derived {
                        bindings: renderer.bind(&gpu),
                        gpu,
                        atom_map: Arc::new(atoms),
                        bond_atoms: Arc::new(Vec::new()),
                        strands: None,
                    },
                    radius: width * DOT_RADIUS_PER_WIDTH,
                }
            }
        };
        Ok(Self {
            key,
            style,
            mesh,
            drawn,
        })
    }

    fn valid_for(&self, key: &Key) -> bool {
        self.key == *key
    }

    /// The sizes a dot piece of `radius` angstroms draws with.
    pub(super) fn dot_sizes(radius: f32) -> AtomSizes {
        AtomSizes {
            radius_scale: 0.0,
            radius_offset: radius,
            bond_radius: 0.0,
        }
    }
}

/// `mesh`'s triangles as a soup: three vertices each, coloured and
/// attributed (`real` maps a mesh atom to the structure's) by their
/// nearest atom.
fn soup(
    mesh: &ExpandedMesh,
    real: impl Fn(u32) -> u32,
    all_colors: &[u32],
) -> (Vec<Vec3>, Vec<Vec3>, Vec<u32>, Vec<u32>) {
    let n = mesh.indices.len();
    let (mut positions, mut normals) = (Vec::with_capacity(n), Vec::with_capacity(n));
    let (mut colors, mut source) = (Vec::with_capacity(n), Vec::with_capacity(n));
    for &v in &mesh.indices {
        let v = v as usize;
        let atom = real(mesh.source[v]);
        positions.push(mesh.positions[v]);
        normals.push(mesh.normals[v]);
        colors.push(all_colors[atom as usize]);
        source.push(atom);
    }
    (positions, normals, colors, source)
}

impl GpuCache {
    /// The mesh display of a surface rep: a finished mesh job is uploaded,
    /// a new one starts when the surface itself changed (the old mesh stays
    /// up meanwhile), and a change of style or colouring alone redraws from
    /// the mesh already built.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn sync_surface_mesh(
        &mut self,
        r: &mut RepEntry,
        loaded: &LoadedStructure,
        rep: &Rep,
        style: Style,
        (frame, recolor): (usize, bool),
        ctx: &GpuContext,
        renderer: &Renderer,
        waker: &Option<Waker>,
        colors: &mut dyn FnMut() -> Vec<u32>,
    ) {
        let structure = &loaded.structure;
        let atoms = match rep.representation {
            SceneRepresentation::Ses => without_water(structure, &r.atoms),
            _ => r.atoms.clone(),
        };
        let key = Key {
            kind: rep.representation,
            parameter: surface_parameter(rep),
            atoms,
            frame,
        };
        if let Some(done) = poll_job(&mut r.mesh_job) {
            match done {
                Ok((built, mesh)) => {
                    let key = Key {
                        frame: built,
                        ..key.clone()
                    };
                    self.finish_mesh(r, ctx, renderer, key, style, Arc::new(mesh), &colors());
                }
                Err(()) => self.messages.push("surface mesh build failed".into()),
            }
        }
        let have = match &r.geometry {
            RepGeometry::Mesh(m) => Some(m),
            _ => None,
        };
        match have {
            Some(m) if m.valid_for(&key) => {
                if m.style != style || recolor {
                    let (key, mesh) = (m.key.clone(), m.mesh.clone());
                    self.finish_mesh(r, ctx, renderer, key, style, mesh, &colors());
                }
            }
            _ if r.mesh_job.is_none() && !r.failed => {
                self.start_mesh(r, loaded, rep, &key, frame, waker, colors);
            }
            _ => {}
        }
    }

    fn start_mesh(
        &mut self,
        r: &mut RepEntry,
        loaded: &LoadedStructure,
        rep: &Rep,
        key: &Key,
        frame: usize,
        waker: &Option<Waker>,
        colors: &mut dyn FnMut() -> Vec<u32>,
    ) {
        let (positions, radii, _) = atom_arrays(&loaded.structure, &key.atoms, frame, &colors());
        if rep.representation == SceneRepresentation::SkinSurface
            && positions.len() > SKIN_MAX_ATOMS
        {
            r.failed = true;
            self.messages.push(format!(
                "skin surface is limited to {SKIN_MAX_ATOMS} atoms for now ({} has {} in this rep)",
                loaded.label,
                positions.len()
            ));
            return;
        }
        if matches!(r.geometry, RepGeometry::Pending) {
            self.infos.push(format!(
                "meshing the surface of {} ({} atoms); it appears when ready",
                loaded.label,
                positions.len()
            ));
        }
        r.mesh_job = Some(spawn_mesh(
            rep.representation,
            positions,
            radii,
            key.parameter,
            frame,
            waker.clone(),
        ));
    }

    fn finish_mesh(
        &mut self,
        r: &mut RepEntry,
        ctx: &GpuContext,
        renderer: &Renderer,
        key: Key,
        style: Style,
        mesh: Arc<ExpandedMesh>,
        all_colors: &[u32],
    ) {
        if mesh.indices.len() / 3 > MAX_TRIANGLES {
            r.failed = true;
            r.geometry = RepGeometry::Pending;
            self.messages.push(format!(
                "the surface mesh has {} triangles, more than the {MAX_TRIANGLES} it can draw; \
                 use solid or a smaller selection",
                mesh.indices.len() / 3
            ));
            return;
        }
        match SurfaceMesh::upload(ctx, renderer, key, style, mesh, all_colors) {
            Ok(mesh) => r.geometry = RepGeometry::Mesh(mesh),
            Err(OutOfGpuMemory) => {
                r.failed = true;
                r.geometry = RepGeometry::Pending;
                self.messages
                    .push("not enough GPU memory for the surface mesh".into());
            }
        }
    }
}
