//! The solvent-excluded surface as a volume: `vv_core::ses::Ses` baked
//! into the same distance-field volume the Gaussian surface's march draws
//! (`shaders/ses_volume.wgsl`), so it costs what that march costs --
//! one texture fetch per step, empty bricks skipped -- instead of a
//! ray-cast billboard per patch. Accurate to about a voxel, where the
//! patch renderer (`ses_surface`) is exact.

use bytemuck::{Pod, Zeroable};
use glam::Vec3;
use wgpu::util::DeviceExt;

use crate::context::GpuContext;
use crate::scene::{
    try_buffer_init, with_voxel_budget, GaussianSurfaceGpu, GaussianSurfaceParams, OutOfGpuMemory,
    VolumeBox, GAUSSIAN_MAX_VOXELS,
};
use crate::style::Material;

/// Layout must match `Params` in `shaders/ses_volume.wgsl`.
#[repr(C)]
#[derive(Clone, Copy, Debug, Pod, Zeroable)]
struct SesVolumeParams {
    vol_origin: [f32; 3],
    voxel: f32,
    dims: [u32; 3],
    probe: f32,
    brick_dims: [u32; 3],
    clamp: f32,
    grid_origin: [f32; 3],
    cell_size: f32,
    grid_dims: [i32; 3],
    _pad: u32,
}

/// Voxel slices per submit: small enough that one submit never runs long
/// enough for the OS to reset an integrated GPU.
const SLAB: u32 = 16;

/// The bake's compute pipeline, built once per device
/// (`GpuContext::ses_volume_bake`).
pub(crate) struct SesVolumeBake {
    pipeline: wgpu::ComputePipeline,
    layout: wgpu::BindGroupLayout,
    slab_layout: wgpu::BindGroupLayout,
}

impl SesVolumeBake {
    pub(crate) fn new(device: &wgpu::Device) -> Self {
        let entry = |binding, ty| wgpu::BindGroupLayoutEntry {
            binding,
            visibility: wgpu::ShaderStages::COMPUTE,
            ty,
            count: None,
        };
        let buffer = |binding, ty| {
            entry(
                binding,
                wgpu::BindingType::Buffer {
                    ty,
                    has_dynamic_offset: false,
                    min_binding_size: None,
                },
            )
        };
        let read = wgpu::BufferBindingType::Storage { read_only: true };
        let layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("ses volume bake"),
            entries: &[
                buffer(0, read),
                buffer(1, read),
                buffer(2, wgpu::BufferBindingType::Uniform),
                buffer(3, read),
                buffer(4, read),
                entry(
                    5,
                    wgpu::BindingType::StorageTexture {
                        access: wgpu::StorageTextureAccess::WriteOnly,
                        format: wgpu::TextureFormat::Rgba16Float,
                        view_dimension: wgpu::TextureViewDimension::D3,
                    },
                ),
                buffer(6, wgpu::BufferBindingType::Storage { read_only: false }),
                buffer(7, read),
                buffer(8, read),
                buffer(9, read),
            ],
        });
        let slab_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("ses volume slab"),
            entries: &[buffer(0, wgpu::BufferBindingType::Uniform)],
        });
        let module = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("ses volume bake"),
            source: wgpu::ShaderSource::Wgsl(include_str!("shaders/ses_volume.wgsl").into()),
        });
        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("ses volume bake"),
            bind_group_layouts: &[Some(&layout), Some(&slab_layout)],
            immediate_size: 0,
        });
        let pipeline = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: Some("ses volume bake"),
            layout: Some(&pipeline_layout),
            module: &module,
            entry_point: Some("bake"),
            compilation_options: Default::default(),
            cache: None,
        });
        Self {
            pipeline,
            layout,
            slab_layout,
        }
    }
}

impl GaussianSurfaceGpu {
    /// `ses` (built from `positions` and van der Waals `radii`) baked into
    /// a volume the Gaussian surface's march draws, each point coloured as
    /// `colors[atom]` of the atom it lies on or nearest. Records and submits
    /// the bake; out of GPU memory it retries with coarser voxels, as
    /// [`GaussianSurfaceGpu::upload`] does.
    pub fn ses(
        ctx: &GpuContext,
        ses: &vv_core::ses::Ses,
        positions: &[Vec3],
        radii: &[f32],
        colors: &[u32],
    ) -> Result<Self, OutOfGpuMemory> {
        assert_eq!(positions.len(), radii.len());
        assert_eq!(positions.len(), colors.len());
        assert!(!positions.is_empty(), "need at least one atom");
        with_voxel_budget(GAUSSIAN_MAX_VOXELS, |budget| {
            Self::ses_with_budget(ctx, ses, positions, radii, colors, budget)
        })
    }

    fn ses_with_budget(
        ctx: &GpuContext,
        ses: &vv_core::ses::Ses,
        positions: &[Vec3],
        radii: &[f32],
        colors: &[u32],
        max_voxels: f64,
    ) -> Result<Self, OutOfGpuMemory> {
        let device = &ctx.device;
        let probe = ses.probe_radius;
        let vdw_max = radii.iter().copied().fold(0.0, f32::max);
        // The excluded volume lies inside the SAS.
        let grid_box = VolumeBox::around(ctx, positions, vdw_max + probe, max_voxels);
        // Exact within two voxels of the surface: what the march's
        // trilinear samples and gradient read around a crossing.
        let clamp = probe + 2.0 * grid_box.voxel;
        let cell_size = vdw_max + probe + clamp;
        let indices: Vec<u32> = (0..positions.len() as u32).collect();
        let grid = vv_core::Grid::build(positions, &indices, cell_size);

        let atoms: Vec<[f32; 4]> = positions
            .iter()
            .zip(radii)
            .map(|(p, &r)| [p.x, p.y, p.z, r])
            .collect();
        let mut torus_starts = vec![0u32; positions.len() + 1];
        for t in &ses.tori {
            torus_starts[t.atoms[0] as usize + 1] += 1;
        }
        for i in 0..positions.len() {
            torus_starts[i + 1] += torus_starts[i];
        }
        // `build` sorts tori by atoms, so each atom's are contiguous.
        let tori: Vec<[u32; 4]> = ses
            .tori
            .iter()
            .map(|t| [t.atoms[1], t.probes[0], t.probes[1], 0])
            .collect();
        let probes: Vec<[f32; 4]> = ses
            .probes
            .iter()
            .map(|p| p.center.extend(0.0).to_array())
            .collect();
        let storage = |label: &str, bytes: &[u8]| {
            let bytes = if bytes.is_empty() {
                &[0u8; 16][..]
            } else {
                bytes
            };
            try_buffer_init(ctx, label, bytes, wgpu::BufferUsages::STORAGE)
        };
        let atoms_buffer = storage("ses volume atoms", bytemuck::cast_slice(&atoms))?;
        let colors_buffer = storage("ses volume colors", bytemuck::cast_slice(colors))?;
        let cell_starts = storage(
            "ses volume cell starts",
            bytemuck::cast_slice(grid.cell_starts()),
        )?;
        let cell_atoms = storage(
            "ses volume cell atoms",
            bytemuck::cast_slice(grid.cell_atoms()),
        )?;
        let torus_starts = storage(
            "ses volume torus starts",
            bytemuck::cast_slice(&torus_starts),
        )?;
        let tori = storage("ses volume tori", bytemuck::cast_slice(&tori))?;
        let probes = storage("ses volume probes", bytemuck::cast_slice(&probes))?;
        let bake_params = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("ses volume params"),
            contents: bytemuck::bytes_of(&SesVolumeParams {
                vol_origin: grid_box.min.to_array(),
                voxel: grid_box.voxel,
                dims: grid_box.dims.to_array(),
                probe,
                brick_dims: grid_box.brick_dims.to_array(),
                clamp,
                grid_origin: grid.origin().to_array(),
                cell_size: grid.cell_size(),
                grid_dims: grid.dims().to_array(),
                _pad: 0,
            }),
            usage: wgpu::BufferUsages::UNIFORM,
        });
        let (volume, volume_view, bricks, sampler) = grid_box.allocate(ctx)?;

        let bake = ctx
            .ses_volume_bake
            .get_or_init(|| SesVolumeBake::new(device));
        let buffers = [
            (0, &atoms_buffer),
            (1, &colors_buffer),
            (2, &bake_params),
            (3, &cell_starts),
            (4, &cell_atoms),
            (6, &bricks),
            (7, &torus_starts),
            (8, &tori),
            (9, &probes),
        ];
        let mut entries: Vec<wgpu::BindGroupEntry> = buffers
            .iter()
            .map(|&(binding, b)| wgpu::BindGroupEntry {
                binding,
                resource: b.as_entire_binding(),
            })
            .collect();
        entries.push(wgpu::BindGroupEntry {
            binding: 5,
            resource: wgpu::BindingResource::TextureView(&volume_view),
        });
        let bind = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("ses volume bake"),
            layout: &bake.layout,
            entries: &entries,
        });
        let dims = grid_box.dims;
        for z0 in (0..dims.z).step_by(SLAB as usize) {
            let slab = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some("ses volume slab"),
                contents: bytemuck::cast_slice(&[z0, 0, 0, 0]),
                usage: wgpu::BufferUsages::UNIFORM,
            });
            let slab_bind = device.create_bind_group(&wgpu::BindGroupDescriptor {
                label: Some("ses volume slab"),
                layout: &bake.slab_layout,
                entries: &[wgpu::BindGroupEntry {
                    binding: 0,
                    resource: slab.as_entire_binding(),
                }],
            });
            let mut encoder = device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
                label: Some("ses volume bake"),
            });
            {
                let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                    label: Some("ses volume bake"),
                    timestamp_writes: None,
                });
                pass.set_pipeline(&bake.pipeline);
                pass.set_bind_group(0, &bind, &[]);
                pass.set_bind_group(1, &slab_bind, &[]);
                pass.dispatch_workgroups(
                    dims.x.div_ceil(4),
                    dims.y.div_ceil(4),
                    SLAB.min(dims.z - z0).div_ceil(4),
                );
            }
            ctx.queue.submit([encoder.finish()]);
        }

        // What the march reads.
        let base_params = GaussianSurfaceParams {
            view_inv: glam::Mat4::IDENTITY.to_cols_array_2d(),
            vol_origin: grid_box.min.to_array(),
            voxel: grid_box.voxel,
            dims: dims.to_array(),
            isovalue: 1.0,
            brick_dims: grid_box.brick_dims.to_array(),
            blob_factor: 0.0,
            material: Material::default(),
        };
        let params = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("ses volume march params"),
            contents: bytemuck::bytes_of(&base_params),
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
        });
        Ok(Self {
            volume,
            volume_view,
            sampler,
            bricks,
            params,
            base_params,
            atom_count: positions.len() as u32,
            bounds_center: grid_box.center(),
            bounds_radius: grid_box.radius(),
            inputs: None,
        })
    }
}
