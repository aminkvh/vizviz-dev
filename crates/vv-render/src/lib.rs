//! GPU renderer built on wgpu.
//!
//! Atoms live in paged storage buffers. Each frame a compute pass culls
//! against the frustum and sorts survivors into two buckets by projected
//! size: ray-cast sphere impostors (quads) and single-pixel points. Both
//! draw through indirect arguments the compute pass filled in, so the CPU
//! never touches per-atom data after upload.

pub mod camera;
pub mod color;
pub mod coloring;
pub mod context;
pub mod lod;
pub mod occlusion;
pub mod path_trace;
pub mod postprocess;
pub mod renderer;
pub mod scene;
pub mod selection;
pub mod ses_surface;
pub mod ses_volume;
pub mod style;
pub mod timing;
mod trace_surfaces;

pub use camera::{Camera, Projection};
pub use color::ColorScheme;
pub use context::{GpuContext, GpuError};
pub use lod::AdaptiveLod;
pub use occlusion::OcclusionVolume;
pub use postprocess::downsample_2x_srgb;
pub use renderer::{
    AtomSizes, CartoonBindings, CartoonItem, CartoonMesh, DrawItem, FrameTimes,
    GaussianSurfaceItem, PatchSurface, PatchSurfaceItem, Pick, RenderSettings, Renderer,
    Representation, SesBindings, BOND_ID_FLAG,
};
pub use scene::{
    colors_for, colors_for_fragments, colors_for_ss, colors_from_scalar, colors_from_scalar_in,
    pack_rgba, scalar_range, CartoonGpu, CartoonParams, CartoonSectionGpu, DrawState,
    GaussianSurfaceGpu, GaussianSurfaceParams, GlycanGpu, GlycanVertexGpu, GpuStructure,
    MeshDisplay, OutOfGpuMemory, PageParams, SkinPatchGpu, SkinSurfaceGpu, SkinSurfaceParams, Soup,
};
pub use selection::ItemSelection;
pub use ses_surface::{SesGpu, SesLayout};
pub use style::{Light, Lighting, LightingPreset, Material, MaterialPreset, StylePreset};
pub use wgpu;

/// Human-readable list of GPU adapters wgpu can see on this machine.
pub fn enumerate_adapters() -> Vec<String> {
    let instance = GpuContext::instance();
    pollster::block_on(instance.enumerate_adapters(wgpu::Backends::all()))
        .into_iter()
        .map(|adapter| {
            let info = adapter.get_info();
            format!("{} [{:?}, {:?}]", info.name, info.backend, info.device_type)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    #[test]
    fn enumerating_adapters_does_not_panic() {
        // CI runners may have zero adapters; that is fine here.
        let _ = super::enumerate_adapters();
    }
}
