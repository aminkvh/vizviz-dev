//! Frame orchestration: cull -> opaque -> depth pyramid.
//!
//! The depth pyramid built at the end of frame N feeds occlusion culling in
//! frame N+1, so hidden atoms are dropped in the compute pass instead of
//! being rasterized and rejected by the late depth test.

use std::sync::Arc;

use bytemuck::{Pod, Zeroable};
use glam::{Vec2, Vec3};
use wgpu::util::DeviceExt;

use crate::camera::{Camera, CameraUniform, KEEP_ALL};
use crate::context::GpuContext;
use crate::scene::{
    CartoonGpu, CartoonParams, DrawState, GaussianSurfaceGpu, GaussianSurfaceParams, GlycanGpu,
    GpuStructure, Page, PageParams, PageState, SkinSurfaceGpu, SkinSurfaceParams,
    BONDS_INDIRECT_OFFSET, CARTOON_INDIRECT_RESET, CLUSTER_ATOMS, INDIRECT_RESET,
    LATE_INDIRECT_OFFSET, LINES_INDIRECT_OFFSET, PAGE_ATOMS, POINTS_INDIRECT_OFFSET,
};
use crate::selection::{ItemSelection, SelectionOutline};
use crate::ses_surface::{SesGpu, SesParams};
use crate::style::{Lighting, Material, StylePreset};
use crate::timing::{GpuTimer, PassTimes};

pub const COLOR_FORMAT: wgpu::TextureFormat = wgpu::TextureFormat::Rgba8UnormSrgb;
pub const NORMAL_FORMAT: wgpu::TextureFormat = wgpu::TextureFormat::Rgba8Unorm;
pub const DEPTH_FORMAT: wgpu::TextureFormat = wgpu::TextureFormat::Depth32Float;
/// Glass pass: weighted premultiplied colour and alpha
/// (`shaders/shading.wgsl`'s `glass()`).
const ACCUM_FORMAT: wgpu::TextureFormat = wgpu::TextureFormat::Rgba16Float;
/// Glass pass: the fraction of light that gets through every layer.
const REVEAL_FORMAT: wgpu::TextureFormat = wgpu::TextureFormat::R8Unorm;
/// Glass prepass: (`id_base`, patch + 1) of the surface patch nearest
/// each pixel (`fs_skin_surface_owner`, `fs_ses_surface_owner`).
const OWNER_FORMAT: wgpu::TextureFormat = wgpu::TextureFormat::Rg32Uint;
/// Picking target: a frame-wide atom id **plus one**, or a bond id plus
/// one with [`BOND_ID_FLAG`] set (see `fs_cylinder`). `0` means "nothing
/// drawn here". Zero rather than `u32::MAX` because the attachment is
/// cleared through a float clear color, and only `0.0` converts to the
/// same integer on every backend: Vulkan turned `u32::MAX as f64` into
/// `0xffffffff`, DX12 did not, and background pixels picked atom 0 there.
/// Resolved back to a [`Pick`] by `Renderer::pick`.
pub const ID_FORMAT: wgpu::TextureFormat = wgpu::TextureFormat::R32Uint;
const NO_PICK: u32 = 0;
/// Set on a pick id that names a bond; the low 31 bits are the frame-wide
/// bond id. Must match `BOND_ID_FLAG` in `shaders/draw.wgsl`.
pub const BOND_ID_FLAG: u32 = 0x8000_0000;
const HIZ_FORMAT: wgpu::TextureFormat = wgpu::TextureFormat::R32Float;
const CULL_WORKGROUP: u32 = 256;
const HIZ_WORKGROUP: u32 = 8;
/// Outline tuning (see `shaders/outline.wgsl`). A depth jump of 0.6 A is
/// about a third of a carbon's van der Waals radius: enough that the depth
/// ramp inside one sphere stays quiet at typical zoom, small enough that
/// one atom in front of another still gets a line.
const OUTLINE_DEPTH_THRESHOLD: f32 = 0.6;
/// Normals more than ~60 degrees apart count as a crease.
const OUTLINE_NORMAL_THRESHOLD: f32 = 0.5;
const OUTLINE_STRENGTH: f32 = 0.85;
/// Ambient-occlusion radius in Angstroms (`shaders/ao.wgsl` clamps it to
/// 10..96 px on screen): a few atoms, so grooves and pockets darken.
const AO_RADIUS: f32 = 24.0;
/// The AO target: one channel.
const AO_FORMAT: wgpu::TextureFormat = wgpu::TextureFormat::Rg8Unorm;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum Representation {
    /// Atoms at full van der Waals radius, no bonds.
    #[default]
    Spacefill,
    /// Small spheres joined by bond cylinders.
    BallAndStick,
    /// A smooth tube through the backbone: spheres and cylinders of
    /// `vv_core::TUBE_RADIUS` over spline samples (`vv_core::backbone`),
    /// uploaded as their own `GpuStructure` via `GpuStructure::from_parts`.
    /// The stored radius already is the tube radius, so no scaling.
    Tube,
    /// The solvent-accessible surface as spheres: every atom grown by the
    /// 1.4 A water probe.
    Sas,
    /// Thin bonds with balls of the same radius (licorice).
    Sticks,
    /// Bonds as 1-px lines.
    Lines,
}

/// The water probe radius the SAS style grows atoms by (Angstrom).
pub const SAS_PROBE_RADIUS: f32 = 1.4;
/// Sticks ball and stick radius (Angstrom).
pub const STICKS_RADIUS: f32 = 0.3;

/// A real-triangle mesh for [`Renderer::render_all`]'s `cartoons` list —
/// not a [`DrawItem`]: there is no meaningful `radius_scale`/`bond_radius`
/// for a mesh, so this was never going to fit the impostor
/// `Representation` enum without a nonsense answer to either. Two kinds
/// share the list because they share everything else about how a mesh
/// composites into the frame (bounds, pick ids, the opaque pass): a
/// ribbon/arrow (`vv_core::cartoon`), built and culled on the GPU from
/// stored cross-sections, and a glycan's SNFG glyphs
/// (`vv_core::glycan::mesh`), a plain triangle soup rebuilt on the CPU
/// every frame and uploaded whole.
#[derive(Clone, Copy)]
pub enum CartoonMesh<'a> {
    Ribbon(&'a CartoonGpu),
    Glycan(&'a GlycanGpu),
}

impl<'a> CartoonMesh<'a> {
    fn bounds(self) -> (Vec3, f32) {
        match self {
            CartoonMesh::Ribbon(g) => (g.bounds_center, g.bounds_radius),
            CartoonMesh::Glycan(g) => (g.bounds_center, g.bounds_radius),
        }
    }

    /// The pick id space this mesh claims: one id per section (ribbon)
    /// or per vertex (glycan, flat-shaded so a vertex names one glyph).
    fn pick_count(self) -> u32 {
        match self {
            CartoonMesh::Ribbon(g) => g.section_count,
            CartoonMesh::Glycan(g) => g.vertex_count,
        }
    }

    /// This frame's [`CartoonParams`] destination.
    fn params(self) -> &'a wgpu::Buffer {
        match self {
            CartoonMesh::Ribbon(g) => &g.params,
            CartoonMesh::Glycan(g) => &g.params,
        }
    }
}

/// One structure's mesh geometry for [`Renderer::render_all`]'s
/// `cartoons` list.
pub struct CartoonItem<'a> {
    pub mesh: CartoonMesh<'a>,
    pub bindings: &'a CartoonBindings,
    pub material: Material,
}

/// A mesh's bind groups for drawing, and — a ribbon only — for its
/// level-of-detail compute pass and its per-frame section rebuild.
pub enum CartoonBindings {
    Ribbon {
        lod: wgpu::BindGroup,
        draw: wgpu::BindGroup,
        frame: wgpu::BindGroup,
    },
    Glycan {
        draw: wgpu::BindGroup,
    },
}

impl CartoonBindings {
    fn draw(&self) -> &wgpu::BindGroup {
        match self {
            CartoonBindings::Ribbon { draw, .. } | CartoonBindings::Glycan { draw } => draw,
        }
    }
}

/// One structure's Gaussian-surface data (`vv_render::scene::
/// GaussianSurfaceGpu`) for [`Renderer::render_all`]'s separate
/// `gaussian_surfaces` list, drawn as a full-screen ray-march pass
/// rather than per-atom impostors -- see `shaders/gaussian_surface.wgsl`.
pub struct GaussianSurfaceItem<'a> {
    pub gpu: &'a GaussianSurfaceGpu,
    pub bindings: &'a wgpu::BindGroup,
    pub material: Material,
}

/// A surface drawn as one ray-cast billboard per patch
/// (`shaders/billboard.wgsl`), with its bind groups: the skin surface's
/// mixed-cell patches (all drawn) or the SES's convex, toroidal and
/// concave ones (culled like atoms, [`SesBindings`]).
#[derive(Clone, Copy)]
pub enum PatchSurface<'a> {
    Skin(&'a SkinSurfaceGpu, &'a wgpu::BindGroup),
    Ses(&'a SesGpu, &'a SesBindings),
}

/// An SES's bind groups: the surface itself, the cull state of its patch
/// bounds, and per page the survivors the vertex stage reads.
pub struct SesBindings {
    surface: wgpu::BindGroup,
    pages: PageBindings,
    views: Vec<wgpu::BindGroup>,
}

impl PatchSurface<'_> {
    pub fn patch_count(self) -> u32 {
        match self {
            PatchSurface::Skin(g, _) => g.patch_count,
            PatchSurface::Ses(g, _) => g.patch_count,
        }
    }

    pub fn atom_count(self) -> u32 {
        match self {
            PatchSurface::Skin(g, _) => g.atom_count,
            PatchSurface::Ses(g, _) => g.atom_count,
        }
    }

    pub fn bounds(self) -> (glam::Vec3, f32) {
        match self {
            PatchSurface::Skin(g, _) => (g.bounds_center, g.bounds_radius),
            PatchSurface::Ses(g, _) => (g.bounds_center, g.bounds_radius),
        }
    }

    fn bind_group(&self) -> &wgpu::BindGroup {
        match self {
            PatchSurface::Skin(_, b) => b,
            PatchSurface::Ses(_, b) => &b.surface,
        }
    }

    /// This frame's inverse view, pick id base and material.
    fn write_params(
        self,
        queue: &wgpu::Queue,
        view_inv: glam::Mat4,
        id_base: u32,
        material: Material,
    ) {
        match self {
            PatchSurface::Skin(g, _) => queue.write_buffer(
                &g.params,
                0,
                bytemuck::bytes_of(&SkinSurfaceParams {
                    view_inv: view_inv.to_cols_array_2d(),
                    shrink: g.shrink,
                    atom_count: g.atom_count,
                    patch_count: g.patch_count,
                    id_base,
                    material,
                }),
            ),
            PatchSurface::Ses(g, _) => queue.write_buffer(
                &g.params,
                0,
                bytemuck::bytes_of(&SesParams {
                    view_inv: view_inv.to_cols_array_2d(),
                    probe_radius: g.probe_radius,
                    atom_count: g.atom_count,
                    patch_count: g.patch_count,
                    id_base,
                    material,
                }),
            ),
        }
    }
}

/// One structure's patch surface for [`Renderer::render_all`]. Pickable:
/// a hit names the atom its patch belongs to (nearest, where it spans
/// several).
pub struct PatchSurfaceItem<'a> {
    pub surface: PatchSurface<'a>,
    pub material: Material,
}

/// The three pipelines a patch surface kind draws with.
struct PatchPipelines {
    opaque: wgpu::RenderPipeline,
    glass: wgpu::RenderPipeline,
    /// Glass prepass: which patch owns each pixel's nearest front face.
    owner: wgpu::RenderPipeline,
    /// The clip plane's solid cross-section, one atom-sphere disc at a
    /// time (`vs_atom_cap`/`fs_atom_cap`): opaque only, see `draw_ses`.
    atom_cap: wgpu::RenderPipeline,
}

impl Representation {
    pub fn radius_scale(self) -> f32 {
        match self {
            Representation::Spacefill | Representation::Tube | Representation::Sas => 1.0,
            Representation::BallAndStick => 0.25,
            Representation::Sticks | Representation::Lines => 0.0,
        }
    }

    /// Added to every drawn atom radius, after `radius_scale`.
    pub fn radius_offset(self) -> f32 {
        match self {
            Representation::Sas => SAS_PROBE_RADIUS,
            Representation::Sticks => STICKS_RADIUS,
            _ => 0.0,
        }
    }

    /// Bond cylinder radius in angstroms; `0` means no bonds are drawn.
    pub fn bond_radius(self) -> f32 {
        match self {
            Representation::Spacefill | Representation::Sas => 0.0,
            Representation::BallAndStick => 0.15,
            Representation::Tube => vv_core::TUBE_RADIUS,
            Representation::Sticks => STICKS_RADIUS,
            // Bounds for culling only: lines are one pixel wide.
            Representation::Lines => 0.1,
        }
    }

    /// The drawn radius of an atom of van der Waals radius `vdw`.
    pub fn atom_radius(self, vdw: f32) -> f32 {
        vdw * self.radius_scale() + self.radius_offset()
    }

    /// The next representation in the cycle spacefill -> ball-and-stick
    /// -> tube.
    pub fn toggle(self) -> Self {
        match self {
            Representation::Spacefill => Representation::BallAndStick,
            Representation::BallAndStick => Representation::Tube,
            _ => Representation::Spacefill,
        }
    }
}

#[derive(Clone, Debug)]
pub struct RenderSettings {
    pub representation: Representation,
    /// Atoms whose projected radius is below this many pixels draw as points.
    pub quad_px_threshold: f32,
    /// Drop atoms hidden behind what the previous frame drew.
    pub occlusion_culling: bool,
    /// Global lighting for every item (`shaders/shading.wgsl`).
    pub lighting: Lighting,
    /// The material of [`Renderer::render`]'s one item; `render_all`
    /// items carry their own.
    pub material: Material,
    pub background: wgpu::Color,
    /// Colour at the top of the frame for a vertical gradient from
    /// `background` (at the bottom); `None` is a solid background.
    pub background_top: Option<wgpu::Color>,
    /// Smooths impostor silhouettes in the *displayed* image
    /// (`Renderer::display_view`) every frame. Doesn't affect
    /// `read_color`/screenshots: SSAA export supersedes it there with
    /// real subpixel detail instead of FXAA's single-sample edge guess.
    pub fxaa: bool,
    /// Dark lines on depth and normal discontinuities (silhouettes and
    /// the creases between overlapping atoms), applied before FXAA and
    /// included in `read_color`/screenshots. `StylePreset::outline`
    /// suggests a default per preset.
    pub outline: bool,
    /// Outline line width in pixels. An export rendered at 2x for SSAA
    /// should double this so the lines downsample to the on-screen width.
    pub outline_width: u32,
    /// Selection halo colour, linear RGBA (see [`Renderer::set_selection`]).
    pub selection_color: [f32; 4],
    /// Selection halo width in target pixels; scale it with any
    /// supersampling so it keeps its on-screen width.
    pub selection_width: u32,
    /// Screen-space ambient occlusion strength, 0 = off
    /// (`shaders/ao.wgsl`). `StylePreset::ao` suggests a default.
    pub ao: f32,
    /// Depth cue: how far the far side of the scene fades toward the
    /// background, 0 = off. `StylePreset::depth_cue` suggests a default.
    pub depth_cue: f32,
    /// Screen-space shadows from the key light, 0 = off, 1 = full
    /// (`shaders/ao.wgsl`). `LightingPreset::shadows` suggests a default.
    pub shadows: f32,
    /// Depth of field, 0 = off, 1 = full: sharp at the rotation centre
    /// (the camera's pivot, else its target), blurring with distance from
    /// it (`shaders/dof.wgsl`).
    pub dof: f32,
    /// World-space clip plane `[nx, ny, nz, d]`: points with `n.p + d < 0`
    /// are cut away, and spheres and the Gaussian surface show their
    /// cross-section there. `n` must be unit length.
    pub clip: Option<[f32; 4]>,
    /// A sphere around everything drawn, for the depth cue's range. When
    /// `None`, `render_all` takes it from the items themselves.
    pub scene_bounds: Option<(glam::Vec3, f32)>,
    /// Glass self-occlusion culling (docs/RENDERING.md's "Transparency"):
    /// an approximation that can visibly darken a densely overlapping
    /// glass scene (dropped deep layers), so it defaults off and is
    /// meant for a camera actively moving, not a still view -- `false`
    /// for screenshots/exports and `Renderer::render`'s callers unless
    /// they opt in. `false` costs nothing beyond the existing opaque
    /// occlusion test; see `culled_by_depth` in `shaders/cull.wgsl`.
    pub fast_glass: bool,
}

impl Default for RenderSettings {
    fn default() -> Self {
        Self {
            representation: Representation::Spacefill,
            quad_px_threshold: 1.0,
            occlusion_culling: true,
            lighting: StylePreset::default().lighting(),
            material: StylePreset::default().material(),
            background: StylePreset::default().background(),
            background_top: None,
            fxaa: true,
            outline: StylePreset::default().outline(),
            outline_width: 1,
            selection_color: [0.0, 1.0, 0.0, 1.0],
            selection_width: 2,
            ao: StylePreset::default().ao(),
            depth_cue: StylePreset::default().depth_cue(),
            shadows: StylePreset::default().lighting_preset().shadows(),
            dof: 0.0,
            clip: None,
            scene_bounds: None,
            fast_glass: false,
        }
    }
}

/// Uniform for the post passes; layout must match `Post` in
/// `shaders/ao.wgsl` and `shaders/outline.wgsl`.
#[repr(C)]
#[derive(Clone, Copy, Debug, Pod, Zeroable)]
struct PostUniform {
    /// Clip -> view, to rebuild view-space positions from depth (both
    /// projections).
    inv_proj: [[f32; 4]; 4],
    outline_enabled: u32,
    outline_width: u32,
    near: f32,
    /// Minimum view-depth jump (Angstroms) that counts as an edge; grows
    /// with distance in the shader.
    depth_threshold: f32,
    /// `dot(n_center, n_neighbour)` below this is a crease.
    normal_threshold: f32,
    /// How much of the pixel's colour the line removes (1 = black).
    outline_strength: f32,
    ao_strength: f32,
    ao_radius: f32,
    fog_strength: f32,
    /// View depths where the depth cue starts and reaches full strength.
    fog_near: f32,
    fog_far: f32,
    /// `proj[1][1]`, to turn a world radius into pixels.
    proj_y: f32,
    /// Linear RGBA the depth cue fades toward.
    background: [f32; 4],
    orthographic: u32,
    shadow_strength: f32,
    _pad: [u32; 2],
    /// Linear RGBA at the top of the frame; the background fades from
    /// `background` at the bottom to this.
    background_top: [f32; 4],
    /// View-space direction toward the key light (w unused).
    light_dir: [f32; 4],
    /// View -> world, to look up the occlusion volume.
    view_inv: [[f32; 4]; 4],
    /// Occlusion volume corner (xyz) and voxel size (w).
    vol_origin: [f32; 4],
    /// Occlusion volume size (xyz); w = 1 when there is one.
    vol_extent: [f32; 4],
    /// `RenderSettings::clip`, or a plane that keeps everything.
    clip_world: [f32; 4],
    /// `Camera::near_cut_plane`, in the same form.
    cut_world: [f32; 4],
}

/// Uniform for the depth-of-field pass; layout must match `Dof` in
/// `shaders/dof.wgsl`.
#[repr(C)]
#[derive(Clone, Copy, Debug, Pod, Zeroable)]
struct DofUniform {
    inv_proj: [[f32; 4]; 4],
    /// View depth that stays sharp.
    focus: f32,
    /// Largest blur radius in pixels.
    max_radius: f32,
    /// Depth distance from `focus` at which the blur is largest.
    range: f32,
    _pad: f32,
}

/// Largest depth-of-field blur at full strength, in pixels of a 1080-line
/// frame (scaled with the target's height).
const DOF_MAX_RADIUS: f32 = 14.0;

pub type FrameTimes = PassTimes;

struct Layouts {
    frame: wgpu::BindGroupLayout,
    cull: wgpu::BindGroupLayout,
    draw: wgpu::BindGroupLayout,
    hiz_copy: wgpu::BindGroupLayout,
    hiz_down: wgpu::BindGroupLayout,
    blit: wgpu::BindGroupLayout,
    outline: wgpu::BindGroupLayout,
    ao: wgpu::BindGroupLayout,
    cartoon: wgpu::BindGroupLayout,
    cartoon_lod: wgpu::BindGroupLayout,
    cartoon_frame: wgpu::BindGroupLayout,
    glycan: wgpu::BindGroupLayout,
    gaussian_surface: wgpu::BindGroupLayout,
    skin_surface: wgpu::BindGroupLayout,
    ses_surface: wgpu::BindGroupLayout,
    ses_page: wgpu::BindGroupLayout,
    glass_resolve: wgpu::BindGroupLayout,
    glass_owner: wgpu::BindGroupLayout,
    occlusion: wgpu::BindGroupLayout,
    dof: wgpu::BindGroupLayout,
}

struct Targets {
    width: u32,
    height: u32,
    color: wgpu::Texture,
    color_view: wgpu::TextureView,
    normal_view: wgpu::TextureView,
    depth_view: wgpu::TextureView,
    id_texture: wgpu::Texture,
    id_view: wgpu::TextureView,
    /// Outline pass output: `color` with edge lines when enabled, a plain
    /// copy otherwise. What `read_color` (screenshots) reads.
    outlined: wgpu::Texture,
    outlined_view: wgpu::TextureView,
    /// Reads `color_view`, `normal_view`, `depth_view`, `ao_view`;
    /// writes `outlined_view`.
    outline_bind_group: wgpu::BindGroup,
    /// Ambient occlusion (`shaders/ao.wgsl`), read by the outline pass.
    ao_view: wgpu::TextureView,
    /// Reads `normal_view`, `depth_view`; writes `ao_view`.
    ao_bind_group: wgpu::BindGroup,
    /// Post-process output: what's actually displayed (`Renderer::
    /// display_view`). Always populated, even when FXAA is off (a plain
    /// copy of `outlined_view`) — see the comment on `display_view` for why.
    fxaa_color: wgpu::Texture,
    fxaa_view: wgpu::TextureView,
    /// `fxaa_color` read as plain `Rgba8Unorm`: egui samples its textures
    /// as already gamma-encoded, so an sRGB view would decode them twice.
    fxaa_gamma_view: wgpu::TextureView,
    /// Reads `outlined_view`, writes `fxaa_view`.
    fxaa_bind_group: wgpu::BindGroup,
    /// Reads `normal_view`, for where caps were drawn.
    cap_depth_bind_group: wgpu::BindGroup,
    /// Depth of field reads this copy of `outlined` and writes `outlined`.
    dof_source: wgpu::Texture,
    dof_bind_group: wgpu::BindGroup,
    /// Glass pass targets; resolved onto `outlined_view`.
    depth: wgpu::Texture,
    /// The opaque depth, then the nearest glass skin surface: the depth
    /// buffer of the glass prepass.
    glass_depth: wgpu::Texture,
    glass_depth_view: wgpu::TextureView,
    /// Reads `glass_depth`, writes the glass depth pyramid `frame_bind_
    /// group` exposes to cull.wgsl's `glass_hidden` -- this frame's copy,
    /// right after the "glass id" pass, before `glass_depth` is
    /// overwritten as Skin/SES scratch: next frame's glass self-occlusion
    /// cull reference.
    glass_hiz_copy_bind_group: wgpu::BindGroup,
    glass_owner_view: wgpu::TextureView,
    glass_owner_bind_group: wgpu::BindGroup,
    accum_view: wgpu::TextureView,
    reveal_view: wgpu::TextureView,
    glass_resolve_bind_group: wgpu::BindGroup,
    hiz_mip_count: u32,
    /// Level 0 copies the depth buffer, level i min-reduces level i-1.
    hiz_bind_groups: Vec<wgpu::BindGroup>,
    /// Camera uniform + full pyramid view, shared by every pass.
    frame_bind_group: wgpu::BindGroup,
}

pub struct Renderer {
    ctx: Arc<GpuContext>,
    layouts: Layouts,
    targets: Targets,
    camera_buffer: wgpu::Buffer,
    lighting_buffer: wgpu::Buffer,
    outline_buffer: wgpu::Buffer,
    dof_buffer: wgpu::Buffer,
    dof_pipeline: wgpu::RenderPipeline,
    cull_pipeline: wgpu::ComputePipeline,
    cull_bonds_pipeline: wgpu::ComputePipeline,
    /// Occlusion phase 2 (`PHASE = 2` in shaders/cull.wgsl).
    cull_late_pipeline: wgpu::ComputePipeline,
    cull_bonds_late_pipeline: wgpu::ComputePipeline,
    cull_clusters_pipeline: wgpu::ComputePipeline,
    cull_clusters_late_pipeline: wgpu::ComputePipeline,
    sphere_pipeline: wgpu::RenderPipeline,
    sphere_glass_pipeline: wgpu::RenderPipeline,
    sphere_pick_pipeline: wgpu::RenderPipeline,
    point_pipeline: wgpu::RenderPipeline,
    point_glass_pipeline: wgpu::RenderPipeline,
    point_pick_pipeline: wgpu::RenderPipeline,
    cylinder_pipeline: wgpu::RenderPipeline,
    cylinder_glass_pipeline: wgpu::RenderPipeline,
    cylinder_pick_pipeline: wgpu::RenderPipeline,
    /// Bonds below half a pixel, as 1-px lines.
    line_pipeline: wgpu::RenderPipeline,
    line_glass_pipeline: wgpu::RenderPipeline,
    line_pick_pipeline: wgpu::RenderPipeline,
    cartoon_pipeline: wgpu::RenderPipeline,
    cartoon_glass_pipeline: wgpu::RenderPipeline,
    cartoon_pick_pipeline: wgpu::RenderPipeline,
    cartoon_lod_pipeline: wgpu::ComputePipeline,
    cartoon_frame_pipeline: wgpu::ComputePipeline,
    /// `scene::cartoon_join_indices`, shared by every cartoon.
    cartoon_indices: wgpu::Buffer,
    glycan_pipeline: wgpu::RenderPipeline,
    glycan_glass_pipeline: wgpu::RenderPipeline,
    glycan_pick_pipeline: wgpu::RenderPipeline,
    gaussian_surface_pipeline: wgpu::RenderPipeline,
    gaussian_glass_pipeline: wgpu::RenderPipeline,
    skin: PatchPipelines,
    ses: PatchPipelines,
    glass_resolve_pipeline: wgpu::RenderPipeline,
    hiz_copy_pipeline: wgpu::ComputePipeline,
    hiz_down_pipeline: wgpu::ComputePipeline,
    outline_pipeline: wgpu::RenderPipeline,
    ao_pipeline: wgpu::RenderPipeline,
    fxaa_pipeline: wgpu::RenderPipeline,
    /// Moves caps' depth onto their cut plane (`shaders/cap_depth.wgsl`).
    cap_depth_pipeline: wgpu::RenderPipeline,
    passthrough_pipeline: wgpu::RenderPipeline,
    blit_sampler: wgpu::Sampler,
    /// Long-range AO and shadows (`set_occlusion_volume`); `None` leaves
    /// only the screen-space terms.
    occlusion: Option<Arc<crate::occlusion::OcclusionVolume>>,
    occlusion_sampler: wgpu::Sampler,
    /// Bound when there is no volume (the shader then skips it).
    empty_volume: crate::occlusion::OcclusionVolume,
    occlusion_bind_group: wgpu::BindGroup,
    timer: Option<GpuTimer>,
    last_times: Option<PassTimes>,
    /// Id ranges the most recent frame assigned to each draw item.
    pick_ranges: Vec<PickRange>,
    selection: SelectionOutline,
}

/// Per-page bind groups, created once per (renderer, structure) pair.
/// One draw state's bind groups over a structure's pages
/// (`Renderer::bind`, `Renderer::bind_state`), plus handles to the
/// per-frame buffers the renderer rewrites for that state.
pub struct PageBindings {
    pages: Vec<PageBinding>,
}

struct PageBinding {
    cull: wgpu::BindGroup,
    draw: wgpu::BindGroup,
    params: wgpu::Buffer,
    indirect: wgpu::Buffer,
    occluded: wgpu::Buffer,
}

impl PageBinding {
    /// This frame's parameters, and empty survivor counts for the cull.
    fn reset(&self, queue: &wgpu::Queue, encoder: &mut wgpu::CommandEncoder, params: &PageParams) {
        queue.write_buffer(&self.params, 0, bytemuck::bytes_of(params));
        queue.write_buffer(&self.indirect, 0, bytemuck::cast_slice(&INDIRECT_RESET));
        encoder.clear_buffer(&self.occluded, 0, None);
    }
}

/// Pick ids at and above this belong to nothing: an SES patch drawn as a
/// point names no atom (its fragments pick through `SesParams::id_base`).
const NO_PICK_ID_BASE: u32 = 1 << 30;

/// One structure in a composited frame; see [`Renderer::render_all`].
/// Every field is a reference or already `Copy`, so this is too, letting
/// a caller reuse one item across several `render_all` calls (e.g. a
/// test comparing a scene with and without an extra item) without
/// juggling borrows.
#[derive(Clone, Copy)]
pub struct DrawItem<'a> {
    pub structure: &'a GpuStructure,
    pub bindings: &'a PageBindings,
    pub representation: Representation,
    pub material: Material,
    /// Atom and bond sizes; `AtomSizes::of(representation)` unless a rep
    /// tunes them.
    pub sizes: AtomSizes,
}

/// How big an item draws its atoms and bonds (see `Representation::
/// atom_radius`): each atom's van der Waals radius times `radius_scale`
/// plus `radius_offset`, bonds `bond_radius` (0: none).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AtomSizes {
    pub radius_scale: f32,
    pub radius_offset: f32,
    pub bond_radius: f32,
}

impl AtomSizes {
    /// What `representation` draws with untouched.
    pub fn of(representation: Representation) -> Self {
        Self {
            radius_scale: representation.radius_scale(),
            radius_offset: representation.radius_offset(),
            bond_radius: representation.bond_radius(),
        }
    }

    pub fn atom_radius(&self, vdw: f32) -> f32 {
        vdw * self.radius_scale + self.radius_offset
    }
}

/// What [`Renderer::pick`] found under a pixel. `item` indexes the list
/// passed to the most recent `render_all` (always `0` after `render`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Pick {
    Atom {
        item: usize,
        atom: u32,
    },
    /// `bond` indexes the `BondTable` the structure was uploaded with.
    /// Only ball-and-stick draws cylinders, so only it can produce this.
    Bond {
        item: usize,
        bond: u32,
    },
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct PickRange {
    pub(crate) atom_base: u32,
    pub(crate) atom_count: u32,
    pub(crate) bond_base: u32,
    pub(crate) bond_count: u32,
}

pub(crate) fn uniform_entry(
    binding: u32,
    visibility: wgpu::ShaderStages,
) -> wgpu::BindGroupLayoutEntry {
    wgpu::BindGroupLayoutEntry {
        binding,
        visibility,
        ty: wgpu::BindingType::Buffer {
            ty: wgpu::BufferBindingType::Uniform,
            has_dynamic_offset: false,
            min_binding_size: None,
        },
        count: None,
    }
}

pub(crate) fn storage_entry(
    binding: u32,
    read_only: bool,
    visibility: wgpu::ShaderStages,
) -> wgpu::BindGroupLayoutEntry {
    wgpu::BindGroupLayoutEntry {
        binding,
        visibility,
        ty: wgpu::BindingType::Buffer {
            ty: wgpu::BufferBindingType::Storage { read_only },
            has_dynamic_offset: false,
            min_binding_size: None,
        },
        count: None,
    }
}

fn occlusion_bind_group(
    device: &wgpu::Device,
    layouts: &Layouts,
    view: &wgpu::TextureView,
    sampler: &wgpu::Sampler,
) -> wgpu::BindGroup {
    device.create_bind_group(&wgpu::BindGroupDescriptor {
        label: Some("occlusion volume"),
        layout: &layouts.occlusion,
        entries: &[
            wgpu::BindGroupEntry {
                binding: 0,
                resource: wgpu::BindingResource::TextureView(view),
            },
            wgpu::BindGroupEntry {
                binding: 1,
                resource: wgpu::BindingResource::Sampler(sampler),
            },
        ],
    })
}

pub(crate) fn texture_entry(
    binding: u32,
    sample_type: wgpu::TextureSampleType,
    visibility: wgpu::ShaderStages,
) -> wgpu::BindGroupLayoutEntry {
    wgpu::BindGroupLayoutEntry {
        binding,
        visibility,
        ty: wgpu::BindingType::Texture {
            sample_type,
            view_dimension: wgpu::TextureViewDimension::D2,
            multisampled: false,
        },
        count: None,
    }
}

fn storage_texture_entry(binding: u32) -> wgpu::BindGroupLayoutEntry {
    wgpu::BindGroupLayoutEntry {
        binding,
        visibility: wgpu::ShaderStages::COMPUTE,
        ty: wgpu::BindingType::StorageTexture {
            access: wgpu::StorageTextureAccess::WriteOnly,
            format: HIZ_FORMAT,
            view_dimension: wgpu::TextureViewDimension::D2,
        },
        count: None,
    }
}

/// A full-screen pass writing only depth, wherever its fragment says.
fn cap_depth_pipeline(
    device: &wgpu::Device,
    layouts: &Layouts,
    source: &str,
) -> wgpu::RenderPipeline {
    let module = shader(device, "cap depth", source.to_string());
    let layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
        label: Some("cap depth"),
        bind_group_layouts: &[Some(&layouts.frame), Some(&layouts.blit)],
        ..Default::default()
    });
    device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
        label: Some("cap depth"),
        layout: Some(&layout),
        vertex: wgpu::VertexState {
            module: &module,
            entry_point: Some("vs_fullscreen"),
            compilation_options: Default::default(),
            buffers: &[],
        },
        primitive: Default::default(),
        depth_stencil: Some(wgpu::DepthStencilState {
            format: DEPTH_FORMAT,
            depth_write_enabled: Some(true),
            depth_compare: Some(wgpu::CompareFunction::Always),
            stencil: Default::default(),
            bias: Default::default(),
        }),
        multisample: Default::default(),
        fragment: Some(wgpu::FragmentState {
            module: &module,
            entry_point: Some("fs_cap_depth"),
            compilation_options: Default::default(),
            targets: &[],
        }),
        multiview_mask: None,
        cache: None,
    })
}

/// `source` with the fragment `entries` marked for the early depth test
/// where the adapter supports it. Only for impostors drawn on the near
/// side of what they ray-cast: their depth is never nearer than the
/// billboard's (`less_equal` in reversed Z), so a hidden fragment is
/// rejected before its ray is solved.
fn with_early_depth(ctx: &GpuContext, mut source: String, entries: &[&str]) -> String {
    if ctx.early_depth {
        for entry in entries {
            source = source.replace(
                &format!("@fragment\nfn {entry}("),
                &format!("@fragment @early_depth_test(less_equal)\nfn {entry}("),
            );
        }
    }
    source
}

fn shader(device: &wgpu::Device, label: &str, source: String) -> wgpu::ShaderModule {
    device.create_shader_module(wgpu::ShaderModuleDescriptor {
        label: Some(label),
        source: wgpu::ShaderSource::Wgsl(source.into()),
    })
}

fn compute_pipeline(
    device: &wgpu::Device,
    label: &str,
    layouts: &[&wgpu::BindGroupLayout],
    module: &wgpu::ShaderModule,
    entry_point: &str,
) -> wgpu::ComputePipeline {
    let groups: Vec<Option<&wgpu::BindGroupLayout>> = layouts.iter().map(|l| Some(*l)).collect();
    let layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
        label: Some(label),
        bind_group_layouts: &groups,
        ..Default::default()
    });
    device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
        label: Some(label),
        layout: Some(&layout),
        module,
        entry_point: Some(entry_point),
        compilation_options: Default::default(),
        cache: None,
    })
}

impl Renderer {
    /// Offscreen targets of `width` x `height`; the app shows the result
    /// through egui (`display_view`), headless users read it back.
    pub fn new(ctx: Arc<GpuContext>, width: u32, height: u32) -> Self {
        let device = &ctx.device;
        let compute = wgpu::ShaderStages::COMPUTE;
        let vertex = wgpu::ShaderStages::VERTEX;

        let camera_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("camera"),
            contents: bytemuck::bytes_of(&CameraUniform::zeroed()),
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
        });
        let lighting_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("lighting"),
            contents: bytemuck::bytes_of(&Lighting::default()),
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
        });
        let dof_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("dof"),
            contents: bytemuck::bytes_of(&DofUniform::zeroed()),
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
        });
        let outline_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("outline"),
            contents: bytemuck::bytes_of(&PostUniform::zeroed()),
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
        });

        let layouts = Layouts {
            frame: device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
                label: Some("frame"),
                entries: &[
                    wgpu::BindGroupLayoutEntry {
                        binding: 0,
                        visibility: compute | wgpu::ShaderStages::VERTEX_FRAGMENT,
                        ty: wgpu::BindingType::Buffer {
                            ty: wgpu::BufferBindingType::Uniform,
                            has_dynamic_offset: false,
                            min_binding_size: None,
                        },
                        count: None,
                    },
                    texture_entry(
                        1,
                        wgpu::TextureSampleType::Float { filterable: false },
                        compute,
                    ),
                    // Style material params (see shaders/draw.wgsl); only
                    // read by the draw fragment shaders, but harmless to
                    // include in the cull/Hi-Z pipelines' layouts too,
                    // since a bind group layout may have entries a given
                    // shader doesn't reference.
                    wgpu::BindGroupLayoutEntry {
                        binding: 2,
                        visibility: wgpu::ShaderStages::FRAGMENT,
                        ty: wgpu::BindingType::Buffer {
                            ty: wgpu::BufferBindingType::Uniform,
                            has_dynamic_offset: false,
                            min_binding_size: None,
                        },
                        count: None,
                    },
                    // Last frame's nearest-glass depth, one texel per
                    // pixel; only cull.wgsl's `glass_hidden` reads it.
                    texture_entry(
                        3,
                        wgpu::TextureSampleType::Float { filterable: false },
                        compute,
                    ),
                ],
            }),
            // Binding 7 in both page layouts is the `PageParams` uniform
            // (shaders/atoms.wgsl); the numbers in between differ per
            // layout because cull and draw see different index buffers.
            cull: device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
                label: Some("cull page"),
                entries: &[
                    storage_entry(0, true, compute),
                    storage_entry(1, true, compute),
                    storage_entry(2, false, compute),
                    storage_entry(3, false, compute),
                    storage_entry(4, false, compute),
                    storage_entry(5, true, compute),
                    storage_entry(6, false, compute),
                    uniform_entry(7, compute),
                    storage_entry(8, false, compute),
                    storage_entry(9, false, compute),
                    storage_entry(10, true, compute),
                ],
            }),
            draw: device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
                label: Some("draw page"),
                entries: &[
                    storage_entry(0, true, vertex),
                    storage_entry(1, true, vertex),
                    storage_entry(2, true, vertex),
                    storage_entry(3, true, vertex),
                    storage_entry(4, true, vertex),
                    storage_entry(5, true, vertex),
                    storage_entry(6, true, vertex),
                    uniform_entry(7, wgpu::ShaderStages::VERTEX_FRAGMENT),
                ],
            }),
            hiz_copy: device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
                label: Some("hiz copy"),
                entries: &[
                    texture_entry(0, wgpu::TextureSampleType::Depth, compute),
                    storage_texture_entry(1),
                ],
            }),
            hiz_down: device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
                label: Some("hiz down"),
                entries: &[
                    texture_entry(
                        0,
                        wgpu::TextureSampleType::Float { filterable: false },
                        compute,
                    ),
                    storage_texture_entry(1),
                ],
            }),
            dof: device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
                label: Some("dof"),
                entries: &[
                    texture_entry(
                        0,
                        wgpu::TextureSampleType::Float { filterable: false },
                        wgpu::ShaderStages::FRAGMENT,
                    ),
                    texture_entry(
                        1,
                        wgpu::TextureSampleType::Depth,
                        wgpu::ShaderStages::FRAGMENT,
                    ),
                    uniform_entry(2, wgpu::ShaderStages::FRAGMENT),
                ],
            }),
            occlusion: device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
                label: Some("occlusion volume"),
                entries: &[
                    wgpu::BindGroupLayoutEntry {
                        binding: 0,
                        visibility: wgpu::ShaderStages::FRAGMENT,
                        ty: wgpu::BindingType::Texture {
                            sample_type: wgpu::TextureSampleType::Float { filterable: true },
                            view_dimension: wgpu::TextureViewDimension::D3,
                            multisampled: false,
                        },
                        count: None,
                    },
                    wgpu::BindGroupLayoutEntry {
                        binding: 1,
                        visibility: wgpu::ShaderStages::FRAGMENT,
                        ty: wgpu::BindingType::Sampler(wgpu::SamplerBindingType::Filtering),
                        count: None,
                    },
                ],
            }),
            glass_owner: device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
                label: Some("glass owner"),
                entries: &[texture_entry(
                    0,
                    wgpu::TextureSampleType::Uint,
                    wgpu::ShaderStages::FRAGMENT,
                )],
            }),
            glass_resolve: device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
                label: Some("glass resolve"),
                entries: &[
                    texture_entry(
                        0,
                        wgpu::TextureSampleType::Float { filterable: false },
                        wgpu::ShaderStages::FRAGMENT,
                    ),
                    texture_entry(
                        1,
                        wgpu::TextureSampleType::Float { filterable: false },
                        wgpu::ShaderStages::FRAGMENT,
                    ),
                ],
            }),
            blit: device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
                label: Some("blit"),
                entries: &[
                    texture_entry(
                        0,
                        wgpu::TextureSampleType::Float { filterable: true },
                        wgpu::ShaderStages::FRAGMENT,
                    ),
                    wgpu::BindGroupLayoutEntry {
                        binding: 1,
                        visibility: wgpu::ShaderStages::FRAGMENT,
                        ty: wgpu::BindingType::Sampler(wgpu::SamplerBindingType::Filtering),
                        count: None,
                    },
                ],
            }),
            outline: device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
                label: Some("outline"),
                entries: &[
                    texture_entry(
                        0,
                        wgpu::TextureSampleType::Float { filterable: false },
                        wgpu::ShaderStages::FRAGMENT,
                    ),
                    texture_entry(
                        1,
                        wgpu::TextureSampleType::Float { filterable: false },
                        wgpu::ShaderStages::FRAGMENT,
                    ),
                    texture_entry(
                        2,
                        wgpu::TextureSampleType::Depth,
                        wgpu::ShaderStages::FRAGMENT,
                    ),
                    uniform_entry(3, wgpu::ShaderStages::FRAGMENT),
                    texture_entry(
                        4,
                        wgpu::TextureSampleType::Float { filterable: false },
                        wgpu::ShaderStages::FRAGMENT,
                    ),
                ],
            }),
            ao: device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
                label: Some("ao"),
                entries: &[
                    texture_entry(
                        0,
                        wgpu::TextureSampleType::Float { filterable: false },
                        wgpu::ShaderStages::FRAGMENT,
                    ),
                    texture_entry(
                        1,
                        wgpu::TextureSampleType::Depth,
                        wgpu::ShaderStages::FRAGMENT,
                    ),
                    uniform_entry(2, wgpu::ShaderStages::FRAGMENT),
                ],
            }),
            cartoon: device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
                label: Some("cartoon"),
                entries: &[
                    uniform_entry(0, wgpu::ShaderStages::VERTEX_FRAGMENT),
                    storage_entry(1, true, wgpu::ShaderStages::VERTEX),
                    storage_entry(2, true, wgpu::ShaderStages::VERTEX),
                ],
            }),
            cartoon_frame: device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
                label: Some("cartoon frame"),
                entries: &[
                    storage_entry(0, true, compute),
                    storage_entry(1, true, compute),
                    storage_entry(2, false, compute),
                ],
            }),
            cartoon_lod: device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
                label: Some("cartoon lod"),
                entries: &[
                    storage_entry(0, true, compute),
                    storage_entry(1, true, compute),
                    storage_entry(2, false, compute),
                    storage_entry(3, false, compute),
                ],
            }),
            // A plain vertex buffer (uniform + one storage array): no
            // level-of-detail or per-frame rebuild passes, unlike
            // `cartoon` above.
            glycan: device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
                label: Some("glycan"),
                entries: &[
                    uniform_entry(0, wgpu::ShaderStages::VERTEX_FRAGMENT),
                    storage_entry(1, true, wgpu::ShaderStages::VERTEX),
                ],
            }),
            gaussian_surface: device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
                label: Some("gaussian surface"),
                entries: &[
                    uniform_entry(0, wgpu::ShaderStages::FRAGMENT),
                    wgpu::BindGroupLayoutEntry {
                        binding: 1,
                        visibility: wgpu::ShaderStages::FRAGMENT,
                        ty: wgpu::BindingType::Texture {
                            sample_type: wgpu::TextureSampleType::Float { filterable: true },
                            view_dimension: wgpu::TextureViewDimension::D3,
                            multisampled: false,
                        },
                        count: None,
                    },
                    wgpu::BindGroupLayoutEntry {
                        binding: 2,
                        visibility: wgpu::ShaderStages::FRAGMENT,
                        ty: wgpu::BindingType::Sampler(wgpu::SamplerBindingType::Filtering),
                        count: None,
                    },
                    storage_entry(3, true, wgpu::ShaderStages::FRAGMENT),
                ],
            }),
            // Patches (binding 3) and params are read by the vertex stage
            // too: each patch's billboard is sized in `vs_skin_surface`.
            // Atoms (binding 0) also, for `vs_atom_cap`'s own billboard.
            skin_surface: device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
                label: Some("skin surface"),
                entries: &[
                    storage_entry(0, true, wgpu::ShaderStages::VERTEX_FRAGMENT),
                    storage_entry(1, true, wgpu::ShaderStages::FRAGMENT),
                    uniform_entry(2, wgpu::ShaderStages::VERTEX_FRAGMENT),
                    storage_entry(3, true, wgpu::ShaderStages::VERTEX_FRAGMENT),
                    storage_entry(4, true, wgpu::ShaderStages::FRAGMENT),
                    storage_entry(5, true, wgpu::ShaderStages::FRAGMENT),
                ],
            }),
            // Atoms, patches and probes size the billboards in
            // `vs_ses_surface`.
            ses_surface: device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
                label: Some("ses surface"),
                entries: &[
                    storage_entry(0, true, wgpu::ShaderStages::VERTEX_FRAGMENT),
                    storage_entry(1, true, wgpu::ShaderStages::FRAGMENT),
                    uniform_entry(2, wgpu::ShaderStages::VERTEX_FRAGMENT),
                    storage_entry(4, true, wgpu::ShaderStages::FRAGMENT),
                    storage_entry(5, true, wgpu::ShaderStages::VERTEX_FRAGMENT),
                    storage_entry(6, true, wgpu::ShaderStages::FRAGMENT),
                ],
            }),
            // One page of an SES's patch bounds, its survivors, first
            // patch and patches (`PageView` in shaders/ses_surface.wgsl).
            ses_page: device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
                label: Some("ses page"),
                entries: &[
                    storage_entry(0, true, wgpu::ShaderStages::VERTEX),
                    storage_entry(1, true, wgpu::ShaderStages::VERTEX),
                    uniform_entry(2, wgpu::ShaderStages::VERTEX_FRAGMENT),
                    storage_entry(3, true, wgpu::ShaderStages::VERTEX_FRAGMENT),
                ],
            }),
        };

        let common = include_str!("shaders/atoms.wgsl");
        let cull_module = shader(
            device,
            "cull",
            format!("{common}\n{}", include_str!("shaders/cull.wgsl")),
        );
        // Spheres and cylinders lie wholly behind their billboards.
        let shading = include_str!("shaders/shading.wgsl");
        let draw_source = with_early_depth(
            &ctx,
            format!("{common}\n{shading}\n{}", include_str!("shaders/draw.wgsl")),
            &[
                "fs_sphere",
                "fs_sphere_glass",
                "fs_sphere_pick",
                "fs_cylinder",
                "fs_cylinder_glass",
                "fs_cylinder_pick",
            ],
        );
        let draw_module = shader(device, "draw", draw_source);
        let hiz_copy_module = shader(
            device,
            "hiz copy",
            include_str!("shaders/hiz_copy.wgsl").to_string(),
        );
        let hiz_down_module = shader(
            device,
            "hiz down",
            include_str!("shaders/hiz_down.wgsl").to_string(),
        );
        let fullscreen = include_str!("shaders/fullscreen.wgsl");
        let blit_module = shader(
            device,
            "blit",
            format!("{fullscreen}\n{}", include_str!("shaders/blit.wgsl")),
        );
        let outline_module = shader(
            device,
            "outline",
            format!("{fullscreen}\n{}", include_str!("shaders/outline.wgsl")),
        );

        let cull_pipeline = compute_pipeline(
            device,
            "cull",
            &[&layouts.frame, &layouts.cull],
            &cull_module,
            "cull",
        );
        let cull_bonds_pipeline = compute_pipeline(
            device,
            "cull bonds",
            &[&layouts.frame, &layouts.cull],
            &cull_module,
            "cull_bonds",
        );
        let late_pipeline = |label: &str, entry: &str| {
            let layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
                label: Some(label),
                bind_group_layouts: &[Some(&layouts.frame), Some(&layouts.cull)],
                ..Default::default()
            });
            device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
                label: Some(label),
                layout: Some(&layout),
                module: &cull_module,
                entry_point: Some(entry),
                compilation_options: wgpu::PipelineCompilationOptions {
                    constants: &[("PHASE", 2.0)],
                    ..Default::default()
                },
                cache: None,
            })
        };
        let cull_late_pipeline = late_pipeline("cull late", "cull");
        let cull_clusters_late_pipeline = late_pipeline("cull clusters late", "cull_clusters");
        let cull_clusters_pipeline = compute_pipeline(
            device,
            "cull clusters",
            &[&layouts.frame, &layouts.cull],
            &cull_module,
            "cull_clusters",
        );
        let cull_bonds_late_pipeline = late_pipeline("cull bonds late", "cull_bonds");
        let hiz_copy_pipeline = compute_pipeline(
            device,
            "hiz copy",
            &[&layouts.hiz_copy],
            &hiz_copy_module,
            "copy_depth",
        );
        let hiz_down_pipeline = compute_pipeline(
            device,
            "hiz down",
            &[&layouts.hiz_down],
            &hiz_down_module,
            "downsample",
        );

        let draw_pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("draw"),
            bind_group_layouts: &[Some(&layouts.frame), Some(&layouts.draw)],
            ..Default::default()
        });
        let targets_desc = [
            Some(wgpu::ColorTargetState {
                format: COLOR_FORMAT,
                blend: None,
                write_mask: wgpu::ColorWrites::ALL,
            }),
            Some(wgpu::ColorTargetState {
                format: NORMAL_FORMAT,
                blend: None,
                write_mask: wgpu::ColorWrites::ALL,
            }),
            Some(wgpu::ColorTargetState {
                format: ID_FORMAT,
                blend: None,
                write_mask: wgpu::ColorWrites::ALL,
            }),
        ];
        let depth_state = wgpu::DepthStencilState {
            format: DEPTH_FORMAT,
            depth_write_enabled: Some(true),
            // Reversed-Z: larger is closer.
            depth_compare: Some(wgpu::CompareFunction::Greater),
            stencil: Default::default(),
            bias: Default::default(),
        };
        let make_draw_pipeline = |label: &str, vs: &str, fs: &str, topology| {
            device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
                label: Some(label),
                layout: Some(&draw_pipeline_layout),
                vertex: wgpu::VertexState {
                    module: &draw_module,
                    entry_point: Some(vs),
                    compilation_options: Default::default(),
                    buffers: &[],
                },
                primitive: wgpu::PrimitiveState {
                    topology,
                    cull_mode: None,
                    ..Default::default()
                },
                depth_stencil: Some(depth_state.clone()),
                multisample: Default::default(),
                fragment: Some(wgpu::FragmentState {
                    module: &draw_module,
                    entry_point: Some(fs),
                    compilation_options: Default::default(),
                    targets: &targets_desc,
                }),
                multiview_mask: None,
                cache: None,
            })
        };
        let sphere_pipeline = make_draw_pipeline(
            "spheres",
            "vs_sphere",
            "fs_sphere",
            wgpu::PrimitiveTopology::TriangleList,
        );
        let point_pipeline = make_draw_pipeline(
            "points",
            "vs_point",
            "fs_point",
            wgpu::PrimitiveTopology::PointList,
        );
        let cylinder_pipeline = make_draw_pipeline(
            "cylinders",
            "vs_cylinder",
            "fs_cylinder",
            wgpu::PrimitiveTopology::TriangleList,
        );
        let line_pipeline = make_draw_pipeline(
            "bond lines",
            "vs_line",
            "fs_line",
            wgpu::PrimitiveTopology::LineList,
        );

        // Cartoon ribbons: real triangles, built in the vertex shader from
        // stored cross-sections (shaders/cartoon.wgsl).
        let cartoon_lod_module = shader(
            device,
            "cartoon lod",
            include_str!("shaders/cartoon_lod.wgsl").to_string(),
        );
        let cartoon_lod_pipeline = compute_pipeline(
            device,
            "cartoon lod",
            &[&layouts.frame, &layouts.cartoon_lod],
            &cartoon_lod_module,
            "lod",
        );
        let cartoon_frame_module = shader(
            device,
            "cartoon frame",
            include_str!("shaders/cartoon_frame.wgsl").to_string(),
        );
        let cartoon_frame_pipeline = compute_pipeline(
            device,
            "cartoon frame",
            &[&layouts.cartoon_frame],
            &cartoon_frame_module,
            "cartoon_frame",
        );
        let cartoon_indices = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("cartoon join indices"),
            contents: bytemuck::cast_slice(&crate::scene::cartoon_join_indices()),
            usage: wgpu::BufferUsages::INDEX,
        });
        let cartoon_module = shader(
            device,
            "cartoon",
            format!("{shading}\n{}", include_str!("shaders/cartoon.wgsl")),
        );
        let cartoon_pipeline_layout =
            device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
                label: Some("cartoon"),
                bind_group_layouts: &[Some(&layouts.frame), Some(&layouts.cartoon)],
                ..Default::default()
            });
        let cartoon_pipeline = device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
            label: Some("cartoon"),
            layout: Some(&cartoon_pipeline_layout),
            vertex: wgpu::VertexState {
                module: &cartoon_module,
                entry_point: Some("vs_cartoon"),
                compilation_options: Default::default(),
                buffers: &[],
            },
            primitive: wgpu::PrimitiveState {
                topology: wgpu::PrimitiveTopology::TriangleList,
                // Back faces too: where a plane cuts a tube open they
                // draw its cap (fs_cartoon). Otherwise the front faces
                // hide them (the mesh is wound counter-clockwise from
                // outside, per `vv_core::cartoon`'s
                // `triangle_winding_agrees_with_outward_vertex_normals`).
                cull_mode: None,
                ..Default::default()
            },
            depth_stencil: Some(depth_state.clone()),
            multisample: Default::default(),
            fragment: Some(wgpu::FragmentState {
                module: &cartoon_module,
                entry_point: Some("fs_cartoon"),
                compilation_options: Default::default(),
                targets: &targets_desc,
            }),
            multiview_mask: None,
            cache: None,
        });

        // Glycan glyphs: a plain triangle soup, no compute passes (see
        // shaders/glycan_mesh.wgsl).
        let glycan_module = shader(
            device,
            "glycan mesh",
            format!("{shading}\n{}", include_str!("shaders/glycan_mesh.wgsl")),
        );
        let glycan_pipeline_layout =
            device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
                label: Some("glycan"),
                bind_group_layouts: &[Some(&layouts.frame), Some(&layouts.glycan)],
                ..Default::default()
            });
        let glycan_pipeline = device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
            label: Some("glycan"),
            layout: Some(&glycan_pipeline_layout),
            vertex: wgpu::VertexState {
                module: &glycan_module,
                entry_point: Some("vs_glycan"),
                compilation_options: Default::default(),
                buffers: &[],
            },
            primitive: wgpu::PrimitiveState {
                topology: wgpu::PrimitiveTopology::TriangleList,
                // No cap to draw where a plane cuts a glyph (fs_glycan
                // discards it whole instead), but double-sided anyway:
                // cheap insurance against a winding slip in any one of
                // the 9 SNFG shapes (`vv_core::glycan::mesh`'s own
                // `convex_shapes_have_outward_facing_normals` test).
                cull_mode: None,
                ..Default::default()
            },
            depth_stencil: Some(depth_state.clone()),
            multisample: Default::default(),
            fragment: Some(wgpu::FragmentState {
                module: &glycan_module,
                entry_point: Some("fs_glycan"),
                compilation_options: Default::default(),
                targets: &targets_desc,
            }),
            multiview_mask: None,
            cache: None,
        });

        // Gaussian surface: a full-screen triangle (no vertex buffer, like
        // the post-process passes below) but writing the same three opaque
        // targets + depth cartoons do (see shaders/gaussian_surface.wgsl),
        // so it composites and depth-tests against every other
        // representation in the same opaque pass rather than needing one
        // of its own.
        let gaussian_surface_module = shader(
            device,
            "gaussian surface",
            format!(
                "{shading}\n{}",
                include_str!("shaders/gaussian_surface.wgsl")
            ),
        );
        let gaussian_surface_pipeline_layout =
            device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
                label: Some("gaussian surface"),
                bind_group_layouts: &[Some(&layouts.frame), Some(&layouts.gaussian_surface)],
                ..Default::default()
            });
        let gaussian_surface_pipeline =
            device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
                label: Some("gaussian surface"),
                layout: Some(&gaussian_surface_pipeline_layout),
                vertex: wgpu::VertexState {
                    module: &gaussian_surface_module,
                    entry_point: Some("vs_gaussian_surface"),
                    compilation_options: Default::default(),
                    buffers: &[],
                },
                primitive: Default::default(),
                depth_stencil: Some(depth_state.clone()),
                multisample: Default::default(),
                fragment: Some(wgpu::FragmentState {
                    module: &gaussian_surface_module,
                    entry_point: Some("fs_gaussian_surface"),
                    compilation_options: Default::default(),
                    targets: &targets_desc,
                }),
                multiview_mask: None,
                cache: None,
            });

        // Skin surface and SES: per-patch billboards (six vertices each,
        // no vertex buffer -- the patch storage buffer is indexed by
        // vertex_index / 6, like the sphere impostors index their atoms)
        // into the same opaque pass, so patches depth-test against each
        // other and every other representation. A patch lies inside its
        // bounding sphere, behind its billboard.
        let billboard = include_str!("shaders/billboard.wgsl");
        let patch_module = |label: &str, stem: &str, source: &str| {
            shader(
                device,
                label,
                with_early_depth(
                    &ctx,
                    format!("{shading}\n{billboard}\n{source}"),
                    &[&format!("fs_{stem}"), &format!("fs_{stem}_owner")],
                ),
            )
        };
        let skin_surface_module = patch_module(
            "skin surface",
            "skin_surface",
            &format!(
                "{}
{}",
                include_str!("shaders/skin_patch.wgsl"),
                include_str!("shaders/skin_surface.wgsl")
            ),
        );
        let ses_surface_module = patch_module(
            "ses surface",
            "ses_surface",
            &format!(
                "{}
{}",
                include_str!("shaders/ses_patch.wgsl"),
                include_str!("shaders/ses_surface.wgsl")
            ),
        );

        // Glass: the same surfaces into the weighted-blended OIT targets,
        // depth-tested against the opaque pass but writing no depth,
        // normal or id, so AO, shadows, outlines and picking all see
        // through them.
        let glass_targets = [
            Some(wgpu::ColorTargetState {
                format: ACCUM_FORMAT,
                blend: Some(wgpu::BlendState {
                    color: wgpu::BlendComponent {
                        src_factor: wgpu::BlendFactor::One,
                        dst_factor: wgpu::BlendFactor::One,
                        operation: wgpu::BlendOperation::Add,
                    },
                    alpha: wgpu::BlendComponent {
                        src_factor: wgpu::BlendFactor::One,
                        dst_factor: wgpu::BlendFactor::One,
                        operation: wgpu::BlendOperation::Add,
                    },
                }),
                write_mask: wgpu::ColorWrites::ALL,
            }),
            Some(wgpu::ColorTargetState {
                format: REVEAL_FORMAT,
                blend: Some(wgpu::BlendState {
                    color: wgpu::BlendComponent {
                        src_factor: wgpu::BlendFactor::Zero,
                        dst_factor: wgpu::BlendFactor::OneMinusSrc,
                        operation: wgpu::BlendOperation::Add,
                    },
                    alpha: wgpu::BlendComponent::REPLACE,
                }),
                write_mask: wgpu::ColorWrites::ALL,
            }),
        ];
        let glass_depth = wgpu::DepthStencilState {
            depth_write_enabled: Some(false),
            ..depth_state.clone()
        };
        let make_glass_pipeline =
            |label: &str,
             layout: &wgpu::PipelineLayout,
             module: &wgpu::ShaderModule,
             vs: &str,
             fs: &str,
             depth_compare: wgpu::CompareFunction| {
                device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
                    label: Some(label),
                    layout: Some(layout),
                    vertex: wgpu::VertexState {
                        module,
                        entry_point: Some(vs),
                        compilation_options: Default::default(),
                        buffers: &[],
                    },
                    primitive: Default::default(),
                    depth_stencil: Some(wgpu::DepthStencilState {
                        depth_compare: Some(depth_compare),
                        ..glass_depth.clone()
                    }),
                    multisample: Default::default(),
                    fragment: Some(wgpu::FragmentState {
                        module,
                        entry_point: Some(fs),
                        compilation_options: Default::default(),
                        targets: &glass_targets,
                    }),
                    multiview_mask: None,
                    cache: None,
                })
            };
        let gaussian_glass_pipeline = make_glass_pipeline(
            "gaussian surface glass",
            &gaussian_surface_pipeline_layout,
            &gaussian_surface_module,
            "vs_gaussian_surface",
            "fs_gaussian_surface_glass",
            wgpu::CompareFunction::Greater,
        );

        // Every impostor/mesh representation's own glass (accum/reveal)
        // and pick (id-only) pipelines: transparent atoms, bonds,
        // cartoons and glycans are unrestricted multi-layer weighted OIT
        // (no owner buffer — each ray-cast or rasterized primitive is
        // already exactly one front-facing layer on its own), unlike the
        // patch surfaces below. Pick uses the *opaque* depth-write state
        // (`depth_state`: write enabled, compare greater) against the
        // dedicated `glass_depth` scratch buffer the "glass id" prepass
        // clears to the far plane (`renderer.rs`'s `render_all`), so the
        // nearest transparent fragment at a pixel wins exactly like an
        // opaque one would.
        let pick_targets = [Some(wgpu::ColorTargetState {
            format: ID_FORMAT,
            blend: None,
            write_mask: wgpu::ColorWrites::ALL,
        })];
        let make_pick_pipeline =
            |label: &str,
             layout: &wgpu::PipelineLayout,
             module: &wgpu::ShaderModule,
             vs: &str,
             fs: &str,
             topology: wgpu::PrimitiveTopology| {
                device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
                    label: Some(label),
                    layout: Some(layout),
                    vertex: wgpu::VertexState {
                        module,
                        entry_point: Some(vs),
                        compilation_options: Default::default(),
                        buffers: &[],
                    },
                    primitive: wgpu::PrimitiveState {
                        topology,
                        cull_mode: None,
                        ..Default::default()
                    },
                    depth_stencil: Some(depth_state.clone()),
                    multisample: Default::default(),
                    fragment: Some(wgpu::FragmentState {
                        module,
                        entry_point: Some(fs),
                        compilation_options: Default::default(),
                        targets: &pick_targets,
                    }),
                    multiview_mask: None,
                    cache: None,
                })
            };
        let sphere_glass_pipeline = make_glass_pipeline(
            "spheres glass",
            &draw_pipeline_layout,
            &draw_module,
            "vs_sphere",
            "fs_sphere_glass",
            wgpu::CompareFunction::Greater,
        );
        let sphere_pick_pipeline = make_pick_pipeline(
            "spheres pick",
            &draw_pipeline_layout,
            &draw_module,
            "vs_sphere",
            "fs_sphere_pick",
            wgpu::PrimitiveTopology::TriangleList,
        );
        let point_glass_pipeline = make_glass_pipeline(
            "points glass",
            &draw_pipeline_layout,
            &draw_module,
            "vs_point",
            "fs_point_glass",
            wgpu::CompareFunction::Greater,
        );
        let point_pick_pipeline = make_pick_pipeline(
            "points pick",
            &draw_pipeline_layout,
            &draw_module,
            "vs_point",
            "fs_point_pick",
            wgpu::PrimitiveTopology::PointList,
        );
        let cylinder_glass_pipeline = make_glass_pipeline(
            "cylinders glass",
            &draw_pipeline_layout,
            &draw_module,
            "vs_cylinder",
            "fs_cylinder_glass",
            wgpu::CompareFunction::Greater,
        );
        let cylinder_pick_pipeline = make_pick_pipeline(
            "cylinders pick",
            &draw_pipeline_layout,
            &draw_module,
            "vs_cylinder",
            "fs_cylinder_pick",
            wgpu::PrimitiveTopology::TriangleList,
        );
        let line_glass_pipeline = make_glass_pipeline(
            "bond lines glass",
            &draw_pipeline_layout,
            &draw_module,
            "vs_line",
            "fs_line_glass",
            wgpu::CompareFunction::Greater,
        );
        let line_pick_pipeline = make_pick_pipeline(
            "bond lines pick",
            &draw_pipeline_layout,
            &draw_module,
            "vs_line",
            "fs_line_pick",
            wgpu::PrimitiveTopology::LineList,
        );
        let cartoon_glass_pipeline = make_glass_pipeline(
            "cartoon glass",
            &cartoon_pipeline_layout,
            &cartoon_module,
            "vs_cartoon",
            "fs_cartoon_glass",
            wgpu::CompareFunction::Greater,
        );
        let cartoon_pick_pipeline = make_pick_pipeline(
            "cartoon pick",
            &cartoon_pipeline_layout,
            &cartoon_module,
            "vs_cartoon",
            "fs_cartoon_pick",
            wgpu::PrimitiveTopology::TriangleList,
        );
        let glycan_glass_pipeline = make_glass_pipeline(
            "glycan glass",
            &glycan_pipeline_layout,
            &glycan_module,
            "vs_glycan",
            "fs_glycan_glass",
            wgpu::CompareFunction::Greater,
        );
        let glycan_pick_pipeline = make_pick_pipeline(
            "glycan pick",
            &glycan_pipeline_layout,
            &glycan_module,
            "vs_glycan",
            "fs_glycan_pick",
            wgpu::PrimitiveTopology::TriangleList,
        );

        // Group 2 is the glass owner texture (glass only); group 3, where
        // given, what the vertex stage walks (the SES's culled pages).
        let patch_pipelines = |label: &str,
                               surface_layout: &wgpu::BindGroupLayout,
                               walk: Option<&wgpu::BindGroupLayout>,
                               module: &wgpu::ShaderModule,
                               stem: &str| {
            let layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
                label: Some(label),
                bind_group_layouts: &[Some(&layouts.frame), Some(surface_layout), None, walk],
                ..Default::default()
            });
            let glass_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
                label: Some(label),
                bind_group_layouts: &[
                    Some(&layouts.frame),
                    Some(surface_layout),
                    Some(&layouts.glass_owner),
                    walk,
                ],
                ..Default::default()
            });
            let vs = format!("vs_{stem}");
            let pipeline = |vs: &str, fs: &str, targets: &[Option<wgpu::ColorTargetState>]| {
                device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
                    label: Some(label),
                    layout: Some(&layout),
                    vertex: wgpu::VertexState {
                        module,
                        entry_point: Some(vs),
                        compilation_options: Default::default(),
                        buffers: &[],
                    },
                    primitive: Default::default(),
                    depth_stencil: Some(depth_state.clone()),
                    multisample: Default::default(),
                    fragment: Some(wgpu::FragmentState {
                        module,
                        entry_point: Some(fs),
                        compilation_options: Default::default(),
                        targets,
                    }),
                    multiview_mask: None,
                    cache: None,
                })
            };
            PatchPipelines {
                opaque: pipeline(&vs, &format!("fs_{stem}"), &targets_desc),
                owner: pipeline(
                    &vs,
                    &format!("fs_{stem}_owner"),
                    &[Some(wgpu::ColorTargetState {
                        format: OWNER_FORMAT,
                        blend: None,
                        write_mask: wgpu::ColorWrites::ALL,
                    })],
                ),
                glass: make_glass_pipeline(
                    label,
                    &glass_layout,
                    module,
                    &vs,
                    &format!("fs_{stem}_glass"),
                    wgpu::CompareFunction::Greater,
                ),
                atom_cap: pipeline("vs_atom_cap", "fs_atom_cap", &targets_desc),
            }
        };
        let skin = patch_pipelines(
            "skin surface",
            &layouts.skin_surface,
            None,
            &skin_surface_module,
            "skin_surface",
        );
        let ses = patch_pipelines(
            "ses surface",
            &layouts.ses_surface,
            Some(&layouts.ses_page),
            &ses_surface_module,
            "ses_surface",
        );

        // Post-process chain, all full-screen triangles over COLOR_FORMAT:
        // `color` -> outline pass -> `outlined` -> FXAA or copy -> `fxaa`
        // (the display). Every pass runs every frame (the outline pass is
        // uniform-gated, FXAA is a pipeline swap), so no texture in the
        // chain ever changes identity -- the display texture is what's
        // registered with egui, and switching it would need
        // re-registration exactly like a resize does.
        let make_post_pipeline =
            |label: &str, layout: &wgpu::BindGroupLayout, module: &wgpu::ShaderModule, fs: &str| {
                let pipeline_layout =
                    device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
                        label: Some(label),
                        bind_group_layouts: &[Some(layout)],
                        ..Default::default()
                    });
                device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
                    label: Some(label),
                    layout: Some(&pipeline_layout),
                    vertex: wgpu::VertexState {
                        module,
                        entry_point: Some("vs_fullscreen"),
                        compilation_options: Default::default(),
                        buffers: &[],
                    },
                    primitive: Default::default(),
                    depth_stencil: None,
                    multisample: Default::default(),
                    fragment: Some(wgpu::FragmentState {
                        module,
                        entry_point: Some(fs),
                        compilation_options: Default::default(),
                        targets: &[Some(wgpu::ColorTargetState {
                            format: COLOR_FORMAT,
                            blend: None,
                            write_mask: wgpu::ColorWrites::ALL,
                        })],
                    }),
                    multiview_mask: None,
                    cache: None,
                })
            };
        let outline_pipeline =
            make_post_pipeline("outline", &layouts.outline, &outline_module, "fs_outline");
        let dof_module = shader(
            device,
            "dof",
            format!("{fullscreen}\n{}", include_str!("shaders/dof.wgsl")),
        );
        let dof_pipeline = make_post_pipeline("dof", &layouts.dof, &dof_module, "fs_dof");
        let ao_module = shader(
            device,
            "ao",
            format!("{fullscreen}\n{}", include_str!("shaders/ao.wgsl")),
        );
        let ao_pipeline = {
            let layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
                label: Some("ao"),
                bind_group_layouts: &[Some(&layouts.ao), Some(&layouts.occlusion)],
                ..Default::default()
            });
            device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
                label: Some("ao"),
                layout: Some(&layout),
                vertex: wgpu::VertexState {
                    module: &ao_module,
                    entry_point: Some("vs_fullscreen"),
                    compilation_options: Default::default(),
                    buffers: &[],
                },
                primitive: Default::default(),
                depth_stencil: None,
                multisample: Default::default(),
                fragment: Some(wgpu::FragmentState {
                    module: &ao_module,
                    entry_point: Some("fs_ao"),
                    compilation_options: Default::default(),
                    targets: &[Some(wgpu::ColorTargetState {
                        format: AO_FORMAT,
                        blend: None,
                        write_mask: wgpu::ColorWrites::ALL,
                    })],
                }),
                multiview_mask: None,
                cache: None,
            })
        };
        let fxaa_pipeline = make_post_pipeline("fxaa", &layouts.blit, &blit_module, "fs_fxaa");
        let cap_depth_pipeline = cap_depth_pipeline(
            device,
            &layouts,
            &format!(
                "{fullscreen}\n{shading}\n{}",
                include_str!("shaders/cap_depth.wgsl")
            ),
        );
        let occlusion_sampler = device.create_sampler(&wgpu::SamplerDescriptor {
            label: Some("occlusion volume"),
            mag_filter: wgpu::FilterMode::Linear,
            min_filter: wgpu::FilterMode::Linear,
            mipmap_filter: wgpu::MipmapFilterMode::Linear,
            ..Default::default()
        });
        let empty_volume = crate::occlusion::OcclusionVolume::empty(&ctx);
        let occlusion_bind_group =
            occlusion_bind_group(device, &layouts, &empty_volume.view, &occlusion_sampler);
        let glass_resolve_module = shader(
            device,
            "glass resolve",
            format!(
                "{fullscreen}\n{}",
                include_str!("shaders/glass_resolve.wgsl")
            ),
        );
        let glass_resolve_pipeline = {
            let layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
                label: Some("glass resolve"),
                bind_group_layouts: &[Some(&layouts.glass_resolve)],
                ..Default::default()
            });
            device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
                label: Some("glass resolve"),
                layout: Some(&layout),
                vertex: wgpu::VertexState {
                    module: &glass_resolve_module,
                    entry_point: Some("vs_fullscreen"),
                    compilation_options: Default::default(),
                    buffers: &[],
                },
                primitive: Default::default(),
                depth_stencil: None,
                multisample: Default::default(),
                fragment: Some(wgpu::FragmentState {
                    module: &glass_resolve_module,
                    entry_point: Some("fs_glass_resolve"),
                    compilation_options: Default::default(),
                    targets: &[Some(wgpu::ColorTargetState {
                        format: COLOR_FORMAT,
                        blend: Some(wgpu::BlendState::ALPHA_BLENDING),
                        write_mask: wgpu::ColorWrites::ALL,
                    })],
                }),
                multiview_mask: None,
                cache: None,
            })
        };
        let passthrough_pipeline =
            make_post_pipeline("post passthrough", &layouts.blit, &blit_module, "fs_blit");
        let blit_sampler = device.create_sampler(&wgpu::SamplerDescriptor {
            label: Some("blit"),
            mag_filter: wgpu::FilterMode::Linear,
            min_filter: wgpu::FilterMode::Linear,
            ..Default::default()
        });

        let timer = GpuTimer::new(&ctx);
        let targets = Targets::new(
            device,
            width,
            height,
            &layouts,
            &camera_buffer,
            &lighting_buffer,
            &outline_buffer,
            &dof_buffer,
            &blit_sampler,
        );
        let selection = SelectionOutline::new(device, fullscreen, &targets.id_view);
        Self {
            ctx,
            layouts,
            targets,
            camera_buffer,
            lighting_buffer,
            outline_buffer,
            dof_buffer,
            dof_pipeline,
            cull_pipeline,
            cull_bonds_pipeline,
            cull_late_pipeline,
            cull_bonds_late_pipeline,
            cull_clusters_pipeline,
            cull_clusters_late_pipeline,
            sphere_pipeline,
            sphere_glass_pipeline,
            sphere_pick_pipeline,
            point_pipeline,
            point_glass_pipeline,
            point_pick_pipeline,
            cylinder_pipeline,
            cylinder_glass_pipeline,
            cylinder_pick_pipeline,
            line_pipeline,
            line_glass_pipeline,
            line_pick_pipeline,
            cartoon_pipeline,
            cartoon_glass_pipeline,
            cartoon_pick_pipeline,
            cartoon_lod_pipeline,
            cartoon_frame_pipeline,
            cartoon_indices,
            glycan_pipeline,
            glycan_glass_pipeline,
            glycan_pick_pipeline,
            gaussian_surface_pipeline,
            hiz_copy_pipeline,
            hiz_down_pipeline,
            outline_pipeline,
            ao_pipeline,
            fxaa_pipeline,
            cap_depth_pipeline,
            passthrough_pipeline,
            gaussian_glass_pipeline,
            skin,
            ses,
            glass_resolve_pipeline,
            blit_sampler,
            occlusion: None,
            occlusion_sampler,
            empty_volume,
            occlusion_bind_group,
            timer,
            last_times: None,
            pick_ranges: Vec::new(),
            selection,
        }
    }

    pub fn context(&self) -> &Arc<GpuContext> {
        &self.ctx
    }

    pub fn size(&self) -> (u32, u32) {
        (self.targets.width, self.targets.height)
    }

    pub fn resize(&mut self, width: u32, height: u32) {
        if width == 0 || height == 0 {
            return;
        }
        if (width, height) != self.size() {
            self.targets = Targets::new(
                &self.ctx.device,
                width,
                height,
                &self.layouts,
                &self.camera_buffer,
                &self.lighting_buffer,
                &self.outline_buffer,
                &self.dof_buffer,
                &self.blit_sampler,
            );
            self.selection
                .rebind_ids(&self.ctx.device, &self.targets.id_view);
        }
    }

    /// What the selection outline marks from the next `render_all` on:
    /// one [`ItemSelection`] per pick range, in pick order (items, then
    /// cartoons, then patch surfaces). The frame-wide bitset is rebuilt
    /// and uploaded only when this or the frame's id ranges change. An
    /// empty list, or no set bit, draws no outline.
    pub fn set_selection(&mut self, items: Vec<ItemSelection>) {
        self.selection.set(items);
    }

    /// The raw opaque-pass output, before the outline pass and FXAA.
    /// `read_color` (and so screenshots) reads the *outlined* image, one
    /// step later; nothing reads this directly except the post passes.
    pub fn color_view(&self) -> &wgpu::TextureView {
        &self.targets.color_view
    }

    pub fn color_texture(&self) -> &wgpu::Texture {
        &self.targets.color
    }

    /// What's actually shown to the user: `color_view` after the
    /// post-process pass (FXAA, or a plain copy when `RenderSettings::fxaa`
    /// is off). Always the same texture regardless of that setting — vv-app
    /// registers it with egui once and only re-registers on resize
    /// (`ensure_viewport_texture`/`maybe_resize_viewport`), so this can't
    /// switch between two different textures without that registration
    /// going stale exactly the way an un-handled resize would.
    pub fn display_view(&self) -> &wgpu::TextureView {
        &self.targets.fxaa_view
    }

    /// `display_view` as egui must sample it (`fxaa_gamma_view`).
    pub fn display_view_for_ui(&self) -> &wgpu::TextureView {
        &self.targets.fxaa_gamma_view
    }

    /// Clears both `color_view` and `display_view` to `background`, for
    /// when there's nothing to render (e.g. no structure loaded). Clears
    /// both targets, not just one, so a screenshot taken in this state
    /// (`color_view`) and what's on screen (`display_view`) agree.
    pub fn clear(&self, encoder: &mut wgpu::CommandEncoder, background: wgpu::Color) {
        fn pass<'a>(
            encoder: &'a mut wgpu::CommandEncoder,
            view: &'a wgpu::TextureView,
            clear: wgpu::Color,
        ) {
            let _pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("clear"),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view,
                    resolve_target: None,
                    ops: wgpu::Operations {
                        load: wgpu::LoadOp::Clear(clear),
                        store: wgpu::StoreOp::Store,
                    },
                    depth_slice: None,
                })],
                ..Default::default()
            });
        }
        pass(encoder, &self.targets.color_view, background);
        pass(encoder, &self.targets.outlined_view, background);
        pass(encoder, &self.targets.fxaa_view, background);
    }

    /// GPU pass times of the most recent frame whose results have arrived.
    pub fn last_times(&self) -> Option<FrameTimes> {
        self.last_times
    }

    /// Bind groups for `structure` drawn with its own draw state.
    pub fn bind(&self, structure: &GpuStructure) -> PageBindings {
        self.bind_state(structure, &structure.state)
    }

    /// Bind groups for `structure`'s geometry drawn with `state` (another
    /// rep of the same atoms: its own colours, parameters and survivors).
    pub fn bind_state(&self, structure: &GpuStructure, state: &DrawState) -> PageBindings {
        let device = &self.ctx.device;
        fn entry(binding: u32, buffer: &wgpu::Buffer) -> wgpu::BindGroupEntry<'_> {
            wgpu::BindGroupEntry {
                binding,
                resource: buffer.as_entire_binding(),
            }
        }
        let pages = structure
            .pages
            .iter()
            .zip(&state.pages)
            .map(|(p, s): (&Page, &PageState)| PageBinding {
                cull: device.create_bind_group(&wgpu::BindGroupDescriptor {
                    label: Some("cull page"),
                    layout: &self.layouts.cull,
                    entries: &[
                        entry(0, &p.atoms),
                        entry(1, &s.colors),
                        entry(2, &s.visible_quads),
                        entry(3, &s.visible_points),
                        entry(4, &s.indirect),
                        entry(5, &p.bonds),
                        entry(6, &s.visible_bonds),
                        entry(7, &s.params),
                        entry(8, &s.occluded),
                        entry(9, &s.clusters),
                        entry(10, &p.cluster_bounds),
                    ],
                }),
                draw: device.create_bind_group(&wgpu::BindGroupDescriptor {
                    label: Some("draw page"),
                    layout: &self.layouts.draw,
                    entries: &[
                        entry(0, &p.atoms),
                        entry(1, &s.colors),
                        entry(2, &s.visible_quads),
                        entry(3, &s.visible_points),
                        entry(4, &p.bonds),
                        entry(5, &s.visible_bonds),
                        entry(6, &p.bond_ids),
                        entry(7, &s.params),
                    ],
                }),
                params: s.params.clone(),
                indirect: s.indirect.clone(),
                occluded: s.occluded.clone(),
            })
            .collect();
        PageBindings { pages }
    }

    /// The bind groups a [`CartoonGpu`] needs — created once per mesh,
    /// like `bind`.
    pub fn bind_cartoon(&self, mesh: &CartoonGpu) -> CartoonBindings {
        let device = &self.ctx.device;
        fn entry(binding: u32, buffer: &wgpu::Buffer) -> wgpu::BindGroupEntry<'_> {
            wgpu::BindGroupEntry {
                binding,
                resource: buffer.as_entire_binding(),
            }
        }
        CartoonBindings::Ribbon {
            lod: device.create_bind_group(&wgpu::BindGroupDescriptor {
                label: Some("cartoon lod"),
                layout: &self.layouts.cartoon_lod,
                entries: &[
                    entry(0, &mesh.sections),
                    entry(1, &mesh.spans),
                    entry(2, &mesh.pairs),
                    entry(3, &mesh.indirect),
                ],
            }),
            draw: device.create_bind_group(&wgpu::BindGroupDescriptor {
                label: Some("cartoon"),
                layout: &self.layouts.cartoon,
                entries: &[
                    entry(0, &mesh.params),
                    entry(1, &mesh.sections),
                    entry(2, &mesh.pairs),
                ],
            }),
            frame: device.create_bind_group(&wgpu::BindGroupDescriptor {
                label: Some("cartoon frame"),
                layout: &self.layouts.cartoon_frame,
                entries: &[
                    entry(0, &mesh.recipes),
                    entry(1, &mesh.spline),
                    entry(2, &mesh.sections),
                ],
            }),
        }
    }

    /// The bind group a [`GlycanGpu`] needs to draw — created once per
    /// upload, like `bind_cartoon`; unlike a ribbon, there is no
    /// level-of-detail or per-frame rebuild pass to bind.
    pub fn bind_glycan(&self, mesh: &GlycanGpu) -> CartoonBindings {
        CartoonBindings::Glycan {
            draw: self
                .ctx
                .device
                .create_bind_group(&wgpu::BindGroupDescriptor {
                    label: Some("glycan"),
                    layout: &self.layouts.glycan,
                    entries: &[
                        wgpu::BindGroupEntry {
                            binding: 0,
                            resource: mesh.params.as_entire_binding(),
                        },
                        wgpu::BindGroupEntry {
                            binding: 1,
                            resource: mesh.vertices.as_entire_binding(),
                        },
                    ],
                }),
        }
    }

    pub fn bind_gaussian_surface(&self, gpu: &GaussianSurfaceGpu) -> wgpu::BindGroup {
        self.ctx
            .device
            .create_bind_group(&wgpu::BindGroupDescriptor {
                label: Some("gaussian surface"),
                layout: &self.layouts.gaussian_surface,
                entries: &[
                    wgpu::BindGroupEntry {
                        binding: 0,
                        resource: gpu.params.as_entire_binding(),
                    },
                    wgpu::BindGroupEntry {
                        binding: 1,
                        resource: wgpu::BindingResource::TextureView(&gpu.volume_view),
                    },
                    wgpu::BindGroupEntry {
                        binding: 2,
                        resource: wgpu::BindingResource::Sampler(&gpu.sampler),
                    },
                    wgpu::BindGroupEntry {
                        binding: 3,
                        resource: gpu.bricks.as_entire_binding(),
                    },
                ],
            })
    }

    pub fn bind_skin_surface(&self, gpu: &SkinSurfaceGpu) -> wgpu::BindGroup {
        self.ctx
            .device
            .create_bind_group(&wgpu::BindGroupDescriptor {
                label: Some("skin surface"),
                layout: &self.layouts.skin_surface,
                entries: &[
                    wgpu::BindGroupEntry {
                        binding: 0,
                        resource: gpu.atoms.as_entire_binding(),
                    },
                    wgpu::BindGroupEntry {
                        binding: 1,
                        resource: gpu.colors.as_entire_binding(),
                    },
                    wgpu::BindGroupEntry {
                        binding: 2,
                        resource: gpu.params.as_entire_binding(),
                    },
                    wgpu::BindGroupEntry {
                        binding: 3,
                        resource: gpu.patches.as_entire_binding(),
                    },
                    wgpu::BindGroupEntry {
                        binding: 4,
                        resource: gpu.competitor_starts.as_entire_binding(),
                    },
                    wgpu::BindGroupEntry {
                        binding: 5,
                        resource: gpu.competitors.as_entire_binding(),
                    },
                ],
            })
    }

    pub fn bind_ses_surface(&self, gpu: &SesGpu) -> SesBindings {
        let device = &self.ctx.device;
        let buffers = [
            (0, &gpu.atoms),
            (1, &gpu.colors),
            (2, &gpu.params),
            (4, &gpu.caps),
            (5, &gpu.probes),
            (6, &gpu.probe_neighbors),
        ];
        let entries: Vec<wgpu::BindGroupEntry> = buffers
            .iter()
            .map(|&(binding, b)| wgpu::BindGroupEntry {
                binding,
                resource: b.as_entire_binding(),
            })
            .collect();
        let surface = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("ses surface"),
            layout: &self.layouts.ses_surface,
            entries: &entries,
        });
        let views = gpu
            .bounds
            .pages
            .iter()
            .zip(&gpu.bounds.state.pages)
            .zip(&gpu.page_bases)
            .zip(&gpu.patches)
            .map(|(((page, state), base), patches)| {
                device.create_bind_group(&wgpu::BindGroupDescriptor {
                    label: Some("ses page"),
                    layout: &self.layouts.ses_page,
                    entries: &[
                        wgpu::BindGroupEntry {
                            binding: 0,
                            resource: page.atoms.as_entire_binding(),
                        },
                        wgpu::BindGroupEntry {
                            binding: 1,
                            resource: state.visible_quads.as_entire_binding(),
                        },
                        wgpu::BindGroupEntry {
                            binding: 2,
                            resource: base.as_entire_binding(),
                        },
                        wgpu::BindGroupEntry {
                            binding: 3,
                            resource: patches.as_entire_binding(),
                        },
                    ],
                })
            })
            .collect();
        SesBindings {
            surface,
            pages: self.bind(&gpu.bounds),
            views,
        }
    }

    fn patch_pipelines(&self, surface: PatchSurface<'_>) -> &PatchPipelines {
        match surface {
            PatchSurface::Skin(..) => &self.skin,
            PatchSurface::Ses(..) => &self.ses,
        }
    }

    /// The scene's occlusion volume for long-range AO and shadows, or
    /// `None` for screen-space only. Kept until replaced; rebuild it when
    /// what is drawn changes, not when the camera moves.
    pub fn set_occlusion_volume(&mut self, volume: Option<Arc<crate::occlusion::OcclusionVolume>>) {
        let view = volume.as_ref().map_or(&self.empty_volume.view, |v| &v.view);
        self.occlusion_bind_group = occlusion_bind_group(
            &self.ctx.device,
            &self.layouts,
            view,
            &self.occlusion_sampler,
        );
        self.occlusion = volume;
    }

    /// Records one frame showing a single structure, drawn with
    /// `settings.representation`. A thin wrapper over [`Renderer::render_all`]
    /// for benchmarks, tests, and headless callers.
    pub fn render(
        &mut self,
        encoder: &mut wgpu::CommandEncoder,
        camera: &Camera,
        structure: &GpuStructure,
        bindings: &PageBindings,
        settings: &RenderSettings,
    ) {
        let item = DrawItem {
            structure,
            bindings,
            representation: settings.representation,
            sizes: AtomSizes::of(settings.representation),
            material: settings.material,
        };
        self.render_all(
            encoder,
            camera,
            std::slice::from_ref(&item),
            &[],
            &[],
            &[],
            settings,
        );
    }

    /// Records one frame compositing every item in `items` and every mesh
    /// in `cartoons` into the same targets: one clear, one cull pass and
    /// one opaque pass covering all of them (so they depth-test against
    /// each other), then post-process and the depth pyramid. Each item
    /// draws with its own representation; `settings.representation` is
    /// ignored here. Pick ids are assigned in `items` order, then
    /// `cartoons` order continuing the same id space (one vertex = one
    /// pick id, resolved back as `Pick::Atom { item, atom: vertex_index }`
    /// — the caller maps `vertex_index` to a real atom the same way it
    /// already does for tube samples, via `CartoonMesh::source`), and
    /// resolved back by [`Renderer::pick`].
    pub fn render_all(
        &mut self,
        encoder: &mut wgpu::CommandEncoder,
        camera: &Camera,
        items: &[DrawItem<'_>],
        cartoons: &[CartoonItem<'_>],
        gaussian_surfaces: &[GaussianSurfaceItem<'_>],
        patch_surfaces: &[PatchSurfaceItem<'_>],
        settings: &RenderSettings,
    ) {
        let viewport = Vec2::new(self.targets.width as f32, self.targets.height as f32);
        let mut uniform = camera.uniform(viewport);
        uniform.quad_px_threshold = settings.quad_px_threshold;
        uniform.hiz_mip_count = self.targets.hiz_mip_count;
        uniform.occlusion = settings.occlusion_culling as u32;
        if let Some(plane) = settings.clip {
            uniform.clip = camera.view_clip(plane);
        }
        // Skin and SES draw their solid cross-section as a separate pass
        // over every substrate atom (see `PatchPipelines::atom_cap`), so it
        // only runs while something can actually be cutting the surface:
        // an explicit clip plane, or a dolly forward. Perspective's near
        // cut is otherwise always technically live (just past the near
        // plane) but never visually reached without one of these.
        let cut_active = settings.clip.is_some() || camera.dolly > 0.0;
        self.ctx
            .queue
            .write_buffer(&self.camera_buffer, 0, bytemuck::bytes_of(&uniform));
        self.ctx.queue.write_buffer(
            &self.lighting_buffer,
            0,
            bytemuck::bytes_of(&settings.lighting),
        );
        // Depth cue range: the scene's own depth span from this eye.
        let bounds = settings.scene_bounds.or_else(|| {
            let spheres = items
                .iter()
                .map(|i| (i.structure.center, i.structure.radius))
                .chain(cartoons.iter().map(|c| c.mesh.bounds()))
                .chain(
                    gaussian_surfaces
                        .iter()
                        .map(|g| (g.gpu.bounds_center, g.gpu.bounds_radius)),
                )
                .chain(patch_surfaces.iter().map(|s| s.surface.bounds()));
            spheres.reduce(|(c1, r1), (c2, r2)| {
                let d = (c2 - c1).length();
                if d + r2 <= r1 {
                    (c1, r1)
                } else if d + r1 <= r2 {
                    (c2, r2)
                } else {
                    let r = (d + r1 + r2) * 0.5;
                    (c1 + (c2 - c1) * ((r - r1) / d.max(1e-6)), r)
                }
            })
        });
        let (fog_near, fog_far, dof_range) = match bounds {
            Some((center, radius)) => {
                let view_dir = (camera.target - camera.eye()).normalize_or_zero();
                let d = (center - camera.eye()).dot(view_dir);
                let (near, far) = fog_range(d, radius);
                (near, far, 0.8 * radius)
            }
            None => (0.0, 0.0, 1.0),
        };
        let proj = camera.proj(viewport.x / viewport.y);
        let bg = settings.background;
        let post = PostUniform {
            inv_proj: proj.inverse().to_cols_array_2d(),
            outline_enabled: settings.outline as u32,
            outline_width: settings.outline_width.max(1),
            near: camera.near,
            depth_threshold: OUTLINE_DEPTH_THRESHOLD,
            normal_threshold: OUTLINE_NORMAL_THRESHOLD,
            outline_strength: OUTLINE_STRENGTH,
            ao_strength: settings.ao.max(0.0),
            ao_radius: AO_RADIUS,
            fog_strength: settings.depth_cue.clamp(0.0, 1.0),
            fog_near,
            fog_far,
            proj_y: proj.y_axis.y,
            background: [bg.r as f32, bg.g as f32, bg.b as f32, bg.a as f32],
            orthographic: (camera.projection == crate::camera::Projection::Orthographic) as u32,
            shadow_strength: settings.shadows.clamp(0.0, 1.0),
            _pad: [0; 2],
            // Screen-space shadows/AO stay single-light: the first live
            // light stands in for "the key light" here.
            light_dir: {
                let [x, y, z] = settings.lighting.lights[0].dir;
                [x, y, z, 0.0]
            },
            view_inv: camera.view().inverse().to_cols_array_2d(),
            vol_origin: match &self.occlusion {
                Some(v) => v.origin.extend(v.voxel).to_array(),
                None => [0.0, 0.0, 0.0, 1.0],
            },
            vol_extent: match &self.occlusion {
                Some(v) => (glam::UVec3::from_array(v.dims).as_vec3() * v.voxel)
                    .extend(1.0)
                    .to_array(),
                None => [1.0, 1.0, 1.0, 0.0],
            },
            clip_world: settings.clip.unwrap_or(KEEP_ALL),
            cut_world: camera.near_cut_plane().unwrap_or(KEEP_ALL),
            background_top: {
                let top = settings.background_top.unwrap_or(bg);
                [top.r as f32, top.g as f32, top.b as f32, top.a as f32]
            },
        };
        self.ctx
            .queue
            .write_buffer(&self.outline_buffer, 0, bytemuck::bytes_of(&post));

        // Per-page parameters and the frame-wide id ranges behind `pick`.
        // Atom ids are page-relative in the shaders, so every page gets its
        // own base; bond ids are already structure-wide (`Page::bond_ids`),
        // so all pages of one item share that item's bond base.
        self.pick_ranges.clear();
        let (mut atom_base, mut bond_base) = (0u32, 0u32);
        for item in items {
            let glass_cull_margin = if settings.fast_glass {
                crate::scene::glass_cull_margin(
                    &item.material,
                    item.structure.atom_count,
                    item.structure.radius,
                )
            } else {
                0.0
            };
            for (p, page) in item.bindings.pages.iter().enumerate() {
                let params = PageParams {
                    radius_scale: item.sizes.radius_scale,
                    bond_radius: item.sizes.bond_radius,
                    radius_offset: item.sizes.radius_offset,
                    lines: (item.representation == Representation::Lines) as u32,
                    min_quad_radius: 0.0,
                    glass_cull_margin,
                    atom_id_base: atom_base + (p * PAGE_ATOMS) as u32,
                    bond_id_base: bond_base,
                    material: item.material,
                };
                page.reset(&self.ctx.queue, encoder, &params);
            }
            let range = PickRange {
                atom_base,
                atom_count: item.structure.atom_count as u32,
                bond_base,
                bond_count: item.structure.bond_table_len as u32,
            };
            self.pick_ranges.push(range);
            atom_base += range.atom_count;
            bond_base += range.bond_count;
        }
        // Cartoons continue the same frame-wide atom id space, one
        // `PickRange` per mesh appended right after `items`' own ranges
        // so their `item` index picks up where `items.len()` left off —
        // exactly the order `vv-app`'s `draw_items` appends their
        // `DrawSource`s in. No bond ids: a cartoon never draws bonds.
        for item in cartoons {
            let params = CartoonParams {
                id_base: atom_base,
                _pad: [0; 3],
                material: item.material,
            };
            self.ctx
                .queue
                .write_buffer(item.mesh.params(), 0, bytemuck::bytes_of(&params));
            let count = item.mesh.pick_count();
            self.pick_ranges.push(PickRange {
                atom_base,
                atom_count: count,
                bond_base,
                bond_count: 0,
            });
            atom_base += count;
        }
        // Gaussian surfaces are not part of the frame-wide pick id space
        // (not pickable yet — see shaders/gaussian_surface.wgsl's phase-1
        // cuts), so this only needs the one field that actually changes
        // frame to frame: the inverse view, letting the shader turn one
        // per-pixel ray into the world space its atoms live in without
        // transforming every atom itself every frame.
        let view_inv = camera.view().inverse();
        for item in gaussian_surfaces {
            let params = GaussianSurfaceParams {
                view_inv: view_inv.to_cols_array_2d(),
                material: item.material,
                ..item.gpu.base_params
            };
            self.ctx
                .queue
                .write_buffer(&item.gpu.params, 0, bytemuck::bytes_of(&params));
        }
        // Patch surfaces *are* pickable -- a hit resolves to an atom of
        // its patch -- so each continues the frame-wide atom id space
        // after the cartoons, one `PickRange` per surface covering its
        // structure's atoms (`vv-app`'s `draw_items` appends their
        // `DrawSource`s in the same order). `view_inv` for the same
        // reason as the Gaussian surfaces.
        for item in patch_surfaces {
            item.surface
                .write_params(&self.ctx.queue, view_inv, atom_base, item.material);
            if let PatchSurface::Ses(g, b) = item.surface {
                // Patch bounds cull as plain spheres.
                let params = PageParams {
                    radius_scale: 1.0,
                    bond_radius: 0.0,
                    radius_offset: 0.0,
                    lines: 0,
                    min_quad_radius: g.probe_radius,
                    // SES glass keeps its own owner-buffer prepass
                    // (docs/RENDERING.md's "Transparency"); untouched by
                    // glass self-occlusion culling.
                    glass_cull_margin: 0.0,
                    atom_id_base: NO_PICK_ID_BASE,
                    bond_id_base: 0,
                    material: item.material,
                };
                for page in &b.pages.pages {
                    page.reset(&self.ctx.queue, encoder, &params);
                }
            }
            self.pick_ranges.push(PickRange {
                atom_base,
                atom_count: item.surface.atom_count(),
                bond_base,
                bond_count: 0,
            });
            atom_base += item.surface.atom_count();
        }
        self.selection.prepare(
            &self.ctx.device,
            &self.ctx.queue,
            &self.pick_ranges,
            settings.selection_color,
            settings.selection_width,
        );
        let draws_bonds =
            |item: &DrawItem<'_>| item.sizes.bond_radius > 0.0 && item.structure.bond_count > 0;

        fn compute_stamps(
            timer: &Option<GpuTimer>,
            begin: u32,
        ) -> Option<wgpu::ComputePassTimestampWrites<'_>> {
            timer.as_ref().map(|t| wgpu::ComputePassTimestampWrites {
                query_set: t.query_set(),
                beginning_of_pass_write_index: Some(begin),
                end_of_pass_write_index: Some(begin + 1),
            })
        }

        // Every SES, culled and drawn page by page like atoms.
        let ses: Vec<(&SesGpu, &SesBindings, Material)> = patch_surfaces
            .iter()
            .filter_map(|i| match i.surface {
                PatchSurface::Ses(g, b) if g.patch_count > 0 => Some((g, b, i.material)),
                _ => None,
            })
            .collect();
        let cull = |encoder: &mut wgpu::CommandEncoder,
                    clusters: &wgpu::ComputePipeline,
                    atoms: &wgpu::ComputePipeline,
                    bonds: &wgpu::ComputePipeline,
                    stamps: Option<wgpu::ComputePassTimestampWrites<'_>>| {
            let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: Some("cull"),
                timestamp_writes: stamps,
            });
            pass.set_bind_group(0, &self.targets.frame_bind_group, &[]);
            for item in items {
                pass.set_pipeline(clusters);
                for (page, binding) in item.structure.pages.iter().zip(&item.bindings.pages) {
                    pass.set_bind_group(1, &binding.cull, &[]);
                    let groups = page
                        .atom_count
                        .div_ceil(CLUSTER_ATOMS)
                        .div_ceil(CULL_WORKGROUP);
                    pass.dispatch_workgroups(groups, 1, 1);
                }
                pass.set_pipeline(atoms);
                for (page, binding) in item.structure.pages.iter().zip(&item.bindings.pages) {
                    pass.set_bind_group(1, &binding.cull, &[]);
                    pass.dispatch_workgroups(page.atom_count.div_ceil(CULL_WORKGROUP), 1, 1);
                }
                if draws_bonds(item) {
                    pass.set_pipeline(bonds);
                    for (page, binding) in item.structure.pages.iter().zip(&item.bindings.pages) {
                        if page.bond_count > 0 {
                            pass.set_bind_group(1, &binding.cull, &[]);
                            pass.dispatch_workgroups(
                                page.bond_count.div_ceil(CULL_WORKGROUP),
                                1,
                                1,
                            );
                        }
                    }
                }
            }
            // SES patch bounds: clusters, then patches; no bonds.
            for (gpu, b, _) in &ses {
                let pages = || gpu.bounds.pages.iter().zip(&b.pages.pages);
                pass.set_pipeline(clusters);
                for (page, binding) in pages() {
                    pass.set_bind_group(1, &binding.cull, &[]);
                    let groups = page
                        .atom_count
                        .div_ceil(CLUSTER_ATOMS)
                        .div_ceil(CULL_WORKGROUP);
                    pass.dispatch_workgroups(groups, 1, 1);
                }
                pass.set_pipeline(atoms);
                for (page, binding) in pages() {
                    pass.set_bind_group(1, &binding.cull, &[]);
                    pass.dispatch_workgroups(page.atom_count.div_ceil(CULL_WORKGROUP), 1, 1);
                }
            }
        };
        let hiz = |encoder: &mut wgpu::CommandEncoder,
                   stamps: Option<wgpu::ComputePassTimestampWrites<'_>>| {
            let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: Some("depth pyramid"),
                timestamp_writes: stamps,
            });
            let (mut w, mut h) = (self.targets.width, self.targets.height);
            for (level, bind_group) in self.targets.hiz_bind_groups.iter().enumerate() {
                if level == 0 {
                    pass.set_pipeline(&self.hiz_copy_pipeline);
                } else {
                    pass.set_pipeline(&self.hiz_down_pipeline);
                    w = (w / 2).max(1);
                    h = (h / 2).max(1);
                }
                pass.set_bind_group(0, bind_group, &[]);
                pass.dispatch_workgroups(w.div_ceil(HIZ_WORKGROUP), h.div_ceil(HIZ_WORKGROUP), 1);
            }
        };
        // Draws every opaque item's culled atoms and bonds into an open
        // pass; a glass item's own atoms and bonds draw only in the glass
        // pass instead (`draw_atoms_glass` below).
        let draw_atoms = |pass: &mut wgpu::RenderPass<'_>, late: bool| {
            let base = if late { LATE_INDIRECT_OFFSET } else { 0 };
            for item in items {
                if item.material.is_glass() {
                    continue;
                }
                pass.set_pipeline(&self.sphere_pipeline);
                for binding in &item.bindings.pages {
                    pass.set_bind_group(1, &binding.draw, &[]);
                    pass.draw_indirect(&binding.indirect, base);
                }
                pass.set_pipeline(&self.point_pipeline);
                for binding in &item.bindings.pages {
                    pass.set_bind_group(1, &binding.draw, &[]);
                    pass.draw_indirect(&binding.indirect, base + POINTS_INDIRECT_OFFSET);
                }
                if draws_bonds(item) {
                    pass.set_pipeline(&self.cylinder_pipeline);
                    for (page, binding) in item.structure.pages.iter().zip(&item.bindings.pages) {
                        if page.bond_count > 0 {
                            pass.set_bind_group(1, &binding.draw, &[]);
                            pass.draw_indirect(&binding.indirect, base + BONDS_INDIRECT_OFFSET);
                        }
                    }
                    pass.set_pipeline(&self.line_pipeline);
                    for (page, binding) in item.structure.pages.iter().zip(&item.bindings.pages) {
                        if page.bond_count > 0 {
                            pass.set_bind_group(1, &binding.draw, &[]);
                            pass.draw_indirect(&binding.indirect, base + LINES_INDIRECT_OFFSET);
                        }
                    }
                }
            }
        };
        // A glass item's culled atoms and bonds, into the accum/reveal
        // targets (`sphere_glass_pipeline` etc.) or, with `pick: true`,
        // into `id_view` alone (`sphere_pick_pipeline` etc. — see the
        // "glass id" pass in `render_all`). `late` draws both phases'
        // survivors: by the time the colour glass pass runs, occlusion
        // phase 2 has already completed, unlike the id-only pass, which
        // runs before phase 1's own draw and so only ever sees `late =
        // false` (documented gap: a glass atom occlusion-culling only
        // recovers in phase 2 keeps last frame's pick for one frame).
        let draw_atoms_glass = |pass: &mut wgpu::RenderPass<'_>, pick: bool, late: bool| {
            let base = if late { LATE_INDIRECT_OFFSET } else { 0 };
            for item in items {
                if !item.material.is_glass() {
                    continue;
                }
                pass.set_pipeline(if pick {
                    &self.sphere_pick_pipeline
                } else {
                    &self.sphere_glass_pipeline
                });
                for binding in &item.bindings.pages {
                    pass.set_bind_group(1, &binding.draw, &[]);
                    pass.draw_indirect(&binding.indirect, base);
                }
                pass.set_pipeline(if pick {
                    &self.point_pick_pipeline
                } else {
                    &self.point_glass_pipeline
                });
                for binding in &item.bindings.pages {
                    pass.set_bind_group(1, &binding.draw, &[]);
                    pass.draw_indirect(&binding.indirect, base + POINTS_INDIRECT_OFFSET);
                }
                if draws_bonds(item) {
                    pass.set_pipeline(if pick {
                        &self.cylinder_pick_pipeline
                    } else {
                        &self.cylinder_glass_pipeline
                    });
                    for (page, binding) in item.structure.pages.iter().zip(&item.bindings.pages) {
                        if page.bond_count > 0 {
                            pass.set_bind_group(1, &binding.draw, &[]);
                            pass.draw_indirect(&binding.indirect, base + BONDS_INDIRECT_OFFSET);
                        }
                    }
                    pass.set_pipeline(if pick {
                        &self.line_pick_pipeline
                    } else {
                        &self.line_glass_pipeline
                    });
                    for (page, binding) in item.structure.pages.iter().zip(&item.bindings.pages) {
                        if page.bond_count > 0 {
                            pass.set_bind_group(1, &binding.draw, &[]);
                            pass.draw_indirect(&binding.indirect, base + LINES_INDIRECT_OFFSET);
                        }
                    }
                }
            }
        };
        // A glass mesh's own triangles, colour or (`pick: true`) id only.
        let draw_cartoons_glass = |pass: &mut wgpu::RenderPass<'_>, pick: bool| {
            for item in cartoons {
                if !item.material.is_glass() {
                    continue;
                }
                pass.set_bind_group(1, item.bindings.draw(), &[]);
                match item.mesh {
                    CartoonMesh::Ribbon(g) => {
                        pass.set_pipeline(if pick {
                            &self.cartoon_pick_pipeline
                        } else {
                            &self.cartoon_glass_pipeline
                        });
                        pass.set_index_buffer(
                            self.cartoon_indices.slice(..),
                            wgpu::IndexFormat::Uint16,
                        );
                        pass.draw_indexed_indirect(&g.indirect, 0);
                    }
                    CartoonMesh::Glycan(g) => {
                        pass.set_pipeline(if pick {
                            &self.glycan_pick_pipeline
                        } else {
                            &self.glycan_glass_pipeline
                        });
                        pass.draw(0..g.vertex_count, 0..1);
                    }
                }
            }
        };
        let any_glass_items = items.iter().any(|i| i.material.is_glass());
        let any_glass_cartoons = cartoons.iter().any(|c| c.material.is_glass());
        let any_pick_glass = any_glass_items || any_glass_cartoons;
        // Draws every opaque SES's culled patches, and those below a pixel
        // as points in the colour of their atom.
        let draw_ses = |pass: &mut wgpu::RenderPass<'_>, late: bool| {
            let base = if late { LATE_INDIRECT_OFFSET } else { 0 };
            for (g, b, material) in &ses {
                if material.is_glass() {
                    continue;
                }
                pass.set_pipeline(&self.ses.opaque);
                pass.set_bind_group(1, &b.surface, &[]);
                for (view, page) in b.views.iter().zip(&b.pages.pages) {
                    pass.set_bind_group(3, view, &[]);
                    pass.draw_indirect(&page.indirect, base);
                }
                pass.set_pipeline(&self.point_pipeline);
                for page in &b.pages.pages {
                    pass.set_bind_group(1, &page.draw, &[]);
                    pass.draw_indirect(&page.indirect, base + POINTS_INDIRECT_OFFSET);
                }
                if !late && cut_active {
                    pass.set_pipeline(&self.ses.atom_cap);
                    pass.set_bind_group(1, &b.surface, &[]);
                    pass.draw(0..g.atom_count * 6, 0..1);
                }
            }
        };

        // Ribbons only: a glycan mesh is rebuilt whole on the CPU every
        // frame and carries no GPU-side spline or level-of-detail state.
        let ribbons: Vec<(&CartoonGpu, &wgpu::BindGroup, &wgpu::BindGroup)> = cartoons
            .iter()
            .filter_map(|item| match (item.mesh, item.bindings) {
                (CartoonMesh::Ribbon(g), CartoonBindings::Ribbon { lod, frame, .. }) => {
                    Some((g, lod, frame))
                }
                _ => None,
            })
            .collect();

        // Cartoons whose atoms moved: sections from the new spline first.
        let moved: Vec<_> = ribbons
            .iter()
            .copied()
            .filter(|(g, ..)| g.take_stale())
            .collect();
        if !moved.is_empty() {
            let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: Some("cartoon frame"),
                timestamp_writes: None,
            });
            pass.set_pipeline(&self.cartoon_frame_pipeline);
            for (g, _, frame) in moved {
                pass.set_bind_group(0, frame, &[]);
                pass.dispatch_workgroups(g.section_count.div_ceil(64), 1, 1);
            }
        }
        if !ribbons.is_empty() {
            for (g, ..) in ribbons.iter().copied() {
                self.ctx.queue.write_buffer(
                    &g.indirect,
                    0,
                    bytemuck::cast_slice(&CARTOON_INDIRECT_RESET),
                );
            }
            let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: Some("cartoon lod"),
                timestamp_writes: None,
            });
            pass.set_pipeline(&self.cartoon_lod_pipeline);
            pass.set_bind_group(0, &self.targets.frame_bind_group, &[]);
            for (g, lod, _) in ribbons.iter().copied() {
                pass.set_bind_group(1, lod, &[]);
                pass.dispatch_workgroups(g.span_count.div_ceil(64), 1, 1);
            }
        }
        cull(
            encoder,
            &self.cull_clusters_pipeline,
            &self.cull_pipeline,
            &self.cull_bonds_pipeline,
            compute_stamps(&self.timer, 0),
        );

        // Every glass item's id, before anything opaque draws: the first
        // (and, here, only) write to `id_view` this frame, so the opaque
        // pass below can `Load` instead of `Clear` and simply draw over
        // it wherever it draws — an opaque hit always wins because it
        // always draws after this, regardless of which is nearer. Depth
        // is `glass_depth`, cleared to the far plane and shared by no one
        // else at this point in the frame, so overlapping glass atoms
        // still resolve to the nearest one. Only phase 1's survivors: see
        // `draw_atoms_glass`'s doc for the resulting gap.
        if any_pick_glass {
            let mut pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("glass id"),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view: &self.targets.id_view,
                    resolve_target: None,
                    ops: wgpu::Operations {
                        load: wgpu::LoadOp::Clear(wgpu::Color::TRANSPARENT),
                        store: wgpu::StoreOp::Store,
                    },
                    depth_slice: None,
                })],
                depth_stencil_attachment: Some(wgpu::RenderPassDepthStencilAttachment {
                    view: &self.targets.glass_depth_view,
                    depth_ops: Some(wgpu::Operations {
                        load: wgpu::LoadOp::Clear(0.0),
                        store: wgpu::StoreOp::Store,
                    }),
                    stencil_ops: None,
                }),
                ..Default::default()
            });
            pass.set_bind_group(0, &self.targets.frame_bind_group, &[]);
            draw_atoms_glass(&mut pass, true, false);
            draw_cartoons_glass(&mut pass, true);
            drop(pass);
            // Only while `fast_glass` is on: this reference is otherwise
            // never read (`glass_cull_margin` stays 0 above), and
            // skipping the copy means turning `fast_glass` back on after
            // a while sees a stale-but-safe zero reference (no culling)
            // for one frame rather than a pose from whenever it was last
            // on, which could be very wrong. Copy out this frame's
            // nearest-glass depth before anything else reuses
            // `glass_depth` as scratch (the Skin/SES owner prepass,
            // below) -- next frame's glass self-occlusion cull reference
            // (cull.wgsl's `glass_hidden`, docs/RENDERING.md's
            // "Transparency").
            if settings.fast_glass {
                let mut pyramid = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                    label: Some("glass hiz"),
                    timestamp_writes: None,
                });
                pyramid.set_pipeline(&self.hiz_copy_pipeline);
                pyramid.set_bind_group(0, &self.targets.glass_hiz_copy_bind_group, &[]);
                pyramid.dispatch_workgroups(
                    self.targets.width.div_ceil(HIZ_WORKGROUP),
                    self.targets.height.div_ceil(HIZ_WORKGROUP),
                    1,
                );
            }
        }

        let draw_stamps = self
            .timer
            .as_ref()
            .map(|t| wgpu::RenderPassTimestampWrites {
                query_set: t.query_set(),
                beginning_of_pass_write_index: Some(2),
                end_of_pass_write_index: Some(3),
            });
        {
            fn attachment(
                view: &wgpu::TextureView,
                clear: wgpu::Color,
            ) -> Option<wgpu::RenderPassColorAttachment<'_>> {
                Some(wgpu::RenderPassColorAttachment {
                    view,
                    resolve_target: None,
                    ops: wgpu::Operations {
                        load: wgpu::LoadOp::Clear(clear),
                        store: wgpu::StoreOp::Store,
                    },
                    depth_slice: None,
                })
            }
            // Cleared to NO_PICK (0) here, the one float clear value every
            // backend converts to the same integer -- unless the "glass
            // id" pass above already populated it with a transparent
            // item's own ids, in which case this pass keeps them wherever
            // it does not itself draw (an opaque hit always overwrites,
            // regardless of depth: see the "glass id" pass's doc).
            let id_attachment = Some(wgpu::RenderPassColorAttachment {
                view: &self.targets.id_view,
                resolve_target: None,
                ops: wgpu::Operations {
                    load: if any_pick_glass {
                        wgpu::LoadOp::Load
                    } else {
                        wgpu::LoadOp::Clear(wgpu::Color::TRANSPARENT)
                    },
                    store: wgpu::StoreOp::Store,
                },
                depth_slice: None,
            });
            let mut pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("opaque"),
                color_attachments: &[
                    attachment(&self.targets.color_view, settings.background),
                    attachment(&self.targets.normal_view, wgpu::Color::BLACK),
                    id_attachment,
                ],
                depth_stencil_attachment: Some(wgpu::RenderPassDepthStencilAttachment {
                    view: &self.targets.depth_view,
                    depth_ops: Some(wgpu::Operations {
                        load: wgpu::LoadOp::Clear(0.0),
                        store: wgpu::StoreOp::Store,
                    }),
                    stencil_ops: None,
                }),
                timestamp_writes: draw_stamps,
                ..Default::default()
            });
            pass.set_bind_group(0, &self.targets.frame_bind_group, &[]);
            draw_atoms(&mut pass, false);
            for item in cartoons {
                if item.material.is_glass() {
                    continue;
                }
                pass.set_bind_group(1, item.bindings.draw(), &[]);
                match item.mesh {
                    CartoonMesh::Ribbon(g) => {
                        pass.set_pipeline(&self.cartoon_pipeline);
                        pass.set_index_buffer(
                            self.cartoon_indices.slice(..),
                            wgpu::IndexFormat::Uint16,
                        );
                        pass.draw_indexed_indirect(&g.indirect, 0);
                    }
                    CartoonMesh::Glycan(g) => {
                        pass.set_pipeline(&self.glycan_pipeline);
                        pass.draw(0..g.vertex_count, 0..1);
                    }
                }
            }
            if !gaussian_surfaces.is_empty() {
                pass.set_pipeline(&self.gaussian_surface_pipeline);
                for item in gaussian_surfaces {
                    if item.material.is_glass() {
                        continue;
                    }
                    pass.set_bind_group(1, item.bindings, &[]);
                    // Full-screen triangle, no vertex/index buffer -- see
                    // shaders/gaussian_surface.wgsl's `vs_gaussian_surface`.
                    pass.draw(0..3, 0..1);
                }
            }
            for item in patch_surfaces {
                if let PatchSurface::Skin(gpu, bindings) = item.surface {
                    if gpu.patch_count == 0 || item.material.is_glass() {
                        continue;
                    }
                    pass.set_pipeline(&self.skin.opaque);
                    pass.set_bind_group(1, bindings, &[]);
                    // Six vertices per patch billboard (shaders/billboard.wgsl).
                    pass.draw(0..gpu.patch_count * 6, 0..1);
                    if cut_active {
                        pass.set_pipeline(&self.skin.atom_cap);
                        pass.set_bind_group(1, bindings, &[]);
                        pass.draw(0..gpu.atom_count * 6, 0..1);
                    }
                }
            }
            draw_ses(&mut pass, false);
        }

        // Occlusion phase 2: pyramid from what phase 1 drew, re-test what
        // it rejected, draw the survivors on top (see shaders/cull.wgsl).
        if settings.occlusion_culling && !(items.is_empty() && ses.is_empty()) {
            hiz(encoder, None);
            cull(
                encoder,
                &self.cull_clusters_late_pipeline,
                &self.cull_late_pipeline,
                &self.cull_bonds_late_pipeline,
                None,
            );
            let load = |view| {
                Some(wgpu::RenderPassColorAttachment {
                    view,
                    resolve_target: None,
                    ops: wgpu::Operations {
                        load: wgpu::LoadOp::Load,
                        store: wgpu::StoreOp::Store,
                    },
                    depth_slice: None,
                })
            };
            let mut pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("opaque late"),
                color_attachments: &[
                    load(&self.targets.color_view),
                    load(&self.targets.normal_view),
                    load(&self.targets.id_view),
                ],
                depth_stencil_attachment: Some(wgpu::RenderPassDepthStencilAttachment {
                    view: &self.targets.depth_view,
                    depth_ops: Some(wgpu::Operations {
                        load: wgpu::LoadOp::Load,
                        store: wgpu::StoreOp::Store,
                    }),
                    stencil_ops: None,
                }),
                ..Default::default()
            });
            pass.set_bind_group(0, &self.targets.frame_bind_group, &[]);
            draw_atoms(&mut pass, true);
            draw_ses(&mut pass, true);
        }
        {
            let mut pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("cap depth"),
                color_attachments: &[],
                depth_stencil_attachment: Some(wgpu::RenderPassDepthStencilAttachment {
                    view: &self.targets.depth_view,
                    depth_ops: Some(wgpu::Operations {
                        load: wgpu::LoadOp::Load,
                        store: wgpu::StoreOp::Store,
                    }),
                    stencil_ops: None,
                }),
                ..Default::default()
            });
            pass.set_pipeline(&self.cap_depth_pipeline);
            pass.set_bind_group(0, &self.targets.frame_bind_group, &[]);
            pass.set_bind_group(1, &self.targets.cap_depth_bind_group, &[]);
            pass.draw(0..3, 0..1);
        }
        // The pyramid next frame's phase 1 tests against.
        hiz(encoder, compute_stamps(&self.timer, 4));

        // Post-process: AO -> `ao`; `color` + `ao` -> composite (outline)
        // pass -> `outlined` -> FXAA or copy -> `fxaa_view` (the display).
        // Every pass always runs, so no texture in the chain ever changes
        // identity — see the comment on `display_view` for why that
        // matters. Timed as one span (timestamps 6..7).
        let post = |encoder: &mut wgpu::CommandEncoder,
                    label: &str,
                    pipeline: &wgpu::RenderPipeline,
                    bind_group: &wgpu::BindGroup,
                    target: &wgpu::TextureView,
                    end_stamp: Option<u32>| {
            let mut pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some(label),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view: target,
                    resolve_target: None,
                    ops: wgpu::Operations {
                        load: wgpu::LoadOp::Clear(wgpu::Color::BLACK),
                        store: wgpu::StoreOp::Store,
                    },
                    depth_slice: None,
                })],
                timestamp_writes: end_stamp.and_then(|i| {
                    self.timer
                        .as_ref()
                        .map(|t| wgpu::RenderPassTimestampWrites {
                            query_set: t.query_set(),
                            beginning_of_pass_write_index: None,
                            end_of_pass_write_index: Some(i),
                        })
                }),
                ..Default::default()
            });
            pass.set_pipeline(pipeline);
            pass.set_bind_group(0, bind_group, &[]);
            pass.draw(0..3, 0..1);
        };
        {
            // AO into its own target; a cleared-white target when off
            // (the shader returns 1 early), so the composite never reads
            // a stale frame.
            let mut pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("ao"),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view: &self.targets.ao_view,
                    resolve_target: None,
                    ops: wgpu::Operations {
                        load: wgpu::LoadOp::Clear(wgpu::Color::WHITE),
                        store: wgpu::StoreOp::Store,
                    },
                    depth_slice: None,
                })],
                timestamp_writes: self
                    .timer
                    .as_ref()
                    .map(|t| wgpu::RenderPassTimestampWrites {
                        query_set: t.query_set(),
                        beginning_of_pass_write_index: Some(6),
                        end_of_pass_write_index: None,
                    }),
                ..Default::default()
            });
            if settings.ao > 0.0 || settings.shadows > 0.0 {
                pass.set_pipeline(&self.ao_pipeline);
                pass.set_bind_group(0, &self.targets.ao_bind_group, &[]);
                pass.set_bind_group(1, &self.occlusion_bind_group, &[]);
                pass.draw(0..3, 0..1);
            }
        }
        post(
            encoder,
            "outline",
            &self.outline_pipeline,
            &self.targets.outline_bind_group,
            &self.targets.outlined_view,
            None,
        );
        let glass_gaussian: Vec<_> = gaussian_surfaces
            .iter()
            .filter(|i| i.material.is_glass())
            .collect();
        // A skin's every patch; an SES's survivors of both cull phases.
        let draw_patches = |pass: &mut wgpu::RenderPass<'_>,
                            surface: PatchSurface<'_>,
                            pipeline: &wgpu::RenderPipeline| {
            pass.set_pipeline(pipeline);
            pass.set_bind_group(1, surface.bind_group(), &[]);
            match surface {
                PatchSurface::Skin(gpu, _) => pass.draw(0..gpu.patch_count * 6, 0..1),
                PatchSurface::Ses(_, b) => {
                    for (view, page) in b.views.iter().zip(&b.pages.pages) {
                        pass.set_bind_group(3, view, &[]);
                        pass.draw_indirect(&page.indirect, 0);
                        pass.draw_indirect(&page.indirect, LATE_INDIRECT_OFFSET);
                    }
                }
            }
        };
        let glass_patches: Vec<_> = patch_surfaces
            .iter()
            .filter(|i| i.material.is_glass() && i.surface.patch_count() > 0)
            .collect();
        if any_pick_glass || !glass_gaussian.is_empty() || !glass_patches.is_empty() {
            encoder.copy_texture_to_texture(
                self.targets.depth.as_image_copy(),
                self.targets.glass_depth.as_image_copy(),
                self.targets.glass_depth.size(),
            );
            if !glass_patches.is_empty() {
                let mut pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                    label: Some("glass owner"),
                    color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                        view: &self.targets.glass_owner_view,
                        resolve_target: None,
                        ops: wgpu::Operations {
                            load: wgpu::LoadOp::Clear(wgpu::Color::TRANSPARENT),
                            store: wgpu::StoreOp::Store,
                        },
                        depth_slice: None,
                    })],
                    depth_stencil_attachment: Some(wgpu::RenderPassDepthStencilAttachment {
                        view: &self.targets.glass_depth_view,
                        depth_ops: Some(wgpu::Operations {
                            load: wgpu::LoadOp::Load,
                            store: wgpu::StoreOp::Discard,
                        }),
                        stencil_ops: None,
                    }),
                    ..Default::default()
                });
                pass.set_bind_group(0, &self.targets.frame_bind_group, &[]);
                for item in &glass_patches {
                    draw_patches(
                        &mut pass,
                        item.surface,
                        &self.patch_pipelines(item.surface).owner,
                    );
                }
            }
            let mut pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("glass"),
                color_attachments: &[
                    Some(wgpu::RenderPassColorAttachment {
                        view: &self.targets.accum_view,
                        resolve_target: None,
                        ops: wgpu::Operations {
                            load: wgpu::LoadOp::Clear(wgpu::Color::TRANSPARENT),
                            store: wgpu::StoreOp::Store,
                        },
                        depth_slice: None,
                    }),
                    Some(wgpu::RenderPassColorAttachment {
                        view: &self.targets.reveal_view,
                        resolve_target: None,
                        ops: wgpu::Operations {
                            load: wgpu::LoadOp::Clear(wgpu::Color::WHITE),
                            store: wgpu::StoreOp::Store,
                        },
                        depth_slice: None,
                    }),
                ],
                depth_stencil_attachment: Some(wgpu::RenderPassDepthStencilAttachment {
                    view: &self.targets.depth_view,
                    depth_ops: Some(wgpu::Operations {
                        load: wgpu::LoadOp::Load,
                        store: wgpu::StoreOp::Store,
                    }),
                    stencil_ops: None,
                }),
                ..Default::default()
            });
            pass.set_bind_group(0, &self.targets.frame_bind_group, &[]);
            draw_atoms_glass(&mut pass, false, false);
            draw_atoms_glass(&mut pass, false, true);
            draw_cartoons_glass(&mut pass, false);
            pass.set_bind_group(2, &self.targets.glass_owner_bind_group, &[]);
            pass.set_pipeline(&self.gaussian_glass_pipeline);
            for item in glass_gaussian {
                pass.set_bind_group(1, item.bindings, &[]);
                pass.draw(0..3, 0..1);
            }
            for item in glass_patches {
                draw_patches(
                    &mut pass,
                    item.surface,
                    &self.patch_pipelines(item.surface).glass,
                );
            }
            drop(pass);
            let mut pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("glass resolve"),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view: &self.targets.outlined_view,
                    resolve_target: None,
                    ops: wgpu::Operations {
                        load: wgpu::LoadOp::Load,
                        store: wgpu::StoreOp::Store,
                    },
                    depth_slice: None,
                })],
                ..Default::default()
            });
            pass.set_pipeline(&self.glass_resolve_pipeline);
            pass.set_bind_group(0, &self.targets.glass_resolve_bind_group, &[]);
            pass.draw(0..3, 0..1);
        }
        if settings.dof > 0.0 {
            let focus = camera.pivot.unwrap_or(camera.target);
            let focus_depth = -(camera.view() * focus.extend(1.0)).z;
            let uniform = DofUniform {
                inv_proj: proj.inverse().to_cols_array_2d(),
                focus: focus_depth.max(1e-3),
                max_radius: DOF_MAX_RADIUS * settings.dof.min(1.0) * self.targets.height as f32
                    / 1080.0,
                range: dof_range.max(1.0),
                _pad: 0.0,
            };
            self.ctx
                .queue
                .write_buffer(&self.dof_buffer, 0, bytemuck::bytes_of(&uniform));
            encoder.copy_texture_to_texture(
                self.targets.outlined.as_image_copy(),
                self.targets.dof_source.as_image_copy(),
                self.targets.outlined.size(),
            );
            post(
                encoder,
                "dof",
                &self.dof_pipeline,
                &self.targets.dof_bind_group,
                &self.targets.outlined_view,
                None,
            );
        }
        let aa = if settings.fxaa {
            &self.fxaa_pipeline
        } else {
            &self.passthrough_pipeline
        };
        post(
            encoder,
            "post",
            aa,
            &self.targets.fxaa_bind_group,
            &self.targets.fxaa_view,
            Some(7),
        );
        // After FXAA so the halo stays crisp; display only, so exports
        // (`read_color`) never show it.
        self.selection.draw(encoder, &self.targets.fxaa_view);

        if let Some(timer) = &mut self.timer {
            if let Some(times) = timer.end_frame(encoder) {
                self.last_times = Some(times);
            }
        }
    }

    /// Must be called after the frame's command buffer is submitted so the
    /// timestamp readback for that frame can be scheduled.
    pub fn after_submit(&mut self) {
        if let Some(timer) = &mut self.timer {
            timer.after_submit(&self.ctx);
        }
    }

    fn read_texture(&self, texture: &wgpu::Texture) -> Vec<u8> {
        read_texture_rgba8(&self.ctx, texture)
    }

    /// Blocking readback of the color target after the outline pass and
    /// before FXAA, as tightly packed RGBA8 rows. For tests and screenshots
    /// (SSAA export wants outlines but not FXAA), not the frame loop.
    pub fn read_color(&self) -> Vec<u8> {
        self.read_texture(&self.targets.outlined)
    }

    /// Blocking readback of the depth target (reversed Z: 0 is the far
    /// plane and the background), row-major. For tests.
    pub fn read_depth(&self) -> Vec<f32> {
        self.read_texture(&self.targets.depth)
            .chunks_exact(4)
            .map(|b| f32::from_le_bytes(b.try_into().unwrap()))
            .collect()
    }

    /// Blocking readback of the pick-id target (see [`ID_FORMAT`]),
    /// row-major. For tests.
    pub fn read_ids(&self) -> Vec<u32> {
        self.read_texture(&self.targets.id_texture)
            .chunks_exact(4)
            .map(|b| u32::from_le_bytes(b.try_into().unwrap()))
            .collect()
    }

    /// Blocking readback of `display_view` (after FXAA, if it ran) as
    /// tightly packed RGBA8 rows. For tests and diagnostics; the frame
    /// loop never reads this back, it just shows it.
    pub fn read_display_color(&self) -> Vec<u8> {
        self.read_texture(&self.targets.fxaa_color)
    }

    /// Blocking readback of what is under pixel `(x, y)` in the id buffer of
    /// the most recently rendered frame: an atom, a bond cylinder (only
    /// ball-and-stick draws any), or `None` for background.
    ///
    /// Reuses that frame's already-culled draw rather than re-running
    /// cull: an atom the occlusion pass dropped this frame is not drawn,
    /// so — correctly — it cannot be picked either. Call this right after
    /// `render()`/`render_all()`, before the next frame overwrites the id
    /// buffer.
    pub fn pick(&self, x: u32, y: u32) -> Option<Pick> {
        let (width, height) = self.size();
        if x >= width || y >= height {
            return None;
        }
        let bytes_per_row = wgpu::COPY_BYTES_PER_ROW_ALIGNMENT;
        let device = &self.ctx.device;
        let buffer = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("pick readback"),
            size: bytes_per_row as u64,
            usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        let mut encoder = device.create_command_encoder(&Default::default());
        encoder.copy_texture_to_buffer(
            wgpu::TexelCopyTextureInfo {
                texture: &self.targets.id_texture,
                mip_level: 0,
                origin: wgpu::Origin3d { x, y, z: 0 },
                aspect: wgpu::TextureAspect::All,
            },
            wgpu::TexelCopyBufferInfo {
                buffer: &buffer,
                layout: wgpu::TexelCopyBufferLayout {
                    offset: 0,
                    bytes_per_row: Some(bytes_per_row),
                    rows_per_image: None,
                },
            },
            wgpu::Extent3d {
                width: 1,
                height: 1,
                depth_or_array_layers: 1,
            },
        );
        self.ctx.queue.submit([encoder.finish()]);

        let slice = buffer.slice(..);
        let (tx, rx) = std::sync::mpsc::channel();
        slice.map_async(wgpu::MapMode::Read, move |r| {
            let _ = tx.send(r);
        });
        let _ = device.poll(wgpu::PollType::wait_indefinitely());
        rx.recv().expect("map callback").expect("map readback");
        let data = slice.get_mapped_range().expect("mapped range");
        let id = u32::from_le_bytes(data[0..4].try_into().unwrap());
        drop(data);
        buffer.unmap();
        self.resolve_pick(id)
    }

    /// Maps a raw id from the picking target back to the draw item and
    /// local index it came from, using the ranges `render_all` assigned.
    fn resolve_pick(&self, id: u32) -> Option<Pick> {
        if id == NO_PICK {
            return None;
        }
        let is_bond = id & BOND_ID_FLAG != 0;
        // Shaders write `id + 1` so that 0 is free for "nothing".
        let global = (id & !BOND_ID_FLAG) - 1;
        self.pick_ranges
            .iter()
            .enumerate()
            .find_map(|(item, range)| {
                let (base, count) = if is_bond {
                    (range.bond_base, range.bond_count)
                } else {
                    (range.atom_base, range.atom_count)
                };
                (global >= base && global < base + count).then(|| {
                    let local = global - base;
                    if is_bond {
                        Pick::Bond { item, bond: local }
                    } else {
                        Pick::Atom { item, atom: local }
                    }
                })
            })
    }
}

impl Targets {
    #[allow(clippy::too_many_arguments)]
    fn new(
        device: &wgpu::Device,
        width: u32,
        height: u32,
        layouts: &Layouts,
        camera_buffer: &wgpu::Buffer,
        lighting_buffer: &wgpu::Buffer,
        outline_buffer: &wgpu::Buffer,
        dof_buffer: &wgpu::Buffer,
        blit_sampler: &wgpu::Sampler,
    ) -> Self {
        let size = wgpu::Extent3d {
            width,
            height,
            depth_or_array_layers: 1,
        };
        let make = |label: &str, format, usage, mip_level_count| {
            device.create_texture(&wgpu::TextureDescriptor {
                label: Some(label),
                size,
                mip_level_count,
                sample_count: 1,
                dimension: wgpu::TextureDimension::D2,
                format,
                usage,
                view_formats: &[],
            })
        };
        let attach = wgpu::TextureUsages::RENDER_ATTACHMENT;
        let color = make(
            "color",
            COLOR_FORMAT,
            attach | wgpu::TextureUsages::TEXTURE_BINDING | wgpu::TextureUsages::COPY_SRC,
            1,
        );
        let normal = make(
            "normal",
            NORMAL_FORMAT,
            attach | wgpu::TextureUsages::TEXTURE_BINDING,
            1,
        );
        let depth = make(
            "depth",
            DEPTH_FORMAT,
            attach | wgpu::TextureUsages::TEXTURE_BINDING | wgpu::TextureUsages::COPY_SRC,
            1,
        );
        let glass_depth = make(
            "glass depth",
            DEPTH_FORMAT,
            attach | wgpu::TextureUsages::COPY_DST | wgpu::TextureUsages::TEXTURE_BINDING,
            1,
        );
        let glass_depth_view = glass_depth.create_view(&Default::default());
        // Single texel, no mip chain (`glass_hidden` samples its own
        // atom's projected pixel only) -- reuses `hiz_copy_pipeline`,
        // level 0 of the opaque pyramid's own build.
        let glass_hiz = make(
            "glass depth pyramid",
            HIZ_FORMAT,
            wgpu::TextureUsages::TEXTURE_BINDING | wgpu::TextureUsages::STORAGE_BINDING,
            1,
        );
        let glass_hiz_view = glass_hiz.create_view(&Default::default());
        let glass_hiz_copy_bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("glass depth pyramid"),
            layout: &layouts.hiz_copy,
            entries: &[
                wgpu::BindGroupEntry {
                    binding: 0,
                    resource: wgpu::BindingResource::TextureView(&glass_depth_view),
                },
                wgpu::BindGroupEntry {
                    binding: 1,
                    resource: wgpu::BindingResource::TextureView(&glass_hiz_view),
                },
            ],
        });
        let glass_owner_view = make(
            "glass owner",
            OWNER_FORMAT,
            attach | wgpu::TextureUsages::TEXTURE_BINDING,
            1,
        )
        .create_view(&Default::default());
        let glass_owner_bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("glass owner"),
            layout: &layouts.glass_owner,
            entries: &[wgpu::BindGroupEntry {
                binding: 0,
                resource: wgpu::BindingResource::TextureView(&glass_owner_view),
            }],
        });
        let id_texture = make(
            "id",
            ID_FORMAT,
            attach | wgpu::TextureUsages::TEXTURE_BINDING | wgpu::TextureUsages::COPY_SRC,
            1,
        );
        // Outline pass output; what `read_color` (screenshots) reads.
        let outlined = make(
            "outlined",
            COLOR_FORMAT,
            attach | wgpu::TextureUsages::TEXTURE_BINDING | wgpu::TextureUsages::COPY_SRC,
            1,
        );
        let outlined_view = outlined.create_view(&Default::default());
        // The FXAA (or passthrough) output: the display. `read_display_color`
        // (tests, diagnostics) needs COPY_SRC here too.
        let gamma_format = COLOR_FORMAT.remove_srgb_suffix();
        let fxaa_color = device.create_texture(&wgpu::TextureDescriptor {
            label: Some("fxaa"),
            size,
            mip_level_count: 1,
            sample_count: 1,
            dimension: wgpu::TextureDimension::D2,
            format: COLOR_FORMAT,
            usage: attach | wgpu::TextureUsages::TEXTURE_BINDING | wgpu::TextureUsages::COPY_SRC,
            view_formats: &[gamma_format],
        });
        let fxaa_view = fxaa_color.create_view(&Default::default());
        let fxaa_gamma_view = fxaa_color.create_view(&wgpu::TextureViewDescriptor {
            format: Some(gamma_format),
            ..Default::default()
        });
        let accum_view = make(
            "glass accum",
            ACCUM_FORMAT,
            attach | wgpu::TextureUsages::TEXTURE_BINDING,
            1,
        )
        .create_view(&Default::default());
        let reveal_view = make(
            "glass reveal",
            REVEAL_FORMAT,
            attach | wgpu::TextureUsages::TEXTURE_BINDING,
            1,
        )
        .create_view(&Default::default());
        let glass_resolve_bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("glass resolve"),
            layout: &layouts.glass_resolve,
            entries: &[
                wgpu::BindGroupEntry {
                    binding: 0,
                    resource: wgpu::BindingResource::TextureView(&accum_view),
                },
                wgpu::BindGroupEntry {
                    binding: 1,
                    resource: wgpu::BindingResource::TextureView(&reveal_view),
                },
            ],
        });

        let hiz_mip_count = (width.max(height) as f32).log2().floor() as u32 + 1;
        let hiz = make(
            "depth pyramid",
            HIZ_FORMAT,
            wgpu::TextureUsages::TEXTURE_BINDING | wgpu::TextureUsages::STORAGE_BINDING,
            hiz_mip_count,
        );
        let level_view = |level: u32| {
            hiz.create_view(&wgpu::TextureViewDescriptor {
                label: Some("depth pyramid level"),
                base_mip_level: level,
                mip_level_count: Some(1),
                ..Default::default()
            })
        };
        let depth_view = depth.create_view(&Default::default());
        let mut hiz_bind_groups = Vec::with_capacity(hiz_mip_count as usize);
        for level in 0..hiz_mip_count {
            let dst = level_view(level);
            let (layout, src) = if level == 0 {
                (&layouts.hiz_copy, depth_view.clone())
            } else {
                (&layouts.hiz_down, level_view(level - 1))
            };
            hiz_bind_groups.push(device.create_bind_group(&wgpu::BindGroupDescriptor {
                label: Some("depth pyramid level"),
                layout,
                entries: &[
                    wgpu::BindGroupEntry {
                        binding: 0,
                        resource: wgpu::BindingResource::TextureView(&src),
                    },
                    wgpu::BindGroupEntry {
                        binding: 1,
                        resource: wgpu::BindingResource::TextureView(&dst),
                    },
                ],
            }));
        }
        let hiz_view = hiz.create_view(&Default::default());
        let frame_bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("frame"),
            layout: &layouts.frame,
            entries: &[
                wgpu::BindGroupEntry {
                    binding: 0,
                    resource: camera_buffer.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 1,
                    resource: wgpu::BindingResource::TextureView(&hiz_view),
                },
                wgpu::BindGroupEntry {
                    binding: 2,
                    resource: lighting_buffer.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 3,
                    resource: wgpu::BindingResource::TextureView(&glass_hiz_view),
                },
            ],
        });

        let color_view = color.create_view(&Default::default());
        let normal_view = normal.create_view(&Default::default());
        // Half resolution: a quarter of the samples; the composite pass's
        // depth-aware blur upsamples it.
        let ao = device.create_texture(&wgpu::TextureDescriptor {
            label: Some("ao"),
            size: wgpu::Extent3d {
                width: width.div_ceil(2),
                height: height.div_ceil(2),
                depth_or_array_layers: 1,
            },
            mip_level_count: 1,
            sample_count: 1,
            dimension: wgpu::TextureDimension::D2,
            format: AO_FORMAT,
            usage: attach | wgpu::TextureUsages::TEXTURE_BINDING,
            view_formats: &[],
        });
        let ao_view = ao.create_view(&Default::default());
        let dof_source = make(
            "dof source",
            COLOR_FORMAT,
            wgpu::TextureUsages::TEXTURE_BINDING | wgpu::TextureUsages::COPY_DST,
            1,
        );
        let dof_bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("dof"),
            layout: &layouts.dof,
            entries: &[
                wgpu::BindGroupEntry {
                    binding: 0,
                    resource: wgpu::BindingResource::TextureView(
                        &dof_source.create_view(&Default::default()),
                    ),
                },
                wgpu::BindGroupEntry {
                    binding: 1,
                    resource: wgpu::BindingResource::TextureView(&depth_view),
                },
                wgpu::BindGroupEntry {
                    binding: 2,
                    resource: dof_buffer.as_entire_binding(),
                },
            ],
        });
        let ao_bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("ao"),
            layout: &layouts.ao,
            entries: &[
                wgpu::BindGroupEntry {
                    binding: 0,
                    resource: wgpu::BindingResource::TextureView(&normal_view),
                },
                wgpu::BindGroupEntry {
                    binding: 1,
                    resource: wgpu::BindingResource::TextureView(&depth_view),
                },
                wgpu::BindGroupEntry {
                    binding: 2,
                    resource: outline_buffer.as_entire_binding(),
                },
            ],
        });
        let outline_bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("outline"),
            layout: &layouts.outline,
            entries: &[
                wgpu::BindGroupEntry {
                    binding: 0,
                    resource: wgpu::BindingResource::TextureView(&color_view),
                },
                wgpu::BindGroupEntry {
                    binding: 1,
                    resource: wgpu::BindingResource::TextureView(&normal_view),
                },
                wgpu::BindGroupEntry {
                    binding: 2,
                    resource: wgpu::BindingResource::TextureView(&depth_view),
                },
                wgpu::BindGroupEntry {
                    binding: 3,
                    resource: outline_buffer.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 4,
                    resource: wgpu::BindingResource::TextureView(&ao_view),
                },
            ],
        });
        let blit_bind_group = |label: &str, view: &wgpu::TextureView| {
            device.create_bind_group(&wgpu::BindGroupDescriptor {
                label: Some(label),
                layout: &layouts.blit,
                entries: &[
                    wgpu::BindGroupEntry {
                        binding: 0,
                        resource: wgpu::BindingResource::TextureView(view),
                    },
                    wgpu::BindGroupEntry {
                        binding: 1,
                        resource: wgpu::BindingResource::Sampler(blit_sampler),
                    },
                ],
            })
        };
        let fxaa_bind_group = blit_bind_group("fxaa", &outlined_view);
        let cap_depth_bind_group = blit_bind_group("cap depth", &normal_view);

        Self {
            width,
            height,
            color_view,
            normal_view,
            depth_view,
            id_view: id_texture.create_view(&Default::default()),
            id_texture,
            color,
            outlined,
            outlined_view,
            outline_bind_group,
            ao_view,
            ao_bind_group,
            fxaa_color,
            fxaa_view,
            fxaa_gamma_view,
            fxaa_bind_group,
            cap_depth_bind_group,
            depth,
            dof_source,
            dof_bind_group,
            glass_depth,
            glass_depth_view,
            glass_hiz_copy_bind_group,
            glass_owner_view,
            glass_owner_bind_group,
            accum_view,
            reveal_view,
            glass_resolve_bind_group,
            hiz_mip_count,
            hiz_bind_groups,
            frame_bind_group,
        }
    }
}

/// Blocking readback of a 4-byte-per-texel 2D texture (it needs
/// `COPY_SRC`) as tightly packed rows, in the texture's own channel order.
/// Diagnostics and exports only, never the frame loop.
pub fn read_texture_rgba8(ctx: &GpuContext, texture: &wgpu::Texture) -> Vec<u8> {
    let (width, height) = (texture.width(), texture.height());
    let bytes_per_row = (width * 4).div_ceil(wgpu::COPY_BYTES_PER_ROW_ALIGNMENT)
        * wgpu::COPY_BYTES_PER_ROW_ALIGNMENT;
    let device = &ctx.device;
    let buffer = device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("color readback"),
        size: bytes_per_row as u64 * height as u64,
        usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
        mapped_at_creation: false,
    });
    let mut encoder = device.create_command_encoder(&Default::default());
    encoder.copy_texture_to_buffer(
        texture.as_image_copy(),
        wgpu::TexelCopyBufferInfo {
            buffer: &buffer,
            layout: wgpu::TexelCopyBufferLayout {
                offset: 0,
                bytes_per_row: Some(bytes_per_row),
                rows_per_image: None,
            },
        },
        wgpu::Extent3d {
            width,
            height,
            depth_or_array_layers: 1,
        },
    );
    ctx.queue.submit([encoder.finish()]);

    let slice = buffer.slice(..);
    let (tx, rx) = std::sync::mpsc::channel();
    slice.map_async(wgpu::MapMode::Read, move |r| {
        let _ = tx.send(r);
    });
    let _ = device.poll(wgpu::PollType::wait_indefinitely());
    rx.recv().expect("map callback").expect("map readback");
    let data = slice.get_mapped_range().expect("mapped range");
    let mut out = Vec::with_capacity((width * height * 4) as usize);
    for row in data.chunks(bytes_per_row as usize) {
        out.extend_from_slice(&row[..(width * 4) as usize]);
    }
    drop(data);
    buffer.unmap();
    out
}

/// View-depth range of the depth cue for a scene sphere whose centre is
/// `center_depth` in front of the eye: clear up to the sphere's front face,
/// fully cued 1.5 radii later. Anchored to the front face (not the centre)
/// so an eye that closes on or enters the structure sees it clearer, and
/// the ramp keeps its width instead of shrinking with the camera distance.
fn fog_range(center_depth: f32, radius: f32) -> (f32, f32) {
    let near = (center_depth - radius).max(0.0);
    (near, near + 1.5 * radius)
}

#[cfg(test)]
mod tests {
    use super::fog_range;

    fn fog_at(depth: f32, (near, far): (f32, f32)) -> f32 {
        ((depth - near) / (far - near)).clamp(0.0, 1.0)
    }

    #[test]
    fn fog_leaves_the_front_face_clear_and_fogs_the_back() {
        let range = fog_range(100.0, 10.0);
        assert_eq!(fog_at(90.0, range), 0.0);
        assert!(fog_at(110.0, range) > 0.5);
    }

    #[test]
    fn fog_lessens_at_a_fixed_point_as_the_eye_closes_in() {
        let point = 60.0_f32;
        let mut last = f32::MAX;
        for eye_to_center in [100.0, 60.0, 30.0, 10.0, 0.0] {
            let depth = point - (100.0 - eye_to_center);
            if depth <= 0.0 {
                continue;
            }
            let fog = fog_at(depth, fog_range(eye_to_center, 10.0));
            assert!(fog <= last + 1e-6, "{fog} > {last} at {eye_to_center}");
            last = fog;
        }
    }

    #[test]
    fn fog_ramp_width_does_not_shrink_inside_the_structure() {
        let (near, far) = fog_range(3.0, 10.0);
        assert_eq!(near, 0.0);
        assert_eq!(far, 15.0);
    }
}
