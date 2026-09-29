//! Device/queue setup shared by the app, the bench, and headless tests.

use std::sync::Arc;

#[derive(thiserror::Error, Debug)]
pub enum GpuError {
    #[error("no compatible GPU adapter found")]
    NoAdapter(#[from] wgpu::RequestAdapterError),
    #[error("no GPU adapter matches `{0}`")]
    NoMatchingAdapter(String),
    #[error("failed to create GPU device")]
    NoDevice(#[from] wgpu::RequestDeviceError),
    #[error("failed to create a window surface")]
    NoSurface(#[from] wgpu::CreateSurfaceError),
}

pub struct GpuContext {
    pub instance: wgpu::Instance,
    pub adapter: wgpu::Adapter,
    pub device: wgpu::Device,
    pub queue: wgpu::Queue,
    /// Whether GPU timestamp queries are available on this device.
    pub timestamps: bool,
    /// Whether fragment shaders may declare conservative depth
    /// (`@early_depth_test`), letting hidden impostor fragments be
    /// rejected before they run.
    pub early_depth: bool,
    /// Compute pipelines built on first use by GPU-side geometry uploads
    /// (they have a context, not a `Renderer`).
    pub(crate) gaussian_bake: std::sync::OnceLock<crate::scene::GaussianBake>,
    pub(crate) ses_volume_bake: std::sync::OnceLock<crate::ses_volume::SesVolumeBake>,
}

/// Backends to try, in order. Each `wgpu::Instance` initializes every
/// backend it is asked for, and on Windows that is where startup time
/// goes: Vulkan alone takes ~0.4 s to enumerate, DX12 ~0.05 s but ~0.5 s
/// to create a device, and OpenGL ~1.3 s for a backend this renderer
/// never picks. So: one native backend first (the one every GPU test has
/// been run on), the rest only if it has no usable adapter.
fn backend_preference() -> [wgpu::Backends; 2] {
    if cfg!(target_os = "macos") {
        [wgpu::Backends::METAL, wgpu::Backends::PRIMARY]
    } else {
        [wgpu::Backends::VULKAN, wgpu::Backends::PRIMARY]
    }
}

/// A surface for a window, tied to the instance that will present to it.
pub type SurfaceFactory<'a> =
    &'a dyn Fn(&wgpu::Instance) -> Result<wgpu::Surface<'static>, wgpu::CreateSurfaceError>;

impl GpuContext {
    /// An instance over `backends`, honouring the `WGPU_*` environment
    /// overrides (`WGPU_BACKEND=dx12`, validation flags).
    pub fn instance_for(backends: wgpu::Backends) -> wgpu::Instance {
        let mut desc = wgpu::InstanceDescriptor::new_without_display_handle();
        desc.backends = backends;
        wgpu::Instance::new(desc.with_env())
    }

    /// The instance a caller should use when it has no window: the
    /// preferred native backend (see `backend_preference`). Headless users
    /// wanting the fallback should go through [`GpuContext::new`].
    pub fn instance() -> wgpu::Instance {
        Self::instance_for(backend_preference()[0])
    }

    /// Headless context: the preferred backend, then everything else.
    pub fn new(adapter_filter: Option<&str>) -> Result<Arc<Self>, GpuError> {
        let mut last = None;
        for backends in backend_preference() {
            match Self::with_instance(Self::instance_for(backends), None, adapter_filter) {
                Ok(ctx) => return Ok(ctx),
                Err(e) => last = Some(e),
            }
        }
        Err(last.expect("at least one backend attempted"))
    }

    /// Context plus a surface for a window: the preferred backend, then
    /// everything else. `make_surface` is called once per attempt because
    /// a surface belongs to the instance that created it.
    pub fn for_window(
        make_surface: SurfaceFactory<'_>,
        adapter_filter: Option<&str>,
    ) -> Result<(Arc<Self>, wgpu::Surface<'static>), GpuError> {
        let mut last = None;
        for backends in backend_preference() {
            let instance = Self::instance_for(backends);
            let surface = match make_surface(&instance) {
                Ok(s) => s,
                Err(e) => {
                    last = Some(GpuError::NoSurface(e));
                    continue;
                }
            };
            match Self::with_instance(instance, Some(&surface), adapter_filter) {
                Ok(ctx) => return Ok((ctx, surface)),
                Err(e) => last = Some(e),
            }
        }
        Err(last.expect("at least one backend attempted"))
    }

    /// `adapter_filter`, when given, selects the first adapter whose name
    /// contains it (case-insensitive) instead of the highest-performance one.
    pub fn with_instance(
        instance: wgpu::Instance,
        surface: Option<&wgpu::Surface<'_>>,
        adapter_filter: Option<&str>,
    ) -> Result<Arc<Self>, GpuError> {
        let adapter = match adapter_filter {
            Some(filter) => {
                let needle = filter.to_lowercase();
                pollster::block_on(instance.enumerate_adapters(wgpu::Backends::all()))
                    .into_iter()
                    .find(|a| {
                        a.get_info().name.to_lowercase().contains(&needle)
                            && surface.is_none_or(|s| a.is_surface_supported(s))
                    })
                    .ok_or_else(|| GpuError::NoMatchingAdapter(filter.to_string()))?
            }
            None => pollster::block_on(instance.request_adapter(&wgpu::RequestAdapterOptions {
                power_preference: wgpu::PowerPreference::HighPerformance,
                compatible_surface: surface,
                force_fallback_adapter: false,
                ..Default::default()
            }))?,
        };

        let timestamps = adapter.features().contains(wgpu::Features::TIMESTAMP_QUERY);
        let mut required_features = wgpu::Features::empty();
        if timestamps {
            required_features |= wgpu::Features::TIMESTAMP_QUERY;
        }
        let early_depth = adapter
            .features()
            .contains(wgpu::Features::SHADER_EARLY_DEPTH_TEST);
        if early_depth {
            required_features |= wgpu::Features::SHADER_EARLY_DEPTH_TEST;
        }

        // wgpu's defaults cap every buffer at 256 MB and every storage
        // binding at 128 MB; a 4M-atom cartoon mesh alone is ~480 MB
        // (8GLV). Ask for what the adapter really has.
        let supported = adapter.limits();
        let mut required_limits = wgpu::Limits::default().using_resolution(supported.clone());
        required_limits.max_buffer_size = supported.max_buffer_size;
        required_limits.max_storage_buffer_binding_size = supported.max_storage_buffer_binding_size;
        required_limits.max_storage_buffers_per_shader_stage =
            supported.max_storage_buffers_per_shader_stage;
        let (device, queue) =
            pollster::block_on(adapter.request_device(&wgpu::DeviceDescriptor {
                label: Some("vizviz"),
                required_features,
                required_limits,
                ..Default::default()
            }))?;

        Ok(Arc::new(Self {
            instance,
            adapter,
            device,
            queue,
            timestamps,
            early_depth,
            gaussian_bake: Default::default(),
            ses_volume_bake: Default::default(),
        }))
    }

    pub fn adapter_name(&self) -> String {
        let info = self.adapter.get_info();
        format!("{} [{:?}]", info.name, info.backend)
    }
}
