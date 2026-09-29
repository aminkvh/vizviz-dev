//! The running application: everything that exists once a window and a GPU
//! context do. Split out of `main.rs`, which is the entry point (CLI
//! parsing, the winit event loop glue) and shouldn't also be the place
//! that knows how a frame gets rendered.

use std::collections::{HashSet, VecDeque};
use std::sync::Arc;
use std::time::{Duration, Instant};

use egui_dock::{DockArea, DockState};
use vv_render::{AdaptiveLod, Camera, GpuContext, RenderSettings, Renderer};
use vv_scene::{Command, CommandHistory, Scene, StructureId};
use winit::window::Window;

use crate::gpu_cache::{frame_structure, DrawSource, GpuCache, Waker};
use crate::layout::{self, Tab};
use crate::prefetch::Prefetch;
use crate::ui::{
    color32_from_wgpu, label_draws, measurement_draws, project_to_pixel, timeline_target, AppUi,
    ConsoleState, DialogState, ExportDialogState, ExportRequest, ExportStructureDialogState,
    FetchDialogState, LayoutRequest, Notice, PaletteState, Readout, TimelineState, ViewSettings,
    WorkspaceDialogState, CAPTION_FONT_PX, CAPTION_TEXT_COLOR, LABEL_DOT_COLOR, LABEL_DOT_RADIUS,
    LABEL_FONT_PX, LABEL_TEXT_COLOR, MEASURE_DASH, MEASURE_LINE_COLOR,
};
use crate::Options;

mod movie_run;

pub(crate) const DEFAULT_VIEWPORT: (u32, u32) = (1024, 768);
const RESIZE_DEBOUNCE_FRAMES: u32 = 4;

pub(crate) struct State {
    pub(crate) window: Arc<Window>,
    pub(crate) surface: wgpu::Surface<'static>,
    pub(crate) config: wgpu::SurfaceConfiguration,
    pub(crate) ctx: Arc<GpuContext>,

    renderer: Renderer,
    gpu_cache: GpuCache,
    timeline: TimelineState,
    movie: crate::movie_state::MovieState,
    movie_export: Option<movie_run::ExportJob>,
    /// Reads a playing trajectory's next frames ahead of the timeline.
    prefetch: Prefetch,
    /// Asks for a redraw from any thread (a worker finished).
    waker: Waker,
    /// Multi-frame structures the timeline has already been opened for,
    /// so closing it stays closed until the next such load.
    timeline_offered: HashSet<StructureId>,
    /// What each item the last frame drew stands for; a `Pick { item }`
    /// indexes this to find its structure and real atom.
    draw_order: Vec<DrawSource>,
    /// What the renderer's selection outline was last built from.
    outline_key: crate::selection_outline::OutlineKey,
    camera: Camera,
    /// `camera` as of the end of the last rendered frame: `camera !=
    /// previous_camera` this frame means the user is actively dragging
    /// or the view is otherwise moving, the signal `RenderSettings::
    /// fast_glass` uses (docs/RENDERING.md's "Transparency").
    previous_camera: Camera,
    lod: AdaptiveLod,
    view: ViewSettings,
    viewport_size: (u32, u32),
    pending_viewport_size: Option<egui::Vec2>,
    last_requested_px: (u32, u32),
    stable_frames: u32,

    scene: Scene,
    history: CommandHistory,
    /// View settings and camera jumps in the same undo
    /// history as scene edits (`crate::app_history`'s own doc).
    app_history: crate::app_history::AppHistory,
    /// The last committed `view` (`AppUi::commit_view_edit`'s baseline):
    /// updated on every commit, and directly on session load/restore
    /// (no edit pushed there -- loading isn't a user edit to undo).
    view_baseline: ViewSettings,
    /// The Structures panel's own "File info" target, when it differs
    /// from the current structure: `info_ui` shows this structure
    /// instead, without changing which one the ribbon edits -- unlike
    /// "current", read-only display has no correctness risk in showing a
    /// structure other than the last-loaded one.
    info_target: Option<StructureId>,
    /// The Structures panel's own click-to-select current structure: an
    /// override on top of `gpu_cache::current_structure`'s "most recently
    /// loaded" default, read by `AppUi::current()`. Unlike `info_target`
    /// this *does* change what the Represent ribbon tab edits: ribbon
    /// clicks that act on "the current structure" (`rep`/`color`/
    /// `material`) carry this structure's id explicitly so the tab and
    /// the row stay in sync (`ribbon::with_current_override`).
    current_override: Option<StructureId>,
    /// The Selections-panel row hovered this frame -- recomputed fresh on every UI pass, so
    /// releasing the hover clears it on its own. `render_scene` reads it
    /// (one frame after the panel sets it, same lag as `draw_order`) to
    /// union its rep's atoms into the viewport's selection outline.
    row_highlight: Option<(StructureId, vv_scene::RepId)>,
    /// The atom the pointer rests over in the viewport, with no button
    /// held (deliverable: hover glow) -- same one-frame lag as
    /// `row_highlight`, into the same outline.
    hover: crate::viewport::HoverPick,
    /// The viewport's pending Ctrl+click measurement chain.
    measure: crate::measure::MeasureChain,
    /// The Inspector's Contacts…/Surface area result card.
    result_card: Option<crate::ui::ResultCard>,

    egui_ctx: egui::Context,
    pub(crate) egui_state: egui_winit::State,
    egui_renderer: egui_wgpu::Renderer,
    viewport_texture: Option<egui::TextureId>,
    pub(crate) dock_state: DockState<Tab>,
    pub(crate) workspaces: Vec<(String, DockState<Tab>)>,
    /// Opened at every launch when set (`startlayout`).
    start_layout: Option<DockState<Tab>>,
    /// Theme and ribbon state, saved with the layout.
    pub(crate) prefs: crate::theme::UiPrefs,
    /// The theme egui currently has, so a change is applied once.
    applied_theme: Option<crate::theme::ThemeMode>,
    ribbon: crate::ribbon::RibbonState,
    dialog: DialogState,
    export_dialog: ExportDialogState,
    export_structure_dialog: ExportStructureDialogState,
    fetch_dialog: FetchDialogState,
    workspace_dialog: WorkspaceDialogState,
    prefs_dialog: crate::ui::PreferencesDialogState,
    palette: PaletteState,
    /// Click order of the active selection, for angles and dihedrals.
    pick_order: Vec<u32>,
    console: ConsoleState,
    /// Remaining `--exec` commands, run one per frame (`next_exec_line`).
    exec_queue: VecDeque<String>,
    /// A path-traced render running on a worker (`render`).
    render_job: Option<crate::render::RenderJob>,
    /// `--listen`: the server, and the request being run.
    live: Option<crate::live::LiveServer>,
    live_request: Option<LiveRun>,
    /// Files dropped on the window since the last frame.
    dropped: Vec<std::path::PathBuf>,
    /// A file is being dragged over the window (draw the drop hint).
    drop_hover: bool,
    /// The most recent failure, shown over the viewport for a few seconds
    /// so a failed open is visible even with the Log tab closed.
    notice: Option<Notice>,
    /// The view cube's in-flight turn, or the double-click residue
    /// zoom's.
    cube_anim: Option<crate::ui::CubeAnim>,
    /// The camera as it stood right before a pending jump, set by
    /// `AppUi::mark_camera_jump` and consumed by `commit_camera_jump`
    /// once `cube_anim` (which can span several frames) settles -- a
    /// `State` field, not a per-frame local, so it survives until then.
    camera_jump: Option<Camera>,
    /// Set when an `--exec` command fails; the rest of the script is
    /// dropped, the app quits, and `main` exits 1 so a shell or CI job
    /// can tell.
    pub(crate) exec_failed: bool,
    quit_requested: bool,
    selection_summary: crate::ui::SelectionSummary,
    /// `history.edits()` when a session was last saved or opened.
    saved_edits: u64,
    quit_dialog: bool,
    quit_confirmed: bool,
    /// The quit dialog's "Don't remind me again" checkbox, reset to unset
    /// each time the dialog opens; only written to `prefs` if the person
    /// then actually quits, not on Cancel.
    quit_dont_ask: bool,
    /// The start card stays off for the rest of this run once dismissed;
    /// starts dismissed for a script, `--listen`, a command-line file, or
    /// `prefs.hide_start_card`. `welcome` clears it.
    start_card_dismissed: bool,
    /// Studio ▸ Frame preview / `framing`: overlays the render's aspect
    /// ratio on the viewport. Not undoable or saved -- a preview aid, not
    /// a look.
    studio_frame: bool,
    sequence: crate::sequence::SequenceState,
    window_request: Option<(u32, u32)>,
    /// `uishot`: capture the whole window once this many more frames
    /// are drawn.
    ui_shot: Option<(std::path::PathBuf, u8)>,
    log: Vec<String>,

    pub(crate) frame_times: VecDeque<f32>,
    last_frame: Instant,
}

impl State {
    /// Brings up the GPU, egui and the scene. `progress` is called with a
    /// short label at each stage, for the splash's status line -- it may
    /// repaint synchronously, since nothing here needs the event loop.
    pub(crate) fn new(
        opening: crate::splash::Opening,
        options: &Options,
        progress: &mut dyn FnMut(&str),
    ) -> State {
        let t0 = Instant::now();
        let lap = |what: &str| {
            if std::env::var_os("VIZVIZ_TIMING").is_some() {
                eprintln!(
                    "startup: {what} at {:.0} ms",
                    t0.elapsed().as_secs_f64() * 1e3
                );
            }
        };
        let crate::splash::Opening {
            window,
            layout:
                layout::SavedLayout {
                    launch: dock_state,
                    workspaces,
                    prefs,
                    start: start_layout,
                },
        } = opening;
        progress("Starting the GPU\u{2026}");
        let instance = GpuContext::instance();
        lap("wgpu instance");
        let surface = instance
            .create_surface(window.clone())
            .expect("create surface");
        let ctx = GpuContext::with_instance(instance, Some(&surface), options.adapter.as_deref())
            .expect("gpu context");
        eprintln!("adapter: {}", ctx.adapter_name());
        lap("adapter + device");

        let caps = surface.get_capabilities(&ctx.adapter);
        let surface_format = caps
            .formats
            .iter()
            .copied()
            .find(|f| f.is_srgb())
            .unwrap_or(caps.formats[0]);
        let size = window.inner_size();
        let config = wgpu::SurfaceConfiguration {
            usage: wgpu::TextureUsages::RENDER_ATTACHMENT,
            format: surface_format,
            width: size.width.max(1),
            height: size.height.max(1),
            present_mode: if options.vsync {
                wgpu::PresentMode::AutoVsync
            } else {
                wgpu::PresentMode::AutoNoVsync
            },
            alpha_mode: caps.alpha_modes[0],
            color_space: wgpu::SurfaceColorSpace::Auto,
            view_formats: vec![],
            desired_maximum_frame_latency: 2,
        };
        surface.configure(&ctx.device, &config);
        lap("surface configured");

        progress("Building the renderer\u{2026}");
        let renderer = Renderer::new(ctx.clone(), DEFAULT_VIEWPORT.0, DEFAULT_VIEWPORT.1);
        lap("renderer (pipelines)");

        progress("Loading the interface\u{2026}");
        let egui_ctx = egui::Context::default();
        crate::theme::install_fonts(&egui_ctx);
        let egui_state = egui_winit::State::new(
            egui_ctx.clone(),
            egui::ViewportId::ROOT,
            &window,
            Some(window.scale_factor() as f32),
            None,
            None,
        );
        let egui_renderer = egui_wgpu::Renderer::new(
            &ctx.device,
            surface_format,
            egui_wgpu::RendererOptions::default(),
        );

        lap("egui");
        let mut scene = Scene::new();
        let mut history = CommandHistory::default();
        let mut log = Vec::new();
        let mut camera = Camera::framing(glam::Vec3::ZERO, 10.0);
        if let Some(path) = &options.file {
            match history.dispatch(&mut scene, Command::LoadStructure { path: path.clone() }) {
                Ok(()) => {
                    log.push(format!("loaded {}", path.display()));
                    let id = scene.structures().next_back().map(|(id, _)| id);
                    if let Some(framed) = id.and_then(|id| frame_structure(&scene, id)) {
                        camera = framed;
                    }
                }
                Err(e) => log.push(format!("error: {e}")),
            }
        }

        // Background work (builds, prefetched frames) wakes the window
        // when it finishes; it may be idle.
        let mut gpu_cache = GpuCache::default();
        let waker_window = window.clone();
        let waker: Waker = Arc::new(move || waker_window.request_redraw());
        gpu_cache.set_waker(waker.clone());
        let live = options.listen.and_then(|port| {
            match crate::live::LiveServer::start(port, waker.clone()) {
                Ok(server) => {
                    eprintln!("listening on 127.0.0.1:{}", server.port);
                    Some(server)
                }
                Err(e) => {
                    eprintln!("--listen: {e}");
                    None
                }
            }
        });
        // Starts dismissed for a script/client (its own reason for
        // launching) or a file already named on the command line.
        let start_card_dismissed =
            options.scripted() || options.file.is_some() || prefs.hide_start_card;
        State {
            window,
            surface,
            config,
            ctx,
            renderer,
            gpu_cache,
            timeline: TimelineState::default(),
            movie: Default::default(),
            movie_export: None,
            prefetch: Prefetch::default(),
            waker,
            timeline_offered: HashSet::new(),
            draw_order: Vec::new(),
            outline_key: Default::default(),
            previous_camera: camera.clone(),
            camera,
            lod: AdaptiveLod::default(),
            view: ViewSettings::default(),
            viewport_size: DEFAULT_VIEWPORT,
            pending_viewport_size: None,
            last_requested_px: DEFAULT_VIEWPORT,
            stable_frames: 0,
            scene,
            history,
            app_history: crate::app_history::AppHistory::default(),
            view_baseline: ViewSettings::default(),
            info_target: None,
            current_override: None,
            row_highlight: None,
            hover: crate::viewport::HoverPick::default(),
            measure: crate::measure::MeasureChain::default(),
            result_card: None,
            egui_ctx,
            egui_state,
            egui_renderer,
            viewport_texture: None,
            dock_state,
            workspaces,
            start_layout,
            prefs,
            applied_theme: None,
            ribbon: Default::default(),
            dialog: DialogState::default(),
            export_dialog: ExportDialogState::default(),
            export_structure_dialog: ExportStructureDialogState::default(),
            fetch_dialog: FetchDialogState::default(),
            workspace_dialog: WorkspaceDialogState::default(),
            prefs_dialog: crate::ui::PreferencesDialogState::default(),
            palette: PaletteState::default(),
            pick_order: Vec::new(),
            console: ConsoleState::default(),
            exec_queue: options.exec.iter().cloned().collect(),
            render_job: None,
            live,
            live_request: None,
            dropped: Vec::new(),
            drop_hover: false,
            notice: None,
            cube_anim: None,
            camera_jump: None,
            exec_failed: false,
            quit_requested: false,
            selection_summary: Default::default(),
            saved_edits: 0,
            quit_dialog: false,
            quit_confirmed: false,
            quit_dont_ask: false,
            start_card_dismissed,
            studio_frame: false,
            sequence: Default::default(),
            window_request: None,
            ui_shot: None,
            log,
            frame_times: VecDeque::with_capacity(120),
            last_frame: Instant::now(),
        }
    }

    pub(crate) fn resize_window(&mut self, width: u32, height: u32) {
        if width == 0 || height == 0 {
            return;
        }
        self.config.width = width;
        self.config.height = height;
        self.surface.configure(&self.ctx.device, &self.config);
    }

    /// Applies a debounced viewport resize: only once the requested size
    /// has been stable for a few frames, so dragging a dock divider
    /// rescales the displayed image smoothly instead of rebuilding the
    /// renderer's targets (and its whole depth-pyramid bind-group chain)
    /// every single frame.
    fn maybe_resize_viewport(&mut self) {
        let requested = self
            .pending_viewport_size
            .map(|v| {
                // Physical pixels times the supersampling factor: the
                // offscreen target is bigger than the tab and egui's
                // linear filter averages it down when drawing the image.
                let scale = self.egui_ctx.pixels_per_point() * self.view.render_scale;
                (
                    (v.x * scale).round().max(1.0) as u32,
                    (v.y * scale).round().max(1.0) as u32,
                )
            })
            .unwrap_or(self.viewport_size);

        if requested != self.last_requested_px {
            self.last_requested_px = requested;
            self.stable_frames = 0;
        } else {
            self.stable_frames += 1;
        }

        if self.stable_frames == RESIZE_DEBOUNCE_FRAMES && requested != self.viewport_size {
            self.renderer.resize(requested.0, requested.1);
            self.viewport_size = requested;
            if let Some(id) = self.viewport_texture {
                self.egui_renderer.update_egui_texture_from_wgpu_texture(
                    &self.ctx.device,
                    self.renderer.display_view_for_ui(),
                    wgpu::FilterMode::Linear,
                    id,
                );
            }
        }
    }

    fn ensure_viewport_texture(&mut self) {
        if self.viewport_texture.is_none() {
            let id = self.egui_renderer.register_native_texture(
                &self.ctx.device,
                self.renderer.display_view_for_ui(),
                wgpu::FilterMode::Linear,
            );
            self.viewport_texture = Some(id);
        }
    }

    fn render_scene(&mut self) {
        self.gpu_cache
            .sync(&self.scene, &self.ctx, &self.renderer, self.playing());
        if let Some(volume) = self.gpu_cache.update_occlusion(&self.scene, &self.ctx) {
            self.renderer.set_occlusion_volume(volume);
        }
        // Build progress: the row's spinner already shows it.
        self.log.extend(self.gpu_cache.take_infos());
        for message in self.gpu_cache.take_messages() {
            self.log.push(message.clone());
            self.notice = Some(Notice::error(message));
        }
        let (order, items, cartoons, gaussian_surfaces, skin_surfaces) =
            self.gpu_cache.draw_items(&self.scene);
        self.draw_order = order;
        crate::selection_outline::sync(
            &mut self.renderer,
            &mut self.outline_key,
            &self.scene,
            &self.draw_order,
            self.row_highlight,
            self.hover.atom,
        );

        // Glass self-occlusion culling can visibly darken a densely
        // overlapping glass scene (docs/RENDERING.md's "Transparency"),
        // so it only runs while the camera is actually moving or a
        // trajectory is playing -- a still frame always renders exact.
        // Compared against last frame's own camera, not a drag/animate
        // flag ui.rs would have to set: `Camera` already carries every
        // knob (pan, zoom, dolly, rotation, projection) a still view
        // needs to hold steady.
        let camera_moving = self.camera != self.previous_camera || self.playing();
        self.previous_camera = self.camera.clone();

        let settings = RenderSettings {
            // Per structure via `DrawItem`; this field is unused by render_all.
            representation: Default::default(),
            quad_px_threshold: self.lod.threshold_px,
            occlusion_culling: self.view.occlusion_culling,
            fast_glass: camera_moving,
            lighting: self.view.lighting(),
            material: Default::default(),
            background: self.view.background,
            background_top: self.view.background_top,
            fxaa: self.view.fxaa,
            outline: self.view.outline,
            outline_width: 1,
            selection_color: RenderSettings::default().selection_color,
            selection_width: self.selection_width(),
            ao: self.view.ao,
            depth_cue: self.view.depth_cue,
            shadows: self.view.shadows,
            dof: self.view.dof,
            clip: self.view.clip.map(|c| c.world(&self.camera)),
            scene_bounds: None,
        };

        let mut encoder = self.ctx.device.create_command_encoder(&Default::default());
        if items.is_empty()
            && cartoons.is_empty()
            && gaussian_surfaces.is_empty()
            && skin_surfaces.is_empty()
        {
            // Nothing loaded: clear the target so the viewport shows the
            // background instead of a stale or uninitialized frame.
            self.renderer.clear(&mut encoder, settings.background);
        } else {
            self.renderer.render_all(
                &mut encoder,
                &self.camera,
                &items,
                &cartoons,
                &gaussian_surfaces,
                &skin_surfaces,
                &settings,
            );
        }
        self.ctx.queue.submit([encoder.finish()]);

        if self.view.adaptive {
            if let Some(times) = self.renderer.last_times() {
                self.lod.update(times.total_ms());
            }
        }
    }

    /// Two logical pixels, in the viewport target's supersampled pixels.
    fn selection_width(&self) -> u32 {
        let scale = self.egui_ctx.pixels_per_point() * self.view.render_scale;
        (2.0 * scale).round().max(1.0) as u32
    }

    /// Renders one frame off to the side, on a temporary, throwaway
    /// `Renderer`, and writes it to a PNG. Never touches the live
    /// `self.renderer`: reusing it would mean either fighting over its GPU
    /// timer's single in-flight slot (a second `render`+`after_submit` this
    /// frame steals the slot `render_scene()` already armed, so the live
    /// frame's timings never get read back — found the hard way) or, for
    /// SSAA, resizing it and re-registering its egui texture mid-frame. A
    /// one-time shader-compile hitch on the temp renderer is a fine trade
    /// for sidestepping both.
    fn export_screenshot(&mut self, request: ExportRequest) {
        let (vw, vh) = self.renderer.size();
        let Some(rgba) = self.render_offscreen(vw, vh, request.ssaa) else {
            return;
        };
        let outcome = write_screenshot(&request.path, &rgba, vw, vh)
            .map(|()| format!("Saved {} ({vw}x{vh})", request.path.display()));
        self.export_outcome(outcome);
    }

    /// The scene through `self.camera` at `vw x vh` as RGBA8, labels and
    /// captions baked in. `None` (logged) when nothing is drawn.
    fn render_offscreen(&mut self, vw: u32, vh: u32, ssaa: bool) -> Option<Vec<u8>> {
        // `render_scene()` already ran earlier this frame, so the cache is
        // current; the same draw list the viewport just used (tubes and
        // ligand sticks included) is what gets exported.
        let (_, live_items, live_cartoons, live_gaussian_surfaces, live_skin_surfaces) =
            self.gpu_cache.draw_items(&self.scene);
        if live_items.is_empty()
            && live_cartoons.is_empty()
            && live_gaussian_surfaces.is_empty()
            && live_skin_surfaces.is_empty()
        {
            self.log.push("export: no structure loaded".into());
            return None;
        }
        let scale = if ssaa { 2 } else { 1 };
        let settings = RenderSettings {
            representation: Default::default(),
            // Force full detail regardless of the live adaptive threshold:
            // this is a one-shot export, not the interactive frame loop
            // AdaptiveLod exists to protect.
            quad_px_threshold: self.lod.min_px,
            occlusion_culling: self.view.occlusion_culling,
            // Exact, not the interactive approximation: a one-shot export
            // is exactly the still view `fast_glass` is off for.
            fast_glass: false,
            lighting: self.view.lighting(),
            material: Default::default(),
            background: self.view.background,
            background_top: self.view.background_top,
            // Irrelevant here: `read_color` below reads the outlined
            // target, before FXAA. `false` just picks the cheaper of the
            // two passes for output nothing reads.
            fxaa: false,
            outline: self.view.outline,
            // Drawn at export resolution, so scale the line width with it
            // or the 2x SSAA downsample halves it.
            outline_width: scale,
            // Never drawn: this renderer is given no selection.
            selection_color: RenderSettings::default().selection_color,
            selection_width: scale,
            ao: self.view.ao,
            depth_cue: self.view.depth_cue,
            shadows: self.view.shadows,
            dof: self.view.dof,
            clip: self.view.clip.map(|c| c.world(&self.camera)),
            scene_bounds: None,
        };
        let (rw, rh) = (vw * scale, vh * scale);
        let mut temp = Renderer::new(self.ctx.clone(), rw, rh);
        temp.set_occlusion_volume(self.gpu_cache.occlusion());
        // The live items' bind groups serve this renderer too (same
        // layouts); each rep's draw state is rewritten per frame anyway.
        let mut encoder = self.ctx.device.create_command_encoder(&Default::default());
        temp.render_all(
            &mut encoder,
            &self.camera,
            &live_items,
            &live_cartoons,
            &live_gaussian_surfaces,
            &live_skin_surfaces,
            &settings,
        );
        self.ctx.queue.submit([encoder.finish()]);
        temp.after_submit();
        let rendered = temp.read_color();
        let mut rgba = if ssaa {
            vv_render::downsample_2x_srgb(&rendered, rw, rh)
        } else {
            rendered
        };

        // Baked in after SSAA downsampling (sharp, not softened by the
        // average) at the same `vw x vh` the file will actually be, using
        // the identical current-frame lookup the live overlay uses so a
        // label doesn't drift to a stale frame's position in the export.
        self.bake_labels_and_captions(&mut rgba, vw, vh);
        Some(rgba)
    }

    /// Starts a path-traced render of the current view on a worker
    /// (`crate::render`); one at a time.
    fn start_render(&mut self, request: crate::render::RenderRequest) {
        if self.render_job.is_some() {
            self.render_outcome(Err("render: one is already running".into()));
            return;
        }
        let input = crate::render::trace_scene(&self.scene, &self.gpu_cache, request.quality);
        if !input.skipped.is_empty() {
            let line = format!("render leaves out: {}", input.skipped.join(", "));
            self.log.push(line.clone());
            self.notice = Some(Notice::info(line));
        }
        if input.is_empty() {
            self.render_outcome(Err("render: nothing it can draw is shown".into()));
            return;
        }
        let linear = |c: wgpu::Color| [c.r as f32, c.g as f32, c.b as f32];
        let settings = vv_render::path_trace::TraceSettings {
            width: request.width,
            height: request.height,
            samples: request.samples,
            lighting: self.view.lighting(),
            background: linear(self.view.background),
            background_top: linear(self.view.background_top.unwrap_or(self.view.background)),
            transparent: request.transparent,
            clip: self.view.clip.map(|c| c.world(&self.camera)),
            unshadowed: false,
            light_spread: request.quality.light_spread(),
            ao: vv_render::path_trace::RENDER_AO,
            direct: vv_render::path_trace::RENDER_DIRECT,
        };
        self.render_job = Some(crate::render::RenderJob::start(
            self.ctx.clone(),
            input,
            self.camera.clone(),
            settings,
            request.path,
            self.waker.clone(),
        ));
    }

    /// Shows a running render's progress; writes the image when it is done.
    fn poll_render(&mut self) {
        let Some(job) = &self.render_job else {
            return;
        };
        let Some(done) = job.poll() else {
            return;
        };
        let job = self.render_job.take().expect("checked above");
        let result = done.and_then(|rgba| {
            write_screenshot(&job.path, &rgba, job.width, job.height).map_err(|e| e.to_string())
        });
        self.render_outcome(
            result
                .map(|()| {
                    format!(
                        "rendered {} ({}x{}, {} samples, {:.1} s)",
                        job.path.display(),
                        job.width,
                        job.height,
                        job.samples,
                        job.started.elapsed().as_secs_f32()
                    )
                })
                .map_err(|e| format!("render failed: {e}")),
        );
    }

    /// Reports how a render went, in the log, the notice bar and on stderr
    /// (as `--exec` output is); a failure stops a running script, as a
    /// failed command does.
    fn render_outcome(&mut self, outcome: Result<String, String>) {
        let (line, failed) = match outcome {
            Ok(line) => (line, false),
            Err(line) => (line, true),
        };
        eprintln!("{line}");
        self.log.push(line.clone());
        if failed {
            self.notice = Some(Notice::error(line));
            if self.live_request.is_some() {
                self.fail_live_request();
            } else if !self.exec_queue.is_empty() {
                self.exec_failed = true;
                self.exec_queue.clear();
                self.quit_requested = true;
            }
        } else {
            self.notice = Some(Notice::info(line));
        }
    }

    /// Vector export, phase 1: every atom (of any structure, any
    /// representation) as a painter's-algorithm-sorted SVG `<circle>`, at
    /// its true van der Waals screen radius, with whole-primitive
    /// back-to-front ordering and no circle-circle clipping. Reuses `project_to_pixel` (shared with the
    /// live overlay and raster export) and `gpu_cache::colors_of` (shared
    /// with the GPU path), so a vector figure agrees with what's on
    /// screen. Deliberately not attempted here, stated as phase-2 cuts:
    /// bonds/cylinders (ball-and-stick, tube) and cartoon ribbon
    /// triangles -- every atom is drawn as if in spacefill regardless of
    /// its actual live representation, and there is no depth-buffer
    /// visibility culling (a fully atom hidden behind a dense core is
    /// still emitted; painter's-algorithm ordering still composites it
    /// correctly, at the cost of a larger file than strictly necessary).
    fn export_svg(&mut self, request: &ExportRequest) {
        let (vw, vh) = self.renderer.size();
        let (document, atom_count) =
            build_svg_document(&self.scene, &self.camera, vw, vh, self.view.background);
        let outcome = svg::save(&request.path, &document).map(|()| {
            format!(
                "Saved {} ({vw}x{vh}, {atom_count} atoms)",
                request.path.display()
            )
        });
        self.export_outcome(outcome);
    }

    /// A finished screenshot in the log and the notice bar, naming the
    /// full path it went to.
    fn export_outcome(&mut self, outcome: std::io::Result<String>) {
        let notice = match outcome {
            Ok(line) => Notice::info(line),
            Err(e) => Notice::error(format!("export failed: {e}")),
        };
        self.log.push(notice.text.clone());
        self.notice = Some(notice);
    }

    /// CPU-rasterizes labels and captions directly into an export buffer
    /// (`crate::text`), since `export_screenshot` never runs an egui paint
    /// pass the way the live viewport's `labels_and_captions_ui` does, so
    /// export needs its own drawing path rather than reusing egui's
    /// renderer.
    fn bake_labels_and_captions(&self, rgba: &mut [u8], width: u32, height: u32) {
        let view_proj = self.camera.proj(width as f32 / height.max(1) as f32) * self.camera.view();
        for draw in label_draws(&self.scene) {
            let Some((x, y, _depth)) =
                project_to_pixel(view_proj, width as f32, height as f32, draw.world)
            else {
                continue;
            };
            crate::text::draw_dot(rgba, width, height, x, y, LABEL_DOT_RADIUS, LABEL_DOT_COLOR);
            crate::text::draw_text(
                rgba,
                width,
                height,
                x + 6.0,
                y - 14.0,
                &draw.text,
                LABEL_FONT_PX,
                LABEL_TEXT_COLOR,
            );
        }
        for m in measurement_draws(&self.scene) {
            let project = |w| {
                project_to_pixel(view_proj, width as f32, height as f32, w).map(|(x, y, _)| (x, y))
            };
            let Some(points) = m
                .path
                .iter()
                .map(|&w| project(w))
                .collect::<Option<Vec<_>>>()
            else {
                continue;
            };
            for pair in points.windows(2) {
                crate::text::draw_dashed_line(
                    rgba,
                    width,
                    height,
                    pair[0],
                    pair[1],
                    MEASURE_DASH,
                    MEASURE_LINE_COLOR,
                );
            }
            if let Some((x, y)) = project(m.anchor) {
                crate::text::draw_text(
                    rgba,
                    width,
                    height,
                    x + 4.0,
                    y - 12.0,
                    &m.text,
                    LABEL_FONT_PX,
                    LABEL_TEXT_COLOR,
                );
            }
        }
        for caption in self.scene.captions() {
            let x = caption.x.clamp(0.0, 1.0) * width as f32;
            let y = caption.y.clamp(0.0, 1.0) * height as f32;
            crate::text::draw_text(
                rgba,
                width,
                height,
                x,
                y,
                &caption.text,
                CAPTION_FONT_PX,
                CAPTION_TEXT_COLOR,
            );
        }
    }

    /// Steps the timeline's frame when playback is on and the frame period
    /// has elapsed. Runs before `render_scene` so the new coordinates are
    /// uploaded in the same frame.
    fn advance_playback(&mut self) {
        if !self.timeline.playing {
            return;
        }
        let Some((id, count)) = timeline_target(&self.scene) else {
            self.timeline.playing = false;
            return;
        };
        let now = Instant::now();
        let period = Duration::from_secs_f32(1.0 / self.timeline.fps.max(0.1));
        if let Some(last) = self.timeline.last_step {
            if now.duration_since(last) < period {
                return;
            }
        }
        let Some(loaded) = self.scene.structure(id) else {
            return;
        };
        let frame = loaded.frame;
        let next = if frame + 1 < count {
            frame + 1
        } else if self.timeline.looping {
            0
        } else {
            self.timeline.playing = false;
            frame
        };
        // A streamed trajectory: read half a second ahead, and step only
        // onto a frame already in memory (retried each redraw), so slow
        // reads slow playback down instead of stalling the window.
        let ahead = (self.timeline.fps * 0.5).ceil().max(4.0) as usize;
        self.prefetch.ahead(
            id,
            &loaded.structure,
            next,
            ahead,
            self.timeline.looping,
            Some(self.waker.clone()),
        );
        if !loaded.structure.has_frame(next) {
            return;
        }
        self.timeline.last_step = Some(now);
        self.scene.set_frame_live(id, next);
    }

    /// Steps the camera's in-flight turn (a view-cube face, or the
    /// double-click residue zoom), if any, and drops it once it reaches
    /// its target.
    fn advance_cube_anim(&mut self) {
        let Some(anim) = &self.cube_anim else {
            return;
        };
        let (orientation, target, distance, finished) = anim.step();
        self.camera.orientation = orientation;
        self.camera.target = target;
        self.camera.distance = distance;
        if finished {
            self.cube_anim = None;
        }
    }

    /// Opens the Timeline tab under the viewport the first time a
    /// structure with more than one frame appears (the Timeline tab
    /// auto-opens on trajectory load).
    fn maybe_open_timeline(&mut self) {
        let newcomers: Vec<StructureId> = self
            .scene
            .structures()
            .filter(|(id, s)| s.structure.frame_count() > 1 && !self.timeline_offered.contains(id))
            .map(|(id, _)| id)
            .collect();
        if newcomers.is_empty() {
            return;
        }
        layout::open_near_viewport(&mut self.dock_state, Tab::Timeline);
        self.timeline_offered.extend(newcomers);
    }

    pub(crate) fn drop_file(&mut self, path: std::path::PathBuf) {
        self.drop_hover = false;
        self.dropped.push(path);
        self.window.request_redraw();
    }

    pub(crate) fn set_drop_hover(&mut self, hovering: bool) {
        self.drop_hover = hovering;
        self.window.request_redraw();
    }

    /// What `layout.json` keeps of this session.
    pub(crate) fn saved_layout(&self) -> layout::SavedLayout {
        layout::SavedLayout {
            launch: self.dock_state.clone(),
            workspaces: self.workspaces.clone(),
            prefs: self.prefs.clone(),
            start: self.start_layout.clone(),
        }
    }

    /// Whether the app has asked to exit (the `quit` command). Reading it
    /// clears it.
    pub(crate) fn take_quit_request(&mut self) -> bool {
        std::mem::take(&mut self.quit_requested)
    }

    /// Whether to exit now. A person with changes no saved session holds
    /// is asked first (the answer comes back through `quit`); a script or
    /// live client never is.
    pub(crate) fn may_quit(&mut self, scripted: bool) -> bool {
        let unsaved =
            self.scene.structures().next().is_some() && self.saved_edits != self.history.edits();
        if scripted || !unsaved || self.quit_confirmed || self.prefs.dont_confirm_quit {
            return true;
        }
        self.quit_dialog = true;
        self.window.request_redraw();
        false
    }

    /// The next `--exec` command to run this frame, if any. One per frame,
    /// and only once the viewport has settled at its docked size, so a
    /// `load` is uploaded and drawn before a following `screenshot` runs,
    /// and the screenshot is the size the user will see, not the startup
    /// placeholder.
    fn next_exec_line(&mut self) -> Option<String> {
        // Also wait out background cartoon/surface builds, so `rep cartoon;
        // screenshot x.png` captures the cartoon, not the atoms drawn
        // meanwhile.
        if self.exec_queue.is_empty() || self.busy() {
            return None;
        }
        self.exec_queue.pop_front()
    }

    /// Whether background work a command started is still running: the
    /// next scripted line waits for it.
    fn busy(&self) -> bool {
        self.stable_frames <= RESIZE_DEBOUNCE_FRAMES
            || self.gpu_cache.building()
            || self.movie_export.is_some()
            || self.render_job.is_some()
            || self.ui_shot.is_some()
    }

    /// The next line of the `--listen` request being run, taking the next
    /// request when none is; one line per frame, as `--exec` runs.
    fn next_live_line(&mut self) -> Option<String> {
        if self.live_request.is_none() {
            let request = self.live.as_ref()?.poll()?;
            self.live_request = Some(LiveRun {
                lines: request.lines.into(),
                log_start: self.log.len(),
                failed: false,
                reply: request.reply,
            });
        }
        if self.busy() {
            return None;
        }
        self.live_request.as_mut()?.lines.pop_front()
    }

    /// Ends the running `--listen` request as failed: its remaining lines
    /// are dropped once what is running finishes.
    fn fail_live_request(&mut self) {
        if let Some(run) = &mut self.live_request {
            run.failed = true;
            run.lines.clear();
        }
    }

    /// Answers the running `--listen` request once its lines have run and
    /// nothing they started is still going.
    fn finish_live_request(&mut self) {
        let done = self
            .live_request
            .as_ref()
            .is_some_and(|run| run.lines.is_empty());
        if !done || self.busy() {
            return;
        }
        let run = self.live_request.take().expect("checked above");
        let _ = run.reply.send(crate::live::Reply {
            ok: !run.failed,
            output: self.log[run.log_start..].to_vec(),
        });
    }

    fn apply_layout_request(&mut self, request: LayoutRequest) {
        match request {
            LayoutRequest::Switch(name) => {
                let built_in = layout::built_in_workspaces()
                    .into_iter()
                    .find(|(n, _)| n.eq_ignore_ascii_case(&name))
                    .map(|(_, s)| s);
                let saved = self
                    .workspaces
                    .iter()
                    .find(|(n, _)| n.eq_ignore_ascii_case(&name))
                    .map(|(_, s)| s.clone());
                if let Some(state) = built_in.or(saved) {
                    self.dock_state = state;
                }
            }
            LayoutRequest::AddPanel(path, tab) => {
                if self.dock_state.find_tab(&tab).is_some() {
                    return;
                }
                self.dock_state.set_focused_node_and_surface(path);
                self.dock_state.push_to_focused_leaf(tab);
                if !layout::viewport_in_place(&self.dock_state) {
                    self.dock_state
                        .remove_tab(self.dock_state.find_tab(&tab).expect("just added"));
                    layout::open_near_viewport(&mut self.dock_state, tab);
                }
            }
            LayoutRequest::OpenPanel(tab) => match self.dock_state.find_tab(&tab) {
                Some(found) => {
                    let _ = self.dock_state.set_active_tab(found);
                }
                None => layout::open_near_viewport(&mut self.dock_state, tab),
            },
            LayoutRequest::ClosePanel(tab) => {
                if let Some(found) = self.dock_state.find_tab(&tab) {
                    self.dock_state.remove_tab(found);
                }
            }
            LayoutRequest::SetStart => self.start_layout = Some(self.dock_state.clone()),
            LayoutRequest::ResetStart => self.start_layout = None,
            LayoutRequest::Save(name) => {
                let snapshot = self.dock_state.clone();
                match self
                    .workspaces
                    .iter_mut()
                    .find(|(n, _)| n.eq_ignore_ascii_case(&name))
                {
                    Some(existing) => existing.1 = snapshot,
                    None => self.workspaces.push((name, snapshot)),
                }
            }
            LayoutRequest::GrowSequence(px) => layout::grow_sequence(&mut self.dock_state, px),
            LayoutRequest::LoadState(state) => {
                if layout::viewport_in_place(&state) {
                    self.dock_state = *state;
                }
            }
        }
    }

    pub(crate) fn redraw(&mut self) {
        self.advance_playback();
        self.step_movie();
        self.advance_cube_anim();
        self.maybe_open_timeline();
        self.render_scene();
        self.ensure_viewport_texture();

        let raw_input = self.egui_state.take_egui_input(&self.window);
        let readout = self.readout();
        // Taken inside the closure: egui reruns it when a pass is
        // discarded (a new grid sizing itself), and these must run once.
        let mut exec_line = self.next_exec_line();
        let mut live_line = if exec_line.is_none() {
            self.next_live_line()
        } else {
            None
        };
        let mut live_failed = false;
        let mut dropped = std::mem::take(&mut self.dropped);
        let mut exec_failed = false;
        let mut requested_viewport_size = None;
        let mut export_request: Option<ExportRequest> = None;
        let mut render_request: Option<crate::render::RenderRequest> = None;
        let mut layout_request: Option<LayoutRequest> = None;
        // A plain clone, not a borrow: `self.dock_state` is borrowed
        // mutably below for `DockArea`, so `savesession` reads this
        // frame-start copy instead (see `AppUi::dock_layout`'s doc).
        let dock_layout_snapshot = self.dock_state.clone();
        let egui_ctx = self.egui_ctx.clone();
        if self.applied_theme != Some(self.prefs.theme) {
            crate::theme::apply(&egui_ctx, self.prefs.theme);
            let was = self.applied_theme.unwrap_or(crate::theme::ThemeMode::Dark);
            self.view.follow_theme(was, self.prefs.theme);
            self.view_baseline.follow_theme(was, self.prefs.theme);
            self.applied_theme = Some(self.prefs.theme);
        }
        // Set again below by a hovered or current selection row, if its
        // panel is open.
        self.row_highlight = None;
        let render_progress = self.render_job.as_ref().map(|job| job.progress());
        let mut full_output = egui_ctx.run_ui(raw_input, |ui| {
            let closed_panels = layout::closed_panels(&self.dock_state);
            let mut app_ui = AppUi {
                scene: &mut self.scene,
                history: &mut self.history,
                app_history: &mut self.app_history,
                camera_jump: &mut self.camera_jump,
                view_baseline: &mut self.view_baseline,
                info_target: &mut self.info_target,
                current_override: &mut self.current_override,
                row_highlight: &mut self.row_highlight,
                hover: &mut self.hover,
                measure: &mut self.measure,
                result_card: &mut self.result_card,
                camera: &mut self.camera,
                view: &mut self.view,
                renderer: &self.renderer,
                gpu_cache: &self.gpu_cache,
                draw_order: &self.draw_order,
                timeline: &mut self.timeline,
                movie: &mut self.movie,
                lod: &mut self.lod,
                viewport_texture: self.viewport_texture,
                requested_viewport_size: &mut requested_viewport_size,
                dialog: &mut self.dialog,
                export_dialog: &mut self.export_dialog,
                export_structure_dialog: &mut self.export_structure_dialog,
                fetch_dialog: &mut self.fetch_dialog,
                export_request: &mut export_request,
                render_request: &mut render_request,
                workspace_dialog: &mut self.workspace_dialog,
                prefs_dialog: &mut self.prefs_dialog,
                workspaces: &mut self.workspaces,
                layout_request: &mut layout_request,
                closed_panels,
                dock_layout: &dock_layout_snapshot,
                palette: &mut self.palette,
                pick_order: &mut self.pick_order,
                drop_hover: self.drop_hover,
                notice: &mut self.notice,
                cube_anim: &mut self.cube_anim,
                render_progress,
                console: &mut self.console,
                quit_requested: &mut self.quit_requested,
                selection_summary: &mut self.selection_summary,
                saved_edits: &mut self.saved_edits,
                quit_dialog: &mut self.quit_dialog,
                quit_confirmed: &mut self.quit_confirmed,
                quit_dont_ask: &mut self.quit_dont_ask,
                start_card_dismissed: &mut self.start_card_dismissed,
                studio_frame: &mut self.studio_frame,
                sequence: &mut self.sequence,
                window_request: &mut self.window_request,
                ui_shot: &mut self.ui_shot,
                prefs: &mut self.prefs,
                ribbon: &mut self.ribbon,
                log: &mut self.log,
                readout: readout.clone(),
            };
            for path in dropped.drain(..) {
                let verb = if path.extension().is_some_and(|e| {
                    e.eq_ignore_ascii_case("vviz") || e.eq_ignore_ascii_case("json")
                }) {
                    "loadsession"
                } else {
                    "open"
                };
                app_ui.run_command_logged(&format!("{verb} {}", path.display()));
            }
            if let Some(line) = exec_line.take() {
                // Scripted runs are often unattended: mirror the console
                // to stderr so the terminal shows what happened, and stop
                // at the first failure (see `exec_failed`) rather than
                // carrying on and exiting 0 with the figure never written.
                let before = app_ui.log.len();
                let ok = app_ui.run_command_logged(&line);
                for entry in &app_ui.log[before..] {
                    eprintln!("{entry}");
                }
                if !ok {
                    exec_failed = true;
                    *app_ui.quit_requested = true;
                }
            }
            if let Some(line) = live_line.take() {
                // A client's failing line ends its request, not the app.
                live_failed |= !app_ui.run_command_logged(&line);
            }
            app_ui.handle_keys(ui.ctx());
            egui::Panel::top("ribbon").show(ui, |ui| app_ui.ribbon_ui(ui));
            egui::Panel::left("tools")
                .resizable(false)
                .exact_size(crate::theme::CONTROL_HEIGHT + 2.0 * crate::theme::space::TIGHT)
                .show(ui, |ui| app_ui.tool_strip_ui(ui));
            egui::CentralPanel::default()
                .frame(egui::Frame::NONE)
                .show(ui, |ui| {
                    // Each tab closes from its own x; no second x per panel.
                    // The viewport keeps no header of its own.
                    layout::hide_viewport_tab_bar(&mut self.dock_state);
                    DockArea::new(&mut self.dock_state)
                        .style(crate::theme::dock_style(ui.style(), app_ui.prefs.theme))
                        .show_leaf_close_all_buttons(false)
                        .show_add_buttons(true)
                        .show_add_popup(true)
                        .hidable_tab_bars(true)
                        .show_inside(ui, &mut app_ui);
                    if !layout::viewport_in_place(&self.dock_state) {
                        self.dock_state = dock_layout_snapshot.clone();
                    }
                });
            app_ui.open_dialog(ui.ctx());
            app_ui.export_dialog_ui(ui.ctx());
            app_ui.export_structure_dialog_ui(ui.ctx());
            app_ui.fetch_dialog_ui(ui.ctx());
            app_ui.workspace_dialog_ui(ui.ctx(), &self.dock_state);
            app_ui.preferences_dialog_ui(ui.ctx());
            app_ui.quit_dialog_ui(ui.ctx());
            app_ui.palette_ui(ui.ctx());
            // Last, so every ribbon/panel/popover mutation this
            // frame is included in what gets committed.
            app_ui.commit_view_edit(ui.ctx());
            app_ui.commit_camera_jump();
        });
        self.pending_viewport_size = requested_viewport_size;
        self.maybe_resize_viewport();
        if exec_failed {
            self.exec_failed = true;
            self.exec_queue.clear();
        }
        if live_failed {
            self.fail_live_request();
        }
        if let Some(request) = layout_request {
            self.apply_layout_request(request);
        }
        if let Some(mut request) = render_request {
            request.path = export_target(&request.path);
            self.start_render(request);
        }
        self.poll_render();
        self.finish_live_request();
        if let Some(mut request) = export_request {
            request.path = export_target(&request.path);
            let is_svg = request
                .path
                .extension()
                .and_then(|e| e.to_str())
                .is_some_and(|e| e.eq_ignore_ascii_case("svg"));
            if is_svg {
                self.export_svg(&request);
            } else {
                self.export_screenshot(request);
            }
        }
        self.take_movie_export_request();
        self.capture_movie_frame();

        self.egui_state
            .handle_platform_output(&self.window, full_output.platform_output.clone());
        let paint_jobs = self
            .egui_ctx
            .tessellate(full_output.shapes.clone(), full_output.pixels_per_point);

        for (id, deltas) in &full_output.textures_delta.set {
            for delta in deltas {
                self.egui_renderer
                    .update_texture(&self.ctx.device, &self.ctx.queue, *id, delta);
            }
        }
        for id in &full_output.textures_delta.free {
            self.egui_renderer.free_texture(id);
        }
        // `TexturesDelta` panics (debug builds) or leaks (release) if
        // dropped non-empty: every entry must be handled, then cleared,
        // including on the early returns below.
        full_output.textures_delta.clear();

        use wgpu::CurrentSurfaceTexture as Current;
        let frame = match self.surface.get_current_texture() {
            Current::Success(frame) | Current::Suboptimal(frame) => frame,
            Current::Outdated | Current::Lost => {
                self.resize_window(self.config.width, self.config.height);
                return;
            }
            Current::Timeout | Current::Occluded => return,
            other => {
                eprintln!("surface: {other:?}");
                return;
            }
        };
        let view = frame.texture.create_view(&Default::default());

        let screen_descriptor = egui_wgpu::ScreenDescriptor {
            size_in_pixels: [self.config.width, self.config.height],
            pixels_per_point: full_output.pixels_per_point,
        };
        let mut encoder = self.ctx.device.create_command_encoder(&Default::default());
        let command_buffers = self.egui_renderer.update_buffers(
            &self.ctx.device,
            &self.ctx.queue,
            &mut encoder,
            &paint_jobs,
            &screen_descriptor,
        );
        {
            let pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("egui"),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view: &view,
                    resolve_target: None,
                    ops: wgpu::Operations {
                        load: wgpu::LoadOp::Clear(wgpu::Color::BLACK),
                        store: wgpu::StoreOp::Store,
                    },
                    depth_slice: None,
                })],
                ..Default::default()
            });
            let mut pass = pass.forget_lifetime();
            self.egui_renderer
                .render(&mut pass, &paint_jobs, &screen_descriptor);
        }
        self.ctx
            .queue
            .submit(command_buffers.into_iter().chain([encoder.finish()]));
        self.renderer.after_submit();
        match self.ui_shot.take() {
            Some((path, 0)) => self.capture_window(&path, &paint_jobs, &screen_descriptor),
            Some((path, frames)) => {
                self.ui_shot = Some((path, frames - 1));
                self.window.request_redraw();
            }
            None => {}
        }
        if let Some((w, h)) = self.window_request.take() {
            self.window.set_maximized(false);
            let _ = self
                .window
                .request_inner_size(winit::dpi::LogicalSize::new(w, h));
        }
        self.ctx.queue.present(frame);

        let now = Instant::now();
        let dt = now.duration_since(self.last_frame).as_secs_f32() * 1e3;
        self.last_frame = now;
        if self.frame_times.len() == 120 {
            self.frame_times.pop_front();
        }
        self.frame_times.push_back(dt);

        // Self-scheduling: the event loop (`ControlFlow::Wait`) only ticks
        // again when something asks it to. Keep asking while there's
        // exec-queue work waiting (including the resize-debounce window
        // it's gated behind, see `next_exec_line`), a trajectory playing,
        // or egui itself wants to redraw soon (a focused cursor blinking,
        // a hover fade); otherwise a static viewport truly goes idle.
        let egui_wants_soon = full_output
            .viewport_output
            .get(&egui::ViewportId::ROOT)
            .is_some_and(|v| v.repaint_delay < Duration::from_secs(1));
        if !self.exec_queue.is_empty()
            || self.stable_frames <= RESIZE_DEBOUNCE_FRAMES
            || self.playing()
            || self.movie.active()
            || self.cube_anim.is_some()
            || self.render_job.is_some()
            || self.live_request.is_some()
            || egui_wants_soon
        {
            self.window.request_redraw();
        }
    }

    /// Draws this frame's UI again into an offscreen texture and writes it
    /// as a PNG: the whole window, for documentation and UI review.
    fn capture_window(
        &mut self,
        path: &std::path::Path,
        paint_jobs: &[egui::ClippedPrimitive],
        screen: &egui_wgpu::ScreenDescriptor,
    ) {
        let (width, height) = (self.config.width, self.config.height);
        let texture = self.ctx.device.create_texture(&wgpu::TextureDescriptor {
            label: Some("window capture"),
            size: wgpu::Extent3d {
                width,
                height,
                depth_or_array_layers: 1,
            },
            mip_level_count: 1,
            sample_count: 1,
            dimension: wgpu::TextureDimension::D2,
            format: self.config.format,
            usage: wgpu::TextureUsages::RENDER_ATTACHMENT | wgpu::TextureUsages::COPY_SRC,
            view_formats: &[],
        });
        let view = texture.create_view(&Default::default());
        let mut encoder = self.ctx.device.create_command_encoder(&Default::default());
        {
            let pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("window capture"),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view: &view,
                    resolve_target: None,
                    ops: wgpu::Operations {
                        load: wgpu::LoadOp::Clear(wgpu::Color::BLACK),
                        store: wgpu::StoreOp::Store,
                    },
                    depth_slice: None,
                })],
                ..Default::default()
            });
            let mut pass = pass.forget_lifetime();
            self.egui_renderer.render(&mut pass, paint_jobs, screen);
        }
        self.ctx.queue.submit([encoder.finish()]);
        let mut rgba = vv_render::renderer::read_texture_rgba8(&self.ctx, &texture);
        if matches!(
            self.config.format,
            wgpu::TextureFormat::Bgra8Unorm | wgpu::TextureFormat::Bgra8UnormSrgb
        ) {
            for px in rgba.chunks_mut(4) {
                px.swap(0, 2);
            }
        }
        match write_png(path, &rgba, width, height) {
            Ok(()) => self.log.push(format!(
                "captured the window to {} ({width}x{height})",
                path.display()
            )),
            Err(e) => self.log.push(format!("window capture failed: {e}")),
        }
    }

    fn readout(&self) -> Readout {
        let avg = if self.frame_times.is_empty() {
            0.0
        } else {
            self.frame_times.iter().sum::<f32>() / self.frame_times.len() as f32
        };
        let times = self.renderer.last_times().unwrap_or_default();
        Readout {
            fps: 1e3 / avg.max(1e-3),
            frame_ms: avg,
            cull_ms: times.cull_ms,
            draw_ms: times.draw_ms,
            hiz_ms: times.hiz_ms,
            quad_px_threshold: self.lod.threshold_px,
        }
    }
}

fn write_png(path: &std::path::Path, rgba: &[u8], width: u32, height: u32) -> std::io::Result<()> {
    let file = std::fs::File::create(path)?;
    let mut encoder = png::Encoder::new(std::io::BufWriter::new(file), width, height);
    encoder.set_color(png::ColorType::Rgba);
    encoder.set_depth(png::BitDepth::Eight);
    let mut writer = encoder
        .write_header()
        .map_err(|e| std::io::Error::other(e.to_string()))?;
    writer
        .write_image_data(rgba)
        .map_err(|e| std::io::Error::other(e.to_string()))
}

/// JPEG has no alpha channel; `rgba` is written as opaque RGB. Quality 90
/// matches what most figure-export tools default to (visually lossless
/// for line art/molecular renders, well under half the file size of 100).
fn write_jpeg(path: &std::path::Path, rgba: &[u8], width: u32, height: u32) -> std::io::Result<()> {
    let file = std::fs::File::create(path)?;
    let encoder = jpeg_encoder::Encoder::new(std::io::BufWriter::new(file), 90);
    encoder
        .encode(
            rgba,
            width as u16,
            height as u16,
            jpeg_encoder::ColorType::Rgba,
        )
        .map_err(std::io::Error::other)
}

/// Chooses PNG or JPEG by `path`'s extension (`.jpg`/`.jpeg`, case
/// insensitive; anything else is PNG) -- the same "infer format from the
/// file name" convention `vv_io::Format::from_path` already uses for
/// loading, so `screenshot out.jpg` just works like `screenshot out.png`.
/// `path` made absolute (against the current directory, as a script
/// expects) with its folder created, so a save always says where it went.
fn export_target(path: &std::path::Path) -> std::path::PathBuf {
    let path = std::path::absolute(path).unwrap_or_else(|_| path.to_owned());
    if let Some(dir) = path.parent() {
        let _ = std::fs::create_dir_all(dir);
    }
    path
}

fn write_screenshot(
    path: &std::path::Path,
    rgba: &[u8],
    width: u32,
    height: u32,
) -> std::io::Result<()> {
    let is_jpeg = path
        .extension()
        .and_then(|e| e.to_str())
        .is_some_and(|e| e.eq_ignore_ascii_case("jpg") || e.eq_ignore_ascii_case("jpeg"));
    if is_jpeg {
        write_jpeg(path, rgba, width, height)
    } else {
        write_png(path, rgba, width, height)
    }
}

/// Builds the SVG document for `export_svg` -- pure CPU math over
/// `Scene`/`Camera` with no GPU dependency at all, split out from the
/// `State` method so it's directly unit-testable without a wgpu device.
/// Returns the document and the number of atoms actually drawn (for the
/// log line). See `State::export_svg`'s doc comment for the algorithm and
/// its stated phase-2 cuts.
struct SvgCircle {
    x: f32,
    y: f32,
    r: f32,
    depth: f32,
    color: [u8; 3],
}

/// Every atom of every loaded structure, projected to screen space at its
/// true van der Waals radius and sorted farthest-first (painter's
/// algorithm: SVG paints in document order, so this is a correctness
/// requirement, not an optimization) -- pure CPU math over
/// `Scene`/`Camera`, no GPU dependency, so it's directly unit-testable.
/// `view_proj`/`proj_scale` come from the caller so this and
/// `build_svg_document` always agree on the exact same projection.
fn atom_circles(
    scene: &Scene,
    view_proj: vv_core::glam::Mat4,
    proj_scale: f32,
    width: f32,
    height: f32,
) -> Vec<SvgCircle> {
    let mut circles = Vec::new();
    for (_id, loaded) in scene.structures() {
        let frame = loaded
            .frame
            .min(loaded.structure.frame_count().saturating_sub(1));
        let coords = loaded.structure.frame(frame);
        let positions = coords.positions();
        let colors = crate::gpu_cache::colors_of(loaded, &loaded.rep().coloring, frame);
        let elements = &loaded.structure.topology.element;
        for (atom, &pos) in positions.iter().enumerate() {
            let Some((x, y, depth)) = project_to_pixel(view_proj, width, height, pos) else {
                continue;
            };
            if depth <= 0.0 {
                continue;
            }
            let r = elements[atom].vdw_radius() * proj_scale / depth;
            if r <= 0.0 {
                continue;
            }
            let bytes = colors[atom].to_le_bytes();
            circles.push(SvgCircle {
                x,
                y,
                r,
                depth,
                color: [bytes[0], bytes[1], bytes[2]],
            });
        }
    }
    // Farthest first, so nearer circles (drawn later, per SVG's own
    // paint-order-is-document-order rule) cover the far ones -- the same
    // convention `cull.wgsl` and every impostor shader use for "closer
    // wins."
    circles.sort_by(|a, b| b.depth.total_cmp(&a.depth));
    circles
}

/// Builds the SVG document for `export_svg` -- pure CPU math over
/// `Scene`/`Camera` with no GPU dependency at all, split out from the
/// `State` method so it's directly unit-testable without a wgpu device.
/// Returns the document and the number of atoms actually drawn (for the
/// log line). See `State::export_svg`'s doc comment for the algorithm and
/// its stated phase-2 cuts.
fn build_svg_document(
    scene: &Scene,
    camera: &Camera,
    width: u32,
    height: u32,
    background: wgpu::Color,
) -> (svg::Document, usize) {
    let (fw, fh) = (width as f32, height.max(1) as f32);
    let aspect = fw / fh;
    let proj = camera.proj(aspect);
    let view_proj = proj * camera.view();
    // Same "pixels per world unit at depth 1" `cull.wgsl` uses for the
    // live GPU cull pass (`CameraUniform::proj_scale`), evaluated at this
    // export's own resolution.
    let proj_scale = proj.y_axis.y * fh * 0.5;
    let circles = atom_circles(scene, view_proj, proj_scale, fw, fh);
    let atom_count = circles.len();

    let hex = |c: [u8; 3]| format!("#{:02x}{:02x}{:02x}", c[0], c[1], c[2]);
    let bg = color32_from_wgpu(background);
    let mut document = svg::Document::new()
        .set("viewBox", (0, 0, width, height))
        .set("width", width)
        .set("height", height)
        .add(
            svg::node::element::Rectangle::new()
                .set("x", 0)
                .set("y", 0)
                .set("width", width)
                .set("height", height)
                .set("fill", hex([bg.r(), bg.g(), bg.b()])),
        );
    for c in &circles {
        document = document.add(
            svg::node::element::Circle::new()
                .set("cx", c.x)
                .set("cy", c.y)
                .set("r", c.r)
                .set("fill", hex(c.color)),
        );
    }
    for draw in label_draws(scene) {
        let Some((x, y, depth)) = project_to_pixel(view_proj, fw, fh, draw.world) else {
            continue;
        };
        if depth <= 0.0 {
            continue;
        }
        document = document.add(
            svg::node::element::Circle::new()
                .set("cx", x)
                .set("cy", y)
                .set("r", LABEL_DOT_RADIUS)
                .set("fill", hex(LABEL_DOT_COLOR)),
        );
        document = document.add(
            svg::node::element::Text::new("")
                .set("x", x + 6.0)
                .set("y", y - 6.0)
                .set("font-family", "sans-serif")
                .set("font-size", LABEL_FONT_PX)
                .set("fill", hex(LABEL_TEXT_COLOR))
                .add(svg::node::Text::new(draw.text)),
        );
    }
    for m in measurement_draws(scene) {
        let project = |w| project_to_pixel(view_proj, fw, fh, w).map(|(x, y, _)| (x, y));
        let Some(points) = m
            .path
            .iter()
            .map(|&w| project(w))
            .collect::<Option<Vec<_>>>()
        else {
            continue;
        };
        let path: Vec<String> = points.iter().map(|(x, y)| format!("{x},{y}")).collect();
        document = document.add(
            svg::node::element::Polyline::new()
                .set("points", path.join(" "))
                .set("fill", "none")
                .set("stroke", hex(MEASURE_LINE_COLOR))
                .set("stroke-width", 1.5)
                .set(
                    "stroke-dasharray",
                    format!("{} {}", MEASURE_DASH[0], MEASURE_DASH[1]),
                ),
        );
        if let Some((x, y)) = project(m.anchor) {
            document = document.add(
                svg::node::element::Text::new("")
                    .set("x", x + 4.0)
                    .set("y", y - 4.0)
                    .set("font-family", "sans-serif")
                    .set("font-size", LABEL_FONT_PX)
                    .set("fill", hex(LABEL_TEXT_COLOR))
                    .add(svg::node::Text::new(m.text)),
            );
        }
    }
    for caption in scene.captions() {
        let x = caption.x.clamp(0.0, 1.0) * fw;
        let y = caption.y.clamp(0.0, 1.0) * fh;
        document = document.add(
            svg::node::element::Text::new("")
                .set("x", x)
                .set("y", y + CAPTION_FONT_PX)
                .set("font-family", "sans-serif")
                .set("font-size", CAPTION_FONT_PX)
                .set("fill", hex(CAPTION_TEXT_COLOR))
                .add(svg::node::Text::new(caption.text.clone())),
        );
    }
    (document, atom_count)
}

/// A `--listen` request being run (`State::next_live_line`).
struct LiveRun {
    lines: VecDeque<String>,
    /// Where its output starts in the log.
    log_start: usize,
    failed: bool,
    reply: std::sync::mpsc::Sender<crate::live::Reply>,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn write_png_round_trips_dimensions_and_pixels() {
        let path = std::env::temp_dir().join("vizviz_write_png_test.png");
        let rgba: Vec<u8> = vec![
            10, 20, 30, 255, 40, 50, 60, 255, // row 0
            70, 80, 90, 255, 100, 110, 120, 255, // row 1
        ];
        write_png(&path, &rgba, 2, 2).expect("write png");

        let file = std::fs::File::open(&path).expect("reopen written png");
        let decoder = png::Decoder::new(file);
        let mut reader = decoder.read_info().expect("read png header");
        assert_eq!(reader.info().width, 2);
        assert_eq!(reader.info().height, 2);
        let mut buf = vec![0u8; reader.output_buffer_size()];
        let frame = reader.next_frame(&mut buf).expect("read png frame");
        assert_eq!(&buf[..frame.buffer_size()], rgba.as_slice());

        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn write_png_reports_an_error_for_an_unwritable_path() {
        let path = std::path::PathBuf::from("/nonexistent-directory/vizviz-test.png");
        assert!(write_png(&path, &[0, 0, 0, 255], 1, 1).is_err());
    }

    #[test]
    fn write_jpeg_produces_a_well_formed_file() {
        let path = std::env::temp_dir().join("vizviz_write_jpeg_test.jpg");
        let rgba: Vec<u8> = vec![
            10, 20, 30, 255, 40, 50, 60, 255, 70, 80, 90, 255, 100, 110, 120, 255,
        ];
        write_jpeg(&path, &rgba, 2, 2).expect("write jpeg");
        let bytes = std::fs::read(&path).expect("reopen written jpeg");
        // SOI (Start Of Image) / EOI (End Of Image) markers: the minimal
        // structural check without pulling in a JPEG decoder just for this.
        assert_eq!(&bytes[0..2], &[0xFF, 0xD8], "missing SOI marker");
        assert_eq!(
            &bytes[bytes.len() - 2..],
            &[0xFF, 0xD9],
            "missing EOI marker"
        );
        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn write_screenshot_picks_the_format_from_the_extension() {
        let dir = std::env::temp_dir();
        let rgba = [10u8, 20, 30, 255];
        let png_path = dir.join("vizviz_write_screenshot_test.png");
        write_screenshot(&png_path, &rgba, 1, 1).expect("write png via write_screenshot");
        assert_eq!(&std::fs::read(&png_path).unwrap()[1..4], b"PNG");
        std::fs::remove_file(&png_path).ok();

        for ext in ["jpg", "JPEG"] {
            let jpg_path = dir.join(format!("vizviz_write_screenshot_test.{ext}"));
            write_screenshot(&jpg_path, &rgba, 1, 1).expect("write jpeg via write_screenshot");
            assert_eq!(&std::fs::read(&jpg_path).unwrap()[0..2], &[0xFF, 0xD8]);
            std::fs::remove_file(&jpg_path).ok();
        }
    }

    fn load_1crn_framed(scene: &mut Scene, history: &mut CommandHistory) -> (StructureId, Camera) {
        let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../fixtures/small/1CRN.pdb");
        history
            .dispatch(scene, Command::LoadStructure { path })
            .unwrap();
        let (id, loaded) = scene.structures().next().unwrap();
        let (center, radius) = loaded.structure.frame(0).bounding_sphere().unwrap();
        (id, Camera::framing(center, radius))
    }

    #[test]
    fn atom_circles_draws_every_atom_with_a_positive_radius() {
        let mut scene = Scene::new();
        let mut history = CommandHistory::new(10);
        let (id, camera) = load_1crn_framed(&mut scene, &mut history);
        let atoms = scene.structure(id).unwrap().structure.atom_count();

        let (w, h) = (400.0, 300.0);
        let proj = camera.proj(w / h);
        let view_proj = proj * camera.view();
        let proj_scale = proj.y_axis.y * h * 0.5;
        let circles = atom_circles(&scene, view_proj, proj_scale, w, h);

        // A structure framed to fit the view should have every atom
        // land inside the frustum with a real, positive screen radius.
        assert_eq!(circles.len(), atoms, "every atom should be drawn");
        assert!(circles.iter().all(|c| c.r > 0.0));
        assert!(circles.iter().all(|c| c.depth > 0.0));
    }

    #[test]
    fn atom_circles_are_sorted_farthest_first() {
        let mut scene = Scene::new();
        let mut history = CommandHistory::new(10);
        let (_id, camera) = load_1crn_framed(&mut scene, &mut history);

        let (w, h) = (400.0, 300.0);
        let proj = camera.proj(w / h);
        let view_proj = proj * camera.view();
        let proj_scale = proj.y_axis.y * h * 0.5;
        let circles = atom_circles(&scene, view_proj, proj_scale, w, h);

        // Painter's algorithm: depth must be non-increasing along the
        // list, so SVG's document-order-is-paint-order rule draws nearer
        // circles (later in the list) on top of farther ones.
        assert!(
            circles.windows(2).all(|w| w[0].depth >= w[1].depth),
            "circles must be sorted farthest-first"
        );
    }

    #[test]
    fn build_svg_document_produces_a_well_formed_document() {
        let mut scene = Scene::new();
        let mut history = CommandHistory::new(10);
        let (id, camera) = load_1crn_framed(&mut scene, &mut history);
        let atoms = scene.structure(id).unwrap().structure.atom_count();

        history
            .dispatch(
                &mut scene,
                Command::SetLabel {
                    id,
                    atom: 0,
                    text: Some("N-term".into()),
                },
            )
            .unwrap();
        history
            .dispatch(
                &mut scene,
                Command::SetCaption {
                    caption: vv_scene::Caption {
                        name: "title".into(),
                        x: 0.5,
                        y: 0.05,
                        text: "Crambin".into(),
                    },
                },
            )
            .unwrap();

        let bg = wgpu::Color {
            r: 0.1,
            g: 0.1,
            b: 0.1,
            a: 1.0,
        };
        let (document, count) = build_svg_document(&scene, &camera, 400, 300, bg);
        assert_eq!(count, atoms);
        let text = document.to_string();
        assert!(text.starts_with("<svg"), "{text}");
        assert_eq!(
            text.matches("<circle").count(),
            atoms + 1,
            "atoms + 1 label dot"
        );
        assert!(text.contains("Crambin"), "caption text missing: {text}");
        assert!(text.contains("N-term"), "label text missing: {text}");
    }
}
