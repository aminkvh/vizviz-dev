//! egui panels: the `TabViewer` for the dock area, the menu bar, and the
//! command palette. All state that changes the document goes through
//! `history.dispatch` — panels never mutate `Scene` directly.

use std::time::Instant;

use egui::{Color32, Context, RichText, Ui, Vec2, WidgetText};
use egui_dock::TabViewer as TabViewerTrait;
use fuzzy_matcher::skim::SkimMatcherV2;
use fuzzy_matcher::FuzzyMatcher;

use vv_core::glam::Vec4;
use vv_core::Element;
use vv_render::{AdaptiveLod, Camera, LightingPreset, Pick, Renderer, StylePreset};
use vv_scene::{empty_mask, Command, CommandHistory, Scene, StructureId};

use crate::layout::Tab;
use crate::widgets::{self, count, Variant};
use crate::{measure, viewport};

/// A clip plane's facing. `View` re-faces the camera every frame, so its
/// cut always looks the same on screen; the others are fixed in the
/// scene, so rotating shows the cut from every side.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ClipOrientation {
    View,
    X,
    Y,
    Z,
    /// The view direction at the moment it was captured (`clip_ui`'s
    /// "Capture view"), fixed in the scene from then on.
    Fixed(vv_core::glam::Vec3),
}

/// `RenderSettings::clip`'s source of truth: [`ClipPlane::world`] derives
/// the world-space plane from it and the current camera every frame,
/// since `View` must track the camera and the others don't move with it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ClipPlane {
    pub orientation: ClipOrientation,
    /// Distance from the rotation centre along the (possibly flipped)
    /// normal, in Å; the normal's side is kept.
    pub depth: f32,
    pub flipped: bool,
}

impl ClipPlane {
    /// `RenderSettings::clip`'s `[nx, ny, nz, d]` form, evaluated against
    /// `camera` fresh every frame (only `View` actually depends on it).
    pub fn world(&self, camera: &Camera) -> [f32; 4] {
        let mut n = match self.orientation {
            ClipOrientation::View => camera.orientation * vv_core::glam::Vec3::NEG_Z,
            ClipOrientation::X => vv_core::glam::Vec3::X,
            ClipOrientation::Y => vv_core::glam::Vec3::Y,
            ClipOrientation::Z => vv_core::glam::Vec3::Z,
            ClipOrientation::Fixed(v) => v,
        };
        if self.flipped {
            n = -n;
        }
        n.extend(-n.dot(camera.target) - self.depth).to_array()
    }
}

/// `Clone + PartialEq` so a whole-struct snapshot can be diffed for
/// undo -- see `State::commit_view_edit`.
#[derive(Clone, PartialEq)]
pub struct ViewSettings {
    pub occlusion_culling: bool,
    pub adaptive: bool,
    /// The look last picked: it set the lights below (`set_style`).
    pub style: StylePreset,
    /// Global lighting; materials are per structure (`material`).
    pub lighting: LightingPreset,
    /// The lights the preset set, as edited since.
    pub lights: LightRig,
    /// A slider in Look ▸ Lights ▾ has moved `lights` away from
    /// `lighting`'s own values, so that popover shows "Custom" instead of
    /// a stale preset name. Set here (not derived by comparing `lights`
    /// against a freshly recomputed `LightRig::of(lighting)`) because
    /// `atan2`/`asin` round-tripping isn't guaranteed bit-identical
    /// across call sites, which made an equality check flag "Custom"
    /// spuriously right after picking a preset.
    pub lights_custom: bool,
    /// ACES tonemapping (`vv_render::style::TONEMAP_ACES`).
    pub tonemap: bool,
    /// Seeded from `style`'s suggested background whenever the user picks
    /// a new preset. Kept as its own field (not derived from `style` every
    /// frame) so a screenshot export can override just this without
    /// touching the shading model — see docs/RENDERING.md.
    pub background: wgpu::Color,
    /// Top colour of a vertical gradient down to `background`; `None` is
    /// solid.
    pub background_top: Option<wgpu::Color>,
    /// Smooths the displayed viewport; never applied to exports, where
    /// SSAA supersedes it — see `RenderSettings::fxaa`.
    pub fxaa: bool,
    /// Dark edge lines (`RenderSettings::outline`). Seeded from the preset
    /// like `background`, and overridable the same way.
    pub outline: bool,
    /// Ambient occlusion strength, 0 = off (`RenderSettings::ao`). Seeded
    /// from the preset like `outline`.
    pub ao: f32,
    /// Depth cue strength, 0 = off (`RenderSettings::depth_cue`).
    pub depth_cue: f32,
    /// Shadow strength, 0 = off (`RenderSettings::shadows`).
    pub shadows: f32,
    /// Depth of field, 0 = off (`RenderSettings::dof`).
    pub dof: f32,
    /// Live supersampling: the viewport renders at this multiple of its
    /// on-screen size and egui's linear filter averages it down. Real
    /// subpixel coverage, so thin tubes and sticks stop shimmering where
    /// FXAA's single-sample guess cannot help. 2x costs 4x the pixels.
    pub render_scale: f32,
    /// What a left drag does: R rotate, T translate, S scale,
    /// C center (click an atom to make it the rotation center).
    pub mouse_mode: MouseMode,
    /// Home ▸ Select: what a left drag selects, when not `Click`. Only
    /// meaningful with `mouse_mode` `Rotate`; either one resets the other.
    pub select_tool: crate::select_tool::SelectTool,
    /// The clip plane, or off. [`ClipPlane::world`] turns this into
    /// `RenderSettings::clip`.
    pub clip: Option<ClipPlane>,
    /// Bookmarked camera positions ("Studio ▸ Views ▾ ▸ Save view");
    /// clicking one animates the camera there (a separate, undoable
    /// camera jump -- see `mark_camera_jump`), not part of this edit.
    pub saved_views: Vec<SavedView>,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum MouseMode {
    #[default]
    Rotate,
    Translate,
    Scale,
    Center,
    /// Analyze ▸ Label: click atoms to label them (sticky, unlike
    /// `Center`'s one-shot revert to `Rotate` -- Escape is what exits).
    Label,
}

impl MouseMode {
    fn label(self) -> &'static str {
        match self {
            MouseMode::Rotate => "Rotate",
            MouseMode::Translate => "Translate",
            MouseMode::Scale => "Scale",
            MouseMode::Center => "Center: click an atom",
            MouseMode::Label => "Label: click atoms, Escape to stop",
        }
    }
}

/// What the Inspector shows of the active selection: its atoms per
/// element. Built when the selection changes, not every frame: a
/// 4M-atom selection is a 4M-bit mask.
#[derive(Default)]
pub struct SelectionSummary {
    /// The selection's mask, by address.
    key: Option<usize>,
    by_element: Vec<(&'static str, usize)>,
}

impl SelectionSummary {
    fn of(
        mask: &vv_core::fixedbitset::FixedBitSet,
        loaded: &vv_scene::LoadedStructure,
        key: usize,
    ) -> Self {
        let mut by_element = std::collections::BTreeMap::<&'static str, usize>::new();
        for atom in mask.ones() {
            *by_element
                .entry(Element::symbol(loaded.structure.topology.element[atom]))
                .or_default() += 1;
        }
        Self {
            key: Some(key),
            by_element: by_element.into_iter().collect(),
        }
    }
}

/// A message shown over the viewport for a few seconds: red for
/// failures, blue for information (e.g. what's already loaded).
#[derive(Clone, Debug)]
pub struct Notice {
    pub text: String,
    pub when: Instant,
    pub error: bool,
    /// UX_PATTERNS "Feedback": "an action when useful (Undo, Show)" --
    /// (button label, the command it runs).
    pub action: Option<(&'static str, &'static str)>,
}

impl Notice {
    pub fn error(text: String) -> Self {
        Self {
            text,
            when: Instant::now(),
            error: true,
            action: None,
        }
    }

    pub fn info(text: String) -> Self {
        Self {
            text,
            when: Instant::now(),
            error: false,
            action: None,
        }
    }

    pub fn with_action(mut self, label: &'static str, command: &'static str) -> Self {
        self.action = Some((label, command));
        self
    }
}

/// A Contacts…/Surface area result (UX_PATTERNS "Feedback": "a result
/// card in the Inspector with Copy, plus a toast"), replaced by the next
/// one or dismissed by hand.
#[derive(Clone)]
pub struct ResultCard {
    pub title: String,
    pub body: String,
}

/// The view cube's click/middle-click turn, eased over
/// [`CUBE_ANIM_SECS`] instead of snapping (`overlays::view_cube`'s click
/// sets this instead of calling `Camera::look_from` directly; also used
/// for the double-click residue zoom, which additionally animates
/// `distance`; `State::advance_cube_anim` steps it each frame before the
/// render).
pub struct CubeAnim {
    pub from: glam::Quat,
    pub to: glam::Quat,
    pub from_target: glam::Vec3,
    pub to_target: glam::Vec3,
    pub from_distance: f32,
    pub to_distance: f32,
    pub start: Instant,
}

pub const CUBE_ANIM_SECS: f32 = 0.3;

impl CubeAnim {
    /// Where a click on `side` (the view cube's face label) ends up:
    /// `look_from` on a clone, so this needs no mutable camera access.
    /// Keeps `distance` (the view cube never zooms).
    pub fn to_face(camera: &Camera, side: &str) -> Option<Self> {
        let (from, up) = crate::overlays::face(side)?;
        let mut target = camera.clone();
        target.look_from(from, up);
        Some(Self::to_camera(camera, &target))
    }

    /// Eases orientation, target and distance from `camera` to `target`.
    pub fn to_camera(camera: &Camera, target: &Camera) -> Self {
        Self {
            from: camera.orientation,
            to: target.orientation,
            from_target: camera.target,
            to_target: target.target,
            from_distance: camera.distance,
            to_distance: target.distance,
            start: Instant::now(),
        }
    }

    /// The eased orientation/target/distance at now, and whether the
    /// animation has finished (the caller then drops it and snaps
    /// exactly to the target values, since easing never quite reaches
    /// 1.0 by clock time alone).
    pub fn step(&self) -> (glam::Quat, glam::Vec3, f32, bool) {
        let t = (self.start.elapsed().as_secs_f32() / CUBE_ANIM_SECS).min(1.0);
        let eased = t * t * (3.0 - 2.0 * t); // smoothstep
        let orientation = self.from.slerp(self.to, eased);
        let target = self.from_target.lerp(self.to_target, eased);
        let distance = self.from_distance + (self.to_distance - self.from_distance) * eased;
        (orientation, target, distance, t >= 1.0)
    }
}

// The light list (`StudioLight`/`LightRig`) and saved camera views
// (`SavedView`) live in `studio.rs`, re-exported here since `ViewSettings`
// carries them.
pub use crate::studio::{LightRig, SavedView};

/// The viewport background a theme starts with: #17171C or #F2F4F7.
pub(crate) fn default_background(theme: crate::theme::ThemeMode) -> wgpu::Color {
    let rgb = match theme {
        crate::theme::ThemeMode::Dark => Color32::from_rgb(0x17, 0x17, 0x1C),
        crate::theme::ThemeMode::Light => Color32::from_rgb(0xF2, 0xF4, 0xF7),
    };
    wgpu_from_color32(rgb)
}

impl ViewSettings {
    /// Swaps `was`'s default background for `now`'s; a colour the user
    /// picked stays.
    pub(crate) fn follow_theme(
        &mut self,
        was: crate::theme::ThemeMode,
        now: crate::theme::ThemeMode,
    ) {
        if self.background_top.is_none() && self.background == default_background(was) {
            self.background = default_background(now);
        }
    }
}

impl Default for ViewSettings {
    fn default() -> Self {
        let style = StylePreset::default();
        Self {
            occlusion_culling: true,
            adaptive: true,
            background: default_background(crate::theme::ThemeMode::Dark),
            background_top: None,
            outline: style.outline(),
            ao: 0.0,
            depth_cue: style.depth_cue(),
            shadows: style.lighting_preset().shadows(),
            dof: 0.0,
            lighting: style.lighting_preset(),
            lights: LightRig::of(style.lighting_preset()),
            lights_custom: false,
            tonemap: style.lighting().tonemap > 0.5,
            style,
            fxaa: true,
            render_scale: 1.0,
            mouse_mode: MouseMode::default(),
            select_tool: Default::default(),
            clip: None,
            saved_views: Vec::new(),
        }
    }
}

/// Precomputed, display-only performance numbers (recomputed each frame in
/// `main.rs`; the UI never touches the renderer's timing internals).
#[derive(Clone, Default)]
pub struct Readout {
    pub fps: f32,
    pub frame_ms: f32,
    pub cull_ms: f32,
    pub draw_ms: f32,
    pub hiz_ms: f32,
    pub quad_px_threshold: f32,
}

#[derive(Default)]
pub struct DialogState {
    pub open: bool,
    pub path: String,
}

pub struct ExportDialogState {
    pub open: bool,
    pub path: String,
    pub ssaa: bool,
}

impl Default for ExportDialogState {
    fn default() -> Self {
        Self {
            open: false,
            path: crate::layout::export_dir()
                .join("screenshot.png")
                .display()
                .to_string(),
            ssaa: true,
        }
    }
}

/// What the user asked for; `main.rs` picks this up after the UI pass and
/// does the actual GPU work, since `AppUi` only ever borrows the live
/// `Renderer` immutably (see `viewport_ui`) and export needs a full
/// off-to-the-side render at a different resolution for SSAA.
pub struct ExportRequest {
    pub path: std::path::PathBuf,
    pub ssaa: bool,
}

/// File ▸ Export ▸ Structure… / a structure row's ⋯ ▸ "Export
/// structure…": format, which atoms and which frame(s), open on
/// `structure` so a non-current row exports itself, not the current one.
#[derive(Default)]
pub struct ExportStructureDialogState {
    pub open: bool,
    pub structure: Option<StructureId>,
    pub format: usize,
    pub current_selection: bool,
    pub all_frames: bool,
}

/// `ExportStructureDialogState::format`'s options: `None` infers the
/// writer from whatever extension is typed into the save dialog (as
/// `savestructure` always does); `Some` forces one, appending its
/// extension when the typed name doesn't already end in a matching one.
const EXPORT_FORMATS: &[(&str, Option<&str>)] = &[
    ("By extension", None),
    ("PDB", Some("pdb")),
    ("mmCIF", Some("cif")),
    ("XYZ", Some("xyz")),
    ("PQR", Some("pqr")),
    ("GRO", Some("gro")),
];

/// File > Fetch from RCSB...
pub struct FetchDialogState {
    pub open: bool,
    pub ids: String,
    pub assembly: bool,
    pub assembly_number: u32,
}

impl Default for FetchDialogState {
    fn default() -> Self {
        Self {
            open: false,
            ids: String::new(),
            assembly: false,
            assembly_number: 1,
        }
    }
}

/// The Preferences dialog, gear in the quick-access
/// bar or `preferences`. Every field it edits (`UiPrefs`, `ViewSettings`,
/// `AdaptiveLod`) applies and saves immediately; this only tracks the
/// dialog's own open state and Danger's "reset all" confirmation ("Delete
/// saved layouts" reuses `WorkspaceDialogState::deleting`).
#[derive(Default)]
pub struct PreferencesDialogState {
    pub open: bool,
    pub confirm_reset: bool,
}

#[derive(Default)]
pub struct WorkspaceDialogState {
    pub open: bool,
    pub name: String,
    /// The saved workspaces, each with Delete (`deletelayout` alone).
    pub deleting: bool,
    /// The workspace whose deletion is waiting for a yes.
    pub confirm: Option<String>,
}

#[derive(Default)]
pub struct PaletteState {
    pub open: bool,
    pub query: String,
}

/// The one-line console at the bottom of the Log tab.
#[derive(Default)]
pub struct ConsoleState {
    pub input: String,
}

/// A workspace change the UI asked for. Applied by `State` after the UI
/// pass, because the dock state is mutably borrowed by the dock area
/// while panels draw (see `menu_bar`).
#[derive(Clone, Debug)]
pub enum LayoutRequest {
    /// Open the current arrangement at every launch (`startlayout`).
    SetStart,
    /// Launch with the last arrangement again (`startlayout reset`).
    ResetStart,
    Switch(String),
    Save(String),
    /// Show a panel, reopening it if closed.
    OpenPanel(Tab),
    /// Close a panel (View ▸ Panels ▾'s checklist; a no-op if already
    /// closed).
    ClosePanel(Tab),
    /// Reopen a closed panel as a tab of that node.
    AddPanel(egui_dock::NodePath, Tab),
    /// Replace the dock layout outright (a loaded session's `dock_layout`).
    LoadState(Box<egui_dock::DockState<Tab>>),
}

/// Playback state of the Timeline tab. View state like `State::frames`,
/// not document state.
pub struct TimelineState {
    pub playing: bool,
    pub fps: f32,
    pub looping: bool,
    /// When the last automatic frame step happened.
    pub last_step: Option<Instant>,
}

impl Default for TimelineState {
    fn default() -> Self {
        Self {
            playing: false,
            fps: 10.0,
            looping: true,
            last_step: None,
        }
    }
}

/// The structure the timeline scrubs, with its frame count: the current
/// one if it has more than one frame, otherwise the first that does.
pub fn timeline_target(scene: &Scene) -> Option<(StructureId, usize)> {
    let multi = |id: StructureId| {
        scene
            .structure(id)
            .map(|s| s.structure.frame_count())
            .filter(|&n| n > 1)
            .map(|n| (id, n))
    };
    crate::gpu_cache::current_structure(scene)
        .and_then(multi)
        .or_else(|| scene.structures().find_map(|(id, _)| multi(id)))
}

/// The mode chip's "?": every mouse binding, the left-drag row
/// naming whichever mode is current.
fn mouse_bindings_popup(ui: &mut Ui, left_drag: String) {
    ui.set_max_width(260.0);
    widgets::section(ui, "Mouse");
    let rows = [
        ("Left drag", left_drag),
        (
            "Middle drag",
            "Roll around the axis out of the screen".to_owned(),
        ),
        ("Right drag", "Pan".to_owned()),
        ("Wheel", "Zoom".to_owned()),
        (
            "Shift+wheel",
            "Dolly: moves the camera, cutting the scene".to_owned(),
        ),
        ("Click an atom/bond", "Select".to_owned()),
        ("Ctrl+click", "Add to the selection".to_owned()),
        (
            "Ctrl+click, again",
            "Extend a measure chain: distance, angle, dihedral".to_owned(),
        ),
        (
            "Double-click an atom",
            "Center, zoom to it, and show its neighborhood".to_owned(),
        ),
        (
            "View cube click/middle-click",
            "Look from that face".to_owned(),
        ),
    ];
    widgets::table(ui, "mouse-bindings", &rows);
}

/// Rotation by drag is a trackball -- the molecule's near surface
/// follows the cursor. `delta` is in egui points (y down); turning the
/// scene by `t` about screen axis `a` moves the near point by `t * (a × z)`,
/// so a drag `(dx, dy)` needs `a = (dy, dx, 0)`.
/// The turn about the screen's out-of-view axis that carries the pointer
/// from `before` to `after` (both relative to the viewport centre, in egui
/// points, y down): positive is counter-clockwise as seen.
fn roll_angle(before: Vec2, after: Vec2) -> f32 {
    if before.length() < 4.0 || after.length() < 4.0 {
        return 0.0;
    }
    -(before.x * after.y - before.y * after.x).atan2(before.dot(after))
}

fn trackball_rotate(camera: &mut Camera, delta: Vec2) {
    let len = delta.length();
    if len < f32::EPSILON {
        return;
    }
    let axis = glam::Vec3::new(delta.y, delta.x, 0.0).normalize();
    camera.rotate_scene(axis, len * 0.005);
}

/// `wgpu::Color` is linear; `Color32` is sRGB. `Rgba` is egui's linear
/// type, so it's the bridge in both directions.
pub(crate) fn color32_from_wgpu(c: wgpu::Color) -> Color32 {
    Color32::from(egui::Rgba::from_rgba_premultiplied(
        c.r as f32, c.g as f32, c.b as f32, c.a as f32,
    ))
}

/// Projects a world position through `view_proj` to a pixel offset across
/// a `width x height` canvas (origin top-left), plus its camera-space
/// depth (view-space distance in front of the camera -- the same
/// quantity `cull.wgsl` calls `-vc.z`, usable both to scale a screen
/// radius via `CameraUniform::proj_scale` and to depth-sort primitives
/// back-to-front). Shared by the live viewport overlay
/// (`labels_and_captions_ui`, where the canvas is the viewport `Rect`),
/// the screenshot exporter (`state::export_screenshot`, where it's the
/// raw export buffer), and SVG export (`state::export_svg`), so a label
/// or atom lands in the same place in all three. `None` when the point is
/// behind the camera or far enough outside the frustum that drawing it
/// would just be screen-edge clutter.
pub(crate) fn project_to_pixel(
    view_proj: vv_core::glam::Mat4,
    width: f32,
    height: f32,
    world: vv_core::glam::Vec3,
) -> Option<(f32, f32, f32)> {
    let clip = view_proj * Vec4::new(world.x, world.y, world.z, 1.0);
    if clip.w <= 1e-4 {
        return None;
    }
    let ndc_x = clip.x / clip.w;
    let ndc_y = clip.y / clip.w;
    if !(-1.5..=1.5).contains(&ndc_x) || !(-1.5..=1.5).contains(&ndc_y) {
        return None;
    }
    Some((
        (ndc_x * 0.5 + 0.5) * width,
        (1.0 - (ndc_y * 0.5 + 0.5)) * height,
        clip.w,
    ))
}

/// `world`'s pixel position across a `width x height` canvas with no
/// frustum bound, or `None` behind the camera: for lines the caller clips.
pub(crate) fn in_front_pixel(
    view_proj: vv_core::glam::Mat4,
    width: f32,
    height: f32,
    world: vv_core::glam::Vec3,
) -> Option<(f32, f32)> {
    let clip = view_proj * Vec4::new(world.x, world.y, world.z, 1.0);
    (clip.w > 1e-4).then(|| {
        (
            (clip.x / clip.w * 0.5 + 0.5) * width,
            (1.0 - (clip.y / clip.w * 0.5 + 0.5)) * height,
        )
    })
}

/// "Rendering 45%" over a thin bar, in the viewport's lower-right corner.
fn render_progress_badge(ui: &Ui, viewport: egui::Rect, progress: f32) {
    let t = crate::theme::Tokens::current(ui.ctx());
    let painter = ui.painter_at(viewport);
    let text = format!("Rendering {:.0}%", progress * 100.0);
    let galley = painter.layout_no_wrap(text, egui::TextStyle::Small.resolve(ui.style()), t.text);
    let pad = crate::theme::space::GAP;
    let size = galley.size() + egui::vec2(2.0 * pad, 2.0 * pad);
    let corner = viewport.right_bottom() - egui::vec2(pad, pad);
    let badge = egui::Rect::from_min_size(corner - size, size);
    painter.rect(
        badge,
        crate::theme::radius::CONTROL,
        t.surface,
        egui::Stroke::new(1.0, t.border),
        egui::StrokeKind::Inside,
    );
    painter.galley(badge.min + egui::vec2(pad, pad - 2.0), galley, t.text);
    let bar = egui::Rect::from_min_max(
        egui::pos2(badge.left() + pad, badge.bottom() - 4.0),
        egui::pos2(
            badge.left() + pad + progress.clamp(0.0, 1.0) * (badge.width() - 2.0 * pad),
            badge.bottom() - 2.0,
        ),
    );
    painter.rect_filled(bar, 1.0, t.primary);
}

/// An annotation's text on a dark pill, legible on any background or
/// atom colour; `at` is the text's bottom-left.
fn annotation_tag(ui: &Ui, painter: &egui::Painter, at: egui::Pos2, text: &str) {
    let color = Color32::from_rgb(
        LABEL_TEXT_COLOR[0],
        LABEL_TEXT_COLOR[1],
        LABEL_TEXT_COLOR[2],
    );
    let galley = painter.layout_no_wrap(
        text.to_owned(),
        egui::TextStyle::Small.resolve(ui.style()),
        color,
    );
    let rect = egui::Rect::from_min_size(at - egui::vec2(0.0, galley.size().y), galley.size());
    painter.rect_filled(
        rect.expand2(egui::vec2(4.0, 2.0)),
        3.0,
        Color32::from_black_alpha(170),
    );
    painter.galley(rect.min, galley, color);
}

/// One label as drawn: which atom it's anchored to (so a right-click
/// hit-test can remove exactly that one), its current-frame world
/// position, and its text.
pub(crate) struct LabelDraw {
    pub id: StructureId,
    pub atom: u32,
    pub world: vv_core::glam::Vec3,
    pub text: String,
}

/// Every structure-anchored label's current-frame world position and
/// text, across all loaded structures -- shared by the live viewport
/// overlay and the screenshot exporter for the same reason
/// `project_to_pixel` is.
pub(crate) fn label_draws(scene: &Scene) -> Vec<LabelDraw> {
    let mut draws = Vec::new();
    for (id, loaded) in scene.structures() {
        if loaded.labels.is_empty() {
            continue;
        }
        let frame = loaded
            .frame
            .min(loaded.structure.frame_count().saturating_sub(1));
        let coords = loaded.structure.frame(frame);
        let positions = coords.positions();
        for (&atom, text) in &loaded.labels {
            if let Some(&world) = positions.get(atom as usize) {
                draws.push(LabelDraw {
                    id,
                    atom,
                    world,
                    text: text.clone(),
                });
            }
        }
    }
    draws
}

/// One measurement as drawn: the atoms' current-frame positions (the
/// dashed path through them), where its value goes, and the value.
pub(crate) struct MeasurementDraw {
    pub path: Vec<vv_core::glam::Vec3>,
    pub anchor: vv_core::glam::Vec3,
    pub text: String,
}

/// Every measurement on screen, measured at its structure's current
/// frame -- so they follow a playing trajectory. The value sits at the
/// middle of a distance, on an angle's vertex, and at the middle of a
/// dihedral's central bond. Shared by the live overlay and both exports,
/// as `label_draws` is.
pub(crate) fn measurement_draws(scene: &Scene) -> Vec<MeasurementDraw> {
    let mut draws = Vec::new();
    for (_id, loaded) in scene.structures() {
        if loaded.measurements.is_empty() {
            continue;
        }
        let coords = loaded.structure.frame(loaded.frame);
        let positions = coords.positions();
        for m in &loaded.measurements {
            let path: Vec<_> = m.atoms().iter().map(|&a| positions[a as usize]).collect();
            let anchor = match path[..] {
                [a, b] | [_, a, b, _] => (a + b) * 0.5,
                [_, b, _] => b,
                _ => unreachable!("2 to 4 atoms"),
            };
            draws.push(MeasurementDraw {
                text: m.text(positions),
                path,
                anchor,
            });
        }
    }
    draws
}

/// Label/caption drawing constants shared between the live egui overlay
/// (`labels_and_captions_ui`) and the CPU-rasterized export path
/// (`state::export_screenshot`, via `crate::text`), so an exported figure
/// matches what the live viewport shows. Pixel sizes are independent
/// approximations of egui's own `TextStyle::Small`/`Heading` (a different
/// font renders them), not an exact match -- exact glyph-for-glyph parity
/// isn't the goal, matching proportions is.
/// The Selections panel's add-a-selection field, a stable id so the tool
/// strip's button can focus it from outside the panel.
pub(crate) const SELECTION_FIELD_ID: &str = "vizviz-selection-expr";

/// The Structures panel's rename-in-progress key (see `inline_rename_ui`).
fn structure_rename_key() -> egui::Id {
    egui::Id::new("structure-rename")
}

/// A selection row's expression-being-edited key (see `inline_rename_ui`);
/// keyed by `(StructureId, RepId)`, distinct from `structure_rename_key`'s
/// own storage so renaming a structure and editing a selection's
/// expression never collide.
fn selection_expr_key() -> egui::Id {
    egui::Id::new("selection-expr-edit")
}

/// Which row `storage_key` is mid-rename, if any (`inline_rename_ui`'s
/// own storage, read without disturbing it).
pub(crate) fn renaming_id<Id: Clone + Send + Sync + 'static>(
    ui: &Ui,
    storage_key: egui::Id,
) -> Option<Id> {
    ui.data_mut(|d| d.get_temp::<(Id, String)>(storage_key))
        .map(|(id, _)| id)
}

/// UX_PATTERNS "Edit: inline": a row mid-rename, its label replaced by a
/// text field -- Enter saves; Escape or clicking elsewhere reverts. `storage_key` scopes both
/// which row and its draft text, so several lists (Structures, Sets)
/// each rename independently; the caller checks `renaming_id` first to
/// decide whether to draw this instead of the row's usual look, and
/// dispatches its own rename command with the name this returns (only on
/// a real, non-empty, changed Enter -- `None` otherwise, including every
/// frame the field is still open).
pub(crate) fn inline_rename_ui<
    Id: Clone + PartialEq + std::hash::Hash + std::fmt::Debug + Send + Sync + 'static,
>(
    ui: &mut Ui,
    storage_key: egui::Id,
    id: Id,
    original: &str,
) -> Option<String> {
    let mut text = ui
        .data_mut(|d| d.get_temp::<(Id, String)>(storage_key))
        .filter(|(row, _)| *row == id)
        .map(|(_, text)| text)
        .unwrap_or_else(|| original.to_owned());
    let field_id = storage_key.with(format!("{id:?}"));
    let opened = storage_key.with("opened");
    let edit = ui
        .add(
            egui::TextEdit::singleline(&mut text)
                .id(field_id)
                .desired_width(f32::INFINITY),
        )
        .on_hover_text("Enter to save; Escape or click elsewhere to cancel");
    if !ui.data(|d| d.get_temp::<bool>(opened).unwrap_or(false)) {
        edit.request_focus();
        ui.data_mut(|d| d.insert_temp(opened, true));
    }
    if !edit.lost_focus() {
        ui.data_mut(|d| d.insert_temp(storage_key, (id, text)));
        return None;
    }
    ui.data_mut(|d| {
        d.remove::<(Id, String)>(storage_key);
        d.remove::<bool>(opened);
    });
    let enter = ui.input(|i| i.key_pressed(egui::Key::Enter));
    let trimmed = text.trim();
    (enter && !trimmed.is_empty() && trimmed != original).then(|| trimmed.to_owned())
}

pub(crate) const LABEL_DOT_RADIUS: f32 = 2.5;
pub(crate) const LABEL_DOT_COLOR: [u8; 3] = [255, 210, 90];
pub(crate) const LABEL_TEXT_COLOR: [u8; 3] = [255, 255, 255];
pub(crate) const LABEL_FONT_PX: f32 = 14.0;
pub(crate) const MEASURE_LINE_COLOR: [u8; 3] = [255, 210, 90];
/// Dash and gap of a measurement's line, in pixels.
pub(crate) const MEASURE_DASH: [f32; 2] = [6.0, 4.0];
pub(crate) const CAPTION_TEXT_COLOR: [u8; 3] = [255, 255, 255];
pub(crate) const CAPTION_FONT_PX: f32 = 26.0;

fn wgpu_from_color32(c: Color32) -> wgpu::Color {
    let rgba = egui::Rgba::from(c);
    wgpu::Color {
        r: rgba.r() as f64,
        g: rgba.g() as f64,
        b: rgba.b() as f64,
        a: rgba.a() as f64,
    }
}

/// Built fresh each frame in `main.rs` from `State`'s fields. Implements
/// `TabViewer` for the dock area and also draws the menu bar and palette,
/// which aren't tabs but need the same borrows.
pub struct AppUi<'a> {
    pub scene: &'a mut Scene,
    pub history: &'a mut CommandHistory,
    /// View settings and camera jumps in the same undo history as
    /// scene edits (`crate::app_history`).
    pub app_history: &'a mut crate::app_history::AppHistory,
    /// Set by `view reset`/`view face`/a view-cube click (the explicit
    /// "camera jump" call sites) to the camera *before* they
    /// move it; `State::redraw` turns this into an `app_history` entry
    /// after the frame, once the jump (and any cube-click animation
    /// start) is done. Ordinary drag/pan/zoom/wheel never set this.
    pub camera_jump: &'a mut Option<Camera>,
    /// `AppUi::commit_view_edit`'s baseline (see `State`'s own doc on
    /// its copy) -- also corrected by `app_undo`/`app_redo` so the
    /// frame right after an undo doesn't read as a brand new edit.
    pub view_baseline: &'a mut ViewSettings,
    /// The Info panel's own target structure, when a row's "File info"
    /// picked one other than the current structure.
    pub info_target: &'a mut Option<StructureId>,
    /// A Structures-panel row click's override of "the current structure"
    /// (`current()`'s doc). `None` once its structure closes.
    pub current_override: &'a mut Option<StructureId>,
    /// A selection row hovered or current this frame (`State::
    /// row_highlight`'s doc): the viewport outlines its atoms.
    pub row_highlight: &'a mut Option<(StructureId, vv_scene::RepId)>,
    /// The atom under the pointer, with no button held (`State::hover`).
    pub hover: &'a mut crate::viewport::HoverPick,
    /// The viewport's pending Ctrl+click measurement chain.
    pub measure: &'a mut crate::measure::MeasureChain,
    /// The Inspector's Contacts…/Surface area result card.
    pub result_card: &'a mut Option<ResultCard>,
    pub camera: &'a mut Camera,
    pub view: &'a mut ViewSettings,
    pub renderer: &'a Renderer,
    /// Read-only, for a rep row's own loading spinner.
    pub gpu_cache: &'a crate::gpu_cache::GpuCache,
    /// What each item the last frame drew stands for (see `State`).
    pub draw_order: &'a [crate::gpu_cache::DrawSource],
    pub timeline: &'a mut TimelineState,
    pub movie: &'a mut crate::movie_state::MovieState,
    pub lod: &'a mut AdaptiveLod,
    pub viewport_texture: Option<egui::TextureId>,
    pub requested_viewport_size: &'a mut Option<Vec2>,
    pub dialog: &'a mut DialogState,
    pub export_dialog: &'a mut ExportDialogState,
    pub export_structure_dialog: &'a mut ExportStructureDialogState,
    pub fetch_dialog: &'a mut FetchDialogState,
    pub export_request: &'a mut Option<ExportRequest>,
    pub render_request: &'a mut Option<crate::render::RenderRequest>,
    pub workspace_dialog: &'a mut WorkspaceDialogState,
    pub prefs_dialog: &'a mut PreferencesDialogState,
    /// User-saved layouts (name -> dock state), persisted alongside the
    /// live layout. Not the built-in ones (`layout::built_in_workspaces`)
    /// — those are generated fresh and never stored here, so they can't
    /// accumulate duplicates across launches; see `layout.rs`.
    pub workspaces: &'a mut Vec<(String, egui_dock::DockState<Tab>)>,
    pub layout_request: &'a mut Option<LayoutRequest>,
    /// Panels not in the dock, offered by each tab bar's add button.
    pub closed_panels: Vec<Tab>,
    /// The dock layout as this frame began, for `savesession` (the live
    /// one is borrowed by the dock area while panels draw).
    pub dock_layout: &'a egui_dock::DockState<Tab>,
    pub palette: &'a mut PaletteState,
    /// Atoms of the active selection in click order (see `select_atoms`).
    pub pick_order: &'a mut Vec<u32>,
    /// A file is being dragged over the window.
    pub drop_hover: bool,
    /// Last notice over the viewport (see `State::notice`).
    pub notice: &'a mut Option<Notice>,
    /// The view cube's in-flight turn, if any.
    pub cube_anim: &'a mut Option<CubeAnim>,
    /// A running path-traced render's progress, 0..1.
    pub render_progress: Option<f32>,
    pub console: &'a mut ConsoleState,
    pub quit_requested: &'a mut bool,
    pub selection_summary: &'a mut SelectionSummary,
    /// The scene's edits when a session was last saved or opened.
    pub saved_edits: &'a mut u64,
    /// Asking whether to save before quitting.
    pub quit_dialog: &'a mut bool,
    /// Quitting was confirmed (saved, or chosen without saving).
    pub quit_confirmed: &'a mut bool,
    /// The quit dialog's own "Don't remind me again" checkbox.
    pub quit_dont_ask: &'a mut bool,
    /// The start card (`start_card_ui`) stays hidden for the rest of this
    /// run once closed, even with nothing loaded. `welcome` clears it.
    pub start_card_dismissed: &'a mut bool,
    /// Studio ▸ Frame preview / `framing` (see `studio_frame`'s own doc
    /// on `State`).
    pub studio_frame: &'a mut bool,
    /// `window WxH`: the window size to take, in logical pixels.
    pub window_request: &'a mut Option<(u32, u32)>,
    /// Where to write a capture of the whole window (`uishot`), and how
    /// many more frames to draw first.
    pub ui_shot: &'a mut Option<(std::path::PathBuf, u8)>,
    pub prefs: &'a mut crate::theme::UiPrefs,
    pub ribbon: &'a mut crate::ribbon::RibbonState,
    pub log: &'a mut Vec<String>,
    pub readout: Readout,
}

impl AppUi<'_> {
    pub(crate) fn dispatch(&mut self, command: Command) {
        match self.history.dispatch(self.scene, command) {
            Ok(()) => self.app_history.push(crate::app_history::Edit::Scene),
            Err(e) => self.fail(e.to_string()),
        }
    }

    /// `undo`/`redo`: across scene edits and view/camera edits
    /// alike. A `Scene` marker replays on `vv_scene::CommandHistory`,
    /// kept in lockstep with it by `dispatch`; `View`/`Camera` just
    /// swap in their own snapshot.
    pub(crate) fn app_undo(&mut self) -> String {
        use crate::app_history::Edit;
        match self.app_history.pop_undo() {
            None => "nothing to undo".into(),
            Some(Edit::Scene) => match self.history.undo(self.scene) {
                Ok(true) => {
                    self.app_history.push_undone(Edit::Scene);
                    "undone".into()
                }
                Ok(false) => "nothing to undo".into(),
                Err(e) => e.to_string(),
            },
            Some(Edit::View(before, after)) => {
                *self.view = (*before).clone();
                *self.view_baseline = (*before).clone();
                self.app_history.push_undone(Edit::View(before, after));
                "undone".into()
            }
            Some(Edit::Camera(before, after)) => {
                *self.camera = (*before).clone();
                self.app_history.push_undone(Edit::Camera(before, after));
                "undone".into()
            }
        }
    }

    pub(crate) fn app_redo(&mut self) -> String {
        use crate::app_history::Edit;
        match self.app_history.pop_redo() {
            None => "nothing to redo".into(),
            Some(Edit::Scene) => match self.history.redo(self.scene) {
                Ok(true) => {
                    self.app_history.push_done(Edit::Scene);
                    "redone".into()
                }
                Ok(false) => "nothing to redo".into(),
                Err(e) => e.to_string(),
            },
            Some(Edit::View(before, after)) => {
                *self.view = (*after).clone();
                *self.view_baseline = (*after).clone();
                self.app_history.push_done(Edit::View(before, after));
                "redone".into()
            }
            Some(Edit::Camera(before, after)) => {
                *self.camera = (*after).clone();
                self.app_history.push_done(Edit::Camera(before, after));
                "redone".into()
            }
        }
    }

    /// Commits a whole-`ViewSettings` edit (lights, background,
    /// effects, clip, look/light presets) when it differs from the last
    /// committed baseline and nothing is still being dragged (a slider,
    /// a color wheel) -- one undo step per committed change, live
    /// preview while dragging, matching UX_PATTERNS "Edit: Saving" for
    /// scene edits. Called once, at the end of every frame's UI pass.
    ///
    /// Performance settings (`occlusion_culling`, `adaptive`, `fxaa`,
    /// `render_scale` -- Preferences ▸ Performance, App-scoped, and
    /// "Reset all preferences" explicitly says it isn't undoable) live
    /// on `ViewSettings` too but aren't undoable, so they're
    /// masked out of the comparison: on their own they never trigger a
    /// step. A change to one of them landing in the same frame as a
    /// real view edit still rides along in that step's "after" (and so
    /// back out on undo) -- accepted as a rare coincidence rather than
    /// engineered around.
    pub(crate) fn commit_view_edit(&mut self, ctx: &Context) {
        if ctx.input(|i| i.pointer.any_down()) {
            return;
        }
        let mut comparable = self.view.clone();
        comparable.occlusion_culling = self.view_baseline.occlusion_culling;
        comparable.adaptive = self.view_baseline.adaptive;
        comparable.fxaa = self.view_baseline.fxaa;
        comparable.render_scale = self.view_baseline.render_scale;
        // Mouse mode is interaction state (R/T/S/C, or Escape exiting
        // Label), not a "look" that belongs in undo history.
        comparable.mouse_mode = self.view_baseline.mouse_mode;
        comparable.select_tool = self.view_baseline.select_tool;
        if comparable == *self.view_baseline {
            return;
        }
        let before = self.view_baseline.clone();
        let after = self.view.clone();
        *self.view_baseline = after.clone();
        self.app_history.push(crate::app_history::Edit::View(
            Box::new(before),
            Box::new(after),
        ));
        self.offer_undo_on_current_toast();
    }

    /// `commit_view_edit`/`commit_camera_jump` run after the frame's own
    /// `run_command_logged` already decided whether to attach Undo to
    /// its toast (checking `CommandHistory::edits()`, which a view/
    /// camera edit never grows) -- so once one of them actually commits,
    /// this reaches back and adds Undo to whatever toast is still
    /// showing, if it doesn't already have an action.
    fn offer_undo_on_current_toast(&mut self) {
        if let Some(notice) = self.notice.as_mut() {
            if notice.action.is_none() {
                notice.action = Some(("Undo", "undo"));
            }
        }
    }

    /// Records the camera's state right before a deliberate "jump"
    /// (`view reset`, `view face`, a view-cube click -- its three call
    /// sites) so `commit_camera_jump` can turn it into an `app_history`
    /// entry; keeps the *first* jump's before-state if more than one
    /// runs in the same frame (a script), rather than overwriting it.
    pub(crate) fn mark_camera_jump(&mut self) {
        if self.camera_jump.is_none() {
            *self.camera_jump = Some(self.camera.clone());
        }
    }

    /// Turns a pending `camera_jump` (view reset/face, a cube click)
    /// into one `app_history` entry, if it actually moved the camera.
    /// Called once, at the end of every frame's UI pass. A cube click's
    /// jump doesn't finish moving the camera until `cube_anim` eases it
    /// there over several frames (`State::advance_cube_anim`, stepped
    /// before this runs each frame), so this waits for that to settle
    /// rather than consuming `camera_jump` -- and committing an empty
    /// step -- on the click's own first frame.
    pub(crate) fn commit_camera_jump(&mut self) {
        if self.cube_anim.is_some() {
            return;
        }
        let Some(before) = self.camera_jump.take() else {
            return;
        };
        if *self.camera == before {
            return;
        }
        self.app_history.push(crate::app_history::Edit::Camera(
            Box::new(before),
            Box::new(self.camera.clone()),
        ));
        self.offer_undo_on_current_toast();
    }

    /// Logs a failure and shows it over the viewport for a few seconds.
    pub(crate) fn fail(&mut self, message: String) {
        self.log.push(format!("error: {message}"));
        *self.notice = Some(Notice::error(message));
    }

    /// The structure ribbon actions and the Structures panel act on:
    /// whichever row was last clicked (`current_override`), if it's still
    /// loaded, else the most recently loaded structure.
    pub(crate) fn current(&self) -> Option<StructureId> {
        self.current_override
            .filter(|&id| self.scene.structure(id).is_some())
            .or_else(|| crate::gpu_cache::current_structure(self.scene))
    }

    /// Any modal dialog is up (Open, Fetch, Export, Export structure,
    /// Preferences, workspace save/delete, quit): the viewport's own
    /// Escape (clear-selection) steps aside for it; the same check
    /// serves the "closes the top popover" priority.
    pub(crate) fn any_dialog_open(&self) -> bool {
        self.dialog.open
            || self.export_dialog.open
            || self.export_structure_dialog.open
            || self.fetch_dialog.open
            || self.workspace_dialog.open
            || self.workspace_dialog.deleting
            || self.workspace_dialog.confirm.is_some()
            || self.prefs_dialog.open
            || self.prefs_dialog.confirm_reset
            || *self.quit_dialog
    }

    /// Applies a look: lighting, effects and every rep's material (one
    /// undo step). The background stays the person's own choice.
    /// A look preset: lights and effects only -- a structure's
    /// material is its own choice (Represent > Material), untouched here.
    /// `open_path` still seeds a newly loaded structure's material from
    /// `style.material_preset()`, which is a load default, not this
    /// coupling.
    pub(crate) fn set_style(&mut self, style: StylePreset) {
        self.view.lights = LightRig::of(style.lighting_preset());
        self.view.lights_custom = false;
        self.view.lighting = style.lighting_preset();
        self.view.tonemap = style.lighting().tonemap > 0.5;
        self.view.style = style;
    }

    /// Picks a lighting preset. Effects are separate (Preferences ▸
    /// Effects), so a preset never turns one back on.
    pub(crate) fn set_lighting(&mut self, lighting: LightingPreset) {
        self.view.lighting = lighting;
        self.view.lights = LightRig::of(lighting);
        self.view.lights_custom = false;
    }

    /// Sets the material of structure `id`'s current rep.
    pub(crate) fn set_material(&mut self, id: StructureId, material: vv_scene::Material) {
        let Some(rep) = self.scene.structure(id).map(|s| s.rep()) else {
            return;
        };
        if rep.material != material {
            let rep = rep.id;
            self.dispatch(Command::SetMaterial { id, rep, material });
        }
    }

    /// Loads a structure from a typed, pasted, or dropped path. Tolerates
    /// what people actually paste: surrounding quotes (Explorer's "Copy as
    /// path"), a leading `~`, stray whitespace.
    pub(crate) fn open_path(&mut self, path: String) {
        let Some(path) = clean_path(&path) else {
            return;
        };
        if !path.exists() {
            self.fail(format!("no such file: {}", path.display()));
            return;
        }
        self.log.push(format!("loading {}…", path.display()));
        let already: Vec<String> = self
            .scene
            .structures()
            .map(|(id, s)| format!("{} (#{})", s.label, id.to_raw()))
            .collect();
        let opened = path.clone();
        self.dispatch(Command::LoadStructure { path });
        if self.scene.structures().count() > already.len() {
            *self.info_target = None;
            self.prefs.add_recent(opened);
            let loaded = self.scene.structures().next_back().map(|(id, _)| id);
            if let Some(id) = loaded {
                self.set_material(id, self.view.style.material_preset());
            }
            // Pivot on what was just opened, keeping the user's projection
            // and viewing angle: a load must not reset how they set up the
            // view, but it must rotate about the structure they just added,
            // not wherever the rest of the scene happens to be.
            self.frame_current(true);
            if !already.is_empty() {
                let text = format!("added to the scene; already open: {}", already.join(", "));
                self.log.push(text.clone());
                *self.notice = Some(Notice::info(text));
            }
        }
    }

    /// The platform's own file picker (blocking; the frame simply waits).
    /// Returns `None` both when the user cancels and when the platform has
    /// no picker to offer (a Linux box with no desktop portal), which is
    /// why "Open Path..." keeps the plain text dialog reachable.
    pub(crate) fn open_file_dialog(&mut self) {
        let picked = rfd::FileDialog::new()
            .set_title("Open Structure")
            .add_filter("Structures", &["cif", "mmcif", "pdbx", "pdb", "ent", "gz"])
            .add_filter("All files", &["*"])
            .pick_file();
        if let Some(path) = picked {
            self.open_path(path.display().to_string());
        }
    }

    /// Frames the current structure (`AppUi::current`: the most recently
    /// loaded one still open), resetting the viewing angle too. Called
    /// after `load`/`fetch`/`loadtrajectory` on the console, and exposed
    /// as `view reset` for when the user has panned/zoomed away and wants
    /// to snap back to it.
    pub(crate) fn reset_camera(&mut self) {
        self.frame_current(false);
    }

    /// Fits the current structure in view, pivoting on its own centre.
    /// Always keeps the projection (a user setting, not framing); with
    /// `keep_orientation`, also the viewing angle, which is what loading
    /// another file wants (`open_path`).
    pub(crate) fn frame_current(&mut self, keep_orientation: bool) {
        match self
            .current()
            .and_then(|id| crate::gpu_cache::frame_structure(self.scene, id))
        {
            Some(mut camera) => {
                camera.projection = self.camera.projection;
                if keep_orientation {
                    camera.orientation = self.camera.orientation;
                }
                *self.camera = camera;
            }
            None => self.log.push("no structure to frame".into()),
        }
    }

    /// Makes `atom` of structure `id` the rotation center without moving
    /// the picture: later rotations pivot on it.
    pub(crate) fn center_on_atom(&mut self, id: StructureId, atom: u32) {
        if let Some(loaded) = self.scene.structure(id) {
            let frame = loaded
                .frame
                .min(loaded.structure.frame_count().saturating_sub(1));
            if let Some(&p) = loaded.structure.frame(frame).positions().get(atom as usize) {
                self.camera.pivot = Some(p);
                self.log.push(format!(
                    "rotating about atom {} of #{}",
                    atom + 1,
                    id.to_raw()
                ));
            }
        }
    }

    /// Analyze ▸ Label's click-atoms mode: labels `atom` with its
    /// own compact tag (`atom_tag`); a re-click on an already-labeled atom
    /// is a no-op rather than re-labeling it with the same text.
    pub(crate) fn label_atom(&mut self, id: StructureId, atom: u32) {
        let Some(loaded) = self.scene.structure(id) else {
            return;
        };
        if loaded.labels.contains_key(&atom) {
            return;
        }
        let text = atom_tag(&loaded.structure.topology, atom as usize);
        self.dispatch(Command::SetLabel {
            id,
            atom,
            text: Some(text.clone()),
        });
        self.undoable(&format!("Labeled {text}."));
    }

    /// Selects `atoms` of structure `id`. With `add`, they join the
    /// existing selection when it is on the same structure (e.g. picking a
    /// second atom for a distance measurement, see `inspector_ui`); a
    /// Ctrl+click that lands on a different structure starts fresh rather
    /// than silently doing nothing.
    pub(crate) fn select_atoms(&mut self, id: StructureId, atoms: &[u32], add: bool) {
        let Some(loaded) = self.scene.structure(id) else {
            return;
        };
        let n = loaded.structure.atom_count();
        let mut bits = match self.scene.active_selection() {
            Some(active) if add && active.structure == id => (*active.mask).clone(),
            _ => (*empty_mask(n)).clone(),
        };
        for &atom in atoms {
            bits.insert(atom as usize);
        }
        // Remember click order: a dihedral is A-B-C-D, not "these four".
        if !(add
            && self
                .scene
                .active_selection()
                .is_some_and(|a| a.structure == id))
        {
            self.pick_order.clear();
        }
        for &atom in atoms {
            if !self.pick_order.contains(&atom) {
                self.pick_order.push(atom);
            }
        }
        self.dispatch(Command::Select {
            id,
            mask: std::sync::Arc::new(bits),
        });
    }

    /// The selected atoms in the order they were clicked, when the
    /// selection is exactly the clicked atoms; otherwise (an expression, a
    /// saved set, an undo) in index order.
    fn ordered_selection(&self, mask: &vv_scene::Mask) -> Vec<usize> {
        let count = mask.count_ones(..);
        if self.pick_order.len() == count && self.pick_order.iter().all(|&a| mask[a as usize]) {
            self.pick_order.iter().map(|&a| a as usize).collect()
        } else {
            mask.ones().collect()
        }
    }
}

/// A usable path from user input, or `None` when it is blank.
fn clean_path(text: &str) -> Option<std::path::PathBuf> {
    let mut t = text.trim();
    for (open, close) in [('"', '"'), ('\'', '\''), ('<', '>')] {
        if t.len() >= 2 && t.starts_with(open) && t.ends_with(close) {
            t = t[1..t.len() - 1].trim();
        }
    }
    if t.is_empty() {
        return None;
    }
    let expanded = if let Some(rest) = t.strip_prefix('~') {
        let home = std::env::var_os("USERPROFILE").or_else(|| std::env::var_os("HOME"));
        match home {
            Some(h) => std::path::PathBuf::from(h).join(rest.trim_start_matches(['/', '\\'])),
            None => std::path::PathBuf::from(t),
        }
    } else {
        std::path::PathBuf::from(t)
    };
    Some(expanded)
}

impl TabViewerTrait for AppUi<'_> {
    type Tab = Tab;

    fn id(&mut self, tab: &mut Self::Tab) -> egui::Id {
        egui::Id::new(*tab as usize)
    }

    fn title(&mut self, tab: &mut Self::Tab) -> WidgetText {
        tab.title().into()
    }

    fn is_closeable(&self, tab: &Self::Tab) -> bool {
        *tab != Tab::Viewport
    }

    fn add_popup(&mut self, ui: &mut Ui, path: egui_dock::NodePath) {
        if self.closed_panels.is_empty() {
            widgets::caption(ui, "Every panel is open");
            return;
        }
        for &tab in &self.closed_panels {
            if ui.button(tab.title()).clicked() {
                *self.layout_request = Some(LayoutRequest::AddPanel(path, tab));
                ui.close();
            }
        }
    }

    fn ui(&mut self, ui: &mut Ui, tab: &mut Self::Tab) {
        if *tab == Tab::Viewport {
            self.viewport_ui(ui);
            return;
        }
        // These scroll themselves: the log (its console line pinned) and
        // the sequence (both ways).
        if matches!(*tab, Tab::Log | Tab::Sequence) {
            egui::Frame::NONE
                .inner_margin(crate::theme::space::INSET)
                .show(ui, |ui| self.panel_ui(ui, *tab));
            return;
        }
        egui::ScrollArea::vertical()
            .id_salt(("panel", *tab as usize))
            .auto_shrink([false, false])
            .show(ui, |ui| {
                egui::Frame::NONE
                    .inner_margin(crate::theme::space::GAP)
                    .show(ui, |ui| self.panel_ui(ui, *tab));
            });
    }
}

impl AppUi<'_> {
    fn panel_ui(&mut self, ui: &mut Ui, tab: Tab) {
        match tab {
            Tab::Viewport => unreachable!("the viewport is not a padded panel"),
            Tab::Structures => self.scene_ui(ui),
            Tab::Selection => self.selection_ui(ui),
            Tab::Inspector => self.inspector_ui(ui),
            Tab::Log => self.log_ui(ui),
            Tab::Sequence => self.sequence_ui(ui),
            // Unreachable in practice (`layout::strip_scene_settings`
            // drops it from every layout/session on load) but not
            // provably so from a hand-edited file, so this stays a
            // graceful message rather than a panic.
            Tab::SceneSettings => {
                widgets::caption(ui, "This panel is gone; see the Look tab.");
            }
            Tab::Timeline => self.timeline_ui(ui),
            Tab::Movie => self.movie_ui(ui),
            Tab::Info => self.info_ui(ui),
        }
    }
}

/// An atom's compact tag for a live measurement, a label's default text,
/// or an Annotations row: "CA·ALA12" (name, residue name, residue
/// number). Shared so the three read the same way.
fn atom_tag(top: &vv_core::Topology, atom: usize) -> String {
    let r = top.residue_index[atom] as usize;
    format!(
        "{}\u{00b7}{}{}",
        top.atom_name(atom),
        top.residue_name(r),
        top.residues[r].auth_seq_id
    )
}

/// One option of `rep_options_ui`: a slider staged in `ui`'s temp storage
/// while dragged and pushed to `set` on release (one undo step), or
/// `radius_by`'s segmented Constant/B-factor choice.
fn rep_option_ui(
    ui: &mut Ui,
    id: StructureId,
    rep: &vv_scene::Rep,
    o: &vv_scene::scene::RepOption,
    set: &mut Vec<(&'static str, Option<f32>)>,
) {
    if o.name == "radius_by" {
        let now = rep.option(o.name).unwrap_or(o.default) >= 0.5;
        ui.label(o.label);
        if let Some(i) = widgets::segmented(ui, &["Constant", "B-factor"], Some(now as usize)) {
            set.push((o.name, Some(i as f32)));
        }
        return;
    }
    let key = egui::Id::new(("rep option", id, rep.id, o.name));
    let mut value = ui
        .data_mut(|d| d.get_temp::<f32>(key))
        .unwrap_or_else(|| rep.option(o.name).unwrap_or(o.default));
    let slider = widgets::slider(ui, o.label, &mut value, o.min..=o.max, o.unit);
    if slider.dragged() {
        ui.data_mut(|d| d.insert_temp(key, value));
    } else if slider.drag_stopped() || slider.changed() {
        ui.data_mut(|d| d.remove::<f32>(key));
        set.push((o.name, Some(value)));
    }
}

/// What a selection row's ⋯ menu or eye asked for.
enum RepPick {
    Style(vv_scene::Representation),
    Coloring(vv_scene::ColorScheme),
    Material(vv_scene::Material),
    ToggleVisible,
    SelectAtoms,
    EditExpression,
    Duplicate,
    Delete,
    Options,
}

/// A selection row's ⋯ menu: Style, Coloring and Material as submenus,
/// then the row's own actions.
fn rep_menu_ui(ui: &mut Ui, rep: &vv_scene::Rep, can_delete: bool) -> Option<RepPick> {
    use crate::ribbon::{MATERIALS, REPS};
    let mut pick = None;
    ui.menu_button("Style", |ui| {
        let now = format!(
            "rep {}",
            vv_scene::session::representation_name(rep.representation)
        );
        pick = choice_list(ui, REPS, &now)
            .and_then(|c| vv_scene::script::parse_representation(c.strip_prefix("rep ")?).ok())
            .map(RepPick::Style);
    });
    let glycan = rep.representation == vv_scene::Representation::Glycan;
    ui.add_enabled_ui(!glycan, |ui| {
        ui.menu_button("Coloring", |ui| {
            if let Some(coloring) = coloring_menu_ui(ui, &rep.coloring) {
                pick = Some(RepPick::Coloring(coloring));
            }
        })
    })
    .response
    .on_disabled_hover_text("SNFG colours are fixed");
    ui.menu_button("Material", |ui| {
        let now = format!("material {}", rep.material.name());
        pick = choice_grid(ui, MATERIALS, &now, 3)
            .and_then(|c| vv_scene::Material::parse(c.strip_prefix("material ")?))
            .map(RepPick::Material);
    });
    ui.separator();
    for (label, action) in [
        ("Select these atoms", RepPick::SelectAtoms),
        ("Edit expression", RepPick::EditExpression),
        ("Duplicate", RepPick::Duplicate),
        ("Options…", RepPick::Options),
    ] {
        if ui.button(label).clicked() {
            pick = Some(action);
        }
    }
    let danger = crate::theme::Tokens::current(ui.ctx()).danger;
    if can_delete && ui.button(RichText::new("Delete").color(danger)).clicked() {
        pick = Some(RepPick::Delete);
    }
    pick
}

/// `options` (a ribbon list, e.g. `ribbon::REPS`) as checkable items, the
/// one whose command is `now` checked; the picked item's command.
fn choice_list(
    ui: &mut Ui,
    options: &[(&'static str, &'static str)],
    now: &str,
) -> Option<&'static str> {
    let mut picked = None;
    for (label, command) in options {
        if ui.selectable_label(*command == now, *label).clicked() {
            picked = Some(*command);
        }
    }
    picked
}

/// `choice_list` laid out `columns` wide, filled row by row.
fn choice_grid(
    ui: &mut Ui,
    options: &[(&'static str, &'static str)],
    now: &str,
    columns: usize,
) -> Option<&'static str> {
    let mut picked = None;
    // Keyed by the slice, so several grids can share one menu.
    egui::Grid::new(ui.id().with(options.as_ptr() as usize))
        .spacing(egui::vec2(
            crate::theme::space::TIGHT,
            crate::theme::space::TIGHT,
        ))
        .show(ui, |ui| {
            for (i, (label, command)) in options.iter().enumerate() {
                if ui.selectable_label(*command == now, *label).clicked() {
                    picked = Some(*command);
                }
                if i % columns == columns - 1 {
                    ui.end_row();
                }
            }
        });
    picked
}

/// The schemes grouped under a caption per `ribbon::COLORING_GROUPS`
/// (Atoms, Chains, Residues, Hydrophobicity), then a swatch per named
/// colour for one solid colour (any other colour: Options… once it is
/// solid).
fn coloring_menu_ui(
    ui: &mut Ui,
    coloring: &vv_scene::ColorScheme,
) -> Option<vv_scene::ColorScheme> {
    let now = format!("color {}", coloring.name());
    use crate::ribbon::{COLORINGS, COLORING_GROUPS};
    let mut picked = None;
    for (g, (heading, start)) in COLORING_GROUPS.iter().enumerate() {
        let end = COLORING_GROUPS
            .get(g + 1)
            .map_or(COLORINGS.len(), |(_, s)| *s);
        if g > 0 {
            ui.separator();
        }
        widgets::caption(ui, heading);
        if let Some(command) = choice_grid(ui, &COLORINGS[*start..end], &now, 3) {
            picked = command
                .strip_prefix("color ")
                .and_then(vv_scene::ColorScheme::parse);
        }
    }
    ui.separator();
    widgets::caption(ui, "Solid color");
    egui::Grid::new("solid-color-swatches")
        .spacing(egui::vec2(4.0, 4.0))
        .show(ui, |ui| {
            for (i, (name, rgb)) in vv_scene::scene::NAMED_COLORS.iter().enumerate() {
                let current = *coloring == vv_scene::ColorScheme::Constant(*rgb);
                if swatch(ui, *rgb, current).on_hover_text(*name).clicked() {
                    picked = Some(vv_scene::ColorScheme::Constant(*rgb));
                }
                if i % 6 == 5 {
                    ui.end_row();
                }
            }
        });
    picked
}

fn swatch(ui: &mut Ui, [r, g, b]: [u8; 3], current: bool) -> egui::Response {
    let (rect, response) = ui.allocate_exact_size(egui::vec2(20.0, 20.0), egui::Sense::click());
    let t = crate::theme::Tokens::current(ui.ctx());
    let stroke = if current || response.hovered() {
        egui::Stroke::new(2.0, t.primary)
    } else {
        egui::Stroke::new(1.0, t.border)
    };
    ui.painter().rect(
        rect,
        crate::theme::radius::CONTROL,
        Color32::from_rgb(r, g, b),
        stroke,
        egui::StrokeKind::Inside,
    );
    response
}

/// A solid colour's full picker, staged while the pointer is held and
/// returned on release (one undo step per pick, not one per frame).
fn solid_color_ui(ui: &mut Ui, key: egui::Id, rgb: [u8; 3]) -> Option<[u8; 3]> {
    let staged = ui.data_mut(|d| d.get_temp::<Color32>(key));
    let mut color = staged.unwrap_or(Color32::from_rgb(rgb[0], rgb[1], rgb[2]));
    let changed = widgets::color_editor(ui, key.with("editor"), &mut color);
    if ui.input(|i| i.pointer.any_down()) {
        if changed || staged.is_some() {
            ui.data_mut(|d| d.insert_temp(key, color));
        }
        return None;
    }
    ui.data_mut(|d| d.remove::<Color32>(key));
    (changed || staged.is_some()).then(|| [color.r(), color.g(), color.b()])
}

impl AppUi<'_> {
    /// The view cube and scale bar over the viewport image at `rect`, or,
    /// with nothing loaded and not dismissed this run, the start card.
    fn overlays_ui(&mut self, ui: &mut Ui, rect: egui::Rect) {
        let tokens = crate::theme::Tokens::of(self.prefs.theme);
        if self.scene.structures().next().is_none() && !*self.start_card_dismissed {
            self.start_card_ui(ui, rect, &tokens);
            return;
        }
        if let Some(side) = crate::overlays::view_cube(ui, rect, self.camera.orientation, &tokens) {
            if let Some(anim) = CubeAnim::to_face(self.camera, side) {
                self.mark_camera_jump();
                *self.cube_anim = Some(anim);
            }
        }
        // World height of the view at the target's depth, either projection.
        let half_height = self.camera.distance * (self.camera.fov_y * 0.5).tan();
        let b = self.view.background;
        let background = egui::Rgba::from_rgb(b.r as f32, b.g as f32, b.b as f32);
        crate::overlays::scale_bar(
            ui,
            rect,
            2.0 * half_height / rect.height().max(1.0),
            background.into(),
        );
        self.render_frame_overlay(ui, rect);
    }

    /// Studio ▸ Frame preview / `framing`, or File ▸ Export ▸ Render's
    /// popover being open: a safe-frame rectangle at the render's W:H
    /// aspect, matching what `render_form_size` reads off its own form
    /// (or the viewport's own size, before the form has ever been
    /// opened -- the overlay then just fills the viewport).
    fn render_frame_overlay(&self, ui: &Ui, rect: egui::Rect) {
        let render_popover_open =
            egui::Popup::is_id_open(ui.ctx(), egui::Id::new(("ribbon-popover", "file.render")));
        if !render_popover_open && !*self.studio_frame {
            return;
        }
        let (w, h) = render_form_size(ui.ctx(), self.renderer);
        crate::overlays::render_frame(ui, rect, w as f32 / h.max(1) as f32);
    }

    /// The ways to start: the logo, Open/Fetch/Open session, up to 4 most
    /// recent files, and "Don't show again". Closable (`welcome` brings it
    /// back for the rest of the run); its own checkbox turns off the
    /// automatic showing at launch (`docs/COMMANDS.md`'s `welcome`).
    fn start_card_ui(&mut self, ui: &mut Ui, rect: egui::Rect, tokens: &crate::theme::Tokens) {
        use crate::theme::space;
        use egui_phosphor::regular as icon;
        let recents: Vec<std::path::PathBuf> =
            self.prefs.recent_files.iter().take(4).cloned().collect();
        let logo = crate::start::logo_texture(ui.ctx(), self.prefs.theme);
        let desired = 320.0
            + if logo.is_some() { 135.0 } else { 0.0 }
            + if recents.is_empty() {
                0.0
            } else {
                36.0 + 32.0 * recents.len() as f32
            };
        // Never taller than the viewport itself: a small docked panel or a
        // long recents list scrolls (one `ScrollArea` around everything,
        // including the header) instead of drawing outside the card.
        let height = desired.min((rect.height() - space::WIDE * 2.0).max(240.0));
        let card = egui::Rect::from_center_size(rect.center(), egui::vec2(380.0, height));
        let mut child = ui.new_child(egui::UiBuilder::new().max_rect(card));
        let mut run = None;
        let mut close = false;
        egui::Frame::NONE
            .fill(tokens.surface)
            .stroke(egui::Stroke::new(1.0, tokens.border))
            .corner_radius(crate::theme::radius::CONTAINER)
            .inner_margin(space::WIDE)
            .show(&mut child, |ui| {
                ui.set_width(card.width() - 2.0 * space::WIDE);
                egui::ScrollArea::vertical()
                    .max_height(height - 2.0 * space::WIDE)
                    .show(ui, |ui| {
                        ui.horizontal(|ui| {
                            if let Some(logo) = &logo {
                                let width = 160.0_f32.min(logo.size_vec2().x);
                                let size = logo.size_vec2() * (width / logo.size_vec2().x);
                                ui.vertical_centered(|ui| ui.image((logo.id(), size)));
                            }
                            ui.with_layout(egui::Layout::right_to_left(egui::Align::TOP), |ui| {
                                if widgets::button(ui, icon::X, "", Variant::Ghost)
                                    .on_hover_text("Close")
                                    .clicked()
                                {
                                    close = true;
                                }
                            });
                        });
                        ui.label(
                            RichText::new("Open a structure")
                                .text_style(crate::theme::title_style()),
                        );
                        widgets::caption(ui, "PDB or mmCIF — or drop a file here.");
                        ui.add_space(space::INSET);
                        for (glyph, label, command) in [
                            (icon::FOLDER_OPEN, "Open a file…  (Ctrl+O)", "open"),
                            (icon::CLOUD_ARROW_DOWN, "Fetch from the PDB…", "fetch"),
                            (icon::FILE_ARROW_UP, "Open a session…", "loadsession"),
                        ] {
                            if widgets::button(ui, glyph, label, Variant::Ghost).clicked() {
                                run = Some(command.to_string());
                            }
                        }
                        if !recents.is_empty() {
                            ui.add_space(space::INSET);
                            widgets::section(ui, "Recent");
                            for path in &recents {
                                let name = path
                                    .file_name()
                                    .map(|n| n.to_string_lossy().into_owned())
                                    .unwrap_or_else(|| path.display().to_string());
                                let recent_verb = if path.extension().is_some_and(|e| {
                                    e.eq_ignore_ascii_case("vviz") || e.eq_ignore_ascii_case("json")
                                }) {
                                    "loadsession"
                                } else {
                                    "open"
                                };
                                if widgets::list_row(ui, &name, false, |_| {})
                                    .on_hover_text(path.display().to_string())
                                    .clicked()
                                {
                                    run = Some(format!("{recent_verb} {}", path.display()));
                                }
                            }
                        }
                        ui.add_space(space::INSET);
                        widgets::caption(
                            ui,
                            "Ctrl+P finds any command. Tap Alt (or F10) for key tips.",
                        );
                        ui.add_space(space::INSET);
                        if widgets::switch(ui, self.prefs.hide_start_card, "Don't show again")
                            .clicked()
                        {
                            self.prefs.hide_start_card = !self.prefs.hide_start_card;
                            close = self.prefs.hide_start_card;
                        }
                    });
            });
        if let Some(command) = run {
            self.run_command_logged(&command);
        }
        if close {
            *self.start_card_dismissed = true;
        }
    }

    /// What the mode chip and the mouse popup say a left drag does.
    fn left_drag_label(&self) -> String {
        if self.view.select_tool.draws() && self.view.mouse_mode == MouseMode::Rotate {
            self.view.select_tool.drag_label()
        } else {
            self.view.mouse_mode.label().to_owned()
        }
    }

    /// Collects the left drag's path in viewport-local pixels, paints it,
    /// and on release selects inside it: Shift adds, Ctrl/Alt subtracts.
    fn select_drag_ui(&mut self, ui: &mut Ui, response: &egui::Response) {
        let path_id = response.id.with("select-path");
        let mut path: Vec<egui::Pos2> = ui.data_mut(|d| d.get_temp(path_id)).unwrap_or_default();
        if let Some(pos) = response.interact_pointer_pos() {
            let local = (pos - response.rect.min).to_pos2();
            if path.last().is_none_or(|last| last.distance(local) >= 2.0) {
                path.push(local);
            }
        }
        let shape = self.view.select_tool.shape;
        if let Some(region) = crate::select_tool::Region::from_path(shape, &path) {
            let primary = crate::theme::Tokens::of(self.prefs.theme).primary;
            let painter = ui.painter_at(response.rect);
            crate::select_tool::paint_region(&painter, response.rect.min, &region, primary);
        }
        ui.data_mut(|d| d.insert_temp(path_id, path));
        ui.ctx().request_repaint();
    }

    /// Runs the drag the release just ended, if any.
    fn select_drag_release(&mut self, ui: &Ui, response: &egui::Response) {
        let path_id = response.id.with("select-path");
        let path: Vec<egui::Pos2> = ui.data_mut(|d| d.remove_temp(path_id)).unwrap_or_default();
        let shape = self.view.select_tool.shape;
        let Some(region) = crate::select_tool::Region::from_path(shape, &path) else {
            return;
        };
        use crate::select_tool::Combine;
        let how = ui.input(
            |i| match (i.modifiers.shift, i.modifiers.command || i.modifiers.alt) {
                (_, true) => Combine::Subtract,
                (true, false) => Combine::Add,
                (false, false) => Combine::Replace,
            },
        );
        let size = (response.rect.width(), response.rect.height());
        self.select_in_region(&region, size, how);
    }

    fn viewport_ui(&mut self, ui: &mut Ui) {
        let available = ui.available_size();
        *self.requested_viewport_size = Some(available);

        let Some(texture) = self.viewport_texture else {
            ui.centered_and_justified(|ui| ui.label("no viewport texture yet"));
            return;
        };
        let (rw, rh) = self.renderer.size();
        let aspect = rw as f32 / rh.max(1) as f32;
        let image = egui::Image::new((texture, available)).sense(egui::Sense::click_and_drag());
        let response = ui.add(image);
        self.overlays_ui(ui, response.rect);
        if let Some(progress) = self.render_progress {
            render_progress_badge(ui, response.rect, progress);
        }

        let selecting = self.view.select_tool.draws() && self.view.mouse_mode == MouseMode::Rotate;
        if selecting && response.dragged_by(egui::PointerButton::Primary) {
            self.select_drag_ui(ui, &response);
        } else if response.dragged_by(egui::PointerButton::Secondary) {
            let delta = response.drag_delta();
            self.camera.pan(delta.x, delta.y);
        } else if response.dragged_by(egui::PointerButton::Primary) {
            let delta = response.drag_delta();
            match self.view.mouse_mode {
                MouseMode::Translate => self.camera.pan(delta.x, delta.y),
                // Drag up to zoom in, same rate as the wheel.
                MouseMode::Scale => self.camera.zoom((delta.y * 0.004).exp()),
                // Label mode still rotates on drag; only a click labels.
                MouseMode::Rotate | MouseMode::Center | MouseMode::Label => {
                    trackball_rotate(self.camera, delta)
                }
            }
        } else if response.dragged_by(egui::PointerButton::Middle) {
            if let Some(pos) = response.interact_pointer_pos() {
                let center = response.rect.center();
                let angle = roll_angle(pos - response.drag_delta() - center, pos - center);
                self.camera.rotate_scene(glam::Vec3::Z, angle);
            }
        }
        if selecting && response.drag_stopped_by(egui::PointerButton::Primary) {
            self.select_drag_release(ui, &response);
        }
        ui.input(|i| {
            // Exponential: every wheel notch (or trackpad distance) is the
            // same fraction closer or farther, at any zoom level.
            let scroll = i.smooth_scroll_delta.y;
            if i.modifiers.shift && response.hovered() {
                // egui turns Shift+wheel into horizontal scroll.
                self.camera.dolly_wheel(scroll + i.smooth_scroll_delta.x);
            } else if scroll != 0.0 && response.hovered() {
                self.camera.zoom((-scroll * 0.0025).exp());
            }
        });

        // Ctrl (Cmd on Mac) held: releasing it or Escape cancels the
        // pending measure chain (not any measurement already shown).
        let ctrl = ui.input(|i| i.modifiers.command);
        if self.measure.is_pending() && (!ctrl || ui.input(|i| i.key_pressed(egui::Key::Escape))) {
            self.measure.cancel();
        }

        // Double-click an atom: center on it, select its residue, and
        // zoom to its ~5 Å neighborhood.
        if response.double_clicked() {
            if let Some(pos) = response.interact_pointer_pos() {
                let (px, py) = viewport::pointer_pixel(pos, response.rect, rw, rh);
                if let Some(Pick::Atom { item, atom }) = self.renderer.pick(px, py) {
                    if let Some(source) = self.draw_order.get(item) {
                        let (id, atom) = (source.id, source.atom(atom));
                        self.zoom_to_residue(id, atom);
                    }
                }
            }
        }

        // A single atom/bond click: Center/Label act specially; Ctrl adds
        // to the selection and extends the measure chain, a plain click
        // replaces it. Guarded against the double-click's own second
        // release, which egui also reports as `clicked()`.
        if response.clicked() && !response.double_clicked() {
            if let Some(pos) = response.interact_pointer_pos() {
                let (px, py) = viewport::pointer_pixel(pos, response.rect, rw, rh);
                let add = ctrl;
                let pick = self.renderer.pick(px, py);
                // `item` indexes the draw list of the frame just rendered,
                // which `draw_order` mirrors (see `State::render_scene`).
                match pick {
                    Some(Pick::Atom { item, atom })
                        if self.view.mouse_mode == MouseMode::Center =>
                    {
                        if let Some(source) = self.draw_order.get(item) {
                            let (id, atom) = (source.id, source.atom(atom));
                            self.center_on_atom(id, atom);
                            // Back to rotating once centered.
                            self.view.mouse_mode = MouseMode::Rotate;
                        }
                    }
                    Some(Pick::Atom { item, atom }) if self.view.mouse_mode == MouseMode::Label => {
                        // Sticky, unlike Center: stays in label mode so
                        // several atoms can be labeled in a row.
                        if let Some(source) = self.draw_order.get(item) {
                            let (id, atom) = (source.id, source.atom(atom));
                            self.label_atom(id, atom);
                        }
                    }
                    Some(Pick::Atom { item, atom }) => {
                        if let Some(source) = self.draw_order.get(item) {
                            let (id, atom) = (source.id, source.atom(atom));
                            self.select_atoms(id, &[atom], add);
                            if add {
                                self.push_measure(id, atom);
                            }
                        }
                    }
                    Some(Pick::Bond { item, bond }) => {
                        // A bond selects both of its atoms; the inspector
                        // then shows its length. Derived geometry (tube
                        // segments, ligand sticks) carries its own pairs.
                        // Kept out of the measure chain (which is atom by
                        // atom): a bond click only ever extends selection.
                        let pair = self.draw_order.get(item).and_then(|source| {
                            let pair = match &source.bond_atoms {
                                Some(pairs) => *pairs.get(bond as usize)?,
                                None => *self
                                    .scene
                                    .structure(source.id)?
                                    .bonds
                                    .pairs
                                    .get(bond as usize)?,
                            };
                            Some((source.id, pair))
                        });
                        if let Some((id, pair)) = pair {
                            self.select_atoms(id, &pair, add);
                        }
                    }
                    None => {
                        if !add {
                            self.dispatch(Command::ClearSelection);
                        }
                    }
                }
            }
        }

        // The mouse mode (what a left drag does); frame timings on hover.
        let stats = format!(
            "{:.0} fps ({:.1} ms): cull {:.2}, draw {:.2}, depth pyramid {:.2} ms\n\
             atoms under {:.1} px draw as points{}{}",
            self.readout.fps,
            self.readout.frame_ms,
            self.readout.cull_ms,
            self.readout.draw_ms,
            self.readout.hiz_ms,
            self.readout.quad_px_threshold,
            if self.view.adaptive {
                " (adaptive)"
            } else {
                ""
            },
            if self.view.occlusion_culling {
                ""
            } else {
                "; occlusion culling off"
            },
        );
        let mode = format!("{} · keys R T S C", self.left_drag_label());
        let (chip, help) = widgets::mode_chip(ui, response.rect.min, &mode);
        chip.on_hover_text(stats);
        let left_drag = self.left_drag_label();
        egui::Popup::from_toggle_button_response(&help)
            .close_behavior(egui::PopupCloseBehavior::CloseOnClickOutside)
            .show(|ui| mouse_bindings_popup(ui, left_drag));

        self.update_hover(ui, &response, rw, rh);

        let right_click = response
            .secondary_clicked()
            .then(|| response.interact_pointer_pos())
            .flatten();
        if let Some((id, atom)) =
            self.labels_and_captions_ui(ui, response.rect, aspect, right_click)
        {
            self.dispatch(Command::SetLabel {
                id,
                atom,
                text: None,
            });
            self.undoable("Removed the label.");
        }

        if ctrl {
            self.draw_pending_measure(ui, response.rect, aspect);
        }

        if let Some(Notice {
            text: message,
            when,
            error,
            action,
        }) = self.notice.clone()
        {
            // Info fades after a few seconds; an error stays until dismissed.
            let expired = !error && when.elapsed().as_secs_f32() > 4.0;
            // Bottom-centered and compact, not a full-width bar; Item 5a:
            // leave the scale bar's own corner clear so both show.
            let width = (response.rect.width() - 2.0 * crate::theme::space::GAP).min(420.0);
            let cx = response.rect.center().x;
            let toast_rect = egui::Rect::from_min_max(
                egui::pos2(cx - width * 0.5, response.rect.min.y),
                egui::pos2(
                    cx + width * 0.5,
                    response.rect.max.y - crate::overlays::SCALE_BAR_CLEARANCE,
                ),
            );
            match widgets::toast(ui, toast_rect, &message, error, action) {
                widgets::ToastOutcome::Dismissed => *self.notice = None,
                widgets::ToastOutcome::Action(command) => {
                    *self.notice = None;
                    self.run_command_logged(command);
                }
                widgets::ToastOutcome::None if expired => *self.notice = None,
                widgets::ToastOutcome::None => {
                    if !error {
                        ui.ctx()
                            .request_repaint_after(std::time::Duration::from_millis(250));
                    }
                }
            }
        }
        if self.drop_hover {
            widgets::drop_target(
                ui,
                response.rect,
                "Drop to open (.cif, .pdb, .gz, or a session .json)",
            );
        }
    }

    /// Hover glow (no button held): re-picks only when the pointer
    /// actually moved (`HoverPick`, a blocking GPU readback otherwise
    /// re-run on every incidental repaint), and nudges one extra repaint
    /// when the hovered atom changes -- `render_scene` (state.rs) reads
    /// `self.hover` a frame later, like `row_highlight`, so without this
    /// the glow would always lag one pointer move behind.
    fn update_hover(&mut self, ui: &Ui, response: &egui::Response, rw: u32, rh: u32) {
        let pixel = (response.hovered() && !ui.input(|i| i.pointer.any_down()))
            .then(|| response.hover_pos())
            .flatten()
            .map(|pos| viewport::pointer_pixel(pos, response.rect, rw, rh));
        if !self.hover.moved(pixel) {
            return;
        }
        let picked = pixel.and_then(|(x, y)| match self.renderer.pick(x, y) {
            Some(Pick::Atom { item, atom }) => {
                self.draw_order.get(item).map(|s| (s.id, s.atom(atom)))
            }
            _ => None,
        });
        if picked != self.hover.atom {
            ui.ctx().request_repaint();
        }
        self.hover.atom = picked;
    }

    /// Adds `atom` to the pending Ctrl+click measure chain, dispatching
    /// whatever `MeasureChain::push` says to show (one undo step).
    fn push_measure(&mut self, id: StructureId, atom: u32) {
        let measure::ChainStep::Show { hide, show } = self.measure.push(id, atom) else {
            return;
        };
        let mut commands: Vec<Command> = hide
            .into_iter()
            .map(|measurement| Command::SetMeasurement {
                id,
                measurement,
                shown: false,
            })
            .collect();
        commands.push(Command::SetMeasurement {
            id,
            measurement: show,
            shown: true,
        });
        self.dispatch(Command::Batch(commands));
        self.undoable("Measured.");
    }

    /// The pending chain's live preview while Ctrl is held: a ring on its
    /// anchor atom and a rubber band from its last atom to the pointer,
    /// with the live distance to whatever's hovered.
    fn draw_pending_measure(&self, ui: &Ui, rect: egui::Rect, aspect: f32) {
        let (Some((id, anchor)), Some((_, tip)), Some(pointer)) = (
            self.measure.anchor(),
            self.measure.tip(),
            ui.input(|i| i.pointer.hover_pos()),
        ) else {
            return;
        };
        let Some(loaded) = self.scene.structure(id) else {
            return;
        };
        let frame = loaded
            .frame
            .min(loaded.structure.frame_count().saturating_sub(1));
        let coords = loaded.structure.frame(frame);
        let positions = coords.positions();
        let view_proj = self.camera.proj(aspect) * self.camera.view();
        let project = |atom: u32| {
            positions.get(atom as usize).and_then(|&world| {
                project_to_pixel(view_proj, rect.width(), rect.height(), world)
                    .map(|(x, y, _)| rect.min + egui::vec2(x, y))
            })
        };
        let (Some(anchor_2d), Some(tip_2d)) = (project(anchor), project(tip)) else {
            return;
        };
        let hover_atom = self
            .hover
            .atom
            .filter(|&(hid, _)| hid == id)
            .map(|(_, a)| a);
        let painter = ui.painter_at(rect);
        measure::draw_pending(
            &painter,
            anchor_2d,
            tip_2d,
            pointer,
            hover_atom.and_then(project),
        );
        if let Some(hover_atom) = hover_atom {
            if let (Some(&tip_pos), Some(&hover_pos)) = (
                positions.get(tip as usize),
                positions.get(hover_atom as usize),
            ) {
                painter.text(
                    pointer + egui::vec2(10.0, -10.0),
                    egui::Align2::LEFT_BOTTOM,
                    format!("{:.2} \u{c5}", tip_pos.distance(hover_pos)),
                    egui::TextStyle::Small.resolve(ui.style()),
                    Color32::from_rgb(
                        LABEL_TEXT_COLOR[0],
                        LABEL_TEXT_COLOR[1],
                        LABEL_TEXT_COLOR[2],
                    ),
                );
            }
        }
    }

    /// Double-click: centers on `atom`, selects its whole
    /// residue, frames it plus `NEIGHBORHOOD_RADIUS` Å around it
    /// (animated), and adds a ball-and-stick layer for that neighborhood.
    /// The scene edit (select + add rep) is one undo step; the camera
    /// jump is a separate one (different undo domains -- see
    /// `mark_camera_jump`'s doc).
    pub(crate) fn zoom_to_residue(&mut self, id: StructureId, atom: u32) {
        const NEIGHBORHOOD_RADIUS: f32 = 5.0;
        let Some(loaded) = self.scene.structure(id) else {
            return;
        };
        let top = loaded.structure.topology.clone();
        let frame = loaded
            .frame
            .min(loaded.structure.frame_count().saturating_sub(1));
        let coords = loaded.structure.frame(frame);
        let current_rep = loaded.rep().clone();
        let next_rep_id = loaded.next_rep_id;
        let rep_count = loaded.reps.len();
        let previous: Vec<vv_scene::RepId> = loaded
            .reps
            .iter()
            .filter(|r| viewport::is_neighborhood_expr(&r.selection))
            .map(|r| r.id)
            .collect();

        let (center, radius) =
            viewport::bounding_sphere(coords.positions(), viewport::residue_atoms(&top, atom));
        let residue_expr = viewport::residue_expr(&top, atom);
        let neighborhood = viewport::neighborhood_expr(&top, atom, NEIGHBORHOOD_RADIUS);
        if let Err(e) = vv_core::select::parse(&neighborhood) {
            self.fail(format!("neighborhood selection `{neighborhood}`: {e}"));
            return;
        }

        self.mark_camera_jump();
        self.center_on_atom(id, atom);
        let mut target = self.camera.clone();
        target.target = center;
        target.distance = Camera::framing(center, radius + NEIGHBORHOOD_RADIUS).distance;
        *self.cube_anim = Some(CubeAnim::to_camera(self.camera, &target));

        let new_rep = vv_scene::Rep {
            id: vv_scene::RepId(next_rep_id),
            selection: neighborhood,
            representation: vv_scene::Representation::BallAndStick,
            ..current_rep
        };
        self.pick_order.clear();
        // Added before the old ones go, so a structure never drops to
        // zero reps mid-batch (`SceneError::LastRep`).
        let mut commands = vec![
            Command::SelectExpr {
                id,
                expr: residue_expr,
            },
            Command::AddRep {
                id,
                index: rep_count,
                rep: new_rep,
            },
        ];
        commands.extend(
            previous
                .into_iter()
                .map(|rep| Command::RemoveRep { id, rep }),
        );
        self.dispatch(Command::Batch(commands));
        self.undoable("Added the residue's neighborhood.");
    }

    /// Brings `selection_summary` up to the active selection.
    fn refresh_selection_summary(&mut self) {
        let Some(active) = self.scene.active_selection() else {
            return;
        };
        let Some(loaded) = self.scene.structure(active.structure) else {
            return;
        };
        let key = std::sync::Arc::as_ptr(&active.mask) as usize;
        if self.selection_summary.key != Some(key) {
            *self.selection_summary = SelectionSummary::of(&active.mask, loaded, key);
        }
    }

    /// Structure-anchored labels (billboarded at each labeled atom's
    /// current-frame world position) and screen-space captions
    /// (fractional viewport position, independent of the camera) --
    /// docs/UI_DESIGN.md's two annotation types. Painted last so they sit
    /// above the rendered scene but below the notice bar/drop overlay.
    /// `right_click`, when a right-click landed this frame, is hit-tested
    /// against each label's dot ("right-click a label to remove it");
    /// the caller dispatches the removal, since this only borrows `&self`.
    fn labels_and_captions_ui(
        &self,
        ui: &Ui,
        rect: egui::Rect,
        aspect: f32,
        right_click: Option<egui::Pos2>,
    ) -> Option<(StructureId, u32)> {
        const HIT_RADIUS: f32 = 8.0;
        let painter = ui.painter_at(rect);
        let view_proj = self.camera.proj(aspect) * self.camera.view();
        let mut removed = None;
        for draw in label_draws(self.scene) {
            let Some((dx, dy, _depth)) =
                project_to_pixel(view_proj, rect.width(), rect.height(), draw.world)
            else {
                continue;
            };
            let anchor = rect.min + egui::vec2(dx, dy);
            if right_click.is_some_and(|pos| pos.distance(anchor) <= HIT_RADIUS) {
                removed = Some((draw.id, draw.atom));
            }
            painter.circle_filled(
                anchor,
                LABEL_DOT_RADIUS,
                Color32::from_rgb(LABEL_DOT_COLOR[0], LABEL_DOT_COLOR[1], LABEL_DOT_COLOR[2]),
            );
            annotation_tag(ui, &painter, anchor + egui::vec2(6.0, -6.0), &draw.text);
        }
        let line_color = Color32::from_rgb(
            MEASURE_LINE_COLOR[0],
            MEASURE_LINE_COLOR[1],
            MEASURE_LINE_COLOR[2],
        );
        for m in measurement_draws(self.scene) {
            // Unbounded: zoomed in, an atom just off screen must not hide
            // the whole measurement; `painter` clips to the viewport.
            let project = |w| {
                in_front_pixel(view_proj, rect.width(), rect.height(), w)
                    .map(|(x, y)| rect.min + egui::vec2(x, y))
            };
            let Some(points) = m
                .path
                .iter()
                .map(|&w| project(w))
                .collect::<Option<Vec<_>>>()
            else {
                continue;
            };
            painter.extend(egui::Shape::dashed_line(
                &points,
                egui::Stroke::new(1.5, line_color),
                MEASURE_DASH[0],
                MEASURE_DASH[1],
            ));
            if let Some(at) = project(m.anchor) {
                // Kept on screen while any part of the measurement is.
                let at = at.clamp(
                    rect.min + egui::vec2(4.0, 20.0),
                    rect.max - egui::vec2(80.0, 4.0),
                );
                annotation_tag(ui, &painter, at + egui::vec2(4.0, -4.0), &m.text);
            }
        }
        for caption in self.scene.captions() {
            let pos = egui::pos2(
                rect.min.x + caption.x.clamp(0.0, 1.0) * rect.width(),
                rect.min.y + caption.y.clamp(0.0, 1.0) * rect.height(),
            );
            painter.text(
                pos,
                egui::Align2::LEFT_TOP,
                &caption.text,
                egui::TextStyle::Heading.resolve(ui.style()),
                Color32::from_rgb(
                    CAPTION_TEXT_COLOR[0],
                    CAPTION_TEXT_COLOR[1],
                    CAPTION_TEXT_COLOR[2],
                ),
            );
        }
        removed
    }

    /// What is loaded, one compact row each ("Structures", formerly
    /// "Scene"): name, atom count, an eye for the whole structure,
    /// ⋯ (rename, close, file info, export). Clicking a row makes it
    /// current -- the Represent ribbon tab follows (`AppUi::current`'s
    /// doc) -- and shows its selections below.
    fn scene_ui(&mut self, ui: &mut Ui) {
        use egui_phosphor::regular as icon;
        let current = self.current();
        let rows: Vec<(StructureId, String, usize, bool)> = self
            .scene
            .structures()
            .map(|(id, s)| (id, s.label.clone(), s.structure.atom_count(), s.visible))
            .collect();
        if rows.is_empty() {
            if widgets::empty_state(ui, icon::ATOM, "No structures yet", "Open a structure") {
                self.open_file_dialog();
            }
            return;
        }
        let key = structure_rename_key();
        let renaming = renaming_id::<StructureId>(ui, key);
        for (id, label, atoms, visible) in rows {
            if renaming == Some(id) {
                if let Some(new_label) = inline_rename_ui(ui, key, id, &label) {
                    self.dispatch(Command::SetStructureLabel {
                        id,
                        label: new_label,
                    });
                }
                continue;
            }
            let row = format!("{label}  ·  {}", count(atoms, "atom"));
            let mut clicked = None;
            let mut toggle = false;
            let pick = widgets::list_row(ui, &row, Some(id) == current, |ui| {
                // Trailing icons draw right-to-left: the menu first puts
                // it rightmost, so the eye (left of it) reads "eye, ⋯".
                let items = [
                    ("Rename", false),
                    ("File info", false),
                    ("Export structure…", false),
                    ("Close", true),
                ];
                clicked = widgets::menu_button(ui, &items).0.map(|i| items[i].0);
                let eye = if visible { icon::EYE } else { icon::EYE_SLASH };
                toggle = widgets::button(ui, eye, "", Variant::Ghost)
                    .on_hover_text(if visible {
                        "Hide this structure (skipped while playing)"
                    } else {
                        "Show this structure"
                    })
                    .clicked();
            });
            if pick.clicked() {
                *self.current_override = Some(id);
            }
            if pick.double_clicked() {
                ui.data_mut(|d| d.insert_temp(key, (id, label.clone())));
            }
            if toggle {
                self.dispatch(Command::ShowStructure {
                    id,
                    visible: !visible,
                });
            }
            match clicked {
                Some("Rename") => {
                    ui.data_mut(|d| d.insert_temp(key, (id, label.clone())));
                }
                Some("File info") => {
                    *self.info_target = Some(id);
                    *self.layout_request = Some(LayoutRequest::OpenPanel(Tab::Info));
                }
                Some("Export structure…") => {
                    self.export_structure_dialog.structure = Some(id);
                    self.export_structure_dialog.open = true;
                }
                Some("Close") => {
                    self.dispatch(Command::CloseStructure { id });
                    self.undoable(&format!("Closed {label}."));
                }
                _ => {}
            }
        }
    }

    /// The current structure's selections: one row per rep -- this *is*
    /// the rep model, not a parallel one -- its expression and style, an
    /// eye, and ⋯ (Style, Coloring, Material, then row actions), then a
    /// field to add another. A row's own atoms outline in the viewport
    /// while it is hovered (`row_highlight`'s doc).
    fn selections_ui(&mut self, ui: &mut Ui, id: StructureId) {
        let Some(loaded) = self.scene.structure(id) else {
            return;
        };
        let reps: Vec<vv_scene::Rep> = loaded.reps.clone();
        let current_rep = loaded.current_rep;
        let key = selection_expr_key();
        let renaming = renaming_id::<(StructureId, vv_scene::RepId)>(ui, key);
        let mut hovered_rep = None;
        for (i, rep) in reps.iter().enumerate() {
            if renaming == Some((id, rep.id)) {
                if let Some(new_expr) = inline_rename_ui(ui, key, (id, rep.id), &rep.selection) {
                    self.dispatch(Command::SetRepSelection {
                        id,
                        rep: rep.id,
                        selection: new_expr,
                    });
                }
                continue;
            }
            let (row, pick) = self.selection_row_ui(ui, id, rep, i == current_rep, reps.len() > 1);
            if ui
                .ctx()
                .input(|i| i.pointer.hover_pos())
                .is_some_and(|p| row.rect.contains(p))
            {
                hovered_rep = Some(rep.id);
            }
            if row.clicked() {
                self.dispatch(Command::SetCurrentRep { id, rep: rep.id });
            }
            if row.double_clicked() {
                self.apply_rep_pick(ui, id, rep, i, RepPick::EditExpression);
            }
            if let Some(pick) = pick {
                self.apply_rep_pick(ui, id, rep, i, pick);
            }
            let popup_id = egui::Id::new(("selection-options-popover", id, rep.id));
            egui::Popup::new(popup_id, ui.ctx().clone(), row.rect, ui.layer_id())
                .open_memory(None)
                .close_behavior(egui::PopupCloseBehavior::CloseOnClickOutside)
                .show(|ui| {
                    ui.set_max_width(240.0);
                    self.rep_options_ui(ui, id, rep);
                });
        }
        *self.row_highlight = hovered_rep.map(|rep| (id, rep));
        ui.add_space(crate::theme::space::TIGHT);
        self.add_selection_ui(ui, id, &reps);
    }

    /// One selection row: its words, then (right to left) ⋯, the eye and a
    /// spinner while it builds. Returns the whole row's response and what
    /// its ⋯ menu or eye asked for.
    fn selection_row_ui(
        &self,
        ui: &mut Ui,
        id: StructureId,
        rep: &vv_scene::Rep,
        is_current: bool,
        can_delete: bool,
    ) -> (egui::Response, Option<RepPick>) {
        use egui_phosphor::regular as icon;
        let building = self.gpu_cache.rep_building(
            id,
            rep.id,
            rep.representation == vv_scene::Representation::Cartoon,
        );
        let mut pick = None;
        let mut trailing = None;
        let words = format!("{} · {}", rep.selection, crate::ribbon::style_label(rep));
        let text = widgets::list_row(ui, &words, is_current, |ui| {
            let menu = widgets::menu(ui, |ui| pick = rep_menu_ui(ui, rep, can_delete));
            let eye = if rep.visible {
                icon::EYE
            } else {
                icon::EYE_SLASH
            };
            let eye = widgets::button(ui, eye, "", Variant::Ghost).on_hover_text(if rep.visible {
                "Hide this selection"
            } else {
                "Show this selection"
            });
            if eye.clicked() {
                pick = Some(RepPick::ToggleVisible);
            }
            if building {
                ui.add(egui::Spinner::new().size(14.0))
                    .on_hover_text("Building…");
            }
            trailing = Some(menu.rect);
        });
        let text = text.on_hover_text(crate::ribbon::rep_summary(rep));
        let rect = trailing.map_or(text.rect, |r| text.rect.union(r));
        (text.with_new_rect(rect), pick)
    }

    fn apply_rep_pick(
        &mut self,
        ui: &Ui,
        id: StructureId,
        rep: &vv_scene::Rep,
        index: usize,
        pick: RepPick,
    ) {
        let edit =
            |command| Command::Batch(vec![Command::SetCurrentRep { id, rep: rep.id }, command]);
        match pick {
            RepPick::Style(representation) => self.dispatch(edit(Command::SetRepresentation {
                id,
                rep: rep.id,
                representation,
            })),
            RepPick::Coloring(coloring) => self.dispatch(edit(Command::SetColoring {
                id,
                rep: rep.id,
                coloring,
            })),
            RepPick::Material(material) => self.dispatch(edit(Command::SetMaterial {
                id,
                rep: rep.id,
                material,
            })),
            RepPick::ToggleVisible => self.dispatch(Command::ShowRep {
                id,
                rep: rep.id,
                visible: !rep.visible,
            }),
            RepPick::SelectAtoms => {
                self.dispatch(Command::SelectExpr {
                    id,
                    expr: rep.selection.clone(),
                });
                let n = self
                    .scene
                    .active_selection()
                    .map_or(0, |a| a.mask.count_ones(..));
                self.undoable(&format!("Selected {}.", count(n, "atom")));
            }
            RepPick::EditExpression => {
                let editing = ((id, rep.id), rep.selection.clone());
                ui.data_mut(|d| d.insert_temp(selection_expr_key(), editing));
            }
            RepPick::Duplicate => {
                if let Some(loaded) = self.scene.structure(id) {
                    let new_rep = vv_scene::Rep {
                        id: vv_scene::RepId(loaded.next_rep_id),
                        ..rep.clone()
                    };
                    self.dispatch(Command::AddRep {
                        id,
                        index: index + 1,
                        rep: new_rep,
                    });
                    self.undoable("Duplicated the selection.");
                }
            }
            RepPick::Delete => {
                self.dispatch(Command::RemoveRep { id, rep: rep.id });
                self.undoable("Deleted the selection.");
            }
            RepPick::Options => {
                self.dispatch(Command::SetCurrentRep { id, rep: rep.id });
                egui::Popup::open_id(
                    ui.ctx(),
                    egui::Id::new(("selection-options-popover", id, rep.id)),
                );
            }
        }
    }

    /// The list's own "+" (UX_PATTERNS "Create"): an expression field,
    /// checked as typed. Enter or the "+" adds a new selection with the current
    /// row's style/coloring/material/options -- matching Duplicate and
    /// Add rep's own defaulting, so there is one rule, not a guess at
    /// what the new row "should" look like.
    fn add_selection_ui(&mut self, ui: &mut Ui, id: StructureId, reps: &[vv_scene::Rep]) {
        use egui_phosphor::regular as icon;
        let field_id = egui::Id::new(("add-selection", id));
        let mut text = ui
            .data_mut(|d| d.get_temp::<String>(field_id))
            .unwrap_or_default();
        let (edit, plus) = ui
            .horizontal(|ui| {
                let edit = ui.add(
                    egui::TextEdit::singleline(&mut text)
                        .id(egui::Id::new(SELECTION_FIELD_ID))
                        .hint_text("chain A and name CA")
                        .desired_width(ui.available_width() - 36.0),
                );
                let plus = widgets::button(ui, icon::PLUS, "", Variant::Secondary)
                    .on_hover_text("Add a selection");
                (edit, plus)
            })
            .inner;
        let t = crate::theme::Tokens::current(ui.ctx());
        match self.add_selection_preview(ui, id, &text) {
            Some(Ok(n)) => {
                ui.label(
                    RichText::new(format!("Matches {}", count(n, "atom")))
                        .small()
                        .color(t.success),
                );
            }
            Some(Err(e)) => {
                ui.label(RichText::new(e).small().color(t.danger));
            }
            None => {}
        }
        let enter = edit.lost_focus() && ui.input(|i| i.key_pressed(egui::Key::Enter));
        if !(enter || plus.clicked()) {
            ui.data_mut(|d| d.insert_temp(field_id, text));
            return;
        }
        let expr = text.trim().to_owned();
        if expr.is_empty() {
            return;
        }
        if let Err(e) = vv_core::select::parse(&expr) {
            self.fail(format!("bad selection `{expr}`: {e}"));
            ui.data_mut(|d| d.insert_temp(field_id, text));
            return;
        }
        let Some(loaded) = self.scene.structure(id) else {
            return;
        };
        let current = &reps[loaded.current_rep];
        let new_rep = vv_scene::Rep {
            id: vv_scene::RepId(loaded.next_rep_id),
            selection: expr.clone(),
            ..current.clone()
        };
        let index = reps.len();
        self.dispatch(Command::AddRep {
            id,
            index,
            rep: new_rep,
        });
        self.undoable(&format!("Added the selection {expr}."));
        ui.data_mut(|d| d.remove::<String>(field_id));
    }

    /// The add-selection field's live match count/error, debounced 300 ms
    /// (a selection over millions of atoms takes a noticeable moment) -- kept in `ui`'s
    /// temp storage, per structure, so it needs no new `AppUi` field.
    fn add_selection_preview(
        &self,
        ui: &Ui,
        id: StructureId,
        text: &str,
    ) -> Option<Result<usize, String>> {
        const PAUSE: std::time::Duration = std::time::Duration::from_millis(300);
        let trimmed = text.trim();
        if trimmed.is_empty() {
            return None;
        }
        type Cache = (String, Instant, Option<Result<usize, String>>);
        let key = egui::Id::new(("add-selection-preview", id));
        let mut cache = ui
            .data(|d| d.get_temp::<Cache>(key))
            .filter(|(seen, ..)| seen == trimmed)
            .unwrap_or_else(|| (trimmed.to_owned(), Instant::now(), None));
        if let Some(result) = cache.2.clone() {
            return Some(result);
        }
        let waited = cache.1.elapsed();
        if waited < PAUSE {
            ui.data_mut(|d| d.insert_temp(key, cache));
            ui.ctx().request_repaint_after(PAUSE - waited);
            return None;
        }
        let loaded = self.scene.structure(id)?;
        let result = vv_core::select(
            &loaded.structure.topology,
            loaded.structure.frame(loaded.frame).positions(),
            trimmed,
        )
        .map(|mask| mask.count_ones(..))
        .map_err(|e| e.to_string());
        cache.2 = Some(result.clone());
        ui.data_mut(|d| d.insert_temp(key, cache));
        Some(result)
    }

    /// The current rep's options (`vv_scene::Representation::options`),
    /// the first two shown directly and the rest (only Tube's three today)
    /// folded under "Advanced". A drag previews and commits (a toast with
    /// Undo, item 5b) as one undo step when released; `radius_by` (Tube)
    /// is a two-way choice stored as `0.0`/`1.0`, shown as a segmented
    /// control rather than a slider.
    fn rep_options_ui(&mut self, ui: &mut Ui, id: StructureId, rep: &vv_scene::Rep) {
        const ALWAYS_SHOWN: usize = 2;
        if let vv_scene::ColorScheme::Constant(rgb) = rep.coloring {
            widgets::section(ui, "Color");
            let key = egui::Id::new(("solid color", id, rep.id));
            if let Some(rgb) = solid_color_ui(ui, key, rgb) {
                self.dispatch(Command::SetColoring {
                    id,
                    rep: rep.id,
                    coloring: vv_scene::ColorScheme::Constant(rgb),
                });
                self.undoable("Set the color.");
            }
        }
        let options = rep.representation.options();
        if options.is_empty() {
            widgets::caption(ui, "This style has no options.");
            return;
        }
        widgets::section(ui, "Options");
        let mut set = Vec::new();
        for o in options.iter().take(ALWAYS_SHOWN) {
            rep_option_ui(ui, id, rep, o, &mut set);
        }
        if options.len() > ALWAYS_SHOWN {
            egui::CollapsingHeader::new("Advanced")
                .id_salt(("rep options advanced", id, rep.id))
                .show(ui, |ui| {
                    for o in &options[ALWAYS_SHOWN..] {
                        rep_option_ui(ui, id, rep, o, &mut set);
                    }
                });
        }
        if !rep.options.is_empty()
            && widgets::button(ui, "", "Defaults", Variant::Ghost)
                .on_hover_text("Every option of this rep back to its default")
                .clicked()
        {
            set.extend(options.iter().map(|o| (o.name, None)));
        }
        if !set.is_empty() {
            let names: Vec<&str> = set.iter().map(|(name, _)| *name).collect();
            let commands = set
                .into_iter()
                .map(|(name, value)| Command::SetRepOption {
                    id,
                    rep: rep.id,
                    name: name.into(),
                    value,
                })
                .collect();
            self.dispatch(Command::Batch(commands));
            self.undoable(&format!("Set {}.", names.join(", ")));
        }
    }

    /// The current structure's selection layers (`selections_ui`).
    fn selection_ui(&mut self, ui: &mut Ui) {
        match self.current() {
            Some(id) => self.selections_ui(ui, id),
            None => {
                let icon = egui_phosphor::regular::SELECTION;
                if widgets::empty_state(ui, icon, "No structure yet", "Open a structure") {
                    self.open_file_dialog();
                }
            }
        }
    }

    /// The selection's details, a Contacts…/Surface area result card,
    /// and every measurement, label and caption in the scene.
    fn inspector_ui(&mut self, ui: &mut Ui) {
        self.selection_details_ui(ui);
        self.result_card_ui(ui);
        self.annotations_ui(ui);
    }

    fn selection_details_ui(&mut self, ui: &mut Ui) {
        let Some(active) = self.scene.active_selection() else {
            ui.weak("No selection.");
            return;
        };
        let Some(loaded) = self.scene.structure(active.structure) else {
            ui.weak("Selected structure was closed.");
            return;
        };
        let count = active.mask.count_ones(..);
        if count == 0 {
            ui.weak("Empty selection.");
            return;
        }
        if count == 1 {
            let atom = active.mask.ones().next().unwrap();
            let top = &loaded.structure.topology;
            let residue_idx = top.residue_index[atom] as usize;
            let chain_idx = top.chain_of_atom(atom) as usize;
            // The current frame, not always the first one.
            let p = loaded.structure.frame(loaded.frame).positions()[atom];
            widgets::table(
                ui,
                "inspector",
                &[
                    ("Atom", top.atom_name(atom).to_owned()),
                    ("Element", Element::symbol(top.element[atom]).to_owned()),
                    (
                        "Residue",
                        format!(
                            "{} {}",
                            top.residue_name(residue_idx),
                            top.residues[residue_idx].auth_seq_id
                        ),
                    ),
                    ("Chain", top.chain_name(chain_idx).to_owned()),
                    (
                        "B-factor",
                        format!("{:.2}", top.b_factor.get(atom).copied().unwrap_or(0.0)),
                    ),
                    ("Position", format!("{:.3}, {:.3}, {:.3} Å", p.x, p.y, p.z)),
                ],
            );
            widgets::caption(
                ui,
                "Ctrl+click more atoms: 2 = distance, 3 = angle, 4 = dihedral.",
            );
        } else {
            widgets::section(
                ui,
                &format!("{} in {}", widgets::count(count, "atom"), loaded.label),
            );
            if (2..=4).contains(&count) {
                let atoms = self.ordered_selection(&active.mask);
                let top = &loaded.structure.topology;
                // The current frame, so this follows a playing
                // trajectory the same way the viewport's own measurements do.
                let coords = loaded.structure.frame(loaded.frame);
                let p = coords.positions();
                let chain = atoms
                    .iter()
                    .map(|&a| atom_tag(top, a))
                    .collect::<Vec<_>>()
                    .join("\u{2013}");
                let line = match atoms[..] {
                    [a, b] => {
                        let what = if loaded.bonds.contains(a as u32, b as u32) {
                            "Bond"
                        } else {
                            "Distance"
                        };
                        format!("{what} {chain}: {:.3} \u{c5}", vv_core::distance(p, a, b))
                    }
                    [a, b, c] => format!(
                        "Angle {chain}: {:.2}\u{b0}",
                        vv_core::analysis::angle(p, a, b, c)
                    ),
                    [a, b, c, d] => format!(
                        "Dihedral {chain}: {:.2}\u{b0}",
                        vv_core::analysis::dihedral(p, a, b, c, d)
                    ),
                    _ => unreachable!(),
                };
                ui.separator();
                ui.label(line);
                if count > 2 {
                    ui.weak("in click order; Ctrl+click adds the next atom");
                }
                ui.separator();
            }
            self.refresh_selection_summary();
            let rows: Vec<(&str, String)> = self
                .selection_summary
                .by_element
                .iter()
                .map(|&(symbol, n)| (symbol, n.to_string()))
                .collect();
            widgets::table(ui, "by element", &rows);
        }
    }

    /// Analyze ▸ Contacts…/Surface area's own result (UX_PATTERNS
    /// "Feedback"): the text `run_command_logged` already stashed, with
    /// Copy, replaced by the next result or dismissed by hand.
    fn result_card_ui(&mut self, ui: &mut Ui) {
        let Some(card) = self.result_card.clone() else {
            return;
        };
        widgets::section(ui, &card.title);
        ui.add(egui::Label::new(&card.body).wrap());
        ui.horizontal(|ui| {
            if widgets::button(ui, "", "Copy", Variant::Secondary).clicked() {
                ui.ctx().copy_text(card.body.clone());
            }
            if widgets::button(ui, "", "Dismiss", Variant::Ghost).clicked() {
                *self.result_card = None;
            }
        });
    }

    /// Every measurement, label and caption in the scene: each row
    /// deletes (undoable) with its "x", and -- for measurements and
    /// labels, which are anchored to atoms -- selects those atoms on
    /// click. A caption has no atoms to select, so its row only deletes.
    fn annotations_ui(&mut self, ui: &mut Ui) {
        use egui_phosphor::regular as icon;
        enum AnnotationDelete {
            Measurement(vv_scene::Measurement),
            Label(u32),
        }
        struct Row {
            id: StructureId,
            words: String,
            select: Option<Vec<u32>>,
            delete: AnnotationDelete,
        }
        let mut rows = Vec::new();
        for (id, loaded) in self.scene.structures() {
            let top = &loaded.structure.topology;
            let coords = loaded.structure.frame(loaded.frame);
            let positions = coords.positions();
            for m in &loaded.measurements {
                let tags: Vec<String> = m
                    .atoms()
                    .iter()
                    .map(|&a| atom_tag(top, a as usize))
                    .collect();
                let kind = m.kind();
                let title = kind[..1].to_uppercase() + &kind[1..];
                rows.push(Row {
                    id,
                    words: format!("{title} {}: {}", tags.join("\u{2013}"), m.text(positions)),
                    select: Some(m.atoms().to_vec()),
                    delete: AnnotationDelete::Measurement(m.clone()),
                });
            }
            for (&atom, text) in &loaded.labels {
                rows.push(Row {
                    id,
                    words: format!("{}: {text}", atom_tag(top, atom as usize)),
                    select: Some(vec![atom]),
                    delete: AnnotationDelete::Label(atom),
                });
            }
        }
        let captions = self.scene.captions().to_vec();
        if rows.is_empty() && captions.is_empty() {
            return;
        }
        widgets::section(ui, "Annotations");
        for row in rows {
            let mut delete = false;
            let clicked = widgets::list_row(ui, &row.words, false, |ui| {
                delete = widgets::button(ui, icon::X, "", Variant::Ghost)
                    .on_hover_text("Remove (Ctrl+Z to undo)")
                    .clicked();
            })
            .clicked();
            if clicked {
                if let Some(atoms) = &row.select {
                    self.select_atoms(row.id, atoms, false);
                }
            }
            if delete {
                match row.delete {
                    AnnotationDelete::Measurement(measurement) => {
                        self.dispatch(Command::SetMeasurement {
                            id: row.id,
                            measurement,
                            shown: false,
                        });
                    }
                    AnnotationDelete::Label(atom) => {
                        self.dispatch(Command::SetLabel {
                            id: row.id,
                            atom,
                            text: None,
                        });
                    }
                }
                self.undoable("Removed.");
            }
        }
        for caption in captions {
            let words = format!("{}: {}", caption.name, caption.text);
            let mut delete = false;
            widgets::list_row(ui, &words, false, |ui| {
                delete = widgets::button(ui, icon::X, "", Variant::Ghost)
                    .on_hover_text("Remove (Ctrl+Z to undo)")
                    .clicked();
            });
            if delete {
                self.dispatch(Command::DeleteCaption {
                    name: caption.name.clone(),
                });
                self.undoable(&format!("Removed the caption {}.", caption.name));
            }
        }
    }

    /// The log, with a console line under it: type a command
    /// (docs/COMMANDS.md), press Enter, and the echo and result appear
    /// above. The same commands `--exec` and the palette run.
    fn log_ui(&mut self, ui: &mut Ui) {
        let mut submitted = false;
        egui::Panel::bottom("console_line")
            .show_separator_line(false)
            .show(ui, |ui| {
                ui.horizontal(|ui| {
                    ui.monospace(">");
                    let resp = ui.add(
                        egui::TextEdit::singleline(&mut self.console.input)
                            .hint_text("command (try `help`)")
                            .font(egui::TextStyle::Monospace)
                            .desired_width(f32::INFINITY),
                    );
                    if resp.lost_focus() && ui.input(|i| i.key_pressed(egui::Key::Enter)) {
                        submitted = true;
                        resp.request_focus();
                    }
                });
            });
        egui::ScrollArea::vertical()
            .stick_to_bottom(true)
            .auto_shrink([false, false])
            .show(ui, |ui| {
                for line in self.log.iter() {
                    ui.label(line);
                }
            });
        if submitted {
            let line = std::mem::take(&mut self.console.input);
            self.run_command_logged(&line);
        }
    }

    /// A structure's annotations, as parsed from its file: the usual
    /// summary first, then every kept category as a table. Shows
    /// `info_target` (a row's "File info") when it names a structure
    /// still loaded, else the current one.
    fn info_ui(&mut self, ui: &mut Ui) {
        let id = self
            .info_target
            .filter(|id| self.scene.structure(*id).is_some())
            .or_else(|| self.current());
        let Some(id) = id else {
            ui.weak("No structure loaded.");
            return;
        };
        let Some(loaded) = self.scene.structure(id) else {
            return;
        };
        let t = &loaded.structure.topology;
        let ann = &t.annotations;
        ui.label(RichText::new(&loaded.label).strong());
        if let Some(path) = &loaded.path {
            ui.weak(path.display().to_string());
        }
        ui.separator();
        let summary: Vec<(&str, Option<String>)> = vec![
            ("Title", ann.title().map(str::to_owned)),
            ("Method", ann.method().map(str::to_owned)),
            (
                "Resolution",
                ann.resolution().map(|r| format!("{r:.2} \u{c5}")),
            ),
            ("Deposited", ann.deposition_date().map(str::to_owned)),
            ("Organism", ann.organism().map(str::to_owned)),
            ("Keywords", ann.keywords().map(str::to_owned)),
            ("Citation", ann.citation_title().map(str::to_owned)),
            ("DOI", ann.doi().map(str::to_owned)),
            (
                "UniProt",
                Some(ann.uniprot_accessions())
                    .filter(|v| !v.is_empty())
                    .map(|v| v.join(", ")),
            ),
            (
                "Counts",
                Some(
                    [
                        count(t.atom_count(), "atom"),
                        count(t.residue_count(), "residue"),
                        count(t.chain_count(), "chain"),
                        count(loaded.structure.frame_count(), "frame"),
                    ]
                    .join(", "),
                ),
            ),
            (
                "Residues",
                vv_core::residue_class::describe(&t.class_counts),
            ),
        ];
        let rows: Vec<(&str, String)> = summary
            .into_iter()
            .filter_map(|(key, value)| Some((key, value?)))
            .collect();
        widgets::table(ui, "info summary", &rows);
        let entities = ann.entities();
        if !entities.is_empty() {
            widgets::section(ui, "Entities");
            for (eid, description) in &entities {
                ui.label(format!("{eid}: {description}"));
            }
        }
        if ann.is_empty() {
            widgets::caption(ui, "This file carries no header annotations.");
            return;
        }
        ui.add_space(crate::theme::space::PAD);
        egui::CollapsingHeader::new(format!(
            "Every category in the file ({})",
            ann.categories.len()
        ))
        .show(ui, |ui| self.raw_categories_ui(ui, &ann.categories));
    }

    /// The file's mmCIF categories as it wrote them, one table each.
    fn raw_categories_ui(
        &self,
        ui: &mut Ui,
        categories: &[vv_core::annotations::AnnotationCategory],
    ) {
        for category in categories {
            egui::CollapsingHeader::new(format!(
                "_{} ({} row{})",
                category.name,
                category.rows.len(),
                if category.rows.len() == 1 { "" } else { "s" }
            ))
            .id_salt(&category.name)
            .show(ui, |ui| {
                egui::Grid::new(("info_cat", &category.name))
                    .num_columns(2)
                    .spacing([12.0, 2.0])
                    .striped(true)
                    .show(ui, |ui| {
                        for (r, row) in category.rows.iter().enumerate().take(50) {
                            for (item, value) in category.items.iter().zip(row) {
                                if value.is_empty() {
                                    continue;
                                }
                                let key = if category.rows.len() > 1 {
                                    format!("[{r}] {item}")
                                } else {
                                    item.clone()
                                };
                                ui.label(RichText::new(key).monospace().weak());
                                ui.add(egui::Label::new(value).wrap());
                                ui.end_row();
                            }
                        }
                    });
            });
        }
    }

    /// Everything about how the scene looks that isn't document state.
    /// A superset of the View menu's quick toggles (style, occlusion,
    /// adaptive quality, FXAA, outline) with room to explain them, plus
    /// background color, supersampling, and the point-threshold slider,
    /// which the View menu doesn't expose.
    /// The scene's look, all in one place: lights, background, effects,
    /// the clip plane and detail. The ribbon's Lighting tab switches the
    /// same things with the same controls.
    /// Frame scrubber for a multi-model structure (an NMR ensemble, a
    /// multi-model PDB). Writes `loaded.frame` directly via
    /// `Scene::set_frame_live` (not `Command::SetFrame`: scrubbing this
    /// slider changes many times a second and would flood undo history
    /// otherwise — see `LoadedStructure::frame`'s doc). `State` uploads
    /// the new coordinates before the next render and advances playback.
    fn timeline_ui(&mut self, ui: &mut Ui) {
        use egui_phosphor::regular as icon;
        let Some((id, count)) = timeline_target(self.scene) else {
            widgets::caption(
                ui,
                "No structure with several frames. Load a trajectory (Select tab) or open an NMR ensemble to play it here.",
            );
            return;
        };
        let loaded = self.scene.structure(id);
        let label = loaded.map(|s| s.label.clone()).unwrap_or_default();
        let last = count - 1;
        let mut frame = loaded.map_or(0, |s| s.frame).min(last);
        ui.horizontal(|ui| {
            ui.label(format!("{label}  ·  {}", widgets::count(count, "frame")));
            if widgets::icon_button(ui, icon::SKIP_BACK, false)
                .on_hover_text("First frame")
                .clicked()
            {
                frame = 0;
            }
            if widgets::icon_button(ui, icon::CARET_LEFT, false)
                .on_hover_text("Previous frame")
                .clicked()
            {
                frame = frame.saturating_sub(1);
            }
            let (glyph, tip) = if self.timeline.playing {
                (icon::PAUSE, "Pause")
            } else {
                (icon::PLAY, "Play")
            };
            if widgets::icon_button(ui, glyph, self.timeline.playing)
                .on_hover_text(tip)
                .clicked()
            {
                self.timeline.playing = !self.timeline.playing;
                self.timeline.last_step = None;
            }
            if widgets::icon_button(ui, icon::CARET_RIGHT, false)
                .on_hover_text("Next frame")
                .clicked()
            {
                frame = (frame + 1).min(last);
            }
            if widgets::icon_button(ui, icon::SKIP_FORWARD, false)
                .on_hover_text("Last frame")
                .clicked()
            {
                frame = last;
            }
            ui.add(
                egui::DragValue::new(&mut self.timeline.fps)
                    .range(0.5..=120.0)
                    .suffix(" fps"),
            )
            .on_hover_text("Playback speed");
            if widgets::switch(ui, self.timeline.looping, "Loop").clicked() {
                self.timeline.looping = !self.timeline.looping;
            }
        });
        widgets::slider_int(ui, "Frame", &mut frame, 0..=last);
        self.scene.set_frame_live(id, frame);
    }

    /// Open by path: for when the system file picker is unavailable or
    /// the path is already on the clipboard.
    /// File ▸ Open…: a typed path or Browse for the platform's own
    /// picker, so there's one Open entry point instead of a separate
    /// "Open path" prompt.
    pub fn open_dialog(&mut self, ctx: &Context) {
        if !self.dialog.open {
            return;
        }
        let path = &mut self.dialog.path;
        let mut browse = false;
        let answer = widgets::dialog(
            ctx,
            "open path",
            "Open a structure",
            "Open",
            Variant::Primary,
            |ui| {
                ui.label("Path");
                ui.add(egui::TextEdit::singleline(path).desired_width(f32::INFINITY))
                    .request_focus();
                if widgets::button(ui, "", "Browse…", Variant::Secondary).clicked() {
                    browse = true;
                }
                widgets::caption(ui, "A .cif, .mmcif or .pdb file, optionally .gz.");
            },
        );
        if browse {
            self.dialog.open = false;
            self.open_file_dialog();
            return;
        }
        match answer {
            widgets::Answer::Confirm => {
                let path = std::mem::take(&mut self.dialog.path);
                self.dialog.open = false;
                self.open_path(path);
            }
            widgets::Answer::Cancel => self.dialog.open = false,
            widgets::Answer::Pending => {}
        }
    }

    /// Fetch by PDB ID: one or more, optionally the biological assembly
    /// instead of the asymmetric unit. Runs the `fetch` command per ID, so
    /// the download is the same code path as the console.
    pub fn fetch_dialog_ui(&mut self, ctx: &Context) {
        if !self.fetch_dialog.open {
            return;
        }
        let state = &mut self.fetch_dialog;
        let answer = widgets::dialog(
            ctx,
            "fetch",
            "Fetch from the PDB",
            "Fetch",
            Variant::Primary,
            |ui| {
                ui.label("PDB IDs");
                ui.add(
                    egui::TextEdit::singleline(&mut state.ids)
                        .hint_text("4HHB 1CRN")
                        .desired_width(f32::INFINITY),
                )
                .request_focus();
                widgets::caption(ui, "Separate several with spaces or commas. Downloads are kept, never fetched twice.");
                ui.add_space(crate::theme::space::GAP);
                ui.horizontal(|ui| {
                    if widgets::switch(ui, state.assembly, "Biological assembly").clicked() {
                        state.assembly = !state.assembly;
                    }
                    ui.add_enabled(
                        state.assembly,
                        egui::DragValue::new(&mut state.assembly_number).range(1..=99),
                    );
                });
            },
        );
        let ids: Vec<String> = self
            .fetch_dialog
            .ids
            .split(|c: char| c.is_whitespace() || c == ',')
            .filter(|w| !w.is_empty())
            .map(str::to_owned)
            .collect();
        match answer {
            widgets::Answer::Confirm if !ids.is_empty() => {
                self.fetch_dialog.open = false;
                let suffix = if self.fetch_dialog.assembly {
                    format!(" assembly {}", self.fetch_dialog.assembly_number)
                } else {
                    String::new()
                };
                for id in ids {
                    if self.run_command_logged(&format!("fetch {id}{suffix}")) {
                        self.reset_camera();
                    }
                }
            }
            widgets::Answer::Cancel => self.fetch_dialog.open = false,
            _ => {}
        }
    }

    /// Screenshot export: a path and 2x supersampling. Submitting records
    /// the request in `export_request`; `State` renders and writes it after
    /// this frame's UI pass, since that needs the renderer.
    pub fn export_dialog_ui(&mut self, ctx: &Context) {
        if !self.export_dialog.open {
            return;
        }
        let state = &mut self.export_dialog;
        let answer = widgets::dialog(
            ctx,
            "export",
            "Save a screenshot",
            "Save",
            Variant::Primary,
            |ui| {
                ui.label("Save to");
                save_path_field(ui, &mut state.path, &["png", "jpg", "svg"]);
                widgets::caption(ui, "A .png, .jpg or .svg file.");
                ui.add_space(crate::theme::space::GAP);
                if widgets::switch(ui, state.ssaa, "2x supersampling").clicked() {
                    state.ssaa = !state.ssaa;
                }
            },
        );
        match answer {
            widgets::Answer::Confirm if !self.export_dialog.path.trim().is_empty() => {
                *self.export_request = Some(ExportRequest {
                    path: self.export_dialog.path.trim().into(),
                    ssaa: self.export_dialog.ssaa,
                });
                self.export_dialog.open = false;
            }
            widgets::Answer::Cancel => self.export_dialog.open = false,
            _ => {}
        }
    }

    /// Format, atoms and frames, then the platform's own save dialog
    /// -- writes with `vv_io::save` directly against `structure` rather
    /// than running `savestructure` (which always acts on whichever
    /// structure loaded last), so exporting from a non-current row's ⋯
    /// menu exports that row, not the current one.
    pub fn export_structure_dialog_ui(&mut self, ctx: &Context) {
        if !self.export_structure_dialog.open {
            return;
        }
        let Some(id) = self
            .export_structure_dialog
            .structure
            .or_else(|| self.current())
        else {
            self.export_structure_dialog.open = false;
            return;
        };
        let Some(label) = self.scene.structure(id).map(|s| s.label.clone()) else {
            self.export_structure_dialog.open = false;
            return;
        };
        let has_selection = self
            .scene
            .active_selection()
            .is_some_and(|a| a.structure == id);
        let names: Vec<&str> = EXPORT_FORMATS.iter().map(|(n, _)| *n).collect();
        let state = &mut self.export_structure_dialog;
        let answer = widgets::dialog(
            ctx,
            "export structure",
            &format!("Export {label}"),
            "Export…",
            Variant::Primary,
            |ui| {
                ui.label("Format");
                if let Some(i) = widgets::segmented(ui, &names, Some(state.format)) {
                    state.format = i;
                }
                ui.label("Atoms");
                if let Some(i) = widgets::segmented(
                    ui,
                    &["All", "Current selection"],
                    Some(state.current_selection as usize),
                ) {
                    state.current_selection = i == 1;
                }
                if state.current_selection && !has_selection {
                    widgets::caption(
                        ui,
                        "No selection on this structure: every atom will be written.",
                    );
                }
                ui.label("Frames");
                if let Some(i) = widgets::segmented(
                    ui,
                    &["Current frame", "All frames"],
                    Some(state.all_frames as usize),
                ) {
                    state.all_frames = i == 1;
                }
            },
        );
        match answer {
            widgets::Answer::Cancel => self.export_structure_dialog.open = false,
            widgets::Answer::Pending => {}
            widgets::Answer::Confirm => {
                self.export_structure_dialog.open = false;
                self.run_export_structure(id, &label, has_selection);
            }
        }
    }

    /// The save dialog and write itself, once the form is confirmed.
    fn run_export_structure(&mut self, id: StructureId, label: &str, has_selection: bool) {
        let (format_name, ext) = EXPORT_FORMATS[self.export_structure_dialog.format];
        let mut picker = rfd::FileDialog::new()
            .set_title("Export structure")
            .set_file_name(format!("{label}.{}", ext.unwrap_or("pdb")));
        picker = match ext {
            Some(ext) => picker.add_filter(format_name, &[ext]),
            None => picker.add_filter(
                "Structure",
                &["pdb", "ent", "cif", "mmcif", "xyz", "pqr", "gro"],
            ),
        };
        let Some(mut path) = picker.save_file() else {
            return;
        };
        if let Some(ext) = ext {
            if !path
                .extension()
                .is_some_and(|e| e.eq_ignore_ascii_case(ext))
            {
                path.set_extension(ext);
            }
        }
        let Some(loaded) = self.scene.structure(id) else {
            return;
        };
        let mask = (self.export_structure_dialog.current_selection && has_selection)
            .then(|| self.scene.active_selection())
            .flatten()
            .map(|a| a.mask.clone());
        let frames: Vec<usize> = if self.export_structure_dialog.all_frames {
            (0..loaded.structure.frame_count()).collect()
        } else {
            vec![loaded.frame]
        };
        let opts = vv_io::SaveOptions {
            atoms: mask.as_deref(),
            frames: &frames,
        };
        match vv_io::save(&loaded.structure, &path, &opts) {
            Ok(warnings) => {
                let atoms = mask
                    .as_ref()
                    .map_or(loaded.structure.atom_count(), |m| m.count_ones(..));
                self.log.push(format!(
                    "saved {} atom(s), {} frame(s) to {}",
                    atoms,
                    frames.len(),
                    path.display()
                ));
                self.log
                    .extend(warnings.iter().map(|w| format!("warning: {w}")));
                let message = if warnings.is_empty() {
                    format!("Saved {}.", path.display())
                } else {
                    format!("Saved {} — {}", path.display(), warnings.join("; "))
                };
                *self.notice = Some(Notice::info(message));
            }
            Err(e) => self.fail(e.to_string()),
        }
    }

    /// Preferences: Only me, Performance, Danger. Every
    /// control applies and saves at once (UX_PATTERNS "Settings screens");
    /// there is no Save/Cancel, only Close.
    pub fn preferences_dialog_ui(&mut self, ctx: &Context) {
        self.confirm_reset_preferences_ui(ctx);
        if !self.prefs_dialog.open {
            return;
        }
        let t = crate::theme::Tokens::current(ctx);
        let modal = egui::Modal::new(egui::Id::new("preferences"))
            .frame(
                egui::Frame::NONE
                    .fill(t.surface)
                    .stroke(egui::Stroke::new(1.0, t.border))
                    .corner_radius(crate::theme::radius::CONTAINER)
                    .inner_margin(crate::theme::space::WIDE),
            )
            .show(ctx, |ui| {
                ui.set_width(440.0);
                ui.label(RichText::new("Preferences").text_style(crate::theme::title_style()));
                ui.add_space(crate::theme::space::INSET);
                egui::ScrollArea::vertical()
                    .max_height(480.0)
                    .show(ui, |ui| {
                        self.prefs_only_me_ui(ui);
                        self.prefs_effects_ui(ui);
                        self.prefs_performance_ui(ui);
                        self.prefs_danger_ui(ui);
                    });
                ui.add_space(crate::theme::space::WIDE);
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    widgets::button(ui, "", "Close", Variant::Secondary).clicked()
                })
                .inner
            });
        let close =
            ctx.input(|i| i.key_pressed(egui::Key::Escape)) || modal.should_close() || modal.inner;
        if close {
            self.prefs_dialog.open = false;
        }
    }

    /// Post effects over the whole view; saved with a session, undoable.
    fn prefs_effects_ui(&mut self, ui: &mut Ui) {
        widgets::section(ui, "Effects");
        widgets::caption(ui, "The whole view. Saved with a session.");
        let lighting = self.view.lighting;
        let view = &mut *self.view;
        effect_ui(
            ui,
            "Ambient occlusion",
            &mut view.ao,
            lighting.ao().max(1.0),
            2.0,
        );
        effect_ui(
            ui,
            "Shadows",
            &mut view.shadows,
            lighting.shadows().max(0.6),
            1.0,
        );
        effect_ui(
            ui,
            "Depth cue",
            &mut view.depth_cue,
            lighting.depth_cue().max(0.5),
            1.0,
        );
        effect_ui(ui, "Depth of field", &mut view.dof, 0.5, 1.0);
        if widgets::switch(ui, view.outline, "Outline").clicked() {
            view.outline = !view.outline;
        }
    }

    /// Saved per user, this computer only: theme, start behavior, recent
    /// files.
    fn prefs_only_me_ui(&mut self, ui: &mut Ui) {
        widgets::section(ui, "Only me");
        widgets::caption(ui, "Saved on this computer, for you.");
        let names = ["Dark", "Light"];
        let selected = match self.prefs.theme {
            crate::theme::ThemeMode::Dark => 0,
            crate::theme::ThemeMode::Light => 1,
        };
        if let Some(i) = widgets::segmented(ui, &names, Some(selected)) {
            self.prefs.theme = if i == 0 {
                crate::theme::ThemeMode::Dark
            } else {
                crate::theme::ThemeMode::Light
            };
        }
        if widgets::switch(ui, !self.prefs.hide_start_card, "Start card at launch")
            .on_hover_text("Open, Fetch, Open session, recent files, over an empty viewport")
            .clicked()
        {
            self.prefs.hide_start_card = !self.prefs.hide_start_card;
        }
        if widgets::switch(
            ui,
            !self.prefs.dont_confirm_quit,
            "Ask to save before quitting",
        )
        .clicked()
        {
            self.prefs.dont_confirm_quit = !self.prefs.dont_confirm_quit;
        }
        ui.horizontal(|ui| {
            if widgets::button(ui, "", "Start with this layout", Variant::Secondary)
                .on_hover_text("Every launch opens the panel arrangement you have now")
                .clicked()
            {
                *self.layout_request = Some(LayoutRequest::SetStart);
            }
            if widgets::button(ui, "", "Use the default at start", Variant::Secondary).clicked() {
                *self.layout_request = Some(LayoutRequest::ResetStart);
            }
        });
        if widgets::button(
            ui,
            "",
            &format!(
                "Clear recent files ({})",
                widgets::count(self.prefs.recent_files.len(), "file")
            ),
            Variant::Secondary,
        )
        .clicked()
        {
            self.prefs.recent_files.clear();
        }
    }

    /// Rendering cost/quality, saved with the layout like everything else
    /// in `UiPrefs` (adaptive quality and point size side by side, so
    /// the one greying the other out is visible).
    fn prefs_performance_ui(&mut self, ui: &mut Ui) {
        widgets::section(ui, "Performance");
        ui.horizontal(|ui| {
            if widgets::switch(ui, self.view.adaptive, "Adaptive quality").clicked() {
                self.view.adaptive = !self.view.adaptive;
            }
            ui.add_enabled_ui(!self.view.adaptive, |ui| {
                let (min, max) = (self.lod.min_px, self.lod.max_px);
                widgets::slider(
                    ui,
                    "Point size",
                    &mut self.lod.threshold_px,
                    min..=max,
                    " px",
                );
            });
        });
        widgets::caption(
            ui,
            "Adaptive quality moves the point-size threshold itself to hold the frame rate.",
        );
        if widgets::switch(ui, self.view.occlusion_culling, "Occlusion culling")
            .on_hover_text("Skip atoms hidden behind the previous frame")
            .clicked()
        {
            self.view.occlusion_culling = !self.view.occlusion_culling;
        }
        if widgets::switch(ui, self.view.fxaa, "Smooth edges (FXAA)").clicked() {
            self.view.fxaa = !self.view.fxaa;
        }
        let names = ["Off", "1.5x", "2x"];
        let selected = match self.view.render_scale {
            s if s >= 1.9 => 2,
            s if s >= 1.4 => 1,
            _ => 0,
        };
        ui.horizontal(|ui| {
            ui.label("Supersample");
            if let Some(i) = widgets::segmented(ui, &names, Some(selected)) {
                self.view.render_scale = [1.0, 1.5, 2.0][i];
            }
        });
    }

    fn prefs_danger_ui(&mut self, ui: &mut Ui) {
        widgets::section(ui, "Danger");
        if widgets::button(ui, "", "Reset all preferences", Variant::Danger).clicked() {
            self.prefs_dialog.confirm_reset = true;
        }
        if widgets::button(ui, "", "Delete saved layouts", Variant::Danger).clicked() {
            self.workspace_dialog.deleting = true;
        }
    }

    /// "Reset all preferences" (UX_PATTERNS "Delete": not undoable, so it
    /// confirms). Resets `UiPrefs` and the Performance group's own
    /// settings only — lights, background, clip and the rest keep their
    /// own per-popover Reset, untouched here.
    fn confirm_reset_preferences_ui(&mut self, ctx: &Context) {
        if !self.prefs_dialog.confirm_reset {
            return;
        }
        let answer = widgets::dialog(
            ctx,
            "reset preferences",
            "Reset all preferences?",
            "Reset",
            Variant::Danger,
            |ui| {
                ui.label(
                    "Theme, start behavior and performance settings go back to their \
                     defaults. This can't be undone.",
                );
            },
        );
        match answer {
            widgets::Answer::Confirm => {
                *self.prefs = crate::theme::UiPrefs::default();
                self.view.adaptive = true;
                self.view.occlusion_culling = true;
                self.view.fxaa = true;
                self.view.render_scale = 1.0;
                *self.lod = AdaptiveLod::default();
                self.prefs_dialog.confirm_reset = false;
            }
            widgets::Answer::Cancel => self.prefs_dialog.confirm_reset = false,
            widgets::Answer::Pending => {}
        }
    }

    /// Save-current-layout-as dialog. `dock_state` is a parameter for the
    /// same reason `menu_bar` takes one: whatever's rendering the dock
    /// area needs its own `&mut` borrow of it elsewhere in the frame, so
    /// it can't live as an `AppUi` field.
    pub fn workspace_dialog_ui(&mut self, ctx: &Context, dock_state: &egui_dock::DockState<Tab>) {
        self.delete_workspace_ui(ctx);
        if !self.workspace_dialog.open {
            return;
        }
        let name = &mut self.workspace_dialog.name;
        let answer = widgets::dialog(
            ctx,
            "save workspace",
            "Save workspace",
            "Save",
            Variant::Primary,
            |ui| {
                ui.label("Name");
                ui.add(egui::TextEdit::singleline(name).desired_width(f32::INFINITY))
                    .request_focus();
                widgets::caption(ui, "Saving under an existing name replaces that workspace.");
            },
        );
        let submit = match answer {
            widgets::Answer::Cancel => {
                self.workspace_dialog.open = false;
                false
            }
            widgets::Answer::Confirm => true,
            widgets::Answer::Pending => false,
        };
        if submit {
            let name = self.workspace_dialog.name.trim();
            if !name.is_empty() {
                if let Some(existing) = self.workspaces.iter_mut().find(|(n, _)| n == name) {
                    existing.1 = dock_state.clone();
                } else {
                    self.workspaces.push((name.to_string(), dock_state.clone()));
                }
                self.workspace_dialog.name.clear();
                self.workspace_dialog.open = false;
            }
        }
    }

    /// Quitting with changes no session holds: save them, or not. The
    /// "Don't remind me again" checkbox only takes effect if this dialog
    /// actually ends in quitting, not on Cancel.
    pub fn quit_dialog_ui(&mut self, ctx: &Context) {
        if !*self.quit_dialog {
            return;
        }
        let dont_ask = &mut *self.quit_dont_ask;
        let answer = widgets::dialog_with(
            ctx,
            "quit",
            "Save the session before quitting?",
            &[
                ("Quit without saving", Variant::Danger),
                ("Save session…", Variant::Primary),
            ],
            |ui| {
                ui.label("The scene has changes that no saved session holds.");
                ui.add_space(crate::theme::space::GAP);
                if widgets::switch(ui, *dont_ask, "Don't remind me again").clicked() {
                    *dont_ask = !*dont_ask;
                }
            },
        );
        let quit = match answer {
            Some(Some(0)) => true,
            Some(Some(_)) => self.run_command_logged("savesession") && self.saved(),
            _ => false,
        };
        if answer.is_some() {
            *self.quit_dialog = false;
        }
        if quit {
            if *self.quit_dont_ask {
                self.prefs.dont_confirm_quit = true;
            }
            *self.quit_confirmed = true;
            *self.quit_requested = true;
        }
        if answer.is_some() {
            *self.quit_dont_ask = false;
        }
    }

    /// Says what an undoable action did, with a clickable Undo --
    /// the button now carries the "how to take it back" the text used to
    /// spell out.
    pub(crate) fn undoable(&mut self, done: &str) {
        *self.notice = Some(Notice::info(done.to_owned()).with_action("Undo", "undo"));
    }

    /// No changes since the last saved or opened session.
    pub(crate) fn saved(&self) -> bool {
        *self.saved_edits == self.history.edits()
    }

    /// The saved workspaces to delete from (`deletelayout` alone), and the
    /// confirmation naming the one about to go.
    fn delete_workspace_ui(&mut self, ctx: &Context) {
        use egui_phosphor::regular as icon;
        if let Some(name) = self.workspace_dialog.confirm.clone() {
            let question = format!("Delete workspace \u{201c}{name}\u{201d}?");
            let answer = widgets::dialog(
                ctx,
                "confirm delete workspace",
                &question,
                "Delete",
                Variant::Danger,
                |ui| {
                    ui.label("Its panel arrangement is forgotten. This can't be undone.");
                },
            );
            match answer {
                widgets::Answer::Confirm => {
                    self.workspace_dialog.confirm = None;
                    self.run_command_logged(&format!("deletelayout {name}"));
                }
                widgets::Answer::Cancel => self.workspace_dialog.confirm = None,
                widgets::Answer::Pending => {}
            }
            return;
        }
        if !self.workspace_dialog.deleting {
            return;
        }
        let names: Vec<String> = self.workspaces.iter().map(|(n, _)| n.clone()).collect();
        let mut picked = None;
        let answer = widgets::dialog(
            ctx,
            "delete workspaces",
            "Delete a workspace",
            "Done",
            Variant::Secondary,
            |ui| {
                if names.is_empty() {
                    ui.label("No saved workspaces. The built-in ones (Default, Trajectory, Analysis, Compact) always stay.");
                }
                for name in &names {
                    widgets::list_row(ui, name, false, |ui| {
                        if widgets::button(ui, icon::TRASH, "Delete", Variant::Danger).clicked() {
                            picked = Some(name.clone());
                        }
                    });
                }
            },
        );
        if answer != widgets::Answer::Pending {
            self.workspace_dialog.deleting = false;
        }
        if picked.is_some() {
            self.workspace_dialog.confirm = picked;
        }
    }

    /// Fuzzy command palette, toggled by Ctrl+P. Lists every registered
    /// command (`commands::entries`); the text box doubles as a command
    /// line, so `select chain A` + Enter runs it directly instead of
    /// fuzzy-matching.
    pub fn palette_ui(&mut self, ctx: &Context) {
        if ctx.input(|i| i.key_pressed(egui::Key::P) && i.modifiers.ctrl) {
            self.palette.open = !self.palette.open;
        }
        if !self.palette.open {
            return;
        }
        let entries = crate::commands::entries();
        let query = self.palette.query.trim().to_owned();
        let is_command_line = crate::commands::validate(&query).is_ok()
            && !query.is_empty()
            && (query.contains(' ') || entries.iter().any(|e| e.id == query));
        let matcher = SkimMatcherV2::default();
        let mut ranked: Vec<(crate::commands::Entry, i64)> = entries
            .iter()
            .filter_map(|e| {
                if query.is_empty() {
                    return Some((*e, 0));
                }
                let haystack = format!("{} {} {}", e.title, e.id, e.keywords.join(" "));
                matcher
                    .fuzzy_match(&haystack, &query)
                    .map(|score| (*e, score))
            })
            .collect();
        ranked.sort_by_key(|(_, score)| std::cmp::Reverse(*score));

        let mut open = self.palette.open;
        let mut chosen: Option<crate::commands::Entry> = None;
        let mut run_line: Option<String> = None;
        egui::Window::new("Find a command")
            .open(&mut open)
            .collapsible(false)
            .show(ctx, |ui| {
                let resp = ui.add(
                    egui::TextEdit::singleline(&mut self.palette.query)
                        .hint_text("search, or type a command")
                        .desired_width(360.0),
                );
                resp.request_focus();
                let enter = ui.input(|i| i.key_pressed(egui::Key::Enter));
                if enter && is_command_line {
                    run_line = Some(query.clone());
                } else if enter {
                    chosen = ranked.first().map(|(e, _)| *e);
                }
                egui::ScrollArea::vertical()
                    .max_height(320.0)
                    .show(ui, |ui| {
                        for (entry, _) in &ranked {
                            // Item 5c: a bare multi-verb command (`mode`,
                            // `panel`) has no exact shortcut of its own,
                            // but its family's keys ("R/T/S/C") are still
                            // worth showing here.
                            let shortcut = crate::keys::shortcut_for(entry.id)
                                .or_else(|| crate::keys::shortcut_family(entry.id));
                            let label = match shortcut {
                                Some(shortcut) => {
                                    format!("{}  ·  {}  ({shortcut})", entry.title, entry.usage)
                                }
                                None => format!("{}  ·  {}", entry.title, entry.usage),
                            };
                            if widgets::list_row(ui, &label, false, |_| {})
                                .on_hover_text(entry.help)
                                .clicked()
                            {
                                chosen = Some(*entry);
                            }
                        }
                    });
            });
        // egui::Window has no built-in Escape-to-close (unlike Modal/
        // Popup): "Escape closes the top popover" needs it here.
        if ctx.input(|i| i.key_pressed(egui::Key::Escape)) {
            open = false;
        }
        self.palette.open = open;
        if let Some(line) = run_line {
            self.palette.open = false;
            self.palette.query.clear();
            self.run_command_logged(&line);
        } else if let Some(entry) = chosen {
            if entry.needs_args() {
                // Leave the palette open with the verb typed, ready for
                // its arguments.
                self.palette.query = format!("{} ", entry.id);
            } else {
                self.palette.open = false;
                self.palette.query.clear();
                self.run_command_logged(entry.id);
            }
        }
    }
}

impl ViewSettings {
    /// The lighting uniform these settings describe.
    pub fn lighting(&self) -> vv_render::Lighting {
        let mut lighting = self.lighting.lighting();
        self.lights.apply(&mut lighting);
        if !self.tonemap {
            lighting.tonemap = vv_render::style::TONEMAP_NONE;
        }
        lighting
    }
}

/// A lighting preset as the interface names it.
pub(crate) fn preset_label(preset: LightingPreset) -> &'static str {
    match preset {
        LightingPreset::Default => "Default",
        LightingPreset::Soft => "Soft",
        LightingPreset::Full => "Full",
        LightingPreset::Flat => "Flat",
        LightingPreset::Illustrative => "Illustrative",
        LightingPreset::Gentle => "Gentle",
        LightingPreset::TwinLight => "Twin light",
        LightingPreset::ThreePoint => "Three-point",
    }
}

/// Look ▸ Background ▾.
pub(crate) fn look_background_popover(app: &mut AppUi<'_>, ui: &mut Ui) {
    ui.set_max_width(240.0);
    let gradient = app.view.background_top.is_some();
    if widgets::switch(ui, gradient, "Gradient").clicked() {
        app.view.background_top = (!gradient).then(|| gradient_top(app.view.background));
    }
    let top_key = egui::Id::new("background-editing-top");
    let mut editing_top = gradient && ui.data(|d| d.get_temp::<bool>(top_key)).unwrap_or(false);
    if gradient {
        if let Some(i) = widgets::segmented(ui, &["Bottom", "Top"], Some(editing_top as usize)) {
            editing_top = i == 1;
            ui.data_mut(|d| d.insert_temp(top_key, editing_top));
        }
    }
    ui.add_space(crate::theme::space::GAP);
    let target = match &mut app.view.background_top {
        Some(top) if editing_top => top,
        _ => &mut app.view.background,
    };
    let mut color = color32_from_wgpu(*target);
    if widgets::color_editor(ui, egui::Id::new("background-color"), &mut color) {
        *target = wgpu_from_color32(color);
    }
    ui.add_space(crate::theme::space::GAP);
    if widgets::button(ui, "", "Reset", Variant::Secondary).clicked() {
        app.view.background = default_background(app.prefs.theme);
        app.view.background_top = None;
    }
}

/// An effect's switch and, while on, its strength. Switching on goes to
/// what the lighting preset suggests (`on`).
fn effect_ui(ui: &mut Ui, label: &str, value: &mut f32, on: f32, max: f32) {
    let enabled = *value > 0.0;
    if widgets::switch(ui, enabled, label).clicked() {
        *value = if enabled { 0.0 } else { on };
    }
    if *value > 0.0 {
        widgets::slider(ui, "Strength", value, 0.0..=max, "");
    }
}

/// View ▸ Camera's "Rotate" form (UX_PATTERNS "Create"): axis and
/// degrees, kept across frames in egui's own per-widget storage since
/// the popover has no dedicated state of its own.
pub(crate) fn view_rotate_popover(app: &mut AppUi<'_>, ui: &mut Ui) {
    ui.set_max_width(200.0);
    let id = egui::Id::new("view-rotate-form");
    let (mut axis, mut degrees) = ui
        .data_mut(|d| d.get_temp::<(usize, f32)>(id))
        .unwrap_or((1, 15.0));
    ui.label("Axis");
    if let Some(i) = widgets::segmented(ui, &["X", "Y", "Z"], Some(axis)) {
        axis = i;
    }
    ui.label("Degrees");
    ui.add(egui::DragValue::new(&mut degrees).suffix("°"));
    ui.data_mut(|d| d.insert_temp(id, (axis, degrees)));
    ui.add_space(crate::theme::space::GAP);
    if let Some(submit) = widgets::form_actions(ui, "Rotate") {
        if submit {
            let axis_word = ["x", "y", "z"][axis];
            app.run_command_logged(&format!("rotate {axis_word} {degrees}"));
        }
        ui.close();
    }
}

/// View ▸ Clip's ▾: everything but on/off, which is the switch
/// itself (`ribbon::Kind::SwitchEffect`).
pub(crate) fn view_clip_popover(app: &mut AppUi<'_>, ui: &mut Ui) {
    use egui_phosphor::regular as icon;
    ui.set_max_width(220.0);
    let Some(clip) = &mut app.view.clip else {
        return;
    };
    let facings = ["View", "X", "Y", "Z"];
    let selected = match clip.orientation {
        ClipOrientation::View => Some(0),
        ClipOrientation::X => Some(1),
        ClipOrientation::Y => Some(2),
        ClipOrientation::Z => Some(3),
        ClipOrientation::Fixed(_) => None,
    };
    ui.horizontal(|ui| {
        if let Some(i) = widgets::segmented(ui, &facings, selected) {
            clip.orientation = [
                ClipOrientation::View,
                ClipOrientation::X,
                ClipOrientation::Y,
                ClipOrientation::Z,
            ][i];
        }
        if widgets::icon_button(
            ui,
            icon::VIDEO_CAMERA,
            matches!(clip.orientation, ClipOrientation::Fixed(_)),
        )
        .on_hover_text("Fix the plane along the direction you're looking now")
        .clicked()
        {
            clip.orientation = ClipOrientation::Fixed(app.camera.orientation * glam::Vec3::NEG_Z);
        }
    });
    if widgets::switch(ui, clip.flipped, "Flip side").clicked() {
        clip.flipped = !clip.flipped;
    }
    widgets::slider(ui, "Depth", &mut clip.depth, -60.0..=60.0, " Å");
}

/// View ▸ Layout ▾: one line of layout chips (a saved layout's own delete
/// is its right-click, confirmed through `workspace_dialog`), then Save…
/// and Start with this as icon buttons.
pub(crate) fn view_layout_popover(app: &mut AppUi<'_>, ui: &mut Ui) {
    use egui_phosphor::regular as icon;
    ui.set_max_width(360.0);
    let mut run = None;
    let mut close = false;
    ui.horizontal_wrapped(|ui| {
        for (name, _) in crate::layout::built_in_workspaces() {
            if widgets::button(ui, "", name, Variant::Secondary).clicked() {
                run = Some(format!("layout {name}"));
            }
        }
        for (name, _) in app.workspaces.clone() {
            let chip = widgets::button(ui, "", &name, Variant::Secondary)
                .on_hover_text("Right-click to delete");
            if chip.clicked() {
                run = Some(format!("layout {name}"));
            }
            if chip.secondary_clicked() {
                app.workspace_dialog.confirm = Some(name.clone());
                close = true;
            }
        }
        if widgets::icon_button(ui, icon::FLOPPY_DISK, false)
            .on_hover_text("Save this layout…")
            .clicked()
        {
            app.workspace_dialog.open = true;
            close = true;
        }
        if widgets::icon_button(ui, icon::FLAG, false)
            .on_hover_text("Start with this layout")
            .clicked()
        {
            run = Some("startlayout".into());
        }
    });
    if let Some(line) = run {
        app.run_command_logged(&line);
        close = true;
    }
    if close {
        ui.close();
    }
}

/// View ▸ Panels ▾: a checklist, Ctrl+1-8 covering the eight built-in
/// tab-bar shortcuts that still mean something (`panel look`, Ctrl+4,
/// covers the Look tab now). The Movie panel has no shortcut.
pub(crate) fn view_panels_popover(app: &mut AppUi<'_>, ui: &mut Ui) {
    ui.set_max_width(180.0);
    for tab in Tab::PANELS {
        let open = !app.closed_panels.contains(&tab);
        let command = format!("panel {}", crate::ribbon::panel_word(tab));
        let row = widgets::switch(ui, open, tab.title());
        // Item 5c: each row's own Ctrl+N, not just the popover's first.
        let row = match crate::keys::shortcut_for(&command) {
            Some(shortcut) => row.on_hover_text(shortcut),
            None => row,
        };
        if row.clicked() {
            *app.layout_request = Some(if open {
                LayoutRequest::ClosePanel(tab)
            } else {
                LayoutRequest::OpenPanel(tab)
            });
        }
    }
}

/// File ▸ Add trajectory's form (UX_PATTERNS "Create"): the topology is
/// the current structure by default (the trajectory is attached to it),
/// or another file (loaded as a new structure); the trajectory
/// path is typed or Browse. `loadtrajectory` parses on whitespace with
/// no quoting, so a path with spaces still won't work -- a pre-existing
/// limit of that command, not new here.
pub(crate) fn file_trajectory_popover(app: &mut AppUi<'_>, ui: &mut Ui) {
    ui.set_max_width(280.0);
    let id = egui::Id::new("file-trajectory-form");
    let (mut use_current, mut topology, mut trajectory) = ui
        .data_mut(|d| d.get_temp::<(bool, String, String)>(id))
        .unwrap_or_else(|| (true, String::new(), String::new()));

    let current_path = app
        .current()
        .and_then(|id| app.scene.structure(id))
        .and_then(|s| s.path.clone());

    ui.label("Topology");
    if let Some(i) = widgets::segmented(
        ui,
        &["Current structure", "Choose…"],
        Some(usize::from(!use_current)),
    ) {
        use_current = i == 0;
    }
    let mut browse_topology = false;
    if use_current {
        let shown = current_path
            .as_ref()
            .map(|p| p.display().to_string())
            .unwrap_or_else(|| "No structure open".into());
        widgets::caption(ui, &shown);
    } else {
        ui.horizontal(|ui| {
            ui.add(egui::TextEdit::singleline(&mut topology).desired_width(160.0));
            if widgets::button(ui, "", "Browse…", Variant::Secondary).clicked() {
                browse_topology = true;
            }
        });
    }
    ui.label("Trajectory file");
    let mut browse_trajectory = false;
    ui.horizontal(|ui| {
        ui.add(egui::TextEdit::singleline(&mut trajectory).desired_width(160.0));
        if widgets::button(ui, "", "Browse…", Variant::Secondary).clicked() {
            browse_trajectory = true;
        }
    });
    if browse_topology {
        if let Some(p) = rfd::FileDialog::new()
            .set_title("Topology file")
            .pick_file()
        {
            topology = p.display().to_string();
        }
    }
    if browse_trajectory {
        if let Some(p) = rfd::FileDialog::new()
            .set_title("Trajectory file")
            .pick_file()
        {
            trajectory = p.display().to_string();
        }
    }
    ui.data_mut(|d| d.insert_temp(id, (use_current, topology.clone(), trajectory.clone())));
    ui.add_space(crate::theme::space::GAP);
    if let Some(submit) = widgets::form_actions(ui, "Add trajectory") {
        if submit {
            let trajectory = trajectory.trim();
            let topo = if use_current {
                current_path.map(|p| p.display().to_string())
            } else {
                Some(topology.clone())
            };
            let trajectory_arg = vv_scene::script::quote_arg(trajectory);
            match (use_current.then(|| app.current()).flatten(), topo) {
                (Some(id), _) if !trajectory.is_empty() => {
                    app.run_command_logged(&format!(
                        "attachtrajectory {trajectory_arg} {}",
                        id.to_raw()
                    ));
                }
                (_, Some(t)) if !trajectory.is_empty() => {
                    let t = vv_scene::script::quote_arg(&t);
                    let trajectory = trajectory_arg;
                    app.run_command_logged(&format!("loadtrajectory {t} {trajectory}"));
                }
                _ => app.fail("choose a topology and a trajectory file".into()),
            }
        }
        ui.data_mut(|d| d.remove::<(bool, String, String)>(id));
        ui.close();
    }
}

/// A path field with a Browse… button opening the system save dialog
/// at the field's folder.
fn save_path_field(ui: &mut Ui, path: &mut String, extensions: &[&str]) {
    ui.horizontal(|ui| {
        let browse = widgets::button(ui, "", "Browse…", Variant::Secondary);
        ui.add(egui::TextEdit::singleline(path).desired_width(ui.available_width()));
        if browse.clicked() {
            let current = std::path::PathBuf::from(path.as_str());
            let mut dialog = rfd::FileDialog::new().add_filter("Image", extensions);
            if let Some(dir) = current.parent().filter(|d| d.is_dir()) {
                dialog = dialog.set_directory(dir);
            }
            if let Some(name) = current.file_name() {
                dialog = dialog.set_file_name(name.to_string_lossy());
            }
            if let Some(picked) = dialog.save_file() {
                *path = picked.display().to_string();
            }
        }
    });
}

/// File ▸ Export ▸ Render's form (UX_PATTERNS "Create"): size, samples,
/// transparent background. `Save to` is a plain typed path, matching
/// the existing Export ▸ Image… dialog (no native save picker there
/// either).
/// `file_render_popover`'s own form width/height, or its default (the
/// current viewport size) before it's ever been opened -- read from its
/// egui temp data, so `render_frame_overlay` can preview the same size
/// without a live borrow of that popover.
fn render_form_size(ctx: &egui::Context, renderer: &vv_render::Renderer) -> (u32, u32) {
    ctx.data(|d| d.get_temp::<RenderForm>(egui::Id::new("file-render-form")))
        .map(|f| (f.width, f.height))
        .unwrap_or_else(|| renderer.size())
}

#[derive(Clone)]
struct RenderForm {
    path: String,
    width: u32,
    height: u32,
    samples: u32,
    transparent: bool,
    quality: usize,
    aspect_linked: bool,
}

const MAX_RENDER_SIDE: u32 = 16384;

/// The size after the user edited one side of `before`: with the aspect
/// linked, the other side follows so `before`'s ratio holds (rounded, at
/// least 1, at most `MAX_RENDER_SIDE`).
fn linked_size(before: (u32, u32), edited: (u32, u32)) -> (u32, u32) {
    let follow = |value: u32, from: u32, to: u32| {
        let scaled = (value as f64 * to as f64 / from.max(1) as f64).round();
        scaled.clamp(1.0, MAX_RENDER_SIDE as f64) as u32
    };
    if edited.0 != before.0 {
        (edited.0, follow(edited.0, before.0, before.1))
    } else if edited.1 != before.1 {
        (follow(edited.1, before.1, before.0), edited.1)
    } else {
        edited
    }
}

fn render_size_row(ui: &mut Ui, form: &mut RenderForm) {
    use egui_phosphor::regular as icon;
    let before = (form.width, form.height);
    ui.horizontal(|ui| {
        ui.label("Size");
        ui.add(egui::DragValue::new(&mut form.width).range(1..=MAX_RENDER_SIDE));
        let glyph = if form.aspect_linked {
            icon::LINK
        } else {
            icon::LINK_BREAK
        };
        if widgets::icon_button(ui, glyph, form.aspect_linked)
            .on_hover_text("Keep the aspect ratio")
            .clicked()
        {
            form.aspect_linked = !form.aspect_linked;
        }
        ui.add(egui::DragValue::new(&mut form.height).range(1..=MAX_RENDER_SIDE));
    });
    if form.aspect_linked {
        (form.width, form.height) = linked_size(before, (form.width, form.height));
    }
}

pub(crate) fn file_render_popover(app: &mut AppUi<'_>, ui: &mut Ui) {
    ui.set_max_width(240.0);
    let id = egui::Id::new("file-render-form");
    let mut form = ui
        .data_mut(|d| d.get_temp::<RenderForm>(id))
        .unwrap_or_else(|| {
            let (width, height) = render_form_size(ui.ctx(), app.renderer);
            let default = crate::render::Quality::default();
            let quality = crate::render::Quality::ALL
                .iter()
                .position(|&q| q == default)
                .expect("Quality::default is in Quality::ALL");
            RenderForm {
                path: crate::layout::export_dir()
                    .join("render.png")
                    .display()
                    .to_string(),
                width,
                height,
                samples: default.default_samples(),
                transparent: false,
                quality,
                aspect_linked: true,
            }
        });
    ui.label("Save to");
    save_path_field(ui, &mut form.path, &["png", "jpg"]);
    render_size_row(ui, &mut form);
    ui.horizontal(|ui| {
        ui.label("Samples");
        ui.add(egui::DragValue::new(&mut form.samples).range(1..=4096));
    });
    ui.label("Quality");
    if let Some(picked) = widgets::segmented(ui, &["Draft", "High", "Ultra"], Some(form.quality)) {
        form.quality = picked;
        form.samples = crate::render::Quality::ALL[picked].default_samples();
    }
    if widgets::switch(ui, form.transparent, "Transparent background").clicked() {
        form.transparent = !form.transparent;
    }
    ui.data_mut(|d| d.insert_temp(id, form.clone()));
    ui.add_space(crate::theme::space::GAP);
    if let Some(submit) = widgets::form_actions(ui, "Render") {
        if submit && !form.path.trim().is_empty() {
            let t = if form.transparent { " transparent" } else { "" };
            let q = crate::render::Quality::ALL[form.quality].name();
            let path = vv_scene::script::quote_arg(form.path.trim());
            let RenderForm {
                width,
                height,
                samples,
                ..
            } = form;
            app.run_command_logged(&format!("render {path} {width}x{height} {samples} {q}{t}"));
        }
        ui.data_mut(|d| d.remove::<RenderForm>(id));
        ui.close();
    }
}

/// Analyze ▸ Caption's form (UX_PATTERNS "Create"): name, fractional
/// viewport position and text. Dispatched directly (`Command::
/// SetCaption`), not run as a `caption` command line, since the text can
/// hold anything a hand-built command string would have to escape.
pub(crate) fn analyze_caption_popover(app: &mut AppUi<'_>, ui: &mut Ui) {
    ui.set_max_width(220.0);
    let id = egui::Id::new("analyze-caption-form");
    let (mut name, mut x, mut y, mut text) = ui
        .data_mut(|d| d.get_temp::<(String, f32, f32, String)>(id))
        .unwrap_or_else(|| ("caption1".to_owned(), 0.5, 0.1, String::new()));
    ui.label("Name");
    ui.add(egui::TextEdit::singleline(&mut name).desired_width(f32::INFINITY));
    ui.horizontal(|ui| {
        ui.label("Position");
        ui.add(egui::DragValue::new(&mut x).range(0.0..=1.0).speed(0.01));
        ui.add(egui::DragValue::new(&mut y).range(0.0..=1.0).speed(0.01));
    });
    ui.label("Text");
    ui.add(egui::TextEdit::singleline(&mut text).desired_width(f32::INFINITY));
    ui.data_mut(|d| d.insert_temp(id, (name.clone(), x, y, text.clone())));
    ui.add_space(crate::theme::space::GAP);
    let Some(submit) = widgets::form_actions(ui, "Add caption") else {
        return;
    };
    if !submit {
        ui.data_mut(|d| d.remove::<(String, f32, f32, String)>(id));
        ui.close();
        return;
    }
    let (trimmed_name, trimmed_text) = (name.trim(), text.trim());
    if trimmed_name.is_empty() || trimmed_text.is_empty() {
        app.fail("a caption needs a name and text".into());
        return;
    }
    app.dispatch(Command::SetCaption {
        caption: vv_scene::Caption {
            name: trimmed_name.to_owned(),
            x,
            y,
            text: trimmed_text.to_owned(),
        },
    });
    app.undoable(&format!("Added the caption {trimmed_name}."));
    ui.data_mut(|d| d.remove::<(String, f32, f32, String)>(id));
    ui.close();
}

/// Analyze ▸ Contacts' form (UX_PATTERNS "Create"): cutoff and two
/// selections, run as `contacts` (its own well-tested pipeline) rather
/// than reimplemented here; `run_command_logged` stashes its output as
/// the Inspector's result card.
pub(crate) fn analyze_contacts_popover(app: &mut AppUi<'_>, ui: &mut Ui) {
    ui.set_max_width(220.0);
    let id = egui::Id::new("analyze-contacts-form");
    let (mut cutoff, mut a, mut b) = ui
        .data_mut(|d| d.get_temp::<(f32, String, String)>(id))
        .unwrap_or_else(|| (4.5, String::new(), String::new()));
    ui.label("Cutoff");
    ui.add(
        egui::DragValue::new(&mut cutoff)
            .range(0.1..=100.0)
            .speed(0.1)
            .suffix(" Å"),
    );
    ui.label("Selection A");
    ui.add(
        egui::TextEdit::singleline(&mut a)
            .hint_text("chain A")
            .desired_width(f32::INFINITY),
    );
    ui.label("Selection B");
    ui.add(
        egui::TextEdit::singleline(&mut b)
            .hint_text("chain B")
            .desired_width(f32::INFINITY),
    );
    ui.data_mut(|d| d.insert_temp(id, (cutoff, a.clone(), b.clone())));
    ui.add_space(crate::theme::space::GAP);
    let Some(submit) = widgets::form_actions(ui, "Find contacts") else {
        return;
    };
    if !submit {
        ui.data_mut(|d| d.remove::<(f32, String, String)>(id));
        ui.close();
        return;
    }
    let (trimmed_a, trimmed_b) = (a.trim(), b.trim());
    if trimmed_a.is_empty() || trimmed_b.is_empty() {
        app.fail("choose two selections".into());
        return;
    }
    app.run_command_logged(&format!("contacts {cutoff} {trimmed_a} | {trimmed_b}"));
    ui.data_mut(|d| d.remove::<(f32, String, String)>(id));
    ui.close();
}

/// A first guess at a gradient's top colour: the background lifted
/// toward white, so the fade reads as light from above.
pub(crate) fn gradient_top(bottom: wgpu::Color) -> wgpu::Color {
    let lift = |c: f64| c + (1.0 - c) * 0.35;
    wgpu::Color {
        r: lift(bottom.r),
        g: lift(bottom.g),
        b: lift(bottom.b),
        a: bottom.a,
    }
}

#[cfg(test)]
mod path_tests {
    use super::clean_path;

    #[test]
    fn pasted_paths_lose_their_quotes_and_tildes() {
        let p = |s: &str| clean_path(s).map(|p| p.display().to_string());
        assert_eq!(p(r"  C:\data\x.cif  "), Some(r"C:\data\x.cif".into()));
        assert_eq!(
            p(r#""C:\My Files\4hhb.cif""#),
            Some(r"C:\My Files\4hhb.cif".into())
        );
        assert_eq!(p("'/tmp/a b.pdb'"), Some("/tmp/a b.pdb".into()));
        assert_eq!(p(""), None);
        assert_eq!(p(r#"  ""  "#), None);
        let home = std::env::var("USERPROFILE")
            .or_else(|_| std::env::var("HOME"))
            .unwrap();
        let expanded = p("~/Downloads/x.pdb").unwrap();
        assert!(
            expanded.starts_with(&home) && expanded.ends_with("x.pdb"),
            "{expanded}"
        );
    }
}

#[cfg(test)]
mod cube_anim_tests {
    use super::{Camera, CubeAnim, CUBE_ANIM_SECS};
    use std::time::{Duration, Instant};

    #[test]
    fn steps_partway_then_reaches_exactly_the_target_face() {
        let camera = Camera::framing(glam::Vec3::ZERO, 1.0);
        let mut anim = CubeAnim::to_face(&camera, "+X").expect("+X is a real face");

        anim.start = Instant::now() - Duration::from_secs_f32(CUBE_ANIM_SECS * 0.5);
        let (mid, _, _, finished) = anim.step();
        assert!(!finished);
        assert_ne!(mid, anim.from, "halfway should have moved");
        assert_ne!(mid, anim.to, "halfway shouldn't already be there");

        anim.start = Instant::now() - Duration::from_secs_f32(CUBE_ANIM_SECS * 2.0);
        let (end, _, _, finished) = anim.step();
        assert!(finished);
        assert_eq!(end, anim.to);
    }

    #[test]
    fn to_camera_also_eases_distance() {
        let from = Camera::framing(glam::Vec3::ZERO, 1.0);
        let to = Camera::framing(glam::Vec3::ZERO, 10.0);
        let mut anim = CubeAnim::to_camera(&from, &to);
        anim.start = Instant::now() - Duration::from_secs_f32(CUBE_ANIM_SECS * 2.0);
        let (_, _, distance, finished) = anim.step();
        assert!(finished);
        assert_eq!(distance, to.distance);
    }

    #[test]
    fn an_unknown_face_yields_no_animation() {
        let camera = Camera::framing(glam::Vec3::ZERO, 1.0);
        assert!(CubeAnim::to_face(&camera, "nonsense").is_none());
    }
}

#[cfg(test)]
mod trackball_tests {
    use super::{trackball_rotate, Camera, Vec2};

    /// The molecule's point that faced the camera moves the same way on
    /// screen as the drag (egui y is down, so an upward drag is `-y`).
    #[test]
    fn the_near_point_follows_the_drag_direction() {
        let on_screen = |drag: Vec2| {
            let mut camera = Camera::framing(glam::Vec3::ZERO, 1.0);
            let near = camera.orientation * glam::Vec3::Z;
            trackball_rotate(&mut camera, drag);
            glam::Vec2::new(
                near.dot(camera.orientation * glam::Vec3::X),
                near.dot(camera.orientation * glam::Vec3::Y),
            )
        };
        let right = on_screen(Vec2::new(50.0, 0.0));
        assert!(right.x > 0.05 && right.y.abs() < 1e-4, "{right:?}");
        let up = on_screen(Vec2::new(0.0, -50.0));
        assert!(up.y > 0.05 && up.x.abs() < 1e-4, "{up:?}");
    }

    #[test]
    fn a_counter_clockwise_drag_is_a_positive_roll() {
        let quarter = super::roll_angle(Vec2::new(10.0, 0.0), Vec2::new(0.0, -10.0));
        assert!(
            (quarter - std::f32::consts::FRAC_PI_2).abs() < 1e-5,
            "{quarter}"
        );
        assert_eq!(
            super::roll_angle(Vec2::new(1.0, 0.0), Vec2::new(0.0, 1.0)),
            0.0
        );
    }

    #[test]
    fn a_zero_drag_leaves_the_orientation_unchanged() {
        let mut camera = Camera::framing(glam::Vec3::ZERO, 1.0);
        let before = camera.orientation;
        trackball_rotate(&mut camera, Vec2::ZERO);
        assert_eq!(camera.orientation, before);
    }
}

#[cfg(test)]
mod linked_size_tests {
    use super::linked_size;

    #[test]
    fn editing_width_scales_height_to_the_old_ratio() {
        assert_eq!(linked_size((1920, 1080), (960, 1080)), (960, 540));
    }

    #[test]
    fn editing_height_scales_width_to_the_old_ratio() {
        assert_eq!(linked_size((1920, 1080), (1920, 540)), (960, 540));
    }

    #[test]
    fn results_round_and_stay_in_range() {
        assert_eq!(linked_size((100, 33), (50, 33)), (50, 17));
        assert_eq!(linked_size((16384, 1), (1, 1)), (1, 1));
        assert_eq!(linked_size((1, 16384), (16384, 16384)), (16384, 16384));
    }

    #[test]
    fn an_unchanged_size_is_left_alone() {
        assert_eq!(linked_size((640, 480), (640, 480)), (640, 480));
    }
}
