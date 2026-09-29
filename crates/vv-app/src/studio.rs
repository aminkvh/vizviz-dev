//! Studio mode: the light list (`StudioLight`/`LightRig`, packed into
//! `vv_render::Lighting`), saved camera views (`SavedView`), and the
//! Studio ribbon tab's popovers.

use egui::{Color32, Ui};
use vv_render::{Camera, LightingPreset, Projection};

use crate::ui::AppUi;
use crate::widgets::{self, Variant};

pub const MAX_LIGHTS: usize = vv_render::style::MAX_LIGHTS;

/// One light, as the scene manager edits it: on/off, colour (sRGB, as
/// picked), intensity and direction in degrees -- azimuth right of the
/// view axis, elevation above it, both full range so a light can sit
/// behind the subject (a rim light) or below it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StudioLight {
    pub enabled: bool,
    pub color: Color32,
    pub intensity: f32,
    pub azimuth: f32,
    pub elevation: f32,
}

impl StudioLight {
    /// Filler rows offered by "+" for the list's 3rd/4th slot: off, back
    /// and to a side, so adding a light doesn't start unaimed.
    const EXTRA: [StudioLight; 2] = [
        StudioLight {
            enabled: false,
            color: Color32::WHITE,
            intensity: 0.6,
            azimuth: -130.0,
            elevation: 20.0,
        },
        StudioLight {
            enabled: false,
            color: Color32::WHITE,
            intensity: 0.5,
            azimuth: 130.0,
            elevation: -10.0,
        },
    ];

    fn of_gpu(light: vv_render::Light) -> Self {
        let [x, y, z] = light.dir;
        Self {
            enabled: true,
            color: Color32::WHITE,
            intensity: light.intensity,
            azimuth: x.atan2(z).to_degrees(),
            elevation: y.clamp(-1.0, 1.0).asin().to_degrees(),
        }
    }

    fn dir(self) -> [f32; 3] {
        let (az, el) = (self.azimuth.to_radians(), self.elevation.to_radians());
        [el.cos() * az.sin(), el.sin(), el.cos() * az.cos()]
    }

    fn to_gpu(self) -> vv_render::Light {
        vv_render::Light {
            dir: self.dir(),
            intensity: self.intensity,
            color: srgb_to_linear3(self.color),
            _pad: 0.0,
        }
    }

    pub fn name(index: usize) -> &'static str {
        match index {
            0 => "Key",
            1 => "Fill",
            2 => "Light 3",
            _ => "Light 4",
        }
    }
}

fn srgb_to_linear3(c: Color32) -> [f32; 3] {
    [c.r(), c.g(), c.b()].map(|v| vv_render::style::srgb_to_linear(v as f32 / 255.0))
}

/// The Studio light list: up to [`MAX_LIGHTS`] [`StudioLight`]s (`count`
/// of them in the list -- "+"/remove; each also has its own on/off
/// switch) plus ambient. A lighting preset fills it ([`LightRig::of`]);
/// [`LightRig::apply`] packs it into the GPU `Lighting` uniform.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LightRig {
    pub lights: [StudioLight; MAX_LIGHTS],
    /// How many rows the list shows (0..=[`MAX_LIGHTS`]); a light past
    /// this keeps its own settings (re-adding it remembers them) but
    /// never lights the scene.
    pub count: usize,
    pub ambient: f32,
    pub ambient_color: Color32,
}

impl LightRig {
    pub fn of(preset: LightingPreset) -> Self {
        let l = preset.lighting();
        Self {
            lights: [
                StudioLight::of_gpu(l.lights[0]),
                StudioLight::of_gpu(l.lights[1]),
                StudioLight::EXTRA[0],
                StudioLight::EXTRA[1],
            ],
            count: 2,
            ambient: l.ambient,
            ambient_color: Color32::WHITE,
        }
    }

    /// Packs the enabled, in-list lights densely into `lighting` from
    /// index 0, so shaders never branch on a per-light flag -- only on
    /// `Lighting::count`.
    pub(crate) fn apply(&self, lighting: &mut vv_render::Lighting) {
        let mut n = 0;
        for light in self.lights.iter().take(self.count).filter(|l| l.enabled) {
            lighting.lights[n] = light.to_gpu();
            n += 1;
        }
        for slot in &mut lighting.lights[n..] {
            *slot = vv_render::Light::default();
        }
        lighting.count = n as f32;
        lighting.ambient = self.ambient;
        lighting.ambient_color = srgb_to_linear3(self.ambient_color);
    }

    /// Adds the next slot (up to [`MAX_LIGHTS`]), enabled: "+" in the
    /// Studio popover.
    pub fn add(&mut self) {
        if self.count < MAX_LIGHTS {
            self.lights[self.count].enabled = true;
            self.count += 1;
        }
    }

    /// Drops list row `i`; the rest shift down, keeping their own
    /// settings (the removed light's settings end up parked past the
    /// new count, so re-adding it with "+" brings them back).
    pub fn remove(&mut self, i: usize) {
        if i < self.count {
            self.lights[i..self.count].rotate_left(1);
            self.count -= 1;
        }
    }
}

/// A bookmarked camera position ("Save view"): everything [`Camera`]
/// needs to return to it, minus interaction-only state (`pivot`).
#[derive(Clone, Debug, PartialEq)]
pub struct SavedView {
    pub name: String,
    pub target: glam::Vec3,
    pub distance: f32,
    pub orientation: glam::Quat,
    pub fov_y: f32,
    pub near: f32,
    pub projection: Projection,
    pub scene_radius: f32,
    pub dolly: f32,
}

impl SavedView {
    pub fn capture(name: String, camera: &Camera) -> Self {
        Self {
            name,
            target: camera.target,
            distance: camera.distance,
            orientation: camera.orientation,
            fov_y: camera.fov_y,
            near: camera.near,
            projection: camera.projection,
            scene_radius: camera.scene_radius,
            dolly: camera.dolly,
        }
    }

    pub fn to_camera(&self) -> Camera {
        Camera {
            target: self.target,
            distance: self.distance,
            orientation: self.orientation,
            fov_y: self.fov_y,
            near: self.near,
            projection: self.projection,
            scene_radius: self.scene_radius,
            pivot: None,
            dolly: self.dolly,
        }
    }
}

/// The next "View N" not already taken by `existing`.
pub fn unique_view_name(existing: &[SavedView]) -> String {
    let mut n = existing.len() + 1;
    loop {
        let name = format!("View {n}");
        if !existing.iter().any(|v| v.name == name) {
            return name;
        }
        n += 1;
    }
}

/// The Studio/Look ▸ Lights ▾ rename-in-progress key (see
/// `crate::ui::inline_rename_ui`).
fn view_rename_key() -> egui::Id {
    egui::Id::new("saved-view-rename")
}

/// A colour preview, painted in place: not `widgets::color_picker`,
/// whose own popup would close this one (see its doc). The colour is
/// edited inline in the detail section below, for whichever light is
/// selected.
fn swatch_preview(ui: &mut Ui, color: Color32) {
    let (rect, _) = ui.allocate_exact_size(egui::Vec2::splat(16.0), egui::Sense::hover());
    if ui.is_rect_visible(rect) {
        let t = crate::theme::Tokens::current(ui.ctx());
        ui.painter().rect(
            rect,
            crate::theme::radius::CONTROL,
            color,
            egui::Stroke::new(1.0, t.border),
            egui::StrokeKind::Outside,
        );
    }
}

/// Look ▸ Lights ▾ and Studio ▸ Lights ▾ (same body -- see ribbon.rs's
/// `TABS` for why both exist): a preset picker, the light list (switch,
/// colour preview, intensity; add up to 4, remove), a direction gizmo
/// plus az/el sliders and a colour editor for whichever light is
/// selected, an ambient row, and Reset.
pub(crate) fn lights_popover(app: &mut AppUi<'_>, ui: &mut Ui) {
    ui.set_max_width(300.0);

    let names: Vec<&str> = LightingPreset::ALL
        .iter()
        .map(|p| crate::ui::preset_label(*p))
        .collect();
    if app.view.lights_custom {
        widgets::caption(ui, "Custom");
    }
    let now = (!app.view.lights_custom)
        .then(|| {
            LightingPreset::ALL
                .iter()
                .position(|p| *p == app.view.lighting)
        })
        .flatten();
    if let Some(i) = widgets::segmented(ui, &names, now) {
        app.set_lighting(LightingPreset::ALL[i]);
    }
    ui.add_space(crate::theme::space::GAP);
    ui.separator();

    // The list, gizmo, two colour editors and ambient easily run past a
    // short window's height; scrolling (not clipping) keeps Reset and
    // the ambient row reachable regardless.
    egui::ScrollArea::vertical()
        .max_height(420.0)
        .show(ui, |ui| lights_popover_body(app, ui));
}

fn lights_popover_body(app: &mut AppUi<'_>, ui: &mut Ui) {
    use egui_phosphor::regular as icon;
    let selected_key = egui::Id::new("studio-light-selected");
    let count = app.view.lights.count;
    let mut selected = ui
        .data(|d| d.get_temp::<usize>(selected_key))
        .unwrap_or(0)
        .min(count.saturating_sub(1));

    let mut moved = false;
    let snapshot = app.view.lights.lights;
    let mut toggled = None;
    let mut removed = None;
    let mut picked = None;
    for (i, light) in snapshot.iter().enumerate().take(count) {
        let label = format!("{}  ·  {:.2}", StudioLight::name(i), light.intensity);
        let row = widgets::list_row(ui, &label, selected == i, |ui| {
            if widgets::icon_button(ui, icon::TRASH, false)
                .on_hover_text("Remove")
                .clicked()
            {
                removed = Some(i);
            }
            swatch_preview(ui, light.color);
            if widgets::switch(ui, light.enabled, "").clicked() {
                toggled = Some(i);
            }
        });
        if row.clicked() {
            picked = Some(i);
        }
    }
    if let Some(i) = picked {
        selected = i;
    }
    if let Some(i) = toggled {
        app.view.lights.lights[i].enabled = !app.view.lights.lights[i].enabled;
        moved = true;
    }
    if let Some(i) = removed {
        app.view.lights.remove(i);
        selected = selected.min(app.view.lights.count.saturating_sub(1));
        moved = true;
    }
    ui.data_mut(|d| d.insert_temp(selected_key, selected));

    if app.view.lights.count == 0 {
        widgets::caption(ui, "No lights -- ambient only");
    } else {
        ui.add_space(crate::theme::space::GAP);
        ui.separator();
        widgets::caption(ui, StudioLight::name(selected));
        // Not `snapshot`: a remove this frame already shifted the real
        // list, and `snapshot` still holds the pre-shift rows.
        let dots: Vec<(f32, f32, Color32, bool)> = app.view.lights.lights[..app.view.lights.count]
            .iter()
            .map(|l| (l.azimuth, l.elevation, l.color, l.enabled))
            .collect();
        if let Some((az, el)) = widgets::direction_gizmo(ui, 140.0, &dots, selected) {
            let light = &mut app.view.lights.lights[selected];
            light.azimuth = az;
            light.elevation = el;
            moved = true;
        }
        let light = &mut app.view.lights.lights[selected];
        moved |= widgets::slider_moved(&widgets::slider(
            ui,
            "Azimuth",
            &mut light.azimuth,
            -180.0..=180.0,
            "°",
        ));
        moved |= widgets::slider_moved(&widgets::slider(
            ui,
            "Elevation",
            &mut light.elevation,
            -90.0..=90.0,
            "°",
        ));
        moved |= widgets::slider_moved(&widgets::slider(
            ui,
            "Intensity",
            &mut light.intensity,
            0.0..=2.0,
            "",
        ));
        let mut color = light.color;
        if widgets::color_editor(
            ui,
            egui::Id::new(("studio-light-color", selected)),
            &mut color,
        ) {
            app.view.lights.lights[selected].color = color;
            moved = true;
        }
    }

    ui.add_space(crate::theme::space::GAP);
    if app.view.lights.count < MAX_LIGHTS
        && widgets::button(ui, icon::PLUS, "Add light", Variant::Secondary).clicked()
    {
        app.view.lights.add();
        moved = true;
    }

    ui.add_space(crate::theme::space::GAP);
    widgets::section(ui, "Ambient");
    let mut ambient_color = app.view.lights.ambient_color;
    if widgets::color_editor(
        ui,
        egui::Id::new("studio-ambient-color"),
        &mut ambient_color,
    ) {
        app.view.lights.ambient_color = ambient_color;
        moved = true;
    }
    moved |= widgets::slider_moved(&widgets::slider(
        ui,
        "Ambient",
        &mut app.view.lights.ambient,
        0.0..=1.5,
        "",
    ));

    if moved {
        app.view.lights_custom = true;
    }
    if widgets::switch(ui, app.view.tonemap, "Tonemap highlights")
        .on_hover_text("Rolls bright highlights off instead of clipping them")
        .clicked()
    {
        app.view.tonemap = !app.view.tonemap;
    }
    ui.add_space(crate::theme::space::GAP);
    if widgets::button(ui, "", "Reset", Variant::Secondary).clicked() {
        app.view.lights = LightRig::of(app.view.lighting);
        app.view.lights_custom = false;
        app.view.tonemap = app.view.lighting.lighting().tonemap > 0.5;
    }
}

/// Studio ▸ Views ▾: "Save view" bookmarks the current camera; the list
/// below flies back to one on click (undoable, like any other camera
/// jump) and renames/deletes through the usual inline/⋯ patterns.
pub(crate) fn views_popover(app: &mut AppUi<'_>, ui: &mut Ui) {
    use egui_phosphor::regular as icon;
    ui.set_max_width(260.0);
    if widgets::button(ui, icon::BOOKMARK_SIMPLE, "Save view", Variant::Secondary).clicked() {
        let name = unique_view_name(&app.view.saved_views);
        app.view
            .saved_views
            .push(SavedView::capture(name, app.camera));
    }
    ui.add_space(crate::theme::space::GAP);

    let views = app.view.saved_views.clone();
    if views.is_empty() {
        widgets::caption(ui, "No saved views yet");
    }
    let key = view_rename_key();
    let renaming = crate::ui::renaming_id::<String>(ui, key);
    let mut jump = None;
    let mut rename_target = None;
    let mut delete = None;
    for (i, view) in views.iter().enumerate() {
        if renaming.as_deref() == Some(view.name.as_str()) {
            if let Some(new_name) =
                crate::ui::inline_rename_ui(ui, key, view.name.clone(), &view.name)
            {
                app.view.saved_views[i].name = new_name;
            }
            continue;
        }
        let items = [("Rename", false), ("Delete", true)];
        let mut picked = None;
        let clicked = widgets::list_row(ui, &view.name, false, |ui| {
            picked = widgets::menu_button(ui, &items).0;
        })
        .clicked();
        match picked.map(|i| items[i].0) {
            Some("Rename") => rename_target = Some(view.name.clone()),
            Some("Delete") => delete = Some(i),
            _ => {}
        }
        if clicked {
            jump = Some(i);
        }
    }
    if let Some(name) = rename_target {
        ui.data_mut(|d| d.insert_temp(key, (name.clone(), name)));
    }
    if let Some(i) = jump {
        let target = views[i].to_camera();
        app.mark_camera_jump();
        *app.cube_anim = Some(crate::ui::CubeAnim::to_camera(app.camera, &target));
    }
    if let Some(i) = delete {
        app.view.saved_views.remove(i);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_preset_maps_key_and_fill_onto_lights_0_and_1_enabled() {
        let rig = LightRig::of(LightingPreset::Default);
        let expected = LightingPreset::Default.lighting();
        assert_eq!(rig.count, 2);
        assert!(rig.lights[0].enabled && rig.lights[1].enabled);
        assert_eq!(rig.lights[0].intensity, expected.lights[0].intensity);
        assert_eq!(rig.lights[1].intensity, expected.lights[1].intensity);
        assert!(!rig.lights[2].enabled && !rig.lights[3].enabled);
    }

    /// Fill sits behind the scene in most presets (negative view-space
    /// z): the ±180° azimuth range must be able to express that, not
    /// just the ±90° a front-only light would need. `Default`'s fill is in
    /// front, so this uses `Full`, whose fill is behind the scene.
    #[test]
    fn fill_light_azimuth_can_exceed_90_degrees() {
        let rig = LightRig::of(LightingPreset::Full);
        assert!(
            rig.lights[1].azimuth.abs() > 90.0,
            "fill azimuth {}",
            rig.lights[1].azimuth
        );
    }

    #[test]
    fn azimuth_elevation_round_trip_through_gpu_direction() {
        for &(az, el) in &[(0.0, 0.0), (45.0, 20.0), (-170.0, -80.0), (179.0, 89.0)] {
            let light = StudioLight {
                enabled: true,
                color: Color32::WHITE,
                intensity: 1.0,
                azimuth: az,
                elevation: el,
            };
            let back = StudioLight::of_gpu(light.to_gpu());
            assert!((back.azimuth - az).abs() < 1e-2, "{az} -> {}", back.azimuth);
            assert!(
                (back.elevation - el).abs() < 1e-2,
                "{el} -> {}",
                back.elevation
            );
        }
    }

    #[test]
    fn apply_packs_only_enabled_in_list_lights_densely() {
        let mut rig = LightRig::of(LightingPreset::Default);
        rig.lights[0].enabled = false; // Key off
        rig.add(); // a 3rd light, enabled
        let mut lighting = vv_render::Lighting::default();
        rig.apply(&mut lighting);
        assert_eq!(lighting.count, 2.0, "Fill + the new 3rd light");
        assert_eq!(lighting.lights[0].intensity, rig.lights[1].intensity);
    }

    #[test]
    fn remove_shifts_the_rest_down_and_keeps_their_settings() {
        let mut rig = LightRig::of(LightingPreset::Default);
        rig.add();
        let (fill, third) = (rig.lights[1], rig.lights[2]);
        rig.remove(0); // drop Key
        assert_eq!(rig.count, 2);
        assert_eq!(rig.lights[0], fill);
        assert_eq!(rig.lights[1], third);
    }

    #[test]
    fn view_names_stay_unique() {
        let existing = vec![
            SavedView::capture("View 1".into(), &Camera::framing(glam::Vec3::ZERO, 1.0)),
            SavedView::capture("View 2".into(), &Camera::framing(glam::Vec3::ZERO, 1.0)),
        ];
        assert_eq!(unique_view_name(&existing), "View 3");
        assert_eq!(unique_view_name(&[]), "View 1");
    }

    #[test]
    fn saved_view_round_trips_camera_fields() {
        let mut camera = Camera::framing(glam::Vec3::new(1.0, 2.0, 3.0), 5.0);
        camera.orientation = glam::Quat::from_rotation_y(0.7);
        camera.dolly = 1.5;
        camera.projection = Projection::Orthographic;
        let saved = SavedView::capture("Front".into(), &camera);
        let back = saved.to_camera();
        assert_eq!(back.target, camera.target);
        assert_eq!(back.distance, camera.distance);
        assert_eq!(back.orientation, camera.orientation);
        assert_eq!(back.projection, camera.projection);
        assert_eq!(back.dolly, camera.dolly);
        assert_eq!(back.pivot, None);
    }
}
