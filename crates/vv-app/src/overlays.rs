//! Viewport overlays drawn with egui's painter over the rendered image:
//! the view cube (click a face to look at the scene from that side) and a
//! live scale bar. Interface only -- never in screenshots -- and no GPU
//! pass of their own.

use egui::{Color32, Pos2, Rect, Stroke, Ui, Vec2};
use glam::{Quat, Vec3};

/// The six faces: outward normal, then two in-face axes (right, up as
/// seen from outside), and a label.
const FACES: [(Vec3, Vec3, Vec3, &str); 6] = [
    (Vec3::Z, Vec3::X, Vec3::Y, "+Z"),
    (Vec3::NEG_Z, Vec3::NEG_X, Vec3::Y, "-Z"),
    (Vec3::X, Vec3::NEG_Z, Vec3::Y, "+X"),
    (Vec3::NEG_X, Vec3::Z, Vec3::Y, "-X"),
    (Vec3::Y, Vec3::X, Vec3::NEG_Z, "+Y"),
    (Vec3::NEG_Y, Vec3::X, Vec3::Z, "-Y"),
];

/// Half the cube's on-screen size, in points.
const CUBE_HALF: f32 = 26.0;

/// The axes drawn through the cube, in the red-green-blue order every
/// 3D tool uses for x, y, z (scene data, not interface chrome).
const AXES: [(Vec3, &str, Color32); 3] = [
    (Vec3::X, "X", Color32::from_rgb(0xE5, 0x48, 0x4D)),
    (Vec3::Y, "Y", Color32::from_rgb(0x46, 0xA7, 0x58)),
    (Vec3::Z, "Z", Color32::from_rgb(0x3E, 0x8E, 0xF7)),
];

/// Draws the view cube in the viewport's top-right corner; returns the
/// label of the face clicked (`view face` takes it).
pub fn view_cube(
    ui: &Ui,
    viewport: Rect,
    orientation: Quat,
    tokens: &crate::theme::Tokens,
) -> Option<&'static str> {
    let centre = Pos2::new(
        viewport.right() - CUBE_HALF * 2.2,
        viewport.top() + CUBE_HALF * 2.2,
    );
    let to_view = orientation.inverse();
    let project = |p: Vec3| {
        let v = to_view * p;
        centre + Vec2::new(v.x, -v.y) * CUBE_HALF * 0.8
    };
    let area = Rect::from_center_size(centre, Vec2::splat(CUBE_HALF * 3.2));
    let response = ui.interact(area, ui.id().with("view cube"), egui::Sense::click());
    let pointer = response.hover_pos();
    let painter = ui.painter_at(viewport);
    // Axes behind the faces: each shows where it leaves the cube, so the
    // orientation reads even looking straight at one face.
    for (axis, name, color) in AXES {
        let tip = project(axis * 1.9);
        painter.line_segment([project(Vec3::ZERO), tip], Stroke::new(2.0, color));
        painter.text(
            tip,
            egui::Align2::CENTER_CENTER,
            name,
            egui::FontId::proportional(crate::theme::text::CAPTION),
            color,
        );
    }
    let mut clicked = None;
    for (normal, right, up, label) in FACES {
        let facing = (to_view * normal).z;
        if facing <= 0.05 {
            continue;
        }
        let corners: Vec<Pos2> = [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)]
            .iter()
            .map(|&(a, b)| project(normal + right * a + up * b))
            .collect();
        let hovered = pointer.is_some_and(|p| inside(&corners, p));
        let fill = if hovered {
            tokens.primary
        } else {
            tokens.surface_raised.gamma_multiply(0.55 + 0.45 * facing)
        };
        painter.add(egui::Shape::convex_polygon(
            corners.clone(),
            fill,
            Stroke::new(1.0, tokens.border),
        ));
        let text = if hovered {
            tokens.on_primary
        } else {
            tokens.text
        };
        painter.text(
            project(normal),
            egui::Align2::CENTER_CENTER,
            label,
            egui::FontId::proportional(crate::theme::text::CAPTION),
            text.gamma_multiply(facing.sqrt()),
        );
        // Left or middle click turns to face `label`; both act the
        // same, so the cube works on trackpads with no distinct middle
        // button and on mice where middle is the natural "look here".
        if hovered && (response.clicked() || response.clicked_by(egui::PointerButton::Middle)) {
            clicked = Some(label);
        }
    }
    clicked
}

fn inside(polygon: &[Pos2], p: Pos2) -> bool {
    let n = polygon.len();
    let mut sign = 0.0f32;
    for i in 0..n {
        let (a, b) = (polygon[i], polygon[(i + 1) % n]);
        let cross = (b - a).x * (p - a).y - (b - a).y * (p - a).x;
        if cross != 0.0 {
            if sign != 0.0 && cross.signum() != sign {
                return false;
            }
            sign = cross.signum();
        }
    }
    true
}

/// The longest 1, 2 or 5 x 10^n Angstrom that fits in `max` Angstrom.
pub fn nice_length(max: f32) -> f32 {
    let power = 10f32.powf(max.log10().floor());
    [5.0, 2.0, 1.0]
        .into_iter()
        .map(|m| m * power)
        .find(|&l| l <= max)
        .unwrap_or(power)
}

/// Item 5a: space the scale bar needs at the viewport's bottom-left --
/// its ticks (±5 either side of the line, 22 px up) and label (another
/// ~16 px above that), plus a small margin. `viewport_ui` shrinks the
/// toast's rect by this much so a toast never covers it.
pub const SCALE_BAR_CLEARANCE: f32 = 48.0;

/// A scale bar in the viewport's bottom-left corner. `angstrom_per_point`
/// is the world size of one interface point at the rotation centre's
/// depth; the text contrasts with the viewport `background`.
pub fn scale_bar(ui: &Ui, viewport: Rect, angstrom_per_point: f32, background: Color32) {
    if !(angstrom_per_point.is_finite() && angstrom_per_point > 0.0) {
        return;
    }
    let length = nice_length(140.0 * angstrom_per_point);
    let width = length / angstrom_per_point;
    let light_background =
        background.r() as u32 + background.g() as u32 + background.b() as u32 > 3 * 128;
    let ink = crate::theme::Tokens::of(if light_background {
        crate::theme::ThemeMode::Light
    } else {
        crate::theme::ThemeMode::Dark
    })
    .text;
    let left = Pos2::new(viewport.left() + 16.0, viewport.bottom() - 22.0);
    let right = left + Vec2::new(width, 0.0);
    let painter = ui.painter_at(viewport);
    let stroke = Stroke::new(2.0, ink);
    painter.line_segment([left, right], stroke);
    for end in [left, right] {
        painter.line_segment(
            [end - Vec2::new(0.0, 5.0), end + Vec2::new(0.0, 5.0)],
            stroke,
        );
    }
    // Angstrom is the unit of the scene; nm only past 10 nm.
    let label = if length >= 100.0 {
        format!("{} nm", length / 10.0)
    } else {
        format!("{length} Å")
    };
    painter.text(
        left + Vec2::new(width * 0.5, -8.0),
        egui::Align2::CENTER_BOTTOM,
        label,
        egui::FontId::proportional(crate::theme::text::CAPTION),
        ink,
    );
}

/// The render's `aspect` (width/height) as a rectangle centred in
/// `viewport`, letterboxed on whichever axis is relatively taller; the
/// area outside it is dimmed so framing matches the exported image
/// (Studio ▸ Frame preview / `framing`, or while File ▸ Export ▸
/// Render's popover is open). Pure overlay -- never in a screenshot or
/// render, which reads the camera and settings directly, not this rect.
pub fn render_frame(ui: &Ui, viewport: Rect, aspect: f32) {
    if !(aspect.is_finite() && aspect > 0.0) {
        return;
    }
    let frame = fit_aspect(viewport, aspect);
    let painter = ui.painter_at(viewport);
    let dim = Color32::from_black_alpha(140);
    // The four bands outside `frame`, clockwise from the top.
    for band in [
        Rect::from_min_max(
            viewport.left_top(),
            Pos2::new(viewport.right(), frame.top()),
        ),
        Rect::from_min_max(
            Pos2::new(frame.right(), frame.top()),
            viewport.right_bottom(),
        ),
        Rect::from_min_max(
            Pos2::new(viewport.left(), frame.bottom()),
            viewport.right_bottom(),
        ),
        Rect::from_min_max(
            viewport.left_top(),
            Pos2::new(frame.left(), viewport.bottom()),
        ),
    ] {
        painter.rect_filled(band, 0.0, dim);
    }
    painter.rect_stroke(
        frame,
        0.0,
        Stroke::new(1.0, Color32::WHITE),
        egui::StrokeKind::Inside,
    );
}

/// The largest `aspect`-shaped rectangle centred in and fully inside
/// `viewport`.
fn fit_aspect(viewport: Rect, aspect: f32) -> Rect {
    let (w, h) = (viewport.width(), viewport.height());
    let size = if w / h > aspect {
        Vec2::new(h * aspect, h)
    } else {
        Vec2::new(w, w / aspect)
    };
    Rect::from_center_size(viewport.center(), size)
}

/// The outward normal and screen-up of the face labelled `face` ("+X"...).
pub fn face(face: &str) -> Option<(Vec3, Vec3)> {
    FACES
        .iter()
        .find(|f| f.3.eq_ignore_ascii_case(face))
        .map(|f| (f.0, f.2))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fit_aspect_letterboxes_the_narrower_axis() {
        let viewport = Rect::from_min_size(Pos2::ZERO, Vec2::new(800.0, 600.0));
        // Wider than the viewport (2:1 in a 4:3 area): full width, bars
        // top and bottom.
        let wide = fit_aspect(viewport, 2.0);
        assert_eq!(wide.width(), 800.0);
        assert_eq!(wide.height(), 400.0);
        assert_eq!(wide.center(), viewport.center());
        // Taller than the viewport (1:2): full height, bars left/right.
        let tall = fit_aspect(viewport, 0.5);
        assert_eq!(tall.height(), 600.0);
        assert_eq!(tall.width(), 300.0);
        assert_eq!(tall.center(), viewport.center());
        // Same aspect: the whole viewport, no bars.
        let same = fit_aspect(viewport, 800.0 / 600.0);
        assert!((same.width() - 800.0).abs() < 1e-3);
        assert!((same.height() - 600.0).abs() < 1e-3);
    }

    #[test]
    fn nice_lengths_are_1_2_5_steps_that_fit() {
        assert_eq!(nice_length(9.9), 5.0);
        assert_eq!(nice_length(10.0), 10.0);
        assert_eq!(nice_length(47.0), 20.0);
        assert_eq!(nice_length(140.0), 100.0);
        assert_eq!(nice_length(0.3), 0.2);
    }

    #[test]
    fn the_default_view_shows_the_plus_z_face_only_head_on() {
        let to_view = Quat::IDENTITY.inverse();
        let visible: Vec<&str> = FACES
            .iter()
            .filter(|(n, ..)| (to_view * *n).z > 0.05)
            .map(|f| f.3)
            .collect();
        assert_eq!(visible, ["+Z"]);
    }

    #[test]
    fn looking_from_a_face_shows_that_face_head_on_upright() {
        for (normal, right, up, label) in FACES {
            let mut camera = vv_render::Camera::framing(Vec3::ZERO, 10.0);
            camera.orbit(0.7, -0.4);
            let (from, face_up) = face(label).unwrap();
            camera.look_from(from, face_up);
            let o = camera.orientation;
            assert!((o * Vec3::Z).distance(normal) < 1e-5, "{label}");
            assert!((o * Vec3::Y).distance(up) < 1e-5, "{label}");
            assert!((o * Vec3::X).distance(right) < 1e-5, "{label}");
        }
    }
}
