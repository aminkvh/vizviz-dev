//! Orbit camera with a reversed-Z, infinite-far perspective projection.

use bytemuck::{Pod, Zeroable};
use glam::{Mat3, Mat4, Quat, Vec2, Vec3, Vec4};

/// How the camera maps view space to clip space. Orthographic drops
/// distance-based size falloff entirely (parallel rays, not converging
/// on a single eye point); the default. The difference isn't confined
/// to [`Camera::proj`]: every shader that reconstructs a per-fragment ray by casting from a single
/// eye point (`draw.wgsl`'s sphere/cylinder impostors, `gaussian_
/// surface.wgsl`/`skin_surface.wgsl`'s full-screen ray-casts) or sizes
/// something on screen by dividing by view depth (`cull.wgsl`'s LOD/
/// occlusion-margin math) is a perspective-specific formula and must
/// branch on this -- see `CameraUniform::projection` and each shader's
/// own comments at the sites that branch on it.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum Projection {
    Perspective,
    #[default]
    Orthographic,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Camera {
    pub target: Vec3,
    pub distance: f32,
    /// Maps camera axes to world: the eye sits along `orientation * +Z`
    /// from `target`, screen-up is `orientation * +Y`. A quaternion, not
    /// yaw/pitch, so rotation has no pole and never stops.
    pub orientation: Quat,
    pub fov_y: f32,
    pub near: f32,
    /// In orthographic mode, `fov_y` still sets the apparent framing (via
    /// [`Camera::proj`]'s half-height derivation from `distance`), so
    /// toggling projection with the same `distance` keeps the subject
    /// roughly the same size on screen instead of jumping.
    pub projection: Projection,
    /// Radius of what was framed (`framing`). Orthographic zoom scales
    /// the view (`distance` sets the half-height) but keeps the eye at
    /// least this far out beyond the target, so zooming in never moves
    /// the eye into the structure and slices it (it did on 8GLV).
    pub scene_radius: f32,
    /// Where rotation pivots, when set: orbiting swings the
    /// view around this point instead of `target`, so picking a center
    /// doesn't move the picture. `None` pivots on `target`.
    pub pivot: Option<Vec3>,
    /// How far (Angstrom, >= 0) the camera has moved forward: a dolly, as
    /// opposed to `zoom`. The scene nearer than the camera is cut away
    /// along a view-facing plane ([`Camera::near_cut`]) and capped.
    pub dolly: f32,
}

impl Camera {
    /// Frames a bounding sphere so it fits the vertical field of view.
    pub fn framing(center: Vec3, radius: f32) -> Self {
        let fov_y = 45f32.to_radians();
        let distance = radius / (fov_y * 0.5).sin() * 1.05;
        Self {
            target: center,
            distance,
            orientation: Quat::IDENTITY,
            fov_y,
            near: 0.1,
            projection: Projection::default(),
            scene_radius: radius,
            pivot: None,
            dolly: 0.0,
        }
    }

    /// The orientation for a yaw about world +Y and an elevation above
    /// the XZ plane (radians): how older sessions stored the camera.
    pub fn angles(yaw: f32, pitch: f32) -> Quat {
        Quat::from_rotation_y(yaw) * Quat::from_rotation_x(-pitch)
    }

    /// How far the eye sits from `target`: `distance` less the dolly in
    /// perspective (negative once it has passed the target); in
    /// orthographic, never inside the scene.
    /// Half the world height on screen at the target plane, in either
    /// projection (orthographic zoom sets it through `distance`).
    pub fn visible_half_height(&self) -> f32 {
        self.distance * (self.fov_y * 0.5).tan()
    }

    pub fn eye_distance(&self) -> f32 {
        match self.projection {
            Projection::Perspective => self.distance - self.dolly,
            Projection::Orthographic => self
                .distance
                .max(self.scene_radius * 2.0 + self.near * 10.0),
        }
    }

    /// The orthographic far plane: past the whole scene from the (pushed
    /// out) eye, with generous slack. Finite: no 1/z falloff to push to a
    /// limit.
    pub fn orthographic_far(&self) -> f32 {
        self.eye_distance() + self.scene_radius * 2.0 + self.distance * 10.0 + self.near * 1000.0
    }

    pub fn eye(&self) -> Vec3 {
        self.target + self.orientation * Vec3::Z * self.eye_distance()
    }

    /// Unit view direction, from the eye into the scene.
    pub fn forward(&self) -> Vec3 {
        self.orientation * Vec3::NEG_Z
    }

    /// Looks along `forward` rather than at `target`, which a dolly can
    /// carry the eye past.
    pub fn view(&self) -> Mat4 {
        Mat4::look_to_rh(self.eye(), self.forward(), self.orientation * Vec3::Y)
    }

    /// View depth of the plane the camera cuts the scene at, if any.
    /// Perspective always cuts just past the near plane, so an eye inside
    /// the scene sees capped cross-sections instead of the rasterizer's
    /// ragged near clip (the margin keeps cap fragments off that plane).
    /// Orthographic has no eye position to move, so its camera starts at
    /// the front of the framed scene and cuts only once dollied in.
    pub fn near_cut(&self) -> Option<f32> {
        match self.projection {
            Projection::Perspective => Some(self.near * 2.0),
            Projection::Orthographic => {
                (self.dolly > 0.0).then(|| self.eye_distance() - self.scene_radius + self.dolly)
            }
        }
    }

    /// [`Camera::near_cut`] as a world-space clip plane (`n.p + d < 0` cut
    /// away), in the form `RenderSettings::clip` takes.
    pub fn near_cut_plane(&self) -> Option<[f32; 4]> {
        let f = self.forward();
        self.near_cut()
            .map(|depth| f.extend(-f.dot(self.eye() + f * depth)).to_array())
    }

    /// Moves the camera `by` Angstrom forward (negative: back), never
    /// behind where it started.
    pub fn dolly_by(&mut self, by: f32) {
        self.dolly = (self.dolly + by).max(0.0);
    }

    /// Dollies for a wheel scroll of `points` (about 50 a notch; positive
    /// is forward), a notch being a twentieth of the framed radius.
    pub fn dolly_wheel(&mut self, points: f32) {
        self.dolly_by(points / 50.0 * (self.scene_radius * 0.05).max(0.5));
    }

    /// A world-space clip plane `[nx, ny, nz, d]` (`n.p + d < 0` cut away)
    /// in this camera's view space, as the shaders test it.
    pub fn view_clip(&self, [x, y, z, d]: [f32; 4]) -> [f32; 4] {
        let view = self.view();
        let n = view.transform_vector3(Vec3::new(x, y, z)).normalize();
        // The view is rigid: a point p maps to R p + t, so the plane's
        // offset moves by -n' . t.
        n.extend(d - n.dot(view.w_axis.truncate())).to_array()
    }

    pub fn proj(&self, aspect: f32) -> Mat4 {
        match self.projection {
            Projection::Perspective => {
                Mat4::perspective_infinite_reverse_rh(self.fov_y, aspect, self.near)
            }
            Projection::Orthographic => {
                // Half-extents at `distance` from `fov_y`, so switching
                // from perspective at the same `distance` keeps the
                // subject the same apparent size instead of jumping --
                // see the struct doc.
                let half_height = self.distance * (self.fov_y * 0.5).tan();
                let half_width = half_height * aspect;
                // No infinite-far trick for orthographic (there's no
                // 1/z falloff to push to a limit): a generous finite far
                // plane, scaled with `distance` so it still comfortably
                // contains the frustum after zooming.
                let far = self.orthographic_far();
                orthographic_reversed_rh(half_width, half_height, self.near, far)
            }
        }
    }

    /// Trackball rotation about the screen's own axes: `delta_yaw` about
    /// screen-up, `delta_pitch` about screen-right. No limit in any
    /// direction.
    pub fn orbit(&mut self, delta_yaw: f32, delta_pitch: f32) {
        self.turn(Quat::from_rotation_y(delta_yaw) * Quat::from_rotation_x(-delta_pitch));
    }

    /// Turns the scene by `angle` radians about a screen axis (x right,
    /// y up, z toward the viewer).
    pub fn rotate_scene(&mut self, axis: Vec3, angle: f32) {
        self.turn(Quat::from_axis_angle(axis.normalize(), -angle));
    }

    /// Puts the eye on the `from` side of the target with `up` screen-up
    /// (both unit and perpendicular), keeping a pivot where it is on screen.
    pub fn look_from(&mut self, from: Vec3, up: Vec3) {
        let goal = Quat::from_mat3(&Mat3::from_cols(up.cross(from), up, from));
        self.turn(self.orientation.inverse() * goal);
    }

    /// Applies `q` (in the camera's own frame) to the orientation.
    fn turn(&mut self, q: Quat) {
        let before = self.orientation;
        self.orientation = (self.orientation * q).normalize();
        if let Some(pivot) = self.pivot {
            // Swing the target around the pivot by the same rotation, so
            // the pivot stays put on screen.
            let delta = self.orientation * before.inverse();
            self.target = pivot + delta * (self.target - pivot);
        }
    }

    pub fn zoom(&mut self, factor: f32) {
        self.distance = (self.distance * factor).max(self.near * 2.0);
    }

    /// Slides `target` (the orbit pivot) across the view plane by
    /// `delta_x`/`delta_y` screen points, dragging the scene under the
    /// cursor rather than rotating around it. Scaled by `distance` so the
    /// same drag distance feels the same whether zoomed in or out.
    pub fn pan(&mut self, delta_x: f32, delta_y: f32) {
        let right = self.orientation * Vec3::X;
        let up = self.orientation * Vec3::Y;
        let scale = self.distance * 0.0015;
        self.target += (right * -delta_x + up * delta_y) * scale;
    }

    /// Camera-derived fields only; the renderer fills in the rest.
    pub fn uniform(&self, viewport: Vec2) -> CameraUniform {
        let view = self.view();
        let proj = self.proj(viewport.x / viewport.y);
        let view_proj = proj * view;
        CameraUniform {
            view_proj: view_proj.to_cols_array_2d(),
            view: view.to_cols_array_2d(),
            proj: proj.to_cols_array_2d(),
            planes: frustum_planes(view_proj),
            viewport: viewport.to_array(),
            quad_px_threshold: 1.0,
            // Pixels per world unit at view depth 1.
            proj_scale: proj.y_axis.y * viewport.y * 0.5,
            hiz_mip_count: 0,
            occlusion: 0,
            near: self.near,
            projection: match self.projection {
                Projection::Perspective => 0,
                Projection::Orthographic => 1,
            },
            clip: KEEP_ALL,
            cut: self
                .near_cut_plane()
                .map_or(KEEP_ALL, |plane| self.view_clip(plane)),
        }
    }
}

/// A right-handed orthographic projection, reversed-Z (near -> depth 1,
/// far -> depth 0, matching [`Camera::proj`]'s perspective branch and
/// `shaders/*.wgsl`'s shared convention): a plain linear remap since
/// orthographic has no 1/z divide, derived directly rather than composed
/// from glam's non-reversed `orthographic_rh` to keep the near/far
/// algebra in one visible place. `clip.w` is always `1.0` (true
/// orthographic: no perspective division).
fn orthographic_reversed_rh(half_width: f32, half_height: f32, near: f32, far: f32) -> Mat4 {
    Mat4::from_cols(
        Vec4::new(1.0 / half_width, 0.0, 0.0, 0.0),
        Vec4::new(0.0, 1.0 / half_height, 0.0, 0.0),
        Vec4::new(0.0, 0.0, 1.0 / (far - near), 0.0),
        Vec4::new(0.0, 0.0, far / (far - near), 1.0),
    )
}

/// Layout must match `Camera` in the WGSL shaders.
#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
pub struct CameraUniform {
    pub view_proj: [[f32; 4]; 4],
    pub view: [[f32; 4]; 4],
    pub proj: [[f32; 4]; 4],
    pub planes: [[f32; 4]; 6],
    pub viewport: [f32; 2],
    pub quad_px_threshold: f32,
    pub proj_scale: f32,
    /// Mip levels in the depth pyramid; the renderer fills this in.
    pub hiz_mip_count: u32,
    /// Non-zero enables occlusion culling against the previous frame.
    pub occlusion: u32,
    pub near: f32,
    /// `0` = perspective, `1` = orthographic -- see [`Projection`]. Every
    /// shader that redeclares `Camera` must branch on this wherever it
    /// casts a ray from a single eye point or sizes something by
    /// dividing by view depth; see [`Projection`]'s own doc for the full
    /// list of sites.
    pub projection: u32,
    /// View-space clip plane (`RenderSettings::clip`): points with
    /// `dot(xyz, p) + w < 0` are cut away. `[0, 0, 0, 1]` keeps all.
    pub clip: [f32; 4],
    /// [`Camera::near_cut`] in view space, in `clip`'s form.
    pub cut: [f32; 4],
}

/// A clip plane that cuts nothing.
pub const KEEP_ALL: [f32; 4] = [0.0, 0.0, 0.0, 1.0];

/// Six normalized planes (left, right, bottom, top, near, far) with inward
/// normals, so `dot(n, p) + d >= -r` means a sphere may be visible. The
/// near/far pair is read directly off the reversed-Z convention
/// (`z_clip == 0` at the far boundary, `z_clip == w_clip` at the near
/// one) both [`Camera::proj`] branches share, so this needs no changes
/// for [`Projection::Orthographic`]'s finite far plane even though it
/// was written against the infinite-far perspective matrix alone.
fn frustum_planes(m: Mat4) -> [[f32; 4]; 6] {
    let row = |i: usize| Vec4::new(m.x_axis[i], m.y_axis[i], m.z_axis[i], m.w_axis[i]);
    let (r0, r1, r2, r3) = (row(0), row(1), row(2), row(3));
    let planes = [r3 + r0, r3 - r0, r3 + r1, r3 - r1, r3 - r2, r2];
    planes.map(|p| {
        let len = p.truncate().length().max(1e-12);
        (p / len).to_array()
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rotating_the_scene_about_y_brings_its_right_side_to_the_front() {
        let mut camera = Camera::framing(Vec3::ZERO, 10.0);
        let right = Vec3::new(5.0, 0.0, 0.0);
        let before = (camera.view() * right.extend(1.0)).z;
        camera.rotate_scene(Vec3::Y, -std::f32::consts::FRAC_PI_2);
        let after = (camera.view() * right.extend(1.0)).z;
        assert!(after > before + 4.0, "{before} -> {after}");
    }

    #[test]
    fn uniform_layout_is_352_bytes() {
        assert_eq!(std::mem::size_of::<CameraUniform>(), 352);
    }

    #[test]
    fn framing_puts_sphere_inside_frustum() {
        let center = Vec3::new(10.0, -5.0, 3.0);
        let radius = 50.0;
        let cam = Camera::framing(center, radius);
        let u = cam.uniform(Vec2::new(1920.0, 1080.0));
        for plane in u.planes {
            let n = Vec3::new(plane[0], plane[1], plane[2]);
            assert!(
                n.dot(center) + plane[3] >= -radius,
                "sphere culled by plane {plane:?}"
            );
        }
        let behind = cam.eye() + (cam.eye() - center).normalize() * 10.0;
        let near = u.planes[4];
        let n = Vec3::new(near[0], near[1], near[2]);
        assert!(
            n.dot(behind) + near[3] < -1.0,
            "point behind camera should be culled"
        );
    }

    #[test]
    fn reversed_z_maps_near_to_one_and_far_to_zero() {
        let cam = Camera::framing(Vec3::ZERO, 10.0);
        let proj = cam.proj(1.0);
        let at = |z: f32| {
            let c = proj * Vec4::new(0.0, 0.0, -z, 1.0);
            c.z / c.w
        };
        assert!((at(cam.near) - 1.0).abs() < 1e-5);
        assert!(at(1.0e6) < 1e-4);
        assert!(at(1.0) > at(2.0));
    }

    /// Rotation never stops: tilting over the top carries the eye past the
    /// pole to the far side, and a full turn comes back to the start.
    #[test]
    fn rotation_goes_over_the_pole_and_all_the_way_round() {
        let mut cam = Camera::framing(Vec3::ZERO, 10.0);
        let start = cam.eye();
        let steps = 72;
        let step = std::f32::consts::TAU / steps as f32;
        for _ in 0..steps / 2 {
            cam.orbit(0.0, step);
        }
        assert!(
            cam.eye().z < -0.99 * start.z,
            "half a turn should reach the far side: {:?}",
            cam.eye()
        );
        for _ in 0..steps / 2 {
            cam.orbit(0.0, step);
        }
        assert!(
            cam.eye().distance(start) < 1e-2 * start.length(),
            "full turn: {:?}",
            cam.eye()
        );
    }

    /// With a pivot set, orbiting keeps the pivot at the
    /// same place on screen (the picture doesn't jump to it).
    #[test]
    fn orbiting_about_a_pivot_keeps_it_still_on_screen() {
        for projection in [Projection::Perspective, Projection::Orthographic] {
            let mut cam = Camera::framing(Vec3::ZERO, 20.0);
            cam.projection = projection;
            let pivot = Vec3::new(6.0, -4.0, 3.0);
            cam.pivot = Some(pivot);
            let screen = |c: &Camera| {
                let clip = c.proj(1.0) * c.view() * pivot.extend(1.0);
                clip.truncate().truncate() / clip.w
            };
            let before = screen(&cam);
            cam.orbit(0.4, 0.25);
            cam.orbit(-0.1, 0.3);
            let after = screen(&cam);
            assert!(
                before.distance(after) < 1e-4,
                "{projection:?}: pivot moved on screen {before:?} -> {after:?}"
            );
        }
    }

    #[test]
    fn orthographic_reversed_z_maps_near_to_one_and_far_to_zero() {
        let mut cam = Camera::framing(Vec3::ZERO, 10.0);
        cam.projection = Projection::Orthographic;
        let proj = cam.proj(1.0);
        let at = |z: f32| {
            let c = proj * Vec4::new(0.0, 0.0, -z, 1.0);
            // True orthographic: `c.w` is always 1, unlike perspective.
            assert!(
                (c.w - 1.0).abs() < 1e-5,
                "clip.w should stay 1, got {}",
                c.w
            );
            c.z
        };
        assert!((at(cam.near) - 1.0).abs() < 1e-5);
        // `far` is finite here (`Camera::proj`'s doc), so depth actually
        // reaches (rather than merely approaches) 0 at it.
        assert!(at(cam.orthographic_far()).abs() < 1e-4);
        assert!(
            at(1.0) > at(2.0),
            "depth should decrease linearly with distance"
        );
    }

    #[test]
    fn orthographic_does_not_scale_apparent_size_with_depth() {
        // The defining property this whole feature is for: a fixed
        // world-space offset projects to the same clip-space size no
        // matter how far along the view axis it's measured, unlike
        // perspective where it shrinks with distance.
        let mut cam = Camera::framing(Vec3::ZERO, 10.0);
        cam.projection = Projection::Orthographic;
        let proj = cam.proj(1.0);
        let width_at = |z: f32| {
            let c = proj * Vec4::new(1.0, 0.0, -z, 1.0);
            c.x / c.w
        };
        assert!((width_at(cam.near) - width_at(cam.distance)).abs() < 1e-5);
        assert!((width_at(cam.distance) - width_at(cam.distance * 5.0)).abs() < 1e-5);
    }

    #[test]
    fn orthographic_frames_the_same_sphere_perspective_does() {
        // Same `distance` (from `framing`), switched to orthographic
        // after the fact: the frustum should still contain the sphere
        // `framing` sized it for -- see `Camera::proj`'s doc on why
        // `fov_y`/`distance` still drive the orthographic half-extents.
        let center = Vec3::new(-3.0, 2.0, 1.0);
        let radius = 25.0;
        let mut cam = Camera::framing(center, radius);
        cam.projection = Projection::Orthographic;
        let u = cam.uniform(Vec2::new(1600.0, 900.0));
        for plane in u.planes {
            let n = Vec3::new(plane[0], plane[1], plane[2]);
            assert!(
                n.dot(center) + plane[3] >= -radius,
                "sphere culled by plane {plane:?}"
            );
        }
    }

    #[test]
    fn pan_moves_target_and_eye_by_the_same_offset() {
        let mut cam = Camera::framing(Vec3::new(1.0, 2.0, 3.0), 20.0);
        let (orientation, distance) = (cam.orientation, cam.distance);
        let eye_before = cam.eye();
        cam.pan(30.0, -15.0);
        // Panning translates the whole rig: orientation/distance (and so the
        // view direction) stay fixed, and the eye moves by the same
        // world-space offset as the target did.
        assert_eq!(cam.orientation, orientation);
        assert_eq!(cam.distance, distance);
        let target_offset = cam.target - Vec3::new(1.0, 2.0, 3.0);
        let eye_offset = cam.eye() - eye_before;
        assert!(
            (target_offset - eye_offset).length() < 1e-4,
            "target and eye should shift together: {target_offset:?} vs {eye_offset:?}"
        );
        assert!(
            target_offset.length() > 0.0,
            "pan should actually move something"
        );
    }

    #[test]
    fn zero_pan_does_not_move_the_target() {
        let mut cam = Camera::framing(Vec3::ZERO, 10.0);
        cam.pan(0.0, 0.0);
        assert_eq!(cam.target, Vec3::ZERO);
    }

    /// Orthographic: no cut until dollied, then one scene radius of dolly
    /// cuts through the target, facing the viewer, without changing the
    /// view. Perspective: the eye itself moves, even past the target.
    #[test]
    fn dolly_cuts_from_the_front_of_the_scene() {
        let mut cam = Camera::framing(Vec3::new(1.0, 2.0, 3.0), 10.0);
        cam.orbit(0.3, -0.2);
        assert_eq!(cam.near_cut_plane(), None);
        let (view, proj) = (cam.view(), cam.proj(1.0));
        cam.dolly_by(10.0);
        let [x, y, z, d] = cam.near_cut_plane().unwrap();
        let n = Vec3::new(x, y, z);
        assert!((n.dot(cam.target) + d).abs() < 1e-3, "through the target");
        assert!(n.dot(cam.forward()) > 0.999, "facing the viewer");
        assert_eq!((cam.view(), cam.proj(1.0)), (view, proj));
        cam.dolly_by(-20.0);
        assert_eq!(cam.dolly, 0.0, "never behind where it started");

        cam.projection = Projection::Perspective;
        let before = cam.eye();
        cam.dolly_by(cam.distance + 5.0);
        let moved = cam.eye() - before;
        assert!((moved - cam.forward() * (cam.distance + 5.0)).length() < 1e-3);
        let ahead = cam.view().transform_point3(cam.eye() + cam.forward());
        assert!(
            (ahead - Vec3::NEG_Z).length() < 1e-4,
            "still looking forward"
        );
    }
}
