// Billboards for surfaces drawn one ray-cast patch at a time (skin, SES):
// six vertices per patch, no vertex buffer, each quad sized to contain the
// patch's bounding sphere; the fragment shader re-derives its ray from the
// interpolated position and solves the real geometry.

struct VsOut {
    @builtin(position) clip: vec4<f32>,
    @location(0) view_pos: vec3<f32>,
    @location(1) @interpolate(flat) patch_index: u32,
};

// Two triangles: (-1,-1) (1,-1) (1,1) / (-1,-1) (1,1) (-1,1).
fn corner(vi: u32) -> vec2<f32> {
    let cx = select(-1.0, 1.0, vi == 1u || vi == 2u || vi == 4u);
    let cy = select(-1.0, 1.0, vi == 2u || vi == 4u || vi == 5u);
    return vec2<f32>(cx, cy);
}

// Same enlargement `draw.wgsl`'s `vs_sphere` uses: a billboard on the
// near tangent plane always contains the perspective silhouette of a
// sphere of radius r at distance d.
// `draw.wgsl`'s `vs_sphere` places its billboard at the sphere's own
// near tangent point (`center + toward_eye * r`, distance `d - r` from
// the camera) sized by `silhouette_enlarge(d, r) = d / sqrt(d^2 - r^2)`
// -- what matters for a fragment on that quad is not the quad's own
// distance but its *direction* from the camera (the fragment shader
// re-derives its own ray from the interpolated position and solves the
// real geometry), so this only needs the corner directions to reach the
// tangent cone's half-angle `asin(r / d)`; using `d` rather than the
// plane's own distance `d - r` in the formula overshoots that on
// purpose, a comfortable margin against float32 rounding right at the
// silhouette (the tight, non-oversized formula `(d-r)*r/sqrt(d^2-r^2)`
// was tried here first and measurably lost silhouette pixels to
// rounding -- caught by `vv-cpu/tests/agreement.rs`, not assumed).
//
// That placement is only valid while the camera is safely outside the
// sphere -- an atom's own VDW radius, what this technique was designed
// for, makes "camera inside the atom" a non-issue in practice. A skin-
// surface patch's *bounding* sphere is routinely much larger (it covers
// a whole mixed cell, not one atom), so at ordinary close-up zoom the
// camera can end up inside or very near one, and the unclamped formula
// then places the plane BEHIND the camera (`toward_eye * r` overshoots
// past the origin once `r > d`) -- hardware clipping silently drops
// real coverage there, a second real bug the same test caught. Fixed by
// clamping the plane's distance to `l = max(d - r, margin)` while
// keeping the *original* generous half-width formula (evaluated at the
// true `d`, not the clamped `l`): since half-width at fixed direction-
// coverage only grows as the plane moves closer to the camera, using
// the same numerator at a closer-than-`d-r` plane is still safely
// sufficient, never tight.
// View-space position of billboard corner `vertex % 6` for a bounding
// sphere at view-space `center`, radius `r`.
fn billboard(center: vec3<f32>, r: f32, vertex: u32) -> vec3<f32> {
    let c = corner(vertex % 6u);
    var pos: vec3<f32>;
    if (cam.projection == 1u) {
        // Orthographic: silhouette size is distance-independent; a
        // square of half-width r exactly contains a circle of radius r,
        // wherever the plane sits along the fixed view direction.
        pos = center + vec3<f32>(0.0, 0.0, r) + vec3<f32>(c.x, c.y, 0.0) * r;
    } else if (length(center) < r + max(cam.near * 2.0, r * 0.02)) {
        // The eye is inside (or touching) the sphere, so any direction may
        // hit it: cover the whole view, just past the near plane.
        let l = cam.near * 2.0;
        pos = vec3<f32>(c.x * l / cam.proj[0][0], c.y * l / cam.proj[1][1], -l);
    } else {
        let d = length(center);
        let margin = max(cam.near * 2.0, r * 0.02);
        let l = max(d - r, margin);
        let half_width = r * d / sqrt(max(d * d - r * r, margin * margin));
        let toward_eye = select(vec3<f32>(0.0, 0.0, 1.0), -center / d, d > 1e-6);
        // Object-facing basis, NOT `vec3(c.x, c.y, 0.0)`'s screen-aligned
        // one (what `vs_sphere` uses): for a patch far off the camera's
        // central axis -- routine for a bound sphere, which is larger
        // and less centered than a single atom typically is -- a
        // screen-aligned quad is foreshortened relative to the true
        // sightline, and the isotropic half-width above (correct for a
        // quad that actually faces the patch) under-covers in the
        // away-from-center direction. This real, measured gap and this
        // fix both came from `vv-cpu/tests/agreement.rs`'s CPU/GPU
        // comparison at a wide, off-center camera angle, not a
        // hypothetical -- an isolated close-up render centered
        // directly on the same patch did NOT reproduce it, which is
        // what pointed at off-axis foreshortening specifically.
        let up_hint = select(vec3<f32>(0.0, 1.0, 0.0), vec3<f32>(1.0, 0.0, 0.0), abs(toward_eye.y) > 0.99);
        let right = normalize(cross(up_hint, toward_eye));
        let up = cross(toward_eye, right);
        pos = center + toward_eye * (d - l) + (right * c.x + up * c.y) * half_width;
    }

    return pos;
}
