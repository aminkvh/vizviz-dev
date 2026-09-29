// Moves each cap's depth onto the plane that cut it. Mesh and patch
// surfaces draw a cap as the back face the cut exposes, at that face's
// own depth, so the nearest surviving hit still wins the depth test;
// they mark it with normal.a = 0. AO, the outline, the depth cue and
// next frame's occlusion culling need it where it shows: on the plane.
// fullscreen.wgsl and shading.wgsl are prepended.

// Must match `atoms.wgsl`'s `Camera`.
struct Camera {
    view_proj: mat4x4<f32>,
    view: mat4x4<f32>,
    proj: mat4x4<f32>,
    planes: array<vec4<f32>, 6>,
    viewport: vec2<f32>,
    quad_px_threshold: f32,
    proj_scale: f32,
    hiz_mip_count: u32,
    occlusion: u32,
    near: f32,
    projection: u32,
    clip: vec4<f32>,
    cut: vec4<f32>,
};
@group(0) @binding(0) var<uniform> cam: Camera;

@group(1) @binding(0) var normals: texture_2d<f32>;

@fragment
fn fs_cap_depth(in: Vs) -> @builtin(frag_depth) f32 {
    if (textureLoad(normals, vec2<i32>(in.clip.xy), 0).a > 0.5) {
        discard;
    }
    let ndc = vec2<f32>(in.uv.x * 2.0 - 1.0, 1.0 - in.uv.y * 2.0);
    let lateral = vec2<f32>(ndc.x / cam.proj[0][0], ndc.y / cam.proj[1][1]);
    var ro = vec3<f32>(0.0);
    var rd = normalize(vec3<f32>(lateral, -1.0));
    if (cam.projection == 1u) {
        ro = vec3<f32>(lateral, 0.0);
        rd = vec3<f32>(0.0, 0.0, -1.0);
    }
    let k = kept_span(ro, rd, 0.0, 3.4e38);
    if (k.start <= 0.0) {
        discard;
    }
    let clip = cam.proj * vec4<f32>(ro + rd * k.start, 1.0);
    return clip.z / clip.w;
}
