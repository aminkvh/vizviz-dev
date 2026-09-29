// Composite pass over the opaque pass's color, normal and depth targets,
// plus ao.wgsl's occlusion. Runs every frame (a plain copy when every
// effect is off) so the texture chain never changes shape -- see the
// comment on `Renderer::display_view`. The vertex shader comes from
// fullscreen.wgsl. In order:
//
// 1. Ambient occlusion and shadows: a depth-aware 4x4 average of the
//    half-resolution AO target (r = occlusion, g = key-light visibility;
//    upsamples it and cancels its 2x2 rotation pattern without bleeding
//    across edges). A shadow removes up to 60% of the pixel's light.
// 2. Depth cue: fade toward the background with view depth across the
//    scene's own depth range (`fog_near`..`fog_far`), on a linear ramp.
//
// Background pixels get the vertical gradient from `background` (bottom)
// to `background_top`; equal colours make it solid.
// 3. Outline: darkens pixels on a depth or normal discontinuity.
//    - depth: the neighbour is more than `depth_threshold` (grown with
//      view depth, so distant structures don't dissolve into noise)
//      nearer or farther. Background counts as infinitely far, which is
//      what draws the silhouette on the object's own rim pixels.
//    - normal: the neighbour faces a clearly different way, which
//      catches the crease where two overlapping atoms meet even though
//      the depth there is continuous.

@group(0) @binding(0) var color: texture_2d<f32>;
@group(0) @binding(1) var normal: texture_2d<f32>;
@group(0) @binding(2) var depth: texture_depth_2d;

// Layout must match `vv_render::renderer::PostUniform` (and `Post` in
// ao.wgsl).
struct Post {
    inv_proj: mat4x4<f32>,
    outline_enabled: u32,
    outline_width: u32,
    near: f32,
    depth_threshold: f32,
    normal_threshold: f32,
    outline_strength: f32,
    ao_strength: f32,
    ao_radius: f32,
    fog_strength: f32,
    fog_near: f32,
    fog_far: f32,
    proj_y: f32,
    background: vec4<f32>,
    orthographic: u32,
    shadow_strength: f32,
    _pad0: u32,
    _pad1: u32,
    background_top: vec4<f32>,
    // View-space direction toward the key light (xyz).
    light_dir: vec4<f32>,
    view_inv: mat4x4<f32>,
    // Occlusion volume: corner (xyz) and voxel size (w); size in world
    // units (xyz) and 1 when there is a volume (w).
    vol_origin: vec4<f32>,
    vol_extent: vec4<f32>,
    // World-space clip plane; the volume is empty on its cut side.
    clip_world: vec4<f32>,
    // The camera's near cut, as `clip_world`.
    cut_world: vec4<f32>,
};
@group(0) @binding(3) var<uniform> post: Post;
@group(0) @binding(4) var ao_target: texture_2d<f32>;

// Positive view depth of a depth-buffer value (0 = background). Both
// projections are centred, so only the z and w rows of the inverse
// projection matter.
fn view_depth(d: f32) -> f32 {
    if (d <= 0.0) {
        return 1e30;
    }
    let z = post.inv_proj[2][2] * d + post.inv_proj[3][2];
    let w = post.inv_proj[2][3] * d + post.inv_proj[3][3];
    return -z / w;
}

fn blurred_ao(p: vec2<i32>, zc: f32, dims: vec2<i32>) -> vec2<f32> {
    var sum = vec2<f32>(0.0);
    var weight = 0.0;
    let tolerance = max(0.5, zc * 0.02);
    for (var y = -2; y < 2; y++) {
        for (var x = -2; x < 2; x++) {
            let q = clamp(p + vec2<i32>(x, y), vec2<i32>(0), dims - 1);
            let zq = view_depth(textureLoad(depth, q, 0));
            if (abs(zq - zc) < tolerance) {
                let ao_dims = vec2<i32>(textureDimensions(ao_target));
                sum += textureLoad(ao_target, min(q / 2, ao_dims - 1), 0).rg;
                weight += 1.0;
            }
        }
    }
    return select(vec2<f32>(1.0), sum / weight, weight > 0.0);
}

@fragment
fn fs_outline(in: Vs) -> @location(0) vec4<f32> {
    let dims = vec2<i32>(textureDimensions(color));
    let p = vec2<i32>(in.clip.xy);
    let c = textureLoad(color, p, 0);
    let dc = textureLoad(depth, p, 0);
    let up = 1.0 - (f32(p.y) + 0.5) / f32(dims.y);
    let bg = mix(post.background, post.background_top, up);
    if (dc <= 0.0) {
        return bg;
    }
    let zc = view_depth(dc);
    var rgb = c.rgb;

    if (post.ao_strength > 0.0 || post.shadow_strength > 0.0) {
        let ao = blurred_ao(p, zc, dims);
        rgb *= ao.x * (1.0 - 0.6 * post.shadow_strength * (1.0 - ao.y));
    }
    if (post.fog_strength > 0.0 && post.fog_far > post.fog_near) {
        let t = clamp((zc - post.fog_near) / (post.fog_far - post.fog_near), 0.0, 1.0);
        rgb = mix(rgb, bg.rgb, post.fog_strength * t);
    }

    if (post.outline_enabled == 0u) {
        return vec4<f32>(rgb, c.a);
    }
    let nc = normalize(textureLoad(normal, p, 0).xyz * 2.0 - 1.0);
    let threshold = max(post.depth_threshold, zc * 0.01);
    let w = i32(post.outline_width);
    let offsets = array<vec2<i32>, 4>(
        vec2<i32>(w, 0), vec2<i32>(-w, 0), vec2<i32>(0, w), vec2<i32>(0, -w),
    );
    var edge = false;
    for (var k = 0; k < 4; k++) {
        let q = clamp(p + offsets[k], vec2<i32>(0), dims - 1);
        let dq = textureLoad(depth, q, 0);
        let zq = view_depth(dq);
        if (abs(zq - zc) > threshold) {
            edge = true;
            break;
        }
        if (dq > 0.0) {
            let nq = normalize(textureLoad(normal, q, 0).xyz * 2.0 - 1.0);
            if (dot(nc, nq) < post.normal_threshold) {
                edge = true;
                break;
            }
        }
    }
    if (!edge) {
        return vec4<f32>(rgb, c.a);
    }
    return vec4<f32>(rgb * (1.0 - post.outline_strength), c.a);
}
