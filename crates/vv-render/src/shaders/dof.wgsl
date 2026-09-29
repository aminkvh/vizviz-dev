// Depth of field over the composited image (a copy of `outlined`, which
// this pass then overwrites, so screenshots include it). The vertex
// shader comes from fullscreen.wgsl.
//
// Each pixel's blur radius (circle of confusion) grows with its view
// depth's distance from the focus depth. A gather over a golden-angle
// spiral: a sample counts when its own blur reaches this pixel, and a
// sample behind this pixel counts only as far as this pixel's own blur
// reaches, so a blurred background never bleeds over a sharp subject.

struct Dof {
    // Clip -> view rows needed for view depth (see outline.wgsl).
    inv_proj: mat4x4<f32>,
    focus: f32,
    // Largest blur radius, in pixels.
    max_radius: f32,
    // Depth distance from the focus at which the blur is largest: half
    // the scene's depth span.
    range: f32,
    _pad0: f32,
};

@group(0) @binding(0) var color: texture_2d<f32>;
@group(0) @binding(1) var depth: texture_depth_2d;
@group(0) @binding(2) var<uniform> dof: Dof;

const SAMPLES: u32 = 32u;

fn view_depth(d: f32) -> f32 {
    if (d <= 0.0) {
        return 1e30;
    }
    let z = dof.inv_proj[2][2] * d + dof.inv_proj[3][2];
    let w = dof.inv_proj[2][3] * d + dof.inv_proj[3][3];
    return -z / w;
}

fn blur_radius(z: f32) -> f32 {
    return dof.max_radius * clamp(abs(z - dof.focus) / dof.range, 0.0, 1.0);
}

@fragment
fn fs_dof(in: Vs) -> @location(0) vec4<f32> {
    let dims = vec2<i32>(textureDimensions(color));
    let p = vec2<i32>(in.clip.xy);
    let zc = view_depth(textureLoad(depth, p, 0));
    let rc = blur_radius(zc);
    var sum = textureLoad(color, p, 0);
    var weight = 1.0;
    for (var i = 0u; i < SAMPLES; i++) {
        let f = (f32(i) + 0.5) / f32(SAMPLES);
        let dist = sqrt(f) * dof.max_radius;
        let angle = f32(i) * 2.39996323;
        let q = clamp(p + vec2<i32>(round(vec2<f32>(cos(angle), sin(angle)) * dist)), vec2<i32>(0), dims - 1);
        let zq = view_depth(textureLoad(depth, q, 0));
        var reach = blur_radius(zq);
        if (zq > zc) {
            reach = min(reach, rc);
        }
        let w = clamp(reach - dist + 1.0, 0.0, 1.0);
        sum += textureLoad(color, q, 0) * w;
        weight += w;
    }
    return sum / weight;
}
