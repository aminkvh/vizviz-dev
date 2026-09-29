// Composites the glass pass (weighted blended OIT, see `glass()` in
// shading.wgsl) over the opaque image: the weighted average colour,
// covering 1 - reveal of the pixel (blend: src alpha, one minus src
// alpha). The vertex shader comes from fullscreen.wgsl.

@group(0) @binding(0) var accum: texture_2d<f32>;
@group(0) @binding(1) var reveal: texture_2d<f32>;

@fragment
fn fs_glass_resolve(in: Vs) -> @location(0) vec4<f32> {
    let p = vec2<i32>(in.clip.xy);
    let r = textureLoad(reveal, p, 0).r;
    if (r >= 1.0) {
        discard;
    }
    let a = textureLoad(accum, p, 0);
    let color = a.rgb / max(a.a, 1e-5);
    return vec4<f32>(color, 1.0 - r);
}
