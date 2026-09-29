// Full-screen copy and FXAA over one source texture. The vertex shader
// comes from fullscreen.wgsl, prepended at pipeline creation.

@group(0) @binding(0) var src: texture_2d<f32>;
@group(0) @binding(1) var src_sampler: sampler;

@fragment
fn fs_blit(in: Vs) -> @location(0) vec4<f32> {
    return textureSample(src, src_sampler, in.uv);
}

// FXAA (luma edge-detect + directional blend along the local gradient),
// following the structure of Timothy Lottes' public NVIDIA FXAA
// whitepaper -- a published technique, not code taken from any viewer
// (see the "clean-room" note in README.md: parsers are ours, cited
// published techniques are fine). Smooths the ray-cast impostor
// silhouettes that MSAA can't touch here (frag_depth disables it; see
// docs/RENDERING.md), at the cost of a full-screen pass every frame.

fn luma(c: vec3<f32>) -> f32 {
    return dot(c, vec3<f32>(0.299, 0.587, 0.114));
}

@fragment
fn fs_fxaa(in: Vs) -> @location(0) vec4<f32> {
    let texel = 1.0 / vec2<f32>(textureDimensions(src));
    let rgbM = textureSample(src, src_sampler, in.uv).rgb;
    let rgbN = textureSample(src, src_sampler, in.uv + vec2<f32>(0.0, -texel.y)).rgb;
    let rgbS = textureSample(src, src_sampler, in.uv + vec2<f32>(0.0, texel.y)).rgb;
    let rgbE = textureSample(src, src_sampler, in.uv + vec2<f32>(texel.x, 0.0)).rgb;
    let rgbW = textureSample(src, src_sampler, in.uv + vec2<f32>(-texel.x, 0.0)).rgb;
    let rgbNW = textureSample(src, src_sampler, in.uv + vec2<f32>(-texel.x, -texel.y)).rgb;
    let rgbNE = textureSample(src, src_sampler, in.uv + vec2<f32>(texel.x, -texel.y)).rgb;
    let rgbSW = textureSample(src, src_sampler, in.uv + vec2<f32>(-texel.x, texel.y)).rgb;
    let rgbSE = textureSample(src, src_sampler, in.uv + vec2<f32>(texel.x, texel.y)).rgb;

    let lM = luma(rgbM);
    let lN = luma(rgbN);
    let lS = luma(rgbS);
    let lE = luma(rgbE);
    let lW = luma(rgbW);
    let lNW = luma(rgbNW);
    let lNE = luma(rgbNE);
    let lSW = luma(rgbSW);
    let lSE = luma(rgbSE);

    let lMin = min(lM, min(min(lN, lS), min(lE, lW)));
    let lMax = max(lM, max(max(lN, lS), max(lE, lW)));
    let range = lMax - lMin;

    // Flat region: no visible edge here, skip the blur to keep detail.
    if (range < max(0.0312, lMax * 0.125)) {
        return vec4<f32>(rgbM, 1.0);
    }

    // Edge direction from the gradient of the four diagonal corners
    // (a cheap Sobel-like estimate), converted to a texel-space step.
    var dir = vec2<f32>(
        -((lNW + lNE) - (lSW + lSE)),
        (lNW + lSW) - (lNE + lSE),
    );
    let dir_reduce = max((lN + lS + lE + lW) * 0.03125, 0.0078125);
    let rcp_dir_min = 1.0 / (min(abs(dir.x), abs(dir.y)) + dir_reduce);
    dir = clamp(dir * rcp_dir_min, vec2<f32>(-8.0), vec2<f32>(8.0)) * texel;

    let rgbA = 0.5 * (
        textureSample(src, src_sampler, in.uv + dir * (1.0 / 3.0 - 0.5)).rgb +
        textureSample(src, src_sampler, in.uv + dir * (2.0 / 3.0 - 0.5)).rgb
    );
    let rgbB = rgbA * 0.5 + 0.25 * (
        textureSample(src, src_sampler, in.uv + dir * -0.5).rgb +
        textureSample(src, src_sampler, in.uv + dir * 0.5).rgb
    );
    let lB = luma(rgbB);
    if (lB < lMin || lB > lMax) {
        return vec4<f32>(rgbA, 1.0);
    }
    return vec4<f32>(rgbB, 1.0);
}
