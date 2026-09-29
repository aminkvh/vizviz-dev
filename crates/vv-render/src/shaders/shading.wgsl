// Shared lighting and material shading, prepended to every shader that
// shades a surface. Ambient tints the object's colour; the outline darkens
// the whole colour toward the silhouette. GPU twin of
// `vv_render::style::shade` and `style::alpha`.

// Layout must match `vv_render::style::Light`.
struct Light {
    dir: vec3<f32>,
    intensity: f32,
    color: vec3<f32>,
    _pad: f32,
};

// Layout must match `vv_render::style::Lighting`.
struct Lighting {
    lights: array<Light, 4>,
    ambient_color: vec3<f32>,
    ambient: f32,
    tonemap: f32,
    // How many of `lights` (from index 0) actually light the scene.
    count: f32,
    _pad0: f32,
    _pad1: f32,
};
@group(0) @binding(2) var<uniform> lighting: Lighting;

// `vv_render::style::Material`, carried in each shader's params as three
// vec4s: (ambient, diffuse, specular, shininess), (toon_bands, opacity,
// outline, outline_width), (transmode, rim, -, -).
struct Material {
    ambient: f32,
    diffuse: f32,
    specular: f32,
    shininess: f32,
    toon_bands: f32,
    opacity: f32,
    outline: f32,
    outline_width: f32,
    transmode: f32,
    rim: f32,
};

fn material(m0: vec4<f32>, m1: vec4<f32>, m2: vec4<f32>) -> Material {
    return Material(m0.x, m0.y, m0.z, m0.w, m1.x, m1.y, m1.z, m1.w, m2.x, m2.y);
}

fn tonemap(c: vec3<f32>) -> vec3<f32> {
    if (lighting.tonemap < 0.5) {
        return c;
    }
    let mapped = (c * (2.51 * c + 0.03)) / (c * (2.43 * c + 0.59) + 0.14);
    return clamp(mapped, vec3<f32>(0.0), vec3<f32>(1.0));
}

fn blinn(n: vec3<f32>, light: vec3<f32>, view_dir: vec3<f32>, shininess: f32) -> f32 {
    return pow(max(dot(n, normalize(light - view_dir)), 0.0), shininess);
}

// Palette colours are sRGB display values; lighting works in linear light, and the sRGB target encodes
// the result back. Shaders decode where they unpack a colour (per
// vertex where they can: per fragment costs ~3 ms on an iGPU), so every
// `base` below is linear.
fn srgb_to_linear(c: vec3<f32>) -> vec3<f32> {
    return select(pow((c + 0.055) / 1.055, vec3<f32>(2.4)), c / 12.92, c <= vec3<f32>(0.04045));
}

fn linear_rgba(c: vec4<f32>) -> vec4<f32> {
    return vec4<f32>(srgb_to_linear(c.rgb), c.a);
}

fn plane_distance(plane: vec4<f32>, p: vec3<f32>) -> f32 {
    return dot(plane.xyz, p) + plane.w;
}

// Signed distance of view-space `p` from the planes that cut the scene
// (the clip plane and the camera's near cut); negative is cut away.
// Always positive with both off.
fn clip_distance(p: vec3<f32>) -> f32 {
    return min(plane_distance(cam.clip, p), plane_distance(cam.cut, p));
}

// The stretch [start, end] of a view-space ray, within [t0, t1], that no
// plane cuts away (empty when start > end), and the outward normal of the
// cut face where it starts: what a solid shows where a plane opened it
// (when start > t0).
struct Kept {
    start: f32,
    end: f32,
    normal: vec3<f32>,
};

fn kept_by(k: Kept, plane: vec4<f32>, ro: vec3<f32>, rd: vec3<f32>) -> Kept {
    var out = k;
    let rate = dot(plane.xyz, rd);
    let d = plane_distance(plane, ro);
    if (abs(rate) < 1e-12) {
        if (d < 0.0) {
            out.end = out.start - 1.0;
        }
        return out;
    }
    let t = -d / rate;
    if (rate > 0.0 && t > out.start) {
        out.start = t;
        out.normal = -plane.xyz;
    } else if (rate < 0.0) {
        out.end = min(out.end, t);
    }
    return out;
}

fn kept_span(ro: vec3<f32>, rd: vec3<f32>, t0: f32, t1: f32) -> Kept {
    let k = Kept(t0, t1, vec3<f32>(0.0, 0.0, 1.0));
    return kept_by(kept_by(k, cam.clip, ro, rd), cam.cut, ro, rd);
}

// Where the cap of a cut-open solid lies on the camera ray through view-
// space `hit`, the ray's nearest surviving hit and a back face; `start` < 0
// when no plane cut the ray before it (a back face seen some other way).
fn cap_on_ray(hit: vec3<f32>, orthographic: bool) -> Kept {
    var ro = vec3<f32>(0.0);
    var rd = normalize(hit);
    if (orthographic) {
        ro = vec3<f32>(hit.xy, 0.0);
        rd = vec3<f32>(0.0, 0.0, -1.0);
    }
    var k = kept_span(ro, rd, 0.0, dot(hit - ro, rd));
    if (k.start <= 0.0 || k.start > k.end) {
        k.start = -1.0;
    }
    return k;
}

// Caps of overlapping primitives share one plane: each is pulled toward
// the viewer by its power at the point (`inside`: r^2 less the squared
// distance to its centre or axis), so the one it is deepest in wins and
// neighbours meet on their radical line, not in draw order: well above
// depth resolution, yet at most halfway back to the billboard (at `t_quad`
// on the ray), which the early depth test needs the cap to stay behind.
fn cap_nudge(inside: f32, t: f32, t_quad: f32) -> f32 {
    return min(1e-2 * inside, 0.5 * (t - t_quad));
}

// A view-space sphere (bounding a patch or an atom) entirely on the cut
// side of an active plane -- `(0, 0, 0, 1)` is `KEEP_ALL`, always
// distance 1 from anywhere, which a real plane's unit normal never is, so
// an inactive plane never reports this. A patch surface's own vertex
// shader (`vs_ses_surface`/`vs_skin_surface`) calls this on the patch's
// *bounding* sphere to degenerate one wholly past the cut before it
// rasterizes: the removed side is always the nearest thing on its own
// rays, so nothing else in the depth buffer can early-Z reject it, and
// without this every one of them still pays a full ray-cast down to its
// own `discard` -- measured turning a 3x slowdown on a 4779-atom SES
// (patches on the cut side dominating draw time) into clipped being as
// fast as, or faster than, whole (RTX 3060 Laptop, every zoom preset).
fn wholly_cut(center: vec3<f32>, r: f32) -> bool {
    return plane_distance(cam.clip, center) <= -r || plane_distance(cam.cut, center) <= -r;
}

// A cheap, direction-independent over-approximation of "no ray through
// this view-space sphere can find a cap": true unless an active plane
// crosses it, or `wholly_cut` removes it outright. An atom-disc cap layer
// (`ses_surface.wgsl`/`skin_surface.wgsl`'s `vs_atom_cap`) calls this per
// atom to degenerate the ones a clip plane cannot possibly touch, so its
// fragment shader only ever runs near the cut.
fn no_cap_possible(center: vec3<f32>, r: f32) -> bool {
    let touches = (dot(cam.clip.xyz, cam.clip.xyz) > 0.5 && abs(plane_distance(cam.clip, center)) < r)
        || (dot(cam.cut.xyz, cam.cut.xyz) > 0.5 && abs(plane_distance(cam.cut, center)) < r);
    return wholly_cut(center, r) || !touches;
}

// Diffuse-lit colour and specular, before the tonemap.
struct Lit {
    diffuse: vec3<f32>,
    specular: f32,
};

// `n` and `view_dir` (eye into the scene) in view space.
fn light_surface(base: vec3<f32>, n: vec3<f32>, view_dir: vec3<f32>, m: Material) -> Lit {
    return light_surface_seen(base, n, view_dir, m, array<f32, 4>(1.0, 1.0, 1.0, 1.0), 1.0);
}

// As `light_surface`, with how much of each light (`seen_lights[i]`) and
// the sky/ambient (`seen_sky`) reaches the point: all of it in the
// viewport; what shadow and occlusion rays find in the path tracer
// (path_trace.wgsl).
fn light_surface_seen(base: vec3<f32>, n: vec3<f32>, view_dir: vec3<f32>, m: Material, seen_lights: array<f32, 4>, seen_sky: f32) -> Lit {
    var diffuse_light = lighting.ambient_color * (lighting.ambient * seen_sky);
    let count = u32(lighting.count);
    for (var i = 0u; i < count; i++) {
        let l = lighting.lights[i];
        var ndotl = max(dot(n, l.dir), 0.0);
        if (m.toon_bands > 0.5) {
            ndotl = floor(ndotl * m.toon_bands) / m.toon_bands;
        }
        diffuse_light += l.color * (l.intensity * ndotl * seen_lights[i]);
    }
    var light = vec3<f32>(m.ambient) + m.diffuse * diffuse_light;
    if (m.outline > 0.0) {
        let facing = -dot(n, view_dir);
        let edge = 1.0 - pow(max(1.0 - facing * facing, 0.0), (1.0 - m.outline_width) * 32.0);
        light *= max(mix(1.0, edge, m.outline), 0.0);
    }
    if (m.rim > 0.0) {
        let facing = clamp(-dot(n, view_dir), 0.0, 1.0);
        light += vec3<f32>(m.rim * pow(1.0 - facing, 3.0));
    }
    // First light only, achromatic: the fill light sits behind the
    // molecule (its highlight lands on rims); extending this to every light and colour costs a `pow` per light
    // per fragment for a highlight that reads the same either way.
    let key = lighting.lights[0];
    let spec = m.specular * key.intensity * seen_lights[0] * blinn(n, key.dir, view_dir, m.shininess);
    return Lit(base * light, spec);
}

fn shade(base: vec3<f32>, n: vec3<f32>, view_dir: vec3<f32>, m: Material) -> vec3<f32> {
    let lit = light_surface(base, n, view_dir, m);
    return tonemap(lit.diffuse + vec3<f32>(lit.specular));
}

// Weighted blended order-independent transparency (McGuire & Bavoil,
// JCGT 2013): every glass fragment adds its premultiplied colour and
// alpha, weighted to favour the nearest layers, to `accum`, and
// multiplies `reveal` (the light that gets through) by 1 - alpha;
// shaders/glass_resolve.wgsl composites the two over the opaque image.
// Used by every representation's own `_glass` fragment entry point
// (spheres, cylinders, cartoons, glycans, Gaussian/skin/SES surfaces),
// so they all deposit into the same pair of targets and blend as one
// unsorted set regardless of draw order or which shader wrote them.
struct GlassColor {
    @location(0) accum: vec4<f32>,
    @location(1) reveal: vec4<f32>,
};

// As `GlassColor`, plus the fragment's own ray-cast depth: for an
// impostor (sphere, cylinder, or a patch/voxel surface), the rasterized
// billboard depth is only a conservative bound, not the true surface hit
// — the glass pass's depth test (against the opaque depth, write
// disabled) needs the real one, exactly as the opaque pass's `frag_depth`
// does. A rasterized triangle mesh (cartoon, glycan) has no such gap and
// returns `GlassColor` directly instead.
struct GlassOut {
    @location(0) accum: vec4<f32>,
    @location(1) reveal: vec4<f32>,
    @builtin(frag_depth) depth: f32,
};

// Depth-favouring weight from McGuire & Bavoil's reference weight
// function: strongly favours nearby layers within camera range while
// never fully zeroing out a distant one.
fn glass_weight(view_depth: f32) -> f32 {
    return clamp(10.0 / (1e-5 + pow(view_depth / 5.0, 2.0) + pow(view_depth / 200.0, 6.0)), 1e-2, 3e3);
}

// Shared by every `_glass` entry point once it has its own premultiplied
// colour and alpha: `glass_resolve.wgsl` reconstructs the average visible
// colour as `accum.rgb / accum.a`, which recovers the lit colour (bounded
// to at most ~1, like an opaque surface) only when the whole numerator —
// specular included — was scaled by `alpha` here, not added at full
// strength. Unscaled specular would push that reconstructed colour above
// 1.0 for any material whose specular exceeds its own alpha (Faint keeps
// specular 1.0 at opacity 0.1), and a fixed-point render target clamps
// colour to 1.0 *before* blending, so the highlight would blow out to
// flat white however transparent the surface is.
fn glass_accum(premultiplied: vec3<f32>, alpha: f32, view_depth: f32) -> GlassColor {
    var out: GlassColor;
    out.accum = vec4<f32>(premultiplied, alpha) * glass_weight(view_depth);
    out.reveal = vec4<f32>(alpha);
    return out;
}

// Only front faces reach here (callers discard the rest), so a closed
// surface is one layer. Alpha is `material_alpha` exactly (no Fresnel
// boost): the flat-opacity glass materials have the same opacity at every
// angle, and this must equal what `path_trace.wgsl` tests its own
// transmission probability against, or a rasterized glass sphere would
// read glossier than the same sphere path-traced.
fn glass_lit(base: vec3<f32>, n: vec3<f32>, view_dir: vec3<f32>, m: Material, view_depth: f32) -> GlassColor {
    let lit = light_surface(base, n, view_dir, m);
    let alpha = material_alpha(n, view_dir, m);
    let premultiplied = tonemap(lit.diffuse + vec3<f32>(lit.specular)) * alpha;
    return glass_accum(premultiplied, alpha, view_depth);
}

// As `glass_lit`, for a ray-cast impostor whose fragment must also
// override its own depth (see `GlassOut`'s doc).
fn glass(base: vec3<f32>, n: vec3<f32>, view_dir: vec3<f32>, m: Material, view_depth: f32, depth: f32) -> GlassOut {
    let c = glass_lit(base, n, view_dir, m, view_depth);
    var out: GlassOut;
    out.accum = c.accum;
    out.reveal = c.reveal;
    out.depth = depth;
    return out;
}

// As `glass`, for a sub-pixel point or line with no usable surface
// normal: no specular, no `transmode` angle term (both need one), just
// `shade_point`'s ambient+diffuse average scaled by opacity. Rasterized
// depth (no ray-cast gap), so it returns `GlassColor`, not `GlassOut`.
fn glass_point(base: vec3<f32>, m: Material, view_depth: f32) -> GlassColor {
    let alpha = m.opacity;
    return glass_accum(shade_point(base, m) * alpha, alpha, view_depth);
}

// Angle-dependent opacity (Merritt & Bacon 1997, Methods Enzymol. 277:505)
// when `transmode` is set: clearer face-on, more opaque at grazing angles.
fn material_alpha(n: vec3<f32>, view_dir: vec3<f32>, m: Material) -> f32 {
    if (m.transmode < 0.5) {
        return m.opacity;
    }
    let a = 1.0 + cos(3.1415926 * (1.0 - m.opacity) * -dot(n, view_dir));
    return a * a * 0.25;
}

// A sub-pixel dot has no normal: roughly what a lit sphere averages, each
// light contributing a flat share (0.4, chosen so two lights -- the old
// key+fill default -- sum to the old formula's 0.6+0.2).
fn shade_point(base: vec3<f32>, m: Material) -> vec3<f32> {
    var diffuse_light = lighting.ambient_color * lighting.ambient;
    let count = u32(lighting.count);
    for (var i = 0u; i < count; i++) {
        diffuse_light += lighting.lights[i].color * (lighting.lights[i].intensity * 0.4);
    }
    let light = vec3<f32>(m.ambient) + m.diffuse * diffuse_light;
    return tonemap(base * light);
}
