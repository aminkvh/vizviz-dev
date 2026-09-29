// GPU ray-caster for the solvent-excluded surface: one billboard per
// patch (billboard.wgsl), each fragment solving only its own patch
// (ses_patch.wgsl); the depth buffer resolves the rest.

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
    // `vv_render::camera::CameraUniform::cut`, in `clip`'s form.
    cut: vec4<f32>,
};
@group(0) @binding(0) var<uniform> cam: Camera;

// xyz = position, w = van der Waals radius.
@group(1) @binding(0) var<storage, read> atoms: array<vec4<f32>>;
@group(1) @binding(1) var<storage, read> colors: array<u32>;

// Layout must match `vv_render::ses_surface::SesParams`.
struct SesParams {
    view_inv: mat4x4<f32>,
    probe_radius: f32,
    atom_count: u32,
    patch_count: u32,
    id_base: u32,
    material0: vec4<f32>,
    material1: vec4<f32>,
    material2: vec4<f32>,
};
@group(1) @binding(2) var<uniform> params: SesParams;

@group(1) @binding(4) var<storage, read> caps: array<u32>;
@group(1) @binding(5) var<storage, read> probes: array<Probe>;
@group(1) @binding(6) var<storage, read> probe_neighbors: array<u32>;

fn unpack_color(c: u32) -> vec3<f32> {
    return vec3<f32>(f32(c & 0xffu), f32((c >> 8u) & 0xffu), f32((c >> 16u) & 0xffu)) / 255.0;
}

fn probe_radius() -> f32 {
    return params.probe_radius;
}

fn atom_color(a: u32) -> vec3<f32> {
    return srgb_to_linear(unpack_color(colors[a]));
}

// One page of the patches' bounding spheres (`SesGpu::bounds`), the
// cull pass's survivors in it, and its patches, all page-relative.
struct PageView {
    base: u32,
};
@group(3) @binding(0) var<storage, read> bounds: array<vec4<f32>>;
@group(3) @binding(1) var<storage, read> visible: array<u32>;
@group(3) @binding(2) var<uniform> page: PageView;
// Top two bits: kind. Convex: (atom, caps start, caps end, -); torus:
// (atom i, atom j, probe a, probe b) with FULL_CIRCLE for a whole circle;
// concave: (probe, -, -, -).
@group(3) @binding(3) var<storage, read> patches: array<vec4<u32>>;

// Below this probe radius on screen concave patches are skipped, and
// convex patches drawn whole: every crevice is then under a pixel or so.
const PROBE_PX: f32 = 2.0;
// Tori whose pair's `fillet_size` is under this on screen are skipped too,
// and convex patches ignore those pairs' caps, so the atom's own sphere
// shows there instead of a hole.
const TORUS_PX: f32 = 3.0;
// Skipping a fillet happens a little below the size at which its atom
// stops cutting the cap, so a fillet a bit farther away than its atom is
// never skipped while the cap is still cut.
const SKIP_MARGIN: f32 = 0.9;

// World units per pixel at view-space `v` (`cull.wgsl`'s screen_px, inverted).
fn world_per_px(v: vec3<f32>) -> f32 {
    if (cam.projection == 1u) {
        return 1.0 / cam.proj_scale;
    }
    return max(-v.z, 1e-3) / cam.proj_scale;
}

fn skip_fillet(p: vec4<u32>, view_center: vec3<f32>) -> bool {
    let px = SKIP_MARGIN * world_per_px(view_center);
    let small = params.probe_radius < PROBE_PX * px;
    if (kind_of(p) == KIND_TORUS) {
        return small || fillet_size(first_of(p), p.y) < TORUS_PX * px;
    }
    return small;
}

@vertex
fn vs_ses_surface(@builtin(vertex_index) vertex: u32) -> VsOut {
    let local = visible[vertex / 6u];
    let b = bounds[local];
    let center = (cam.view * vec4<f32>(b.xyz, 1.0)).xyz;
    var out: VsOut;
    out.patch_index = page.base + local;
    let p = patches[local];
    if (kind_of(p) != KIND_CONVEX && skip_fillet(p, center)) {
        // A fillet under a pixel: a degenerate triangle outside the view.
        out.clip = vec4<f32>(2.0, 2.0, 2.0, 1.0);
        out.view_pos = center;
        return out;
    }
    if (wholly_cut(center, b.w)) {
        // Past the cut: it would only discard, and being the nearest
        // thing on its own rays, nothing else can early-Z reject it --
        // skip rasterizing it at all.
        out.clip = vec4<f32>(2.0, 2.0, 2.0, 1.0);
        out.view_pos = center;
        return out;
    }
    let pos = billboard(center, b.w, vertex);
    out.clip = cam.proj * vec4<f32>(pos, 1.0);
    out.view_pos = pos;
    return out;
}


struct SesHit {
    depth: f32,
    view_depth: f32,
    normal_view: vec3<f32>,
    dir_view: vec3<f32>,
    color: vec3<f32>,
    atom: u32,
};

// This patch's own surface, cut open where a plane crosses it: the solid
// cross-section is drawn separately, by every atom's own disc
// (`vs_atom_cap`/`fs_atom_cap` below) -- a patch only ever covers the
// solvent-exposed shell, so a buried atom has none; recasting just this
// one patch past the cut would leave buried interior pixels an
// inconsistent patchwork of hits and holes.
fn ses_hit(in: VsOut) -> SesHit {
    var ro_v: vec3<f32>;
    var rd_v: vec3<f32>;
    if (cam.projection == 1u) {
        ro_v = in.view_pos;
        rd_v = vec3<f32>(0.0, 0.0, -1.0);
    } else {
        ro_v = vec3<f32>(0.0, 0.0, 0.0);
        rd_v = normalize(in.view_pos);
    }
    let ro = (params.view_inv * vec4<f32>(ro_v, 1.0)).xyz;
    let rd = normalize((params.view_inv * vec4<f32>(rd_v, 0.0)).xyz);
    // An orthographic ray starts on the billboard; step back so the
    // patch's own bounding sphere lies wholly in front.
    let p = patches[in.patch_index - page.base];
    let bound = ses_bound(p);
    var o = ro;
    if (cam.projection == 1u) {
        o = ro - rd * (dot(ro - bound.xyz, -rd) + bound.w + 1.0);
    }
    let view_center = (cam.view * vec4<f32>(bound.xyz, 1.0)).xyz;
    let px = world_per_px(view_center);
    let min_fillet = TORUS_PX * px;
    let whole = params.probe_radius < PROBE_PX * px;
    let s = cast_patch(p, o, rd, 0.0, min_fillet, whole);
    let dir_view = normalize((cam.view * vec4<f32>(rd, 0.0)).xyz);
    if (!s.ok) {
        discard;
    }
    let hit_view = (cam.view * vec4<f32>(o + rd * s.t, 1.0)).xyz;
    if (clip_distance(hit_view) < 0.0) {
        discard;
    }
    let clip = cam.proj * vec4<f32>(hit_view, 1.0);
    let normal_view = normalize((cam.view * vec4<f32>(s.normal, 0.0)).xyz);
    return SesHit(clip.z / clip.w, -hit_view.z, normal_view, dir_view, s.color, s.atom);
}

fn cast_patch(
    p: vec4<u32>,
    o: vec3<f32>,
    rd: vec3<f32>,
    t_min: f32,
    min_fillet: f32,
    whole: bool,
) -> Surf {
    switch (kind_of(p)) {
        case KIND_CONVEX: { return hit_convex(p, o, rd, t_min, min_fillet, whole); }
        case KIND_TORUS: { return hit_torus(p, o, rd, t_min); }
        default: { return hit_concave(p, o, rd, t_min); }
    }
}

struct FsOut {
    @location(0) color: vec4<f32>,
    @location(1) normal: vec4<f32>,
    @location(2) id: u32,
    @builtin(frag_depth) depth: f32,
};

@fragment
fn fs_ses_surface(in: VsOut) -> FsOut {
    let h = ses_hit(in);
    let m = material(params.material0, params.material1, params.material2);
    var out: FsOut;
    out.depth = h.depth;
    out.color = vec4<f32>(shade(h.color, h.normal_view, h.dir_view, m), 1.0);
    out.normal = vec4<f32>(h.normal_view * 0.5 + 0.5, 1.0);
    out.id = params.id_base + h.atom + 1u;
    return out;
}

// Glass keeps one layer: the prepass records which patch owns each
// pixel's nearest front-facing hit (as skin_surface.wgsl does).
struct OwnerOut {
    @location(0) owner: vec2<u32>,
    @builtin(frag_depth) depth: f32,
};

@fragment
fn fs_ses_surface_owner(in: VsOut) -> OwnerOut {
    let h = ses_hit(in);
    if (dot(h.normal_view, h.dir_view) > 0.0) {
        discard;
    }
    return OwnerOut(vec2<u32>(params.id_base, in.patch_index + 1u), h.depth);
}

@group(2) @binding(0) var glass_owner: texture_2d<u32>;

@fragment
fn fs_ses_surface_glass(in: VsOut) -> GlassOut {
    let owner = textureLoad(glass_owner, vec2<i32>(in.clip.xy), 0).xy;
    if (any(owner != vec2<u32>(params.id_base, in.patch_index + 1u))) {
        discard;
    }
    let h = ses_hit(in);
    let m = material(params.material0, params.material1, params.material2);
    return glass(h.color, h.normal_view, h.dir_view, m, h.view_depth, h.depth);
}

// The solid interior a clip plane exposes, one atom at a time: a billboard
// per atom (not per patch, so every atom has one, including buried atoms a
// patch never covers), discarded except where the plane crosses that
// atom's own van der Waals sphere -- the same analytic disc `draw.wgsl`'s
// `fs_sphere` cuts spacefill with, reused here instead of a per-patch
// re-cast (see `ses_hit`'s comment). Cavities and the gaps between atoms
// stay open: nothing is drawn there.
struct AtomCapVs {
    @builtin(position) clip: vec4<f32>,
    @location(0) view_pos: vec3<f32>,
    @location(1) @interpolate(flat) center: vec3<f32>,
    @location(2) @interpolate(flat) radius: f32,
    @location(3) @interpolate(flat) atom: u32,
};

@vertex
fn vs_atom_cap(@builtin(vertex_index) vertex: u32) -> AtomCapVs {
    let idx = vertex / 6u;
    let a = atoms[idx];
    let center = (cam.view * vec4<f32>(a.xyz, 1.0)).xyz;
    var out: AtomCapVs;
    out.center = center;
    out.radius = a.w;
    out.atom = idx;
    if (no_cap_possible(center, a.w)) {
        out.clip = vec4<f32>(2.0, 2.0, 2.0, 1.0); // outside the frustum
        out.view_pos = center;
        return out;
    }
    out.view_pos = billboard(center, a.w, vertex);
    out.clip = cam.proj * vec4<f32>(out.view_pos, 1.0);
    return out;
}

@fragment
fn fs_atom_cap(in: AtomCapVs) -> FsOut {
    var ro: vec3<f32>;
    var rd: vec3<f32>;
    if (cam.projection == 1u) {
        ro = in.view_pos;
        rd = vec3<f32>(0.0, 0.0, -1.0);
    } else {
        ro = vec3<f32>(0.0, 0.0, 0.0);
        rd = normalize(in.view_pos);
    }
    let c = in.center;
    let r = in.radius;
    let oc = ro - c;
    let b = dot(rd, oc);
    let disc = b * b - dot(oc, oc) + r * r;
    if (disc < 0.0) {
        discard;
    }
    let t_near = -b - sqrt(disc);
    let kept = kept_span(ro, rd, t_near, -b + sqrt(disc));
    if (kept.start > kept.end || kept.start <= t_near) {
        // Not cut here: no cap, and the patch surface already draws the
        // uncut part of this atom's contribution.
        discard;
    }
    var t = kept.start;
    let q = ro + rd * t - c;
    t -= cap_nudge(r * r - dot(q, q), t, dot(in.view_pos - ro, rd));
    if (t <= 0.0) {
        discard;
    }
    let p = ro + rd * t;
    let clip = cam.proj * vec4<f32>(p, 1.0);
    var m = material(params.material0, params.material1, params.material2);
    m.specular = 0.0;
    var out: FsOut;
    out.depth = clip.z / clip.w;
    out.color = vec4<f32>(shade(atom_color(in.atom), kept.normal, rd, m), 1.0);
    out.normal = vec4<f32>(kept.normal * 0.5 + 0.5, 1.0);
    out.id = params.id_base + in.atom + 1u;
    return out;
}
