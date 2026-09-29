// Frustum + occlusion cull and size classification. One thread per atom
// (or bond); survivors take a slot from a workgroup-local atomic counter,
// then one thread claims the workgroup's range in the visible-index
// buffers with a single global atomic add per bucket. (A shared-memory
// prefix scan was used first and silently dropped wave tails on AMD;
// atomics are robust.)

// Two-phase occlusion culling. Phase 1 tests against the depth pyramid
// of the previous frame and records what it rejects; the renderer then
// rebuilds the pyramid from what phase 1 drew, and phase 2 re-tests only
// the rejected atoms/bonds against it, drawing the ones that are visible
// after all (newly disoccluded by camera motion). Without phase 2,
// rotating culls visible atoms and leaves holes until the camera stops.
override PHASE: u32 = 1u;

// Depth pyramid (reversed-Z, min-reduced), so texel = farthest depth in
// its footprint: last frame's in phase 1, this frame's in phase 2.
@group(0) @binding(1) var hiz: texture_2d<f32>;
// Last frame's nearest-glass depth (`Renderer::render_all`'s "glass id"
// pass, copied out before it is reused as Skin/SES scratch), one texel
// per pixel, no mip chain -- see `glass_hidden`.
@group(0) @binding(3) var glass_hiz: texture_2d<f32>;

@group(1) @binding(2) var<storage, read_write> visible_quads: array<u32>;
@group(1) @binding(3) var<storage, read_write> visible_points: array<u32>;
// Four DrawIndirectArgs per occlusion phase: quads [0..4), points
// [4..8), bond cylinders [8..12), bond lines [12..16) for phase 1, the
// same +16 for phase 2 (whose entries go after phase 1's in the visible
// lists; its `first_vertex` points past them). Quads and cylinders
// accumulate 6 vertices each, lines 2, points 1.
@group(1) @binding(4) var<storage, read_write> indirect: array<atomic<u32>, 32>;
@group(1) @binding(5) var<storage, read> bonds: array<vec2<u32>>;
@group(1) @binding(6) var<storage, read_write> visible_bonds: array<u32>;
// One bit per atom, then one per bond (from word `bond_bits_start()`),
// that phase 1 rejected as occluded; cleared by the renderer every frame.
@group(1) @binding(8) var<storage, read_write> occluded_bits: array<atomic<u32>>;

// Cull state per CLUSTER_ATOMS consecutive atoms: 0 outside the frustum,
// 1 visible, 2 hidden (phase 1), 3 revealed by phase 2. Atoms in clusters
// that are outside or hidden skip their own tests.
@group(1) @binding(9) var<storage, read_write> cluster_state: array<u32>;
// xyz centre, w radius, at full van der Waals radius.
@group(1) @binding(10) var<storage, read> cluster_bounds: array<vec4<f32>>;

const CLUSTER_ATOMS: u32 = 64u;
const CLUSTER_OUT: u32 = 0u;
const CLUSTER_VISIBLE: u32 = 1u;
const CLUSTER_HIDDEN: u32 = 2u;
const CLUSTER_REVEALED: u32 = 3u;

fn bond_bits_start() -> u32 {
    return (arrayLength(&atoms) + 31u) / 32u;
}

fn was_occluded(word: u32, i: u32) -> bool {
    return ((atomicLoad(&occluded_bits[word]) >> (i & 31u)) & 1u) == 1u;
}

var<workgroup> local_q: atomic<u32>;
var<workgroup> local_p: atomic<u32>;
var<workgroup> base_q: u32;
var<workgroup> base_p: u32;

// Wang hash -> [0, 1). Used to jitter the quad/point threshold per atom so
// a structure whose atoms all project to the same size does not flip
// between "all quads" and "all points" as the threshold moves.
fn hash01(x: u32) -> f32 {
    var h = x;
    h = (h ^ 61u) ^ (h >> 16u);
    h = h * 9u;
    h = h ^ (h >> 4u);
    h = h * 0x27d4eb2du;
    h = h ^ (h >> 15u);
    return f32(h & 0xffffu) / 65536.0;
}

// Reads the projected on-screen size a world-space `size` at `view_depth`
// (a positive distance, `-vc.z`) actually occupies, in pixels -- for
// perspective that shrinks with distance (`cam.proj_scale` is "pixels
// per world unit at view depth 1"), for orthographic it's constant, so
// `cam.projection` decides which. Every LOD/occlusion-margin site below
// used to divide by `view_depth` unconditionally; that's simply wrong
// once `view_depth` no longer changes apparent size (see `cam.
// projection`'s doc).
fn screen_px(size: f32, view_depth: f32) -> f32 {
    if (cam.projection == 1u) {
        return size * cam.proj_scale;
    }
    return size * cam.proj_scale / max(view_depth, 1e-3);
}

// The reversed-Z depth value at view-space `z` (negative, in front of the
// camera) for whichever of `cam.proj`'s two forms is active -- read
// directly off the matrix instead of hardcoding either form's own
// closed-form shortcut, so this needs no `cam.projection` branch at all
// and stays correct if either matrix ever changes shape. Both this
// codebase's projections make `clip.x`/`clip.y` irrelevant to depth, so
// plugging in `x = y = 0` is exact, not an approximation.
fn depth_from_view_z(z: f32) -> f32 {
    let num = cam.proj[2][2] * z + cam.proj[3][2];
    let den = cam.proj[2][3] * z + cam.proj[3][3];
    return num / den;
}

fn in_frustum(c: vec3<f32>, r: f32) -> bool {
    for (var k = 0u; k < 6u; k++) {
        let pl = cam.planes[k];
        if (dot(pl.xyz, c) + pl.w < -r) {
            return false;
        }
    }
    return true;
}

// True when the whole sphere lies behind what the previous frame drew.
// `vc` is the view-space centre. Conservative: the sphere is treated as
// half a radius nearer than it is, so small camera motion between frames
// does not cull visible atoms.
fn occluded(vc: vec3<f32>, r: f32) -> bool {
    let d = length(vc);
    if (d <= r * 1.05) {
        return false;
    }
    let view_depth = -vc.z;
    let near_z = view_depth - r * 1.5;
    if (near_z <= cam.near) {
        return false;
    }
    let sphere_depth = depth_from_view_z(-near_z);

    // Perspective-only silhouette enlarge (see `draw.wgsl`'s
    // `silhouette_enlarge`, same reasoning): orthographic's is always 1.
    var enlarge = 1.0;
    if (cam.projection == 0u) {
        enlarge = d / sqrt(max(d * d - r * r, 1e-6));
    }
    let half_px = screen_px(r * enlarge, view_depth);
    let clip = cam.proj * vec4<f32>(vc, 1.0);
    let ndc = clip.xy / clip.w;
    let center_px = vec2<f32>(
        (ndc.x * 0.5 + 0.5) * cam.viewport.x,
        (0.5 - ndc.y * 0.5) * cam.viewport.y,
    );
    let lo = max(center_px - half_px, vec2<f32>(0.0));
    let hi = min(center_px + half_px, cam.viewport - 1.0);
    if (lo.x > hi.x || lo.y > hi.y) {
        return false;
    }
    // Pick the level where the footprint spans at most 2x2 texels.
    let size_px = max(hi.x - lo.x, hi.y - lo.y);
    let level = min(u32(ceil(log2(max(size_px, 1.0)))), cam.hiz_mip_count - 1u);
    let dims = textureDimensions(hiz, level);
    let scale = 1.0 / f32(1u << level);
    let t0 = min(vec2<u32>(lo * scale), dims - 1u);
    let t1 = min(vec2<u32>(hi * scale), dims - 1u);
    let a = textureLoad(hiz, vec2<u32>(t0.x, t0.y), level).r;
    let b = textureLoad(hiz, vec2<u32>(t1.x, t0.y), level).r;
    let c = textureLoad(hiz, vec2<u32>(t0.x, t1.y), level).r;
    let e = textureLoad(hiz, vec2<u32>(t1.x, t1.y), level).r;
    let farthest_occluder = min(min(a, b), min(c, e));
    return sphere_depth < farthest_occluder;
}

// True when `vc` lies more than `margin` world units behind the nearest
// glass surface the previous frame's "glass id" pass recorded at `vc`'s
// own screen pixel -- past the point where weighted-blended OIT's
// contribution rounds to zero (see `vv_render::scene::glass_cull_margin`
// for `margin` and docs/RENDERING.md's "Transparency" for the
// derivation). Single centre texel, not a footprint query like
// `occluded`: cheaper, and exact for one atom, but `cull_clusters` also
// calls this on a whole cluster's bounding sphere, whose on-screen
// footprint can span many pixels -- a cluster straddling a glass
// silhouette reads only its centre's reference, which is the source of
// the pixel error the quality-bound test in `tests/headless.rs` measures
// (not a footprint query bug to fix, since that test bounds it). A stale
// reference (no glass drawn there last frame, or the scene just gained
// its first glass item) reads 0 and never culls.
fn glass_hidden(vc: vec3<f32>, margin: f32) -> bool {
    let view_depth = -vc.z;
    let near_z = view_depth - margin;
    if (near_z <= cam.near) {
        return false;
    }
    let probe_depth = depth_from_view_z(-near_z);
    let clip = cam.proj * vec4<f32>(vc, 1.0);
    let ndc = clip.xy / clip.w;
    let dims = vec2<i32>(textureDimensions(glass_hiz));
    let px = vec2<i32>(vec2<f32>(
        (ndc.x * 0.5 + 0.5) * cam.viewport.x,
        (0.5 - ndc.y * 0.5) * cam.viewport.y,
    ));
    if (any(px < vec2<i32>(0)) || any(px >= dims)) {
        return false;
    }
    let nearest_glass = textureLoad(glass_hiz, px, 0).r;
    return nearest_glass > 0.0 && probe_depth < nearest_glass;
}

// The glass half of `culled_by_depth`, on its own: `cull`/`cull_bonds`
// need it separately from `occluded` so phase 2 only re-tests atoms
// `occluded` actually rejected (see their own `occluded_bits` use) --
// re-testing a glass-culled atom would just read the same, unchanged
// `glass_hiz` again for the same answer.
fn glass_only_hidden(vc: vec3<f32>) -> bool {
    return cam.occlusion != 0u
        && params.glass_cull_margin > 0.0
        && glass_hidden(vc, params.glass_cull_margin);
}

// True when either the opaque Hi-Z or, for a glass page, the glass
// self-occlusion margin says `vc`/`r` cannot contribute a visible pixel
// (see `occluded` and `glass_hidden`).
fn culled_by_depth(vc: vec3<f32>, r: f32) -> bool {
    if (cam.occlusion == 0u) {
        return false;
    }
    return occluded(vc, r) || glass_only_hidden(vc);
}

// One thread per cluster: frustum and occlusion on its bounding sphere.
@compute @workgroup_size(256)
fn cull_clusters(@builtin(global_invocation_id) gid: vec3<u32>) {
    let c = gid.x;
    if (c >= arrayLength(&cluster_bounds) || (PHASE == 2u && cluster_state[c] != CLUSTER_HIDDEN)) {
        return;
    }
    let center = cluster_bounds[c].xyz;
    let r = cluster_bounds[c].w;
    var state = CLUSTER_VISIBLE;
    if (!in_frustum(center, r)) {
        state = CLUSTER_OUT;
    } else if (culled_by_depth((cam.view * vec4<f32>(center, 1.0)).xyz, r)) {
        state = CLUSTER_HIDDEN;
    }
    if (PHASE == 1u) {
        cluster_state[c] = state;
    } else if (state == CLUSTER_VISIBLE) {
        cluster_state[c] = CLUSTER_REVEALED;
    }
}

@compute @workgroup_size(256)
fn cull(
    @builtin(global_invocation_id) gid: vec3<u32>,
    @builtin(local_invocation_index) lid: u32,
) {
    if (lid == 0u) {
        atomicStore(&local_q, 0u);
        atomicStore(&local_p, 0u);
    }
    workgroupBarrier();

    let n = arrayLength(&atoms);
    let i = gid.x;
    var bucket = 0u;
    // Phase 1: atoms of visible clusters. Phase 2: atoms of revealed
    // clusters, and atoms phase 1 rejected inside visible ones.
    var candidate = false;
    if (i < n) {
        let cs = cluster_state[i / CLUSTER_ATOMS];
        if (PHASE == 1u) {
            candidate = cs == CLUSTER_VISIBLE;
        } else {
            candidate = cs == CLUSTER_REVEALED
                || (cs == CLUSTER_VISIBLE && was_occluded(i >> 5u, i));
        }
    }
    if (candidate) {
        let a = atoms[i];
        let c = a.xyz;
        let r = atom_radius(a);
        if (in_frustum(c, r)) {
            let vc = (cam.view * vec4<f32>(c, 1.0)).xyz;
            let opaque_hidden = cam.occlusion != 0u && occluded(vc, r);
            let hidden = opaque_hidden || glass_only_hidden(vc);
            if (opaque_hidden && PHASE == 1u) {
                atomicOr(&occluded_bits[i >> 5u], 1u << (i & 31u));
            }
            if (!hidden) {
                let px = screen_px(max(r, params.min_quad_radius), -vc.z);
                let threshold = cam.quad_px_threshold * (0.6 + 0.8 * hash01(i));
                bucket = select(2u, 1u, px >= threshold);
            }
        }
    }
    var slot = 0u;
    if (bucket == 1u) {
        slot = atomicAdd(&local_q, 1u);
    } else if (bucket == 2u) {
        slot = atomicAdd(&local_p, 1u);
    }
    workgroupBarrier();

    if (lid == 0u) {
        // Quads are drawn non-instanced: 6 vertices per atom.
        if (PHASE == 1u) {
            base_q = atomicAdd(&indirect[0], atomicLoad(&local_q) * 6u) / 6u;
            base_p = atomicAdd(&indirect[4], atomicLoad(&local_p));
        } else {
            let q1 = atomicLoad(&indirect[0]);
            let p1 = atomicLoad(&indirect[4]);
            atomicStore(&indirect[18], q1);
            atomicStore(&indirect[22], p1);
            base_q = q1 / 6u + atomicAdd(&indirect[16], atomicLoad(&local_q) * 6u) / 6u;
            base_p = p1 + atomicAdd(&indirect[20], atomicLoad(&local_p));
        }
    }
    workgroupBarrier();

    if (bucket == 1u) {
        visible_quads[base_q + slot] = i;
    } else if (bucket == 2u) {
        visible_points[base_p + slot] = i;
    }
}

var<workgroup> local_b: atomic<u32>;
var<workgroup> local_l: atomic<u32>;
var<workgroup> base_b: u32;
var<workgroup> base_l: u32;

// Bonds too thin for a cylinder are drawn as 1-px lines instead of being
// dropped, so ball-and-stick and tube keep their connectivity at any zoom
// (dropping them left big structures as clouds of dots). Lines fill
// `visible_bonds` from the end, cylinders from the start.
@compute @workgroup_size(256)
fn cull_bonds(
    @builtin(global_invocation_id) gid: vec3<u32>,
    @builtin(local_invocation_index) lid: u32,
) {
    if (lid == 0u) {
        atomicStore(&local_b, 0u);
        atomicStore(&local_l, 0u);
    }
    workgroupBarrier();

    let n = arrayLength(&bonds);
    let i = gid.x;
    // 0 = culled, 1 = cylinder, 2 = line.
    var kind = 0u;
    // A bond lies inside its two atoms' cluster spheres, so it is skipped
    // when both clusters are off-screen or hidden, and phase 2 tests it
    // again only if one of them was revealed.
    var candidate = false;
    var bnd = vec2<u32>(0u);
    if (i < n) {
        bnd = bonds[i];
        let ca = cluster_state[bnd.x / CLUSTER_ATOMS];
        let cb = cluster_state[bnd.y / CLUSTER_ATOMS];
        let any_visible = ca == CLUSTER_VISIBLE || cb == CLUSTER_VISIBLE;
        if (PHASE == 1u) {
            candidate = any_visible;
        } else {
            candidate = was_occluded(bond_bits_start() + (i >> 5u), i)
                || (!any_visible && (ca == CLUSTER_REVEALED || cb == CLUSTER_REVEALED));
        }
    }
    if (candidate && params.bond_radius > 0.0) {
        let a = atoms[bnd.x].xyz;
        let b = atoms[bnd.y].xyz;
        let m = (a + b) * 0.5;
        let rad = length(b - a) * 0.5 + params.bond_radius;
        if (in_frustum(m, rad)) {
            let vc = (cam.view * vec4<f32>(m, 1.0)).xyz;
            let opaque_hidden = cam.occlusion != 0u && occluded(vc, rad);
            let hidden = opaque_hidden || glass_only_hidden(vc);
            if (opaque_hidden && PHASE == 1u) {
                atomicOr(&occluded_bits[bond_bits_start() + (i >> 5u)], 1u << (i & 31u));
            }
            if (!hidden) {
                // Cylinder when at least half a pixel thick; else a line,
                // but only if it spans a couple of pixels (shorter ones
                // are covered by their atoms' dots).
                let depth = -vc.z;
                if (params.lines != 0u) {
                    kind = 2u;
                } else if (screen_px(params.bond_radius, depth) >= 0.5) {
                    kind = 1u;
                } else if (screen_px(length(b - a), depth) >= 3.0) {
                    kind = 2u;
                }
            }
        }
    }
    var slot = 0u;
    if (kind == 1u) {
        slot = atomicAdd(&local_b, 1u);
    } else if (kind == 2u) {
        slot = atomicAdd(&local_l, 1u);
    }
    workgroupBarrier();
    if (lid == 0u) {
        // Cylinders: 6 vertices each; lines: 2.
        if (PHASE == 1u) {
            base_b = atomicAdd(&indirect[8], atomicLoad(&local_b) * 6u) / 6u;
            base_l = atomicAdd(&indirect[12], atomicLoad(&local_l) * 2u) / 2u;
        } else {
            let b1 = atomicLoad(&indirect[8]);
            let l1 = atomicLoad(&indirect[12]);
            atomicStore(&indirect[26], b1);
            atomicStore(&indirect[30], l1);
            base_b = b1 / 6u + atomicAdd(&indirect[24], atomicLoad(&local_b) * 6u) / 6u;
            base_l = l1 / 2u + atomicAdd(&indirect[28], atomicLoad(&local_l) * 2u) / 2u;
        }
    }
    workgroupBarrier();
    if (kind == 1u) {
        visible_bonds[base_b + slot] = i;
    } else if (kind == 2u) {
        visible_bonds[n - 1u - (base_l + slot)] = i;
    }
}
