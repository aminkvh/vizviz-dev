// Sphere impostors (ray-cast quads), single-pixel points, and bond
// cylinders (ray-cast capped cylinders on oriented quads).

@group(1) @binding(2) var<storage, read> visible_quads: array<u32>;
@group(1) @binding(3) var<storage, read> visible_points: array<u32>;
@group(1) @binding(4) var<storage, read> bonds: array<vec2<u32>>;
@group(1) @binding(5) var<storage, read> visible_bonds: array<u32>;
// Index of each bond in the source BondTable, for picking (see fs_cylinder).
@group(1) @binding(6) var<storage, read> bond_ids: array<u32>;

// Pick ids are written as index + 1 so that 0 (the cleared value) means
// "nothing here" on every backend. Set on a picked id when it names a bond
// rather than an atom; the remaining 31 bits are the frame-wide bond id +
// 1. Must match `vv_render::renderer::BOND_ID_FLAG`.
const BOND_ID_FLAG: u32 = 0x80000000u;

struct SphereVs {
    @builtin(position) clip: vec4<f32>,
    @location(0) view_pos: vec3<f32>,
    @location(1) @interpolate(flat) center: vec3<f32>,
    @location(2) @interpolate(flat) radius: f32,
    @location(3) @interpolate(flat) color: vec4<f32>,
    // Frame-wide atom id, for the picking target; see fs_sphere.
    @location(4) @interpolate(flat) id: u32,
};

struct SphereFs {
    @location(0) color: vec4<f32>,
    @location(1) normal: vec4<f32>,
    // R32Uint picking target; 0 = "nothing here", otherwise id + 1. Read
    // back on click only (Renderer::pick), never every frame.
    @location(2) id: u32,
    @builtin(frag_depth) depth: f32,
};

// Two triangles: (-1,-1) (1,-1) (1,1) / (-1,-1) (1,1) (-1,1).
fn corner(vi: u32) -> vec2<f32> {
    let cx = select(-1.0, 1.0, vi == 1u || vi == 2u || vi == 4u);
    let cy = select(-1.0, 1.0, vi == 2u || vi == 4u || vi == 5u);
    return vec2<f32>(cx, cy);
}

// Enlargement so a billboard on the near tangent plane always contains
// the perspective silhouette of a sphere of radius r at distance d.
// Perspective only: an orthographic silhouette is always the sphere's
// true equatorial circle, at any distance, so callers use `1.0` instead
// in that mode rather than calling this at all.
fn silhouette_enlarge(d: f32, r: f32) -> f32 {
    return select(d / sqrt(max(d * d - r * r, 1e-6)), 4.0, d <= r * 1.001);
}

@vertex
fn vs_sphere(@builtin(vertex_index) vertex: u32) -> SphereVs {
    // Non-instanced: six consecutive vertices make one atom's quad.
    let idx = visible_quads[vertex / 6u];
    let a = atoms[idx];
    let center = (cam.view * vec4<f32>(a.xyz, 1.0)).xyz;
    let r = atom_radius(a);
    let c = corner(vertex % 6u);

    // "Toward the eye" and the silhouette enlarge factor are both
    // perspective concepts (a billboard facing a *point*, sized for that
    // point's foreshortening); orthographic has no such point, only a
    // fixed view direction, and no foreshortening at any distance -- see
    // `cam.projection`'s doc.
    var p: vec3<f32>;
    if (cam.projection == 1u) {
        p = center + vec3<f32>(c.x, c.y, 1.0) * r;
    } else {
        let d = length(center);
        let toward_eye = -center / max(d, 1e-6);
        // At most the tangent plane, but never nearer than 1.5 near:
        // past the near plane (never clipped) yet before the camera's cut
        // at 2 near (a cap stays behind it, as early depth needs). A
        // screen-filling quad when the eye is inside the sphere.
        let front = min(r, d - cam.near * 1.5);
        p = center + toward_eye * front + vec3<f32>(c.x, c.y, 0.0) * (r * silhouette_enlarge(d, r));
        if (d <= r * 1.001) {
            p = vec3<f32>(c.x * 1e3, c.y * 1e3, -1.0) * (cam.near * 1.5);
        }
    }

    var out: SphereVs;
    out.clip = cam.proj * vec4<f32>(p, 1.0);
    out.view_pos = p;
    out.center = center;
    out.radius = r;
    out.color = linear_rgba(unpack_color(colors[idx]));
    out.id = params.atom_id_base + idx + 1u;
    return out;
}

// Perspective rays fan out from the eye at the view-space origin
// (`ro = 0`); orthographic rays are parallel (`rd` fixed), each starting
// from its own fragment's already-correct 3D position on the billboard
// (perspective-correct interpolation reduces to plain linear
// interpolation when `clip.w` is always 1, as the orthographic
// projection makes it -- so `view_pos` is already the right point on the
// right ray, no reconstruction needed). See `cam.projection`'s doc.
fn eye_ray(view_pos: vec3<f32>) -> array<vec3<f32>, 2> {
    if (cam.projection == 1u) {
        return array<vec3<f32>, 2>(view_pos, vec3<f32>(0.0, 0.0, -1.0));
    }
    return array<vec3<f32>, 2>(vec3<f32>(0.0), normalize(view_pos));
}

struct SphereHit {
    hit: bool,
    t: f32,
    n: vec3<f32>,
    m: Material,
};

// Shared by `fs_sphere`, `fs_sphere_glass` and `fs_sphere_pick`: the
// ray-cast hit (possibly a clip-plane cap) of one atom's sphere, or
// `hit = false` where the caller must discard.
fn sphere_hit(ro: vec3<f32>, rd: vec3<f32>, view_pos: vec3<f32>, c: vec3<f32>, r: f32, m_in: Material) -> SphereHit {
    var out: SphereHit;
    out.hit = false;
    // General ray-sphere intersection (reduces exactly to the old
    // `ro = 0`-specialized form when `ro` is zero: `oc = -c`, `b = -dot
    // (rd, c)`, same discriminant, same `t`).
    let oc = ro - c;
    let b = dot(rd, oc);
    let disc = b * b - dot(oc, oc) + r * r;
    if (disc < 0.0) {
        return out;
    }
    var t = -b - sqrt(disc);
    var n = normalize(ro + rd * t - c);
    var m = m_in;
    // Cut: where the ray first reaches the kept side inside the sphere,
    // draw the flat cross-section.
    let kept = kept_span(ro, rd, t, -b + sqrt(disc));
    if (kept.start > kept.end) {
        return out;
    }
    if (kept.start > t) {
        t = kept.start;
        n = kept.normal;
        // A flat cut facing the light would mirror it: keep caps matte.
        m.specular = 0.0;
        let q = ro + rd * t - c;
        t -= cap_nudge(r * r - dot(q, q), t, dot(view_pos - ro, rd));
    }
    if (t <= 0.0) {
        return out;
    }
    out.hit = true;
    out.t = t;
    out.n = n;
    out.m = m;
    return out;
}

@fragment
fn fs_sphere(in: SphereVs) -> SphereFs {
    let eye = eye_ray(in.view_pos);
    let m0 = material(params.material0, params.material1, params.material2);
    let h = sphere_hit(eye[0], eye[1], in.view_pos, in.center, in.radius, m0);
    if (!h.hit) {
        discard;
    }
    let p = eye[0] + eye[1] * h.t;
    let clip = cam.proj * vec4<f32>(p, 1.0);

    var out: SphereFs;
    out.depth = clip.z / clip.w;
    out.color = vec4<f32>(shade(in.color.rgb, h.n, eye[1], h.m), 1.0);
    out.normal = vec4<f32>(h.n * 0.5 + 0.5, 1.0);
    out.id = in.id;
    return out;
}

// Transparent: weighted-blended OIT into `accum`/`reveal` instead of the
// opaque targets, depth-tested against the opaque pass (see
// `make_glass_pipeline` in renderer.rs) but writing no depth of its own.
@fragment
fn fs_sphere_glass(in: SphereVs) -> GlassOut {
    let eye = eye_ray(in.view_pos);
    let m0 = material(params.material0, params.material1, params.material2);
    let h = sphere_hit(eye[0], eye[1], in.view_pos, in.center, in.radius, m0);
    if (!h.hit) {
        discard;
    }
    let p = eye[0] + eye[1] * h.t;
    let clip = cam.proj * vec4<f32>(p, 1.0);
    return glass(in.color.rgb, h.n, eye[1], h.m, -p.z, clip.z / clip.w);
}

struct PickDepthOut {
    @location(0) id: u32,
    @builtin(frag_depth) depth: f32,
};

// A transparent atom's id, for the pick/selection prepass
// (`renderer.rs`'s "glass id" pass): only reached where no opaque
// fragment already claimed the pixel this frame (see `RenderSettings`'
// picking policy in docs/RENDERING.md), and only the nearest such atom
// wins, exactly as an opaque pick would.
@fragment
fn fs_sphere_pick(in: SphereVs) -> PickDepthOut {
    let eye = eye_ray(in.view_pos);
    let m0 = material(params.material0, params.material1, params.material2);
    let h = sphere_hit(eye[0], eye[1], in.view_pos, in.center, in.radius, m0);
    if (!h.hit) {
        discard;
    }
    let p = eye[0] + eye[1] * h.t;
    let clip = cam.proj * vec4<f32>(p, 1.0);
    var out: PickDepthOut;
    out.id = in.id;
    out.depth = clip.z / clip.w;
    return out;
}

struct PointVs {
    @builtin(position) clip: vec4<f32>,
    @location(0) @interpolate(flat) color: vec4<f32>,
    @location(1) @interpolate(flat) id: u32,
    @location(2) @interpolate(flat) view_depth: f32,
};

struct PointFs {
    @location(0) color: vec4<f32>,
    @location(1) normal: vec4<f32>,
    @location(2) id: u32,
};

@vertex
fn vs_point(@builtin(vertex_index) vi: u32) -> PointVs {
    let idx = visible_points[vi];
    let a = atoms[idx];
    let view_pos = (cam.view * vec4<f32>(a.xyz, 1.0)).xyz;
    var out: PointVs;
    out.clip = cam.view_proj * vec4<f32>(a.xyz, 1.0);
    if (clip_distance(view_pos) < 0.0) {
        out.clip = vec4<f32>(2.0, 2.0, 2.0, 1.0); // outside the frustum
    }
    out.color = linear_rgba(unpack_color(colors[idx]));
    out.id = params.atom_id_base + idx + 1u;
    out.view_depth = -view_pos.z;
    return out;
}

@fragment
fn fs_point(in: PointVs) -> PointFs {
    var out: PointFs;
    out.color = vec4<f32>(shade_point(in.color.rgb, material(params.material0, params.material1, params.material2)), 1.0);
    out.normal = vec4<f32>(0.5, 0.5, 1.0, 1.0);
    out.id = in.id;
    return out;
}

@fragment
fn fs_point_glass(in: PointVs) -> GlassColor {
    return glass_point(in.color.rgb, material(params.material0, params.material1, params.material2), in.view_depth);
}

@fragment
fn fs_point_pick(in: PointVs) -> @location(0) u32 {
    return in.id;
}

struct CylinderVs {
    @builtin(position) clip: vec4<f32>,
    @location(0) view_pos: vec3<f32>,
    @location(1) @interpolate(flat) start: vec3<f32>,
    @location(2) @interpolate(flat) axis: vec3<f32>,
    @location(3) @interpolate(flat) length: f32,
    @location(4) @interpolate(flat) color_a: vec4<f32>,
    @location(5) @interpolate(flat) color_b: vec4<f32>,
    @location(6) @interpolate(flat) id: u32,
};

@vertex
fn vs_cylinder(@builtin(vertex_index) vertex: u32) -> CylinderVs {
    let bi = visible_bonds[vertex / 6u];
    let bnd = bonds[bi];
    let pa = (cam.view * vec4<f32>(atoms[bnd.x].xyz, 1.0)).xyz;
    let pb = (cam.view * vec4<f32>(atoms[bnd.y].xyz, 1.0)).xyz;
    let r = params.bond_radius;
    let m = (pa + pb) * 0.5;
    let d = pb - pa;
    let len = length(d);
    let e = d / max(len, 1e-6);

    // Billboard facing the eye, aligned with the bond on screen. `v` points
    // from the eye toward the bond.
    var v: vec3<f32>;
    var enlarge: f32;
    if (cam.projection == 1u) {
        v = vec3<f32>(0.0, 0.0, -1.0);
        enlarge = 1.0;
    } else {
        let dist = length(m);
        v = m / max(dist, 1e-6);
        let bound = len * 0.5 + r;
        enlarge = silhouette_enlarge(dist, bound);
    }
    var along = e - v * dot(e, v);
    let along_len = length(along);
    var across: vec3<f32>;
    if (along_len < 1e-4) {
        // Looking straight down the bond: any perpendicular frame works.
        across = cross(v, vec3<f32>(0.0, 1.0, 0.0));
        if (length(across) < 1e-4) {
            across = cross(v, vec3<f32>(1.0, 0.0, 0.0));
        }
        across = normalize(across);
        along = normalize(cross(across, v));
    } else {
        along = along / along_len;
        across = normalize(cross(along, v));
    }
    let half_along = (len * 0.5 * along_len + r) * enlarge;
    let half_across = r * enlarge;
    let c = corner(vertex % 6u);
    // In front of the whole cylinder (tilted ends included), so its
    // ray-cast depth is never nearer than the billboard's: conservative
    // depth. A nearer quad only covers more of the screen.
    let front = min(r + len * 0.5 * abs(dot(e, v)), max(dot(m, v) - cam.near * 1.5, 0.0));
    let p = m - v * front + along * (c.x * half_along) + across * (c.y * half_across);

    var out: CylinderVs;
    out.clip = cam.proj * vec4<f32>(p, 1.0);
    out.view_pos = p;
    out.start = pa;
    out.axis = e;
    out.length = len;
    out.color_a = linear_rgba(unpack_color(colors[bnd.x]));
    out.color_b = linear_rgba(unpack_color(colors[bnd.y]));
    out.id = BOND_ID_FLAG | (params.bond_id_base + bond_ids[bi] + 1u);
    return out;
}

struct CylinderHit {
    hit: bool,
    t: f32,
    n: vec3<f32>,
    m: Material,
};

// Shared by `fs_cylinder`, `fs_cylinder_glass` and `fs_cylinder_pick`:
// the ray-cast hit (body, end cap, or a clip-plane cap) of one bond's
// finite cylinder, or `hit = false` where the caller must discard.
fn cylinder_hit(ro: vec3<f32>, rd: vec3<f32>, view_pos: vec3<f32>, a: vec3<f32>, e: vec3<f32>, len: f32, r: f32, m_in: Material) -> CylinderHit {
    var out: CylinderHit;
    out.hit = false;
    // Infinite cylinder: components perpendicular to the axis. General
    // ray form (reduces exactly to the old `ro = 0` case: `oc = -a`,
    // matching before).
    let oc = ro - a;
    let rd_p = rd - e * dot(rd, e);
    let oc_p = oc - e * dot(oc, e);
    let qa = dot(rd_p, rd_p);
    let qb = 2.0 * dot(rd_p, oc_p);
    let qc = dot(oc_p, oc_p) - r * r;
    var t = -1.0;
    var n = vec3<f32>(0.0);
    if (qa > 1e-8) {
        let disc = qb * qb - 4.0 * qa * qc;
        if (disc >= 0.0) {
            let tt = (-qb - sqrt(disc)) / (2.0 * qa);
            let p = ro + rd * tt;
            let s = dot(p - a, e);
            if (tt > 0.0 && s >= 0.0 && s <= len) {
                t = tt;
                n = normalize((p - a) - e * s);
            }
        }
    }
    if (t < 0.0) {
        // End caps: the nearer disk hit, if any. `ta`/`tb` solve
        // `dot((ro + rd*t) - a, e) == 0` generally (matches the old
        // `dot(a, e) / denom` exactly when `ro = 0`).
        let denom = dot(rd, e);
        if (abs(denom) > 1e-6) {
            let b = a + e * len;
            let ta = dot(a - ro, e) / denom;
            let tb = dot(b - ro, e) / denom;
            let pa_hit = ro + rd * ta;
            let pb_hit = ro + rd * tb;
            let hit_a = ta > 0.0 && dot(pa_hit - a, pa_hit - a) <= r * r;
            let hit_b = tb > 0.0 && dot(pb_hit - b, pb_hit - b) <= r * r;
            if (hit_a && (!hit_b || ta < tb)) {
                t = ta;
                n = -e;
            } else if (hit_b) {
                t = tb;
                n = e;
            }
            if (dot(n, rd) > 0.0) {
                n = -n;
            }
        }
    }
    if (t < 0.0) {
        return out;
    }
    var m = m_in;
    if (clip_distance(ro + rd * t) < 0.0) {
        // Capped where the ray reaches the kept side, if that is still
        // inside the (convex) cylinder.
        let kept = kept_span(ro, rd, t, 3.4e38);
        let q = ro + rd * kept.start;
        let along = dot(q - a, e);
        let off_axis = (q - a) - e * along;
        if (kept.start > kept.end || along < 0.0 || along > len || dot(off_axis, off_axis) > r * r) {
            return out;
        }
        t = kept.start - cap_nudge(r * r - dot(off_axis, off_axis), kept.start, dot(view_pos - ro, rd));
        n = kept.normal;
        m.specular = 0.0;
    }
    out.hit = true;
    out.t = t;
    out.n = n;
    out.m = m;
    return out;
}

// Same `ro`/`rd` split as `fs_sphere` -- see its comment and
// `cam.projection`'s doc.
fn cylinder_eye_ray(in: CylinderVs) -> array<vec3<f32>, 2> {
    if (cam.projection == 1u) {
        // Start the parallel ray in front of the whole cylinder, not at
        // the billboard: the billboard sits one radius in front of the
        // bond's *middle*, so a bond tilted toward the viewer has its near
        // half in front of it, and starting there cut that half off
        // (pencil-tip sticks, notched tubes).
        let ro = vec3<f32>(in.view_pos.xy, in.start.z + in.length + params.bond_radius * 2.0 + 1.0);
        return array<vec3<f32>, 2>(ro, vec3<f32>(0.0, 0.0, -1.0));
    }
    return array<vec3<f32>, 2>(vec3<f32>(0.0), normalize(in.view_pos));
}

@fragment
fn fs_cylinder(in: CylinderVs) -> SphereFs {
    let eye = cylinder_eye_ray(in);
    let m0 = material(params.material0, params.material1, params.material2);
    let h = cylinder_hit(eye[0], eye[1], in.view_pos, in.start, in.axis, in.length, params.bond_radius, m0);
    if (!h.hit) {
        discard;
    }
    let p = eye[0] + eye[1] * h.t;
    let s = dot(p - in.start, in.axis);
    let color = select(in.color_b, in.color_a, s < in.length * 0.5);
    let clip = cam.proj * vec4<f32>(p, 1.0);

    var out: SphereFs;
    out.depth = clip.z / clip.w;
    out.color = vec4<f32>(shade(color.rgb, h.n, eye[1], h.m), 1.0);
    out.normal = vec4<f32>(h.n * 0.5 + 0.5, 1.0);
    out.id = in.id;
    return out;
}

@fragment
fn fs_cylinder_glass(in: CylinderVs) -> GlassOut {
    let eye = cylinder_eye_ray(in);
    let m0 = material(params.material0, params.material1, params.material2);
    let h = cylinder_hit(eye[0], eye[1], in.view_pos, in.start, in.axis, in.length, params.bond_radius, m0);
    if (!h.hit) {
        discard;
    }
    let p = eye[0] + eye[1] * h.t;
    let s = dot(p - in.start, in.axis);
    let color = select(in.color_b, in.color_a, s < in.length * 0.5);
    let clip = cam.proj * vec4<f32>(p, 1.0);
    return glass(color.rgb, h.n, eye[1], h.m, -p.z, clip.z / clip.w);
}

@fragment
fn fs_cylinder_pick(in: CylinderVs) -> PickDepthOut {
    let eye = cylinder_eye_ray(in);
    let m0 = material(params.material0, params.material1, params.material2);
    let h = cylinder_hit(eye[0], eye[1], in.view_pos, in.start, in.axis, in.length, params.bond_radius, m0);
    if (!h.hit) {
        discard;
    }
    let p = eye[0] + eye[1] * h.t;
    let clip = cam.proj * vec4<f32>(p, 1.0);
    var out: PickDepthOut;
    out.id = in.id;
    out.depth = clip.z / clip.w;
    return out;
}

// A bond too thin for a cylinder (see `cull_bonds`): a 1-px line between
// its atoms, filled from the end of `visible_bonds`. Two vertices per
// bond; color blends from one atom's to the other's.
struct LineVs {
    @builtin(position) clip: vec4<f32>,
    @location(0) color: vec4<f32>,
    @location(1) @interpolate(flat) id: u32,
    @location(2) view_depth: f32,
};

@vertex
fn vs_line(@builtin(vertex_index) vertex: u32) -> LineVs {
    let n = arrayLength(&visible_bonds);
    let bi = visible_bonds[n - 1u - vertex / 2u];
    let bnd = bonds[bi];
    let atom = select(bnd.x, bnd.y, vertex % 2u == 1u);
    let view_pos = (cam.view * vec4<f32>(atoms[atom].xyz, 1.0)).xyz;
    var out: LineVs;
    out.clip = cam.view_proj * vec4<f32>(atoms[atom].xyz, 1.0);
    let ends_view = array<vec3<f32>, 2>(
        (cam.view * vec4<f32>(atoms[bnd.x].xyz, 1.0)).xyz,
        (cam.view * vec4<f32>(atoms[bnd.y].xyz, 1.0)).xyz,
    );
    if (min(clip_distance(ends_view[0]), clip_distance(ends_view[1])) < 0.0) {
        out.clip = vec4<f32>(2.0, 2.0, 2.0, 1.0); // outside the frustum
    }
    out.color = linear_rgba(unpack_color(colors[atom]));
    out.id = BOND_ID_FLAG | (params.bond_id_base + bond_ids[bi] + 1u);
    out.view_depth = -view_pos.z;
    return out;
}

@fragment
fn fs_line(in: LineVs) -> PointFs {
    var out: PointFs;
    out.color = vec4<f32>(shade_point(in.color.rgb, material(params.material0, params.material1, params.material2)), 1.0);
    out.normal = vec4<f32>(0.5, 0.5, 1.0, 1.0);
    out.id = in.id;
    return out;
}

@fragment
fn fs_line_glass(in: LineVs) -> GlassColor {
    return glass_point(in.color.rgb, material(params.material0, params.material1, params.material2), in.view_depth);
}

@fragment
fn fs_line_pick(in: LineVs) -> @location(0) u32 {
    return in.id;
}
