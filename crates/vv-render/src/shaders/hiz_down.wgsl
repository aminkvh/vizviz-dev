// Depth pyramid level i from level i-1: min of the 2x2 footprint, folding
// the leftover column/row of odd-sized sources into the last texel so the
// pyramid never claims a region is nearer than it is (reversed-Z: smaller
// is farther).

@group(0) @binding(0) var src: texture_2d<f32>;
@group(0) @binding(1) var dst: texture_storage_2d<r32float, write>;

fn load(c: vec2<u32>, size: vec2<u32>) -> f32 {
    return textureLoad(src, min(c, size - 1u), 0).r;
}

@compute @workgroup_size(8, 8)
fn downsample(@builtin(global_invocation_id) id: vec3<u32>) {
    let dst_size = textureDimensions(dst);
    if (id.x >= dst_size.x || id.y >= dst_size.y) {
        return;
    }
    let src_size = textureDimensions(src);
    let base = id.xy * 2u;
    var m = load(base, src_size);
    m = min(m, load(base + vec2<u32>(1u, 0u), src_size));
    m = min(m, load(base + vec2<u32>(0u, 1u), src_size));
    m = min(m, load(base + vec2<u32>(1u, 1u), src_size));

    let extra_x = (src_size.x & 1u) == 1u && id.x == dst_size.x - 1u;
    let extra_y = (src_size.y & 1u) == 1u && id.y == dst_size.y - 1u;
    if (extra_x) {
        let x = src_size.x - 1u;
        m = min(m, load(vec2<u32>(x, base.y), src_size));
        m = min(m, load(vec2<u32>(x, base.y + 1u), src_size));
    }
    if (extra_y) {
        let y = src_size.y - 1u;
        m = min(m, load(vec2<u32>(base.x, y), src_size));
        m = min(m, load(vec2<u32>(base.x + 1u, y), src_size));
    }
    if (extra_x && extra_y) {
        m = min(m, load(src_size - 1u, src_size));
    }
    textureStore(dst, id.xy, vec4<f32>(m, 0.0, 0.0, 0.0));
}
