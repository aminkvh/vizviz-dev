// Selection outline (`vv_render::selection`): colours every pixel that is
// not selected itself but has a selected pixel within `width` pixels (a
// rounded disc), so the halo sits just outside the selection and never
// covers it. The vertex shader comes from fullscreen.wgsl.

@group(0) @binding(0) var ids: texture_2d<u32>;

// Frame-wide bitsets over pick ids minus one; see `renderer::ID_FORMAT`.
@group(1) @binding(0) var<storage, read> atom_bits: array<u32>;
@group(1) @binding(1) var<storage, read> bond_bits: array<u32>;

// Layout must match `StyleUniform` in selection.rs.
struct Style {
    // Linear RGBA.
    color: vec4<f32>,
    // Halo width in target pixels.
    width: u32,
};
@group(1) @binding(2) var<uniform> style: Style;

// Must match `renderer::BOND_ID_FLAG`.
const BOND_ID_FLAG: u32 = 0x80000000u;

// Ids past the bitsets (the `NO_PICK_ID_BASE` range) are never selected.
fn is_selected(id: u32) -> bool {
    if (id == 0u) {
        return false;
    }
    let i = (id & ~BOND_ID_FLAG) - 1u;
    let word = i >> 5u;
    let bit = 1u << (i & 31u);
    if ((id & BOND_ID_FLAG) != 0u) {
        return word < arrayLength(&bond_bits) && (bond_bits[word] & bit) != 0u;
    }
    return word < arrayLength(&atom_bits) && (atom_bits[word] & bit) != 0u;
}

@fragment
fn fs_selection_outline(in: Vs) -> @location(0) vec4<f32> {
    let dims = vec2<i32>(textureDimensions(ids));
    let p = vec2<i32>(in.clip.xy);
    if (is_selected(textureLoad(ids, p, 0).r)) {
        discard;
    }
    let w = i32(style.width);
    let reach = w * w + w;
    for (var dy = -w; dy <= w; dy++) {
        for (var dx = -w; dx <= w; dx++) {
            let q = p + vec2<i32>(dx, dy);
            if (dx * dx + dy * dy > reach || any(q < vec2<i32>(0)) || any(q >= dims)) {
                continue;
            }
            if (is_selected(textureLoad(ids, q, 0).r)) {
                return style.color;
            }
        }
    }
    discard;
}
