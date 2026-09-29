// Bakes the Gaussian ("blobby") density field into a 3D texture, so
// `gaussian_surface.wgsl` ray-marches one filtered texture fetch per step
// instead of summing every nearby atom at every step (Krone,
// Stone, Ertl & Schulten, EuroVis 2012). Two passes:
//
// - `scatter`, one thread per atom: adds the atom's Gaussian
//   (`vv_core::gaussian_surface`, 1.0 at its van der Waals surface) to
//   every voxel within its cutoff, and bids for the voxel's colour with
//   how much it contributes. Few voxels per atom, so this is far cheaper
//   than every voxel summing every atom of its neighbourhood.
// - `resolve`, one thread per voxel: rgb = colour of the atom that
//   contributed most, a = density; flags every 8^3 brick holding a voxel
//   at or above the isovalue, so the march can skip the rest.
//
// WGSL has no float atomics: density adds up in fixed point, and the
// colour bid packs (contribution, atom) into one u32 for atomicMax.

struct Params {
    vol_origin: vec3<f32>,
    voxel: f32,
    dims: vec3<u32>,
    isovalue: f32,
    brick_dims: vec3<u32>,
    blob_factor: f32,
    // An atom's cutoff over its radius (`vv_core::gaussian_surface::
    // cutoff_radius`).
    cutoff_ratio: f32,
    atom_count: u32,
    // Bits of the colour bid holding the atom index.
    atom_bits: u32,
    _pad: u32,
};

// xyz = position, w = van der Waals radius.
@group(0) @binding(0) var<storage, read> atoms: array<vec4<f32>>;
@group(0) @binding(1) var<storage, read> colors: array<u32>;
@group(0) @binding(2) var<uniform> params: Params;
// Per voxel: density in fixed point, and the colour bid (0: no atom).
@group(0) @binding(3) var<storage, read_write> density: array<atomic<u32>>;
@group(0) @binding(4) var<storage, read_write> owner: array<atomic<u32>>;
@group(0) @binding(5) var volume: texture_storage_3d<rgba16float, write>;
@group(0) @binding(6) var<storage, read_write> bricks: array<atomic<u32>>;

const BRICK: u32 = 8u;
// Density units per 1.0: fine enough that rounding moves no surface.
const FIXED: f32 = 1048576.0;
// A contribution counts up to this in a colour bid.
const BID_MAX: f32 = 8.0;

// Linear colour: the palette's sRGB value decoded (as `srgb_to_linear`
// in shading.wgsl), so the volume blends and the march shades in linear
// light.
fn unpack_color(c: u32) -> vec3<f32> {
    let s = vec3<f32>(f32(c & 0xffu), f32((c >> 8u) & 0xffu), f32((c >> 16u) & 0xffu)) / 255.0;
    return select(pow((s + 0.055) / 1.055, vec3<f32>(2.4)), s / 12.92, s <= vec3<f32>(0.04045));
}

fn voxel_index(v: vec3<u32>) -> u32 {
    return v.x + v.y * params.dims.x + v.z * params.dims.x * params.dims.y;
}

@compute @workgroup_size(64)
fn scatter(@builtin(global_invocation_id) gid: vec3<u32>) {
    let i = gid.x;
    if (i >= params.atom_count) {
        return;
    }
    let a = atoms[i];
    let cutoff = a.w * params.cutoff_ratio;
    let lo = vec3<i32>(ceil((a.xyz - cutoff - params.vol_origin) / params.voxel));
    let hi = vec3<i32>(floor((a.xyz + cutoff - params.vol_origin) / params.voxel));
    let first = vec3<u32>(max(lo, vec3<i32>(0)));
    let last = vec3<u32>(min(hi, vec3<i32>(params.dims) - 1));
    let levels = f32((1u << (32u - params.atom_bits)) - 1u);
    for (var z = first.z; z <= last.z; z++) {
        for (var y = first.y; y <= last.y; y++) {
            for (var x = first.x; x <= last.x; x++) {
                let p = params.vol_origin + vec3<f32>(f32(x), f32(y), f32(z)) * params.voxel;
                let d = p - a.xyz;
                let d2 = dot(d, d);
                if (d2 > cutoff * cutoff) {
                    continue;
                }
                let g = exp(-params.blob_factor * (d2 / (a.w * a.w) - 1.0));
                let v = voxel_index(vec3<u32>(x, y, z));
                atomicAdd(&density[v], u32(g * FIXED));
                // At least 1, so a bid is never "no atom".
                let bid = max(u32(min(g / BID_MAX, 1.0) * levels), 1u);
                atomicMax(&owner[v], (bid << params.atom_bits) | i);
            }
        }
    }
}

fn flag_brick(b: vec3<u32>) {
    let i = b.x + b.y * params.brick_dims.x + b.z * params.brick_dims.x * params.brick_dims.y;
    atomicOr(&bricks[i], 1u);
}

@compute @workgroup_size(4, 4, 4)
fn resolve(@builtin(global_invocation_id) gid: vec3<u32>) {
    if (any(gid >= params.dims)) {
        return;
    }
    let v = voxel_index(gid);
    let sum = f32(atomicLoad(&density[v])) / FIXED;
    let bid = atomicLoad(&owner[v]);
    let atom = bid & ((1u << params.atom_bits) - 1u);
    let color = select(vec3<f32>(0.0), unpack_color(colors[atom]), bid != 0u);
    textureStore(volume, gid, vec4<f32>(color, min(sum, 60000.0)));

    // Trilinear filtering between voxel v and v + 1 can reach the
    // isovalue only if one of them does, and the world-space region of
    // brick b interpolates voxels 8b ..= 8b + 8. So a voxel at 8b + 0
    // also belongs to brick b - 1 along that axis.
    if (sum >= params.isovalue) {
        let own = gid / BRICK;
        let edge = (gid % BRICK == vec3<u32>(0u)) & (gid > vec3<u32>(0u));
        let hi_b = min(own, params.brick_dims - 1u);
        let lo_b = min(select(own, own - 1u, edge), hi_b);
        for (var z = lo_b.z; z <= hi_b.z; z++) {
            for (var y = lo_b.y; y <= hi_b.y; y++) {
                for (var x = lo_b.x; x <= hi_b.x; x++) {
                    flag_brick(vec3<u32>(x, y, z));
                }
            }
        }
    }
}
