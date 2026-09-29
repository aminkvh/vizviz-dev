// Depth pyramid level 0: copy the depth buffer into an R32Float texture.

@group(0) @binding(0) var src_depth: texture_depth_2d;
@group(0) @binding(1) var dst: texture_storage_2d<r32float, write>;

@compute @workgroup_size(8, 8)
fn copy_depth(@builtin(global_invocation_id) id: vec3<u32>) {
    let size = textureDimensions(dst);
    if (id.x >= size.x || id.y >= size.y) {
        return;
    }
    let d = textureLoad(src_depth, id.xy, 0);
    textureStore(dst, id.xy, vec4<f32>(d, 0.0, 0.0, 0.0));
}
