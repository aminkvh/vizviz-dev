//! CPU-side image post-processing for exports: 2x supersampling downsample.
//! Pure pixel math, no GPU state, so it's tested without a GPU adapter.

/// `COLOR_FORMAT` is `Rgba8UnormSrgb`: the RGB channels are sRGB-encoded,
/// alpha is linear Unorm (wgpu's `*Srgb` formats never encode alpha). A
/// plain average of sRGB bytes darkens edges — this decodes to linear,
/// averages, and re-encodes, which is what makes 2x SSAA actually look
/// smoother instead of just blurrier.
fn srgb_to_linear(c: u8) -> f32 {
    let c = c as f32 / 255.0;
    if c <= 0.04045 {
        c / 12.92
    } else {
        ((c + 0.055) / 1.055).powf(2.4)
    }
}

fn linear_to_srgb(c: f32) -> u8 {
    let c = c.clamp(0.0, 1.0);
    let s = if c <= 0.0031308 {
        c * 12.92
    } else {
        1.055 * c.powf(1.0 / 2.4) - 0.055
    };
    (s * 255.0).round() as u8
}

/// Box-filters a tightly packed RGBA8 image (as from `Renderer::read_color`)
/// from `width`x`height` down to half that in each dimension. Both must be
/// even.
pub fn downsample_2x_srgb(pixels: &[u8], width: u32, height: u32) -> Vec<u8> {
    assert_eq!(width % 2, 0, "width must be even to downsample by 2x");
    assert_eq!(height % 2, 0, "height must be even to downsample by 2x");
    assert_eq!(pixels.len(), (width * height * 4) as usize);
    let (out_w, out_h) = (width / 2, height / 2);
    let at = |x: u32, y: u32, c: u32| pixels[((y * width + x) * 4 + c) as usize];

    let mut out = vec![0u8; (out_w * out_h * 4) as usize];
    for oy in 0..out_h {
        for ox in 0..out_w {
            let (x0, y0) = (ox * 2, oy * 2);
            let mut linear_rgb = [0f32; 3];
            let mut alpha_sum = 0u32;
            for dy in 0..2 {
                for dx in 0..2 {
                    let (x, y) = (x0 + dx, y0 + dy);
                    for (c, acc) in linear_rgb.iter_mut().enumerate() {
                        *acc += srgb_to_linear(at(x, y, c as u32));
                    }
                    alpha_sum += at(x, y, 3) as u32;
                }
            }
            let o = ((oy * out_w + ox) * 4) as usize;
            for (c, sum) in linear_rgb.into_iter().enumerate() {
                out[o + c] = linear_to_srgb(sum / 4.0);
            }
            out[o + 3] = (alpha_sum / 4) as u8;
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn uniform_image_downsamples_to_the_same_color() {
        let px = [10u8, 20, 30, 255].repeat(4); // 2x2, one color everywhere
        let out = downsample_2x_srgb(&px, 2, 2);
        assert_eq!(out, [10, 20, 30, 255]);
    }

    #[test]
    fn averaging_is_gamma_correct_not_naive() {
        // Two black, two white -> linear average 0.5, not sRGB average 0.5.
        let black = [0u8, 0, 0, 255];
        let white = [255u8, 255, 255, 255];
        let mut px = Vec::new();
        px.extend_from_slice(&black);
        px.extend_from_slice(&white);
        px.extend_from_slice(&white);
        px.extend_from_slice(&black);
        let out = downsample_2x_srgb(&px, 2, 2);
        // Linear-correct mid-gray from 50% black/white is ~188, well above
        // the naive sRGB-space average of 128.
        assert!(
            out[0] > 170 && out[0] < 200,
            "expected a gamma-correct mid-gray, got {}",
            out[0]
        );
        assert_eq!(out[3], 255, "alpha averages linearly, no gamma curve");
    }
}
