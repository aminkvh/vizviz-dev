//! RGB <-> planar YUV 4:2:0, BT.709 limited range (Y 16..235, CbCr 16..240).
//!
//! Down: Y per pixel; chroma from the mean RGB of each 2x2 block.
//! Up: chroma sampled bilinearly, sample centres at the middle of each 2x2
//! block (JPEG/MPEG-1 siting; H.264's default is left-cosited, so a real
//! player may differ by half a chroma sample horizontally).

pub struct Rgb {
    pub w: usize,
    pub h: usize,
    pub px: Vec<u8>,
}

pub struct Yuv {
    pub w: usize,
    pub h: usize,
    pub y: Vec<u8>,
    pub u: Vec<u8>,
    pub v: Vec<u8>,
}

fn clamp8(v: f32) -> u8 {
    v.round().clamp(0.0, 255.0) as u8
}

fn luma(r: f32, g: f32, b: f32) -> f32 {
    16.0 + 219.0 * (0.2126 * r + 0.7152 * g + 0.0722 * b) / 255.0
}

fn chroma(r: f32, g: f32, b: f32) -> (f32, f32) {
    let y = 0.2126 * r + 0.7152 * g + 0.0722 * b;
    let cb = 128.0 + 224.0 * (b - y) / (255.0 * 1.8556);
    let cr = 128.0 + 224.0 * (r - y) / (255.0 * 1.5748);
    (cb, cr)
}

impl Yuv {
    pub fn from_rgb(rgb: &Rgb) -> Yuv {
        let (w, h) = (rgb.w, rgb.h);
        let (cw, ch) = (w / 2, h / 2);
        let mut out = Yuv {
            w,
            h,
            y: vec![0; w * h],
            u: vec![0; cw * ch],
            v: vec![0; cw * ch],
        };
        for (i, p) in rgb.px.chunks_exact(3).enumerate() {
            out.y[i] = clamp8(luma(p[0] as f32, p[1] as f32, p[2] as f32));
        }
        for cy in 0..ch {
            for cx in 0..cw {
                let mut acc = [0.0f32; 3];
                for (dx, dy) in [(0, 0), (1, 0), (0, 1), (1, 1)] {
                    let o = ((cy * 2 + dy) * w + cx * 2 + dx) * 3;
                    for k in 0..3 {
                        acc[k] += rgb.px[o + k] as f32 * 0.25;
                    }
                }
                let (cb, cr) = chroma(acc[0], acc[1], acc[2]);
                out.u[cy * cw + cx] = clamp8(cb);
                out.v[cy * cw + cx] = clamp8(cr);
            }
        }
        out
    }

    pub fn to_rgb(&self) -> Rgb {
        let (w, h) = (self.w, self.h);
        let (cw, ch) = (w / 2, h / 2);
        let mut px = vec![0u8; w * h * 3];
        for y in 0..h {
            let (y0, y1, fy) = taps(y, ch);
            for x in 0..w {
                let (x0, x1, fx) = taps(x, cw);
                let cb = bilinear(&self.u, cw, (x0, x1, fx), (y0, y1, fy));
                let cr = bilinear(&self.v, cw, (x0, x1, fx), (y0, y1, fy));
                let yy = (self.y[y * w + x] as f32 - 16.0) * 255.0 / 219.0;
                let (pb, pr) = ((cb - 128.0) * 255.0 / 224.0, (cr - 128.0) * 255.0 / 224.0);
                let r = yy + 1.5748 * pr;
                let b = yy + 1.8556 * pb;
                let g = (yy - 0.2126 * r - 0.0722 * b) / 0.7152;
                let o = (y * w + x) * 3;
                px[o] = clamp8(r);
                px[o + 1] = clamp8(g);
                px[o + 2] = clamp8(b);
            }
        }
        Rgb { w, h, px }
    }
}

/// The two chroma samples around full-resolution position `i` and the
/// weight of the second.
fn taps(i: usize, n: usize) -> (usize, usize, f32) {
    let pos = (i as f32 - 0.5) / 2.0;
    let lo = pos.floor().max(0.0);
    let frac = (pos - lo).clamp(0.0, 1.0);
    let lo = lo as usize;
    (lo.min(n - 1), (lo + 1).min(n - 1), frac)
}

fn bilinear(
    plane: &[u8],
    stride: usize,
    (x0, x1, fx): (usize, usize, f32),
    (y0, y1, fy): (usize, usize, f32),
) -> f32 {
    let at = |x: usize, y: usize| plane[y * stride + x] as f32;
    let top = at(x0, y0) * (1.0 - fx) + at(x1, y0) * fx;
    let bottom = at(x0, y1) * (1.0 - fx) + at(x1, y1) * fx;
    top * (1.0 - fy) + bottom * fy
}
