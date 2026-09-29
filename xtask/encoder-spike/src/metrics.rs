//! Quality metrics against the source frames.
//!
//! PSNR = 10 log10(255^2 / MSE), MSE over every sample of the frame
//! (Y: luma plane; RGB: all three channels), capped at 100 dB; the table
//! reports the mean of per-frame PSNR.
//!
//! SSIM: standard formula with C1 = (0.01*255)^2, C2 = (0.03*255)^2 on
//! 8x8 uniform windows placed every 4 pixels, averaged over the frame
//! (RGB: mean of the three channels' SSIM).
//!
//! Edge PSNR: same MSE but only over "edge pixels" of the source: a pixel
//! whose RGB differs by more than 32 in any channel from its right or lower
//! neighbour, plus those neighbours (thin coloured lines and their sides).

use crate::yuv::{Rgb, Yuv};

#[derive(Default, Clone, Copy)]
pub struct Score {
    pub psnr_y: f64,
    pub psnr_rgb: f64,
    pub ssim_y: f64,
    pub ssim_rgb: f64,
    pub psnr_edge: f64,
}

pub fn psnr(mse: f64) -> f64 {
    if mse <= 0.0 {
        100.0
    } else {
        (10.0 * (255.0f64 * 255.0 / mse).log10()).min(100.0)
    }
}

fn mse(a: &[u8], b: &[u8]) -> f64 {
    let sum: u64 = a
        .iter()
        .zip(b)
        .map(|(&x, &y)| (x as i32 - y as i32).pow(2) as u64)
        .sum();
    sum as f64 / a.len() as f64
}

pub fn edge_mask(src: &Rgb) -> Vec<bool> {
    let (w, h) = (src.w, src.h);
    let differs =
        |i: usize, j: usize| (0..3).any(|k| src.px[i * 3 + k].abs_diff(src.px[j * 3 + k]) > 32);
    let mut mask = vec![false; w * h];
    for y in 0..h {
        for x in 0..w {
            let i = y * w + x;
            for j in [(x + 1 < w).then_some(i + 1), (y + 1 < h).then_some(i + w)]
                .into_iter()
                .flatten()
            {
                if differs(i, j) {
                    mask[i] = true;
                    mask[j] = true;
                }
            }
        }
    }
    mask
}

fn masked_mse(a: &Rgb, b: &Rgb, mask: &[bool]) -> Option<f64> {
    let (mut sum, mut n) = (0u64, 0u64);
    for (i, &m) in mask.iter().enumerate().filter(|(_, m)| **m) {
        let _ = m;
        for k in 0..3 {
            sum += (a.px[i * 3 + k] as i32 - b.px[i * 3 + k] as i32).pow(2) as u64;
            n += 1;
        }
    }
    (n > 0).then(|| sum as f64 / n as f64)
}

fn ssim_plane(a: &[u8], b: &[u8], w: usize, h: usize, stride3: usize, off: usize) -> f64 {
    const C1: f64 = (0.01 * 255.0) * (0.01 * 255.0);
    const C2: f64 = (0.03 * 255.0) * (0.03 * 255.0);
    let at = |p: &[u8], x: usize, y: usize| p[(y * w + x) * stride3 + off] as f64;
    let (mut total, mut count) = (0.0, 0u64);
    for wy in (0..=h - 8).step_by(4) {
        for wx in (0..=w - 8).step_by(4) {
            let (mut sa, mut sb, mut saa, mut sbb, mut sab) = (0.0, 0.0, 0.0, 0.0, 0.0);
            for y in wy..wy + 8 {
                for x in wx..wx + 8 {
                    let (p, q) = (at(a, x, y), at(b, x, y));
                    sa += p;
                    sb += q;
                    saa += p * p;
                    sbb += q * q;
                    sab += p * q;
                }
            }
            let n = 64.0;
            let (ma, mb) = (sa / n, sb / n);
            let (va, vb, cov) = (saa / n - ma * ma, sbb / n - mb * mb, sab / n - ma * mb);
            total += ((2.0 * ma * mb + C1) * (2.0 * cov + C2))
                / ((ma * ma + mb * mb + C1) * (va + vb + C2));
            count += 1;
        }
    }
    total / count as f64
}

/// `src_yuv` is `Yuv::from_rgb(src)`: the luma the encoder was given.
pub fn score(src: &Rgb, src_yuv: &Yuv, mask: &[bool], out: &Yuv) -> Score {
    let back = out.to_rgb();
    let (w, h) = (src.w, src.h);
    let ssim_rgb = (0..3)
        .map(|k| ssim_plane(&src.px, &back.px, w, h, 3, k))
        .sum::<f64>()
        / 3.0;
    Score {
        psnr_y: psnr(mse(&src_yuv.y, &out.y)),
        psnr_rgb: psnr(mse(&src.px, &back.px)),
        ssim_y: ssim_plane(&src_yuv.y, &out.y, w, h, 1, 0),
        ssim_rgb,
        psnr_edge: masked_mse(src, &back, mask).map_or(100.0, psnr),
    }
}

pub fn mean(scores: &[Score]) -> Score {
    let n = scores.len() as f64;
    let sum = |f: fn(&Score) -> f64| scores.iter().map(f).sum::<f64>() / n;
    Score {
        psnr_y: sum(|s| s.psnr_y),
        psnr_rgb: sum(|s| s.psnr_rgb),
        ssim_y: sum(|s| s.ssim_y),
        ssim_rgb: sum(|s| s.ssim_rgb),
        psnr_edge: sum(|s| s.psnr_edge),
    }
}
