//! `gifinfo FILE.gif [FRAMES_DIR]`: what an exported GIF holds and how close
//! its composited frames are to the source PNGs (area-averaged to the GIF
//! size in sRGB space, an implementation independent of the app's).

use std::path::Path;

use crate::metrics::psnr;
use crate::yuv::Rgb;

pub fn run(path: &Path, sources: Option<Vec<Rgb>>) {
    let mut options = gif::DecodeOptions::new();
    options.set_color_output(gif::ColorOutput::RGBA);
    let mut decoder = options
        .read_info(std::fs::File::open(path).unwrap())
        .unwrap();
    let (w, h) = (decoder.width() as usize, decoder.height() as usize);
    println!("canvas {w}x{h}, repeat {:?}", decoder.repeat());
    let mut canvas = vec![0u8; w * h * 3];
    let (mut n, mut partial, mut ticks) = (0usize, 0usize, 0u32);
    let mut scores = Vec::new();
    while let Some(frame) = decoder.read_next_frame().unwrap() {
        partial += ((frame.width as usize, frame.height as usize) != (w, h)) as usize;
        paint(&mut canvas, w, frame);
        ticks += u32::from(frame.delay);
        if let Some(src) = sources.as_ref().and_then(|s| s.get(n)) {
            scores.push(vs_source(&canvas, src, w, h));
        }
        n += 1;
    }
    println!("frames {n}, partial-rect frames {partial}, total delay {ticks} cs");
    if !scores.is_empty() {
        let mean = scores.iter().sum::<f64>() / scores.len() as f64;
        let worst = scores.iter().cloned().fold(f64::MAX, f64::min);
        println!("RGB PSNR vs source: mean {mean:.2} dB, worst frame {worst:.2} dB");
    }
}

fn paint(canvas: &mut [u8], w: usize, frame: &gif::Frame) {
    for (i, px) in frame.buffer.chunks_exact(4).enumerate() {
        if px[3] == 0 {
            continue;
        }
        let x = frame.left as usize + i % frame.width as usize;
        let y = frame.top as usize + i / frame.width as usize;
        canvas[(y * w + x) * 3..(y * w + x) * 3 + 3].copy_from_slice(&px[..3]);
    }
}

fn vs_source(canvas: &[u8], src: &Rgb, w: usize, h: usize) -> f64 {
    let (sx, sy) = (src.w as f64 / w as f64, src.h as f64 / h as f64);
    let mut sq = 0.0;
    for y in 0..h {
        for x in 0..w {
            for k in 0..3 {
                let mean = area_mean(
                    src,
                    (x as f64 * sx, (x + 1) as f64 * sx),
                    (y as f64 * sy, (y + 1) as f64 * sy),
                    k,
                );
                sq += (mean - canvas[(y * w + x) * 3 + k] as f64).powi(2);
            }
        }
    }
    psnr(sq / (w * h * 3) as f64)
}

fn area_mean(src: &Rgb, (x0, x1): (f64, f64), (y0, y1): (f64, f64), k: usize) -> f64 {
    let (mut acc, mut total) = (0.0, 0.0);
    for j in y0.floor() as usize..(y1.ceil() as usize).min(src.h) {
        let wy = y1.min((j + 1) as f64) - y0.max(j as f64);
        for i in x0.floor() as usize..(x1.ceil() as usize).min(src.w) {
            let wx = x1.min((i + 1) as f64) - x0.max(i as f64);
            acc += wx * wy * src.px[(j * src.w + i) * 3 + k] as f64;
            total += wx * wy;
        }
    }
    acc / total
}
