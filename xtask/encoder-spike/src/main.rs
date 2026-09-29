//! Throwaway harness behind docs/MOVIE_ENCODER_SPIKE.md.
//!
//!   encoder-spike FRAMES_DIR OUT_DIR [floor] [h264] [av1]
//!
//! FRAMES_DIR holds `frame_NNNN.png` from `movie export`. Prints one
//! markdown table row per encoder setting.

mod av1;
mod gifinfo;
mod h264;
mod metrics;
mod yuv;

use std::path::{Path, PathBuf};

use metrics::{edge_mask, mean, score, Score};
use openh264::encoder::{RateControlMode, UsageType};
use yuv::{Rgb, Yuv};

const FPS: u32 = 30;

fn load_frames(dir: &Path) -> Vec<Rgb> {
    let mut paths: Vec<PathBuf> = std::fs::read_dir(dir)
        .unwrap()
        .map(|e| e.unwrap().path())
        .filter(|p| p.extension().is_some_and(|e| e == "png"))
        .collect();
    paths.sort();
    paths.iter().map(|p| load_png(p, true)).collect()
}

fn load_png(path: &Path, require_opaque: bool) -> Rgb {
    let decoder = png::Decoder::new(std::fs::File::open(path).unwrap());
    let mut reader = decoder.read_info().unwrap();
    let mut buf = vec![0; reader.output_buffer_size()];
    let info = reader.next_frame(&mut buf).unwrap();
    assert_eq!(info.color_type, png::ColorType::Rgba, "expected RGBA PNG");
    let px = buf[..info.buffer_size()]
        .chunks_exact(4)
        .flat_map(|p| {
            assert!(!require_opaque || p[3] == 255, "frame is not opaque");
            [p[0], p[1], p[2]]
        })
        .collect();
    Rgb {
        w: info.width as usize,
        h: info.height as usize,
        px,
    }
}

struct Reference {
    rgb: Vec<Rgb>,
    yuv: Vec<Yuv>,
    masks: Vec<Vec<bool>>,
}

fn scores(r: &Reference, decoded: &[Yuv]) -> Score {
    assert_eq!(decoded.len(), r.rgb.len(), "decoded frame count differs");
    let threads = std::thread::available_parallelism().map_or(4, |n| n.get());
    let per: Vec<Score> = std::thread::scope(|s| {
        let handles: Vec<_> = (0..threads)
            .map(|t| {
                s.spawn(move || {
                    (t..decoded.len())
                        .step_by(threads)
                        .map(|i| score(&r.rgb[i], &r.yuv[i], &r.masks[i], &decoded[i]))
                        .collect::<Vec<_>>()
                })
            })
            .collect();
        handles
            .into_iter()
            .flat_map(|h| h.join().unwrap())
            .collect()
    });
    mean(&per)
}

fn row(name: &str, target: &str, fps_enc: f64, bytes: u64, frames: usize, s: Score) {
    let mbps = bytes as f64 * 8.0 / (frames as f64 / FPS as f64) / 1e6;
    println!(
        "| {name} | {target} | {fps_enc:.1} | {:.2} | {mbps:.1} | {:.2} | {:.2} | {:.4} | {:.4} | {:.2} |",
        bytes as f64 / 1e6, s.psnr_y, s.psnr_rgb, s.ssim_y, s.ssim_rgb, s.psnr_edge
    );
}

fn run_floor(r: &Reference) {
    let s = scores(r, &r.yuv);
    row("4:2:0 round trip, no codec", "-", 0.0, 0, r.rgb.len(), s);
}

fn run_h264(r: &Reference, out: &Path) {
    let cases = [
        (
            UsageType::CameraVideoRealTime,
            RateControlMode::Bitrate,
            true,
            "camera bitrate+skip",
        ),
        (
            UsageType::CameraVideoRealTime,
            RateControlMode::Quality,
            false,
            "camera quality",
        ),
        (
            UsageType::ScreenContentRealTime,
            RateControlMode::Bitrate,
            true,
            "screen bitrate+skip",
        ),
    ];
    for (usage, mode, skip, label) in cases {
        for mbps in [2u32, 4, 8, 20] {
            let path = out.join(format!("h264_{label}_{mbps}.mp4").replace([' ', '+'], "_"));
            let settings = h264::Settings {
                bitrate_bps: mbps * 1_000_000,
                fps: FPS as f32,
                usage,
                mode,
                skip,
            };
            let e = h264::encode_to_mp4(&r.yuv, &settings, &path);
            eprintln!("mux {:.3}s", e.mux_secs);
            let decoded = h264::decode_mp4(&path, FPS as f32);
            let s = scores(r, &decoded);
            let name = format!("openh264 {label}");
            row(
                &name,
                &format!("{mbps} Mbps"),
                r.rgb.len() as f64 / e.encode_secs,
                e.bytes,
                r.rgb.len(),
                s,
            );
        }
    }
}

fn run_av1(r: &Reference, out: &Path) {
    for speed in [10u8, 6] {
        for mbps in [2i32, 4, 8, 20] {
            let path = out.join(format!("av1_s{speed}_{mbps}.ivf"));
            let e = av1::encode_to_ivf(&r.yuv, FPS, mbps * 1_000_000, speed, &path);
            let s = scores(r, &e.recon);
            let label = format!("rav1e speed {speed}");
            row(
                &label,
                &format!("{mbps} Mbps"),
                r.rgb.len() as f64 / e.encode_secs,
                e.bytes,
                r.rgb.len(),
                s,
            );
        }
    }
}

/// `compare A.png B.png`: RGB PSNR and mean signed difference per channel.
fn compare(a: &Path, b: &Path) {
    let (x, y) = (load_png(a, false), load_png(b, false));
    let n = x.px.len() as f64;
    let sq: f64 =
        x.px.iter()
            .zip(&y.px)
            .map(|(&p, &q)| (p as f64 - q as f64).powi(2))
            .sum();
    let bias: Vec<f64> = (0..3)
        .map(|k| {
            x.px.iter()
                .zip(&y.px)
                .skip(k)
                .step_by(3)
                .map(|(&p, &q)| q as f64 - p as f64)
                .sum::<f64>()
                / (n / 3.0)
        })
        .collect();
    println!(
        "PSNR {:.2} dB, mean (B-A) per channel {:?}",
        metrics::psnr(sq / n),
        bias
    );
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args[1] == "gifinfo" {
        let sources = args.get(3).map(|d| load_frames(Path::new(d)));
        return gifinfo::run(Path::new(&args[2]), sources);
    }
    if args[1] == "compare" {
        return compare(Path::new(&args[2]), Path::new(&args[3]));
    }
    let (frames_dir, out) = (Path::new(&args[1]), Path::new(&args[2]));
    let which: Vec<&str> = args[3..].iter().map(String::as_str).collect();
    std::fs::create_dir_all(out).unwrap();
    let rgb = load_frames(frames_dir);
    let yuv: Vec<Yuv> = rgb.iter().map(Yuv::from_rgb).collect();
    let masks = rgb.iter().map(edge_mask).collect::<Vec<_>>();
    let edge_share = masks
        .iter()
        .map(|m| m.iter().filter(|b| **b).count())
        .sum::<usize>() as f64
        / (masks.len() * masks[0].len()) as f64;
    eprintln!(
        "{} frames {}x{}, edge pixels {:.2}%",
        rgb.len(),
        rgb[0].w,
        rgb[0].h,
        edge_share * 100.0
    );
    let r = Reference { rgb, yuv, masks };
    println!("| encoder | target | enc fps | MB | actual Mbps | PSNR Y | PSNR RGB | SSIM Y | SSIM RGB | edge PSNR RGB |");
    println!("|---|---|---|---|---|---|---|---|---|---|");
    if which.contains(&"floor") {
        run_floor(&r);
    }
    if which.contains(&"h264") {
        run_h264(&r, out);
    }
    if which.contains(&"av1") {
        run_av1(&r, out);
    }
}
