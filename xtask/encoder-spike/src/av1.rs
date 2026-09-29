//! rav1e encode into an IVF file. rav1e has no decoder, so quality is
//! measured on the encoder's own reconstruction (`Packet::rec`), which a
//! conforming decoder reproduces bit-exactly; the IVF itself is not
//! independently decoded here.

use std::io::Write;
use std::path::Path;
use std::time::Instant;

use rav1e::prelude::*;

use crate::yuv::Yuv;

pub struct Encoded {
    pub encode_secs: f64,
    pub bytes: u64,
    pub recon: Vec<Yuv>,
}

fn config(w: usize, h: usize, fps: u32, bitrate_bps: i32, speed: u8) -> Config {
    let mut enc = EncoderConfig::with_speed_preset(speed);
    enc.width = w;
    enc.height = h;
    enc.bit_depth = 8;
    enc.chroma_sampling = ChromaSampling::Cs420;
    enc.time_base = Rational::new(1, fps as u64);
    enc.bitrate = bitrate_bps;
    enc.min_key_frame_interval = 0;
    enc.max_key_frame_interval = 240;
    enc.pixel_range = PixelRange::Limited;
    enc.color_description = Some(ColorDescription {
        color_primaries: ColorPrimaries::BT709,
        transfer_characteristics: TransferCharacteristics::BT709,
        matrix_coefficients: MatrixCoefficients::BT709,
    });
    Config::new().with_encoder_config(enc)
}

fn plane_bytes(p: &Plane<u8>) -> Vec<u8> {
    let (w, h) = (p.cfg.width, p.cfg.height);
    p.rows_iter()
        .take(h)
        .flat_map(|r| r[..w].to_vec())
        .collect()
}

pub fn encode_to_ivf(frames: &[Yuv], fps: u32, bitrate_bps: i32, speed: u8, out: &Path) -> Encoded {
    let (w, h) = (frames[0].w, frames[0].h);
    let mut ctx: Context<u8> = config(w, h, fps, bitrate_bps, speed).new_context().unwrap();
    let mut packets: Vec<(u64, Vec<u8>)> = Vec::new();
    let mut recon: Vec<Option<Yuv>> = (0..frames.len()).map(|_| None).collect();
    let t = Instant::now();
    let mut collect = |ctx: &mut Context<u8>, packets: &mut Vec<(u64, Vec<u8>)>| loop {
        match ctx.receive_packet() {
            Ok(pkt) => {
                if let Some(rec) = &pkt.rec {
                    recon[pkt.input_frameno as usize] = Some(Yuv {
                        w,
                        h,
                        y: plane_bytes(&rec.planes[0]),
                        u: plane_bytes(&rec.planes[1]),
                        v: plane_bytes(&rec.planes[2]),
                    });
                }
                packets.push((pkt.input_frameno, pkt.data));
            }
            Err(EncoderStatus::Encoded) => {}
            Err(EncoderStatus::NeedMoreData | EncoderStatus::LimitReached) => break,
            Err(e) => panic!("rav1e: {e}"),
        }
    };
    for f in frames {
        let mut frame = ctx.new_frame();
        frame.planes[0].copy_from_raw_u8(&f.y, w, 1);
        frame.planes[1].copy_from_raw_u8(&f.u, w / 2, 1);
        frame.planes[2].copy_from_raw_u8(&f.v, w / 2, 1);
        ctx.send_frame(frame).unwrap();
        collect(&mut ctx, &mut packets);
    }
    ctx.flush();
    collect(&mut ctx, &mut packets);
    let encode_secs = t.elapsed().as_secs_f64();
    write_ivf(out, w, h, fps, &packets);
    Encoded {
        encode_secs,
        bytes: std::fs::metadata(out).unwrap().len(),
        recon: recon
            .into_iter()
            .map(|r| r.expect("no reconstruction returned"))
            .collect(),
    }
}

fn write_ivf(path: &Path, w: usize, h: usize, fps: u32, packets: &[(u64, Vec<u8>)]) {
    let mut f = std::io::BufWriter::new(std::fs::File::create(path).unwrap());
    f.write_all(b"DKIF").unwrap();
    f.write_all(&0u16.to_le_bytes()).unwrap();
    f.write_all(&32u16.to_le_bytes()).unwrap();
    f.write_all(b"AV01").unwrap();
    f.write_all(&(w as u16).to_le_bytes()).unwrap();
    f.write_all(&(h as u16).to_le_bytes()).unwrap();
    f.write_all(&fps.to_le_bytes()).unwrap();
    f.write_all(&1u32.to_le_bytes()).unwrap();
    f.write_all(&(packets.len() as u32).to_le_bytes()).unwrap();
    f.write_all(&0u32.to_le_bytes()).unwrap();
    for (pts, data) in packets {
        f.write_all(&(data.len() as u32).to_le_bytes()).unwrap();
        f.write_all(&pts.to_le_bytes()).unwrap();
        f.write_all(data).unwrap();
    }
}
