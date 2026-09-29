//! openh264 encode, MP4 mux with the `mp4` crate, and decode-back through
//! the same MP4 (which checks the mux as well as the codec).

use std::fs::File;
use std::io::BufReader;
use std::path::Path;
use std::time::Instant;

use mp4::{
    AvcConfig, FourCC, MediaConfig, Mp4Config, Mp4Reader, Mp4Sample, Mp4Writer, TrackConfig,
    TrackType,
};
use openh264::decoder::Decoder;
use openh264::encoder::{
    BitRate, Complexity, Encoder, EncoderConfig, FrameRate, FrameType, Profile, RateControlMode,
    UsageType, VuiConfig,
};
use openh264::formats::{YUVBuffer, YUVSource};
use openh264::OpenH264API;

use crate::yuv::Yuv;

pub struct Settings {
    pub bitrate_bps: u32,
    pub fps: f32,
    pub usage: UsageType,
    pub mode: RateControlMode,
    pub skip: bool,
}

pub struct Encoded {
    pub encode_secs: f64,
    pub mux_secs: f64,
    pub bytes: u64,
}

#[cfg(feature = "libloading")]
fn api() -> OpenH264API {
    let path = std::env::var("OPENH264_DLL").expect("set OPENH264_DLL to Cisco's library");
    OpenH264API::from_blob_path(path).expect("Cisco library rejected")
}

#[cfg(not(feature = "libloading"))]
fn api() -> OpenH264API {
    OpenH264API::from_source()
}

pub fn encode_to_mp4(frames: &[Yuv], s: &Settings, out: &Path) -> Encoded {
    let config = EncoderConfig::new()
        .bitrate(BitRate::from_bps(s.bitrate_bps))
        .max_frame_rate(FrameRate::from_hz(s.fps))
        .rate_control_mode(s.mode)
        .usage_type(s.usage)
        .profile(Profile::High)
        .complexity(Complexity::High)
        .skip_frames(s.skip)
        .vui(VuiConfig::bt709());
    let mut encoder = Encoder::with_api_config(api(), config).expect("encoder init");
    let (mut stream, mut sync, mut skipped) = (Vec::new(), Vec::new(), 0usize);
    let t = Instant::now();
    for f in frames {
        let mut planes = Vec::with_capacity(f.y.len() * 3 / 2);
        planes.extend_from_slice(&f.y);
        planes.extend_from_slice(&f.u);
        planes.extend_from_slice(&f.v);
        let buf = YUVBuffer::from_vec(planes, f.w, f.h);
        let bits = encoder.encode(&buf).expect("encode");
        skipped += matches!(bits.frame_type(), FrameType::Skip) as usize;
        sync.push(matches!(bits.frame_type(), FrameType::IDR | FrameType::I));
        let mut packet = Vec::new();
        bits.write_vec(&mut packet);
        stream.push(packet);
    }
    let encode_secs = t.elapsed().as_secs_f64();
    eprintln!("skipped {skipped} of {}", frames.len());
    let t = Instant::now();
    mux(&stream, &sync, frames[0].w, frames[0].h, s.fps, out);
    Encoded {
        encode_secs,
        mux_secs: t.elapsed().as_secs_f64(),
        bytes: std::fs::metadata(out).unwrap().len(),
    }
}

/// NAL units of an Annex B buffer, without start codes.
fn nals(data: &[u8]) -> Vec<&[u8]> {
    let mut starts = Vec::new();
    let mut i = 0;
    while i + 3 <= data.len() {
        if data[i..i + 3] == [0, 0, 1] {
            starts.push(i + 3);
            i += 3;
        } else {
            i += 1;
        }
    }
    (0..starts.len())
        .map(|k| {
            let end = starts.get(k + 1).map_or(data.len(), |&s| s - 3);
            let mut nal = &data[starts[k]..end];
            while nal.last() == Some(&0) {
                nal = &nal[..nal.len() - 1];
            }
            nal
        })
        .collect()
}

fn mux(packets: &[Vec<u8>], sync: &[bool], w: usize, h: usize, fps: f32, out: &Path) {
    let (mut sps, mut pps) = (Vec::new(), Vec::new());
    let mut samples: Vec<(Vec<u8>, bool, u32)> = Vec::new();
    for packet in packets {
        let mut avcc = Vec::new();
        for nal in nals(packet) {
            match nal[0] & 0x1f {
                7 => sps = nal.to_vec(),
                8 => pps = nal.to_vec(),
                1 | 5 => {
                    avcc.extend_from_slice(&(nal.len() as u32).to_be_bytes());
                    avcc.extend_from_slice(nal);
                }
                _ => {}
            }
        }
        if avcc.is_empty() {
            if let Some(last) = samples.last_mut() {
                last.2 += 1;
            }
        } else {
            samples.push((avcc, false, 1u32));
        }
    }
    for (sample, is_sync) in samples.iter_mut().zip(sync_of_nonempty(packets, sync)) {
        sample.1 = is_sync;
    }
    let timescale = 90_000u32;
    let duration = (timescale as f32 / fps).round() as u32;
    let mut writer = Mp4Writer::write_start(
        std::io::BufWriter::new(File::create(out).unwrap()),
        &Mp4Config {
            major_brand: FourCC::from(*b"isom"),
            minor_version: 512,
            compatible_brands: vec![FourCC::from(*b"isom"), FourCC::from(*b"avc1")],
            timescale,
        },
    )
    .unwrap();
    writer
        .add_track(&TrackConfig {
            track_type: TrackType::Video,
            timescale,
            language: "und".into(),
            media_conf: MediaConfig::AvcConfig(AvcConfig {
                width: w as u16,
                height: h as u16,
                seq_param_set: sps,
                pic_param_set: pps,
            }),
        })
        .unwrap();
    let mut start = 0u64;
    for (bytes, is_sync, repeat) in samples {
        let ticks = duration * repeat;
        writer
            .write_sample(
                1,
                &Mp4Sample {
                    start_time: start,
                    duration: ticks,
                    rendering_offset: 0,
                    is_sync,
                    bytes: bytes.into(),
                },
            )
            .unwrap();
        start += ticks as u64;
    }
    writer.write_end().unwrap();
}

fn sync_of_nonempty(packets: &[Vec<u8>], sync: &[bool]) -> Vec<bool> {
    packets
        .iter()
        .zip(sync)
        .filter(|(p, _)| !p.is_empty())
        .map(|(_, s)| *s)
        .collect()
}

fn annexb(avcc: &[u8]) -> Vec<u8> {
    let mut out = Vec::new();
    let mut i = 0;
    while i + 4 <= avcc.len() {
        let len = u32::from_be_bytes(avcc[i..i + 4].try_into().unwrap()) as usize;
        out.extend_from_slice(&[0, 0, 0, 1]);
        out.extend_from_slice(&avcc[i + 4..i + 4 + len]);
        i += 4 + len;
    }
    out
}

/// Reads the MP4 back and decodes every sample with openh264's decoder.
pub fn decode_mp4(path: &Path, fps: f32) -> Vec<Yuv> {
    let size = std::fs::metadata(path).unwrap().len();
    let mut mp4 = Mp4Reader::read_header(BufReader::new(File::open(path).unwrap()), size)
        .expect("mp4 header");
    let (id, track) = mp4.tracks().iter().next().map(|(k, v)| (*k, v)).unwrap();
    let mut decoder = Decoder::with_api_config(api(), Default::default()).expect("decoder");
    let mut header = Vec::new();
    for ps in [
        track.sequence_parameter_set(),
        track.picture_parameter_set(),
    ] {
        header.extend_from_slice(&[0, 0, 0, 1]);
        header.extend_from_slice(ps.unwrap());
    }
    decoder.decode(&header).expect("parameter sets");
    let unit = (90_000.0 / fps).round() as u32;
    let mut out = Vec::new();
    for n in 1..=mp4.sample_count(id).unwrap() {
        let sample = mp4.read_sample(id, n).unwrap().unwrap();
        if let Some(yuv) = decoder.decode(&annexb(&sample.bytes)).expect("decode") {
            let frame = copy_planes(&yuv);
            let repeat = (sample.duration / unit).max(1) as usize;
            for _ in 1..repeat {
                out.push(copy_planes(&yuv));
            }
            out.push(frame);
        }
    }
    for yuv in decoder.flush_remaining().unwrap() {
        out.push(copy_planes(&yuv));
    }
    out
}

fn copy_planes(d: &impl YUVSource) -> Yuv {
    let (w, h) = d.dimensions();
    let (sy, su, sv) = d.strides();
    let rows = |plane: &[u8], stride: usize, pw: usize, ph: usize| -> Vec<u8> {
        (0..ph)
            .flat_map(|r| plane[r * stride..r * stride + pw].iter().copied())
            .collect()
    };
    Yuv {
        w,
        h,
        y: rows(d.y(), sy, w, h),
        u: rows(d.u(), su, w / 2, h / 2),
        v: rows(d.v(), sv, w / 2, h / 2),
    }
}
