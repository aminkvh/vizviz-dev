//! A movie as an animated GIF: frames are downscaled by area averaging in
//! linear light, quantised to a palette of their own, and written as
//! changed-region deltas. A worker thread does the work so the export loop
//! only ever holds a few frames.

use std::borrow::Cow;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::sync::mpsc::{sync_channel, SyncSender};
use std::thread::JoinHandle;

use gif::{DisposalMethod, Encoder, EncodingError, Frame, Repeat};

use crate::gif_quant::{build_palette, map_to_palette, Pixels, Rgb};

/// GIF's delay unit is 1/100 s and browsers play 0 and 1 as 10, so 2 (50
/// fps) is the fastest reliable rate.
pub const MAX_FPS: f64 = 50.0;
pub const DEFAULT_WIDTH: u32 = 800;
/// Frames waiting for the worker; the export blocks past this.
const QUEUE: usize = 3;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct GifSettings {
    /// Widest the GIF may be; frames are never enlarged.
    pub width: u32,
    pub dither: bool,
}

impl Default for GifSettings {
    fn default() -> Self {
        Self {
            width: DEFAULT_WIDTH,
            dither: true,
        }
    }
}

/// The GIF's size for `src` frames at no more than `width`, same aspect.
pub fn output_size((sw, sh): (u32, u32), width: u32) -> (u16, u16) {
    let w = width.clamp(1, sw.max(1)).min(u32::from(u16::MAX));
    let h = (u64::from(sh) * u64::from(w) + u64::from(sw) / 2) / u64::from(sw.max(1));
    (w as u16, h.clamp(1, u64::from(u16::MAX)) as u16)
}

/// Which movie frames become GIF frames, and how long each is shown.
/// Movie rates above `MAX_FPS` drop frames evenly; the rest keep them all.
/// Delays are whole hundredths that add up to the true elapsed time, so
/// 30 fps plays as 3, 4, 3, 3, 4, 3 ... rather than a steady, too-fast 3.
pub struct Timing {
    movie_fps: f64,
    gif_fps: f64,
    seen: u64,
    shown: u64,
}

impl Timing {
    pub fn new(fps: f32) -> Self {
        let movie_fps = f64::from(fps).max(1.0);
        Self {
            movie_fps,
            gif_fps: movie_fps.min(MAX_FPS),
            seen: 0,
            shown: 0,
        }
    }

    /// The delay for the next movie frame, in hundredths of a second, or
    /// `None` when that frame is dropped.
    pub fn next_delay(&mut self) -> Option<u16> {
        let ratio = self.gif_fps / self.movie_fps;
        let kept = ((self.seen + 1) as f64 * ratio).floor() > (self.seen as f64 * ratio).floor();
        self.seen += 1;
        if !kept {
            return None;
        }
        let at = |k: u64| (k as f64 * 100.0 / self.gif_fps).round() as u64;
        let delay = at(self.shown + 1) - at(self.shown);
        self.shown += 1;
        Some(delay as u16)
    }
}

fn srgb_to_linear(v: u8) -> f32 {
    let c = f32::from(v) / 255.0;
    if c <= 0.04045 {
        c / 12.92
    } else {
        ((c + 0.055) / 1.055).powf(2.4)
    }
}

fn linear_to_srgb(v: f32) -> u8 {
    let c = v.clamp(0.0, 1.0);
    let s = if c <= 0.003_130_8 {
        c * 12.92
    } else {
        1.055 * c.powf(1.0 / 2.4) - 0.055
    };
    (s * 255.0).round() as u8
}

/// For each of `dst` output cells, the `(source index, weight)` pairs
/// covering it; weights sum to 1.
fn coverage(src: usize, dst: usize) -> Vec<Vec<(usize, f32)>> {
    let scale = src as f64 / dst as f64;
    (0..dst)
        .map(|i| {
            let (from, to) = (i as f64 * scale, (i + 1) as f64 * scale);
            (from.floor() as usize..(to.ceil() as usize).min(src))
                .map(|j| {
                    let overlap = to.min((j + 1) as f64) - from.max(j as f64);
                    (j, (overlap / scale) as f32)
                })
                .collect()
        })
        .collect()
}

/// RGBA8 `src` scaled to `dst` (both `(width, height)`) as RGB8, each
/// output pixel the area-weighted mean of the pixels under it.
pub fn downscale(rgba: &[u8], src: (usize, usize), dst: (usize, usize)) -> Vec<u8> {
    let linear: Vec<f32> = (0..256).map(|v| srgb_to_linear(v as u8)).collect();
    let across = coverage(src.0, dst.0);
    let mut rows = vec![0.0f32; dst.0 * src.1 * 3];
    for y in 0..src.1 {
        for (x, taps) in across.iter().enumerate() {
            for k in 0..3 {
                rows[(y * dst.0 + x) * 3 + k] = taps
                    .iter()
                    .map(|&(j, w)| w * linear[usize::from(rgba[(y * src.0 + j) * 4 + k])])
                    .sum();
            }
        }
    }
    let down = coverage(src.1, dst.1);
    let mut out = vec![0u8; dst.0 * dst.1 * 3];
    for (y, taps) in down.iter().enumerate() {
        for i in 0..dst.0 * 3 {
            let v: f32 = taps.iter().map(|&(j, w)| w * rows[j * dst.0 * 3 + i]).sum();
            out[y * dst.0 * 3 + i] = linear_to_srgb(v);
        }
    }
    out
}

struct Rect {
    left: usize,
    top: usize,
    width: usize,
    height: usize,
}

/// The smallest rectangle holding every pixel that differs, if any.
fn changed_rect(before: &[u8], after: &[u8], width: usize) -> Option<Rect> {
    let (mut min_x, mut min_y, mut max_x, mut max_y) = (usize::MAX, usize::MAX, 0, 0);
    for (i, (a, b)) in before
        .chunks_exact(3)
        .zip(after.chunks_exact(3))
        .enumerate()
    {
        if a != b {
            let (x, y) = (i % width, i / width);
            (min_x, max_x) = (min_x.min(x), max_x.max(x));
            (min_y, max_y) = (min_y.min(y), max_y.max(y));
        }
    }
    (min_x != usize::MAX).then(|| Rect {
        left: min_x,
        top: min_y,
        width: max_x - min_x + 1,
        height: max_y - min_y + 1,
    })
}

/// `rect`'s pixels of `rgb`, and whether each differs from `before`.
fn crop(rgb: &[u8], before: Option<&[u8]>, width: usize, rect: &Rect) -> (Vec<u8>, Vec<bool>) {
    let (mut pixels, mut changed) = (Vec::new(), Vec::new());
    for y in rect.top..rect.top + rect.height {
        for x in rect.left..rect.left + rect.width {
            let at = (y * width + x) * 3;
            pixels.extend_from_slice(&rgb[at..at + 3]);
            changed.push(before.is_none_or(|b| b[at..at + 3] != rgb[at..at + 3]));
        }
    }
    (pixels, changed)
}

fn flat(palette: &[Rgb]) -> Vec<u8> {
    palette.iter().flatten().copied().collect()
}

/// Writes frames as they arrive.
pub struct GifWriter<W: Write> {
    encoder: Encoder<W>,
    source: (usize, usize),
    size: (usize, usize),
    dither: bool,
    timing: Timing,
    previous: Option<Vec<u8>>,
}

impl<W: Write> GifWriter<W> {
    pub fn new(
        out: W,
        source: (u32, u32),
        fps: f32,
        settings: GifSettings,
    ) -> Result<Self, EncodingError> {
        let (w, h) = output_size(source, settings.width);
        let mut encoder = Encoder::new(out, w, h, &[])?;
        encoder.set_repeat(Repeat::Infinite)?;
        Ok(Self {
            encoder,
            source: (source.0 as usize, source.1 as usize),
            size: (usize::from(w), usize::from(h)),
            dither: settings.dither,
            timing: Timing::new(fps),
            previous: None,
        })
    }

    /// Adds the next movie frame (RGBA8, `source` size).
    pub fn push(&mut self, rgba: &[u8]) -> Result<(), EncodingError> {
        let Some(delay) = self.timing.next_delay() else {
            return Ok(());
        };
        let rgb = downscale(rgba, self.source, self.size);
        let frame = match self.changes(&rgb) {
            Some((rect, before)) => self.quantized(&rgb, before, &rect, delay),
            None => still_frame(delay),
        };
        self.encoder.write_frame(&frame)?;
        self.previous = Some(rgb);
        Ok(())
    }

    /// What to draw for `rgb`: the whole frame the first time, then the
    /// changed rectangle against the last frame (`None` if nothing changed).
    fn changes(&self, rgb: &[u8]) -> Option<(Rect, Option<&[u8]>)> {
        let Some(before) = &self.previous else {
            let all = Rect {
                left: 0,
                top: 0,
                width: self.size.0,
                height: self.size.1,
            };
            return Some((all, None));
        };
        changed_rect(before, rgb, self.size.0).map(|rect| (rect, Some(before.as_slice())))
    }

    /// `rect` of `rgb` on its own palette; with `before`, pixels that did
    /// not change are left transparent so the last frame shows through.
    fn quantized(
        &self,
        rgb: &[u8],
        before: Option<&[u8]>,
        rect: &Rect,
        delay: u16,
    ) -> Frame<'static> {
        let (pixels, changed) = crop(rgb, before, self.size.0, rect);
        let transparent = changed.contains(&false);
        let slots = if transparent { 255 } else { 256 };
        let mut palette = build_palette(&pixels, &changed, slots);
        let clear = palette.len() as u8;
        if transparent {
            palette.push([0, 0, 0]);
        }
        let all = Pixels {
            rgb: &pixels,
            active: &changed,
            width: rect.width,
        };
        let indices = map_to_palette(&all, &palette, clear, self.dither);
        Frame {
            delay,
            dispose: DisposalMethod::Keep,
            transparent: transparent.then_some(clear),
            left: rect.left as u16,
            top: rect.top as u16,
            width: rect.width as u16,
            height: rect.height as u16,
            palette: Some(flat(&palette)),
            buffer: Cow::Owned(indices),
            ..Frame::default()
        }
    }

    pub fn finish(self) -> Result<W, EncodingError> {
        self.encoder.into_inner()
    }
}

/// A frame identical to the last: one transparent pixel that only adds time.
fn still_frame(delay: u16) -> Frame<'static> {
    Frame {
        delay,
        dispose: DisposalMethod::Keep,
        transparent: Some(0),
        width: 1,
        height: 1,
        palette: Some(vec![0, 0, 0]),
        buffer: Cow::Owned(vec![0]),
        ..Frame::default()
    }
}

/// A `GifWriter` on a worker thread, fed through a small queue.
pub struct GifJob {
    frames: Option<SyncSender<Vec<u8>>>,
    worker: Option<JoinHandle<Result<(), String>>>,
    /// How the worker ended, once it has.
    result: Option<Result<(), String>>,
    path: PathBuf,
}

impl GifJob {
    pub fn start(
        path: PathBuf,
        source: (u32, u32),
        fps: f32,
        settings: GifSettings,
    ) -> Result<Self, String> {
        let file =
            File::create(&path).map_err(|e| format!("can't create {}: {e}", path.display()))?;
        let writer = GifWriter::new(BufWriter::new(file), source, fps, settings)
            .map_err(|e| e.to_string())?;
        let (frames, queue) = sync_channel::<Vec<u8>>(QUEUE);
        let worker = std::thread::spawn(move || run(writer, queue));
        Ok(Self {
            frames: Some(frames),
            worker: Some(worker),
            result: None,
            path,
        })
    }

    /// Queues a frame; blocks while the worker is `QUEUE` frames behind.
    pub fn push(&mut self, rgba: Vec<u8>) -> Result<(), String> {
        let sent = self.frames.as_ref().is_some_and(|tx| tx.send(rgba).is_ok());
        if sent {
            return Ok(());
        }
        Err(self
            .join()
            .err()
            .unwrap_or_else(|| "the GIF writer stopped".into()))
    }

    /// No more frames are coming.
    pub fn close(&mut self) {
        self.frames = None;
    }

    /// `Some` once the worker has finished, with how it went.
    pub fn poll(&mut self) -> Option<Result<(), String>> {
        if self.worker.as_ref().is_some_and(JoinHandle::is_finished) {
            let _ = self.join();
        }
        self.result.clone()
    }

    /// Waits for the worker, which ends once the queue is closed.
    fn join(&mut self) -> Result<(), String> {
        self.frames = None;
        let Some(worker) = self.worker.take() else {
            return self.result.clone().unwrap_or(Ok(()));
        };
        let result = worker
            .join()
            .unwrap_or_else(|_| Err("the GIF writer crashed".into()));
        self.result = Some(result.clone());
        result
    }

    /// Stops the worker and deletes the file unless it was already complete.
    pub fn abort(mut self) {
        let complete = self.poll() == Some(Ok(()));
        let _ = self.join();
        if !complete {
            let _ = std::fs::remove_file(&self.path);
        }
    }
}

fn run(
    mut writer: GifWriter<BufWriter<File>>,
    queue: std::sync::mpsc::Receiver<Vec<u8>>,
) -> Result<(), String> {
    for rgba in queue {
        writer.push(&rgba).map_err(|e| e.to_string())?;
    }
    let mut out = writer.finish().map_err(|e| e.to_string())?;
    out.flush().map_err(|e| e.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn delays(fps: f32, frames: usize) -> Vec<u16> {
        let mut timing = Timing::new(fps);
        (0..frames).filter_map(|_| timing.next_delay()).collect()
    }

    fn solid(w: usize, h: usize, rgb: [u8; 3]) -> Vec<u8> {
        (0..w * h)
            .flat_map(|_| [rgb[0], rgb[1], rgb[2], 255])
            .collect()
    }

    fn painted(w: usize, h: usize, t: usize) -> Vec<u8> {
        let mut px = solid(w, h, [243, 245, 249]);
        for y in 8..24 {
            for x in 8 + t..24 + t {
                let at = (y * w + x) * 4;
                px[at..at + 3].copy_from_slice(&[(x * 5) as u8, (y * 9) as u8, 90]);
            }
        }
        px
    }

    fn decode(bytes: &[u8]) -> (gif::Decoder<&[u8]>, Vec<gif::Frame<'static>>) {
        let mut opts = gif::DecodeOptions::new();
        opts.set_color_output(gif::ColorOutput::Indexed);
        let mut decoder = opts.read_info(bytes).unwrap();
        let mut frames = Vec::new();
        while let Some(frame) = decoder.read_next_frame().unwrap() {
            frames.push(frame.clone());
        }
        (decoder, frames)
    }

    #[test]
    fn thirty_fps_alternates_three_and_four_hundredths() {
        assert_eq!(delays(30.0, 6), [3, 4, 3, 3, 4, 3]);
    }

    #[test]
    fn delays_add_up_to_the_movie_length_and_never_drop_below_two() {
        for fps in [10.0, 24.0, 25.0, 29.97, 30.0, 50.0] {
            let d = delays(fps, 300);
            let total: u32 = d.iter().map(|&x| u32::from(x)).sum();
            let exact = 300.0 * 100.0 / f64::from(fps);
            assert!(
                (f64::from(total) - exact).abs() <= 1.0,
                "{fps}: {total} vs {exact}"
            );
            assert!(d.iter().all(|&x| x >= 2), "{fps}");
        }
    }

    #[test]
    fn rates_above_fifty_drop_frames_evenly() {
        let d = delays(60.0, 60);
        assert_eq!(d.len(), 50);
        assert_eq!(d.iter().map(|&x| u32::from(x)).sum::<u32>(), 100);
    }

    #[test]
    fn size_keeps_the_aspect_and_never_enlarges() {
        assert_eq!(output_size((1280, 720), 800), (800, 450));
        assert_eq!(output_size((640, 480), 800), (640, 480));
        assert_eq!(output_size((1280, 720), 1), (1, 1));
    }

    #[test]
    fn downscaling_a_flat_colour_keeps_it_exactly() {
        let out = downscale(&solid(37, 23, [243, 245, 249]), (37, 23), (11, 7));
        assert!(out.chunks_exact(3).all(|p| p == [243, 245, 249]));
    }

    #[test]
    fn downscaling_by_two_averages_in_linear_light() {
        let mut px = solid(2, 2, [0, 0, 0]);
        px[0..3].copy_from_slice(&[255, 255, 255]);
        px[4..7].copy_from_slice(&[255, 255, 255]);
        let out = downscale(&px, (2, 2), (1, 1));
        assert_eq!(
            out[0], 188,
            "half white in linear light is sRGB 188, not 128"
        );
    }

    #[test]
    fn the_file_decodes_to_the_frames_written_and_loops_forever() {
        let (w, h) = (64, 40);
        let settings = GifSettings {
            width: 32,
            dither: true,
        };
        let mut writer = GifWriter::new(Vec::new(), (w as u32, h as u32), 30.0, settings).unwrap();
        for t in 0..5 {
            writer.push(&painted(w, h, t)).unwrap();
        }
        let bytes = writer.finish().unwrap();
        let (decoder, frames) = decode(&bytes);
        assert_eq!(frames.len(), 5);
        assert_eq!((decoder.width(), decoder.height()), (32, 20));
        assert_eq!(decoder.repeat(), Repeat::Infinite);
        assert_eq!(frames.iter().map(|f| u32::from(f.delay)).sum::<u32>(), 17);
    }

    #[test]
    fn later_frames_only_cover_what_moved() {
        let (w, h) = (64, 40);
        let settings = GifSettings {
            width: 64,
            dither: true,
        };
        let mut writer = GifWriter::new(Vec::new(), (w as u32, h as u32), 30.0, settings).unwrap();
        writer.push(&painted(w, h, 0)).unwrap();
        writer.push(&painted(w, h, 2)).unwrap();
        writer.push(&painted(w, h, 2)).unwrap();
        let (_, frames) = decode(&writer.finish().unwrap());
        assert_eq!((frames[0].width, frames[0].height), (64, 40));
        assert!(frames[1].width < 64 && frames[1].left >= 8);
        assert_eq!((frames[2].width, frames[2].height), (1, 1));
    }
}
