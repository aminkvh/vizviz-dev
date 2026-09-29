//! The Movie as `State` runs it: following the playhead in the live view,
//! and exporting frame by frame. An export sets the camera and frames for
//! one movie frame, lets the normal `render_scene` upload them, then reads
//! the frame back on the throwaway renderer (`render_offscreen`) -- so the
//! GPU always holds this frame's coordinates when it is captured.

use std::path::PathBuf;
use std::process::Child;

use super::*;
use crate::movie::{apply_at, Frames};
use crate::movie_export;
use crate::movie_gif::GifJob;
use crate::movie_state::{frame_counts, ExportAsk};

/// Frames a streamed trajectory reads ahead of the one being shown.
const READ_AHEAD: usize = 8;

/// A running export, and what the view was before it took over.
pub(super) struct ExportJob {
    dir: PathBuf,
    total: usize,
    next: usize,
    ssaa: bool,
    png: bool,
    mp4: bool,
    /// The GIF being written on its worker thread.
    gif: Option<GifJob>,
    /// Every frame is rendered; only ffmpeg and the GIF worker are left.
    finishing: bool,
    /// The camera at time 0.
    base: Camera,
    saved_camera: Camera,
    saved_frames: Vec<(StructureId, usize)>,
    /// This frame's camera once its structure frames are set (and any
    /// streamed frame has arrived); the capture re-applies it, so a drag
    /// during the UI pass cannot move the frame.
    posed: Option<Camera>,
    /// ffmpeg, once every PNG is written.
    encoder: Option<Child>,
}

/// How one of an export's slow tails (ffmpeg, the GIF worker) stands.
enum Tail {
    Running,
    Done,
    Failed(String),
}

impl ExportJob {
    fn ffmpeg_tail(&mut self) -> Tail {
        let Some(child) = self.encoder.as_mut() else {
            return Tail::Done;
        };
        match child.try_wait() {
            Ok(Some(status)) if status.success() => Tail::Done,
            Ok(Some(_)) => Tail::Failed("Wrote the frames; ffmpeg failed to encode them".into()),
            Ok(None) => Tail::Running,
            Err(e) => Tail::Failed(format!("Wrote the frames; ffmpeg: {e}")),
        }
    }

    fn gif_tail(&mut self) -> Tail {
        let Some(gif) = self.gif.as_mut() else {
            return Tail::Done;
        };
        match gif.poll() {
            None => Tail::Running,
            Some(Ok(())) => Tail::Done,
            Some(Err(e)) => Tail::Failed(format!("GIF export failed: {e}")),
        }
    }

    /// Still running, all done, or the first failure.
    fn tails(&mut self) -> Tail {
        match (self.ffmpeg_tail(), self.gif_tail()) {
            (Tail::Failed(e), _) | (_, Tail::Failed(e)) => Tail::Failed(e),
            (Tail::Done, Tail::Done) => Tail::Done,
            _ => Tail::Running,
        }
    }

    fn summary(&self) -> String {
        let mut parts = Vec::new();
        if self.png {
            parts.push(format!("{} frames", self.total));
        }
        if self.mp4 {
            parts.push(movie_export::MP4_NAME.to_string());
        }
        if self.gif.is_some() {
            parts.push(movie_export::GIF_NAME.to_string());
        }
        format!("Wrote {} to {}", parts.join(", "), self.dir.display())
    }

    /// Stops what still runs and removes an unfinished GIF.
    fn abandon(self) {
        if let Some(mut child) = self.encoder {
            let _ = child.kill();
        }
        if let Some(gif) = self.gif {
            gif.abort();
        }
    }
}

impl State {
    pub(super) fn playing(&self) -> bool {
        self.timeline.playing || self.movie.playing || self.movie_export.is_some()
    }

    /// Start of every redraw: follows the playhead, or poses the next
    /// export frame.
    pub(super) fn step_movie(&mut self) {
        if self.movie_export.is_some() {
            self.pose_export_frame();
            return;
        }
        self.movie.adopt_manual_camera(&self.camera);
        self.movie.tick(Instant::now());
        if self.movie.active() {
            self.show_movie_time();
        }
    }

    fn show_movie_time(&mut self) {
        if self.movie.playing {
            self.timeline.playing = false;
        }
        let frames = self
            .movie
            .apply(&mut self.camera, &frame_counts(&self.scene));
        if !self.show_frames(&frames) {
            self.movie.touch();
        }
    }

    /// Puts each structure on its frame. False while a streamed frame is
    /// still being read, so the caller retries.
    fn show_frames(&mut self, frames: &Frames) -> bool {
        let mut ready = true;
        for (&raw, &frame) in frames {
            let id = StructureId::from_raw(raw);
            let Some(loaded) = self.scene.structure(id) else {
                continue;
            };
            let waker = Some(self.waker.clone());
            self.prefetch
                .ahead(id, &loaded.structure, frame, READ_AHEAD, false, waker);
            let (present, changed) = (loaded.structure.has_frame(frame), loaded.frame != frame);
            match (present, changed) {
                (false, _) => ready = false,
                (true, true) => self.scene.set_frame_live(id, frame),
                (true, false) => {}
            }
        }
        ready
    }

    /// After the UI pass: a `movie export` or the panel's Export button.
    pub(super) fn take_movie_export_request(&mut self) {
        let Some(ask) = self.movie.export_request.take() else {
            return;
        };
        if let Err(error) = self.start_movie_export(ask) {
            self.end_export(error, true);
        }
    }

    fn start_movie_export(&mut self, ask: ExportAsk) -> Result<(), String> {
        if self.movie_export.is_some() {
            return Err("an export is already running".into());
        }
        if self.movie.movie.is_empty() {
            return Err("the movie has no clips yet".into());
        }
        if self.scene.structures().next().is_none() {
            return Err("no structure loaded".into());
        }
        if !ask.png && ask.gif.is_none() {
            return Err("no output chosen: switch on PNG frames or GIF".into());
        }
        std::fs::create_dir_all(&ask.dir)
            .map_err(|e| format!("can't create {}: {e}", ask.dir.display()))?;
        let gif = self.start_gif(&ask)?;
        self.movie_export = Some(ExportJob {
            total: self.movie.movie.frame_total(),
            next: 0,
            ssaa: ask.ssaa,
            png: ask.png,
            mp4: ask.mp4 && ask.png && movie_export::ffmpeg_available(),
            gif,
            finishing: false,
            base: self.movie.base_for(&self.camera),
            saved_camera: self.camera.clone(),
            saved_frames: self
                .scene
                .structures()
                .map(|(id, s)| (id, s.frame))
                .collect(),
            posed: None,
            encoder: None,
            dir: ask.dir,
        });
        self.movie.playing = false;
        self.movie.exporting = true;
        self.timeline.playing = false;
        Ok(())
    }

    fn start_gif(&self, ask: &ExportAsk) -> Result<Option<GifJob>, String> {
        let movie = &self.movie.movie;
        ask.gif
            .map(|settings| {
                let path = ask.dir.join(movie_export::GIF_NAME);
                GifJob::start(path, (movie.width, movie.height), movie.fps, settings)
            })
            .transpose()
    }

    fn pose_export_frame(&mut self) {
        let Some(job) = &self.movie_export else {
            return;
        };
        let finishing = job.finishing;
        if self.movie.cancel_export {
            if !finishing {
                self.restore_after_export();
            }
            self.end_export("Export cancelled".into(), false);
            return;
        }
        if finishing {
            return;
        }
        let t = self.movie.movie.time_of_frame(job.next);
        let (camera, frames) =
            apply_at(&self.movie.movie, t, &job.base, &frame_counts(&self.scene));
        self.camera = camera.clone();
        let ready = self.show_frames(&frames);
        if let Some(job) = &mut self.movie_export {
            job.posed = ready.then_some(camera);
        }
    }

    /// After the UI pass (and so after `render_scene`): writes the posed
    /// frame, or checks on ffmpeg and the GIF once all are rendered.
    pub(super) fn capture_movie_frame(&mut self) {
        let Some(job) = &self.movie_export else {
            return;
        };
        if job.finishing {
            self.poll_tails();
        } else if job.posed.is_some() && !self.gpu_cache.building() {
            self.write_export_frame();
        }
    }

    fn write_export_frame(&mut self) {
        let Some(job) = &self.movie_export else {
            return;
        };
        let Some(posed) = job.posed.clone() else {
            return;
        };
        self.camera = posed;
        let (index, total, ssaa) = (job.next, job.total, job.ssaa);
        let (w, h) = (self.movie.movie.width, self.movie.movie.height);
        let written = match self.render_offscreen(w, h, ssaa) {
            Some(rgba) => self.emit_frame(rgba),
            None => Err("nothing to draw".into()),
        };
        if let Err(e) = written {
            self.restore_after_export();
            self.end_export(format!("Export failed at frame {}: {e}", index + 1), true);
            return;
        }
        self.movie.status = Some(format!("Exporting frame {}/{total}", index + 1));
        if let Some(job) = &mut self.movie_export {
            job.next += 1;
            job.posed = None;
        }
        if index + 1 == total {
            self.frames_done();
        }
    }

    /// Writes this frame's PNG and hands it to the GIF worker, as asked.
    fn emit_frame(&mut self, rgba: Vec<u8>) -> Result<(), String> {
        let (w, h) = (self.movie.movie.width, self.movie.movie.height);
        let Some(job) = &mut self.movie_export else {
            return Ok(());
        };
        if job.png {
            let path = movie_export::frame_path(&job.dir, job.next, job.total);
            write_screenshot(&path, &rgba, w, h).map_err(|e| e.to_string())?;
        }
        job.gif.as_mut().map_or(Ok(()), |gif| gif.push(rgba))
    }

    /// Every frame is rendered: give the view back, then let ffmpeg and
    /// the GIF worker finish.
    fn frames_done(&mut self) {
        self.restore_after_export();
        let fps = self.movie.movie.fps;
        let Some(job) = &mut self.movie_export else {
            return;
        };
        job.finishing = true;
        if let Some(gif) = &mut job.gif {
            gif.close();
        }
        if job.mp4 {
            match movie_export::start_encoding(&job.dir, fps, job.total) {
                Ok(child) => job.encoder = Some(child),
                Err(e) => {
                    let line = format!("Wrote the frames; ffmpeg didn't start: {e}");
                    self.end_export(line, true);
                    return;
                }
            }
        }
        self.movie.status = Some("Encoding".into());
        self.poll_tails();
    }

    fn poll_tails(&mut self) {
        let Some(job) = &mut self.movie_export else {
            return;
        };
        match job.tails() {
            Tail::Running => {}
            Tail::Done => {
                let line = job.summary();
                self.end_export(line, false);
            }
            Tail::Failed(line) => self.end_export(line, true),
        }
    }

    fn restore_after_export(&mut self) {
        let Some(job) = &self.movie_export else {
            return;
        };
        self.camera = job.saved_camera.clone();
        for (id, frame) in job.saved_frames.clone() {
            self.scene.set_frame_live(id, frame);
        }
    }

    fn end_export(&mut self, line: String, failed: bool) {
        if let Some(job) = self.movie_export.take() {
            job.abandon();
        }
        self.movie.exporting = false;
        self.movie.cancel_export = false;
        self.log.push(line.clone());
        self.notice = Some(if failed {
            Notice::error(line.clone())
        } else {
            Notice::info(line.clone())
        });
        self.movie.status = Some(line);
    }
}
