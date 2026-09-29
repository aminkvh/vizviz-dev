//! The Movie panel's live state: where the playhead is, whether it plays,
//! and how the movie's camera is kept apart from the camera the user
//! moves by hand. Only `movie` itself is saved in sessions.

use std::collections::BTreeMap;
use std::path::PathBuf;
use std::time::Instant;

use vv_render::Camera;
use vv_scene::Scene;

use crate::movie::{apply_at, Frames, Movie};
use crate::movie_gif::GifSettings;

/// A request to export the movie, taken by `State` after the UI pass.
pub struct ExportAsk {
    pub dir: PathBuf,
    pub ssaa: bool,
    /// The numbered PNG sequence (the MP4 is encoded from it).
    pub png: bool,
    pub mp4: bool,
    pub gif: Option<GifSettings>,
}

pub struct MovieState {
    pub movie: Movie,
    /// The playhead, in seconds.
    pub time: f32,
    pub playing: bool,
    /// The clip being edited: (track, index).
    pub selected: Option<(usize, usize)>,
    pub export_request: Option<ExportAsk>,
    pub cancel_export: bool,
    /// An export is running (set by `State`).
    pub exporting: bool,
    /// Which outputs the export row asks for.
    pub ssaa: bool,
    pub png: bool,
    pub mp4: bool,
    pub gif: bool,
    pub gif_settings: GifSettings,
    /// The last export line: progress, the folder written, or why it failed.
    pub status: Option<String>,
    /// The time changed or a clip was edited, so the view needs redoing.
    dirty: bool,
    last_tick: Option<Instant>,
    /// The camera at time 0, taken from the live camera the first time the
    /// movie is applied after the user last moved it.
    base: Option<Camera>,
    /// The camera the movie last set: any other live camera was moved by
    /// hand, which makes it the new starting view.
    applied: Option<Camera>,
}

impl Default for MovieState {
    fn default() -> Self {
        Self {
            movie: Movie::default(),
            time: 0.0,
            playing: false,
            selected: None,
            export_request: None,
            cancel_export: false,
            exporting: false,
            ssaa: true,
            png: true,
            mp4: true,
            gif: false,
            gif_settings: GifSettings::default(),
            status: None,
            dirty: false,
            last_tick: None,
            base: None,
            applied: None,
        }
    }
}

impl MovieState {
    /// An export into `dir` with the outputs the panel has switched on.
    pub fn export_ask(&self, dir: PathBuf) -> ExportAsk {
        ExportAsk {
            dir,
            ssaa: self.ssaa,
            png: self.png,
            mp4: self.mp4 && self.png,
            gif: self.gif.then_some(self.gif_settings),
        }
    }

    pub fn set_time(&mut self, seconds: f32) {
        self.time = seconds.clamp(0.0, self.movie.duration);
        self.dirty = true;
    }

    /// Something in the movie changed: show it at the current time again.
    pub fn touch(&mut self) {
        self.dirty = true;
    }

    /// True while the view needs a redraw to follow the movie.
    pub fn active(&self) -> bool {
        self.playing || self.dirty
    }

    pub fn toggle_play(&mut self) {
        self.playing = !self.playing;
        self.last_tick = None;
        if self.playing && self.time >= self.movie.duration {
            self.time = 0.0;
        }
        self.dirty = true;
    }

    /// Moves the playhead by the real time since the last call.
    pub fn tick(&mut self, now: Instant) {
        if !self.playing {
            return;
        }
        let dt = self
            .last_tick
            .map_or(0.0, |t| now.duration_since(t).as_secs_f32());
        self.last_tick = Some(now);
        self.time += dt;
        if self.time >= self.movie.duration {
            self.time = self.movie.duration;
            self.playing = false;
        }
        self.dirty = true;
    }

    /// If `live` isn't the camera the movie last set, the user moved it:
    /// that view becomes the movie's start and the playhead goes back to 0.
    pub fn adopt_manual_camera(&mut self, live: &Camera) {
        if self.applied.as_ref().is_some_and(|a| a != live) {
            self.applied = None;
            self.playing = false;
            self.time = 0.0;
        }
    }

    /// The camera time 0 looks through, for `live` if the movie hasn't
    /// taken the camera over yet.
    pub fn base_for(&self, live: &Camera) -> Camera {
        match (&self.applied, &self.base) {
            (Some(_), Some(base)) => base.clone(),
            _ => live.clone(),
        }
    }

    /// Sets `camera` to the movie's view at the playhead and returns the
    /// frames to show, clearing the pending redraw.
    pub fn apply(&mut self, camera: &mut Camera, counts: &BTreeMap<u32, usize>) -> Frames {
        let base = self.base_for(camera);
        let (view, frames) = apply_at(&self.movie, self.time, &base, counts);
        self.base = Some(base);
        *camera = view.clone();
        self.applied = Some(view);
        self.dirty = false;
        frames
    }

    /// Drops the movie's hold on the camera (a new movie, a loaded session).
    pub fn release_camera(&mut self) {
        self.applied = None;
        self.base = None;
        self.playing = false;
    }
}

/// Frame count per structure (raw id), as `apply_at` takes it.
pub fn frame_counts(scene: &Scene) -> BTreeMap<u32, usize> {
    scene
        .structures()
        .map(|(id, s)| (id.to_raw(), s.structure.frame_count()))
        .collect()
}

/// `1:03.20`-style time for a readout.
pub fn clock(seconds: f32) -> String {
    let whole = seconds.max(0.0);
    let minutes = (whole / 60.0).floor() as u32;
    format!("{minutes}:{:05.2}", whole - minutes as f32 * 60.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::movie::{Axis, Clip, Easing, Effect, TrackKind};
    use glam::Vec3;

    fn spinning() -> MovieState {
        let mut state = MovieState::default();
        let t = state.movie.add_track(TrackKind::Rotate);
        let clip = Clip {
            start: 0.0,
            length: 2.0,
            easing: Easing::Linear,
            effect: Effect::Rotate {
                axis: Axis::Y,
                degrees: 90.0,
            },
        };
        state.movie.add_clip(t, clip).unwrap();
        state
    }

    fn camera() -> Camera {
        Camera::framing(Vec3::ZERO, 10.0)
    }

    #[test]
    fn scrubbing_back_to_zero_restores_the_starting_camera() {
        let mut state = spinning();
        let mut cam = camera();
        let counts = BTreeMap::new();
        state.set_time(2.0);
        state.apply(&mut cam, &counts);
        assert_ne!(cam, camera());
        state.set_time(0.0);
        state.apply(&mut cam, &counts);
        assert_eq!(cam.orientation, camera().orientation);
    }

    #[test]
    fn moving_the_camera_by_hand_makes_it_the_new_start_and_rewinds() {
        let mut state = spinning();
        let mut cam = camera();
        let counts = BTreeMap::new();
        state.set_time(1.0);
        state.apply(&mut cam, &counts);
        cam.zoom(2.0);
        state.adopt_manual_camera(&cam);
        assert_eq!(state.time, 0.0);
        state.set_time(0.0);
        let hand_moved = cam.clone();
        state.apply(&mut cam, &counts);
        assert_eq!(cam.distance, hand_moved.distance);
    }

    #[test]
    fn our_own_camera_is_not_mistaken_for_a_hand_move() {
        let mut state = spinning();
        let mut cam = camera();
        state.set_time(1.0);
        state.apply(&mut cam, &BTreeMap::new());
        state.adopt_manual_camera(&cam);
        assert_eq!(state.time, 1.0);
    }

    #[test]
    fn playing_advances_then_stops_at_the_end() {
        let mut state = spinning();
        state.movie.set_duration(2.0).unwrap();
        state.toggle_play();
        let start = Instant::now();
        state.tick(start);
        state.tick(start + std::time::Duration::from_secs_f32(1.0));
        assert!((state.time - 1.0).abs() < 1e-3);
        state.tick(start + std::time::Duration::from_secs_f32(5.0));
        assert_eq!(state.time, 2.0);
        assert!(!state.playing);
    }

    #[test]
    fn playing_from_the_end_starts_over() {
        let mut state = spinning();
        state.set_time(state.movie.duration);
        state.toggle_play();
        assert_eq!(state.time, 0.0);
    }

    #[test]
    fn the_clock_reads_minutes_and_hundredths() {
        assert_eq!(clock(3.2), "0:03.20");
        assert_eq!(clock(63.25), "1:03.25");
    }
}
