//! The movie model: tracks of clips on one shared time axis, and
//! `apply_at`, the single function that says what the scene looks like at
//! any time. Plain data and math, no GPU or UI, so it saves in sessions
//! and every rule here is unit-tested.
//!
//! Camera clips compose in track order onto the camera at time 0 (a
//! finished clip stays applied); trajectory clips pick a frame per
//! structure. Times are seconds, angles degrees, frames zero-based.

use std::collections::BTreeMap;

use glam::{Quat, Vec3};
use serde::{Deserialize, Serialize};
use vv_render::Camera;

/// Shortest clip, so a clip is always grabbable and never zero-length.
pub const MIN_CLIP: f32 = 0.05;
/// Tolerance for "touching" clips and a clip ending exactly at the end.
const EPS: f32 = 1e-4;

#[derive(Clone, Copy, PartialEq, Eq, Debug, Serialize, Deserialize)]
pub enum Easing {
    Linear,
    EaseInOut,
}

impl Easing {
    pub const ALL: [Easing; 2] = [Easing::Linear, Easing::EaseInOut];

    pub fn name(self) -> &'static str {
        match self {
            Easing::Linear => "linear",
            Easing::EaseInOut => "ease",
        }
    }

    pub fn parse(word: &str) -> Option<Easing> {
        Self::ALL.into_iter().find(|e| e.name() == word)
    }

    /// Progress 0..=1 to eased progress 0..=1.
    pub fn apply(self, p: f32) -> f32 {
        let p = p.clamp(0.0, 1.0);
        match self {
            Easing::Linear => p,
            Easing::EaseInOut => p * p * (3.0 - 2.0 * p),
        }
    }
}

/// A screen axis, as `Camera::rotate_scene` takes it.
#[derive(Clone, Copy, PartialEq, Eq, Debug, Serialize, Deserialize)]
pub enum Axis {
    X,
    Y,
    Z,
}

impl Axis {
    pub const ALL: [Axis; 3] = [Axis::X, Axis::Y, Axis::Z];

    pub fn name(self) -> &'static str {
        match self {
            Axis::X => "x",
            Axis::Y => "y",
            Axis::Z => "z",
        }
    }

    pub fn parse(word: &str) -> Option<Axis> {
        Self::ALL.into_iter().find(|a| a.name() == word)
    }

    fn vector(self) -> Vec3 {
        match self {
            Axis::X => Vec3::X,
            Axis::Y => Vec3::Y,
            Axis::Z => Vec3::Z,
        }
    }
}

/// Where a camera looks from: what a "move to view" clip flies to.
#[derive(Clone, Copy, PartialEq, Debug, Serialize, Deserialize)]
pub struct ViewPose {
    pub target: [f32; 3],
    pub distance: f32,
    pub orientation: [f32; 4],
}

impl ViewPose {
    pub fn of(camera: &Camera) -> Self {
        Self {
            target: camera.target.to_array(),
            distance: camera.distance,
            orientation: camera.orientation.to_array(),
        }
    }

    fn is_finite(&self) -> bool {
        self.target
            .iter()
            .chain(&self.orientation)
            .all(|v| v.is_finite())
            && self.distance.is_finite()
            && self.distance > 0.0
    }

    /// `camera` moved fraction `p` of the way here; exactly here at 1.
    fn blend_into(&self, camera: &mut Camera, p: f32) {
        let (target, orientation) = (
            Vec3::from_array(self.target),
            Quat::from_array(self.orientation).normalize(),
        );
        if p >= 1.0 {
            camera.target = target;
            camera.distance = self.distance;
            camera.orientation = orientation;
            return;
        }
        camera.target = camera.target.lerp(target, p);
        camera.distance += (self.distance - camera.distance) * p;
        camera.orientation = camera.orientation.slerp(orientation, p).normalize();
    }
}

/// What a track holds, and so what its clips can do.
#[derive(Clone, Copy, PartialEq, Eq, Debug, Serialize, Deserialize)]
pub enum TrackKind {
    Rotate,
    Zoom,
    MoveTo,
    Play,
}

impl TrackKind {
    pub const ALL: [TrackKind; 4] = [
        TrackKind::Rotate,
        TrackKind::Zoom,
        TrackKind::MoveTo,
        TrackKind::Play,
    ];

    /// The word `movie addtrack` takes.
    pub fn name(self) -> &'static str {
        match self {
            TrackKind::Rotate => "rotate",
            TrackKind::Zoom => "zoom",
            TrackKind::MoveTo => "moveto",
            TrackKind::Play => "play",
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            TrackKind::Rotate => "Rotate",
            TrackKind::Zoom => "Zoom",
            TrackKind::MoveTo => "Move to view",
            TrackKind::Play => "Play trajectory",
        }
    }

    pub fn parse(word: &str) -> Option<TrackKind> {
        Self::ALL.into_iter().find(|k| k.name() == word)
    }
}

/// One clip's effect. `Play` frames are zero-based and inclusive.
#[derive(Clone, Copy, PartialEq, Debug, Serialize, Deserialize)]
pub enum Effect {
    /// Turns the scene `degrees` about a screen axis over the clip.
    Rotate { axis: Axis, degrees: f32 },
    /// Multiplies the camera distance by a factor running `from` to `to`
    /// (1 = unchanged, above 1 = away), geometrically.
    Zoom { from: f32, to: f32 },
    /// Flies the camera to a pose, replacing whatever earlier tracks did.
    MoveTo(ViewPose),
    /// Steps `structure` through `first..=last` at `fps` frames a second.
    Play {
        structure: u32,
        first: usize,
        last: usize,
        fps: f32,
        looping: bool,
    },
}

impl Effect {
    pub fn kind(&self) -> TrackKind {
        match self {
            Effect::Rotate { .. } => TrackKind::Rotate,
            Effect::Zoom { .. } => TrackKind::Zoom,
            Effect::MoveTo(_) => TrackKind::MoveTo,
            Effect::Play { .. } => TrackKind::Play,
        }
    }

    /// A few words for a clip bar or a listing.
    pub fn label(&self) -> String {
        match self {
            Effect::Rotate { axis, degrees } => {
                format!("{} {degrees:+}\u{b0}", axis.name().to_uppercase())
            }
            Effect::Zoom { from, to } => format!("\u{d7}{from} to \u{d7}{to}"),
            Effect::MoveTo(_) => "to view".into(),
            Effect::Play {
                structure,
                first,
                last,
                ..
            } => format!("#{structure} frames {first}-{last}"),
        }
    }

    fn validate(&self) -> Result<(), String> {
        let ok = match self {
            Effect::Rotate { degrees, .. } => degrees.is_finite(),
            Effect::Zoom { from, to } => [from, to].iter().all(|v| v.is_finite() && **v > 0.0),
            Effect::MoveTo(pose) => pose.is_finite(),
            Effect::Play {
                first, last, fps, ..
            } => first <= last && fps.is_finite() && *fps > 0.0,
        };
        ok.then_some(())
            .ok_or_else(|| format!("a {} clip's settings are out of range", self.kind().name()))
    }
}

#[derive(Clone, Copy, PartialEq, Debug, Serialize, Deserialize)]
pub struct Clip {
    pub start: f32,
    pub length: f32,
    pub easing: Easing,
    pub effect: Effect,
}

impl Clip {
    pub fn end(&self) -> f32 {
        self.start + self.length
    }

    /// Eased progress at `t`: 0 before the clip, 1 from its end on.
    fn progress(&self, t: f32) -> f32 {
        if t >= self.end() {
            1.0
        } else {
            self.easing.apply((t - self.start) / self.length)
        }
    }
}

/// The part of a clip bar being dragged.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Grab {
    Body,
    Left,
    Right,
}

/// Clips of one kind, sorted by start and never overlapping.
#[derive(Clone, PartialEq, Debug, Serialize, Deserialize)]
pub struct Track {
    pub kind: TrackKind,
    pub clips: Vec<Clip>,
}

#[derive(Clone, PartialEq, Debug, Serialize, Deserialize)]
#[serde(default)]
pub struct Movie {
    pub fps: f32,
    pub duration: f32,
    pub width: u32,
    pub height: u32,
    pub tracks: Vec<Track>,
}

impl Default for Movie {
    fn default() -> Self {
        Self {
            fps: 30.0,
            duration: 10.0,
            width: 1280,
            height: 720,
            tracks: Vec::new(),
        }
    }
}

/// Frames a movie is rendered as, and where each sits in time.
impl Movie {
    pub fn frame_total(&self) -> usize {
        ((self.duration * self.fps).ceil() as usize).max(1)
    }

    pub fn time_of_frame(&self, index: usize) -> f32 {
        index as f32 / self.fps
    }

    pub fn is_empty(&self) -> bool {
        self.tracks.iter().all(|t| t.clips.is_empty())
    }

    /// Renames every structure a play clip names (a session stores them by
    /// position, which loading maps to new ids); an unmapped one becomes
    /// `u32::MAX`, which no structure has, so its clip does nothing.
    pub fn remap_structures(&mut self, map: impl Fn(u32) -> Option<u32>) {
        let clips = self.tracks.iter_mut().flat_map(|t| &mut t.clips);
        for clip in clips {
            if let Effect::Play { structure, .. } = &mut clip.effect {
                *structure = map(*structure).unwrap_or(u32::MAX);
            }
        }
    }

    pub fn set_fps(&mut self, fps: f32) -> Result<(), String> {
        if !(1.0..=120.0).contains(&fps) {
            return Err("fps must be between 1 and 120".into());
        }
        self.fps = fps;
        Ok(())
    }

    /// Never shorter than the last clip's end.
    pub fn set_duration(&mut self, seconds: f32) -> Result<(), String> {
        let last = self
            .tracks
            .iter()
            .flat_map(|t| &t.clips)
            .map(Clip::end)
            .fold(MIN_CLIP, f32::max);
        if !seconds.is_finite() || seconds < last - EPS {
            return Err(format!(
                "duration can't be under {last:.2} s (a clip ends there)"
            ));
        }
        self.duration = seconds.min(3600.0);
        Ok(())
    }

    /// Even sizes only, which video encoders need.
    pub fn set_size(&mut self, width: u32, height: u32) {
        let even = |v: u32| (v.clamp(16, 7680)) & !1;
        self.width = even(width);
        self.height = even(height);
    }
}

/// Editing tracks and clips.
impl Movie {
    pub fn add_track(&mut self, kind: TrackKind) -> usize {
        self.tracks.push(Track {
            kind,
            clips: Vec::new(),
        });
        self.tracks.len() - 1
    }

    pub fn remove_track(&mut self, track: usize) -> bool {
        (track < self.tracks.len())
            .then(|| self.tracks.remove(track))
            .is_some()
    }

    pub fn add_clip(&mut self, track: usize, clip: Clip) -> Result<usize, String> {
        self.check(track, None, &clip)?;
        Ok(insert_sorted(&mut self.tracks[track].clips, clip))
    }

    /// Swaps clip `index` for `clip` (a move, resize or edit); the old one
    /// stays if the new one is refused. Returns the clip's new index.
    pub fn replace_clip(
        &mut self,
        track: usize,
        index: usize,
        clip: Clip,
    ) -> Result<usize, String> {
        self.check(track, Some(index), &clip)?;
        let clips = &mut self.tracks[track].clips;
        clips.remove(index);
        Ok(insert_sorted(clips, clip))
    }

    pub fn remove_clip(&mut self, track: usize, index: usize) -> bool {
        match self.tracks.get_mut(track) {
            Some(t) if index < t.clips.len() => {
                t.clips.remove(index);
                true
            }
            _ => false,
        }
    }

    /// Why `clip` can't sit on `track` (replacing clip `skip`), if it can't.
    fn check(&self, track: usize, skip: Option<usize>, clip: &Clip) -> Result<(), String> {
        let t = self.tracks.get(track).ok_or("no such track")?;
        if clip.effect.kind() != t.kind {
            return Err(format!(
                "a {} track holds {} clips",
                t.kind.name(),
                t.kind.name()
            ));
        }
        clip.effect.validate()?;
        if !(clip.start.is_finite() && clip.length.is_finite()) || clip.start < -EPS {
            return Err("a clip can't start before 0".into());
        }
        if clip.length < MIN_CLIP - EPS {
            return Err(format!("a clip must be at least {MIN_CLIP} s long"));
        }
        if clip.end() > self.duration + EPS {
            return Err(format!("the clip runs past the end ({} s)", self.duration));
        }
        let clash = t.clips.iter().enumerate().any(|(i, c)| {
            Some(i) != skip && clip.start < c.end() - EPS && c.start < clip.end() - EPS
        });
        if clash {
            return Err("clips on one track can't overlap".into());
        }
        Ok(())
    }

    /// The free stretch `(from, to)` around clip `index`: how far it may
    /// move or grow without touching a neighbour or the ends.
    pub fn free_range(&self, track: usize, index: usize) -> (f32, f32) {
        let clips = &self.tracks[track].clips;
        let before = index.checked_sub(1).map_or(0.0, |i| clips[i].end());
        let after = clips.get(index + 1).map_or(self.duration, |c| c.start);
        (before, after)
    }

    /// Drags part of a clip by `dt` seconds, stopping at its neighbours,
    /// the ends and the minimum length. Returns the clip's index (a drag
    /// never re-sorts it).
    pub fn drag_clip(
        &mut self,
        track: usize,
        index: usize,
        grab: Grab,
        dt: f32,
    ) -> Result<usize, String> {
        let (lo, hi) = self.free_range(track, index);
        let mut clip = self.tracks[track].clips[index];
        match grab {
            Grab::Body => clip.start = (clip.start + dt).clamp(lo, hi - clip.length),
            Grab::Left => {
                let end = clip.end();
                clip.start = (clip.start + dt).clamp(lo, end - MIN_CLIP);
                clip.length = end - clip.start;
            }
            Grab::Right => clip.length = (clip.length + dt).clamp(MIN_CLIP, hi - clip.start),
        }
        self.replace_clip(track, index, clip)
    }

    /// Where a clip of up to `want` seconds fits on `track`: the earliest
    /// gap at or after `at`, shortened to the gap. `None` if none is big
    /// enough.
    pub fn place(&self, track: usize, at: f32, want: f32) -> Option<(f32, f32)> {
        let mut start = at.max(0.0);
        for c in &self.tracks.get(track)?.clips {
            if c.end() <= start + EPS {
                continue;
            }
            if c.start - start >= MIN_CLIP {
                break;
            }
            start = start.max(c.end());
        }
        let next = self.tracks[track]
            .clips
            .iter()
            .map(|c| c.start)
            .find(|s| *s > start + EPS)
            .unwrap_or(self.duration);
        let length = want.min(next.min(self.duration) - start);
        (length >= MIN_CLIP).then_some((start, length))
    }
}

fn insert_sorted(clips: &mut Vec<Clip>, clip: Clip) -> usize {
    let at = clips.partition_point(|c| c.start <= clip.start);
    clips.insert(at, clip);
    at
}

/// Frames per structure (raw id), only for structures a clip plays.
pub type Frames = BTreeMap<u32, usize>;

/// The camera and per-structure frame at time `t`. Camera clips apply in
/// track order to `base`, the camera at time 0; `frame_counts` (raw id to
/// frame count) keeps each frame in range and drops structures that are
/// gone.
pub fn apply_at(
    movie: &Movie,
    t: f32,
    base: &Camera,
    frame_counts: &BTreeMap<u32, usize>,
) -> (Camera, Frames) {
    let mut camera = base.clone();
    let mut frames = Frames::new();
    for clip in movie.tracks.iter().flat_map(|track| &track.clips) {
        match clip.effect {
            Effect::Play { .. } => play_frame(clip, t, frame_counts, &mut frames),
            _ if t >= clip.start => apply_camera(&mut camera, clip, clip.progress(t)),
            _ => {}
        }
    }
    (camera, frames)
}

fn apply_camera(camera: &mut Camera, clip: &Clip, p: f32) {
    match clip.effect {
        Effect::Rotate { axis, degrees } => {
            camera.rotate_scene(axis.vector(), degrees.to_radians() * p)
        }
        Effect::Zoom { from, to } => camera.zoom(from.powf(1.0 - p) * to.powf(p)),
        Effect::MoveTo(pose) => pose.blend_into(camera, p),
        Effect::Play { .. } => {}
    }
}

/// Before a play clip starts its structure shows the clip's first frame,
/// unless an earlier clip already set one; after, it keeps the last frame
/// reached.
fn play_frame(clip: &Clip, t: f32, counts: &BTreeMap<u32, usize>, frames: &mut Frames) {
    let Effect::Play {
        structure,
        first,
        last,
        fps,
        looping,
    } = clip.effect
    else {
        return;
    };
    let Some(count) = counts.get(&structure) else {
        return;
    };
    let top = count.saturating_sub(1);
    if t < clip.start {
        frames.entry(structure).or_insert(first.min(top));
        return;
    }
    let steps = ((t.min(clip.end()) - clip.start) * fps).floor() as usize;
    let span = last - first + 1;
    let step = if looping {
        steps % span
    } else {
        steps.min(span - 1)
    };
    frames.insert(structure, (first + step).min(top));
}

#[cfg(test)]
mod tests {
    use super::*;

    fn camera() -> Camera {
        Camera::framing(Vec3::ZERO, 10.0)
    }

    fn rotate(start: f32, length: f32, degrees: f32) -> Clip {
        Clip {
            start,
            length,
            easing: Easing::Linear,
            effect: Effect::Rotate {
                axis: Axis::Y,
                degrees,
            },
        }
    }

    fn zoom(start: f32, length: f32, from: f32, to: f32) -> Clip {
        Clip {
            start,
            length,
            easing: Easing::Linear,
            effect: Effect::Zoom { from, to },
        }
    }

    fn play(start: f32, length: f32, first: usize, last: usize, fps: f32, looping: bool) -> Clip {
        Clip {
            start,
            length,
            easing: Easing::Linear,
            effect: Effect::Play {
                structure: 7,
                first,
                last,
                fps,
                looping,
            },
        }
    }

    fn movie_with(kind: TrackKind, clips: &[Clip]) -> Movie {
        let mut movie = Movie::default();
        let track = movie.add_track(kind);
        for clip in clips {
            movie.add_clip(track, *clip).unwrap();
        }
        movie
    }

    fn counts() -> BTreeMap<u32, usize> {
        BTreeMap::from([(7, 100)])
    }

    fn same_view(a: &Camera, b: &Camera) -> bool {
        a.orientation.angle_between(b.orientation) < 1e-4
            && (a.distance - b.distance).abs() < 1e-3
            && a.target.distance(b.target) < 1e-3
    }

    #[test]
    fn easing_hits_both_ends_and_eases_the_middle() {
        for easing in Easing::ALL {
            assert_eq!(easing.apply(0.0), 0.0);
            assert_eq!(easing.apply(1.0), 1.0);
        }
        assert_eq!(Easing::Linear.apply(0.25), 0.25);
        assert!(Easing::EaseInOut.apply(0.25) < 0.25);
        assert_eq!(Easing::EaseInOut.apply(0.5), 0.5);
    }

    #[test]
    fn overlapping_clips_on_one_track_are_refused_but_touching_ones_fit() {
        let mut movie = movie_with(TrackKind::Rotate, &[rotate(1.0, 2.0, 90.0)]);
        assert!(movie.add_clip(0, rotate(2.0, 2.0, 90.0)).is_err());
        assert!(movie.add_clip(0, rotate(0.0, 1.5, 90.0)).is_err());
        assert!(movie.add_clip(0, rotate(3.0, 1.0, 90.0)).is_ok());
        assert!(movie.add_clip(0, rotate(0.0, 1.0, 90.0)).is_ok());
        let starts: Vec<f32> = movie.tracks[0].clips.iter().map(|c| c.start).collect();
        assert_eq!(starts, [0.0, 1.0, 3.0]);
    }

    #[test]
    fn different_tracks_may_overlap_freely() {
        let mut movie = movie_with(TrackKind::Rotate, &[rotate(0.0, 4.0, 90.0)]);
        let zoom_track = movie.add_track(TrackKind::Zoom);
        assert!(movie.add_clip(zoom_track, zoom(1.0, 2.0, 1.0, 2.0)).is_ok());
    }

    #[test]
    fn a_clip_must_match_its_track_and_fit_the_movie() {
        let mut movie = movie_with(TrackKind::Rotate, &[]);
        assert!(movie.add_clip(0, zoom(0.0, 1.0, 1.0, 2.0)).is_err());
        assert!(movie.add_clip(0, rotate(-1.0, 1.0, 90.0)).is_err());
        assert!(movie.add_clip(0, rotate(0.0, 0.0, 90.0)).is_err());
        assert!(movie.add_clip(0, rotate(9.5, 1.0, 90.0)).is_err());
        assert!(movie.add_clip(0, rotate(9.0, 1.0, f32::NAN)).is_err());
        assert!(movie.add_clip(3, rotate(0.0, 1.0, 90.0)).is_err());
    }

    #[test]
    fn a_refused_replace_keeps_the_old_clip() {
        let mut movie = movie_with(
            TrackKind::Rotate,
            &[rotate(0.0, 1.0, 90.0), rotate(2.0, 1.0, 90.0)],
        );
        assert!(movie.replace_clip(0, 0, rotate(1.5, 1.0, 45.0)).is_err());
        assert_eq!(movie.tracks[0].clips[0], rotate(0.0, 1.0, 90.0));
        let moved = movie.replace_clip(0, 0, rotate(4.0, 1.0, 45.0)).unwrap();
        assert_eq!(moved, 1, "moving past a neighbour re-sorts it");
    }

    #[test]
    fn duration_cannot_cut_off_a_clip() {
        let mut movie = movie_with(TrackKind::Rotate, &[rotate(6.0, 2.0, 90.0)]);
        assert!(movie.set_duration(7.0).is_err());
        assert!(movie.set_duration(8.0).is_ok());
    }

    #[test]
    fn place_finds_the_first_gap_and_shortens_to_fit() {
        let movie = movie_with(
            TrackKind::Rotate,
            &[rotate(0.0, 2.0, 90.0), rotate(3.0, 2.0, 90.0)],
        );
        assert_eq!(movie.place(0, 0.0, 5.0), Some((2.0, 1.0)));
        assert_eq!(movie.place(0, 5.0, 3.0), Some((5.0, 3.0)));
        assert_eq!(movie.place(0, 9.99, 3.0), None);
    }

    #[test]
    fn free_range_is_bounded_by_neighbours_and_the_ends() {
        let movie = movie_with(
            TrackKind::Rotate,
            &[rotate(1.0, 1.0, 90.0), rotate(4.0, 1.0, 90.0)],
        );
        assert_eq!(movie.free_range(0, 0), (0.0, 4.0));
        assert_eq!(movie.free_range(0, 1), (2.0, 10.0));
    }

    #[test]
    fn dragging_a_clip_stops_at_its_neighbours_and_the_ends() {
        let mut movie = movie_with(
            TrackKind::Rotate,
            &[rotate(2.0, 1.0, 90.0), rotate(5.0, 1.0, 90.0)],
        );
        movie.drag_clip(0, 0, Grab::Body, -9.0).unwrap();
        assert_eq!(movie.tracks[0].clips[0].start, 0.0);
        movie.drag_clip(0, 0, Grab::Body, 9.0).unwrap();
        assert_eq!(
            movie.tracks[0].clips[0].end(),
            5.0,
            "stops against the next clip"
        );
        movie.drag_clip(0, 1, Grab::Right, 99.0).unwrap();
        assert_eq!(movie.tracks[0].clips[1].end(), 10.0);
        movie.drag_clip(0, 1, Grab::Left, 99.0).unwrap();
        assert!((movie.tracks[0].clips[1].length - MIN_CLIP).abs() < 1e-4);
        assert!(
            (movie.tracks[0].clips[1].end() - 10.0).abs() < 1e-4,
            "the right edge stays put"
        );
    }

    #[test]
    fn time_zero_is_the_base_camera_when_clips_start_later() {
        let movie = movie_with(TrackKind::Rotate, &[rotate(1.0, 2.0, 90.0)]);
        let (cam, frames) = apply_at(&movie, 0.0, &camera(), &counts());
        assert!(same_view(&cam, &camera()));
        assert!(frames.is_empty());
    }

    #[test]
    fn a_rotation_is_partway_in_the_middle_and_stays_done_afterwards() {
        let movie = movie_with(TrackKind::Rotate, &[rotate(1.0, 2.0, 90.0)]);
        let turned = |degrees: f32| {
            let mut c = camera();
            c.rotate_scene(Vec3::Y, degrees.to_radians());
            c
        };
        let (mid, _) = apply_at(&movie, 2.0, &camera(), &counts());
        let (done, _) = apply_at(&movie, 3.0, &camera(), &counts());
        let (later, _) = apply_at(&movie, 9.0, &camera(), &counts());
        assert!(same_view(&mid, &turned(45.0)));
        assert!(same_view(&done, &turned(90.0)));
        assert!(same_view(&later, &turned(90.0)));
    }

    #[test]
    fn a_zoom_runs_geometrically_from_to() {
        let movie = movie_with(TrackKind::Zoom, &[zoom(0.0, 2.0, 1.0, 4.0)]);
        let (mid, _) = apply_at(&movie, 1.0, &camera(), &counts());
        let (end, _) = apply_at(&movie, 2.0, &camera(), &counts());
        assert!((mid.distance / camera().distance - 2.0).abs() < 1e-3);
        assert!((end.distance / camera().distance - 4.0).abs() < 1e-3);
    }

    #[test]
    fn clips_of_different_tracks_run_at_once() {
        let mut movie = movie_with(TrackKind::Rotate, &[rotate(0.0, 2.0, 90.0)]);
        let z = movie.add_track(TrackKind::Zoom);
        movie.add_clip(z, zoom(0.0, 2.0, 1.0, 2.0)).unwrap();
        let (cam, _) = apply_at(&movie, 2.0, &camera(), &counts());
        let mut expected = camera();
        expected.rotate_scene(Vec3::Y, 90f32.to_radians());
        expected.zoom(2.0);
        assert!(same_view(&cam, &expected));
    }

    #[test]
    fn a_move_to_view_lands_exactly_on_its_pose_and_blends_on_the_way() {
        let mut goal = camera();
        goal.rotate_scene(Vec3::X, 60f32.to_radians());
        goal.zoom(0.5);
        goal.target = Vec3::new(3.0, 0.0, 0.0);
        let clip = Clip {
            start: 0.0,
            length: 2.0,
            easing: Easing::Linear,
            effect: Effect::MoveTo(ViewPose::of(&goal)),
        };
        let movie = movie_with(TrackKind::MoveTo, &[clip]);
        let (end, _) = apply_at(&movie, 2.0, &camera(), &counts());
        let (mid, _) = apply_at(&movie, 1.0, &camera(), &counts());
        assert_eq!(end.target, goal.target);
        assert_eq!(end.distance, goal.distance);
        assert!(end.orientation.angle_between(goal.orientation) < 1e-5);
        assert!((mid.target.x - 1.5).abs() < 1e-4);
    }

    #[test]
    fn track_order_decides_how_effects_compose() {
        let mut goal = camera();
        goal.rotate_scene(Vec3::X, 90f32.to_radians());
        let go = Clip {
            start: 0.0,
            length: 1.0,
            easing: Easing::Linear,
            effect: Effect::MoveTo(ViewPose::of(&goal)),
        };
        let mut spin_last = movie_with(TrackKind::MoveTo, &[go]);
        let t = spin_last.add_track(TrackKind::Rotate);
        spin_last.add_clip(t, rotate(0.0, 1.0, 90.0)).unwrap();
        let mut spin_first = movie_with(TrackKind::Rotate, &[rotate(0.0, 1.0, 90.0)]);
        let t = spin_first.add_track(TrackKind::MoveTo);
        spin_first.add_clip(t, go).unwrap();
        let (a, _) = apply_at(&spin_last, 1.0, &camera(), &counts());
        let (b, _) = apply_at(&spin_first, 1.0, &camera(), &counts());
        assert!(
            !same_view(&a, &b),
            "a move-to replaces earlier rotations only"
        );
        assert!(same_view(&b, &goal));
    }

    #[test]
    fn a_trajectory_clip_steps_frames_at_its_speed() {
        let movie = movie_with(TrackKind::Play, &[play(1.0, 4.0, 10, 19, 5.0, false)]);
        let at = |t: f32| apply_at(&movie, t, &camera(), &counts()).1[&7];
        assert_eq!(at(0.0), 10, "before the clip: its first frame");
        assert_eq!(at(1.0), 10);
        assert_eq!(at(2.0), 15);
        assert_eq!(at(4.9), 19, "clamps at the last frame");
        assert_eq!(at(9.0), 19, "keeps the last frame after the clip");
    }

    #[test]
    fn a_looping_trajectory_wraps_inside_its_range() {
        let movie = movie_with(TrackKind::Play, &[play(0.0, 5.0, 10, 14, 5.0, true)]);
        let at = |t: f32| apply_at(&movie, t, &camera(), &counts()).1[&7];
        assert_eq!(at(0.9), 14);
        assert_eq!(at(1.0), 10);
        assert_eq!(at(1.5), 12);
    }

    #[test]
    fn frames_stay_in_range_and_missing_structures_are_skipped() {
        let movie = movie_with(TrackKind::Play, &[play(0.0, 5.0, 0, 500, 100.0, false)]);
        let (_, frames) = apply_at(&movie, 4.0, &camera(), &counts());
        assert_eq!(frames[&7], 99);
        let (_, none) = apply_at(&movie, 4.0, &camera(), &BTreeMap::new());
        assert!(none.is_empty());
    }

    #[test]
    fn later_tracks_win_a_shared_structure() {
        let mut movie = movie_with(TrackKind::Play, &[play(0.0, 5.0, 0, 9, 1.0, false)]);
        let t = movie.add_track(TrackKind::Play);
        movie
            .add_clip(t, play(2.0, 3.0, 50, 59, 1.0, false))
            .unwrap();
        let at = |t: f32| apply_at(&movie, t, &camera(), &counts()).1[&7];
        assert_eq!(at(1.0), 1);
        assert_eq!(at(3.0), 51);
    }

    #[test]
    fn frames_are_numbered_from_the_duration_and_fps() {
        let mut movie = Movie::default();
        movie.duration = 2.0;
        movie.fps = 24.0;
        assert_eq!(movie.frame_total(), 48);
        assert_eq!(movie.time_of_frame(24), 1.0);
    }

    #[test]
    fn sizes_are_even_and_bounded() {
        let mut movie = Movie::default();
        movie.set_size(1921, 1);
        assert_eq!((movie.width, movie.height), (1920, 16));
    }

    #[test]
    fn a_movie_round_trips_through_json() {
        let mut movie = movie_with(TrackKind::Rotate, &[rotate(1.0, 2.0, 90.0)]);
        let p = movie.add_track(TrackKind::Play);
        movie.add_clip(p, play(0.0, 3.0, 2, 8, 12.0, true)).unwrap();
        let json = serde_json::to_string(&movie).unwrap();
        assert_eq!(serde_json::from_str::<Movie>(&json).unwrap(), movie);
    }

    #[test]
    fn remapping_renames_played_structures_and_orphans_unknown_ones() {
        let mut movie = movie_with(TrackKind::Play, &[play(0.0, 2.0, 0, 4, 5.0, false)]);
        movie.remap_structures(|id| (id == 7).then_some(2));
        assert!(matches!(
            movie.tracks[0].clips[0].effect,
            Effect::Play { structure: 2, .. }
        ));
        movie.remap_structures(|_| None);
        assert!(matches!(
            movie.tracks[0].clips[0].effect,
            Effect::Play {
                structure: u32::MAX,
                ..
            }
        ));
    }

    #[test]
    fn a_partial_saved_movie_fills_in_defaults() {
        let movie: Movie = serde_json::from_str(r#"{"fps": 24.0}"#).unwrap();
        assert_eq!(movie.fps, 24.0);
        assert_eq!(movie.duration, Movie::default().duration);
    }
}
