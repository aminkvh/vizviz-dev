//! The Movie panel: transport, movie settings, the track timeline
//! (`movie_timeline`), the selected clip's properties and Export.

use egui::{DragValue, Ui};
use egui_phosphor::regular as icon;

use crate::movie::{Axis, Clip, Easing, Effect, TrackKind, ViewPose, MIN_CLIP};
use crate::movie_export;
use crate::movie_state::clock;
use crate::theme::space;
use crate::ui::{timeline_target, AppUi};
use crate::widgets::{self, Variant};

/// A structure a play clip can name.
struct Playable {
    raw: u32,
    label: String,
    frames: usize,
}

impl AppUi<'_> {
    pub(crate) fn movie_ui(&mut self, ui: &mut Ui) {
        egui::ScrollArea::vertical().show(ui, |ui| {
            self.movie_transport(ui);
            self.movie_settings(ui);
            ui.add_space(space::GAP);
            self.movie_timeline(ui);
            self.movie_add_track_menu(ui);
            ui.add_space(space::GAP);
            self.movie_clip_properties(ui);
            ui.add_space(space::GAP);
            self.movie_export_row(ui);
        });
    }

    fn movie_transport(&mut self, ui: &mut Ui) {
        ui.horizontal(|ui| {
            if widgets::icon_button(ui, icon::SKIP_BACK, false)
                .on_hover_text("Back to the start")
                .clicked()
            {
                self.movie.playing = false;
                self.movie.set_time(0.0);
            }
            let playing = self.movie.playing;
            let (glyph, tip) = if playing {
                (icon::PAUSE, "Pause")
            } else {
                (icon::PLAY, "Play")
            };
            if widgets::icon_button(ui, glyph, playing)
                .on_hover_text(tip)
                .clicked()
            {
                self.movie.toggle_play();
            }
            let (time, duration) = (self.movie.time, self.movie.movie.duration);
            ui.label(format!("{} / {}", clock(time), clock(duration)));
        });
    }

    fn movie_settings(&mut self, ui: &mut Ui) {
        let mut error = None;
        ui.horizontal_wrapped(|ui| {
            let movie = &mut self.movie.movie;
            let mut fps = movie.fps;
            if ui
                .add(DragValue::new(&mut fps).range(1.0..=120.0).suffix(" fps"))
                .changed()
            {
                error = movie.set_fps(fps).err();
            }
            let mut duration = movie.duration;
            let drag = DragValue::new(&mut duration)
                .range(1.0..=600.0)
                .speed(0.1)
                .suffix(" s");
            if ui.add(drag).on_hover_text("Length of the movie").changed() {
                error = movie.set_duration(duration).err();
            }
            let (mut width, mut height) = (movie.width, movie.height);
            let w_changed = size_field(ui, &mut width);
            ui.label("\u{d7}");
            let changed = w_changed | size_field(ui, &mut height);
            if changed {
                movie.set_size(width, height);
            }
        });
        if let Some(e) = error {
            self.fail(e);
        }
    }

    fn movie_add_track_menu(&mut self, ui: &mut Ui) {
        let button = widgets::button(ui, icon::PLUS, "Add track", Variant::Secondary);
        let mut picked = None;
        egui::Popup::from_toggle_button_response(&button)
            .close_behavior(egui::PopupCloseBehavior::CloseOnClick)
            .show(|ui| {
                ui.set_min_width(160.0);
                for kind in TrackKind::ALL {
                    if widgets::button(ui, "", kind.label(), Variant::Ghost).clicked() {
                        picked = Some(kind);
                    }
                }
            });
        if let Some(kind) = picked {
            self.movie.movie.add_track(kind);
        }
        if self.movie.movie.tracks.is_empty() {
            widgets::caption(
                ui,
                "Add a track: rotate, zoom, move to a view or play a trajectory. Clips on different tracks run together.",
            );
        }
    }

    /// Adds a clip to `track` at the playhead (or the next free stretch).
    pub(crate) fn add_clip_at_playhead(&mut self, track: usize) {
        match self.new_clip(track) {
            Ok(clip) => match self.movie.movie.add_clip(track, clip) {
                Ok(index) => self.movie.selected = Some((track, index)),
                Err(e) => self.fail(e),
            },
            Err(e) => self.fail(e),
        }
        self.movie.touch();
    }

    fn new_clip(&self, track: usize) -> Result<Clip, String> {
        let kind = self.movie.movie.tracks[track].kind;
        let (start, length) = self
            .movie
            .movie
            .place(track, self.movie.time, 2.0)
            .ok_or("no room for a clip here: move the playhead or shorten a clip")?;
        Ok(Clip {
            start,
            length,
            easing: Easing::EaseInOut,
            effect: self.default_effect(kind)?,
        })
    }

    fn default_effect(&self, kind: TrackKind) -> Result<Effect, String> {
        Ok(match kind {
            TrackKind::Rotate => Effect::Rotate {
                axis: Axis::Y,
                degrees: 360.0,
            },
            TrackKind::Zoom => Effect::Zoom { from: 1.0, to: 0.5 },
            TrackKind::MoveTo => Effect::MoveTo(ViewPose::of(self.camera)),
            TrackKind::Play => {
                let (id, frames) = timeline_target(self.scene)
                    .ok_or("load a trajectory (or an NMR ensemble) to play it")?;
                Effect::Play {
                    structure: id.to_raw(),
                    first: 0,
                    last: frames - 1,
                    fps: self.timeline.fps,
                    looping: false,
                }
            }
        })
    }

    fn playable_structures(&self) -> Vec<Playable> {
        self.scene
            .structures()
            .filter(|(_, s)| s.structure.frame_count() > 1)
            .map(|(id, s)| Playable {
                raw: id.to_raw(),
                label: format!("#{} {}", id.to_raw(), s.label),
                frames: s.structure.frame_count(),
            })
            .collect()
    }

    fn movie_clip_properties(&mut self, ui: &mut Ui) {
        widgets::section(ui, "Clip");
        let Some((track, index)) = self.movie.selected else {
            widgets::caption(
                ui,
                "Click a clip to edit it. Drag its bar to move it, its ends to resize.",
            );
            return;
        };
        let Some(clip) = self
            .movie
            .movie
            .tracks
            .get(track)
            .and_then(|t| t.clips.get(index))
            .copied()
        else {
            self.movie.selected = None;
            return;
        };
        let mut edited = clip;
        let playable = self.playable_structures();
        let pose = ViewPose::of(self.camera);
        let mut remove = false;
        ui.horizontal_wrapped(|ui| {
            let free = self.movie.movie.free_range(track, index);
            timing_controls(ui, &mut edited, free);
            effect_controls(ui, &mut edited.effect, &playable, pose);
            remove = widgets::button(ui, icon::TRASH, "Delete", Variant::Danger).clicked();
        });
        let changed = remove || edited != clip;
        if remove {
            self.movie.movie.remove_clip(track, index);
            self.movie.selected = None;
        } else if edited != clip {
            match self.movie.movie.replace_clip(track, index, edited) {
                Ok(new_index) => self.movie.selected = Some((track, new_index)),
                Err(e) => self.fail(e),
            }
        }
        if changed {
            self.movie.touch();
        }
    }

    fn movie_export_row(&mut self, ui: &mut Ui) {
        widgets::section(ui, "Export");
        if self.movie.exporting {
            self.export_progress(ui);
            return;
        }
        let ffmpeg = movie_export::ffmpeg_available();
        ui.horizontal_wrapped(|ui| {
            self.export_switches(ui, ffmpeg);
            let ready = !self.movie.movie.is_empty() && (self.movie.png || self.movie.gif);
            let export = ui.add_enabled_ui(ready, |ui| {
                widgets::button(ui, icon::EXPORT, "Export\u{2026}", Variant::Primary)
            });
            if export.inner.clicked() {
                self.pick_export_folder();
            }
            export
                .response
                .on_disabled_hover_text("Add a clip and switch on an output first");
        });
        if self.movie.gif {
            self.gif_options(ui);
        }
        widgets::caption(ui, &self.export_note(ffmpeg));
        if let Some(status) = &self.movie.status {
            widgets::caption(ui, status);
        }
    }

    fn export_switches(&mut self, ui: &mut Ui, ffmpeg: bool) {
        let movie = &mut self.movie;
        if widgets::switch(ui, movie.ssaa, "Smooth edges").clicked() {
            movie.ssaa = !movie.ssaa;
        }
        if widgets::switch(ui, movie.png, "PNG frames").clicked() {
            movie.png = !movie.png;
        }
        if ffmpeg && movie.png && widgets::switch(ui, movie.mp4, "Also MP4").clicked() {
            movie.mp4 = !movie.mp4;
        }
        if widgets::switch(ui, movie.gif, "Also GIF").clicked() {
            movie.gif = !movie.gif;
        }
    }

    fn gif_options(&mut self, ui: &mut Ui) {
        let movie = &mut self.movie;
        ui.horizontal_wrapped(|ui| {
            ui.label("GIF width");
            let widest = movie.movie.width.max(16);
            let drag = DragValue::new(&mut movie.gif_settings.width)
                .range(16..=widest)
                .speed(8.0)
                .suffix(" px");
            ui.add(drag)
                .on_hover_text("Frames are shrunk to this width; never enlarged");
            if widgets::switch(ui, movie.gif_settings.dither, "Dither").clicked() {
                movie.gif_settings.dither = !movie.gif_settings.dither;
            }
        });
    }

    fn export_note(&self, ffmpeg: bool) -> String {
        let mut parts = Vec::new();
        if self.movie.png {
            parts.push("numbered PNGs".to_string());
        }
        if ffmpeg && self.movie.png && self.movie.mp4 {
            parts.push(movie_export::MP4_NAME.to_string());
        }
        if self.movie.gif {
            parts.push(movie_export::GIF_NAME.to_string());
        }
        let list = if parts.is_empty() {
            "Nothing".to_string()
        } else {
            parts.join(", ")
        };
        let hint = if ffmpeg || !self.movie.png {
            ""
        } else {
            " Install ffmpeg (on your PATH) to get an MP4 as well."
        };
        format!("Writes {list} into the folder.{hint}")
    }

    fn export_progress(&mut self, ui: &mut Ui) {
        ui.horizontal(|ui| {
            ui.label(
                self.movie
                    .status
                    .clone()
                    .unwrap_or_else(|| "Exporting".into()),
            );
            if widgets::button(ui, "", "Cancel", Variant::Secondary).clicked() {
                self.movie.cancel_export = true;
            }
        });
    }

    fn pick_export_folder(&mut self) {
        let Some(dir) = rfd::FileDialog::new()
            .set_title("Export the movie into a folder")
            .pick_folder()
        else {
            return;
        };
        self.movie.export_request = Some(self.movie.export_ask(dir));
    }
}

fn timing_controls(ui: &mut Ui, clip: &mut Clip, (lo, hi): (f32, f32)) {
    ui.label("Start");
    seconds_field(ui, &mut clip.start, lo, (hi - clip.length).max(lo));
    ui.label("Length");
    seconds_field(
        ui,
        &mut clip.length,
        MIN_CLIP,
        (hi - clip.start).max(MIN_CLIP),
    );
    let names: Vec<&str> = Easing::ALL.iter().map(|e| e.name()).collect();
    let current = Easing::ALL.iter().position(|e| *e == clip.easing);
    if let Some(i) = widgets::segmented(ui, &names, current) {
        clip.easing = Easing::ALL[i];
    }
}

fn effect_controls(ui: &mut Ui, effect: &mut Effect, playable: &[Playable], pose: ViewPose) {
    match effect {
        Effect::Rotate { axis, degrees } => {
            let names: Vec<&str> = Axis::ALL.iter().map(|a| a.name()).collect();
            let current = Axis::ALL.iter().position(|a| a == axis);
            if let Some(i) = widgets::segmented(ui, &names, current) {
                *axis = Axis::ALL[i];
            }
            ui.add(DragValue::new(degrees).speed(1.0).suffix("\u{b0}"));
        }
        Effect::Zoom { from, to } => {
            factor_field(ui, from);
            ui.label("to");
            factor_field(ui, to);
        }
        Effect::MoveTo(target) => {
            if widgets::button(ui, icon::CROSSHAIR, "Use current view", Variant::Secondary)
                .clicked()
            {
                *target = pose;
            }
        }
        Effect::Play {
            structure,
            first,
            last,
            fps,
            looping,
        } => {
            play_controls(ui, playable, structure, [first, last], fps, looping);
        }
    }
}

fn play_controls(
    ui: &mut Ui,
    playable: &[Playable],
    structure: &mut u32,
    [first, last]: [&mut usize; 2],
    fps: &mut f32,
    looping: &mut bool,
) {
    let names: Vec<&str> = playable.iter().map(|p| p.label.as_str()).collect();
    let current = playable.iter().position(|p| p.raw == *structure);
    if let Some(i) = widgets::select(ui, "movie-play-structure", &names, current, "structure") {
        *structure = playable[i].raw;
        *first = 0;
        *last = playable[i].frames - 1;
    }
    let top = current.map_or(usize::MAX, |i| playable[i].frames - 1);
    ui.label("frames");
    ui.add(DragValue::new(first).range(0..=top.min(*last)));
    ui.label("to");
    ui.add(DragValue::new(last).range(*first..=top));
    ui.add(DragValue::new(fps).range(0.5..=240.0).suffix(" fps"));
    if widgets::switch(ui, *looping, "Loop").clicked() {
        *looping = !*looping;
    }
}

fn seconds_field(ui: &mut Ui, value: &mut f32, min: f32, max: f32) -> bool {
    let field = DragValue::new(value)
        .range(min..=max)
        .speed(0.05)
        .suffix(" s");
    ui.add(field).changed()
}

fn factor_field(ui: &mut Ui, value: &mut f32) -> bool {
    let field = DragValue::new(value)
        .range(0.05..=20.0)
        .speed(0.01)
        .prefix("\u{d7}");
    ui.add(field).changed()
}

fn size_field(ui: &mut Ui, value: &mut u32) -> bool {
    ui.add(DragValue::new(value).range(16..=7680).speed(8.0))
        .changed()
}
