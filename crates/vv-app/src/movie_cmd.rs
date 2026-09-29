//! The `movie` verb: every Movie panel action as a command, for scripts,
//! `--exec` and the palette. Tracks and clips are numbered from 1 here.

use std::collections::BTreeMap;
use std::path::PathBuf;

use crate::layout::Tab;
use crate::movie::{Axis, Clip, Easing, Effect, TrackKind, ViewPose};
use crate::movie_gif::GifSettings;
use crate::movie_state::{clock, frame_counts, ExportAsk};
use crate::ui::{AppUi, LayoutRequest};

pub const USAGE: &str = "movie [addtrack KIND | clip TRACK START LENGTH ARGS [linear|ease] | removetrack TRACK | removeclip TRACK CLIP | set fps|duration|size VALUE | time SECONDS | play [on|off] | list | clear | export DIR [nossaa] [nomp4] [nopng] [gif] [nodither] [gifwidth=N]]";

pub const HELP: &str = "The movie maker. Bare `movie` opens the panel. Tracks (KIND rotate, zoom, moveto, play) hold clips on a shared time axis, in seconds; tracks and clips count from 1. `clip` takes: rotate AXIS DEGREES; zoom FROM TO (factors on the distance); moveto (flies to the view on screen now); play STRUCTURE FIRST LAST FPS [loop] (STRUCTURE is `#N` from `structures` or `current`). `time` shows that moment in the viewport; `export` writes numbered PNGs (and movie.mp4 when ffmpeg is installed) into DIR, one movie frame at a time; `gif` also writes movie.gif (256 colours per frame, dithered unless `nodither`, at most 800 px wide or `gifwidth=N`), `nopng` skips the PNGs and the MP4 built from them.";

/// What `movie clip` needs to know about the app to build an effect.
pub struct ClipContext<'a> {
    pub pose: ViewPose,
    pub current: Option<u32>,
    pub frame_counts: &'a BTreeMap<u32, usize>,
}

/// The effect a `clip` line's arguments describe, and its easing.
pub fn parse_clip(
    kind: TrackKind,
    args: &[&str],
    ctx: &ClipContext,
) -> Result<(Effect, Easing), String> {
    let easing = args
        .iter()
        .find_map(|w| Easing::parse(w))
        .unwrap_or(Easing::EaseInOut);
    let words: Vec<&str> = args
        .iter()
        .copied()
        .filter(|w| Easing::parse(w).is_none())
        .collect();
    let effect = match (kind, words.as_slice()) {
        (TrackKind::Rotate, [axis, degrees]) => Effect::Rotate {
            axis: Axis::parse(axis).ok_or("the axis is x, y or z")?,
            degrees: number(degrees, "degrees")?,
        },
        (TrackKind::Zoom, [from, to]) => Effect::Zoom {
            from: number(from, "FROM")?,
            to: number(to, "TO")?,
        },
        (TrackKind::MoveTo, []) => Effect::MoveTo(ctx.pose),
        (TrackKind::Play, [structure, first, last, fps, rest @ ..]) if rest.len() <= 1 => {
            play_effect(ctx, structure, [first, last, fps], rest.first().copied())?
        }
        _ => {
            return Err(format!(
                "usage for a {} clip: {}",
                kind.name(),
                clip_usage(kind)
            ))
        }
    };
    Ok((effect, easing))
}

fn clip_usage(kind: TrackKind) -> &'static str {
    match kind {
        TrackKind::Rotate => "movie clip TRACK START LENGTH x|y|z DEGREES [linear|ease]",
        TrackKind::Zoom => "movie clip TRACK START LENGTH FROM TO [linear|ease]",
        TrackKind::MoveTo => "movie clip TRACK START LENGTH [linear|ease]",
        TrackKind::Play => "movie clip TRACK START LENGTH #N|current FIRST LAST FPS [loop]",
    }
}

fn play_effect(
    ctx: &ClipContext,
    structure: &str,
    [first, last, fps]: [&str; 3],
    looping: Option<&str>,
) -> Result<Effect, String> {
    let structure = match structure {
        "current" => ctx.current.ok_or("no structure loaded")?,
        other => other
            .trim_start_matches('#')
            .parse()
            .map_err(|_| format!("`{other}` is not a structure (use #N or current)"))?,
    };
    let count = ctx
        .frame_counts
        .get(&structure)
        .ok_or_else(|| format!("no structure #{structure}"))?;
    let (first, last): (usize, usize) = (
        first.parse().map_err(|_| "FIRST is a frame number")?,
        last.parse().map_err(|_| "LAST is a frame number")?,
    );
    if last >= *count {
        return Err(format!("#{structure} has frames 0-{}", count - 1));
    }
    let looping = match looping {
        None => false,
        Some("loop") => true,
        Some(other) => return Err(format!("expected `loop`, got `{other}`")),
    };
    Ok(Effect::Play {
        structure,
        first,
        last,
        fps: number(fps, "FPS")?,
        looping,
    })
}

fn number(word: &str, what: &str) -> Result<f32, String> {
    word.parse::<f32>()
        .ok()
        .filter(|v| v.is_finite())
        .ok_or_else(|| format!("{what} must be a number, got `{word}`"))
}

/// A 1-based track or clip number as an index.
fn index(word: &str, what: &str) -> Result<usize, String> {
    word.parse::<usize>()
        .ok()
        .and_then(|n| n.checked_sub(1))
        .ok_or_else(|| format!("{what} numbers start at 1, got `{word}`"))
}

fn parse_size(word: &str) -> Result<(u32, u32), String> {
    let bad = || format!("size is WIDTHxHEIGHT, got `{word}`");
    let (w, h) = word.split_once(['x', 'X']).ok_or_else(bad)?;
    Ok((w.parse().map_err(|_| bad())?, h.parse().map_err(|_| bad())?))
}

impl AppUi<'_> {
    pub(crate) fn run_movie(&mut self, rest: &str) -> Result<String, String> {
        let words: Vec<&str> = rest.split_whitespace().collect();
        let line = match words.as_slice() {
            [] => {
                *self.layout_request = Some(LayoutRequest::OpenPanel(Tab::Movie));
                "showing Movie".into()
            }
            ["addtrack", kind] => self.movie_add_track(kind)?,
            ["clip", track, start, length, args @ ..] => {
                self.movie_add_clip(track, start, length, args)?
            }
            ["removetrack", track] => self.movie_remove_track(track)?,
            ["removeclip", track, clip] => self.movie_remove_clip(track, clip)?,
            ["set", what, value] => self.movie_set(what, value)?,
            ["time", seconds] => {
                self.movie.set_time(number(seconds, "SECONDS")?);
                format!("time {}", clock(self.movie.time))
            }
            ["play"] => self.movie_play(None)?,
            ["play", state] => self.movie_play(Some(state))?,
            ["list"] => self.movie_list(),
            ["clear"] => {
                self.movie.movie.tracks.clear();
                self.movie.selected = None;
                "movie cleared".into()
            }
            ["export", dir, flags @ ..] => self.movie_export(dir, flags)?,
            _ => return Err(format!("usage: {USAGE}")),
        };
        self.movie.touch();
        Ok(line)
    }

    fn movie_add_track(&mut self, kind: &str) -> Result<String, String> {
        let kind = TrackKind::parse(kind).ok_or("KIND is rotate, zoom, moveto or play")?;
        let track = self.movie.movie.add_track(kind);
        Ok(format!("track {}: {}", track + 1, kind.label()))
    }

    fn movie_add_clip(
        &mut self,
        track: &str,
        start: &str,
        length: &str,
        args: &[&str],
    ) -> Result<String, String> {
        let track = index(track, "track")?;
        let kind = self
            .movie
            .movie
            .tracks
            .get(track)
            .ok_or("no such track")?
            .kind;
        let counts = frame_counts(self.scene);
        let ctx = ClipContext {
            pose: ViewPose::of(self.camera),
            current: self.current().map(|id| id.to_raw()),
            frame_counts: &counts,
        };
        let (effect, easing) = parse_clip(kind, args, &ctx)?;
        let clip = Clip {
            start: number(start, "START")?,
            length: number(length, "LENGTH")?,
            easing,
            effect,
        };
        let at = self.movie.movie.add_clip(track, clip)?;
        self.movie.selected = Some((track, at));
        Ok(format!(
            "track {} clip {}: {}",
            track + 1,
            at + 1,
            clip.effect.label()
        ))
    }

    fn movie_remove_track(&mut self, track: &str) -> Result<String, String> {
        let track = index(track, "track")?;
        if !self.movie.movie.remove_track(track) {
            return Err("no such track".into());
        }
        self.movie.selected = None;
        Ok(format!("removed track {}", track + 1))
    }

    fn movie_remove_clip(&mut self, track: &str, clip: &str) -> Result<String, String> {
        let (track, clip) = (index(track, "track")?, index(clip, "clip")?);
        if !self.movie.movie.remove_clip(track, clip) {
            return Err("no such clip".into());
        }
        self.movie.selected = None;
        Ok(format!("removed track {} clip {}", track + 1, clip + 1))
    }

    fn movie_set(&mut self, what: &str, value: &str) -> Result<String, String> {
        let movie = &mut self.movie.movie;
        match what {
            "fps" => movie.set_fps(number(value, "fps")?)?,
            "duration" => movie.set_duration(number(value, "duration")?)?,
            "size" => {
                let (w, h) = parse_size(value)?;
                movie.set_size(w, h);
            }
            _ => return Err("set takes fps, duration or size".into()),
        }
        Ok(format!(
            "{} fps, {} s, {}x{}",
            movie.fps, movie.duration, movie.width, movie.height
        ))
    }

    fn movie_play(&mut self, state: Option<&str>) -> Result<String, String> {
        let want = match state {
            None => !self.movie.playing,
            Some("on") => true,
            Some("off") => false,
            Some(other) => return Err(format!("expected on or off, got `{other}`")),
        };
        if want != self.movie.playing {
            self.movie.toggle_play();
        }
        Ok(if want { "playing" } else { "paused" }.into())
    }

    fn movie_list(&self) -> String {
        let movie = &self.movie.movie;
        let mut lines = vec![format!(
            "{} fps, {} s, {}x{}",
            movie.fps, movie.duration, movie.width, movie.height
        )];
        for (t, track) in movie.tracks.iter().enumerate() {
            lines.push(format!("track {}: {}", t + 1, track.kind.label()));
            lines.extend(track.clips.iter().enumerate().map(|(c, clip)| {
                format!(
                    "  clip {}: {:.2}-{:.2} s, {}, {}",
                    c + 1,
                    clip.start,
                    clip.end(),
                    clip.easing.name(),
                    clip.effect.label()
                )
            }));
        }
        lines.join("\n")
    }

    fn movie_export(&mut self, dir: &str, flags: &[&str]) -> Result<String, String> {
        self.movie.export_request = Some(export_flags(dir.into(), flags)?);
        Ok(format!("exporting the movie to {dir}"))
    }
}

/// The export `movie export DIR FLAGS` asks for. The verb's defaults are
/// the panel's: PNGs and an MP4 (when ffmpeg exists), no GIF.
fn export_flags(dir: PathBuf, flags: &[&str]) -> Result<ExportAsk, String> {
    let has = |name: &str| flags.contains(&name);
    let width = gif_width_flag(flags)?;
    let known = |f: &&str| {
        ["nossaa", "nomp4", "nopng", "gif", "nodither"].contains(f) || f.starts_with("gifwidth=")
    };
    if let Some(bad) = flags.iter().find(|f| !known(f)) {
        return Err(format!(
            "expected nossaa, nomp4, nopng, gif, nodither or gifwidth=N, got `{bad}`"
        ));
    }
    if !has("gif") && (has("nodither") || width.is_some()) {
        return Err("nodither and gifwidth=N only apply with gif".into());
    }
    if has("nopng") && !has("gif") {
        return Err("nopng with no gif would write nothing".into());
    }
    let defaults = GifSettings::default();
    Ok(ExportAsk {
        dir,
        ssaa: !has("nossaa"),
        png: !has("nopng"),
        mp4: !has("nomp4") && !has("nopng"),
        gif: has("gif").then_some(GifSettings {
            width: width.unwrap_or(defaults.width),
            dither: !has("nodither"),
        }),
    })
}

fn gif_width_flag(flags: &[&str]) -> Result<Option<u32>, String> {
    let Some(value) = flags.iter().find_map(|f| f.strip_prefix("gifwidth=")) else {
        return Ok(None);
    };
    match value.parse::<u32>() {
        Ok(w) if w >= 16 => Ok(Some(w)),
        _ => Err(format!(
            "gifwidth=N needs a width of 16 or more, got `{value}`"
        )),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use glam::Vec3;
    use vv_render::Camera;

    fn parse(kind: TrackKind, line: &str) -> Result<(Effect, Easing), String> {
        let counts = BTreeMap::from([(3, 50)]);
        let ctx = ClipContext {
            pose: ViewPose::of(&Camera::framing(Vec3::ZERO, 5.0)),
            current: Some(3),
            frame_counts: &counts,
        };
        let args: Vec<&str> = line.split_whitespace().collect();
        parse_clip(kind, &args, &ctx)
    }

    #[test]
    fn a_rotate_clip_reads_axis_degrees_and_easing() {
        let (effect, easing) = parse(TrackKind::Rotate, "y 360 linear").unwrap();
        assert_eq!(
            effect,
            Effect::Rotate {
                axis: Axis::Y,
                degrees: 360.0
            }
        );
        assert_eq!(easing, Easing::Linear);
        assert_eq!(
            parse(TrackKind::Rotate, "y 90").unwrap().1,
            Easing::EaseInOut
        );
        assert!(parse(TrackKind::Rotate, "w 90").is_err());
        assert!(parse(TrackKind::Rotate, "y").is_err());
    }

    #[test]
    fn a_zoom_clip_reads_two_factors() {
        let (effect, _) = parse(TrackKind::Zoom, "1 0.5").unwrap();
        assert_eq!(effect, Effect::Zoom { from: 1.0, to: 0.5 });
        assert!(parse(TrackKind::Zoom, "1").is_err());
    }

    #[test]
    fn a_move_clip_takes_the_view_on_screen() {
        let (effect, _) = parse(TrackKind::MoveTo, "").unwrap();
        assert!(matches!(effect, Effect::MoveTo(_)));
        assert!(parse(TrackKind::MoveTo, "extra").is_err());
    }

    #[test]
    fn a_play_clip_names_its_structure_and_checks_the_frames() {
        let (effect, _) = parse(TrackKind::Play, "#3 0 49 24 loop").unwrap();
        assert_eq!(
            effect,
            Effect::Play {
                structure: 3,
                first: 0,
                last: 49,
                fps: 24.0,
                looping: true
            }
        );
        assert!(matches!(
            parse(TrackKind::Play, "current 5 10 12").unwrap().0,
            Effect::Play { structure: 3, .. }
        ));
        assert!(
            parse(TrackKind::Play, "#3 0 50 24").is_err(),
            "past the last frame"
        );
        assert!(
            parse(TrackKind::Play, "#9 0 1 24").is_err(),
            "no such structure"
        );
        assert!(parse(TrackKind::Play, "#3 0 1 24 spin").is_err());
    }

    #[test]
    fn numbers_and_sizes_are_checked() {
        assert_eq!(index("1", "track"), Ok(0));
        assert!(index("0", "track").is_err());
        assert_eq!(parse_size("1920x1080"), Ok((1920, 1080)));
        assert!(parse_size("1920").is_err());
        assert!(number("nan", "x").is_err());
    }

    fn ask(flags: &[&str]) -> Result<ExportAsk, String> {
        export_flags(PathBuf::from("out"), flags)
    }

    #[test]
    fn a_plain_export_keeps_its_old_meaning() {
        let plain = ask(&[]).unwrap();
        assert!(plain.ssaa && plain.png && plain.mp4 && plain.gif.is_none());
        let trimmed = ask(&["nossaa", "nomp4"]).unwrap();
        assert!(!trimmed.ssaa && trimmed.png && !trimmed.mp4);
    }

    #[test]
    fn gif_takes_a_width_and_can_turn_dithering_off() {
        let on = ask(&["gif"]).unwrap().gif.unwrap();
        assert_eq!(on, GifSettings::default());
        let custom = ask(&["gif", "gifwidth=480", "nodither"])
            .unwrap()
            .gif
            .unwrap();
        assert_eq!(
            custom,
            GifSettings {
                width: 480,
                dither: false
            }
        );
    }

    #[test]
    fn nopng_drops_the_mp4_built_from_the_pngs() {
        let gif_only = ask(&["gif", "nopng"]).unwrap();
        assert!(!gif_only.png && !gif_only.mp4 && gif_only.gif.is_some());
    }

    #[test]
    fn contradictory_or_unknown_flags_are_refused() {
        assert!(ask(&["nopng"]).is_err());
        assert!(ask(&["nodither"]).is_err());
        assert!(ask(&["gifwidth=640"]).is_err());
        assert!(ask(&["gif", "gifwidth=abc"]).is_err());
        assert!(ask(&["gif", "gifwidth=4"]).is_err());
        assert!(ask(&["mkv"]).is_err());
    }
}
