//! `sequence color|track|tracks|legend`: the Sequence panel header's
//! choices as commands.

use super::color::{SeqColor, SCHEMES};
use super::tracks::{provider, PROVIDERS};
use crate::commands::parse_on_off;
use crate::ui::AppUi;

impl AppUi<'_> {
    pub(crate) fn sequence_command(&mut self, rest: &str) -> Result<String, String> {
        let (verb, arg) = rest.split_once(char::is_whitespace).unwrap_or((rest, ""));
        let arg = arg.trim();
        match verb {
            "" => Ok(self.sequence_summary()),
            "color" | "colour" => self.sequence_color(arg),
            "track" => self.sequence_track(arg),
            "tracks" => self.sequence_tracks(arg),
            "legend" | "key" => self.sequence_legend(arg),
            other => Err(format!("unknown sequence option `{other}`")),
        }
    }

    fn sequence_summary(&self) -> String {
        format!(
            "sequence color {}; tracks: {}",
            self.sequence.color.word(),
            self.sequence_track_list()
        )
    }

    fn sequence_track_list(&self) -> String {
        match self.sequence.tracks.len() {
            0 => "none".to_string(),
            _ => PROVIDERS
                .iter()
                .filter(|p| self.sequence.track_on(p.id()))
                .map(|p| p.id())
                .collect::<Vec<_>>()
                .join(", "),
        }
    }

    fn sequence_color(&mut self, arg: &str) -> Result<String, String> {
        if arg.is_empty() {
            let words: Vec<&str> = SCHEMES.iter().map(|s| s.1).collect();
            return Ok(format!("sequence color {}", words.join("|")));
        }
        let scheme = SeqColor::parse(arg).ok_or_else(|| {
            let words: Vec<&str> = SCHEMES.iter().map(|s| s.1).collect();
            format!("unknown color `{arg}`; expected {}", words.join(", "))
        })?;
        self.sequence.color = scheme;
        Ok(format!("sequence color {}", scheme.word()))
    }

    fn sequence_track(&mut self, arg: &str) -> Result<String, String> {
        let (name, state) = arg.split_once(char::is_whitespace).unwrap_or((arg, ""));
        let p = provider(name).ok_or_else(|| self.unknown_track(name))?;
        let on = parse_on_off(state.trim(), self.sequence.track_on(p.id()))?;
        self.sequence.set_track(p.id(), on);
        Ok(format!(
            "sequence track {} {}",
            p.id(),
            if on { "on" } else { "off" }
        ))
    }

    fn unknown_track(&self, name: &str) -> String {
        let ids: Vec<&str> = PROVIDERS.iter().map(|p| p.id()).collect();
        format!("unknown track `{name}`; expected {}", ids.join(", "))
    }

    fn sequence_tracks(&mut self, arg: &str) -> Result<String, String> {
        match arg {
            "" => Ok(format!("tracks on: {}", self.sequence_track_list())),
            "all" | "none" => {
                for p in PROVIDERS {
                    self.sequence.set_track(p.id(), arg == "all");
                }
                Ok(format!("tracks: {}", self.sequence_track_list()))
            }
            _ => Err("expected `sequence tracks [all|none]`".into()),
        }
    }

    fn sequence_legend(&mut self, arg: &str) -> Result<String, String> {
        self.sequence.legend = parse_on_off(arg, self.sequence.legend)?;
        Ok(format!(
            "sequence key {}",
            if self.sequence.legend { "on" } else { "off" }
        ))
    }
}
