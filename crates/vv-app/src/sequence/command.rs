//! `sequence color|track|tracks|legend`: the Sequence panel header's
//! choices as commands.

use std::path::PathBuf;

use vv_core::antibody::{CdrDefinition, Scheme};

use super::anarci::Backend;
use super::chain_props;
use super::color::{SeqColor, SCHEMES};
use super::rows::rows_of;
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
            "antibody" => self.sequence_antibody(arg),
            "uniprot" => self.sequence_uniprot(arg),
            "props" => self.sequence_props(),
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
                for p in PROVIDERS.iter().filter(|p| arg == "none" || !p.online()) {
                    self.sequence.set_track(p.id(), arg == "all");
                }
                Ok(format!("tracks: {}", self.sequence_track_list()))
            }
            _ => Err("expected `sequence tracks [all|none]`".into()),
        }
    }

    fn sequence_antibody(&mut self, arg: &str) -> Result<String, String> {
        let (what, name) = arg.split_once(char::is_whitespace).unwrap_or((arg, ""));
        let name = name.trim();
        let settings = &mut self.sequence.antibody;
        match (what, name) {
            ("", _) => {}
            ("scheme", _) => {
                settings.scheme = Scheme::parse(name).ok_or_else(|| {
                    let all: Vec<&str> = Scheme::ALL.iter().map(|s| s.name()).collect();
                    format!("unknown numbering `{name}`; expected {}", all.join(", "))
                })?;
            }
            ("cdr", _) => {
                settings.cdr = CdrDefinition::parse(name).ok_or_else(|| {
                    let all: Vec<&str> = CdrDefinition::ALL.iter().map(|d| d.name()).collect();
                    format!(
                        "unknown CDR definition `{name}`; expected {}",
                        all.join(", ")
                    )
                })?;
            }
            ("backend", _) => {
                settings.backend = Backend::parse(name).ok_or_else(|| {
                    let all: Vec<&str> = Backend::ALL.iter().map(|b| b.name()).collect();
                    format!("unknown backend `{name}`; expected {}", all.join(", "))
                })?;
            }
            ("exe", _) => {
                let path = (!name.is_empty()).then(|| PathBuf::from(name));
                self.sequence.cache.set_anarci_exe(path);
            }
            _ => {
                return Err(
                    "expected `sequence antibody [scheme NAME | cdr NAME | backend native|anarci|abnum | exe PATH]`"
                        .into(),
                )
            }
        }
        let settings = self.sequence.antibody;
        Ok(format!(
            "sequence antibody: numbering {}, CDRs {}, numbers from {}",
            settings.scheme.name(),
            settings.cdr.name(),
            settings.backend.name()
        ))
    }

    fn sequence_uniprot(&mut self, arg: &str) -> Result<String, String> {
        let (track, state) = match arg.strip_prefix("variants") {
            Some(rest) => ("variants", rest.trim()),
            None => ("uniprot", arg),
        };
        let on = parse_on_off(state, self.sequence.track_on(track))?;
        self.sequence.set_track(track, on);
        let what = if track == "variants" {
            "uniprot variants"
        } else {
            "uniprot"
        };
        Ok(format!("sequence {what} {}", if on { "on" } else { "off" }))
    }

    /// One line per protein chain: size, pI, charge and extinction.
    fn sequence_props(&self) -> Result<String, String> {
        let mut lines = Vec::new();
        for row in rows_of(self.scene, |_, _| Default::default()) {
            let Some(loaded) = self.scene.structure(row.structure) else {
                continue;
            };
            let top = &loaded.structure.topology;
            if let Some(p) = chain_props::of_chain(top, row.residues.clone()) {
                lines.push(format!("{}: {}", row.label, chain_props::one_line(&p)));
            }
        }
        match lines.is_empty() {
            true => Err("no protein chain loaded".into()),
            false => Ok(lines.join("\n")),
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
