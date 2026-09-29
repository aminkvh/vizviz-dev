//! The Sequence panel's choices as JSON, for the session file's `view`
//! block. Names are the command words, so a file survives new tracks or
//! schemes: an unknown word is dropped, a missing key changes nothing.

use serde_json::{json, Value};
use vv_core::antibody::{CdrDefinition, Scheme};

use super::color::SeqColor;
use super::tracks::{provider, PROVIDERS};
use super::SequenceState;

impl SequenceState {
    pub fn snapshot(&self) -> Value {
        let tracks: Vec<&str> = PROVIDERS
            .iter()
            .filter(|p| self.track_on(p.id()))
            .map(|p| p.id())
            .collect();
        json!({
            "color": self.color.word(),
            "tracks": tracks,
            "legend": self.legend,
            "antibody": {
                "scheme": self.antibody.scheme.name(),
                "cdr": self.antibody.cdr.name(),
            },
        })
    }

    pub fn restore(&mut self, saved: &Value) {
        if let Some(color) = saved["color"].as_str().and_then(SeqColor::parse) {
            self.color = color;
        }
        if let Some(ids) = saved["tracks"].as_array() {
            self.tracks.clear();
            for p in ids.iter().filter_map(|v| provider(v.as_str()?)) {
                self.set_track(p.id(), true);
            }
        }
        if let Some(legend) = saved["legend"].as_bool() {
            self.legend = legend;
        }
        if let Some(s) = saved["antibody"]["scheme"].as_str().and_then(Scheme::parse) {
            self.antibody.scheme = s;
        }
        if let Some(d) = saved["antibody"]["cdr"]
            .as_str()
            .and_then(CdrDefinition::parse)
        {
            self.antibody.cdr = d;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn choices_round_trip_through_json() {
        let mut state = SequenceState {
            color: SeqColor::Chemistry,
            ..Default::default()
        };
        state.set_track("antibody", true);
        state.set_track("ss", true);
        state.legend = true;
        state.antibody.scheme = Scheme::Imgt;
        state.antibody.cdr = CdrDefinition::North;

        let text = serde_json::to_string(&state.snapshot()).unwrap();
        let mut back = SequenceState::default();
        back.restore(&serde_json::from_str(&text).unwrap());

        assert_eq!(back.color, SeqColor::Chemistry);
        assert!(back.track_on("antibody") && back.track_on("ss") && back.tracks.len() == 2);
        assert!(back.legend);
        assert_eq!(back.antibody, state.antibody);
    }

    #[test]
    fn unknown_words_are_dropped_and_missing_keys_change_nothing() {
        let mut state = SequenceState {
            color: SeqColor::Taylor,
            ..Default::default()
        };
        state.restore(&json!({"tracks": ["ss", "nonsense"], "antibody": {"scheme": "bogus"}}));
        assert_eq!(state.color, SeqColor::Taylor);
        assert_eq!(state.tracks, ["ss"]);
        assert_eq!(state.antibody.scheme, Scheme::Kabat);
    }
}
