//! The keyboard, in one place: shortcuts (each runs a command line, as the
//! ribbon does), key tips (tap Alt or press F10, then the letters shown on
//! the ribbon), the tool strip down the left edge (the viewport's mouse
//! modes and tools), and Escape's close-the-top-thing
//! chain. Nothing here fires while a text field has focus, and
//! single-letter shortcuts wait while key tips show.

use egui::{Key, Modifiers, Ui};
use egui_phosphor::regular as icon;

use crate::ribbon::{KeyTips, TABS};
use crate::ui::{AppUi, MouseMode};
use vv_scene::Command;

pub struct Shortcut {
    pub ctrl: bool,
    pub shift: bool,
    pub key: Key,
    pub command: &'static str,
}

const fn ctrl(key: Key, command: &'static str) -> Shortcut {
    Shortcut {
        ctrl: true,
        shift: false,
        key,
        command,
    }
}

const fn plain(key: Key, command: &'static str) -> Shortcut {
    Shortcut {
        ctrl: false,
        shift: false,
        key,
        command,
    }
}

pub const SHORTCUTS: &[Shortcut] = &[
    ctrl(Key::Z, "undo"),
    ctrl(Key::Y, "redo"),
    Shortcut {
        ctrl: true,
        shift: true,
        key: Key::Z,
        command: "redo",
    },
    ctrl(Key::O, "open"),
    ctrl(Key::S, "savesession"),
    ctrl(Key::Comma, "preferences"),
    Shortcut {
        ctrl: true,
        shift: true,
        key: Key::S,
        command: "screenshot",
    },
    ctrl(Key::Num1, "panel structures"),
    ctrl(Key::Num2, "panel selection"),
    ctrl(Key::Num3, "panel inspector"),
    ctrl(Key::Num4, "panel look"),
    ctrl(Key::Num5, "panel sequence"),
    ctrl(Key::Num6, "panel log"),
    ctrl(Key::Num7, "panel info"),
    ctrl(Key::Num8, "panel timeline"),
    // No Ctrl+9: the viewport is always visible, nothing to show/hide.
    // Mouse-mode keys.
    plain(Key::R, "mode rotate"),
    plain(Key::T, "mode translate"),
    plain(Key::S, "mode scale"),
    plain(Key::C, "mode center"),
    plain(Key::F, "view reset"),
    plain(Key::Q, "selectmode shape next"),
    plain(Key::Space, "play"),
];

/// What a key tip letter reached.
#[derive(Debug, PartialEq, Eq)]
pub enum TipHit {
    Tab(usize),
    /// Action `action` of group `group` of tab `tab`.
    Action {
        tab: usize,
        group: usize,
        action: usize,
    },
}

/// The key tip `typed` names in `mode`, with `tab` open.
pub fn key_tip(mode: KeyTips, tab: usize, typed: &str) -> Option<TipHit> {
    match mode {
        KeyTips::Off => None,
        KeyTips::Tabs => TABS
            .iter()
            .position(|t| t.tip.eq_ignore_ascii_case(typed))
            .map(TipHit::Tab),
        KeyTips::Actions => TABS[tab].groups.iter().enumerate().find_map(|(g, group)| {
            group
                .actions
                .iter()
                .position(|a| a.tip.eq_ignore_ascii_case(typed))
                .map(|a| TipHit::Action {
                    tab,
                    group: g,
                    action: a,
                })
        }),
    }
}

fn letter(key: Key) -> Option<&'static str> {
    Some(match key {
        Key::A => "A",
        Key::B => "B",
        Key::C => "C",
        Key::D => "D",
        Key::E => "E",
        Key::F => "F",
        Key::G => "G",
        Key::H => "H",
        Key::I => "I",
        Key::J => "J",
        Key::K => "K",
        Key::L => "L",
        Key::M => "M",
        Key::N => "N",
        Key::O => "O",
        Key::P => "P",
        Key::Q => "Q",
        Key::R => "R",
        Key::S => "S",
        Key::T => "T",
        Key::U => "U",
        Key::V => "V",
        Key::W => "W",
        Key::X => "X",
        Key::Y => "Y",
        Key::Z => "Z",
        Key::Num1 => "1",
        Key::Num2 => "2",
        Key::Num3 => "3",
        Key::Num4 => "4",
        Key::Num5 => "5",
        Key::Num6 => "6",
        Key::Num7 => "7",
        Key::Num8 => "8",
        Key::Num9 => "9",
        _ => return None,
    })
}

/// A key's display name in a shortcut string ("Z", "1", ",", "Space").
fn key_name(key: Key) -> &'static str {
    match key {
        Key::Comma => ",",
        Key::Space => "Space",
        other => letter(other).unwrap_or("?"),
    }
}

/// A shortcut formatted for a tooltip or the palette ("Ctrl+Z",
/// "Ctrl+Shift+S").
fn format_shortcut(s: &Shortcut) -> String {
    let mut parts = Vec::new();
    if s.ctrl {
        parts.push("Ctrl");
    }
    if s.shift {
        parts.push("Shift");
    }
    parts.push(key_name(s.key));
    parts.join("+")
}

/// The shortcut that runs `command` exactly, formatted for a tooltip or
/// the palette ("Ctrl+Z", "Ctrl+Shift+S"). Only an exact match: a
/// shortcut for one arm of a multi-arg command (e.g. `mode rotate`)
/// would mislabel that command's whole family.
pub fn shortcut_for(command: &str) -> Option<String> {
    SHORTCUTS
        .iter()
        .find(|s| s.command == command)
        .map(format_shortcut)
}

/// Every shortcut bound to one arm of the `verb` family ("mode rotate",
/// "mode translate", ... for `verb` "mode"), joined for a palette row
/// whose own id is the bare verb -- item 5c: `shortcut_for("mode")` is
/// `None` (correctly: no single key names the whole family), but the
/// family itself is still worth showing somewhere, e.g. "R/T/S/C".
pub fn shortcut_family(verb: &str) -> Option<String> {
    let mut keys: Vec<String> = SHORTCUTS
        .iter()
        .filter(|s| {
            s.command
                .strip_prefix(verb)
                .is_some_and(|rest| rest.starts_with(' '))
        })
        .map(format_shortcut)
        .collect();
    keys.dedup();
    (!keys.is_empty()).then(|| keys.join("/"))
}

impl AppUi<'_> {
    /// Runs this frame's shortcuts and key tips. Call before the panels,
    /// so a key that a key tip consumed never reaches the viewport.
    pub fn handle_keys(&mut self, ctx: &egui::Context) {
        if ctx.egui_wants_keyboard_input() {
            self.ribbon.key_tips = KeyTips::Off;
            return;
        }
        // Escape closes the top overlay first (key tips, a popover,
        // the palette -- each closes itself on this same unconsumed key,
        // later this frame); only then does it exit Label mode, else
        // clear the selection. Must run before `ribbon_ui`/`palette_ui`
        // render, so a popover/the palette reads as open before they can
        // self-close it this same frame.
        if ctx.input(|i| i.key_pressed(Key::Escape))
            && !self.any_dialog_open()
            && self.ribbon.key_tips == KeyTips::Off
            && !egui::Popup::is_any_open(ctx)
            && !self.palette.open
        {
            if self.view.mouse_mode == MouseMode::Label {
                self.view.mouse_mode = MouseMode::Rotate;
            } else if self.view.select_tool.draws() {
                self.view.select_tool.shape = crate::select_tool::SelectShape::Click;
            } else if self.scene.active_selection().is_some() {
                self.dispatch(Command::ClearSelection);
            }
        }
        let (alt_now, pressed): (bool, Vec<(Key, Modifiers)>) = ctx.input(|i| {
            let keys = i
                .events
                .iter()
                .filter_map(|e| match e {
                    egui::Event::Key {
                        key,
                        pressed: true,
                        modifiers,
                        ..
                    } => Some((*key, *modifiers)),
                    _ => None,
                })
                .collect();
            (i.modifiers.alt, keys)
        });
        if self.ribbon.alt_tapped(alt_now, !pressed.is_empty()) {
            self.ribbon.toggle_key_tips();
        }
        for (key, modifiers) in pressed {
            // The other common key-tips key, for when Alt goes to the system.
            if key == Key::F10 {
                self.ribbon.toggle_key_tips();
                continue;
            }
            if key == Key::Escape && self.ribbon.key_tips != KeyTips::Off {
                self.ribbon.key_tips = KeyTips::Off;
                continue;
            }
            if self.ribbon.key_tips != KeyTips::Off {
                if let Some(typed) = letter(key) {
                    self.consume(ctx, key, modifiers);
                    self.follow_key_tip(typed);
                }
                continue;
            }
            let shortcut = SHORTCUTS.iter().find(|s| {
                s.key == key
                    && s.ctrl == modifiers.command
                    && s.shift == modifiers.shift
                    && !modifiers.alt
            });
            if let Some(s) = shortcut {
                self.consume(ctx, key, modifiers);
                self.run_command_logged(s.command);
            }
        }
    }

    fn consume(&self, ctx: &egui::Context, key: Key, modifiers: Modifiers) {
        ctx.input_mut(|i| i.consume_key(modifiers, key));
    }

    fn follow_key_tip(&mut self, typed: &str) {
        match key_tip(self.ribbon.key_tips, self.ribbon.tab, typed) {
            Some(TipHit::Tab(tab)) => {
                self.ribbon.tab = tab;
                self.prefs.ribbon_collapsed = false;
                self.ribbon.key_tips = KeyTips::Actions;
            }
            Some(TipHit::Action { tab, group, action }) => {
                self.ribbon.key_tips = KeyTips::Off;
                self.press(&TABS[tab].groups[group].actions[action]);
            }
            None => {}
        }
    }

    /// The tool strip: mouse modes (R T S C) and view tools only
    /// (measure and label are Analyze tab actions, clip is View's; every
    /// other button is a command).
    pub fn tool_strip_ui(&mut self, ui: &mut Ui) {
        let mode = self.view.mouse_mode;
        // Item 5c: the shortcut in each tip comes from `SHORTCUTS` itself
        // (via `shortcut_for`), not a hand-typed "(R)" that could drift.
        let tools: [(&str, &str, &str, Option<MouseMode>); 5] = [
            (
                icon::ARROWS_CLOCKWISE,
                "Rotate",
                "mode rotate",
                Some(MouseMode::Rotate),
            ),
            (
                icon::HAND_GRABBING,
                "Translate",
                "mode translate",
                Some(MouseMode::Translate),
            ),
            (
                icon::MAGNIFYING_GLASS_PLUS,
                "Scale",
                "mode scale",
                Some(MouseMode::Scale),
            ),
            (
                icon::CROSSHAIR,
                "Center: click an atom",
                "mode center",
                Some(MouseMode::Center),
            ),
            (icon::FRAME_CORNERS, "Reset view", "view reset", None),
        ];
        let mut run: Option<&str> = None;
        ui.vertical_centered(|ui| {
            for (glyph, tip, command, tool_mode) in tools {
                let tip = match shortcut_for(command) {
                    Some(shortcut) => format!("{tip} ({shortcut})"),
                    None => tip.into(),
                };
                let button = crate::widgets::icon_button(ui, glyph, tool_mode == Some(mode));
                if button.on_hover_text(tip).clicked() {
                    run = Some(command);
                }
            }
        });
        if let Some(command) = run {
            self.run_command_logged(command);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;

    #[test]
    fn shortcuts_are_known_commands_and_distinct() {
        let mut seen = HashSet::new();
        for s in SHORTCUTS {
            crate::commands::validate(s.command).unwrap();
            assert!(
                seen.insert((s.ctrl, s.shift, s.key)),
                "{:?} bound twice",
                s.key
            );
        }
    }

    #[test]
    fn a_lone_alt_tap_toggles_key_tips_and_alt_with_a_key_does_not() {
        let mut r = crate::ribbon::RibbonState::default();
        // Down, up: a tap.
        assert!(!r.alt_tapped(true, false));
        assert!(r.alt_tapped(false, false));
        // Down, a key, up: a chord, not a tap.
        assert!(!r.alt_tapped(true, false));
        assert!(!r.alt_tapped(true, true));
        assert!(!r.alt_tapped(false, false));
        // No Alt at all.
        assert!(!r.alt_tapped(false, true));
    }

    #[test]
    fn key_tips_reach_a_tab_then_its_actions() {
        let file = TABS.iter().position(|t| t.label == "File").unwrap();
        let look = TABS.iter().position(|t| t.label == "Look").unwrap();
        assert_eq!(key_tip(KeyTips::Tabs, file, "l"), Some(TipHit::Tab(look)));
        let Some(TipHit::Action { tab, group, action }) = key_tip(KeyTips::Actions, look, "L")
        else {
            panic!("no action L on the Look tab");
        };
        assert_eq!(TABS[tab].groups[group].actions[action].label, "Lights");
        assert_eq!(key_tip(KeyTips::Off, file, "L"), None);
        assert_eq!(key_tip(KeyTips::Actions, file, "?"), None);
    }

    #[test]
    fn shortcut_for_formats_modifiers_and_key_and_is_exact() {
        assert_eq!(shortcut_for("undo"), Some("Ctrl+Z".into()));
        assert_eq!(shortcut_for("redo"), Some("Ctrl+Y".into()));
        assert_eq!(shortcut_for("preferences"), Some("Ctrl+,".into()));
        // "mode" (bare) matches no shortcut: only "mode rotate" etc. do.
        assert_eq!(shortcut_for("mode"), None);
    }

    #[test]
    fn shortcut_family_joins_a_verbs_arms_and_ignores_unrelated_prefixes() {
        assert_eq!(
            shortcut_family("mode"),
            Some("R/T/S/C".into()),
            "rotate, translate, scale, center -- in SHORTCUTS' own order"
        );
        assert_eq!(
            shortcut_family("panel"),
            Some("Ctrl+1/Ctrl+2/Ctrl+3/Ctrl+4/Ctrl+5/Ctrl+6/Ctrl+7/Ctrl+8".into())
        );
        // "mo" is a prefix of "mode rotate" but not itself followed by a
        // space in that command, so it must not falsely match.
        assert_eq!(shortcut_family("mo"), None);
        assert_eq!(shortcut_family("nonexistent"), None);
    }

    #[test]
    fn every_shortcut_key_has_a_display_name() {
        // Guards `key_name`: iterates every entry directly (not through
        // `shortcut_for`, which only returns one match per command and
        // would skip "redo"'s second binding, Ctrl+Shift+Z).
        for s in SHORTCUTS {
            assert_ne!(key_name(s.key), "?", "{}: unnamed key", s.command);
        }
    }
}
