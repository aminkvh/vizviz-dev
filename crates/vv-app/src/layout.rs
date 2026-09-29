//! Panel layout: which tabs exist, their default arrangement, and saving
//! the arrangement the user actually ends up with.

use crate::theme::UiPrefs;
use std::path::PathBuf;

use egui_dock::{DockState, NodeIndex};
use serde::{Deserialize, Serialize};

#[derive(Clone, Copy, PartialEq, Eq, Debug, Serialize, Deserialize)]
pub enum Tab {
    Viewport,
    /// What's loaded and how each piece is drawn. Formerly named "Scene";
    /// the alias keeps a layout or session saved under that name loading.
    #[serde(alias = "Scene")]
    Structures,
    Selection,
    Inspector,
    Log,
    /// One-letter residue codes per chain, docked above the viewport.
    Sequence,
    /// Removed: its panel is now the Look ribbon tab. Kept only so an
    /// older `layout.json` or session `dock_layout` still deserializes;
    /// `strip_scene_settings` then drops the tab.
    SceneSettings,
    /// Frame scrubber for multi-model files; auto-opens on load.
    Timeline,
    /// Multi-track movie editor: camera moves and trajectory play on a
    /// shared time axis, exported as frames.
    Movie,
    /// The file's own annotations: title, method, resolution, entities.
    Info,
}

impl Tab {
    /// Every tab but the viewport, which never closes.
    pub const PANELS: [Tab; 8] = [
        Tab::Structures,
        Tab::Selection,
        Tab::Inspector,
        Tab::Sequence,
        Tab::Log,
        Tab::Info,
        Tab::Timeline,
        Tab::Movie,
    ];

    pub fn title(self) -> &'static str {
        match self {
            Tab::Viewport => "Viewport",
            Tab::Structures => "Structures",
            Tab::Selection => "Selections",
            Tab::Inspector => "Inspector",
            Tab::Log => "Log",
            Tab::Sequence => "Sequence",
            Tab::SceneSettings => "Look",
            Tab::Timeline => "Timeline",
            Tab::Movie => "Movie",
            Tab::Info => "Info",
        }
    }
}

/// Drops a leftover `Tab::SceneSettings` tab (its panel is gone; see
/// the variant's doc) from a layout read from disk or a session, so
/// it never tries to show. A no-op once nothing on this machine still
/// has one.
pub fn strip_scene_settings(state: &mut DockState<Tab>) {
    while let Some(found) = state.find_tab(&Tab::SceneSettings) {
        state.remove_tab(found);
    }
}

/// The default layout: the viewport in the middle, the Structures panel and
/// selection left, inspector, info and look right, the sequence strip
/// above the viewport and the log across the bottom, those two collapsed.
/// The timeline is not here: it opens itself when a multi-frame structure
/// loads (`open_near_viewport`).
///
/// `Tree::split_*` always returns `[node_with_original_content,
/// node_with_new_tabs]` regardless of direction, so the viewport's (and
/// the Structures panel's) node index has to be re-threaded through every
/// split that touches it — splitting `NodeIndex::root()` a second time
/// would target whatever the first split left behind at that index, not
/// the tab we actually want to keep splitting.
pub fn default_layout() -> DockState<Tab> {
    let mut state = DockState::new(vec![Tab::Viewport]);
    let tree = state.main_surface_mut();
    let root = NodeIndex::root();

    let [viewport, log] = tree.split_below(root, 0.8, vec![Tab::Log]);
    let [viewport, _right] = tree.split_right(viewport, 0.75, vec![Tab::Inspector, Tab::Info]);
    let [viewport, scene] = tree.split_left(viewport, 0.30, vec![Tab::Structures]);
    // Structures is a short list; Selections holds the layers.
    let [_scene, _selection] = tree.split_below(scene, 0.3, vec![Tab::Selection]);
    let [_viewport, sequence] = tree.split_above(viewport, 0.14, vec![Tab::Sequence]);
    for node in [log, sequence] {
        if let egui_dock::Node::Leaf(leaf) = &mut tree[node] {
            leaf.collapsed = true;
        }
    }
    state
}

/// The viewport sits alone in its own open leaf of the main surface: never
/// tabbed with a panel, collapsed or floated. The dock undoes any move that
/// breaks this.
pub fn viewport_in_place(state: &DockState<Tab>) -> bool {
    let tree = state.main_surface();
    tree.find_tab(&Tab::Viewport)
        .and_then(|(node, _)| tree.leaf(node).ok())
        .is_some_and(|leaf| leaf.tabs.len() == 1 && !leaf.collapsed)
}

/// Forces the viewport's own tab bar hidden every frame: no header, no ▼, no +. Cheap enough to call
/// unconditionally; a no-op once already hidden. Needs `DockArea::
/// hidable_tab_bars(true)` and `theme::dock_style`'s zeroed reveal
/// affordances to actually leave nothing in its place.
pub fn hide_viewport_tab_bar(state: &mut DockState<Tab>) {
    if let Some((node, _)) = state.main_surface().find_tab(&Tab::Viewport) {
        if let egui_dock::Node::Leaf(leaf) = &mut state.main_surface_mut()[node] {
            leaf.tab_bar_hidden = true;
        }
    }
}

/// Just the viewport and inspector, for screen recording or a second
/// monitor.
pub fn compact_layout() -> DockState<Tab> {
    let mut state = DockState::new(vec![Tab::Viewport]);
    let tree = state.main_surface_mut();
    tree.split_right(NodeIndex::root(), 0.7, vec![Tab::Inspector]);
    state
}

/// The default plus a timeline under the viewport.
pub fn trajectory_layout() -> DockState<Tab> {
    let mut state = default_layout();
    open_near_viewport(&mut state, Tab::Timeline);
    state
}

/// Selection and sequence front and centre: selection panel left (full
/// height), sequence strip above the viewport, inspector right with the
/// Structures panel behind it, log below.
pub fn analysis_layout() -> DockState<Tab> {
    let mut state = DockState::new(vec![Tab::Viewport]);
    let tree = state.main_surface_mut();
    let root = NodeIndex::root();
    let [viewport, _right] = tree.split_right(root, 0.72, vec![Tab::Inspector, Tab::Structures]);
    let [viewport, _selection] = tree.split_left(viewport, 0.26, vec![Tab::Selection]);
    let [viewport, _log] = tree.split_below(viewport, 0.8, vec![Tab::Log]);
    tree.split_above(viewport, 0.16, vec![Tab::Sequence]);
    state
}

/// The panels not in `state`, in menu order.
pub fn closed_panels(state: &DockState<Tab>) -> Vec<Tab> {
    Tab::PANELS
        .into_iter()
        .filter(|tab| state.find_tab(tab).is_none())
        .collect()
}

/// How much of the viewport's height stays with it when `tab` opens under
/// it: the movie editor needs a taller strip than the timeline.
fn strip_fraction(tab: Tab) -> f32 {
    match tab {
        Tab::Movie => 0.5,
        _ => 0.82,
    }
}

/// Opens `tab` in a strip under the viewport, if it isn't open already.
/// Used for the timeline's auto-open; falls back to the first leaf when
/// the viewport itself has been closed.
pub fn open_near_viewport(state: &mut DockState<Tab>, tab: Tab) {
    if state.find_tab(&tab).is_some() {
        return;
    }
    let tree = state.main_surface_mut();
    match tree.find_tab(&Tab::Viewport) {
        Some((node, _)) => {
            tree.split_below(node, strip_fraction(tab), vec![tab]);
        }
        None => state.push_to_first_leaf(tab),
    }
}

/// Generated fresh every time, never persisted or user-deletable — see
/// `built_in_workspaces` at the call site in `ui.rs` for why these are
/// kept separate from the user-saved list `load_layout`/`save_layout`
/// round-trip (persisting them would duplicate on every launch).
pub fn built_in_workspaces() -> Vec<(&'static str, DockState<Tab>)> {
    vec![
        ("Default", default_layout()),
        ("Trajectory", trajectory_layout()),
        ("Analysis", analysis_layout()),
        ("Compact", compact_layout()),
    ]
}

/// Per-user configuration directory: `%APPDATA%\vizviz` on Windows,
/// `~/Library/Application Support/vizviz` on macOS, `$XDG_CONFIG_HOME/
/// vizviz` (or `~/.config/vizviz`) elsewhere. `None` when no home can be
/// found, in which case state lands in the working directory as before.
pub fn config_dir() -> Option<PathBuf> {
    let var = |name: &str| std::env::var_os(name).filter(|v| !v.is_empty());
    let base = if cfg!(windows) {
        var("APPDATA").map(PathBuf::from)
    } else if cfg!(target_os = "macos") {
        var("HOME").map(|h| PathBuf::from(h).join("Library/Application Support"))
    } else {
        var("XDG_CONFIG_HOME")
            .map(PathBuf::from)
            .or_else(|| var("HOME").map(|h| PathBuf::from(h).join(".config")))
    }?;
    Some(base.join("vizviz"))
}

/// Where the Image and Render dialogs save by default: Pictures\vizviz
/// (the current directory if the home folder is unknown).
pub fn export_dir() -> PathBuf {
    let home = std::env::var_os(if cfg!(windows) { "USERPROFILE" } else { "HOME" })
        .filter(|v| !v.is_empty());
    home.map_or_else(PathBuf::new, |h| {
        PathBuf::from(h).join("Pictures").join("vizviz")
    })
}

/// Where the layout lives. Used to be a bare relative path, which meant a
/// `vizviz` on PATH scattered `vizviz_layout.json` into every directory it
/// was launched from and never saw the same layout twice.
fn layout_path() -> PathBuf {
    config_dir()
        .map(|dir| dir.join("layout.json"))
        .unwrap_or_else(|| PathBuf::from("vizviz_layout.json"))
}

/// The pre-config-dir location, read once as a fallback so an existing
/// layout survives the move.
fn legacy_layout_path() -> PathBuf {
    PathBuf::from("vizviz_layout.json")
}

/// Everything `layout.json` keeps: the saved workspaces, the start layout
/// and the UI preferences.
#[derive(Serialize, Deserialize)]
pub struct SavedLayout {
    /// The arrangement to open: the start layout, else the built-in
    /// default. Never saved, so a launch never inherits the last session's.
    #[serde(skip, default = "default_layout")]
    pub launch: DockState<Tab>,
    #[serde(default)]
    pub workspaces: Vec<(String, DockState<Tab>)>,
    #[serde(default)]
    pub prefs: UiPrefs,
    /// Set (`startlayout`): every launch opens this arrangement, whatever
    /// the last one was.
    #[serde(default)]
    pub start: Option<DockState<Tab>>,
}

impl Default for SavedLayout {
    fn default() -> Self {
        Self {
            launch: default_layout(),
            workspaces: Vec::new(),
            prefs: UiPrefs::default(),
            start: None,
        }
    }
}

/// The saved state, launching with the start layout when one is set. A
/// missing or unreadable file
/// (one from before workspaces existed no longer parses) gives the
/// built-in default: per-user local state, not worth a migration.
pub fn load_layout() -> SavedLayout {
    let mut saved = [layout_path(), legacy_layout_path()]
        .iter()
        .find_map(|path| {
            let text = std::fs::read_to_string(path).ok()?;
            serde_json::from_str::<SavedLayout>(&text).ok()
        })
        .unwrap_or_default();
    strip_scene_settings(&mut saved.launch);
    for (_, workspace) in &mut saved.workspaces {
        strip_scene_settings(workspace);
    }
    if let Some(start) = &mut saved.start {
        strip_scene_settings(start);
    }
    if let Some(start) = saved.start.as_ref().filter(|s| viewport_in_place(s)) {
        saved.launch = start.clone();
    }
    saved
}

pub fn save_layout(saved: &SavedLayout) {
    let path = layout_path();
    let Ok(json) = serde_json::to_string_pretty(saved) else {
        return;
    };
    if let Some(dir) = path.parent() {
        let _ = std::fs::create_dir_all(dir);
    }
    if let Err(e) = std::fs::write(&path, json) {
        eprintln!("could not save layout to {}: {e}", path.display());
    }
}

/// A layout fresh from `default_layout()` (and friends), or one just
/// switched to (`layout NAME`), leaves every node at `Rect::NOTHING`
/// (infinite coordinates) until its first `DockArea` pass; serde_json
/// encodes non-finite floats as `null`, which an f32 field can't read
/// back. `commands.rs`'s `capture_view` reads `AppUi::dock_layout` -- a
/// snapshot taken at the very start of a frame, before that frame's own
/// pass -- so a session saved right after such a switch would otherwise
/// capture a `dock_layout` that fails to deserialize on load, silently
/// dropping it (`serde_json::from_value(..).ok()`). Gives every node a
/// real rect instead; where it puts things does not matter once loaded,
/// since egui_dock re-lays-out from the tree and its fractions anyway.
pub(crate) fn give_every_node_a_real_rect(state: &mut DockState<Tab>) {
    let rect = egui::Rect::from_min_size(egui::Pos2::ZERO, egui::vec2(800.0, 600.0));
    for node in state.main_surface_mut().iter_mut() {
        node.set_rect(rect);
        // `set_rect` only touches `LeafNode::rect`; `viewport` (the tab
        // body area) is otherwise only ever assigned by egui_dock's own
        // layout pass, so it has to be poked directly here too.
        if let egui_dock::Node::Leaf(leaf) = node {
            leaf.viewport = rect;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sorted_tabs(state: &DockState<Tab>) -> Vec<Tab> {
        let mut tabs: Vec<Tab> = state.iter_all_tabs().map(|(_, t)| *t).collect();
        tabs.sort_by_key(|t| *t as usize);
        tabs
    }

    #[test]
    fn default_layout_contains_every_tab_but_the_timeline_exactly_once() {
        assert_eq!(
            sorted_tabs(&default_layout()),
            [
                Tab::Viewport,
                Tab::Structures,
                Tab::Selection,
                Tab::Inspector,
                Tab::Log,
                Tab::Sequence,
                Tab::Info,
            ]
        );
    }

    #[test]
    fn default_layout_starts_log_and_sequence_collapsed() {
        let state = default_layout();
        let collapsed = |tab: Tab| {
            state
                .iter_all_nodes()
                .any(|(_, n)| n.tabs().is_some_and(|t| t.contains(&tab)) && n.is_collapsed())
        };
        assert!(collapsed(Tab::Log) && collapsed(Tab::Sequence));
        assert!(!collapsed(Tab::Viewport) && !collapsed(Tab::Structures));
    }

    #[test]
    fn the_viewport_stays_alone_open_and_docked() {
        for (name, state) in built_in_workspaces() {
            assert!(viewport_in_place(&state), "{name}");
        }
        let mut shared = default_layout();
        let (viewport, _) = shared.main_surface().find_tab(&Tab::Viewport).unwrap();
        shared.set_focused_node_and_surface(egui_dock::NodePath {
            surface: egui_dock::SurfaceIndex::main(),
            node: viewport,
        });
        shared.push_to_focused_leaf(Tab::Timeline);
        assert!(!viewport_in_place(&shared));
        assert!(!viewport_in_place(&DockState::new(vec![Tab::Log])));
    }

    #[test]
    fn trajectory_layout_adds_the_timeline_once() {
        let state = trajectory_layout();
        assert_eq!(
            state
                .iter_all_tabs()
                .filter(|(_, t)| **t == Tab::Timeline)
                .count(),
            1
        );
        // Opening it again is a no-op.
        let mut again = state.clone();
        open_near_viewport(&mut again, Tab::Timeline);
        assert_eq!(again.iter_all_tabs().count(), state.iter_all_tabs().count());
    }

    #[test]
    fn open_near_viewport_falls_back_when_there_is_no_viewport() {
        let mut state = DockState::new(vec![Tab::Log]);
        open_near_viewport(&mut state, Tab::Timeline);
        assert!(state.find_tab(&Tab::Timeline).is_some());
    }

    /// Collapsed panels stay collapsed across a save.
    #[test]
    fn collapsed_panels_survive_a_save() {
        let mut state = default_layout();
        for (_, node) in state.iter_all_nodes_mut() {
            if let egui_dock::Node::Leaf(leaf) = node {
                if leaf.tabs.contains(&Tab::Log) {
                    leaf.collapsed = true;
                }
            }
        }
        give_every_node_a_real_rect(&mut state);
        let json = serde_json::to_string(&state).unwrap();
        let loaded: DockState<Tab> = serde_json::from_str(&json).unwrap();
        assert!(loaded
            .iter_all_nodes()
            .any(|(_, node)| node.is_leaf() && node.is_collapsed()));
    }

    /// An older layout still names the `SceneSettings` tab; it deserializes (the variant stays in the
    /// enum) and `strip_scene_settings` drops it, leaving everything else.
    #[test]
    fn a_pre_a3_layout_naming_scene_settings_loads_and_strips_it() {
        let mut state = default_layout();
        give_every_node_a_real_rect(&mut state);
        state.push_to_first_leaf(Tab::SceneSettings);
        let json = serde_json::to_string(&state).unwrap();
        let mut loaded: DockState<Tab> = serde_json::from_str(&json).unwrap();
        assert!(loaded.find_tab(&Tab::SceneSettings).is_some());
        strip_scene_settings(&mut loaded);
        assert!(loaded.find_tab(&Tab::SceneSettings).is_none());
        assert_eq!(sorted_tabs(&loaded), sorted_tabs(&default_layout()));
    }

    /// An older layout or session names the tab "Scene"; the `#[serde(alias = "Scene")]` on `Tab::Structures`
    /// deserializes it under its new name with no strip step needed.
    #[test]
    fn a_pre_b1_tab_naming_scene_loads_as_structures() {
        let loaded: Tab = serde_json::from_str("\"Scene\"").unwrap();
        assert_eq!(loaded, Tab::Structures);
    }

    #[test]
    fn layout_round_trips_through_json() {
        let mut state = default_layout();
        give_every_node_a_real_rect(&mut state);
        let json = serde_json::to_string(&state).unwrap();
        let restored: DockState<Tab> = serde_json::from_str(&json).unwrap();
        assert_eq!(sorted_tabs(&state), sorted_tabs(&restored));
    }

    #[test]
    fn saved_workspaces_round_trip_through_the_saved_format() {
        let mut current = default_layout();
        give_every_node_a_real_rect(&mut current);
        let mut compact = compact_layout();
        give_every_node_a_real_rect(&mut compact);
        let workspaces = vec![("My Layout".to_string(), compact)];

        let mut prefs = UiPrefs {
            theme: crate::theme::ThemeMode::Light,
            ribbon_collapsed: true,
            dont_confirm_quit: true,
            hide_start_card: true,
            ..Default::default()
        };
        prefs.add_recent(PathBuf::from("4HHB.cif"));
        let mut start = trajectory_layout();
        give_every_node_a_real_rect(&mut start);
        let saved = SavedLayout {
            launch: current,
            workspaces,
            prefs: prefs.clone(),
            start: Some(start),
        };
        let json = serde_json::to_string(&saved).unwrap();
        let restored: SavedLayout = serde_json::from_str(&json).unwrap();
        assert!(restored.start.unwrap().find_tab(&Tab::Timeline).is_some());
        assert_eq!(restored.prefs, prefs);

        assert_eq!(restored.workspaces.len(), 1);
        assert_eq!(restored.workspaces[0].0, "My Layout");
        assert_eq!(
            restored.workspaces[0].1.iter_all_tabs().count(),
            compact_layout().iter_all_tabs().count()
        );
    }

    #[test]
    fn a_launch_ignores_the_last_arrangement() {
        let json = serde_json::to_string(&SavedLayout {
            launch: compact_layout(),
            ..SavedLayout::default()
        })
        .unwrap();
        let restored: SavedLayout = serde_json::from_str(&json).unwrap();
        assert_eq!(
            sorted_tabs(&restored.launch),
            sorted_tabs(&default_layout())
        );
    }

    #[test]
    fn closed_panels_lists_what_is_missing_and_never_the_viewport() {
        assert_eq!(
            closed_panels(&default_layout()),
            [Tab::Timeline, Tab::Movie]
        );
        let closed = closed_panels(&compact_layout());
        assert!(!closed.contains(&Tab::Viewport) && !closed.contains(&Tab::Inspector));
        assert_eq!(closed.len(), Tab::PANELS.len() - 1);
    }

    #[test]
    fn built_in_workspaces_are_named_and_each_has_a_viewport() {
        let workspaces = built_in_workspaces();
        for name in ["Default", "Trajectory", "Analysis", "Compact"] {
            assert!(workspaces.iter().any(|(n, _)| *n == name), "{name}");
        }
        for (name, state) in &workspaces {
            assert!(state.find_tab(&Tab::Viewport).is_some(), "{name}");
        }
    }
}

#[cfg(test)]
mod config_dir_tests {
    use super::*;

    #[test]
    fn config_dir_is_per_user_and_named_vizviz() {
        // On every CI runner and developer machine a home is set; the
        // fallback branch is only for stripped-down containers.
        let dir = config_dir().expect("a home directory");
        assert!(dir.is_absolute(), "{}", dir.display());
        assert_eq!(dir.file_name().unwrap(), "vizviz");
        assert!(layout_path().starts_with(&dir));
    }
}
