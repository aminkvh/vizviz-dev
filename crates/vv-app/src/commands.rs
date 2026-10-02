//! The app's command line: `vv_scene::script`'s document verbs plus the
//! view verbs only a running window has (style, camera, screenshot,
//! quit). One runner serves `--exec`, the command palette, and the Log
//! console, so anything you can click you can also script.
//!
//! View verbs are not commands on the bus and are not undoable, on
//! purpose: a style preset or a camera move is how you look at the
//! document, not part of it (see docs/RENDERING.md).

use egui_dock::DockState;
use vv_render::{Camera, LightingPreset, Projection, StylePreset};
use vv_scene::script::{self, split_verb};
use vv_scene::StructureId;

use crate::layout::Tab;
use crate::ui::{AppUi, ExportRequest, LayoutRequest, TimelineState, ViewSettings};

/// A palette/`help` entry, from either registry.
#[derive(Clone, Copy, Debug)]
pub struct Entry {
    pub id: &'static str,
    pub title: &'static str,
    pub keywords: &'static [&'static str],
    pub usage: &'static str,
    pub help: &'static str,
}

impl Entry {
    /// Whether the command needs something typed after it (its usage has
    /// a placeholder that isn't optional), so the palette should prefill
    /// the verb instead of running it bare.
    pub fn needs_args(&self) -> bool {
        self.usage
            .split_whitespace()
            .nth(1)
            .is_some_and(|arg| !arg.starts_with('['))
    }
}

/// View verbs. Document verbs come from `vv_scene::script::SPECS`.
pub const VIEW_ENTRIES: &[Entry] = &[
    Entry {
        id: "open",
        title: "Open a structure…",
        keywords: &["file", "dialog", "browse", "load"],
        usage: "open [PATH]",
        help: "Load a file and frame it. With no PATH, shows the file picker.",
    },
    Entry {
        id: "style",
        title: "Set style preset",
        keywords: &[
            "preset", "dark", "white", "publication", "glossy", "flat", "cel", "look", "lighting",
        ],
        usage: "style dark|white|glossy|flat",
        help: "A look: sets the lights. Effects (Preferences ▸ Effects), materials \
               (`material`) and the background (`background`) are each their own choice, \
               left as they are.",
    },
    Entry {
        id: "lighting",
        title: "Set lighting",
        keywords: &[
            "light",
            "lights",
            "soft",
            "full",
            "flat",
            "illustrative",
            "gentle",
            "twin",
            "twinlight",
            "studio",
            "three-point",
            "threepoint",
            "ambient",
            "key",
        ],
        usage: "lighting default|soft|full|flat|illustrative|gentle|twinlight|threepoint",
        help: "Global lighting preset. Effects (Preferences ▸ Effects) and materials \
               (`material`, per structure) are left as they are.",
    },
    Entry {
        id: "tonemap",
        title: "Toggle tonemapping",
        keywords: &["aces", "hdr", "highlights", "filmic"],
        usage: "tonemap [on|off]",
        help: "ACES tonemapping of bright highlights.",
    },
    Entry {
        id: "background",
        title: "Set background",
        keywords: &["color", "colour", "bg", "look", "lighting"],
        usage: "background R G B | #RRGGBB | default",
        help: "Viewport background as 0-255 components, a hex color, or the style's default.",
    },
    Entry {
        id: "mode",
        title: "Mouse mode",
        keywords: &["rotate", "translate", "scale", "center", "label", "drag", "tool"],
        usage: "mode rotate|translate|scale|center|label",
        help: "What a left drag in the viewport does, with the R T S C keys: rotate, \
               translate, scale, or center (click an atom to rotate about it). `label` \
               (Analyze > Label) labels each atom clicked until Escape or another mode.",
    },
    Entry {
        id: "selectmode",
        title: "Selection tool",
        keywords: &["select", "box", "circle", "lasso", "atom", "residue", "chain", "molecule", "drag", "home", "tool", "level"],
        usage: "selectmode shape click|box|circle|lasso|toggle|next | level atom|residue|chain|molecule|next | flyout shapes|more|interface",
        help: "Home > Select: what a left drag in the viewport selects. `shape` draws a box, \
               circle or freehand lasso and selects the atoms of the visible reps inside it \
               on release (Shift adds, Ctrl or Alt subtracts); `click` goes back to \
               orbiting, `toggle` switches between `click` and the last shape used, and \
               `next` (the Q key) cycles the shapes. `level` grows the hits to their whole \
               residues, chains or molecules (connected atoms); `level next` cycles it. \
               `flyout shapes`, `flyout more` and `flyout interface` open the tool's shape \
               flyout, the By type \"More…\" flyout and the Interface form.",
    },
    Entry {
        id: "ribbon",
        title: "Ribbon tab",
        keywords: &["tab", "toolbar", "collapse", "expand", "menu"],
        usage: "ribbon TAB|collapse|expand",
        help: "Open a ribbon tab by name (home, select, represent, ...), or fold the ribbon to its tab names and back.",
    },
    Entry {
        id: "panel",
        title: "Show panel",
        keywords: &["window", "tab", "dock", "open", "inspector", "log", "timeline", "look", "settings"],
        usage: "panel viewport|scene|selection|inspector|sequence|log|info|timeline|movie [collapse|expand]",
        help: "Show a dock panel (reopening it if it was closed) and bring it to the front. \
               `collapse` folds a panel to its tab bar and `expand` opens it again. \
               `settings`/`look` open the Look ribbon tab instead (its panel is gone).",
    },
    Entry {
        id: "play",
        title: "Play trajectory",
        keywords: &["timeline", "animate", "pause", "frames"],
        usage: "play [on|off|FPS]",
        help: "Play or pause the frames of the current multi-frame structure; a number plays at that many frames per second.",
    },
    Entry {
        id: "deletelayout",
        title: "Delete workspace",
        keywords: &["layout", "remove", "workspace"],
        usage: "deletelayout [NAME]",
        help: "Forget a saved workspace (the built-in ones stay); alone, lists them to \
               choose from.",
    },
    Entry {
        id: "startlayout",
        title: "Start layout",
        keywords: &["default", "layout", "workspace", "launch", "startup"],
        usage: "startlayout [reset]",
        help: "Open the current panel arrangement at every launch; `reset` goes back to \
               the built-in default.",
    },
    Entry {
        id: "preferences",
        title: "Preferences",
        keywords: &[
            "settings",
            "prefs",
            "options",
            "theme",
            "performance",
            "reset",
        ],
        usage: "preferences",
        help: "This computer's settings: theme and startup behavior, performance, and \
               resetting them.",
    },
    Entry {
        id: "look",
        title: "Open a look popover",
        keywords: &["lights", "background", "ao", "shadows", "depthcue", "dof", "effects"],
        usage: "look lights|background|ao|shadows|depthcue|dof",
        help: "Switches to the Look tab and opens one of its ▾ popovers, for scripting a \
               screenshot of it or jumping straight to it.",
    },
    Entry {
        id: "theme",
        title: "Set theme",
        keywords: &["dark", "light", "colors", "appearance", "mode", "preferences", "view"],
        usage: "theme dark|light",
        help: "The interface's colors (the viewport's background is `background`). Saved \
               with the layout.",
    },
    Entry {
        id: "welcome",
        title: "Show the start card",
        keywords: &["start", "welcome", "onboarding", "recent", "splash"],
        usage: "welcome [on|off]",
        help: "Shows the viewport's start card now (Open, Open session, Fetch, recent files), \
               closed for the rest of the run once dismissed. `on`/`off` set whether it offers \
               itself automatically over an empty viewport (the default: on, unless a file was \
               named on the command line); its own \"Don't show again\" checkbox is `welcome \
               off`.",
    },
    Entry {
        id: "sequence",
        title: "Sequence coloring and tracks",
        keywords: &[
            "annotation", "annotations", "kabat", "numbering", "glycosylation", "disulfide",
            "liability", "liabilities", "sasa", "hydrophobicity", "track", "legend", "strip",
        ],
        usage: "sequence color SCHEME | track NAME [on|off] | tracks [all|none] | legend [on|off] | antibody scheme|cdr|backend NAME | antibody exe PATH | uniprot [variants] [on|off] | props",
        help: "The Sequence panel's header as commands. `color` picks how residue letters \
               are colored: none, view (as the 3D view), ss, chemistry, hydrophobicity, \
               bfactor, sasa, charge, clustal, zappo, taylor. `track` shows or hides an \
               annotation row under each chain: ss, numbering, missing, disulfide, glycan, \
               liability, ligand, interface, altloc, modified, antibody, burial, conservation, \
               uniprot. `tracks all|none` sets every offline track; `legend` shows the key of \
               the enabled ones. `antibody scheme` picks the numbering (kabat, chothia, imgt, \
               martin) and `antibody cdr` the CDR definition (kabat, chothia, imgt, contact, \
               north), independently; `antibody backend anarci` takes the numbers from the \
               external ANARCI program instead of the built-in profiles (`antibody exe PATH` \
               says where it is), and `antibody backend abnum` from the scheme authors' public \
               web service (Kabat, Chothia and Martin only; sends sequences over plain HTTP; \
               `cdr` stays native). `uniprot on|off` draws UniProt features and `uniprot variants on|off` its natural variants (needs the \
               network). `props` prints each protein chain's mass, pI, charge and extinction.",
    },
    Entry {
        id: "confirmquit",
        title: "Confirm quit",
        keywords: &["quit", "confirm", "remind", "unsaved", "dialog"],
        usage: "confirmquit [on|off]",
        help: "Ask to save before quitting when there are unsaved changes (the default). \
               `off` is what the quit dialog's \"Don't remind me again\" checkbox sets. No \
               argument toggles.",
    },
    Entry {
        id: "uishot",
        title: "Capture window",
        keywords: &["screenshot", "window", "ui", "capture", "docs"],
        usage: "uishot PATH",
        help: "Write the whole window, panels and all, as a PNG (for docs and UI review; \
               `screenshot` writes the viewport alone).",
    },
    Entry {
        id: "window",
        title: "Window size",
        keywords: &["resize", "size", "width", "height", "layout"],
        usage: "window WxH",
        help: "Resize the window to W x H logical pixels (at least 1024 x 640), for checking the layout at a size.",
    },
    Entry {
        id: "dof",
        title: "Set depth of field",
        keywords: &["focus", "blur", "bokeh", "depth"],
        usage: "dof on|off|STRENGTH",
        help: "Depth of field, strength 0-1: sharp at the rotation centre (set it with C or a \
               double-click), blurring with distance from it.",
    },
    Entry {
        id: "rotate",
        title: "Rotate view",
        keywords: &["turn", "spin", "orbit", "angle"],
        usage: "rotate x|y|z DEGREES",
        help: "Turn the scene about a screen axis (x right, y up, z toward you).",
    },
    Entry {
        id: "clip",
        title: "Clip plane",
        keywords: &["cut", "section", "slab", "cross-section", "cutplane", "render", "surface"],
        usage: "clip on|off|view|x|y|z|flip|DEPTH",
        help: "Cut the scene with a plane through the rotation centre: every solid (a sphere, a \
               bond, a tube, a closed surface) shows a flat cap of its own colour where it is \
               cut, and the space between them stays open. `view` faces the camera and tracks \
               it as you rotate; `x`, `y`, `z` fix it to that world axis instead; `flip` keeps \
               the other side. DEPTH (Angstrom) moves the plane along its normal, away from \
               what it keeps.",
    },
    Entry {
        id: "shadows",
        title: "Set shadows",
        keywords: &["shadow", "light", "cast", "look", "lighting"],
        usage: "shadows on|off|STRENGTH",
        help: "Screen-space shadows from the key light, strength 0-1.",
    },
    Entry {
        id: "gradient",
        title: "Set background gradient",
        keywords: &["background", "fade", "vertical", "bg", "look", "lighting"],
        usage: "gradient R G B | #RRGGBB | on | off",
        help: "Fade the background from its color at the bottom to this one at the top; \
               `on` picks a lighter shade, `off` makes it solid.",
    },
    Entry {
        id: "fxaa",
        title: "FXAA on/off",
        keywords: &["antialias", "smooth", "edges", "preferences", "render"],
        usage: "fxaa [on|off]",
        help: "Edge smoothing in the viewport. No argument toggles.",
    },
    Entry {
        id: "outline",
        title: "Outlines on/off",
        keywords: &["edges", "cel", "toon", "lines", "silhouette", "look", "lighting"],
        usage: "outline [on|off]",
        help: "Dark lines on silhouettes and creases (the Flat/Cel preset turns them on). No argument toggles.",
    },
    Entry {
        id: "ao",
        title: "Ambient occlusion",
        keywords: &["ssao", "shadow", "depth", "material", "lighting", "occlusion", "look"],
        usage: "ao [on|off|STRENGTH]",
        help: "Darkens grooves and pockets so shape reads at any zoom. STRENGTH 0-2 (on = 1). No argument toggles.",
    },
    Entry {
        id: "depthcue",
        title: "Depth cue",
        keywords: &["fog", "depth", "haze", "material", "look", "lighting"],
        usage: "depthcue [on|off|STRENGTH]",
        help: "Fades toward the background with distance from the scene's front, most visible near its silhouette. STRENGTH 0-1 (on = 0.5). No argument toggles.",
    },
    Entry {
        id: "supersample",
        title: "Supersampling",
        keywords: &["antialias", "ssaa", "quality", "sharp", "render scale", "preferences", "render"],
        usage: "supersample 1|1.5|2",
        help: "Render the viewport at 1x, 1.5x or 2x its size and average down: sharper thin geometry, 4x the pixels at 2x.",
    },
    Entry {
        id: "occlusion",
        title: "Occlusion culling on/off",
        keywords: &["cull", "hidden", "performance", "preferences", "render"],
        usage: "occlusion [on|off]",
        help: "Skip atoms hidden behind the previous frame. No argument toggles.",
    },
    Entry {
        id: "adaptive",
        title: "Adaptive quality on/off",
        keywords: &["lod", "performance", "fps", "preferences", "render"],
        usage: "adaptive [on|off]",
        help: "Let the point/sphere threshold move to hold the frame rate. No argument toggles.",
    },
    Entry {
        id: "view",
        title: "Camera",
        keywords: &["camera", "reset", "orbit", "zoom", "dolly", "frame", "projection", "orthographic", "perspective"],
        usage: "view reset | orbit YAW PITCH | zoom FACTOR | dolly ANGSTROM | face +X|-X|+Y|-Y|+Z|-Z | projection perspective|orthographic",
        help: "Reset frames every structure (and undoes the dolly); orbit is in degrees; zoom > 1 moves away; dolly moves the camera itself forward (negative: back), cutting the scene flat where it passes, as Shift+wheel does; face looks from that side, as clicking the view cube; projection switches perspective/orthographic.",
    },
    Entry {
        id: "layout",
        title: "Switch workspace",
        keywords: &["workspace", "panels", "default", "trajectory", "analysis", "compact"],
        usage: "layout NAME",
        help: "Switch to a built-in (default, trajectory, analysis, compact) or saved workspace.",
    },
    Entry {
        id: "savelayout",
        title: "Save the workspace as…",
        keywords: &["workspace", "panels"],
        usage: "savelayout [NAME]",
        help: "Save the panel arrangement as a workspace. With no NAME, asks for one.",
    },
    Entry {
        id: "screenshot",
        title: "Save a screenshot…",
        keywords: &["png", "jpg", "jpeg", "svg", "vector", "export", "render", "image", "figure"],
        usage: "screenshot [PATH] [nossaa]",
        help: "Write the viewport as a PNG, a JPEG if PATH ends in .jpg/.jpeg (2x supersampled unless `nossaa`), or a vector SVG if PATH ends in .svg (atoms as painter's-algorithm-sorted circles, plus interaction dashes as lines; no bonds/cartoon yet). With no PATH, opens the dialog.",
    },
    Entry {
        id: "render",
        title: "Render a path-traced image…",
        keywords: &["path", "trace", "ray", "shadows", "occlusion", "figure", "4k", "publication", "quality", "draft", "ultra"],
        usage: "render PATH [WIDTHxHEIGHT] [SAMPLES] [draft|high|ultra] [transparent]",
        help: "Path-trace the view into a PNG (or JPEG): real shadows and ambient occlusion, antialiased, at any size (default the viewport's) with SAMPLES per pixel (default by quality). Quality (default high) sets the sample count and cartoon/tube tessellation together; draft is fast and coarse, ultra the slowest and smoothest. `transparent` leaves the background out. Spacefill, ball-and-stick, SAS, sticks, lines and cartoon reps; the others are named and left out. Runs in the background.",
    },
    Entry {
        id: "framing",
        title: "Render frame preview on/off",
        keywords: &["render", "safe area", "letterbox", "aspect", "crop", "studio"],
        usage: "framing [on|off]",
        help: "Overlay a safe-frame rectangle on the viewport at the render form's width:height, dimmed outside it, so framing matches the exported image. Also shows while File ▸ Export ▸ Render's popover is open. No argument toggles.",
    },
    Entry {
        id: "movie",
        title: "Movie maker",
        keywords: &["animation", "video", "tracks", "clips", "timeline", "mp4", "frames", "record"],
        usage: crate::movie_cmd::USAGE,
        help: crate::movie_cmd::HELP,
    },
    Entry {
        id: "quit",
        title: "Quit",
        keywords: &["exit", "close"],
        usage: "quit",
        help: "Exit the application (saves the panel layout).",
    },
];

fn from_spec(s: &script::Spec) -> Entry {
    Entry {
        id: s.id,
        title: s.title,
        keywords: s.keywords,
        usage: s.usage,
        help: s.help,
    }
}

/// Every command the app understands: document verbs first, then view
/// verbs. This list is the palette, `help`, and the CLI's validator.
pub fn entries() -> Vec<Entry> {
    script::SPECS
        .iter()
        .map(from_spec)
        .chain(VIEW_ENTRIES.iter().copied())
        .collect()
}

/// What a command line's command does, for tooltips.
pub fn help_for(line: &str) -> Option<&'static str> {
    let (verb, _) = split_verb(line.trim());
    entries().into_iter().find(|e| e.id == verb).map(|e| e.help)
}

fn view_entry(verb: &str) -> Option<&'static Entry> {
    let id = match verb {
        "exit" => "quit",
        "bg" => "background",
        "camera" => "view",
        other => other,
    };
    VIEW_ENTRIES.iter().find(|e| e.id == id)
}

/// Checks that a command line names a known verb, without running it.
/// Lets `--exec` typos fail before a window opens.
pub fn validate(line: &str) -> Result<(), String> {
    let (verb, _) = split_verb(line);
    if verb.is_empty() || script::spec(verb).is_some() || view_entry(verb).is_some() {
        Ok(())
    } else {
        Err(format!(
            "unknown command `{verb}`; try `vizviz --exec help`"
        ))
    }
}

/// `validate`, and refuses the forms that open a dialog (a command with
/// no path): an unattended script (`--exec`, a `--listen` client) has
/// nobody to answer it.
pub fn validate_unattended(line: &str) -> Result<(), String> {
    validate(line)?;
    let (verb, rest) = split_verb(line);
    let asks = [
        "open",
        "fetch",
        "savesession",
        "loadsession",
        "savestructure",
        "screenshot",
        "savelayout",
        "deletelayout",
    ];
    if asks.contains(&verb) && rest.trim().is_empty() {
        return Err(format!(
            "`{verb}` needs its argument in a script (alone it opens a dialog)"
        ));
    }
    Ok(())
}

/// `on`/`off`/a number for a 0..`max` strength; no argument toggles
/// between off and `on`.
fn parse_strength(word: &str, current: f32, on: f32, max: f32) -> Result<f32, String> {
    match word {
        "" | "toggle" => Ok(if current > 0.0 { 0.0 } else { on }),
        "on" | "true" | "yes" => Ok(on),
        "off" | "false" | "no" => Ok(0.0),
        other => other
            .parse::<f32>()
            .ok()
            .filter(|v| v.is_finite())
            .map(|v| v.clamp(0.0, max))
            .ok_or_else(|| format!("expected on, off or a number, got `{other}`")),
    }
}

pub(crate) fn parse_on_off(word: &str, current: bool) -> Result<bool, String> {
    match word {
        "" | "toggle" => Ok(!current),
        "on" | "true" | "1" | "yes" => Ok(true),
        "off" | "false" | "0" | "no" => Ok(false),
        other => Err(format!("expected on or off, got `{other}`")),
    }
}

/// The word `style` accepts for a preset; what a session file stores.
pub(crate) fn style_name(style: StylePreset) -> &'static str {
    match style {
        StylePreset::DarkPresentation => "dark",
        StylePreset::PublicationWhite => "white",
        StylePreset::Glossy => "glossy",
        StylePreset::FlatCel => "flat",
    }
}

fn parse_style(word: &str) -> Result<StylePreset, String> {
    Ok(match word {
        "dark" | "dark_presentation" | "presentation" => StylePreset::DarkPresentation,
        "white" | "publication" | "publication_white" => StylePreset::PublicationWhite,
        "glossy" => StylePreset::Glossy,
        "flat" | "cel" | "flat_cel" | "toon" => StylePreset::FlatCel,
        other => {
            return Err(format!(
                "unknown style `{other}`; expected dark, white, glossy, or flat"
            ))
        }
    })
}

fn projection_name(projection: Projection) -> &'static str {
    match projection {
        Projection::Perspective => "perspective",
        Projection::Orthographic => "orthographic",
    }
}

fn parse_projection(word: &str) -> Result<Projection, String> {
    Ok(match word {
        "perspective" | "persp" => Projection::Perspective,
        "orthographic" | "ortho" => Projection::Orthographic,
        other => {
            return Err(format!(
                "unknown projection `{other}`; expected perspective or orthographic"
            ))
        }
    })
}

/// A picked colour as `[r, g, b]` 0..1, matching how `background`/
/// `background_top` already store colour in a session.
fn color32_floats(c: egui::Color32) -> [f32; 3] {
    [c.r(), c.g(), c.b()].map(|v| v as f32 / 255.0)
}

fn parse_color32_floats(v: &serde_json::Value) -> Option<egui::Color32> {
    let n: Vec<f32> = v
        .as_array()?
        .iter()
        .filter_map(|v| v.as_f64())
        .map(|v| v as f32)
        .collect();
    let [r, g, b] = n[..] else { return None };
    let to_u8 = |c: f32| (c.clamp(0.0, 1.0) * 255.0).round() as u8;
    Some(egui::Color32::from_rgb(to_u8(r), to_u8(g), to_u8(b)))
}

fn saved_view_snapshot(v: &crate::studio::SavedView) -> serde_json::Value {
    serde_json::json!({
        "name": v.name,
        "target": [v.target.x, v.target.y, v.target.z],
        "distance": v.distance,
        "orientation": v.orientation.to_array(),
        "fov_y": v.fov_y, "near": v.near,
        "projection": projection_name(v.projection),
        "scene_radius": v.scene_radius, "dolly": v.dolly,
    })
}

/// Reads a `saved_view_snapshot`; `None` skips a malformed entry rather
/// than failing the whole session load.
fn parse_saved_view(v: &serde_json::Value) -> Option<crate::studio::SavedView> {
    let name = v["name"].as_str()?.to_owned();
    let f = |x: &serde_json::Value| x.as_f64().map(|n| n as f32);
    let t: Vec<f32> = v["target"].as_array()?.iter().filter_map(f).collect();
    let [tx, ty, tz] = t[..] else { return None };
    let q: Vec<f32> = v["orientation"].as_array()?.iter().filter_map(f).collect();
    let [qx, qy, qz, qw] = q[..] else { return None };
    Some(crate::studio::SavedView {
        name,
        target: glam::Vec3::new(tx, ty, tz),
        distance: f(&v["distance"])?,
        orientation: glam::Quat::from_xyzw(qx, qy, qz, qw).normalize(),
        fov_y: f(&v["fov_y"])?,
        near: f(&v["near"])?,
        projection: v["projection"]
            .as_str()
            .and_then(|s| parse_projection(s).ok())
            .unwrap_or_default(),
        scene_radius: f(&v["scene_radius"]).unwrap_or(1.0),
        dolly: f(&v["dolly"]).unwrap_or(0.0).max(0.0),
    })
}

fn clip_snapshot(plane: crate::ui::ClipPlane) -> serde_json::Value {
    use crate::ui::ClipOrientation;
    let orientation = match plane.orientation {
        ClipOrientation::View => serde_json::json!("view"),
        ClipOrientation::X => serde_json::json!("x"),
        ClipOrientation::Y => serde_json::json!("y"),
        ClipOrientation::Z => serde_json::json!("z"),
        ClipOrientation::Fixed(v) => serde_json::json!({ "fixed": [v.x, v.y, v.z] }),
    };
    serde_json::json!({
        "orientation": orientation,
        "depth": plane.depth,
        "flipped": plane.flipped,
    })
}

/// Reads a `clip_snapshot`, or an older session (a raw
/// `[x, y, z, d]` world-space plane; converted to a fixed orientation, its
/// depth measured from the given rotation centre).
fn parse_clip_snapshot(v: &serde_json::Value, target: glam::Vec3) -> Option<crate::ui::ClipPlane> {
    use crate::ui::{ClipOrientation, ClipPlane};
    if let Some(plane) = v.as_array() {
        let n: Vec<f32> = plane
            .iter()
            .filter_map(|v| v.as_f64())
            .map(|v| v as f32)
            .collect();
        let [x, y, z, d] = n[..] else { return None };
        let normal = glam::Vec3::new(x, y, z);
        return Some(ClipPlane {
            orientation: ClipOrientation::Fixed(normal),
            depth: -d - normal.dot(target),
            flipped: false,
        });
    }
    let orientation = match &v["orientation"] {
        serde_json::Value::String(s) => match s.as_str() {
            "view" => ClipOrientation::View,
            "x" => ClipOrientation::X,
            "y" => ClipOrientation::Y,
            "z" => ClipOrientation::Z,
            _ => return None,
        },
        serde_json::Value::Object(_) => {
            let f = v["orientation"]["fixed"].as_array()?;
            let n: Vec<f32> = f
                .iter()
                .filter_map(|v| v.as_f64())
                .map(|v| v as f32)
                .collect();
            let [x, y, z] = n[..] else { return None };
            ClipOrientation::Fixed(glam::Vec3::new(x, y, z))
        }
        _ => return None,
    };
    Some(ClipPlane {
        orientation,
        depth: v["depth"].as_f64()? as f32,
        flipped: v["flipped"].as_bool().unwrap_or(false),
    })
}

fn clip_description(plane: crate::ui::ClipPlane) -> String {
    use crate::ui::ClipOrientation;
    let facing = match plane.orientation {
        ClipOrientation::View => "view",
        ClipOrientation::X => "x",
        ClipOrientation::Y => "y",
        ClipOrientation::Z => "z",
        ClipOrientation::Fixed(_) => "captured view",
    };
    let flipped = if plane.flipped { ", flipped" } else { "" };
    format!("clip at {:.1} A, {facing}{flipped}", plane.depth)
}

/// `R G B [A]` in 0-255, or `#RRGGBB`/`#RRGGBBAA`, as a linear-space
/// `wgpu::Color` (sRGB bytes decoded the same way the style presets'
/// hand-picked backgrounds are expressed).
fn parse_color(text: &str) -> Result<wgpu::Color, String> {
    let bad = || format!("expected `R G B`, `#RRGGBB`, or `default`, got `{text}`");
    let bytes: Vec<u8> = if let Some(hex) = text.strip_prefix('#') {
        if hex.len() != 6 && hex.len() != 8 {
            return Err(bad());
        }
        (0..hex.len())
            .step_by(2)
            .map(|i| u8::from_str_radix(&hex[i..i + 2], 16).map_err(|_| bad()))
            .collect::<Result<_, _>>()?
    } else {
        let parts: Vec<u8> = text
            .split_whitespace()
            .map(|p| p.parse::<u8>().map_err(|_| bad()))
            .collect::<Result<_, _>>()?;
        if parts.len() != 3 && parts.len() != 4 {
            return Err(bad());
        }
        parts
    };
    let linear = |b: u8| {
        let c = b as f64 / 255.0;
        if c <= 0.04045 {
            c / 12.92
        } else {
            ((c + 0.055) / 1.055).powf(2.4)
        }
    };
    Ok(wgpu::Color {
        r: linear(bytes[0]),
        g: linear(bytes[1]),
        b: linear(bytes[2]),
        a: bytes.get(3).map_or(1.0, |&a| a as f64 / 255.0),
    })
}

fn parse_f32(word: &str, what: &str) -> Result<f32, String> {
    word.parse::<f32>()
        .ok()
        .filter(|v| v.is_finite())
        .ok_or_else(|| format!("expected a number for {what}, got `{word}`"))
}

/// Everything the app's `view` half of a session file captures: render
/// settings, the camera, the dock layout at save time, and playback --
/// opaque JSON `vv_scene::session` carries without looking inside (see
/// that module's doc). A free function, not an `AppUi` method, so it
/// (and `restore_view`) can be unit-tested without a GPU-backed `State`.
fn capture_view(
    view: &ViewSettings,
    camera: &Camera,
    dock_layout: &DockState<Tab>,
    timeline: &TimelineState,
    movie: &crate::movie::Movie,
) -> serde_json::Value {
    // A layout switched to only shortly before saving may not have been
    // through a `DockArea` pass yet, leaving non-finite rects that
    // `restore_view` could never deserialize back — see
    // `layout::give_every_node_a_real_rect`'s doc.
    let mut dock_layout = dock_layout.clone();
    crate::layout::give_every_node_a_real_rect(&mut dock_layout);
    serde_json::json!({
        "style": style_name(view.style),
        "lighting": view.lighting.name(),
        "lights": {
            "count": view.lights.count,
            "ambient": view.lights.ambient,
            "ambient_color": color32_floats(view.lights.ambient_color),
            // Every slot, not just the `count` in the list: a light
            // dropped by "-" keeps its settings if "+" brings it back.
            "list": view.lights.lights.iter().map(|l| serde_json::json!({
                "enabled": l.enabled,
                "color": color32_floats(l.color),
                "intensity": l.intensity,
                "azimuth": l.azimuth,
                "elevation": l.elevation,
            })).collect::<Vec<_>>(),
        },
        "saved_views": view.saved_views.iter().map(saved_view_snapshot).collect::<Vec<_>>(),
        // Without this, a reloaded session with hand-tuned lights shows the
        // last preset name in Look ▸ Lights ▾ instead of "Custom".
        "lights_custom": view.lights_custom,
        "tonemap": view.tonemap,
        "background": [view.background.r, view.background.g,
                       view.background.b, view.background.a],
        "background_top": view.background_top.map(|c| [c.r, c.g, c.b, c.a]),
        "fxaa": view.fxaa,
        "outline": view.outline,
        "ao": view.ao,
        "depth_cue": view.depth_cue,
        "shadows": view.shadows,
        "dof": view.dof,
        "clip": view.clip.map(clip_snapshot),
        // The `supersample` command's setting; keeps its original field
        // name for on-disk compatibility with session files saved
        // before that command existed.
        "render_scale": view.render_scale,
        "occlusion_culling": view.occlusion_culling,
        "adaptive": view.adaptive,
        "camera": {
            "target": [camera.target.x, camera.target.y, camera.target.z],
            "distance": camera.distance, "orientation": camera.orientation.to_array(),
            "fov_y": camera.fov_y, "near": camera.near,
            "projection": projection_name(camera.projection),
            "scene_radius": camera.scene_radius, "dolly": camera.dolly,
        },
        // Opaque to vv-scene: a schema-3 file simply lacks these two, and
        // a load leaves the dock layout and playback as they already are.
        "dock_layout": serde_json::to_value(&dock_layout).ok(),
        "playback": {
            "playing": timeline.playing,
            "fps": timeline.fps,
            "looping": timeline.looping,
        },
        "movie": movie,
    })
}

/// Applies a `view` block written by `capture_view`. Missing or
/// malformed fields are skipped one by one, so an older or hand-edited
/// session still restores what it can. Returns the saved dock layout, if
/// any and well-formed: the caller threads it through `LayoutRequest`
/// (`dock_state` is borrowed elsewhere while panels draw, so this
/// function can't set it directly -- see `AppUi::dock_layout`'s doc).
fn restore_view(
    view: &serde_json::Value,
    settings: &mut ViewSettings,
    camera: &mut Camera,
    timeline: &mut TimelineState,
    movie: &mut crate::movie::Movie,
) -> Option<DockState<Tab>> {
    if let Some(style) = view["style"].as_str().and_then(|s| parse_style(s).ok()) {
        settings.style = style;
        settings.lighting = style.lighting_preset();
        settings.tonemap = style.lighting().tonemap > 0.5;
    }
    if let Some(lighting) = view["lighting"].as_str().and_then(LightingPreset::parse) {
        settings.lighting = lighting;
        settings.lights = crate::studio::LightRig::of(lighting);
    }
    match view["lights"]["list"].as_array() {
        Some(list) => {
            for (i, entry) in list.iter().take(vv_render::style::MAX_LIGHTS).enumerate() {
                let light = &mut settings.lights.lights[i];
                if let Some(enabled) = entry["enabled"].as_bool() {
                    light.enabled = enabled;
                }
                if let Some(color) = parse_color32_floats(&entry["color"]) {
                    light.color = color;
                }
                if let Some(v) = entry["intensity"].as_f64() {
                    light.intensity = v as f32;
                }
                if let Some(v) = entry["azimuth"].as_f64() {
                    light.azimuth = v as f32;
                }
                if let Some(v) = entry["elevation"].as_f64() {
                    light.elevation = v as f32;
                }
            }
            if let Some(count) = view["lights"]["count"].as_u64() {
                settings.lights.count = (count as usize).min(vv_render::style::MAX_LIGHTS);
            }
        }
        // A session saved before the light list existed: one scalar
        // key/fill pair, mapped onto lights 0/1 (still enabled, still
        // in the list, as the preset just set them).
        None => {
            if let Some(v) = view["lights"]["key"].as_f64() {
                settings.lights.lights[0].intensity = v as f32;
            }
            if let Some(v) = view["lights"]["key_azimuth"].as_f64() {
                settings.lights.lights[0].azimuth = v as f32;
            }
            if let Some(v) = view["lights"]["key_elevation"].as_f64() {
                settings.lights.lights[0].elevation = v as f32;
            }
            if let Some(v) = view["lights"]["fill"].as_f64() {
                settings.lights.lights[1].intensity = v as f32;
            }
        }
    }
    if let Some(v) = view["lights"]["ambient"].as_f64() {
        settings.lights.ambient = v as f32;
    }
    if let Some(c) = parse_color32_floats(&view["lights"]["ambient_color"]) {
        settings.lights.ambient_color = c;
    }
    if let Some(list) = view["saved_views"].as_array() {
        settings.saved_views = list.iter().filter_map(parse_saved_view).collect();
    }
    if let Some(tonemap) = view["tonemap"].as_bool() {
        settings.tonemap = tonemap;
    }
    if let Some(custom) = view["lights_custom"].as_bool() {
        settings.lights_custom = custom;
    }
    if let Some(top) = view["background_top"].as_array() {
        let n: Vec<f64> = top.iter().filter_map(|v| v.as_f64()).collect();
        if let [r, g, b, a] = n[..] {
            settings.background_top = Some(wgpu::Color { r, g, b, a });
        }
    }
    if let Some(bg) = view["background"].as_array() {
        let n: Vec<f64> = bg.iter().filter_map(|v| v.as_f64()).collect();
        if let [r, g, b, a] = n[..] {
            settings.background = wgpu::Color { r, g, b, a };
        }
    }
    for (key, field) in [
        ("fxaa", &mut settings.fxaa),
        ("outline", &mut settings.outline),
        ("occlusion_culling", &mut settings.occlusion_culling),
        ("adaptive", &mut settings.adaptive),
    ] {
        if let Some(b) = view[key].as_bool() {
            *field = b;
        }
    }
    if let Some(v) = view["ao"].as_f64() {
        settings.ao = (v as f32).clamp(0.0, 2.0);
    }
    if let Some(v) = view["dof"].as_f64() {
        settings.dof = (v as f32).clamp(0.0, 1.0);
    }
    if let Some(v) = view["shadows"].as_f64() {
        settings.shadows = v as f32;
    }
    if let Some(v) = view["depth_cue"].as_f64() {
        settings.depth_cue = (v as f32).clamp(0.0, 1.0);
    }
    if let Some(scale) = view["render_scale"].as_f64() {
        if (1.0..=2.0).contains(&scale) {
            settings.render_scale = scale as f32;
        }
    }
    let cam = &view["camera"];
    let f = |v: &serde_json::Value| v.as_f64().map(|x| x as f32);
    // Older sessions stored yaw/pitch instead.
    let orientation = match cam["orientation"].as_array() {
        Some(q) => match q.iter().filter_map(f).collect::<Vec<_>>()[..] {
            [x, y, z, w] => Some(glam::Quat::from_xyzw(x, y, z, w).normalize()),
            _ => None,
        },
        None => f(&cam["yaw"])
            .zip(f(&cam["pitch"]))
            .map(|(yaw, pitch)| Camera::angles(yaw, pitch)),
    };
    if let (Some(t), Some(distance), Some(orientation), Some(fov_y), Some(near)) = (
        cam["target"].as_array(),
        f(&cam["distance"]),
        orientation,
        f(&cam["fov_y"]),
        f(&cam["near"]),
    ) {
        let t: Vec<f32> = t.iter().filter_map(f).collect();
        if let [x, y, z] = t[..] {
            // Missing/unrecognized `projection` (older session files
            // never wrote it) keeps the projection already in effect,
            // rather than silently resetting to perspective.
            let projection = cam["projection"]
                .as_str()
                .and_then(|s| parse_projection(s).ok())
                .unwrap_or(camera.projection);
            *camera = Camera {
                target: glam::Vec3::new(x, y, z),
                distance,
                orientation,
                fov_y,
                near,
                projection,
                scene_radius: f(&cam["scene_radius"]).unwrap_or(camera.scene_radius),
                pivot: None,
                dolly: f(&cam["dolly"]).unwrap_or(0.0).max(0.0),
            };
        }
    }
    // After the camera, so an old raw-plane session's depth is measured
    // from the target it actually had.
    if view["clip"].is_object() || view["clip"].is_array() {
        settings.clip = parse_clip_snapshot(&view["clip"], camera.target);
    }
    if let Some(playing) = view["playback"]["playing"].as_bool() {
        timeline.playing = playing;
    }
    if let Some(fps) = view["playback"]["fps"].as_f64() {
        if fps > 0.0 && fps <= 240.0 {
            timeline.fps = fps as f32;
        }
    }
    if let Some(looping) = view["playback"]["looping"].as_bool() {
        timeline.looping = looping;
    }
    *movie = serde_json::from_value(view["movie"].clone()).unwrap_or_default();
    let mut dock: DockState<Tab> = serde_json::from_value(view["dock_layout"].clone()).ok()?;
    crate::layout::strip_scene_settings(&mut dock);
    Some(dock)
}

impl AppUi<'_> {
    /// Switches to ribbon tab `tab` and asks it to force one of its ▾
    /// popovers open next time it draws (`ribbon::RibbonState::
    /// pending_popover`) -- `clip`/`layout`/`panel`/`look`'s
    /// no-argument forms, for scripting a screenshot of a popover or
    /// jumping straight to it.
    fn open_popover(&mut self, tab: &str, key: &'static str) {
        self.ribbon.tab = crate::ribbon::tab_named(tab).unwrap_or(self.ribbon.tab);
        self.prefs.ribbon_collapsed = false;
        self.ribbon.pending_popover = Some(key);
    }

    /// Runs one command line. `Ok` is a short result for the log.
    pub fn run_command(&mut self, line: &str) -> Result<String, String> {
        let (verb, rest) = split_verb(line);
        if verb.is_empty() {
            return Ok(String::new());
        }
        if let Some(spec) = script::spec(verb) {
            match spec.id {
                "help" => return Ok(self.help(rest)),
                // Undo/redo cross scene edits and view/camera edits
                // alike now (`app_history`), so these no longer go
                // straight to `script::run_line`/`vv_scene::
                // CommandHistory` -- `app_undo`/`app_redo` call into it
                // themselves for a `Scene` marker.
                "undo" => return Ok(self.app_undo()),
                "redo" => return Ok(self.app_redo()),
                // Sessions carry the app's view state too, so the app
                // handles these two itself (with dialogs when no path).
                "fetch" if rest.is_empty() => {
                    self.fetch_dialog.open = true;
                    return Ok(String::new());
                }
                "structures" if rest == "menu" => {
                    self.ribbon.open_structure_menu = true;
                    return Ok(String::new());
                }
                "loadtrajectory" if rest.is_empty() => {
                    self.open_popover("file", "file.trajectory");
                    return Ok(String::new());
                }
                // Represent ▸ Add rep's form (style + "Draws").
                "addrep" if rest.is_empty() => {
                    self.open_popover("represent", "represent.addrep");
                    return Ok(String::new());
                }
                // Represent ▸ Coloring's popover.
                "color" if rest.is_empty() => {
                    self.open_popover("represent", "represent.coloring");
                    return Ok(String::new());
                }
                // Analyze ▸ Caption's/Contacts' forms.
                "caption" if rest.is_empty() => {
                    self.open_popover("analyze", "analyze.caption");
                    return Ok(String::new());
                }
                "contacts" if rest.is_empty() => {
                    self.open_popover("analyze", "analyze.contacts");
                    return Ok(String::new());
                }
                // File ▸ Export ▸ Structure… and a structure row's ⋯
                // Export both run this bare; the dialog then writes with
                // `vv_io::save` directly (`export_structure_dialog_ui`),
                // not by re-running `savestructure` itself.
                "savestructure" if rest.is_empty() => {
                    self.export_structure_dialog.structure = None;
                    self.export_structure_dialog.open = true;
                    return Ok(String::new());
                }
                "savesession" => return self.save_session(rest),
                "loadsession" => return self.load_session(rest),
                _ => {}
            }
            let out = script::run_line(self.scene, self.history, line).map_err(|e| e.0)?;
            if spec.id == "load" || spec.id == "fetch" || spec.id == "loadtrajectory" {
                self.reset_camera();
                // A row's "File info" can point Info at a structure
                // other than the current one; a fresh load should still
                // follow the current one again, as it always used to.
                *self.info_target = None;
            }
            return Ok(out);
        }
        let Some(entry) = view_entry(verb) else {
            return Err(format!("unknown command `{verb}`; try `help`"));
        };
        let usage = || format!("usage: {}", entry.usage);
        match entry.id {
            "open" => {
                if rest.is_empty() {
                    self.dialog.open = true;
                    Ok(String::new())
                } else {
                    self.open_path(rest.to_owned());
                    Ok(format!("opened {rest}"))
                }
            }
            "style" => {
                let style = parse_style(rest)?;
                self.set_style(style);
                Ok(format!("style: {}", style.label()))
            }
            "lighting" => {
                let lighting = LightingPreset::parse(rest).ok_or_else(|| {
                    "expected default, soft, full, flat or illustrative".to_string()
                })?;
                self.set_lighting(lighting);
                Ok(format!("lighting: {rest}"))
            }
            "tonemap" => {
                self.view.tonemap = parse_on_off(rest, self.view.tonemap)?;
                Ok(format!(
                    "tonemap {}",
                    if self.view.tonemap { "on" } else { "off" }
                ))
            }
            "background" => {
                self.view.background = match rest {
                    "" => return Err(usage()),
                    "default" | "style" => self.view.style.background(),
                    other => parse_color(other)?,
                };
                Ok("background set".into())
            }
            "dof" => {
                self.view.dof = match rest {
                    "on" => 0.6,
                    "off" => 0.0,
                    other => other.parse::<f32>().map_err(|_| usage())?.clamp(0.0, 1.0),
                };
                Ok(format!("dof {}", self.view.dof))
            }
            "rotate" if rest.is_empty() => {
                self.open_popover("view", "view.rotate");
                Ok("rotate".into())
            }
            "rotate" => {
                let (axis_word, degrees) = rest.split_once(' ').ok_or_else(usage)?;
                let axis = match axis_word {
                    "x" => glam::Vec3::X,
                    "y" => glam::Vec3::Y,
                    "z" => glam::Vec3::Z,
                    _ => return Err(usage()),
                };
                let degrees: f32 = degrees.trim().parse().map_err(|_| usage())?;
                self.camera.rotate_scene(axis, degrees.to_radians());
                Ok(format!("rotated {degrees} about {axis_word}"))
            }
            "clip" if rest.is_empty() => {
                self.open_popover("view", "view.clip");
                Ok("clip".into())
            }
            "clip" => {
                let default = || crate::ui::ClipPlane {
                    orientation: crate::ui::ClipOrientation::View,
                    depth: 0.0,
                    flipped: false,
                };
                self.view.clip = match rest {
                    "off" => None,
                    "on" => Some(default()),
                    "view" | "x" | "y" | "z" => {
                        let orientation = match rest {
                            "view" => crate::ui::ClipOrientation::View,
                            "x" => crate::ui::ClipOrientation::X,
                            "y" => crate::ui::ClipOrientation::Y,
                            _ => crate::ui::ClipOrientation::Z,
                        };
                        let mut plane = self.view.clip.unwrap_or_else(default);
                        plane.orientation = orientation;
                        Some(plane)
                    }
                    "flip" => {
                        let mut plane = self.view.clip.unwrap_or_else(default);
                        plane.flipped = !plane.flipped;
                        Some(plane)
                    }
                    other => {
                        let depth: f32 = other.parse().map_err(|_| usage())?;
                        let mut plane = self.view.clip.unwrap_or_else(default);
                        plane.depth = depth;
                        Some(plane)
                    }
                };
                Ok(match self.view.clip {
                    Some(plane) => clip_description(plane),
                    None => "clip off".into(),
                })
            }
            "shadows" => {
                self.view.shadows = match rest {
                    "on" => 1.0,
                    "off" => 0.0,
                    other => other.parse::<f32>().map_err(|_| usage())?.clamp(0.0, 1.0),
                };
                Ok(format!("shadows {}", self.view.shadows))
            }
            "gradient" => {
                self.view.background_top = match rest {
                    "" => return Err(usage()),
                    "off" => None,
                    "on" => Some(crate::ui::gradient_top(self.view.background)),
                    other => Some(parse_color(other)?),
                };
                Ok("gradient set".into())
            }
            "fxaa" => {
                self.view.fxaa = parse_on_off(rest, self.view.fxaa)?;
                Ok(format!(
                    "fxaa {}",
                    if self.view.fxaa { "on" } else { "off" }
                ))
            }
            "framing" => {
                *self.studio_frame = parse_on_off(rest, *self.studio_frame)?;
                Ok(format!(
                    "framing {}",
                    if *self.studio_frame { "on" } else { "off" }
                ))
            }
            "ao" => {
                self.view.ao = parse_strength(rest, self.view.ao, 1.0, 2.0)?;
                Ok(format!("ambient occlusion {:.2}", self.view.ao))
            }
            "depthcue" => {
                self.view.depth_cue = parse_strength(rest, self.view.depth_cue, 0.5, 1.0)?;
                Ok(format!("depth cue {:.2}", self.view.depth_cue))
            }
            "outline" => {
                self.view.outline = parse_on_off(rest, self.view.outline)?;
                Ok(format!(
                    "outlines {}",
                    if self.view.outline { "on" } else { "off" }
                ))
            }
            "supersample" => {
                let scale = match rest {
                    "1" | "off" | "1x" => 1.0,
                    "1.5" | "1.5x" => 1.5,
                    "2" | "2x" | "on" => 2.0,
                    other => {
                        return Err(format!("expected 1, 1.5 or 2, got `{other}`"));
                    }
                };
                self.view.render_scale = scale;
                Ok(format!("supersampling {scale}x"))
            }
            "occlusion" => {
                self.view.occlusion_culling = parse_on_off(rest, self.view.occlusion_culling)?;
                Ok(format!(
                    "occlusion culling {}",
                    if self.view.occlusion_culling {
                        "on"
                    } else {
                        "off"
                    }
                ))
            }
            "adaptive" => {
                self.view.adaptive = parse_on_off(rest, self.view.adaptive)?;
                Ok(format!(
                    "adaptive quality {}",
                    if self.view.adaptive { "on" } else { "off" }
                ))
            }
            "view" => {
                let (what, args) = split_verb(rest);
                let args: Vec<&str> = args.split_whitespace().collect();
                match (what, args.as_slice()) {
                    ("reset", []) => {
                        self.mark_camera_jump();
                        self.reset_camera();
                        Ok("view reset".into())
                    }
                    ("orbit", [yaw, pitch]) => {
                        let yaw = parse_f32(yaw, "yaw")?.to_radians();
                        let pitch = parse_f32(pitch, "pitch")?.to_radians();
                        self.camera.orbit(yaw, pitch);
                        Ok("orbited".into())
                    }
                    ("zoom", [factor]) => {
                        let factor = parse_f32(factor, "zoom")?;
                        if factor <= 0.0 {
                            return Err("zoom factor must be positive".into());
                        }
                        self.camera.zoom(factor);
                        Ok("zoomed".into())
                    }
                    ("dolly", [by]) => {
                        self.camera.dolly_by(parse_f32(by, "dolly")?);
                        Ok(format!("dolly {:.1} A", self.camera.dolly))
                    }
                    ("face", [side]) => {
                        let (from, up) = crate::overlays::face(side)
                            .ok_or_else(|| format!("no face `{side}`: +X -X +Y -Y +Z -Z"))?;
                        self.mark_camera_jump();
                        self.camera.look_from(from, up);
                        Ok(format!("looking from {}", side.to_uppercase()))
                    }
                    ("projection", [mode]) => {
                        self.camera.projection = parse_projection(mode)?;
                        Ok(format!(
                            "projection {}",
                            projection_name(self.camera.projection)
                        ))
                    }
                    _ => Err(usage()),
                }
            }
            "layout" if rest.is_empty() => {
                self.open_popover("view", "view.layout");
                Ok("layout".into())
            }
            "layout" => {
                let known = crate::layout::built_in_workspaces()
                    .iter()
                    .any(|(n, _)| n.eq_ignore_ascii_case(rest))
                    || self
                        .workspaces
                        .iter()
                        .any(|(n, _)| n.eq_ignore_ascii_case(rest));
                if !known {
                    return Err(format!("no workspace named `{rest}`"));
                }
                *self.layout_request = Some(LayoutRequest::Switch(rest.to_owned()));
                Ok(format!("layout: {rest}"))
            }
            "startlayout" => match rest {
                "" => {
                    *self.layout_request = Some(LayoutRequest::SetStart);
                    Ok("every launch opens this layout".into())
                }
                "reset" => {
                    *self.layout_request = Some(LayoutRequest::ResetStart);
                    Ok("launches open the default layout again".into())
                }
                _ => Err(usage()),
            },
            "savelayout" => {
                if rest.is_empty() {
                    self.workspace_dialog.open = true;
                    Ok(String::new())
                } else {
                    *self.layout_request = Some(LayoutRequest::Save(rest.to_owned()));
                    Ok(format!("saved layout `{rest}`"))
                }
            }
            "screenshot" => {
                let mut words: Vec<&str> = rest.split_whitespace().collect();
                let ssaa = !words.contains(&"nossaa");
                words.retain(|w| *w != "nossaa");
                if words.is_empty() {
                    self.export_dialog.open = true;
                    return Ok(String::new());
                }
                let path = words.join(" ");
                *self.export_request = Some(ExportRequest {
                    path: path.clone().into(),
                    ssaa,
                });
                Ok(format!("exporting {path}"))
            }
            "render" if rest.is_empty() => {
                self.open_popover("file", "file.render");
                Ok("render".into())
            }
            "render" => {
                let mut words = script::split_args(rest);
                let transparent = words.iter().any(|w| w == "transparent");
                words.retain(|w| w != "transparent");
                let quality = words
                    .iter()
                    .find_map(|w| crate::render::Quality::parse(w))
                    .unwrap_or_default();
                words.retain(|w| crate::render::Quality::parse(w).is_none());
                let (vw, vh) = self.renderer.size();
                let mut request = crate::render::RenderRequest {
                    path: std::path::PathBuf::new(),
                    width: vw,
                    height: vh,
                    samples: quality.default_samples(),
                    transparent,
                    quality,
                };
                let mut rest_words = words.into_iter();
                let path = rest_words.next().ok_or_else(usage)?;
                request.path = (&path).into();
                for word in rest_words {
                    if let Some((w, h)) = word.split_once(['x', 'X']) {
                        request.width = w.parse().map_err(|_| usage())?;
                        request.height = h.parse().map_err(|_| usage())?;
                    } else {
                        request.samples = word.parse().map_err(|_| usage())?;
                    }
                }
                if request.width == 0 || request.height == 0 || request.samples == 0 {
                    return Err(usage());
                }
                let line = format!(
                    "rendering {path} ({}x{}, {} samples, {} quality)",
                    request.width,
                    request.height,
                    request.samples,
                    quality.name()
                );
                *self.render_request = Some(request);
                Ok(line)
            }
            "selectmode" if rest.starts_with("flyout") => {
                let key = match rest.strip_prefix("flyout").map(str::trim) {
                    Some("shapes") => crate::home::SHAPES_FLYOUT,
                    Some("more") => crate::home::MORE_FLYOUT,
                    Some("interface") => crate::home_interface::FORM_FLYOUT,
                    _ => return Err(usage()),
                };
                self.open_popover("home", key);
                Ok(format!("opened the {rest} flyout"))
            }
            "selectmode" => {
                self.view.select_tool.set(rest).map_err(|()| usage())?;
                if rest.starts_with("level") {
                    self.snap_selection_to_level();
                }
                if self.view.select_tool.draws() {
                    self.view.mouse_mode = crate::ui::MouseMode::Rotate;
                }
                Ok(format!("selection tool {rest}"))
            }
            "movie" => self.run_movie(rest),
            "mode" => {
                self.view.mouse_mode = match rest {
                    "rotate" => crate::ui::MouseMode::Rotate,
                    "translate" => crate::ui::MouseMode::Translate,
                    "scale" => crate::ui::MouseMode::Scale,
                    "center" => crate::ui::MouseMode::Center,
                    "label" => crate::ui::MouseMode::Label,
                    _ => return Err(usage()),
                };
                self.view.select_tool.shape = crate::select_tool::SelectShape::Click;
                Ok(format!("mouse mode {rest}"))
            }
            "ribbon" => {
                match rest {
                    "collapse" => self.prefs.ribbon_collapsed = true,
                    "expand" => self.prefs.ribbon_collapsed = false,
                    name => {
                        let tab = crate::ribbon::tab_named(name).ok_or_else(usage)?;
                        self.ribbon.tab = tab;
                        self.prefs.ribbon_collapsed = false;
                    }
                }
                Ok(format!("ribbon {rest}"))
            }
            // "settings"/"look" named the Look panel before it was folded
            // into the Look ribbon tab; the deep link still works, just
            // as a tab switch instead of a panel, with a toast saying so.
            "panel" if rest == "settings" || rest == "look" => {
                self.ribbon.tab = crate::ribbon::tab_named("look").unwrap_or(self.ribbon.tab);
                self.prefs.ribbon_collapsed = false;
                *self.notice = Some(crate::ui::Notice::info(
                    "Look settings are now in the Look tab.".into(),
                ));
                Ok("showing the Look tab".into())
            }
            "panel" if rest.is_empty() => {
                self.open_popover("view", "view.panels");
                Ok("panel".into())
            }
            "panel" if rest.ends_with(" collapse") || rest.ends_with(" expand") => {
                let (name, how) = rest.rsplit_once(' ').ok_or_else(usage)?;
                let tab = crate::ribbon::panel_named(name).ok_or_else(usage)?;
                *self.layout_request = Some(LayoutRequest::CollapsePanel(tab, how == "collapse"));
                Ok(format!("{how} {}", tab.title()))
            }
            "panel" => {
                let tab = crate::ribbon::panel_named(rest).ok_or_else(usage)?;
                self.ribbon.focus_expression = tab == Tab::Selection;
                *self.layout_request = Some(LayoutRequest::OpenPanel(tab));
                Ok(format!("showing {}", tab.title()))
            }
            "play" => {
                if crate::ui::timeline_target(self.scene).is_none() {
                    return Err("no structure with several frames is loaded".into());
                }
                if let Ok(fps) = rest.parse::<f32>() {
                    if !(fps > 0.0 && fps <= 240.0) {
                        return Err("frames per second must be in (0, 240]".into());
                    }
                    self.timeline.fps = fps;
                    self.timeline.playing = true;
                } else {
                    self.timeline.playing = parse_on_off(rest, self.timeline.playing)?;
                }
                Ok(if self.timeline.playing {
                    "playing"
                } else {
                    "paused"
                }
                .into())
            }
            "deletelayout" if rest.is_empty() => {
                self.workspace_dialog.deleting = true;
                Ok(String::new())
            }
            "deletelayout" => {
                let before = self.workspaces.len();
                self.workspaces
                    .retain(|(n, _)| !n.eq_ignore_ascii_case(rest));
                if self.workspaces.len() == before {
                    return Err(format!("no saved workspace named `{rest}`"));
                }
                Ok(format!("deleted workspace {rest}"))
            }
            "sequence" => self.sequence_command(rest),
            "preferences" => {
                self.prefs_dialog.open = true;
                Ok(String::new())
            }
            "look" => {
                let key = match rest {
                    "lights" => "look.lights",
                    "background" => "look.background",
                    "ao" | "shadows" | "depthcue" | "dof" => {
                        self.prefs_dialog.open = true;
                        return Ok(format!("look {rest}"));
                    }
                    _ => return Err(usage()),
                };
                self.open_popover("look", key);
                Ok(format!("look {rest}"))
            }
            "theme" => {
                let mode = crate::theme::ThemeMode::parse(rest).ok_or_else(usage)?;
                self.prefs.theme = mode;
                Ok(format!("theme {}", mode.name()))
            }
            "welcome" => match rest {
                "" => {
                    *self.start_card_dismissed = false;
                    Ok("showing the start card".into())
                }
                _ => {
                    let show = parse_on_off(rest, !self.prefs.hide_start_card)?;
                    self.prefs.hide_start_card = !show;
                    *self.start_card_dismissed = !show;
                    Ok(format!(
                        "start card at launch: {}",
                        if show { "on" } else { "off" }
                    ))
                }
            },
            "confirmquit" => {
                let confirm = parse_on_off(rest, !self.prefs.dont_confirm_quit)?;
                self.prefs.dont_confirm_quit = !confirm;
                Ok(format!(
                    "confirm quit {}",
                    if confirm { "on" } else { "off" }
                ))
            }
            "uishot" => {
                if rest.is_empty() {
                    return Err(usage());
                }
                // A theme, layout or popover change made just before
                // lands (and finishes fading in) a few frames later.
                *self.ui_shot = Some((rest.into(), 12));
                Ok(format!("capturing the window to {rest}"))
            }
            "window" => {
                let size = rest
                    .split_once('x')
                    .and_then(|(w, h)| Some((w.trim().parse().ok()?, h.trim().parse().ok()?)));
                let Some((w, h)) = size else {
                    return Err(usage());
                };
                *self.window_request = Some((w, h));
                Ok(format!("window {w}x{h}"))
            }
            "quit" => {
                *self.quit_requested = true;
                Ok("quitting".into())
            }
            other => unreachable!("every view entry is handled: {other}"),
        }
    }

    /// The view half of a session file: what `vv_scene::session` carries
    /// as opaque JSON. Each structure's current frame is document state
    /// now (`vv_scene::session::SavedStructure::frame`, restored through
    /// the document half via `Command::SetFrame`), not part of this block
    /// — see that module's doc for why an old file's `"frames"` key here
    /// (from before that existed) is deliberately no longer read.
    fn view_snapshot(&self) -> serde_json::Value {
        // Play clips are saved by position in the scene, as the session's
        // own structures are, since ids differ once loaded again.
        let order: Vec<u32> = self.scene.structures().map(|(id, _)| id.to_raw()).collect();
        let mut movie = self.movie.movie.clone();
        movie.remap_structures(|raw| order.iter().position(|&id| id == raw).map(|i| i as u32));
        let mut view = capture_view(
            self.view,
            self.camera,
            self.dock_layout,
            self.timeline,
            &movie,
        );
        view["sequence"] = self.sequence.snapshot();
        view
    }

    /// Applies a `view` block written by `view_snapshot`. `ids` maps a
    /// saved structure's position to the id it got on load, which the
    /// movie's play clips need.
    fn apply_view_snapshot(&mut self, view: &serde_json::Value, ids: &[Option<StructureId>]) {
        self.movie.release_camera();
        self.movie.selected = None;
        self.sequence.restore(&view["sequence"]);
        if let Some(dock) = restore_view(
            view,
            self.view,
            self.camera,
            self.timeline,
            &mut self.movie.movie,
        ) {
            *self.layout_request = Some(LayoutRequest::LoadState(Box::new(dock)));
        }
        self.movie.movie.remap_structures(|saved| {
            ids.get(saved as usize)
                .copied()
                .flatten()
                .map(|id| id.to_raw())
        });
        // After the camera, so an old raw-plane session's depth is measured
        // from the target it actually had.
        if view["clip"].is_object() || view["clip"].is_array() {
            self.view.clip = parse_clip_snapshot(&view["clip"], self.camera.target);
        }
        // Restoring a session isn't a user edit to undo -- rebases
        // the view-edit baseline instead of letting the frame's own
        // `commit_view_edit` record it as one.
        *self.view_baseline = self.view.clone();
    }

    fn save_session(&mut self, path: &str) -> Result<String, String> {
        let path = if path.is_empty() {
            match rfd::FileDialog::new()
                .set_title("Save Session")
                .set_file_name("session.vviz")
                .add_filter("vizviz session", &["vviz"])
                .save_file()
            {
                Some(p) => p,
                None => return Ok(String::new()),
            }
        } else {
            path.into()
        };
        let view = self.view_snapshot();
        vv_scene::session::save(self.scene, &path, Some(view)).map_err(|e| e.to_string())?;
        *self.saved_edits = self.history.edits();
        Ok(format!("saved session to {}", path.display()))
    }

    fn load_session(&mut self, path: &str) -> Result<String, String> {
        let path = if path.is_empty() {
            match rfd::FileDialog::new()
                .set_title("Open Session")
                .add_filter("vizviz session", &["vviz", "json"])
                .pick_file()
            {
                Some(p) => p,
                None => return Ok(String::new()),
            }
        } else {
            path.into()
        };
        let file = vv_scene::session::read(&path).map_err(|e| e.to_string())?;
        self.prefs.add_recent(path.clone());
        let applied = vv_scene::session::apply(&file, self.scene, self.history);
        *self.saved_edits = self.history.edits();
        *self.info_target = None;
        match &file.view {
            Some(view) => self.apply_view_snapshot(view, &applied.ids),
            None => self.reset_camera(),
        }
        let mut lines = vec![format!(
            "loaded session {}: {} structure(s), {} set(s)",
            path.display(),
            applied.ids.iter().flatten().count(),
            self.scene.selection_sets().len()
        )];
        lines.extend(applied.warnings.iter().map(|w| format!("warning: {w}")));
        Ok(lines.join("\n"))
    }

    /// `help` over both registries.
    fn help(&self, word: &str) -> String {
        if word.is_empty() {
            return entries()
                .iter()
                .map(|e| format!("{:<44} {}", e.usage, e.title))
                .collect::<Vec<_>>()
                .join("\n");
        }
        match entries().into_iter().find(|e| e.id == word) {
            Some(e) => format!("{}\n  {}", e.usage, e.help),
            None => match view_entry(word) {
                Some(e) => format!("{}\n  {}", e.usage, e.help),
                None => script::help(word),
            },
        }
    }

    /// Runs a command and puts the echo and its result (or error) in the
    /// log, one line each, the way a console does. Returns whether it
    /// succeeded, so `--exec` can stop at the first failure.
    pub fn run_command_logged(&mut self, line: &str) -> bool {
        let line = line.trim();
        if line.is_empty() {
            return true;
        }
        self.log.push(format!("> {line}"));
        // Every state-changing command toasts (UX_PATTERNS
        // "Feedback"), using its own first line -- already sentence-
        // shaped ("selected 214 atom(s)"). Compared against `notice`
        // afterward so a handler with something more specific to say
        // (`open_path`'s "already open") isn't clobbered; an Undo action
        // is attached exactly when this command actually entered undo
        // history (`edits()` grew) -- e.g. not true of interaction-only view state.
        let notice_before = self.notice.as_ref().map(|n| n.when);
        let edits_before = self.history.edits();
        match self.run_command(line) {
            Ok(out) => {
                self.log.extend(out.lines().map(str::to_owned));
                let verb = split_verb(line).0;
                // Analyze > Contacts.../Surface area show their
                // result in the Inspector, not just the (often-collapsed)
                // Log.
                if let Some(title) = result_card_title(verb) {
                    *self.result_card = Some(crate::ui::ResultCard {
                        title: title.into(),
                        body: out.clone(),
                    });
                }
                let notice_already_set = self.notice.as_ref().map(|n| n.when) != notice_before;
                if !out.is_empty() && !notice_already_set && !TOAST_SKIP.contains(&verb) {
                    let message = out.lines().next().unwrap_or_default().to_owned();
                    let mut notice = crate::ui::Notice::info(message);
                    if self.history.edits() > edits_before {
                        notice = notice.with_action("Undo", "undo");
                    }
                    *self.notice = Some(notice);
                }
                true
            }
            Err(e) => {
                self.fail(e);
                false
            }
        }
    }
}

/// Which commands' output becomes an Inspector result card, and its title.
fn result_card_title(verb: &str) -> Option<&'static str> {
    match verb {
        "sasa" => Some("Surface area"),
        "contacts" => Some("Contacts"),
        _ => None,
    }
}

/// Verbs the generic toast (`run_command_logged`) skips: pure reads
/// (their own panel, not a toast, is the point; toasts surface results, not
/// every list) and interaction-state changes
/// frequent enough that toasting them would be noise, not signal (`mode`
/// especially -- R T S C are single-key, pressed rapidly).
const TOAST_SKIP: &[&str] = &[
    "help",
    "version",
    "uishot",
    "window",
    "structures",
    "reps",
    "sets",
    "labels",
    "captions",
    "measurements",
    "info",
    "mode",
    "selectmode",
    "sequence",
    "ribbon",
    "panel",
    "undo",
    "redo",
    "repopt",
];

#[cfg(test)]
mod tests {
    use super::*;

    /// DESIGN.md's wording rules, on everything the palette shows:
    /// sentence case, "…" not "...", and no run of spaces (the mark of a
    /// string continued without its backslash).
    #[test]
    fn command_text_follows_the_wording_rules() {
        const ACRONYMS: &[&str] = &["FXAA", "PDB", "SES", "SAS", "SASA"];
        for e in entries() {
            for text in [e.title, e.help, e.usage] {
                assert!(!text.contains("  "), "{}: run of spaces in {text:?}", e.id);
            }
            assert!(!e.title.contains("..."), "{}: {:?}", e.id, e.title);
            for word in e.title.split_whitespace().skip(1) {
                let capital = word.chars().next().is_some_and(char::is_uppercase);
                assert!(
                    !capital || ACRONYMS.contains(&word),
                    "{}: {:?} is not sentence case",
                    e.id,
                    e.title
                );
            }
        }
    }

    #[test]
    fn validate_accepts_both_registries_and_rejects_typos() {
        assert!(validate("load x.cif").is_ok());
        assert!(validate("style white").is_ok());
        assert!(validate("exit").is_ok());
        assert!(validate("").is_ok());
        assert!(validate("laod x.cif")
            .unwrap_err()
            .contains("unknown command `laod`"));
    }

    #[test]
    fn scripts_refuse_the_forms_that_open_a_dialog() {
        for verb in [
            "open",
            "fetch",
            "savesession",
            "loadsession",
            "savestructure",
            "screenshot",
            "savelayout",
        ] {
            let refused = validate_unattended(verb).unwrap_err();
            assert!(refused.contains("opens a dialog"), "{verb}: {refused}");
            assert!(
                validate_unattended(&format!("{verb}   ")).is_err(),
                "{verb}"
            );
            assert!(
                validate_unattended(&format!("{verb} x")).is_ok(),
                "{verb} x"
            );
        }
        assert!(validate_unattended("render out.png").is_ok());
        assert!(validate_unattended("laod x.cif").is_err());
    }

    #[test]
    fn entries_have_unique_ids_and_sensible_arg_flags() {
        let all = entries();
        let mut ids: Vec<&str> = all.iter().map(|e| e.id).collect();
        ids.sort_unstable();
        ids.dedup();
        assert_eq!(ids.len(), all.len(), "duplicate command id");
        let by_id = |id: &str| all.iter().find(|e| e.id == id).unwrap();
        assert!(by_id("load").needs_args());
        assert!(by_id("style").needs_args());
        assert!(!by_id("open").needs_args());
        assert!(!by_id("quit").needs_args());
        assert!(!by_id("close").needs_args());
    }

    #[test]
    fn colors_parse_from_components_and_hex() {
        let white = parse_color("255 255 255").unwrap();
        assert!((white.r - 1.0).abs() < 1e-9 && (white.a - 1.0).abs() < 1e-9);
        let black = parse_color("#000000").unwrap();
        assert_eq!((black.r, black.g, black.b), (0.0, 0.0, 0.0));
        // sRGB 128 decodes to ~0.216 linear, not 0.5.
        let mid = parse_color("#808080").unwrap();
        assert!((mid.r - 0.2158).abs() < 1e-3, "{}", mid.r);
        let translucent = parse_color("#ff000080").unwrap();
        assert!((translucent.a - 128.0 / 255.0).abs() < 1e-9);
        assert!(parse_color("1 2").is_err());
        assert!(parse_color("#12345").is_err());
        assert!(parse_color("red").is_err());
    }

    #[test]
    fn on_off_words_and_toggle() {
        assert_eq!(parse_on_off("on", false), Ok(true));
        assert_eq!(parse_on_off("off", true), Ok(false));
        assert_eq!(parse_on_off("", true), Ok(false));
        assert!(parse_on_off("maybe", true).is_err());
    }

    /// A rich, non-default view (style, lights, effects, clip, camera
    /// dolly/projection, a non-default dock layout, and playback) survives
    /// `capture_view`/`restore_view` -- the pieces of a session's opaque
    /// `view` block that a GPU-backed `State` would otherwise be needed to
    /// exercise (item 5's audit: dock layout and playback were the two
    /// actually missing from what a session restored).
    #[test]
    fn capture_and_restore_round_trip_a_rich_view() {
        let mut lights = crate::studio::LightRig::of(LightingPreset::Full);
        lights.lights[0] = crate::studio::StudioLight {
            enabled: true,
            color: egui::Color32::from_rgb(0xFF, 0x40, 0x40),
            intensity: 1.7,
            azimuth: 20.0,
            elevation: 35.0,
        };
        lights.lights[1].intensity = 0.4;
        lights.ambient = 0.2;
        lights.ambient_color = egui::Color32::from_rgb(0x20, 0x20, 0xE0);
        lights.add(); // a 3rd light, then switched off below
        lights.lights[2].enabled = false;
        let view = ViewSettings {
            style: StylePreset::Glossy,
            lighting: LightingPreset::Full,
            lights,
            lights_custom: true,
            tonemap: true,
            saved_views: vec![crate::studio::SavedView {
                name: "Front".into(),
                target: glam::Vec3::new(4.0, 5.0, 6.0),
                distance: 12.0,
                orientation: glam::Quat::from_rotation_x(0.3),
                fov_y: 0.8,
                near: 0.2,
                projection: Projection::Perspective,
                scene_radius: 9.0,
                dolly: 0.5,
            }],
            background: wgpu::Color {
                r: 0.1,
                g: 0.2,
                b: 0.3,
                a: 1.0,
            },
            background_top: Some(wgpu::Color {
                r: 0.4,
                g: 0.5,
                b: 0.6,
                a: 1.0,
            }),
            fxaa: false,
            outline: true,
            ao: 1.5,
            depth_cue: 0.6,
            shadows: 0.8,
            dof: 0.3,
            clip: Some(crate::ui::ClipPlane {
                orientation: crate::ui::ClipOrientation::Y,
                depth: -5.0,
                flipped: true,
            }),
            render_scale: 1.5,
            occlusion_culling: false,
            adaptive: false,
            ..ViewSettings::default()
        };

        let camera = Camera {
            target: glam::Vec3::new(1.0, 2.0, 3.0),
            distance: 42.0,
            orientation: glam::Quat::from_rotation_y(0.4) * glam::Quat::from_rotation_x(0.2),
            fov_y: 0.9,
            near: 0.05,
            projection: Projection::Orthographic,
            scene_radius: 17.5,
            pivot: None,
            dolly: 3.25,
        };

        let mut dock_layout = crate::layout::compact_layout();
        crate::layout::give_every_node_a_real_rect(&mut dock_layout);

        let timeline = TimelineState {
            playing: true,
            fps: 42.0,
            looping: false,
            ..TimelineState::default()
        };

        let mut movie = crate::movie::Movie::default();
        movie.set_fps(24.0).unwrap();
        let json = capture_view(&view, &camera, &dock_layout, &timeline, &movie);

        let mut view2 = ViewSettings::default();
        let mut camera2 = Camera::framing(glam::Vec3::ZERO, 1.0);
        let mut timeline2 = TimelineState::default();
        let mut movie2 = crate::movie::Movie::default();
        let restored_dock =
            restore_view(&json, &mut view2, &mut camera2, &mut timeline2, &mut movie2)
                .expect("dock layout");
        assert_eq!(movie2, movie);

        assert_eq!(view2.style, StylePreset::Glossy);
        assert_eq!(view2.lighting, LightingPreset::Full);
        assert_eq!(view2.lights.count, 3);
        assert_eq!(view2.lights.lights[0].intensity, 1.7);
        assert_eq!(view2.lights.lights[0].azimuth, 20.0);
        assert_eq!(view2.lights.lights[0].elevation, 35.0);
        assert_eq!(
            view2.lights.lights[0].color,
            egui::Color32::from_rgb(0xFF, 0x40, 0x40)
        );
        assert_eq!(view2.lights.lights[1].intensity, 0.4);
        assert!(!view2.lights.lights[2].enabled);
        assert_eq!(view2.lights.ambient, 0.2);
        assert_eq!(
            view2.lights.ambient_color,
            egui::Color32::from_rgb(0x20, 0x20, 0xE0)
        );
        assert_eq!(view2.saved_views.len(), 1);
        assert_eq!(view2.saved_views[0].name, "Front");
        assert_eq!(view2.saved_views[0].target, glam::Vec3::new(4.0, 5.0, 6.0));
        assert_eq!(view2.saved_views[0].distance, 12.0);
        assert_eq!(view2.saved_views[0].projection, Projection::Perspective);
        assert!(view2.lights_custom, "lights_custom must round-trip");
        assert!(view2.tonemap);
        assert!(!view2.fxaa);
        assert!(view2.outline);
        assert_eq!(view2.ao, 1.5);
        assert_eq!(view2.depth_cue, 0.6);
        assert_eq!(view2.shadows, 0.8);
        assert_eq!(view2.dof, 0.3);
        assert_eq!(view2.clip, view.clip);
        assert_eq!(view2.render_scale, 1.5);
        assert!(!view2.occlusion_culling);
        assert!(!view2.adaptive);
        assert_eq!(
            (view2.background.r, view2.background.g, view2.background.b),
            (0.1, 0.2, 0.3)
        );
        assert_eq!(
            view2.background_top.map(|c| (c.r, c.g, c.b)),
            Some((0.4, 0.5, 0.6))
        );

        assert_eq!(camera2.target, camera.target);
        assert_eq!(camera2.distance, camera.distance);
        assert_eq!(camera2.projection, Projection::Orthographic);
        assert_eq!(camera2.scene_radius, camera.scene_radius);
        assert_eq!(camera2.dolly, camera.dolly);
        assert!(
            camera2.orientation.dot(camera.orientation).abs() > 0.9999,
            "{:?} vs {:?}",
            camera2.orientation,
            camera.orientation
        );

        assert!(timeline2.playing);
        assert_eq!(timeline2.fps, 42.0);
        assert!(!timeline2.looping);

        // The dock layout that comes back has the same tabs as what was
        // saved (the actual regression this test guards: a session used
        // to drop the layout entirely).
        let mut original_tabs: Vec<_> = dock_layout.iter_all_tabs().map(|(_, t)| *t).collect();
        let mut restored_tabs: Vec<_> = restored_dock.iter_all_tabs().map(|(_, t)| *t).collect();
        original_tabs.sort_by_key(|t| *t as usize);
        restored_tabs.sort_by_key(|t| *t as usize);
        assert_eq!(original_tabs, restored_tabs);
        assert!(restored_dock.find_tab(&Tab::Timeline).is_none());
        assert!(restored_dock.find_tab(&Tab::Inspector).is_some());
    }

    /// A schema-3 session (no `dock_layout`/`playback` keys) leaves both
    /// as they already are, rather than failing to restore anything else.
    #[test]
    fn restoring_a_schema_3_view_leaves_dock_layout_and_playback_alone() {
        let old_view = serde_json::json!({
            "style": "white",
            "camera": {
                "target": [0.0, 0.0, 0.0],
                "distance": 10.0,
                "orientation": [0.0, 0.0, 0.0, 1.0],
                "fov_y": 0.7,
                "near": 0.1,
            },
        });
        let mut view = ViewSettings::default();
        let mut camera = Camera::framing(glam::Vec3::ZERO, 1.0);
        let mut timeline = TimelineState {
            playing: true,
            fps: 24.0,
            ..TimelineState::default()
        };
        let mut movie = crate::movie::Movie::default();
        movie.set_fps(60.0).unwrap();
        let dock = restore_view(&old_view, &mut view, &mut camera, &mut timeline, &mut movie);
        assert!(dock.is_none());
        assert_eq!(
            movie,
            crate::movie::Movie::default(),
            "a session without a movie clears it"
        );
        assert!(timeline.playing);
        assert_eq!(timeline.fps, 24.0);
    }

    /// A session saved before the light list existed stored one scalar
    /// key/fill pair (`capture_view`'s old shape); it must still load,
    /// landing on lights 0/1 -- both still enabled and in the list, since
    /// a preset (or the defaults) always fills those two slots first.
    #[test]
    fn old_session_scalar_lights_still_load() {
        let old_view = serde_json::json!({
            "lights": {
                "key": 1.2,
                "key_azimuth": 10.0,
                "key_elevation": 20.0,
                "fill": 0.3,
                "ambient": 0.5,
            },
        });
        let mut view = ViewSettings::default();
        let mut camera = Camera::framing(glam::Vec3::ZERO, 1.0);
        let mut timeline = TimelineState::default();
        let mut movie = crate::movie::Movie::default();
        restore_view(&old_view, &mut view, &mut camera, &mut timeline, &mut movie);
        assert!(view.lights.lights[0].enabled);
        assert_eq!(view.lights.lights[0].intensity, 1.2);
        assert_eq!(view.lights.lights[0].azimuth, 10.0);
        assert_eq!(view.lights.lights[0].elevation, 20.0);
        assert!(view.lights.lights[1].enabled);
        assert_eq!(view.lights.lights[1].intensity, 0.3);
        assert_eq!(view.lights.ambient, 0.5);
        assert_eq!(view.lights.count, 2);
    }
}
