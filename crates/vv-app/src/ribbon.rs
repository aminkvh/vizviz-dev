//! The ribbon: tabs by task, labelled groups, big and small buttons
//! (docs/UI_DESIGN.md). Every action runs a
//! command line through the same runner as `--exec`, the palette and the
//! console, so anything clicked here can also be typed or scripted; the
//! `tests` hold every command to being reachable from here (or listed as
//! console-only, with a reason).
//!
//! An action is a plain command, a prompt (the palette opens with a verb
//! that needs typed input), a toggle (a command pair plus the state it
//! reflects), or a choice (options, each a command, plus which one is
//! current). Each has a key tip, a letter unique within its tab, for
//! reaching it from the keyboard.

use egui::Ui;
use egui_phosphor::regular as icon;

use crate::layout::Tab;
use crate::ui::AppUi;
use crate::widgets::{self, Variant};

pub type Query = fn(&AppUi<'_>) -> bool;
pub type Current = fn(&AppUi<'_>) -> Option<String>;

pub enum Kind {
    Run(&'static str),
    /// Opens the palette with this text, for a verb that needs input.
    Prompt(&'static str),
    Toggle {
        on: &'static str,
        off: &'static str,
        state: Query,
    },
    /// A color picker: `get` reads the color shown, and picking one runs
    /// `verb` followed by its hex (`background #FFFFFF`).
    Color {
        get: fn(&AppUi<'_>) -> egui::Color32,
        verb: &'static str,
    },
    /// `(label, command)` options; `current` names the current command.
    Choice {
        options: &'static [(&'static str, &'static str)],
        current: Current,
        /// A drop-down instead of a gallery of buttons (long lists).
        dropdown: bool,
        /// Disabled, with a tooltip saying why, while this holds --
        /// Coloring on a glycan rep ("SNFG colours are fixed").
        disabled_when: Option<(Query, &'static str)>,
    },
    /// A `NAME ▾` popover (UX_PATTERNS "Edit: popover form"): `key`
    /// identifies it for a deep link to force it open (`RibbonState`'s
    /// `pending_popover`); `commands` are the deep-link commands it holds,
    /// for `help`/tests/tooltips; `body` draws its contents. `ellipsis`
    /// reads it as "NAME…" instead --
    /// DESIGN.md's form trigger, for a popover that asks for input
    /// before it acts (Rotate) rather than a settings disclosure
    /// (Lights, Clip, Layout, ...).
    Popover {
        key: &'static str,
        commands: &'static [&'static str],
        body: PopoverBody,
        ellipsis: bool,
    },
    /// A row of compact controls drawn by `body`, which returns the
    /// command line a click chose; `commands` lists every command it can
    /// run (the first is what its key tip runs).
    Inline {
        commands: fn() -> Vec<&'static str>,
        body: InlineBody,
    },
    /// A disabled action with a tooltip saying why (e.g. "coming soon").
    Disabled {
        reason: &'static str,
    },
    /// An on/off effect (a switch, `on`/`off` its commands) with its
    /// strength one click away in a ▾ popover, shown only while it's on.
    SwitchEffect {
        on: &'static str,
        off: &'static str,
        state: Query,
        key: &'static str,
        body: PopoverBody,
    },
}

pub type PopoverBody = fn(&mut AppUi<'_>, &mut Ui);
pub type InlineBody = fn(&mut AppUi<'_>, &mut Ui) -> Option<String>;

pub struct Action {
    pub label: &'static str,
    pub icon: &'static str,
    /// Key tip: the letter(s) that reach this action while its tab's key
    /// tips show.
    pub tip: &'static str,
    pub kind: Kind,
    /// A large icon-over-label button (the group's main action).
    pub big: bool,
}

pub struct Group {
    pub label: &'static str,
    pub actions: &'static [Action],
}

pub struct RibbonTab {
    pub label: &'static str,
    pub tip: &'static str,
    pub groups: &'static [Group],
    /// Contextual tabs show only while this holds.
    pub shown: Option<Query>,
}

const fn run(
    label: &'static str,
    icon: &'static str,
    tip: &'static str,
    cmd: &'static str,
) -> Action {
    Action {
        label,
        icon,
        tip,
        kind: Kind::Run(cmd),
        big: false,
    }
}

const fn big(
    label: &'static str,
    icon: &'static str,
    tip: &'static str,
    cmd: &'static str,
) -> Action {
    Action {
        label,
        icon,
        tip,
        kind: Kind::Run(cmd),
        big: true,
    }
}

const fn inline(
    label: &'static str,
    icon: &'static str,
    tip: &'static str,
    commands: fn() -> Vec<&'static str>,
    body: InlineBody,
) -> Action {
    Action {
        label,
        icon,
        tip,
        kind: Kind::Inline { commands, body },
        big: false,
    }
}

const fn prompt(
    label: &'static str,
    icon: &'static str,
    tip: &'static str,
    text: &'static str,
) -> Action {
    Action {
        label,
        icon,
        tip,
        kind: Kind::Prompt(text),
        big: false,
    }
}

const fn toggle(
    label: &'static str,
    icon: &'static str,
    tip: &'static str,
    on: &'static str,
    off: &'static str,
    state: Query,
) -> Action {
    Action {
        label,
        icon,
        tip,
        kind: Kind::Toggle { on, off, state },
        big: false,
    }
}

const fn color(
    label: &'static str,
    icon: &'static str,
    tip: &'static str,
    get: fn(&AppUi<'_>) -> egui::Color32,
    verb: &'static str,
) -> Action {
    Action {
        label,
        icon,
        tip,
        kind: Kind::Color { get, verb },
        big: false,
    }
}

const fn choice(
    label: &'static str,
    icon: &'static str,
    tip: &'static str,
    options: &'static [(&'static str, &'static str)],
    current: Current,
    dropdown: bool,
) -> Action {
    Action {
        label,
        icon,
        tip,
        kind: Kind::Choice {
            options,
            current,
            dropdown,
            disabled_when: None,
        },
        big: false,
    }
}

/// `choice`, disabled (with a tooltip saying why) while `disabled_when`
/// holds -- e.g. Coloring on a glycan rep.
#[allow(clippy::too_many_arguments)]
const fn choice_disabled_when(
    label: &'static str,
    icon: &'static str,
    tip: &'static str,
    options: &'static [(&'static str, &'static str)],
    current: Current,
    dropdown: bool,
    disabled_when: Query,
    reason: &'static str,
) -> Action {
    Action {
        label,
        icon,
        tip,
        kind: Kind::Choice {
            options,
            current,
            dropdown,
            disabled_when: Some((disabled_when, reason)),
        },
        big: false,
    }
}

const fn popover(
    label: &'static str,
    icon: &'static str,
    tip: &'static str,
    key: &'static str,
    commands: &'static [&'static str],
    body: PopoverBody,
) -> Action {
    Action {
        label,
        icon,
        tip,
        kind: Kind::Popover {
            key,
            commands,
            body,
            ellipsis: false,
        },
        big: false,
    }
}

/// `popover`, but reading "NAME…": for a
/// popover that's a form asking for input before it acts (Rotate),
/// DESIGN.md's trigger for that, as opposed to a settings disclosure.
const fn popover_form(
    label: &'static str,
    icon: &'static str,
    tip: &'static str,
    key: &'static str,
    commands: &'static [&'static str],
    body: PopoverBody,
) -> Action {
    Action {
        label,
        icon,
        tip,
        kind: Kind::Popover {
            key,
            commands,
            body,
            ellipsis: true,
        },
        big: false,
    }
}

// Currently unused; kept for the next disabled ribbon button.
#[allow(dead_code)]
const fn disabled(
    label: &'static str,
    icon: &'static str,
    tip: &'static str,
    reason: &'static str,
) -> Action {
    Action {
        label,
        icon,
        tip,
        kind: Kind::Disabled { reason },
        big: false,
    }
}

#[allow(clippy::too_many_arguments)]
const fn switch_effect(
    label: &'static str,
    icon: &'static str,
    tip: &'static str,
    on: &'static str,
    off: &'static str,
    state: Query,
    key: &'static str,
    body: PopoverBody,
) -> Action {
    Action {
        label,
        icon,
        tip,
        kind: Kind::SwitchEffect {
            on,
            off,
            state,
            key,
            body,
        },
        big: false,
    }
}

// ---- state queries ----------------------------------------------------------

fn current_rep(a: &AppUi<'_>) -> Option<vv_scene::Rep> {
    Some(a.scene.structure(a.current()?)?.rep().clone())
}
fn rep_style(a: &AppUi<'_>) -> Option<String> {
    let r = current_rep(a)?;
    Some(format!(
        "rep {}",
        vv_scene::session::representation_name(r.representation)
    ))
}
/// The current rep is glycan, so its own Coloring choice is
/// disabled -- SNFG's colours come from the recognized monosaccharide,
/// not the usual coloring schemes, which this rep silently ignores.
fn rep_is_glycan(a: &AppUi<'_>) -> bool {
    current_rep(a).is_some_and(|r| r.representation == vv_scene::Representation::Glycan)
}
fn rep_coloring(a: &AppUi<'_>) -> Option<String> {
    Some(format!("color {}", current_rep(a)?.coloring.name()))
}
fn rep_material(a: &AppUi<'_>) -> Option<String> {
    Some(format!("material {}", current_rep(a)?.material.name()))
}
fn look(a: &AppUi<'_>) -> Option<String> {
    Some(format!(
        "style {}",
        crate::commands::style_name(a.view.style)
    ))
}
fn projection(a: &AppUi<'_>) -> Option<String> {
    Some(match a.camera.projection {
        vv_render::Projection::Perspective => "view projection perspective".into(),
        vv_render::Projection::Orthographic => "view projection orthographic".into(),
    })
}
/// The current rep's single color, or mid gray when it is colored by
/// something else.
fn rep_color(a: &AppUi<'_>) -> egui::Color32 {
    match current_rep(a).map(|r| r.coloring) {
        Some(vv_scene::ColorScheme::Constant([r, g, b])) => egui::Color32::from_rgb(r, g, b),
        _ => egui::Color32::GRAY,
    }
}
fn clip_on(a: &AppUi<'_>) -> bool {
    a.view.clip.is_some()
}
fn outline_on(a: &AppUi<'_>) -> bool {
    a.view.outline
}
fn label_mode_on(a: &AppUi<'_>) -> bool {
    a.view.mouse_mode == crate::ui::MouseMode::Label
}
fn interaction_on(a: &AppUi<'_>, kind: vv_core::interactions::InteractionKind) -> bool {
    a.current()
        .and_then(|id| a.scene.structure(id))
        .is_some_and(|s| s.interactions.contains(&kind))
}
fn hbond_on(a: &AppUi<'_>) -> bool {
    interaction_on(a, vv_core::interactions::InteractionKind::Hbond)
}
fn metal_on(a: &AppUi<'_>) -> bool {
    interaction_on(a, vv_core::interactions::InteractionKind::Metal)
}
fn saltbridge_on(a: &AppUi<'_>) -> bool {
    interaction_on(a, vv_core::interactions::InteractionKind::SaltBridge)
}
fn framing_on(a: &AppUi<'_>) -> bool {
    *a.studio_frame
}
/// File ▸ Recent's popover body: the recently opened files, most recent
/// first (`open`/`loadsession` per extension, as the start card does).
fn recent_popover(app: &mut AppUi<'_>, ui: &mut Ui) {
    ui.set_max_width(280.0);
    let recents = app.prefs.recent_files.clone();
    if recents.is_empty() {
        widgets::caption(ui, "Nothing opened yet");
        return;
    }
    let mut run = None;
    for path in recents {
        let name = path
            .file_name()
            .map(|n| n.to_string_lossy().into_owned())
            .unwrap_or_else(|| path.display().to_string());
        let verb = if path
            .extension()
            .is_some_and(|e| e.eq_ignore_ascii_case("vviz") || e.eq_ignore_ascii_case("json"))
        {
            "loadsession"
        } else {
            "open"
        };
        if widgets::list_row(ui, &name, false, |_| {})
            .on_hover_text(path.display().to_string())
            .clicked()
        {
            run = Some(format!("{verb} {}", path.display()));
            ui.close();
        }
    }
    if let Some(line) = run {
        app.run_command_logged(&line);
    }
}

// ---- the tabs -----------------------------------------------------------------

pub(crate) const REPS: &[(&str, &str)] = &[
    ("Spacefill", "rep spacefill"),
    ("Ball & stick", "rep ball_and_stick"),
    ("Sticks", "rep sticks"),
    ("Lines", "rep lines"),
    ("SAS", "rep sas"),
    ("Tube", "rep tube"),
    ("Cartoon", "rep cartoon"),
    ("Gaussian surface", "rep gaussian_surface"),
    ("Skin surface", "rep skin_surface"),
    ("Molecular surface (SES)", "rep ses"),
    ("Glycan (3D-SNFG)", "rep glycan"),
];

/// Represent ▸ Add rep's form (UX_PATTERNS "Create"): style and
/// which atoms it draws, in one step and one undo entry -- replaces the
/// old Add-then-Atoms two-step flow. Coloring, material and options start
/// as the current rep's own (as the Structures panel's "Add rep" used
/// to), so only style and atoms need deciding up front.
fn represent_add_rep_popover(app: &mut AppUi<'_>, ui: &mut Ui) {
    ui.set_max_width(240.0);
    let id = egui::Id::new("represent-add-rep-form");
    let (mut style, mut atoms) = ui
        .data_mut(|d| d.get_temp::<(usize, String)>(id))
        .unwrap_or_else(|| (0, "all".to_owned()));
    ui.label("Style");
    let names: Vec<&str> = REPS.iter().map(|(name, _)| *name).collect();
    if let Some(i) = widgets::select(ui, "add-rep-style", &names, Some(style), "Style") {
        style = i;
    }
    ui.label("Draws");
    ui.add(
        egui::TextEdit::singleline(&mut atoms)
            .hint_text("all")
            .desired_width(f32::INFINITY),
    );
    ui.data_mut(|d| d.insert_temp(id, (style, atoms.clone())));
    ui.add_space(crate::theme::space::GAP);
    let Some(submit) = widgets::form_actions(ui, "Add rep") else {
        return;
    };
    if !submit {
        ui.data_mut(|d| d.remove::<(usize, String)>(id));
        ui.close();
        return;
    }
    let selection = if atoms.trim().is_empty() {
        "all".to_owned()
    } else {
        atoms.trim().to_owned()
    };
    if let Err(e) = vv_core::select::parse(&selection) {
        app.fail(format!("bad selection `{selection}`: {e}"));
        return;
    }
    let Some(id_) = app.current() else {
        app.fail("no structure to add a rep to".into());
        return;
    };
    let word = REPS[style]
        .1
        .strip_prefix("rep ")
        .expect("REPS commands start with `rep `");
    let representation =
        vv_scene::script::parse_representation(word).expect("REPS names are valid representations");
    let loaded = app.scene.structure(id_).expect("current");
    let rep = vv_scene::Rep {
        id: vv_scene::RepId(loaded.next_rep_id),
        representation,
        selection,
        ..loaded.rep().clone()
    };
    let index = loaded.reps.len();
    let label = REPS[style].0;
    app.dispatch(vv_scene::Command::AddRep {
        id: id_,
        index,
        rep,
    });
    app.undoable(&format!("Added a {label} rep."));
    ui.data_mut(|d| d.remove::<(usize, String)>(id));
    ui.close();
}

pub(crate) const COLORINGS: &[(&str, &str)] = &[
    // Atoms
    ("Element", "color element"),
    ("Carbons by chain", "color hetero"),
    ("B-factor", "color b_factor"),
    ("Occupancy", "color occupancy"),
    ("Molecule class", "color class"),
    // Chains
    ("Chain", "color chain"),
    ("Segment", "color segname"),
    ("Fragment", "color fragment"),
    ("Structure", "color structure"),
    ("Rainbow", "color rainbow"),
    // Residues
    ("Residue type", "color restype"),
    ("Residue name", "color resname"),
    ("Clustal", "color clustal"),
    ("Zappo", "color zappo"),
    ("Taylor", "color taylor"),
    ("Helix propensity", "color helix"),
    ("Strand propensity", "color strand"),
    ("Turn propensity", "color turn"),
    ("Buried index", "color buried"),
    ("Nucleotide", "color nucleotide"),
    ("Purine/pyrimidine", "color purinepyrimidine"),
    // Hydrophobicity
    ("Kyte-Doolittle", "color hydrophobicity"),
    ("Wimley-White", "color ww"),
];

/// Where each caption in [`COLORINGS`] starts (`coloring_menu_ui`'s group
/// headings); kept in sync with it by the `every_coloring_round_trips_
/// and_groups_cover_the_whole_list` test rather than by construction,
/// since a `const` can't easily carry both a flat list (every other
/// `COLORINGS` consumer wants that) and its grouping at once.
pub(crate) const COLORING_GROUPS: &[(&str, usize)] = &[
    ("Atoms", 0),
    ("Chains", 5),
    ("Residues", 10),
    ("Hydrophobicity", 21),
];

pub(crate) const MATERIALS: &[(&str, &str)] = &[
    ("Opaque", "material opaque"),
    ("Flat", "material flat"),
    ("Transparent", "material transparent"),
    ("Brushed metal", "material brushedmetal"),
    ("Diffuse", "material diffuse"),
    ("Faint", "material faint"),
    ("Clear glass", "material clearglass"),
    ("Tinted glass", "material tintedglass"),
    ("Frosted glass", "material frostedglass"),
    ("Glossy", "material glossy"),
    ("Hard plastic", "material hardplastic"),
    ("Soft metal", "material softmetal"),
    ("Steel", "material steel"),
    ("Translucent", "material translucent"),
    ("Inked", "material inked"),
    ("Inked gloss", "material inkedgloss"),
    ("Inked glass", "material inkedglass"),
    ("Textbook", "material textbook"),
    ("Polished", "material polished"),
    ("Chalk", "material chalk"),
    ("Chalk edge", "material chalkedge"),
    ("Bubbled glass", "material bubbledglass"),
    ("Hollow glass", "material hollowglass"),
    ("Mirror chrome", "material mirrorchrome"),
    ("Toon", "material toon"),
    ("Clay", "material clay"),
    ("Rim", "material rim"),
    ("Studio matte", "material studiomatte"),
    ("Lab default", "material labdefault"),
];

const LOOKS: &[(&str, &str)] = &[
    ("Standard", "style dark"),
    ("Matte", "style white"),
    ("Glossy", "style glossy"),
    ("Cel", "style flat"),
];

const PROJECTIONS: &[(&str, &str)] = &[
    ("Orthographic", "view projection orthographic"),
    ("Perspective", "view projection perspective"),
];

/// The ribbon's fixed top-right chrome: global
/// actions, not part of any tab, always in the same place. "Find
/// commands" isn't a scriptable command (it opens the mechanism commands
/// run through) so it isn't here; `quick_access_ui` draws it alongside.
pub const QUICK_ACCESS: &[Action] = &[
    run("Undo", icon::ARROW_COUNTER_CLOCKWISE, "Z", "undo"),
    run("Redo", icon::ARROW_CLOCKWISE, "Y", "redo"),
    run("Save session", icon::FLOPPY_DISK, "S", "savesession"),
    run("Preferences", icon::GEAR, "P", "preferences"),
];

pub const TABS: &[RibbonTab] = &[
    RibbonTab {
        label: "Home",
        tip: "H",
        shown: None,
        groups: &[
            Group {
                label: "Select",
                actions: &[inline(
                    "Select tool and pick level",
                    icon::SELECTION,
                    "S",
                    crate::home::select_commands,
                    crate::home::select_row,
                )],
            },
            Group {
                label: "By type",
                actions: &[inline(
                    "Quick select by type",
                    icon::SPIRAL,
                    "T",
                    crate::home::quick_commands,
                    crate::home::quick_row,
                )],
            },
            Group {
                label: "Modify",
                actions: &[inline(
                    "Modify the selection",
                    icon::SELECTION_INVERSE,
                    "M",
                    crate::home::modify_commands,
                    crate::home::modify_row,
                )],
            },
            Group {
                label: "Interface",
                actions: &[inline(
                    "Interface residues",
                    crate::home_interface::CHAINS_ICON,
                    "I",
                    crate::home_interface::commands,
                    crate::home_interface::interface_row,
                )],
            },
            Group {
                label: "View",
                actions: &[inline(
                    "Reset and screenshot",
                    icon::FRAME_CORNERS,
                    "R",
                    crate::home::view_commands,
                    crate::home::view_row,
                )],
            },
        ],
    },
    RibbonTab {
        label: "File",
        tip: "F",
        shown: None,
        groups: &[
            Group {
                label: "File",
                actions: &[
                    big("Open", icon::FOLDER_OPEN, "O", "open"),
                    run("Fetch from PDB…", icon::CLOUD_ARROW_DOWN, "F", "fetch"),
                    popover(
                        "Recent",
                        icon::CLOCK_COUNTER_CLOCKWISE,
                        "E",
                        "file.recent",
                        &["open"],
                        recent_popover,
                    ),
                    run("Open session…", icon::FILE_ARROW_UP, "L", "loadsession"),
                    run("Save session…", icon::FLOPPY_DISK, "S", "savesession"),
                    popover(
                        "Add trajectory",
                        icon::FILM_STRIP,
                        "T",
                        "file.trajectory",
                        &["loadtrajectory ", "attachtrajectory "],
                        crate::ui::file_trajectory_popover,
                    ),
                    run("Close structure", icon::X_SQUARE, "C", "close"),
                ],
            },
            Group {
                label: "Export",
                actions: &[
                    run("Image…", icon::CAMERA, "I", "screenshot"),
                    run("Movie…", icon::FILM_SLATE, "M", "panel movie"),
                    popover(
                        "Render",
                        icon::APERTURE,
                        "N",
                        "file.render",
                        &["render image.png 3840x2160 64"],
                        crate::ui::file_render_popover,
                    ),
                    run("Structure…", icon::FLOPPY_DISK_BACK, "U", "savestructure"),
                ],
            },
        ],
    },
    RibbonTab {
        label: "Represent",
        tip: "R",
        shown: None,
        groups: &[
            Group {
                label: "Style of the current rep",
                actions: &[choice("Style", icon::ATOM, "S", REPS, rep_style, false)],
            },
            Group {
                label: "Reps",
                actions: &[Action {
                    big: true,
                    ..popover(
                        "Add rep",
                        icon::STACK,
                        "A",
                        "represent.addrep",
                        &["addrep "],
                        represent_add_rep_popover,
                    )
                }],
            },
            Group {
                label: "Coloring",
                actions: &[choice_disabled_when(
                    "Coloring",
                    icon::PALETTE,
                    "G",
                    COLORINGS,
                    rep_coloring,
                    true,
                    rep_is_glycan,
                    "SNFG colours are fixed",
                )],
            },
            Group {
                label: "Custom",
                actions: &[
                    color("One color", icon::EYEDROPPER, "1", rep_color, "color "),
                    prompt("Values", icon::GRAPH, "V", "values "),
                ],
            },
            Group {
                label: "Material",
                actions: &[choice(
                    "Material",
                    icon::DROP,
                    "M",
                    MATERIALS,
                    rep_material,
                    true,
                )],
            },
        ],
    },
    RibbonTab {
        // Every look setting lives here (there is no Look panel), decoupled from materials and, for
        // the four effects with a strength, a switch + its own ▾.
        label: "Look",
        tip: "L",
        shown: None,
        groups: &[
            Group {
                label: "Preset",
                actions: &[choice("Preset", icon::SPARKLE, "Y", LOOKS, look, false)],
            },
            Group {
                label: "Lights",
                actions: &[popover(
                    "Lights",
                    icon::SUN,
                    "L",
                    "look.lights",
                    &["lighting default", "tonemap on"],
                    crate::studio::lights_popover,
                )],
            },
            Group {
                label: "Background",
                actions: &[popover(
                    "Background",
                    icon::PAINT_BRUSH,
                    "B",
                    "look.background",
                    &["background ", "gradient on"],
                    crate::ui::look_background_popover,
                )],
            },
            Group {
                label: "Outline",
                // `ViewSettings::outline` is a plain bool: a switch, no ▾.
                actions: &[toggle(
                    "Outline",
                    icon::POLYGON,
                    "O",
                    "outline on",
                    "outline off",
                    outline_on,
                )],
            },
        ],
    },
    RibbonTab {
        // List/Remove prompts and File info are gone: the Inspector's Annotations list (each row's own "x") and a
        // structure row's ⋯ cover them now.
        label: "Analyze",
        tip: "A",
        shown: None,
        groups: &[
            Group {
                label: "Measure",
                actions: &[Action {
                    big: true,
                    ..run("Measure", icon::RULER, "M", "measure")
                }],
            },
            Group {
                label: "Annotate",
                actions: &[
                    toggle(
                        "Label",
                        icon::TAG,
                        "L",
                        "mode label",
                        "mode rotate",
                        label_mode_on,
                    ),
                    popover(
                        "Caption",
                        icon::TEXT_T,
                        "T",
                        "analyze.caption",
                        &["caption "],
                        crate::ui::analyze_caption_popover,
                    ),
                ],
            },
            Group {
                label: "Interactions",
                actions: &[
                    toggle(
                        "H-bonds",
                        icon::LINK_SIMPLE,
                        "H",
                        "interactions hbond on",
                        "interactions hbond off",
                        hbond_on,
                    ),
                    toggle(
                        "Metal",
                        icon::ATOM,
                        "E",
                        "interactions metal on",
                        "interactions metal off",
                        metal_on,
                    ),
                    toggle(
                        "Salt bridges",
                        icon::LIGHTNING,
                        "B",
                        "interactions saltbridge on",
                        "interactions saltbridge off",
                        saltbridge_on,
                    ),
                ],
            },
            Group {
                label: "Reports",
                actions: &[
                    popover(
                        "Contacts",
                        icon::LINE_SEGMENTS,
                        "C",
                        "analyze.contacts",
                        &["contacts "],
                        crate::ui::analyze_contacts_popover,
                    ),
                    run("Surface area", icon::DROP, "S", "sasa"),
                ],
            },
        ],
    },
    RibbonTab {
        // Clip is a switch + ▾ (a `SwitchEffect`); Workspace and Panels are single ▾ popovers now.
        // Ctrl+4 covers `panel look`; Ctrl+9 is gone.
        label: "View",
        tip: "V",
        shown: None,
        groups: &[
            Group {
                label: "Camera",
                actions: &[
                    big("Reset view", icon::FRAME_CORNERS, "R", "view reset"),
                    choice(
                        "Projection",
                        icon::CUBE,
                        "P",
                        PROJECTIONS,
                        projection,
                        false,
                    ),
                    popover_form(
                        "Rotate",
                        icon::ARROWS_CLOCKWISE,
                        "O",
                        "view.rotate",
                        &["rotate x 15"],
                        crate::ui::view_rotate_popover,
                    ),
                ],
            },
            Group {
                label: "Clip",
                actions: &[switch_effect(
                    "Clip plane",
                    icon::KNIFE,
                    "C",
                    "clip on",
                    "clip off",
                    clip_on,
                    "view.clip",
                    crate::ui::view_clip_popover,
                )],
            },
            Group {
                label: "Workspace",
                actions: &[popover(
                    "Layout",
                    icon::SQUARES_FOUR,
                    "L",
                    "view.layout",
                    &[
                        "layout Default",
                        "savelayout ",
                        "startlayout",
                        "deletelayout",
                    ],
                    crate::ui::view_layout_popover,
                )],
            },
            Group {
                label: "Panels",
                actions: &[popover(
                    "Panels",
                    icon::DOTS_NINE,
                    "N",
                    "view.panels",
                    &["panel structures"],
                    crate::ui::view_panels_popover,
                )],
            },
        ],
    },
    RibbonTab {
        // Scene setup: the same Lights popover
        // as Look ▸ Lights (kept there too -- widening the light list is
        // a data-model change that touches that popover regardless, and
        // deleting it would strand `look lights` and existing muscle
        // memory), plus saved camera views and a render-framing preview,
        // neither of which fits Look's "how it's shaded" scope.
        label: "Studio",
        tip: "S",
        shown: None,
        groups: &[
            Group {
                label: "Lighting",
                actions: &[popover(
                    "Lights",
                    icon::SUN,
                    "L",
                    "studio.lights",
                    &["lighting default", "tonemap on"],
                    crate::studio::lights_popover,
                )],
            },
            Group {
                label: "Camera",
                actions: &[popover(
                    "Views",
                    icon::BOOKMARK_SIMPLE,
                    "V",
                    "studio.views",
                    &[],
                    crate::studio::views_popover,
                )],
            },
            Group {
                label: "Preview",
                actions: &[toggle(
                    "Frame preview",
                    icon::CROP,
                    "F",
                    "framing on",
                    "framing off",
                    framing_on,
                )],
            },
        ],
    },
];

/// Height of a tab's body: three rows of small actions.
pub(crate) const BODY_HEIGHT: f32 = 3.0 * crate::theme::CONTROL_HEIGHT;

/// Which tab is open, and the key-tip state.
#[derive(Default)]
pub struct RibbonState {
    pub tab: usize,
    /// Alt pressed: tab letters show; after one, that tab's action letters.
    pub key_tips: KeyTips,
    /// Alt was down last frame, and whether a key was pressed with it (a
    /// tap of Alt alone toggles key tips).
    pub alt_down: bool,
    pub alt_used: bool,
    /// A folded ribbon showing one tab until its next action runs.
    pub peek: bool,
    /// A popover a deep link (e.g. `look lights`) asked to force open, by
    /// its `Kind::Popover` key; consumed the next time that action draws.
    pub pending_popover: Option<&'static str>,
    /// Set by `panel selection`: the Selections panel's expression field
    /// takes keyboard focus the next time it draws.
    pub focus_expression: bool,
    /// Set by `structures menu`: the current structure row's ⋯ menu opens.
    pub open_structure_menu: bool,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum KeyTips {
    #[default]
    Off,
    Tabs,
    Actions,
}

impl RibbonState {
    /// Feeds one frame's Alt state and whether any key was pressed; true
    /// when Alt was just tapped (down and up with no key between).
    pub fn alt_tapped(&mut self, alt_now: bool, key_pressed: bool) -> bool {
        if alt_now && !self.alt_down {
            self.alt_used = false;
        }
        if alt_now && key_pressed {
            self.alt_used = true;
        }
        let tapped = !alt_now && self.alt_down && !self.alt_used;
        self.alt_down = alt_now;
        tapped
    }

    pub fn toggle_key_tips(&mut self) {
        self.key_tips = match self.key_tips {
            KeyTips::Off => KeyTips::Tabs,
            _ => KeyTips::Off,
        };
    }
}

/// What clicking an action does: a command line, or opening the palette.
enum Effect {
    Run(String),
    Prompt(&'static str),
}

/// Verbs that act on "the current structure" (`vv_scene::script::
/// current`, the most recently loaded one) when given no id.
const CURRENT_STRUCTURE_VERBS: [&str; 5] =
    ["rep", "representation", "color", "material", "interactions"];

/// Appends `id` (a structure's raw slot index, `StructureId::to_raw`) to
/// `line`'s end when its verb acts on "the current structure" and an
/// override is given, so a script command line built from a ribbon click
/// names the structure explicitly instead of falling through to whichever
/// one loaded last. Pure (no `AppUi`) so the verb family and formatting
/// are checkable directly.
fn append_current_structure(line: String, id: Option<u32>) -> String {
    let verb = line.split_whitespace().next().unwrap_or("");
    match (CURRENT_STRUCTURE_VERBS.contains(&verb), id) {
        (true, Some(id)) => format!("{line} {id}"),
        _ => line,
    }
}

/// Represent actions act on "the current structure"; this carries the
/// Structures panel's click-to-select override through to them, so the
/// tab and the row it acts on stay in sync (`AppUi::current`'s doc)
/// instead of the ribbon silently editing whichever structure loaded last
/// while the panel shows a different one as current.
fn with_current_override(line: String, app: &AppUi<'_>) -> String {
    let id = app
        .current_override
        .filter(|&id| app.scene.structure(id).is_some())
        .map(|id| id.to_raw());
    append_current_structure(line, id)
}

fn activate(action: &Action, app: &AppUi<'_>, option: Option<usize>) -> Option<Effect> {
    match &action.kind {
        Kind::Run(cmd) => Some(Effect::Run((*cmd).into())),
        Kind::Inline { commands, .. } => commands().first().map(|c| Effect::Run((*c).into())),
        Kind::Prompt(text) => Some(Effect::Prompt(text)),
        Kind::Toggle { on, off, state } | Kind::SwitchEffect { on, off, state, .. } => {
            Some(Effect::Run(if state(app) {
                (*off).into()
            } else {
                (*on).into()
            }))
        }
        Kind::Color { verb, .. } => Some(Effect::Prompt(verb)),
        Kind::Choice { options, .. } => option.map(|i| Effect::Run(options[i].1.into())),
        // Opened directly from the render loop (it needs a `Ui`), and
        // inert: neither reachable via a key tip.
        Kind::Popover { .. } | Kind::Disabled { .. } => None,
    }
}

/// The ribbon's label for the option whose command is `command`.
fn option_label(options: &[(&'static str, &str)], command: String) -> String {
    options
        .iter()
        .find(|(_, c)| *c == command)
        .map_or(command.clone(), |(label, _)| (*label).to_owned())
}

/// A rep's style as the ribbon names it: "Cartoon".
pub fn style_label(rep: &vv_scene::Rep) -> String {
    let command = format!(
        "rep {}",
        vv_scene::session::representation_name(rep.representation)
    );
    option_label(REPS, command)
}

/// A rep in words, as the ribbon names its style, coloring and material:
/// "Spacefill · Element · Opaque".
pub fn rep_summary(rep: &vv_scene::Rep) -> String {
    [
        style_label(rep),
        coloring_label(&rep.coloring),
        option_label(MATERIALS, format!("material {}", rep.material.name())),
    ]
    .join(" · ")
}

/// "Chain", or "Solid red" / "Solid #ff8000" for one colour.
pub fn coloring_label(coloring: &vv_scene::ColorScheme) -> String {
    match coloring {
        vv_scene::ColorScheme::Constant(_) => format!("Solid {}", coloring.name()),
        _ => option_label(COLORINGS, format!("color {}", coloring.name())),
    }
}

impl AppUi<'_> {
    /// Does what clicking `action` does; a choice steps to its next option.
    pub fn press(&mut self, action: &Action) {
        let option = match &action.kind {
            Kind::Choice {
                options, current, ..
            } => {
                let now = current(self);
                let at = options.iter().position(|(_, c)| Some(*c) == now.as_deref());
                Some(at.map_or(0, |i| (i + 1) % options.len()))
            }
            _ => None,
        };
        if let Some(effect) = activate(action, self, option) {
            self.apply_effect(effect);
        }
    }

    fn apply_effect(&mut self, effect: Effect) {
        match effect {
            Effect::Run(line) => {
                let line = with_current_override(line, self);
                self.run_command_logged(&line);
            }
            Effect::Prompt(text) => {
                self.palette.open = true;
                self.palette.query = text.into();
            }
        }
    }

    /// The tabs showing now (contextual ones only while they apply).
    fn visible_tabs(&self) -> Vec<usize> {
        (0..TABS.len())
            .filter(|&i| TABS[i].shown.is_none_or(|q| q(self)))
            .collect()
    }

    pub fn ribbon_ui(&mut self, ui: &mut Ui) {
        // Below this width the ribbon folds to its tab row, so the panels
        // and viewport keep their room.
        const NARROW: f32 = 1280.0;
        let narrow = ui.ctx().content_rect().width() < NARROW;
        let folded = (self.prefs.ribbon_collapsed || narrow) && !self.ribbon.peek;
        let visible = self.visible_tabs();
        if !visible.contains(&self.ribbon.tab) {
            self.ribbon.tab = 0;
        }
        let tips = self.ribbon.key_tips;
        ui.horizontal(|ui| {
            for &i in &visible {
                let tab = &TABS[i];
                let label = if tips == KeyTips::Tabs {
                    format!("{}  [{}]", tab.label, tab.tip)
                } else {
                    tab.label.into()
                };
                if widgets::tab(ui, &label, self.ribbon.tab == i, tab.shown.is_some()).clicked() {
                    // Folded, a click shows the tab until an action runs;
                    // the same tab again folds it.
                    self.ribbon.peek = !(self.ribbon.peek && self.ribbon.tab == i);
                    self.ribbon.tab = i;
                    if !narrow {
                        self.prefs.ribbon_collapsed = false;
                    }
                }
            }
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                let (chevron, tip) = if folded {
                    (icon::CARET_DOWN, "Show the ribbon")
                } else {
                    (icon::CARET_UP, "Fold the ribbon to its tab names")
                };
                if widgets::icon_button(ui, chevron, false)
                    .on_hover_text(tip)
                    .clicked()
                {
                    if narrow {
                        self.ribbon.peek = folded;
                    } else {
                        self.prefs.ribbon_collapsed = !self.prefs.ribbon_collapsed;
                    }
                }
                ui.separator();
                self.quick_access_ui(ui);
            });
        });
        if folded {
            return;
        }
        ui.separator();
        let mut effect: Option<Effect> = None;
        let tab = &TABS[self.ribbon.tab];
        egui::ScrollArea::horizontal()
            .id_salt("ribbon groups")
            .show(ui, |ui| {
                ui.horizontal(|ui| {
                    for group in tab.groups {
                        effect = self.group_ui(ui, group, tips).or(effect.take());
                        ui.separator();
                    }
                });
            });
        if let Some(e) = effect {
            self.apply_effect(e);
            self.ribbon.peek = false;
        }
    }

    /// One ribbon group: its big buttons, then its small actions in
    /// columns of three (or, for compact rows, one centered column of
    /// rows), with the group's name under them.
    fn group_ui(&mut self, ui: &mut Ui, group: &Group, tips: KeyTips) -> Option<Effect> {
        let mut effect = None;
        ui.vertical(|ui| {
            ui.horizontal(|ui| {
                // One body height for every tab, so switching tabs never
                // moves the panels below.
                ui.set_min_height(BODY_HEIGHT);
                for action in group.actions.iter().filter(|a| a.big) {
                    effect = self.big_button(ui, action, tips).or(effect.take());
                }
                let (rows, small): (Vec<&Action>, Vec<&Action>) = group
                    .actions
                    .iter()
                    .filter(|a| !a.big)
                    .partition(|a| is_row(a));
                for column in small.chunks(3) {
                    ui.vertical(|ui| {
                        ui.spacing_mut().item_spacing.y = 0.0;
                        for action in column {
                            // A choice among other actions says what it
                            // chooses ("Supersample: Off 1.5x 2x").
                            let named = group.actions.len() > 1;
                            effect = self.small_action(ui, action, tips, named).or(effect.take());
                        }
                    });
                }
                if !rows.is_empty() {
                    effect = self.row_column(ui, &rows, tips).or(effect.take());
                }
            });
            widgets::caption(ui, group.label);
        });
        effect
    }

    /// Compact rows (`is_row`) one above the other.
    fn row_column(&mut self, ui: &mut Ui, rows: &[&Action], tips: KeyTips) -> Option<Effect> {
        let mut effect = None;
        ui.vertical(|ui| {
            ui.spacing_mut().item_spacing.y = crate::theme::space::TIGHT;
            for action in rows {
                effect = self.small_action(ui, action, tips, false).or(effect.take());
            }
        });
        effect
    }

    /// The fixed top-right chrome: Undo, Redo, Save
    /// session, Find commands, Preferences, left to right, icon-only so
    /// they never crowd out the tabs. Drawn under a right-to-left layout
    /// (the ribbon's outer row), which places each added button to the
    /// *left* of the one before it -- so this walks `QUICK_ACCESS`
    /// reversed (Preferences first) to land in the order above, splicing
    /// in Find commands (not a scriptable `Action`, so not in the array)
    /// right after Preferences, its neighbour in that order.
    fn quick_access_ui(&mut self, ui: &mut Ui) {
        let mut effect = None;
        let mut find_clicked = false;
        for (i, action) in QUICK_ACCESS.iter().enumerate().rev() {
            let enabled = self.action_enabled(action);
            let r = ui
                .add_enabled_ui(enabled, |ui| widgets::icon_button(ui, action.icon, false))
                .inner;
            if with_help(r, action).clicked() {
                effect = activate(action, self, None);
            }
            if i == QUICK_ACCESS.len() - 1
                && widgets::icon_button(ui, icon::MAGNIFYING_GLASS, false)
                    .on_hover_text("Find commands (Ctrl+P)")
                    .clicked()
            {
                find_clicked = true;
            }
        }
        if find_clicked {
            self.palette.open = true;
        }
        if let Some(e) = effect {
            self.apply_effect(e);
        }
    }

    /// Undo/Redo grey out once neither history (scene or view/camera) has
    /// anything left in that direction.
    fn action_enabled(&self, action: &Action) -> bool {
        match action.label {
            "Undo" => self.history.can_undo() || self.app_history.can_undo(),
            "Redo" => self.history.can_redo() || self.app_history.can_redo(),
            _ => true,
        }
    }

    /// Opens/anchors a `Kind::Popover`'s body on `toggle`'s response --
    /// shared by `big_button` (a group's one big action, e.g. Represent's
    /// Add rep), `small_action`'s own `Kind::Popover` and `SwitchEffect`'s
    /// ▾.
    fn popover_ui(
        &mut self,
        ui: &Ui,
        toggle: &egui::Response,
        key: &'static str,
        body: PopoverBody,
    ) {
        let popup_id = egui::Id::new(("ribbon-popover", key));
        if self.ribbon.pending_popover == Some(key) {
            egui::Popup::open_id(ui.ctx(), popup_id);
            self.ribbon.pending_popover = None;
        }
        egui::Popup::from_toggle_button_response(toggle)
            .id(popup_id)
            .close_behavior(egui::PopupCloseBehavior::CloseOnClickOutside)
            .show(|ui| body(self, ui));
    }

    fn big_button(&mut self, ui: &mut Ui, action: &Action, tips: KeyTips) -> Option<Effect> {
        let label = tip_label(action, tips);
        if let Kind::Popover { key, body, .. } = &action.kind {
            let (key, body) = (*key, *body);
            let r = with_help(
                widgets::large_button(ui, action.icon, &label, false),
                action,
            );
            self.popover_ui(ui, &r, key, body);
            return None;
        }
        let selected = matches!(&action.kind, Kind::Toggle { state, .. } if state(self));
        with_help(
            widgets::large_button(ui, action.icon, &label, selected),
            action,
        )
        .clicked()
        .then(|| activate(action, self, None))
        .flatten()
    }

    fn small_action(
        &mut self,
        ui: &mut Ui,
        action: &Action,
        tips: KeyTips,
        named: bool,
    ) -> Option<Effect> {
        let label = tip_label(action, tips);
        match &action.kind {
            Kind::Choice {
                options,
                current,
                dropdown,
                disabled_when,
            } => {
                if let Some((disabled, reason)) = disabled_when {
                    if disabled(self) {
                        ui.add_enabled_ui(false, |ui| {
                            widgets::button(ui, action.icon, &label, Variant::Ghost)
                        })
                        .inner
                        .on_disabled_hover_text(*reason);
                        return None;
                    }
                }
                let now = current(self);
                let selected = options.iter().position(|(_, c)| Some(*c) == now.as_deref());
                let names: Vec<&str> = options.iter().map(|o| o.0).collect();
                let picked = if *dropdown {
                    widgets::select(ui, action.label, &names, selected, action.label)
                } else {
                    ui.horizontal(|ui| {
                        ui.set_max_width(440.0);
                        if named || tips == KeyTips::Actions {
                            widgets::caption(ui, &label);
                        }
                        widgets::segmented(ui, &names, selected)
                    })
                    .inner
                };
                picked.and_then(|i| activate(action, self, Some(i)))
            }
            Kind::Toggle { state, .. } => {
                let on = state(self);
                with_help(widgets::switch(ui, on, &label), action)
                    .clicked()
                    .then(|| activate(action, self, None))
                    .flatten()
            }
            Kind::Color { get, verb } => {
                let mut picked = get(self);
                ui.horizontal(|ui| {
                    let changed = widgets::color_picker(ui, &mut picked);
                    ui.label(&label);
                    changed.then(|| Effect::Run(format!("{verb}{}", widgets::hex(picked))))
                })
                .inner
            }
            Kind::Run(_) => {
                let r = widgets::button(ui, action.icon, &label, Variant::Ghost);
                with_help(r, action)
                    .clicked()
                    .then(|| activate(action, self, None))
                    .flatten()
            }
            // An ellipsis: it asks for more input before it acts.
            Kind::Prompt(_) => {
                let r = widgets::button(ui, action.icon, &format!("{label}…"), Variant::Ghost);
                with_help(r, action)
                    .clicked()
                    .then(|| activate(action, self, None))
                    .flatten()
            }
            Kind::Popover {
                key,
                body,
                ellipsis,
                ..
            } => {
                let (key, body) = (*key, *body);
                let r = if *ellipsis {
                    widgets::button(ui, action.icon, &format!("{label}…"), Variant::Ghost)
                } else {
                    widgets::disclosure_button(ui, action.icon, &label)
                };
                with_help(r.clone(), action);
                self.popover_ui(ui, &r, key, body);
                None
            }
            Kind::Inline { body, .. } => {
                let body = *body;
                ui.horizontal(|ui| {
                    if tips == KeyTips::Actions {
                        widgets::caption(ui, &format!("[{}]", action.tip));
                    }
                    body(self, ui)
                })
                .inner
                .map(Effect::Run)
            }
            Kind::Disabled { reason } => {
                ui.add_enabled_ui(false, |ui| {
                    widgets::button(ui, action.icon, &label, Variant::Ghost)
                })
                .inner
                .on_disabled_hover_text(*reason);
                None
            }
            Kind::SwitchEffect { key, body, .. } => {
                let (key, body) = (*key, *body);
                let on = matches!(&action.kind, Kind::SwitchEffect { state, .. } if state(self));
                let (sw, more) = widgets::switch_disclosure(ui, on, &label);
                let effect = with_help(sw, action)
                    .clicked()
                    .then(|| activate(action, self, None))
                    .flatten();
                if let Some(more) = more {
                    self.popover_ui(ui, &more, key, body);
                }
                effect
            }
        }
    }
}

/// Whether an action draws as a compact row instead of a stacked button.
fn is_row(action: &Action) -> bool {
    matches!(action.kind, Kind::Inline { .. })
}

/// The command an action runs first, whose help is its tooltip.
fn first_command(action: &Action) -> &'static str {
    match &action.kind {
        Kind::Run(c) | Kind::Prompt(c) | Kind::Color { verb: c, .. } => c,
        Kind::Inline { commands, .. } => commands().first().copied().unwrap_or(""),
        Kind::Toggle { on, .. } | Kind::SwitchEffect { on, .. } => on,
        Kind::Choice { options, .. } => options[0].1,
        // A popover with no scriptable equivalent (Studio ▸ Views ▾:
        // camera bookmarks have no console verb) falls back to its own
        // label, same as `Kind::Disabled`.
        Kind::Popover { commands, .. } => commands.first().copied().unwrap_or(""),
        Kind::Disabled { .. } => "",
    }
}

/// Hovering an action shows what its command does, and its shortcut if
/// it has one.
fn with_help(response: egui::Response, action: &Action) -> egui::Response {
    let command = first_command(action);
    response.on_hover_ui(|ui| {
        ui.set_max_width(320.0);
        let help = crate::commands::help_for(command).unwrap_or(action.label);
        match crate::keys::shortcut_for(command) {
            Some(shortcut) => ui.label(format!("{help}  ({shortcut})")),
            None => ui.label(help),
        };
    })
}

fn tip_label(action: &Action, tips: KeyTips) -> String {
    if tips == KeyTips::Actions {
        format!("{}  [{}]", action.label, action.tip)
    } else {
        action.label.into()
    }
}

/// The tab a `ribbon NAME` word names: its label, lowercase, from its
/// first word (`color` for "Color & Material").
pub fn tab_named(word: &str) -> Option<usize> {
    TABS.iter().position(|t| {
        t.label
            .split_whitespace()
            .next()
            .is_some_and(|first| first.eq_ignore_ascii_case(word))
    })
}

/// The dock panel a `panel NAME` word names.
/// A real, dockable panel `panel NAME` can open; "settings"/"look" are
/// handled before this is called (`commands::run_command`) since they now
/// name a ribbon tab, not a panel. "scene" is kept as an alias of
/// "structures" (the panel's current name) for old scripts.
pub fn panel_named(word: &str) -> Option<Tab> {
    Some(match word {
        "viewport" => Tab::Viewport,
        "structures" | "scene" => Tab::Structures,
        "selection" => Tab::Selection,
        "inspector" => Tab::Inspector,
        "log" => Tab::Log,
        "sequence" => Tab::Sequence,
        "timeline" => Tab::Timeline,
        "movie" => Tab::Movie,
        "info" => Tab::Info,
        _ => return None,
    })
}

/// Inverse of `panel_named`'s canonical word (never the "scene" alias):
/// item 5c, so View ▸ Panels ▾'s rows can show each one's own `Ctrl+N`
/// (`keys::shortcut_for(&format!("panel {word}"))`).
pub fn panel_word(tab: Tab) -> &'static str {
    match tab {
        Tab::Viewport => "viewport",
        Tab::Structures => "structures",
        Tab::Selection => "selection",
        Tab::Inspector => "inspector",
        Tab::Log => "log",
        Tab::Sequence => "sequence",
        Tab::SceneSettings => "settings",
        Tab::Timeline => "timeline",
        Tab::Movie => "movie",
        Tab::Info => "info",
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::{BTreeSet, HashSet};

    /// Commands reached some other way than the ribbon, and why.
    const CONSOLE_ONLY: &[(&str, &str)] = &[
        ("help", "the palette and console are its interface"),
        (
            "version",
            "for scripts and clients checking what they talk to",
        ),
        ("quit", "closing the window"),
        (
            "confirmquit",
            "the quit dialog's own checkbox sets it; rare enough not to need a ribbon button",
        ),
        (
            "welcome",
            "the start card's own checkbox sets it; `open`/`fetch`/`loadsession` cover the rest",
        ),
        (
            "info",
            "a structure row's ⋯ ▸ File info (Info panel) is its clickable form",
        ),
        (
            "label",
            "Analyze ▸ Label's click-atoms mode is its clickable form; a specific atom and \
             text (scripting) stay console/`--exec`/Python-only",
        ),
        (
            "unlabel",
            "the Inspector's Annotations rows' own x is its clickable form",
        ),
        (
            "labels",
            "the Inspector's Annotations list is its clickable form",
        ),
        (
            "uncaption",
            "the Inspector's Annotations rows' own x is its clickable form",
        ),
        (
            "captions",
            "the Inspector's Annotations list is its clickable form",
        ),
        (
            "unmeasure",
            "the Inspector's Annotations rows' own x is its clickable form",
        ),
        (
            "measurements",
            "the Inspector's Annotations list is its clickable form",
        ),
        (
            "movie",
            "the Movie panel (File ▸ Export ▸ Movie…) is its clickable form",
        ),
        (
            "sequence",
            "the Sequence panel's own header (Color, Tracks)",
        ),
        ("uishot", "a documentation and review tool"),
        ("window", "sets the window size for scripted layout checks"),
        ("repopt", "a Selections row's ⋯ ▸ Options…"),
        (
            "ribbon",
            "the ribbon's own tabs and chevron are its clickable form",
        ),
        (
            "look",
            "the Look tab's own ▾ buttons are its clickable form",
        ),
        ("mode", "the tool strip's mouse-mode buttons and R T S C"),
        ("ao", "Preferences ▸ Effects"),
        ("shadows", "Preferences ▸ Effects"),
        ("depthcue", "Preferences ▸ Effects"),
        ("dof", "Preferences ▸ Effects"),
        ("load", "Open covers it (`open PATH` loads and frames)"),
        (
            "select",
            "clicking atoms, Home's By type, Modify and Interface groups, or a Selections row's ⋯ ▸ Select these atoms",
        ),
        ("clear", "Escape in the viewport"),
        ("saveset", "the Selections panel's Save set"),
        ("useset", "a saved set's row in the Selections panel"),
        ("deleteset", "a saved set row's own x in the Selections panel"),
        ("sets", "the Selections panel lists the saved sets"),
        (
            "structures",
            "the Structures panel's rows are its clickable form",
        ),
        (
            "showstructure",
            "the Structures panel row's own eye is its clickable form",
        ),
        ("altloc", "a structure row's ⋯ ▸ Conformer"),
        ("selrep", "a Selections row's ⋯ ▸ Edit expression"),
        ("currep", "clicking a row in the Selections panel"),
        ("showrep", "a selection row's own eye is its clickable form"),
        (
            "delrep",
            "a selection row's own ⋯ ▸ Delete is its clickable form",
        ),
        ("reps", "the Selections panel's rows"),
        (
            "play",
            "the Timeline panel's play button and Space are its clickable form",
        ),
        (
            "frame",
            "the Timeline panel's frame slider is its clickable form",
        ),
        (
            "theme",
            "Preferences > Only me sets it directly (a live field edit, not a run command)",
        ),
        (
            "fxaa",
            "Preferences > Performance sets it directly (a live field edit, not a run command)",
        ),
        (
            "supersample",
            "Preferences > Performance sets it directly (a live field edit, not a run command)",
        ),
        (
            "occlusion",
            "Preferences > Performance sets it directly (a live field edit, not a run command)",
        ),
        (
            "adaptive",
            "Preferences > Performance sets it directly (a live field edit, not a run command)",
        ),
    ];

    fn commands_of(action: &Action) -> Vec<&'static str> {
        match &action.kind {
            Kind::Run(c) | Kind::Prompt(c) | Kind::Color { verb: c, .. } => {
                vec![c]
            }
            Kind::Inline { commands, .. } => commands(),
            Kind::Toggle { on, off, .. } | Kind::SwitchEffect { on, off, .. } => vec![on, off],
            Kind::Choice { options, .. } => options.iter().map(|(_, c)| *c).collect(),
            Kind::Popover { commands, .. } => commands.to_vec(),
            Kind::Disabled { .. } => vec![],
        }
    }

    fn all_actions() -> Vec<&'static Action> {
        TABS.iter()
            .flat_map(|t| t.groups)
            .flat_map(|g| g.actions)
            .chain(QUICK_ACCESS)
            .collect()
    }

    fn verb(line: &str) -> &str {
        line.split_whitespace().next().unwrap_or("")
    }

    #[test]
    fn append_current_structure_only_touches_represent_verbs() {
        assert_eq!(
            append_current_structure("rep cartoon".into(), Some(3)),
            "rep cartoon 3"
        );
        assert_eq!(
            append_current_structure("color chain".into(), Some(7)),
            "color chain 7"
        );
        assert_eq!(
            append_current_structure("material glass1".into(), None),
            "material glass1",
            "no override: line unchanged"
        );
        assert_eq!(
            append_current_structure("view reset".into(), Some(2)),
            "view reset",
            "a verb that isn't structure-scoped is never touched"
        );
        assert_eq!(append_current_structure("undo".into(), Some(2)), "undo");
    }

    #[test]
    fn every_ribbon_action_is_a_known_command() {
        for action in all_actions() {
            for line in commands_of(action) {
                crate::commands::validate(line).unwrap_or_else(|e| panic!("{}: {e}", action.label));
            }
        }
    }

    #[test]
    fn every_material_and_rep_command_round_trips_through_its_name() {
        for &(label, command) in MATERIALS {
            let name = command
                .strip_prefix("material ")
                .expect("a material command");
            let preset = vv_scene::Material::parse(name)
                .unwrap_or_else(|| panic!("{label}: {command}: not a known material"));
            assert_eq!(
                option_label(MATERIALS, format!("material {}", preset.name())),
                label,
                "{command}: name() should map back to its own ribbon label"
            );
        }
        for &(label, command) in REPS {
            let name = command.strip_prefix("rep ").expect("a rep command");
            let rep = vv_scene::script::parse_representation(name)
                .unwrap_or_else(|_| panic!("{label}: {command}: not a known representation"));
            assert_eq!(
                option_label(
                    REPS,
                    format!("rep {}", vv_scene::session::representation_name(rep))
                ),
                label,
                "{command}: representation_name() should map back to its own ribbon label"
            );
        }
    }

    #[test]
    fn every_coloring_round_trips_and_groups_cover_the_whole_list() {
        for &(_, command) in COLORINGS {
            let name = command.strip_prefix("color ").expect("a color command");
            let scheme = vv_scene::ColorScheme::parse(name)
                .unwrap_or_else(|| panic!("{command}: not a known color scheme"));
            assert_eq!(
                vv_scene::ColorScheme::parse(&scheme.name()),
                Some(scheme),
                "{command}: name() should round-trip through parse()"
            );
        }
        let starts: Vec<usize> = COLORING_GROUPS.iter().map(|(_, i)| *i).collect();
        assert_eq!(starts.first(), Some(&0), "the first group covers index 0");
        assert!(
            starts.windows(2).all(|w| w[0] < w[1]),
            "group starts must strictly increase: {starts:?}"
        );
        assert!(
            *starts.last().unwrap() < COLORINGS.len(),
            "every group start must be a real COLORINGS index"
        );
    }

    #[test]
    fn key_tips_are_unique_within_each_tab_and_across_tabs() {
        let mut tab_tips = HashSet::new();
        for tab in TABS {
            assert!(tab_tips.insert(tab.tip), "tab tip {} repeats", tab.tip);
            let mut tips = HashSet::new();
            for action in tab.groups.iter().flat_map(|g| g.actions) {
                assert!(
                    tips.insert(action.tip),
                    "{}: key tip {} repeats",
                    tab.label,
                    action.tip
                );
            }
        }
    }

    #[test]
    fn every_command_is_on_the_ribbon_or_explained() {
        let on_ribbon: BTreeSet<&str> = all_actions()
            .into_iter()
            .flat_map(commands_of)
            .map(verb)
            .collect();
        let explained: HashSet<&str> = CONSOLE_ONLY.iter().map(|(id, _)| *id).collect();
        let missing: Vec<&str> = crate::commands::entries()
            .iter()
            .map(|e| e.id)
            .filter(|id| {
                let alias = match *id {
                    "representation" => "rep",
                    other => other,
                };
                !on_ribbon.contains(id) && !on_ribbon.contains(alias) && !explained.contains(id)
            })
            .collect();
        assert!(missing.is_empty(), "not on the ribbon: {missing:?}");
    }
}
