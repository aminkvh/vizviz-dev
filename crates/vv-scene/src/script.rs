//! Text commands over the command bus: `load 4HHB.cif`, `select chain A
//! and name CA`, `color chain`, `undo`. One parser serves the app's
//! `--exec` flag, its command palette and Log console, and Python's
//! `Session.exec()`, so a script means the same thing everywhere.
//!
//! Only *document* verbs live here: things that change the `Scene` and
//! go through undo history. View settings (style, FXAA, camera,
//! screenshots) are the app's business and are layered on top in
//! `vv-app`; a headless `Session` can't offer them, so it shouldn't
//! advertise them. `SPECS` is the registry: it feeds the palette, `help`,
//! and later an agent-readable manifest (docs/COMMANDS.md).

use std::path::Path;
use std::sync::Arc;

use crate::command::{Command, SceneError};
use crate::history::CommandHistory;
use crate::scene::{
    Caption, ColorScheme, LoadedStructure, Material, Measurement, Rep, RepId, Representation,
    Scene, StructureId,
};
use crate::values::ValueChannel;

/// One registered command. `usage` is what `help` prints; `keywords`
/// widen fuzzy matching in the palette.
#[derive(Clone, Copy, Debug)]
pub struct Spec {
    pub id: &'static str,
    pub title: &'static str,
    pub keywords: &'static [&'static str],
    pub usage: &'static str,
    pub help: &'static str,
}

pub const SPECS: &[Spec] = &[
    Spec {
        id: "load",
        title: "Load structure",
        keywords: &["open", "file", "read"],
        usage: "load PATH",
        help: "Load an mmCIF or PDB file (optionally .gz). It becomes the current structure.",
    },
    Spec {
        id: "loadtrajectory",
        title: "Load trajectory",
        keywords: &[
            "dcd", "xtc", "trr", "netcdf", "psf", "prmtop", "charmm", "amber", "trajectory",
            "trj", "frames", "md",
        ],
        usage: "loadtrajectory TOPOLOGY TRAJECTORY",
        help: "Load a topology file (mmCIF, PDB, CHARMM/NAMD/X-PLOR PSF, or AMBER PRMTOP) \
               plus a DCD, XTC, TRR, AMBER NetCDF, or PDB/mmCIF (one frame per model) \
               trajectory as one multi-frame structure; scrub frames in the Timeline tab. A \
               PSF/PRMTOP's own bond list is used as-is. If TOPOLOGY is already open, the \
               trajectory is attached to that structure instead (see attachtrajectory). \
               Neither path may contain spaces.",
    },
    Spec {
        id: "attachtrajectory",
        title: "Attach trajectory",
        keywords: &["dcd", "xtc", "trr", "netcdf", "trajectory", "trj", "frames", "md"],
        usage: "attachtrajectory TRAJECTORY [STRUCTURE]",
        help: "Give an open structure (default: the current one) the frames of a DCD, XTC, TRR, \
               AMBER NetCDF, or PDB/mmCIF trajectory, keeping its reps, selections and \
               labels. The trajectory must have the structure's atom count; undoable. The \
               path may not contain spaces.",
    },
    Spec {
        id: "fetch",
        title: "Fetch from the PDB",
        keywords: &["download", "pdb", "rcsb", "id", "assembly", "biological unit"],
        usage: "fetch ID [bcif | assembly N]",
        help: "Download an entry by PDB ID into the per-user cache and load it: the asymmetric unit (as BinaryCIF with `bcif`), or biological assembly N.",
    },
    Spec {
        id: "close",
        title: "Close structure",
        keywords: &["remove", "unload"],
        usage: "close [ID]",
        help: "Close a structure (default: the current one). Undoable.",
    },
    Spec {
        id: "structures",
        title: "List structures",
        keywords: &["list", "loaded", "ids"],
        usage: "structures",
        help: "List loaded structures with their ids.",
    },
    Spec {
        id: "showstructure",
        title: "Show or hide structure",
        keywords: &["hide", "visible", "toggle", "eye"],
        usage: "showstructure on|off [ID]",
        help: "Show or hide a whole structure (default: the current one), independent of any \
               rep's own `showrep`. A hidden structure draws nothing and is skipped during \
               trajectory playback: showing it again catches up to the current frame at once.",
    },
    Spec {
        id: "altloc",
        title: "Alternate locations",
        keywords: &["conformer", "alt", "disorder", "occupancy"],
        usage: "altloc first|all|LABEL [ID]",
        help: "Choose which alternate location of each residue reps draw (and build bonds and \
               surfaces from): `first`, the default, the conformer with the largest occupancy; \
               `all`; or one label such as `A` (a residue without it shows its `first`). Every \
               conformer stays in the data and is still reachable with `altloc A` selections.",
    },
    Spec {
        id: "select",
        title: "Select atoms",
        keywords: &["expression", "query", "atoms"],
        usage: "select EXPR | add EXPR | remove EXPR | invert | expand residue|chain|molecule|within N | grow [N] | shrink [N|within N] | interface A [to|with B] [within N]",
        help: "Select atoms of the current structure matching EXPR (docs/SELECTION.md). \
               `add` keeps the current selection and adds EXPR, `remove` takes EXPR out of it; `invert` selects what is not selected. \
               `expand` grows the selection to whole residues, chains or molecules, or to everything \
               within N Å; `grow`/`shrink` move each selected run N residues along its chain; \
               `shrink within N` peels off atoms within N Å of the unselected. `interface A to B` \
               selects the residues of A within N Å of B (`with` adds B's side; no B means any \
               other polymer; N defaults to 5).",
    },
    Spec {
        id: "clear",
        title: "Clear selection",
        keywords: &["deselect", "none"],
        usage: "clear",
        help: "Clear the active selection.",
    },
    Spec {
        id: "saveset",
        title: "Save selection set",
        keywords: &["named", "selection", "store"],
        usage: "saveset NAME",
        help: "Save the active selection under NAME (overwrites an existing set).",
    },
    Spec {
        id: "useset",
        title: "Use selection set",
        keywords: &["recall", "named", "selection"],
        usage: "useset NAME",
        help: "Make the saved set NAME the active selection.",
    },
    Spec {
        id: "deleteset",
        title: "Delete selection set",
        keywords: &["remove", "named", "selection"],
        usage: "deleteset NAME",
        help: "Delete the saved set NAME. Undoable.",
    },
    Spec {
        id: "sets",
        title: "List selection sets",
        keywords: &["named", "selections"],
        usage: "sets",
        help: "List saved selection sets.",
    },
    Spec {
        id: "representation",
        title: "Set representation",
        keywords: &[
            "rep", "spacefill", "ballstick", "ball-and-stick", "sticks", "tube", "backbone",
            "cartoon", "ribbon", "surface", "gaussian", "blob", "blobby", "skin", "skinsurface",
            "ses", "glycan", "snfg", "carbohydrate",
        ],
        usage: "representation spacefill|ballstick|sticks|lines|sas|tube|cartoon|gaussiansurface|skinsurface|ses|glycan [ID]",
        help: "Draw the current rep of a structure (default: current) as spacefill, \
               ball-and-stick, a backbone tube, DSSP-shaped cartoon ribbons (`vv_core::dssp`/`cartoon`), a ray-marched \
               Gaussian surface (`vv_core::gaussian_surface`), a skin surface \
               (`vv_core::skin_surface`, Edelsbrunner's mixed complex), the analytic solvent-excluded \
               surface (`ses`), or 3D-SNFG glycan symbols (`glycan`, `vv_core::glycan`). `rep` works too.",
    },
    Spec {
        id: "addrep",
        title: "Add representation",
        keywords: &["rep", "representation", "layer", "new", "ligand", "pocket"],
        usage: "addrep STYLE [SELECTION]",
        help: "Draw part of the current structure another way: a new rep with STYLE (as \
               `representation`) for the atoms SELECTION picks (default all; \
               docs/SELECTION.md), coloured and made like the current rep, and made current. \
               e.g. `addrep ballstick hetero and not water`.",
    },
    Spec {
        id: "delrep",
        title: "Delete representation",
        keywords: &["remove", "rep", "layer"],
        usage: "delrep [N]",
        help: "Delete rep N (default: the current one) of the current structure. A structure \
               keeps at least one rep; `showrep N off` hides it instead. Undoable.",
    },
    Spec {
        id: "selrep",
        title: "Set representation selection",
        keywords: &["select", "rep", "atoms", "subset"],
        usage: "selrep SELECTION",
        help: "Which atoms the current rep draws (docs/SELECTION.md; `all` for every atom).",
    },
    Spec {
        id: "currep",
        title: "Current representation",
        keywords: &["rep", "edit", "which"],
        usage: "currep N",
        help: "Make rep N the one `representation`, `color`, `material` and `selrep` edit.",
    },
    Spec {
        id: "repopt",
        title: "Representation option",
        keywords: &[
            "size", "radius", "probe", "scale", "bond", "blob", "shrink", "tune", "bases",
            "ligands", "ions", "water", "glycans", "lipids", "additives",
        ],
        usage: "repopt [NAME VALUE|default]",
        help: "Tune the current rep: atom `scale`, `bond` or tube `radius` in Angstrom, the \
               SAS/SES `probe` radius, Gaussian `blob`, skin `shrink`. A cartoon or tube \
               also takes `ligands`, `ions`, `glycans`, `lipids`, `water`, `additives` \
               (on|off), and a cartoon `bases` (plate|ladder|stick). Alone, lists the \
               options of the current rep's style.",
    },
    Spec {
        id: "showrep",
        title: "Show or hide representation",
        keywords: &["hide", "visible", "rep", "toggle"],
        usage: "showrep N on|off",
        help: "Show or hide rep N of the current structure.",
    },
    Spec {
        id: "reps",
        title: "List representations",
        keywords: &["rep", "layers", "list"],
        usage: "reps",
        help: "List the current structure's reps: index, style, coloring, material, \
               selection; `*` marks the current one.",
    },
    Spec {
        id: "material",
        title: "Set material",
        keywords: &[
            "opaque", "flat", "glossy", "chalk", "textbook", "glass", "transparent", "inked",
            "steel", "metal", "toon", "shiny", "specular", "clay", "rim", "studiomatte",
            "labdefault",
        ],
        usage: "material NAME [ID]",
        help: "Set what the current rep of a structure (default: current) is made of: how it reflects the \
               scene's lighting (`lighting`). Opaque, flat (no shading -- the shading term is \
               exactly the base colour; the viewport's tonemap/ao/shadows/depthcue still apply \
               on top, as for any material), transparent, brushedmetal, diffuse, faint, \
               clearglass, tintedglass, frostedglass, glossy, hardplastic, softmetal, steel, \
               translucent, inked, inkedgloss, inkedglass, textbook, polished, chalk, chalkedge, \
               bubbledglass, hollowglass, mirrorchrome; plus toon, clay (a soft, shadow-free \
               studio look), rim (a bright edge glow instead of an outline's darkened one), \
               studiomatte and labdefault (a matte studio look and a flat, broadly lit look). Old names \
               (glass1-3, aochalky, goodsell, ...) still parse; docs/COMMANDS.md \
               lists them. Undoable.",
    },
    Spec {
        id: "frame",
        title: "Set frame",
        keywords: &["timeline", "model", "trajectory", "coordinate", "scrub"],
        usage: "frame N [ID]",
        help: "Show coordinate set N (0-based) of a structure (default: the current structure \
               if it has more than one frame, else the first loaded structure that does). \
               Undoable; the frame a saved session persists. Playing back or scrubbing the \
               Timeline tab does not go through this (it would flood undo history).",
    },
    Spec {
        id: "color",
        title: "Set coloring",
        keywords: &["colour", "element", "chain", "bfactor", "values", "structure", "secondary", "restype", "rainbow", "hetero", "resname", "occupancy", "hydrophobicity", "red", "blue"],
        usage: "color element|chain|structure|restype|resname|rainbow|hetero|bfactor|occupancy|hydrophobicity|NAME|#RRGGBB|values CHANNEL [ID]",
        help: "Color a structure (default: current) by element, chain, secondary structure, residue type (acidic/basic/polar/nonpolar), rainbow N->C, carbons-by-chain (hetero), B-factor, or an attached value channel (see `values`).",
    },
    Spec {
        id: "values",
        title: "Per-atom values",
        keywords: &["channel", "npy", "sasa", "score", "column", "attach"],
        usage: "values [NAME PATH | remove NAME] [ID]",
        help: "Attach a column from a .npy or text file (one number per atom or per residue, optionally per frame) under NAME, list the attached channels, or remove one. `color values NAME` colors by it.",
    },
    Spec {
        id: "label",
        title: "Label atom",
        keywords: &["text", "tag", "structure-anchored", "residue"],
        usage: "label ATOM TEXT",
        help: "Attach a text label to atom index ATOM of the current structure, billboarded in the viewport and following it every frame. Overwrites an existing label on that atom.",
    },
    Spec {
        id: "unlabel",
        title: "Remove label",
        keywords: &["text", "tag"],
        usage: "unlabel ATOM",
        help: "Remove the label from atom index ATOM of the current structure, if any.",
    },
    Spec {
        id: "labels",
        title: "List labels",
        keywords: &["text", "tags"],
        usage: "labels",
        help: "List every atom label across all loaded structures.",
    },
    Spec {
        id: "caption",
        title: "Add caption",
        keywords: &["text", "screen-space", "figure", "title"],
        usage: "caption NAME X Y TEXT",
        help: "Add (or replace) a free-floating text caption named NAME at fractional viewport position X Y (0..1 from the top-left), not tied to any structure.",
    },
    Spec {
        id: "uncaption",
        title: "Remove caption",
        keywords: &["text", "screen-space", "figure"],
        usage: "uncaption NAME",
        help: "Remove the caption named NAME.",
    },
    Spec {
        id: "captions",
        title: "List captions",
        keywords: &["text", "screen-space", "figure"],
        usage: "captions",
        help: "List every screen-space caption.",
    },
    Spec {
        id: "undo",
        title: "Undo",
        keywords: &["back", "revert"],
        usage: "undo",
        help: "Undo the last command.",
    },
    Spec {
        id: "redo",
        title: "Redo",
        keywords: &["forward", "again"],
        usage: "redo",
        help: "Redo the last undone command.",
    },
    Spec {
        id: "info",
        title: "Structure info",
        keywords: &["annotations", "header", "title", "resolution", "method", "metadata"],
        usage: "info [ID]",
        help: "The file's annotations for a structure (default: current): title, method, resolution, organism, entities, citation.",
    },
    Spec {
        id: "measure",
        title: "Measure",
        keywords: &["distance", "angle", "dihedral", "torsion", "geometry"],
        usage: "measure [A B [C [D]]]",
        help: "Distance (2 atom indices), angle (3), or dihedral (4, IUPAC sign) in the current structure at its current frame, kept on screen and measured again as the frame changes. No atoms: the selection's, when it has 2 to 4.",
    },
    Spec {
        id: "unmeasure",
        title: "Remove measurement",
        keywords: &["distance", "angle", "dihedral", "clear"],
        usage: "unmeasure A B [C [D]] | all",
        help: "Take a measurement of the current structure off the screen, or all of them.",
    },
    Spec {
        id: "measurements",
        title: "List measurements",
        keywords: &["distance", "angle", "dihedral"],
        usage: "measurements",
        help: "Every measurement on screen, with its value at the current frame.",
    },
    Spec {
        id: "interactions",
        title: "Interactions",
        keywords: &["hbond", "hydrogen", "metal", "coordination", "salt", "bridge", "dashes"],
        usage: "interactions [hbond|metal|saltbridge on|off] [id]",
        help: "Draw dashed lines for hydrogen bonds, metal coordination or salt bridges of the current structure. Alone, lists which are on.",
    },
    Spec {
        id: "contacts",
        title: "Contacts",
        keywords: &["interface", "neighbors", "within", "pairs"],
        usage: "contacts CUTOFF EXPR_A | EXPR_B",
        help: "Atom pairs within CUTOFF A between two selections of the current structure: counts and the residue pairs.",
    },
    Spec {
        id: "sasa",
        title: "Solvent-accessible surface area",
        keywords: &["area", "exposure", "accessible", "fastsasa", "buried", "surface"],
        usage: "sasa [FILTER] [| SELECT]",
        help: "Per-atom solvent-accessible area of the current structure's current frame, as FastSASA computes it (Shrake-Rupley, 100 points, 1.4 A probe, ProtOr radii, element radii for the rest). FILTER is the atoms in the calculation (default `not hydrogen and not hetero`, FastSASA's default; name hydrogens or hetero atoms to include them, e.g. `sasa all`); SELECT reports a subset with the rest still occluding, like FastSASA's --select. Attached as the value channel `sasa`: `color values sasa`.",
    },
    Spec {
        id: "savestructure",
        title: "Export structure…",
        keywords: &["export", "write", "pdb", "mmcif", "cif", "xyz", "pqr", "gro", "save"],
        usage: "savestructure PATH [SELECTION] [frame N|allframes]",
        help: "Write the current structure to PATH: format from its extension (.pdb/.ent, .cif/.mmcif/.pdbx, .xyz, .pqr, .gro; any of them optionally .gz). SELECTION restricts which atoms (docs/SELECTION.md; default all); frame N picks one coordinate set (default the current frame), allframes writes every one. PATH may not contain spaces.",
    },
    Spec {
        id: "savesession",
        title: "Save session…",
        keywords: &["session", "save", "project", "workspace"],
        usage: "savesession [PATH.vviz]",
        help: "Write loaded structures (as paths), value channels (as .npy sidecars), selection sets, and the active selection to a .vviz session file (JSON). In the app, no PATH opens a save dialog.",
    },
    Spec {
        id: "loadsession",
        title: "Open session…",
        keywords: &["session", "open", "load", "restore", "project"],
        usage: "loadsession [PATH.vviz]",
        help: "Replace what is open with a saved session (.vviz, or a legacy .json); every step is undoable. In the app, no PATH opens a file picker.",
    },
    Spec {
        id: "version",
        title: "Version",
        keywords: &["about", "semver", "language"],
        usage: "version",
        help: "The command language's version (semver; docs/VERSIONING.md).",
    },
    Spec {
        id: "help",
        title: "Help",
        keywords: &["commands", "usage", "?"],
        usage: "help [COMMAND]",
        help: "List commands, or show one command's usage.",
    },
];

/// A command that could not be parsed or failed to apply. Plain text;
/// callers decide where it goes (a log line, a Python exception, stderr).
#[derive(Clone, Debug, PartialEq, Eq, thiserror::Error)]
#[error("{0}")]
pub struct ScriptError(pub String);

impl From<SceneError> for ScriptError {
    fn from(e: SceneError) -> Self {
        ScriptError(e.to_string())
    }
}

/// The command spec whose id (or alias) is `word`.
pub fn spec(word: &str) -> Option<&'static Spec> {
    let id = match word {
        "rep" => "representation",
        "colour" => "color",
        "?" => "help",
        other => other,
    };
    SPECS.iter().find(|s| s.id == id)
}

/// Splits script text into command lines: one per line or `;`, trimmed,
/// with blank lines and `#` comments dropped. Quoting is deliberately not
/// supported: paths with spaces work because a command's tail is taken
/// whole, and `;` inside a path is rare enough to live with.
pub fn split_script(text: &str) -> Vec<String> {
    text.lines()
        .flat_map(|line| line.split(';'))
        .map(|s| s.trim())
        .filter(|s| !s.is_empty() && !s.starts_with('#'))
        .map(str::to_owned)
        .collect()
}

/// `(verb, rest)` of a command line; `rest` is trimmed and may be empty.
pub fn split_verb(line: &str) -> (&str, &str) {
    let line = line.trim();
    match line.split_once(char::is_whitespace) {
        Some((verb, rest)) => (verb, rest.trim()),
        None => (line, ""),
    }
}

/// `rest` split on whitespace, with `"..."` keeping a path that
/// contains spaces as one argument.
pub fn split_args(rest: &str) -> Vec<String> {
    let mut args = Vec::new();
    let mut chars = rest.trim().chars().peekable();
    while let Some(&c) = chars.peek() {
        if c.is_whitespace() {
            chars.next();
        } else if c == '"' {
            chars.next();
            args.push(chars.by_ref().take_while(|&c| c != '"').collect());
        } else {
            let mut word = String::new();
            while let Some(&c) = chars.peek().filter(|c| !c.is_whitespace()) {
                word.push(c);
                chars.next();
            }
            args.push(word);
        }
    }
    args
}

/// `path` as one command argument: quoted when it holds a space.
pub fn quote_arg(path: &str) -> String {
    if path.contains(char::is_whitespace) {
        format!("\"{path}\"")
    } else {
        path.to_owned()
    }
}

/// The structure named by `word` (a raw id from `structures`), or the
/// current one when `word` is empty.
/// The rep that `representation`, `color` and `material` edit.
fn current_rep(scene: &Scene, id: StructureId) -> RepId {
    scene.structure(id).expect("resolved").rep().id
}

/// Rep number `word` (an index, as `reps` lists them) of structure `id`;
/// the current rep when `word` is empty.
fn rep_at(scene: &Scene, id: StructureId, word: &str) -> Result<RepId, ScriptError> {
    let loaded = scene.structure(id).expect("resolved");
    if word.is_empty() {
        return Ok(loaded.rep().id);
    }
    let n: usize = word
        .parse()
        .map_err(|_| ScriptError(format!("expected a rep number, got `{word}`")))?;
    loaded
        .reps
        .get(n)
        .map(|r| r.id)
        .ok_or_else(|| ScriptError(format!("no rep {n} (`reps` lists them)")))
}

/// The `[frame N|allframes]` tail of `savestructure`'s arguments.
enum FrameSpec {
    Current,
    One(usize),
    All,
}

/// Splits `savestructure`'s trailing `[SELECTION] [frame N|allframes]`
/// into the selection text (joined back with single spaces; empty means
/// every atom) and the frame spec, by peeling a recognized frame suffix
/// off the end of the whitespace-separated words.
fn split_frame_spec(rest: &str) -> (String, FrameSpec) {
    let words: Vec<&str> = rest.split_whitespace().collect();
    let (selection, spec) = match words.as_slice() {
        [head @ .., "allframes"] => (head, FrameSpec::All),
        [head @ .., "frame", n] if n.parse::<usize>().is_ok() => {
            (head, FrameSpec::One(n.parse().expect("checked")))
        }
        all => (all, FrameSpec::Current),
    };
    (selection.join(" "), spec)
}

/// The open structure loaded from `path`, if any.
fn open_structure_at(scene: &Scene, path: &Path) -> Option<StructureId> {
    let wanted = std::fs::canonicalize(path).unwrap_or_else(|_| path.to_path_buf());
    scene.structures().find_map(|(id, s)| {
        let open = s.path.as_ref()?;
        let open = std::fs::canonicalize(open).unwrap_or_else(|_| open.clone());
        (open == wanted).then_some(id)
    })
}

fn attach_trajectory(
    scene: &mut Scene,
    history: &mut CommandHistory,
    id: StructureId,
    trajectory: &str,
) -> Result<String, ScriptError> {
    history.dispatch(
        scene,
        Command::AttachTrajectory {
            id,
            trajectory: trajectory.into(),
        },
    )?;
    let loaded = scene.structure(id).expect("just attached");
    Ok(format!(
        "attached {} frames to {} (#{})",
        loaded.structure.frame_count(),
        loaded.label,
        id.to_raw()
    ))
}

fn resolve_structure(scene: &Scene, word: &str) -> Result<StructureId, ScriptError> {
    if word.is_empty() {
        return current(scene);
    }
    let raw: u32 = word
        .parse()
        .map_err(|_| ScriptError(format!("expected a structure id, got `{word}`")))?;
    let id = StructureId::from_raw(raw);
    if scene.structure(id).is_none() {
        return Err(ScriptError(format!("no structure with id {raw}")));
    }
    Ok(id)
}

/// The most recently loaded structure still open: what verbs act on when
/// no id is given. Mirrors the app's notion of "current".
/// Atom indices of structure `id`, whitespace-separated in `words`.
fn atom_indices(scene: &Scene, id: StructureId, words: &str) -> Result<Vec<u32>, ScriptError> {
    let n = scene
        .structure(id)
        .expect("current exists")
        .structure
        .atom_count();
    words
        .split_whitespace()
        .map(|w| {
            w.parse::<u32>()
                .ok()
                .filter(|&a| (a as usize) < n)
                .ok_or_else(|| ScriptError(format!("`{w}` is not an atom index below {n}")))
        })
        .collect()
}

/// The atoms of the active selection on structure `id`, in index order,
/// when there are 2 to 4 of them.
fn selected_atoms(scene: &Scene, id: StructureId) -> Result<Vec<u32>, ScriptError> {
    let atoms: Vec<u32> = scene
        .active_selection()
        .filter(|s| s.structure == id)
        .map(|s| s.mask.ones().map(|a| a as u32).collect())
        .unwrap_or_default();
    if (2..=4).contains(&atoms.len()) {
        Ok(atoms)
    } else {
        Err(ScriptError(format!(
            "select 2 to 4 atoms to measure (the selection has {})",
            atoms.len()
        )))
    }
}

/// `distance CA(ALA 5)-CB(ALA 5): 1.530 A`, at the current frame.
fn describe_measurement(loaded: &LoadedStructure, m: &Measurement) -> String {
    let t = &loaded.structure.topology;
    let coords = loaded.structure.frame(loaded.frame);
    let names: Vec<String> = m
        .atoms()
        .iter()
        .map(|&a| {
            let r = t.residue_index[a as usize] as usize;
            format!(
                "{}({} {})",
                t.atom_name(a as usize),
                t.residue_name(r),
                t.residues[r].auth_seq_id
            )
        })
        .collect();
    format!(
        "{} {}: {}",
        m.kind(),
        names.join("-"),
        m.text(coords.positions())
    )
}

/// The forms of `select`: a new selection, or one combined with the
/// active selection of the same structure.
#[derive(Debug, PartialEq)]
enum SelectForm<'a> {
    Replace(&'a str),
    Add(&'a str),
    Remove(&'a str),
    Invert,
}

fn select_form(rest: &str) -> SelectForm<'_> {
    if rest == "invert" {
        return SelectForm::Invert;
    }
    match rest.split_once(char::is_whitespace) {
        Some(("add", expr)) if !expr.trim().is_empty() => SelectForm::Add(expr.trim()),
        Some(("remove", expr)) if !expr.trim().is_empty() => SelectForm::Remove(expr.trim()),
        _ => SelectForm::Replace(rest),
    }
}

/// The mask `Add` (union with the active selection) or `Invert` (its
/// complement; everything when nothing is selected) leaves selected.
fn combined_selection(
    scene: &Scene,
    id: StructureId,
    form: SelectForm<'_>,
) -> Result<crate::Mask, ScriptError> {
    let loaded = scene
        .structure(id)
        .ok_or_else(|| ScriptError("no structure loaded".into()))?;
    let structure = &loaded.structure;
    let same = scene.active_selection().filter(|a| a.structure == id);
    let mut mask = same.map_or_else(
        || vv_core::fixedbitset::FixedBitSet::with_capacity(structure.atom_count()),
        |a| (*a.mask).clone(),
    );
    let hits_of = |expr: &str| {
        loaded
            .select(expr, 0)
            .map_err(|e| ScriptError(format!("{expr}: {e}")))
    };
    match form {
        SelectForm::Add(expr) => mask.union_with(&hits_of(expr)?),
        SelectForm::Remove(expr) => mask.difference_with(&hits_of(expr)?),
        SelectForm::Invert => mask.toggle_range(..),
        SelectForm::Replace(_) => unreachable!("Replace goes through SelectExpr"),
    }
    Ok(Arc::new(mask))
}

pub fn current(scene: &Scene) -> Result<StructureId, ScriptError> {
    scene
        .structures()
        .next_back()
        .map(|(id, _)| id)
        .ok_or_else(|| ScriptError("no structure loaded".into()))
}

/// The structure `frame` acts on when no id is given: the current
/// structure if it has more than one frame, else the first loaded
/// structure that does. Matches the app's Timeline tab's own notion of
/// "what am I scrubbing" (`vv-app`'s `timeline_target`) — this is the
/// document-level source of truth now, so that function defers to this
/// one rather than the two independently searching the same way.
pub fn frame_target(scene: &Scene, word: &str) -> Result<StructureId, ScriptError> {
    if !word.is_empty() {
        return resolve_structure(scene, word);
    }
    let multi = |id: StructureId| {
        scene
            .structure(id)
            .is_some_and(|s| s.structure.frame_count() > 1)
    };
    current(scene)
        .ok()
        .filter(|&id| multi(id))
        .or_else(|| {
            scene
                .structures()
                .find(|&(id, _)| multi(id))
                .map(|(id, _)| id)
        })
        .ok_or_else(|| ScriptError("no multi-frame structure loaded".into()))
}

/// Inverse of `Representation`'s own canonical names, plus the aliases the
/// console/palette/`--exec`/`Session.exec()` accept. `Session.set_representation`
/// (`vv-py`) delegates here too, so a name accepted by one is accepted by all.
pub fn parse_representation(word: &str) -> Result<Representation, ScriptError> {
    match word {
        "spacefill" | "cpk" | "spheres" => Ok(Representation::Spacefill),
        "ballstick" | "ball_and_stick" | "ball-and-stick" => Ok(Representation::BallAndStick),
        "tube" | "backbone" => Ok(Representation::Tube),
        "cartoon" | "ribbon" => Ok(Representation::Cartoon),
        "gaussiansurface" | "gaussian_surface" | "surface" | "gaussian" | "blob" | "blobby" => {
            Ok(Representation::GaussianSurface)
        }
        "skinsurface" | "skin_surface" | "skin" => Ok(Representation::SkinSurface),
        "ses" | "molecular" | "msms" | "connolly" => Ok(Representation::Ses),
        "sas" | "accessible" => Ok(Representation::Sas),
        // "sticks" used to mean ball-and-stick here; it now names its own
        // rep. `licorice` and `sticks_only` are older names that still parse.
        "sticks" | "licorice" | "sticks_only" => Ok(Representation::Sticks),
        "lines" | "wire" | "wireframe" => Ok(Representation::Lines),
        "glycan" | "snfg" | "3d-snfg" | "carbohydrate" => Ok(Representation::Glycan),
        other => Err(ScriptError(format!(
            "unknown representation `{other}`; expected spacefill, ballstick, sticks, lines, sas, tube, cartoon, gaussiansurface, skinsurface, ses, or glycan"
        ))),
    }
}

fn parse_coloring(word: &str) -> Result<ColorScheme, ScriptError> {
    ColorScheme::parse(word).ok_or_else(|| {
        ScriptError(format!(
            "unknown coloring `{word}`; expected element, chain, bfactor, or values NAME"
        ))
    })
}

fn usage(id: &str) -> ScriptError {
    let usage = spec(id).map(|s| s.usage).unwrap_or(id);
    ScriptError(format!("usage: {usage}"))
}

/// Runs one command line against the scene. `Ok` carries a short
/// human-readable result (what happened, or what was listed).
pub fn run_line(
    scene: &mut Scene,
    history: &mut CommandHistory,
    line: &str,
) -> Result<String, ScriptError> {
    let (verb, rest) = split_verb(line);
    if verb.is_empty() {
        return Ok(String::new());
    }
    let Some(spec) = spec(verb) else {
        return Err(ScriptError(format!("unknown command `{verb}`; try `help`")));
    };
    match spec.id {
        "load" => {
            if rest.is_empty() {
                return Err(usage("load"));
            }
            history.dispatch(scene, Command::LoadStructure { path: rest.into() })?;
            let (id, loaded) = scene.structures().next_back().expect("just loaded");
            Ok(format!(
                "loaded {} as #{} ({} atoms)",
                loaded.label,
                id.to_raw(),
                loaded.structure.atom_count()
            ))
        }
        "loadtrajectory" => {
            let [topology, trajectory] =
                <[String; 2]>::try_from(split_args(rest)).map_err(|_| usage("loadtrajectory"))?;
            if let Some(id) = open_structure_at(scene, Path::new(&topology)) {
                return attach_trajectory(scene, history, id, &trajectory);
            }
            history.dispatch(
                scene,
                Command::LoadTrajectory {
                    topology: topology.into(),
                    trajectory: trajectory.into(),
                },
            )?;
            let (id, loaded) = scene.structures().next_back().expect("just loaded");
            Ok(format!(
                "loaded {} as #{} ({} atoms, {} frames)",
                loaded.label,
                id.to_raw(),
                loaded.structure.atom_count(),
                loaded.structure.frame_count()
            ))
        }
        "attachtrajectory" => {
            let args = split_args(rest);
            let (trajectory, structure) = match args.as_slice() {
                [t] => (t, ""),
                [t, s] => (t, s.as_str()),
                _ => return Err(usage("attachtrajectory")),
            };
            let id = resolve_structure(scene, structure)?;
            attach_trajectory(scene, history, id, trajectory)
        }
        "fetch" => {
            let mut words = rest.split_whitespace();
            let Some(pdb_id) = words.next() else {
                return Err(usage("fetch"));
            };
            let cache = vv_io::fetch::cache_dir()
                .ok_or_else(|| ScriptError("no cache directory (no home directory?)".into()))?;
            let fetched = match (words.next(), words.next(), words.next()) {
                (None, _, _) => {
                    vv_io::fetch::fetch(pdb_id, vv_io::fetch::Assembly::AsymmetricUnit, &cache)
                }
                (Some(w), None, _) if w.eq_ignore_ascii_case("bcif") => {
                    vv_io::fetch::fetch_bcif(pdb_id, &cache)
                }
                (Some(w), Some(n), None) if w.eq_ignore_ascii_case("assembly") => {
                    let n: u32 = n.parse().map_err(|_| usage("fetch"))?;
                    vv_io::fetch::fetch(
                        pdb_id,
                        vv_io::fetch::Assembly::Biological(n.max(1)),
                        &cache,
                    )
                }
                _ => return Err(usage("fetch")),
            };
            let path = fetched.map_err(|e| ScriptError(e.to_string()))?;
            history.dispatch(scene, Command::LoadStructure { path: path.clone() })?;
            let (id, loaded) = scene.structures().next_back().expect("just loaded");
            Ok(format!(
                "fetched {} as #{} ({} atoms) -> {}",
                loaded.label,
                id.to_raw(),
                loaded.structure.atom_count(),
                path.display()
            ))
        }
        "close" => {
            let id = resolve_structure(scene, rest)?;
            history.dispatch(scene, Command::CloseStructure { id })?;
            Ok(format!("closed #{}", id.to_raw()))
        }
        "showstructure" => {
            let (state, id_word) = split_verb(rest);
            let visible = match state {
                "on" => true,
                "off" => false,
                _ => return Err(usage("showstructure")),
            };
            let id = resolve_structure(scene, id_word)?;
            history.dispatch(scene, Command::ShowStructure { id, visible })?;
            Ok(format!("#{} {state}", id.to_raw()))
        }
        "altloc" => {
            let (policy, id_word) = split_verb(rest);
            let policy: vv_core::altloc::AltlocPolicy =
                policy.parse().map_err(|_| usage("altloc"))?;
            let id = resolve_structure(scene, id_word)?;
            history.dispatch(scene, Command::SetAltloc { id, policy })?;
            Ok(format!("#{} altloc {policy}", id.to_raw()))
        }
        "structures" => {
            let mut lines: Vec<String> = scene
                .structures()
                .map(|(id, s)| {
                    format!(
                        "#{} {} ({} atoms, {} frames){}",
                        id.to_raw(),
                        s.label,
                        s.structure.atom_count(),
                        s.structure.frame_count(),
                        if s.visible { "" } else { ", hidden" }
                    )
                })
                .collect();
            if lines.is_empty() {
                lines.push("no structures loaded".into());
            }
            Ok(lines.join("\n"))
        }
        "select" => {
            if rest.is_empty() {
                return Err(usage("select"));
            }
            if let Some(done) = crate::select_ops::run(scene, history, rest)? {
                return Ok(done);
            }
            let id = current(scene)?;
            match select_form(rest) {
                SelectForm::Replace(expr) => history.dispatch(
                    scene,
                    Command::SelectExpr {
                        id,
                        expr: expr.to_owned(),
                    },
                )?,
                form => {
                    let mask = combined_selection(scene, id, form)?;
                    history.dispatch(scene, Command::Select { id, mask })?;
                }
            }
            let n = scene
                .active_selection()
                .map_or(0, |a| a.mask.count_ones(..));
            Ok(format!("selected {n} atom(s)"))
        }
        "clear" => {
            history.dispatch(scene, Command::ClearSelection)?;
            Ok("selection cleared".into())
        }
        "saveset" => {
            if rest.is_empty() {
                return Err(usage("saveset"));
            }
            history.dispatch(scene, Command::SaveSelectionSet { name: rest.into() })?;
            Ok(format!("saved selection set `{rest}`"))
        }
        "useset" => {
            if rest.is_empty() {
                return Err(usage("useset"));
            }
            let set = scene
                .selection_set(rest)
                .ok_or_else(|| ScriptError(format!("no selection set named `{rest}`")))?;
            let (id, mask) = (set.structure, Arc::clone(&set.mask));
            let n = mask.count_ones(..);
            history.dispatch(scene, Command::Select { id, mask })?;
            Ok(format!("selected {n} atom(s) from `{rest}`"))
        }
        "deleteset" => {
            if rest.is_empty() {
                return Err(usage("deleteset"));
            }
            history.dispatch(scene, Command::DeleteSelectionSet { name: rest.into() })?;
            Ok(format!("deleted selection set `{rest}`"))
        }
        "sets" => {
            let mut lines: Vec<String> = scene
                .selection_sets()
                .iter()
                .map(|s| match &s.expr {
                    Some(expr) => format!(
                        "{} (#{}, {} atoms) = {expr}",
                        s.name,
                        s.structure.to_raw(),
                        s.mask.count_ones(..)
                    ),
                    None => format!(
                        "{} (#{}, {} atoms)",
                        s.name,
                        s.structure.to_raw(),
                        s.mask.count_ones(..)
                    ),
                })
                .collect();
            if lines.is_empty() {
                lines.push("no selection sets".into());
            }
            Ok(lines.join("\n"))
        }
        "representation" => {
            let (what, id_word) = split_verb(rest);
            if what.is_empty() {
                return Err(usage("representation"));
            }
            let representation = parse_representation(what)?;
            let id = resolve_structure(scene, id_word)?;
            let rep = current_rep(scene, id);
            history.dispatch(
                scene,
                Command::SetRepresentation {
                    id,
                    rep,
                    representation,
                },
            )?;
            Ok(format!("#{} drawn as {what}", id.to_raw()))
        }
        "material" => {
            let (what, id_word) = split_verb(rest);
            let material = Material::parse(what).ok_or_else(|| usage("material"))?;
            let id = resolve_structure(scene, id_word)?;
            let rep = current_rep(scene, id);
            history.dispatch(scene, Command::SetMaterial { id, rep, material })?;
            Ok(format!("#{} made {what}", id.to_raw()))
        }
        "addrep" => {
            let (what, selection) = split_verb(rest);
            if what.is_empty() {
                return Err(usage("addrep"));
            }
            let representation = parse_representation(what)?;
            let id = current(scene)?;
            let loaded = scene.structure(id).expect("current");
            let rep = Rep {
                id: RepId(loaded.next_rep_id),
                selection: if selection.is_empty() {
                    "all".into()
                } else {
                    selection.into()
                },
                representation,
                ..loaded.rep().clone()
            };
            vv_core::select::parse(&rep.selection)
                .map_err(|e| ScriptError(format!("bad selection `{}`: {e}", rep.selection)))?;
            let index = loaded.reps.len();
            history.dispatch(scene, Command::AddRep { id, index, rep })?;
            Ok(format!(
                "#{} rep {index}: {what} of {}",
                id.to_raw(),
                scene.structure(id).expect("current").rep().selection
            ))
        }
        "delrep" => {
            let id = current(scene)?;
            let rep = rep_at(scene, id, rest)?;
            history.dispatch(scene, Command::RemoveRep { id, rep })?;
            Ok(format!("#{} rep deleted", id.to_raw()))
        }
        "selrep" => {
            if rest.is_empty() {
                return Err(usage("selrep"));
            }
            let id = current(scene)?;
            let rep = current_rep(scene, id);
            history.dispatch(
                scene,
                Command::SetRepSelection {
                    id,
                    rep,
                    selection: rest.to_string(),
                },
            )?;
            Ok(format!("#{} rep draws {rest}", id.to_raw()))
        }
        "currep" => {
            let id = current(scene)?;
            let rep = rep_at(scene, id, rest)?;
            history.dispatch(scene, Command::SetCurrentRep { id, rep })?;
            Ok(format!("#{} editing rep {rest}", id.to_raw()))
        }
        "repopt" => {
            let id = current(scene)?;
            let loaded = scene.structure(id).expect("current");
            let rep = loaded.rep().clone();
            let options = rep.representation.options();
            if rest.is_empty() {
                if options.is_empty() {
                    return Ok(format!(
                        "{} has no options",
                        crate::session::representation_name(rep.representation)
                    ));
                }
                return Ok(options
                    .iter()
                    .map(|o| {
                        let now = rep.option(o.name).unwrap_or(o.default);
                        if o.choices.is_empty() {
                            format!(
                                "{} = {now}{} ({} to {}; {})",
                                o.name, o.unit, o.min, o.max, o.label
                            )
                        } else {
                            format!(
                                "{} = {} ({}; {})",
                                o.name,
                                o.show(now),
                                o.choices.join("|"),
                                o.label
                            )
                        }
                    })
                    .collect::<Vec<_>>()
                    .join("\n"));
            }
            let (name, value) = split_verb(rest);
            let Some(spec) = options.iter().find(|o| o.name == name) else {
                return Err(ScriptError(format!(
                    "{} has no option `{name}` (try `repopt`)",
                    crate::session::representation_name(rep.representation)
                )));
            };
            let value = match value {
                "default" => None,
                v => {
                    let v = spec.parse(v).ok_or_else(|| usage("repopt"))?;
                    if !(spec.min..=spec.max).contains(&v) {
                        return Err(ScriptError(format!(
                            "{name} must be {} to {}",
                            spec.min, spec.max
                        )));
                    }
                    Some(v)
                }
            };
            history.dispatch(
                scene,
                Command::SetRepOption {
                    id,
                    rep: rep.id,
                    name: name.into(),
                    value,
                },
            )?;
            Ok(format!(
                "{name} = {}",
                spec.show(value.unwrap_or(spec.default))
            ))
        }
        "showrep" => {
            let (n, state) = split_verb(rest);
            let visible = match state {
                "on" => true,
                "off" => false,
                _ => return Err(usage("showrep")),
            };
            let id = current(scene)?;
            let rep = rep_at(scene, id, n)?;
            history.dispatch(scene, Command::ShowRep { id, rep, visible })?;
            Ok(format!("#{} rep {n} {state}", id.to_raw()))
        }
        "reps" => {
            let id = current(scene)?;
            let loaded = scene.structure(id).expect("current");
            let lines: Vec<String> = loaded
                .reps
                .iter()
                .enumerate()
                .map(|(i, r)| {
                    format!(
                        "{}{i}: {} {} {} | {}{}",
                        if i == loaded.current_rep { "*" } else { " " },
                        crate::session::representation_name(r.representation),
                        r.coloring.name(),
                        r.material.name(),
                        r.selection,
                        if r.visible { "" } else { " (hidden)" },
                    )
                })
                .collect();
            Ok(lines.join("\n"))
        }
        "frame" => {
            let (n_word, id_word) = split_verb(rest);
            if n_word.is_empty() {
                return Err(usage("frame"));
            }
            let n: usize = n_word.parse().map_err(|_| usage("frame"))?;
            let id = frame_target(scene, id_word)?;
            history.dispatch(scene, Command::SetFrame { id, frame: n })?;
            Ok(format!("#{} frame {n}", id.to_raw()))
        }
        "color" => {
            let (what, mut id_word) = split_verb(rest);
            if what.is_empty() {
                return Err(usage("color"));
            }
            let coloring = if what == "values" {
                let (name, after) = split_verb(id_word);
                if name.is_empty() {
                    return Err(usage("color"));
                }
                id_word = after;
                ColorScheme::Values(name.to_string())
            } else {
                parse_coloring(what)?
            };
            let id = resolve_structure(scene, id_word)?;
            let missing = match &coloring {
                ColorScheme::Values(name) => !scene
                    .structure(id)
                    .expect("resolved")
                    .values
                    .contains_key(name),
                _ => false,
            };
            let label = coloring.name();
            let rep = current_rep(scene, id);
            history.dispatch(scene, Command::SetColoring { id, rep, coloring })?;
            Ok(if missing {
                format!(
                    "#{} colored by {label} (no such channel yet; drawn by element until `values` attaches it)",
                    id.to_raw()
                )
            } else {
                format!("#{} colored by {label}", id.to_raw())
            })
        }
        "values" => {
            let words: Vec<&str> = rest.split_whitespace().collect();
            match words.as_slice() {
                [] | [_] if words.first().is_none_or(|w| w.parse::<u32>().is_ok()) => {
                    let id = resolve_structure(scene, words.first().copied().unwrap_or(""))?;
                    let loaded = scene.structure(id).expect("resolved");
                    if loaded.values.is_empty() {
                        return Ok(format!("#{}: no value channels", id.to_raw()));
                    }
                    let lines: Vec<String> = loaded
                        .values
                        .iter()
                        .map(|(name, c)| {
                            let (lo, hi) = c.range();
                            let frames = if c.frames() > 1 {
                                format!(", {} frames", c.frames())
                            } else {
                                String::new()
                            };
                            format!("{name}: {lo:.3} .. {hi:.3}{frames}")
                        })
                        .collect();
                    Ok(lines.join("\n"))
                }
                ["remove", name, rest @ ..] if rest.len() <= 1 => {
                    let id = resolve_structure(scene, rest.first().copied().unwrap_or(""))?;
                    history.dispatch(
                        scene,
                        Command::SetValues {
                            id,
                            name: name.to_string(),
                            channel: None,
                        },
                    )?;
                    Ok(format!("removed values `{name}` from #{}", id.to_raw()))
                }
                [name, path, rest @ ..] if rest.len() <= 1 => {
                    let id = resolve_structure(scene, rest.first().copied().unwrap_or(""))?;
                    let loaded = scene.structure(id).expect("resolved");
                    let channel = ValueChannel::read(std::path::Path::new(path), &loaded.structure)
                        .map_err(|e| ScriptError(e.to_string()))?;
                    let (lo, hi) = channel.range();
                    let frames = channel.frames();
                    history.dispatch(
                        scene,
                        Command::SetValues {
                            id,
                            name: name.to_string(),
                            channel: Some(channel),
                        },
                    )?;
                    Ok(format!(
                        "attached values `{name}` to #{} ({lo:.3} .. {hi:.3}{}); `color values {name}` to see it",
                        id.to_raw(),
                        if frames > 1 {
                            format!(", {frames} frames")
                        } else {
                            String::new()
                        }
                    ))
                }
                _ => Err(usage("values")),
            }
        }
        "label" => {
            let (atom_word, text) = split_verb(rest);
            if atom_word.is_empty() || text.is_empty() {
                return Err(usage("label"));
            }
            let atom: u32 = atom_word.parse().map_err(|_| usage("label"))?;
            let id = current(scene)?;
            history.dispatch(
                scene,
                Command::SetLabel {
                    id,
                    atom,
                    text: Some(text.to_string()),
                },
            )?;
            Ok(format!("labeled atom {atom} of #{}", id.to_raw()))
        }
        "unlabel" => {
            if rest.is_empty() {
                return Err(usage("unlabel"));
            }
            let atom: u32 = rest.parse().map_err(|_| usage("unlabel"))?;
            let id = current(scene)?;
            history.dispatch(
                scene,
                Command::SetLabel {
                    id,
                    atom,
                    text: None,
                },
            )?;
            Ok(format!(
                "removed label from atom {atom} of #{}",
                id.to_raw()
            ))
        }
        "labels" => {
            let mut lines: Vec<String> = scene
                .structures()
                .flat_map(|(id, s)| {
                    s.labels
                        .iter()
                        .map(move |(atom, text)| format!("#{} atom {atom}: {text}", id.to_raw()))
                })
                .collect();
            if lines.is_empty() {
                lines.push("no labels".into());
            }
            Ok(lines.join("\n"))
        }
        "caption" => {
            let (name, after) = split_verb(rest);
            let (x_word, after) = split_verb(after);
            let (y_word, text) = split_verb(after);
            if name.is_empty() || text.is_empty() {
                return Err(usage("caption"));
            }
            let x: f32 = x_word.parse().map_err(|_| usage("caption"))?;
            let y: f32 = y_word.parse().map_err(|_| usage("caption"))?;
            history.dispatch(
                scene,
                Command::SetCaption {
                    caption: Caption {
                        name: name.to_string(),
                        x,
                        y,
                        text: text.to_string(),
                    },
                },
            )?;
            Ok(format!("caption `{name}`"))
        }
        "uncaption" => {
            if rest.is_empty() {
                return Err(usage("uncaption"));
            }
            history.dispatch(scene, Command::DeleteCaption { name: rest.into() })?;
            Ok(format!("removed caption `{rest}`"))
        }
        "captions" => {
            let mut lines: Vec<String> = scene
                .captions()
                .iter()
                .map(|c| format!("{} ({:.2}, {:.2}): {}", c.name, c.x, c.y, c.text))
                .collect();
            if lines.is_empty() {
                lines.push("no captions".into());
            }
            Ok(lines.join("\n"))
        }
        "undo" => Ok(if history.undo(scene)? {
            "undone".into()
        } else {
            "nothing to undo".into()
        }),
        "redo" => Ok(if history.redo(scene)? {
            "redone".into()
        } else {
            "nothing to redo".into()
        }),
        "info" => {
            let id = resolve_structure(scene, rest)?;
            let loaded = scene.structure(id).expect("resolved");
            let t = &loaded.structure.topology;
            let a = &t.annotations;
            let mut lines = vec![format!(
                "#{} {}: {} atoms, {} residues, {} chains, {} frame(s)",
                id.to_raw(),
                loaded.label,
                t.atom_count(),
                t.residue_count(),
                t.chain_count(),
                loaded.structure.frame_count()
            )];
            let mut push = |key: &str, value: Option<String>| {
                if let Some(v) = value {
                    lines.push(format!("{key:<12}{v}"));
                }
            };
            push("title", a.title().map(str::to_owned));
            push("method", a.method().map(str::to_owned));
            push("resolution", a.resolution().map(|r| format!("{r:.2} A")));
            push("deposited", a.deposition_date().map(str::to_owned));
            push("organism", a.organism().map(str::to_owned));
            push("keywords", a.keywords().map(str::to_owned));
            push("citation", a.citation_title().map(str::to_owned));
            push("doi", a.doi().map(str::to_owned));
            let uniprot = a.uniprot_accessions();
            push("uniprot", (!uniprot.is_empty()).then(|| uniprot.join(", ")));
            for (eid, description) in a.entities() {
                push("entity", Some(format!("{eid}: {description}")));
            }
            Ok(lines.join("\n"))
        }
        "measure" => {
            let id = current(scene)?;
            let atoms = if rest.is_empty() {
                selected_atoms(scene, id)?
            } else {
                atom_indices(scene, id, rest)?
            };
            let measurement = Measurement::new(atoms).ok_or_else(|| usage("measure"))?;
            let line =
                describe_measurement(scene.structure(id).expect("current exists"), &measurement);
            history.dispatch(
                scene,
                Command::SetMeasurement {
                    id,
                    measurement,
                    shown: true,
                },
            )?;
            Ok(line)
        }
        "unmeasure" => {
            let id = current(scene)?;
            let doomed = if rest == "all" {
                scene
                    .structure(id)
                    .expect("current exists")
                    .measurements
                    .clone()
            } else {
                let atoms = atom_indices(scene, id, rest)?;
                vec![Measurement::new(atoms).ok_or_else(|| usage("unmeasure"))?]
            };
            let count = doomed.len();
            let edits = doomed
                .into_iter()
                .map(|measurement| Command::SetMeasurement {
                    id,
                    measurement,
                    shown: false,
                })
                .collect();
            history.dispatch(scene, Command::Batch(edits))?;
            Ok(format!(
                "removed {count} measurement(s) from #{}",
                id.to_raw()
            ))
        }
        "measurements" => {
            let mut lines: Vec<String> = scene
                .structures()
                .flat_map(|(id, s)| {
                    s.measurements
                        .iter()
                        .map(move |m| format!("#{} {}", id.to_raw(), describe_measurement(s, m)))
                })
                .collect();
            if lines.is_empty() {
                lines.push("no measurements".into());
            }
            Ok(lines.join("\n"))
        }
        "interactions" => {
            if rest.is_empty() {
                let id = current(scene)?;
                let on = &scene.structure(id).expect("current exists").interactions;
                let names: Vec<&str> = on.iter().map(|k| k.name()).collect();
                return Ok(if names.is_empty() {
                    "no interactions shown".into()
                } else {
                    names.join(", ")
                });
            }
            let (kind, tail) = split_verb(rest);
            let (state, id_word) = split_verb(tail);
            let id = resolve_structure(scene, id_word)?;
            let kind = vv_core::interactions::InteractionKind::parse(kind)
                .ok_or_else(|| usage("interactions"))?;
            let on = match state {
                "on" => true,
                "off" => false,
                _ => return Err(usage("interactions")),
            };
            history.dispatch(scene, Command::SetInteraction { id, kind, on })?;
            Ok(format!("#{} {} {state}", id.to_raw(), kind.name()))
        }
        "sasa" => {
            let (filter, report) = match rest.split_once('|') {
                Some((filter, report)) => (filter.trim(), Some(report.trim())),
                None => (rest, None),
            };
            let filter = if filter.is_empty() {
                "not hydrogen and not hetero"
            } else {
                filter
            };
            let id = current(scene)?;
            let loaded = scene.structure(id).expect("current exists");
            let structure = &loaded.structure;
            let coords = structure.frame(loaded.frame);
            let all = coords.positions();
            let pick = |expr: &str| -> Result<Vec<bool>, ScriptError> {
                let mask = loaded
                    .select(expr, loaded.frame)
                    .map_err(|e| ScriptError(format!("bad selection `{expr}`: {e}")))?;
                let mut picked = vec![false; all.len()];
                for a in mask.ones() {
                    picked[a] = true;
                }
                Ok(picked)
            };
            let universe = pick(filter)?;
            let conformer = vv_core::sasa::first_conformer(&structure.topology);
            let atoms: Vec<usize> = (0..all.len())
                .filter(|&a| universe[a] && conformer[a])
                .collect();
            if atoms.is_empty() {
                return Err(ScriptError(format!("`{filter}` selects no atoms")));
            }
            let reported = match report {
                Some(expr) => pick(expr)?,
                None => universe,
            };
            let radii = vv_core::sasa::radii(&structure.topology);
            let positions: Vec<_> = atoms.iter().map(|&a| all[a]).collect();
            let picked: Vec<f32> = atoms.iter().map(|&a| radii[a]).collect();
            let area = vv_core::sasa::shrake_rupley(
                &positions,
                &picked,
                vv_core::ses::WATER_PROBE,
                vv_core::sasa::DEFAULT_POINTS,
            );
            // Atoms not reported have no value (drawn mid-scale).
            let mut values = vec![f32::NAN; structure.atom_count()];
            let mut total = 0.0;
            let mut count = 0;
            for (&a, &v) in atoms.iter().zip(&area) {
                if reported[a] {
                    values[a] = v;
                    total += v;
                    count += 1;
                }
            }
            if count == 0 {
                return Err(ScriptError(format!(
                    "`{}` selects none of the {} atoms in `{filter}`",
                    report.unwrap_or(filter),
                    atoms.len()
                )));
            }
            let n = structure.atom_count();
            let channel = ValueChannel::new(values, &[n], n, structure.frame_count())
                .map_err(|e| ScriptError(e.to_string()))?;
            history.dispatch(
                scene,
                Command::SetValues {
                    id,
                    name: "sasa".into(),
                    channel: Some(channel),
                },
            )?;
            Ok(format!(
                "SASA of {count} atoms (of {} in context): {total:.1} A^2; attached as `sasa` \
                 (`color values sasa`)",
                atoms.len()
            ))
        }
        "contacts" => {
            let (cutoff_word, exprs) = split_verb(rest);
            let cutoff: f32 = cutoff_word
                .parse()
                .ok()
                .filter(|c: &f32| c.is_finite() && *c > 0.0)
                .ok_or_else(|| usage("contacts"))?;
            let Some((expr_a, expr_b)) = exprs.split_once('|') else {
                return Err(usage("contacts"));
            };
            let id = current(scene)?;
            let loaded = scene.structure(id).expect("current exists");
            let t = &loaded.structure.topology;
            let coords = loaded.structure.frame(0);
            let p = coords.positions();
            let group = |expr: &str| -> Result<Vec<u32>, ScriptError> {
                let expr = expr.trim();
                loaded
                    .select(expr, 0)
                    .map(|bits| bits.ones().map(|i| i as u32).collect())
                    .map_err(|e| ScriptError(format!("bad selection `{expr}`: {e}")))
            };
            let (a, b) = (group(expr_a)?, group(expr_b)?);
            let mut found = Vec::new();
            vv_core::analysis::contacts_into(p, &a, &b, cutoff, &mut found);
            let residues = vv_core::analysis::residue_pairs(t, &found);
            let mut lines = vec![format!(
                "{} atom pairs within {cutoff} A ({} vs {} atoms), {} residue pairs",
                found.len(),
                a.len(),
                b.len(),
                residues.len()
            )];
            let name = |r: u32| {
                let rec = &t.residues[r as usize];
                format!(
                    "{}{}:{}",
                    t.chain_name(rec.chain as usize),
                    t.residue_name(r as usize),
                    rec.auth_seq_id
                )
            };
            for pair in residues.iter().take(30) {
                lines.push(format!("  {} - {}", name(pair[0]), name(pair[1])));
            }
            if residues.len() > 30 {
                lines.push(format!("  ... {} more", residues.len() - 30));
            }
            Ok(lines.join(
                "
",
            ))
        }
        "savestructure" => {
            let (path, rest) = split_verb(rest);
            if path.is_empty() {
                return Err(usage("savestructure"));
            }
            let id = current(scene)?;
            let loaded = scene.structure(id).expect("resolved");
            let (selection, frame_spec) = split_frame_spec(rest);

            let mask = if selection.is_empty() {
                None
            } else {
                let m = loaded
                    .select(&selection, loaded.frame)
                    .map_err(|e| ScriptError(format!("bad selection `{selection}`: {e}")))?;
                if m.count_ones(..) == 0 {
                    return Err(ScriptError(format!("`{selection}` selects no atoms")));
                }
                Some(m)
            };
            let frames = match frame_spec {
                FrameSpec::Current => vec![loaded.frame],
                FrameSpec::All => (0..loaded.structure.frame_count()).collect(),
                FrameSpec::One(n) => {
                    if n >= loaded.structure.frame_count() {
                        return Err(ScriptError(format!(
                            "frame {n} out of range: structure has {} frame(s)",
                            loaded.structure.frame_count()
                        )));
                    }
                    vec![n]
                }
            };
            let opts = vv_io::SaveOptions {
                atoms: mask.as_ref(),
                frames: &frames,
            };
            let warnings = vv_io::save(&loaded.structure, path, &opts)
                .map_err(|e| ScriptError(e.to_string()))?;
            let mut lines = vec![format!(
                "saved {} atom(s), {} frame(s) to {path}",
                mask.as_ref()
                    .map_or(loaded.structure.atom_count(), |m| m.count_ones(..)),
                frames.len(),
            )];
            lines.extend(warnings.iter().map(|w| format!("warning: {w}")));
            Ok(lines.join("\n"))
        }
        "savesession" => {
            if rest.is_empty() {
                return Err(usage("savesession"));
            }
            crate::session::save(scene, std::path::Path::new(rest), None)
                .map_err(|e| ScriptError(e.to_string()))?;
            Ok(format!("saved session to {rest}"))
        }
        "loadsession" => {
            if rest.is_empty() {
                return Err(usage("loadsession"));
            }
            let file = crate::session::read(std::path::Path::new(rest))
                .map_err(|e| ScriptError(e.to_string()))?;
            let applied = crate::session::apply(&file, scene, history);
            let mut lines = vec![format!(
                "loaded session {rest}: {} structure(s), {} set(s)",
                applied.ids.iter().flatten().count(),
                scene.selection_sets().len()
            )];
            lines.extend(applied.warnings.iter().map(|w| format!("warning: {w}")));
            Ok(lines.join(
                "
",
            ))
        }
        "help" => Ok(help(rest)),
        "version" => Ok(format!("command language {}", crate::LANGUAGE_VERSION)),
        other => unreachable!("every spec id is handled: {other}"),
    }
}

/// `help` output: every command's usage, or one command's usage and help.
pub fn help(word: &str) -> String {
    if word.is_empty() {
        return SPECS
            .iter()
            .map(|s| format!("{:<44} {}", s.usage, s.title))
            .collect::<Vec<_>>()
            .join("\n");
    }
    match spec(word) {
        Some(s) => format!("{}\n  {}", s.usage, s.help),
        None => format!("unknown command `{word}`"),
    }
}

/// Runs a whole script (see [`split_script`]), stopping at the first
/// error. Returns each successful line's result, then the error if any.
pub fn run_script(
    scene: &mut Scene,
    history: &mut CommandHistory,
    text: &str,
) -> (Vec<String>, Option<(String, ScriptError)>) {
    let mut outputs = Vec::new();
    for line in split_script(text) {
        match run_line(scene, history, &line) {
            Ok(out) => outputs.push(out),
            Err(e) => return (outputs, Some((line, e))),
        }
    }
    (outputs, None)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture(name: &str) -> String {
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../fixtures/small")
            .join(name)
            .display()
            .to_string()
    }

    #[test]
    fn splits_on_semicolons_newlines_and_drops_comments() {
        let lines = split_script("load a.cif; select chain A\n# comment\n\n  color chain ;");
        assert_eq!(lines, ["load a.cif", "select chain A", "color chain"]);
    }

    #[test]
    fn every_spec_id_parses_and_aliases_resolve() {
        for s in SPECS {
            assert!(spec(s.id).is_some(), "{}", s.id);
        }
        assert_eq!(spec("rep").unwrap().id, "representation");
        assert_eq!(spec("colour").unwrap().id, "color");
        assert!(spec("bogus").is_none());
    }

    #[test]
    fn quoted_arguments_keep_their_spaces() {
        assert_eq!(
            split_args(r#""C:/Coarse grain/a.pdb" b.xtc  "x y" z"#),
            ["C:/Coarse grain/a.pdb", "b.xtc", "x y", "z"]
        );
        assert_eq!(
            quote_arg("C:/Coarse grain/a.pdb"),
            r#""C:/Coarse grain/a.pdb""#
        );
        assert_eq!(quote_arg("a.pdb"), "a.pdb");
    }

    #[test]
    fn sticks_is_the_new_name_and_licorice_still_parses() {
        assert_eq!(parse_representation("sticks"), Ok(Representation::Sticks));
        assert_eq!(parse_representation("licorice"), Ok(Representation::Sticks));
        assert_eq!(
            parse_representation("sticks_only"),
            Ok(Representation::Sticks)
        );
        assert_eq!(
            parse_representation("ballstick"),
            Ok(Representation::BallAndStick),
            "sticks no longer means ball-and-stick"
        );
    }

    /// Options are checked against the style's list and range, undo, and
    /// survive a session.
    #[test]
    fn rep_options_set_check_undo_and_save() {
        let mut scene = Scene::new();
        let mut history = CommandHistory::default();
        let mut run = |line: &str| run_line(&mut scene, &mut history, line).map_err(|e| e.0);
        run(&format!("load {}", fixture("1CRN.pdb"))).unwrap();
        run("rep ses").unwrap();
        assert!(run("repopt").unwrap().contains("probe = 1.4"));
        assert_eq!(run("repopt probe 2").unwrap(), "probe = 2");
        assert!(run("repopt probe 9").unwrap_err().contains("0.5 to 3"));
        assert!(run("repopt scale 1")
            .unwrap_err()
            .contains("no option `scale`"));
        let id = scene.structures().next().unwrap().0;
        assert_eq!(
            scene.structure(id).unwrap().rep().option("probe"),
            Some(2.0)
        );

        let path = std::env::temp_dir().join("vizviz_rep_options.json");
        crate::session::save(&scene, &path, None).unwrap();
        let file = crate::session::read(&path).unwrap();
        let mut restored = Scene::new();
        let applied = crate::session::apply(&file, &mut restored, &mut CommandHistory::default());
        let loaded = restored.structure(applied.ids[0].unwrap()).unwrap();
        assert_eq!(loaded.rep().option("probe"), Some(2.0));

        history.undo(&mut scene).unwrap();
        assert_eq!(
            scene.structure(id).unwrap().rep().option("probe"),
            Some(1.4)
        );
    }

    /// `rep glycan` (and its aliases), its `size`/`radius` options, and a
    /// session round trip -- the wiring this crate owns for the 3D-SNFG
    /// glycan representation (`vv_core::glycan` does the detection/geometry).
    #[test]
    fn rep_glycan_aliases_options_and_session_round_trip() {
        let mut scene = Scene::new();
        let mut history = CommandHistory::default();
        run_line(
            &mut scene,
            &mut history,
            &format!("load {}", fixture("1CRN.pdb")),
        )
        .unwrap();

        for word in ["glycan", "snfg", "3d-snfg", "carbohydrate"] {
            run_line(&mut scene, &mut history, &format!("rep {word}")).unwrap();
            assert_eq!(
                scene
                    .structure(current(&scene).unwrap())
                    .unwrap()
                    .rep()
                    .representation,
                Representation::Glycan,
                "`rep {word}` should select the glycan representation"
            );
        }

        let mut run = |line: &str| run_line(&mut scene, &mut history, line).map_err(|e| e.0);
        assert!(run("repopt").unwrap().contains("size = 4"));
        assert!(run("repopt").unwrap().contains("radius = 0.5"));
        assert_eq!(run("repopt size 1.5").unwrap(), "size = 1.5");
        assert_eq!(run("repopt radius 0").unwrap(), "radius = 0");
        assert!(run("repopt size 20").unwrap_err().contains("1.5 to 8"));

        let path = std::env::temp_dir().join("vizviz_rep_glycan.json");
        crate::session::save(&scene, &path, None).unwrap();
        let file = crate::session::read(&path).unwrap();
        let mut restored = Scene::new();
        let applied = crate::session::apply(&file, &mut restored, &mut CommandHistory::default());
        let loaded = restored.structure(applied.ids[0].unwrap()).unwrap();
        assert_eq!(loaded.rep().representation, Representation::Glycan);
        assert_eq!(loaded.rep().option("size"), Some(1.5));
        assert_eq!(loaded.rep().option("radius"), Some(0.0));
    }

    #[test]
    fn interactions_toggle_undo_and_survive_a_session() {
        let mut scene = Scene::new();
        let mut history = CommandHistory::default();
        run_line(
            &mut scene,
            &mut history,
            &format!("load {}", fixture("1CRN.pdb")),
        )
        .unwrap();
        let run = |scene: &mut Scene, history: &mut CommandHistory, line: &str| {
            run_line(scene, history, line).map_err(|e| e.0)
        };
        assert_eq!(
            run(&mut scene, &mut history, "interactions").unwrap(),
            "no interactions shown"
        );
        run(&mut scene, &mut history, "interactions hbond on").unwrap();
        run(&mut scene, &mut history, "interactions metal on").unwrap();
        assert_eq!(
            run(&mut scene, &mut history, "interactions").unwrap(),
            "hbond, metal"
        );
        assert!(run(&mut scene, &mut history, "interactions ionic on").is_err());
        history.undo(&mut scene).unwrap();
        assert_eq!(
            run(&mut scene, &mut history, "interactions").unwrap(),
            "hbond"
        );

        let path = std::env::temp_dir().join("vizviz_interactions.json");
        crate::session::save(&scene, &path, None).unwrap();
        let file = crate::session::read(&path).unwrap();
        let mut restored = Scene::new();
        crate::session::apply(&file, &mut restored, &mut CommandHistory::default());
        let loaded = restored.structures().next().unwrap().1;
        assert_eq!(loaded.interactions.len(), 1);
    }

    #[test]
    fn a_cartoons_choice_options_take_names_and_survive_a_session() {
        let mut scene = Scene::new();
        let mut history = CommandHistory::default();
        run_line(
            &mut scene,
            &mut history,
            &format!("load {}", fixture("1CRN.pdb")),
        )
        .unwrap();
        let mut run = |line: &str| run_line(&mut scene, &mut history, line).map_err(|e| e.0);
        run("rep cartoon").unwrap();
        assert!(run("repopt")
            .unwrap()
            .contains("bases = plate (stick|plate|ladder;"));
        assert_eq!(run("repopt bases ladder").unwrap(), "bases = ladder");
        run("rep tube").unwrap();
        assert_eq!(run("repopt bases ladder").unwrap(), "bases = ladder");
        run("rep cartoon").unwrap();
        assert_eq!(run("repopt ions off").unwrap(), "ions = off");
        assert_eq!(run("repopt water 1").unwrap(), "water = on");
        assert!(run("repopt bases sideways").is_err());
        assert!(run("repopt ions 2").unwrap_err().contains("0 to 1"));

        let path = std::env::temp_dir().join("vizviz_rep_choices.json");
        crate::session::save(&scene, &path, None).unwrap();
        let file = crate::session::read(&path).unwrap();
        let mut restored = Scene::new();
        let applied = crate::session::apply(&file, &mut restored, &mut CommandHistory::default());
        let rep = restored
            .structure(applied.ids[0].unwrap())
            .unwrap()
            .rep()
            .clone();
        assert_eq!(rep.option("bases"), Some(2.0));
        assert_eq!(rep.option("ions"), Some(0.0));
        assert_eq!(rep.option("water"), Some(1.0));
    }

    #[test]
    fn reps_add_edit_list_delete_and_undo() {
        let mut scene = Scene::new();
        let mut history = CommandHistory::default();
        let script = format!(
            "load {}\nrep cartoon\naddrep ballstick hetero and not water\nmaterial glossy",
            fixture("4HHB.cif")
        );
        let (out, err) = run_script(&mut scene, &mut history, &script);
        assert_eq!(err, None, "{out:?}");
        let id = current(&scene).unwrap();
        let loaded = scene.structure(id).unwrap();
        assert_eq!(loaded.reps.len(), 2);
        assert_eq!(loaded.reps[0].representation, Representation::Cartoon);
        assert_eq!(loaded.reps[0].material, Material::Opaque);
        assert_eq!(loaded.rep().representation, Representation::BallAndStick);
        assert_eq!(loaded.rep().selection, "hetero and not water");
        assert_eq!(loaded.rep().material, Material::Glossy);
        let second = loaded.rep().id;

        let listing = run_line(&mut scene, &mut history, "reps").unwrap();
        assert!(listing.contains("*1: ball_and_stick"), "{listing}");
        assert!(run_line(&mut scene, &mut history, "selrep within 5 of").is_err());
        assert!(run_line(&mut scene, &mut history, "currep 7").is_err());

        run_line(&mut scene, &mut history, "delrep").unwrap();
        assert_eq!(scene.structure(id).unwrap().reps.len(), 1);
        assert!(
            run_line(&mut scene, &mut history, "delrep").is_err(),
            "last rep stays"
        );
        run_line(&mut scene, &mut history, "undo").unwrap();
        let loaded = scene.structure(id).unwrap();
        assert_eq!(loaded.reps.len(), 2);
        assert_eq!(loaded.reps[1].id, second, "undo brings back the same rep");
        run_line(&mut scene, &mut history, "showrep 0 off").unwrap();
        assert!(!scene.structure(id).unwrap().reps[0].visible);
    }

    #[test]
    fn showstructure_hides_the_whole_structure_independent_of_reps_and_undoes() {
        let mut scene = Scene::new();
        let mut history = CommandHistory::default();
        run_line(
            &mut scene,
            &mut history,
            &format!("load {}", fixture("4HHB.cif")),
        )
        .unwrap();
        let id = current(&scene).unwrap();
        assert!(scene.structure(id).unwrap().visible);

        run_line(&mut scene, &mut history, "showstructure off").unwrap();
        assert!(!scene.structure(id).unwrap().visible);
        assert!(scene.structure(id).unwrap().rep().visible, "reps untouched");
        let listing = run_line(&mut scene, &mut history, "structures").unwrap();
        assert!(listing.contains("hidden"), "{listing}");

        run_line(&mut scene, &mut history, "undo").unwrap();
        assert!(scene.structure(id).unwrap().visible);
        assert!(run_line(&mut scene, &mut history, "showstructure").is_err());
    }

    #[test]
    fn altloc_sets_the_display_policy_rejects_nonsense_and_undoes() {
        use vv_core::altloc::AltlocPolicy;
        let mut scene = Scene::new();
        let mut history = CommandHistory::default();
        run_line(
            &mut scene,
            &mut history,
            &format!("load {}", fixture("1AKE.pdb")),
        )
        .unwrap();
        let id = current(&scene).unwrap();
        assert_eq!(scene.structure(id).unwrap().altloc, AltlocPolicy::First);

        run_line(&mut scene, &mut history, "altloc B").unwrap();
        assert_eq!(
            scene.structure(id).unwrap().altloc,
            AltlocPolicy::Label(b'B')
        );
        run_line(&mut scene, &mut history, "altloc all").unwrap();
        assert_eq!(scene.structure(id).unwrap().altloc, AltlocPolicy::All);
        assert!(run_line(&mut scene, &mut history, "altloc AB").is_err());
        assert!(run_line(&mut scene, &mut history, "altloc").is_err());

        run_line(&mut scene, &mut history, "undo").unwrap();
        assert_eq!(
            scene.structure(id).unwrap().altloc,
            AltlocPolicy::Label(b'B')
        );
    }

    #[test]
    fn select_forms_are_told_apart() {
        assert_eq!(select_form("protein"), SelectForm::Replace("protein"));
        assert_eq!(select_form("invert"), SelectForm::Invert);
        assert_eq!(select_form("add  water "), SelectForm::Add("water"));
        assert_eq!(select_form("add"), SelectForm::Replace("add"));
    }

    #[test]
    fn select_add_and_invert_combine_with_the_active_selection() {
        let mut scene = Scene::new();
        let mut history = CommandHistory::default();
        let load = format!("load {}", fixture("1CRN.pdb"));
        run_line(&mut scene, &mut history, &load).unwrap();
        let total = scene.structures().next().unwrap().1.structure.atom_count();
        let selected = |scene: &Scene| {
            scene
                .active_selection()
                .map_or(0, |a| a.mask.count_ones(..))
        };

        run_line(&mut scene, &mut history, "select index 0-9").unwrap();
        run_line(&mut scene, &mut history, "select add index 5-19").unwrap();
        assert_eq!(selected(&scene), 20);
        run_line(&mut scene, &mut history, "select invert").unwrap();
        assert_eq!(selected(&scene), total - 20);
        history.undo(&mut scene).unwrap();
        assert_eq!(selected(&scene), 20);

        run_line(&mut scene, &mut history, "clear").unwrap();
        run_line(&mut scene, &mut history, "select add index 0-2").unwrap();
        assert_eq!(
            selected(&scene),
            3,
            "add with nothing selected just selects"
        );
        assert!(run_line(&mut scene, &mut history, "select add nonsense_word").is_err());
    }

    #[test]
    fn a_full_session_from_text() {
        let mut scene = Scene::new();
        let mut history = CommandHistory::default();
        let script = format!(
            "load {}\nselect chain A and name CA\nsaveset ca\nclear\nuseset ca\ncolor chain\nmaterial aochalky\nrep ballstick",
            fixture("4HHB.cif")
        );
        let (out, err) = run_script(&mut scene, &mut history, &script);
        assert_eq!(err, None, "{out:?}");
        assert_eq!(out.len(), 8);
        assert!(out[0].contains("4HHB") && out[0].contains("4779 atoms"));
        assert_eq!(out[1], "selected 141 atom(s)");
        assert_eq!(out[4], "selected 141 atom(s) from `ca`");

        let id = current(&scene).unwrap();
        let loaded = scene.structure(id).unwrap();
        assert_eq!(loaded.rep().coloring, ColorScheme::Chain);
        assert_eq!(loaded.rep().representation, Representation::BallAndStick);
        assert_eq!(
            loaded.rep().material,
            Material::Chalk,
            "aochalky is now an alias for chalk"
        );
        assert!(run_line(&mut scene, &mut history, "material velvet").is_err());
        assert_eq!(
            scene.selection_set("ca").unwrap().expr.as_deref(),
            Some("chain A and name CA")
        );

        assert!(run_line(&mut scene, &mut history, "structures")
            .unwrap()
            .contains("#0 4HHB"));
        assert!(run_line(&mut scene, &mut history, "sets")
            .unwrap()
            .contains("ca (#0, 141 atoms) = chain A and name CA"));
        assert_eq!(
            run_line(&mut scene, &mut history, "undo").unwrap(),
            "undone"
        );
        assert_eq!(
            scene.structure(id).unwrap().rep().representation,
            Representation::Lines
        );
        assert_eq!(
            run_line(&mut scene, &mut history, "redo").unwrap(),
            "redone"
        );
        assert_eq!(
            run_line(&mut scene, &mut history, "close").unwrap(),
            "closed #0"
        );
        assert_eq!(scene.structures().count(), 0);
        assert_eq!(
            run_line(&mut scene, &mut history, "redo").unwrap(),
            "nothing to redo"
        );
    }

    #[test]
    fn errors_name_the_problem_and_stop_the_script() {
        let mut scene = Scene::new();
        let mut history = CommandHistory::default();
        let err = |line: &str, scene: &mut Scene, history: &mut CommandHistory| {
            run_line(scene, history, line).unwrap_err().0
        };
        assert!(err("bogus", &mut scene, &mut history).contains("unknown command `bogus`"));
        assert_eq!(
            err("select chain A", &mut scene, &mut history),
            "no structure loaded"
        );
        assert_eq!(err("load", &mut scene, &mut history), "usage: load PATH");
        assert!(err("load /no/such/file.cif", &mut scene, &mut history).contains("failed to load"));
        run_line(
            &mut scene,
            &mut history,
            &format!("load {}", fixture("1CRN.cif")),
        )
        .unwrap();
        assert!(err("select chain A and", &mut scene, &mut history).contains("bad selection"));
        assert!(err("color plaid", &mut scene, &mut history).contains("unknown coloring"));
        assert!(err("rep sticks 7", &mut scene, &mut history).contains("no structure with id 7"));
        assert!(err("close x", &mut scene, &mut history).contains("expected a structure id"));
        assert!(err("useset nope", &mut scene, &mut history).contains("no selection set"));

        let (out, failure) = run_script(&mut scene, &mut history, "clear; bogus; clear");
        assert_eq!(out, ["selection cleared"]);
        let (line, e) = failure.unwrap();
        assert_eq!(line, "bogus");
        assert!(e.0.contains("unknown command"));
    }

    #[test]
    fn values_attach_color_list_and_remove() {
        let mut scene = Scene::new();
        let mut history = CommandHistory::default();
        run_line(
            &mut scene,
            &mut history,
            &format!("load {}", fixture("1CRN.cif")),
        )
        .unwrap();
        let atoms = scene.structures().next().unwrap().1.structure.atom_count();
        let dir = std::env::temp_dir().join(format!("vizviz_script_values_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("score.txt");
        let text: String = (0..atoms)
            .map(|i| format!("{}\n", i as f32 * 0.5))
            .collect();
        std::fs::write(&path, text).unwrap();

        assert!(run_line(&mut scene, &mut history, "values")
            .unwrap()
            .contains("no value"));
        let out = run_line(
            &mut scene,
            &mut history,
            &format!("values score {}", path.display()),
        )
        .unwrap();
        assert!(out.contains("attached values `score`"), "{out}");
        assert!(run_line(&mut scene, &mut history, "values")
            .unwrap()
            .starts_with("score: 0.000"));

        let out = run_line(&mut scene, &mut history, "color values score").unwrap();
        assert!(!out.contains("no such channel"), "{out}");
        let loaded = scene.structures().next().unwrap().1;
        assert_eq!(loaded.rep().coloring, ColorScheme::Values("score".into()));
        assert_eq!(
            loaded
                .values_for(&loaded.rep().coloring)
                .unwrap()
                .1
                .frame(0)[2],
            1.0
        );

        // Wrong length is a clear error and attaches nothing.
        std::fs::write(&path, "1\n2\n").unwrap();
        let err = run_line(
            &mut scene,
            &mut history,
            &format!("values bad {}", path.display()),
        )
        .unwrap_err();
        assert!(err.0.contains("expected one value per atom"), "{err}");

        run_line(&mut scene, &mut history, "values remove score").unwrap();
        assert!(scene.structures().next().unwrap().1.values.is_empty());
        history.undo(&mut scene).unwrap();
        assert!(scene
            .structures()
            .next()
            .unwrap()
            .1
            .values
            .contains_key("score"));
        assert!(run_line(&mut scene, &mut history, "values remove nope").is_err());
        assert!(run_line(&mut scene, &mut history, "color values").is_err());
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn fetch_rejects_bad_arguments_before_touching_the_network() {
        let mut scene = Scene::new();
        let mut history = CommandHistory::default();
        assert!(run_line(&mut scene, &mut history, "fetch").is_err());
        let err = run_line(&mut scene, &mut history, "fetch notanid").unwrap_err();
        assert!(err.0.contains("not a PDB ID"), "{err}");
        assert!(run_line(&mut scene, &mut history, "fetch 4HHB assembly x").is_err());
        assert!(run_line(&mut scene, &mut history, "fetch 4HHB unit 1").is_err());
        assert_eq!(scene.structures().count(), 0);
    }

    #[test]
    fn info_prints_the_file_annotations() {
        let mut scene = Scene::new();
        let mut history = CommandHistory::default();
        run_line(
            &mut scene,
            &mut history,
            &format!("load {}", fixture("4HHB.cif")),
        )
        .unwrap();
        let out = run_line(&mut scene, &mut history, "info").unwrap();
        assert!(out.contains("4779 atoms"), "{out}");
        assert!(out.to_lowercase().contains("haemoglobin"), "{out}");
        assert!(
            out.contains("X-RAY DIFFRACTION") && out.contains("1.74 A"),
            "{out}"
        );
        assert!(out.contains("P69905"), "{out}");
        assert!(run_line(&mut scene, &mut history, "info 7").is_err());
    }

    #[test]
    fn attachtrajectory_and_loadtrajectory_of_an_open_topology_attach() {
        let mut scene = Scene::new();
        let mut history = CommandHistory::default();
        let xtc = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../fixtures/traj/1crn.xtc");
        let pdb = fixture("1CRN.pdb");
        let mut run =
            |scene: &mut Scene, line: String| run_line(scene, &mut history, &line).map_err(|e| e.0);
        run(&mut scene, format!("load {pdb}")).unwrap();
        let out = run(&mut scene, format!("attachtrajectory {}", xtc.display())).unwrap();
        assert!(out.starts_with("attached 10 frames to 1CRN"), "{out}");
        run(
            &mut scene,
            format!("loadtrajectory {pdb} {}", xtc.display()),
        )
        .unwrap();
        assert_eq!(scene.structures().count(), 1);
        let frames = scene.structures().next().unwrap().1.structure.frame_count();
        assert_eq!(frames, 10);

        let mut other = Scene::new();
        run(&mut other, format!("load {}", fixture("1AKE.pdb"))).unwrap();
        assert!(run(&mut other, format!("attachtrajectory {}", xtc.display())).is_err());
    }

    #[test]
    fn measurements_stay_on_screen_and_follow_the_frame() {
        let mut scene = Scene::new();
        let mut history = CommandHistory::default();
        // A random walk per atom (fixtures/traj/make.py): distances change.
        let xtc = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../fixtures/traj/1crn.xtc");
        let line = format!("loadtrajectory {} {}", fixture("1CRN.pdb"), xtc.display());
        run_line(&mut scene, &mut history, &line).unwrap();
        run_line(&mut scene, &mut history, "measure 0 40").unwrap();
        let listed = run_line(&mut scene, &mut history, "measurements").unwrap();
        assert!(listed.starts_with("#0 distance N(THR 1)-"), "{listed}");
        let loaded = scene.structures().next().unwrap().1;
        let m = &loaded.measurements[0];
        let at = |f: usize| m.value(loaded.structure.frame(f).positions());
        assert_ne!(at(0), at(2), "measured again at each frame");

        // The selection's atoms, in index order.
        run_line(&mut scene, &mut history, "select index 40 0 5").unwrap();
        let angle = run_line(&mut scene, &mut history, "measure").unwrap();
        assert!(angle.starts_with("angle N(THR 1)-"), "{angle}");
        let loaded = scene.structures().next().unwrap().1;
        assert_eq!(loaded.measurements.len(), 2);
        assert_eq!(loaded.measurements[1].atoms(), [0, 5, 40]);

        run_line(&mut scene, &mut history, "unmeasure 0 40").unwrap();
        run_line(&mut scene, &mut history, "unmeasure all").unwrap();
        let none = run_line(&mut scene, &mut history, "measurements").unwrap();
        assert_eq!(none, "no measurements");
        history.undo(&mut scene).unwrap();
        let back = scene.structures().next().unwrap().1;
        assert_eq!(back.measurements.len(), 1);
        run_line(&mut scene, &mut history, "clear").unwrap();
        assert!(run_line(&mut scene, &mut history, "measure").is_err());
    }

    #[test]
    fn measure_and_contacts_report_geometry() {
        let mut scene = Scene::new();
        let mut history = CommandHistory::default();
        run_line(
            &mut scene,
            &mut history,
            &format!("load {}", fixture("4HHB.cif")),
        )
        .unwrap();
        let d = run_line(&mut scene, &mut history, "measure 0 1").unwrap();
        // 4HHB's N-terminal N-CA is a badly modelled 1.16 A in the deposited
        // coordinates, so check the label and a plausible range, not 1.46.
        assert!(d.starts_with("distance N(VAL 1)-CA(VAL 1): 1."), "{d}");
        let a = run_line(&mut scene, &mut history, "measure 0 1 2").unwrap();
        assert!(a.starts_with("angle ") && a.contains(" deg"), "{a}");
        let t = run_line(&mut scene, &mut history, "measure 0 1 2 3").unwrap();
        assert!(t.starts_with("dihedral "), "{t}");
        assert!(run_line(&mut scene, &mut history, "measure 0").is_err());
        assert!(run_line(&mut scene, &mut history, "measure 0 99999").is_err());

        let c = run_line(
            &mut scene,
            &mut history,
            "contacts 4 chain A and protein | chain B and protein",
        )
        .unwrap();
        assert!(c.starts_with("102 atom pairs within 4 A"), "{c}");
        assert!(c.contains("residue pairs"));
        assert!(run_line(&mut scene, &mut history, "contacts 4 chain A").is_err());
        assert!(run_line(&mut scene, &mut history, "contacts x chain A | chain B").is_err());
    }

    #[test]
    fn savestructure_writes_a_selection_and_reports_warnings_and_errors() {
        let mut scene = Scene::new();
        let mut history = CommandHistory::default();
        run_line(
            &mut scene,
            &mut history,
            &format!("load {}", fixture("1CRN.pdb")),
        )
        .unwrap();

        let dir = std::env::temp_dir().join("vizviz-scene-tests");
        std::fs::create_dir_all(&dir).unwrap();

        let all_out = dir.join("savestructure_all.cif");
        let out = run_line(
            &mut scene,
            &mut history,
            &format!("savestructure {}", all_out.display()),
        )
        .unwrap();
        assert!(out.starts_with("saved 327 atom(s), 1 frame(s)"), "{out}");
        assert_eq!(vv_io::load(&all_out).unwrap().atom_count(), 327);

        let ca_out = dir.join("savestructure_ca.pdb");
        let out = run_line(
            &mut scene,
            &mut history,
            &format!("savestructure {} name CA", ca_out.display()),
        )
        .unwrap();
        assert!(out.starts_with("saved 46 atom(s)"), "{out}");
        assert_eq!(vv_io::load(&ca_out).unwrap().atom_count(), 46);

        let frame_out = dir.join("savestructure_frame.pdb");
        let out = run_line(
            &mut scene,
            &mut history,
            &format!("savestructure {} frame 0", frame_out.display()),
        )
        .unwrap();
        assert!(out.contains("1 frame(s)"), "{out}");

        assert!(run_line(&mut scene, &mut history, "savestructure").is_err());
        assert!(run_line(
            &mut scene,
            &mut history,
            &format!("savestructure {} resname XYZ", dir.join("x.pdb").display())
        )
        .unwrap_err()
        .0
        .contains("selects no atoms"));
        assert!(run_line(
            &mut scene,
            &mut history,
            &format!("savestructure {} frame 9", dir.join("x.pdb").display())
        )
        .unwrap_err()
        .0
        .contains("out of range"));
        assert!(run_line(
            &mut scene,
            &mut history,
            &format!("savestructure {}", dir.join("x.unknownext").display())
        )
        .is_err());
    }

    #[test]
    fn labels_and_captions_from_text() {
        let mut scene = Scene::new();
        let mut history = CommandHistory::default();
        run_line(
            &mut scene,
            &mut history,
            &format!("load {}", fixture("1CRN.cif")),
        )
        .unwrap();

        assert!(run_line(&mut scene, &mut history, "labels")
            .unwrap()
            .contains("no labels"));
        let out = run_line(&mut scene, &mut history, "label 0 Thr1 N").unwrap();
        assert!(out.contains("labeled atom 0"), "{out}");
        let loaded = scene.structures().next().unwrap().1;
        assert_eq!(loaded.labels.get(&0).map(String::as_str), Some("Thr1 N"));
        assert!(run_line(&mut scene, &mut history, "labels")
            .unwrap()
            .contains("atom 0: Thr1 N"));

        run_line(&mut scene, &mut history, "unlabel 0").unwrap();
        assert!(run_line(&mut scene, &mut history, "labels")
            .unwrap()
            .contains("no labels"));

        assert!(run_line(&mut scene, &mut history, "label 99999 x").is_err());
        assert!(run_line(&mut scene, &mut history, "label").is_err());
        assert!(run_line(&mut scene, &mut history, "label 0").is_err());

        assert!(run_line(&mut scene, &mut history, "captions")
            .unwrap()
            .contains("no captions"));
        let out = run_line(&mut scene, &mut history, "caption title 0.5 0.05 Crambin").unwrap();
        assert!(out.contains("caption `title`"), "{out}");
        assert_eq!(scene.captions().len(), 1);
        assert!(run_line(&mut scene, &mut history, "captions")
            .unwrap()
            .contains("title (0.50, 0.05): Crambin"));

        run_line(&mut scene, &mut history, "uncaption title").unwrap();
        assert!(scene.captions().is_empty());
        assert!(run_line(&mut scene, &mut history, "uncaption title").is_err());
        assert!(run_line(&mut scene, &mut history, "caption title x 0.05 text").is_err());
    }

    #[test]
    fn frame_dispatches_set_frame_defaults_to_the_multiframe_structure_and_validates() {
        let mut scene = Scene::new();
        let mut history = CommandHistory::default();
        // A single-frame structure loaded first, then the trajectory:
        // bare `frame N` (no id) must still find the multi-frame one,
        // not "current" (the single-frame structure, loaded last).
        run_line(
            &mut scene,
            &mut history,
            &format!("load {}", fixture("1CRN.cif")),
        )
        .unwrap();
        run_line(
            &mut scene,
            &mut history,
            &format!(
                "loadtrajectory {} {}",
                fixture("1CRN.pdb"),
                fixture("1CRN_traj.dcd")
            ),
        )
        .unwrap();
        let single_id = scene.structures().next().unwrap().0;
        let traj_id = scene.structures().next_back().unwrap().0;
        assert_ne!(single_id, traj_id);

        let out = run_line(&mut scene, &mut history, "frame 2").unwrap();
        assert_eq!(out, format!("#{} frame 2", traj_id.to_raw()));
        assert_eq!(scene.structure(traj_id).unwrap().frame, 2);
        assert_eq!(scene.structure(single_id).unwrap().frame, 0);

        // Explicit id still works and is undoable like every other verb.
        run_line(
            &mut scene,
            &mut history,
            &format!("frame 0 {}", traj_id.to_raw()),
        )
        .unwrap();
        assert_eq!(scene.structure(traj_id).unwrap().frame, 0);
        assert_eq!(
            run_line(&mut scene, &mut history, "undo").unwrap(),
            "undone"
        );
        assert_eq!(scene.structure(traj_id).unwrap().frame, 2);

        // Out of range and non-numeric are clear errors, not a panic.
        assert!(err_of(&mut scene, &mut history, "frame 99").contains("out of range"));
        assert!(err_of(&mut scene, &mut history, "frame x").contains("usage"));

        // No multi-frame structure at all: a clear error, not a panic or
        // a silent fall-through to the single-frame one.
        let mut scene2 = Scene::new();
        let mut history2 = CommandHistory::default();
        run_line(
            &mut scene2,
            &mut history2,
            &format!("load {}", fixture("1CRN.cif")),
        )
        .unwrap();
        assert!(err_of(&mut scene2, &mut history2, "frame 0").contains("no multi-frame structure"));
    }

    fn err_of(scene: &mut Scene, history: &mut CommandHistory, line: &str) -> String {
        run_line(scene, history, line).unwrap_err().0
    }

    #[test]
    fn sasa_matches_fastsasa_and_attaches_a_channel() {
        let mut scene = Scene::new();
        let mut history = CommandHistory::default();
        run_line(
            &mut scene,
            &mut history,
            &format!("load {}", fixture("1CRN.pdb")),
        )
        .unwrap();
        let total = |out: String| -> f32 {
            out.split(": ")
                .nth(1)
                .and_then(|s| s.split(' ').next())
                .unwrap()
                .parse()
                .unwrap()
        };
        // FastSASA on 1CRN.pdb: 3001.18 A^2 (fixtures/small/1CRN.sasa.pdb).
        let whole = total(run_line(&mut scene, &mut history, "sasa").unwrap());
        assert!((whole - 3001.18).abs() < 3.0, "{whole}");
        let loaded = scene.structures().next().unwrap().1;
        assert!(loaded.values.contains_key("sasa"));
        // A reported part, with the rest still in the calculation.
        let part = total(run_line(&mut scene, &mut history, "sasa | resid 1-10").unwrap());
        assert!(part > 0.0 && part < whole);
        assert!(run_line(&mut scene, &mut history, "sasa resname XYZ").is_err());
        assert!(run_line(&mut scene, &mut history, "sasa | resname XYZ").is_err());

        // Hydrogens count when the selection names them: FastSASA
        // `--hydrogen` on 1D3Z's first model gives 5035.61 A^2.
        run_line(
            &mut scene,
            &mut history,
            &format!("load {}", fixture("1D3Z.pdb")),
        )
        .unwrap();
        let with_h = total(run_line(&mut scene, &mut history, "sasa all").unwrap());
        assert!((with_h - 5035.61).abs() < 5.0, "{with_h}");
    }

    #[test]
    fn help_lists_every_command() {
        let all = help("");
        for s in SPECS {
            assert!(all.contains(s.usage), "{}", s.usage);
        }
        assert!(help("select").contains("docs/SELECTION.md"));
        assert!(help("nope").contains("unknown command"));
    }
}
