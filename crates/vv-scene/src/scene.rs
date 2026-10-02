//! The document: loaded structures, selections, and their commands.

use std::collections::BTreeMap;
use std::path::PathBuf;
use std::sync::Arc;

use vv_core::altloc::AltlocPolicy;
use vv_core::fixedbitset::FixedBitSet;
use vv_core::glam::Vec3;
use vv_core::{BondTable, CoordSet, ResidueClass, Structure};

use crate::coloring::{ColorOverride, Property, PropertyKind};
use crate::selection::{Mask, SelectionSet};
use crate::slotmap::{Id, SlotMap};
use crate::values::ValueChannel;

pub type StructureId = Id<LoadedStructure>;

/// Mirrors `vv_render::Representation`; kept separate so this crate has no
/// GPU dependency. `vv-app` maps between the two.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum Representation {
    #[default]
    Spacefill,
    BallAndStick,
    /// Backbone tube for polymers, sticks for ligands, no water
    /// (cartoon-lite; see `vv_core::backbone`).
    Tube,
    /// Ribbon/arrow geometry from DSSP secondary structure
    /// (`vv_core::cartoon`, `vv_core::dssp`): helices and strands as flat
    /// ribbons, everything else as a thin near-round one. Nucleotides get
    /// base plates; ligands, ions, glycans and lipids are drawn alongside
    /// (see docs/RENDERING.md), tuned by `repopt`.
    Cartoon,
    /// Ray-marched Gaussian ("blobby") molecular surface
    /// (`vv_core::gaussian_surface`, `vv_render::GaussianSurfaceGpu`): a
    /// full-screen pass, not per-atom impostors, drawn through
    /// `vv-app`'s own `gaussian_surfaces` list the same way `Cartoon` is
    /// drawn through `cartoons`.
    GaussianSurface,
    /// Skin surface (Edelsbrunner 1999; `vv_core::skin_surface`,
    /// `vv_core::feasible_cell`): the full mixed complex -- sphere,
    /// hyperboloid, and void patches for every vertex, edge, triangle and
    /// tetrahedron -- each ray-cast on its own impostor billboard and
    /// pickable by atom, drawn through `vv-app`'s own `skin_surfaces`
    /// list the same way `GaussianSurface` is drawn through
    /// `gaussian_surfaces`.
    SkinSurface,
    /// The solvent-excluded (molecular, Connolly) surface for a 1.4 A
    /// water probe, analytic (`vv_core::ses`), ray-cast per patch like
    /// `SkinSurface`.
    Ses,
    /// The solvent-accessible surface as spheres: atoms grown by the
    /// 1.4 A water probe.
    Sas,
    /// Thin bonds with balls of the same radius (licorice).
    Sticks,
    /// Bonds as 1-px lines.
    Lines,
    /// 3D-SNFG glyphs (`vv_core::glycan`) for glycan residues only,
    /// regardless of `selrep` -- like `Tube` restricting itself to
    /// backbone and ligand sticks. Port of the `3D-SNFG.tcl` script v1
    /// (Thieker, Hadden, Schulten & Woods, 2016, Glycobiology 26(8):
    /// 786-787, doi:10.1093/glycob/cww076); symbols: Varki et al. (2015)
    /// Glycobiology 25(12):1323-1324, doi:10.1093/glycob/cwv091.
    Glycan,
}

/// What a structure is made of: a material preset.
pub use vv_core::MaterialPreset as Material;

/// How a rep's atoms are colored. `vv_render::coloring::colors_of`
/// resolves it to one color per atom.
#[derive(Clone, Debug, PartialEq, Default)]
pub enum ColorScheme {
    #[default]
    Element,
    Chain,
    /// Secondary structure: the file's records for a
    /// single structure, DSSP per frame for a trajectory.
    SecondaryStructure,
    /// Acidic / basic / polar / nonpolar residues.
    ResidueType,
    /// Blue to red from N- to C-terminus along each chain.
    Rainbow,
    /// Carbons by chain, other atoms by element.
    Hetero,
    /// One hue per residue name.
    ResidueName,
    /// One hue per segment id (`Topology::segid`), or per `auth_asym_id`
    /// in a structure with no segment ids.
    SegmentName,
    /// One hue per connected component of the bond graph (a
    /// "fragment"); needs a `BondTable`, so unlike every other scheme here
    /// it has no `vv_render::color::ColorScheme` counterpart -- see
    /// `vv_render::colors_for_fragments`.
    Fragment,
    /// Zappo physicochemical grouping (Livingstone & Barton 1993).
    Zappo,
    /// Taylor (1997) physicochemical colour wheel.
    Taylor,
    /// Simplified (non-alignment) Clustal colouring (Thompson et al. 1997).
    Clustal,
    /// Nucleotide identity (A/C/G/T/U).
    Nucleotide,
    /// Purine (A, G) vs pyrimidine (C, T, U).
    PurinePyrimidine,
    /// One fixed colour per residue molecule class (`vv_core::ResidueClass`).
    MoleculeClass,
    /// Every atom one colour (sRGB bytes).
    Constant([u8; 3]),
    /// A continuous property on a ramp: B-factor, hydrophobicity, charge,
    /// a value channel, ...
    Property(Property),
}

/// The standard named colours, for `color NAME`.
pub const NAMED_COLORS: [(&str, [u8; 3]); 18] = [
    ("blue", [0, 0, 255]),
    ("red", [255, 0, 0]),
    ("gray", [89, 89, 89]),
    ("orange", [255, 128, 0]),
    ("yellow", [255, 255, 0]),
    ("tan", [128, 128, 89]),
    ("silver", [153, 153, 153]),
    ("green", [0, 255, 0]),
    ("white", [255, 255, 255]),
    ("pink", [255, 153, 153]),
    ("cyan", [64, 191, 191]),
    ("purple", [166, 0, 166]),
    ("lime", [128, 255, 0]),
    ("mauve", [255, 145, 200]),
    ("ochre", [128, 77, 13]),
    ("iceblue", [128, 153, 255]),
    ("black", [0, 0, 0]),
    ("magenta", [255, 0, 255]),
];

fn named_color(name: &str) -> Option<[u8; 3]> {
    let name = if name == "grey" { "gray" } else { name };
    NAMED_COLORS
        .iter()
        .find(|(n, _)| *n == name)
        .map(|(_, c)| *c)
}

impl From<PropertyKind> for ColorScheme {
    fn from(kind: PropertyKind) -> Self {
        ColorScheme::Property(Property::new(kind))
    }
}

impl ColorScheme {
    /// Every scheme that takes no parameters and is no property or
    /// solid color, for menus and tests.
    pub fn plain() -> Vec<ColorScheme> {
        use ColorScheme::*;
        vec![
            Element,
            Chain,
            SecondaryStructure,
            ResidueType,
            Rainbow,
            Hetero,
            ResidueName,
            SegmentName,
            Fragment,
            Zappo,
            Taylor,
            Clustal,
            Nucleotide,
            PurinePyrimidine,
            MoleculeClass,
        ]
    }

    /// The name the command language and session files use.
    pub fn name(&self) -> String {
        match self {
            ColorScheme::Element => "element".into(),
            ColorScheme::Chain => "chain".into(),
            ColorScheme::SecondaryStructure => "structure".into(),
            ColorScheme::ResidueType => "restype".into(),
            ColorScheme::Rainbow => "rainbow".into(),
            ColorScheme::Hetero => "hetero".into(),
            ColorScheme::ResidueName => "resname".into(),
            ColorScheme::SegmentName => "segname".into(),
            ColorScheme::Fragment => "fragment".into(),
            ColorScheme::Zappo => "zappo".into(),
            ColorScheme::Taylor => "taylor".into(),
            ColorScheme::Clustal => "clustal".into(),
            ColorScheme::Nucleotide => "nucleotide".into(),
            ColorScheme::PurinePyrimidine => "purinepyrimidine".into(),
            ColorScheme::MoleculeClass => "class".into(),
            ColorScheme::Constant(rgb) => match NAMED_COLORS.iter().find(|(_, c)| c == rgb) {
                Some((name, _)) => (*name).into(),
                None => format!("#{:02x}{:02x}{:02x}", rgb[0], rgb[1], rgb[2]),
            },
            ColorScheme::Property(p) => p.name(),
        }
    }

    /// Inverse of `name`, accepting every spelling a surface takes (the
    /// console's `values NAME` as well as the session file's `values:NAME`),
    /// so Python and session files can take exactly what the console does.
    pub fn parse(text: &str) -> Option<ColorScheme> {
        let words: Vec<&str> = text.split_whitespace().collect();
        match ColorScheme::parse_words(&words)? {
            (scheme, used) if used == words.len() => Some(scheme),
            _ => None,
        }
    }

    /// A coloring from the front of `words` and how many words it used;
    /// a property takes its `scale`/`range`/`ramp` words with it, and
    /// whatever follows (a structure id) is left alone.
    pub fn parse_words(words: &[&str]) -> Option<(ColorScheme, usize)> {
        if let Some((property, used)) = Property::parse_words(words) {
            return Some((ColorScheme::Property(property), used));
        }
        let scheme = match *words.first()? {
            "element" => ColorScheme::Element,
            "chain" => ColorScheme::Chain,
            "structure" | "ss" | "secondary" => ColorScheme::SecondaryStructure,
            "restype" | "residue_type" => ColorScheme::ResidueType,
            "rainbow" | "index" | "spectrum" => ColorScheme::Rainbow,
            "hetero" | "byhetero" | "cbc" => ColorScheme::Hetero,
            "resname" | "residue_name" => ColorScheme::ResidueName,
            "segname" | "segid" | "segment" => ColorScheme::SegmentName,
            "fragment" => ColorScheme::Fragment,
            "zappo" => ColorScheme::Zappo,
            "taylor" => ColorScheme::Taylor,
            "clustal" => ColorScheme::Clustal,
            "nucleotide" | "nuc" => ColorScheme::Nucleotide,
            "purinepyrimidine" | "purpyr" => ColorScheme::PurinePyrimidine,
            "class" | "moleculeclass" | "molclass" => ColorScheme::MoleculeClass,
            other => ColorScheme::Constant(parse_solid(other)?),
        };
        Some((scheme, 1))
    }
}

/// `#rrggbb` or a name from `NAMED_COLORS`.
pub fn parse_solid(word: &str) -> Option<[u8; 3]> {
    match word.strip_prefix('#') {
        Some(hex) if hex.len() == 6 => {
            let byte = |i: usize| u8::from_str_radix(hex.get(i..i + 2)?, 16).ok();
            Some([byte(0)?, byte(2)?, byte(4)?])
        }
        Some(_) => None,
        None => named_color(word),
    }
}

/// A rep's identity within its structure. Survives reordering, removal
/// and undo, so render caches follow the rep rather than its index.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct RepId(pub u32);

/// One way of drawing part of a structure (a "rep"): which atoms (a
/// selection expression, docs/SELECTION.md), drawn how, coloured how,
/// made of what.
#[derive(Clone, Debug, PartialEq)]
pub struct Rep {
    pub id: RepId,
    pub selection: String,
    pub representation: Representation,
    pub coloring: ColorScheme,
    pub material: Material,
    pub visible: bool,
    /// Options set away from their default (`Representation::options`),
    /// by name. Kept across a change of representation: an option name
    /// means the same wherever it appears (`probe` in SAS and SES).
    pub options: std::collections::BTreeMap<String, f32>,
}

/// A size or quality a representation is tuned by: its name in the
/// command language, the label the interface shows, and its range.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RepOption {
    pub name: &'static str,
    pub label: &'static str,
    pub default: f32,
    pub min: f32,
    pub max: f32,
    pub unit: &'static str,
    /// Names for the integer values `0..choices.len()` when the option is
    /// a choice rather than a size: `repopt NAME <choice>` and the
    /// interface use them; the stored value stays a plain float.
    pub choices: &'static [&'static str],
}

impl RepOption {
    /// The value a typed `text` means: a choice name, else a number.
    pub fn parse(&self, text: &str) -> Option<f32> {
        match self.choices.iter().position(|c| *c == text) {
            Some(i) => Some(i as f32),
            None => text.parse().ok(),
        }
    }

    /// `value` as the user would type it.
    pub fn show(&self, value: f32) -> String {
        match self.choices.get(value as usize) {
            Some(name) => (*name).to_string(),
            None => format!("{value}"),
        }
    }
}

const fn option(
    name: &'static str,
    label: &'static str,
    default: f32,
    min: f32,
    max: f32,
    unit: &'static str,
) -> RepOption {
    RepOption {
        name,
        label,
        default,
        min,
        max,
        unit,
        choices: &[],
    }
}

const fn choice(
    name: &'static str,
    label: &'static str,
    default: usize,
    choices: &'static [&'static str],
) -> RepOption {
    RepOption {
        name,
        label,
        default: default as f32,
        min: 0.0,
        max: (choices.len() - 1) as f32,
        unit: "",
        choices,
    }
}

const ON_OFF: &[&str] = &["off", "on"];

const SCALE: RepOption = option("scale", "Atom size", 1.0, 0.2, 1.5, "× vdW");
const BALL: RepOption = option("scale", "Ball size", 0.25, 0.1, 0.6, "× vdW");
const BOND: RepOption = option("bond", "Bond radius", 0.15, 0.05, 0.5, " Å");
const STICK: RepOption = option("bond", "Stick radius", 0.3, 0.1, 0.8, " Å");
const TUBE: RepOption = option("radius", "Tube radius", 0.3, 0.1, 1.0, " Å");
/// Putty's thin end: the drawn radius at the lowest B-factor of the
/// tube's atoms; `radius` is the thick end, at the highest.
const TUBE_RADIUS_MIN: RepOption =
    option("radius_min", "Min radius (thin end)", 0.1, 0.05, 1.0, " Å");
/// `0`: constant `radius`. `1`: putty -- `radius_min` to `radius` over the
/// B-factor range of the tube's atoms (`vv_core::backbone::putty_radius`).
/// A plain float like every other option, so `repopt`/sessions need no
/// special case; the UI shows it as a two-way choice, not a slider.
const RADIUS_BY: RepOption = option("radius_by", "Radius by B-factor", 0.0, 0.0, 1.0, "");
/// Scales putty's spread: the thick end becomes `radius_min + putty *
/// (radius - radius_min)`; 0 flattens it to `radius_min`.
const PUTTY: RepOption = option("putty", "Putty strength", 1.0, 0.0, 4.0, "×");
const PROBE: RepOption = option("probe", "Probe radius", 1.4, 0.5, 3.0, " Å");
const BLOB: RepOption = option("blob", "Blobbiness", 2.0, 0.5, 4.0, "");
const SHRINK: RepOption = option("shrink", "Shrink", 0.5, 0.2, 0.8, "");
/// A surface drawn solid, as the edges of its triangle mesh, or as the
/// mesh's vertices (`vv_core::{gaussian_mesh, ses_mesh, skin_mesh}`).
const SURFACE: RepOption = choice("surface", "Display", 0, &["solid", "mesh", "dots"]);
/// Mesh line width in pixels; a dot is `line / 10` Å in radius.
const MESH_LINE: RepOption = option("line", "Mesh line width", 1.5, 0.5, 6.0, " px");
/// 1.5 is the source script's "icon" preset, 4.0 its "full" preset; a
/// continuous option covers both and anything between or beyond.
const GLYCAN_SIZE: RepOption = option("size", "Shape size", 4.0, 1.5, 8.0, " Å");
/// Shares its name with `TUBE`'s, so switching a rep between `Tube` and
/// `Glycan` keeps a sensible cylinder radius; 0 hides linkages, as the
/// source script's icon preset does by zeroing its own cylinder radius.
const GLYCAN_RADIUS: RepOption = option("radius", "Linkage radius", 0.5, 0.0, 1.5, " Å");

/// How a cartoon or tube draws each nucleotide's base: a stick to the pairing
/// atom, a plate on its ring atoms (default), or one rung per base pair.
const BASES: RepOption = choice("bases", "Bases", 1, &["stick", "plate", "ladder"]);
/// A cartoon's helices as flat ribbons, or as straight cylinders along the
/// helix axis (`vv_core::cartoon::CartoonPlan::with_cylinder_helices`).
const HELIX: RepOption = choice("helix", "Helices", 0, &["ribbon", "cylinder"]);
const LIGANDS: RepOption = choice("ligands", "Ligands, cofactors", 1, ON_OFF);
const IONS: RepOption = choice("ions", "Ions", 1, ON_OFF);
const GLYCANS: RepOption = choice("glycans", "Glycans", 1, ON_OFF);
const LIPIDS: RepOption = choice("lipids", "Lipids", 1, ON_OFF);
const WATER: RepOption = choice("water", "Water", 0, ON_OFF);
const ADDITIVES: RepOption = choice("additives", "Crystallization additives", 0, ON_OFF);

impl Representation {
    /// Option `name` from the options set in `set`, or its default; `None`
    /// when this representation has no such option.
    pub fn option(self, set: &std::collections::BTreeMap<String, f32>, name: &str) -> Option<f32> {
        let spec = self.options().iter().find(|o| o.name == name)?;
        Some(set.get(name).copied().unwrap_or(spec.default))
    }

    /// What this representation can be tuned by; the defaults are what
    /// it draws with untouched.
    pub fn options(self) -> &'static [RepOption] {
        match self {
            Representation::Spacefill => &[SCALE],
            Representation::Sas => &[PROBE],
            Representation::Ses => &[PROBE, SURFACE, MESH_LINE],
            Representation::BallAndStick => &[BALL, BOND],
            Representation::Sticks => &[STICK],
            // The first two show outright (`rep_options_ui`): putty's
            // switch and strength, the tuning people reach for.
            Representation::Tube => &[
                RADIUS_BY,
                PUTTY,
                TUBE,
                TUBE_RADIUS_MIN,
                BASES,
                LIGANDS,
                IONS,
                GLYCANS,
                LIPIDS,
                WATER,
                ADDITIVES,
            ],
            Representation::Cartoon => &[
                HELIX, BASES, LIGANDS, IONS, GLYCANS, LIPIDS, WATER, ADDITIVES,
            ],
            Representation::GaussianSurface => &[BLOB, SURFACE, MESH_LINE],
            Representation::SkinSurface => &[SHRINK, SURFACE, MESH_LINE],
            Representation::Glycan => &[GLYCAN_SIZE, GLYCAN_RADIUS],
            Representation::Lines => &[],
        }
    }
}

impl Rep {
    /// Option `name` of this rep's representation: as set, or its
    /// default. `None` when the representation has no such option.
    pub fn option(&self, name: &str) -> Option<f32> {
        self.representation.option(&self.options, name)
    }

    /// Every atom, drawn as `representation`, default colouring and
    /// material.
    pub fn new(id: RepId, representation: Representation) -> Rep {
        Rep {
            id,
            selection: "all".into(),
            representation,
            coloring: ColorScheme::default(),
            material: Material::default(),
            visible: true,
            options: Default::default(),
        }
    }

    /// Whether the selection is every atom (no mask needed).
    pub fn selects_all(&self) -> bool {
        self.selection.trim() == "all"
    }
}

/// A cartoon where there is a polymer to trace, else lines (the cheapest
/// style that shows every bond at any size). The cartoon draws what it
/// cannot trace as licorice, so nothing a protein structure holds is lost.
fn default_representation(structure: &Structure) -> Representation {
    let counts = &structure.topology.class_counts;
    let polymer = counts[ResidueClass::Protein.index()] + counts[ResidueClass::Nucleic.index()];
    if polymer > 0 {
        Representation::Cartoon
    } else {
        Representation::Lines
    }
}

#[derive(Clone, Debug)]
pub struct LoadedStructure {
    pub structure: Structure,
    pub path: Option<PathBuf>,
    /// The trajectory file its frames stream from (`Command::
    /// LoadTrajectory`), for a saved session to load again.
    pub trajectory: Option<PathBuf>,
    pub label: String,
    /// Whole-structure show/hide (Structures panel row's eye), independent
    /// of any rep's own `visible`: hidden draws nothing and, during
    /// playback, is skipped by `GpuCache::sync` entirely (no coordinate
    /// upload, spline/surface/glycan rebuild or occlusion contribution) --
    /// `frame` below still advances, so showing it again catches up in one
    /// sync rather than replaying every missed frame.
    pub visible: bool,
    /// Which alternate location of each residue reps draw
    /// (`vv_core::altloc`); every conformer stays in the data.
    pub altloc: AltlocPolicy,
    /// How the structure is drawn, bottom to top; never empty.
    pub reps: Vec<Rep>,
    /// The rep that `representation`, `color` and `material` edit.
    pub current_rep: usize,
    /// The next `RepId` to hand out.
    pub next_rep_id: u32,
    /// Perceived eagerly at load time; see docs/RENDERING.md for why this
    /// is cheap enough not to bother making lazy yet.
    pub bonds: BondTable,
    /// Per-atom value channels from outside (`Command::SetValues`), by
    /// name. Persisted as `.npy` sidecars next to a saved session
    /// (`vv_scene::session`).
    pub values: BTreeMap<String, ValueChannel>,
    /// Colors laid over every rep's coloring, in order, later winning
    /// (`Command::SetColorOverrides`).
    pub color_overrides: Vec<ColorOverride>,
    /// Text labels anchored to specific atoms (`Command::SetLabel`),
    /// billboarded in the viewport and following that atom's current
    /// frame position. Keyed by atom index. Persisted in session files
    /// (`vv_scene::session::SavedStructure::labels`).
    pub labels: BTreeMap<u32, String>,
    /// Distances, angles and dihedrals kept on screen
    /// (`Command::SetMeasurement`), re-measured at the current frame.
    pub measurements: Vec<Measurement>,
    /// Contact overlays switched on (`Command::SetInteraction`).
    pub interactions: std::collections::BTreeSet<vv_core::interactions::InteractionKind>,
    /// Which coordinate set (of `structure.frame_count()`) is current.
    /// Always `0` for a single-frame structure. Two ways to change it,
    /// deliberately: `Command::SetFrame` (undoable, for a script/console/
    /// Python's deliberate one-shot set — this is the value a saved
    /// session persists) and [`Scene::set_frame_live`] (not undoable, for
    /// playback and scrubbing, which change many times a second and would
    /// flood undo history if each step were a `Command` — the same
    /// reasoning that already keeps camera moves off the command bus).
    pub frame: usize,
}

impl LoadedStructure {
    /// A freshly loaded structure: one rep over every atom, drawn as
    /// `default_representation`.
    pub fn new(
        structure: Structure,
        path: Option<PathBuf>,
        label: String,
        bonds: BondTable,
    ) -> Self {
        let representation = default_representation(&structure);
        LoadedStructure {
            structure,
            path,
            trajectory: None,
            label,
            visible: true,
            altloc: AltlocPolicy::default(),
            reps: vec![Rep::new(RepId(0), representation)],
            current_rep: 0,
            next_rep_id: 1,
            bonds,
            values: Default::default(),
            color_overrides: Vec::new(),
            labels: Default::default(),
            measurements: Vec::new(),
            interactions: Default::default(),
            frame: 0,
        }
    }

    /// The atoms the altloc policy displays, or `None` when all are. The
    /// one mask every drawing, picking and selecting path shares.
    pub fn shown_atoms(&self) -> Option<FixedBitSet> {
        vv_core::altloc::visible_atoms(&self.structure.topology, self.altloc)
    }

    /// `expr` over `frame`, its `shown` keyword following the altloc policy.
    pub fn select(
        &self,
        expr: &str,
        frame: usize,
    ) -> Result<vv_core::fixedbitset::FixedBitSet, vv_core::SelectError> {
        vv_core::select_under(
            &self.structure.topology,
            self.structure.frame(frame).positions(),
            expr,
            self.altloc,
        )
    }

    /// `atoms` without the ones the altloc policy hides.
    pub fn only_shown(&self, atoms: &[u32]) -> Vec<u32> {
        let shown = self.shown_atoms();
        atoms
            .iter()
            .copied()
            .filter(|&a| shown.as_ref().is_none_or(|s| s.contains(a as usize)))
            .collect()
    }

    /// Coordinates of `frame` with each hidden conformer atom moved onto
    /// the shown one it stands in for (`vv_core::altloc::stand_ins`), so
    /// geometry found by atom name follows the displayed conformer.
    pub fn drawn_positions(&self, frame: usize) -> Arc<CoordSet> {
        let coords = self.structure.frame(frame);
        let Some(shown) = self.shown_atoms() else {
            return coords;
        };
        let mut out: Vec<Vec3> = coords.positions().to_vec();
        for (hidden, twin) in vv_core::altloc::stand_ins(&self.structure.topology, &shown) {
            out[hidden as usize] = out[twin as usize];
        }
        Arc::new(CoordSet::new(out))
    }

    /// `keep` (per-atom, `None` for all) with each hidden trace atom
    /// answering for its stand-in, for filters over name-found atoms.
    pub fn trace_keep(&self, keep: &Option<Vec<bool>>) -> Option<Vec<bool>> {
        let mut keep = keep.clone()?;
        if let Some(shown) = self.shown_atoms() {
            for (hidden, twin) in vv_core::altloc::stand_ins(&self.structure.topology, &shown) {
                keep[hidden as usize] = keep[twin as usize];
            }
        }
        Some(keep)
    }

    /// `pairs` with each hidden conformer atom replaced by the shown atom
    /// it stands in for, for bonds between name-found atoms.
    pub fn as_shown(&self, pairs: &[[u32; 2]]) -> Vec<[u32; 2]> {
        let Some(shown) = self.shown_atoms() else {
            return pairs.to_vec();
        };
        let twins: std::collections::HashMap<u32, u32> =
            vv_core::altloc::stand_ins(&self.structure.topology, &shown)
                .into_iter()
                .collect();
        let swap = |a: u32| twins.get(&a).copied().unwrap_or(a);
        pairs.iter().map(|&[a, b]| [swap(a), swap(b)]).collect()
    }

    /// The rep that edits apply to.
    pub fn rep(&self) -> &Rep {
        &self.reps[self.current_rep]
    }

    pub fn rep_index(&self, id: RepId) -> Option<usize> {
        self.reps.iter().position(|r| r.id == id)
    }

    /// The value channel `coloring` refers to, if it is attached.
    pub fn values_for<'a>(
        &'a self,
        coloring: &'a ColorScheme,
    ) -> Option<(&'a str, &'a ValueChannel)> {
        match coloring {
            ColorScheme::Property(Property {
                kind: PropertyKind::Values(name),
                ..
            }) => self.values.get(name).map(|c| (name.as_str(), c)),
            _ => None,
        }
    }

    /// Each atom's override color at `frame` (`None` where no override
    /// reaches it), later overrides winning; empty when there are none.
    /// A target that no longer selects (a stale expression) is skipped.
    pub fn override_colors(&self, frame: usize) -> Vec<Option<[u8; 3]>> {
        if self.color_overrides.is_empty() {
            return Vec::new();
        }
        let mut out = vec![None; self.structure.atom_count()];
        for o in &self.color_overrides {
            let Ok(atoms) = self.select(&o.target.expression(), frame) else {
                continue;
            };
            for a in atoms.ones() {
                out[a] = Some(o.color);
            }
        }
        out
    }
}

/// A distance (2 atoms), angle (3) or dihedral (4, IUPAC sign) between
/// atoms of one structure: shown in
/// the viewport and measured again at whatever frame is current.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Measurement(Vec<u32>);

impl Measurement {
    /// `None` unless there are 2, 3 or 4 atoms.
    pub fn new(atoms: Vec<u32>) -> Option<Self> {
        (2..=4).contains(&atoms.len()).then_some(Self(atoms))
    }

    pub fn atoms(&self) -> &[u32] {
        &self.0
    }

    /// "distance", "angle" or "dihedral".
    pub fn kind(&self) -> &'static str {
        ["distance", "angle", "dihedral"][self.0.len() - 2]
    }

    /// Its value at `positions`: Angstrom for a distance, degrees
    /// otherwise.
    pub fn value(&self, positions: &[vv_core::glam::Vec3]) -> f32 {
        use vv_core::analysis::{angle, dihedral, distance};
        match self.0.iter().map(|&a| a as usize).collect::<Vec<_>>()[..] {
            [a, b] => distance(positions, a, b),
            [a, b, c] => angle(positions, a, b, c),
            [a, b, c, d] => dihedral(positions, a, b, c, d),
            _ => unreachable!("2 to 4 atoms"),
        }
    }

    /// The value with its unit, as the viewport shows it.
    pub fn text(&self, positions: &[vv_core::glam::Vec3]) -> String {
        match self.0.len() {
            2 => format!("{:.2} A", self.value(positions)),
            _ => format!("{:.1} deg", self.value(positions)),
        }
    }
}

/// The transient "current" selection driving the inspector, as opposed to
/// a saved, named `SelectionSet`.
#[derive(Clone, Debug, PartialEq)]
pub struct ActiveSelection {
    pub structure: StructureId,
    pub mask: Mask,
    /// The expression this selection was made from (`Command::SelectExpr`),
    /// or `None` for a pick. Carried into a saved set so it can be shown,
    /// edited, and re-applied.
    pub expr: Option<String>,
}

/// A free-floating text caption for figure composition, not tied to any
/// structure coordinate (docs/UI_DESIGN.md's "screen-space" annotation
/// type, as opposed to the structure-anchored kind on
/// `LoadedStructure::labels`). Position is fractional (0..1 from the
/// viewport's top-left) so it lands in the same place at any window size.
#[derive(Clone, Debug, PartialEq)]
pub struct Caption {
    pub name: String,
    pub x: f32,
    pub y: f32,
    pub text: String,
}

#[derive(Clone, Debug, Default)]
pub struct Scene {
    structures: SlotMap<LoadedStructure>,
    selection_sets: Vec<SelectionSet>,
    captions: Vec<Caption>,
    active: Option<ActiveSelection>,
    /// Bumped on every applied command; panels compare against their last
    /// seen value to know whether to recompute cached UI state.
    version: u64,
}

impl Scene {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn version(&self) -> u64 {
        self.version
    }

    pub fn structure(&self, id: StructureId) -> Option<&LoadedStructure> {
        self.structures.get(id)
    }

    pub fn structures(&self) -> impl DoubleEndedIterator<Item = (StructureId, &LoadedStructure)> {
        self.structures.iter()
    }

    pub fn active_selection(&self) -> Option<&ActiveSelection> {
        self.active.as_ref()
    }

    pub fn selection_sets(&self) -> &[SelectionSet] {
        &self.selection_sets
    }

    pub fn selection_set(&self, name: &str) -> Option<&SelectionSet> {
        self.selection_sets.iter().find(|s| s.name == name)
    }

    pub fn captions(&self) -> &[Caption] {
        &self.captions
    }

    pub(crate) fn structures_mut(&mut self) -> &mut SlotMap<LoadedStructure> {
        &mut self.structures
    }

    pub(crate) fn captions_mut(&mut self) -> &mut Vec<Caption> {
        &mut self.captions
    }

    pub(crate) fn set_active(
        &mut self,
        active: Option<ActiveSelection>,
    ) -> Option<ActiveSelection> {
        std::mem::replace(&mut self.active, active)
    }

    pub(crate) fn selection_sets_mut(&mut self) -> &mut Vec<SelectionSet> {
        &mut self.selection_sets
    }

    pub(crate) fn touch(&mut self) {
        self.version += 1;
    }

    /// Sets `id`'s current frame directly, bypassing the command bus (no
    /// undo entry) — see [`LoadedStructure::frame`]'s doc for why this
    /// exists alongside `Command::SetFrame` rather than instead of it.
    /// Clamped to a valid index; a missing `id` or a single-frame
    /// structure is silently a no-op, matching how a camera move on a
    /// closed/nonexistent structure would have nothing to act on either.
    pub fn set_frame_live(&mut self, id: StructureId, frame: usize) {
        if let Some(loaded) = self.structures.get_mut(id) {
            let last = loaded.structure.frame_count().saturating_sub(1);
            loaded.frame = frame.min(last);
        }
        self.touch();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn colorings_round_trip_by_name() {
        let tuned = Property {
            kind: PropertyKind::Charge,
            ramp: Some(crate::coloring::Ramp::Gray),
            range: Some([-2.0, 0.5]),
        };
        let mut schemes = ColorScheme::plain();
        schemes.extend(PropertyKind::builtin().into_iter().map(ColorScheme::from));
        schemes.extend([
            ColorScheme::Constant([255, 0, 0]),
            ColorScheme::Constant([1, 2, 3]),
            ColorScheme::from(PropertyKind::Values("sasa".into())),
            ColorScheme::Property(tuned),
        ]);
        for scheme in schemes {
            assert_eq!(
                ColorScheme::parse(&scheme.name()),
                Some(scheme.clone()),
                "{scheme:?}"
            );
        }
        assert_eq!(ColorScheme::Constant([255, 0, 0]).name(), "red");
        assert_eq!(ColorScheme::parse("grey"), ColorScheme::parse("gray"));
        assert_eq!(ColorScheme::parse("#zzzzzz"), None);
    }
}
