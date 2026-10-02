//! One packed color per atom of a scene structure: a rep's coloring (a
//! scheme, or a property on a ramp over a range), then the structure's
//! per-target overrides on top. Shared by the app and Python so both draw
//! the same colors, and the source of each property coloring's legend.

use std::borrow::Cow;

use vv_scene::{ColorScheme as Scheme, HydroScale, LoadedStructure, Property, PropertyKind, Ramp};

use crate::color::{self, ramp_color, rgba, ColorScheme as Builtin};
use crate::scene::{colors_for, colors_for_fragments, colors_for_ss, scalar_range};

/// One packed color per atom of `loaded` under `coloring` at `frame`, with
/// the structure's color overrides laid over it. A `Values` property whose
/// channel is not attached draws by element.
pub fn colors_of(loaded: &LoadedStructure, coloring: &Scheme, frame: usize) -> Vec<u32> {
    let mut colors = scheme_colors(loaded, coloring, frame);
    let overrides = loaded.override_colors(frame);
    for (color, over) in colors.iter_mut().zip(overrides) {
        if let Some([r, g, b]) = over {
            *color = rgba(r, g, b);
        }
    }
    colors
}

fn scheme_colors(loaded: &LoadedStructure, coloring: &Scheme, frame: usize) -> Vec<u32> {
    let topology = &loaded.structure.topology;
    match coloring {
        Scheme::Property(p) => property_colors(loaded, p, frame),
        Scheme::SecondaryStructure => {
            let positions = loaded.structure.frame(frame);
            let single = loaded.structure.frame_count() == 1;
            let codes =
                vv_core::cartoon::secondary_structure(topology, positions.positions(), single);
            colors_for_ss(topology, &codes)
        }
        Scheme::Fragment => colors_for_fragments(topology, &loaded.bonds),
        plain => colors_for(builtin(plain), topology),
    }
}

/// The renderer's scheme for a plain (parameterless, whole-structure)
/// coloring; the ones needing more input are handled by the caller.
fn builtin(scheme: &Scheme) -> Builtin {
    match scheme {
        Scheme::Chain => Builtin::Chain,
        Scheme::ResidueType => Builtin::ResidueType,
        Scheme::Rainbow => Builtin::Rainbow,
        Scheme::Hetero => Builtin::Hetero,
        Scheme::ResidueName => Builtin::ResidueName,
        Scheme::SegmentName => Builtin::SegmentName,
        Scheme::Zappo => Builtin::Zappo,
        Scheme::Taylor => Builtin::Taylor,
        Scheme::Clustal => Builtin::Clustal,
        Scheme::Nucleotide => Builtin::Nucleotide,
        Scheme::PurinePyrimidine => Builtin::PurinePyrimidine,
        Scheme::MoleculeClass => Builtin::Class,
        Scheme::Constant([r, g, b]) => Builtin::Constant(rgba(*r, *g, *b)),
        Scheme::Element | Scheme::SecondaryStructure | Scheme::Fragment | Scheme::Property(_) => {
            Builtin::Element
        }
    }
}

type Table = fn(&str) -> Option<f32>;

/// Higher is more hydrophobic, so Wimley-White's transfer free energy is
/// negated.
fn wimley_white_hydrophobicity(residue: &str) -> Option<f32> {
    color::wimley_white(residue).map(|v| -v)
}

/// The per-residue table behind a property, if it has one.
fn residue_table(kind: &PropertyKind) -> Option<Table> {
    Some(match kind {
        PropertyKind::Hydrophobicity(HydroScale::KyteDoolittle) => color::kyte_doolittle,
        PropertyKind::Hydrophobicity(HydroScale::WimleyWhite) => wimley_white_hydrophobicity,
        PropertyKind::Hydrophobicity(HydroScale::Eisenberg) => color::eisenberg,
        PropertyKind::Charge => color::charge,
        PropertyKind::HelixPropensity => color::helix_propensity,
        PropertyKind::StrandPropensity => color::strand_propensity,
        PropertyKind::TurnPropensity => color::turn_propensity,
        PropertyKind::BuriedIndex => color::buried_index,
        PropertyKind::BFactor | PropertyKind::Occupancy | PropertyKind::Values(_) => return None,
    })
}

/// A property's values for every atom; NaN where it has none (a residue
/// outside its table, or a missing number). `None` for a value channel
/// that is not attached.
fn atom_values<'a>(
    loaded: &'a LoadedStructure,
    kind: &PropertyKind,
    frame: usize,
) -> Option<Cow<'a, [f32]>> {
    let topology = &loaded.structure.topology;
    if let Some(table) = residue_table(kind) {
        let per_residue: Vec<f32> = (0..topology.residues.len())
            .map(|r| table(topology.residue_name(r)).unwrap_or(f32::NAN))
            .collect();
        let per_atom = topology.residue_index.iter();
        return Some(per_atom.map(|&r| per_residue[r as usize]).collect());
    }
    match kind {
        PropertyKind::BFactor => Some(Cow::Borrowed(&topology.b_factor[..])),
        PropertyKind::Occupancy => Some(Cow::Borrowed(&topology.occupancy[..])),
        PropertyKind::Values(name) => loaded
            .values
            .get(name)
            .map(|channel| Cow::Borrowed(channel.frame(frame))),
        _ => None,
    }
}

/// The range a property spans unless told otherwise: a published scale's
/// own, or the data's. `None` for a value channel that is not attached.
fn default_range(loaded: &LoadedStructure, kind: &PropertyKind) -> Option<[f32; 2]> {
    let topology = &loaded.structure.topology;
    let symmetric = |m: f32| [-m, m];
    Some(match kind {
        PropertyKind::Hydrophobicity(HydroScale::KyteDoolittle) => symmetric(4.5),
        PropertyKind::Hydrophobicity(HydroScale::WimleyWhite) => symmetric(3.64),
        PropertyKind::Hydrophobicity(HydroScale::Eisenberg) => symmetric(2.53),
        PropertyKind::Charge => symmetric(1.0),
        PropertyKind::HelixPropensity => [0.57, 1.51],
        PropertyKind::StrandPropensity => [0.37, 1.70],
        PropertyKind::TurnPropensity => [0.47, 1.56],
        PropertyKind::BuriedIndex => [0.05, 4.6],
        PropertyKind::BFactor => scalar_range(&topology.b_factor).into(),
        PropertyKind::Occupancy => scalar_range(&topology.occupancy).into(),
        PropertyKind::Values(name) => loaded.values.get(name)?.range().into(),
    })
}

/// Diverging ramps where zero or the midpoint means something, a
/// sequential one for scales that only run from low to high.
fn default_ramp(kind: &PropertyKind) -> Ramp {
    match kind {
        PropertyKind::Hydrophobicity(_) => Ramp::TealWhiteGold,
        PropertyKind::Charge => Ramp::RedWhiteBlue,
        PropertyKind::BFactor | PropertyKind::Occupancy | PropertyKind::Values(_) => {
            Ramp::BlueWhiteRed
        }
        PropertyKind::HelixPropensity
        | PropertyKind::StrandPropensity
        | PropertyKind::TurnPropensity
        | PropertyKind::BuriedIndex => Ramp::Viridis,
    }
}

fn property_colors(loaded: &LoadedStructure, p: &Property, frame: usize) -> Vec<u32> {
    let topology = &loaded.structure.topology;
    let (Some(values), Some(range)) = (
        atom_values(loaded, &p.kind, frame),
        default_range(loaded, &p.kind),
    ) else {
        return colors_for(Builtin::Element, topology);
    };
    let [lo, hi] = p.range.unwrap_or(range);
    let ramp = p.ramp.unwrap_or_else(|| default_ramp(&p.kind));
    let by_residue = residue_table(&p.kind).is_some();
    values
        .iter()
        .enumerate()
        .map(|(atom, &v)| match v.is_nan() && by_residue {
            true => color::by_element(topology.element[atom]),
            false => ramp_color(ramp, (v - lo) / (hi - lo)),
        })
        .collect()
}

/// What a continuous coloring's legend shows: its ramp from `min` to `max`.
#[derive(Clone, Debug, PartialEq)]
pub struct Legend {
    pub title: String,
    pub ramp: Ramp,
    pub min: f32,
    pub max: f32,
    pub low: &'static str,
    pub high: &'static str,
}

/// The legend of `coloring`, or `None` for anything that is not a
/// property with its values available (a value channel not yet attached).
pub fn legend(loaded: &LoadedStructure, coloring: &Scheme) -> Option<Legend> {
    let Scheme::Property(p) = coloring else {
        return None;
    };
    let [min, max] = p.range.or(default_range(loaded, &p.kind))?;
    let (low, high) = match p.kind {
        PropertyKind::Hydrophobicity(_) => ("hydrophilic", "hydrophobic"),
        PropertyKind::Charge => ("negative", "positive"),
        PropertyKind::BuriedIndex => ("exposed", "buried"),
        _ => ("low", "high"),
    };
    Some(Legend {
        title: p.kind.label(),
        ramp: p.ramp.unwrap_or_else(|| default_ramp(&p.kind)),
        min,
        max,
        low,
        high,
    })
}

/// The range and ramp a property uses when none is given, for the
/// interface's reset buttons.
pub fn property_defaults(
    loaded: &LoadedStructure,
    kind: &PropertyKind,
) -> Option<([f32; 2], Ramp)> {
    Some((default_range(loaded, kind)?, default_ramp(kind)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use vv_scene::{CommandHistory, Scene};

    fn load_1crn() -> Scene {
        let mut scene = Scene::new();
        let mut history = CommandHistory::default();
        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../../fixtures/small/1CRN.cif");
        vv_scene::run_line(&mut scene, &mut history, &format!("load {path}")).unwrap();
        scene
    }

    fn run(scene: &mut Scene, line: &str) {
        let mut history = CommandHistory::default();
        vv_scene::run_line(scene, &mut history, line).unwrap();
    }

    fn first(scene: &Scene) -> &LoadedStructure {
        scene.structures().next().unwrap().1
    }

    fn property(text: &str) -> Scheme {
        Scheme::parse(text).unwrap()
    }

    #[test]
    fn hydrophobicity_colors_follow_the_chosen_scale_and_ramp() {
        let scene = load_1crn();
        let loaded = first(&scene);
        let topology = &loaded.structure.topology;
        let residue_color = |scheme: &str, name: &str| {
            let colors = colors_of(loaded, &property(scheme), 0);
            let atom = (0..topology.atom_count())
                .find(|&a| topology.residue_name(topology.residue_index[a] as usize) == name)
                .unwrap();
            colors[atom]
        };
        // Arg is the most hydrophilic and Ile the most hydrophobic residue
        // of Kyte & Doolittle's table: the ends of the default ramp.
        assert_eq!(
            residue_color("hydrophobicity", "ARG"),
            ramp_color(Ramp::TealWhiteGold, 0.0)
        );
        assert_eq!(
            residue_color("hydrophobicity", "ILE"),
            ramp_color(Ramp::TealWhiteGold, 1.0)
        );
        assert_eq!(
            residue_color("hydrophobicity ramp viridis", "ILE"),
            ramp_color(Ramp::Viridis, 1.0)
        );
        // A narrower range saturates sooner.
        assert_eq!(
            residue_color("hydrophobicity range -1 1", "LEU"),
            ramp_color(Ramp::TealWhiteGold, 1.0)
        );
        // Wimley-White runs the other way round, so Phe is hydrophobic.
        assert_eq!(
            residue_color("hydrophobicity scale ww", "PHE"),
            ramp_color(Ramp::TealWhiteGold, (1.71 + 3.64) / 7.28)
        );
        // A residue outside the table keeps its element color.
        let water = (0..topology.atom_count())
            .find(|&a| topology.residue_name(topology.residue_index[a] as usize) == "HOH");
        if let Some(a) = water {
            let colors = colors_of(loaded, &property("hydrophobicity"), 0);
            assert_eq!(colors[a], color::by_element(topology.element[a]));
        }
    }

    #[test]
    fn the_legend_reports_the_range_in_force() {
        let scene = load_1crn();
        let legend_of = |scene: &Scene, text: &str| legend(first(scene), &property(text)).unwrap();

        let kd = legend_of(&scene, "hydrophobicity");
        assert_eq!((kd.min, kd.max, kd.ramp), (-4.5, 4.5, Ramp::TealWhiteGold));
        let ww = legend_of(&scene, "hydrophobicity scale ww range -2 3 ramp gray");
        assert_eq!((ww.min, ww.max, ww.ramp), (-2.0, 3.0, Ramp::Gray));

        let topology = &first(&scene).structure.topology;
        let (lo, hi) = scalar_range(&topology.b_factor);
        let b = legend_of(&scene, "b_factor");
        assert_eq!(
            (b.min, b.max),
            (lo, hi),
            "B-factor spans the structure's own"
        );
        assert!(hi > lo);

        assert!(legend(first(&scene), &Scheme::Chain).is_none());
        assert!(
            legend(first(&scene), &property("values:nothing")).is_none(),
            "no channel, no legend"
        );
    }

    #[test]
    fn overrides_win_over_every_scheme_and_leave_the_rest() {
        let mut scene = load_1crn();
        run(&mut scene, "color set name CA #ff00ff");
        let loaded = first(&scene);
        let ca = loaded.select("name CA", 0).unwrap();
        for scheme in [Scheme::Element, Scheme::Chain, property("hydrophobicity")] {
            let base = scheme_colors(loaded, &scheme, 0);
            let colors = colors_of(loaded, &scheme, 0);
            for (a, (&c, &b)) in colors.iter().zip(&base).enumerate() {
                let expected = if ca.contains(a) { rgba(255, 0, 255) } else { b };
                assert_eq!(c, expected, "{scheme:?} atom {a}");
            }
        }
    }
}
