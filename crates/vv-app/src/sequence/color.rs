//! Residue coloring of the sequence: a scheme picks one color per residue
//! of a structure (`Color32::TRANSPARENT` for none).

use egui::Color32;
use vv_core::dssp::DsspCode;
use vv_render::color as vc;
use vv_scene::LoadedStructure;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum SeqColor {
    #[default]
    None,
    /// The 3D view's coloring of the structure's current representation.
    View,
    SecondaryStructure,
    Chemistry,
    Hydrophobicity,
    BFactor,
    Sasa,
    Charge,
    Clustal,
    Zappo,
    Taylor,
}

/// `(scheme, command word, menu label)`.
pub const SCHEMES: [(SeqColor, &str, &str); 11] = [
    (SeqColor::None, "none", "None"),
    (SeqColor::View, "view", "Same as 3D view"),
    (SeqColor::SecondaryStructure, "ss", "Secondary structure"),
    (SeqColor::Chemistry, "chemistry", "Residue chemistry"),
    (SeqColor::Hydrophobicity, "hydrophobicity", "Hydrophobicity"),
    (SeqColor::BFactor, "bfactor", "B-factor"),
    (SeqColor::Sasa, "sasa", "Burial (SASA)"),
    (SeqColor::Charge, "charge", "Charge"),
    (SeqColor::Clustal, "clustal", "Clustal"),
    (SeqColor::Zappo, "zappo", "Zappo"),
    (SeqColor::Taylor, "taylor", "Taylor"),
];

impl SeqColor {
    pub fn word(self) -> &'static str {
        SCHEMES.iter().find(|s| s.0 == self).unwrap().1
    }

    pub fn parse(word: &str) -> Option<SeqColor> {
        let word = word.to_ascii_lowercase();
        let alias = match word.as_str() {
            "3d" | "same" | "structure" => "view",
            "secondary" => "ss",
            "type" | "residue" => "chemistry",
            "hydro" | "kd" => "hydrophobicity",
            "b" => "bfactor",
            "burial" | "exposure" => "sasa",
            other => other,
        };
        SCHEMES.iter().find(|s| s.1 == alias).map(|s| s.0)
    }

    /// Whether the colors name discrete classes (drawn as soft tints)
    /// rather than sample a continuous ramp.
    pub fn categorical(self) -> bool {
        !matches!(
            self,
            SeqColor::None | SeqColor::Hydrophobicity | SeqColor::BFactor | SeqColor::Sasa
        )
    }

    /// Whether the colors change with the frame.
    pub fn per_frame(self) -> bool {
        matches!(
            self,
            SeqColor::View | SeqColor::SecondaryStructure | SeqColor::Sasa
        )
    }
}

/// What a scheme reads.
pub struct ColorInput<'a> {
    pub loaded: &'a LoadedStructure,
    pub frame: usize,
    pub positions: &'a [vv_core::glam::Vec3],
    /// Relative SASA per residue, when the scheme needs it.
    pub rel_sasa: Option<&'a [f32]>,
}

pub fn from_packed(packed: u32) -> Color32 {
    Color32::from_rgb(packed as u8, (packed >> 8) as u8, (packed >> 16) as u8)
}

/// Share of a categorical color kept when it is mixed into the panel.
const TINT: f32 = 0.38;
/// WCAG AA for body text.
const MIN_CONTRAST: f32 = 4.5;

/// `fill` as a soft tint over `panel`.
pub fn tint(fill: Color32, panel: Color32) -> Color32 {
    let mix = |a: u8, b: u8| (b as f32 + (a as f32 - b as f32) * TINT).round() as u8;
    Color32::from_rgb(
        mix(fill.r(), panel.r()),
        mix(fill.g(), panel.g()),
        mix(fill.b(), panel.b()),
    )
}

/// WCAG 2 relative luminance of an sRGB color.
fn luminance(c: Color32) -> f32 {
    let lin = |v: u8| {
        let v = v as f32 / 255.0;
        if v <= 0.04045 {
            v / 12.92
        } else {
            ((v + 0.055) / 1.055).powf(2.4)
        }
    };
    0.2126 * lin(c.r()) + 0.7152 * lin(c.g()) + 0.0722 * lin(c.b())
}

pub fn contrast_ratio(a: Color32, b: Color32) -> f32 {
    let (la, lb) = (luminance(a), luminance(b));
    (la.max(lb) + 0.05) / (la.min(lb) + 0.05)
}

/// `preferred` text if it reads on `fill` (WCAG AA), else black or white,
/// whichever contrasts more.
pub fn readable_on(fill: Color32, preferred: Color32) -> Color32 {
    if contrast_ratio(preferred, fill) >= MIN_CONTRAST {
        return preferred;
    }
    if contrast_ratio(Color32::BLACK, fill) >= contrast_ratio(Color32::WHITE, fill) {
        Color32::BLACK
    } else {
        Color32::WHITE
    }
}

pub fn residue_colors(scheme: SeqColor, input: &ColorInput) -> Vec<Color32> {
    let top = &input.loaded.structure.topology;
    match scheme {
        SeqColor::None => vec![Color32::TRANSPARENT; top.residue_count()],
        SeqColor::View => view_colors(input),
        SeqColor::SecondaryStructure => secondary_colors(input),
        SeqColor::BFactor => b_factor_colors(input),
        SeqColor::Sasa => sasa_colors(input),
        _ => by_name(scheme, top),
    }
}

/// The 3D view's color of each residue's CA (else first) atom.
fn view_colors(input: &ColorInput) -> Vec<Color32> {
    let loaded = input.loaded;
    let top = &loaded.structure.topology;
    let coloring = &loaded.reps[loaded.current_rep.min(loaded.reps.len() - 1)].coloring;
    let atoms = crate::gpu_cache::colors_of(loaded, coloring, input.frame);
    top.residues
        .iter()
        .map(|rec| {
            let ca = rec
                .atoms
                .clone()
                .find(|&a| top.atom_name(a as usize) == "CA");
            from_packed(atoms[ca.unwrap_or(rec.atoms.start) as usize])
        })
        .collect()
}

/// Helix, strand and turn get their cartoon colors; coil stays plain.
fn secondary_colors(input: &ColorInput) -> Vec<Color32> {
    let top = &input.loaded.structure.topology;
    let single = input.loaded.structure.frame_count() == 1;
    vv_core::cartoon::secondary_structure(top, input.positions, single)
        .into_iter()
        .map(|code| match code {
            DsspCode::Bend | DsspCode::Coil | DsspCode::None => Color32::TRANSPARENT,
            _ => from_packed(vc::by_secondary_structure(code)),
        })
        .collect()
}

fn b_factor_colors(input: &ColorInput) -> Vec<Color32> {
    let top = &input.loaded.structure.topology;
    let means: Vec<Option<f32>> = top
        .residues
        .iter()
        .enumerate()
        .map(|(r, rec)| {
            let polymer = crate::sequence::rows::is_polymer(top, r as u32);
            let sum: f32 = rec.atoms.clone().map(|a| top.b_factor[a as usize]).sum();
            (polymer && !rec.atoms.is_empty()).then(|| sum / rec.atoms.len() as f32)
        })
        .collect();
    let known = means.iter().flatten().copied();
    let lo = known.clone().fold(f32::MAX, f32::min);
    let hi = known.fold(f32::MIN, f32::max);
    means
        .iter()
        .map(|m| {
            m.map_or(Color32::TRANSPARENT, |b| {
                from_packed(vc::by_b_factor(b, lo, hi))
            })
        })
        .collect()
}

/// White (fully exposed) to blue (buried).
fn sasa_colors(input: &ColorInput) -> Vec<Color32> {
    let n = input.loaded.structure.topology.residue_count();
    let Some(rel) = input.rel_sasa else {
        return vec![Color32::TRANSPARENT; n];
    };
    rel.iter()
        .map(|&v| {
            if v.is_nan() {
                return Color32::TRANSPARENT;
            }
            let buried = (1.0 - v).clamp(0.0, 1.0);
            let mix =
                |from: u8, to: u8| (from as f32 + (to as f32 - from as f32) * buried).round() as u8;
            Color32::from_rgb(mix(255, 0x2E), mix(255, 0x5E), mix(255, 0xAA))
        })
        .collect()
}

fn name_color(scheme: SeqColor, name: &str) -> Option<u32> {
    match scheme {
        SeqColor::Chemistry => vc::by_residue_type(name).or_else(|| vc::by_nucleotide(name)),
        SeqColor::Hydrophobicity => {
            vc::kyte_doolittle(name).map(|v| vc::by_hydropathy(v, vc::KD_SCALE))
        }
        SeqColor::Charge => charge_color(name),
        SeqColor::Clustal => vc::by_clustal(name),
        SeqColor::Zappo => vc::by_zappo(name),
        SeqColor::Taylor => vc::by_taylor(name),
        _ => None,
    }
}

fn by_name(scheme: SeqColor, top: &vv_core::Topology) -> Vec<Color32> {
    (0..top.residue_count())
        .map(|r| name_color(scheme, top.residue_name(r)).map_or(Color32::TRANSPARENT, from_packed))
        .collect()
}

/// Acidic red, basic blue, histidine pale blue (mostly neutral at pH 7.4).
fn charge_color(name: &str) -> Option<u32> {
    match name {
        "ASP" | "GLU" => Some(vc::rgba(0xE0, 0x55, 0x55)),
        "LYS" | "ARG" => Some(vc::rgba(0x55, 0x77, 0xE0)),
        "HIS" => Some(vc::rgba(0xB4, 0xC8, 0xF0)),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn scheme_words_and_aliases_parse() {
        for (scheme, word, _) in SCHEMES {
            assert_eq!(SeqColor::parse(word), Some(scheme));
        }
        assert_eq!(SeqColor::parse("3D"), Some(SeqColor::View));
        assert_eq!(SeqColor::parse("burial"), Some(SeqColor::Sasa));
        assert_eq!(SeqColor::parse("nope"), None);
    }

    const RESIDUES: [&str; 25] = [
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET",
        "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "MSE", "A", "C", "G", "U",
    ];

    fn chip_colors(scheme: SeqColor) -> Vec<Color32> {
        let mut out: Vec<Color32> = RESIDUES
            .iter()
            .filter_map(|n| name_color(scheme, n))
            .map(from_packed)
            .collect();
        if scheme == SeqColor::SecondaryStructure {
            use vv_core::dssp::DsspCode::*;
            out.extend(
                [AlphaHelix, Helix3_10, HelixPi, Strand, Bridge, Turn]
                    .map(|c| from_packed(vc::by_secondary_structure(c))),
            );
        }
        out
    }

    #[test]
    fn every_categorical_chip_is_readable_in_both_themes() {
        use crate::theme::{ThemeMode, Tokens};
        for mode in [ThemeMode::Light, ThemeMode::Dark] {
            let t = Tokens::of(mode);
            for (scheme, word, _) in SCHEMES.iter().filter(|s| s.0.categorical()) {
                for fill in chip_colors(*scheme) {
                    let chip = tint(fill, t.surface);
                    let text = readable_on(chip, t.text);
                    assert!(
                        contrast_ratio(text, chip) >= MIN_CONTRAST,
                        "{word} {fill:?} in {mode:?}"
                    );
                }
            }
        }
    }

    #[test]
    fn a_tint_is_softer_than_its_color_and_ramps_stay_readable() {
        let panel = Color32::from_rgb(0xF6, 0xF8, 0xFC);
        let soft = tint(Color32::from_rgb(255, 0, 0), panel);
        assert!(soft.g() > 100 && soft.r() > 200, "{soft:?}");
        for v in (0..=255u8).step_by(15) {
            let fill = Color32::from_rgb(v, 255 - v, 128);
            let text = readable_on(fill, Color32::from_gray(40));
            assert!(contrast_ratio(text, fill) >= MIN_CONTRAST, "{fill:?}");
        }
    }
}
