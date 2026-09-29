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

/// Text that reads on `fill`: near-black on light chips, white on dark.
pub fn text_on(fill: Color32) -> Color32 {
    let lum = 0.299 * fill.r() as f32 + 0.587 * fill.g() as f32 + 0.114 * fill.b() as f32;
    if lum > 150.0 {
        Color32::from_gray(25)
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

fn by_name(scheme: SeqColor, top: &vv_core::Topology) -> Vec<Color32> {
    (0..top.residue_count())
        .map(|r| {
            let name = top.residue_name(r);
            let packed = match scheme {
                SeqColor::Chemistry => {
                    vc::by_residue_type(name).or_else(|| vc::by_nucleotide(name))
                }
                SeqColor::Hydrophobicity => {
                    vc::kyte_doolittle(name).map(|v| vc::by_hydropathy(v, vc::KD_SCALE))
                }
                SeqColor::Charge => charge_color(name),
                SeqColor::Clustal => vc::by_clustal(name),
                SeqColor::Zappo => vc::by_zappo(name),
                SeqColor::Taylor => vc::by_taylor(name),
                _ => None,
            };
            packed.map_or(Color32::TRANSPARENT, from_packed)
        })
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

    #[test]
    fn text_contrasts_with_its_chip() {
        assert_eq!(text_on(Color32::WHITE), Color32::from_gray(25));
        assert_eq!(text_on(Color32::from_rgb(0x2E, 0x5E, 0xAA)), Color32::WHITE);
    }
}
