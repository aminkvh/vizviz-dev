//! Color schemes producing packed RGBA8 (little-endian `0xAABBGGRR`).

use vv_core::Element;

pub const fn rgba(r: u8, g: u8, b: u8) -> u32 {
    (r as u32) | ((g as u32) << 8) | ((b as u32) << 16) | (0xFF << 24)
}

/// Which per-atom attribute a structure's colors are derived from. Mirrors
/// `vv_scene::ColorScheme`; kept separate so this crate stays the only one
/// that knows what an RGBA8 color looks like.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum ColorScheme {
    /// Conventional CPK-style element colors.
    #[default]
    Element,
    /// One well-separated hue per chain, independent of chain count.
    Chain,
    /// Blue (low) -> white (mid) -> red (high) over the structure's own
    /// b-factor range.
    BFactor,
    /// Secondary structure, one colour per class. `colors_for` uses
    /// the file's HELIX/SHEET records; `colors_for_ss` takes an
    /// assignment (e.g. DSSP for a trajectory frame).
    SecondaryStructure,
    /// Residue type: acidic red, basic blue, polar green, nonpolar white.
    ResidueType,
    /// Blue (N-terminus) to red (C-terminus) along each chain.
    Rainbow,
    /// Carbons by chain, every other element by element.
    Hetero,
    /// One well-separated hue per residue name.
    ResidueName,
    /// The shared scalar ramp over the structure's occupancy range.
    Occupancy,
    /// Kyte-Doolittle hydropathy of amino acids (Kyte & Doolittle 1982,
    /// J. Mol. Biol. 157:105), teal (hydrophilic) to white to goldenrod
    /// (hydrophobic); other atoms by element.
    Hydrophobicity,
    /// Wimley-White whole-residue octanol hydropathy; same ramp as
    /// `Hydrophobicity`.
    WimleyWhite,
    /// One hue per `auth_asym_id` (see
    /// `vv_scene::ColorScheme::SegmentName`).
    SegmentName,
    /// Zappo physicochemical grouping (Livingstone & Barton 1993).
    Zappo,
    /// Taylor (1997) physicochemical colour wheel.
    Taylor,
    /// Simplified (non-alignment) Clustal colouring (Thompson et al. 1997).
    Clustal,
    /// Chou-Fasman alpha-helix propensity.
    HelixPropensity,
    /// Chou-Fasman beta-strand propensity.
    StrandPropensity,
    /// Chou-Fasman turn propensity.
    TurnPropensity,
    /// Fraction of residues of this type typically buried.
    BuriedIndex,
    /// Nucleotide identity (A/C/G/T/U).
    Nucleotide,
    /// Purine (A, G) vs pyrimidine (C, T, U).
    PurinePyrimidine,
    /// One fixed colour per residue molecule class.
    Class,
    /// Every atom one packed colour.
    Constant(u32),
}

/// The fixed display colour of a molecule class: a display convention,
/// picked to keep the eight classes distinguishable, not data.
pub fn by_class(class: vv_core::ResidueClass) -> u32 {
    use vv_core::ResidueClass as C;
    match class {
        C::Protein => rgba(0x6F, 0x8F, 0xD0),
        C::Nucleic => rgba(0xE8, 0xA3, 0x3D),
        C::Lipid => rgba(0xF2, 0xD1, 0x6B),
        C::Glycan => rgba(0x4F, 0xB8, 0xA8),
        C::Water => rgba(0x9A, 0xD3, 0xF0),
        C::Ion => rgba(0x6F, 0xCF, 0x6F),
        C::SmallMolecule => rgba(0xD4, 0x55, 0xA8),
        C::Other => rgba(0x9A, 0x9A, 0x9A),
    }
}

/// Kyte & Doolittle (1982) hydropathy of a standard amino acid.
pub fn kyte_doolittle(residue: &str) -> Option<f32> {
    Some(match residue {
        "ILE" => 4.5,
        "VAL" => 4.2,
        "LEU" => 3.8,
        "PHE" => 2.8,
        "CYS" => 2.5,
        "MET" | "MSE" => 1.9,
        "ALA" => 1.8,
        "GLY" => -0.4,
        "THR" => -0.7,
        "SER" => -0.8,
        "TRP" => -0.9,
        "TYR" => -1.3,
        "PRO" => -1.6,
        "HIS" => -3.2,
        "GLU" | "GLN" | "ASP" | "ASN" => -3.5,
        "LYS" => -3.9,
        "ARG" => -4.5,
        _ => return None,
    })
}

/// Wimley & White (1996) Nat. Struct. Biol. 3:842-848, Table 1: whole-
/// residue water -> n-octanol transfer free energy (kcal/mol). Positive is
/// unfavorable to octanol (hydrophilic) -- the opposite sign convention to
/// `kyte_doolittle`, so callers feeding `by_hydropathy` negate it.
/// Charged residues use their physiological form (Asp-, Glu-, Lys+, Arg+);
/// His uses its mostly-neutral form at pH 7, as `kyte_doolittle` does.
pub fn wimley_white(residue: &str) -> Option<f32> {
    Some(match residue {
        "TRP" => -2.09,
        "PHE" => -1.71,
        "LEU" => -1.25,
        "ILE" => -1.12,
        "MET" | "MSE" => -0.67,
        "TYR" => -0.71,
        "VAL" => -0.46,
        "CYS" => -0.02,
        "PRO" => 0.14,
        "THR" => 0.25,
        "SER" => 0.46,
        "HIS" => 0.11,
        "GLN" => 0.77,
        "ALA" => 0.50,
        "GLY" => 1.15,
        "ASN" => 0.85,
        "ARG" => 1.81,
        "LYS" => 2.80,
        "GLU" => 3.63,
        "ASP" => 3.64,
        _ => return None,
    })
}

/// Kyte & Doolittle's own most extreme magnitude (Arg -4.5, Ile +4.5):
/// `by_hydropathy`'s normalizing half-range for that table.
pub const KD_SCALE: f32 = 4.5;
/// Wimley-White's own most extreme magnitude (Asp- +3.64): `by_hydropathy`'s
/// normalizing half-range for that table.
pub const WW_SCALE: f32 = 3.64;

/// Teal (most hydrophilic) -> white (0) -> goldenrod (most hydrophobic),
/// `value` normalized by `scale` (the table's own most extreme magnitude:
/// `KD_SCALE` or `WW_SCALE`).
pub fn by_hydropathy(value: f32, scale: f32) -> u32 {
    let t = (value / scale).clamp(-1.0, 1.0);
    let (end, k) = if t < 0.0 {
        ((0x00, 0x8B, 0x8B), -t)
    } else {
        ((0xDA, 0xA5, 0x20), t)
    };
    let mix = |a: u8, b: u8| (a as f32 + (b as f32 - a as f32) * k).round() as u8;
    rgba(mix(0xFF, end.0), mix(0xFF, end.1), mix(0xFF, end.2))
}

/// A hue per residue name, stable across structures (from the name, not
/// the order it appears in).
pub fn by_residue_name(name: &str) -> u32 {
    let hash = name
        .bytes()
        .fold(2166136261u32, |h, b| (h ^ b as u32).wrapping_mul(16777619));
    by_chain(hash % 64)
}

/// Secondary-structure colours: purple, blue, red, yellow, tan, cyan,
/// white.
pub fn by_secondary_structure(code: vv_core::dssp::DsspCode) -> u32 {
    use vv_core::dssp::DsspCode::*;
    match code {
        AlphaHelix => rgba(0xA6, 0x00, 0xA6),
        Helix3_10 => rgba(0x00, 0x00, 0xFF),
        HelixPi => rgba(0xFF, 0x00, 0x00),
        Strand => rgba(0xFF, 0xFF, 0x00),
        Bridge => rgba(0x80, 0x80, 0x33),
        Turn => rgba(0x40, 0xBF, 0xBF),
        Bend | Coil | None => rgba(0xFF, 0xFF, 0xFF),
    }
}

/// Residue-type colours; `None` for anything that isn't a standard amino
/// acid (the caller falls back to element colors).
pub fn by_residue_type(name: &str) -> Option<u32> {
    match name {
        "ASP" | "GLU" => Some(rgba(0xFF, 0x00, 0x00)),
        "ARG" | "LYS" | "HIS" | "HID" | "HIE" | "HIP" => Some(rgba(0x00, 0x00, 0xFF)),
        "SER" | "THR" | "ASN" | "GLN" | "TYR" | "CYS" => Some(rgba(0x00, 0xCC, 0x00)),
        "ALA" | "GLY" | "VAL" | "LEU" | "ILE" | "MET" | "PHE" | "TRP" | "PRO" => {
            Some(rgba(0xFF, 0xFF, 0xFF))
        }
        _ => None,
    }
}

/// Zappo scheme (Livingstone & Barton 1993 physicochemical grouping): hydrophobic pink, aromatic orange,
/// negative red, positive blue, hydrophilic green, Gly/Pro magenta.
pub fn by_zappo(name: &str) -> Option<u32> {
    const PINK: u32 = rgba(0xFF, 0xAF, 0xAF);
    const MID_BLUE: u32 = rgba(0x64, 0x64, 0xFF);
    const GREEN: u32 = rgba(0x00, 0xFF, 0x00);
    const RED: u32 = rgba(0xFF, 0x00, 0x00);
    const YELLOW: u32 = rgba(0xFF, 0xFF, 0x00);
    const MAGENTA: u32 = rgba(0xFF, 0x00, 0xFF);
    const ORANGE: u32 = rgba(0xFF, 0xC8, 0x00);
    Some(match name {
        "ALA" | "ILE" | "LEU" | "MET" | "MSE" | "VAL" => PINK,
        "ARG" | "LYS" | "HIS" => MID_BLUE,
        "ASN" | "GLN" | "SER" | "THR" => GREEN,
        "ASP" | "GLU" => RED,
        "CYS" => YELLOW,
        "GLY" | "PRO" => MAGENTA,
        "PHE" | "TRP" | "TYR" => ORANGE,
        _ => return None,
    })
}

/// Taylor (1997) Protein Eng. 10(7):743-746, "Residual colours: a proposal
/// for aminochromography".
pub fn by_taylor(name: &str) -> Option<u32> {
    Some(match name {
        "ALA" => rgba(0xCC, 0xFF, 0x00),
        "ARG" => rgba(0x00, 0x00, 0xFF),
        "ASN" => rgba(0xCC, 0x00, 0xFF),
        "ASP" => rgba(0xFF, 0x00, 0x00),
        "CYS" => rgba(0xFF, 0xFF, 0x00),
        "GLN" => rgba(0xFF, 0x00, 0xCC),
        "GLU" => rgba(0xFF, 0x00, 0x66),
        "GLY" => rgba(0xFF, 0x99, 0x00),
        "HIS" => rgba(0x00, 0x66, 0xFF),
        "ILE" => rgba(0x66, 0xFF, 0x00),
        "LEU" => rgba(0x33, 0xFF, 0x00),
        "LYS" => rgba(0x66, 0x00, 0xFF),
        "MET" | "MSE" => rgba(0x00, 0xFF, 0x00),
        "PHE" => rgba(0x00, 0xFF, 0x66),
        "PRO" => rgba(0xFF, 0xCC, 0x00),
        "SER" => rgba(0xFF, 0x33, 0x00),
        "THR" => rgba(0xFF, 0x66, 0x00),
        "TRP" => rgba(0x00, 0xCC, 0xFF),
        "TYR" => rgba(0x00, 0xFF, 0xCC),
        "VAL" => rgba(0x99, 0xFF, 0x00),
        _ => return None,
    })
}

/// ClustalX's default per-residue-type colours (Thompson et al. 1997,
/// Nucleic Acids Res. 25:4876), simplified: the original also raises or drops a residue's colour by how
/// conserved its alignment column is (e.g. Cys only turns pink at >=85%
/// column identity); this colours each residue by its own type alone, as
/// if every column were 100% that residue -- so every threshold the real
/// scheme has for that type, including Cys's, applies unconditionally.
pub fn by_clustal(name: &str) -> Option<u32> {
    const BLUE: u32 = rgba(0x80, 0xB3, 0xE6);
    const RED: u32 = rgba(0xE6, 0x33, 0x1A);
    const GREEN: u32 = rgba(0x1A, 0xCC, 0x1A);
    const PINK: u32 = rgba(0xE6, 0x80, 0x80);
    const MAGENTA: u32 = rgba(0xCC, 0x4D, 0xCC);
    const ORANGE: u32 = rgba(0xE6, 0x99, 0x4D);
    const CYAN: u32 = rgba(0x1A, 0xB3, 0xB3);
    const YELLOW: u32 = rgba(0xCC, 0xCC, 0x00);
    Some(match name {
        "ALA" | "ILE" | "LEU" | "MET" | "MSE" | "PHE" | "TRP" | "VAL" => BLUE,
        "ARG" | "LYS" => RED,
        "ASN" | "GLN" | "SER" | "THR" => GREEN,
        "ASP" | "GLU" => MAGENTA,
        "CYS" => PINK,
        "GLY" => ORANGE,
        "HIS" | "TYR" => CYAN,
        "PRO" => YELLOW,
        _ => return None,
    })
}

/// Chou & Fasman (1978) Adv. Enzymol. 47:45-148 alpha-helix propensity
/// (values, and the 0.57..1.51 range `by_helix_propensity` normalizes by).
pub fn helix_propensity(name: &str) -> Option<f32> {
    Some(match name {
        "GLU" => 1.51,
        "MET" | "MSE" => 1.45,
        "ALA" => 1.42,
        "LEU" => 1.21,
        "LYS" => 1.16,
        "PHE" => 1.13,
        "GLN" => 1.11,
        "TRP" => 1.08,
        "ILE" => 1.08,
        "VAL" => 1.06,
        "ASP" => 1.01,
        "HIS" => 1.00,
        "ARG" => 0.98,
        "THR" => 0.83,
        "SER" => 0.77,
        "CYS" => 0.70,
        "TYR" => 0.69,
        "ASN" => 0.67,
        "PRO" => 0.57,
        "GLY" => 0.57,
        _ => return None,
    })
}

/// Green (least helix-forming) -> magenta (most).
pub fn by_helix_propensity(value: f32) -> u32 {
    let t = ((value - 0.57) / (1.51 - 0.57)).clamp(0.0, 1.0);
    let c = (t * 255.0).round() as u8;
    rgba(c, 255 - c, c)
}

/// Chou & Fasman (1978) beta-strand propensity (values and the
/// 0.37..1.70 range `by_strand_propensity` normalizes by).
pub fn strand_propensity(name: &str) -> Option<f32> {
    Some(match name {
        "VAL" => 1.70,
        "ILE" => 1.60,
        "TYR" => 1.47,
        "CYS" => 1.19,
        "THR" => 1.19,
        "TRP" => 1.37,
        "PHE" => 1.38,
        "LEU" => 1.30,
        "MET" | "MSE" => 1.05,
        "GLN" => 1.10,
        "ARG" => 0.93,
        "ASN" => 0.89,
        "HIS" => 0.87,
        "ALA" => 0.83,
        "SER" => 0.75,
        "GLY" => 0.75,
        "LYS" => 0.74,
        "PRO" => 0.55,
        "ASP" => 0.54,
        "GLU" => 0.37,
        _ => return None,
    })
}

/// Blue (least strand-forming) -> yellow (most).
pub fn by_strand_propensity(value: f32) -> u32 {
    let t = ((value - 0.37) / (1.70 - 0.37)).clamp(0.0, 1.0);
    let c = (t * 255.0).round() as u8;
    rgba(c, c, 255 - c)
}

/// Chou & Fasman (1978) turn propensity (values and the 0.47..1.56 range
/// `by_turn_propensity` normalizes by).
pub fn turn_propensity(name: &str) -> Option<f32> {
    Some(match name {
        "ASN" => 1.56,
        "GLY" => 1.56,
        "PRO" => 1.52,
        "ASP" => 1.46,
        "SER" => 1.43,
        "CYS" => 1.19,
        "TYR" => 1.14,
        "LYS" => 1.01,
        "GLN" => 0.98,
        "THR" => 0.96,
        "TRP" => 0.96,
        "ARG" => 0.95,
        "HIS" => 0.95,
        "GLU" => 0.74,
        "ALA" => 0.66,
        "MET" | "MSE" => 0.60,
        "PHE" => 0.60,
        "LEU" => 0.59,
        "VAL" => 0.50,
        "ILE" => 0.47,
        _ => return None,
    })
}

/// Cyan (least turn-forming) -> red (most).
pub fn by_turn_propensity(value: f32) -> u32 {
    let t = ((value - 0.47) / (1.56 - 0.47)).clamp(0.0, 1.0);
    let c = (t * 255.0).round() as u8;
    rgba(c, 255 - c, 255 - c)
}

/// Fraction of residues of this type found buried (values and the
/// 0.05..4.6 range `by_buried_index` normalizes by).
pub fn buried_index(name: &str) -> Option<f32> {
    Some(match name {
        "CYS" => 4.6,
        "ILE" => 3.1,
        "VAL" => 2.9,
        "LEU" => 2.4,
        "PHE" => 2.2,
        "MET" | "MSE" => 1.9,
        "GLY" => 1.8,
        "ALA" => 1.7,
        "TRP" => 1.6,
        "HIS" => 0.8,
        "SER" => 0.8,
        "THR" => 0.7,
        "TYR" => 0.5,
        "ASN" => 0.4,
        "ASP" => 0.4,
        "GLN" => 0.3,
        "GLU" => 0.3,
        "ARG" => 0.1,
        "LYS" => 0.05,
        _ => return None,
    })
}

/// Green (least buried) -> blue (most).
pub fn by_buried_index(value: f32) -> u32 {
    let t = ((value - 0.05) / (4.6 - 0.05)).clamp(0.0, 1.0);
    let c = (t * 255.0).round() as u8;
    rgba(0, 255 - c, c)
}

/// A DNA residue name's base letter ("DA" -> "A"), so the nucleotide
/// tables need one entry per base regardless of DNA/RNA naming.
fn base_of(name: &str) -> &str {
    match name {
        "DA" | "DC" | "DG" | "DT" | "DI" => &name[1..],
        _ => name,
    }
}

/// Nucleotide colours: A green, C orange, G red, T/U
/// blue.
pub fn by_nucleotide(name: &str) -> Option<u32> {
    Some(match base_of(name) {
        "A" => rgba(0x64, 0xF7, 0x3F),
        "C" => rgba(0xFF, 0xB3, 0x40),
        "G" => rgba(0xEB, 0x41, 0x3C),
        "T" | "U" => rgba(0x3C, 0x88, 0xEE),
        _ => return None,
    })
}

/// Purines (A, G) orchid,
/// pyrimidines (C, T, U) turquoise.
pub fn by_purine_pyrimidine(name: &str) -> Option<u32> {
    Some(match base_of(name) {
        "A" | "G" => rgba(0xFF, 0x83, 0xFA),
        "C" | "T" | "U" => rgba(0x40, 0xE0, 0xD0),
        _ => return None,
    })
}

/// Blue at `t = 0` through green to red at `t = 1`.
pub fn rainbow(t: f32) -> u32 {
    hsv_to_rgb(240.0 * (1.0 - t.clamp(0.0, 1.0)), 0.85, 0.95)
}

/// Conventional CPK-style element colors (Corey & Pauling 1953, Rev. Sci.
/// Instrum. 24:621); carbon is lighter than the usual 0x909090, which
/// reads near black once shaded. Display convention, not data.
pub fn by_element(e: Element) -> u32 {
    match e.atomic_number() {
        1 => rgba(0xFF, 0xFF, 0xFF),
        6 => rgba(0xB4, 0xB4, 0xB4),
        7 => rgba(0x30, 0x50, 0xF8),
        8 => rgba(0xFF, 0x0D, 0x0D),
        9 => rgba(0x90, 0xE0, 0x50),
        11 => rgba(0xAB, 0x5C, 0xF2),
        12 => rgba(0x8A, 0xFF, 0x00),
        15 => rgba(0xFF, 0x80, 0x00),
        16 => rgba(0xFF, 0xFF, 0x30),
        17 => rgba(0x1F, 0xF0, 0x1F),
        19 => rgba(0x8F, 0x40, 0xD4),
        20 => rgba(0x3D, 0xFF, 0x00),
        26 => rgba(0xE0, 0x66, 0x33),
        30 => rgba(0x7D, 0x80, 0xB0),
        35 => rgba(0xA6, 0x29, 0x29),
        53 => rgba(0x94, 0x00, 0x94),
        _ => rgba(0xFF, 0x14, 0x93),
    }
}

fn lerp_u8(a: u8, b: u8, t: f32) -> u8 {
    (a as f32 + (b as f32 - a as f32) * t).round() as u8
}

fn hsv_to_rgb(h: f32, s: f32, v: f32) -> u32 {
    let c = v * s;
    let hp = h / 60.0;
    let x = c * (1.0 - (hp % 2.0 - 1.0).abs());
    let (r1, g1, b1) = match hp as u32 {
        0 => (c, x, 0.0),
        1 => (x, c, 0.0),
        2 => (0.0, c, x),
        3 => (0.0, x, c),
        4 => (x, 0.0, c),
        _ => (c, 0.0, x),
    };
    let m = v - c;
    rgba(
        ((r1 + m) * 255.0).round() as u8,
        ((g1 + m) * 255.0).round() as u8,
        ((b1 + m) * 255.0).round() as u8,
    )
}

/// A well-separated color per chain index, spacing hues by the golden angle
/// so it looks reasonable whether a structure has 2 chains or 60 — no need
/// to know the total chain count up front.
pub fn by_chain(chain_index: u32) -> u32 {
    let hue = (chain_index as f32 * 137.508) % 360.0;
    hsv_to_rgb(hue, 0.65, 0.95)
}

const B_FACTOR_LOW: (u8, u8, u8) = (0x21, 0x4C, 0xE0);
const B_FACTOR_MID: (u8, u8, u8) = (0xFF, 0xFF, 0xFF);
const B_FACTOR_HIGH: (u8, u8, u8) = (0xE0, 0x30, 0x30);

/// Blue (`min`) -> white (midpoint) -> red (`max`). A structure where every
/// atom shares one b-factor (or has none recorded, i.e. `min == max`) comes
/// back white rather than dividing by zero.
pub fn by_b_factor(value: f32, min: f32, max: f32) -> u32 {
    if max <= min {
        return rgba(B_FACTOR_MID.0, B_FACTOR_MID.1, B_FACTOR_MID.2);
    }
    let t = ((value - min) / (max - min)).clamp(0.0, 1.0);
    let (from, to, k) = if t < 0.5 {
        (B_FACTOR_LOW, B_FACTOR_MID, t * 2.0)
    } else {
        (B_FACTOR_MID, B_FACTOR_HIGH, (t - 0.5) * 2.0)
    };
    rgba(
        lerp_u8(from.0, to.0, k),
        lerp_u8(from.1, to.1, k),
        lerp_u8(from.2, to.2, k),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn packs_little_endian_rgba() {
        assert_eq!(rgba(0x11, 0x22, 0x33), 0xFF33_2211);
        assert_eq!(by_element(Element::CARBON), 0xFFB4_B4B4);
    }

    #[test]
    fn chain_colors_are_well_separated_and_deterministic() {
        assert_eq!(by_chain(0), by_chain(0));
        let colors: Vec<u32> = (0..8).map(by_chain).collect();
        for i in 0..colors.len() {
            for j in (i + 1)..colors.len() {
                assert_ne!(colors[i], colors[j], "chains {i} and {j} collided");
            }
        }
    }

    #[test]
    fn b_factor_ramps_blue_to_red_through_white() {
        assert_eq!(by_b_factor(0.0, 0.0, 100.0), rgba(0x21, 0x4C, 0xE0));
        assert_eq!(by_b_factor(50.0, 0.0, 100.0), rgba(0xFF, 0xFF, 0xFF));
        assert_eq!(by_b_factor(100.0, 0.0, 100.0), rgba(0xE0, 0x30, 0x30));
        // Out-of-range values clamp instead of extrapolating.
        assert_eq!(by_b_factor(-10.0, 0.0, 100.0), by_b_factor(0.0, 0.0, 100.0));
        assert_eq!(
            by_b_factor(200.0, 0.0, 100.0),
            by_b_factor(100.0, 0.0, 100.0)
        );
    }

    #[test]
    fn b_factor_with_no_spread_is_white() {
        assert_eq!(by_b_factor(5.0, 5.0, 5.0), rgba(0xFF, 0xFF, 0xFF));
    }

    #[test]
    fn kd_and_ww_hydropathy_use_their_own_scale_and_agree_on_direction() {
        assert_eq!(kyte_doolittle("ILE"), Some(4.5));
        assert_eq!(wimley_white("TRP"), Some(-2.09));
        assert_eq!(wimley_white("ASP"), Some(3.64));
        assert_eq!(wimley_white("XXX"), None);
        // Each table's own most extreme magnitude saturates the ramp.
        assert_eq!(
            by_hydropathy(KD_SCALE, KD_SCALE),
            by_hydropathy(100.0, KD_SCALE)
        );
        assert_eq!(
            by_hydropathy(-WW_SCALE, WW_SCALE),
            by_hydropathy(-100.0, WW_SCALE)
        );
        // Wimley-White's sign is opposite Kyte-Doolittle's (positive =
        // hydrophilic, not hydrophobic): callers negate it, which must
        // actually flip which end of the ramp a value lands on.
        let trp = wimley_white("TRP").unwrap();
        assert_ne!(
            by_hydropathy(trp, WW_SCALE),
            by_hydropathy(-trp, WW_SCALE),
            "negating a nonzero value must change its colour"
        );
    }

    #[test]
    fn zappo_taylor_and_clustal_group_known_residues() {
        assert_eq!(by_zappo("ASP"), by_zappo("GLU"), "both acidic, both red");
        assert_ne!(by_zappo("ASP"), by_zappo("LYS"));
        assert_eq!(by_taylor("ALA"), Some(rgba(0xCC, 0xFF, 0x00)));
        assert_eq!(by_taylor("XXX"), None);
        // Clustal's Cys-specific override applies unconditionally under
        // the simplified, non-alignment reduction.
        assert_eq!(by_clustal("CYS"), Some(rgba(0xE6, 0x80, 0x80)));
        assert_eq!(by_clustal("ASP"), by_clustal("GLU"));
    }

    #[test]
    fn propensity_scales_span_their_own_min_to_max() {
        assert_eq!(helix_propensity("GLU"), Some(1.51));
        assert_eq!(by_helix_propensity(1.51), by_helix_propensity(100.0));
        assert_eq!(by_helix_propensity(0.57), by_helix_propensity(-100.0));

        assert_eq!(strand_propensity("VAL"), Some(1.70));
        assert_eq!(by_strand_propensity(1.70), by_strand_propensity(100.0));

        assert_eq!(turn_propensity("ASN"), Some(1.56));
        assert_eq!(by_turn_propensity(1.56), by_turn_propensity(100.0));

        assert_eq!(buried_index("CYS"), Some(4.6));
        assert_eq!(by_buried_index(4.6), by_buried_index(100.0));
    }

    #[test]
    fn nucleotide_schemes_cover_dna_and_rna_spellings() {
        assert_eq!(by_nucleotide("A"), by_nucleotide("DA"));
        assert_ne!(by_nucleotide("A"), by_nucleotide("G"));
        assert_eq!(by_purine_pyrimidine("A"), by_purine_pyrimidine("G"));
        assert_eq!(by_purine_pyrimidine("DT"), by_purine_pyrimidine("U"));
        assert_ne!(by_purine_pyrimidine("A"), by_purine_pyrimidine("C"));
        assert_eq!(by_nucleotide("XXX"), None);
    }
}
