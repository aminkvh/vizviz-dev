//! Material presets: 23 classic named looks with fixed coefficients, plus
//! a handful of cheap extra looks. Plain data so the scene model can name
//! a material without a GPU dependency; `vv_render::Material::from` turns
//! one into what the shaders read. Old identifiers still `parse()` as
//! hidden aliases; `name()` only ever returns the current one.

/// A preset's Phong-style coefficients, plus two fields only the extra
/// looks use (`toon_bands`, `rim`).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialSettings {
    /// Light the surface has regardless of the lights.
    pub ambient: f32,
    pub specular: f32,
    pub diffuse: f32,
    /// 0..1; the specular exponent is `10^(3 * shininess)`.
    pub shininess: f32,
    pub opacity: f32,
    /// How much diffuse light fades toward the silhouette (can exceed 1).
    pub outline: f32,
    pub outline_width: f32,
    /// 1 = angle-dependent opacity (Merritt & Bacon 1997, Methods Enzymol.
    /// 277:505).
    pub transmode: f32,
    /// Diffuse quantized into this many bands; 0 is smooth.
    pub toon_bands: f32,
    /// Extra light added at grazing angles, brightening the silhouette
    /// instead of `outline` darkening it; 0 is off.
    pub rim: f32,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum MaterialPreset {
    #[default]
    Opaque,
    /// No shading: ambient 1, diffuse 0, specular 0, so the per-pixel
    /// shading term is exactly the base colour, unaffected by the key/
    /// fill lights or the lighting's own ambient term (it is additive,
    /// never multiplied by them). The viewport's own post effects --
    /// tonemap, ambient occlusion, shadows, depth cue -- still apply on
    /// top, same as any other material; `tonemap off; ao off; shadows
    /// off; depthcue off` gets the literal, untouched colour.
    Flat,
    Transparent,
    BrushedMetal,
    Diffuse,
    /// Faint and see-through: full specular at low opacity.
    Faint,
    ClearGlass,
    TintedGlass,
    FrostedGlass,
    Glossy,
    HardPlastic,
    /// A dim, pastel-tinted metal.
    SoftMetal,
    Steel,
    Translucent,
    /// Ink-outlined, matte.
    Inked,
    /// Ink-outlined, shiny.
    InkedGloss,
    /// Ink-outlined, transparent.
    InkedGlass,
    /// Flat, high-contrast illustrative look, after David Goodsell's
    /// textbook molecular illustrations.
    Textbook,
    /// Shiny with strong ambient occlusion.
    Polished,
    /// Matte with strong ambient occlusion.
    Chalk,
    /// Ink-outlined, with strong ambient occlusion.
    ChalkEdge,
    /// A blown-glass bubble: thin, angle-dependent opacity.
    BubbledGlass,
    /// Thinner and more transparent than `BubbledGlass`.
    HollowGlass,
    /// A mirror-bright chrome look (drawn as `Opaque` here; its mirror
    /// reflection is ray-tracing only).
    MirrorChrome,
    /// Three-band cel shading.
    Toon,
    /// Flat, shadow-soft studio look: high ambient, gentle diffuse,
    /// almost no specular -- a cheap stand-in for a matcap (lit-sphere)
    /// texture.
    Clay,
    /// A Fresnel-edge outline, brightened instead of darkened: a thin
    /// glow at grazing angles.
    Rim,
    /// A matte studio look: low ambient (0.14), soft diffuse and
    /// specular (0.45 each), specular exponent 55.0.
    StudioMatte,
    /// A flat, broadly-lit default look: no ambient term, diffuse 0.8,
    /// specular 0.3, specular exponent 30.0.
    LabDefault,
}

/// The classic named materials: ambient, specular, diffuse, shininess
/// (0..1), opacity, outline, outline width, transmode. A mirror term
/// (ray tracing only) is left out, so MirrorChrome draws as Opaque here.
#[rustfmt::skip]
const BASE: [(MaterialPreset, &str, [f32; 8]); 23] = {
    use MaterialPreset as M;
    [
        (M::Opaque,       "opaque",       [0.0,  0.5,  0.65, 0.534, 1.0,  0.0,  0.0,  0.0]),
        (M::Transparent,  "transparent",  [0.0,  0.5,  0.65, 0.534, 0.3,  0.0,  0.0,  0.0]),
        (M::BrushedMetal, "brushedmetal", [0.08, 0.34, 0.39, 0.15,  1.0,  0.0,  0.0,  0.0]),
        (M::Diffuse,      "diffuse",      [0.0,  0.0,  0.62, 0.53,  1.0,  0.0,  0.0,  0.0]),
        (M::Faint,        "faint",        [0.0,  1.0,  0.0,  0.23,  0.1,  0.0,  0.0,  0.0]),
        (M::ClearGlass,   "clearglass",   [0.0,  0.65, 0.5,  0.53,  0.15, 0.0,  0.0,  0.0]),
        (M::TintedGlass,  "tintedglass",  [0.52, 0.22, 0.76, 0.59,  0.68, 0.0,  0.0,  0.0]),
        (M::FrostedGlass, "frostedglass", [0.15, 0.75, 0.25, 0.8,   0.5,  0.0,  0.0,  0.0]),
        (M::Glossy,       "glossy",       [0.0,  1.0,  0.65, 0.88,  1.0,  0.0,  0.0,  0.0]),
        (M::HardPlastic,  "hardplastic",  [0.0,  0.28, 0.56, 0.69,  1.0,  0.0,  0.0,  0.0]),
        (M::SoftMetal,    "softmetal",    [0.0,  0.55, 0.26, 0.19,  1.0,  0.0,  0.0,  0.0]),
        (M::Steel,        "steel",        [0.25, 0.38, 0.0,  0.32,  1.0,  0.0,  0.0,  0.0]),
        (M::Translucent,  "translucent",  [0.0,  0.6,  0.7,  0.3,   0.8,  0.0,  0.0,  0.0]),
        (M::Inked,        "inked",        [0.0,  0.0,  0.66, 0.75,  1.0,  0.62, 0.94, 0.0]),
        (M::InkedGloss,   "inkedgloss",   [0.0,  0.96, 0.66, 0.75,  1.0,  0.76, 0.94, 0.0]),
        (M::InkedGlass,   "inkedglass",   [0.0,  0.5,  0.66, 0.75,  0.62, 0.62, 0.94, 0.0]),
        (M::Textbook,     "textbook",     [0.52, 0.0,  1.0,  0.0,   1.0,  4.0,  0.9,  0.0]),
        (M::Polished,     "polished",     [0.0,  0.2,  0.85, 0.53,  1.0,  0.0,  0.0,  0.0]),
        (M::Chalk,        "chalk",        [0.0,  0.0,  0.85, 0.53,  1.0,  0.0,  0.0,  0.0]),
        (M::ChalkEdge,    "chalkedge",    [0.0,  0.2,  0.9,  0.53,  1.0,  0.62, 0.93, 0.0]),
        (M::BubbledGlass, "bubbledglass", [0.04, 1.0,  0.34, 1.0,   0.1,  0.0,  0.0,  1.0]),
        (M::HollowGlass,  "hollowglass",  [0.25, 1.0,  0.34, 1.0,   0.04, 0.0,  0.0,  1.0]),
        (M::MirrorChrome, "mirrorchrome", [0.0,  0.5,  0.65, 0.53,  1.0,  0.0,  0.0,  0.0]),
    ]
};

/// The extra looks: ambient, specular, diffuse, shininess, opacity,
/// outline, outline width, transmode, toon bands, rim -- same order as
/// `BASE`'s columns, with the two extra fields appended.
#[rustfmt::skip]
const EXTRAS: [(MaterialPreset, &str, [f32; 10]); 6] = {
    use MaterialPreset as M;
    [
        (M::Flat,        "flat",        [1.0,  0.0,  0.0,  0.0,  1.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
        (M::Toon,        "toon",        [0.0,  0.0,  0.9,  0.0,  1.0, 0.0, 0.0, 0.0, 3.0, 0.0]),
        (M::Clay,        "clay",        [0.35, 0.05, 0.65, 0.2,  1.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
        (M::Rim,         "rim",         [0.0,  0.3,  0.6,  0.5,  1.0, 0.0, 0.0, 0.0, 0.0, 1.0]),
        (M::StudioMatte, "studiomatte", [0.14, 0.45, 0.45, 0.58, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
        (M::LabDefault,  "labdefault",  [0.0,  0.3,  0.8,  0.49, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
    ]
};

/// Old identifiers that `parse()` still accepts so old scripts and saved
/// sessions keep working; never returned by `name()`.
#[rustfmt::skip]
const ALIASES: [(&str, MaterialPreset); 15] = {
    use MaterialPreset as M;
    [
        ("ghost",           M::Faint),
        ("glass1",          M::ClearGlass),
        ("glass2",          M::TintedGlass),
        ("glass3",          M::FrostedGlass),
        ("metallicpastel",  M::SoftMetal),
        ("edgy",            M::Inked),
        ("edgyshiny",       M::InkedGloss),
        ("edgyglass",       M::InkedGlass),
        ("goodsell",        M::Textbook),
        ("aoshiny",         M::Polished),
        ("aochalky",        M::Chalk),
        ("aoedgy",          M::ChalkEdge),
        ("blownglass",      M::BubbledGlass),
        ("glassbubble",     M::HollowGlass),
        ("rtchrome",        M::MirrorChrome),
    ]
};

impl MaterialPreset {
    pub const ALL: [MaterialPreset; 29] = {
        let mut all = [MaterialPreset::Toon; 29];
        let mut i = 0;
        while i < BASE.len() {
            all[i] = BASE[i].0;
            i += 1;
        }
        let mut j = 0;
        while j < EXTRAS.len() {
            all[BASE.len() + j] = EXTRAS[j].0;
            j += 1;
        }
        all
    };

    pub fn name(self) -> &'static str {
        if let Some((_, name, _)) = BASE.iter().find(|(m, ..)| *m == self) {
            return name;
        }
        match EXTRAS.iter().find(|(m, ..)| *m == self) {
            Some((_, name, _)) => name,
            None => "toon",
        }
    }

    pub fn parse(name: &str) -> Option<MaterialPreset> {
        MaterialPreset::ALL
            .into_iter()
            .find(|p| p.name() == name)
            .or_else(|| ALIASES.iter().find(|(n, _)| *n == name).map(|(_, p)| *p))
    }

    /// This preset's settings, from `BASE` or `EXTRAS`.
    pub fn settings(self) -> MaterialSettings {
        if let Some((.., v)) = BASE.iter().find(|(m, ..)| *m == self) {
            let [ambient, specular, diffuse, shininess, opacity, outline, outline_width, transmode] =
                *v;
            return MaterialSettings {
                ambient,
                specular,
                diffuse,
                shininess,
                opacity,
                outline,
                outline_width,
                transmode,
                toon_bands: 0.0,
                rim: 0.0,
            };
        }
        let v = match EXTRAS.iter().find(|(m, ..)| *m == self) {
            Some((.., v)) => *v,
            None => EXTRAS[0].2,
        };
        let [ambient, specular, diffuse, shininess, opacity, outline, outline_width, transmode, toon_bands, rim] =
            v;
        MaterialSettings {
            ambient,
            specular,
            diffuse,
            shininess,
            opacity,
            outline,
            outline_width,
            transmode,
            toon_bands,
            rim,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn names_round_trip_and_are_unique() {
        for (i, a) in MaterialPreset::ALL.iter().enumerate() {
            assert_eq!(MaterialPreset::parse(a.name()), Some(*a));
            for b in &MaterialPreset::ALL[i + 1..] {
                assert_ne!(a.name(), b.name());
            }
        }
    }

    #[test]
    fn old_names_still_parse_and_never_shadow_a_current_name() {
        for (alias, preset) in ALIASES {
            assert_eq!(
                MaterialPreset::parse(alias),
                Some(preset),
                "alias {alias:?}"
            );
            assert!(
                MaterialPreset::ALL.iter().all(|p| p.name() != alias),
                "{alias:?} is also a current name"
            );
        }
    }
}
