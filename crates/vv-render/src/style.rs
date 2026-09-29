//! How things are lit and what they are made of, kept apart:
//!
//! - [`Lighting`] is global (one per view): key and fill lights, ambient,
//!   tonemapping. Presets: [`LightingPreset`].
//! - [`Material`] belongs to each drawn object: how it reflects that
//!   light. Presets: [`MaterialPreset`].
//! - [`StylePreset`] is a one-click look: a lighting preset, a material,
//!   a background and post-effect defaults.
//!
//! Shading is a Phong-style model (`shaders/shading.wgsl`, and
//! `shade` below for the CPU): `base * (m.ambient + m.diffuse *
//! (L.ambient + key n.l + fill n.l2))`, darkened toward the silhouette by
//! the material's outline, plus Blinn specular from the key light, then the
//! tonemap.

use bytemuck::{Pod, Zeroable};
pub use vv_core::material::MaterialPreset;

/// Up to this many directional lights can be live at once (`Lighting::
/// count`); the app widens/narrows the list, GPU code just loops
/// `0..count` and never branches on a per-light enabled flag.
pub const MAX_LIGHTS: usize = 4;

/// One directional light. Layout is part of `Lighting`'s (`shaders/
/// shading.wgsl` must match).
#[repr(C)]
#[derive(Clone, Copy, Debug, PartialEq, Default, Pod, Zeroable)]
pub struct Light {
    /// Unit direction toward the light, in view space.
    pub dir: [f32; 3],
    pub intensity: f32,
    /// Linear light colour -- decoded from the sRGB colour it was set
    /// with (`srgb_to_linear`) once by the caller, not per fragment.
    pub color: [f32; 3],
    pub _pad: f32,
}

impl Light {
    pub const WHITE: [f32; 3] = [1.0, 1.0, 1.0];

    fn new(dir: [f32; 3], intensity: f32) -> Self {
        Light {
            dir,
            intensity,
            color: Self::WHITE,
            _pad: 0.0,
        }
    }
}

/// Global lighting: up to [`MAX_LIGHTS`] directional lights (packed
/// densely from index 0, `count` of them live) plus ambient. Layout must
/// match `Lighting` in `shaders/shading.wgsl`.
#[repr(C)]
#[derive(Clone, Copy, Debug, PartialEq, Pod, Zeroable)]
pub struct Lighting {
    pub lights: [Light; MAX_LIGHTS],
    pub ambient_color: [f32; 3],
    pub ambient: f32,
    /// One of the `TONEMAP_*` selectors.
    pub tonemap: f32,
    /// How many of `lights` (from index 0) actually light the scene.
    pub count: f32,
    pub _pad: [f32; 2],
}

impl Default for Lighting {
    fn default() -> Self {
        LightingPreset::default().lighting()
    }
}

/// What an object is made of. Layout must match the three `vec4`s every
/// shader's params carry (`shaders/shading.wgsl`'s `material()`).
#[repr(C)]
#[derive(Clone, Copy, Debug, PartialEq, Pod, Zeroable)]
pub struct Material {
    /// Light the surface has regardless of the lights.
    pub ambient: f32,
    pub diffuse: f32,
    pub specular: f32,
    /// Specular exponent.
    pub shininess: f32,
    /// Diffuse quantized into this many bands (cel look); 0 is smooth.
    pub toon_bands: f32,
    pub opacity: f32,
    /// How much diffuse light fades toward the silhouette (can
    /// exceed 1).
    pub outline: f32,
    /// 0..1; wider toward 1.
    pub outline_width: f32,
    /// 1 = angle-dependent opacity (Merritt & Bacon 1997, Methods Enzymol.
    /// 277:505): clearer face-on, more opaque at grazing angles.
    pub transmode: f32,
    /// Extra light added at grazing angles, brightening the silhouette
    /// (Fresnel edge light); 0 is off.
    pub rim: f32,
    pub _pad: [f32; 2],
}

impl Material {
    /// Drawn in the transparent (glass) pass rather than the opaque one.
    pub fn is_glass(&self) -> bool {
        self.opacity < 1.0 || self.transmode > 0.5
    }
}

impl Default for Material {
    fn default() -> Self {
        MaterialPreset::default().into()
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum LightingPreset {
    /// A key light from the upper left front, a fill from the upper
    /// right front and an ambient floor (key:fill:ambient 0.65:0.35:0.25).
    /// Both lights are off-axis, not near the eye, so the surface shows a
    /// shading gradient. Brighter than `TwinLight`, whose two
    /// full-intensity lights and zero ambient leave anything facing away
    /// from both lights black.
    #[default]
    Default,
    /// Ambient only, strong occlusion: shape from shadowed crevices
    /// (key 0, fill 0, ambient 1.5).
    Soft,
    /// Strong key light with shadows and occlusion (key 0.7, fill 0.3,
    /// ambient 0.8, key-light shadows on).
    Full,
    /// Ambient only, with silhouettes. Ambient 1.5 with Opaque's diffuse
    /// of 0.65 puts ambient x diffuse near 1, so surface colour shows
    /// unshaded.
    Flat,
    /// Near-flat light, thick silhouettes, white background (Goodsell).
    Illustrative,
    /// `Soft`'s lights with softer occlusion: lower ambient-occlusion
    /// strength for a less sooty look.
    Gentle,
    /// Two full-intensity lights, no ambient term, depth cueing off.
    /// `Default` is deliberately brighter than this.
    TwinLight,
    /// A bright three-point-style look: a strong key, a fill opened up
    /// further than `Default`'s, and ambient standing in for a third,
    /// back/rim light. vizviz has only two directional lights plus
    /// ambient, so this approximates a three-point rig rather than
    /// adding a third light channel. (Named `Studio` in older sessions;
    /// renamed to leave that word for the Studio-mode feature.)
    ThreePoint,
}

impl LightingPreset {
    pub const ALL: [LightingPreset; 8] = [
        LightingPreset::Default,
        LightingPreset::Soft,
        LightingPreset::Full,
        LightingPreset::Flat,
        LightingPreset::Illustrative,
        LightingPreset::Gentle,
        LightingPreset::TwinLight,
        LightingPreset::ThreePoint,
    ];

    /// Old identifiers `parse()` still accepts; never returned by `name()`.
    const ALIASES: [(&'static str, LightingPreset); 1] = [("studio", LightingPreset::ThreePoint)];

    pub fn name(self) -> &'static str {
        match self {
            LightingPreset::Default => "default",
            LightingPreset::Soft => "soft",
            LightingPreset::Full => "full",
            LightingPreset::Flat => "flat",
            LightingPreset::Illustrative => "illustrative",
            LightingPreset::Gentle => "gentle",
            LightingPreset::TwinLight => "twinlight",
            LightingPreset::ThreePoint => "threepoint",
        }
    }

    pub fn parse(name: &str) -> Option<LightingPreset> {
        LightingPreset::ALL
            .into_iter()
            .find(|p| p.name() == name)
            .or_else(|| {
                LightingPreset::ALIASES
                    .iter()
                    .find(|(n, _)| *n == name)
                    .map(|(_, p)| *p)
            })
    }

    pub fn lighting(self) -> Lighting {
        let (key, fill, ambient) = match self {
            LightingPreset::Default => (0.65, 0.35, 0.25),
            LightingPreset::Soft | LightingPreset::Gentle => (0.0, 0.0, 1.5),
            LightingPreset::Full => (0.7, 0.3, 0.8),
            LightingPreset::Flat => (0.0, 0.0, 1.5),
            LightingPreset::Illustrative => (0.2, 0.0, 1.3),
            LightingPreset::TwinLight => (1.0, 1.0, 0.0),
            LightingPreset::ThreePoint => (1.0, 0.7, 0.5),
        };
        // Key nearly from the eye, fill from above right. Full's key comes
        // from higher and to the side instead, so its shadows fall visibly
        // beside their casters. `Default`'s key sits upper-left-front
        // rather than at the eye, so facing normals actually vary in
        // brightness. TwinLight and `Default` put the fill in front
        // (+0.5); every other preset's sits behind the scene (-0.5).
        let key_dir = match self {
            LightingPreset::Default => normalize([-0.5, 0.6, 0.6]),
            LightingPreset::Full => normalize([0.55, 0.65, 0.55]),
            _ => normalize([-0.1, 0.1, 1.0]),
        };
        let fill_dir = match self {
            LightingPreset::Default | LightingPreset::TwinLight => normalize([1.0, 2.0, 0.5]),
            _ => normalize([1.0, 2.0, -0.5]),
        };
        Lighting {
            lights: [
                Light::new(key_dir, key),
                Light::new(fill_dir, fill),
                Light::default(),
                Light::default(),
            ],
            ambient_color: Light::WHITE,
            ambient,
            tonemap: TONEMAP_ACES,
            count: 2.0,
            _pad: [0.0; 2],
        }
    }

    /// Ambient-occlusion strength (`RenderSettings::ao`).
    pub fn ao(self) -> f32 {
        match self {
            LightingPreset::Default | LightingPreset::Full | LightingPreset::ThreePoint => 1.0,
            LightingPreset::Soft => 1.6,
            LightingPreset::Gentle => 1.0,
            LightingPreset::Flat | LightingPreset::TwinLight => 0.0,
            LightingPreset::Illustrative => 1.2,
        }
    }

    /// Depth-cue strength (`RenderSettings::depth_cue`). Off for
    /// `TwinLight`, `Flat` and `Illustrative`; every other preset keeps
    /// it visible, strongest in `Full`'s dramatic look and weakest in `ThreePoint`'s clean
    /// product-shot look. `Default` is light, not `Full`'s 0.7: on a
    /// near-white background, 0.6 faded the far half of a molecule into
    /// the paper.
    pub fn depth_cue(self) -> f32 {
        match self {
            LightingPreset::Flat | LightingPreset::Illustrative | LightingPreset::TwinLight => 0.0,
            LightingPreset::Default => 0.3,
            LightingPreset::Full => 0.7,
            LightingPreset::ThreePoint => 0.4,
            _ => 0.6,
        }
    }

    /// Shadow strength (`RenderSettings::shadows`).
    pub fn shadows(self) -> f32 {
        match self {
            LightingPreset::Full => 1.0,
            LightingPreset::ThreePoint => 0.6,
            _ => 0.0,
        }
    }

    /// Whether silhouettes are on (`RenderSettings::outline`).
    pub fn outline(self) -> bool {
        matches!(self, LightingPreset::Flat | LightingPreset::Illustrative)
    }
}

impl From<MaterialPreset> for Material {
    fn from(preset: MaterialPreset) -> Material {
        let v = preset.settings();
        Material {
            ambient: v.ambient,
            diffuse: v.diffuse,
            specular: v.specular,
            // 0..1 shininess maps to the exponent 10^(3s): 0.53 -> 38.9,
            // 0.8 -> 251.2.
            shininess: 10f32.powf(3.0 * v.shininess),
            toon_bands: v.toon_bands,
            opacity: v.opacity,
            outline: v.outline,
            outline_width: v.outline_width,
            transmode: v.transmode,
            rim: v.rim,
            _pad: [0.0; 2],
        }
    }
}

/// A one-click look: lighting, material, background and post effects.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum StylePreset {
    #[default]
    DarkPresentation,
    PublicationWhite,
    Glossy,
    FlatCel,
}

impl StylePreset {
    pub const ALL: [StylePreset; 4] = [
        StylePreset::DarkPresentation,
        StylePreset::PublicationWhite,
        StylePreset::Glossy,
        StylePreset::FlatCel,
    ];

    pub fn label(self) -> &'static str {
        match self {
            StylePreset::DarkPresentation => "Dark Presentation",
            StylePreset::PublicationWhite => "Publication White",
            StylePreset::Glossy => "Glossy",
            StylePreset::FlatCel => "Flat / Cel",
        }
    }

    pub fn lighting_preset(self) -> LightingPreset {
        match self {
            StylePreset::FlatCel => LightingPreset::Flat,
            _ => LightingPreset::Default,
        }
    }

    pub fn material_preset(self) -> MaterialPreset {
        match self {
            StylePreset::DarkPresentation => MaterialPreset::Opaque,
            StylePreset::PublicationWhite => MaterialPreset::Diffuse,
            StylePreset::Glossy => MaterialPreset::Glossy,
            StylePreset::FlatCel => MaterialPreset::Toon,
        }
    }

    /// The preset's lighting. Tonemapping compresses highlights that
    /// would clip in the 8-bit target; the white and cel looks stay
    /// linear so their colors print as picked.
    pub fn lighting(self) -> Lighting {
        let mut lighting = self.lighting_preset().lighting();
        if matches!(self, StylePreset::PublicationWhite | StylePreset::FlatCel) {
            lighting.tonemap = TONEMAP_NONE;
        }
        lighting
    }

    pub fn material(self) -> Material {
        self.material_preset().into()
    }

    pub fn outline(self) -> bool {
        self.lighting_preset().outline()
    }

    pub fn ao(self) -> f32 {
        self.lighting_preset().ao()
    }

    pub fn depth_cue(self) -> f32 {
        self.lighting_preset().depth_cue()
    }

    /// The background the look implies; seeds `RenderSettings::background`
    /// when the look is picked (an export may still override it).
    pub fn background(self) -> wgpu::Color {
        let (r, g, b) = match self {
            StylePreset::DarkPresentation | StylePreset::Glossy => (0.09, 0.09, 0.11),
            StylePreset::PublicationWhite => (1.0, 1.0, 1.0),
            StylePreset::FlatCel => (0.82, 0.85, 0.90),
        };
        wgpu::Color { r, g, b, a: 1.0 }
    }
}

fn normalize(v: [f32; 3]) -> [f32; 3] {
    let l = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    [v[0] / l, v[1] / l, v[2] / l]
}

fn add3(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn scale3(a: [f32; 3], s: f32) -> [f32; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

fn mul3(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [a[0] * b[0], a[1] * b[1], a[2] * b[2]]
}

/// `Lighting::tonemap` values.
pub const TONEMAP_NONE: f32 = 0.0;
/// The ACES filmic fit (Narkowicz 2015): `x(2.51x+0.03) / (x(2.43x+0.59)+0.14)`.
pub const TONEMAP_ACES: f32 = 1.0;

/// Applies the curve `selector` names to one linear color channel. CPU
/// twin of `tonemap()` in `shaders/shading.wgsl`.
pub fn tonemap_channel(x: f32, selector: f32) -> f32 {
    if selector < 0.5 {
        return x;
    }
    ((x * (2.51 * x + 0.03)) / (x * (2.43 * x + 0.59) + 0.14)).clamp(0.0, 1.0)
}

/// CPU twin of `shade()` in `shaders/shading.wgsl` (`vv-cpu` uses it so
/// both backends agree). `base` is a palette (sRGB) colour, decoded here;
/// the shaders decode it where they unpack it. `n`, `view_dir` and the light directions are in
/// the same space; `view_dir` points from the eye into the scene.
pub fn shade(
    base: [f32; 3],
    n: [f32; 3],
    view_dir: [f32; 3],
    l: &Lighting,
    m: &Material,
) -> [f32; 3] {
    let dot = |a: [f32; 3], b: [f32; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    let mut diffuse_light = scale3(l.ambient_color, l.ambient);
    let count = (l.count.round() as usize).min(MAX_LIGHTS);
    for light in &l.lights[..count] {
        let mut ndotl = dot(n, light.dir).max(0.0);
        if m.toon_bands > 0.5 {
            ndotl = (ndotl * m.toon_bands).floor() / m.toon_bands;
        }
        diffuse_light = add3(diffuse_light, scale3(light.color, light.intensity * ndotl));
    }
    let mut light_v = add3([m.ambient; 3], mul3(diffuse_light, [m.diffuse; 3]));
    if m.outline > 0.0 {
        let facing = -dot(n, view_dir);
        let edge = 1.0
            - (1.0 - facing * facing)
                .max(0.0)
                .powf((1.0 - m.outline_width) * 32.0);
        light_v = scale3(light_v, (1.0 + (edge - 1.0) * m.outline).max(0.0));
    }
    if m.rim > 0.0 {
        let facing = (-dot(n, view_dir)).clamp(0.0, 1.0);
        light_v = add3(light_v, [m.rim * (1.0 - facing).powf(3.0); 3]);
    }
    // Specular stays achromatic and tied to the first light only, as
    // before (extending it to every light and colour costs a `pow` per
    // light per fragment for a highlight that reads the same either way).
    let key = l.lights[0];
    let h = normalize([
        key.dir[0] - view_dir[0],
        key.dir[1] - view_dir[1],
        key.dir[2] - view_dir[2],
    ]);
    let spec = m.specular * key.intensity * dot(n, h).max(0.0).powf(m.shininess);
    let mut out = [0.0; 3];
    for c in 0..3 {
        out[c] = tonemap_channel(srgb_to_linear(base[c]) * light_v[c] + spec, l.tonemap);
    }
    out
}

/// sRGB display value to linear light; CPU twin of `srgb_to_linear()` in
/// `shaders/shading.wgsl`. Palette colours are display values.
pub fn srgb_to_linear(c: f32) -> f32 {
    if c <= 0.04045 {
        c / 12.92
    } else {
        ((c + 0.055) / 1.055).powf(2.4)
    }
}

/// CPU twin of `material_alpha()` in `shaders/shading.wgsl`.
pub fn alpha(n: [f32; 3], view_dir: [f32; 3], m: &Material) -> f32 {
    if m.transmode < 0.5 {
        return m.opacity;
    }
    let facing = -(n[0] * view_dir[0] + n[1] * view_dir[1] + n[2] * view_dir[2]);
    let a = 1.0 + (std::f32::consts::PI * (1.0 - m.opacity) * facing).cos();
    a * a * 0.25
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn uniform_layouts() {
        assert_eq!(std::mem::size_of::<Light>(), 32);
        assert_eq!(std::mem::size_of::<Lighting>(), 160);
        assert_eq!(std::mem::size_of::<Material>(), 48);
    }

    #[test]
    fn presets_are_distinct_and_round_trip_by_name() {
        for (i, a) in MaterialPreset::ALL.iter().enumerate() {
            assert_eq!(MaterialPreset::parse(a.name()), Some(*a));
            for b in &MaterialPreset::ALL[i + 1..] {
                assert_ne!(Material::from(*a), Material::from(*b), "{a:?} vs {b:?}");
            }
        }
        // Distinct as a whole look, not just the lights: Soft and Gentle
        // share the same key/fill/ambient and differ only in `ao`.
        let signature = |p: LightingPreset| {
            (
                p.lighting(),
                p.ao(),
                p.depth_cue(),
                p.shadows(),
                p.outline(),
            )
        };
        for (i, a) in LightingPreset::ALL.iter().enumerate() {
            assert_eq!(LightingPreset::parse(a.name()), Some(*a));
            for b in &LightingPreset::ALL[i + 1..] {
                assert_ne!(signature(*a), signature(*b), "{a:?} vs {b:?}");
            }
        }
    }

    #[test]
    fn old_lighting_names_still_parse() {
        for (alias, preset) in LightingPreset::ALIASES {
            assert_eq!(
                LightingPreset::parse(alias),
                Some(preset),
                "alias {alias:?}"
            );
            assert!(
                LightingPreset::ALL.iter().all(|p| p.name() != alias),
                "{alias:?} is also a current lighting name"
            );
        }
    }

    #[test]
    fn toon_material_quantizes_the_key_light_into_bands() {
        let lighting = Lighting {
            count: 1.0, // the first light only
            ambient: 0.0,
            tonemap: TONEMAP_NONE,
            ..LightingPreset::Default.lighting()
        };
        let toon = Material::from(MaterialPreset::Toon);
        let mut seen = std::collections::BTreeSet::new();
        for i in 0..20 {
            let angle = (i as f32 / 19.0) * std::f32::consts::PI - std::f32::consts::FRAC_PI_2;
            let n = [angle.sin(), 0.0, angle.cos()];
            let shaded = shade([1.0; 3], n, [0.0, 0.0, -1.0], &lighting, &toon);
            seen.insert((shaded[0] * 1000.0).round() as i32);
        }
        assert!(seen.len() <= 4, "3 bands plus unlit, saw {seen:?}");
    }

    #[test]
    fn a_coloured_light_tints_the_shaded_result() {
        let mut lighting = Lighting {
            ambient: 0.0,
            tonemap: TONEMAP_NONE,
            count: 1.0,
            ..LightingPreset::Default.lighting()
        };
        lighting.lights[0].color = [1.0, 0.2, 0.2]; // red
        let m = Material::from(MaterialPreset::Diffuse);
        let n = [0.0, 0.0, 1.0];
        let shaded = shade([1.0; 3], n, [0.0, 0.0, -1.0], &lighting, &m);
        assert!(
            shaded[0] > shaded[1] && shaded[0] > shaded[2],
            "a red light should tint red brighter than green/blue: {shaded:?}"
        );
    }

    /// A disabled light (past `count`) contributes nothing, even though
    /// its own fields are still whatever the caller left them at.
    #[test]
    fn lights_past_count_are_ignored() {
        let mut one = LightingPreset::Default.lighting();
        one.count = 1.0;
        let mut two = one;
        two.lights[1] = two.lights[0]; // a second, identical light
        two.count = 1.0; // ...but still not counted
        let m = Material::from(MaterialPreset::Diffuse);
        let n = [0.0, 0.0, 1.0];
        assert_eq!(
            shade([1.0; 3], n, [0.0, 0.0, -1.0], &one, &m),
            shade([1.0; 3], n, [0.0, 0.0, -1.0], &two, &m),
        );
    }

    #[test]
    fn aces_compresses_highlights_and_keeps_midtones_close() {
        assert_eq!(tonemap_channel(0.5, TONEMAP_NONE), 0.5);
        let mid = tonemap_channel(0.18, TONEMAP_ACES);
        assert!((0.18..0.30).contains(&mid), "midtone 0.18 -> {mid}");
        let hot = tonemap_channel(3.0, TONEMAP_ACES);
        assert!(hot < 1.0 && hot > 0.9, "3.0 -> {hot}");
        let mut last = 0.0;
        for i in 1..=100 {
            let y = tonemap_channel(i as f32 * 0.05, TONEMAP_ACES);
            assert!(y >= last);
            last = y;
        }
    }
}
