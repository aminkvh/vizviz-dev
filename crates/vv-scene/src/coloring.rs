//! What colorings are made of beyond a plain scheme name: continuous
//! properties with a chosen ramp and range, and per-target overrides laid
//! over any scheme.

use std::fmt;

/// A color ramp's identity; `vv_render::color::ramp_color` knows its colors.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Ramp {
    BlueWhiteRed,
    RedWhiteBlue,
    TealWhiteGold,
    Viridis,
    Rainbow,
    Gray,
}

impl Ramp {
    pub const ALL: [Ramp; 6] = [
        Ramp::BlueWhiteRed,
        Ramp::RedWhiteBlue,
        Ramp::TealWhiteGold,
        Ramp::Viridis,
        Ramp::Rainbow,
        Ramp::Gray,
    ];

    pub fn name(self) -> &'static str {
        match self {
            Ramp::BlueWhiteRed => "bwr",
            Ramp::RedWhiteBlue => "rwb",
            Ramp::TealWhiteGold => "teal-gold",
            Ramp::Viridis => "viridis",
            Ramp::Rainbow => "rainbow",
            Ramp::Gray => "gray",
        }
    }

    /// The label the interface shows.
    pub fn label(self) -> &'static str {
        match self {
            Ramp::BlueWhiteRed => "Blue-white-red",
            Ramp::RedWhiteBlue => "Red-white-blue",
            Ramp::TealWhiteGold => "Teal-white-gold",
            Ramp::Viridis => "Viridis",
            Ramp::Rainbow => "Rainbow",
            Ramp::Gray => "Gray",
        }
    }

    pub fn parse(name: &str) -> Option<Ramp> {
        let name = name.to_ascii_lowercase();
        Ramp::ALL.into_iter().find(|r| r.name() == name)
    }
}

/// The hydrophobicity tables a `Hydrophobicity` property can use.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum HydroScale {
    /// Kyte & Doolittle (1982).
    KyteDoolittle,
    /// Wimley & White (1996), octanol.
    WimleyWhite,
    /// Eisenberg, Schwarz, Komaromy & Wall (1984) consensus.
    Eisenberg,
}

impl HydroScale {
    pub const ALL: [HydroScale; 3] = [
        HydroScale::KyteDoolittle,
        HydroScale::WimleyWhite,
        HydroScale::Eisenberg,
    ];

    pub fn name(self) -> &'static str {
        match self {
            HydroScale::KyteDoolittle => "kd",
            HydroScale::WimleyWhite => "ww",
            HydroScale::Eisenberg => "eisenberg",
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            HydroScale::KyteDoolittle => "Kyte-Doolittle",
            HydroScale::WimleyWhite => "Wimley-White",
            HydroScale::Eisenberg => "Eisenberg",
        }
    }

    pub fn parse(name: &str) -> Option<HydroScale> {
        let name = name.to_ascii_lowercase();
        let name = match name.as_str() {
            "kyte-doolittle" | "kytedoolittle" => "kd",
            "wimley-white" | "wimleywhite" | "wimley_white" => "ww",
            other => other,
        };
        HydroScale::ALL.into_iter().find(|s| s.name() == name)
    }
}

/// What a continuous coloring measures.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum PropertyKind {
    Hydrophobicity(HydroScale),
    /// Formal charge of the residue at pH 7.
    Charge,
    BFactor,
    Occupancy,
    /// Chou-Fasman alpha-helix propensity.
    HelixPropensity,
    /// Chou-Fasman beta-strand propensity.
    StrandPropensity,
    /// Chou-Fasman turn propensity.
    TurnPropensity,
    /// Fraction of residues of this type typically buried.
    BuriedIndex,
    /// The named per-atom value channel (`LoadedStructure::values`); falls
    /// back to element colors while none of that name is attached.
    Values(String),
}

impl PropertyKind {
    /// The parameterless kinds, for menus and tests.
    pub fn builtin() -> Vec<PropertyKind> {
        let mut all: Vec<PropertyKind> = HydroScale::ALL
            .into_iter()
            .map(PropertyKind::Hydrophobicity)
            .collect();
        all.extend([
            PropertyKind::Charge,
            PropertyKind::BFactor,
            PropertyKind::Occupancy,
            PropertyKind::HelixPropensity,
            PropertyKind::StrandPropensity,
            PropertyKind::TurnPropensity,
            PropertyKind::BuriedIndex,
        ]);
        all
    }

    /// The command word(s); the inverse of [`Property::parse_words`].
    pub fn name(&self) -> String {
        match self {
            PropertyKind::Hydrophobicity(HydroScale::KyteDoolittle) => "hydrophobicity".into(),
            PropertyKind::Hydrophobicity(s) => format!("hydrophobicity scale {}", s.name()),
            PropertyKind::Charge => "charge".into(),
            PropertyKind::BFactor => "b_factor".into(),
            PropertyKind::Occupancy => "occupancy".into(),
            PropertyKind::HelixPropensity => "helix".into(),
            PropertyKind::StrandPropensity => "strand".into(),
            PropertyKind::TurnPropensity => "turn".into(),
            PropertyKind::BuriedIndex => "buried".into(),
            PropertyKind::Values(name) => format!("values:{name}"),
        }
    }

    /// The label the interface shows.
    pub fn label(&self) -> String {
        match self {
            PropertyKind::Hydrophobicity(s) => format!("Hydrophobicity ({})", s.label()),
            PropertyKind::Charge => "Charge".into(),
            PropertyKind::BFactor => "B-factor".into(),
            PropertyKind::Occupancy => "Occupancy".into(),
            PropertyKind::HelixPropensity => "Helix propensity".into(),
            PropertyKind::StrandPropensity => "Strand propensity".into(),
            PropertyKind::TurnPropensity => "Turn propensity".into(),
            PropertyKind::BuriedIndex => "Buried index".into(),
            PropertyKind::Values(name) => format!("Values: {name}"),
        }
    }
}

/// A continuous coloring: a property mapped onto a ramp over a range.
/// `None` ramp or range means the property's own default.
#[derive(Clone, Debug, PartialEq)]
pub struct Property {
    pub kind: PropertyKind,
    pub ramp: Option<Ramp>,
    pub range: Option<[f32; 2]>,
}

impl Property {
    pub fn new(kind: PropertyKind) -> Property {
        Property {
            kind,
            ramp: None,
            range: None,
        }
    }

    /// `name()`: the kind, then ` range LO HI` and ` ramp NAME` when set.
    pub fn name(&self) -> String {
        let mut out = self.kind.name();
        if let Some([lo, hi]) = self.range {
            out.push_str(&format!(" range {lo} {hi}"));
        }
        if let Some(ramp) = self.ramp {
            out.push_str(&format!(" ramp {}", ramp.name()));
        }
        out
    }

    /// A property from the front of `words` and how many words it used:
    /// a kind word (`hydrophobicity`, `ww`, `b_factor`, `values NAME`,
    /// ...), then any of `scale NAME`, `range LO HI`, `ramp NAME`.
    pub fn parse_words(words: &[&str]) -> Option<(Property, usize)> {
        let (kind, mut used) = parse_kind(words)?;
        let mut property = Property::new(kind);
        while let Some(&word) = words.get(used) {
            match (word, &words[used + 1..]) {
                ("scale", [name, ..]) => match (&mut property.kind, HydroScale::parse(name)) {
                    (PropertyKind::Hydrophobicity(scale), Some(parsed)) => *scale = parsed,
                    _ => return None,
                },
                ("ramp", [name, ..]) => property.ramp = Some(Ramp::parse(name)?),
                ("range", [lo, hi, ..]) => {
                    let (lo, hi) = (lo.parse().ok()?, hi.parse().ok()?);
                    property.range = Some([lo, hi]);
                    used += 1;
                }
                _ => break,
            }
            used += 2;
        }
        Some((property, used))
    }
}

fn parse_kind(words: &[&str]) -> Option<(PropertyKind, usize)> {
    let first = *words.first()?;
    let one = |kind| Some((kind, 1));
    match first {
        "hydrophobicity" | "hydropathy" | "kd" => {
            one(PropertyKind::Hydrophobicity(HydroScale::KyteDoolittle))
        }
        "ww" | "wimleywhite" | "wimley_white" => {
            one(PropertyKind::Hydrophobicity(HydroScale::WimleyWhite))
        }
        "eisenberg" => one(PropertyKind::Hydrophobicity(HydroScale::Eisenberg)),
        "charge" => one(PropertyKind::Charge),
        "b_factor" | "bfactor" | "b-factor" => one(PropertyKind::BFactor),
        "occupancy" => one(PropertyKind::Occupancy),
        "helix" | "helixpropensity" => one(PropertyKind::HelixPropensity),
        "strand" | "strandpropensity" => one(PropertyKind::StrandPropensity),
        "turn" | "turnpropensity" => one(PropertyKind::TurnPropensity),
        "buried" | "buriedindex" => one(PropertyKind::BuriedIndex),
        "values" => {
            let name = words.get(1).filter(|n| !n.is_empty())?;
            Some((PropertyKind::Values((*name).to_string()), 2))
        }
        other => other
            .strip_prefix("values:")
            .filter(|n| !n.is_empty())
            .map(|n| (PropertyKind::Values(n.to_string()), 1)),
    }
}

/// A set of atoms a color override names.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ColorTarget {
    Element(String),
    /// An atom name (`CA`).
    Name(String),
    Resname(String),
    /// A 0-based atom index.
    Atom(u32),
    /// Any selection expression.
    Selection(String),
}

impl ColorTarget {
    /// `kind` (`element`, `name`, `resname`, `atom`, `sel`) and its argument.
    pub fn parse(kind: &str, arg: &str) -> Result<ColorTarget, String> {
        let arg = arg.trim();
        if arg.is_empty() {
            return Err(format!("`{kind}` needs a value"));
        }
        Ok(match kind {
            "element" => ColorTarget::Element(arg.to_string()),
            "name" => ColorTarget::Name(arg.to_string()),
            "resname" => ColorTarget::Resname(arg.to_string()),
            "atom" => ColorTarget::Atom(
                arg.parse()
                    .map_err(|_| format!("atom needs a 0-based index, not `{arg}`"))?,
            ),
            "sel" => ColorTarget::Selection(arg.to_string()),
            other => {
                return Err(format!(
                    "unknown target `{other}`; expected element, name, resname, atom or sel"
                ))
            }
        })
    }

    pub fn kind(&self) -> &'static str {
        match self {
            ColorTarget::Element(_) => "element",
            ColorTarget::Name(_) => "name",
            ColorTarget::Resname(_) => "resname",
            ColorTarget::Atom(_) => "atom",
            ColorTarget::Selection(_) => "sel",
        }
    }

    /// The selection expression that picks this target's atoms.
    pub fn expression(&self) -> String {
        match self {
            ColorTarget::Element(e) => format!("element {e}"),
            ColorTarget::Name(n) => format!("name {n}"),
            ColorTarget::Resname(r) => format!("resname {r}"),
            ColorTarget::Atom(i) => format!("index {i}"),
            ColorTarget::Selection(expr) => expr.clone(),
        }
    }
}

/// `name CA`, as the command language and session files write it.
impl fmt::Display for ColorTarget {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ColorTarget::Element(v) | ColorTarget::Name(v) | ColorTarget::Resname(v) => {
                write!(f, "{} {v}", self.kind())
            }
            ColorTarget::Atom(i) => write!(f, "atom {i}"),
            ColorTarget::Selection(expr) => write!(f, "sel {expr}"),
        }
    }
}

impl std::str::FromStr for ColorTarget {
    type Err = String;

    fn from_str(text: &str) -> Result<Self, String> {
        let (kind, arg) = text
            .trim()
            .split_once(char::is_whitespace)
            .unwrap_or((text, ""));
        ColorTarget::parse(kind, arg)
    }
}

/// `color` laid over a structure's base coloring on the atoms `target`
/// names; later overrides win where they overlap.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ColorOverride {
    pub target: ColorTarget,
    pub color: [u8; 3],
}

#[cfg(test)]
mod tests {
    use super::*;

    fn words(text: &str) -> Vec<&str> {
        text.split_whitespace().collect()
    }

    #[test]
    fn property_parameters_round_trip_through_their_name() {
        let mut property = Property::new(PropertyKind::Hydrophobicity(HydroScale::WimleyWhite));
        property.range = Some([-2.5, 2.5]);
        property.ramp = Some(Ramp::Viridis);
        let name = property.name();
        assert_eq!(name, "hydrophobicity scale ww range -2.5 2.5 ramp viridis");
        let (parsed, used) = Property::parse_words(&words(&name)).unwrap();
        assert_eq!((parsed, used), (property, 8));
    }

    #[test]
    fn old_spellings_and_trailing_words_are_left_alone() {
        let (p, used) = Property::parse_words(&words("ww 3")).unwrap();
        assert_eq!(
            p.kind,
            PropertyKind::Hydrophobicity(HydroScale::WimleyWhite)
        );
        assert_eq!(used, 1, "a trailing structure id is not a parameter");
        let (p, used) = Property::parse_words(&words("values sasa ramp gray 2")).unwrap();
        assert_eq!(p.kind, PropertyKind::Values("sasa".into()));
        assert_eq!(used, 4);
        assert!(Property::parse_words(&words("chain")).is_none());
        assert!(Property::parse_words(&words("charge scale kd")).is_none());
        assert!(Property::parse_words(&words("charge ramp nope")).is_none());
    }

    #[test]
    fn targets_parse_print_and_select() {
        for text in [
            "element C",
            "name CA",
            "resname HEM",
            "atom 123",
            "sel chain A and helix",
        ] {
            let target: ColorTarget = text.parse().unwrap();
            assert_eq!(target.to_string(), text);
        }
        assert_eq!(
            "atom 7".parse::<ColorTarget>().unwrap().expression(),
            "index 7"
        );
        assert!("atom x".parse::<ColorTarget>().is_err());
        assert!("bogus 1".parse::<ColorTarget>().is_err());
        assert!("name".parse::<ColorTarget>().is_err());
    }
}
