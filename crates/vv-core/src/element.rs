//! Chemical elements: symbol <-> atomic number.
//!
//! Radii and colors are added later with cited sources (see docs/VALIDATION.md).

/// A chemical element identified by atomic number. `0` means unknown.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug, Default, PartialOrd, Ord)]
// `Pod` lets `&[Element]` be viewed as `&[u8]` without a copy (the Python
// bindings expose it as a NumPy array); sound because of `repr(transparent)`.
#[derive(bytemuck::Pod, bytemuck::Zeroable)]
#[repr(transparent)]
pub struct Element(u8);

/// Index = atomic number; index 0 is the unknown placeholder.
const SYMBOLS: [&str; 119] = [
    "X", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S",
    "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge",
    "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd",
    "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm",
    "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn",
    "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
];

const FALLBACK_RADIUS: f32 = 2.0;

/// `Element::vdw_radius`, indexed by atomic number (0 = unknown).
const VDW_RADII: [f32; 101] = [
    FALLBACK_RADIUS,
    1.10,
    1.40,
    1.82,
    1.53,
    1.92,
    1.70, // H He Li Be B C
    1.55,
    1.52,
    1.47,
    1.54,
    2.27,
    1.73, // N O F Ne Na Mg
    1.84,
    2.10,
    1.80,
    1.80,
    1.75,
    1.88, // Al Si P S Cl Ar
    2.75,
    2.31,
    2.30,
    2.15,
    2.05,
    2.05, // K Ca Sc Ti V Cr
    2.05,
    2.05,
    2.00,
    1.63,
    1.40,
    1.39, // Mn Fe Co Ni Cu Zn
    1.87,
    2.11,
    1.85,
    1.90,
    1.85,
    2.02, // Ga Ge As Se Br Kr
    3.03,
    2.49,
    2.40,
    2.30,
    2.15,
    2.10, // Rb Sr Y Zr Nb Mo
    2.05,
    2.05,
    2.00,
    1.63,
    1.72,
    1.58, // Tc Ru Rh Pd Ag Cd
    1.93,
    2.17,
    2.06,
    2.06,
    1.98,
    2.16, // In Sn Sb Te I Xe
    3.43,
    2.68,
    2.98,
    2.88,
    2.92,
    2.95, // Cs Ba La Ce Pr Nd
    2.00,
    2.90,
    2.87,
    2.83,
    2.79,
    2.87, // Pm Sm Eu Gd Tb Dy
    2.81,
    2.83,
    2.79,
    2.80,
    2.74,
    2.25, // Ho Er Tm Yb Lu Hf
    2.20,
    2.10,
    2.05,
    2.00,
    2.00,
    1.72, // Ta W Re Os Ir Pt
    1.66,
    1.55,
    1.96,
    2.02,
    2.07,
    1.97, // Au Hg Tl Pb Bi Po
    2.02,
    2.20,
    3.48,
    2.83,
    2.80,
    2.93, // At Rn Fr Ra Ac Th
    2.88,
    1.86,
    2.82,
    2.81,
    2.83,
    3.05, // Pa U Np Pu Am Cm
    3.40,
    3.05,
    2.70,
    2.00, // Bk Cf Es Fm
];

impl Element {
    pub const UNKNOWN: Element = Element(0);
    pub const HYDROGEN: Element = Element(1);
    pub const CARBON: Element = Element(6);
    pub const NITROGEN: Element = Element(7);
    pub const OXYGEN: Element = Element(8);
    pub const SULFUR: Element = Element(16);
    pub const MAX_ATOMIC_NUMBER: u8 = 118;

    pub const fn from_atomic_number(z: u8) -> Option<Element> {
        if z >= 1 && z <= Self::MAX_ATOMIC_NUMBER {
            Some(Element(z))
        } else {
            None
        }
    }

    pub const fn atomic_number(self) -> u8 {
        self.0
    }

    pub const fn is_unknown(self) -> bool {
        self.0 == 0
    }

    pub const fn is_hydrogen(self) -> bool {
        self.0 == 1
    }

    pub fn symbol(self) -> &'static str {
        SYMBOLS[self.0 as usize]
    }

    /// Single-bond covalent radius in angstroms, used for bond perception.
    ///
    /// Source: Cordero et al., Dalton Trans. 2008, 2832 (low-spin values
    /// for Mn/Fe/Co; sp3 for carbon). Elements not listed fall back to
    /// 1.5 Å, a display convention, not a measured value.
    pub fn covalent_radius(self) -> f32 {
        match self.0 {
            1 => 0.31,
            2 => 0.28,
            3 => 1.28,
            4 => 0.96,
            5 => 0.84,
            6 => 0.76,
            7 => 0.71,
            8 => 0.66,
            9 => 0.57,
            10 => 0.58,
            11 => 1.66,
            12 => 1.41,
            13 => 1.21,
            14 => 1.11,
            15 => 1.07,
            16 => 1.05,
            17 => 1.02,
            18 => 1.06,
            19 => 2.03,
            20 => 1.76,
            21 => 1.70,
            22 => 1.60,
            23 => 1.53,
            24 => 1.39,
            25 => 1.39,
            26 => 1.32,
            27 => 1.26,
            28 => 1.24,
            29 => 1.32,
            30 => 1.22,
            31 => 1.22,
            32 => 1.20,
            33 => 1.19,
            34 => 1.20,
            35 => 1.20,
            36 => 1.16,
            37 => 2.20,
            38 => 1.95,
            42 => 1.54,
            44 => 1.46,
            45 => 1.42,
            46 => 1.39,
            47 => 1.45,
            48 => 1.44,
            50 => 1.39,
            53 => 1.39,
            55 => 2.44,
            56 => 2.15,
            74 => 1.62,
            77 => 1.41,
            78 => 1.36,
            79 => 1.36,
            80 => 1.32,
            82 => 1.46,
            83 => 1.48,
            _ => 1.50,
        }
    }

    /// Metals get looser bond-perception rules (coordination, not covalent).
    pub fn is_metal(self) -> bool {
        matches!(
            self.0,
            3 | 4 | 11..=13 | 19..=31 | 37..=50 | 55..=83
        )
    }

    /// Van der Waals radius in angstroms, for every element up to Fm.
    ///
    /// Bondi 1964 (with its own Ni-U metals), then Mantina et al. 2009 for
    /// main-group gaps, Batsanov 2001
    /// for the d block, Alvarez 2013 for the f block; H is Rowland & Taylor
    /// 1996. Per-element sources: docs/RADII_RESEARCH.md. Pm, Z >= 100 and
    /// unknown atoms have no measured value: 2.0 Å is a display convention.
    pub fn vdw_radius(self) -> f32 {
        VDW_RADII
            .get(self.0 as usize)
            .copied()
            .unwrap_or(FALLBACK_RADIUS)
    }

    /// Case-insensitive lookup ("FE", "Fe", "fe" all work). Surrounding ASCII
    /// whitespace is ignored. Anything that is not an element symbol maps to
    /// [`Element::UNKNOWN`].
    pub fn from_symbol(symbol: &[u8]) -> Element {
        let s = symbol.trim_ascii();
        // Deuterium (PDB element column, mmCIF `type_symbol`) is hydrogen.
        if matches!(s, [b'D' | b'd']) {
            return Element::HYDROGEN;
        }
        let (a, b) = match s {
            [a] => (a.to_ascii_uppercase(), 0u8),
            [a, b] => (a.to_ascii_uppercase(), b.to_ascii_lowercase()),
            _ => return Element::UNKNOWN,
        };
        for (z, sym) in SYMBOLS.iter().enumerate().skip(1) {
            let sb = sym.as_bytes();
            let matches = match sb {
                [x] => *x == a && b == 0,
                [x, y] => *x == a && *y == b,
                _ => false,
            };
            if matches {
                return Element(z as u8);
            }
        }
        Element::UNKNOWN
    }

    /// Best guess from an atom name when a file (or an array from another
    /// package) carries no element column: `CA` is carbon, `FE` iron,
    /// `1HB` hydrogen. Two-letter symbols are rare inside polymers, so the
    /// first letter wins unless the pair is a known two-letter ion.
    ///
    /// A bare two-letter name (`HG`, `FE`, `ZN`, ...) is read as that ion,
    /// but one with more name left over after it (`HG11`, `HE21`: Val/Ile's
    /// and Gln's branched hydrogens) is hydrogen with a numeric locant, not
    /// mercury or helium with one - no residue in a biomolecule has an
    /// atom named that closely after a bonded hydrogen, so the extra
    /// characters are the tell.
    pub fn from_atom_name(name: &[u8]) -> Element {
        let trimmed = name.trim_ascii();
        let letters: Vec<u8> = trimmed
            .iter()
            .copied()
            .filter(u8::is_ascii_alphabetic)
            .take(2)
            .collect();
        if letters.is_empty() {
            return Element::UNKNOWN;
        }
        let two = Element::from_symbol(&letters);
        let one = Element::from_symbol(&letters[..1]);
        let hydrogen_with_locant =
            one == Element::HYDROGEN && !two.is_unknown() && trimmed.len() > 2;
        if letters.len() == 2
            && !two.is_unknown()
            && !hydrogen_with_locant
            && (one.is_unknown() || is_two_letter_ion(&letters))
        {
            two
        } else {
            one
        }
    }
}

/// Standard atomic weight in daltons, index = atomic number (index 0
/// unused). Source: CIAAW/IUPAC 2021 conventional values (radioactive
/// elements with no stable isotope, 43/61/84-86, use their most common
/// isotope's mass instead). Covers H through Rn (1-86); an element beyond
/// that range, or one this repository has no plausible use for, is not
/// tabulated and [`Element::from_mass`] never returns it.
const STANDARD_WEIGHT: [f32; 87] = [
    0.0, 1.008, 4.0026, 6.94, 9.0122, 10.81, 12.011, 14.007, 15.999, 18.998, 20.180, 22.990,
    24.305, 26.982, 28.085, 30.974, 32.06, 35.45, 39.948, 39.098, 40.078, 44.956, 47.867, 50.942,
    51.996, 54.938, 55.845, 58.933, 58.693, 63.546, 65.38, 69.723, 72.630, 74.922, 78.971, 79.904,
    83.798, 85.468, 87.62, 88.906, 91.224, 92.906, 95.95, 98.0, 101.07, 102.91, 106.42, 107.87,
    112.41, 114.82, 118.71, 121.76, 127.60, 126.90, 131.29, 132.91, 137.33, 138.91, 140.12, 140.91,
    144.24, 145.0, 150.36, 151.96, 157.25, 158.93, 162.50, 164.93, 167.26, 168.93, 173.05, 174.97,
    178.49, 180.95, 183.84, 186.21, 190.23, 192.22, 195.08, 196.97, 200.59, 204.38, 207.2, 208.98,
    209.0, 210.0, 222.0,
];

/// A mass this far from any tabulated standard atomic weight is not a
/// confident match (e.g. hydrogen mass-repartitioned to ~3 Da, or a
/// virtual site/lone pair at ~0): callers should fall back to the atom
/// name instead of trusting a nearest-neighbor guess.
const MASS_TOLERANCE: f32 = 0.5;

impl Element {
    /// Best guess from an atomic mass in daltons, for formats that carry a
    /// mass but not an atomic number (PSF, and PRMTOP without
    /// `ATOMIC_NUMBER`): the tabulated element ([`STANDARD_WEIGHT`])
    /// nearest `mass`, or [`Element::UNKNOWN`] when none is within
    /// [`MASS_TOLERANCE`]. Hydrogen mass repartitioning (heavy-hydrogen MD
    /// schemes shift mass from a heavy atom to its bonded hydrogens) can
    /// defeat this for the affected atoms; the atom name is the fallback.
    pub fn from_mass(mass: f32) -> Element {
        STANDARD_WEIGHT
            .iter()
            .enumerate()
            .skip(1)
            .map(|(z, &w)| (z, (mass - w).abs()))
            .min_by(|a, b| a.1.total_cmp(&b.1))
            .filter(|&(_, diff)| diff <= MASS_TOLERANCE)
            .map(|(z, _)| Element(z as u8))
            .unwrap_or(Element::UNKNOWN)
    }
}

/// Monatomic ion residue names, mapped straight to their element: PDB CCD
/// codes (`NA`, `CL`, `ZN`, ...), CHARMM's own names (`SOD`, `CLA`, `POT`,
/// `CAL`), and Amber's charge-suffixed `atomic_ions.lib` names (`Na+`,
/// `Cl-`, `MG2+`, ...). Deliberately keyed by *residue* name, not atom
/// name: [`Element::from_atom_name`] cannot tell a calcium ion's atom
/// named `CA` from an alpha carbon named the same, but a residue named
/// `CA` with one atom is unambiguous. See [`ion_element`].
const ION_ALIASES: &[(&str, Element)] = &[
    // Plain element-symbol residue names (PDB CCD and GROMACS ion .itp
    // files); robust to a source file's column alignment, unlike the
    // legacy-PDB fixed-column convention `vv_io::pdb` otherwise relies on.
    ("LI", Element(3)),
    ("NA", Element(11)),
    ("MG", Element(12)),
    ("AL", Element(13)),
    ("CL", Element(17)),
    ("K", Element(19)),
    ("CA", Element(20)),
    ("MN", Element(25)),
    ("FE", Element(26)),
    ("FE2", Element(26)),
    ("CO", Element(27)),
    ("NI", Element(28)),
    ("CU", Element(29)),
    ("CU1", Element(29)),
    ("ZN", Element(30)),
    ("BR", Element(35)),
    ("RB", Element(37)),
    ("SR", Element(38)),
    ("AG", Element(47)),
    ("CD", Element(48)),
    ("IOD", Element(53)),
    ("CS", Element(55)),
    ("BA", Element(56)),
    ("PT", Element(78)),
    ("AU", Element(79)),
    ("HG", Element(80)),
    ("PB", Element(82)),
    // CHARMM
    ("SOD", Element(11)),
    ("POT", Element(19)),
    ("CAL", Element(20)),
    ("CLA", Element(17)),
    ("ZN2", Element(30)),
    ("CES", Element(55)),
    ("RUB", Element(37)),
    ("LIT", Element(3)),
    ("BAR", Element(56)),
    ("CD2", Element(48)),
    ("F", Element(9)),
    // Amber `atomic_ions.lib`
    ("LI+", Element(3)),
    ("NA+", Element(11)),
    ("MG2+", Element(12)),
    ("AL3+", Element(13)),
    ("CL-", Element(17)),
    ("K+", Element(19)),
    ("CA2+", Element(20)),
    ("MN2+", Element(25)),
    ("FE2+", Element(26)),
    ("FE3+", Element(26)),
    ("CO2+", Element(27)),
    ("NI2+", Element(28)),
    ("CU+", Element(29)),
    ("CU2+", Element(29)),
    ("ZN2+", Element(30)),
    ("BR-", Element(35)),
    ("RB+", Element(37)),
    ("SR2+", Element(38)),
    ("AG+", Element(47)),
    ("CD2+", Element(48)),
    ("I-", Element(53)),
    ("CS+", Element(55)),
    ("BA2+", Element(56)),
    ("AU+", Element(79)),
    ("AU3+", Element(79)),
    ("HG2+", Element(80)),
    ("PB2+", Element(82)),
    ("F-", Element(9)),
];

/// The element a monatomic ion residue named `resname` is, or `None` when
/// it is not one of [`ION_ALIASES`]. Case-insensitive (Amber's own files
/// are usually mixed-case, e.g. `Na+`).
///
/// Callers must also check that the residue really has exactly one atom
/// before trusting this: `CO` and `NI` are also plausible ligand or side
/// chain atom names, and only a single-atom residue makes the ion
/// reading unambiguous (`vv_core::builder::fix_ion_elements`).
pub fn ion_element(resname: &str) -> Option<Element> {
    ION_ALIASES
        .iter()
        .find(|(name, _)| name.eq_ignore_ascii_case(resname))
        .map(|&(_, e)| e)
}

fn is_two_letter_ion(letters: &[u8]) -> bool {
    matches!(
        &letters.to_ascii_uppercase()[..],
        b"FE"
            | b"ZN"
            | b"MG"
            | b"MN"
            | b"CU"
            | b"NI"
            | b"CO"
            | b"CL"
            | b"BR"
            | b"NA"
            | b"SE"
            | b"CD"
            | b"HG"
            | b"PT"
            | b"AU"
            | b"AG"
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn deuterium_reads_as_hydrogen() {
        assert_eq!(Element::from_symbol(b"D"), Element::HYDROGEN);
        assert_eq!(Element::from_symbol(b" d"), Element::HYDROGEN);
        assert_eq!(Element::from_atom_name(b" D1 "), Element::HYDROGEN);
    }

    #[test]
    fn every_element_with_published_data_has_its_own_radius() {
        let sym = |s: &str| Element::from_symbol(s.as_bytes());
        assert_eq!(sym("C").vdw_radius(), 1.70);
        assert_eq!(sym("H").vdw_radius(), 1.10);
        assert_eq!(sym("Fe").vdw_radius(), 2.05);
        assert_eq!(sym("Zn").vdw_radius(), 1.39);
        assert_eq!(sym("Gd").vdw_radius(), 2.83);
        for z in 1..=100u8 {
            let r = Element(z).vdw_radius();
            assert!((1.0..3.6).contains(&r), "Z={z}: {r}");
        }
        assert_eq!(Element::UNKNOWN.vdw_radius(), FALLBACK_RADIUS);
    }

    #[test]
    fn element_from_atom_name() {
        assert_eq!(Element::from_atom_name(b"CA"), Element::CARBON);
        assert_eq!(Element::from_atom_name(b"N"), Element::NITROGEN);
        assert_eq!(Element::from_atom_name(b"FE"), Element::from_symbol(b"Fe"));
        assert_eq!(Element::from_atom_name(b"C1'"), Element::CARBON);
        assert_eq!(Element::from_atom_name(b"1HB"), Element::HYDROGEN);
        assert_eq!(Element::from_atom_name(b"ZN"), Element::from_symbol(b"Zn"));
        assert_eq!(Element::from_atom_name(b"OXT"), Element::OXYGEN);
        assert_eq!(Element::from_atom_name(b"123"), Element::UNKNOWN);
    }

    #[test]
    fn symbol_roundtrip_all_elements() {
        for z in 1..=Element::MAX_ATOMIC_NUMBER {
            let e = Element::from_atomic_number(z).unwrap();
            assert_eq!(Element::from_symbol(e.symbol().as_bytes()), e, "z = {z}");
        }
    }

    #[test]
    fn symbol_lookup_is_case_insensitive_and_trims() {
        assert_eq!(
            Element::from_symbol(b"FE"),
            Element::from_atomic_number(26).unwrap()
        );
        assert_eq!(
            Element::from_symbol(b"fe"),
            Element::from_atomic_number(26).unwrap()
        );
        assert_eq!(Element::from_symbol(b" C "), Element::CARBON);
        assert_eq!(Element::from_symbol(b"c"), Element::CARBON);
    }

    #[test]
    fn unknown_symbols_map_to_unknown() {
        assert_eq!(Element::from_symbol(b""), Element::UNKNOWN);
        assert_eq!(Element::from_symbol(b"Xx"), Element::UNKNOWN);
        assert_eq!(Element::from_symbol(b"Fee"), Element::UNKNOWN);
        assert_eq!(Element::from_symbol(b"X"), Element::UNKNOWN);
        assert!(Element::from_atomic_number(0).is_none());
        assert!(Element::from_atomic_number(119).is_none());
    }

    #[test]
    fn branched_hydrogen_locants_are_not_read_as_helium_or_mercury() {
        // Gln's amide hydrogens HE21/HE22, Val/Ile's branched methyl
        // hydrogens HG11/HG12/HG13/HG21/HG22/HG23: "HE"/"HG" spell real
        // elements, but the trailing locant digits mark these as hydrogen.
        for name in [
            "HE21", "HE22", "HG11", "HG12", "HG13", "HG21", "HG22", "HG23",
        ] {
            assert_eq!(
                Element::from_atom_name(name.as_bytes()),
                Element::HYDROGEN,
                "{name}"
            );
        }
        // A bare two-letter name is still read as the ion: no trailing
        // locant to say otherwise.
        assert_eq!(Element::from_atom_name(b"HG"), Element::from_symbol(b"Hg"));
        assert_eq!(Element::from_atom_name(b"FE"), Element::from_symbol(b"Fe"));
    }

    #[test]
    fn ion_element_covers_pdb_charmm_and_amber_names() {
        assert_eq!(
            ion_element("CA"),
            Some(Element::from_atomic_number(20).unwrap())
        ); // PDB CCD calcium
        assert_eq!(
            ion_element("ca"),
            Some(Element::from_atomic_number(20).unwrap())
        );
        assert_eq!(
            ion_element("SOD"),
            Some(Element::from_atomic_number(11).unwrap())
        ); // CHARMM sodium
        assert_eq!(
            ion_element("CLA"),
            Some(Element::from_atomic_number(17).unwrap())
        ); // CHARMM chloride
        assert_eq!(
            ion_element("POT"),
            Some(Element::from_atomic_number(19).unwrap())
        ); // CHARMM potassium
        assert_eq!(
            ion_element("Na+"),
            Some(Element::from_atomic_number(11).unwrap())
        ); // Amber
        assert_eq!(
            ion_element("MG2+"),
            Some(Element::from_atomic_number(12).unwrap())
        );
        assert_eq!(ion_element("ALA"), None);
    }

    #[test]
    fn from_mass_finds_the_nearest_standard_atomic_weight() {
        assert_eq!(Element::from_mass(1.008), Element::HYDROGEN);
        assert_eq!(Element::from_mass(12.011), Element::CARBON);
        assert_eq!(Element::from_mass(14.007), Element::NITROGEN);
        assert_eq!(Element::from_mass(15.9994), Element::OXYGEN);
        assert_eq!(Element::from_mass(32.06), Element::SULFUR);
        assert_eq!(
            Element::from_mass(30.974),
            Element::from_atomic_number(15).unwrap()
        ); // P
        assert_eq!(
            Element::from_mass(22.98977),
            Element::from_atomic_number(11).unwrap()
        ); // Na
        assert_eq!(
            Element::from_mass(35.45),
            Element::from_atomic_number(17).unwrap()
        ); // Cl
        assert_eq!(
            Element::from_mass(65.38),
            Element::from_atomic_number(30).unwrap()
        ); // Zn
    }

    #[test]
    fn from_mass_rejects_a_mass_too_far_from_any_element() {
        // Hydrogen mass-repartitioned to carry a heavy atom's shifted mass,
        // and a virtual site / lone pair with ~0 mass: neither is close
        // enough to any tabulated weight to be a confident element guess.
        assert_eq!(Element::from_mass(3.024), Element::UNKNOWN);
        assert_eq!(Element::from_mass(0.0), Element::UNKNOWN);
    }

    #[test]
    fn single_letter_symbol_does_not_match_two_letter_prefix() {
        // "C" must be carbon, not the first two-letter element starting with C.
        assert_eq!(Element::from_symbol(b"C"), Element::CARBON);
        assert_eq!(Element::from_symbol(b"N"), Element::NITROGEN);
        assert_eq!(Element::from_symbol(b"Na").symbol(), "Na");
    }
}
