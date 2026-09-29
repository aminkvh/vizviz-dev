//! Residue-name tables, one entry per name, grouped by the source that
//! defines them. Lookup is case-insensitive; `super::of_name` applies the
//! groups in the order [`GROUPS`] lists them, so the first class wins.
//! Monosaccharides live in `crate::glycan`'s own table, monatomic ions in
//! `crate::element::ion_element`; neither is repeated here.

use super::ResidueClass;

// ---- protein ---------------------------------------------------------

/// The 20 standard amino acids plus the selenocysteine/pyrrolysine codes
/// and ambiguity codes of IUPAC-IUB 1983 (Eur. J. Biochem. 138:9).
const AMINO_ACIDS: &[&str] = &[
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET",
    "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "SEC", "PYL", "ASX", "GLX", "UNK",
];

/// Protonation-state and disulfide variants: CHARMM (MacKerell et al.
/// 1998, J Phys Chem B 102:3586), Amber (`leaprc.protein.*`), and the
/// GROMACS force-field residue databases (`aminoacids.rtp`).
const PROTEIN_VARIANTS: &[&str] = &[
    // CHARMM: neutral His delta/epsilon, doubly protonated His, protonated
    // Asp/Glu, neutral Lys, and the `CYS2` spelling of a disulfide Cys.
    "HSD", "HSE", "HSP", "ASPP", "GLUP", "LSN", "CYS2",
    // Amber: His tautomers/protonated, disulfide Cys, thiolate Cys,
    // neutral Asp/Glu, deprotonated Lys.
    "HID", "HIE", "HIP", "CYX", "CYM", "ASH", "GLH", "LYN",
    // GROMACS (GROMOS and OPLS-AA databases).
    "HISA", "HISB", "HISH", "HISD", "HISE", "HISP", "HIS1", "HIS2", "CYSH", "LYSH", "ARGN", "ASPH",
    "GLUH", "ASP1", "GLU1",
    // GLYCAM glycoprotein-linkage residues: the amino acid a glycan is
    // attached to (the sugars themselves are in `crate::glycan`).
    "NLN", "OLS", "OLT", "ZOLS", "ZOLT",
];

/// Modified and non-standard residues that the wwPDB Chemical Component
/// Dictionary (Westbrook et al. 2015, Nucleic Acids Res 43:D364) lists as
/// peptide-linking, plus D-amino acids.
const PROTEIN_MODIFIED: &[&str] = &[
    "MSE", "HYP", "SEP", "TPO", "PTR", "CSO", "CSD", "OCS", "CSS", "CME", "SMC", "CAS", "CSX",
    "MLY", "M3L", "MLZ", "ALY", "KCX", "LLP", "PCA", "TYS", "HIC", "ABA", "NLE", "ORN", "FME",
    "AIB", "MHO", "SAC", "SNN", "TPQ", "OMT", "LYZ", "NEP", "MEN", "BMT", "DAL", "DAR", "DAS",
    "DCY", "DGL", "DHI", "DIL", "DLE", "DLY", "DPN", "DPR", "DSG", "DSN", "DTH", "DTR", "DTY",
    "DVA", "MED", "MVA", "STA", "CGU", "TRQ", "TYQ",
];

// ---- nucleic acids ---------------------------------------------------

/// PDB nucleotide codes: `A C G U I N` (RNA/any), `D`-prefixed DNA
/// (wwPDB CCD), and modified bases from the CCD.
const NUCLEOTIDES: &[&str] = &[
    "A", "C", "G", "U", "I", "N", "T", "DA", "DC", "DG", "DT", "DI", "DU", "DN", "PSU", "5MC",
    "7MG", "1MA", "2MG", "M2G", "OMC", "OMG", "H2U", "5MU", "4SU", "5BU", "BRU", "CBR", "6MA",
    "MA6", "OMU", "UR3", "1MG", "5CM", "A2M", "G7M", "8OG", "5IU", "CFZ",
];

/// Amber (`leaprc.DNA.*`, `leaprc.RNA.*`) and CHARMM/GROMACS nucleotide
/// variants: 5'/3'-terminal templates, `R`-prefixed RNA, and the
/// three-letter base names of the CHARMM27 nucleic-acid set (Foloppe &
/// MacKerell 2000, J Comput Chem 21:86).
const NUCLEOTIDE_VARIANTS: &[&str] = &[
    "DA5", "DA3", "DAN", "DC5", "DC3", "DCN", "DG5", "DG3", "DT5", "DT3", "DTN", "A5", "A3", "AN",
    "C5", "C3", "CN", "G5", "G3", "GN", "U5", "U3", "UN", "RA", "RC", "RG", "RU", "RA5", "RA3",
    "RAN", "RC5", "RC3", "RCN", "RG5", "RG3", "RGN", "RU5", "RU3", "RUN", "ADE", "CYT", "GUA",
    "THY", "URA", "DADE", "DCYT", "DGUA", "DTHY", "RADE", "RCYT", "RGUA", "RURA",
];

// ---- water -----------------------------------------------------------

/// Water models: SPC (Berendsen et al. 1981), SPC/E (Berendsen, Grigera &
/// Straatsma 1987, J Phys Chem 91:6269), TIP3P/TIP4P (Jorgensen et al.
/// 1983, J Chem Phys 79:926), TIP5P (Mahoney & Jorgensen 2000), OPC
/// (Izadi, Anandakrishnan & Onufriev 2014), Drude SWM4 (Lamoureux et al.
/// 2006); plus PDB and package spellings (`HOH`, `WAT`, `SOL`, `T3P`...).
const WATER: &[&str] = &[
    "HOH", "WAT", "H2O", "DOD", "D2O", "SOL", "TIP", "TIP3", "TIP3P", "TP3", "TP3M", "TIP4",
    "TIP4P", "TIP4PEW", "TP4", "TIP5", "TIP5P", "TP5", "SPC", "SPCE", "OPC", "OPC3", "T3P", "T4P",
    "T5P", "SWM4", "SWM6", "WATER",
    // Martini coarse-grained water beads (Marrink et al. 2007, J Phys Chem
    // B 111:7812): standard, antifreeze, polarizable.
    "W", "WN", "PW",
];

// ---- lipids ----------------------------------------------------------

/// CHARMM36 lipids (Klauda et al. 2010, J Phys Chem B 114:7830) and the
/// extended CHARMM36 set: phosphatidylcholines,
/// -ethanolamines, -serines, -glycerols, -acids, lyso forms, inositol
/// lipids, sterols, sphingolipids and cardiolipins.
const LIPIDS_CHARMM: &[&str] = &[
    "DLPC", "DMPC", "DPPC", "DSPC", "DOPC", "DEPC", "DAPC", "DGPC", "DUPC", "DIPC", "DBPC", "DVPC",
    "DNPC", "POPC", "PLPC", "PAPC", "SOPC", "SDPC", "PDPC", "PYPC", "LPPC", "LMPC", "LSPC", "DLPE",
    "DMPE", "DPPE", "DSPE", "DOPE", "DEPE", "DAPE", "POPE", "PLPE", "PAPE", "SOPE", "SDPE", "DLPS",
    "DMPS", "DPPS", "DSPS", "DOPS", "POPS", "SOPS", "DLPG", "DMPG", "DPPG", "DSPG", "DOPG", "POPG",
    "PLPG", "DLPA", "DMPA", "DPPA", "DOPA", "POPA",
    // Inositol phospholipids and PIP2/PIP3.
    "POPI", "SAPI", "SAPI2", "SAPI24", "SAPI25", "PIP2", "PIP3", "PI2A", "PI3P",
    // Cholesterol and other sterols (Lim et al. 2012, J Phys Chem B 116:203).
    "CHL1", "CHOL", "ERG", "ERGO", "CHSD", "STIG", "SITO",
    // Sphingomyelins and ceramides (Venable et al. 2014, Biophys J 107:134).
    "PSM", "SSM", "NSM", "LSM", "ESM", "CER", "CER160", "CER180", "CER181", "CER240", "CER241",
    // Cardiolipins.
    "TOCL", "TLCL", "TMCL", "TYCL", "PVCL", "CDL1", "CDL2",
];

/// Amber Lipid21 (Dickson et al. 2022, J Chem Theory Comput 18:1726) and
/// Lipid17: each lipid is split into a head-group residue plus one residue
/// per acyl tail, so every piece must classify. `PA` names both the
/// palmitoyl tail and phosphatidic acid (both lipid, so no conflict);
/// `LA`, `AR`, `PA` are also element symbols and are lipid only for
/// multi-atom residues (see [`MULTI_ATOM_ONLY`]).
const LIPIDS_AMBER: &[&str] = &[
    // Head groups: phosphatidylcholine, -ethanolamine, -serine, -glycerol
    // (charged and neutral), phosphatidic acid, cholesterol.
    "PC", "PE", "PS", "PGR", "PGS", "PH-", "PA", "CHL",
    // Acyl tails: palmitoyl (PA above), oleoyl, stearoyl, lauroyl,
    // myristoyl, arachidonoyl, linoleoyl (docosahexaenoyl `DHA` is below).
    "OL", "ST", "LA", "MY", "AR", "LIN",
];

/// Martini lipids (Marrink et al. 2007, J Phys Chem B 111:7812; Marrink &
/// Tieleman 2013, Chem Soc Rev 42:6801): the letters match CHARMM's for
/// the common lipids; these are the coarse-grained-only names. The same
/// spellings serve Slipids (Jambeck & Lyubartsev 2012, J Phys Chem B
/// 116:3164), GROMOS 54A7 (Poger, van Gunsteren & Mark 2010, J Comput
/// Chem 31:1117) and Berger (Berger et al. 1997, Biophys J 72:2002)
/// lipids, which are covered by the CHARMM group above.
const LIPIDS_MARTINI: &[&str] = &["DPSM", "DXCE", "DPCE", "CARD", "CDL0"];

/// wwPDB Chemical Component Dictionary lipid, fatty-acid, sterol and
/// detergent codes seen in membrane-protein co-crystals: the alkyl
/// detergents stand in for the membrane.
const LIPIDS_CCD: &[&str] = &[
    // Phospholipids.
    "PEE", "PC1", "PGT", "PGV", "PTY", "3PE", "POV", "LHG", "LPE", "EPH", "CDL",
    // Glycolipids.
    "LMG", "DGD", "SQD", // Fatty acids and monoacylglycerols.
    "OLC", "OLA", "OLB", "PLM", "MYR", "STE", // Sterols and bile acids.
    "CLR", "Y01", "CHD", // Alkyl detergents.
    "LMT", "LDA", "BOG", "DMU", "UDM", "HTG", "BNG", "C8E", "SDS", // Alkanes.
    "D10", "D12", "OCT", "HEX", "UND",
];

/// Structural placeholders that carry no chemistry.
const OTHER: &[&str] = &["DUM", "DUMM", "LP", "EP", "MW"];

pub(super) const GROUPS: &[(ResidueClass, &[&[&str]])] = &[
    (
        ResidueClass::Protein,
        &[AMINO_ACIDS, PROTEIN_VARIANTS, PROTEIN_MODIFIED],
    ),
    (ResidueClass::Nucleic, &[NUCLEOTIDES, NUCLEOTIDE_VARIANTS]),
    (ResidueClass::Water, &[WATER]),
    (
        ResidueClass::Lipid,
        &[LIPIDS_CHARMM, LIPIDS_AMBER, LIPIDS_MARTINI, LIPIDS_CCD],
    ),
    (ResidueClass::Other, &[OTHER]),
];

/// Lipid names that are also element symbols (lanthanum, argon,
/// protactinium): lipid only when the residue has more than one atom.
pub(super) const MULTI_ATOM_ONLY: &[&str] = &["LA", "AR", "PA"];

/// Lipid names the glycan table also uses. Lipid21's docosahexaenoyl tail
/// `DHA` (about 54 atoms) shares its code with a 24-atom ketoaldonic sugar
/// acid; the lipid reading wins from [`LIPID_OVER_GLYCAN_MIN_ATOMS`] atoms.
pub(super) const LIPID_OVER_GLYCAN: &[&str] = &["DHA"];
pub(super) const LIPID_OVER_GLYCAN_MIN_ATOMS: usize = 30;
