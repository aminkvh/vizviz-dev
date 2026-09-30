//! Legacy PDB reader, written from the wwPDB format description (v3.3).
//! Fixed columns, hybrid-36 serial/residue numbers, MODEL/ENDMDL frames,
//! CONECT and LINK explicit bonds, HELIX/SHEET secondary structure.

use std::collections::HashMap;

use vv_core::glam::Vec3;
use vv_core::{
    AnnotationCategory, Annotations, AtomExtra, AtomRow, Element, ExplicitBond, ExplicitBondKind,
    SecondaryStructure, Structure, Topology, TopologyBuilder,
};

use crate::float::parse_f32;
use crate::pdb_seqres::Seqres;
use crate::ss_range::{self, SsRange};
use crate::ParseError;

/// End of the z coordinate (column 54): the shortest record we accept.
const COORD_END: usize = 54;

fn field(line: &[u8], start: usize, end: usize) -> &[u8] {
    let end = end.min(line.len());
    if start >= end {
        return b"";
    }
    let f = &line[start..end];
    let s = f.iter().position(|b| *b != b' ').unwrap_or(f.len());
    let e = f.iter().rposition(|b| *b != b' ').map_or(s, |p| p + 1);
    &f[s..e]
}

/// Decodes a hybrid-36 number of the given width (plain decimal first).
pub fn hybrid36(s: &[u8], width: usize) -> Option<i32> {
    if s.is_empty() {
        return None;
    }
    if let Some(v) = crate::float::parse_i32(s) {
        return Some(v);
    }
    if s.len() != width {
        return None;
    }
    let digit = |c: u8| -> Option<i32> {
        match c {
            b'0'..=b'9' => Some((c - b'0') as i32),
            b'A'..=b'Z' => Some((c - b'A') as i32 + 10),
            b'a'..=b'z' => Some((c - b'a') as i32 + 10),
            _ => None,
        }
    };
    let upper = s[0].is_ascii_uppercase();
    let lower = s[0].is_ascii_lowercase();
    if !upper && !lower {
        return None;
    }
    let mut v = 0i32;
    for &c in s {
        v = v.checked_mul(36)?.checked_add(digit(c)?)?;
    }
    // "A000.." starts right after the decimal range; "a000.." starts after
    // the 26 upper-case blocks.
    let base = 10i32.pow(width as u32);
    let a_block = 10 * 36i32.pow(width as u32 - 1);
    let upper_span = 26 * 36i32.pow(width as u32 - 1);
    Some(if upper {
        v - a_block + base
    } else {
        v - a_block + base + upper_span
    })
}

/// Inverse of [`hybrid36`]: the plain decimal string when `n` fits in
/// `width` characters, else the hybrid-36 extension (uppercase for the
/// first overflow block, lowercase for the second); `None` past that
/// (rare: needs `n >= 10^width + 2 * 26 * 36^(width-1)`).
pub fn encode_hybrid36(n: i64, width: usize) -> Option<String> {
    let plain = n.to_string();
    if plain.len() <= width {
        return Some(plain);
    }
    if n < 0 || width == 0 {
        return None;
    }
    let base = 10i64.checked_pow(width as u32)?;
    let a_block = 10i64 * 36i64.pow(width as u32 - 1);
    let upper_span = 26i64 * 36i64.pow(width as u32 - 1);
    let (v, lower) = if n < base + upper_span {
        (n - base + a_block, false)
    } else if n < base + 2 * upper_span {
        (n - base - upper_span + a_block, true)
    } else {
        return None;
    };
    let mut v = v as u64;
    let mut buf = vec![0u8; width];
    for slot in buf.iter_mut().rev() {
        let digit = (v % 36) as u8;
        *slot = if digit < 10 {
            b'0' + digit
        } else {
            b'A' + (digit - 10)
        };
        v /= 36;
    }
    if lower {
        buf.iter_mut().for_each(u8::make_ascii_lowercase);
    }
    String::from_utf8(buf).ok()
}

fn element_of(line: &[u8], name: &[u8]) -> Element {
    let symbol = field(line, 76, 78);
    if !symbol.is_empty() {
        let e = Element::from_symbol(symbol);
        if !e.is_unknown() {
            return e;
        }
    }
    // Name columns 13-16: a two-letter element is left-justified at 13,
    // e.g. "FE  " or "ZN1 ". A one-letter element's name is padded with a
    // leading space instead (" CA " is an alpha carbon, "CA  " a calcium
    // ion) - except hydrogen's, which often needs all 4 columns for a
    // branched locant ("HG11", "HE21": Val/Ile's and Gln's own hydrogens,
    // not mercury or helium) and so cannot be told apart by padding alone;
    // columns 3-4 being letters, not the locant's digits, is the tell.
    let two = if name.len() == 4 && name[0].is_ascii_alphabetic() && name[1].is_ascii_alphabetic() {
        Element::from_symbol(&name[..2])
    } else {
        Element::UNKNOWN
    };
    let is_he_or_hg = two == Element::from_symbol(b"He") || two == Element::from_symbol(b"Hg");
    let has_locant = name[2] != b' ' || name[3] != b' ';
    let hydrogen_with_locant = is_he_or_hg && has_locant;
    if !two.is_unknown() && !hydrogen_with_locant {
        return two;
    }
    name.iter()
        .find(|b| b.is_ascii_alphabetic())
        .map_or(Element::UNKNOWN, |&b| Element::from_symbol(&[b]))
}

/// Deuterium from the element column, else (column blank) from the name
/// field. Only alignments no heavy atom uses count: a lone `D` right-justified
/// at column 14 (` D  `, ` DA `), or a four-character hydrogen locant
/// (`DD21`, `DG12`). A two-letter element at column 13 (`DY  `, dysprosium)
/// is never read as deuterium.
fn is_deuterium(symbol: &[u8], name: &[u8; 4]) -> bool {
    match symbol {
        [] => {
            let one_letter = name[0] == b' ' && name[1] == b'D';
            let locant = name[0] == b'D'
                && name[1].is_ascii_alphabetic()
                && name[2].is_ascii_digit()
                && name[3].is_ascii_digit();
            one_letter || locant
        }
        s => matches!(s, b"D" | b"d"),
    }
}

/// A LINK record names its two atoms by (chain, residue number, insertion
/// code, atom name) rather than serial number, so resolving one needs a
/// lookup built alongside `serial_to_atom` while atoms are read.
type AtomKey = (u8, i32, u8, [u8; 4]);

fn atom_key(
    line: &[u8],
    name_start: usize,
    chain_col: usize,
    seq_start: usize,
    icode_col: usize,
) -> Option<AtomKey> {
    let name_raw = &line[name_start.min(line.len())..(name_start + 4).min(line.len())];
    let mut name = [b' '; 4];
    name[..name_raw.len()].copy_from_slice(name_raw);
    let chain = field(line, chain_col, chain_col + 1)
        .first()
        .copied()
        .unwrap_or(b' ');
    let seq = hybrid36(field(line, seq_start, seq_start + 4), 4)?;
    let icode = field(line, icode_col, icode_col + 1)
        .first()
        .copied()
        .unwrap_or(0);
    Some((chain, seq, icode, name))
}

/// The two atoms a LINK record names (wwPDB v3.3 columns 13-27 and 43-57).
fn link_atoms(line: &[u8]) -> Option<(AtomKey, AtomKey)> {
    Some((
        atom_key(line, 12, 21, 22, 26)?,
        atom_key(line, 42, 51, 52, 56)?,
    ))
}

/// Columns 79-80: `2+` per the format, and the sign-first `+2` some
/// writers emit. A bare sign is +-1.
/// SSBOND (v3.3 columns 16-22 and 30-36) names residues, not atoms; a
/// disulfide is always between the two `SG` atoms.
fn ssbond_atoms(line: &[u8]) -> Option<(AtomKey, AtomKey)> {
    let key = |chain_col: usize, seq_col: usize| {
        let chain = field(line, chain_col, chain_col + 1).first().copied();
        let seq = hybrid36(field(line, seq_col, seq_col + 4), 4)?;
        let icode = field(line, seq_col + 4, seq_col + 5).first().copied();
        Some((chain.unwrap_or(b' '), seq, icode.unwrap_or(0), *b" SG "))
    };
    Some((key(15, 17)?, key(29, 31)?))
}

/// `chain_col` and `seq_col` locate the initial residue (its insertion code
/// follows the 4-column number); the terminal residue is at columns 34-38.
fn ss_range(line: &[u8], chain_col: usize, seq_col: usize) -> Option<SsRange> {
    let residue = |seq_col: usize| {
        let seq = hybrid36(field(line, seq_col, seq_col + 4), 4)?;
        let icode = field(line, seq_col + 4, seq_col + 5).first().copied();
        Some((seq, icode.unwrap_or(0)))
    };
    Some(SsRange {
        chain: text(line, chain_col, chain_col + 1).to_string(),
        first: residue(seq_col)?,
        last: residue(33)?,
    })
}

fn charge_of(line: &[u8]) -> i8 {
    let (digit, sign) = match field(line, 78, 80) {
        [d, s] if d.is_ascii_digit() => (*d, *s),
        [s, d] if d.is_ascii_digit() => (*d, *s),
        [s] => (b'1', *s),
        _ => return 0,
    };
    let n = (digit - b'0') as i8;
    match sign {
        b'+' => n,
        b'-' => -n,
        _ => 0,
    }
}

/// Chain name: column 22, else the legacy segment id (columns 73-76) that
/// MD tools write instead when they run out of one-letter chains.
fn chain_name(line: &[u8]) -> &str {
    let chain = text(line, 21, 22);
    if chain.is_empty() {
        text(line, 72, 76)
    } else {
        chain
    }
}

/// Header records collected while scanning, then mapped onto mmCIF
/// category/item names so consumers see one vocabulary for both formats.
#[derive(Default)]
struct PdbHeader {
    classification: String,
    deposition_date: String,
    id: String,
    keywords: String,
    method: String,
    /// The resolution token as written (`"1.74"`), not a reformatted float.
    resolution: String,
    journal_title: String,
    doi: String,
    compnd: String,
    source: String,
    /// (database, accession) from DBREF.
    dbref: Vec<(String, String)>,
    /// `a b c alpha beta gamma` from CRYST1, as written.
    cell: Vec<String>,
    space_group: String,
    z_pdb: String,
    /// Inside the `REMARK 465` residue table (after its column header).
    in_missing_table: bool,
    /// `(residue name, chain, sequence number, insertion code)` of every
    /// `REMARK 465` residue, model 1 only.
    missing: Vec<(String, String, String, String)>,
}

fn text(line: &[u8], start: usize, end: usize) -> &str {
    std::str::from_utf8(field(line, start, end)).unwrap_or("")
}

/// Appends a continuation line, separating it from the previous one.
fn join(into: &mut String, part: &str) {
    if part.is_empty() {
        return;
    }
    if !into.is_empty() {
        into.push(' ');
    }
    into.push_str(part);
}

impl PdbHeader {
    fn record(&mut self, record: &[u8], line: &[u8]) {
        match record {
            b"HEADER" => {
                self.classification = text(line, 10, 50).to_string();
                self.deposition_date = text(line, 50, 59).to_string();
                self.id = text(line, 62, 66).to_string();
            }
            b"CRYST1" => {
                let cols = [(6, 15), (15, 24), (24, 33), (33, 40), (40, 47), (47, 54)];
                self.cell = cols
                    .iter()
                    .map(|&(a, b)| text(line, a, b).to_string())
                    .collect();
                self.space_group = text(line, 55, 66).to_string();
                self.z_pdb = text(line, 66, 70).to_string();
            }
            b"COMPND" => join(&mut self.compnd, text(line, 10, 80)),
            b"SOURCE" => join(&mut self.source, text(line, 10, 80)),
            b"KEYWDS" => join(&mut self.keywords, text(line, 10, 80)),
            b"EXPDTA" => join(&mut self.method, text(line, 10, 80)),
            b"REMARK" if field(line, 7, 10) == b"2" => {
                // `REMARK   2 RESOLUTION.    1.74 ANGSTROMS.`; NMR entries say
                // `NOT APPLICABLE`, which parse_f32 rejects.
                let rest = text(line, 11, 80);
                if let Some(after) = rest.split_once("RESOLUTION.") {
                    if let Some(tok) = after.1.split_whitespace().next() {
                        if parse_f32(tok.as_bytes()).is_some() {
                            self.resolution = tok.to_string();
                        }
                    }
                }
            }
            b"REMARK" if field(line, 7, 10) == b"465" => self.missing_residue(line),
            b"JRNL  " => match field(line, 12, 16) {
                b"TITL" => join(&mut self.journal_title, text(line, 19, 80)),
                b"DOI" => join(&mut self.doi, text(line, 19, 80)),
                _ => {}
            },
            b"DBREF " => {
                let db = text(line, 26, 32);
                let acc = text(line, 33, 41);
                if !db.is_empty() && !acc.is_empty() {
                    self.dbref.push((db.to_string(), acc.to_string()));
                }
            }
            _ => {}
        }
    }

    /// One `REMARK 465` line: the table header switches row parsing on;
    /// a row is residue name (cols 16-18), chain (20), number (22-26) and
    /// insertion code (27). Rows of models after the first are skipped.
    fn missing_residue(&mut self, line: &[u8]) {
        let words: Vec<&str> = text(line, 11, 80).split_whitespace().collect();
        if words == ["M", "RES", "C", "SSSEQI"] {
            self.in_missing_table = true;
            return;
        }
        let model = text(line, 10, 14);
        if !self.in_missing_table || !(model.is_empty() || model == "1") {
            return;
        }
        let (name, chain) = (text(line, 15, 18), text(line, 19, 20));
        let (seq, ins) = (text(line, 21, 26), text(line, 26, 27));
        if !name.is_empty() && !seq.is_empty() {
            self.missing
                .push((name.into(), chain.into(), seq.into(), ins.into()));
        }
    }

    fn into_annotations(self, title: &str) -> Annotations {
        let mut out = Annotations::default();
        let mut single = |name: &str, items: &[(&str, &str)]| {
            let present: Vec<(&str, &str)> = items
                .iter()
                .copied()
                .filter(|(_, v)| !v.is_empty())
                .collect();
            if !present.is_empty() {
                out.categories.push(AnnotationCategory {
                    name: name.to_string(),
                    items: present.iter().map(|(k, _)| k.to_string()).collect(),
                    rows: vec![present.iter().map(|(_, v)| v.to_string()).collect()],
                });
            }
        };
        single("entry", &[("id", &self.id)]);
        single("struct", &[("title", title)]);
        single(
            "struct_keywords",
            &[
                ("pdbx_keywords", &self.classification),
                ("text", &self.keywords),
            ],
        );
        single("exptl", &[("method", &self.method)]);
        single("refine", &[("ls_d_res_high", &self.resolution)]);
        single(
            "citation",
            &[
                ("title", &self.journal_title),
                ("pdbx_database_id_DOI", &self.doi),
            ],
        );
        single(
            "pdbx_database_status",
            &[("recvd_initial_deposition_date", &self.deposition_date)],
        );

        if self.cell.iter().all(|v| !v.is_empty()) && self.cell.len() == 6 {
            let c = &self.cell;
            single(
                "cell",
                &[
                    ("length_a", &c[0]),
                    ("length_b", &c[1]),
                    ("length_c", &c[2]),
                    ("angle_alpha", &c[3]),
                    ("angle_beta", &c[4]),
                    ("angle_gamma", &c[5]),
                    ("Z_PDB", &self.z_pdb),
                ],
            );
            single("symmetry", &[("space_group_name_H-M", &self.space_group)]);
        }

        let entities = spec_per_mol_id(&self.compnd, "MOLECULE");
        if !entities.is_empty() {
            out.categories.push(AnnotationCategory {
                name: "entity".to_string(),
                items: vec!["id".to_string(), "pdbx_description".to_string()],
                rows: entities.into_iter().map(|(id, v)| vec![id, v]).collect(),
            });
        }
        let strands = spec_per_mol_id(&self.compnd, "CHAIN");
        if !strands.is_empty() {
            out.categories.push(AnnotationCategory {
                name: "entity_poly".to_string(),
                items: vec!["entity_id".to_string(), "pdbx_strand_id".to_string()],
                rows: strands
                    .into_iter()
                    .map(|(id, chains)| vec![id, chains.replace(' ', "")])
                    .collect(),
            });
        }
        let organisms = spec_per_mol_id(&self.source, "ORGANISM_SCIENTIFIC");
        if !organisms.is_empty() {
            out.categories.push(AnnotationCategory {
                name: "entity_src_gen".to_string(),
                items: vec![
                    "entity_id".to_string(),
                    "pdbx_gene_src_scientific_name".to_string(),
                ],
                rows: organisms.into_iter().map(|(id, v)| vec![id, v]).collect(),
            });
        }
        if !self.missing.is_empty() {
            out.categories
                .push(missing_residues_category(&self.missing));
        }
        let unp: Vec<Vec<String>> = self
            .dbref
            .into_iter()
            .filter(|(db, _)| db == "UNP")
            .map(|(db, acc)| vec![db, acc])
            .collect();
        if !unp.is_empty() {
            out.categories.push(AnnotationCategory {
                name: "struct_ref".to_string(),
                items: vec!["db_name".to_string(), "pdbx_db_accession".to_string()],
                rows: unp,
            });
        }
        out
    }
}

/// `REMARK 465` rows as the mmCIF `pdbx_unobs_or_zero_occ_residues`
/// category (all unobserved: `occupancy_flag` 0).
fn missing_residues_category(rows: &[(String, String, String, String)]) -> AnnotationCategory {
    let items = [
        "polymer_flag",
        "occupancy_flag",
        "auth_asym_id",
        "auth_comp_id",
        "auth_seq_id",
        "PDB_ins_code",
    ];
    AnnotationCategory {
        name: "pdbx_unobs_or_zero_occ_residues".to_string(),
        items: items.map(String::from).to_vec(),
        rows: rows
            .iter()
            .map(|(name, chain, seq, ins)| {
                ["Y", "0", chain, name, seq, ins].map(String::from).to_vec()
            })
            .collect(),
    }
}

/// Reads a COMPND/SOURCE specification list (`KEY: value;` tokens) and
/// returns `(MOL_ID, value of key)` for every MOL_ID block that has the key.
fn spec_per_mol_id(spec: &str, key: &str) -> Vec<(String, String)> {
    let mut out = Vec::new();
    let mut current: Option<String> = None;
    for token in spec.split(';') {
        let Some((k, v)) = token.split_once(':') else {
            continue;
        };
        let (k, v) = (k.trim(), v.trim());
        if k == "MOL_ID" {
            current = Some(v.to_string());
        } else if k == key && !v.is_empty() {
            if let Some(id) = &current {
                out.push((id.clone(), v.to_string()));
            }
        }
    }
    out
}

pub fn parse(src: &[u8]) -> Result<Structure, ParseError> {
    let mut builder = TopologyBuilder::with_capacity(src.len() / 81 + 1);
    let mut header = PdbHeader::default();
    let mut extra_models: Vec<Vec<Vec3>> = Vec::new();
    let mut in_first_model = true;
    let mut seen_model = false;
    let mut serial_to_atom: HashMap<u32, u32> = HashMap::new();
    let mut link_lookup: HashMap<AtomKey, u32> = HashMap::new();
    let mut conect: Vec<(u32, u32)> = Vec::new();
    let mut link: Vec<(AtomKey, AtomKey)> = Vec::new();
    let mut ssbond: Vec<(AtomKey, AtomKey)> = Vec::new();
    let mut helices: Vec<SsRange> = Vec::new();
    let mut strands: Vec<SsRange> = Vec::new();
    let mut title = String::new();
    let mut id = String::new();
    // Blocks after the first `END` (frames separated by END alone).
    let mut after_end = false;
    let mut block: Vec<(Vec3, [u8; 4])> = Vec::new();
    let mut seqres = Seqres::default();
    let mut last_segid: &[u8] = b"";

    for (line_no, raw) in src.split(|&b| b == b'\n').enumerate() {
        let line = match raw.last() {
            Some(b'\r') => &raw[..raw.len() - 1],
            _ => raw,
        };
        if line.len() < 3 {
            continue;
        }
        let mut padded = *b"      ";
        let head = &line[..line.len().min(6)];
        padded[..head.len()].copy_from_slice(head);
        let record = &padded[..];
        if after_end && !matches!(record, b"ATOM  " | b"HETATM" | b"END   ") {
            continue;
        }
        match record {
            b"ATOM  " | b"HETATM" => {
                if line.len() < COORD_END {
                    return Err(ParseError::Malformed {
                        line: line_no + 1,
                        message: "atom record ends before the z coordinate (column 54)".into(),
                    });
                }
                let mut name = [b' '; 4];
                name.copy_from_slice(&line[12..16]);
                let coord = |a, b| parse_f32(field(line, a, b));
                let (Some(x), Some(y), Some(z)) = (coord(30, 38), coord(38, 46), coord(46, 54))
                else {
                    return Err(ParseError::Malformed {
                        line: line_no + 1,
                        message: "unreadable coordinates".into(),
                    });
                };
                let position = Vec3::new(x, y, z);
                if after_end {
                    block.push((position, name));
                    continue;
                }
                if !in_first_model {
                    extra_models.last_mut().unwrap().push(position);
                    continue;
                }
                let serial = hybrid36(field(line, 6, 11), 5).unwrap_or(0).max(0) as u32;
                let seq = hybrid36(field(line, 22, 26), 4).unwrap_or(0);
                let chain = chain_name(line);
                // Columns 18-21: MD packages use column 21 for a fourth
                // residue-name character; the format leaves it blank.
                let comp = text(line, 17, 21);
                let segid = field(line, 72, 76);
                if segid != last_segid {
                    builder.end_chain();
                    last_segid = segid;
                }
                let extra = AtomExtra {
                    segid: std::str::from_utf8(segid).unwrap_or(""),
                    deuterium: is_deuterium(field(line, 76, 78), &name),
                    ..Default::default()
                };
                let row = AtomRow {
                    element: element_of(line, &name),
                    name,
                    serial,
                    alt_loc: field(line, 16, 17).first().copied().unwrap_or(0),
                    comp,
                    asym: chain,
                    auth_asym: chain,
                    seq_id: seq,
                    auth_seq_id: seq,
                    ins_code: field(line, 26, 27).first().copied().unwrap_or(0),
                    entity: 0,
                    position,
                    occupancy: coord(54, 60).unwrap_or(1.0),
                    b_factor: coord(60, 66).unwrap_or(0.0),
                    charge: charge_of(line),
                    hetero: record == b"HETATM",
                };
                serial_to_atom
                    .entry(serial)
                    .or_insert(builder.atom_count() as u32);
                if let Some(key) = atom_key(line, 12, 21, 22, 26) {
                    link_lookup
                        .entry(key)
                        .or_insert(builder.atom_count() as u32);
                }
                builder.push_with(&row, &extra);
            }
            b"SEQRES" => seqres.record(line),
            b"MODEL " => {
                if seen_model {
                    in_first_model = false;
                    extra_models.push(Vec::new());
                }
                seen_model = true;
            }
            b"ENDMDL" => {}
            b"TER   " => builder.end_chain(),
            b"END   " => {
                if after_end && !accept_block(&builder.topology.name, &mut block, &mut extra_models)
                {
                    break;
                }
                after_end = true;
            }
            b"SSBOND" => {
                if let Some(pair) = ssbond_atoms(line) {
                    ssbond.push(pair);
                }
            }
            b"CONECT" => {
                let Some(a) = hybrid36(field(line, 6, 11), 5) else {
                    continue;
                };
                for start in [11, 16, 21, 26] {
                    if let Some(b) = hybrid36(field(line, start, start + 5), 5) {
                        conect.push((a as u32, b as u32));
                    }
                }
            }
            b"LINK  " => {
                if let Some(pair) = link_atoms(line) {
                    link.push(pair);
                }
            }
            b"HELIX " => {
                if let Some(r) = ss_range(line, 19, 21) {
                    helices.push(r);
                }
            }
            b"SHEET " => {
                if let Some(r) = ss_range(line, 21, 22) {
                    strands.push(r);
                }
            }
            b"TITLE " => {
                let t = std::str::from_utf8(field(line, 10, 80)).unwrap_or("");
                if !title.is_empty() {
                    title.push(' ');
                }
                title.push_str(t);
            }
            b"HEADER" => {
                id = std::str::from_utf8(field(line, 62, 66))
                    .unwrap_or("")
                    .to_string();
                header.record(record, line);
            }
            _ => header.record(record, line),
        }
    }

    accept_block(&builder.topology.name, &mut block, &mut extra_models);
    if builder.atom_count() == 0 {
        return Err(ParseError::NoAtoms);
    }
    let atom_count = builder.atom_count();
    let topology = &mut builder.topology;
    topology.polymer_hint = seqres.hints(topology);
    topology.annotations = header.into_annotations(&title);
    topology.title = title;
    topology.id = id;
    // The old convention for a double/triple bond: an atom's own CONECT
    // line repeats the partner's serial once per bond order. Count each
    // directed (from atom's line) occurrence, then take each normalized
    // pair's higher side -- many files only encode the repeat on one
    // atom's line, not both.
    let mut directed: HashMap<(u32, u32), u32> = HashMap::new();
    for (a, b) in &conect {
        if let (Some(&i), Some(&j)) = (serial_to_atom.get(a), serial_to_atom.get(b)) {
            if i != j {
                *directed.entry((i, j)).or_insert(0) += 1;
            }
        }
    }
    let mut repeats: HashMap<[u32; 2], u32> = HashMap::new();
    for (&(i, j), &count) in &directed {
        let pair = if i < j { [i, j] } else { [j, i] };
        let entry = repeats.entry(pair).or_insert(0);
        *entry = (*entry).max(count);
    }
    // Sorted rather than left in `HashMap` iteration order, so a written
    // file's CONECT lines come out the same on every run.
    let mut repeats: Vec<([u32; 2], u32)> = repeats.into_iter().collect();
    repeats.sort_unstable_by_key(|&(pair, _)| pair);
    for (pair, count) in repeats {
        let order = match count {
            0 | 1 => vv_core::BondOrder::Single,
            2 => vv_core::BondOrder::Double,
            _ => vv_core::BondOrder::Triple,
        };
        topology.explicit_bonds.push(ExplicitBond::with_order(
            pair[0],
            pair[1],
            ExplicitBondKind::Covalent,
            order,
        ));
    }
    let resolve = |pairs: &[(AtomKey, AtomKey)]| -> Vec<[u32; 2]> {
        pairs
            .iter()
            .filter_map(|(a, b)| Some([*link_lookup.get(a)?, *link_lookup.get(b)?]))
            .filter(|[i, j]| i != j)
            .collect()
    };
    for [i, j] in resolve(&link) {
        let metal =
            topology.element[i as usize].is_metal() || topology.element[j as usize].is_metal();
        let kind = if metal {
            ExplicitBondKind::Metal
        } else {
            ExplicitBondKind::Covalent
        };
        add_bond(topology, i, j, kind);
    }
    for [i, j] in resolve(&ssbond) {
        add_bond(topology, i, j, ExplicitBondKind::Disulfide);
    }
    ss_range::apply(topology, &helices, SecondaryStructure::Helix);
    ss_range::apply(topology, &strands, SecondaryStructure::Strand);

    let extra: Vec<Vec<Vec3>> = extra_models
        .into_iter()
        .filter(|m| m.len() == atom_count)
        .collect();
    builder.finish_with_frames(extra).map_err(ParseError::from)
}

/// Adds a LINK/SSBOND bond, upgrading the kind of a CONECT bond between the
/// same atoms instead of duplicating it.
fn add_bond(topology: &mut Topology, i: u32, j: u32, kind: ExplicitBondKind) {
    let same = |b: &ExplicitBond| b.atoms == [i, j] || b.atoms == [j, i];
    match topology.explicit_bonds.iter_mut().find(|b| same(b)) {
        Some(existing) => existing.kind = kind,
        None => topology.explicit_bonds.push(ExplicitBond::new(i, j, kind)),
    }
}

/// Takes the atoms collected after an `END` as another frame when their
/// names match the first block's, one for one (MD output that separates
/// frames with `END` alone). An empty block is fine; a mismatch means a
/// second, concatenated entry, so the caller stops reading.
fn accept_block(
    names: &[[u8; 4]],
    block: &mut Vec<(Vec3, [u8; 4])>,
    frames: &mut Vec<Vec<Vec3>>,
) -> bool {
    let taken = std::mem::take(block);
    if taken.is_empty() {
        return true;
    }
    let matches = taken.len() == names.len() && taken.iter().zip(names).all(|(a, n)| a.1 == *n);
    if matches {
        frames.push(taken.into_iter().map(|(p, _)| p).collect());
    }
    matches
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hybrid36_decodes_decimal_and_extended_ranges() {
        assert_eq!(hybrid36(b"12345", 5), Some(12345));
        assert_eq!(hybrid36(b"99999", 5), Some(99999));
        assert_eq!(hybrid36(b"A0000", 5), Some(100000));
        assert_eq!(hybrid36(b"A0001", 5), Some(100001));
        assert_eq!(hybrid36(b"a0000", 5), Some(100000 + 26 * 36i32.pow(4)));
        assert_eq!(hybrid36(b"9999", 4), Some(9999));
        assert_eq!(hybrid36(b"A000", 4), Some(10000));
        assert_eq!(hybrid36(b"", 4), None);
    }

    #[test]
    fn hybrid36_encode_round_trips_through_decode() {
        for width in [4usize, 5] {
            let base = 10i64.pow(width as u32);
            let upper_span = 26i64 * 36i64.pow(width as u32 - 1);
            for n in [
                0,
                1,
                base - 1,
                base,
                base + 1,
                base + upper_span - 1,
                base + upper_span,
                base + 2 * upper_span - 1,
            ] {
                let s = encode_hybrid36(n, width).unwrap_or_else(|| panic!("{n} w{width}"));
                assert!(s.len() <= width, "{s:?} too wide for {width}");
                assert_eq!(hybrid36(s.as_bytes(), width), Some(n as i32), "{n} -> {s}");
            }
        }
        assert_eq!(encode_hybrid36(-3, 4), Some("-3".to_string()));
        assert_eq!(
            encode_hybrid36(10i64.pow(5) + 2 * 26 * 36i64.pow(4), 5),
            None
        );
    }

    const SMALL: &str = "HEADER    PLANT PROTEIN                           30-APR-81   1CRN
TITLE     WATER STRUCTURE OF A HYDRATED CRYSTAL
HELIX    1  H1 ILE A    7  PRO A   19  1                                  13
SHEET    1  S1 2 THR A   1  CYS A   4  0
MODEL        1
ATOM      1  N   THR A   1      17.047  14.099   3.625  1.00 13.79           N
ATOM      2  CA  THR A   1      16.967  12.784   4.338  1.00 10.80           C
ATOM      3  C   THR A   1      15.685  12.755   5.133  1.00  9.19           C
ATOM      4  N   THR A   2      15.115  11.555   5.265  1.00  9.62           N
HETATM    5 CA    CA A 101      10.000  10.000  10.000  1.00 20.00          CA2+
HETATM    6  O   HOH B 201      20.000  20.000  20.000  1.00 30.00           O
ENDMDL
MODEL        2
ATOM      1  N   THR A   1      17.147  14.099   3.625  1.00 13.79           N
ATOM      2  CA  THR A   1      17.067  12.784   4.338  1.00 10.80           C
ATOM      3  C   THR A   1      15.785  12.755   5.133  1.00  9.19           C
ATOM      4  N   THR A   2      15.215  11.555   5.265  1.00  9.62           N
HETATM    5 CA    CA A 101      10.100  10.000  10.000  1.00 20.00          CA2+
HETATM    6  O   HOH B 201      20.100  20.000  20.000  1.00 30.00           O
ENDMDL
CONECT    3    4
END
";

    #[test]
    fn parses_records_models_and_connectivity() {
        let s = parse(SMALL.as_bytes()).unwrap();
        let t = &s.topology;
        assert_eq!(t.validate(), Ok(()));
        assert_eq!(t.id, "1CRN");
        assert_eq!(t.title, "WATER STRUCTURE OF A HYDRATED CRYSTAL");
        assert_eq!(t.atom_count(), 6);
        assert_eq!(t.residue_count(), 4);
        assert_eq!(t.chain_count(), 2);
        assert_eq!(t.atom_name(1), "CA");
        assert_eq!(t.element[1], Element::CARBON);
        assert_eq!(t.element[4], Element::from_atomic_number(20).unwrap());
        assert_eq!(t.charge[4], 2);
        assert_eq!(t.flags[4], vv_core::flags::HETERO);
        assert_eq!(t.residues[3].seq_id, 201);
        assert_eq!(t.residues[0].ss, SecondaryStructure::Strand);
        assert_eq!(t.residues[1].ss, SecondaryStructure::Strand);
        assert_eq!(t.explicit_bonds.len(), 1);
        assert_eq!(t.explicit_bonds[0].atoms, [2, 3]);
        assert_eq!(s.frame_count(), 2);
        assert_eq!(s.frame(1).positions()[0].x, 17.147);
    }

    /// Builds one fixed-column PDB record, 80 columns, space-padded.
    /// `name` is written exactly as given at columns 13-16, so a test
    /// controls its own leading-space convention (a 2-character name like
    /// `" CA "` for an alpha carbon vs `"CA  "` for a calcium ion).
    fn atom_line(
        record: &str,
        serial: u32,
        name: &str,
        comp: &str,
        chain: char,
        resseq: i32,
    ) -> Vec<u8> {
        let mut line = vec![b' '; 80];
        line[..record.len()].copy_from_slice(record.as_bytes());
        let put = |line: &mut [u8], start: usize, s: &str| {
            line[start..start + s.len()].copy_from_slice(s.as_bytes())
        };
        put(&mut line, 6, &format!("{serial:>5}"));
        put(&mut line, 12, name);
        put(&mut line, 17, comp);
        line[21] = chain as u8;
        put(&mut line, 22, &format!("{resseq:>4}"));
        put(&mut line, 30, "   0.000");
        put(&mut line, 38, "   0.000");
        put(&mut line, 46, "   0.000");
        put(&mut line, 54, "  1.00");
        put(&mut line, 60, "  0.00");
        line
    }

    fn link_line(
        name1: &str,
        comp1: &str,
        chain1: char,
        resseq1: i32,
        name2: &str,
        comp2: &str,
        chain2: char,
        resseq2: i32,
    ) -> Vec<u8> {
        let mut line = vec![b' '; 80];
        line[..6].copy_from_slice(b"LINK  ");
        let put = |line: &mut [u8], start: usize, s: &str| {
            line[start..start + s.len()].copy_from_slice(s.as_bytes())
        };
        put(&mut line, 12, name1);
        put(&mut line, 17, comp1);
        line[21] = chain1 as u8;
        put(&mut line, 22, &format!("{resseq1:>4}"));
        put(&mut line, 42, name2);
        put(&mut line, 47, comp2);
        line[51] = chain2 as u8;
        put(&mut line, 52, &format!("{resseq2:>4}"));
        line
    }

    #[test]
    fn link_records_bond_atoms_named_by_residue_not_serial() {
        // Two residues far enough apart that distance-based bond
        // perception alone would not connect them: only the LINK record
        // should produce the explicit bond.
        let mut src = Vec::new();
        src.extend(atom_line("ATOM  ", 1, "C1  ", "NAG", 'A', 1));
        src.extend(b"\n");
        src.extend({
            let mut l = atom_line("ATOM  ", 2, "ND2 ", "ASN", 'A', 2);
            l[30..38].copy_from_slice(b"  50.000");
            l
        });
        src.extend(b"\n");
        src.extend(link_line("C1  ", "NAG", 'A', 1, "ND2 ", "ASN", 'A', 2));
        src.extend(b"\nEND\n");

        let s = parse(&src).unwrap();
        let t = &s.topology;
        assert_eq!(t.explicit_bonds.len(), 1);
        assert_eq!(t.explicit_bonds[0].atoms, [0, 1]);
    }

    /// A CONECT record naming a serial once per bond order (the old
    /// convention this format used before `_chem_comp_bond`): repeats of
    /// the same neighbour within one line count as extra order.
    fn conect_line(from: u32, others: &[u32]) -> Vec<u8> {
        let mut line = vec![b' '; 80];
        line[..6].copy_from_slice(b"CONECT");
        let put = |line: &mut [u8], start: usize, v: u32| {
            let s = format!("{v:>5}");
            line[start..start + 5].copy_from_slice(s.as_bytes())
        };
        put(&mut line, 6, from);
        for (i, &o) in others.iter().enumerate() {
            put(&mut line, 11 + i * 5, o);
        }
        line
    }

    #[test]
    fn a_conect_record_repeating_a_neighbor_is_a_double_bond() {
        let mut src = Vec::new();
        src.extend(atom_line("ATOM  ", 1, " C1 ", "LIG", 'A', 1));
        src.extend(b"\n");
        src.extend(atom_line("ATOM  ", 2, " C2 ", "LIG", 'A', 1));
        src.extend(b"\n");
        src.extend(conect_line(1, &[2, 2]));
        src.extend(b"\nEND\n");

        let s = parse(&src).unwrap();
        let t = &s.topology;
        assert_eq!(t.explicit_bonds.len(), 1);
        assert_eq!(t.explicit_bonds[0].atoms, [0, 1]);
        assert_eq!(t.explicit_bonds[0].order, vv_core::BondOrder::Double);
    }

    #[test]
    fn a_plain_conect_record_is_a_single_bond() {
        let mut src = Vec::new();
        src.extend(atom_line("ATOM  ", 1, " C1 ", "LIG", 'A', 1));
        src.extend(b"\n");
        src.extend(atom_line("ATOM  ", 2, " C2 ", "LIG", 'A', 1));
        src.extend(b"\n");
        src.extend(conect_line(1, &[2]));
        src.extend(b"\nEND\n");

        let s = parse(&src).unwrap();
        assert_eq!(
            s.topology.explicit_bonds[0].order,
            vv_core::BondOrder::Single
        );
    }

    #[test]
    fn unresolvable_link_record_is_ignored_not_an_error() {
        let mut src = Vec::new();
        src.extend(atom_line("ATOM  ", 1, "C1  ", "NAG", 'A', 1));
        src.extend(b"\n");
        // Names a residue that was never read.
        src.extend(link_line("C1  ", "NAG", 'A', 1, "ND2 ", "ASN", 'A', 99));
        src.extend(b"\nEND\n");

        let s = parse(&src).unwrap();
        assert!(s.topology.explicit_bonds.is_empty());
    }

    #[test]
    fn a_blank_element_column_still_tells_a_calcium_ion_from_an_alpha_carbon() {
        let mut src = Vec::new();
        src.extend(atom_line("ATOM  ", 1, " CA ", "ALA", 'A', 1)); // alpha carbon
        src.extend(b"\n");
        src.extend(atom_line("HETATM", 2, "CA  ", "CA", 'A', 2)); // calcium ion
        src.extend(b"\nEND\n");

        let s = parse(&src).unwrap();
        let t = &s.topology;
        assert_eq!(t.element[0], Element::CARBON);
        assert_eq!(t.element[1], Element::from_atomic_number(20).unwrap());
    }

    #[test]
    fn a_blank_element_column_does_not_read_branched_hydrogens_as_helium_or_mercury() {
        let mut src = Vec::new();
        src.extend(atom_line("ATOM  ", 1, "HG11", "VAL", 'A', 1));
        src.extend(b"\n");
        src.extend(atom_line("ATOM  ", 2, "HE21", "GLN", 'A', 2));
        src.extend(b"\nEND\n");

        let s = parse(&src).unwrap();
        let t = &s.topology;
        assert_eq!(t.element[0], Element::HYDROGEN);
        assert_eq!(t.element[1], Element::HYDROGEN);
    }
}
