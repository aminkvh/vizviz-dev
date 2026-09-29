//! mmCIF reader, written from the PDBx/mmCIF dictionary.
//!
//! The `_atom_site` table is parsed in parallel chunks straight into the
//! columnar model (one row per line, no per-atom allocation); every other
//! category goes through the generic CIF reader in `cif.rs`.

use std::collections::HashMap;

use rayon::prelude::*;
use vv_core::glam::Vec3;
use vv_core::{
    AnnotationCategory, Annotations, AtomExtra, AtomRow, BondOrder, Element, ExplicitBond,
    ExplicitBondKind, SecondaryStructure, Structure, Topology, TopologyBuilder,
};

use crate::cif::{self, is_null, Category, LoopLocation};
use crate::float::{parse_f32, parse_i32};
use crate::mmcif_entity::EntityInfo;
use crate::ss_range::{self, SsRange};
use crate::ParseError;

/// The PDBx dictionary has no segment id. This local item carries a
/// PDB/MD segid so it survives a round trip through mmCIF.
pub(crate) const SEGID_ITEM: &str = "vizviz_segid";

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Role {
    Group,
    Id,
    Element,
    Name,
    AltLoc,
    Comp,
    Asym,
    Entity,
    SeqId,
    InsCode,
    X,
    Y,
    Z,
    Occupancy,
    BFactor,
    Charge,
    AuthSeqId,
    AuthAsym,
    AuthAtomName,
    Model,
    Segid,
}

impl Role {
    fn from_item(item: &str) -> Option<Role> {
        Some(match item {
            "group_PDB" => Role::Group,
            "id" => Role::Id,
            "type_symbol" => Role::Element,
            "label_atom_id" => Role::Name,
            "label_alt_id" => Role::AltLoc,
            "label_comp_id" => Role::Comp,
            "label_asym_id" => Role::Asym,
            "label_entity_id" => Role::Entity,
            "label_seq_id" => Role::SeqId,
            "pdbx_PDB_ins_code" => Role::InsCode,
            "Cartn_x" => Role::X,
            "Cartn_y" => Role::Y,
            "Cartn_z" => Role::Z,
            "occupancy" => Role::Occupancy,
            "B_iso_or_equiv" => Role::BFactor,
            "pdbx_formal_charge" => Role::Charge,
            "auth_seq_id" => Role::AuthSeqId,
            "auth_asym_id" => Role::AuthAsym,
            "auth_atom_id" => Role::AuthAtomName,
            "pdbx_PDB_model_num" => Role::Model,
            SEGID_ITEM => Role::Segid,
            _ => return None,
        })
    }
}

/// Column index of each role, or `usize::MAX` when absent.
struct Columns {
    index: [usize; 21],
    count: usize,
}

impl Columns {
    fn new(items: &[&str]) -> Result<Self, ParseError> {
        let mut index = [usize::MAX; 21];
        for (i, item) in items.iter().enumerate() {
            if let Some(role) = Role::from_item(item) {
                index[role as usize] = i;
            }
        }
        let cols = Self {
            index,
            count: items.len(),
        };
        for (role, name) in [
            (Role::X, "Cartn_x"),
            (Role::Y, "Cartn_y"),
            (Role::Z, "Cartn_z"),
            (Role::Comp, "label_comp_id"),
            (Role::Asym, "label_asym_id"),
        ] {
            if !cols.has(role) {
                return Err(ParseError::MissingColumn(name));
            }
        }
        if !cols.has(Role::Name) && !cols.has(Role::AuthAtomName) {
            return Err(ParseError::MissingColumn("label_atom_id"));
        }
        Ok(cols)
    }

    fn has(&self, role: Role) -> bool {
        self.index[role as usize] != usize::MAX
    }

    fn get<'a>(&self, fields: &[&'a [u8]], role: Role) -> Option<&'a [u8]> {
        let i = self.index[role as usize];
        if i == usize::MAX {
            None
        } else {
            fields.get(i).copied().filter(|v| !is_null(v))
        }
    }
}

/// Splits one row into whitespace-separated fields, honoring quotes.
fn split_fields<'a>(line: &'a [u8], out: &mut Vec<&'a [u8]>) {
    out.clear();
    let mut i = 0;
    let n = line.len();
    while i < n {
        while i < n && line[i].is_ascii_whitespace() {
            i += 1;
        }
        if i >= n {
            break;
        }
        let c = line[i];
        if c == b'\'' || c == b'"' {
            let start = i + 1;
            let mut j = start;
            let mut closed = false;
            while j < n {
                if line[j] == c && (j + 1 >= n || line[j + 1].is_ascii_whitespace()) {
                    closed = true;
                    break;
                }
                j += 1;
            }
            if closed {
                out.push(&line[start..j]);
                i = j + 1;
                continue;
            }
            i = start - 1;
        }
        let start = i;
        while i < n && !line[i].is_ascii_whitespace() {
            i += 1;
        }
        out.push(&line[start..i]);
    }
}

fn atom_name(bytes: &[u8]) -> [u8; 4] {
    let mut name = [b' '; 4];
    let n = bytes.len().min(4);
    name[..n].copy_from_slice(&bytes[..n]);
    name
}

fn parse_charge(v: Option<&[u8]>) -> i8 {
    v.and_then(parse_i32).map_or(0, |c| c.clamp(-9, 9) as i8)
}

fn as_str(v: Option<&[u8]>) -> &str {
    v.and_then(|b| std::str::from_utf8(b).ok()).unwrap_or("")
}

struct ChunkResult {
    builder: TopologyBuilder,
    /// Coordinates of rows belonging to other models, keyed by model number.
    other_models: Vec<(i32, Vec<Vec3>)>,
    error: Option<ParseError>,
}

fn parse_chunk(chunk: &[u8], cols: &Columns, first_model: i32, line_offset: usize) -> ChunkResult {
    let estimate = chunk.len() / 80 + 1;
    let mut builder = TopologyBuilder::with_capacity(estimate);
    let mut other_models: Vec<(i32, Vec<Vec3>)> = Vec::new();
    let mut fields: Vec<&[u8]> = Vec::with_capacity(cols.count);
    for (line_no, line) in (line_offset + 1..).zip(chunk.split(|&b| b == b'\n')) {
        let line = match line.last() {
            Some(b'\r') => &line[..line.len() - 1],
            _ => line,
        };
        if line.is_empty() || line.iter().all(u8::is_ascii_whitespace) || line[0] == b'#' {
            continue;
        }
        split_fields(line, &mut fields);
        if fields.len() != cols.count {
            return ChunkResult {
                builder,
                other_models,
                error: Some(ParseError::Malformed {
                    line: line_no,
                    message: format!(
                        "{ROW_FIELD_MISMATCH} {} fields, header has {}",
                        fields.len(),
                        cols.count
                    ),
                }),
            };
        }
        let coord = |role| cols.get(&fields, role).and_then(parse_f32);
        let (Some(x), Some(y), Some(z)) = (coord(Role::X), coord(Role::Y), coord(Role::Z)) else {
            return ChunkResult {
                builder,
                other_models,
                error: Some(ParseError::Malformed {
                    line: line_no,
                    message: "unreadable coordinates".into(),
                }),
            };
        };
        let position = Vec3::new(x, y, z);
        let model = cols
            .get(&fields, Role::Model)
            .and_then(parse_i32)
            .unwrap_or(first_model);
        if model != first_model {
            match other_models.iter_mut().find(|(m, _)| *m == model) {
                Some((_, v)) => v.push(position),
                None => other_models.push((model, vec![position])),
            }
            continue;
        }

        let name_bytes = cols
            .get(&fields, Role::Name)
            .or_else(|| cols.get(&fields, Role::AuthAtomName))
            .unwrap_or(b"");
        let symbol = cols.get(&fields, Role::Element);
        let element = symbol
            .map(Element::from_symbol)
            .filter(|e| !e.is_unknown())
            .unwrap_or_else(|| Element::from_atom_name(name_bytes));
        let extra = AtomExtra {
            long_name: std::str::from_utf8(name_bytes).ok().filter(|n| n.len() > 4),
            segid: as_str(cols.get(&fields, Role::Segid)),
            deuterium: matches!(symbol, Some(b"D" | b"d")),
        };
        let asym = as_str(cols.get(&fields, Role::Asym));
        let auth_asym = cols
            .get(&fields, Role::AuthAsym)
            .map_or(asym, |a| as_str(Some(a)));
        let auth_seq_id = cols.get(&fields, Role::AuthSeqId).and_then(parse_i32);
        let seq_id = cols.get(&fields, Role::SeqId).and_then(parse_i32);
        let row = AtomRow {
            element,
            name: atom_name(name_bytes),
            serial: cols
                .get(&fields, Role::Id)
                .and_then(parse_i32)
                .unwrap_or(0)
                .max(0) as u32,
            alt_loc: cols.get(&fields, Role::AltLoc).map_or(0, |a| a[0]),
            comp: as_str(cols.get(&fields, Role::Comp)),
            asym,
            auth_asym,
            seq_id: seq_id.or(auth_seq_id).unwrap_or(0),
            auth_seq_id: auth_seq_id.or(seq_id).unwrap_or(0),
            ins_code: cols.get(&fields, Role::InsCode).map_or(0, |c| c[0]),
            entity: cols
                .get(&fields, Role::Entity)
                .and_then(parse_i32)
                .unwrap_or(0) as u16,
            position,
            occupancy: cols
                .get(&fields, Role::Occupancy)
                .and_then(parse_f32)
                .unwrap_or(1.0),
            b_factor: cols
                .get(&fields, Role::BFactor)
                .and_then(parse_f32)
                .unwrap_or(0.0),
            charge: parse_charge(cols.get(&fields, Role::Charge)),
            hetero: cols
                .get(&fields, Role::Group)
                .is_some_and(|g| g == b"HETATM"),
        };
        builder.push_with(&row, &extra);
    }
    ChunkResult {
        builder,
        other_models,
        error: None,
    }
}

/// Splits `body` into roughly equal chunks on line boundaries.
fn chunk_ranges(body: &[u8], target: usize) -> Vec<std::ops::Range<usize>> {
    let mut ranges = Vec::new();
    let mut start = 0;
    while start < body.len() {
        let mut end = (start + target).min(body.len());
        if end < body.len() {
            end = memchr::memchr(b'\n', &body[end..]).map_or(body.len(), |n| end + n + 1);
        }
        ranges.push(start..end);
        start = end;
    }
    ranges
}

const ROW_FIELD_MISMATCH: &str = "atom_site row has";

/// Re-emits `values` one row of `ncols` per line, quoted where a value
/// would not survive whitespace splitting. Lets the line-oriented fast path
/// read tables whose rows wrap over several lines or come as `_atom_site.x
/// value` pairs, both of which CIF allows.
fn one_row_per_line<'a>(values: impl Iterator<Item = &'a [u8]>, ncols: usize) -> Vec<u8> {
    let mut out = Vec::new();
    for (i, v) in values.enumerate() {
        let plain = !v.is_empty()
            && !v.iter().any(u8::is_ascii_whitespace)
            && !matches!(v[0], b'\'' | b'"' | b'#' | b';' | b'_');
        let quote = if v.contains(&b'\'') { b'"' } else { b'\'' };
        if !plain {
            out.push(quote);
        }
        out.extend_from_slice(v);
        if !plain {
            out.push(quote);
        }
        out.push(if (i + 1) % ncols == 0 { b'\n' } else { b' ' });
    }
    out
}

/// The loop body re-flowed, when its token count is a whole number of rows.
fn reflow_rows(body: &[u8], ncols: usize) -> Option<Vec<u8>> {
    let values: Vec<&[u8]> = cif::Tokenizer::new(body)
        .take_while(|t| t.kind == cif::TokenKind::Value)
        .map(|t| t.text)
        .collect();
    (ncols > 0 && values.len() % ncols == 0).then(|| one_row_per_line(values.into_iter(), ncols))
}

fn parse_atom_site(
    src: &[u8],
    loc: &LoopLocation<'_>,
    entity: &EntityInfo,
) -> Result<Structure, ParseError> {
    let body = &src[loc.body.clone()];
    let line_base = memchr::memchr_iter(b'\n', &src[..loc.body.start]).count();
    match parse_rows(&loc.items, body, line_base, entity) {
        Err(ParseError::Malformed { message, .. }) if message.starts_with(ROW_FIELD_MISMATCH) => {
            let flat = reflow_rows(body, loc.items.len());
            flat.map_or_else(
                || parse_rows(&loc.items, body, line_base, entity),
                |flat| parse_rows(&loc.items, &flat, 0, entity),
            )
        }
        other => other,
    }
}

fn parse_rows(
    items: &[&str],
    body: &[u8],
    line_base: usize,
    entity: &EntityInfo,
) -> Result<Structure, ParseError> {
    let cols = Columns::new(items)?;

    let first_model = {
        let mut fields = Vec::new();
        let first_line = body
            .split(|&b| b == b'\n')
            .find(|l| !l.is_empty() && l[0] != b'#' && !l.iter().all(u8::is_ascii_whitespace));
        first_line
            .map(|l| {
                split_fields(l, &mut fields);
                cols.get(&fields, Role::Model)
                    .and_then(parse_i32)
                    .unwrap_or(1)
            })
            .unwrap_or(1)
    };

    let threads = rayon::current_num_threads().max(1);
    let target = (body.len() / (threads * 4)).clamp(1 << 20, 64 << 20);
    let ranges = chunk_ranges(body, target);
    let line_starts: Vec<usize> = {
        let mut acc = line_base;
        ranges
            .iter()
            .map(|r| {
                let s = acc;
                acc += memchr::memchr_iter(b'\n', &body[r.clone()]).count();
                s
            })
            .collect()
    };

    let results: Vec<ChunkResult> = ranges
        .par_iter()
        .zip(line_starts.par_iter())
        .map(|(r, &line_start)| parse_chunk(&body[r.clone()], &cols, first_model, line_start))
        .collect();

    let mut builder = TopologyBuilder::new();
    let mut models: Vec<(i32, Vec<Vec3>)> = Vec::new();
    for result in results {
        if let Some(e) = result.error {
            return Err(e);
        }
        builder.append(result.builder);
        for (model, mut coords) in result.other_models {
            match models.iter_mut().find(|(m, _)| *m == model) {
                Some((_, v)) => v.append(&mut coords),
                None => models.push((model, coords)),
            }
        }
    }
    if builder.atom_count() == 0 {
        return Err(ParseError::NoAtoms);
    }
    let atom_count = builder.atom_count();
    models.sort_by_key(|(m, _)| *m);
    let extra: Vec<Vec<Vec3>> = models
        .into_iter()
        .filter(|(_, v)| v.len() == atom_count)
        .map(|(_, v)| v)
        .collect();
    builder.topology.polymer_hint = entity.hints(&builder.topology);
    builder.finish_with_frames(extra).map_err(ParseError::from)
}

/// One row's range, by label chain and label sequence number; the
/// insertion code (`pdbx_beg_PDB_ins_code`) tells apart residues that
/// share a number when a writer had no unique `label_seq_id` to give.
fn ss_range_of(cat: &Category<'_>, row: usize) -> Option<SsRange> {
    let residue = |p: &str| {
        let seq = cat
            .get(row, &format!("{p}_label_seq_id"))
            .and_then(parse_i32)?;
        let ins = cat
            .get(row, &format!("pdbx_{p}_PDB_ins_code"))
            .map_or(0, |c| c[0]);
        Some((seq, ins))
    };
    Some(SsRange {
        chain: cat.get_str(row, "beg_label_asym_id")?.to_string(),
        first: residue("beg")?,
        last: residue("end")?,
    })
}

fn apply_secondary_structure(topology: &mut Topology, cats: &HashMap<String, Category<'_>>) {
    if let Some(conf) = cats.get("struct_conf") {
        for row in 0..conf.rows.len() {
            let kind = conf.get_str(row, "conf_type_id").unwrap_or("");
            let ss = if kind.starts_with("HELX") {
                SecondaryStructure::Helix
            } else if kind.starts_with("STRN") {
                SecondaryStructure::Strand
            } else {
                SecondaryStructure::Coil
            };
            if let Some(range) = ss_range_of(conf, row) {
                ss_range::apply(topology, &[range], ss);
            }
        }
    }
    if let Some(sheet) = cats.get("struct_sheet_range") {
        let ranges: Vec<SsRange> = (0..sheet.rows.len())
            .filter_map(|row| ss_range_of(sheet, row))
            .collect();
        ss_range::apply(topology, &ranges, SecondaryStructure::Strand);
    }
}

/// Finds an atom by residue key and atom name.
fn find_atom(
    topology: &Topology,
    lookup: &HashMap<(vv_core::InternId, i32, u8), u32>,
    asym: &str,
    seq_id: i32,
    ins: u8,
    name: &str,
    alt: Option<u8>,
) -> Option<u32> {
    let asym_id = topology.names.lookup(asym)?;
    let residue = *lookup.get(&(asym_id, seq_id, ins))?;
    let range = topology.residues[residue as usize].atoms.clone();
    range.into_iter().find(|&a| {
        topology.atom_name(a as usize) == name
            && alt.is_none_or(|alt| {
                topology.alt_loc[a as usize] == 0 || topology.alt_loc[a as usize] == alt
            })
    })
}

fn apply_struct_conn(topology: &mut Topology, cats: &HashMap<String, Category<'_>>) {
    let Some(conn) = cats.get("struct_conn") else {
        return;
    };
    let mut lookup: HashMap<(vv_core::InternId, i32, u8), u32> =
        HashMap::with_capacity(topology.residues.len());
    for (i, res) in topology.residues.iter().enumerate() {
        let asym = topology.chains[res.chain as usize].label_asym;
        lookup
            .entry((asym, res.seq_id, res.ins_code))
            .or_insert(i as u32);
        if res.auth_seq_id != res.seq_id {
            lookup
                .entry((asym, res.auth_seq_id, res.ins_code))
                .or_insert(i as u32);
        }
    }
    let mut bonds = Vec::new();
    for row in 0..conn.rows.len() {
        let kind = match conn.get_str(row, "conn_type_id").unwrap_or("") {
            "disulf" => ExplicitBondKind::Disulfide,
            "metalc" => ExplicitBondKind::Metal,
            t if t.starts_with("covale") => ExplicitBondKind::Covalent,
            _ => continue,
        };
        let partner = |p: &str| {
            let asym = conn.get_str(row, &format!("{p}_label_asym_id"))?;
            let seq = conn
                .get(row, &format!("{p}_label_seq_id"))
                .or_else(|| conn.get(row, &format!("{p}_auth_seq_id")))
                .and_then(parse_i32)?;
            let ins = conn
                .get(row, &format!("pdbx_{p}_PDB_ins_code"))
                .map_or(0, |c| c[0]);
            let name = conn.get_str(row, &format!("{p}_label_atom_id"))?;
            let alt = conn
                .get(row, &format!("pdbx_{p}_label_alt_id"))
                .map(|a| a[0]);
            find_atom(topology, &lookup, asym, seq, ins, name, alt)
        };
        if let (Some(a), Some(b)) = (partner("ptnr1"), partner("ptnr2")) {
            bonds.push(ExplicitBond::new(a, b, kind));
        }
    }
    topology.explicit_bonds = bonds;
}

/// `_chem_comp_bond.value_order` -> [`BondOrder`] (case-insensitive: the
/// wwPDB CCD's own ligand exports use `SING`/`DOUB`, a full entry's inline
/// copy -- as our `fixtures/small/*.cif` carry -- uses `sing`/`doub`).
/// Unrecognized text (or `SING` itself) is `Single`, so only a real
/// non-single order costs an entry in the sparse override list.
fn parse_value_order(s: &str) -> BondOrder {
    match s.to_ascii_uppercase().as_str() {
        "DOUB" => BondOrder::Double,
        "TRIP" | "QUAD" => BondOrder::Triple,
        "AROM" | "DELO" => BondOrder::Aromatic,
        _ => BondOrder::Single,
    }
}

/// The atom named `name` within `res`'s own atom range, or `None` (a
/// nonstandard naming convention the file's own chemistry can't resolve
/// against, same fallback as a residue template's coverage check).
fn find_in_residue(topology: &Topology, res: &vv_core::ResidueRec, name: &str) -> Option<u32> {
    res.atoms
        .clone()
        .find(|&a| topology.atom_name(a as usize) == name)
}

/// Source #1 (see `vv_core::bonds` module doc): a file's own
/// `_chem_comp_bond` loop, when it carries one -- a "full" mmCIF entry's
/// inline copy for every component present (as `fixtures/small/*.cif`
/// do) or a standalone CCD ligand export. Grouped by `comp_id` once (a
/// handful of distinct residue types even in a huge structure -- CCD
/// bonds are per component type, not per instance), then resolved by
/// atom name within each matching residue's own range, non-`Single` rows
/// only: `perceive` still decides *whether* two atoms are bonded, this
/// only overrides the order once it agrees they are.
fn apply_chem_comp_bond(topology: &mut Topology, cats: &HashMap<String, Category<'_>>) {
    let Some(cat) = cats.get("chem_comp_bond") else {
        return;
    };
    // Keyed by `InternId`, not `&str`: every residue looks itself up here,
    // so hashing a short string (or scanning the intern table by name)
    // per residue -- rather than per distinct component type -- would
    // dominate this pass on a huge structure (4M atoms measured ~150ms
    // of it, `docs/VALIDATION.md`-style scratch timing on 8GLV).
    let mut by_comp: HashMap<vv_core::InternId, Vec<(&str, &str, BondOrder)>> = HashMap::new();
    for row in 0..cat.rows.len() {
        let Some(comp) = cat.get_str(row, "comp_id") else {
            continue;
        };
        let order = parse_value_order(cat.get_str(row, "value_order").unwrap_or(""));
        if order == BondOrder::Single {
            continue; // nothing to override
        }
        let (Some(a1), Some(a2)) = (cat.get_str(row, "atom_id_1"), cat.get_str(row, "atom_id_2"))
        else {
            continue;
        };
        // Only a component type some residue actually has needs an id
        // (`lookup` never interns): a ligand-only CCD block otherwise
        // adds ids for compounds this structure never uses.
        let Some(id) = topology.names.lookup(comp) else {
            continue;
        };
        by_comp.entry(id).or_default().push((a1, a2, order));
    }
    if by_comp.is_empty() {
        return;
    }
    let mut out = Vec::new();
    for res in &topology.residues {
        let Some(bonds) = by_comp.get(&res.comp) else {
            continue;
        };
        for &(n1, n2, order) in bonds {
            if let (Some(a), Some(b)) = (
                find_in_residue(topology, res, n1),
                find_in_residue(topology, res, n2),
            ) {
                out.push((if a < b { [a, b] } else { [b, a] }, order));
            }
        }
    }
    topology.chem_comp_bond_order = out;
}

/// Header categories worth keeping, with an optional item whitelist and row
/// cap so that verbose tables (author lists, revision logs) stay small.
/// `None` items means "every item the file has".
const ANNOTATION_CATEGORIES: &[(&str, Option<&[&str]>, usize)] = &[
    ("entry", None, usize::MAX),
    ("struct", None, usize::MAX),
    ("struct_keywords", None, usize::MAX),
    ("exptl", None, usize::MAX),
    (
        "refine",
        Some(&["ls_d_res_high", "ls_R_factor_R_work", "ls_R_factor_R_free"]),
        usize::MAX,
    ),
    ("cell", None, usize::MAX),
    ("symmetry", None, usize::MAX),
    ("citation", None, usize::MAX),
    ("citation_author", None, 50),
    ("pdbx_database_status", None, usize::MAX),
    ("em_3d_reconstruction", None, usize::MAX),
    ("entity", None, usize::MAX),
    ("entity_src_gen", None, usize::MAX),
    ("entity_src_nat", None, usize::MAX),
    ("struct_ref", None, usize::MAX),
    ("pdbx_struct_assembly", None, usize::MAX),
    ("pdbx_struct_assembly_gen", None, usize::MAX),
    ("database_2", None, usize::MAX),
    ("audit_author", None, 50),
    (
        "pdbx_audit_revision_history",
        Some(&["revision_date", "major_revision", "minor_revision"]),
        50,
    ),
];

/// Copies the header categories out of the parsed CIF. Values keep the
/// file's text (quotes already stripped); `?`/`.` become empty strings.
/// Nothing here can fail: absent categories and items are simply skipped.
fn collect_annotations(cats: &HashMap<String, Category<'_>>) -> Annotations {
    let mut out = Annotations::default();
    for &(name, wanted, max_rows) in ANNOTATION_CATEGORIES {
        let Some(cat) = cats.get(name) else {
            continue;
        };
        // Column indices to keep, in file order.
        let keep: Vec<usize> = (0..cat.items.len())
            .filter(|&c| {
                wanted.is_none_or(|w| w.iter().any(|i| i.eq_ignore_ascii_case(cat.items[c])))
            })
            .collect();
        if keep.is_empty() {
            continue;
        }
        let rows: Vec<Vec<String>> = cat
            .rows
            .iter()
            .take(max_rows)
            .map(|row| {
                keep.iter()
                    .map(|&c| match row.get(c) {
                        Some(v) if !is_null(v) => String::from_utf8_lossy(v).trim().to_string(),
                        _ => String::new(),
                    })
                    .collect()
            })
            .collect();
        out.categories.push(AnnotationCategory {
            name: name.to_string(),
            items: keep.iter().map(|&c| cat.items[c].to_string()).collect(),
            rows,
        });
    }
    out
}

/// The first data block holding an `_atom_site` table, as a structure
/// plus that block's other categories.
fn read_first_atom_block(
    src: &[u8],
) -> Result<(Structure, HashMap<String, Category<'_>>), ParseError> {
    for block in 0.. {
        let read = cif::read_block(src, &["atom_site"], block);
        let entity = EntityInfo::collect(&read.categories);
        if let Some(loc) = read.skipped.get("atom_site") {
            return Ok((parse_atom_site(src, loc, &entity)?, read.categories));
        }
        if let Some(single) = read.categories.get("atom_site") {
            let ncols = single.items.len();
            let body = one_row_per_line(single.rows[0].iter().copied(), ncols);
            return Ok((
                parse_rows(&single.items, &body, 0, &entity)?,
                read.categories,
            ));
        }
        if read.blocks_seen <= block {
            break;
        }
    }
    Err(ParseError::NoAtoms)
}

pub fn parse(src: &[u8]) -> Result<Structure, ParseError> {
    let (mut structure, cats) = read_first_atom_block(src)?;
    let topology = std::sync::Arc::get_mut(&mut structure.topology).expect("fresh structure");
    if let Some(entry) = cats.get("entry") {
        topology.id = entry.get_str(0, "id").unwrap_or("").to_string();
    }
    if let Some(s) = cats.get("struct") {
        topology.title = s.get_str(0, "title").unwrap_or("").trim().to_string();
    }
    topology.annotations = collect_annotations(&cats);
    apply_secondary_structure(topology, &cats);
    apply_struct_conn(topology, &cats);
    apply_chem_comp_bond(topology, &cats);
    Ok(structure)
}

#[cfg(test)]
mod tests {
    use super::*;

    const SMALL: &str = "data_test
_entry.id TEST
_struct.title 'Two residues'
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_asym_id
_atom_site.pdbx_PDB_model_num
ATOM 1 N N . ALA A 1 1 ? 0.000 0.000 0.000 1.00 10.0 1 A 1
ATOM 2 C CA . ALA A 1 1 ? 1.458 0.000 0.000 1.00 10.0 1 A 1
ATOM 3 C C . ALA A 1 1 ? 2.009 1.420 0.000 1.00 10.0 1 A 1
ATOM 4 N N . GLY A 1 2 ? 3.332 1.536 0.000 1.00 10.0 2 A 1
HETATM 5 O O . HOH B 2 . ? 9.000 9.000 9.000 1.00 30.0 101 B 1
ATOM 1 N N . ALA A 1 1 ? 0.100 0.000 0.000 1.00 10.0 1 A 2
ATOM 2 C CA . ALA A 1 1 ? 1.558 0.000 0.000 1.00 10.0 1 A 2
ATOM 3 C C . ALA A 1 1 ? 2.109 1.420 0.000 1.00 10.0 1 A 2
ATOM 4 N N . GLY A 1 2 ? 3.432 1.536 0.000 1.00 10.0 2 A 2
HETATM 5 O O . HOH B 2 . ? 9.100 9.000 9.000 1.00 30.0 101 B 2
#
loop_
_struct_conf.conf_type_id
_struct_conf.beg_label_asym_id
_struct_conf.beg_label_seq_id
_struct_conf.end_label_seq_id
HELX_P A 1 2
#
";

    #[test]
    fn splits_quoted_fields() {
        let mut out = Vec::new();
        split_fields(b"ATOM 1 C \"C1'\" 'A B' . 1.0", &mut out);
        assert_eq!(
            out,
            vec![&b"ATOM"[..], b"1", b"C", b"C1'", b"A B", b".", b"1.0"]
        );
    }

    #[test]
    fn parses_a_small_two_model_file() {
        let s = parse(SMALL.as_bytes()).unwrap();
        let t = &s.topology;
        assert_eq!(t.validate(), Ok(()));
        assert_eq!(t.id, "TEST");
        assert_eq!(t.title, "Two residues");
        assert_eq!(t.atom_count(), 5);
        assert_eq!(t.residue_count(), 3);
        assert_eq!(t.chain_count(), 2);
        assert_eq!(t.atom_name(1), "CA");
        assert_eq!(t.residue_name(2), "HOH");
        assert_eq!(t.residues[2].auth_seq_id, 101);
        assert_eq!(t.flags[4], vv_core::flags::HETERO);
        assert_eq!(t.residues[0].ss, SecondaryStructure::Helix);
        assert_eq!(t.residues[1].ss, SecondaryStructure::Helix);
        assert_eq!(t.residues[2].ss, SecondaryStructure::Unknown);
        assert_eq!(s.frame_count(), 2);
        assert_eq!(s.frame(0).positions()[1].x, 1.458);
        assert_eq!(s.frame(1).positions()[1].x, 1.558);
    }

    #[test]
    fn rejects_files_without_coordinates() {
        let err = parse(b"data_x\n_entry.id X\n").unwrap_err();
        assert!(matches!(err, ParseError::NoAtoms));
        let err = parse(b"data_x\nloop_\n_atom_site.id\n_atom_site.Cartn_x\n1 0.0\n").unwrap_err();
        assert!(matches!(err, ParseError::MissingColumn(_)));
    }

    #[test]
    fn reports_malformed_rows_with_line_numbers() {
        let bad = SMALL.replace("ATOM 4 N N . GLY A 1 2 ? 3.332", "ATOM 4 N N . GLY A 1 2 ?");
        let err = parse(bad.as_bytes()).unwrap_err();
        match err {
            ParseError::Malformed { line, .. } => assert_eq!(line, 26),
            other => panic!("unexpected {other:?}"),
        }
    }
}
