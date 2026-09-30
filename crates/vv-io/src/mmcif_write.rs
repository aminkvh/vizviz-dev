//! Minimal mmCIF (PDBx) writer: one `_atom_site` loop, optionally
//! secondary-structure loops and one `_struct_conn` loop for explicit bonds, enough for round-trip tests
//! and for generating large benchmark files. No column-width limits, so
//! (unlike `pdb_write`) atom counts and chain names never need remapping.

use std::io::{self, Write};

use vv_core::fixedbitset::FixedBitSet;
use vv_core::{flags, ExplicitBondKind, SecondaryStructure, Structure, Topology};

use crate::atom_site::SEGID_ITEM;
use crate::mmcif_entity_write::write_entity_tables;
use crate::polymer_layout::Layout;
use crate::write::{read_frames, residue_atoms, ss_runs, SsRun};

const ATOM_SITE_ITEMS: [&str; 19] = [
    "group_PDB",
    "id",
    "type_symbol",
    "label_atom_id",
    "label_alt_id",
    "label_comp_id",
    "label_asym_id",
    "label_entity_id",
    "label_seq_id",
    "pdbx_PDB_ins_code",
    "Cartn_x",
    "Cartn_y",
    "Cartn_z",
    "occupancy",
    "B_iso_or_equiv",
    "pdbx_formal_charge",
    "auth_seq_id",
    "auth_asym_id",
    "pdbx_PDB_model_num",
];

const STRUCT_CONN_ITEMS: [&str; 11] = [
    "conn_type_id",
    "ptnr1_label_asym_id",
    "ptnr1_label_seq_id",
    "pdbx_ptnr1_PDB_ins_code",
    "ptnr1_label_atom_id",
    "pdbx_ptnr1_label_alt_id",
    "ptnr2_label_asym_id",
    "ptnr2_label_seq_id",
    "pdbx_ptnr2_PDB_ins_code",
    "ptnr2_label_atom_id",
    "pdbx_ptnr2_label_alt_id",
];

/// Quotes a token that contains an apostrophe (CIF's own quoting
/// convention); every other token is written bare.
pub(crate) fn token(s: &str) -> String {
    if s.contains('\'') {
        format!("\"{s}\"")
    } else {
        s.to_string()
    }
}

/// Any value as one CIF cell: bare when that reads back the same, else
/// quoted, or a text field when it holds both quote characters or a
/// line break. An empty value is the null `?`.
pub(crate) fn cell(s: &str) -> String {
    let bare = !s.is_empty()
        && !s.starts_with(['_', '#', '$', '\'', '"', ';', '[', ']', '?', '.'])
        && !s.contains(char::is_whitespace);
    if bare {
        return s.to_string();
    }
    match s {
        "" => "?".to_string(),
        s if s.contains('\n') || (s.contains('\'') && s.contains('"')) => {
            format!("\n;{s}\n;\n")
        }
        s if s.contains('\'') => format!("\"{s}\""),
        s => format!("'{s}'"),
    }
}

fn conn_type(kind: ExplicitBondKind) -> &'static str {
    match kind {
        ExplicitBondKind::Disulfide => "disulf",
        ExplicitBondKind::Metal => "metalc",
        ExplicitBondKind::Covalent | ExplicitBondKind::Other => "covale",
    }
}

/// One loop of secondary-structure ranges: `category` and `items` name the
/// category, `lead` gives each row's leading columns (its index is `i`).
fn write_ss_loop(
    t: &Topology,
    labels: &Labels<'_>,
    out: &mut impl Write,
    runs: &[&SsRun],
    category: &str,
    items: &[&str],
    lead: impl Fn(usize) -> String,
) -> io::Result<()> {
    if runs.is_empty() {
        return Ok(());
    }
    writeln!(out, "loop_")?;
    for item in items {
        writeln!(out, "_{category}.{item}")?;
    }
    let ins = |c: u8| if c == 0 { '?' } else { c as char };
    for (i, run) in runs.iter().enumerate() {
        let (a, b) = (
            &t.residues[run.first as usize],
            &t.residues[run.last as usize],
        );
        let asym = labels.of(run.first as usize);
        writeln!(
            out,
            "{} {asym} {} {} {} {}",
            lead(i),
            a.seq_id,
            ins(a.ins_code),
            b.seq_id,
            ins(b.ins_code)
        )?;
    }
    writeln!(out, "#")
}

/// Helix runs as `_struct_conf` and strand runs as `_struct_sheet_range`,
/// by label chain and label sequence number; the insertion-code columns
/// keep residues that share a number apart (a PDB-derived structure has no
/// unique `label_seq_id`).
fn write_secondary_structure(
    t: &Topology,
    labels: &Labels<'_>,
    atoms: Option<&FixedBitSet>,
    out: &mut impl Write,
) -> io::Result<()> {
    let runs = ss_runs(t, atoms);
    let of = |ss| runs.iter().filter(|r| r.ss == ss).collect::<Vec<_>>();
    let range = [
        "beg_label_asym_id",
        "beg_label_seq_id",
        "pdbx_beg_PDB_ins_code",
    ];
    let end = ["end_label_seq_id", "pdbx_end_PDB_ins_code"];
    let items = |lead: &[&'static str]| [lead, &range[..], &end[..]].concat();
    let helices = of(SecondaryStructure::Helix);
    write_ss_loop(
        t,
        labels,
        out,
        &helices,
        "struct_conf",
        &items(&["conf_type_id"]),
        |_| "HELX_P".to_string(),
    )?;
    let strands = of(SecondaryStructure::Strand);
    let lead = ["sheet_id", "id"];
    write_ss_loop(
        t,
        labels,
        out,
        &strands,
        "struct_sheet_range",
        &items(&lead),
        |i| format!("S{} 1", i + 1),
    )
}

/// The trailing segid cell of an atom row (`""` when the column is not
/// written, `.` for a chain record without one).
fn segid_cell(t: &Topology, chain: u32, with_segid: bool) -> String {
    match (with_segid, t.segid(chain as usize)) {
        (false, _) => String::new(),
        (true, "") => " .".to_string(),
        (true, segid) => format!(" {}", token(segid)),
    }
}

/// The `label_asym_id` of each written residue: one per segment of a chain
/// record (`Layout::asym_ids`), so a polymer and the ligands and waters
/// sharing its chain record stay separate. `auth_asym_id` stays the chain
/// name.
struct Labels<'a> {
    layout: &'a Layout,
    asym: Vec<String>,
}

impl Labels<'_> {
    fn of(&self, residue: usize) -> &str {
        self.layout
            .segment(residue)
            .map_or("?", |s| self.asym[s].as_str())
    }
}

/// Writes `structure` as mmCIF. `atoms` selects a subset (`None` is every
/// atom); `frames` are the coordinate sets to emit, each as its own
/// `pdbx_PDB_model_num` (the original frame index + 1). No warnings are
/// possible in this format, but the signature matches the other writers.
pub fn write(
    structure: &Structure,
    atoms: Option<&FixedBitSet>,
    frames: &[usize],
    out: &mut impl Write,
) -> io::Result<Vec<String>> {
    let t = &structure.topology;
    let id = if t.id.is_empty() {
        "vizviz"
    } else {
        t.id.as_str()
    };
    writeln!(out, "data_{id}")?;
    writeln!(out, "#")?;
    writeln!(out, "_entry.id {id}")?;
    if !t.title.is_empty() {
        writeln!(out, "_struct.title '{}'", t.title.replace('\'', ""))?;
    }
    writeln!(out, "#")?;
    let layout = Layout::new(t, atoms);
    write_entity_tables(t, &layout, out)?;
    writeln!(out, "loop_")?;
    for item in ATOM_SITE_ITEMS {
        writeln!(out, "_atom_site.{item}")?;
    }
    let with_segid = t.segids.iter().any(|&s| !t.names.get(s).is_empty());
    if with_segid {
        writeln!(out, "_atom_site.{SEGID_ITEM}")?;
    }
    let labels = Labels {
        layout: &layout,
        asym: layout.asym_ids(t),
    };
    let held = read_frames(structure, frames)?;
    for (&frame, coords) in frames.iter().zip(&held) {
        let positions = coords.positions();
        for (ri, res) in t.residues.iter().enumerate() {
            let chain = &t.chains[res.chain as usize];
            let asym = labels.of(ri);
            let entity = layout.segment(ri).map_or(0, |s| layout.entity_of[s] + 1);
            let auth_asym = t.names.get(chain.auth_asym);
            let comp = t.names.get(res.comp);
            let segid = segid_cell(t, res.chain, with_segid);
            for a in residue_atoms(res, atoms) {
                let a = a as usize;
                let name = token(t.atom_name(a));
                let group = if t.flags.get(a).is_some_and(|f| f & flags::HETERO != 0) {
                    "HETATM"
                } else {
                    "ATOM"
                };
                let alt = t.alt_loc.get(a).copied().filter(|&c| c != 0);
                let ins = if res.ins_code == 0 {
                    '?'
                } else {
                    res.ins_code as char
                };
                let p = positions[a];
                writeln!(
                    out,
                    "{group} {} {} {name} {} {comp} {asym} {} {} {ins} {:.3} {:.3} {:.3} {:.2} {:.2} {} {} {auth_asym} {}{segid}",
                    t.serial.get(a).copied().unwrap_or(a as u32 + 1),
                    if t.is_deuterium(a) { "D" } else { t.element[a].symbol() },
                    alt.map_or(".".to_string(), |c| (c as char).to_string()),
                    entity,
                    res.seq_id,
                    p.x,
                    p.y,
                    p.z,
                    t.occupancy.get(a).copied().unwrap_or(1.0),
                    t.b_factor.get(a).copied().unwrap_or(0.0),
                    t.charge.get(a).copied().unwrap_or(0),
                    res.auth_seq_id,
                    frame + 1,
                )?;
            }
        }
    }
    writeln!(out, "#")?;
    write_secondary_structure(t, &labels, atoms, out)?;

    let bonds: Vec<_> = t
        .explicit_bonds
        .iter()
        .filter(|b| {
            b.atoms
                .iter()
                .all(|&a| atoms.is_none_or(|m| m.contains(a as usize)))
        })
        .collect();
    if !bonds.is_empty() {
        writeln!(out, "loop_")?;
        for item in STRUCT_CONN_ITEMS {
            writeln!(out, "_struct_conn.{item}")?;
        }
        for bond in bonds {
            let partner = |a: u32| {
                let a = a as usize;
                let ri = t.residue_index[a] as usize;
                let res = &t.residues[ri];
                let asym = labels.of(ri);
                let ins = if res.ins_code == 0 {
                    '?'
                } else {
                    res.ins_code as char
                };
                let alt = t
                    .alt_loc
                    .get(a)
                    .copied()
                    .filter(|&c| c != 0)
                    .map_or(".".to_string(), |c| (c as char).to_string());
                format!(
                    "{asym} {} {ins} {} {alt}",
                    res.seq_id,
                    token(t.atom_name(a))
                )
            };
            writeln!(
                out,
                "{} {} {}",
                conn_type(bond.kind),
                partner(bond.atoms[0]),
                partner(bond.atoms[1]),
            )?;
        }
        writeln!(out, "#")?;
    }
    Ok(Vec::new())
}
