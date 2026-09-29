//! Legacy PDB writer, wwPDB format v3.3 fixed columns.
//!
//! Serial numbers and residue sequence numbers beyond 5/4 decimal digits
//! use the hybrid-36 extension (`pdb::encode_hybrid36`) rather than
//! wrapping, matching what `pdb::hybrid36` already decodes on read.
//! Chains are written by author name, one character each; multi-character
//! names map onto a free character from `A-Za-z0-9`, with a warning.

use std::io::{self, Write};

use vv_core::fixedbitset::FixedBitSet;
use vv_core::{
    flags, Element, ExplicitBondKind, ResidueRec, SecondaryStructure, Structure, Topology,
};

use crate::pdb::encode_hybrid36;
use crate::write::{read_frames, residue_atoms, ss_runs, SsRun};

/// PDB columns 13-16 for one atom: the element symbol occupies columns
/// 13-14 (right-justified for a 1-letter symbol, left-justified for a
/// 2-letter one), and any remaining name characters follow in 15-16. A
/// name that already fills all 4 columns (branched hydrogen locants like
/// `HG11`) is written verbatim. Inverse of `pdb::element_of`'s heuristic.
fn atom_name_field(name: &str, element: Element) -> [u8; 4] {
    let mut field = [b' '; 4];
    let bytes = name.as_bytes();
    let n = bytes.len().min(4);
    if n == 4 {
        field.copy_from_slice(&bytes[..4]);
        return field;
    }
    let start = if element.symbol().len() == 2 { 0 } else { 1 };
    let end = (start + n).min(4);
    field[start..end].copy_from_slice(&bytes[..end - start]);
    field
}

const CHAIN_POOL: &[u8; 62] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";

/// One PDB chain character per chain, from its author name: every chain
/// sharing a name (a protein chain and its ligands and waters) shares the
/// character, as the PDB archive writes them. Names longer than one
/// character take a pool character no single-character name uses (`?`
/// once the pool runs out), with one warning per such name. Only chains
/// `included` gets a real character.
fn assign_chain_ids(
    topology: &Topology,
    included: impl Fn(usize) -> bool,
) -> (Vec<u8>, Vec<String>) {
    let name_of = |c: &vv_core::ChainRec| topology.names.get(c.auth_asym);
    let literal = |name: &str| match name.as_bytes() {
        [c] if c.is_ascii_alphanumeric() => Some(*c),
        _ => None,
    };
    let mut used = [false; 256];
    for (i, chain) in topology.chains.iter().enumerate() {
        if let Some(c) = literal(name_of(chain)).filter(|_| included(i)) {
            used[c as usize] = true;
        }
    }
    let mut remapped: std::collections::HashMap<&str, u8> = Default::default();
    let mut warnings = Vec::new();
    let ids = topology
        .chains
        .iter()
        .enumerate()
        .map(|(i, chain)| {
            let name = name_of(chain);
            if !included(i) {
                return b'?';
            }
            if let Some(c) = literal(name) {
                return c;
            }
            *remapped.entry(name).or_insert_with(|| {
                let id = CHAIN_POOL
                    .iter()
                    .copied()
                    .find(|&c| !used[c as usize])
                    .unwrap_or(b'?');
                used[id as usize] = true;
                warnings.push(format!(
                    "chain `{name}` written as `{}` (PDB chain ids are one character)",
                    id as char
                ));
                id
            })
        })
        .collect();
    (ids, warnings)
}

/// Whether chain `chain` has any atom `atoms` keeps.
fn chain_included(topology: &Topology, atoms: Option<&FixedBitSet>, chain: usize) -> bool {
    let Some(mask) = atoms else { return true };
    topology.chains[chain].residues.clone().any(|r| {
        topology.residues[r as usize]
            .atoms
            .clone()
            .any(|a| mask.contains(a as usize))
    })
}

/// A field that never exceeds `width` decimal digits (with an optional
/// leading `-`) is written as-is; wider values switch to hybrid-36.
/// Falls back to `?` fill only when even hybrid-36 overflows (a residue
/// count in the tens of millions).
fn number_field(n: i64, width: usize) -> String {
    encode_hybrid36(n, width).unwrap_or_else(|| "?".repeat(width))
}

fn charge_field(charge: i8) -> String {
    if charge == 0 {
        "  ".to_string()
    } else {
        format!(
            "{}{}",
            charge.unsigned_abs(),
            if charge < 0 { '-' } else { '+' }
        )
    }
}

/// One step of the plan built once from the topology (independent of
/// frame): an atom to place, or a TER closing the polymer run that ended
/// at `residue`. Replayed for every model.
enum PlanItem {
    Atom(u32),
    Ter {
        serial: u32,
        chain_byte: u8,
        residue: u32,
    },
}

/// The write order and serial numbering: residues in file order, a TER
/// right after the last non-HETATM residue of a run before a HETATM
/// residue or a chain change, and one more if the structure ends on one.
fn build_plan(
    topology: &Topology,
    mask: Option<&FixedBitSet>,
    chain_ids: &[u8],
) -> (Vec<PlanItem>, Vec<Option<u32>>) {
    let mut plan = Vec::new();
    let mut serial_of = vec![None; topology.atom_count()];
    let mut serial = 1u32;
    let mut open_polymer: Option<(u32, u32)> = None; // (chain, residue) awaiting TER
    let is_hetero = |atom: u32| {
        topology
            .flags
            .get(atom as usize)
            .is_some_and(|f| f & flags::HETERO != 0)
    };
    for (ri, res) in topology.residues.iter().enumerate() {
        let atoms = residue_atoms(res, mask);
        let Some(&first) = atoms.first() else {
            continue;
        };
        let hetero = is_hetero(first);
        if let Some((chain, residue)) = open_polymer {
            if chain != res.chain || hetero {
                plan.push(PlanItem::Ter {
                    serial,
                    chain_byte: chain_ids[chain as usize],
                    residue,
                });
                serial += 1;
                open_polymer = None;
            }
        }
        for a in atoms {
            serial_of[a as usize] = Some(serial);
            plan.push(PlanItem::Atom(a));
            serial += 1;
        }
        if !hetero {
            open_polymer = Some((res.chain, ri as u32));
        }
    }
    if let Some((chain, residue)) = open_polymer {
        plan.push(PlanItem::Ter {
            serial,
            chain_byte: chain_ids[chain as usize],
            residue,
        });
    }
    (plan, serial_of)
}

fn write_padded(out: &mut impl Write, line: &str) -> io::Result<()> {
    writeln!(out, "{line:<80}")
}

/// Columns 18-21. The format reserves 18-20; column 21 (blank there) holds
/// a fourth character, as `pdb::parse` accepts. Truncated past that.
fn residue_name_field(topology: &Topology, residue: &ResidueRec) -> String {
    let name = topology.names.get(residue.comp);
    let n = name.len().min(4);
    format!("{:<4}", &name[..n])
}

/// The 3-column residue name of records whose column 21 is not part of it
/// (`SSBOND`, `HELIX`, `SHEET`).
fn residue_name3(topology: &Topology, residue: &ResidueRec) -> String {
    let name = topology.names.get(residue.comp);
    format!("{:<3}", &name[..name.len().min(3)])
}

fn icode_char(residue: &ResidueRec) -> char {
    match residue.ins_code {
        0 => ' ',
        c => c as char,
    }
}

fn atom_line(
    topology: &Topology,
    position: vv_core::glam::Vec3,
    atom: u32,
    serial: u32,
    chain_byte: u8,
) -> String {
    let a = atom as usize;
    let residue = &topology.residues[topology.residue_index[a] as usize];
    let hetero = topology
        .flags
        .get(a)
        .is_some_and(|f| f & flags::HETERO != 0);
    let name = atom_name_field(topology.atom_name(a), topology.element[a]);
    let name = std::str::from_utf8(&name).expect("ASCII field");
    let alt = topology
        .alt_loc
        .get(a)
        .copied()
        .filter(|&c| c != 0)
        .unwrap_or(b' ') as char;
    let icode = icode_char(residue);
    format!(
        "{rec:<6}{serial:>5} {name}{alt}{resname}{chain}{resseq:>4}{icode}   {x:>8.3}{y:>8.3}{z:>8.3}{occ:>6.2}{bf:>6.2}          {element:>2}{charge}",
        rec = if hetero { "HETATM" } else { "ATOM" },
        serial = number_field(serial as i64, 5),
        resname = residue_name_field(topology, residue),
        chain = chain_byte as char,
        resseq = number_field(residue.auth_seq_id as i64, 4),
        x = position.x,
        y = position.y,
        z = position.z,
        occ = topology.occupancy.get(a).copied().unwrap_or(1.0),
        bf = topology.b_factor.get(a).copied().unwrap_or(0.0),
        element = topology.element[a].symbol(),
        charge = charge_field(topology.charge.get(a).copied().unwrap_or(0)),
    )
}

fn ter_line(topology: &Topology, residue: u32, serial: u32, chain_byte: u8) -> String {
    let res = &topology.residues[residue as usize];
    let icode = icode_char(res);
    format!(
        "{rec:<6}{serial:>5}      {resname}{chain}{resseq:>4}{icode}",
        rec = "TER",
        serial = number_field(serial as i64, 5),
        resname = residue_name_field(topology, res),
        chain = chain_byte as char,
        resseq = number_field(res.auth_seq_id as i64, 4),
    )
}

/// `CRYST1` from `_cell`/`_symmetry`-shaped annotations (mmCIF's `cell`
/// category kept verbatim by the reader), only when all six geometry
/// values parse; `None` when no cell is known.
fn cryst1_line(topology: &Topology) -> Option<String> {
    let cell = topology.annotations.category("cell")?;
    let get = |item: &str| cell.get(item, 0)?.trim().parse::<f32>().ok();
    let (a, b, c) = (get("length_a")?, get("length_b")?, get("length_c")?);
    let (alpha, beta, gamma) = (get("angle_alpha")?, get("angle_beta")?, get("angle_gamma")?);
    let space_group = topology
        .annotations
        .get("symmetry", "space_group_name_H-M")
        .unwrap_or("P 1");
    let z = cell
        .get("Z_PDB", 0)
        .and_then(|v| v.trim().parse::<u32>().ok())
        .unwrap_or(1);
    Some(format!(
        "{rec:<6}{a:>9.3}{b:>9.3}{c:>9.3}{alpha:>7.2}{beta:>7.2}{gamma:>7.2} {sg:<11}{z:>4}",
        rec = "CRYST1",
        sg = space_group,
    ))
}

/// A disulfide's `SSBOND` record (bonus: `CONECT` already round-trips the
/// bond on its own, since the reader parses only `CONECT`/`LINK`, not
/// `SSBOND`).
fn ssbond_line(topology: &Topology, atoms: [u32; 2], serial: u32, chain_ids: &[u8]) -> String {
    let residue = |a: u32| &topology.residues[topology.residue_index[a as usize] as usize];
    let (r1, r2) = (residue(atoms[0]), residue(atoms[1]));
    let icode = icode_char;
    format!(
        "{rec:<6} {serial:>3} {n1} {c1} {s1:>4}{i1}   {n2} {c2} {s2:>4}{i2}",
        rec = "SSBOND",
        n1 = residue_name3(topology, r1),
        c1 = chain_ids[r1.chain as usize] as char,
        s1 = number_field(r1.auth_seq_id as i64, 4),
        i1 = icode(r1),
        n2 = residue_name3(topology, r2),
        c2 = chain_ids[r2.chain as usize] as char,
        s2 = number_field(r2.auth_seq_id as i64, 4),
        i2 = icode(r2),
    )
}

/// `HELIX` and `SHEET` records (v3.3 columns); every strand is written as
/// a one-strand sheet, which is all the reader needs.
fn ss_lines(topology: &Topology, runs: &[SsRun], chain_ids: &[u8]) -> Vec<String> {
    let (mut helix_no, mut strand_no) = (0, 0);
    runs.iter()
        .map(|run| {
            let (a, b) = (
                &topology.residues[run.first as usize],
                &topology.residues[run.last as usize],
            );
            let chain = chain_ids[a.chain as usize] as char;
            let (n1, n2) = (residue_name3(topology, a), residue_name3(topology, b));
            let (s1, s2) = (number_field(a.auth_seq_id as i64, 4), number_field(b.auth_seq_id as i64, 4));
            let (i1, i2) = (icode_char(a), icode_char(b));
            if run.ss == SecondaryStructure::Helix {
                helix_no += 1;
                let len = run.last - run.first + 1;
                format!(
                    "HELIX  {helix_no:>3} {helix_no:>3} {n1} {chain} {s1:>4}{i1} {n2} {chain} {s2:>4}{i2} 1{cmt:<30} {len:>5}",
                    cmt = ""
                )
            } else {
                strand_no += 1;
                format!(
                    "SHEET  {strand_no:>3} {strand_no:>3} 1 {n1} {chain}{s1:>4}{i1} {n2} {chain}{s2:>4}{i2} 0"
                )
            }
        })
        .collect()
}

fn conect_lines(new_serial: &[Option<u32>], bonds: &[vv_core::ExplicitBond]) -> Vec<String> {
    // The source file may already list a bond in both directions (as
    // separate CONECT records); normalize to one unordered pair each so
    // writing both directions below doesn't quadruple it.
    let mut pairs: Vec<[u32; 2]> = bonds
        .iter()
        .filter_map(|b| {
            let [a, c] = b.atoms;
            match (new_serial[a as usize], new_serial[c as usize]) {
                (Some(_), Some(_)) if a < c => Some([a, c]),
                (Some(_), Some(_)) => Some([c, a]),
                _ => None,
            }
        })
        .collect();
    pairs.sort_unstable();
    pairs.dedup();

    let n = new_serial.len();
    let mut adjacency: Vec<Vec<u32>> = vec![Vec::new(); n];
    for [a, b] in pairs {
        adjacency[a as usize].push(new_serial[b as usize].expect("paired"));
        adjacency[b as usize].push(new_serial[a as usize].expect("paired"));
    }
    let mut lines = Vec::new();
    for (atom, partners) in adjacency.iter().enumerate() {
        if partners.is_empty() {
            continue;
        }
        let base = new_serial[atom].expect("has partners, so was selected");
        for chunk in partners.chunks(4) {
            let tail: String = chunk
                .iter()
                .map(|&p| format!("{:>5}", number_field(p as i64, 5)))
                .collect();
            lines.push(format!(
                "{rec:<6}{base:>5}{tail}",
                rec = "CONECT",
                base = number_field(base as i64, 5),
            ));
        }
    }
    lines
}

/// Writes `structure` as legacy PDB. `atoms` selects a subset (`None` is
/// every atom); `frames` are the coordinate sets to emit, wrapped in
/// `MODEL`/`ENDMDL` when there is more than one. Returns one warning per
/// chain whose id had to be remapped.
pub fn write(
    structure: &Structure,
    atoms: Option<&FixedBitSet>,
    frames: &[usize],
    out: &mut impl Write,
) -> io::Result<Vec<String>> {
    let topology = &structure.topology;
    let (chain_ids, warnings) = assign_chain_ids(topology, |c| chain_included(topology, atoms, c));
    let (plan, new_serial) = build_plan(topology, atoms, &chain_ids);

    for line in ss_lines(topology, &ss_runs(topology, atoms), &chain_ids) {
        write_padded(out, &line)?;
    }
    if let Some(cryst1) = cryst1_line(topology) {
        write_padded(out, &cryst1)?;
    }

    let mut ssbond_serial = 1u32;
    for bond in &topology.explicit_bonds {
        if bond.kind != ExplicitBondKind::Disulfide {
            continue;
        }
        let [a, b] = bond.atoms;
        if new_serial[a as usize].is_some() && new_serial[b as usize].is_some() {
            write_padded(
                out,
                &ssbond_line(topology, [a, b], ssbond_serial, &chain_ids),
            )?;
            ssbond_serial += 1;
        }
    }

    let held = read_frames(structure, frames)?;
    let multi = frames.len() > 1;
    for (model, coords) in frames.iter().zip(&held) {
        if multi {
            write_padded(out, &format!("{:<6}    {:>4}", "MODEL", model + 1))?;
        }
        let positions = coords.positions();
        for item in &plan {
            match *item {
                PlanItem::Atom(a) => {
                    let serial = new_serial[a as usize].expect("every planned atom got a serial");
                    let chain_byte = chain_ids[topology.chain_of_atom(a as usize) as usize];
                    let line = atom_line(topology, positions[a as usize], a, serial, chain_byte);
                    write_padded(out, &line)?;
                }
                PlanItem::Ter {
                    serial,
                    chain_byte,
                    residue,
                } => {
                    write_padded(out, &ter_line(topology, residue, serial, chain_byte))?;
                }
            }
        }
        if multi {
            write_padded(out, "ENDMDL")?;
        }
    }

    for line in conect_lines(&new_serial, &topology.explicit_bonds) {
        write_padded(out, &line)?;
    }
    write_padded(out, "END")?;
    Ok(warnings)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pdb;

    #[test]
    fn atom_name_field_places_by_element_symbol_length() {
        assert_eq!(&atom_name_field("CA", Element::CARBON), b" CA ");
        assert_eq!(
            &atom_name_field("CA", Element::from_atomic_number(20).unwrap()),
            b"CA  "
        );
        assert_eq!(
            &atom_name_field("FE", Element::from_atomic_number(26).unwrap()),
            b"FE  "
        );
        assert_eq!(&atom_name_field("OXT", Element::OXYGEN), b" OXT");
        assert_eq!(&atom_name_field("HG11", Element::HYDROGEN), b"HG11");
    }

    #[test]
    fn chain_ids_prefer_literal_single_characters_and_flag_remaps() {
        let mut names = vv_core::Interner::new();
        let a = names.intern("A");
        let aa = names.intern("AA");
        let topology = Topology {
            chains: vec![
                vv_core::ChainRec {
                    residues: 0..0,
                    label_asym: a,
                    auth_asym: a,
                    entity: 1,
                },
                vv_core::ChainRec {
                    residues: 0..0,
                    label_asym: aa,
                    auth_asym: aa,
                    entity: 1,
                },
            ],
            names,
            ..Default::default()
        };
        let (ids, warnings) = assign_chain_ids(&topology, |_| true);
        assert_eq!(ids[0], b'A');
        assert_eq!(ids[1], b'B');
        assert_eq!(warnings.len(), 1);
        assert!(warnings[0].contains("AA"));
        let (ids, warnings) = assign_chain_ids(&topology, |c| c == 0);
        assert_eq!(
            (ids[0], warnings.len()),
            (b'A', 0),
            "an excluded chain never warns"
        );
    }

    #[test]
    fn number_field_switches_to_hybrid36_past_width() {
        assert_eq!(number_field(99999, 5), "99999");
        assert_eq!(number_field(100000, 5), "A0000");
        assert_eq!(
            pdb::hybrid36(number_field(100000, 5).as_bytes(), 5),
            Some(100000)
        );
    }
}
