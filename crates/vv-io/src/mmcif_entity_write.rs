//! The entity tables of the mmCIF writer: `_entity`, `_entity_poly`,
//! `_pdbx_entity_nonpoly` and `_chem_comp`, which `mmcif_entity` reads back
//! as each residue's polymer status.

use std::collections::HashMap;
use std::io::{self, Write};

use vv_core::{InternId, Topology};

use crate::entity_notes::{EntityNotes, ENTITY_TABLES};
use crate::mmcif_write::{cell, token};
use crate::polymer_layout::{Entity, Kind, Layout};

fn entity_type(kind: Kind) -> &'static str {
    match kind {
        Kind::Protein | Kind::Nucleic => "polymer",
        Kind::Branched => "branched",
        Kind::NonPolymer => "non-polymer",
        Kind::Water => "water",
    }
}

/// A deoxynucleotide name: `DA`, `DC`, `DG`, `DT`, `DU`, `DI`.
fn is_deoxy(comp: &str) -> bool {
    let b = comp.as_bytes();
    b.len() == 2 && b[0] == b'D' && b"ACGTUI".contains(&b[1])
}

fn poly_type(t: &Topology, entity: &Entity) -> &'static str {
    if entity.kind == Kind::Protein {
        return "polypeptide(L)";
    }
    let deoxy = entity
        .comps
        .iter()
        .filter(|&&c| is_deoxy(t.names.get(c)))
        .count();
    match deoxy {
        0 => "polyribonucleotide",
        n if n == entity.comps.len() => "polydeoxyribonucleotide",
        _ => "polydeoxyribonucleotide/polyribonucleotide hybrid",
    }
}

/// One-letter code of a monomer, the parent letter for a common variant,
/// `X` (protein) or `N` (nucleic) otherwise.
fn one_letter(comp: &str, kind: Kind) -> char {
    let upper = comp.to_ascii_uppercase();
    if kind == Kind::Protein {
        const CODES: [(&str, char); 24] = [
            ("ALA", 'A'),
            ("ARG", 'R'),
            ("ASN", 'N'),
            ("ASP", 'D'),
            ("CYS", 'C'),
            ("GLN", 'Q'),
            ("GLU", 'E'),
            ("GLY", 'G'),
            ("HIS", 'H'),
            ("ILE", 'I'),
            ("LEU", 'L'),
            ("LYS", 'K'),
            ("MET", 'M'),
            ("PHE", 'F'),
            ("PRO", 'P'),
            ("SER", 'S'),
            ("THR", 'T'),
            ("TRP", 'W'),
            ("TYR", 'Y'),
            ("VAL", 'V'),
            ("MSE", 'M'),
            ("SEC", 'U'),
            ("PYL", 'O'),
            ("UNK", 'X'),
        ];
        return CODES
            .iter()
            .find(|(n, _)| *n == upper)
            .map_or('X', |&(_, c)| c);
    }
    let base = upper
        .strip_prefix('D')
        .filter(|_| is_deoxy(&upper))
        .unwrap_or(&upper);
    match base {
        "A" | "C" | "G" | "T" | "U" | "I" => base.chars().next().expect("non-empty"),
        _ => 'N',
    }
}

fn strand_ids(t: &Topology, layout: &Layout, entity: &Entity) -> String {
    let mut ids: Vec<&str> = Vec::new();
    for &s in &entity.segments {
        let chain = &t.chains[layout.segments[s].chain as usize];
        let id = t.names.get(chain.auth_asym);
        if !ids.contains(&id) {
            ids.push(id);
        }
    }
    ids.join(",")
}

fn loop_header(out: &mut impl Write, category: &str, items: &[&str]) -> io::Result<()> {
    writeln!(out, "loop_")?;
    for item in items {
        writeln!(out, "_{category}.{item}")?;
    }
    Ok(())
}

/// `_entity` items about the molecule itself, carried over when the file
/// gave them (counts and weights would go stale with an edited structure).
const CARRIED_ENTITY_ITEMS: [&str; 4] =
    ["pdbx_description", "src_method", "pdbx_ec", "pdbx_mutation"];

fn write_entity(layout: &Layout, notes: &EntityNotes, out: &mut impl Write) -> io::Result<()> {
    let count = layout.entities.len();
    let carried: Vec<&str> = CARRIED_ENTITY_ITEMS
        .into_iter()
        .filter(|item| (0..count).any(|e| notes.entity_value(e, item).is_some()))
        .collect();
    let items: Vec<&str> = ["id", "type"].into_iter().chain(carried.clone()).collect();
    loop_header(out, "entity", &items)?;
    for (i, e) in layout.entities.iter().enumerate() {
        write!(out, "{} {}", i + 1, entity_type(e.kind))?;
        for item in &carried {
            write!(out, " {}", cell(notes.entity_value(i, item).unwrap_or("")))?;
        }
        writeln!(out)?;
    }
    writeln!(out, "#")
}

/// The file's source and branch tables, rows renumbered to the written
/// entities they describe.
fn write_entity_sources(
    layout: &Layout,
    notes: &EntityNotes,
    out: &mut impl Write,
) -> io::Result<()> {
    for category in ENTITY_TABLES {
        let mut header_written = false;
        for entity in 0..layout.entities.len() {
            let Some((cat, rows)) = notes.rows(category, entity) else {
                continue;
            };
            let id_col = cat.column("entity_id").expect("rows matched on it");
            for row in rows {
                if !header_written {
                    let items: Vec<&str> = cat.items.iter().map(String::as_str).collect();
                    loop_header(out, category, &items)?;
                    header_written = true;
                }
                let cells: Vec<String> = cat.rows[row]
                    .iter()
                    .enumerate()
                    .map(|(c, v)| match c == id_col {
                        true => (entity + 1).to_string(),
                        false => cell(v),
                    })
                    .collect();
                writeln!(out, "{}", cells.join(" "))?;
            }
        }
        if header_written {
            writeln!(out, "#")?;
        }
    }
    Ok(())
}

fn write_entity_poly(t: &Topology, layout: &Layout, out: &mut impl Write) -> io::Result<()> {
    let polymers: Vec<_> = layout
        .entities
        .iter()
        .enumerate()
        .filter(|(_, e)| e.kind.is_polymer())
        .collect();
    if polymers.is_empty() {
        return Ok(());
    }
    let items = [
        "entity_id",
        "type",
        "pdbx_seq_one_letter_code_can",
        "pdbx_strand_id",
    ];
    loop_header(out, "entity_poly", &items)?;
    for (i, e) in polymers {
        let code: String = e
            .comps
            .iter()
            .map(|&c| one_letter(t.names.get(c), e.kind))
            .collect();
        let strands = token(&strand_ids(t, layout, e));
        writeln!(out, "{} '{}' {code} {strands}", i + 1, poly_type(t, e))?;
    }
    writeln!(out, "#")
}

fn write_entity_nonpoly(
    t: &Topology,
    layout: &Layout,
    notes: &EntityNotes,
    out: &mut impl Write,
) -> io::Result<()> {
    let single: Vec<_> = layout
        .entities
        .iter()
        .enumerate()
        .filter(|(_, e)| matches!(e.kind, Kind::NonPolymer | Kind::Water))
        .collect();
    if single.is_empty() {
        return Ok(());
    }
    loop_header(
        out,
        "pdbx_entity_nonpoly",
        &["entity_id", "name", "comp_id"],
    )?;
    for (i, e) in single {
        let comp = token(t.names.get(e.comps[0]));
        let name = match (e.kind, notes.entity_value(i, "pdbx_description")) {
            (Kind::Water, _) => "water".to_string(),
            (_, Some(description)) => cell(description),
            _ => comp.clone(),
        };
        writeln!(out, "{} {name} {comp}", i + 1)?;
    }
    writeln!(out, "#")
}

fn chem_comp_type(t: &Topology, kind: Kind, comp: InternId) -> &'static str {
    match kind {
        Kind::Protein => "L-peptide linking",
        Kind::Nucleic if is_deoxy(t.names.get(comp)) => "DNA linking",
        Kind::Nucleic => "RNA linking",
        Kind::Branched => "saccharide",
        Kind::NonPolymer | Kind::Water => "non-polymer",
    }
}

/// Each residue name once, in file order; a name used both in a polymer
/// and as a free ligand takes the polymer type, as the dictionary has one
/// entry per component.
fn write_chem_comp(t: &Topology, layout: &Layout, out: &mut impl Write) -> io::Result<()> {
    let mut order: Vec<InternId> = Vec::new();
    let mut types: HashMap<InternId, (Kind, &'static str)> = HashMap::new();
    for segment in &layout.segments {
        let range = segment.residues.start as usize..segment.residues.end as usize;
        for res in &t.residues[range] {
            let kind = segment.kind;
            let ty = chem_comp_type(t, kind, res.comp);
            match types.get(&res.comp) {
                None => {
                    order.push(res.comp);
                    types.insert(res.comp, (kind, ty));
                }
                Some(&(seen, _)) if !seen.is_polymer() && kind.is_polymer() => {
                    types.insert(res.comp, (kind, ty));
                }
                Some(_) => {}
            }
        }
    }
    loop_header(out, "chem_comp", &["id", "type"])?;
    for comp in order {
        writeln!(out, "{} '{}'", token(t.names.get(comp)), types[&comp].1)?;
    }
    writeln!(out, "#")
}

pub(crate) fn write_entity_tables(
    t: &Topology,
    layout: &Layout,
    out: &mut impl Write,
) -> io::Result<()> {
    if layout.entities.is_empty() {
        return Ok(());
    }
    let notes = EntityNotes::new(t, layout);
    write_entity(layout, &notes, out)?;
    write_entity_sources(layout, &notes, out)?;
    write_entity_poly(t, layout, out)?;
    write_entity_nonpoly(t, layout, &notes, out)?;
    write_chem_comp(t, layout, out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn monomers_map_to_parent_letters() {
        assert_eq!(one_letter("ala", Kind::Protein), 'A');
        assert_eq!(one_letter("MSE", Kind::Protein), 'M');
        assert_eq!(one_letter("ZZZ", Kind::Protein), 'X');
        assert_eq!(one_letter("DT", Kind::Nucleic), 'T');
        assert_eq!(one_letter("G", Kind::Nucleic), 'G');
        assert_eq!(one_letter("5BU", Kind::Nucleic), 'N');
    }
}
