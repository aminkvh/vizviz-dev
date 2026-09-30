//! `COMPND` and `SOURCE` (wwPDB v3.3): per polymer, the molecule name, the
//! chains that carry it and its source organism, from the entity notes the
//! file was read with.

use crate::entity_notes::EntityNotes;
use crate::polymer_layout::Layout;
use vv_core::Topology;

/// Characters of specification text a record line holds after its prefix.
const LINE_TEXT: usize = 68;

/// One record's specification, wrapped at word boundaries into lines with
/// the continuation number in columns 9-10.
fn wrap(record: &str, spec: &str) -> Vec<String> {
    let mut lines: Vec<String> = Vec::new();
    for word in spec.split(' ') {
        match lines.last_mut() {
            Some(line) if line.len() + 1 + word.len() <= LINE_TEXT => {
                line.push(' ');
                line.push_str(word);
            }
            _ => lines.push(word.to_string()),
        }
    }
    lines
        .iter()
        .enumerate()
        .map(|(i, text)| match i {
            0 => format!("{record:<6}    {text}"),
            _ => format!("{record:<6}  {:>2} {text}", i + 1),
        })
        .collect()
}

fn chains_of(layout: &Layout, entity: usize, chain_ids: &[u8]) -> String {
    let mut ids: Vec<char> = Vec::new();
    for &s in &layout.entities[entity].segments {
        let id = chain_ids[layout.segments[s].chain as usize] as char;
        if !ids.contains(&id) {
            ids.push(id);
        }
    }
    ids.iter()
        .map(char::to_string)
        .collect::<Vec<_>>()
        .join(", ")
}

/// `COMPND` then `SOURCE` lines for the polymers that have a name or an
/// organism; `MOL_ID`s count those polymers from 1.
pub(crate) fn compnd_source_lines(t: &Topology, layout: &Layout, chain_ids: &[u8]) -> Vec<String> {
    let notes = EntityNotes::new(t, layout);
    let (mut compnd, mut source) = (Vec::new(), Vec::new());
    let mut mol_id = 0;
    for (i, entity) in layout.entities.iter().enumerate() {
        let name = notes.entity_value(i, "pdbx_description");
        let organism = notes.organism(i);
        if !entity.kind.is_polymer() || (name.is_none() && organism.is_none()) {
            continue;
        }
        mol_id += 1;
        let mut spec = vec![format!("MOL_ID: {mol_id};")];
        if let Some(name) = name {
            spec.push(format!("MOLECULE: {name};"));
        }
        spec.push(format!("CHAIN: {};", chains_of(layout, i, chain_ids)));
        compnd.extend(wrap("COMPND", &spec.join(" ")));
        if let Some(organism) = organism {
            let spec = format!("MOL_ID: {mol_id}; ORGANISM_SCIENTIFIC: {organism};");
            source.extend(wrap("SOURCE", &spec));
        }
    }
    compnd.extend(source);
    compnd
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_long_specification_wraps_with_numbered_continuations() {
        let spec = format!("MOL_ID: 1; MOLECULE: {};", "WORD ".repeat(30).trim_end());
        let lines = wrap("COMPND", &spec);
        assert!(lines.len() > 1);
        assert!(lines[0].starts_with("COMPND    MOL_ID: 1;"));
        assert!(lines[1].starts_with("COMPND   2 "));
        assert!(lines.iter().all(|l| l.len() <= 80));
    }
}
