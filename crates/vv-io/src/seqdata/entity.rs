//! Each polymer chain's full sequence as a structure file states it: the
//! mmCIF `_pdbx_poly_seq_scheme` table (which also carries author numbers)
//! or the PDB `SEQRES` records. BinaryCIF is not read.

use std::collections::HashMap;
use std::io::Read;
use std::path::Path;

use vv_core::seqfeat::EntityChain;

use crate::cif;

/// Chains of `path` by chain id: the `label_asym_id` for mmCIF, the chain
/// character for PDB. `None` when the file has no such table or is in a
/// format this does not read.
pub fn read(path: &Path) -> Option<HashMap<String, EntityChain>> {
    let text = read_text(path)?;
    let name = path.file_name()?.to_str()?.to_ascii_lowercase();
    let name = name.strip_suffix(".gz").unwrap_or(&name);
    let chains = if name.ends_with(".pdb") || name.ends_with(".ent") {
        seqres(&text)
    } else if name.ends_with(".cif") || name.ends_with(".mmcif") || name.ends_with(".pdbx") {
        poly_seq_scheme(text.as_bytes())
    } else {
        return None;
    };
    (!chains.is_empty()).then_some(chains)
}

fn read_text(path: &Path) -> Option<String> {
    let bytes = std::fs::read(path).ok()?;
    if !bytes.starts_with(&[0x1f, 0x8b]) {
        return String::from_utf8(bytes).ok();
    }
    let mut text = String::new();
    flate2::read::GzDecoder::new(&bytes[..])
        .read_to_string(&mut text)
        .ok()?;
    Some(text)
}

fn poly_seq_scheme(src: &[u8]) -> HashMap<String, EntityChain> {
    let (cats, _) = cif::read_categories(src, &["atom_site"]);
    let mut chains: HashMap<String, EntityChain> = HashMap::new();
    let Some(cat) = cats.get("pdbx_poly_seq_scheme") else {
        return chains;
    };
    let mut last: HashMap<String, i32> = HashMap::new();
    for row in 0..cat.rows.len() {
        let (Some(asym), Some(mon), Some(seq)) = (
            cat.get_str(row, "asym_id"),
            cat.get_str(row, "mon_id"),
            cat.get_str(row, "seq_id")
                .and_then(|s| s.parse::<i32>().ok()),
        ) else {
            continue;
        };
        // Microheterogeneity repeats a position with another monomer.
        if last.insert(asym.to_string(), seq) == Some(seq) {
            continue;
        }
        let number = cat
            .get_str(row, "pdb_seq_num")
            .and_then(|s| s.parse::<i32>().ok())
            .map(|n| (n, cat.get(row, "pdb_ins_code").map_or(0, |s| s[0])));
        let chain = chains.entry(asym.to_string()).or_default();
        chain.names.push(mon.to_string());
        chain.numbers.push(number);
    }
    chains
}

/// `SEQRES` names sit in 13 columns of 3 characters from column 20.
fn seqres(text: &str) -> HashMap<String, EntityChain> {
    let mut chains: HashMap<String, EntityChain> = HashMap::new();
    for line in text.lines().filter(|l| l.starts_with("SEQRES")) {
        let Some(chain) = line.get(11..12) else {
            continue;
        };
        let entry = chains.entry(chain.to_string()).or_default();
        let mut column = 19;
        while let Some(name) = line.get(column..column + 3) {
            if !name.trim().is_empty() {
                entry.names.push(name.trim().to_string());
            }
            column += 4;
        }
    }
    chains
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture(name: &str) -> std::path::PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../fixtures/small")
            .join(name)
    }

    #[test]
    fn mmcif_gives_names_and_author_numbers() {
        let chains = read(&fixture("4HHB.cif")).unwrap();
        let a = &chains["A"];
        assert_eq!(a.names.len(), 141);
        assert_eq!(a.names[0], "VAL");
        assert_eq!(a.numbers[0], Some((1, 0)));
        assert_eq!(chains["B"].names.len(), 146);
    }

    #[test]
    fn pdb_seqres_gives_names_in_order() {
        let chains = read(&fixture("4HHB.pdb")).unwrap();
        assert_eq!(chains["A"].names.len(), 141);
        assert_eq!(chains["B"].names[..2], ["VAL", "HIS"]);
        assert!(chains["A"].numbers.is_empty());
    }

    #[test]
    fn a_file_without_a_table_or_of_another_format_gives_none() {
        assert!(read(&fixture("1CRN_traj.dcd")).is_none());
        assert!(read(&fixture("4HHB.bcif")).is_none());
        assert!(read(Path::new("no/such/file.pdb")).is_none());
    }
}
