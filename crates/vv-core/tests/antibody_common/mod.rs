//! Shared reader for the git-ignored `fixtures/real/` PDB entries.

#![allow(dead_code)]

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use vv_core::antibody::{find_domains, ChainType, Domain, Label, Scheme};

pub struct Residue {
    pub label: Label,
    pub aa: char,
}

pub struct Chain {
    pub id: char,
    pub residues: Vec<Residue>,
}

impl Chain {
    pub fn sequence(&self) -> String {
        self.residues.iter().map(|r| r.aa).collect()
    }

    pub fn author(&self, d: &Domain) -> Vec<Label> {
        self.residues[d.start..d.end]
            .iter()
            .map(|r| r.label)
            .collect()
    }
}

pub fn one_letter(name: &str) -> char {
    const AA: [(&str, char); 22] = [
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
        ("PCA", 'E'),
    ];
    AA.iter().find(|(n, _)| *n == name).map_or('X', |(_, c)| *c)
}

/// Polymer chains of a PDB file from their CA atoms (first altloc only).
pub fn read_chains(path: &Path) -> Vec<Chain> {
    let text = std::fs::read_to_string(path).unwrap_or_default();
    let mut chains: BTreeMap<char, Vec<Residue>> = BTreeMap::new();
    for line in text
        .lines()
        .filter(|l| l.len() > 26 && (l.starts_with("ATOM") || l.starts_with("HETATM")))
    {
        let (name, res) = (&line[12..16], &line[17..20]);
        let hetero_ok = !line.starts_with("HETATM") || matches!(res, "MSE" | "PCA");
        if name != " CA " || !hetero_ok || !matches!(line.as_bytes()[16], b' ' | b'A') {
            continue;
        }
        let number: u16 = line[22..26].trim().parse().unwrap_or(0);
        let label = match line.as_bytes()[26] {
            b' ' => Label::new(number),
            c => Label::with_insertion(number, c as char),
        };
        let chain = chains.entry(line.as_bytes()[21] as char).or_default();
        if chain.last().is_none_or(|r| r.label != label) {
            chain.push(Residue {
                label,
                aa: one_letter(res),
            });
        }
    }
    chains
        .into_iter()
        .map(|(id, residues)| Chain { id, residues })
        .collect()
}

/// `fixtures/real`, or `$VIZVIZ_FIXTURES_REAL` (worktrees keep the data in the main checkout).
pub fn fixtures_root() -> PathBuf {
    std::env::var_os("VIZVIZ_FIXTURES_REAL").map_or_else(
        || Path::new(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/real"),
        PathBuf::from,
    )
}

pub fn fixtures(dir: &str) -> Vec<(String, PathBuf)> {
    let root = fixtures_root().join(dir);
    let mut files: Vec<(String, PathBuf)> = std::fs::read_dir(root)
        .map(|entries| {
            entries
                .filter_map(Result::ok)
                .map(|e| e.path())
                .filter(|p| p.extension().is_some_and(|e| e == "pdb"))
                .map(|p| (p.file_stem().unwrap().to_string_lossy().to_string(), p))
                .collect()
        })
        .unwrap_or_default();
    files.sort();
    files
}

/// Everything `find_domains` reports for a sequence, scores as raw bits so
/// that two runs compare exactly.
pub type Fingerprint = Vec<(ChainType, usize, usize, u32, Vec<Vec<(usize, Label)>>)>;

pub fn fingerprint(sequence: &str) -> Fingerprint {
    find_domains(sequence)
        .iter()
        .map(|d| {
            let numbered = Scheme::ALL.iter().map(|&s| d.numbering(s)).collect();
            (d.chain, d.start, d.end, d.score.to_bits(), numbered)
        })
        .collect()
}

/// Runs `work` on a rayon pool of exactly `threads` threads.
pub fn in_pool<T: Send>(threads: usize, work: impl FnOnce() -> T + Send) -> T {
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .expect("thread pool")
        .install(work)
}

/// Polymer name of each chain from the `COMPND` records.
pub fn molecule_names(path: &Path) -> BTreeMap<char, String> {
    let text = std::fs::read_to_string(path).unwrap_or_default();
    let mut names = BTreeMap::new();
    let mut current = String::new();
    for line in text.lines().take_while(|l| !l.starts_with("ATOM")) {
        let Some(body) = line
            .strip_prefix("COMPND")
            .map(|b| b.get(4..).unwrap_or("").trim())
        else {
            continue;
        };
        if let Some(name) = body.strip_prefix("MOLECULE:") {
            current = name.trim().trim_end_matches(';').to_string();
        } else if let Some(ids) = body.strip_prefix("CHAIN:") {
            for id in ids.trim().trim_end_matches(';').split(',') {
                if let Some(c) = id.trim().chars().next() {
                    names.insert(c, current.clone());
                }
            }
        }
    }
    names
}
