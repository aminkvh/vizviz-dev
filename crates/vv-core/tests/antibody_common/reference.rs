//! The scheme authors' own renumbering (`fixtures/real/reference/`) as
//! records, the comparison set drawn from it, and its even/odd halves.

#![allow(dead_code)]

use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

use vv_core::antibody::{find_domains, ChainType, Domain, Label, Scheme};

use super::{fixtures, read_chains};

const SEED_IDS: &str = include_str!("../../src/antibody/seed_ids.txt");

pub const SCHEMES: [(&str, Scheme, char); 3] = [
    ("kabat", Scheme::Kabat, 'K'),
    ("chothia", Scheme::Chothia, 'C'),
    ("martin", Scheme::Martin, 'M'),
];

pub struct Record {
    /// `<file stem>:<chain>`, e.g. `1CIC_1:B`.
    pub id: String,
    pub file: PathBuf,
    /// Four-character PDB id of the file.
    pub entry: String,
    pub chain: ChainType,
    pub domain: Domain,
    pub sequence: String,
    pub chain_sequence: String,
    pub ours: Vec<Label>,
    pub reference: Vec<Label>,
    /// `REMARK 950` types the chain as an antigen.
    pub antigen: bool,
}

impl Record {
    pub fn agrees(&self) -> bool {
        self.ours == self.reference
    }
}

fn seed_entries() -> BTreeSet<&'static str> {
    SEED_IDS.split_whitespace().collect()
}

/// Chain type letter of each chain label in `REMARK 950 CHAIN <type> <label> <original>`.
fn chain_types(file: &Path) -> BTreeMap<char, char> {
    let text = std::fs::read_to_string(file).unwrap_or_default();
    text.lines()
        .filter_map(|l| l.strip_prefix("REMARK 950 CHAIN "))
        .filter_map(|l| {
            let f: Vec<&str> = l.split_whitespace().collect();
            Some((f.get(1)?.chars().next()?, f[0].chars().next()?))
        })
        .collect()
}

/// Every non-seed domain of one scheme's reference files.
pub fn records(dir: &str, scheme: Scheme) -> Vec<Record> {
    let seeds = seed_entries();
    let mut out = Vec::new();
    for (stem, path) in fixtures(&format!("reference/{dir}")) {
        let entry = stem[..4].to_ascii_uppercase();
        if seeds.contains(entry.as_str()) {
            continue;
        }
        let types = chain_types(&path);
        for chain in read_chains(&path) {
            let sequence = chain.sequence();
            for d in find_domains(&sequence) {
                out.push(Record {
                    id: format!("{stem}:{}", chain.id),
                    file: path.clone(),
                    entry: entry.clone(),
                    chain: d.chain,
                    domain: d.clone(),
                    sequence: sequence[d.start..d.end].to_string(),
                    chain_sequence: sequence.clone(),
                    ours: d.numbering(scheme).iter().map(|(_, l)| *l).collect(),
                    reference: chain.author(&d),
                    antigen: types.get(&chain.id) == Some(&'A'),
                });
            }
        }
    }
    out
}

/// 0 for the even-indexed entries of the sorted non-seed PDB ids of the
/// three reference bundles together, 1 for the odd-indexed.
pub fn halves() -> BTreeMap<String, usize> {
    let seeds = seed_entries();
    let entries: BTreeSet<String> = SCHEMES
        .iter()
        .flat_map(|(dir, _, _)| fixtures(&format!("reference/{dir}")))
        .map(|(stem, _)| stem[..4].to_ascii_uppercase())
        .filter(|e| !seeds.contains(e.as_str()))
        .collect();
    entries
        .into_iter()
        .zip([0, 1].into_iter().cycle())
        .collect()
}

/// Domains the reference numbered twice with different results although
/// the program saw the same chain: identical chain sequence and domain
/// sequence. The numbering most copies follow is kept; with no majority all
/// copies are dropped. Keyed by record id.
pub fn inconsistent_twins(recs: &[Record]) -> BTreeMap<String, String> {
    let mut groups: BTreeMap<(&str, &str), Vec<&Record>> = BTreeMap::new();
    for r in recs.iter().filter(|r| !r.antigen) {
        groups
            .entry((&r.chain_sequence, &r.sequence))
            .or_default()
            .push(r);
    }
    let mut out = BTreeMap::new();
    for group in groups.values() {
        let mut votes: BTreeMap<&[Label], usize> = BTreeMap::new();
        for r in group {
            *votes.entry(&r.reference).or_default() += 1;
        }
        if votes.len() < 2 {
            continue;
        }
        let top = *votes.values().max().expect("votes");
        let leaders: Vec<&[Label]> = votes
            .iter()
            .filter(|(_, n)| **n == top)
            .map(|(l, _)| *l)
            .collect();
        for r in group {
            let loses = leaders.len() > 1 || !leaders.contains(&r.reference.as_slice());
            if loses {
                let note = format!(
                    "{} of {} identical copies agree with it",
                    votes[r.reference.as_slice()],
                    group.len()
                );
                out.insert(r.id.clone(), note);
            }
        }
    }
    out
}

/// Indices into `recs` of the comparison set: domains the reference
/// scheme-numbered and numbered consistently.
pub fn comparison_set(recs: &[Record]) -> Vec<usize> {
    let twins = inconsistent_twins(recs);
    (0..recs.len())
        .filter(|&i| !recs[i].antigen && !twins.contains_key(&recs[i].id))
        .collect()
}
