//! Kabat, Chothia and Martin numbering must equal the schemes' authors' own
//! renumbering of the non-redundant Fv set (`fixtures/real/reference/`, see
//! docs/ANTIBODY.md; the scheme authors' antibody structure database,
//! Ferdous & Martin 2018, Database, http://www.abybank.org/abdb/Data,
//! `NR_LH_Combined_*.tar.bz2`) on every domain, except the ones listed in
//! [`EXCEPTIONS`], and must not vary from run to run or with thread count.
//!
//! `cargo test --release -p vv-core --test antibody_reference_exact -- --ignored --nocapture`

use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

use rayon::prelude::*;
use vv_core::antibody::{find_domains, Label, Scheme};

mod antibody_common;
use antibody_common::{fingerprint, fixtures, in_pool, read_chains, Chain, Fingerprint};

const SEED_IDS: &str = include_str!("../src/antibody/seed_ids.txt");

/// Why a reference domain may differ from ours.
#[derive(Clone, Copy)]
enum Evidence {
    /// The reference treats the chain as an antigen (`REMARK 950` type `A`)
    /// and keeps the deposited numbers: no scheme was applied.
    Antigen,
    /// Another domain of the same scheme has this exact sequence and is
    /// numbered differently.
    SameDomain(&'static str),
    /// Another domain holds the same residues (given) and numbers them
    /// differently: the flanking chain, not the window, decided.
    SameWindow(&'static str, &'static str),
    /// No rule, and no twin in the set, accounts for it; the string names
    /// where the numbering departs. Not shown to be irreducible.
    Open(&'static str),
}

struct Known {
    id: &'static str,
    /// Schemes it differs in: `K`abat, `C`hothia, `M`artin.
    schemes: &'static str,
    evidence: Evidence,
}

const fn known(id: &'static str, schemes: &'static str, evidence: Evidence) -> Known {
    Known {
        id,
        schemes,
        evidence,
    }
}

use Evidence::{Antigen, Open};

const FR2: &str = "heavy FR2 deletion site";
const FR3: &str = "kappa FR3 gap site";
const START: &str = "N-terminal start";
const LOOP: &str = "CDR/framework boundary";

const EXCEPTIONS: &[Known] = &[
    known("1CIC_1:B", "KCM", Antigen),
    known("1CIC_2:D", "KCM", Antigen),
    known("1DVF_1:B", "KCM", Antigen),
    known("1IAI_2:I", "KCM", Antigen),
    known("1IAI_2:M", "KCM", Antigen),
    known("5XAJ_2:D", "KCM", Antigen),
    known("5XAJ_2:F", "KCM", Antigen),
    known("5WOB_4:Q", "KCM", Antigen),
    known("4K7P_2:X", "CM", Antigen),
    known("4K7P_2:Y", "M", Antigen),
    known("4XCF_1:H", "KCM", Evidence::SameDomain("4XAW_1:H")),
    known(
        "5CEY_1:L",
        "KCM",
        Evidence::SameWindow("5CEY_2:L", "YVRPLSVA"),
    ),
    known(
        "6NNJ_1:L",
        "KC",
        Evidence::SameWindow("5CEY_2:L", "YVRPLSVA"),
    ),
    known("1MFE_1:H", "KCM", Open(FR2)),
    known("4LLV_3:H", "KCM", Open(FR2)),
    known("3UTZ_1:L", "KCM", Open(FR3)),
    known("4LLV_3:L", "KCM", Open(FR3)),
    known("5EOC_2:L", "KCM", Open(FR3)),
    known("5VTA_2:L", "KCM", Open(FR3)),
    known("6BPC_1:L", "KCM", Open(FR3)),
    known("1OAY_2:L", "KC", Open(START)),
    known("3GK8_1:H", "KCM", Open(START)),
    known("1QFW_1:H", "KCM", Open(LOOP)),
    known("1QFW_1:L", "KCM", Open(LOOP)),
    known("4YDL_1:H", "KCM", Open(LOOP)),
];

const SCHEMES: [(&str, Scheme, char); 3] = [
    ("kabat", Scheme::Kabat, 'K'),
    ("chothia", Scheme::Chothia, 'C'),
    ("martin", Scheme::Martin, 'M'),
];

struct Record {
    /// `<file stem>:<chain>`, e.g. `1CIC_1:B`.
    id: String,
    file: PathBuf,
    /// Four-character PDB id of the file.
    entry: String,
    sequence: String,
    ours: Vec<Label>,
    reference: Vec<Label>,
}

impl Record {
    fn agrees(&self) -> bool {
        self.ours == self.reference
    }
}

/// Every non-seed domain of one scheme's reference files.
fn records(dir: &str, scheme: Scheme) -> Vec<Record> {
    let seeds: BTreeSet<&str> = SEED_IDS.split_whitespace().collect();
    let mut out = Vec::new();
    for (stem, path) in fixtures(&format!("reference/{dir}")) {
        let entry = stem[..4].to_ascii_uppercase();
        if seeds.contains(entry.as_str()) {
            continue;
        }
        for chain in read_chains(&path) {
            let sequence = chain.sequence();
            for d in find_domains(&sequence) {
                out.push(Record {
                    id: format!("{stem}:{}", chain.id),
                    file: path.clone(),
                    entry: entry.clone(),
                    sequence: sequence[d.start..d.end].to_string(),
                    ours: d.numbering(scheme).iter().map(|(_, l)| *l).collect(),
                    reference: chain.author(&d),
                });
            }
        }
    }
    out
}

/// Chain type letters of `REMARK 950 CHAIN <type> <label> <original>`.
fn chain_type(file: &Path, label: char) -> Option<char> {
    std::fs::read_to_string(file)
        .ok()?
        .lines()
        .filter_map(|l| l.strip_prefix("REMARK 950 CHAIN "))
        .map(|l| l.split_whitespace().collect::<Vec<_>>())
        .find(|f| f.get(1).is_some_and(|c| c.starts_with(label)))
        .and_then(|f| f[0].chars().next())
}

fn check_evidence(recs: &[Record], record: &Record, evidence: Evidence) {
    let other = |id: &str| {
        recs.iter()
            .find(|r| r.id == id)
            .unwrap_or_else(|| panic!("{id} missing from the reference set"))
    };
    match evidence {
        Evidence::Antigen => {
            let label = record.id.chars().last().expect("chain id");
            assert_eq!(chain_type(&record.file, label), Some('A'), "{}", record.id);
        }
        Evidence::SameDomain(id) => {
            let twin = other(id);
            assert_eq!(twin.sequence, record.sequence, "{id} is not a twin");
            assert_ne!(twin.reference, record.reference, "{id} agrees with it");
        }
        Evidence::SameWindow(id, window) => {
            let twin = other(id);
            let at = |r: &Record| r.sequence.find(window).expect("window in both");
            let (a, b) = (at(record), at(twin));
            let span = window.len();
            assert_ne!(
                record.reference[a..a + span],
                twin.reference[b..b + span],
                "{id} numbers the window alike"
            );
        }
        Evidence::Open(cause) => assert!(!cause.is_empty()),
    }
}

/// Domains of the even-indexed and odd-indexed entries (PDB ids in sorted
/// order), as `(domains, identical)` each.
fn parity_split(recs: &[Record]) -> [(usize, usize); 2] {
    let entries: BTreeSet<&str> = recs.iter().map(|r| r.entry.as_str()).collect();
    let half: BTreeMap<&str, usize> = entries
        .into_iter()
        .zip([0, 1].into_iter().cycle())
        .collect();
    let mut out = [(0, 0); 2];
    for r in recs {
        let slot = &mut out[half[r.entry.as_str()]];
        slot.0 += 1;
        slot.1 += usize::from(r.agrees());
    }
    out
}

#[test]
#[ignore = "needs fixtures/real/reference; use --release"]
fn reference_numbering_is_exact_up_to_documented_exceptions() {
    for (dir, scheme, letter) in SCHEMES {
        let recs = records(dir, scheme);
        if recs.is_empty() {
            println!("{dir}: no reference fixtures, skipped");
            return;
        }
        let listed: BTreeMap<&str, &Known> = EXCEPTIONS
            .iter()
            .filter(|k| k.schemes.contains(letter))
            .map(|k| (k.id, k))
            .collect();
        let stray: Vec<&str> = recs
            .iter()
            .filter(|r| !r.agrees() && !listed.contains_key(r.id.as_str()))
            .map(|r| r.id.as_str())
            .collect();
        assert!(stray.is_empty(), "{dir}: unlisted differences {stray:?}");
        for (id, k) in &listed {
            let record = recs.iter().find(|r| r.id == *id).expect("listed domain");
            assert!(
                !record.agrees(),
                "{dir}: {id} agrees now, drop it from EXCEPTIONS"
            );
            check_evidence(&recs, record, k.evidence);
        }
        let [even, odd] = parity_split(&recs);
        let identical = recs.iter().filter(|r| r.agrees()).count();
        println!(
            "{dir}: {identical}/{} identical (even half {}/{}, odd half {}/{}); {} listed exceptions",
            recs.len(),
            even.1,
            even.0,
            odd.1,
            odd.0,
            listed.len()
        );
    }
}

fn sequences_of(dirs: &[&str]) -> Vec<String> {
    dirs.iter()
        .flat_map(|dir| fixtures(dir))
        .flat_map(|(_, path)| {
            read_chains(&path)
                .iter()
                .map(Chain::sequence)
                .collect::<Vec<_>>()
        })
        .collect()
}

fn fingerprints(sequences: &[String]) -> Vec<Fingerprint> {
    sequences.par_iter().map(|s| fingerprint(s)).collect()
}

#[test]
#[ignore = "needs fixtures/real; use --release"]
fn real_chains_number_identically_every_run_and_thread_count() {
    let sequences = sequences_of(&["reference/kabat", "holdout"]);
    if sequences.is_empty() {
        return;
    }
    let first: Vec<Fingerprint> = sequences.iter().map(|s| fingerprint(s)).collect();
    let again: Vec<Fingerprint> = sequences.iter().map(|s| fingerprint(s)).collect();
    assert!(first == again, "second serial run differs");
    for threads in [1, 4] {
        let pooled = in_pool(threads, || fingerprints(&sequences));
        assert!(
            first == pooled,
            "{threads} threads differ from the serial run"
        );
    }
}
