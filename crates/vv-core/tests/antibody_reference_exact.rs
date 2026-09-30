//! Kabat, Chothia and Martin numbering must equal the schemes' authors' own
//! renumbering of the non-redundant Fv set (`fixtures/real/reference/`, see
//! docs/ANTIBODY.md; the scheme authors' antibody structure database,
//! Ferdous & Martin 2018, Database, http://www.abybank.org/abdb/Data,
//! `NR_LH_Combined_*.tar.bz2`) on every domain, except the ones listed in
//! [`EXCEPTIONS`], and must not vary from run to run or with thread count.
//!
//! `cargo test --release -p vv-core --test antibody_reference_exact -- --ignored --nocapture`

use std::collections::{BTreeMap, BTreeSet};

use rayon::prelude::*;
use vv_core::antibody::{find_domains, Label, Scheme};

mod antibody_common;
use antibody_common::{fingerprint, fixtures, in_pool, read_chains, Chain, Fingerprint};

const SEED_IDS: &str = include_str!("../src/antibody/seed_ids.txt");

/// Why a reference domain may differ from ours.
#[derive(Clone, Copy)]
enum Evidence {
    /// The reference wrote plain consecutive numbers: no scheme was applied
    /// (its Cys 92 or 88 lands on 95, 96 or 106).
    Sequential,
    /// Another domain of the same scheme has this exact sequence and is
    /// numbered differently.
    SameDomain(&'static str),
    /// Another domain starts with the same `n` residues and is numbered
    /// differently over them.
    SamePrefix(&'static str, usize),
    /// A gap or chain start the reference placed where no length rule, and
    /// no twin in the set, accounts for it. Not shown to be irreducible.
    Unexplained,
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

use Evidence::{Sequential, Unexplained};

const EXCEPTIONS: &[Known] = &[
    known("1CIC_1:B", "KCM", Sequential),
    known("1CIC_2:D", "KCM", Sequential),
    known("1DVF_1:B", "KCM", Sequential),
    known("1IAI_2:I", "KCM", Sequential),
    known("1IAI_2:M", "KCM", Sequential),
    known("5XAJ_2:D", "KCM", Sequential),
    known("5XAJ_2:F", "KCM", Sequential),
    known("1T2Q_1:H", "KCM", Evidence::SameDomain("2D03_1:H")),
    known("4G6A_1:H", "KCM", Evidence::SameDomain("4G6A_2:H")),
    known("4R26_1:L", "KC", Evidence::SameDomain("6MCO_1:L")),
    known("6MCO_1:L", "KC", Evidence::SameDomain("4R26_1:L")),
    known("4XCF_1:H", "KCM", Evidence::SameDomain("4XAW_1:H")),
    known("4K7P_2:X", "CM", Evidence::SameDomain("4K7P_1:L")),
    known("4K7P_2:Y", "M", Evidence::SameDomain("4K7P_1:H")),
    known("3U6R_1:H", "KCM", Evidence::SamePrefix("1R70_1:H", 12)),
    known("4N0Y_1:H", "KCM", Evidence::SamePrefix("1R70_1:H", 9)),
    known("2HH0_1:H", "KCM", Evidence::SamePrefix("1R70_1:H", 5)),
    known("1MFE_1:H", "KCM", Unexplained),
    known("1OAY_2:L", "KC", Unexplained),
    known("1QFW_1:H", "KCM", Unexplained),
    known("1QFW_1:L", "KCM", Unexplained),
    known("3UTZ_1:L", "KCM", Unexplained),
    known("4JY6_1:L", "KCM", Unexplained),
    known("4LLV_3:H", "KCM", Unexplained),
    known("4LLV_3:L", "KCM", Unexplained),
    known("4UOM_1:L", "KC", Unexplained),
    known("4YDL_1:H", "KCM", Unexplained),
    known("5EOC_2:L", "KCM", Unexplained),
    known("5FYL_1:L", "KC", Unexplained),
    known("5VTA_2:L", "KCM", Unexplained),
    known("5WB9_1:L", "KCM", Unexplained),
    known("5WOB_4:Q", "KCM", Unexplained),
    known("6AOD_1:L", "KC", Unexplained),
    known("6BPC_1:L", "KCM", Unexplained),
];

const SCHEMES: [(&str, Scheme, char); 3] = [
    ("kabat", Scheme::Kabat, 'K'),
    ("chothia", Scheme::Chothia, 'C'),
    ("martin", Scheme::Martin, 'M'),
];

struct Record {
    /// `<file stem>:<chain>`, e.g. `1CIC_1:B`.
    id: String,
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

fn check_evidence(recs: &[Record], record: &Record, evidence: Evidence) {
    let other = |id: &str| {
        recs.iter()
            .find(|r| r.id == id)
            .unwrap_or_else(|| panic!("{id} missing from the reference set"))
    };
    match evidence {
        Evidence::Sequential => {
            let numbers: Vec<u16> = record.reference.iter().map(|l| l.number).collect();
            let plain = record.reference.iter().all(|l| l.insertion().is_none());
            assert!(plain && numbers.windows(2).all(|w| w[1] == w[0] + 1));
        }
        Evidence::SameDomain(id) => {
            let twin = other(id);
            assert_eq!(twin.sequence, record.sequence, "{id} is not a twin");
            assert_ne!(twin.reference, record.reference, "{id} agrees with it");
        }
        Evidence::SamePrefix(id, n) => {
            let twin = other(id);
            assert_eq!(twin.sequence[..n], record.sequence[..n], "{id} prefix");
            assert_ne!(twin.reference[..n], record.reference[..n], "{id} agrees");
        }
        Evidence::Unexplained => {}
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
