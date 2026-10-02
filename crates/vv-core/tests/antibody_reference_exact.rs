//! Kabat, Chothia and Martin numbering against the schemes' authors' own
//! renumbering of the non-redundant Fv set (`fixtures/real/reference/`, see
//! docs/ANTIBODY.md; the scheme authors' antibody structure database,
//! Ferdous & Martin 2018, Database, http://www.abybank.org/abdb/Data,
//! `NR_LH_Combined_*.tar.bz2`).
//!
//! Two sets are reported. The raw set is every non-seed reference domain.
//! The comparison set drops the chains the reference types as antigen
//! (`REMARK 950` type `A`, deposited numbers kept) and the copies of a
//! domain that the reference numbers inconsistently (same chain and domain
//! sequence, different labels; majority kept, no majority drops all).
//! The test fails on any comparison-set difference not listed in
//! [`RESIDUAL`], on a listed one that starts to agree, and on any raw
//! difference outside the antigen, twin and residual groups. It also
//! requires numbering not to vary from run to run or with thread count.
//!
//! `cargo test --release -p vv-core --test antibody_reference_exact -- --ignored --nocapture`
//!
//! With `VIZVIZ_EVIDENCE=1` the first test also prints, for each residual
//! domain, the rule's, our and the reference's labels over the differing
//! stretch and every placement of that stretch's residues onto the base
//! labels between its ends, scored by the even-half consensus.

use std::collections::{BTreeMap, BTreeSet};

use rayon::prelude::*;
use vv_core::antibody::{ChainType, Label, Scheme};

mod antibody_common;
use antibody_common::reference::{
    comparison_set, halves, inconsistent_twins, records, Record, SCHEMES,
};
use antibody_common::{fingerprint, fixtures, in_pool, read_chains, Chain, Fingerprint};

/// Why a comparison-set domain may differ from ours.
#[derive(Clone, Copy)]
enum Evidence {
    /// Other domains hold the same residues (given) and the reference
    /// numbers them differently: the window does not decide the labels, so
    /// no rule over it can match both.
    SameWindow(&'static [&'static str], &'static str),
    /// The loop between these labels has one pattern per length across the
    /// reference set; this domain is the only one to deviate from it.
    LoopException(u16, u16),
    /// No identical window in the set and no rule found; the string names
    /// the site. The evidence output shows the consensus does not single
    /// out the reference's placement.
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

use Evidence::{LoopException, Open, SameWindow};

const FR3: &str = "kappa FR3 gap site";
const LAMBDA_START: &str = "lambda N-terminal start";
const HEAVY_START: &str = "heavy N-terminal start";

const RESIDUAL: &[Known] = &[
    known(
        "1MFE_1:H",
        "KCM",
        SameWindow(&["1MEX_1:H", "6EAY_1:H"], "PGLEWIG"),
    ),
    known("1QFW_1:L", "KCM", LoopException(24, 34)),
    known("3UTZ_1:L", "KCM", Open(FR3)),
    known("4LLV_3:L", "KCM", Open(FR3)),
    known("5EOC_2:L", "KCM", Open(FR3)),
    known("5VTA_2:L", "KCM", Open(FR3)),
    known("6BPC_1:L", "KCM", Open(FR3)),
    known("5CEY_1:L", "KCM", Open(LAMBDA_START)),
    known("6NNJ_1:L", "KC", Open(LAMBDA_START)),
    known("1OAY_2:L", "KC", Open(LAMBDA_START)),
    known("3GK8_1:H", "KCM", Open(HEAVY_START)),
];

/// Labels of the stretch of `r` whose reference labels lie in `lo..=hi`.
fn loop_pattern(r: &Record, lo: u16, hi: u16) -> Vec<Label> {
    r.reference
        .iter()
        .copied()
        .filter(|l| (lo..=hi).contains(&l.number))
        .collect()
}

fn check_loop_exception(recs: &[Record], record: &Record, (lo, hi): (u16, u16)) {
    let mine = loop_pattern(record, lo, hi);
    let mut counts: BTreeMap<Vec<Label>, usize> = BTreeMap::new();
    for r in recs
        .iter()
        .filter(|r| r.chain == record.chain && !r.antigen && r.id != record.id)
    {
        *counts.entry(loop_pattern(r, lo, hi)).or_default() += 1;
    }
    counts.retain(|p, _| p.len() == mine.len());
    assert!(!counts.contains_key(&mine), "{} is not alone", record.id);
    let usual = counts.values().max().copied().unwrap_or(0);
    assert!(usual >= 50, "{}: no common pattern at this length", record.id);
}

fn check_evidence(recs: &[Record], record: &Record, evidence: Evidence) {
    match evidence {
        SameWindow(ids, window) => {
            let span = window.len();
            let at = |r: &Record| r.sequence.find(window).expect("window in both");
            let mine = &record.reference[at(record)..][..span];
            for id in ids {
                let other = recs
                    .iter()
                    .find(|r| r.id == *id)
                    .unwrap_or_else(|| panic!("{id} missing from the reference set"));
                assert_ne!(mine, &other.reference[at(other)..][..span], "{id} agrees");
            }
        }
        LoopException(lo, hi) => check_loop_exception(recs, record, (lo, hi)),
        Open(cause) => assert!(!cause.is_empty()),
    }
}

/// Every loop of the reference set takes one label pattern per length and
/// chain type, except the listed domains (sole exceptions, see
/// docs/ANTIBODY.md). The length alone fixes the labels, so a numbering
/// that follows the pattern cannot be wrong about a loop it was not told
/// about.
const LOOPS: [(&str, ChainType, u16, u16); 6] = [
    ("H1", ChainType::Heavy, 26, 35),
    ("H3", ChainType::Heavy, 93, 102),
    ("L1 kappa", ChainType::Kappa, 24, 34),
    ("L1 lambda", ChainType::Lambda, 24, 34),
    ("L3 kappa", ChainType::Kappa, 89, 97),
    ("L3 lambda", ChainType::Lambda, 89, 97),
];

/// Domains whose loop deviates from the pattern of every other domain of the
/// same length: 4XCF_1:H against 121 H3 of that length (its SEQRES and
/// observed sequence equal 4XAW_1:H's, which follows the pattern), 1QFW_1:L
/// L1.
const LOOP_EXCEPTIONS: &[&str] = &["4XCF_1:H", "1QFW_1:L"];

#[test]
#[ignore = "needs fixtures/real/reference; use --release"]
fn reference_loops_follow_one_pattern_per_length() {
    for (dir, scheme, _) in SCHEMES {
        let recs = records(dir, scheme);
        for (name, chain, lo, hi) in LOOPS {
            let mut by_length: BTreeMap<usize, BTreeMap<Vec<Label>, Vec<&Record>>> =
                BTreeMap::new();
            for r in recs.iter().filter(|r| r.chain == chain && !r.antigen) {
                let p = loop_pattern(r, lo, hi);
                by_length.entry(p.len()).or_default().entry(p).or_default().push(r);
            }
            for (len, patterns) in by_length {
                let commonest = patterns.values().map(Vec::len).max().expect("a pattern");
                for (pattern, domains) in &patterns {
                    let strays: Vec<&str> = domains
                        .iter()
                        .map(|r| r.id.as_str())
                        .filter(|id| !LOOP_EXCEPTIONS.contains(id))
                        .collect();
                    let minority = domains.len() < commonest;
                    assert!(
                        !minority || strays.is_empty(),
                        "{dir} {name} length {len}: {strays:?} deviate ({})",
                        show(pattern)
                    );
                }
            }
            let wrong: Vec<&str> = comparison_set(&recs)
                .into_iter()
                .map(|i| &recs[i])
                .filter(|r| r.chain == chain && !LOOP_EXCEPTIONS.contains(&r.id.as_str()))
                .filter(|r| our_loop(r, lo, hi) != loop_pattern(r, lo, hi))
                .map(|r| r.id.as_str())
                .collect();
            assert!(wrong.is_empty(), "{dir} {name}: ours differ in {wrong:?}");
        }
    }
}

fn our_loop(r: &Record, lo: u16, hi: u16) -> Vec<Label> {
    let inside = |l: &&Label| (lo..=hi).contains(&l.number);
    r.ours.iter().filter(inside).copied().collect()
}

fn pct(agree: usize, total: usize) -> String {
    format!("{agree}/{total}")
}

fn tally(recs: &[Record], idx: &[usize], half: &BTreeMap<String, usize>) -> [(usize, usize); 4] {
    let mut out = [(0, 0); 4];
    for &i in idx {
        let r = &recs[i];
        for slot in [0, 1 + half[&r.entry]] {
            out[slot].0 += usize::from(r.agrees());
            out[slot].1 += 1;
        }
    }
    out
}

#[test]
#[ignore = "needs fixtures/real/reference; use --release"]
fn reference_numbering_is_exact_up_to_documented_residuals() {
    let half = halves();
    for (dir, scheme, letter) in SCHEMES {
        let recs = records(dir, scheme);
        if recs.is_empty() {
            println!("{dir}: no reference fixtures, skipped");
            return;
        }
        let twins = inconsistent_twins(&recs);
        let cmp = comparison_set(&recs);
        let in_cmp: BTreeSet<usize> = cmp.iter().copied().collect();
        let listed: BTreeMap<&str, &Known> = RESIDUAL
            .iter()
            .filter(|k| k.schemes.contains(letter))
            .map(|k| (k.id, k))
            .collect();

        let misses: BTreeSet<&str> = cmp
            .iter()
            .map(|&i| &recs[i])
            .filter(|r| !r.agrees())
            .map(|r| r.id.as_str())
            .collect();
        let expected: BTreeSet<&str> = listed.keys().copied().collect();
        assert_eq!(misses, expected, "{dir}: comparison-set differences");

        for r in recs
            .iter()
            .enumerate()
            .filter(|(i, r)| !in_cmp.contains(i) && !r.agrees())
        {
            assert!(
                r.1.antigen || twins.contains_key(&r.1.id),
                "{dir}: {} differs outside the comparison set unexplained",
                r.1.id
            );
        }
        for (id, k) in &listed {
            let record = recs.iter().find(|r| r.id == *id).expect("listed domain");
            check_evidence(&recs, record, k.evidence);
        }

        let all: Vec<usize> = (0..recs.len()).collect();
        let raw = tally(&recs, &all, &half);
        let set = tally(&recs, &cmp, &half);
        println!(
            "{dir}: raw {} (even {}, odd {}); comparison {} (even {}, odd {}); {} antigen, {} twin copies excluded",
            pct(raw[0].0, raw[0].1),
            pct(raw[1].0, raw[1].1),
            pct(raw[2].0, raw[2].1),
            pct(set[0].0, set[0].1),
            pct(set[1].0, set[1].1),
            pct(set[2].0, set[2].1),
            recs.iter().filter(|r| r.antigen).count(),
            twins.len(),
        );
        if std::env::var_os("VIZVIZ_EVIDENCE").is_some() {
            print_evidence(dir, scheme, &recs, &listed, &half);
        }
    }
}

const AA: &[u8; 20] = b"ARNDCQEGHILKMFPSTWYV";
const BACKGROUND: [f64; 20] = [
    0.078, 0.051, 0.045, 0.054, 0.019, 0.043, 0.063, 0.074, 0.022, 0.051, 0.091, 0.057, 0.022,
    0.039, 0.052, 0.071, 0.058, 0.013, 0.032, 0.064,
];
const PSEUDOCOUNT: f64 = 5.0;
/// Largest number of placements printed and scored per domain.
const MAX_PLACEMENTS: usize = 200_000;

type Counts = BTreeMap<Label, [f64; 20]>;

/// Residue counts per label over the even half of the comparison set, for
/// one chain type.
fn even_counts(recs: &[Record], half: &BTreeMap<String, usize>, chain: ChainType) -> Counts {
    let mut counts = Counts::new();
    for &i in &comparison_set(recs) {
        let r = &recs[i];
        if r.chain != chain || half[&r.entry] != 0 {
            continue;
        }
        for (aa, label) in r.sequence.bytes().zip(&r.reference) {
            if let Some(a) = AA.iter().position(|&x| x == aa) {
                counts.entry(*label).or_insert([0.0; 20])[a] += 1.0;
            }
        }
    }
    counts
}

fn log_odds(counts: &Counts, label: Label, aa: u8) -> f64 {
    let Some(a) = AA.iter().position(|&x| x == aa) else {
        return 0.0;
    };
    let Some(row) = counts.get(&label) else {
        return 0.0;
    };
    let total: f64 = row.iter().sum();
    let p = (row[a] + PSEUDOCOUNT * BACKGROUND[a]) / (total + PSEUDOCOUNT);
    (p / BACKGROUND[a]).log2()
}

fn fit(counts: &Counts, labels: &[Label], seq: &[u8]) -> f64 {
    labels
        .iter()
        .zip(seq)
        .map(|(l, &aa)| log_odds(counts, *l, aa))
        .sum()
}

fn show(labels: &[Label]) -> String {
    labels
        .iter()
        .map(Label::to_string)
        .collect::<Vec<_>>()
        .join(" ")
}

/// Calls `visit` with every strictly rising choice of `n` of `slots`.
fn choose(
    slots: &[Label],
    n: usize,
    from: usize,
    picked: &mut Vec<Label>,
    visit: &mut dyn FnMut(&[Label]),
) {
    if picked.len() == n {
        visit(picked);
        return;
    }
    let need = n - picked.len();
    for at in from..=slots.len().saturating_sub(need) {
        picked.push(slots[at]);
        choose(slots, n, at + 1, picked, visit);
        picked.pop();
    }
}

fn binomial(n: usize, k: usize) -> f64 {
    (0..k).fold(1.0, |acc, i| acc * (n - i) as f64 / (i + 1) as f64)
}

fn print_evidence(
    dir: &str,
    scheme: Scheme,
    recs: &[Record],
    listed: &BTreeMap<&str, &Known>,
    half: &BTreeMap<String, usize>,
) {
    let mut tables: BTreeMap<String, Counts> = BTreeMap::new();
    for id in listed.keys() {
        let r = recs.iter().find(|r| r.id == *id).expect("listed domain");
        let counts = tables
            .entry(r.chain.name().to_string())
            .or_insert_with(|| even_counts(recs, half, r.chain));
        let seq: Vec<u8> = r.sequence.bytes().collect();
        let rule = r.domain.rule_labels(scheme);
        let differs = |k: usize| r.ours[k] != r.reference[k] || rule[k] != r.reference[k];
        let first = (0..seq.len()).find(|&k| differs(k)).expect("a difference");
        let last = (0..seq.len())
            .rev()
            .find(|&k| differs(k))
            .expect("a difference");
        let span = first..=last;
        let (res, ours, refl, rl) = (
            &seq[span.clone()],
            &r.ours[span.clone()],
            &r.reference[span.clone()],
            &rule[span.clone()],
        );
        println!(
            "\n[{dir}] {id} {} residues {}..={}: {}",
            r.chain.name(),
            first,
            last,
            String::from_utf8_lossy(res)
        );
        println!("  rule      {:>7.1}  {}", fit(counts, rl, res), show(rl));
        println!(
            "  ours      {:>7.1}  {}",
            fit(counts, ours, res),
            show(ours)
        );
        println!(
            "  reference {:>7.1}  {}",
            fit(counts, refl, res),
            show(refl)
        );
        let lo = refl
            .iter()
            .chain(ours)
            .map(|l| l.number)
            .min()
            .expect("labels");
        let hi = refl
            .iter()
            .chain(ours)
            .map(|l| l.number)
            .max()
            .expect("labels");
        let lettered = refl.iter().chain(ours).any(|l| l.insertion().is_some());
        let slots: Vec<Label> = (lo..=hi).map(Label::new).collect();
        if lettered || binomial(slots.len(), res.len()) > MAX_PLACEMENTS as f64 {
            println!("  placements not enumerated (insertion letters or too many)");
            continue;
        }
        let reference_fit = fit(counts, refl, res);
        let mut all: Vec<(f64, Vec<Label>)> = Vec::new();
        choose(&slots, res.len(), 0, &mut Vec::new(), &mut |p| {
            all.push((fit(counts, p, res), p.to_vec()));
        });
        all.sort_by(|a, b| b.0.total_cmp(&a.0));
        let rank = 1 + all
            .iter()
            .filter(|(f, _)| *f > reference_fit + 1e-9)
            .count();
        let present = all.iter().any(|(_, p)| p.as_slice() == refl);
        println!(
            "  {} placements onto {lo}..{hi}; reference among them: {present}, rank {rank}; best:",
            all.len()
        );
        for (f, p) in all.iter().take(3) {
            println!("    {f:>7.1}  {}", show(p));
        }
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
