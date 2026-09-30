//! T-cell receptor variable domains against wwPDB entries that did not seed
//! the receptor profiles: `fixtures/real/tcr/holdout/` (receptor complexes,
//! `tcr_seed_ids.txt` lists what was left out of it) and the antibody entries
//! in `fixtures/real/holdout/`. Every test skips when a fixture is absent.
//!
//! `cargo test --release -p vv-core --test antibody_tcr -- --ignored --nocapture`
//! prints the confusion table and the IMGT agreement.

use std::collections::{BTreeMap, BTreeSet};

use vv_core::antibody::{
    find_variable_domains, find_variable_domains_with, ChainType, Domain, Label, Scheme,
};

mod antibody_common;
use antibody_common::{fixtures, molecule_names, read_chains, Chain};

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug)]
enum Truth {
    Alpha,
    Beta,
    /// Gamma or delta receptor chains: no profile yet, so not scored.
    GammaDelta,
    Antibody,
    Other,
}

/// Constant-domain motifs of receptor chains, which name the chain type even
/// where a deposited name does not (human and mouse C-alpha, C-beta).
fn constant_domain(seq: &str) -> Option<Truth> {
    if seq.contains("PAVYQL") {
        Some(Truth::Alpha)
    } else if seq.contains("VFPPEVAV") {
        Some(Truth::Beta)
    } else {
        None
    }
}

fn truth_of(name: &str, seq: &str) -> Truth {
    let n = name.to_ascii_lowercase();
    let word = |w: &str| n.split(|c: char| !c.is_alphanumeric()).any(|x| x == w);
    let receptor = [
        "t cell receptor",
        "t-cell receptor",
        "tcr",
        "trav",
        "trbv",
        "tra protein",
    ]
    .iter()
    .any(|w| n.contains(w));
    if let Some(kind) = constant_domain(seq) {
        return kind;
    }
    if receptor && (n.contains("gamma") || n.contains("delta")) {
        Truth::GammaDelta
    } else if receptor && (n.contains("alpha") || word("a")) {
        Truth::Alpha
    } else if receptor && (n.contains("beta") || word("b")) {
        Truth::Beta
    } else if [
        "fab",
        "antibody",
        "immunoglobulin",
        "igg",
        "heavy",
        "light",
        "nanobody",
        "scfv",
    ]
    .iter()
    .any(|w| n.contains(w))
    {
        Truth::Antibody
    } else {
        Truth::Other
    }
}

fn label_of(kind: Option<ChainType>) -> &'static str {
    match kind {
        Some(ChainType::TcrAlpha) => "TCR alpha",
        Some(ChainType::TcrBeta) => "TCR beta",
        Some(ChainType::Heavy | ChainType::Kappa | ChainType::Lambda) => "antibody",
        None => "none",
    }
}

type Confusion = BTreeMap<(Truth, &'static str), usize>;

fn tally_chain(confusion: &mut Confusion, truth: Truth, chain: &Chain) -> Vec<Domain> {
    let domains = find_variable_domains(&chain.sequence());
    if domains.is_empty() {
        *confusion.entry((truth, label_of(None))).or_default() += 1;
    }
    for d in &domains {
        *confusion
            .entry((truth, label_of(Some(d.chain))))
            .or_default() += 1;
    }
    domains
}

/// The depositor put the receptor anchors (Cys, Trp, Cys, FR4 Phe and Gly)
/// at IMGT 23, 41, 104, 118 and 119, so the rest of the numbering is
/// comparable.
fn author_used_imgt(d: &Domain, author: &[Label]) -> bool {
    let ours = d.numbering(Scheme::Imgt);
    [23u16, 41, 104, 118, 119].iter().all(|&a| {
        ours.iter()
            .zip(author)
            .any(|((_, l), au)| l.number == a && l.insertion().is_none() && l == au)
    })
}

/// IMGT numbers of the three CDR loops.
fn in_loop(l: &Label) -> bool {
    matches!(l.number, 27..=38 | 56..=65 | 105..=117)
}

#[derive(Default)]
struct Agreement {
    domains: usize,
    exact: usize,
    residues: usize,
    agreeing: usize,
    framework: (usize, usize),
    loops: (usize, usize),
    /// Domains and how many have every framework residue right, by whether
    /// the depositor numbered the extra FR3 residues 84A to 84C.
    framework_exact: BTreeMap<bool, (usize, usize)>,
    /// Author framework numbers our numbering disagrees with, by number.
    framework_misses: BTreeMap<u16, usize>,
    framework_wrong: Vec<String>,
    misses: Vec<String>,
}

fn record(a: &mut Agreement, id: &str, chain: &Chain, d: &Domain) {
    let author = chain.author(d);
    let ours = d.numbering(Scheme::Imgt);
    let agree = ours
        .iter()
        .zip(&author)
        .filter(|((_, l), au)| l == *au)
        .count();
    for ((_, l), au) in ours.iter().zip(&author) {
        let part = if in_loop(l) {
            &mut a.loops
        } else {
            &mut a.framework
        };
        part.0 += usize::from(l == au);
        part.1 += 1;
        if l != au && !in_loop(l) {
            *a.framework_misses.entry(au.number).or_default() += 1;
        }
    }
    let framework_right = ours
        .iter()
        .zip(&author)
        .filter(|((_, l), _)| !in_loop(l))
        .all(|((_, l), au)| l == au);
    let extras = author.contains(&Label::with_insertion(84, 'A'));
    let entry = a.framework_exact.entry(extras).or_default();
    entry.0 += 1;
    entry.1 += usize::from(framework_right);
    if !framework_right && !extras && a.framework_wrong.len() < 20 {
        a.framework_wrong.push(format!("{id}:{}", chain.id));
    }
    a.domains += 1;
    a.residues += author.len();
    a.agreeing += agree;
    if agree == author.len() {
        a.exact += 1;
    } else if a.misses.len() < 10 {
        a.misses.push(format!("{id}:{}", chain.id));
    }
}

struct Report {
    /// Entries with at least one receptor chain, and those with an
    /// IMGT-framed receptor domain.
    entries: (BTreeSet<String>, BTreeSet<String>),
    confusion: Confusion,
    receptor_chains: BTreeMap<Truth, usize>,
    agreement: BTreeMap<Truth, Agreement>,
}

fn receptor_entries() -> Report {
    let mut r = Report {
        entries: Default::default(),
        confusion: Confusion::new(),
        receptor_chains: BTreeMap::new(),
        agreement: BTreeMap::new(),
    };
    for (id, path) in fixtures("tcr/holdout") {
        let names = molecule_names(&path);
        for chain in read_chains(&path) {
            let truth = names
                .get(&chain.id)
                .map_or(Truth::Other, |n| truth_of(n, &chain.sequence()));
            let domains = tally_chain(&mut r.confusion, truth, &chain);
            if matches!(truth, Truth::Alpha | Truth::Beta) {
                r.entries.0.insert(id.clone());
                *r.receptor_chains.entry(truth).or_default() += 1;
                let expected = if truth == Truth::Alpha {
                    ChainType::TcrAlpha
                } else {
                    ChainType::TcrBeta
                };
                for d in domains.iter().filter(|d| d.chain == expected) {
                    if author_used_imgt(d, &chain.author(d)) {
                        r.entries.1.insert(id.clone());
                        record(r.agreement.entry(truth).or_default(), &id, &chain, d);
                    }
                }
            }
        }
    }
    r
}

fn antibody_entries() -> Confusion {
    let mut confusion = Confusion::new();
    for (_, path) in fixtures("holdout") {
        let names = molecule_names(&path);
        for chain in read_chains(&path) {
            let name = names.get(&chain.id).map_or("", String::as_str);
            let truth = match truth_of(name, &chain.sequence()) {
                receptor @ (Truth::Alpha | Truth::Beta | Truth::GammaDelta) => receptor,
                _ => Truth::Antibody,
            };
            tally_chain(&mut confusion, truth, &chain);
        }
    }
    confusion
}

fn count(c: &Confusion, truth: Truth, found: &str) -> usize {
    c.get(&(truth, found)).copied().unwrap_or(0)
}

#[test]
fn receptors_are_detected_as_receptors_and_never_as_antibodies() {
    let r = receptor_entries();
    let chains: usize = r.receptor_chains.values().sum();
    if chains == 0 {
        return;
    }
    for truth in [Truth::Alpha, Truth::Beta] {
        assert_eq!(count(&r.confusion, truth, "antibody"), 0, "{truth:?}");
    }
    let (alpha, beta) = (Truth::Alpha, Truth::Beta);
    assert!(count(&r.confusion, alpha, "TCR alpha") * 100 >= r.receptor_chains[&alpha] * 90);
    assert!(count(&r.confusion, beta, "TCR beta") * 100 >= r.receptor_chains[&beta] * 90);
    assert_eq!(count(&r.confusion, alpha, "TCR beta"), 0);
    assert_eq!(count(&r.confusion, beta, "TCR alpha"), 0);
}

#[test]
fn antibodies_are_never_reported_as_receptors() {
    let c = antibody_entries();
    assert_eq!(count(&c, Truth::Antibody, "TCR alpha"), 0);
    assert_eq!(count(&c, Truth::Antibody, "TCR beta"), 0);
}

#[test]
fn receptor_numbering_matches_imgt_framed_depositions() {
    let r = receptor_entries();
    for (truth, a) in &r.agreement {
        if a.domains < 20 {
            continue;
        }
        assert!(a.agreeing * 100 >= a.residues * 97, "{truth:?}");
    }
}

fn print_confusion(name: &str, c: &Confusion) {
    println!("{name}:");
    for truth in [
        Truth::Alpha,
        Truth::Beta,
        Truth::GammaDelta,
        Truth::Antibody,
        Truth::Other,
    ] {
        let row: Vec<String> = ["TCR alpha", "TCR beta", "antibody", "none"]
            .iter()
            .map(|f| format!("{f} {}", count(c, truth, f)))
            .collect();
        println!("  {truth:?}: {}", row.join(", "));
    }
}

#[test]
#[ignore = "prints the confusion table and IMGT agreement; needs fixtures/real"]
fn receptor_report() {
    let r = receptor_entries();
    println!(
        "{} held-out receptor entries, {} with an IMGT-framed receptor domain",
        r.entries.0.len(),
        r.entries.1.len()
    );
    print_confusion(
        "receptor entries (chains named by the depositor)",
        &r.confusion,
    );
    print_confusion("antibody entries", &antibody_entries());
    for (truth, chains) in &r.receptor_chains {
        println!("{truth:?}: {chains} chains named as receptor chains");
    }
    for (truth, a) in &r.agreement {
        println!(
            "{truth:?} with IMGT-framed author numbering: {} domains, whole-domain identical {}, residues {}/{} ({:.2}%), framework {}/{} ({:.2}%), loops {}/{} ({:.2}%); misses {}",
            a.domains,
            a.exact,
            a.agreeing,
            a.residues,
            100.0 * a.agreeing as f64 / a.residues.max(1) as f64,
            a.framework.0,
            a.framework.1,
            100.0 * a.framework.0 as f64 / a.framework.1.max(1) as f64,
            a.loops.0,
            a.loops.1,
            100.0 * a.loops.0 as f64 / a.loops.1.max(1) as f64,
            a.misses.join(" ")
        );
        println!(
            "  every framework residue right, by author 84A-84C present: {:?}",
            a.framework_exact
        );
        println!("  framework numbers we miss: {:?}", a.framework_misses);
        println!(
            "  framework wrong without 84A: {}",
            a.framework_wrong.join(" ")
        );
    }
}

#[test]
#[ignore = "TCR_ID=3DX9 TCR_CHAIN=A; prints where the IMGT numbering differs"]
fn receptor_numbering_diff() {
    let (Ok(id), Ok(chain_id)) = (std::env::var("TCR_ID"), std::env::var("TCR_CHAIN")) else {
        return;
    };
    let (_, path) = fixtures("tcr/holdout")
        .into_iter()
        .find(|(n, _)| n.eq_ignore_ascii_case(&id))
        .unwrap();
    let chain = read_chains(&path)
        .into_iter()
        .find(|c| c.id.to_string() == chain_id)
        .unwrap();
    let seq = chain.sequence();
    for d in find_variable_domains(&seq) {
        println!("{:?} {}..{}", d.chain, d.start, d.end);
        for ((i, ours), theirs) in d.numbering(Scheme::Imgt).iter().zip(chain.author(&d)) {
            let mark = if *ours == theirs { ' ' } else { '*' };
            println!(
                "{mark}{i} {} ours {ours} theirs {theirs}",
                seq.as_bytes()[*i] as char
            );
        }
    }
}
#[test]
#[ignore = "lists the chains the confusion table counts as errors"]
fn receptor_problems() {
    for (dir, is_receptor_dir) in [("tcr/holdout", true), ("holdout", false)] {
        for (id, path) in fixtures(dir) {
            let names = molecule_names(&path);
            for chain in read_chains(&path) {
                let name = names.get(&chain.id).cloned().unwrap_or_default();
                let truth = truth_of(&name, &chain.sequence());
                let found: Vec<_> = find_variable_domains(&chain.sequence())
                    .iter()
                    .map(|d| (d.chain, format!("{:.2}", d.confidence)))
                    .collect();
                let wanted = match truth {
                    Truth::Alpha => Some(ChainType::TcrAlpha),
                    Truth::Beta => Some(ChainType::TcrBeta),
                    _ => None,
                };
                if truth == Truth::GammaDelta {
                    continue;
                }
                let wrong = match wanted {
                    Some(w) => found.is_empty() || found.iter().any(|(c, _)| *c != w),
                    None => found.iter().any(|(c, _)| !c.is_antibody()),
                };
                if wrong && (is_receptor_dir || found.iter().any(|(c, _)| !c.is_antibody())) {
                    println!(
                        "{dir} {id}:{} {truth:?} len {} name `{name}` found {found:?}",
                        chain.id,
                        chain.residues.len()
                    );
                }
            }
        }
    }
}
/// Best receptor-profile confidence of a chain, 0 when nothing aligns.
fn receptor_confidence(chain: &Chain) -> f32 {
    find_variable_domains_with(&chain.sequence(), 0.0)
        .iter()
        .filter(|d| !d.chain.is_antibody())
        .map(|d| d.confidence)
        .fold(0.0, f32::max)
}

#[test]
#[ignore = "prints receptor-profile confidence by kind of chain; needs fixtures/real"]
fn receptor_confidence_report() {
    let mut by_kind: BTreeMap<&str, Vec<f32>> = BTreeMap::new();
    for (dir, receptor_dir) in [("tcr/holdout", true), ("holdout", false)] {
        for (_, path) in fixtures(dir) {
            let names = molecule_names(&path);
            for chain in read_chains(&path) {
                if chain.residues.len() < 60 {
                    continue;
                }
                let truth = names
                    .get(&chain.id)
                    .map_or(Truth::Other, |n| truth_of(n, &chain.sequence()));
                let kind = match (receptor_dir, truth) {
                    (true, Truth::Alpha | Truth::Beta) => "receptor chain",
                    (true, Truth::Other) => "other chain in a receptor entry",
                    (true, Truth::GammaDelta) => "gamma or delta receptor chain",
                    (true, Truth::Antibody) => "antibody chain in a receptor entry",
                    (false, _) => "chain of an antibody entry",
                };
                by_kind
                    .entry(kind)
                    .or_default()
                    .push(receptor_confidence(&chain));
            }
        }
    }
    for (kind, mut v) in by_kind {
        v.sort_by(f32::total_cmp);
        let at = |p: f64| v[((v.len() - 1) as f64 * p) as usize];
        println!(
            "{kind}: n={} min {:.2} p5 {:.2} median {:.2} p95 {:.2} max {:.2}",
            v.len(),
            v[0],
            at(0.05),
            at(0.5),
            at(0.95),
            v[v.len() - 1]
        );
    }
}
