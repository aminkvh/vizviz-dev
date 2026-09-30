//! Kabat, Chothia and Martin numbering against structures renumbered by the
//! schemes' authors, kept in the git-ignored `fixtures/real/reference/`
//! (`kabat`, `chothia`, `martin`, one PDB file per Fv). Skips when absent.
//!
//! `cargo test --release -p vv-core --test antibody_reference -- --ignored --nocapture`
//! prints the agreement report.

use std::collections::HashSet;

use vv_core::antibody::{find_domains, ChainType, Domain, Label, Scheme};

mod antibody_common;
use antibody_common::{fixtures, read_chains, Chain};

const SEED_IDS: &str = include_str!("../src/antibody/seed_ids.txt");

fn entries(dir: &str) -> Vec<(String, std::path::PathBuf)> {
    fixtures(&format!("reference/{dir}"))
}

#[derive(Default)]
struct Agreement {
    domains: usize,
    exact: usize,
    residues: usize,
    agreeing: usize,
    by_chain: [(usize, usize); 3],
    differing: Vec<String>,
    first_difference: std::collections::BTreeMap<String, (usize, String)>,
}

fn chain_slot(c: ChainType) -> usize {
    match c {
        ChainType::Heavy => 0,
        ChainType::Kappa => 1,
        _ => 2,
    }
}

fn compare(d: &Domain, scheme: Scheme, chain: &Chain, id: &str, out: &mut Agreement) {
    let author: Vec<Label> = chain.author(d);
    let ours = d.numbering(scheme);
    let agree = ours
        .iter()
        .zip(&author)
        .filter(|((_, l), a)| l == *a)
        .count();
    out.domains += 1;
    out.residues += author.len();
    out.agreeing += agree;
    let slot = &mut out.by_chain[chain_slot(d.chain)];
    slot.0 += 1;
    if agree == author.len() {
        out.exact += 1;
        slot.1 += 1;
    } else {
        if out.differing.len() < 12 {
            out.differing.push(format!("{id}:{}", chain.id));
        }
        let seen = out
            .first_difference
            .entry(first_difference(&ours, &author, d.chain))
            .or_insert((0, format!("{id}:{}", chain.id)));
        seen.0 += 1;
    }
}

/// `chain ours>theirs` at the first residue whose labels disagree.
fn first_difference(ours: &[(usize, Label)], author: &[Label], chain: ChainType) -> String {
    let at = ours
        .iter()
        .zip(author)
        .position(|((_, l), a)| l != a)
        .unwrap_or(0);
    format!("{chain:?} {}>{}", ours[at].1, author[at])
}

fn evaluate(scheme: Scheme, dir: &str, skip_seeds: bool) -> Agreement {
    let seeds: HashSet<&str> = SEED_IDS.split_whitespace().collect();
    let mut out = Agreement::default();
    for (name, path) in entries(dir) {
        let id = &name[..4];
        if skip_seeds && seeds.contains(id.to_ascii_uppercase().as_str()) {
            continue;
        }
        for chain in read_chains(&path) {
            for d in find_domains(&chain.sequence()) {
                compare(&d, scheme, &chain, id, &mut out);
            }
        }
    }
    out
}

fn print(name: &str, a: &Agreement) {
    let mut firsts: Vec<_> = a.first_difference.iter().collect();
    firsts.sort_by_key(|(_, n)| std::cmp::Reverse(n.0));
    let firsts: Vec<String> = firsts
        .iter()
        .take(14)
        .map(|(k, (n, ex))| format!("{k}: {n} ({ex})"))
        .collect();
    println!("{name}: first differences {}", firsts.join("; "));
    println!(
        "{name}: {} domains, whole-domain identical {} ({:.1}%), residues {}/{} ({:.2}%)",
        a.domains,
        a.exact,
        100.0 * a.exact as f64 / a.domains.max(1) as f64,
        a.agreeing,
        a.residues,
        100.0 * a.agreeing as f64 / a.residues.max(1) as f64
    );
    println!(
        "  heavy {}/{}, kappa {}/{}, lambda {}/{}; first differing: {}",
        a.by_chain[0].1,
        a.by_chain[0].0,
        a.by_chain[1].1,
        a.by_chain[1].0,
        a.by_chain[2].1,
        a.by_chain[2].0,
        a.differing.join(" ")
    );
}

const SCHEMES: [(&str, Scheme); 3] = [
    ("kabat", Scheme::Kabat),
    ("chothia", Scheme::Chothia),
    ("martin", Scheme::Martin),
];

#[test]
#[ignore = "prints the agreement report; needs fixtures/real/reference"]
fn reference_numbering_report() {
    for (dir, scheme) in SCHEMES {
        print(
            &format!("{dir}, not in seed_ids.txt"),
            &evaluate(scheme, dir, true),
        );
    }
}

#[test]
#[ignore = "REF_SCHEME=martin REF_ID=1A6V REF_CHAIN=L; prints where the numbering differs"]
fn reference_numbering_diff() {
    let (Ok(dir), Ok(id), Ok(chain_id)) = (
        std::env::var("REF_SCHEME"),
        std::env::var("REF_ID"),
        std::env::var("REF_CHAIN"),
    ) else {
        return;
    };
    let scheme = SCHEMES.iter().find(|(d, _)| *d == dir).unwrap().1;
    let (_, path) = entries(&dir)
        .into_iter()
        .find(|(n, _)| n[..4].eq_ignore_ascii_case(&id))
        .unwrap();
    let chain = read_chains(&path)
        .into_iter()
        .find(|c| c.id.to_string() == chain_id)
        .unwrap();
    let seq = chain.sequence();
    for d in find_domains(&seq) {
        let imgt = d.numbering(Scheme::Imgt);
        for (((i, ours), theirs), (_, imgt)) in
            d.numbering(scheme).iter().zip(chain.author(&d)).zip(&imgt)
        {
            let mark = if *ours == theirs { ' ' } else { '*' };
            println!(
                "{}{} {} ours {} theirs {} imgt {}",
                mark,
                i,
                seq.as_bytes()[*i] as char,
                ours,
                theirs,
                imgt
            );
        }
    }
}
#[test]
#[ignore = "REF_DUMP=path; writes one TSV row per reference domain for offline analysis"]
fn reference_numbering_dump() {
    let Ok(path) = std::env::var("REF_DUMP") else {
        return;
    };
    let mut out = String::new();
    for (dir, scheme) in SCHEMES {
        for (name, file) in entries(dir) {
            for chain in read_chains(&file) {
                let seq = chain.sequence();
                for d in find_domains(&seq) {
                    let join = |l: Vec<Label>| {
                        l.iter().map(Label::to_string).collect::<Vec<_>>().join(",")
                    };
                    let ours: Vec<Label> = d.numbering(scheme).iter().map(|(_, l)| *l).collect();
                    out += &format!(
                        "{dir}\t{name}\t{}\t{:?}\t{}\t{}\t{}\t{}\n",
                        chain.id,
                        d.chain,
                        &seq[d.start..d.end],
                        join(ours),
                        join(chain.author(&d)),
                        join(d.numbering(Scheme::Imgt).iter().map(|(_, l)| *l).collect())
                    );
                }
            }
        }
    }
    std::fs::write(path, out).unwrap();
}

/// Base label range of each loop in Kabat/Chothia/Martin numbering.
fn loop_bases(heavy: bool) -> Vec<(&'static str, u16, u16)> {
    match heavy {
        true => vec![
            ("H1", 26, 35),
            ("H2", 50, 65),
            ("H3", 93, 102),
            ("HFR1", 1, 25),
            ("HFR2", 36, 49),
            ("HFR4", 103, 113),
        ],
        false => vec![
            ("L1", 24, 34),
            ("L2", 50, 56),
            ("L3", 89, 97),
            ("LFR1", 1, 23),
            ("LFR2", 35, 49),
            ("LFR3", 57, 88),
            ("LFR4", 98, 108),
        ],
    }
}

type DeletedSets = std::collections::BTreeMap<(String, usize, Vec<u16>), usize>;

fn tally_deletions(d: &Domain, scheme: Scheme, chain: &Chain, out: &mut DeletedSets) {
    let author = chain.author(d);
    let ours = d.numbering(scheme);
    for (name, lo, hi) in loop_bases(d.chain == ChainType::Heavy) {
        let theirs: Vec<Label> = ours
            .iter()
            .zip(&author)
            .filter(|((_, l), _)| (lo..=hi).contains(&l.number))
            .map(|(_, a)| *a)
            .collect();
        let plain: Vec<u16> = theirs
            .iter()
            .filter(|l| l.insertion().is_none())
            .map(|l| l.number)
            .collect();
        let deleted: Vec<u16> = (lo..=hi).filter(|n| !plain.contains(n)).collect();
        if theirs.len() < usize::from(hi - lo) + 1 && plain.len() == theirs.len() {
            *out.entry((name.to_string(), theirs.len(), deleted))
                .or_default() += 1;
        }
    }
}

#[test]
#[ignore = "prints which labels the reference numbering leaves out of short loops"]
fn reference_deletion_report() {
    for (dir, scheme) in SCHEMES {
        let mut sets = DeletedSets::new();
        for (_, path) in entries(dir) {
            for chain in read_chains(&path) {
                for d in find_domains(&chain.sequence()) {
                    tally_deletions(&d, scheme, &chain, &mut sets);
                }
            }
        }
        println!("== {dir}");
        let mut rows: Vec<_> = sets.into_iter().collect();
        rows.sort_by_key(|((name, len, _), n)| (name.clone(), *len, std::cmp::Reverse(*n)));
        for ((name, len, deleted), n) in rows
            .into_iter()
            .filter(|(k, n)| *n >= 3 || (k.0.contains('2') && *n >= 1))
        {
            println!("{name} length {len}: deleted {deleted:?} x{n}");
        }
    }
}
