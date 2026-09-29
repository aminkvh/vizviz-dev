//! Antibody numbering against wwPDB structures whose depositors numbered the
//! chains themselves. The entries live in the git-ignored `fixtures/real/`
//! (see docs/ANTIBODY.md for how they were chosen); every test skips when
//! they are absent, so CI stays green without a download.
//!
//! `cargo test --release -p vv-core --test antibody_pdb -- --ignored --nocapture`
//! prints the accuracy report and the timing run.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use vv_core::antibody::{find_domains, ChainType, Domain, Label, Scheme};

struct Residue {
    label: Label,
    aa: char,
}

struct Chain {
    id: char,
    residues: Vec<Residue>,
}

impl Chain {
    fn sequence(&self) -> String {
        self.residues.iter().map(|r| r.aa).collect()
    }

    fn author(&self, d: &Domain) -> Vec<Label> {
        self.residues[d.start..d.end]
            .iter()
            .map(|r| r.label)
            .collect()
    }
}

fn one_letter(name: &str) -> char {
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
fn read_chains(path: &Path) -> Vec<Chain> {
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

fn fixtures(dir: &str) -> Vec<(String, PathBuf)> {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../fixtures/real")
        .join(dir);
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

fn matches(d: &Domain, scheme: Scheme, author: &[Label]) -> usize {
    d.numbering(scheme)
        .iter()
        .zip(author)
        .filter(|((_, l), a)| l == *a)
        .count()
}

/// The depositor put the four conserved anchors (Cys, Trp, Cys, FR4 Trp/Phe)
/// at their Kabat/Chothia numbers, so the rest of the numbering is comparable.
fn author_used_kabat_frame(d: &Domain, author: &[Label]) -> bool {
    let anchors: [u16; 4] = if d.chain == ChainType::Heavy {
        [22, 36, 92, 103]
    } else {
        [23, 35, 88, 98]
    };
    let ours = d.numbering(Scheme::Kabat);
    anchors.iter().all(|&a| {
        ours.iter()
            .zip(author)
            .any(|((_, l), au)| l.number == a && l.insertion().is_none() && l == au)
    })
}

#[derive(Default)]
struct Tally {
    chains: usize,
    domains: usize,
    anchored: usize,
    exact: [usize; 3],
    exact_any: usize,
    residues: usize,
    residues_agreeing: usize,
    by_chain: BTreeMap<String, (usize, usize)>,
}

const SCHEMES: [Scheme; 3] = [Scheme::Kabat, Scheme::Chothia, Scheme::Martin];

fn evaluate(dir: &str) -> Tally {
    let mut t = Tally::default();
    for (_, path) in fixtures(dir) {
        for chain in read_chains(&path) {
            t.chains += 1;
            for d in find_domains(&chain.sequence()) {
                t.domains += 1;
                let author = chain.author(&d);
                if !author_used_kabat_frame(&d, &author) {
                    continue;
                }
                t.anchored += 1;
                let scores = SCHEMES.map(|s| matches(&d, s, &author));
                let best = scores.into_iter().max().unwrap_or(0);
                for (slot, m) in t.exact.iter_mut().zip(scores) {
                    *slot += usize::from(m == author.len());
                }
                t.residues += author.len();
                t.residues_agreeing += best;
                let row = t.by_chain.entry(format!("{:?}", d.chain)).or_default();
                row.0 += 1;
                row.1 += usize::from(best == author.len());
                t.exact_any += usize::from(best == author.len());
            }
        }
    }
    t
}

fn print_tally(name: &str, t: &Tally) {
    println!(
        "{name}: {} chains, {} domains found, {} with Kabat/Chothia-framed author numbering",
        t.chains, t.domains, t.anchored
    );
    println!(
        "  exact: Kabat {} Chothia {} Martin {} any {} of {}",
        t.exact[0], t.exact[1], t.exact[2], t.exact_any, t.anchored
    );
    println!(
        "  residue agreement (best scheme per domain): {}/{} = {:.2}%",
        t.residues_agreeing,
        t.residues,
        100.0 * t.residues_agreeing as f64 / t.residues.max(1) as f64
    );
    for (chain, (n, exact)) in &t.by_chain {
        println!("  {chain}: {exact}/{n} exact");
    }
}

#[test]
fn author_kabat_numbering_is_reproduced_for_reference_fabs() {
    // Chains whose depositors used Kabat numbering throughout.
    let exact = [
        ("1IGY", "ABCD"),
        ("1HZH", "HK"),
        ("2AAB", "H"),
        ("1N8Z", "A"),
        ("1FVC", "AC"),
        ("1MLC", "AC"),
        ("3HFM", "L"),
        ("1BJ1", "JL"),
    ];
    for (id, chains) in exact {
        let Some((_, path)) = fixtures("").into_iter().find(|(n, _)| n == id) else {
            return;
        };
        for chain in read_chains(&path).iter().filter(|c| chains.contains(c.id)) {
            let found = find_domains(&chain.sequence());
            assert_eq!(found.len(), 1, "{id}:{}", chain.id);
            let author = chain.author(&found[0]);
            assert_eq!(
                matches(&found[0], Scheme::Kabat, &author),
                author.len(),
                "{id}:{}",
                chain.id
            );
        }
    }
}

#[test]
fn trastuzumab_heavy_agrees_with_the_depositor_up_to_cdr_h2() {
    // 1N8Z numbers the heavy chain in Kabat through H52 and sequentially after.
    let Some((_, path)) = fixtures("").into_iter().find(|(n, _)| n == "1N8Z") else {
        return;
    };
    let chain = read_chains(&path)
        .into_iter()
        .find(|c| c.id == 'B')
        .unwrap();
    let d = &find_domains(&chain.sequence())[0];
    let author = chain.author(d);
    let ours = d.numbering(Scheme::Kabat);
    assert!(ours.iter().zip(&author).take(52).all(|((_, l), a)| l == a));
    assert_eq!(ours[52].1.to_string(), "52A");
}

#[test]
fn real_entries_meet_the_accuracy_floor() {
    for dir in ["", "holdout"] {
        let t = evaluate(dir);
        if t.anchored < 50 {
            continue;
        }
        assert!(
            t.residues_agreeing * 100 >= t.residues * 98,
            "{dir}: residue agreement below 98%"
        );
        assert!(
            t.exact_any * 100 >= t.anchored * 70,
            "{dir}: exact domains below 70%"
        );
    }
}

#[test]
fn antigens_in_the_reference_entries_yield_no_domain() {
    let antigens = [
        ("1MLC", "EF"),
        ("3HFM", "Y"),
        ("1MEL", "LM"),
        ("1BJ1", "VW"),
        ("1N8Z", "C"),
        ("4KRL", "A"),
    ];
    for (id, chains) in antigens {
        let Some((_, path)) = fixtures("").into_iter().find(|(n, _)| n == id) else {
            return;
        };
        for chain in read_chains(&path).iter().filter(|c| chains.contains(c.id)) {
            assert!(
                find_domains(&chain.sequence()).is_empty(),
                "{id}:{}",
                chain.id
            );
        }
    }
}

#[test]
#[ignore = "prints the accuracy report; needs fixtures/real"]
fn accuracy_report() {
    print_tally("entries that also seeded the profiles", &evaluate(""));
    print_tally(
        "held-out entries (not in seed_ids.txt)",
        &evaluate("holdout"),
    );
}

#[test]
#[ignore = "timing run; use --release"]
fn timing_10k_sequences() {
    let mut sequences: Vec<String> = Vec::new();
    for dir in ["", "holdout"] {
        for (_, path) in fixtures(dir) {
            sequences.extend(read_chains(&path).iter().map(Chain::sequence));
        }
    }
    if sequences.is_empty() {
        return;
    }
    let batch: Vec<&String> = sequences.iter().cycle().take(10_000).collect();
    let residues: usize = batch.iter().map(|s| s.len()).sum();
    let start = std::time::Instant::now();
    let domains: usize = batch.iter().map(|s| find_domains(s).len()).sum();
    let elapsed = start.elapsed();
    println!(
        "10000 sequences ({residues} residues, {domains} domains): {:.2} s, {:.0} us per sequence, single thread",
        elapsed.as_secs_f64(),
        elapsed.as_micros() as f64 / 10_000.0
    );
}
