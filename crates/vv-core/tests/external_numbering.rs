//! Externally numbered domains: reading ANARCI's text output and using the
//! labels where native ones would go. The samples in `tests/data/` are real
//! ANARCI output (see docs/ANTIBODY.md, "External backend").
//!
//! `cargo test --release -p vv-core --test external_numbering -- --ignored --nocapture`
//! numbers trastuzumab and the reference set with ANARCI when it can be
//! found (`VIZVIZ_ANARCI`, or `ANARCI` on the PATH) and prints the
//! agreement with the native numbering.

use vv_core::antibody::external::{
    cdr_residues_external, parse_anarci, ExternalDomain, NumberedDomain,
};
use vv_core::antibody::{find_domains, CdrDefinition, ChainType, Label, Scheme};

mod antibody_common;

const KABAT: &str = include_str!("data/anarci_kabat.txt");
const SCFV: &str = include_str!("data/anarci_scfv_kabat.txt");
const IMGT: &str = include_str!("data/anarci_imgt.txt");

const HEAVY: &str = "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS";
const KAPPA: &str = "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK";

fn labels(domain: &NumberedDomain) -> Vec<String> {
    domain.labels.iter().map(ToString::to_string).collect()
}

fn external(d: &NumberedDomain, scheme: Scheme) -> ExternalDomain {
    ExternalDomain::new(d.chain, d.range.clone())
        .with_labels(scheme, d.labels.clone())
        .unwrap()
}

#[test]
fn the_kabat_sample_reads_as_two_domains_with_their_insertions() {
    let domains = parse_anarci(KABAT).unwrap();
    assert_eq!(domains.len(), 2, "c2 and c3 were not antibody chains");
    let (h, k) = (&domains[0], &domains[1]);
    assert_eq!(
        (h.name.as_str(), h.chain, h.range.clone()),
        ("c0", ChainType::Heavy, 0..120)
    );
    assert_eq!(
        (k.name.as_str(), k.chain, k.range.clone()),
        ("c1", ChainType::Kappa, 0..107)
    );
    let h = labels(h);
    assert_eq!(&h[51..55], ["52", "52A", "53", "54"]);
    assert!(h.contains(&"82C".to_string()) && h.contains(&"100C".to_string()));
    assert_eq!(h.len(), 120);
}

#[test]
fn a_chain_with_two_domains_keeps_its_name_and_each_range() {
    let domains = parse_anarci(SCFV).unwrap();
    assert_eq!(domains.len(), 2);
    assert!(domains.iter().all(|d| d.name == "s0"));
    assert_eq!(domains[0].chain, ChainType::Heavy);
    assert_eq!(domains[1].chain, ChainType::Kappa);
    assert_eq!(domains[1].range, 135..242);
}

#[test]
fn a_receptor_domain_reads_with_imgt_labels_under_any_scheme() {
    let domains = parse_anarci(IMGT).unwrap();
    assert_eq!(domains.len(), 3, "c2 was not a receptor or antibody chain");
    let d = &domains[2];
    assert_eq!((d.name.as_str(), d.chain), ("c3", ChainType::TcrBeta));
    let e = external(d, Scheme::Imgt);
    let imgt = e.annotate(Scheme::Imgt, CdrDefinition::Imgt).unwrap();
    let kabat = e.annotate(Scheme::Kabat, CdrDefinition::Kabat).unwrap();
    assert_eq!(imgt, kabat);
    assert!(imgt.iter().any(|a| a.region.cdr() == Some(3)));
}

#[test]
fn trastuzumab_cdrs_under_kabat_are_the_published_ones() {
    let domains = parse_anarci(KABAT).unwrap();
    let h = external(&domains[0], Scheme::Kabat);
    let notes = h.annotate(Scheme::Kabat, CdrDefinition::Kabat).unwrap();
    let cdr = |n: u8| -> String {
        notes
            .iter()
            .filter(|a| a.region.cdr() == Some(n))
            .map(|a| HEAVY.as_bytes()[a.index] as char)
            .collect()
    };
    assert_eq!(cdr(1), "DTYIH");
    assert_eq!(cdr(2), "RIYPTNGYTRYADSVKG");
    assert_eq!(cdr(3), "WGGDGFYAMDY");
}

#[test]
fn labels_for_a_missing_scheme_are_refused() {
    let domains = parse_anarci(KABAT).unwrap();
    let h = external(&domains[0], Scheme::Kabat);
    assert!(h.covers(Scheme::Kabat, CdrDefinition::Kabat));
    assert!(!h.covers(Scheme::Imgt, CdrDefinition::Kabat));
    assert!(h.annotate(Scheme::Kabat, CdrDefinition::Chothia).is_none());
    assert!(ExternalDomain::new(ChainType::Heavy, 0..3)
        .with_labels(Scheme::Kabat, vec![Label::new(1)])
        .is_none());
}

#[test]
fn cdr_residues_offset_by_the_chain_start() {
    let domains = parse_anarci(KABAT).unwrap();
    let h = [external(&domains[0], Scheme::Kabat)];
    let found = cdr_residues_external([(1000u32, &h[..])], CdrDefinition::Kabat);
    let h3: Vec<u32> = found
        .iter()
        .filter(|r| r.cdr == 3)
        .map(|r| r.residue)
        .collect();
    assert_eq!(h3.len(), 11);
    assert_eq!(h3[0], 1000 + HEAVY.find("WGGDGFYAMDY").unwrap() as u32);
    assert!(found.iter().all(|r| r.chain == ChainType::Heavy));
}

#[test]
fn malformed_output_is_an_error_not_a_guess() {
    let bad = KABAT.replacen("H 1       E", "H 1", 1);
    assert!(parse_anarci(&bad).is_err());
    let short = KABAT.replacen("|0|119|", "|0|200|", 1);
    assert!(parse_anarci(&short)
        .unwrap_err()
        .contains("numbered residues"));
    assert!(parse_anarci("").unwrap().is_empty());
}

#[test]
fn anarci_and_native_agree_on_the_sample_sequences() {
    let domains = parse_anarci(KABAT).unwrap();
    for (d, seq) in domains.iter().zip([HEAVY, KAPPA]) {
        let native = find_domains(seq);
        let n = native.iter().find(|n| n.chain == d.chain).unwrap();
        assert_eq!(n.start..n.end, d.range);
        let ours: Vec<Label> = n
            .numbering(Scheme::Kabat)
            .into_iter()
            .map(|(_, l)| l)
            .collect();
        let agree = ours.iter().zip(&d.labels).filter(|(a, b)| a == b).count();
        println!("{}: {agree}/{} residues agree", d.name, d.labels.len());
        assert!(agree * 100 >= d.labels.len() * 95, "{}: {agree}", d.name);
    }
}

// ---------------------------------------------------------------------------
// Live comparison with the real ANARCI (ignored: needs the program).

use std::collections::{BTreeMap, HashSet};
use std::io::Write;
use std::process::{Command, Stdio};

use antibody_common::{fixtures, read_chains};

const SEED_IDS: &str = include_str!("../src/antibody/seed_ids.txt");

/// Same stdin-to-tempfile wrapper the application uses inside WSL.
const WSL_SCRIPT: &str = "exe=$1; shift; case $exe in /*) PATH=${exe%/*}:$PATH;; esac; f=$(mktemp) || exit 1; cat >$f; $exe -i $f $@; st=$?; unlink $f; exit $st";

/// ANARCI's text output for `fasta`, from `$VIZVIZ_ANARCI` (`wsl:/path`
/// runs inside WSL) or `ANARCI` on the PATH.
fn run_anarci(fasta: &str, scheme: &str) -> Result<String, String> {
    let setting = std::env::var("VIZVIZ_ANARCI").unwrap_or_else(|_| "ANARCI".into());
    let mut args = vec!["-s", scheme];
    if scheme != "imgt" {
        args.extend(["-r", "ig"]);
    }
    let output = match setting.strip_prefix("wsl:") {
        Some(exe) => {
            let mut child = Command::new("wsl.exe")
                .args(["-e", "sh", "-lc", WSL_SCRIPT, "sh", exe])
                .args(&args)
                .stdin(Stdio::piped())
                .stdout(Stdio::piped())
                .stderr(Stdio::null())
                .spawn()
                .map_err(|e| e.to_string())?;
            child
                .stdin
                .take()
                .unwrap()
                .write_all(fasta.as_bytes())
                .unwrap();
            child.wait_with_output()
        }
        None => {
            let file = std::env::temp_dir().join(format!("vizviz-ref-{}.fa", std::process::id()));
            std::fs::write(&file, fasta).unwrap();
            let out = Command::new(&setting)
                .arg("-i")
                .arg(&file)
                .args(&args)
                .output();
            std::fs::remove_file(&file).ok();
            out
        }
    };
    let output = output.map_err(|e| format!("cannot run `{setting}`: {e}"))?;
    match output.status.success() {
        true => Ok(String::from_utf8_lossy(&output.stdout).into_owned()),
        false => Err(format!("`{setting}` exited with {}", output.status)),
    }
}

#[derive(Default, Debug)]
struct Live {
    both: usize,
    native_only: usize,
    anarci_only: usize,
    native_exact: usize,
    anarci_exact: usize,
    same_exact: usize,
    residues: usize,
    native_ok: usize,
    anarci_ok: usize,
    same_ok: usize,
}

fn pct(a: usize, b: usize) -> f64 {
    100.0 * a as f64 / b.max(1) as f64
}

fn compare_live(scheme: Scheme, dir: &str) -> Result<Live, String> {
    let seeds: HashSet<&str> = SEED_IDS.split_whitespace().collect();
    let mut chains = BTreeMap::new();
    for (name, path) in fixtures(&format!("reference/{dir}")) {
        if seeds.contains(name[..4].to_ascii_uppercase().as_str()) {
            continue;
        }
        for c in read_chains(&path) {
            chains.insert(format!("{name}:{}", c.id), c);
        }
    }
    let fasta: String = chains
        .iter()
        .map(|(k, c)| format!(">{k}\n{}\n", c.sequence()))
        .collect();
    let text = run_anarci(&fasta, dir)?;
    let mut theirs: BTreeMap<String, Vec<NumberedDomain>> = BTreeMap::new();
    for d in parse_anarci(&text)? {
        theirs.entry(d.name.clone()).or_default().push(d);
    }
    let mut out = Live::default();
    for (key, chain) in &chains {
        let native = find_domains(&chain.sequence());
        let others = theirs.remove(key).unwrap_or_default();
        for n in &native {
            let Some(a) = others.iter().find(|a| a.chain == n.chain) else {
                out.native_only += 1;
                continue;
            };
            out.both += 1;
            let ours: Vec<Label> = n.numbering(scheme).into_iter().map(|(_, l)| l).collect();
            let author = |i: usize| chain.residues[i].label;
            let native_exact = ours
                .iter()
                .enumerate()
                .all(|(k, l)| *l == author(n.start + k));
            let anarci_exact = a
                .labels
                .iter()
                .enumerate()
                .all(|(k, l)| *l == author(a.range.start + k));
            out.native_exact += usize::from(native_exact);
            out.anarci_exact += usize::from(anarci_exact);
            out.same_exact += usize::from((n.start..n.end) == a.range && ours == a.labels);
            for i in n.start.max(a.range.start)..n.end.min(a.range.end) {
                let (x, y) = (ours[i - n.start], a.labels[i - a.range.start]);
                out.residues += 1;
                out.native_ok += usize::from(x == author(i));
                out.anarci_ok += usize::from(y == author(i));
                out.same_ok += usize::from(x == y);
            }
        }
        out.anarci_only += others
            .iter()
            .filter(|a| !native.iter().any(|n| n.chain == a.chain))
            .count();
    }
    Ok(out)
}

#[test]
#[ignore = "needs ANARCI and fixtures/real/reference; prints native vs ANARCI agreement"]
fn reference_set_native_vs_anarci() {
    for (dir, scheme) in [
        ("kabat", Scheme::Kabat),
        ("chothia", Scheme::Chothia),
        ("martin", Scheme::Martin),
    ] {
        let live = match compare_live(scheme, dir) {
            Ok(l) => l,
            Err(why) => {
                println!("{dir}: skipped ({why})");
                return;
            }
        };
        println!(
            "{dir}: {} domains found by both ({} native only, {} ANARCI only)\n  \
             whole domain equal to the reference: native {} ({:.1}%), ANARCI {} ({:.1}%); \
             native == ANARCI {} ({:.1}%)\n  \
             residues in both ranges {}: native {:.2}%, ANARCI {:.2}%, native == ANARCI {:.2}%",
            live.both,
            live.native_only,
            live.anarci_only,
            live.native_exact,
            pct(live.native_exact, live.both),
            live.anarci_exact,
            pct(live.anarci_exact, live.both),
            live.same_exact,
            pct(live.same_exact, live.both),
            live.residues,
            pct(live.native_ok, live.residues),
            pct(live.anarci_ok, live.residues),
            pct(live.same_ok, live.residues),
        );
    }
}
