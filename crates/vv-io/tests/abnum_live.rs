//! Abnum against native numbering and the scheme authors' own renumbered
//! reference set. Both tests send sequences to the public server
//! (http://www.bioinf.org.uk/abs/abnum/, plain HTTP, replies cached in the
//! vizviz cache directory, at most one request per half second).
//!
//! `cargo test --release -p vv-io --test abnum_live -- --ignored --nocapture`
//!
//! The reference test needs `fixtures/real/reference` (set
//! `VIZVIZ_FIXTURES_REAL` from a worktree) and prints the table that
//! docs/ANTIBODY.md ("Scheme authors' program as a backend") quotes.

use std::collections::BTreeSet;

use vv_core::antibody::abnum::{flag, parse_abnum, SCHEMES};
use vv_core::antibody::{find_domains, Label, Scheme};

#[path = "../../vv-core/tests/antibody_common/mod.rs"]
mod antibody_common;
use antibody_common::reference::{comparison_set, records, Record, SCHEMES as REFERENCE};

const HEAVY: &str = "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS";
const KAPPA: &str = "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK";

/// `RESIDUAL` of `vv-core/tests/antibody_reference_exact.rs`: the domains
/// where native numbering differs from the reference.
const RESIDUAL: [&str; 9] = [
    "1OAY_2:L", "3GK8_1:H", "3UTZ_1:L", "4LLV_3:L", "5CEY_1:L", "5EOC_2:L", "5VTA_2:L", "6BPC_1:L",
    "6NNJ_1:L",
];
const RANDOM_DOMAINS: usize = 40;

/// Abnum's label for each residue of `sequence` under `scheme`, `None`
/// outside the domain it numbered or when it found none.
fn abnum_labels(sequence: &str, scheme: Scheme) -> Vec<Option<Label>> {
    let field = flag(scheme).expect("an offered scheme");
    let reply = vv_io::seqdata::abnum::number(sequence, field).expect("Abnum answered");
    let mut out = vec![None; sequence.len()];
    if let Some(d) = parse_abnum(sequence, &reply).expect("reply parses") {
        for (slot, label) in out[d.range].iter_mut().zip(d.labels) {
            *slot = Some(label);
        }
    }
    out
}

#[test]
#[ignore = "sends sequences to the public Abnum server"]
fn trastuzumab_matches_native_in_every_scheme() {
    for sequence in [HEAVY, KAPPA] {
        let native = &find_domains(sequence)[0];
        for scheme in SCHEMES {
            let theirs = abnum_labels(sequence, scheme);
            let ours: Vec<Option<Label>> = (0..sequence.len())
                .map(|i| {
                    let at = native.numbering(scheme);
                    at.iter().find(|(index, _)| *index == i).map(|(_, l)| *l)
                })
                .collect();
            assert_eq!(theirs, ours, "{} {}", &sequence[..8], scheme.name());
        }
    }
}

#[derive(Default)]
struct Tally {
    domains: usize,
    exact: usize,
    native_equal: usize,
    residues: usize,
    residues_agree: usize,
}

impl Tally {
    fn add(&mut self, abnum: &[Option<Label>], r: &Record) {
        let span = r.domain.start..r.domain.end;
        let theirs = &abnum[span];
        let agree = theirs
            .iter()
            .zip(&r.reference)
            .filter(|(a, b)| **a == Some(**b))
            .count();
        self.domains += 1;
        self.exact += usize::from(agree == r.reference.len());
        self.native_equal +=
            usize::from(theirs.iter().copied().eq(r.ours.iter().copied().map(Some)));
        self.residues += r.reference.len();
        self.residues_agree += agree;
    }

    fn row(&self, what: &str, scheme: Scheme) -> String {
        format!(
            "| {what} | {} | {}/{} | {}/{} | {}/{} |",
            scheme.name(),
            self.exact,
            self.domains,
            self.native_equal,
            self.domains,
            self.residues_agree,
            self.residues
        )
    }
}

fn first_difference(abnum: &[Option<Label>], r: &Record) -> String {
    let span = r.domain.start..r.domain.end;
    let at = abnum[span]
        .iter()
        .zip(&r.reference)
        .position(|(a, b)| *a != Some(*b));
    match at {
        None => "agrees".to_string(),
        Some(i) => {
            let theirs = abnum[r.domain.start + i].map_or("none".to_string(), |l| l.to_string());
            format!(
                "residue {i} ({}): Abnum {theirs}, reference {}, native {}",
                r.sequence.as_bytes()[i] as char,
                r.reference[i],
                r.ours[i]
            )
        }
    }
}

/// Domains that agree with the reference in every scheme, by a fixed
/// pseudo-random order so the selection repeats.
fn random_agreeing(per_scheme: &[Vec<Record>]) -> Vec<String> {
    let sets: Vec<BTreeSet<String>> = per_scheme
        .iter()
        .map(|recs| {
            comparison_set(recs)
                .into_iter()
                .map(|i| &recs[i])
                .filter(|r| r.agrees())
                .map(|r| r.id.clone())
                .collect()
        })
        .collect();
    let mut ids: Vec<String> = sets[0]
        .iter()
        .filter(|id| sets.iter().all(|s| s.contains(*id)))
        .filter(|id| !RESIDUAL.contains(&id.as_str()))
        .cloned()
        .collect();
    let mut state = 0x9e37_79b9_7f4a_7c15_u64;
    for i in (1..ids.len()).rev() {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        ids.swap(i, (state >> 33) as usize % (i + 1));
    }
    ids.truncate(RANDOM_DOMAINS);
    ids
}

#[test]
#[ignore = "sends reference-set chains to the public Abnum server"]
fn abnum_against_the_reference_residuals_and_a_random_forty() {
    let per_scheme: Vec<Vec<Record>> = REFERENCE.iter().map(|(d, s, _)| records(d, *s)).collect();
    if per_scheme.iter().any(Vec::is_empty) {
        println!("no reference fixtures, skipped");
        return;
    }
    let forty = random_agreeing(&per_scheme);
    assert_eq!(forty.len(), RANDOM_DOMAINS);
    let mut rows = Vec::new();
    for (k, (_, scheme, _)) in REFERENCE.iter().enumerate() {
        let (mut nine, mut random) = (Tally::default(), Tally::default());
        for r in &per_scheme[k] {
            let residual = RESIDUAL.contains(&r.id.as_str());
            if !residual && !forty.contains(&r.id) {
                continue;
            }
            let abnum = abnum_labels(&r.chain_sequence, *scheme);
            if residual {
                nine.add(&abnum, r);
                println!("{} {} {}", scheme.name(), r.id, first_difference(&abnum, r));
            } else {
                random.add(&abnum, r);
            }
        }
        rows.push(nine.row("nine residual domains", *scheme));
        rows.push(random.row("40 random agreeing domains", *scheme));
    }
    println!("| Set | Scheme | Abnum = reference | Abnum = native | Residues Abnum = reference |");
    println!("|---|---|---|---|---|");
    for row in rows {
        println!("{row}");
    }
}
