//! The optional ANARCI numbering backend (Dunbar & Deane 2016,
//! Bioinformatics 32:298): every protein chain of a structure goes through
//! the executable once per scheme, and the labels come back as
//! [`ExternalDomain`]s the antibody track reads in place of native ones.
//! Nothing is bundled; `launch.rs` finds the user's install.

mod launch;

use std::fmt;
use std::ops::Range;
use std::path::Path;

pub use launch::ENV_VAR;
use vv_core::antibody::external::{parse_anarci, ExternalDomain, NumberedDomain};
use vv_core::antibody::Scheme;
use vv_core::residue_class::ResidueClass;
use vv_core::seqfeat::one_letter;
use vv_core::Topology;

use launch::Launcher;

/// Where the antibody track gets its numbers.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Backend {
    #[default]
    Native,
    Anarci,
}

impl Backend {
    pub const ALL: [Backend; 2] = [Backend::Native, Backend::Anarci];

    pub fn name(self) -> &'static str {
        match self {
            Backend::Native => "native",
            Backend::Anarci => "ANARCI",
        }
    }

    /// Case-insensitive [`Self::name`].
    pub fn parse(word: &str) -> Option<Backend> {
        Self::ALL
            .into_iter()
            .find(|b| b.name().eq_ignore_ascii_case(word))
    }
}

#[derive(Debug)]
pub enum Error {
    NotFound,
    Failed(String),
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::NotFound => write!(f, "ANARCI not found: install it or set {ENV_VAR}"),
            Error::Failed(why) => write!(f, "ANARCI failed: {why}"),
        }
    }
}

/// Domains of each chain row, in the order of `chain_rows`; empty for a
/// row that is not a protein chain.
pub type Outcome = Result<Vec<Vec<ExternalDomain>>, Error>;

/// Shortest chain worth sending (a V domain is at least ~90 residues).
const MIN_RESIDUES: usize = 80;

/// The one-letter sequence of each row that ANARCI should see, `None` for
/// the others.
pub fn chain_sequences(top: &Topology, rows: &[(String, Range<u32>)]) -> Vec<Option<String>> {
    rows.iter()
        .map(|(_, residues)| {
            let letters: String = residues.clone().map(|r| letter(top, r)).collect();
            let protein = top.residue_class(residues.start as usize) == ResidueClass::Protein;
            (protein && letters.len() >= MIN_RESIDUES).then_some(letters)
        })
        .collect()
}

fn letter(top: &Topology, residue: u32) -> char {
    let r = residue as usize;
    match top.residue_class(r) {
        ResidueClass::Protein => one_letter(top.residue_name(r))
            .filter(char::is_ascii_uppercase)
            .unwrap_or('X'),
        _ => 'X',
    }
}

/// Numbers the chains with the ANARCI found via `configured`, the
/// environment, the PATH or WSL, for every scheme.
pub fn number(configured: Option<&Path>, chains: &[Option<String>]) -> Outcome {
    if chains.iter().all(Option::is_none) {
        return Ok(vec![Vec::new(); chains.len()]);
    }
    let launcher = Launcher::locate(configured)?;
    number_with(|fasta, scheme| launcher.run(fasta, scheme), chains)
}

/// [`number`] with `run` standing in for the executable: FASTA and scheme
/// in, ANARCI's text out.
pub fn number_with(
    run: impl Fn(&str, Scheme) -> Result<String, Error>,
    chains: &[Option<String>],
) -> Outcome {
    let fasta = fasta_of(chains);
    let mut domains = vec![Vec::new(); chains.len()];
    // IMGT first: it finds every domain, including receptors, and the
    // other schemes add their labels to those.
    let order =
        std::iter::once(Scheme::Imgt).chain(Scheme::ALL.into_iter().filter(|s| *s != Scheme::Imgt));
    for scheme in order {
        let text = run(&fasta, scheme)?;
        let found = parse_anarci(&text).map_err(Error::Failed)?;
        merge(&mut domains, found, scheme);
    }
    Ok(domains)
}

fn fasta_of(chains: &[Option<String>]) -> String {
    let mut out = String::new();
    for (i, seq) in chains.iter().enumerate() {
        if let Some(seq) = seq {
            out.push_str(&format!(">c{i}\n{seq}\n"));
        }
    }
    out
}

fn chain_index(name: &str) -> Option<usize> {
    name.strip_prefix('c')?.parse().ok()
}

/// Adds `found` labels in `scheme` to the domain of the same chain type
/// that lies inside the found range; under IMGT, where every domain is
/// found, the domain is created.
fn merge(domains: &mut [Vec<ExternalDomain>], found: Vec<NumberedDomain>, scheme: Scheme) {
    for d in found {
        let Some(slot) = chain_index(&d.name).and_then(|i| domains.get_mut(i)) else {
            continue;
        };
        let inside = |e: &ExternalDomain| {
            e.chain == d.chain && d.range.start <= e.start && e.end <= d.range.end
        };
        match slot.iter_mut().find(|e| inside(e)) {
            Some(e) => {
                let own = (e.start - d.range.start)..(e.end - d.range.start);
                e.set_labels(scheme, d.labels[own].to_vec());
            }
            None if scheme == Scheme::Imgt => {
                let mut e = ExternalDomain::new(d.chain, d.range);
                e.set_labels(scheme, d.labels);
                slot.push(e);
            }
            None => {}
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use vv_core::antibody::{CdrDefinition, ChainType};

    const IMGT: &str = include_str!("../../../../vv-core/tests/data/anarci_imgt.txt");
    const KABAT: &str = include_str!("../../../../vv-core/tests/data/anarci_kabat.txt");

    /// What the executable printed for the same four chains: the Kabat
    /// sample stands in for every other antibody scheme.
    fn canned(_: &str, scheme: Scheme) -> Result<String, Error> {
        Ok(match scheme {
            Scheme::Imgt => IMGT.to_string(),
            _ => KABAT.to_string(),
        })
    }

    #[test]
    fn labels_from_each_scheme_land_on_the_domain_they_belong_to() {
        let chains = vec![Some("H".repeat(120)), Some("K".repeat(107)), None, None];
        let out = number_with(canned, &chains).unwrap();
        assert_eq!(out.len(), 4);
        let heavy = &out[0];
        assert_eq!(heavy.len(), 1);
        assert_eq!(heavy[0].chain, ChainType::Heavy);
        for scheme in Scheme::ALL {
            assert!(
                heavy[0].covers(scheme, CdrDefinition::Chothia),
                "{scheme:?}"
            );
        }
        let receptor = &out[3];
        assert_eq!(receptor.len(), 1);
        assert_eq!(receptor[0].chain, ChainType::TcrBeta);
        assert!(receptor[0].covers(Scheme::Kabat, CdrDefinition::Kabat));
    }

    #[test]
    fn a_scheme_whose_range_runs_past_imgts_contributes_the_shared_residues() {
        let longer = KABAT
            .replacen("|0|119|", "|0|120|", 1)
            .replacen("//", "H 149       S\n//", 1);
        let run = |_: &str, scheme: Scheme| {
            Ok(match scheme {
                Scheme::Imgt => IMGT.to_string(),
                Scheme::Aho => longer.clone(),
                _ => KABAT.to_string(),
            })
        };
        let chains = vec![Some("H".repeat(120)), Some("K".repeat(107)), None, None];
        let out = number_with(run, &chains).unwrap();
        let aho = out[0][0]
            .annotate(Scheme::Aho, CdrDefinition::Kabat)
            .unwrap();
        assert_eq!(aho.len(), 120);
        assert_eq!(out[0][0].end, 120);
    }

    #[test]
    fn a_failed_run_is_an_error_with_its_reason() {
        let chains = vec![Some("A".repeat(100))];
        let out = number_with(|_, _| Err(Error::Failed("boom".into())), &chains);
        assert_eq!(out.unwrap_err().to_string(), "ANARCI failed: boom");
    }

    #[test]
    fn chains_without_a_long_protein_sequence_skip_the_run() {
        let out = number(Some(Path::new("never-run")), &[None, None]).unwrap();
        assert_eq!(out.len(), 2);
    }

    #[test]
    fn the_not_found_message_says_how_to_fix_it() {
        assert_eq!(
            Error::NotFound.to_string(),
            "ANARCI not found: install it or set VIZVIZ_ANARCI"
        );
    }

    #[test]
    fn backend_words_round_trip() {
        for b in Backend::ALL {
            assert_eq!(Backend::parse(b.name()), Some(b));
        }
        assert_eq!(Backend::parse("anarci"), Some(Backend::Anarci));
        assert_eq!(Backend::parse("other"), None);
    }
}
