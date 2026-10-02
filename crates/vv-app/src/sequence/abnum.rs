//! The optional Abnum numbering backend: the public web service of the
//! scheme authors' program (Abhinandan & Martin 2008, Mol Immunol 45:3832).
//! Each long protein chain goes out once per scheme, one request at a time,
//! and the labels come back as [`ExternalDomain`]s like ANARCI's. The
//! service numbers one domain per sequence.

use vv_core::antibody::abnum::{flag, parse_abnum, SCHEMES};
use vv_core::antibody::external::ExternalDomain;
use vv_io::fetch::FetchError;

use super::anarci::{Error, Outcome};

/// Shown while the backend is selected: the site's own caveat.
pub const PLAIN_HTTP_NOTICE: &str =
    "Sequences are sent over plain HTTP to the public server at bioinf.org.uk";

/// Numbers the chains (see [`super::anarci::chain_sequences`]) with the
/// public server, or with `server` (a URL) in its place.
pub fn number(server: Option<&str>, chains: &[Option<String>]) -> Outcome {
    let Some(url) = server else {
        return number_with(vv_io::seqdata::abnum::number, chains);
    };
    let cache = std::env::temp_dir().join("vizviz-abnum-other-server");
    number_with(
        |sequence, field| vv_io::seqdata::abnum::number_at(url, &cache, sequence, field),
        chains,
    )
}

/// [`number`] with `run` standing in for the server: sequence and scheme
/// field in, the plain-text reply out.
pub fn number_with(
    run: impl Fn(&str, &str) -> Result<String, FetchError>,
    chains: &[Option<String>],
) -> Outcome {
    chains
        .iter()
        .map(|chain| match chain {
            Some(sequence) => number_chain(&run, sequence),
            None => Ok(Vec::new()),
        })
        .collect()
}

/// The chain's domain in every scheme; none if the first reply finds no
/// domain, which spares the other two requests.
fn number_chain(
    run: &impl Fn(&str, &str) -> Result<String, FetchError>,
    sequence: &str,
) -> Result<Vec<ExternalDomain>, Error> {
    let mut domain: Option<ExternalDomain> = None;
    for scheme in SCHEMES {
        let field = flag(scheme).expect("SCHEMES are the offered ones");
        let reply = run(sequence, field).map_err(|_| Error::Unreachable)?;
        let Some(found) = parse_abnum(sequence, &reply).map_err(Error::Failed)? else {
            return Ok(Vec::new());
        };
        let slot = domain.get_or_insert_with(|| ExternalDomain::new(found.chain, found.range));
        let fits = slot.set_labels(scheme, found.labels);
        if !fits {
            return Err(Error::Failed(format!("{} domain differs", scheme.name())));
        }
    }
    Ok(domain.into_iter().collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::cell::RefCell;
    use vv_core::antibody::{CdrDefinition, ChainType, Scheme};

    const HEAVY: &str = "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS";
    const REPLY: &str = include_str!("../../../vv-core/tests/data/abnum_heavy_kabat.txt");
    const ERROR: &str = include_str!("../../../vv-core/tests/data/abnum_error.txt");

    #[test]
    fn every_offered_scheme_is_requested_once_per_chain() {
        let asked = RefCell::new(Vec::new());
        let run = |_: &str, field: &str| {
            asked.borrow_mut().push(field.to_string());
            Ok(REPLY.to_string())
        };
        let chains = vec![Some(format!("GS{HEAVY}")), None];
        let out = number_with(run, &chains).unwrap();
        assert_eq!(*asked.borrow(), ["-k", "-c", "-m"]);
        assert_eq!(out[0].len(), 1);
        assert!(out[1].is_empty());
        let d = &out[0][0];
        assert_eq!((d.chain, d.start, d.end), (ChainType::Heavy, 2, 122));
        for scheme in [Scheme::Kabat, Scheme::Chothia, Scheme::Martin] {
            assert!(d.covers(scheme, CdrDefinition::Kabat), "{scheme:?}");
        }
        assert!(!d.covers(Scheme::Imgt, CdrDefinition::Kabat));
    }

    #[test]
    fn a_chain_without_a_domain_stops_after_one_request() {
        let asked = RefCell::new(0);
        let run = |_: &str, _: &str| {
            *asked.borrow_mut() += 1;
            Ok(ERROR.to_string())
        };
        let out = number_with(run, &[Some("A".repeat(100))]).unwrap();
        assert!(out[0].is_empty());
        assert_eq!(*asked.borrow(), 1);
    }

    #[test]
    fn an_unanswered_request_is_unreachable_not_a_partial_result() {
        let run = |_: &str, _: &str| -> Result<String, FetchError> {
            Err(FetchError::BadRequest("down".into()))
        };
        let out = number_with(run, &[Some(HEAVY.into()), Some(HEAVY.into())]);
        assert_eq!(out.unwrap_err().to_string(), "Abnum unreachable");
    }
}
