//! Reading Abnum's plain-text reply. The samples in `tests/data/` are the
//! service's real replies for trastuzumab (docs/ANTIBODY.md, "Scheme
//! authors' program as a backend").

use vv_core::antibody::abnum::{flag, parse_abnum, SCHEMES};
use vv_core::antibody::{find_domains, ChainType, Scheme};

const HEAVY_KABAT: &str = include_str!("data/abnum_heavy_kabat.txt");
const LIGHT_CHOTHIA: &str = include_str!("data/abnum_light_chothia.txt");
const NOT_ANTIBODY: &str = include_str!("data/abnum_error.txt");

const HEAVY: &str = "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS";
const KAPPA: &str = "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK";

fn labels(text: &str, query: &str) -> Vec<String> {
    let d = parse_abnum(query, text).unwrap().expect("a domain");
    d.labels.iter().map(ToString::to_string).collect()
}

#[test]
fn the_heavy_sample_reads_with_its_insertions() {
    let d = parse_abnum(HEAVY, HEAVY_KABAT).unwrap().unwrap();
    assert_eq!((d.chain, d.range.clone()), (ChainType::Heavy, 0..120));
    let l = labels(HEAVY_KABAT, HEAVY);
    assert_eq!(&l[50..53], ["51", "52", "52A"]);
    assert!(l.contains(&"82C".to_string()) && l.contains(&"100C".to_string()));
    assert_eq!(l.len(), 120);
}

#[test]
fn a_light_chain_is_kappa_or_lambda_by_the_native_profiles() {
    let d = parse_abnum(KAPPA, LIGHT_CHOTHIA).unwrap().unwrap();
    assert_eq!((d.chain, d.range.clone()), (ChainType::Kappa, 0..107));
    let native = find_domains(KAPPA);
    assert_eq!(native[0].chain, ChainType::Kappa);
}

#[test]
fn residues_before_the_domain_are_skipped_in_the_range() {
    let tagged = format!("GSHMAAA{HEAVY}ASTKGPSVFPL");
    let d = parse_abnum(&tagged, HEAVY_KABAT).unwrap().unwrap();
    assert_eq!(d.range, 7..127);
}

#[test]
fn the_error_comment_means_no_domain() {
    assert!(parse_abnum("AAAA", NOT_ANTIBODY).unwrap().is_none());
    assert!(parse_abnum("AAAA", "").unwrap().is_none());
}

#[test]
fn a_warning_comment_among_the_rows_is_ignored() {
    let with = format!("# Warning: unusual framework\n{HEAVY_KABAT}");
    assert_eq!(labels(&with, HEAVY).len(), 120);
}

#[test]
fn damaged_replies_are_errors() {
    assert!(parse_abnum(HEAVY, "<html>busy</html>\n").is_err());
    assert!(parse_abnum("DIQM", HEAVY_KABAT).is_err(), "not in the query");
    assert!(parse_abnum(HEAVY, "H1 E\nL2 V\n").is_err(), "mixed chains");
}

#[test]
fn only_kabat_chothia_and_martin_have_a_request_field() {
    let flags: Vec<_> = SCHEMES.iter().map(|s| flag(*s)).collect();
    assert_eq!(flags, [Some("-k"), Some("-c"), Some("-m")]);
    assert_eq!((flag(Scheme::Imgt), flag(Scheme::Aho)), (None, None));
}
