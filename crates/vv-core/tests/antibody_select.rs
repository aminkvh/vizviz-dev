//! The `cdr` selection keyword on a hand-made trastuzumab heavy + light pair.

use vv_core::antibody::{find_in_residues, CdrDefinition};
use vv_core::glam::Vec3;
use vv_core::{select, AtomRow, Element, Topology, TopologyBuilder};

const HEAVY: &str = "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS";
const LIGHT: &str = "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK";

fn three(letter: char) -> &'static str {
    match letter {
        'A' => "ALA",
        'C' => "CYS",
        'D' => "ASP",
        'E' => "GLU",
        'F' => "PHE",
        'G' => "GLY",
        'H' => "HIS",
        'I' => "ILE",
        'K' => "LYS",
        'L' => "LEU",
        'M' => "MET",
        'N' => "ASN",
        'P' => "PRO",
        'Q' => "GLN",
        'R' => "ARG",
        'S' => "SER",
        'T' => "THR",
        'V' => "VAL",
        'W' => "TRP",
        'Y' => "TYR",
        _ => "UNK",
    }
}

fn fab() -> (Topology, Vec<Vec3>) {
    topology_of(&[("H", HEAVY), ("L", LIGHT)])
}

/// One CA per residue, one chain per `(id, sequence)`.
fn topology_of(chains: &[(&str, &str)]) -> (Topology, Vec<Vec3>) {
    let mut b = TopologyBuilder::new();
    let mut serial = 0;
    for &(chain, seq) in chains {
        for (i, letter) in seq.chars().enumerate() {
            serial += 1;
            b.push(&AtomRow {
                element: Element::from_symbol(b"C"),
                name: *b" CA ",
                serial,
                alt_loc: 0,
                comp: three(letter),
                asym: chain,
                auth_asym: chain,
                seq_id: i as i32 + 1,
                auth_seq_id: i as i32 + 1,
                ins_code: 0,
                entity: 1,
                position: Vec3::new(serial as f32 * 3.8, 0.0, 0.0),
                occupancy: 1.0,
                b_factor: 0.0,
                charge: 0,
                hetero: false,
            });
        }
    }
    (b.topology, b.positions)
}

fn picked(expr: &str) -> String {
    let (t, p) = fab();
    let mask = select(&t, &p, expr).unwrap();
    let all: String = HEAVY.chars().chain(LIGHT.chars()).collect();
    mask.ones()
        .map(|atom| all.as_bytes()[atom] as char)
        .collect()
}

#[test]
fn each_cdr_is_the_published_kabat_sequence() {
    for (expr, seq) in [
        ("cdr h1", "DTYIH"),
        ("cdr h2", "RIYPTNGYTRYADSVKG"),
        ("cdr h3", "WGGDGFYAMDY"),
        ("cdr l1", "RASQDVNTAVA"),
        ("cdr l2", "SASFLYS"),
        ("cdr l3", "QQHYTTPPT"),
    ] {
        assert_eq!(picked(expr), seq, "{expr}");
    }
}

#[test]
fn bare_cdr_and_chain_letters_pick_groups() {
    assert_eq!(
        picked("cdr h"),
        ["DTYIH", "RIYPTNGYTRYADSVKG", "WGGDGFYAMDY"].concat()
    );
    assert_eq!(
        picked("cdr"),
        [
            "DTYIH",
            "RIYPTNGYTRYADSVKG",
            "WGGDGFYAMDY",
            "RASQDVNTAVA",
            "SASFLYS",
            "QQHYTTPPT"
        ]
        .concat()
    );
    assert_eq!(picked("cdr h3 l3"), ["WGGDGFYAMDY", "QQHYTTPPT"].concat());
}

#[test]
fn a_definition_word_changes_the_loops_and_composes_with_other_terms() {
    assert_eq!(picked("cdr chothia h1"), "GFNIKDT");
    assert_eq!(picked("cdr imgt l3"), "QQHYTTPPT");
    assert_eq!(picked("cdr h3 and chain L"), "");
    assert_eq!(picked("cdr h3 and not resname TYR"), "WGGDGFAMD");
}

#[test]
fn a_bad_cdr_option_is_an_error_naming_it() {
    let (t, p) = fab();
    let err = select(&t, &p, "cdr h4").unwrap_err();
    assert!(err.message.contains("h4"), "{}", err.message);
}

#[test]
fn the_domains_are_found_once_per_topology_and_not_carried_by_a_clone() {
    let (t, p) = fab();
    assert!(format!("{:?}", t.antibody).contains("empty"));
    let first = select(&t, &p, "cdr h3").unwrap();
    assert!(format!("{:?}", t.antibody).contains("filled"));
    assert_eq!(select(&t, &p, "cdr h3").unwrap(), first);
    assert!(format!("{:?}", t.clone().antibody).contains("empty"));
}

fn bench_chains(id: &str) -> Vec<(String, String)> {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../fixtures/real/bench")
        .join(format!("{id}.tsv"));
    std::fs::read_to_string(path)
        .unwrap_or_default()
        .lines()
        .filter_map(|l| l.split_once('\t'))
        .map(|(c, s)| (c.to_string(), s.to_string()))
        .collect()
}

/// Detection re-run per query, as `cdr` did before the domains were cached.
fn uncached_cdr_query(t: &Topology) -> usize {
    use rayon::prelude::*;
    t.chains
        .par_iter()
        .map(|chain| {
            find_in_residues(t, chain.residues.clone())
                .iter()
                .map(|d| {
                    d.annotate(CdrDefinition::Kabat.native_scheme(), CdrDefinition::Kabat)
                        .iter()
                        .filter(|a| a.region.cdr().is_some())
                        .count()
                })
                .sum::<usize>()
        })
        .sum()
}

#[test]
#[ignore = "timing; needs fixtures/real/bench (Fab-decorated capsids); use --release"]
fn cdr_selection_timing() {
    for id in ["12US", "4UDF"] {
        let chains = bench_chains(id);
        if chains.is_empty() {
            continue;
        }
        let refs: Vec<(&str, &str)> = chains
            .iter()
            .map(|(c, s)| (c.as_str(), s.as_str()))
            .collect();
        let (t, p) = topology_of(&refs);
        let ms = |f: &mut dyn FnMut()| {
            let start = std::time::Instant::now();
            f();
            start.elapsed().as_secs_f64() * 1e3
        };
        let uncached = ms(&mut || {
            uncached_cdr_query(&t);
        });
        let mut picked = 0;
        let first = ms(&mut || picked = select(&t, &p, "cdr h3").unwrap().count_ones(..));
        let repeat = ms(&mut || {
            select(&t, &p, "cdr h3").unwrap();
        });
        println!(
            "{id}: {} chains, {} residues, {picked} atoms in CDR-H3; per query: no cache {uncached:.1} ms, first cached {first:.1} ms, repeat {repeat:.2} ms",
            chains.len(),
            t.residue_count()
        );
    }
}
