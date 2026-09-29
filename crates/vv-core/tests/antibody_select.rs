//! The `cdr` selection keyword on a hand-made trastuzumab heavy + light pair.

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
        _ => unreachable!(),
    }
}

fn fab() -> (Topology, Vec<Vec3>) {
    let mut b = TopologyBuilder::new();
    let mut serial = 0;
    for (chain, seq) in [("H", HEAVY), ("L", LIGHT)] {
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
