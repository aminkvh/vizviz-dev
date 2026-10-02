//! Antibody numbering on constructed and embedded sequences. Real-structure
//! validation lives in `antibody_pdb.rs`.

use std::collections::HashSet;

mod antibody_common;
use antibody_common::{fingerprint, in_pool};

use proptest::prelude::*;
use rayon::prelude::*;
use vv_core::antibody::{
    find_domains, find_domains_in_chain, find_variable_domains, Annotation, CdrDefinition,
    ChainType, Domain, Label, Region, Scheme,
};

// Trastuzumab Fab (wwPDB 1N8Z): framework/CDR split at the IMGT boundaries.
const H_FR1: &str = "EVQLVESGGGLVQPGGSLRLSCAAS";
const H_CDR1: &str = "GFNIKDTY";
const H_FR2: &str = "IHWVRQAPGKGLEWVAR";
const H_CDR2: &str = "IYPTNGYT";
const H_FR3: &str = "RYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYC";
const H_CDR3: &str = "SRWGGDGFYAMDY";
const H_FR4: &str = "WGQGTLVTVSS";
const L_FR1: &str = "DIQMTQSPSSLSASVGDRVTITCRAS";
const L_CDR1: &str = "QDVNTA";
const L_FR2: &str = "VAWYQQKPGKAPKLLIY";
const L_CDR2: &str = "SAS";
const L_FR3: &str = "FLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYC";
const L_CDR3: &str = "QQHYTTPPT";
const L_FR4: &str = "FGQGTKVEIK";

const CH1: &str = "ASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSC";
const CK: &str = "RTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC";

fn heavy(cdr1: &str, cdr2: &str, cdr3: &str) -> String {
    [H_FR1, cdr1, H_FR2, cdr2, H_FR3, cdr3, H_FR4].concat()
}

fn light(cdr1: &str, cdr2: &str, cdr3: &str) -> String {
    [L_FR1, cdr1, L_FR2, cdr2, L_FR3, cdr3, L_FR4].concat()
}

fn trastuzumab_heavy() -> String {
    heavy(H_CDR1, H_CDR2, H_CDR3)
}

fn trastuzumab_light() -> String {
    light(L_CDR1, L_CDR2, L_CDR3)
}

#[test]
fn numbering_is_identical_every_run_and_thread_count() {
    let mut sequences = vec![trastuzumab_heavy(), trastuzumab_light()];
    sequences.push([trastuzumab_light(), CK.to_string(), trastuzumab_heavy()].concat());
    sequences.push(heavy("GFNIKDTYGGGG", "IYPTNGYT", "SRWGGDGFYAMDYAAAA"));
    sequences.push(light("QDVNTAAAA", "SAS", "QQHYTTPPTGG"));
    let serial: Vec<_> = sequences.iter().map(|s| fingerprint(s)).collect();
    assert!(serial == sequences.iter().map(|s| fingerprint(s)).collect::<Vec<_>>());
    for threads in [1, 4] {
        let pooled = in_pool(threads, || {
            sequences
                .par_iter()
                .map(|s| fingerprint(s))
                .collect::<Vec<_>>()
        });
        assert!(serial == pooled, "{threads} threads");
    }
}

#[test]
fn deposited_numbering_keeps_holes_and_skips_tags() {
    let full = ["HHHHHH", &trastuzumab_heavy()].concat();
    let whole = &find_domains(&full)[0];
    let keep: Vec<usize> = (0..full.len()).filter(|i| !(65..69).contains(i)).collect();
    let observed: String = keep.iter().map(|&i| full.as_bytes()[i] as char).collect();
    let found = find_domains_in_chain(&full, &observed);
    assert_eq!(found.len(), 1);
    let expected: Vec<(usize, Label)> = whole
        .numbering(Scheme::Kabat)
        .into_iter()
        .filter_map(|(i, l)| Some((keep.iter().position(|&k| k == i)?, l)))
        .collect();
    assert_eq!(found[0].numbering(Scheme::Kabat), expected);
    assert!(found[0].start >= 6, "the tag is not part of the domain");
}

#[test]
fn deposited_numbering_falls_back_when_the_sequences_disagree() {
    let heavy = trastuzumab_heavy();
    let found = find_domains_in_chain("MKTAYIAKQR", &heavy);
    assert_eq!(found.len(), 1);
    assert_eq!(found[0].end, find_domains(&heavy)[0].end);
}

fn only_domain(seq: &str) -> Domain {
    let mut found = find_domains(seq);
    assert_eq!(found.len(), 1, "expected one domain in {seq}");
    found.remove(0)
}

fn residues(seq: &str, notes: &[Annotation], region: Region) -> String {
    notes
        .iter()
        .filter(|a| a.region == region)
        .map(|a| seq.as_bytes()[a.index] as char)
        .collect()
}

fn cdrs(seq: &str, scheme: Scheme, def: CdrDefinition) -> [String; 3] {
    let notes = only_domain(seq).annotate(scheme, def);
    [Region::Cdr1, Region::Cdr2, Region::Cdr3].map(|r| residues(seq, &notes, r))
}

fn labels(seq: &str, scheme: Scheme) -> Vec<String> {
    only_domain(seq)
        .numbering(scheme)
        .iter()
        .map(|(_, l)| l.to_string())
        .collect()
}

/// Labels of the residues at `start..start + len`, space separated.
fn slice(seq: &str, scheme: Scheme, start: usize, len: usize) -> String {
    labels(seq, scheme)[start..start + len].join(" ")
}

#[test]
fn trastuzumab_kabat_cdrs_match_the_published_sequences() {
    let h = trastuzumab_heavy();
    assert_eq!(
        cdrs(&h, Scheme::Kabat, CdrDefinition::Kabat),
        ["DTYIH", "RIYPTNGYTRYADSVKG", "WGGDGFYAMDY"]
    );
    let l = trastuzumab_light();
    assert_eq!(
        cdrs(&l, Scheme::Kabat, CdrDefinition::Kabat),
        ["RASQDVNTAVA", "SASFLYS", "QQHYTTPPT"]
    );
}

#[test]
fn trastuzumab_cdrs_under_every_definition() {
    let h = trastuzumab_heavy();
    let cases = [
        (CdrDefinition::Chothia, ["GFNIKDT", "YPTNGY", "WGGDGFYAMDY"]),
        (
            CdrDefinition::Imgt,
            ["GFNIKDTY", "IYPTNGYT", "SRWGGDGFYAMDY"],
        ),
        (
            CdrDefinition::Contact,
            ["KDTYIH", "WVARIYPTNGYTR", "SRWGGDGFYAMD"],
        ),
        (
            CdrDefinition::North,
            ["AASGFNIKDTYIH", "RIYPTNGYTR", "WGGDGFYAMDY"],
        ),
    ];
    for (def, expected) in cases {
        assert_eq!(cdrs(&h, Scheme::Kabat, def), expected, "heavy {def:?}");
    }
    let l = trastuzumab_light();
    let cases = [
        (
            CdrDefinition::Chothia,
            ["RASQDVNTAVA", "SASFLYS", "QQHYTTPPT"],
        ),
        (CdrDefinition::Imgt, ["QDVNTA", "SAS", "QQHYTTPPT"]),
        (
            CdrDefinition::Contact,
            ["NTAVAWY", "LLIYSASFLY", "QQHYTTPP"],
        ),
        (
            CdrDefinition::North,
            ["RASQDVNTAVA", "SASFLYS", "QQHYTTPPT"],
        ),
    ];
    for (def, expected) in cases {
        assert_eq!(cdrs(&l, Scheme::Kabat, def), expected, "light {def:?}");
    }
}

#[test]
fn trastuzumab_framework_anchors_in_every_scheme() {
    let h = trastuzumab_heavy();
    let d = only_domain(&h);
    assert_eq!((d.chain, d.start, d.end), (ChainType::Heavy, 0, 120));
    // Residues C22 (Kabat) = IMGT 23, W36 = IMGT 41, C92 = IMGT 104, W103 = IMGT 118.
    let anchors = [
        h.find('C').unwrap(),
        h.find("WVRQ").unwrap(),
        h.find("YYC").unwrap() + 2,
        h.rfind("WGQG").unwrap(),
    ];
    for (scheme, expected) in [
        (Scheme::Kabat, ["22", "36", "92", "103"]),
        (Scheme::Chothia, ["22", "36", "92", "103"]),
        (Scheme::Martin, ["22", "36", "92", "103"]),
        (Scheme::Imgt, ["23", "41", "104", "118"]),
    ] {
        let n = d.numbering(scheme);
        let got = anchors.map(|i| n[i].1.to_string());
        assert_eq!(got, expected, "{scheme:?}");
    }
}

#[test]
fn imgt_cdr_loops_gap_and_grow_symmetrically() {
    let cases = [
        // CDR1 (27..=38): odd lengths keep the extra residue on the N side.
        (0, 7, "27 28 29 30 36 37 38"),
        (0, 8, "27 28 29 30 35 36 37 38"),
        (0, 12, "27 28 29 30 31 32 33 34 35 36 37 38"),
        // CDR2 (56..=65).
        (1, 7, "56 57 58 59 63 64 65"),
        (1, 10, "56 57 58 59 60 61 62 63 64 65"),
        // CDR3 (105..=117): gaps 111, 112, 110, 113; insertions 112.1, 111.1.
        (2, 13, "105 106 107 108 109 110 111 112 113 114 115 116 117"),
        (2, 12, "105 106 107 108 109 110 112 113 114 115 116 117"),
        (2, 11, "105 106 107 108 109 110 113 114 115 116 117"),
        (
            2,
            14,
            "105 106 107 108 109 110 111 112A 112 113 114 115 116 117",
        ),
        (
            2,
            15,
            "105 106 107 108 109 110 111 111A 112A 112 113 114 115 116 117",
        ),
    ];
    for (loop_id, n, expected) in cases {
        let (cdr1, cdr2, cdr3) = match loop_id {
            0 => ("A".repeat(n), H_CDR2.to_string(), H_CDR3.to_string()),
            1 => (H_CDR1.to_string(), "A".repeat(n), H_CDR3.to_string()),
            _ => (H_CDR1.to_string(), H_CDR2.to_string(), "A".repeat(n)),
        };
        let seq = heavy(&cdr1, &cdr2, &cdr3);
        let start = match loop_id {
            0 => H_FR1.len(),
            1 => H_FR1.len() + H_CDR1.len() + H_FR2.len(),
            _ => seq.len() - H_FR4.len() - n,
        };
        assert_eq!(
            slice(&seq, Scheme::Imgt, start, n),
            expected,
            "loop {loop_id} length {n}"
        );
    }
}

#[test]
fn imgt_numbers_skip_the_germline_gap_columns() {
    let h = labels(&trastuzumab_heavy(), Scheme::Imgt);
    assert_eq!(&h[8..11], ["9", "11", "12"], "heavy FR1 has no 10");
    let at_gap = h.iter().position(|l| l == "72").unwrap();
    assert_eq!(&h[at_gap..at_gap + 2], ["72", "74"], "heavy FR3 has no 73");
    let k = labels(&trastuzumab_light(), Scheme::Imgt);
    let at_gap = k.iter().position(|l| l == "72").unwrap();
    assert_eq!(&k[at_gap..at_gap + 2], ["72", "74"]);
    let at_gap = k.iter().position(|l| l == "80").unwrap();
    assert_eq!(
        &k[at_gap..at_gap + 2],
        ["80", "83"],
        "kappa FR3 has no 81, 82"
    );
}

fn cdr3_labels(seq: &str, scheme: Scheme, n: usize) -> Vec<String> {
    let all = labels(seq, scheme);
    let end = all.len() - H_FR4.len();
    all[end - n..end].to_vec()
}

#[test]
fn kabat_heavy_cdr3_insertions_and_deletions() {
    let long = cdr3_labels(&heavy(H_CDR1, H_CDR2, &"A".repeat(30)), Scheme::Kabat, 30);
    assert_eq!(long[..8].join(" "), "93 94 95 96 97 98 99 100");
    assert_eq!(long[8..28].first().map(String::as_str), Some("100A"));
    assert_eq!(long[27], "100T", "twenty insertions run 100A..100T");
    assert_eq!(long[28..].join(" "), "101 102");
    let short = cdr3_labels(&heavy(H_CDR1, H_CDR2, "ARDYAF"), Scheme::Kabat, 6);
    assert_eq!(short.join(" "), "93 94 95 96 101 102");
}

#[test]
fn kabat_heavy_cdr2_and_framework_insertions() {
    let seq_16 = labels(&heavy(H_CDR1, "IYYSGST", H_CDR3), Scheme::Kabat);
    assert!(
        !seq_16.contains(&"52A".to_string()),
        "sixteen residues fill 50..=65"
    );
    assert!(labels(&trastuzumab_heavy(), Scheme::Kabat).contains(&"52A".to_string()));
    let seq_19 = labels(&heavy(H_CDR1, "IYPTNGYTAAA", H_CDR3), Scheme::Kabat);
    assert!(["52A", "52B", "52C"]
        .iter()
        .all(|x| seq_19.contains(&x.to_string())));
    let fr3_82 = seq_19.iter().filter(|x| x.starts_with("82")).count();
    assert_eq!(fr3_82, 4, "82, 82A, 82B, 82C");
}

#[test]
fn chothia_and_kabat_place_h1_and_l1_insertions_differently() {
    let h = heavy("GFNIKDTYSS", H_CDR2, H_CDR3);
    let kabat = labels(&h, Scheme::Kabat);
    let chothia = labels(&h, Scheme::Chothia);
    for x in ["35A", "35B"] {
        assert!(kabat.contains(&x.to_string()) && !chothia.contains(&x.to_string()));
    }
    for x in ["31A", "31B"] {
        assert!(chothia.contains(&x.to_string()) && !kabat.contains(&x.to_string()));
    }
    let l = light("QSLLYSSNQKNYLA", L_CDR2, L_CDR3);
    let kabat = labels(&l, Scheme::Kabat);
    let chothia = labels(&l, Scheme::Chothia);
    assert!(
        ["27A", "27F"]
            .iter()
            .all(|x| kabat.contains(&x.to_string()))
            && !kabat.contains(&"30A".to_string())
    );
    assert!(
        ["30A", "30F"]
            .iter()
            .all(|x| chothia.contains(&x.to_string()))
            && !chothia.contains(&"27A".to_string())
    );
}

#[test]
fn short_light_cdr1_deletes_position_28() {
    let l = light("QSVSY", L_CDR2, L_CDR3);
    let n = labels(&l, Scheme::Kabat);
    let cdr1 = &n[L_FR1.len() - 3..L_FR1.len() + 5 + 2];
    assert_eq!(cdr1.join(" "), "24 25 26 27 29 30 31 32 33 34");
}

#[test]
fn martin_moves_the_light_cdr2_insertion_site() {
    let l = light(L_CDR1, "SASFLYSAA", L_CDR3);
    let chothia = labels(&l, Scheme::Chothia);
    let martin = labels(&l, Scheme::Martin);
    assert!(chothia.contains(&"54A".to_string()));
    assert!(martin.contains(&"52A".to_string()));
}

/// Labels from `first` to `last` inclusive, as a space-separated string.
fn between(all: &[String], first: &str, last: &str) -> String {
    let a = all.iter().position(|x| x == first).unwrap();
    let b = all.iter().position(|x| x == last).unwrap();
    all[a..=b].join(" ")
}

#[test]
fn martin_heavy_framework_three_inserts_at_h72_where_chothia_inserts_at_h82() {
    let h = trastuzumab_heavy();
    let martin = labels(&h, Scheme::Martin);
    let chothia = labels(&h, Scheme::Chothia);
    let expected = |ins: &str, at: &str| {
        let mut out: Vec<String> = (66..=92).map(|n| n.to_string()).collect();
        let pos = out.iter().position(|x| x == at).unwrap();
        for (k, letter) in ["A", "B", "C"].iter().enumerate() {
            out.insert(pos + 1 + k, format!("{ins}{letter}"));
        }
        out.join(" ")
    };
    assert_eq!(between(&martin, "66", "92"), expected("72", "72"));
    assert_eq!(between(&chothia, "66", "92"), expected("82", "82"));
}

#[test]
fn martin_light_framework_sites() {
    let kappa = trastuzumab_light();
    assert_eq!(
        labels(&kappa, Scheme::Martin),
        labels(&kappa, Scheme::Chothia)
    );
    let lambda = "QSVLTQPPSASGTPGQRVTISCSGSSSNIGSNTVNWYQQLPGTAPKLLIYSNNQRPSGVPDRFSGSKSGTSASLAISGLQSEDEADYYCAAWDDSLNGVVFGGGTKLTVLG";
    let (martin, chothia) = (
        labels(lambda, Scheme::Martin),
        labels(lambda, Scheme::Chothia),
    );
    assert_eq!(
        between(&martin, "1", "11"),
        "1 2 3 4 5 6 8 9 10 11",
        "lambda FR1 gap at L7"
    );
    assert_eq!(
        between(&chothia, "1", "11"),
        "1 2 3 4 5 6 7 8 9 11",
        "Chothia gap at L10"
    );
    assert_eq!(between(&martin, "105", "107"), "105 106 106A 107");
}

#[test]
fn deleted_loop_labels_follow_each_schemes_own_site() {
    let dropped = |seq: &str, scheme: Scheme, lo: u16, hi: u16| -> Vec<u16> {
        let l = labels(seq, scheme);
        (lo..=hi).filter(|n| !l.contains(&n.to_string())).collect()
    };
    let h = heavy("GFNIKDT", H_CDR2, H_CDR3);
    assert_eq!(dropped(&h, Scheme::Kabat, 26, 35), [35]);
    assert_eq!(dropped(&h, Scheme::Chothia, 26, 35), [31]);
    assert_eq!(dropped(&h, Scheme::Martin, 26, 35), [31]);
    let l = light("QDVNA", L_CDR2, L_CDR3);
    assert_eq!(dropped(&l, Scheme::Kabat, 24, 34), [28]);
    assert_eq!(dropped(&l, Scheme::Chothia, 24, 34), [31]);
    assert_eq!(dropped(&l, Scheme::Martin, 24, 34), [30]);
    let l2 = light(L_CDR1, "SA", L_CDR3);
    assert_eq!(dropped(&l2, Scheme::Chothia, 50, 56), [54]);
    assert_eq!(dropped(&l2, Scheme::Martin, 50, 56), [52]);
}

// Human T-cell receptors A6 (wwPDB 1AO7) and 1MI5: variable domain and the
// start of the constant domain.
const A6_ALPHA: &str = "KEVEQNSGPLSVPEGAIASLNCTYSDRGSQSFFWYRQYSGKSPELIMSIYSNGDKEDGRFTAQLNKASQYVSLLIRDSQPSDSATYLCAVTTDSWGKLQFGAGTQVVVTPDIQNPDPAVYQLRD";
const A6_BETA: &str = "NAGVTQTPKFQVLKTGQSMTLQCAQDMNHEYMSWYRQDPGMGLRLIHYSVGAGITDQGEVPNGYNVSRSTTEDFPLRLLSAAPSQTSVYFCASRPGLAGGRPEQYFGPGTRLTVTEDLKNVFPPEVAVFEPSE";
const MI5_ALPHA: &str = "KTTQPNSMESNEEEPVHLPCNHSTISGTDYIHWYRQLPSQGPEYVIHGLTSNVNNRMASLAIAEDRKSSTLILHRATLRDAAVYYCILPLAGGTSYGKLTFGQGTILTVHPNIQNPDPAVYQLRDSKSSDKSVCL";
const MI5_BETA: &str = "GVSQSPRYKVAKRGQDVALRCDPISGHVSLFWYQQALGQGPEFLTYFQNEAQLDKSGLPSDRFFAERPEGSVSTLKIQRTQQEDSAVYLCASSLGQAYEQYFGPGTRLTVTEDLKNVFPPEVAVFEPSE";

fn only_receptor(seq: &str) -> Domain {
    let mut found = find_variable_domains(seq);
    assert_eq!(found.len(), 1, "expected one domain in {seq}");
    found.remove(0)
}

#[test]
fn receptor_chains_are_found_as_receptors_and_not_as_antibodies() {
    for (seq, chain) in [
        (A6_ALPHA, ChainType::TcrAlpha),
        (A6_BETA, ChainType::TcrBeta),
        (MI5_ALPHA, ChainType::TcrAlpha),
        (MI5_BETA, ChainType::TcrBeta),
    ] {
        assert_eq!(only_receptor(seq).chain, chain);
        assert!(find_domains(seq).is_empty());
    }
    for (seq, chain) in [
        (trastuzumab_heavy(), ChainType::Heavy),
        (trastuzumab_light(), ChainType::Kappa),
    ] {
        assert_eq!(only_receptor(&seq).chain, chain);
    }
}

#[test]
fn receptor_domains_carry_imgt_anchors_and_the_published_cdrs() {
    let at = |seq: &str, n: u16| {
        let labels = only_receptor(seq).numbering(Scheme::Imgt);
        let (i, _) = labels
            .iter()
            .find(|(_, l)| l.number == n && l.insertion().is_none())
            .unwrap();
        seq.as_bytes()[*i] as char
    };
    for seq in [A6_ALPHA, A6_BETA, MI5_ALPHA, MI5_BETA] {
        let anchors: String = [23, 41, 104, 118, 119]
            .iter()
            .map(|&n| at(seq, n))
            .collect();
        assert_eq!(anchors, "CWCFG", "{seq}");
    }
    let imgt = |seq: &str| {
        let notes = only_receptor(seq).annotate(Scheme::Imgt, CdrDefinition::Imgt);
        [Region::Cdr1, Region::Cdr2, Region::Cdr3].map(|r| residues(seq, &notes, r))
    };
    assert_eq!(imgt(A6_ALPHA), ["DRGSQS", "IYSNGD", "AVTTDSWGKLQ"]);
    assert_eq!(imgt(A6_BETA), ["MNHEY", "SVGAGI", "ASRPGLAGGRPEQY"]);
}

#[test]
fn lambda_chains_are_classified_and_numbered() {
    // Human IGLV1-44 with a J segment and the start of the constant domain.
    let seq = "QSVLTQPPSASGTPGQRVTISCSGSSSNIGSNTVNWYQQLPGTAPKLLIYSNNQRPSGVPDRFSGSKSGTSASLAISGLQSEDEADYYCAAWDDSLNGWVFGGGTKLTVLGQPKAAPSVTLFPPSS";
    let d = only_domain(seq);
    assert_eq!(d.chain, ChainType::Lambda);
    let k = labels(seq, Scheme::Kabat);
    assert_eq!(
        &k[8..11],
        ["9", "11", "12"],
        "lambda numbering has no position 10"
    );
    let tail = &k[d.end - 3..];
    assert_eq!(
        tail.join(" "),
        "106 106A 107",
        "lambda FR4 ends V106 L106A G107"
    );
    let cdrs = cdrs(seq, Scheme::Kabat, CdrDefinition::Kabat);
    assert_eq!(cdrs, ["SGSSSNIGSNTVN", "SNNQRPS", "AAWDDSLNGWV"]);
    let kabat_l1: Vec<String> = k[22..22 + 13].to_vec();
    assert_eq!(
        kabat_l1.join(" "),
        "24 25 26 27 27A 27B 28 29 30 31 32 33 34"
    );
}

#[test]
fn nanobody_domain_is_heavy() {
    // Camelid VHH, wwPDB 1MEL chain A.
    let seq = "DVQLQASGGGSVQAGGSLRLSCAASGYTIGPYCMGWFRQAPGKEREGVAAINMGGGITYYADSVKGRFTISQDNAKNTVYLLMNSLEPEDTAIYYCAADSTIYASYYECGHGLSTGGYGYDSWGQGTQVTVSSGR";
    let d = only_domain(seq);
    assert_eq!(d.chain, ChainType::Heavy);
    let cdr3 = &cdrs(seq, Scheme::Kabat, CdrDefinition::Kabat)[2];
    assert_eq!(cdr3, "DSTIYASYYECGHGLSTGGYGYDS");
}

#[test]
fn fab_chains_stop_at_the_constant_domain() {
    let heavy_chain = [trastuzumab_heavy(), CH1.to_string()].concat();
    let d = only_domain(&heavy_chain);
    assert_eq!(d.end, trastuzumab_heavy().len());
    let light_chain = [trastuzumab_light(), CK.to_string()].concat();
    let d = only_domain(&light_chain);
    assert_eq!(d.chain, ChainType::Kappa);
    // IMGT FR4 has eleven columns, so kappa takes the first constant residue.
    assert_eq!(d.end, trastuzumab_light().len() + 1);
}

#[test]
fn truncated_ends_are_still_numbered() {
    let h = trastuzumab_heavy();
    let cut = &h[4..h.len() - 4];
    let d = only_domain(cut);
    assert_eq!((d.start, d.end), (0, cut.len()));
    let n = d.numbering(Scheme::Kabat);
    assert_eq!(n[0].1.to_string(), "5");
    assert_eq!(n.last().unwrap().1.to_string(), "109");
}

#[test]
fn scfv_yields_both_domains_in_order() {
    let scfv = [
        trastuzumab_heavy(),
        "GGGGSGGGGSGGGGS".to_string(),
        trastuzumab_light(),
    ]
    .concat();
    let found = find_domains(&scfv);
    assert_eq!(found.len(), 2);
    assert_eq!(
        (found[0].chain, found[0].start, found[0].end),
        (ChainType::Heavy, 0, 120)
    );
    assert_eq!(found[1].chain, ChainType::Kappa);
    assert_eq!(found[1].start, 135);
    assert_eq!(found[1].end, 135 + 107);
}

#[test]
fn non_antibody_chains_return_nothing() {
    let negatives = [
        // Hen egg-white lysozyme
        "KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL",
        // Human beta-2 microglobulin
        "MIQRTPKIQVYSRHPAENGKSNFLNCYVSGFHPSDIEVDLLKNGERIEKVEHSDLSFSKDWSFYLLYYTEFTPTEKDEYACRVNHVTLSQPKIVKWDRDM",
        // TCR alpha and beta variable domains with constant regions (1AO7)
        "KEVEQNSGPLSVPEGAIASLNCTYSDRGSQSFFWYRQYSGKSPELIMSIYSNGDKEDGRFTAQLNKASQYVSLLIRDSQPSDSATYLCAVTTDSWGKLQFGAGTQVVVTPDIQNPDPAVYQLRDSKSSDKSVCLFTDFDSQTNVSQSKDSDVYITDKCVLDMRSMDFKSNSAVAWSNKSDFACANAFNNSIIPEDTFFPS",
        "NAGVTQTPKFQVLKTGQSMTLQCAQDMNHEYMSWYRQDPGMGLRLIHYSVGAGITDQGEVPNGYNVSRSTTEDFPLRLLSAAPSQTSVYFCASRPGLAGGRPEQYFGPGTRLTVTEDLKNVFPPEVAVFEPSEAEISHTQKATLVCLATGFYPDHVELSWWVNGKEVHSGVCTDPQPLKEQPALNDSRYALSSRLRVSATFWQNPRNHFRCQVQFYGLSENDEWTQDRAKPVTQIVSAEAWGRAD",
        CH1,
        CK,
        // Human IgG1 Fc CH2-CH3
        "APELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSRDELTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGK",
        // Human PD-1 IgV domain
        "NPPTFSPALLVVTEGDNATFTCSFSNTSESFVLNWYRMSPSNQTDKLAAFPEDRSQPGQDCRFRVTQLPNGRDFHMSVVRARRNDSGTYLCGAISLAPKAQIKESLRAELRVTERRAE",
    ];
    for s in negatives {
        assert!(find_domains(s).is_empty(), "false positive on {s}");
    }
}

#[test]
fn shuffled_antibodies_are_rejected() {
    let mut seq: Vec<u8> = trastuzumab_heavy().into_bytes();
    let mut state = 12345u64;
    for i in (1..seq.len()).rev() {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        seq.swap(i, (state >> 33) as usize % (i + 1));
    }
    assert!(find_domains(&String::from_utf8(seq).unwrap()).is_empty());
}

fn strictly_increasing(labels: &[Label]) -> bool {
    labels.windows(2).all(|w| w[0] < w[1])
}

fn regions_never_go_back(notes: &[Annotation]) -> bool {
    notes
        .windows(2)
        .all(|w| (w[0].region as u8) <= (w[1].region as u8))
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(64))]

    #[test]
    fn labels_are_ordered_and_unique_for_any_cdr_lengths(
        h1 in 4usize..14, h2 in 4usize..16, h3 in 3usize..40,
        l1 in 4usize..13, l2 in 3usize..8, l3 in 5usize..16,
    ) {
        let h = heavy(&"G".repeat(h1), &"Y".repeat(h2), &"D".repeat(h3));
        let l = light(&"S".repeat(l1), &"A".repeat(l2), &"Q".repeat(l3));
        for seq in [h, l] {
            let d = only_domain(&seq);
            for scheme in [Scheme::Kabat, Scheme::Chothia, Scheme::Martin] {
                let n: Vec<Label> = d.numbering(scheme).into_iter().map(|(_, l)| l).collect();
                prop_assert!(strictly_increasing(&n), "{scheme:?}: {n:?}");
            }
            let imgt: Vec<Label> = d.numbering(Scheme::Imgt).into_iter().map(|(_, l)| l).collect();
            prop_assert_eq!(imgt.iter().collect::<HashSet<_>>().len(), imgt.len());
            for def in [CdrDefinition::Kabat, CdrDefinition::Chothia, CdrDefinition::Imgt, CdrDefinition::Contact, CdrDefinition::North] {
                prop_assert!(regions_never_go_back(&d.annotate(Scheme::Kabat, def)), "{def:?}");
            }
        }
    }

    #[test]
    fn arbitrary_text_never_panics(seq in r"[A-Za-z*. -]{0,400}") {
        for d in find_domains(&seq) {
            prop_assert!(d.start < d.end && d.end <= seq.len());
        }
    }

    #[test]
    fn mutated_antibodies_keep_ordered_labels(
        edits in proptest::collection::vec((0usize..240, 0usize..3, 0usize..20), 0..12),
    ) {
        let aa = b"ACDEFGHIKLMNPQRSTVWY";
        let mut seq = [trastuzumab_heavy(), CH1.to_string()].concat().into_bytes();
        for (at, kind, letter) in edits {
            let at = at.min(seq.len() - 1);
            match kind {
                0 => seq[at] = aa[letter],
                1 => seq.insert(at, aa[letter]),
                _ => { seq.remove(at); }
            }
        }
        let seq = String::from_utf8(seq).unwrap();
        for d in find_domains(&seq) {
            for scheme in [Scheme::Kabat, Scheme::Chothia, Scheme::Martin] {
                let n: Vec<Label> = d.numbering(scheme).into_iter().map(|(_, l)| l).collect();
                prop_assert!(strictly_increasing(&n), "{scheme:?}: {n:?}");
            }
        }
    }
}
