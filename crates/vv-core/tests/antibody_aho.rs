//! AHo numbering (Honegger & Plückthun 2001) against the published alignment
//! of Figure 1(b) and, when ixtures/real/holdout exists, against real
//! wwPDB domains.

use vv_core::antibody::{find_domains, Domain, Scheme};

mod antibody_common;
use antibody_common::{fixtures, read_chains};

const FIRST_COLUMN: usize = 3;

fn aho(d: &Domain) -> Vec<(usize, u16, bool)> {
    d.numbering(Scheme::Aho)
        .into_iter()
        .map(|(i, l)| (i, l.number, l.insertion().is_some()))
        .collect()
}

fn occupancy(d: &Domain) -> String {
    let mut cols = vec!['.'; 144];
    for (_, number, _) in aho(d) {
        if let Some(slot) = usize::from(number).checked_sub(FIRST_COLUMN) {
            if slot < cols.len() {
                cols[slot] = '#';
            }
        }
    }
    cols.into_iter().collect()
}

#[test]
fn aho_reproduces_the_published_alignment_of_all_fourteen_domains() {
    for (name, seq, figure) in FIGURE_1B {
        let found = find_domains(seq);
        assert_eq!(found.len(), 1, "{name}");
        assert_eq!(occupancy(&found[0]), figure, "{name}");
    }
}

#[test]
fn published_anchor_residues_sit_at_their_aho_positions() {
    for (name, seq, _) in FIGURE_1B {
        let d = &find_domains(seq)[0];
        let at = |n: u16| {
            aho(d)
                .into_iter()
                .find(|&(_, number, ins)| number == n && !ins)
                .map(|(i, _, _)| seq.as_bytes()[i] as char)
        };
        assert_eq!(at(23), Some('C'), "{name}");
        assert_eq!(at(106), Some('C'), "{name}");
        assert!(matches!(at(43), Some('W')), "{name}");
        assert!(matches!(at(139), Some('W' | 'F')), "{name}");
        assert_eq!(at(140), Some('G'), "{name}");
    }
}

#[test]
fn aho_numbers_rise_strictly_within_one_to_one_forty_nine() {
    for (name, seq, _) in FIGURE_1B {
        let labels = find_domains(seq)[0].numbering(Scheme::Aho);
        assert!(labels.windows(2).all(|w| w[0].1 < w[1].1), "{name}");
        assert!(labels.iter().all(|(_, l)| (1..=149).contains(&l.number)));
    }
}

#[derive(Default)]
struct Anchors {
    domains: usize,
    in_range: usize,
    rising: usize,
    cys23: usize,
    trp43: usize,
    cys106: usize,
    fr4_139: usize,
    gly140: usize,
    inserted: usize,
    /// Domains whose FR4 start is not `[WF]G` at 139-140.
    odd_fr4: Vec<String>,
}

fn residue_at(d: &Domain, seq: &[u8], n: u16) -> Option<u8> {
    aho(d)
        .into_iter()
        .find(|&(_, number, ins)| number == n && !ins)
        .map(|(i, _, _)| seq[i])
}

fn tally(a: &mut Anchors, id: &str, seq: &str, d: &Domain) {
    let labels = d.numbering(Scheme::Aho);
    let bytes = seq.as_bytes();
    a.domains += 1;
    a.in_range += usize::from(labels.iter().all(|(_, l)| (1..=149).contains(&l.number)));
    a.rising += usize::from(labels.windows(2).all(|w| w[0].1 < w[1].1));
    a.inserted += usize::from(labels.iter().any(|(_, l)| l.insertion().is_some()));
    a.cys23 += usize::from(residue_at(d, bytes, 23) == Some(b'C'));
    a.trp43 += usize::from(residue_at(d, bytes, 43) == Some(b'W'));
    a.cys106 += usize::from(residue_at(d, bytes, 106) == Some(b'C'));
    a.fr4_139 += usize::from(matches!(residue_at(d, bytes, 139), Some(b'W' | b'F')));
    a.gly140 += usize::from(residue_at(d, bytes, 140) == Some(b'G'));
    let fr4: String = [139, 140]
        .iter()
        .map(|&n| residue_at(d, bytes, n).map_or('-', char::from))
        .collect();
    if !matches!(fr4.as_bytes(), [b'W' | b'F', b'G']) {
        a.odd_fr4.push(format!("{id} {fr4}"));
    }
}

fn held_out_anchors() -> Anchors {
    let mut a = Anchors::default();
    for (id, path) in fixtures("holdout") {
        for chain in read_chains(&path) {
            let seq = chain.sequence();
            for d in find_domains(&seq) {
                tally(&mut a, &format!("{id}:{}", chain.id), &seq, &d);
            }
        }
    }
    a
}

#[test]
fn held_out_domains_are_numbered_within_one_to_one_forty_nine_on_the_conserved_cysteines() {
    let a = held_out_anchors();
    if a.domains == 0 {
        return;
    }
    assert_eq!(a.in_range, a.domains);
    assert_eq!(a.rising, a.domains);
    assert_eq!(a.cys23, a.domains);
    assert_eq!(a.cys106, a.domains);
    assert_eq!(a.trp43, a.domains);
}

#[test]
#[ignore = "prints the AHo anchor report; needs fixtures/real/holdout"]
fn aho_anchor_report() {
    let a = held_out_anchors();
    println!(
        "{} domains: 1..149 {} rising {} C23 {} W43 {} C106 {} W/F139 {} G140 {} with insertion letters {}",
        a.domains, a.in_range, a.rising, a.cys23, a.trp43, a.cys106, a.fr4_139, a.gly140, a.inserted
    );
    println!("FR4 exceptions: {:?}", a.odd_fr4);
}

/// Empty columns of the AHo stretch `lo..=hi` of a domain.
fn gap_columns(d: &Domain, lo: u16, hi: u16) -> Vec<u16> {
    let used: Vec<u16> = aho(d).into_iter().map(|(_, n, _)| n).collect();
    (lo..=hi).filter(|c| !used.contains(c)).collect()
}

#[test]
#[ignore = "prints the AHo gap statistics; needs fixtures/real/holdout"]
fn aho_gap_report() {
    use std::collections::BTreeMap;
    let mut cdr2: BTreeMap<(bool, usize), usize> = BTreeMap::new();
    let mut vl_ok = (0, 0);
    for (_, path) in fixtures("holdout") {
        for chain in read_chains(&path) {
            for d in find_domains(&chain.sequence()) {
                let heavy = d.chain == vv_core::antibody::ChainType::Heavy;
                let gaps = gap_columns(&d, 58, 68);
                *cdr2.entry((heavy, gaps.len())).or_default() += 1;
                if !heavy && gaps.len() == 8 {
                    vl_ok.1 += usize::from(gaps == (59..=66).collect::<Vec<_>>());
                    vl_ok.0 += 1;
                }
            }
        }
    }
    for ((heavy, gap), n) in cdr2 {
        println!(
            "CDR2 gap of {gap} columns, {}: {n}",
            if heavy { "heavy" } else { "light" }
        );
    }
    println!(
        "light chains with an 8-column CDR2 gap that is exactly 59-66: {}/{}",
        vl_ok.1, vl_ok.0
    );
}

/// `(entry chain, sequence, occupied AHo columns 3..=146)` for the 14 immunoglobulin
/// domains of Honegger & Plückthun (2001) Figure 1(b); `#` is a residue, `.` a gap.
const FIGURE_1B: [(&str, &str, &str); 14] = [
    (
        "1MFA lambda",
        "QIVVTQESALTTSPGETVTLTCRSSTGTVTSGNHANWVQEKPDHLFTGLIGDTNNRAPGVPARFSGSLIGDKAALTITGAQPEDEAIYFCALWSNNHWIFGGGTKLTVLGQPKSSPSVTLFPPSSEG",
        "#####.###################.#####....#####################........##################..#########################.......................############",
    ),
    (
        "2FB4 lambda",
        "QSVLTQPPSASGTPGQRVTISCSGTSSNIGSSTVNWYQQLPGMAPKLLIYRDAMRPSGVPDRFSGSKSGASASLAIGGLQSEDETDYYCAAWDVSLNAYVFGTGTKVTVLGQPKANPTVTLFPPSSEELQANKATLVCLISDFYPGAVTVAWKADGSPVKAGVETTKPSKQSNNKYAASSYLSLTPEQWKSHRSYSCQVTHEGSTVEKTVAPTECS",
        "#####.###################.#####.....####################........##################..##########################.....................#############",
    ),
    (
        "8FAB lambda",
        "SYELTQPPSVSVSPGQTARITCSANALPNQYAYWYQQKPGRAPVMVIYKDTQRPSGIPQRFSSSTSGTTVTLTISGVQAEDEADYYCQAWDNSASIFGGGTKLTVLGQPKAAPSVTLFPPSSEELQANKATLVCLISDFYPGAVTVAWKADSSPIKAGVETTTPSKQSNNKYAASSYLSLTPEQWKSHRSYSCQVTHEGSTVEKTVAPTECS",
        "#####.##################...####.....####################........##################..#########################.......................############",
    ),
    (
        "2FBJ kappa",
        "EIVLTQSPAITAASLGQKVTITCSASSSVSSLHWYQQKSGTSPKPWIYEISKLASGVPARFSGSGSGTSYSLTINTMEAEDAAIYYCQQWTYPLITFGAGTKLELKRADAAPTVSIFPPSSEQLTSGGASVVCFLNNFYPKDINVKWKIDGSERQNGVLNSWTDQDSKDSTYSMSSTLTLTKDEYERHNSYTCEATHKTSTSPIVKSFNRNEC",
        "########################..####.......###################........##################..#########################.......................############",
    ),
    (
        "1A2Y kappa",
        "DIVLTQSPASLSASVGETVTITCRASGNIHNYLAWYQQKQGKSPQLLVYYTTTLADGVPSRFSGSGSGTQYSLKINSLQPEDFGSYYCQHFWSTPRTFGGGTKLEIK",
        "########################..####......####################........##################..#########################.......................############",
    ),
    (
        "25C8 kappa",
        "DIVLTQSPAIMSASLGERVTMTCTASSSVSSSNLHWYQQKPGSSPKLWIYSTSNLASGVPARFSGSGSGTSYSLTISSMEAEDAATYYCHQYHRSPYTFGGGTKLEIKRADAAPTVSIFPPSSEQLTSGGASVVCFLNNFYPKDINVKWKIDGSERQNGVLNSWTDQDSKDSTYSMSSTLTLTKDEYERHNSYTCEATHKTSTSPIVKSFNR",
        "########################..#####.....####################........##################..#########################.......................############",
    ),
    (
        "1F58 kappa",
        "DIVLTQSPASLAVSLGQRATISCKASQGVDFDGASFMNWYQQKPGQPPKLLIFAASTLESGIPARFSGRGSGTDFTLNIHPVEEEDAATYYCQQSHEDPLTFGAGTKLELKRADAAPTVSIFPPSSEQLTSGGASVVCFLNNFYPKDINVKWKIDGSERQNGVLNSWTDQDSKDSTYSMSSTLTLTKDEYERHNSYTCEATHKTSTSPIVKSFNRA",
        "########################..######..######################........##################..#########################.......................############",
    ),
    (
        "1FLR kappa",
        "DVVMTQTPLSLPVSLGDQASISCRSSQSLVHSNGNTYLRWYLQKPGQSPKVLIYKVSNRFSGVPDRFSGSGSGTDFTLKISRVEAEDLGVYFCSQSTHVPWTFGGGTKLEIKRADAAPTVSIFPPSSEQLTSGGASVVCFLNNFYPKDINVKWKIDGSERQNGVLNSWTDQDSKDSTYSMSSTLTLTKDEYERHNSYTCEATHKTSTSPIVKSFNRNEC",
        "########################..#######.######################........##################..#########################.......................############",
    ),
    (
        "1HIL kappa",
        "DIVMTQSPSSLTVTAGEKVTMSCTSSQSLFNSGKQKNYLTWYQQKPGQPPKVLIYWASTRESGVPDRFTGSGSGTDFTLTISSVQAEDLAVYYCQNDYSNPLTFGGGTKLELKRADAAPTVSIFPPSSEQLTSGGASVVCFLNNFYPKDINVKWKIDGSERQNGVLNSWTDQDSKDSTYSMSSTLTLTKDEYERHNSYTCEATHKTSTSPIVKSFNR",
        "########################..##############################........##################..#########################.......................############",
    ),
    (
        "2HMI heavy",
        "QITLKESGPGIVQPSQPFRLTCTFSGFSLSTSGIGVTWIRQPSGKGLEWLATIWWDDDNRYNPSLKSRLTVSKDTSNNQAFLNMMTVETADTAIYYCAQSAITSVTDSAMDHWGQGTSVTVSSAATTPPSVYPLAPGSAAQTNSMVTLGCLVKGYFPEPVTVTWNSGSLSSGVHTFPAVLQSDLYTLSSSVTVPSSTWPSETVTCNVAHPASSTKVDKKI",
        "#####.###################.######...#######################....##################################################.................###############",
    ),
    (
        "1F58 heavy",
        "DVQLQQSGPDLVKPSQSLSLTCTVTGYSITSGYSWHWIRQFPGNKLEWMGYIHYSAGTNYNPSLKSRISITRDTSKNQFFLQLNSVTTEDTATYYCAREEAMPYGNQAYYYAMDCWGQGTTVTVSSAKTTPPSVYPLAPGSAAQTNSMVTLGCLVKGYFPEPVTVTWNSGSLSSGVHTFPAVLQSDLYTLSSSVTVPSSPRPSETVTCNVAHPASSTKVDKKIVPRDC",
        "#####.###################.#####....#######################....####################################################.............#################",
    ),
    (
        "1A2Y heavy",
        "QVQLQESGPGLVAPSQSLSITCTVSGFSLTGYGVNWVRQPPGKGLEWLGMIWGDGNTDYNSALKSRLSISKDNSKSQVFLKMNSLHTDDTARYYCARERDYRLDYWGQGTTLTVSS",
        "#####.###################.#####.....######################....################################################......................############",
    ),
    (
        "1MFA heavy",
        "EVQVQQSGTVVARPGASVKMSCKASGYTFTNYWMHWIKQRPGQGLEWIGAIYPGNSATFYNHKFRAKTKLTAVTSTTTAYMELSSLTSEDSAVYYCTRGGHGYYGDYWGQGASLTVSSAK",
        "#####.###################.#####.....#######################...################################################.....................#############",
    ),
    (
        "1A6V heavy",
        "QVQLQQPGAELVKPGASVKLSCKASGYTFTSYWMHWVKQRPGRGLEWIGRIDPNSGGTKYNEKFKSKATLTVDKPSSTAYMQLSSLTSEDSAVYYCARYDYYGSSYFDYWGQGTTVTVSS",
        "#####.###################.#####.....#######################...#################################################...................##############",
    ),
];
