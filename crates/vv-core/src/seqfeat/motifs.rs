//! Sequence motifs on one-letter strings: N-glycosylation sequons and the
//! chemical-liability motifs antibody engineers screen for (deamidation,
//! isomerization, oxidation, ...; Yang et al. 2019, mAbs 11:106;
//! Sydow et al. 2014, PLoS ONE 9:e100736).
//!
//! Scanners take a `breaks` mask (`true` where a numbering gap precedes
//! that residue) so a motif never straddles an unobserved stretch.

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Liability {
    Deamidation,
    Isomerization,
    Fragmentation,
    Oxidation,
    /// A cysteine the caller has not found in a disulfide.
    FreeCysteine,
    PyroGlutamate,
    Integrin,
}

impl Liability {
    pub const ALL: [Liability; 7] = [
        Liability::Deamidation,
        Liability::Isomerization,
        Liability::Fragmentation,
        Liability::Integrin,
        Liability::FreeCysteine,
        Liability::PyroGlutamate,
        Liability::Oxidation,
    ];

    pub fn label(self) -> &'static str {
        match self {
            Liability::Deamidation => "Deamidation",
            Liability::Isomerization => "Isomerization",
            Liability::Fragmentation => "Fragmentation",
            Liability::Oxidation => "Oxidation",
            Liability::FreeCysteine => "Free cysteine",
            Liability::PyroGlutamate => "N-terminal cyclization",
            Liability::Integrin => "Integrin binding",
        }
    }
}

/// One motif occurrence covering `start..start + len` of the sequence.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Hit {
    pub start: usize,
    pub len: usize,
    pub kind: Liability,
    pub motif: &'static str,
}

const PAIRS: [(&[u8; 2], Liability, &str); 7] = [
    (b"NG", Liability::Deamidation, "NG"),
    (b"NS", Liability::Deamidation, "NS"),
    (b"NT", Liability::Deamidation, "NT"),
    (b"DG", Liability::Isomerization, "DG"),
    (b"DS", Liability::Isomerization, "DS"),
    (b"DT", Liability::Isomerization, "DT"),
    (b"DP", Liability::Fragmentation, "DP"),
];

fn contiguous(breaks: &[bool], start: usize, len: usize) -> bool {
    breaks[start + 1..start + len].iter().all(|&b| !b)
}

/// Every liability motif in `seq`, in sequence order. A residue in two
/// motifs (`NGS`) appears in two hits. Cysteines are all reported as
/// `FreeCysteine`; the caller drops the ones that are disulfide-bonded.
pub fn liabilities(seq: &[u8], breaks: &[bool]) -> Vec<Hit> {
    let mut hits = Vec::new();
    for i in 0..seq.len() {
        single_residue_hit(seq, breaks, i, &mut hits);
        if i + 2 <= seq.len() && contiguous(breaks, i, 2) {
            pair_hit(seq, i, &mut hits);
        }
        if i + 3 <= seq.len() && contiguous(breaks, i, 3) && &seq[i..i + 3] == b"RGD" {
            hits.push(hit(i, 3, Liability::Integrin, "RGD"));
        }
    }
    hits
}

fn hit(start: usize, len: usize, kind: Liability, motif: &'static str) -> Hit {
    Hit {
        start,
        len,
        kind,
        motif,
    }
}

fn single_residue_hit(seq: &[u8], breaks: &[bool], i: usize, hits: &mut Vec<Hit>) {
    match seq[i] {
        b'M' => hits.push(hit(i, 1, Liability::Oxidation, "Met")),
        b'W' => hits.push(hit(i, 1, Liability::Oxidation, "Trp")),
        b'C' => hits.push(hit(i, 1, Liability::FreeCysteine, "Cys")),
        b'Q' if i == 0 && !breaks[0] => {
            hits.push(hit(0, 1, Liability::PyroGlutamate, "N-term Gln"))
        }
        b'E' if i == 0 && !breaks[0] => {
            hits.push(hit(0, 1, Liability::PyroGlutamate, "N-term Glu"))
        }
        _ => {}
    }
}

fn pair_hit(seq: &[u8], i: usize, hits: &mut Vec<Hit>) {
    let pair = &seq[i..i + 2];
    if let Some(&(_, kind, motif)) = PAIRS.iter().find(|(p, _, _)| p.as_slice() == pair) {
        hits.push(hit(i, 2, kind, motif));
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SequonKind {
    /// N-X-S/T, X not proline.
    Canonical,
    /// N-X-C, observed but far less often occupied.
    Rare,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Sequon {
    /// Index of the asparagine.
    pub asn: usize,
    pub kind: SequonKind,
}

/// N-linked glycosylation sequons in `seq`: N-X-S/T and N-X-C with X != P.
pub fn sequons(seq: &[u8], breaks: &[bool]) -> Vec<Sequon> {
    (0..seq.len().saturating_sub(2))
        .filter(|&i| seq[i] == b'N' && seq[i + 1] != b'P' && contiguous(breaks, i, 3))
        .filter_map(|i| {
            let kind = match seq[i + 2] {
                b'S' | b'T' => SequonKind::Canonical,
                b'C' => SequonKind::Rare,
                _ => return None,
            };
            Some(Sequon { asn: i, kind })
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn no_breaks(seq: &[u8]) -> Vec<bool> {
        vec![false; seq.len()]
    }

    fn motifs(seq: &[u8]) -> Vec<(&'static str, usize)> {
        liabilities(seq, &no_breaks(seq))
            .into_iter()
            .filter(|h| h.kind != Liability::FreeCysteine)
            .map(|h| (h.motif, h.start))
            .collect()
    }

    #[test]
    fn deamidation_and_isomerization_motifs() {
        assert_eq!(motifs(b"ANGSDGT"), [("NG", 1), ("DG", 4)]);
        assert_eq!(motifs(b"AANTAA"), [("NT", 2)]);
        assert_eq!(motifs(b"GNSD"), [("NS", 1)]);
        assert_eq!(motifs(b"ADSA"), [("DS", 1)]);
    }

    #[test]
    fn fragmentation_integrin_and_oxidation() {
        assert_eq!(motifs(b"ADPA"), [("DP", 1)]);
        assert_eq!(motifs(b"ARGDA"), [("RGD", 1)]);
        assert_eq!(motifs(b"AMAWA"), [("Met", 1), ("Trp", 3)]);
    }

    #[test]
    fn free_cysteine_and_n_terminal_cyclization() {
        let seq = b"QVCC";
        let hits = liabilities(seq, &no_breaks(seq));
        let cys: Vec<usize> = hits
            .iter()
            .filter(|h| h.kind == Liability::FreeCysteine)
            .map(|h| h.start)
            .collect();
        assert_eq!(cys, [2, 3]);
        assert!(hits
            .iter()
            .any(|h| h.kind == Liability::PyroGlutamate && h.start == 0));
        assert!(!motifs(b"AQ").iter().any(|m| m.0.starts_with("N-term")));
    }

    #[test]
    fn a_numbering_gap_splits_a_motif() {
        let seq = b"ANGA";
        let mut breaks = no_breaks(seq);
        breaks[2] = true;
        assert!(liabilities(seq, &breaks).is_empty());
        assert!(sequons(b"ANAT", &[false, false, true, false]).is_empty());
    }

    #[test]
    fn a_gap_before_the_first_residue_is_not_a_terminus() {
        let mut breaks = no_breaks(b"QA");
        breaks[0] = true;
        assert!(liabilities(b"QA", &breaks).is_empty());
    }

    #[test]
    fn sequons_need_ser_thr_and_no_proline() {
        let found = |s: &[u8]| sequons(s, &no_breaks(s));
        assert_eq!(found(b"AANASA")[0].asn, 2);
        assert_eq!(found(b"NGT")[0].kind, SequonKind::Canonical);
        assert!(found(b"NPS").is_empty(), "N-P-S is never glycosylated");
        assert!(found(b"NAA").is_empty());
        assert_eq!(found(b"NAC")[0].kind, SequonKind::Rare);
        let overlapping = found(b"NNSS");
        assert_eq!(
            overlapping.iter().map(|s| s.asn).collect::<Vec<_>>(),
            [0, 1]
        );
    }
}
