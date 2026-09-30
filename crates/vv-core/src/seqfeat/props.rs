//! Bulk properties of a protein chain from its sequence: molecular weight,
//! net charge and isoelectric point, and extinction coefficient at 280 nm.
//!
//! The charge model is Henderson-Hasselbalch over the ionizable groups,
//! with the pK set of Bjellqvist et al. (1993, Electrophoresis 14:1023;
//! 1994, Electrophoresis 15:529), including their residue-specific pK for
//! the N-terminal amino group and C-terminal carboxyl. The extinction
//! coefficient is the sum over Trp, Tyr and cystine of Pace et al. (1995,
//! Protein Sci 4:2411). Average residue masses are the IUPAC atomic
//! weights summed over each residue's formula.

/// Mass (Da) of the water added when residues join into a chain.
const WATER: f64 = 18.01528;

/// Average residue mass (Da); `None` for a letter without one.
fn residue_mass(c: u8) -> Option<f64> {
    Some(match c {
        b'A' => 71.0788,
        b'R' => 156.1875,
        b'N' => 114.1038,
        b'D' => 115.0886,
        b'C' => 103.1388,
        b'E' => 129.1155,
        b'Q' => 128.1307,
        b'G' => 57.0519,
        b'H' => 137.1411,
        b'I' | b'L' => 113.1594,
        b'K' => 128.1741,
        b'M' => 131.1926,
        b'F' => 147.1766,
        b'P' => 97.1167,
        b'S' => 87.0782,
        b'T' => 101.1051,
        b'W' => 186.2132,
        b'Y' => 163.1760,
        b'V' => 99.1326,
        b'U' => 150.0379,
        b'O' => 237.3018,
        _ => return None,
    })
}

/// Side-chain groups that are positive when protonated: residue, pK.
const BASIC: [(u8, f64); 3] = [(b'K', 10.0), (b'R', 12.0), (b'H', 5.98)];

/// Side-chain groups that are negative when deprotonated: residue, pK.
const ACIDIC: [(u8, f64); 4] = [(b'D', 4.05), (b'E', 4.45), (b'C', 9.0), (b'Y', 10.0)];

fn n_terminal_pk(first: u8) -> f64 {
    match first {
        b'A' => 7.59,
        b'M' => 7.0,
        b'S' => 6.93,
        b'P' => 8.36,
        b'T' => 6.82,
        b'V' => 7.44,
        b'E' => 7.7,
        _ => 7.5,
    }
}

fn c_terminal_pk(last: u8) -> f64 {
    match last {
        b'D' => 4.55,
        b'E' => 4.75,
        _ => 3.55,
    }
}

/// What [`properties`] finds for one chain.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Props {
    /// Residues with an amino-acid letter.
    pub residues: usize,
    /// Letters that are not amino acids and were left out.
    pub skipped: usize,
    /// Average mass in Da.
    pub mass: f64,
    /// Isoelectric point.
    pub pi: f64,
    /// Net charge at pH 7.
    pub charge_ph7: f64,
    /// M^-1 cm^-1 at 280 nm with every Cys reduced.
    pub extinction_reduced: f64,
    /// The same with every pair of Cys as a cystine.
    pub extinction_oxidized: f64,
}

struct Composition {
    counts: [u32; 256],
    first: u8,
    last: u8,
}

impl Composition {
    fn of(sequence: &[u8]) -> Self {
        let mut counts = [0; 256];
        for &c in sequence {
            counts[c as usize] += 1;
        }
        Composition {
            counts,
            first: sequence[0],
            last: sequence[sequence.len() - 1],
        }
    }

    fn count(&self, c: u8) -> f64 {
        f64::from(self.counts[c as usize])
    }

    fn charge_at(&self, ph: f64) -> f64 {
        let base = |pk: f64| 1.0 / (1.0 + 10f64.powf(ph - pk));
        let acid = |pk: f64| -1.0 / (1.0 + 10f64.powf(pk - ph));
        let mut charge = base(n_terminal_pk(self.first)) + acid(c_terminal_pk(self.last));
        charge += BASIC
            .iter()
            .map(|&(c, pk)| self.count(c) * base(pk))
            .sum::<f64>();
        charge += ACIDIC
            .iter()
            .map(|&(c, pk)| self.count(c) * acid(pk))
            .sum::<f64>();
        charge
    }

    /// The pH where the net charge is zero; charge falls with pH, so bisect.
    fn isoelectric_point(&self) -> f64 {
        let (mut low, mut high) = (0.0, 14.0);
        for _ in 0..50 {
            let mid = (low + high) / 2.0;
            if self.charge_at(mid) > 0.0 {
                low = mid;
            } else {
                high = mid;
            }
        }
        (low + high) / 2.0
    }
}

/// Properties of `sequence`, one upper-case letter per residue. Letters
/// with no amino-acid mass (`x`, nucleotides) are skipped; `None` if
/// nothing is left.
pub fn properties(sequence: &[u8]) -> Option<Props> {
    let kept: Vec<u8> = sequence
        .iter()
        .copied()
        .filter(|&c| residue_mass(c).is_some())
        .collect();
    if kept.is_empty() {
        return None;
    }
    let comp = Composition::of(&kept);
    let mass = kept.iter().filter_map(|&c| residue_mass(c)).sum::<f64>() + WATER;
    let aromatic = 5500.0 * comp.count(b'W') + 1490.0 * comp.count(b'Y');
    let cystines = (comp.count(b'C') / 2.0).floor();
    Some(Props {
        residues: kept.len(),
        skipped: sequence.len() - kept.len(),
        mass,
        pi: comp.isoelectric_point(),
        charge_ph7: comp.charge_at(7.0),
        extinction_reduced: aromatic,
        extinction_oxidized: aromatic + 125.0 * cystines,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Hen egg-white lysozyme, mature chain (129 residues).
    const LYSOZYME: &[u8] = b"KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL";

    #[test]
    fn a_dipeptide_has_the_mass_of_its_residues_and_a_water() {
        let p = properties(b"GG").unwrap();
        assert!((p.mass - (2.0 * 57.0519 + WATER)).abs() < 1e-9);
    }

    const HB_ALPHA: &[u8] = b"VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR";
    const HB_BETA: &[u8] = b"VHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH";
    const UBIQUITIN: &[u8] =
        b"MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG";

    /// Mass (Da) of the mature chains as UniProt lists them (P00698,
    /// P69905, P68871, ubiquitin), and pI from the Bjellqvist pK set as an
    /// independent implementation of the published method computes it.
    #[test]
    fn well_known_proteins_match_published_mass_and_pi() {
        let cases: [(&str, &[u8], f64, f64); 4] = [
            ("lysozyme", LYSOZYME, 14_313.0, 9.32),
            ("hemoglobin alpha", HB_ALPHA, 15_126.2, 8.73),
            ("hemoglobin beta", HB_BETA, 15_867.0, 6.81),
            ("ubiquitin", UBIQUITIN, 8_564.7, 6.56),
        ];
        for (name, seq, mass, pi) in cases {
            let p = properties(seq).unwrap();
            assert!((p.mass - mass).abs() < 0.5, "{name} mass {}", p.mass);
            assert!((p.pi - pi).abs() < 0.02, "{name} pI {}", p.pi);
        }
    }

    #[test]
    fn extinction_follows_trp_tyr_and_cystine_counts() {
        // Lysozyme: 6 Trp, 3 Tyr, 8 Cys (4 cystines); the measured
        // A280 of a 1 mg/mL solution is 2.63.
        let p = properties(LYSOZYME).unwrap();
        assert_eq!(p.extinction_reduced, 6.0 * 5500.0 + 3.0 * 1490.0);
        assert_eq!(p.extinction_oxidized, 37_970.0);
        assert!((p.extinction_oxidized / p.mass - 2.63).abs() < 0.05);
        // Ubiquitin's only aromatic residue is Tyr59: 1490.
        assert_eq!(properties(UBIQUITIN).unwrap().extinction_reduced, 1490.0);
    }

    #[test]
    fn charge_falls_monotonically_and_crosses_zero_at_the_pi() {
        let comp = Composition::of(LYSOZYME);
        assert!(comp.charge_at(2.0) > comp.charge_at(7.0));
        assert!(comp.charge_at(7.0) > comp.charge_at(12.0));
        assert!(comp.charge_at(comp.isoelectric_point()).abs() < 1e-6);
    }

    #[test]
    fn non_amino_letters_are_skipped() {
        let p = properties(b"xxAGxx").unwrap();
        assert_eq!((p.residues, p.skipped), (2, 4));
        assert!(properties(b"acgu").is_none());
    }
}
