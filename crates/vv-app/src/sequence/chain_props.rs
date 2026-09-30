//! Bulk properties of a chain's modeled sequence (`vv_core::seqfeat::props`)
//! as the chain tooltip and `sequence props` show them.

use std::ops::Range;

use vv_core::seqfeat::{properties, Props};
use vv_core::Topology;

use super::rows::letters_and_breaks;

/// Properties of the modeled residues of a protein chain, `None` for a
/// nucleic acid or a chain with no amino-acid residue.
pub fn of_chain(top: &Topology, residues: Range<u32>) -> Option<Props> {
    let (letters, _) = letters_and_breaks(top, residues);
    properties(&letters)
}

fn thousands(value: f64) -> String {
    let digits = format!("{:.0}", value);
    let mut out = String::new();
    for (i, c) in digits.chars().enumerate() {
        if i > 0 && (digits.len() - i) % 3 == 0 {
            out.push(',');
        }
        out.push(c);
    }
    out
}

/// `\u{3b5}280`: one number when Cys pairing changes nothing, else both.
fn extinction(p: &Props) -> String {
    if p.extinction_oxidized == p.extinction_reduced {
        format!(
            "{} M\u{207B}\u{B9}cm\u{207B}\u{B9}",
            thousands(p.extinction_reduced)
        )
    } else {
        format!(
            "{} with cystines, {} reduced M\u{207B}\u{B9}cm\u{207B}\u{B9}",
            thousands(p.extinction_oxidized),
            thousands(p.extinction_reduced)
        )
    }
}

/// The lines of the chain tooltip.
pub fn lines(p: &Props) -> Vec<String> {
    let mut lines = vec![
        format!(
            "{} residues modeled, {:.2} kDa",
            p.residues,
            p.mass / 1000.0
        ),
        format!("pI {:.2}, charge {:+.1} at pH 7", p.pi, p.charge_ph7),
        format!("\u{3B5}280 {}", extinction(p)),
    ];
    if p.skipped > 0 {
        lines[0].push_str(&format!(" ({} unknown left out)", p.skipped));
    }
    lines
}

/// One line for `sequence props`.
pub fn one_line(p: &Props) -> String {
    lines(p).join("; ")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn thousands_separators() {
        assert_eq!(thousands(37970.0), "37,970");
        assert_eq!(thousands(970.4), "970");
        assert_eq!(thousands(1_000_000.0), "1,000,000");
    }

    #[test]
    fn a_dipeptide_reads_out_mass_charge_and_extinction() {
        let p = properties(b"WY").unwrap();
        let text = one_line(&p);
        assert!(text.contains("2 residues modeled"), "{text}");
        assert!(text.contains("\u{3B5}280 6,990"), "{text}");
    }
}
