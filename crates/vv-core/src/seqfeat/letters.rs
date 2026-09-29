/// One-letter code for a residue name: upper case for amino acids (a
/// modified residue takes its parent's letter), lower case for
/// nucleotides, `None` for anything else (water, ligands, ions).
pub fn one_letter(name: &str) -> Option<char> {
    Some(match name {
        "ALA" => 'A',
        "ARG" => 'R',
        "ASN" => 'N',
        "ASP" => 'D',
        "CYS" | "CSO" | "CME" | "CSD" | "OCS" | "CAS" => 'C',
        "GLN" | "PCA" => 'Q',
        "GLU" => 'E',
        "GLY" => 'G',
        "HIS" => 'H',
        "ILE" => 'I',
        "LEU" => 'L',
        "LYS" | "MLY" | "M3L" | "KCX" | "LLP" | "ALY" => 'K',
        "MET" | "MSE" => 'M',
        "PHE" => 'F',
        "PRO" | "HYP" => 'P',
        "SER" | "SEP" => 'S',
        "THR" | "TPO" => 'T',
        "TRP" => 'W',
        "TYR" | "PTR" => 'Y',
        "VAL" => 'V',
        "SEC" => 'U',
        "PYL" => 'O',
        "A" | "DA" => 'a',
        "C" | "DC" => 'c',
        "G" | "DG" => 'g',
        "U" => 'u',
        "DT" => 't',
        "I" | "DI" => 'i',
        _ => return None,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn one_letter_codes_cover_the_standard_residues() {
        assert_eq!(one_letter("ALA"), Some('A'));
        assert_eq!(one_letter("TRP"), Some('W'));
        assert_eq!(one_letter("MSE"), Some('M'));
        assert_eq!(one_letter("SEP"), Some('S'));
        assert_eq!(one_letter("DA"), Some('a'));
        assert_eq!(one_letter("U"), Some('u'));
        assert_eq!(one_letter("HOH"), None);
        assert_eq!(one_letter("HEM"), None);
    }
}
