//! Residue-name lists behind the `additive` and `cofactor` roles. Both
//! apply only to small molecules, so a one-atom `CLA` (chloride) never
//! reaches them while a multi-atom `CLA` (chlorophyll a) does.

/// Chemical Component Dictionary codes of the reagents that dominate the
/// non-biological ligands of crystal structures: precipitants, buffers,
/// cryoprotectants and salts. Selected from the artifact list of Yang, Roy
/// & Zhang 2013 (Nucleic Acids Res 41:D1096) and the common conditions of
/// Jancarik & Kim 1991 (J Appl Cryst 24:409) and McPherson & Gavira 2014
/// (Acta Cryst F70:2). Name-only, so a glycerol that truly binds in an
/// active site is still listed.
#[rustfmt::skip]
pub(super) const ADDITIVES: &[&str] = &[
    // Sulfate, phosphate and other inorganic anions.
    "SO4", "SO3", "PO4", "NO3", "CO3", "BCT", "SCN", "AZI", "CAC",
    // Glycerol, glycols, polyethylene glycol fragments, MPD.
    "GOL", "EDO", "PEG", "PGE", "PG4", "P6G", "1PE", "2PE", "PE4", "PE3", "PGO", "MPD", "MRD",
    "15P", "TFP",
    // Carboxylates and small alcohols.
    "ACT", "ACY", "FMT", "CIT", "TAR", "MLI", "OXL", "IPA", "EOH", "MOH", "DMS", "BME", "DTT",
    // Buffers.
    "TRS", "EPE", "MES", "MPO", "IMD", "BTB", "CXS", "HEZ",
];

/// Enzyme cofactors and prosthetic groups by their dictionary codes: heme
/// variants, flavins, nicotinamides, thiamine, pyridoxal phosphate,
/// coenzyme A and B12, S-adenosyl compounds, iron-sulfur clusters,
/// molybdopterin, biotin, lipoate, chlorophylls. After the classes of the
/// CoFactor database (Fischer et al. 2010, Nucleic Acids Res 38:D805).
pub(super) const COFACTORS: &[&str] = &[
    "HEM", "HEC", "HEA", "HEB", "HDD", "DHE", "1HE", "SRM", "FAD", "FMN", "RBF", "NAD", "NAP",
    "NDP", "NAI", "NDC", "TPP", "THD", "PLP", "PMP", "COA", "ACO", "B12", "COB", "SAM", "SAH",
    "FES", "F3S", "SF4", "CLF", "ICS", "MGD", "MTE", "BTN", "LPA", "CLA", "CL0", "BCL", "BCB",
    "PHO", "BPH", "THF", "FOL", "H4B", "UQ1", "PQN", "PL9", "MQ7",
];

pub(super) fn is_additive(name: &str) -> bool {
    ADDITIVES.iter().any(|n| n.eq_ignore_ascii_case(name))
}

pub(super) fn is_cofactor(name: &str) -> bool {
    COFACTORS.iter().any(|n| n.eq_ignore_ascii_case(name))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lists_are_disjoint_and_upper_case() {
        for name in ADDITIVES.iter().chain(COFACTORS) {
            assert_eq!(*name, name.to_ascii_uppercase());
        }
        assert!(ADDITIVES.iter().all(|n| !is_cofactor(n)));
    }
}
