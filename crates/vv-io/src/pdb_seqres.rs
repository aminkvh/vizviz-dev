//! `SEQRES` (wwPDB v3.3): the residue names of each chain's polymer,
//! read as the file's own statement of which residues are polymer.

use std::collections::{HashMap, HashSet};

use vv_core::residue_class::{is_standard_polymer, of_name};
use vv_core::{flags, PolymerHint, ResidueClass, Topology};

/// Monomer names per chain letter, in the residue-name columns of each
/// `SEQRES` record (13 names of 3 characters from column 20, 4 apart).
#[derive(Default)]
pub(crate) struct Seqres {
    chains: HashMap<String, HashSet<String>>,
}

impl Seqres {
    pub(crate) fn record(&mut self, line: &[u8]) {
        let Some(&chain) = line.get(11) else { return };
        let names = self.chains.entry((chain as char).to_string()).or_default();
        let mut col = 19;
        while col + 3 <= line.len() {
            let name = String::from_utf8_lossy(&line[col..col + 3]);
            let name = name.trim();
            if !name.is_empty() {
                names.insert(name.to_string());
            }
            col += 4;
        }
    }

    /// One hint per residue, or empty when the file has no `SEQRES`.
    ///
    /// A residue named in its chain's `SEQRES` is polymer, protein or
    /// nucleic by the majority of that chain's names. A `HETATM` residue is
    /// not, unless it is a modified residue the `SEQRES` names: a standard
    /// residue is written `ATOM` inside a polymer, so a `HETATM` one is a
    /// free ligand.
    pub(crate) fn hints(&self, topology: &Topology) -> Vec<PolymerHint> {
        if self.chains.is_empty() {
            return Vec::new();
        }
        let kinds: HashMap<&str, PolymerHint> = self
            .chains
            .iter()
            .map(|(chain, names)| (chain.as_str(), chain_kind(names)))
            .collect();
        topology
            .residues
            .iter()
            .map(|res| {
                let chain = topology
                    .names
                    .get(topology.chains[res.chain as usize].label_asym);
                let Some(names) = self.chains.get(chain) else {
                    return PolymerHint::Unknown;
                };
                let comp = topology.names.get(res.comp);
                let hetero = is_hetero(topology, res.atoms.start as usize);
                match (names.contains(comp), hetero) {
                    (true, true) if is_standard_polymer(comp) => PolymerHint::NonPolymer,
                    (true, _) => kinds[chain],
                    (false, true) => PolymerHint::NonPolymer,
                    (false, false) => PolymerHint::Unknown,
                }
            })
            .collect()
    }
}

fn is_hetero(topology: &Topology, atom: usize) -> bool {
    topology.flags[atom] & flags::HETERO != 0
}

fn chain_kind(names: &HashSet<String>) -> PolymerHint {
    let count = |class| {
        names
            .iter()
            .filter(|n| of_name(n, usize::MAX) == Some(class))
            .count()
    };
    let (protein, nucleic) = (count(ResidueClass::Protein), count(ResidueClass::Nucleic));
    match (protein, nucleic) {
        (0, 0) => PolymerHint::Unknown,
        (p, n) if p >= n => PolymerHint::Protein,
        _ => PolymerHint::Nucleic,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn record_reads_thirteen_names_per_line() {
        let mut s = Seqres::default();
        s.record(b"SEQRES   1 A   21  GLY ILE VAL GLU GLN CYS CYS ALA SER VAL CYS SER LEU");
        let names = &s.chains["A"];
        assert!(names.contains("GLY") && names.contains("LEU"));
        assert_eq!(names.len(), 9);
    }

    #[test]
    fn majority_of_names_picks_the_polymer_kind() {
        let names = |v: &[&str]| v.iter().map(|s| s.to_string()).collect::<HashSet<_>>();
        assert_eq!(
            chain_kind(&names(&["ALA", "MSE", "GLY"])),
            PolymerHint::Protein
        );
        assert_eq!(
            chain_kind(&names(&["DA", "DT", "DG"])),
            PolymerHint::Nucleic
        );
        assert_eq!(chain_kind(&names(&["XYZ"])), PolymerHint::Unknown);
    }
}
