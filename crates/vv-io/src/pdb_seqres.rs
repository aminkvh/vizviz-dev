//! `SEQRES` (wwPDB v3.3): the residue names of each chain's polymer,
//! read as the file's own statement of which residues are polymer.

use std::collections::{HashMap, HashSet};

use vv_core::residue_class::{is_standard_polymer, of_name};
use vv_core::seqfeat::protein_letter;
use vv_core::{flags, PolymerHint, ResidueClass, Topology};

use crate::polymer_layout::{sequence, Layout};

/// Monomer names per chain letter, in the residue-name columns of each
/// `SEQRES` record (13 names of 3 characters from column 20, 4 apart).
#[derive(Default)]
pub(crate) struct Seqres {
    chains: HashMap<String, HashSet<String>>,
    order: HashMap<String, Vec<String>>,
}

impl Seqres {
    pub(crate) fn record(&mut self, line: &[u8]) {
        let Some(&chain) = line.get(11) else { return };
        let chain = (chain as char).to_string();
        let names = self.chains.entry(chain.clone()).or_default();
        let order = self.order.entry(chain).or_default();
        let mut col = 19;
        while col + 3 <= line.len() {
            let name = String::from_utf8_lossy(&line[col..col + 3]);
            let name = name.trim();
            if !name.is_empty() {
                names.insert(name.to_string());
                order.push(name.to_string());
            }
            col += 4;
        }
    }

    /// The one-letter sequence of each chain record whose `SEQRES` is a
    /// protein, `""` for the others; empty when there is none.
    pub(crate) fn sequences(&self, topology: &Topology) -> Vec<String> {
        let sequences: Vec<String> = topology
            .chains
            .iter()
            .map(|chain| {
                let id = topology.names.get(chain.label_asym);
                match self.chains.get(id).map(chain_kind) {
                    Some(PolymerHint::Protein) => {
                        self.order[id].iter().map(|n| protein_letter(n)).collect()
                    }
                    _ => String::new(),
                }
            })
            .collect();
        match sequences.iter().any(|s| !s.is_empty()) {
            true => sequences,
            false => Vec::new(),
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
                let comp = topology.names.get(res.comp);
                let hetero = is_hetero(topology, res.atoms.start as usize);
                let Some(names) = self.chains.get(chain) else {
                    return if hetero && is_standard_polymer(comp) {
                        PolymerHint::NonPolymer
                    } else {
                        PolymerHint::Unknown
                    };
                };
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

/// `SEQRES` records for the polymer segments of `layout`: per PDB chain
/// character, the residue names in order, 13 to a line (columns: serial
/// 8-10, chain 12, count 14-17, names from 20).
pub(crate) fn seqres_lines(topology: &Topology, layout: &Layout, chain_ids: &[u8]) -> Vec<String> {
    let mut chains: Vec<(u8, Vec<&str>)> = Vec::new();
    for segment in layout.segments.iter().filter(|s| s.kind.is_polymer()) {
        let byte = chain_ids[segment.chain as usize];
        let at = chains
            .iter()
            .position(|(b, _)| *b == byte)
            .unwrap_or_else(|| {
                chains.push((byte, Vec::new()));
                chains.len() - 1
            });
        for comp in sequence(topology, segment) {
            let name = topology.names.get(comp);
            chains[at].1.push(name.get(..3).unwrap_or(name));
        }
    }
    chains
        .iter()
        .flat_map(|(byte, names)| {
            names.chunks(13).enumerate().map(move |(i, chunk)| {
                let listed: Vec<String> = chunk.iter().map(|n| format!("{n:<3}")).collect();
                format!(
                    "SEQRES {:>3} {} {:>4}  {}",
                    i + 1,
                    *byte as char,
                    names.len(),
                    listed.join(" ")
                )
            })
        })
        .collect()
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
