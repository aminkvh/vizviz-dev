//! Static per-atom attributes and the chain -> residue -> atom hierarchy.
//!
//! Atoms are stored chain-major and residue-contiguous, so the hierarchy is
//! just index ranges: a residue owns `atoms`, a chain owns `residues`.

use std::ops::Range;

use glam::Vec3;

use crate::antibody::AntibodyCache;
use crate::residue_class::{self, ClassCounts, PolymerHint, ResidueClass, Roles};
use crate::{Annotations, BondOrder, BondTable, Element, InternId, Interner};

#[derive(Clone, Copy, PartialEq, Eq, Debug, Default)]
#[repr(u8)]
pub enum SecondaryStructure {
    #[default]
    Unknown = 0,
    Coil = 1,
    Helix = 2,
    Strand = 3,
}

/// Per-atom bit flags (`Topology::flags`).
pub mod flags {
    /// HETATM record (ligand, water, ion) rather than polymer ATOM.
    pub const HETERO: u8 = 1 << 0;
    /// Deuterium: the element is hydrogen (every hydrogen rule applies);
    /// the flag only keeps the isotope so a writer can emit `D` again.
    pub const DEUTERIUM: u8 = 1 << 1;
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ResidueRec {
    pub atoms: Range<u32>,
    pub chain: u32,
    /// Residue type name (e.g. "ALA"), interned in `Topology::names`.
    pub comp: InternId,
    pub seq_id: i32,
    pub auth_seq_id: i32,
    pub ins_code: u8,
    pub ss: SecondaryStructure,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ChainRec {
    pub residues: Range<u32>,
    pub label_asym: InternId,
    pub auth_asym: InternId,
    pub entity: u16,
}

/// A bond given explicitly by the file (mmCIF `_struct_conn`, PDB CONECT).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ExplicitBond {
    pub atoms: [u32; 2],
    pub kind: ExplicitBondKind,
    /// Non-`Single` only for a legacy PDB CONECT pair repeated per the old
    /// double/triple-bond convention (`vv_io::pdb`); `struct_conn`/LINK
    /// bonds are always `Single` here (order comes from `_chem_comp_bond`
    /// or the residue templates instead, see `bonds` module doc).
    pub order: BondOrder,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ExplicitBondKind {
    Covalent,
    Disulfide,
    Metal,
    Other,
}

/// Everything about a structure except coordinates. Columnar per atom.
#[derive(Clone, Debug, Default)]
pub struct Topology {
    pub element: Vec<Element>,
    /// Atom name, space padded (" CA ", " OXT").
    pub name: Vec<[u8; 4]>,
    pub serial: Vec<u32>,
    pub b_factor: Vec<f32>,
    pub occupancy: Vec<f32>,
    /// Alternate location id, `0` when none.
    pub alt_loc: Vec<u8>,
    pub charge: Vec<i8>,
    pub flags: Vec<u8>,
    pub residue_index: Vec<u32>,
    pub residues: Vec<ResidueRec>,
    /// One class per residue, or empty until `assign_residue_classes`
    /// runs (`Structure` construction does).
    pub residue_class: Vec<ResidueClass>,
    /// One role set per residue, assigned with `residue_class`.
    pub residue_roles: Vec<Roles>,
    /// Residues per class; all zero while `residue_class` is empty.
    pub class_counts: ClassCounts,
    pub chains: Vec<ChainRec>,
    /// Segment id (PDB columns 73-76) of each chain record; empty when no
    /// producer set any, `""` for a record without one.
    pub segids: Vec<InternId>,
    /// Atoms whose name exceeds the 4-byte `name` column (mmCIF only),
    /// ascending by atom index; `name` then holds the first four bytes.
    pub long_names: Vec<(u32, Box<str>)>,
    /// What the file itself says about each residue's polymer status
    /// (mmCIF entity tables, PDB `SEQRES`); empty when it says nothing.
    /// Overrides the name tables in `assign_residue_classes`.
    pub polymer_hint: Vec<PolymerHint>,
    pub names: Interner,
    pub explicit_bonds: Vec<ExplicitBond>,
    /// Non-`Single` bond orders the file itself named by atom (mmCIF
    /// `_chem_comp_bond.value_order`), atom indices ordered `a < b`:
    /// highest priority of the three order sources (see `bonds` module
    /// doc), resolved onto whatever pair perception/templates produce for
    /// those two atoms.
    pub chem_comp_bond_order: Vec<([u32; 2], BondOrder)>,
    /// Bonds from perception merged with explicit ones; `None` until computed.
    pub bonds: Option<BondTable>,
    /// Bonds an MD topology file (PSF, PRMTOP) supplied directly. When
    /// set, `bonds::perceive` returns these verbatim instead of running
    /// perception (see `bonds` module doc).
    pub md_bonds: Option<Vec<[u32; 2]>>,
    pub title: String,
    pub id: String,
    /// Header metadata the file carried (method, resolution, citation, ...).
    pub annotations: Annotations,
    /// Antibody domains per chain, found on first query.
    pub antibody: AntibodyCache,
}

#[derive(thiserror::Error, Debug, PartialEq, Eq)]
pub enum TopologyError {
    #[error("per-atom column `{column}` has {actual} entries, expected {expected}")]
    ColumnLength {
        column: &'static str,
        expected: usize,
        actual: usize,
    },
    #[error(
        "residue {residue} atom range {start}..{end} is not contiguous with the previous residue"
    )]
    ResidueRange {
        residue: usize,
        start: u32,
        end: u32,
    },
    #[error("residues cover {covered} atoms but there are {atoms}")]
    ResidueCoverage { covered: u32, atoms: usize },
    #[error(
        "chain {chain} residue range {start}..{end} is not contiguous with the previous chain"
    )]
    ChainRange { chain: usize, start: u32, end: u32 },
    #[error("chains cover {covered} residues but there are {residues}")]
    ChainCoverage { covered: u32, residues: usize },
    #[error("atom {atom} has residue_index {actual}, expected {expected}")]
    ResidueIndex {
        atom: u32,
        expected: u32,
        actual: u32,
    },
    #[error("residue {residue} points at chain {actual}, expected {expected}")]
    ResidueChain {
        residue: usize,
        expected: u32,
        actual: u32,
    },
    #[error("bond references atom {atom} but there are {atoms}")]
    BondAtom { atom: u32, atoms: usize },
    #[error("long name given for atom {atom} but there are {atoms}")]
    LongNameAtom { atom: u32, atoms: usize },
}

impl Topology {
    pub fn atom_count(&self) -> usize {
        self.element.len()
    }

    pub fn residue_count(&self) -> usize {
        self.residues.len()
    }

    pub fn chain_count(&self) -> usize {
        self.chains.len()
    }

    /// Trimmed atom name, or `""` when the structure has no name column.
    pub fn atom_name(&self, atom: usize) -> &str {
        if let Ok(i) = self
            .long_names
            .binary_search_by_key(&(atom as u32), |(a, _)| *a)
        {
            return &self.long_names[i].1;
        }
        self.name
            .get(atom)
            .and_then(|n| std::str::from_utf8(n).ok())
            .unwrap_or("")
            .trim()
    }

    /// The chain record's segment id, `""` when it has none.
    pub fn segid(&self, chain: usize) -> &str {
        self.segids.get(chain).map_or("", |&id| self.names.get(id))
    }

    pub fn is_deuterium(&self, atom: usize) -> bool {
        self.flags
            .get(atom)
            .is_some_and(|f| f & flags::DEUTERIUM != 0)
    }

    pub fn residue_name(&self, residue: usize) -> &str {
        self.names.get(self.residues[residue].comp)
    }

    pub fn chain_name(&self, chain: usize) -> &str {
        self.names.get(self.chains[chain].label_asym)
    }

    /// The residue's molecule class: the cached one once classified,
    /// otherwise decided from its name and elements alone.
    pub fn residue_class(&self, residue: usize) -> ResidueClass {
        match self.residue_class.get(residue) {
            Some(&class) => class,
            None => residue_class::classify_alone(self, residue),
        }
    }

    /// The residue's roles: the cached ones once assigned, otherwise
    /// decided from its class and name alone.
    pub fn residue_roles(&self, residue: usize) -> Roles {
        match self.residue_roles.get(residue) {
            Some(&roles) => roles,
            None => {
                residue_class::roles::alone(self.residue_class(residue), self.residue_name(residue))
            }
        }
    }

    /// Classifies every residue, assigns its roles, and caches both with
    /// the per-class counts. `positions` are frame-0 coordinates, used to
    /// perceive the bonds of residues whose names no table knows and to
    /// find membranes and attached glycans.
    pub fn assign_residue_classes(&mut self, positions: &[Vec3]) {
        let classes = residue_class::classify(self, Some(positions));
        self.residue_roles = residue_class::roles::assign(self, &classes, positions);
        self.class_counts = residue_class::counts(&classes);
        self.residue_class = classes;
    }

    /// Per chain record, the index of its chain name (`label_asym`) in
    /// first-seen order. Records that share a name (a PDB `TER` splits one
    /// chain letter into polymer and water records) share an index, so
    /// per-chain colours agree.
    pub fn chain_name_index(&self) -> Vec<u32> {
        let mut seen: Vec<InternId> = Vec::new();
        self.chains
            .iter()
            .map(|c| {
                seen.iter()
                    .position(|&s| s == c.label_asym)
                    .unwrap_or_else(|| {
                        seen.push(c.label_asym);
                        seen.len() - 1
                    }) as u32
            })
            .collect()
    }

    /// Chain index of an atom.
    pub fn chain_of_atom(&self, atom: usize) -> u32 {
        self.residues[self.residue_index[atom] as usize].chain
    }

    /// Checks the structural invariants every producer (parser, generator)
    /// must uphold. O(atoms); meant for tests and debug builds.
    pub fn validate(&self) -> Result<(), TopologyError> {
        let atoms = self.atom_count();
        let column = |column: &'static str, len: usize| {
            if len == atoms {
                Ok(())
            } else {
                Err(TopologyError::ColumnLength {
                    column,
                    expected: atoms,
                    actual: len,
                })
            }
        };
        column("residue_index", self.residue_index.len())?;
        if !self.residue_class.is_empty() && self.residue_class.len() != self.residues.len() {
            return Err(TopologyError::ColumnLength {
                column: "residue_class",
                expected: self.residues.len(),
                actual: self.residue_class.len(),
            });
        }
        for (column, len) in [
            ("residue_roles", self.residue_roles.len()),
            ("polymer_hint", self.polymer_hint.len()),
        ] {
            if len != 0 && len != self.residues.len() {
                return Err(TopologyError::ColumnLength {
                    column,
                    expected: self.residues.len(),
                    actual: len,
                });
            }
        }
        if !self.segids.is_empty() && self.segids.len() != self.chains.len() {
            return Err(TopologyError::ColumnLength {
                column: "segids",
                expected: self.chains.len(),
                actual: self.segids.len(),
            });
        }
        if let Some((a, _)) = self.long_names.iter().find(|(a, _)| *a as usize >= atoms) {
            return Err(TopologyError::LongNameAtom { atom: *a, atoms });
        }
        for (name, len) in [
            ("name", self.name.len()),
            ("serial", self.serial.len()),
            ("b_factor", self.b_factor.len()),
            ("occupancy", self.occupancy.len()),
            ("alt_loc", self.alt_loc.len()),
            ("charge", self.charge.len()),
            ("flags", self.flags.len()),
        ] {
            if len != 0 {
                column(name, len)?;
            }
        }

        let mut next_atom = 0u32;
        for (r, res) in self.residues.iter().enumerate() {
            if res.atoms.start != next_atom || res.atoms.end < res.atoms.start {
                return Err(TopologyError::ResidueRange {
                    residue: r,
                    start: res.atoms.start,
                    end: res.atoms.end,
                });
            }
            for a in res.atoms.clone() {
                let actual = self.residue_index[a as usize];
                if actual != r as u32 {
                    return Err(TopologyError::ResidueIndex {
                        atom: a,
                        expected: r as u32,
                        actual,
                    });
                }
            }
            next_atom = res.atoms.end;
        }
        if next_atom as usize != atoms {
            return Err(TopologyError::ResidueCoverage {
                covered: next_atom,
                atoms,
            });
        }

        let mut next_res = 0u32;
        for (c, chain) in self.chains.iter().enumerate() {
            if chain.residues.start != next_res || chain.residues.end < chain.residues.start {
                return Err(TopologyError::ChainRange {
                    chain: c,
                    start: chain.residues.start,
                    end: chain.residues.end,
                });
            }
            for r in chain.residues.clone() {
                let actual = self.residues[r as usize].chain;
                if actual != c as u32 {
                    return Err(TopologyError::ResidueChain {
                        residue: r as usize,
                        expected: c as u32,
                        actual,
                    });
                }
            }
            next_res = chain.residues.end;
        }
        if next_res as usize != self.residues.len() {
            return Err(TopologyError::ChainCoverage {
                covered: next_res,
                residues: self.residues.len(),
            });
        }

        let bond_atoms = self
            .explicit_bonds
            .iter()
            .flat_map(|b| b.atoms)
            .chain(
                self.bonds
                    .iter()
                    .flat_map(|b| b.pairs.iter().copied().flatten()),
            )
            .chain(self.md_bonds.iter().flatten().copied().flatten());
        for atom in bond_atoms {
            if atom as usize >= atoms {
                return Err(TopologyError::BondAtom { atom, atoms });
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn two_residue_topology() -> Topology {
        let mut names = Interner::new();
        let ala = names.intern("ALA");
        let a = names.intern("A");
        Topology {
            element: vec![
                Element::NITROGEN,
                Element::CARBON,
                Element::CARBON,
                Element::OXYGEN,
            ],
            residue_index: vec![0, 0, 1, 1],
            residues: vec![
                ResidueRec {
                    atoms: 0..2,
                    chain: 0,
                    comp: ala,
                    seq_id: 1,
                    auth_seq_id: 1,
                    ins_code: 0,
                    ss: SecondaryStructure::Unknown,
                },
                ResidueRec {
                    atoms: 2..4,
                    chain: 0,
                    comp: ala,
                    seq_id: 2,
                    auth_seq_id: 2,
                    ins_code: 0,
                    ss: SecondaryStructure::Unknown,
                },
            ],
            chains: vec![ChainRec {
                residues: 0..2,
                label_asym: a,
                auth_asym: a,
                entity: 1,
            }],
            names,
            ..Default::default()
        }
    }

    #[test]
    fn valid_topology_passes() {
        let t = two_residue_topology();
        assert_eq!(t.validate(), Ok(()));
        assert_eq!(t.atom_count(), 4);
        assert_eq!(t.residue_count(), 2);
        assert_eq!(t.chain_count(), 1);
        assert_eq!(t.residue_name(1), "ALA");
        assert_eq!(t.chain_name(0), "A");
        assert_eq!(t.chain_of_atom(3), 0);
    }

    #[test]
    fn empty_topology_is_valid() {
        assert_eq!(Topology::default().validate(), Ok(()));
    }

    #[test]
    fn gaps_and_mismatches_are_reported() {
        let mut t = two_residue_topology();
        t.residues[1].atoms = 3..4;
        assert!(matches!(
            t.validate(),
            Err(TopologyError::ResidueRange { residue: 1, .. })
        ));

        let mut t = two_residue_topology();
        t.residue_index[3] = 0;
        assert!(matches!(
            t.validate(),
            Err(TopologyError::ResidueIndex { atom: 3, .. })
        ));

        let mut t = two_residue_topology();
        t.chains[0].residues = 0..1;
        assert!(matches!(
            t.validate(),
            Err(TopologyError::ChainCoverage { .. })
        ));

        let mut t = two_residue_topology();
        t.residue_index.pop();
        assert!(matches!(
            t.validate(),
            Err(TopologyError::ColumnLength { .. })
        ));

        let mut t = two_residue_topology();
        t.explicit_bonds.push(ExplicitBond {
            atoms: [0, 9],
            kind: ExplicitBondKind::Covalent,
            order: BondOrder::Single,
        });
        assert!(matches!(
            t.validate(),
            Err(TopologyError::BondAtom { atom: 9, .. })
        ));
    }
}
