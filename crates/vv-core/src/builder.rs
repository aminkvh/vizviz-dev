//! Incremental construction of a `Topology` + coordinates from atom rows.
//!
//! Parsers feed rows in file order; residue and chain boundaries are
//! detected from the row keys. Several builders (one per parallel chunk)
//! can be merged in order, which stitches a residue split across a chunk
//! seam back together.

use glam::Vec3;

use crate::{
    flags, ChainRec, CoordSet, Element, InternId, ResidueRec, SecondaryStructure, Structure,
    StructureError, Topology,
};

/// One atom record as read from a file. Borrowed strings avoid allocation.
#[derive(Clone, Copy, Debug)]
pub struct AtomRow<'a> {
    pub element: Element,
    pub name: [u8; 4],
    pub serial: u32,
    pub alt_loc: u8,
    pub comp: &'a str,
    pub asym: &'a str,
    pub auth_asym: &'a str,
    pub seq_id: i32,
    pub auth_seq_id: i32,
    pub ins_code: u8,
    pub entity: u16,
    pub position: Vec3,
    pub occupancy: f32,
    pub b_factor: f32,
    pub charge: i8,
    pub hetero: bool,
}

/// Per-atom fields most rows leave at their default: a name longer than
/// the 4-byte column, the chain's segment id, the deuterium isotope.
#[derive(Clone, Copy, Debug, Default)]
pub struct AtomExtra<'a> {
    pub long_name: Option<&'a str>,
    /// Read when the row opens a new chain record.
    pub segid: &'a str,
    pub deuterium: bool,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct ResidueKey {
    asym: InternId,
    seq_id: i32,
    ins_code: u8,
    comp: InternId,
}

#[derive(Debug, Default)]
pub struct TopologyBuilder {
    pub topology: Topology,
    pub positions: Vec<Vec3>,
    current: Option<ResidueKey>,
}

impl TopologyBuilder {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn with_capacity(atoms: usize) -> Self {
        let mut b = Self::new();
        let t = &mut b.topology;
        t.element.reserve(atoms);
        t.name.reserve(atoms);
        t.serial.reserve(atoms);
        t.b_factor.reserve(atoms);
        t.occupancy.reserve(atoms);
        t.alt_loc.reserve(atoms);
        t.charge.reserve(atoms);
        t.flags.reserve(atoms);
        t.residue_index.reserve(atoms);
        t.residues.reserve(atoms / 8);
        b.positions.reserve(atoms);
        b
    }

    pub fn atom_count(&self) -> usize {
        self.topology.element.len()
    }

    /// Starts a new chain record with the next pushed row even when its
    /// chain name repeats (a PDB `TER` between same-lettered molecules).
    pub fn end_chain(&mut self) {
        self.current = None;
    }

    pub fn push(&mut self, row: &AtomRow<'_>) {
        self.push_with(row, &AtomExtra::default());
    }

    pub fn push_with(&mut self, row: &AtomRow<'_>, extra: &AtomExtra<'_>) {
        let t = &mut self.topology;
        let asym = t.names.intern(row.asym);
        let comp = t.names.intern(row.comp);
        let key = ResidueKey {
            asym,
            seq_id: row.seq_id,
            ins_code: row.ins_code,
            comp,
        };
        let atom = t.element.len() as u32;

        if self.current.as_ref() != Some(&key) {
            let chain_changed = self.current.as_ref().is_none_or(|c| c.asym != asym);
            if chain_changed {
                let auth_asym = t.names.intern(row.auth_asym);
                let start = t.residues.len() as u32;
                let segid = t.names.intern(extra.segid);
                t.segids.push(segid);
                t.chains.push(ChainRec {
                    residues: start..start,
                    label_asym: asym,
                    auth_asym,
                    entity: row.entity,
                });
            }
            let chain = t.chains.len() as u32 - 1;
            t.residues.push(ResidueRec {
                atoms: atom..atom,
                chain,
                comp,
                seq_id: row.seq_id,
                auth_seq_id: row.auth_seq_id,
                ins_code: row.ins_code,
                ss: SecondaryStructure::Unknown,
            });
            t.chains.last_mut().unwrap().residues.end += 1;
            self.current = Some(key);
        }

        let residue = t.residues.len() as u32 - 1;
        t.residues[residue as usize].atoms.end += 1;
        t.element.push(row.element);
        t.name.push(row.name);
        t.serial.push(row.serial);
        t.b_factor.push(row.b_factor);
        t.occupancy.push(row.occupancy);
        t.alt_loc.push(row.alt_loc);
        t.charge.push(row.charge);
        let mut atom_flags = if row.hetero { flags::HETERO } else { 0 };
        if extra.deuterium {
            atom_flags |= flags::DEUTERIUM;
        }
        t.flags.push(atom_flags);
        if let Some(long) = extra.long_name {
            t.long_names.push((atom, long.into()));
        }
        t.residue_index.push(residue);
        self.positions.push(row.position);
    }

    /// Appends `other` (built from the rows that followed this builder's
    /// rows in the file), stitching a residue split across the seam.
    pub fn append(&mut self, mut other: TopologyBuilder) {
        if other.atom_count() == 0 {
            return;
        }
        if self.atom_count() == 0 {
            *self = other;
            return;
        }
        let t = &mut self.topology;
        let o = &mut other.topology;

        // Remap the other builder's interned names into ours.
        let remap: Vec<InternId> = (0..o.names.len())
            .map(|i| t.names.intern(o.names.get(InternId(i as u32))))
            .collect();
        let atom_base = t.element.len() as u32;

        let first_key = {
            let r = &o.residues[0];
            ResidueKey {
                asym: remap[r.chain_asym(o).0 as usize],
                seq_id: r.seq_id,
                ins_code: r.ins_code,
                comp: remap[r.comp.0 as usize],
            }
        };
        let same_residue = self.current.as_ref() == Some(&first_key);
        let same_chain = self
            .current
            .as_ref()
            .is_some_and(|c| c.asym == first_key.asym);

        let residue_base = t.residues.len() as u32 - u32::from(same_residue);
        let chain_base = t.chains.len() as u32 - u32::from(same_chain);

        for (i, r) in o.residues.iter().enumerate() {
            if i == 0 && same_residue {
                t.residues.last_mut().unwrap().atoms.end += r.atoms.len() as u32;
                continue;
            }
            t.residues.push(ResidueRec {
                atoms: r.atoms.start + atom_base..r.atoms.end + atom_base,
                chain: r.chain + chain_base,
                comp: remap[r.comp.0 as usize],
                seq_id: r.seq_id,
                auth_seq_id: r.auth_seq_id,
                ins_code: r.ins_code,
                ss: r.ss,
            });
        }
        for (i, c) in o.chains.iter().enumerate() {
            let residues = c.residues.start + residue_base..c.residues.end + residue_base;
            if i == 0 && same_chain {
                t.chains.last_mut().unwrap().residues.end = residues.end;
                continue;
            }
            t.segids.push(remap[o.segids[i].0 as usize]);
            t.chains.push(ChainRec {
                residues,
                label_asym: remap[c.label_asym.0 as usize],
                auth_asym: remap[c.auth_asym.0 as usize],
                entity: c.entity,
            });
        }

        t.element.append(&mut o.element);
        t.name.append(&mut o.name);
        t.serial.append(&mut o.serial);
        t.b_factor.append(&mut o.b_factor);
        t.occupancy.append(&mut o.occupancy);
        t.alt_loc.append(&mut o.alt_loc);
        t.charge.append(&mut o.charge);
        t.flags.append(&mut o.flags);
        t.long_names.extend(
            o.long_names
                .drain(..)
                .map(|(atom, name)| (atom + atom_base, name)),
        );
        t.residue_index
            .extend(o.residue_index.iter().map(|r| r + residue_base));
        self.positions.append(&mut other.positions);
        self.current = other.current.map(|c| ResidueKey {
            asym: remap[c.asym.0 as usize],
            seq_id: c.seq_id,
            ins_code: c.ins_code,
            comp: remap[c.comp.0 as usize],
        });
    }

    pub fn finish(self) -> Result<Structure, StructureError> {
        self.finish_with_frames(Vec::new())
    }

    /// `extra_frames` are additional coordinate sets (models) with the same
    /// atom count as the built topology.
    pub fn finish_with_frames(
        mut self,
        extra_frames: Vec<Vec<Vec3>>,
    ) -> Result<Structure, StructureError> {
        fix_ion_elements(&mut self.topology);
        let mut frames = Vec::with_capacity(1 + extra_frames.len());
        frames.push(CoordSet::new(self.positions));
        frames.extend(extra_frames.into_iter().map(CoordSet::new));
        Structure::with_frames(self.topology, frames)
    }
}

/// Overrides a single-atom residue's element from its residue name
/// (`Element::ion_element`) when the two disagree: a source file's atom
/// name/column convention can misidentify a monatomic ion (a calcium ion
/// named `CA`, column-aligned like an alpha carbon, reads as carbon; see
/// `Element::from_atom_name`'s docs), but its residue name cannot, once
/// it is known to be the residue's only atom. A multi-atom residue is
/// left alone: `CO` and `NI` are also plausible ligand/fragment names,
/// and only a lone atom makes the ion reading unambiguous.
fn fix_ion_elements(topology: &mut Topology) {
    for residue in &topology.residues {
        if residue.atoms.len() != 1 {
            continue;
        }
        if let Some(element) = crate::element::ion_element(topology.names.get(residue.comp)) {
            topology.element[residue.atoms.start as usize] = element;
        }
    }
}

impl ResidueRec {
    fn chain_asym(&self, t: &Topology) -> InternId {
        t.chains[self.chain as usize].label_asym
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn row<'a>(asym: &'a str, seq: i32, comp: &'a str, name: &str, serial: u32) -> AtomRow<'a> {
        let mut n = [b' '; 4];
        n[..name.len()].copy_from_slice(name.as_bytes());
        AtomRow {
            element: Element::CARBON,
            name: n,
            serial,
            alt_loc: 0,
            comp,
            asym,
            auth_asym: asym,
            seq_id: seq,
            auth_seq_id: seq,
            ins_code: 0,
            entity: 1,
            position: Vec3::new(serial as f32, 0.0, 0.0),
            occupancy: 1.0,
            b_factor: 10.0,
            charge: 0,
            hetero: false,
        }
    }

    #[test]
    fn fixes_a_misread_ion_element_but_leaves_a_real_alpha_carbon_alone() {
        let mut b = TopologyBuilder::new();
        // A calcium ion, its single atom's element misread as carbon (as
        // `Element::from_atom_name` would for a non-canonically-aligned
        // " CA " atom name with no element column).
        b.push(&row("A", 1, "CA", "CA", 1));
        // An alanine's own alpha carbon: same atom name, same (correct)
        // misread-as-carbon element, but not a single-atom residue.
        b.push(&row("A", 2, "ALA", "N", 2));
        b.push(&row("A", 2, "ALA", "CA", 3));
        let s = b.finish().unwrap();
        let t = &s.topology;
        assert_eq!(t.element[0], Element::from_atomic_number(20).unwrap());
        assert_eq!(t.element[2], Element::CARBON);
    }

    #[test]
    fn detects_residue_and_chain_boundaries() {
        let mut b = TopologyBuilder::new();
        b.push(&row("A", 1, "ALA", "N", 1));
        b.push(&row("A", 1, "ALA", "CA", 2));
        b.push(&row("A", 2, "GLY", "N", 3));
        b.push(&row("B", 1, "SER", "N", 4));
        let s = b.finish().unwrap();
        let t = &s.topology;
        assert_eq!(t.validate(), Ok(()));
        assert_eq!(t.residue_count(), 3);
        assert_eq!(t.chain_count(), 2);
        assert_eq!(t.residues[0].atoms, 0..2);
        assert_eq!(t.residues[2].chain, 1);
        assert_eq!(t.chains[1].residues, 2..3);
        assert_eq!(t.atom_name(1), "CA");
        assert_eq!(s.frame(0).positions()[3].x, 4.0);
    }

    #[test]
    fn append_stitches_a_residue_split_across_the_seam() {
        let mut a = TopologyBuilder::new();
        a.push(&row("A", 1, "ALA", "N", 1));
        a.push(&row("A", 2, "GLY", "N", 2));
        a.push(&row("A", 2, "GLY", "CA", 3));
        let mut b = TopologyBuilder::new();
        b.push(&row("A", 2, "GLY", "C", 4));
        b.push(&row("A", 3, "SER", "N", 5));
        b.push(&row("B", 1, "SER", "N", 6));
        a.append(b);
        let s = a.finish().unwrap();
        let t = &s.topology;
        assert_eq!(t.validate(), Ok(()));
        assert_eq!(t.residue_count(), 4);
        assert_eq!(t.residues[1].atoms, 1..4);
        assert_eq!(t.residue_name(1), "GLY");
        assert_eq!(t.chain_count(), 2);
        assert_eq!(t.chains[0].residues, 0..3);
        assert_eq!(t.chains[1].residues, 3..4);
        assert_eq!(t.residue_index, vec![0, 1, 1, 1, 2, 3]);
        assert_eq!(t.serial, vec![1, 2, 3, 4, 5, 6]);
    }

    #[test]
    fn append_across_a_chain_boundary_starts_a_new_chain() {
        let mut a = TopologyBuilder::new();
        a.push(&row("A", 1, "ALA", "N", 1));
        let mut b = TopologyBuilder::new();
        b.push(&row("B", 1, "ALA", "N", 2));
        a.append(b);
        let t = a.finish().unwrap().topology;
        assert_eq!(t.validate(), Ok(()));
        assert_eq!(t.chain_count(), 2);
        assert_eq!(t.chain_name(1), "B");
    }

    #[test]
    fn extras_are_kept_and_shifted_when_builders_are_appended() {
        let mut a = TopologyBuilder::new();
        let extra = AtomExtra {
            segid: "SEG1",
            ..Default::default()
        };
        a.push_with(&row("A", 1, "ALA", "N", 1), &extra);
        let mut b = TopologyBuilder::new();
        let extra = AtomExtra {
            long_name: Some("N12345"),
            segid: "SEG2",
            deuterium: true,
        };
        b.push_with(&row("B", 1, "ALA", "N123", 2), &extra);
        a.append(b);
        let t = a.finish().unwrap().topology;
        assert_eq!(t.validate(), Ok(()));
        assert_eq!((t.segid(0), t.segid(1)), ("SEG1", "SEG2"));
        assert_eq!(t.long_names, vec![(1, "N12345".into())]);
        assert_eq!((t.atom_name(0), t.atom_name(1)), ("N", "N12345"));
        assert!(!t.is_deuterium(0) && t.is_deuterium(1));
    }

    #[test]
    fn append_into_empty_builder() {
        let mut a = TopologyBuilder::new();
        let mut b = TopologyBuilder::new();
        b.push(&row("A", 1, "ALA", "N", 1));
        a.append(b);
        assert_eq!(a.atom_count(), 1);
        let mut c = TopologyBuilder::new();
        c.push(&row("A", 1, "ALA", "CA", 2));
        a.append(c);
        assert_eq!(a.topology.residue_count(), 1);
        assert_eq!(a.topology.residues[0].atoms, 0..2);
    }
}
