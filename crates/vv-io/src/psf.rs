//! Full topology from a CHARMM/NAMD/X-PLOR PSF file (standard, EXT and
//! CHEQ variants).
//!
//! Fields are whitespace-separated, not column-fixed: EXT (wider serial
//! and name columns, for >99,999 atoms) and CHEQ (two extra charge-
//! equilibration columns per atom) both parse the same way, since neither
//! changes the field order or introduces internal whitespace. Read per
//! atom: segid (-> chain), resid, resname, atom name, charge, mass;
//! element from mass (`crate::guess_element`, name as fallback), since a
//! PSF's own "atom type" field is a force-field type code, not an
//! element and is otherwise unused. The `!NBOND` section becomes
//! `Topology::md_bonds`, used verbatim (`vv_core::bonds` module doc).
//!
//! A PSF alone has no coordinates: `parse` returns just a `Topology`, to
//! be paired with a coordinate file or trajectory of the same atom count
//! (`vv_io::load_topology`).

use vv_core::glam::Vec3;
use vv_core::{AtomRow, Topology, TopologyBuilder};

use crate::ParseError;

/// Reads a whole PSF file: atoms, residues and chains (from segid), plus
/// the bond list ([`read_bonds`]) as `Topology::md_bonds`.
pub fn parse(text: &str) -> Result<Topology, ParseError> {
    let (header_line, count) = natom_header(text)?;
    if count == 0 {
        return Err(ParseError::NoAtoms);
    }

    let mut builder = TopologyBuilder::with_capacity(count);
    for (offset, line) in text.lines().skip(header_line + 1).take(count).enumerate() {
        let lineno = header_line + 2 + offset;
        let atom = AtomFields::parse(line, lineno)?;
        builder.push(&atom.as_row());
    }
    if builder.atom_count() != count {
        return Err(malformed(
            header_line,
            "!NATOM has fewer atom lines than its count",
        ));
    }

    let mut topology = builder.topology;
    topology.md_bonds = Some(read_bonds(text)?);
    Ok(topology)
}

fn natom_header(text: &str) -> Result<(usize, usize), ParseError> {
    text.lines()
        .enumerate()
        .find_map(|(i, line)| {
            let count = line.split_whitespace().next()?.parse::<usize>().ok()?;
            line.contains("!NATOM").then_some((i, count))
        })
        .ok_or_else(|| malformed(0, "no !NATOM section"))
}

/// One `!NATOM` line's fields, borrowed from the source line. Field order
/// (index, segid, resid, resname, name, type, charge, mass, ...) is the
/// same for standard, EXT and CHEQ PSFs; only field width differs, which
/// whitespace splitting does not care about.
struct AtomFields<'a> {
    serial: u32,
    segid: &'a str,
    seq_id: i32,
    ins_code: u8,
    resname: &'a str,
    name: &'a str,
    charge: f32,
    mass: f32,
}

impl<'a> AtomFields<'a> {
    fn parse(line: &'a str, lineno: usize) -> Result<Self, ParseError> {
        let f: Vec<&str> = line.split_whitespace().collect();
        if f.len() < 8 {
            return Err(malformed(lineno, "!NATOM line has fewer than 8 fields"));
        }
        let int_err = |_| malformed(lineno, "non-integer atom index in !NATOM");
        let num_err = |_| malformed(lineno, "unreadable charge or mass in !NATOM");
        let (seq_id, ins_code) = parse_resid(f[2])
            .ok_or_else(|| malformed(lineno, "unreadable residue id in !NATOM"))?;
        Ok(AtomFields {
            serial: f[0].parse().map_err(int_err)?,
            segid: f[1],
            seq_id,
            ins_code,
            resname: f[3],
            name: f[4],
            // f[5] is the force-field atom type, not needed for anything
            // this data model tracks.
            charge: f[6].parse().map_err(num_err)?,
            mass: f[7].parse().map_err(num_err)?,
        })
    }

    fn as_row(&self) -> AtomRow<'a> {
        AtomRow {
            element: crate::guess_element(self.mass, self.name),
            name: crate::pad_name4(self.name),
            serial: self.serial,
            alt_loc: 0,
            comp: self.resname,
            asym: self.segid,
            auth_asym: self.segid,
            seq_id: self.seq_id,
            auth_seq_id: self.seq_id,
            ins_code: self.ins_code,
            entity: 0,
            // A PSF has no coordinates; the caller pairs this topology
            // with a coordinate file or trajectory of the same atom count.
            position: Vec3::ZERO,
            // A PSF has neither column; occupancy defaults as if fully
            // present, b-factor as if unmeasured.
            occupancy: 1.0,
            b_factor: 0.0,
            charge: crate::round_partial_charge(self.charge),
            hetero: false,
        }
    }
}

/// CHARMM's `resid` is normally numeric, but some CHARMM-GUI/PDB-derived
/// segments carry a trailing insertion letter (`82A`); split it off
/// rather than fail the line.
fn parse_resid(tok: &str) -> Option<(i32, u8)> {
    if let Ok(v) = tok.parse::<i32>() {
        return Some((v, 0));
    }
    let digits = tok.trim_end_matches(|c: char| c.is_ascii_alphabetic());
    if digits.is_empty() || digits.len() == tok.len() {
        return None;
    }
    let letter = tok.as_bytes()[digits.len()];
    digits.parse::<i32>().ok().map(|v| (v, letter))
}

/// Bond list from a PSF file's `!NBOND` section, as 0-based atom index
/// pairs in file order (the PSF format numbers atoms from 1).
pub fn read_bonds(text: &str) -> Result<Vec<[u32; 2]>, ParseError> {
    let (header_line, count) = text
        .lines()
        .enumerate()
        .find_map(|(i, line)| {
            let count = line.split_whitespace().next()?.parse::<usize>().ok()?;
            line.contains("!NBOND").then_some((i, count))
        })
        .ok_or_else(|| malformed(0, "no !NBOND section"))?;

    let mut values = text
        .lines()
        .skip(header_line + 1)
        .take_while(|line| !line.trim_start().starts_with('!') && !line.trim().is_empty())
        .flat_map(str::split_whitespace)
        .map(|tok| {
            tok.parse::<u32>()
                .map_err(|_| malformed(header_line, "non-integer atom serial in !NBOND"))
        });

    let mut bonds = Vec::with_capacity(count);
    for _ in 0..count {
        let a = values
            .next()
            .ok_or_else(|| malformed(header_line, "!NBOND has fewer pairs than its count"))??;
        let b = values
            .next()
            .ok_or_else(|| malformed(header_line, "!NBOND has fewer pairs than its count"))??;
        if a == 0 || b == 0 {
            return Err(malformed(header_line, "!NBOND atom serial is not 1-based"));
        }
        bonds.push([a - 1, b - 1]);
    }
    Ok(bonds)
}

fn malformed(line: usize, message: &str) -> ParseError {
    ParseError::Malformed {
        line: line + 1,
        message: message.to_string(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use vv_core::Element;

    const SAMPLE: &str = "PSF EXT\n\
\n\
       1 !NTITLE\n\
 REMARKS example\n\
\n\
       3 !NATOM\n\
       1 A    1    ALA  N    NH3   -0.300000       14.0070           0\n\
       2 A    1    ALA  CA   CT1    0.100000       12.0110           0\n\
       3 A    1    ALA  C    C      0.500000       12.0110           0\n\
\n\
       2 !NBOND: bonds\n\
       1       2       2       3\n";

    #[test]
    fn reads_bond_pairs_zero_based() {
        let bonds = read_bonds(SAMPLE).unwrap();
        assert_eq!(bonds, vec![[0, 1], [1, 2]]);
    }

    #[test]
    fn bonds_wrap_across_lines() {
        // Real PSF files wrap at 4 pairs (8 integers) per line; make sure
        // a count that spans multiple lines is read correctly.
        let mut text =
            String::from("PSF\n\n       0 !NTITLE\n\n       0 !NATOM\n\n       5 !NBOND: bonds\n");
        text.push_str("       1       2       2       3       3       4       4       5\n");
        text.push_str("       5       6\n");
        let bonds = read_bonds(&text).unwrap();
        assert_eq!(bonds, vec![[0, 1], [1, 2], [2, 3], [3, 4], [4, 5]]);
    }

    #[test]
    fn missing_section_is_an_error() {
        assert!(read_bonds("PSF\n\n0 !NATOM\n").is_err());
    }

    #[test]
    fn truncated_section_is_an_error() {
        assert!(read_bonds("PSF\n\n       2 !NBOND: bonds\n       1       2\n").is_err());
    }

    #[test]
    fn parses_atoms_chains_residues_elements_and_charges() {
        let t = parse(SAMPLE).unwrap();
        assert_eq!(t.validate(), Ok(()));
        assert_eq!(t.atom_count(), 3);
        assert_eq!(t.residue_count(), 1);
        assert_eq!(t.chain_count(), 1);
        assert_eq!(t.chain_name(0), "A");
        assert_eq!(t.residue_name(0), "ALA");
        assert_eq!(t.residues[0].seq_id, 1);
        assert_eq!(t.atom_name(0), "N");
        assert_eq!(t.element[0], Element::NITROGEN); // mass 14.007
        assert_eq!(t.element[1], Element::CARBON); // mass 12.011
        assert_eq!(t.charge[0], 0, "-0.3 e is not within 0.1 of an integer");
        assert_eq!(
            t.charge[2], 0,
            "0.5 e is the least confident case: kept neutral"
        );
        assert_eq!(t.md_bonds, Some(vec![[0, 1], [1, 2]]));
    }

    #[test]
    fn a_psf_alone_has_no_coordinates() {
        // `parse` returns a bare `Topology`, not a `Structure`: nothing to
        // check for coordinates, which is the point.
        let t = parse(SAMPLE).unwrap();
        let _: vv_core::Topology = t;
    }

    #[test]
    fn empty_natom_is_no_atoms_error() {
        let text = "PSF\n\n       0 !NTITLE\n\n       0 !NATOM\n\n       0 !NBOND: bonds\n";
        assert!(matches!(parse(text), Err(ParseError::NoAtoms)));
    }

    #[test]
    fn declared_bonds_are_used_verbatim_even_when_geometry_would_disagree() {
        // Two atoms placed at a clash distance the PSF does not bond, and
        // two placed far apart that it does: `md_bonds` must reflect the
        // PSF exactly, since `vv_core::bonds::perceive` trusts it verbatim.
        let text = "PSF\n\
\n\
       1 !NTITLE\n\
 REMARKS\n\
\n\
       3 !NATOM\n\
       1 A    1    LIG  C1   C      0.000000       12.0110           0\n\
       2 A    1    LIG  C2   C      0.000000       12.0110           0\n\
       3 A    1    LIG  C3   C      0.000000       12.0110           0\n\
\n\
       1 !NBOND: bonds\n\
       2       3\n";
        let t = parse(text).unwrap();
        let p = [
            vv_core::glam::Vec3::ZERO,
            vv_core::glam::Vec3::new(0.6, 0.0, 0.0), // clash distance from atom 0
            vv_core::glam::Vec3::new(20.0, 0.0, 0.0),
        ];
        let bonds = vv_core::bonds::perceive(&t, &p);
        assert_eq!(bonds.pairs, vec![[1, 2]]);
    }

    #[test]
    fn resid_with_an_insertion_letter_does_not_fail_the_line() {
        let text = "PSF\n\n       0 !NTITLE\n\n       1 !NATOM\n\
       1 A    82A  ALA  N    NH3   -0.300000       14.0070           0\n\
\n       0 !NBOND: bonds\n";
        let t = parse(text).unwrap();
        assert_eq!(t.residues[0].seq_id, 82);
    }

    #[test]
    fn selection_keywords_recognize_charmm_residue_names() {
        let text = "PSF\n\n       0 !NTITLE\n\n       4 !NATOM\n\
       1 A    1    ALA  N    NH3   -0.300000       14.0070           0\n\
       2 B    1    TIP3 OH2  OT    -0.834000       15.9994           0\n\
       3 C    1    SOD  SOD   SOD   1.000000       22.9898           0\n\
       4 D    1    CLA  CLA   CLA  -1.000000       35.4500           0\n\
\n       0 !NBOND: bonds\n";
        let t = parse(text).unwrap();
        let p = vec![vv_core::glam::Vec3::ZERO; 4];
        let protein = vv_core::select(&t, &p, "protein").unwrap();
        let water = vv_core::select(&t, &p, "water").unwrap();
        let ion = vv_core::select(&t, &p, "ion").unwrap();
        assert!(protein[0] && !protein[1] && !protein[2] && !protein[3]);
        assert!(!water[0] && water[1] && !water[2] && !water[3]);
        assert!(!ion[0] && !ion[1] && ion[2] && ion[3]);
        // Ions' near-integer partial charges round to a formal charge.
        assert_eq!(t.charge[2], 1);
        assert_eq!(t.charge[3], -1);
    }
}
