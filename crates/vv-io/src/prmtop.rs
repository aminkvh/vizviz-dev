//! Full topology from an AMBER PRMTOP file (`%FLAG`/`%FORMAT` sections).
//!
//! A section's FORTRAN `%FORMAT` line (`20a4`, `5E16.8`, `10I8`, ...) gives
//! the fixed width of each of its fields; the line-wrapping itself is not
//! meaningful, so every data line up to the next `%FLAG`/`%FORMAT` is
//! concatenated before slicing into `width`-wide fields ([`section`]). A
//! `%COMMENT` line, which some AmberTools versions insert between `%FLAG`
//! and `%FORMAT`, is skipped.
//!
//! Read per atom: `ATOM_NAME`, `CHARGE` (÷18.2223 to elementary charge),
//! `MASS`, `RESIDUE_LABEL`/`RESIDUE_POINTER` for the residue, `ATOMIC_NUMBER`
//! when present (else element from mass, atom name as further fallback --
//! `crate::guess_element`). Chains come from `ATOMS_PER_MOLECULE`, each
//! solute molecule (before `SOLVENT_POINTERS`' first solvent molecule) its
//! own chain and every solvent molecule lumped into one trailing chain --
//! otherwise a topology with thousands of individually numbered waters
//! would get thousands of chains. Neither flag present: one chain. Both
//! bond sections (`BONDS_WITHOUT_HYDROGEN`, `BONDS_INC_HYDROGEN`) become
//! `Topology::md_bonds`, used verbatim (`vv_core::bonds` module doc).
//!
//! A PRMTOP alone has no coordinates: `parse` returns just a `Topology`,
//! to be paired with a coordinate file or trajectory of the same atom
//! count (`vv_io::load_topology`), typically the matching `.inpcrd`/
//! `.rst7` or a NetCDF/DCD trajectory from the same run.

use vv_core::glam::Vec3;
use vv_core::{AtomRow, Element, Topology, TopologyBuilder};

use crate::ParseError;

/// AMBER stores charge in an internal unit; dividing by this gives
/// elementary charge e (AMBER file format specification: 18.2223 =
/// sqrt(332.0636), Coulomb's constant in AMBER's kcal/mol, angstrom, e
/// units).
const AMBER_CHARGE_SCALE: f32 = 18.2223;

/// Reads a whole PRMTOP file: atoms, residues, chains (from
/// `ATOMS_PER_MOLECULE`/`SOLVENT_POINTERS` when present), plus the bond
/// lists ([`read_bonds`]) as `Topology::md_bonds`.
pub fn parse(text: &str) -> Result<Topology, ParseError> {
    let pointers = ints(text, "POINTERS")?;
    let natom = *pointers
        .first()
        .ok_or_else(|| malformed(0, "POINTERS is empty"))? as usize;
    let nres = *pointers
        .get(11)
        .ok_or_else(|| malformed(0, "POINTERS has fewer than 12 entries"))? as usize;
    if natom == 0 {
        return Err(ParseError::NoAtoms);
    }

    let names = strings(text, "ATOM_NAME")?;
    let charges = floats(text, "CHARGE")?;
    let masses = floats(text, "MASS")?;
    let atomic_numbers = optional_ints(text, "ATOMIC_NUMBER")?;
    let res_labels = strings(text, "RESIDUE_LABEL")?;
    let res_pointer = ints(text, "RESIDUE_POINTER")?;
    check_len("ATOM_NAME", names.len(), natom)?;
    check_len("CHARGE", charges.len(), natom)?;
    check_len("MASS", masses.len(), natom)?;
    check_len("RESIDUE_LABEL", res_labels.len(), nres)?;
    check_len("RESIDUE_POINTER", res_pointer.len(), nres)?;

    let residue_of = residue_of_atom(natom, &res_pointer);
    let atoms_per_molecule = optional_ints(text, "ATOMS_PER_MOLECULE")?;
    let solvent_pointers = optional_ints(text, "SOLVENT_POINTERS")?;
    let chain_of = chain_labels(
        natom,
        atoms_per_molecule.as_deref(),
        solvent_pointers.as_deref(),
    );

    let mut builder = TopologyBuilder::with_capacity(natom);
    for a in 0..natom {
        let r = residue_of[a];
        let element = match atomic_numbers.as_ref().and_then(|z| z.get(a)).copied() {
            Some(z) if z > 0 => Element::from_atomic_number(z as u8).unwrap_or(Element::UNKNOWN),
            _ => crate::guess_element(masses[a], &names[a]),
        };
        let row = AtomRow {
            element,
            name: crate::pad_name4(&names[a]),
            serial: a as u32 + 1,
            alt_loc: 0,
            comp: &res_labels[r],
            asym: &chain_of[a],
            auth_asym: &chain_of[a],
            seq_id: r as i32 + 1,
            auth_seq_id: r as i32 + 1,
            ins_code: 0,
            entity: 0,
            // A PRMTOP has no coordinates; the caller pairs this topology
            // with a coordinate file or trajectory of the same atom count.
            position: Vec3::ZERO,
            // A PRMTOP has neither column; occupancy defaults as if
            // fully present, b-factor as if unmeasured.
            occupancy: 1.0,
            b_factor: 0.0,
            charge: crate::round_partial_charge(charges[a] / AMBER_CHARGE_SCALE),
            hetero: false,
        };
        builder.push(&row);
    }

    let mut topology = builder.topology;
    topology.title = strings(text, "TITLE")
        .map(|v| v.concat().trim().to_string())
        .unwrap_or_default();
    topology.md_bonds = Some(read_bonds(text)?);
    Ok(topology)
}

/// Residue index (0-based) of every atom, from 1-based `RESIDUE_POINTER`
/// (each residue's first atom; a residue runs to the next one's, or to
/// `natom` for the last).
fn residue_of_atom(natom: usize, residue_pointer: &[i64]) -> Vec<usize> {
    let starts: Vec<usize> = residue_pointer
        .iter()
        .map(|&p| (p.max(1) - 1) as usize)
        .collect();
    let mut out = vec![0usize; natom];
    for r in 0..starts.len() {
        let start = starts[r].min(natom);
        let end = starts.get(r + 1).copied().unwrap_or(natom).min(natom);
        // A well-formed RESIDUE_POINTER is non-decreasing; guard against a
        // malformed one instead of panicking on a reversed range.
        if start < end {
            out[start..end].fill(r);
        }
    }
    out
}

/// One chain label per atom. `atoms_per_molecule` (`%FLAG
/// ATOMS_PER_MOLECULE`) splits atoms into molecules; `solvent_pointers`
/// (`%FLAG SOLVENT_POINTERS`, `[IPTRES, NSPM, NSPSOL]`) marks molecules
/// `NSPSOL..NSPM` (1-based) as solvent. Each solute molecule gets its own
/// chain letter; every solvent molecule shares one trailing chain, since
/// a solvated system's thousands of individual waters are not thousands
/// of meaningful chains. Without `SOLVENT_POINTERS`, every molecule is
/// treated as solute. Without `ATOMS_PER_MOLECULE`, everything is one
/// chain.
fn chain_labels(
    natom: usize,
    atoms_per_molecule: Option<&[i64]>,
    solvent_pointers: Option<&[i64]>,
) -> Vec<String> {
    let Some(apm) = atoms_per_molecule else {
        return vec!["A".to_string(); natom];
    };
    let total: i64 = apm.iter().sum();
    if total != natom as i64 || apm.iter().any(|&c| c < 0) {
        return vec!["A".to_string(); natom];
    }

    let first_solvent = solvent_pointers.and_then(|sp| sp.get(2).copied());
    let mut molecule_label = Vec::with_capacity(apm.len());
    let mut next = 0usize;
    let mut solvent_label: Option<String> = None;
    for m in 0..apm.len() {
        let is_solvent = first_solvent.is_some_and(|fs| (m as i64 + 1) >= fs);
        let label = if is_solvent {
            solvent_label
                .get_or_insert_with(|| {
                    let l = chain_letters(next);
                    next += 1;
                    l
                })
                .clone()
        } else {
            let l = chain_letters(next);
            next += 1;
            l
        };
        molecule_label.push(label);
    }

    let mut out = Vec::with_capacity(natom);
    for (m, &count) in apm.iter().enumerate() {
        out.extend(std::iter::repeat_n(
            molecule_label[m].clone(),
            count as usize,
        ));
    }
    out
}

/// Spreadsheet-style chain labels: 0 -> "A", 25 -> "Z", 26 -> "AA", ...
fn chain_letters(mut n: usize) -> String {
    let mut s = Vec::new();
    loop {
        s.push(b'A' + (n % 26) as u8);
        if n < 26 {
            break;
        }
        n = n / 26 - 1;
    }
    s.reverse();
    String::from_utf8(s).unwrap()
}

// ---------------------------------------------------------------------------
// `%FLAG`/`%FORMAT` section reading

/// A section's start line (0-based, for error messages), concatenated
/// data (line-wrapping removed) and per-field width from its `%FORMAT`.
fn section(text: &str, flag: &str) -> Option<(usize, String, usize)> {
    let lines: Vec<&str> = text.lines().collect();
    let start = lines.iter().position(|line| {
        let line = line.trim_start();
        line.starts_with("%FLAG") && line.split_whitespace().nth(1) == Some(flag)
    })?;
    let mut i = start + 1;
    while lines
        .get(i)
        .is_some_and(|l| l.trim_start().starts_with("%COMMENT"))
    {
        i += 1;
    }
    let width = field_width(lines.get(i)?)?;
    let data: String = lines[i + 1..]
        .iter()
        .take_while(|l| !l.trim_start().starts_with('%'))
        .flat_map(|l| l.chars())
        .collect();
    Some((start, data, width))
}

/// Field width from a FORTRAN edit descriptor (`%FORMAT(20a4)` -> 4,
/// `%FORMAT(5E16.8)` -> 16, `%FORMAT(10I8)` -> 8): skip the repeat count,
/// the type letter, then read the width digits.
fn field_width(format_line: &str) -> Option<usize> {
    let inner = format_line
        .trim_start()
        .strip_prefix("%FORMAT(")?
        .strip_suffix(')')?;
    let after_repeat = inner.trim_start_matches(|c: char| c.is_ascii_digit());
    let after_type = after_repeat.trim_start_matches(|c: char| c.is_ascii_alphabetic());
    after_type
        .chars()
        .take_while(|c| c.is_ascii_digit())
        .collect::<String>()
        .parse()
        .ok()
}

fn chunks(data: &str, width: usize) -> impl Iterator<Item = &str> {
    let bytes = data.as_bytes();
    let width = width.max(1);
    (0..bytes.len()).step_by(width).map(move |i| {
        let end = (i + width).min(bytes.len());
        std::str::from_utf8(&bytes[i..end]).unwrap_or("").trim()
    })
}

fn require(text: &str, flag: &'static str) -> Result<(usize, String, usize), ParseError> {
    section(text, flag).ok_or(ParseError::MissingSection(flag))
}

fn ints(text: &str, flag: &'static str) -> Result<Vec<i64>, ParseError> {
    let (start, data, width) = require(text, flag)?;
    chunks(&data, width)
        .map(|s| {
            s.parse::<i64>()
                .map_err(|_| malformed(start, &format!("non-integer value in %FLAG {flag}")))
        })
        .collect()
}

fn floats(text: &str, flag: &'static str) -> Result<Vec<f32>, ParseError> {
    let (start, data, width) = require(text, flag)?;
    chunks(&data, width)
        .map(|s| {
            s.parse::<f32>()
                .map_err(|_| malformed(start, &format!("non-numeric value in %FLAG {flag}")))
        })
        .collect()
}

fn strings(text: &str, flag: &'static str) -> Result<Vec<String>, ParseError> {
    let (_, data, width) = require(text, flag)?;
    Ok(chunks(&data, width).map(str::to_string).collect())
}

fn optional_ints(text: &str, flag: &'static str) -> Result<Option<Vec<i64>>, ParseError> {
    match section(text, flag) {
        None => Ok(None),
        Some((start, data, width)) => chunks(&data, width)
            .map(|s| {
                s.parse::<i64>()
                    .map_err(|_| malformed(start, &format!("non-integer value in %FLAG {flag}")))
            })
            .collect::<Result<Vec<_>, _>>()
            .map(Some),
    }
}

fn check_len(flag: &'static str, actual: usize, expected: usize) -> Result<(), ParseError> {
    if actual < expected {
        Err(ParseError::Malformed {
            line: 0,
            message: format!("%FLAG {flag} has {actual} entries, expected at least {expected}"),
        })
    } else {
        Ok(())
    }
}

fn malformed(line: usize, message: &str) -> ParseError {
    ParseError::Malformed {
        line: line + 1,
        message: message.to_string(),
    }
}

// ---------------------------------------------------------------------------
// Bonds

const SECTIONS: [&str; 2] = ["BONDS_WITHOUT_HYDROGEN", "BONDS_INC_HYDROGEN"];

/// Bond list from a PRMTOP file's bond sections, as 0-based atom index
/// pairs. Each entry is 3 integers (atom1*3, atom2*3, bond-type index);
/// the third is discarded, as vizviz does not model bond order/type.
pub fn read_bonds(text: &str) -> Result<Vec<[u32; 2]>, ParseError> {
    let mut bonds = Vec::new();
    for section in SECTIONS {
        bonds.extend(read_bond_section(text, section)?);
    }
    if bonds.is_empty() {
        return Err(ParseError::Malformed {
            line: 0,
            message: "no %FLAG BONDS_WITHOUT_HYDROGEN or BONDS_INC_HYDROGEN section".to_string(),
        });
    }
    Ok(bonds)
}

fn read_bond_section(text: &str, flag: &'static str) -> Result<Vec<[u32; 2]>, ParseError> {
    let values = match optional_ints(text, flag)? {
        Some(v) => v,
        None => return Ok(Vec::new()),
    };
    if values.len() % 3 != 0 {
        return Err(ParseError::Malformed {
            line: 0,
            message: format!("{flag} length is not a multiple of 3"),
        });
    }
    Ok(values
        .chunks_exact(3)
        .map(|c| [(c[0] / 3) as u32, (c[1] / 3) as u32])
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Builds a `%FLAG`/`%FORMAT` section by section, so a test never has
    /// to hand-count fixed-width columns: each field is formatted to
    /// exactly `width` characters (the only thing that has to line up),
    /// and line-wrapping -- which `section` ignores -- is not this
    /// builder's problem either.
    #[derive(Default)]
    struct Prmtop(String);

    impl Prmtop {
        fn ints(mut self, flag: &str, width: usize, values: &[i64]) -> Self {
            self.0
                .push_str(&format!("%FLAG {flag}\n%FORMAT(I{width})\n"));
            for v in values {
                self.0.push_str(&format!("{v:>width$}"));
            }
            self.0.push('\n');
            self
        }

        fn floats(mut self, flag: &str, width: usize, values: &[f32]) -> Self {
            self.0
                .push_str(&format!("%FLAG {flag}\n%FORMAT(E{width}.8)\n"));
            for v in values {
                self.0.push_str(&format!("{v:>width$.8}"));
            }
            self.0.push('\n');
            self
        }

        fn strings(mut self, flag: &str, width: usize, values: &[&str]) -> Self {
            self.0
                .push_str(&format!("%FLAG {flag}\n%FORMAT(a{width})\n"));
            for v in values {
                self.0.push_str(&format!("{v:<width$}"));
            }
            self.0.push('\n');
            self
        }

        fn build(self) -> String {
            self.0
        }
    }

    /// 5 atoms: an ALA backbone (N, CA, C, O) plus a lone Na+ ion, 2
    /// residues. Bonds: N-CA, CA-C, C-O (no hydrogens).
    fn sample() -> String {
        Prmtop::default()
            .ints("POINTERS", 8, &[5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2])
            .strings("ATOM_NAME", 4, &["N", "CA", "C", "O", "NA"])
            .floats(
                "CHARGE",
                16,
                &[-5.46, 1.82, 9.10, -10.34, AMBER_CHARGE_SCALE],
            )
            .floats("MASS", 16, &[14.007, 12.011, 12.011, 15.999, 22.9898])
            .strings("RESIDUE_LABEL", 4, &["ALA", "SOD"])
            .ints("RESIDUE_POINTER", 8, &[1, 5])
            .ints("BONDS_WITHOUT_HYDROGEN", 8, &[0, 3, 1, 3, 6, 1, 6, 9, 1])
            .ints("BONDS_INC_HYDROGEN", 8, &[])
            .build()
    }

    #[test]
    fn reads_both_bond_sections_as_zero_based_pairs() {
        // atom*3 encoding: 0->0, 3->1, 6->2, 9->3.
        assert_eq!(read_bonds(&sample()).unwrap(), vec![[0, 1], [1, 2], [2, 3]]);
    }

    #[test]
    fn missing_sections_is_an_error() {
        let text = Prmtop::default().ints("POINTERS", 8, &[0]).build();
        assert!(read_bonds(&text).is_err());
    }

    #[test]
    fn parses_atoms_residues_charges_and_masses() {
        let t = parse(&sample()).unwrap();
        assert_eq!(t.validate(), Ok(()));
        assert_eq!(t.atom_count(), 5);
        assert_eq!(t.residue_count(), 2);
        assert_eq!(t.residue_name(0), "ALA");
        assert_eq!(t.residue_name(1), "SOD");
        assert_eq!(t.atom_name(0), "N");
        assert_eq!(t.atom_name(4), "NA");
        assert_eq!(t.element[0], Element::NITROGEN); // mass 14.007
        assert_eq!(t.element[1], Element::CARBON); // mass 12.011
        assert_eq!(t.element[3], Element::OXYGEN); // mass 15.999
        assert_eq!(t.element[4], Element::from_atomic_number(11).unwrap()); // Na, mass 22.9898
                                                                            // -5.46 / 18.2223 = -0.2997 e: not within 0.1 of an integer.
        assert_eq!(t.charge[0], 0);
        // AMBER_CHARGE_SCALE / AMBER_CHARGE_SCALE = 1.0 e exactly: Na+'s
        // formal charge.
        assert_eq!(t.charge[4], 1);
        assert_eq!(t.md_bonds, Some(vec![[0, 1], [1, 2], [2, 3]]));
    }

    #[test]
    fn atomic_number_takes_priority_over_mass() {
        // Atom 0's mass would otherwise read as carbon (12.011), tagged
        // nitrogen by ATOMIC_NUMBER: the column, not the mass table, wins.
        let text = Prmtop::default()
            .ints("POINTERS", 8, &[2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1])
            .strings("ATOM_NAME", 4, &["X", "Y"])
            .floats("CHARGE", 16, &[0.0, 0.0])
            .floats("MASS", 16, &[12.011, 1.008])
            .ints("ATOMIC_NUMBER", 8, &[7, 1])
            .strings("RESIDUE_LABEL", 4, &["LIG"])
            .ints("RESIDUE_POINTER", 8, &[1])
            .ints("BONDS_WITHOUT_HYDROGEN", 8, &[])
            .ints("BONDS_INC_HYDROGEN", 8, &[0, 3, 1])
            .build();
        let t = parse(&text).unwrap();
        assert_eq!(t.element[0], Element::NITROGEN);
        assert_eq!(t.element[1], Element::HYDROGEN);
    }

    #[test]
    fn empty_natom_is_no_atoms_error() {
        let text = Prmtop::default()
            .ints("POINTERS", 8, &[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
            .build();
        assert!(matches!(parse(&text), Err(ParseError::NoAtoms)));
    }

    #[test]
    fn declared_bonds_are_used_verbatim_even_when_geometry_would_disagree() {
        let t = parse(&sample()).unwrap();
        // Atoms 3 (O) and 4 (Na) sit at a clash distance the topology
        // does not bond; atoms 0-1 (N-CA) are declared but far apart.
        let p = [
            Vec3::ZERO,
            Vec3::new(20.0, 0.0, 0.0),
            Vec3::new(40.0, 0.0, 0.0),
            Vec3::new(60.0, 0.0, 0.0),
            Vec3::new(60.6, 0.0, 0.0),
        ];
        let bonds = vv_core::bonds::perceive(&t, &p);
        assert!(bonds.contains(0, 1));
        assert!(!bonds.contains(3, 4));
    }

    #[test]
    fn solute_molecules_each_get_a_chain_and_solvent_is_lumped_into_one() {
        // 5 single-atom "molecules": 2 solute, 3 solvent.
        let text = Prmtop::default()
            .ints("POINTERS", 8, &[5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5])
            .strings("ATOM_NAME", 4, &["A1", "A2", "A3", "A4", "A5"])
            .floats("CHARGE", 16, &[0.0; 5])
            .floats("MASS", 16, &[12.011, 12.011, 15.999, 15.999, 15.999])
            .strings("RESIDUE_LABEL", 4, &["R1", "R2", "R3", "R4", "R5"])
            .ints("RESIDUE_POINTER", 8, &[1, 2, 3, 4, 5])
            .ints("ATOMS_PER_MOLECULE", 8, &[1, 1, 1, 1, 1])
            .ints("SOLVENT_POINTERS", 8, &[2, 5, 3])
            .ints("BONDS_WITHOUT_HYDROGEN", 8, &[0, 3, 1])
            .ints("BONDS_INC_HYDROGEN", 8, &[])
            .build();
        let t = parse(&text).unwrap();
        assert_eq!(
            t.chain_count(),
            3,
            "2 solute chains + 1 lumped solvent chain"
        );
        assert_eq!(t.chain_name(0), "A");
        assert_eq!(t.chain_name(1), "B");
        assert_eq!(t.chain_name(2), "C");
        assert_eq!(
            t.chains[2].residues.len(),
            3,
            "all 3 solvent molecules share the trailing chain"
        );
    }

    #[test]
    fn without_atoms_per_molecule_everything_is_one_chain() {
        assert_eq!(parse(&sample()).unwrap().chain_count(), 1);
    }

    #[test]
    fn a_comment_line_between_flag_and_format_is_skipped() {
        // Newer AmberTools versions insert `%COMMENT` between `%FLAG` and
        // `%FORMAT`; the section reader must look past it for the format.
        let text = format!("%FLAG POINTERS\n%COMMENT metadata\n%FORMAT(I8)\n{:>8}\n", 5);
        assert_eq!(ints(&text, "POINTERS").unwrap(), vec![5]);
    }
}
