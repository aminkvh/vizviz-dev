//! Solvent-accessible surface area per atom, the way FastSASA computes it
//! (Shrake & Rupley 1973; FastSASA's defaults): each atom's sphere grown
//! by the probe radius carries test points on a golden-section spiral,
//! and the accessible fraction is the share of points no other grown
//! sphere covers. Radii are FastSASA's default: ProtOr (Tsai et al. 1999)
//! by residue and atom name, else by element.
//! `vv-io/tests/sasa.rs` checks it against FastSASA's own output.

use std::collections::HashMap;
use std::sync::OnceLock;

use glam::Vec3;
use rayon::prelude::*;

use crate::spatial::Grid;
use crate::topology::Topology;

/// FastSASA's default number of test points per atom.
pub const DEFAULT_POINTS: usize = 100;

/// FastSASA's radius for an element it knows (`fastsasa_element_radius`).
fn element_radius(symbol: &str) -> Option<f32> {
    Some(match symbol.to_ascii_uppercase().as_str() {
        "H" | "D" => 1.10,
        "C" => 1.70,
        "N" => 1.55,
        "O" => 1.52,
        "F" => 1.47,
        "P" => 1.80,
        "S" => 1.80,
        "CL" => 1.75,
        "BR" => 1.85,
        "I" => 1.98,
        "NA" => 2.27,
        "MG" => 1.73,
        "K" => 2.75,
        "CA" => 2.31,
        "FE" => 1.80,
        "ZN" => 1.39,
        _ => return None,
    })
}

/// ProtOr's radius per `"RESIDUE:ATOM"` (`protor.config`, FreeSASA's
/// format: a `types:` section of `TYPE RADIUS CLASS`, an `atoms:` section
/// of `RESIDUE ATOM TYPE`).
fn protor() -> &'static HashMap<String, f32> {
    static TABLE: OnceLock<HashMap<String, f32>> = OnceLock::new();
    TABLE.get_or_init(|| {
        let mut types = HashMap::new();
        let mut atoms = HashMap::new();
        let mut section = "";
        for line in include_str!("protor.config").lines() {
            let line = line.split('#').next().unwrap_or("").trim();
            if line.ends_with(':') {
                section = line;
                continue;
            }
            let words: Vec<&str> = line.split_whitespace().collect();
            match (section, words.as_slice()) {
                ("types:", [name, radius, ..]) => {
                    types.insert(*name, radius.parse::<f32>().expect("protor radius"));
                }
                ("atoms:", [residue, atom, kind]) => {
                    atoms.insert(format!("{residue}:{atom}"), types[kind]);
                }
                _ => {}
            }
        }
        atoms
    })
}

/// Each atom's radius for SASA, in FastSASA's order: ProtOr for its
/// residue and name, else for its name in any residue, else FastSASA's
/// radius for its element, else its van der Waals radius. Hydrogens have
/// no ProtOr entry and take the element radius.
pub fn radii(topology: &Topology) -> Vec<f32> {
    let table = protor();
    (0..topology.atom_count())
        .map(|atom| {
            let name = topology.atom_name(atom).to_ascii_uppercase();
            let residue = topology
                .residue_name(topology.residue_index[atom] as usize)
                .to_ascii_uppercase();
            let element = topology.element[atom];
            table
                .get(&format!("{residue}:{name}"))
                .or_else(|| table.get(&format!("ANY:{name}")))
                .copied()
                .or_else(|| element_radius(element.symbol()))
                .unwrap_or_else(|| element.vdw_radius())
        })
        .collect()
}

/// Whether each atom belongs to its residue's first conformer: atoms
/// without an alternate location, and those with the first one the
/// residue lists. FastSASA reads files this way, whatever it selects.
pub fn first_conformer(topology: &Topology) -> Vec<bool> {
    let mut chosen = vec![0u8; topology.residue_count()];
    (0..topology.atom_count())
        .map(|atom| {
            let alt = topology.alt_loc[atom];
            let residue = &mut chosen[topology.residue_index[atom] as usize];
            if alt == 0 {
                return true;
            }
            if *residue == 0 {
                *residue = alt;
            }
            *residue == alt
        })
        .collect()
}

/// The test points: `n` on the unit sphere along a golden-section spiral,
/// from the north pole down (FastSASA's `cpu_test_points`).
pub fn test_points(n: usize) -> Vec<Vec3> {
    let step = std::f64::consts::PI * (3.0 - 5f64.sqrt());
    let dz = 2.0 / n as f64;
    (0..n)
        .map(|k| {
            let z = 1.0 - dz / 2.0 - dz * k as f64;
            let r = (1.0 - z * z).sqrt();
            let longitude = step * k as f64;
            Vec3::new(
                (longitude.cos() * r) as f32,
                (longitude.sin() * r) as f32,
                z as f32,
            )
        })
        .collect()
}

/// Accessible area of each atom in Å², for spheres at `positions` with
/// `radii`, a solvent probe of `probe` Å and `points` test points per
/// atom. Only these atoms occlude each other: pass the ones that count
/// (FastSASA leaves out hydrogens and HETATM records by default).
pub fn shrake_rupley(positions: &[Vec3], radii: &[f32], probe: f32, points: usize) -> Vec<f32> {
    assert_eq!(positions.len(), radii.len());
    if positions.is_empty() {
        return Vec::new();
    }
    let grown: Vec<f32> = radii.iter().map(|r| r + probe).collect();
    let largest = grown.iter().copied().fold(0.0, f32::max);
    let indices: Vec<u32> = (0..positions.len() as u32).collect();
    let grid = Grid::build(positions, &indices, 2.0 * largest);
    let sphere = test_points(points);
    (0..positions.len())
        .into_par_iter()
        .map_init(Vec::new, |near: &mut Vec<u32>, atom| {
            let (c, r) = (positions[atom], grown[atom]);
            near.clear();
            grid.for_each_within(positions, c, r + largest, |other, d2| {
                let reach = r + grown[other as usize];
                if other as usize != atom && d2 < reach * reach {
                    near.push(other);
                }
            });
            let exposed = sphere
                .iter()
                .filter(|&&u| {
                    let p = c + u * r;
                    !near.iter().any(|&o| {
                        let g = grown[o as usize];
                        p.distance_squared(positions[o as usize]) < g * g
                    })
                })
                .count();
            4.0 * std::f32::consts::PI * r * r * exposed as f32 / points as f32
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_lone_atom_is_all_exposed_and_a_buried_one_not() {
        let alone = shrake_rupley(&[Vec3::ZERO], &[1.7], 1.4, DEFAULT_POINTS);
        let full = 4.0 * std::f32::consts::PI * 3.1 * 3.1;
        assert!((alone[0] - full).abs() < 1e-3);
        // An atom at the centre of a tight shell of six.
        let mut positions = vec![Vec3::ZERO];
        positions.extend([
            Vec3::X,
            Vec3::NEG_X,
            Vec3::Y,
            Vec3::NEG_Y,
            Vec3::Z,
            Vec3::NEG_Z,
        ]);
        let area = shrake_rupley(&positions, &[1.7; 7], 1.4, DEFAULT_POINTS);
        assert_eq!(area[0], 0.0);
    }

    #[test]
    fn protor_types_resolve_and_hydrogens_are_left_to_elements() {
        let table = protor();
        assert_eq!(table["ALA:CB"], 1.88);
        assert_eq!(table["ALA:C"], 1.61);
        assert_eq!(table["HOH:O"], 1.46);
        assert!(!table.contains_key("ALA:H"));
    }

    #[test]
    fn test_points_lie_on_the_unit_sphere_and_balance() {
        let points = test_points(DEFAULT_POINTS);
        assert!(points.iter().all(|p| (p.length() - 1.0).abs() < 1e-5));
        let mean = points.iter().copied().sum::<Vec3>() / points.len() as f32;
        assert!(mean.length() < 0.02);
    }
}
