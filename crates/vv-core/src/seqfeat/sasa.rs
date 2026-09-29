//! Per-residue relative solvent accessibility: burial for the sequence
//! strip's coloring.

use glam::Vec3;

use crate::{sasa, ResidueClass, Topology};

/// Theoretical maximum residue accessible area (Å²), extended Gly-X-Gly
/// conformations: Tien et al. 2013, PLoS ONE 8:e80635, Table 1.
fn max_area(name: &str) -> Option<f32> {
    Some(match name {
        "ALA" => 129.0,
        "ARG" => 274.0,
        "ASN" => 195.0,
        "ASP" => 193.0,
        "CYS" => 167.0,
        "GLN" => 225.0,
        "GLU" => 223.0,
        "GLY" => 104.0,
        "HIS" => 224.0,
        "ILE" => 197.0,
        "LEU" => 201.0,
        "LYS" => 236.0,
        "MET" | "MSE" => 224.0,
        "PHE" => 240.0,
        "PRO" => 159.0,
        "SER" => 155.0,
        "THR" => 172.0,
        "TRP" => 285.0,
        "TYR" => 263.0,
        "VAL" => 174.0,
        _ => return None,
    })
}

/// Accessible area of each residue as a fraction of its maximum (1 =
/// fully exposed), `NaN` for residues without a reference maximum. The
/// polymer heavy atoms (first conformer) occlude each other; ligands and
/// water are left out, as for a bare-protein surface.
pub fn relative_sasa(top: &Topology, positions: &[Vec3], probe: f32) -> Vec<f32> {
    let conformer = sasa::first_conformer(top);
    let radii = sasa::radii(top);
    let atoms: Vec<usize> = (0..top.atom_count())
        .filter(|&a| conformer[a] && !top.element[a].is_hydrogen())
        .filter(|&a| in_polymer(top, top.residue_index[a] as usize))
        .collect();
    let coords: Vec<Vec3> = atoms.iter().map(|&a| positions[a]).collect();
    let sphere: Vec<f32> = atoms.iter().map(|&a| radii[a]).collect();
    let area = sasa::shrake_rupley(&coords, &sphere, probe, sasa::DEFAULT_POINTS);
    let mut per_residue = vec![0.0f32; top.residue_count()];
    for (&a, &v) in atoms.iter().zip(&area) {
        per_residue[top.residue_index[a] as usize] += v;
    }
    per_residue
        .iter()
        .enumerate()
        .map(|(r, &v)| max_area(top.residue_name(r)).map_or(f32::NAN, |m| v / m))
        .collect()
}

fn in_polymer(top: &Topology, residue: usize) -> bool {
    matches!(
        top.residue_class(residue),
        ResidueClass::Protein | ResidueClass::Nucleic
    )
}
