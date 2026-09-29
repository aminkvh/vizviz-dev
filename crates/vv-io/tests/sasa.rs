//! SASA against FastSASA. Each `*.sasa.pdb` is FastSASA's own output
//! (`fastsasa --cpu --config-file share/protor.config --format pdb`:
//! Shrake-Rupley, 100 points, 1.4 A probe, ProtOr radii), holding the
//! atoms it computed, each radius as occupancy and SASA as B-factor:
//! - `1CRN.sasa.pdb`: the defaults (no hydrogens, no HETATM);
//! - `1AKE.sasa.pdb`: the defaults, on alternate locations;
//! - `1D3Z.H.sasa.pdb`: `--hydrogen`, on the first NMR model of 1D3Z;
//! - `1UBQ.het.sasa.pdb`: `--hetatm`, waters included.

use std::path::PathBuf;

fn small(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../fixtures/small")
        .join(name)
}

/// Computes SASA over the atoms `filter` selects in `structure` and
/// checks radii, per-atom areas and the total against `reference`.
fn matches_fastsasa(structure: &str, filter: &str, reference: &str) {
    let structure = vv_io::load(small(structure)).unwrap();
    let text = std::fs::read_to_string(small(reference)).unwrap();
    let (radius_ref, sasa_ref): (Vec<f32>, Vec<f32>) = text
        .lines()
        .filter(|l| l.starts_with("ATOM") || l.starts_with("HETATM"))
        .map(|l| {
            let f = |a: usize, b: usize| l[a..b].trim().parse::<f32>().unwrap();
            (f(54, 60), f(60, 66))
        })
        .unzip();
    let coords = structure.frame(0);
    let all = coords.positions();
    let conformer = vv_core::sasa::first_conformer(&structure.topology);
    let atoms: Vec<usize> = vv_core::select(&structure.topology, all, filter)
        .unwrap()
        .ones()
        .filter(|&a| conformer[a])
        .collect();
    assert_eq!(atoms.len(), radius_ref.len(), "{reference}: atom count");

    let every = vv_core::sasa::radii(&structure.topology);
    let radii: Vec<f32> = atoms.iter().map(|&a| every[a]).collect();
    for (a, (&ours, &theirs)) in radii.iter().zip(&radius_ref).enumerate() {
        assert!(
            (ours - theirs).abs() < 1e-3,
            "{reference} atom {a}: radius {ours} vs {theirs}"
        );
    }
    let positions: Vec<_> = atoms.iter().map(|&a| all[a]).collect();
    let sasa = vv_core::sasa::shrake_rupley(&positions, &radii, 1.4, 100);
    let (ours, theirs): (f32, f32) = (sasa.iter().sum(), sasa_ref.iter().sum());
    assert!(
        (ours - theirs).abs() < 0.001 * theirs,
        "{reference}: total {ours} vs {theirs}"
    );
    // f32 against FastSASA's f64: a test point exactly on a sphere may
    // flip, one point in a hundred of that atom's area.
    for (a, (&o, &t)) in sasa.iter().zip(&sasa_ref).enumerate() {
        let one_point = 4.0 * std::f32::consts::PI * (radii[a] + 1.4).powi(2) / 100.0;
        assert!(
            (o - t).abs() <= one_point + 0.01,
            "{reference} atom {a}: {o} vs {t}"
        );
    }
}

#[test]
fn crambin_matches_fastsasa_defaults() {
    matches_fastsasa("1CRN.pdb", "not hydrogen and not hetero", "1CRN.sasa.pdb");
}

#[test]
fn only_the_first_conformer_counts() {
    matches_fastsasa("1AKE.pdb", "not hydrogen and not hetero", "1AKE.sasa.pdb");
}

#[test]
fn hydrogens_match_fastsasa_hydrogen() {
    matches_fastsasa("1D3Z.pdb", "all", "1D3Z.H.sasa.pdb");
}

#[test]
fn waters_match_fastsasa_hetatm() {
    matches_fastsasa("1UBQ.pdb", "all", "1UBQ.het.sasa.pdb");
}
