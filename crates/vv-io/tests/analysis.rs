//! Analysis kernels against real structures: the dihedral sign pinned to
//! known secondary structure, contacts pinned to chemistry, and timings at
//! the scale the kernels are designed for (ignored by default; run with
//! `cargo test -p vv-io --test analysis --release -- --ignored --nocapture`).

use std::time::Instant;

use vv_core::analysis::{
    angle, contact_counts_frames, contacts_into, dihedral, dihedrals_frames, neighbor_pairs_into,
    residue_pairs,
};
use vv_core::glam::Vec3;
use vv_core::{select, Structure};

fn fixture(name: &str) -> Structure {
    let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../fixtures/small")
        .join(name);
    vv_io::load(path).unwrap()
}

/// Atom index of `name` in residue `r` (0-based residue index).
fn atom(s: &Structure, r: usize, name: &str) -> usize {
    let t = &s.topology;
    t.residues[r]
        .atoms
        .clone()
        .map(|a| a as usize)
        .find(|&a| t.atom_name(a) == name)
        .unwrap_or_else(|| panic!("residue {r} has no atom {name}"))
}

fn phi_psi(s: &Structure, r: usize) -> (f32, f32) {
    let coords = s.frame(0);
    let p = coords.positions();
    let phi = dihedral(
        p,
        atom(s, r - 1, "C"),
        atom(s, r, "N"),
        atom(s, r, "CA"),
        atom(s, r, "C"),
    );
    let psi = dihedral(
        p,
        atom(s, r, "N"),
        atom(s, r, "CA"),
        atom(s, r, "C"),
        atom(s, r + 1, "N"),
    );
    (phi, psi)
}

/// Crambin residues 7-19 are an alpha helix. A right-handed helix has phi
/// around -60 and psi around -45; a sign error would put both positive.
#[test]
fn crambin_helix_has_negative_phi_and_psi() {
    let s = fixture("1CRN.cif");
    for r in 8..=17 {
        let (phi, psi) = phi_psi(&s, r);
        assert!((-100.0..=-40.0).contains(&phi), "residue {r} phi {phi}");
        assert!((-70.0..=-10.0).contains(&psi), "residue {r} psi {psi}");
    }
    // Backbone bond angles are chemistry, not geometry of our choosing:
    // N-CA-C is ~111 degrees in every protein.
    let coords = s.frame(0);
    let p = coords.positions();
    for r in 1..40 {
        let tau = angle(p, atom(&s, r, "N"), atom(&s, r, "CA"), atom(&s, r, "C"));
        assert!((100.0..=122.0).contains(&tau), "residue {r} tau {tau}");
    }
}

#[test]
fn hemoglobin_interface_contacts_are_between_the_right_chains() {
    let s = fixture("4HHB.cif");
    let t = &s.topology;
    let coords = s.frame(0);
    let p = coords.positions();
    let ones = |expr: &str| -> Vec<u32> {
        select(t, p, expr)
            .unwrap()
            .ones()
            .map(|i| i as u32)
            .collect()
    };
    let (a, b) = (ones("chain A and protein"), ones("chain B and protein"));
    let mut contacts = Vec::new();
    contacts_into(p, &a, &b, 4.0, &mut contacts);
    // Brute force is the reference: every A/B pair within 4 A, once.
    let expected = a
        .iter()
        .flat_map(|&i| b.iter().map(move |&j| (i, j)))
        .filter(|&(i, j)| p[i as usize].distance(p[j as usize]) <= 4.0)
        .count();
    assert_eq!(contacts.len(), expected);
    assert!(contacts.len() > 50, "the alpha1/beta1 interface is real");
    assert!(contacts.iter().all(|c| c.distance <= 4.0));
    assert!(contacts
        .iter()
        .all(|c| t.chain_name(t.chain_of_atom(c.a as usize) as usize) == "A"));
    assert!(contacts
        .iter()
        .all(|c| t.chain_name(t.chain_of_atom(c.b as usize) as usize) == "B"));
    let residues = residue_pairs(t, &contacts);
    // Reference count from the brute-force pairs, collapsed to residues.
    let mut expected_residues: Vec<[u32; 2]> = a
        .iter()
        .flat_map(|&i| b.iter().map(move |&j| (i, j)))
        .filter(|&(i, j)| p[i as usize].distance(p[j as usize]) <= 4.0)
        .map(|(i, j)| {
            let (ri, rj) = (t.residue_index[i as usize], t.residue_index[j as usize]);
            [ri.min(rj), ri.max(rj)]
        })
        .collect();
    expected_residues.sort_unstable();
    expected_residues.dedup();
    assert_eq!(residues, expected_residues);
    assert!(residues.len() > 20 && residues.len() < contacts.len());

    // Heme iron coordinates the proximal histidine NE2 at about 2.1 A.
    let fe = ones("element FE");
    let ne2 = ones("resname HIS and name NE2");
    let mut coordination = Vec::new();
    contacts_into(p, &fe, &ne2, 2.6, &mut coordination);
    assert_eq!(coordination.len(), 4, "one proximal His per heme");

    // Neighbor pairs within 1.9 A of protein heavy atoms are covalent bonds.
    let protein = ones("protein and not hydrogen");
    let mut pairs = Vec::new();
    neighbor_pairs_into(p, &protein, 1.9, &mut pairs);
    let bonds = vv_core::bonds::perceive(t, p);
    let bonded = pairs.iter().filter(|c| bonds.contains(c.a, c.b)).count();
    assert!(
        bonded as f32 > 0.95 * pairs.len() as f32,
        "{bonded}/{}",
        pairs.len()
    );
}

/// Timing at design scale. Numbers are recorded in docs/ANALYSIS.md.
#[test]
#[ignore]
fn timings_at_scale() {
    let s = fixture("4HHB.cif");
    let t = &s.topology;
    let coords = s.frame(0);
    let base = coords.positions();
    // 40,000 "frames": the same coordinates jittered per frame.
    let frames_owned: Vec<Vec<Vec3>> = (0..40_000)
        .map(|f| {
            let d = (f % 7) as f32 * 0.01;
            base.iter().map(|p| *p + Vec3::splat(d)).collect()
        })
        .collect();
    let frames: Vec<&[Vec3]> = frames_owned.iter().map(|v| v.as_slice()).collect();
    let n = t.residue_count();
    let quads: Vec<[u32; 4]> = (1..n - 1)
        .filter(|&r| {
            ["N", "CA", "C"].iter().all(|nm| {
                t.residues[r]
                    .atoms
                    .clone()
                    .any(|a| t.atom_name(a as usize) == *nm)
            }) && t.residues[r - 1]
                .atoms
                .clone()
                .any(|a| t.atom_name(a as usize) == "C")
        })
        .map(|r| {
            [
                atom(&s, r - 1, "C") as u32,
                atom(&s, r, "N") as u32,
                atom(&s, r, "CA") as u32,
                atom(&s, r, "C") as u32,
            ]
        })
        .collect();
    let mut out = vec![0.0f32; frames.len() * quads.len()];
    let t0 = Instant::now();
    dihedrals_frames(&frames, &quads, &mut out);
    let dt = t0.elapsed();
    eprintln!(
        "dihedrals: {} frames x {} phi = {} values in {:.0} ms ({:.1} M/s)",
        frames.len(),
        quads.len(),
        out.len(),
        dt.as_secs_f64() * 1e3,
        out.len() as f64 / dt.as_secs_f64() / 1e6
    );

    let p = base;
    let ones = |expr: &str| -> Vec<u32> {
        select(t, p, expr)
            .unwrap()
            .ones()
            .map(|i| i as u32)
            .collect()
    };
    let (a, b) = (ones("chain A"), ones("chain B"));
    let t0 = Instant::now();
    let counts = contact_counts_frames(&frames, &a, &b, 4.0);
    let dt = t0.elapsed();
    eprintln!(
        "contact counts: {} frames, chain A ({}) vs B ({}) at 4 A: {} contacts/frame in {:.0} ms ({:.1} us/frame)",
        frames.len(),
        a.len(),
        b.len(),
        counts[0],
        dt.as_secs_f64() * 1e3,
        dt.as_secs_f64() * 1e6 / frames.len() as f64
    );

    let all: Vec<u32> = (0..t.atom_count() as u32).collect();
    let mut pairs = Vec::new();
    let t0 = Instant::now();
    neighbor_pairs_into(p, &all, 4.0, &mut pairs);
    eprintln!(
        "neighbor list: {} atoms at 4 A: {} pairs in {:.1} ms",
        all.len(),
        pairs.len(),
        t0.elapsed().as_secs_f64() * 1e3
    );

    let large = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../fixtures/large/3J3Q.cif.gz");
    if large.exists() {
        let big = vv_io::load(large).unwrap();
        let coords = big.frame(0);
        let p = coords.positions();
        let all: Vec<u32> = (0..big.atom_count() as u32).collect();
        let mut pairs = Vec::new();
        let t0 = Instant::now();
        neighbor_pairs_into(p, &all, 4.0, &mut pairs);
        eprintln!(
            "neighbor list: 3J3Q {} atoms at 4 A: {} pairs in {:.0} ms",
            all.len(),
            pairs.len(),
            t0.elapsed().as_secs_f64() * 1e3
        );
    }
}
