//! Parser checks against real PDB entries in `fixtures/small/` and against
//! the synthetic generator through the writer.

use std::path::PathBuf;
use std::time::{Duration, Instant};

use vv_core::bonds;
use vv_core::glam::Vec3;
use vv_core::{Element, SecondaryStructure, Structure};
use vv_io::synth::{protein_like, SynthParams};

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../fixtures/small")
        .join(name)
}

fn load(name: &str) -> Structure {
    let s = vv_io::load(fixture(name)).unwrap_or_else(|e| panic!("{name}: {e}"));
    s.topology.validate().unwrap();
    s
}

#[test]
fn crambin_mmcif_hierarchy() {
    let s = load("1CRN.cif");
    let t = &s.topology;
    assert_eq!(t.id, "1CRN");
    assert_eq!(t.atom_count(), 327);
    assert_eq!(t.residue_count(), 46);
    assert_eq!(t.chain_count(), 1);
    assert_eq!(t.chain_name(0), "A");
    assert_eq!(t.residue_name(0), "THR");
    assert_eq!(t.atom_name(0), "N");
    assert_eq!(t.element[0], Element::NITROGEN);
    assert_eq!(s.frame(0).positions()[0].x, 17.047);
    assert_eq!(t.residue_name(45), "ASN");
    let helices = t
        .residues
        .iter()
        .filter(|r| r.ss == SecondaryStructure::Helix)
        .count();
    let strands = t
        .residues
        .iter()
        .filter(|r| r.ss == SecondaryStructure::Strand)
        .count();
    assert!(helices >= 15, "helix residues: {helices}");
    assert!(strands >= 4, "strand residues: {strands}");
    // Three disulfides declared in _struct_conn.
    assert_eq!(t.explicit_bonds.len(), 3);
    for b in &t.explicit_bonds {
        assert!(b
            .atoms
            .iter()
            .all(|&a| t.element[a as usize] == Element::SULFUR));
    }
}

#[test]
fn crambin_pdb_matches_mmcif() {
    let cif = load("1CRN.cif");
    let pdb = load("1CRN.pdb");
    assert_eq!(pdb.atom_count(), cif.atom_count());
    assert_eq!(pdb.topology.residue_count(), cif.topology.residue_count());
    for a in 0..cif.atom_count() {
        assert_eq!(pdb.topology.element[a], cif.topology.element[a], "atom {a}");
        assert_eq!(
            pdb.topology.atom_name(a),
            cif.topology.atom_name(a),
            "atom {a}"
        );
        assert_eq!(
            pdb.frame(0).positions()[a],
            cif.frame(0).positions()[a],
            "atom {a}"
        );
        assert_eq!(pdb.topology.serial[a], cif.topology.serial[a]);
    }
    assert_eq!(
        pdb.topology
            .residues
            .iter()
            .map(|r| r.ss)
            .collect::<Vec<_>>(),
        cif.topology
            .residues
            .iter()
            .map(|r| r.ss)
            .collect::<Vec<_>>()
    );
}

#[test]
fn crambin_bond_perception() {
    let s = load("1CRN.cif");
    let t = &s.topology;
    let b = bonds::perceive(t, s.frame(0).positions());
    // Every residue has N-CA, CA-C, C-O; consecutive residues are peptide bonded.
    let find = |r: usize, name: &str| {
        t.residues[r]
            .atoms
            .clone()
            .find(|&a| t.atom_name(a as usize) == name)
            .unwrap()
    };
    for r in 0..t.residue_count() {
        assert!(
            b.contains(find(r, "N"), find(r, "CA")),
            "N-CA in residue {r}"
        );
        assert!(
            b.contains(find(r, "CA"), find(r, "C")),
            "CA-C in residue {r}"
        );
        assert!(b.contains(find(r, "C"), find(r, "O")), "C-O in residue {r}");
        if r + 1 < t.residue_count() {
            assert!(
                b.contains(find(r, "C"), find(r + 1, "N")),
                "peptide bond {r}"
            );
        }
    }
    let sulfurs: Vec<u32> = (0..t.atom_count() as u32)
        .filter(|&a| t.element[a as usize] == Element::SULFUR)
        .collect();
    let bonds = &b;
    let ss_bonds = sulfurs
        .iter()
        .flat_map(|&a| {
            sulfurs
                .iter()
                .filter(move |&&c| c > a && bonds.contains(a, c))
        })
        .count();
    assert_eq!(ss_bonds, 3, "disulfide bonds");
    let degrees = b.degrees(t.atom_count());
    assert!(
        degrees.iter().all(|&d| d >= 1),
        "every heavy atom in crambin is bonded"
    );
    assert!(
        degrees.iter().all(|&d| d <= 4),
        "no atom has more than 4 bonds"
    );
    assert!(
        b.len() > t.atom_count(),
        "proteins have slightly more bonds than atoms"
    );
}

#[test]
fn hemoglobin_chains_and_heme_iron() {
    let s = load("4HHB.cif");
    let t = &s.topology;
    assert_eq!(t.atom_count(), 4779);
    let polymer_chains = t.chains.iter().filter(|c| {
        let first = &t.residues[c.residues.start as usize];
        t.flags[first.atoms.start as usize] & vv_core::flags::HETERO == 0
    });
    assert_eq!(polymer_chains.count(), 4);
    let b = bonds::perceive(t, s.frame(0).positions());
    let irons: Vec<u32> = (0..t.atom_count() as u32)
        .filter(|&a| t.element[a as usize] == Element::from_atomic_number(26).unwrap())
        .collect();
    assert_eq!(irons.len(), 4);
    let degrees = b.degrees(t.atom_count());
    for fe in irons {
        // Four porphyrin nitrogens (the HEM template's Fe-N) plus the
        // proximal histidine NE2 (an explicit `_struct_conn` metalc bond).
        assert_eq!(degrees[fe as usize], 5, "Fe {fe}");
    }
}

const RING_COUNT: &[(&str, usize)] = &[("PRO", 1), ("PHE", 1), ("TYR", 1), ("HIS", 1), ("TRP", 2)];
const STANDARD_AMINO_ACIDS: &[&str] = &[
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET",
    "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
];

/// Bond count for one amino-acid residue's heavy atoms, computed without
/// any reference to `vv_core::bonds` or the CCD templates it uses: a
/// residue's heavy atoms form a spanning tree (`heavy - 1` bonds) plus
/// one extra bond per ring (its cyclomatic number). Counts only the
/// first alternate location, so altloc'd residues aren't double-counted.
fn independent_residue_bond_count(t: &vv_core::Topology, residue: usize) -> usize {
    let res = &t.residues[residue];
    let mut heavy = 0usize;
    let mut first_alt = None;
    for a in res.atoms.clone() {
        let alt = t.alt_loc[a as usize];
        if alt != 0 {
            if first_alt.is_none() {
                first_alt = Some(alt);
            }
            if Some(alt) != first_alt {
                continue;
            }
        }
        heavy += 1;
    }
    let rings = RING_COUNT
        .iter()
        .find(|(name, _)| *name == t.residue_name(residue))
        .map_or(0, |(_, n)| *n);
    heavy.saturating_sub(1) + rings
}

/// Every standard-residue chain's bond count against a from-scratch
/// count independent of `vv_core::bonds` and its CCD templates: the sum
/// of each residue's spanning-tree-plus-rings count, plus one peptide
/// bond per adjacent pair, must equal exactly the perceived bonds among
/// those residues' heavy atoms (no more, no fewer -- this also checks
/// that no spurious extra bond appears between non-adjacent residues).
fn assert_independent_protein_bond_count(name: &str) {
    let s = load(name);
    let t = &s.topology;
    let b = bonds::perceive(t, s.frame(0).positions());
    let is_standard = |a: u32| {
        STANDARD_AMINO_ACIDS.contains(&t.residue_name(t.residue_index[a as usize] as usize))
    };

    let mut expected = 0usize;
    for chain in &t.chains {
        let residues: Vec<usize> = chain
            .residues
            .clone()
            .map(|r| r as usize)
            .filter(|&r| STANDARD_AMINO_ACIDS.contains(&t.residue_name(r)))
            .collect();
        for &r in &residues {
            expected += independent_residue_bond_count(t, r);
        }
        expected += residues.len().saturating_sub(1);
    }
    // Disulfides are neither a residue-internal bond nor an adjacent-pair
    // peptide bond, so the spanning-tree sum above doesn't see them.
    let disulfides = b
        .pairs
        .iter()
        .filter(|&&[a, b]| t.atom_name(a as usize) == "SG" && t.atom_name(b as usize) == "SG")
        .count();
    expected += disulfides;

    let actual = b
        .pairs
        .iter()
        .filter(|&&[a, b]| is_standard(a) && is_standard(b))
        .count();
    // 4HHB (a 1984 deposit) has a handful of exposed, genuinely disordered
    // side-chain tips (Lys NZ, Glu OE2, Asp OD2) modeled 2.5-3.1 A from
    // their spanning-tree parent -- past any reasonable bond distance, so
    // `perceive` correctly leaves them unbonded rather than drawing a
    // visibly wrong long bond. `actual` must never *exceed* the ideal
    // count (that would be a spurious extra bond), and the shortfall from
    // real disorder must stay small.
    assert!(
        actual <= expected,
        "{name}: {actual} > {expected} (extra bond)"
    );
    let shortfall = expected - actual;
    assert!(
        (shortfall as f32) <= 0.02 * expected as f32,
        "{name}: {shortfall}/{expected} bonds missing, more than expected disorder"
    );
}

#[test]
fn independent_bond_count_matches_on_every_protein_fixture() {
    for name in ["1CRN.cif", "1UBQ.cif", "4HHB.cif"] {
        assert_independent_protein_bond_count(name);
    }
}

/// No extra inter-residue bonds: the number of C(i)-N(i+1) bonds between
/// different residues equals exactly the number of adjacent standard-
/// residue pairs within [`vv_core::bonds`]'s peptide-bond window, on
/// every protein fixture -- not just that every real peptide bond is
/// present (`crambin_bond_perception` already checks that for 1CRN) but
/// that nothing else got connected the same way.
#[test]
fn peptide_bonds_have_no_extras() {
    for name in ["1CRN.cif", "1UBQ.cif", "4HHB.cif"] {
        let s = load(name);
        let t = &s.topology;
        let coords = s.frame(0);
        let p = coords.positions();
        let b = bonds::perceive(t, p);

        let find = |r: usize, atom_name: &str| {
            t.residues[r]
                .atoms
                .clone()
                .find(|&a| t.atom_name(a as usize) == atom_name)
        };
        let mut expected_links = 0;
        for chain in &t.chains {
            for r in chain.residues.start..chain.residues.end.saturating_sub(1) {
                let (r, next) = (r as usize, r as usize + 1);
                if let (Some(c), Some(n)) = (find(r, "C"), find(next, "N")) {
                    if p[c as usize].distance(p[n as usize]) <= 1.6 {
                        expected_links += 1;
                    }
                }
            }
        }

        let inter_residue_cn = b
            .pairs
            .iter()
            .filter(|&&[a, b]| {
                t.atom_name(a as usize) == "C"
                    && t.atom_name(b as usize) == "N"
                    && t.residue_index[a as usize] != t.residue_index[b as usize]
            })
            .count();
        assert_eq!(inter_residue_cn, expected_links, "{name}");
    }
}

#[test]
fn ubiquitin_and_adenylate_kinase_load() {
    let u = load("1UBQ.cif");
    assert_eq!(u.atom_count(), 660);
    assert_eq!(u.topology.chain_count(), 2); // protein + water
    let a = load("1AKE.cif");
    assert_eq!(
        a.topology.chains.iter().filter(|c| c.entity == 1).count(),
        2
    );
    assert!(a.atom_count() > 3000);
    assert!(a.topology.title.to_ascii_lowercase().contains("adenylate"));
}

#[test]
fn pdb_and_mmcif_agree_for_every_fixture() {
    for id in ["1UBQ", "4HHB", "1AKE"] {
        let cif = load(&format!("{id}.cif"));
        let pdb = load(&format!("{id}.pdb"));
        assert_eq!(pdb.atom_count(), cif.atom_count(), "{id}");
        for a in 0..cif.atom_count() {
            assert_eq!(
                pdb.frame(0).positions()[a],
                cif.frame(0).positions()[a],
                "{id} atom {a}"
            );
            assert_eq!(
                pdb.topology.element[a], cif.topology.element[a],
                "{id} atom {a}"
            );
        }
    }
}

#[test]
fn synthetic_structure_round_trips_through_the_writer() {
    let original = protein_like(&SynthParams {
        residues_per_chain: 40,
        ..SynthParams::new(20_000)
    });
    let mut bytes = Vec::new();
    vv_io::mmcif_write::write(&original, None, &[0], &mut bytes).unwrap();
    let parsed = vv_io::parse(&bytes, vv_io::Format::Mmcif).unwrap();
    parsed.topology.validate().unwrap();
    let (o, p) = (&original.topology, &parsed.topology);
    assert_eq!(p.atom_count(), o.atom_count());
    assert_eq!(p.residue_count(), o.residue_count());
    assert_eq!(p.chain_count(), o.chain_count());
    assert_eq!(p.element, o.element);
    assert_eq!(p.residue_index, o.residue_index);
    for (a, b) in p.residues.iter().zip(&o.residues) {
        assert_eq!(a.atoms, b.atoms);
        assert_eq!(a.seq_id, b.seq_id);
        assert_eq!(p.names.get(a.comp), o.names.get(b.comp));
    }
    for (a, b) in p.chains.iter().zip(&o.chains) {
        assert_eq!(a.residues, b.residues);
        assert_eq!(p.names.get(a.label_asym), o.names.get(b.label_asym));
    }
    for (a, b) in parsed
        .frame(0)
        .positions()
        .iter()
        .zip(original.frame(0).positions())
    {
        assert!((*a - *b).abs().max_element() <= 0.0006, "{a} vs {b}");
    }
}

/// Runs `perceive` 5 times and returns the minimum duration and the bond
/// count (constant across runs) -- min-of-5 is far less sensitive than a
/// mean to a transient scheduling hiccup from another process on a shared
/// machine.
fn min_of_5(topology: &vv_core::Topology, positions: &[Vec3]) -> (Duration, usize) {
    let mut best = Duration::MAX;
    let mut count = 0;
    for _ in 0..5 {
        let t0 = Instant::now();
        let b = bonds::perceive(topology, positions);
        best = best.min(t0.elapsed());
        count = b.len();
    }
    (best, count)
}

/// Bond perception timing on real structures and a 1M-atom synthetic
/// one. Numbers are recorded in docs/VALIDATION.md. Ignored by default;
/// run with `cargo test -p vv-io --test fixtures --release -- --ignored
/// --nocapture bond_perception_timing`.
#[test]
#[ignore]
fn bond_perception_timing() {
    // Warms up the process-wide template index (built once, lazily, on
    // the first `perceive` call ever) so it isn't charged to 4HHB below.
    let warm = load("1CRN.cif");
    bonds::perceive(&warm.topology, warm.frame(0).positions());

    let s = load("4HHB.cif");
    let (dt, n) = min_of_5(&s.topology, s.frame(0).positions());
    eprintln!(
        "4HHB: {} atoms, {n} bonds, min of 5 = {dt:.2?}",
        s.atom_count()
    );

    // 6X3Z: 17K atoms, real chemistry, almost entirely templated
    // (protein + water + glycans) -- the representative real-structure-
    // at-scale number, unlike the synthetic system below.
    let g = vv_io::load(
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/glycan/6X3Z.pdb"),
    )
    .unwrap();
    let (dt, n) = min_of_5(&g.topology, g.frame(0).positions());
    eprintln!(
        "6X3Z: {} atoms, {n} bonds, min of 5 = {dt:.2?}",
        g.atom_count()
    );

    let big = protein_like(&SynthParams::new(1_000_000));
    let (dt, n) = min_of_5(&big.topology, big.frame(0).positions());
    eprintln!(
        "synthetic 1M: {} atoms, {n} bonds, min of 5 = {dt:.2?}",
        big.atom_count()
    );
}

#[test]
fn gzip_input_is_transparent() {
    let raw = std::fs::read(fixture("1CRN.cif")).unwrap();
    let mut gz = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
    std::io::Write::write_all(&mut gz, &raw).unwrap();
    let compressed = gz.finish().unwrap();
    let dir = std::env::temp_dir().join("vizviz-tests");
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("1CRN.cif.gz");
    std::fs::write(&path, compressed).unwrap();
    let s = vv_io::load(&path).unwrap();
    assert_eq!(s.atom_count(), 327);
}
