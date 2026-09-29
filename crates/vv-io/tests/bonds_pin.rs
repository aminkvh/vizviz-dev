//! Pins `vv_core::bonds::perceive`'s output on real structures and a
//! mid-size synthetic system, so a perception-performance change can be
//! checked for zero effect on what actually gets bonded. Not a science
//! check (docs/VALIDATION.md's fixtures.rs/glycan.rs own that) -- a
//! byte-for-byte regression pin.

use std::path::PathBuf;

use vv_core::bonds::{self, BondTable};
use vv_io::synth::{protein_like, SynthParams};

fn small_fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../fixtures/small")
        .join(name)
}

fn glycan_fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../fixtures/glycan")
        .join(name)
}

/// FNV-1a over the sorted `[a, b]` pairs' little-endian bytes: cheap,
/// deterministic, and any change in which pairs bond changes it.
fn hash_bonds(bonds: &BondTable) -> u64 {
    let mut h: u64 = 0xcbf29ce484222325;
    for [a, b] in &bonds.pairs {
        for byte in a.to_le_bytes().into_iter().chain(b.to_le_bytes()) {
            h ^= byte as u64;
            h = h.wrapping_mul(0x100000001b3);
        }
    }
    h
}

struct Pin {
    count: usize,
    hash: u64,
    first: [u32; 2],
    last: [u32; 2],
}

fn pin_of(bonds: &BondTable) -> Pin {
    Pin {
        count: bonds.len(),
        hash: hash_bonds(bonds),
        first: *bonds.pairs.first().expect("non-empty"),
        last: *bonds.pairs.last().expect("non-empty"),
    }
}

fn assert_pin(bonds: &BondTable, expected: Pin, label: &str) {
    let actual = pin_of(bonds);
    assert_eq!(actual.count, expected.count, "{label}: bond count");
    assert_eq!(actual.first, expected.first, "{label}: first pair");
    assert_eq!(actual.last, expected.last, "{label}: last pair");
    assert_eq!(actual.hash, expected.hash, "{label}: bond set hash");
}

#[test]
fn crambin_1crn_bonds_are_pinned() {
    let s = vv_io::load(small_fixture("1CRN.cif")).unwrap();
    let bonds = bonds::perceive(&s.topology, s.frame(0).positions());
    assert_pin(
        &bonds,
        Pin {
            count: 337,
            hash: 0x690f4b1089a7e082,
            first: [0, 1],
            last: [323, 325],
        },
        "1CRN",
    );
}

#[test]
fn ubiquitin_1ubq_bonds_are_pinned() {
    let s = vv_io::load(small_fixture("1UBQ.cif")).unwrap();
    let bonds = bonds::perceive(&s.topology, s.frame(0).positions());
    assert_pin(
        &bonds,
        Pin {
            count: 608,
            hash: 0xfa6f4bbc23532d22,
            first: [0, 1],
            last: [599, 601],
        },
        "1UBQ",
    );
}

#[test]
fn hemoglobin_4hhb_bonds_are_pinned() {
    let s = vv_io::load(small_fixture("4HHB.cif")).unwrap();
    let bonds = bonds::perceive(&s.topology, s.frame(0).positions());
    assert_pin(
        &bonds,
        Pin {
            count: 4647,
            hash: 0xd9eb0bf5e66d8579,
            first: [0, 1],
            last: [4556, 4557],
        },
        "4HHB",
    );
}

#[test]
fn glycoprotein_6x3z_bonds_are_pinned() {
    let s = vv_io::load(glycan_fixture("6X3Z.pdb")).unwrap();
    let bonds = bonds::perceive(&s.topology, s.frame(0).positions());
    assert_pin(
        &bonds,
        Pin {
            count: 17829,
            hash: 0x6e88eec8e74d716d,
            first: [0, 1],
            last: [17357, 17364],
        },
        "6X3Z",
    );
}

/// 20k-atom synthetic system: small enough for a non-`#[ignore]`-d test in
/// debug mode, but large enough (many residues, packed atoms) to exercise
/// the same over-cap "priority to shortest" contention path the 1M-atom
/// timing fixture does.
#[test]
fn synthetic_20k_bonds_are_pinned() {
    let s = protein_like(&SynthParams::new(20_000));
    let bonds = bonds::perceive(&s.topology, s.frame(0).positions());
    assert_pin(
        &bonds,
        Pin {
            count: 29303,
            hash: 0x8e8242f71d747c87,
            first: [0, 1],
            last: [19998, 19999],
        },
        "synthetic 20k",
    );
}
