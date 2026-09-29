//! BinaryCIF against the same entries as mmCIF and PDB (`fixtures/small`).

use std::collections::HashMap;
use std::path::PathBuf;

use vv_core::glam::Vec3;
use vv_core::{SecondaryStructure, Structure, Topology};
use vv_io::Format;

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

type AtomKey = (String, i32, u8, String, String, u8);

/// Atoms by (author chain, author number, insertion code, residue name,
/// atom name, altloc) -> frame-0 position.
fn keyed(s: &Structure) -> HashMap<AtomKey, Vec3> {
    let t = &s.topology;
    let mut out = HashMap::new();
    for (ri, res) in t.residues.iter().enumerate() {
        let chain = t.names.get(t.chains[res.chain as usize].auth_asym);
        for a in res.atoms.clone().map(|a| a as usize) {
            let key = (
                chain.to_string(),
                res.auth_seq_id,
                res.ins_code,
                t.residue_name(ri).to_string(),
                t.atom_name(a).to_string(),
                t.alt_loc[a],
            );
            assert!(out.insert(key.clone(), s.frame(0).positions()[a]).is_none());
        }
    }
    out
}

fn ss_keys(t: &Topology, ss: SecondaryStructure) -> Vec<(String, i32, u8)> {
    let mut v: Vec<_> = t
        .residues
        .iter()
        .filter(|r| r.ss == ss)
        .map(|r| {
            let chain = t.names.get(t.chains[r.chain as usize].auth_asym);
            (chain.to_string(), r.auth_seq_id, r.ins_code)
        })
        .collect();
    v.sort();
    v
}

fn assert_same_atoms(id: &str, a: &Structure, b: &Structure) {
    let (ka, kb) = (keyed(a), keyed(b));
    assert_eq!(ka.len(), kb.len(), "{id} atoms");
    for (key, pos) in &ka {
        let other = kb
            .get(key)
            .unwrap_or_else(|| panic!("{id}: {key:?} missing"));
        assert!((*pos - *other).abs().max_element() < 1e-3, "{id} {key:?}");
    }
    for ss in [SecondaryStructure::Helix, SecondaryStructure::Strand] {
        assert_eq!(
            ss_keys(&a.topology, ss),
            ss_keys(&b.topology, ss),
            "{id} {ss:?}"
        );
    }
}

#[test]
fn agrees_with_mmcif_and_pdb_on_the_fixture_entries() {
    for id in ["1CRN", "1AKE", "4HHB"] {
        let b = load(&format!("{id}.bcif"));
        let c = load(&format!("{id}.cif"));
        let (tb, tc) = (&b.topology, &c.topology);
        assert_eq!(tb.chain_count(), tc.chain_count(), "{id} chains");
        assert_eq!(tb.element, tc.element, "{id} elements");
        assert_eq!(tb.charge, tc.charge, "{id} charges");
        assert_eq!(tb.residue_class, tc.residue_class, "{id} classes");
        assert_eq!((&tb.id, &tb.title), (&tc.id, &tc.title), "{id} header");
        assert_eq!(b.frame_count(), c.frame_count(), "{id} models");
        assert_same_atoms(id, &b, &c);
        assert_same_atoms(id, &b, &load(&format!("{id}.pdb")));
    }
}

#[test]
fn is_recognised_by_extension_and_by_content() {
    let bytes = std::fs::read(fixture("1CRN.bcif")).unwrap();
    assert_eq!(Format::sniff(&bytes), Some(Format::Bcif));
    assert_eq!(Format::from_path(&fixture("x.bcif.gz")), Some(Format::Bcif));
    let s = vv_io::parse(&bytes, Format::Bcif).unwrap();
    assert_eq!(s.topology.atom_count(), 327);
}

#[test]
fn a_gzipped_file_loads_like_the_plain_one() {
    use std::io::Write;
    let bytes = std::fs::read(fixture("1CRN.bcif")).unwrap();
    let mut enc = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
    enc.write_all(&bytes).unwrap();
    let path = std::env::temp_dir().join(format!("vv_bcif_{}.bcif.gz", std::process::id()));
    std::fs::write(&path, enc.finish().unwrap()).unwrap();
    let gz = vv_io::load(&path).unwrap();
    std::fs::remove_file(&path).ok();
    assert_same_atoms("1CRN.gz", &gz, &load("1CRN.bcif"));
}

#[test]
fn a_truncated_file_is_an_error() {
    let bytes = std::fs::read(fixture("1CRN.bcif")).unwrap();
    for cut in [10, bytes.len() / 2, bytes.len() - 1] {
        assert!(
            vv_io::parse(&bytes[..cut], Format::Bcif).is_err(),
            "cut {cut}"
        );
    }
}
