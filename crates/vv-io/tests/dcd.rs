//! `fixtures/small/sample.dcd` is a real DCD file written by MDAnalysis's
//! own `DCDWriter` (not hand-built, and not written by our own code) --
//! generated with 5 atoms and 4 frames, where atom `i` in frame `f` sits
//! at `(i*10 + f, i*10 + f + 0.5, i*10 + f + 0.25)`. Its exact byte layout
//! was decoded by hand to write `vv_io::dcd`, then these tests check our
//! reader reproduces the known values -- validation against a real
//! external writer, not just internal self-consistency.

use std::path::PathBuf;

use vv_core::glam::Vec3;
use vv_io::dcd::DcdReader;

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../fixtures/small")
        .join(name)
}

fn expected(atom: usize, frame: usize) -> Vec3 {
    let base = (atom * 10 + frame) as f32;
    Vec3::new(base, base + 0.5, base + 0.25)
}

#[test]
fn reports_atom_and_frame_count() {
    let r = DcdReader::open(fixture("sample.dcd")).unwrap();
    assert_eq!(r.atom_count(), 5);
    assert_eq!(r.frame_count(), 4);
}

#[test]
fn reads_known_positions_for_every_frame() {
    let r = DcdReader::open(fixture("sample.dcd")).unwrap();
    for frame in 0..r.frame_count() {
        let positions = r.read_frame(frame).unwrap();
        assert_eq!(positions.len(), 5);
        for (atom, &p) in positions.iter().enumerate() {
            assert_eq!(p, expected(atom, frame), "atom {atom} frame {frame}");
        }
    }
}

#[test]
fn seeking_out_of_order_gives_the_same_answer_as_in_order() {
    let r = DcdReader::open(fixture("sample.dcd")).unwrap();
    // Read frame 3 with nothing read before it: proves a frame is decoded
    // from its own byte offset, not by replaying frames 0..index.
    let last = r.read_frame(3).unwrap();
    assert_eq!(last.len(), 5);
    for (atom, &p) in last.iter().enumerate() {
        assert_eq!(p, expected(atom, 3), "atom {atom}");
    }
    // Then jump backward to frame 0 in the same reader.
    let first = r.read_frame(0).unwrap();
    assert_eq!(first.len(), 5);
    for (atom, &p) in first.iter().enumerate() {
        assert_eq!(p, expected(atom, 0), "atom {atom}");
    }
}

#[test]
fn frame_out_of_range_is_a_clear_error() {
    let r = DcdReader::open(fixture("sample.dcd")).unwrap();
    let err = r.read_frame(4).unwrap_err();
    assert!(matches!(
        err,
        vv_io::dcd::DcdError::FrameOutOfRange {
            index: 4,
            frame_count: 4
        }
    ));
}

#[test]
fn rejects_a_non_dcd_file() {
    assert!(DcdReader::open(fixture("1CRN.pdb")).is_err());
}

#[test]
fn matches_a_real_structures_positions_plus_a_known_offset() {
    // fixtures/small/1CRN_traj.dcd was generated from the real 1CRN.pdb
    // positions (via MDAnalysis's DCDWriter): frame f = base positions +
    // f * (1, 0, 0). Cross-checks the reader against our own PDB parser's
    // output, not just against values we made up.
    let base = vv_io::load(fixture("1CRN.pdb")).unwrap();
    let coords = base.frame(0);
    let base_positions = coords.positions();
    let r = DcdReader::open(fixture("1CRN_traj.dcd")).unwrap();
    assert_eq!(r.atom_count(), base.atom_count());
    assert_eq!(r.frame_count(), 3);
    for frame in 0..3 {
        let positions = r.read_frame(frame).unwrap();
        for atom in 0..base.atom_count() {
            let expected = base_positions[atom] + Vec3::new(frame as f32, 0.0, 0.0);
            let delta = (positions[atom] - expected).abs().max_element();
            assert!(
                delta < 1e-4,
                "atom {atom} frame {frame}: {} vs {expected}",
                positions[atom]
            );
        }
    }
}

#[test]
fn missing_file_is_an_io_error() {
    let err = DcdReader::open(fixture("does-not-exist.dcd"))
        .err()
        .unwrap();
    assert!(matches!(err, vv_io::dcd::DcdError::Io { .. }));
}
