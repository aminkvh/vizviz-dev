//! Trajectory readers against MDAnalysis: `fixtures/traj/make.py` wrote
//! each file and, next to it, the coordinates MDAnalysis reads back from
//! it (raw little-endian f32, frames x atoms x 3).

use std::path::PathBuf;

use vv_core::glam::Vec3;

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../fixtures/traj")
        .join(name)
}

fn reference(name: &str, atoms: usize) -> Vec<Vec<Vec3>> {
    let bytes = std::fs::read(fixture(&format!("{name}.f32"))).unwrap();
    let values: Vec<f32> = bytes
        .chunks_exact(4)
        .map(|b| f32::from_le_bytes([b[0], b[1], b[2], b[3]]))
        .collect();
    values
        .chunks_exact(atoms * 3)
        .map(|f| {
            f.chunks_exact(3)
                .map(|p| Vec3::new(p[0], p[1], p[2]))
                .collect()
        })
        .collect()
}

fn assert_frames_match(read: impl Fn(usize) -> Vec<Vec3>, expected: &[Vec<Vec3>], tolerance: f32) {
    for (k, frame) in expected.iter().enumerate() {
        let got = read(k);
        assert_eq!(got.len(), frame.len());
        let worst = got
            .iter()
            .zip(frame)
            .map(|(a, b)| a.distance(*b))
            .fold(0.0f32, f32::max);
        assert!(worst < tolerance, "frame {k}: off by {worst} A");
    }
}

#[test]
fn xtc_matches_mdanalysis() {
    let reader = vv_io::xtc::XtcReader::open(fixture("1crn.xtc")).unwrap();
    assert_eq!(reader.atom_count(), 327);
    assert_eq!(reader.frame_count(), 10);
    let expected = reference("1crn.xtc", 327);
    // Both decode the same integers; only the float scaling may differ.
    assert_frames_match(|k| reader.read_frame(k).unwrap(), &expected, 1e-4);
    // Random access: a late frame first.
    let last = reader.read_frame(9).unwrap();
    assert!(last[0].distance(expected[9][0]) < 1e-4);
}

#[test]
fn trr_matches_mdanalysis() {
    let reader = vv_io::trr::TrrReader::open(fixture("1crn.trr")).unwrap();
    assert_eq!((reader.atom_count(), reader.frame_count()), (327, 10));
    let expected = reference("1crn.trr", 327);
    assert_frames_match(|k| reader.read_frame(k).unwrap(), &expected, 1e-4);
}

#[test]
fn netcdf_matches_mdanalysis() {
    let reader = vv_io::netcdf::NetCdfReader::open(fixture("1crn.nc")).unwrap();
    assert_eq!((reader.atom_count(), reader.frame_count()), (327, 10));
    let expected = reference("1crn.nc", 327);
    assert_frames_match(|k| reader.read_frame(k).unwrap(), &expected, 1e-4);
}

#[test]
fn every_format_opens_through_one_reader() {
    for name in ["1crn.xtc", "1crn.trr", "1crn.nc"] {
        let t = vv_io::trajectory::Trajectory::open(fixture(name)).unwrap();
        assert_eq!((t.atom_count(), t.frame_count()), (327, 10), "{name}");
        let expected = reference(name, 327);
        assert_frames_match(|k| t.read_frame(k).unwrap(), &expected, 1e-4);
    }
}
