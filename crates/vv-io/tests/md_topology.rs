//! `vv_io::load_topology`'s format dispatch: `.psf`/`.prmtop`/`.parm7`
//! read directly (no coordinates), everything else falls back to
//! `vv_io::load` (`Topology` only, coordinates dropped).

use std::path::PathBuf;

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../fixtures/small")
        .join(name)
}

fn write_temp(name: &str, content: &str) -> PathBuf {
    let dir = std::env::temp_dir().join("vizviz-tests");
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join(name);
    std::fs::write(&path, content).unwrap();
    path
}

const PSF: &str = "PSF\n\n       0 !NTITLE\n\n       2 !NATOM\n\
       1 A    1    LIG  C1   C      0.000000       12.0110           0\n\
       2 A    1    LIG  C2   C      0.000000       12.0110           0\n\
\n       1 !NBOND: bonds\n       1       2\n";

/// Built with `format!`, not hand-counted whitespace: a PRMTOP field's
/// width is exact (no separator between fields), so a mistyped column
/// count silently shifts every later value instead of erroring.
fn prmtop_text() -> String {
    let mut s = String::new();
    s.push_str("%FLAG POINTERS\n%FORMAT(10I8)\n");
    for v in [2i64, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1] {
        s.push_str(&format!("{v:>8}"));
    }
    s.push('\n');
    s.push_str("%FLAG ATOM_NAME\n%FORMAT(20a4)\n");
    for v in ["C1", "C2"] {
        s.push_str(&format!("{v:<4}"));
    }
    s.push('\n');
    s.push_str("%FLAG CHARGE\n%FORMAT(5E16.8)\n");
    for v in [0.0f32, 0.0] {
        s.push_str(&format!("{v:>16.8}"));
    }
    s.push('\n');
    s.push_str("%FLAG MASS\n%FORMAT(5E16.8)\n");
    for v in [12.011f32, 12.011] {
        s.push_str(&format!("{v:>16.8}"));
    }
    s.push('\n');
    s.push_str("%FLAG RESIDUE_LABEL\n%FORMAT(20a4)\n");
    s.push_str(&format!("{:<4}", "LIG"));
    s.push('\n');
    s.push_str("%FLAG RESIDUE_POINTER\n%FORMAT(10I8)\n");
    s.push_str(&format!("{:>8}", 1));
    s.push('\n');
    s.push_str("%FLAG BONDS_WITHOUT_HYDROGEN\n%FORMAT(10I8)\n");
    for v in [0i64, 3, 1] {
        s.push_str(&format!("{v:>8}"));
    }
    s.push('\n');
    s.push_str("%FLAG BONDS_INC_HYDROGEN\n%FORMAT(10I8)\n");
    s.push('\n');
    s
}

#[test]
fn psf_extension_is_read_as_a_topology() {
    let path = write_temp("dispatch.psf", PSF);
    let t = vv_io::load_topology(&path).unwrap();
    assert_eq!(t.atom_count(), 2);
    assert_eq!(t.md_bonds, Some(vec![[0, 1]]));
}

#[test]
fn prmtop_and_parm7_extensions_are_read_as_a_topology() {
    for ext in ["prmtop", "parm7"] {
        let path = write_temp(&format!("dispatch.{ext}"), &prmtop_text());
        let t = vv_io::load_topology(&path).unwrap();
        assert_eq!(t.atom_count(), 2, "{ext}");
        assert_eq!(t.md_bonds, Some(vec![[0, 1]]), "{ext}");
    }
}

#[test]
fn a_pdb_topology_falls_back_to_the_generic_loader() {
    let t = vv_io::load_topology(fixture("1CRN.pdb")).unwrap();
    assert_eq!(t.atom_count(), 327);
    assert_eq!(t.md_bonds, None);
}

#[test]
fn an_unrecognized_extension_is_a_clear_error() {
    let path = write_temp("dispatch.xyz.notaformat", "not a real structure file");
    let err = vv_io::load_topology(&path).unwrap_err();
    assert!(matches!(err, vv_io::ParseError::UnknownFormat(_)));
}
