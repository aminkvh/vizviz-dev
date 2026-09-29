//! Exercises the built `vv-dssp` binary as a subprocess, the way an
//! external user would run it.

use std::path::PathBuf;
use std::process::Command;

fn bin() -> Command {
    Command::new(env!("CARGO_BIN_EXE_vv-dssp"))
}

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../fixtures/small")
        .join(name)
}

#[test]
fn assigns_crambins_helix() {
    let out = bin().arg(fixture("1CRN.cif")).output().unwrap();
    assert!(out.status.success(), "{:?}", out);
    let text = String::from_utf8(out.stdout).unwrap();
    let mut lines = text.lines();
    assert_eq!(lines.next(), Some("CHAIN\tSEQID\tRES\tSS"));
    // Residue 10 (ARG) sits inside the alpha helix the fixture's own HELIX
    // record and vv-io/tests/dssp.rs both pin to residues 7-19.
    assert!(text.lines().any(|l| l == "A\t10\tARG\tH"), "{text}");
}

#[test]
fn help_exits_zero_without_a_file() {
    let out = bin().arg("--help").output().unwrap();
    assert!(out.status.success());
    assert!(String::from_utf8_lossy(&out.stdout).contains("usage: vv-dssp"));
}

#[test]
fn a_missing_file_is_a_clean_error_not_a_panic() {
    let out = bin().arg("no-such-file.cif").output().unwrap();
    assert_eq!(out.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&out.stderr).contains("no-such-file.cif"));
}

#[test]
fn an_unknown_flag_exits_two() {
    let out = bin().arg("--bogus").output().unwrap();
    assert_eq!(out.status.code(), Some(2));
}

#[test]
fn no_arguments_prints_usage_and_exits_two() {
    let out = bin().output().unwrap();
    assert_eq!(out.status.code(), Some(2));
    assert!(String::from_utf8_lossy(&out.stderr).contains("usage: vv-dssp"));
}
