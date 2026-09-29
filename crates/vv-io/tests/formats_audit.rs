//! One test per rule of `docs/FORMATS.md`'s PDB / mmCIF checklist, on small
//! hand-written snippets, plus round-trip and cross-format agreement checks
//! on `fixtures/small`. Real-entry checks (downloaded, not committed) are
//! `#[ignore]`d at the bottom.

use std::collections::HashMap;
use std::path::PathBuf;

use vv_core::glam::Vec3;
use vv_core::{flags, Element, ExplicitBondKind, SecondaryStructure, Structure, Topology};
use vv_io::{Format, SaveOptions};

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

fn pdb(text: &str) -> Structure {
    let s = vv_io::pdb::parse(text.as_bytes()).unwrap();
    s.topology.validate().unwrap();
    s
}

fn cif(text: &str) -> Structure {
    let s = vv_io::mmcif::parse(text.as_bytes()).unwrap();
    s.topology.validate().unwrap();
    s
}

/// One 80-column ATOM/HETATM record; `res` is columns 18-21.
#[allow(clippy::too_many_arguments)]
fn atom(
    rec: &str,
    serial: &str,
    name: &str,
    alt: char,
    res: &str,
    chain: char,
    seq: &str,
    icode: char,
    x: f32,
    element: &str,
    charge: &str,
) -> String {
    format!(
        "{rec:<6}{serial:>5} {name:<4}{alt}{res:<4}{chain}{seq:>4}{icode}   {x:8.3}{y:8.3}{z:8.3}{occ:6.2}{b:6.2}          {element:>2}{charge:<2}\n",
        y = 0.0,
        z = 0.0,
        occ = 1.0,
        b = 10.0,
    )
}

fn simple(name: &str, res: &str, chain: char, seq: &str, icode: char, x: f32) -> String {
    atom("ATOM", "1", name, ' ', res, chain, seq, icode, x, "", "")
}

// ---- PDB: column layout ---------------------------------------------------

#[test]
fn pdb_four_character_residue_name_in_columns_18_to_21() {
    let s = pdb(&atom(
        "ATOM", "1", " CA ", ' ', "ALA2", 'A', "1", ' ', 1.0, "C", "",
    ));
    assert_eq!(s.topology.residue_name(0), "ALA2");
    assert_eq!(s.topology.chain_name(0), "A");
}

#[test]
fn pdb_insertion_codes_make_distinct_residues() {
    let text = simple(" CA ", "SER", 'H', "52", ' ', 0.0)
        + &simple(" CA ", "SER", 'H', "52", 'A', 4.0)
        + &simple(" CA ", "SER", 'H', "52", 'B', 8.0)
        + &simple(" CA ", "SER", 'H', "53", ' ', 12.0);
    let t = pdb(&text).topology;
    assert_eq!(t.residue_count(), 4);
    let codes: Vec<u8> = t.residues.iter().map(|r| r.ins_code).collect();
    assert_eq!(codes, [0, b'A', b'B', 0]);
}

#[test]
fn pdb_alternate_locations_are_kept_and_tagged() {
    let text = atom("ATOM", "1", " CB ", 'A', "SER", 'A', "5", ' ', 1.0, "C", "")
        + &atom("ATOM", "2", " CB ", 'B', "SER", 'A', "5", ' ', 2.0, "C", "");
    let s = pdb(&text);
    assert_eq!(s.topology.atom_count(), 2);
    assert_eq!(s.topology.alt_loc, [b'A', b'B']);
    assert_eq!(s.topology.residue_count(), 1);
}

#[test]
fn pdb_microheterogeneity_splits_residues_by_name() {
    let text = atom("ATOM", "1", " CA ", 'A', "SER", 'A', "5", ' ', 1.0, "C", "")
        + &atom("ATOM", "2", " CA ", 'B', "THR", 'A', "5", ' ', 1.0, "C", "");
    assert_eq!(pdb(&text).topology.residue_count(), 2);
}

#[test]
fn pdb_charge_column_accepts_digit_sign_and_sign_digit() {
    let text = atom(
        "HETATM", "1", "ZN  ", ' ', "ZN", 'A', "1", ' ', 0.0, "ZN", "2+",
    ) + &atom(
        "HETATM", "2", "CL  ", ' ', "CL", 'A', "2", ' ', 5.0, "CL", "1-",
    ) + &atom(
        "HETATM", "3", "FE  ", ' ', "FE", 'A', "3", ' ', 9.0, "FE", "+3",
    ) + &atom(
        "HETATM", "4", "NA  ", ' ', "NA", 'A', "4", ' ', 12.0, "NA", "",
    );
    assert_eq!(pdb(&text).topology.charge, [2, -1, 3, 0]);
}

#[test]
fn pdb_element_column_wins_deuterium_is_hydrogen_and_blank_falls_back_to_name() {
    let text = atom("ATOM", "1", " CA ", ' ', "ALA", 'A', "1", ' ', 0.0, "C", "")
        + &atom("ATOM", "2", "CA  ", ' ', "CA", 'A', "2", ' ', 5.0, "CA", "")
        + &atom("ATOM", "3", " D  ", ' ', "HOH", 'A', "3", ' ', 9.0, "D", "")
        + &atom("ATOM", "4", " OG ", ' ', "SER", 'A', "4", ' ', 12.0, "", "");
    let t = pdb(&text).topology;
    assert_eq!(t.element[0], Element::CARBON);
    assert_eq!(t.element[1], Element::from_atomic_number(20).unwrap());
    assert_eq!(t.element[2], Element::HYDROGEN);
    assert_eq!(t.element[3], Element::OXYGEN);
}

#[test]
fn pdb_truncated_atom_record_is_an_error_not_a_panic() {
    let err = vv_io::pdb::parse(b"ATOM      1  N   ALA A   1\n").unwrap_err();
    assert!(matches!(err, vv_io::ParseError::Malformed { line: 1, .. }));
    assert!(vv_io::pdb::parse(b"ATOM  \n").is_err());
}

#[test]
fn pdb_segment_id_names_the_chain_when_column_22_is_blank() {
    let mut line = simple(" CA ", "ALA", ' ', "1", ' ', 0.0);
    line.replace_range(72..76, "PROA");
    assert_eq!(pdb(&line).topology.chain_name(0), "PROA");
}

#[test]
fn pdb_hetatm_flag_and_serial_hybrid36_past_99999() {
    let mut text = String::new();
    for serial in ["99999", "A0000", "A0001"] {
        text += &atom(
            "HETATM", serial, " O  ", ' ', "HOH", 'A', "1", ' ', 0.0, "O", "",
        );
    }
    let t = pdb(&text).topology;
    assert_eq!(t.serial, [99999, 100000, 100001]);
    assert_eq!(t.flags[0] & flags::HETERO, flags::HETERO);
}

#[test]
fn pdb_residue_number_hybrid36_past_9999() {
    let text = simple(" CA ", "ALA", 'A', "9999", ' ', 0.0)
        + &simple(" CA ", "ALA", 'A', "A000", ' ', 4.0)
        + &simple(" CA ", "ALA", 'A', "A001", ' ', 8.0);
    let t = pdb(&text).topology;
    let seq: Vec<i32> = t.residues.iter().map(|r| r.auth_seq_id).collect();
    assert_eq!(seq, [9999, 10000, 10001]);
}

#[test]
fn pdb_ter_starts_a_new_chain_even_with_a_repeated_chain_letter() {
    let text = simple(" CA ", "ALA", 'A', "1", ' ', 0.0)
        + "TER\n"
        + &simple(" CA ", "ALA", 'A', "1", ' ', 9.0);
    let t = pdb(&text).topology;
    assert_eq!((t.chain_count(), t.residue_count()), (2, 2));
}

#[test]
fn ter_split_chain_records_still_answer_to_one_chain_name() {
    let text = simple(" CA ", "ALA", 'A', "1", ' ', 0.0)
        + "TER
" + &atom(
        "HETATM", "9", " O  ", ' ', "HOH", 'A', "2", ' ', 9.0, "O", "",
    );
    let s = pdb(&text);
    assert_eq!(s.topology.chain_count(), 2);
    assert_eq!(s.topology.chain_name_index(), [0, 0]);
    let positions = s.frame(0).positions().to_vec();
    let hit = vv_core::select(&s.topology, &positions, "chain A").unwrap();
    assert!(hit[0] && hit[1]);
}

#[test]
fn pdb_end_separated_frames_with_matching_atoms_are_frames() {
    let text = simple(" CA ", "ALA", 'A', "1", ' ', 0.0)
        + "END
" + &simple(" CA ", "ALA", 'A', "1", ' ', 1.0)
        + "END
" + &simple(" CA ", "ALA", 'A', "1", ' ', 2.0);
    let s = pdb(&text);
    assert_eq!(s.topology.atom_count(), 1);
    let xs: Vec<f32> = (0..s.frame_count())
        .map(|f| s.frame(f).positions()[0].x)
        .collect();
    assert_eq!(xs, [0.0, 1.0, 2.0]);
}

#[test]
fn pdb_end_then_a_different_entry_keeps_only_the_first() {
    let text = simple(" CA ", "ALA", 'A', "1", ' ', 0.0)
        + "END
" + &simple(" CB ", "ALA", 'A', "1", ' ', 1.0)
        + &simple(" CA ", "GLY", 'A', "2", ' ', 5.0);
    let s = pdb(&text);
    assert_eq!((s.topology.atom_count(), s.frame_count()), (1, 1));
}

#[test]
fn pdb_models_become_frames_and_short_models_are_dropped() {
    let text = format!(
        "MODEL        1\n{}ENDMDL\nMODEL        2\n{}ENDMDL\nMODEL        3\n{}{}ENDMDL\n",
        simple(" CA ", "ALA", 'A', "1", ' ', 0.0),
        simple(" CA ", "ALA", 'A', "1", ' ', 1.0),
        simple(" CA ", "ALA", 'A', "1", ' ', 2.0),
        simple(" CB ", "ALA", 'A', "1", ' ', 3.0),
    );
    let s = pdb(&text);
    assert_eq!(s.frame_count(), 2);
    assert_eq!(s.frame(1).positions()[0].x, 1.0);
}

#[test]
fn pdb_anisou_seqres_and_remarks_are_not_atoms() {
    let text = "REMARK   3 SOMETHING\nSEQRES   1 A    3  ALA GLY SER\nANISOU    1  CA  ALA A   1     1000   1000   1000      0      0      0       C\n".to_string()
        + &simple(" CA ", "ALA", 'A', "1", ' ', 0.0);
    assert_eq!(pdb(&text).topology.atom_count(), 1);
}

#[test]
fn pdb_cryst1_populates_cell_and_symmetry() {
    let text = "CRYST1   40.960   18.650   22.520  90.00  90.77  90.00 P 21          2\n"
        .to_string()
        + &simple(" CA ", "ALA", 'A', "1", ' ', 0.0);
    let t = pdb(&text).topology;
    assert_eq!(t.annotations.get("cell", "length_a"), Some("40.960"));
    assert_eq!(t.annotations.get("cell", "angle_beta"), Some("90.77"));
    assert_eq!(
        t.annotations.get("symmetry", "space_group_name_H-M"),
        Some("P 21")
    );
}

// ---- PDB: secondary structure and links -----------------------------------

#[test]
fn pdb_helix_range_covers_insertion_code_residues_but_not_ligands_reusing_numbers() {
    let mut text = String::from(
        "HELIX    1   1 SER A   52  THR A   53  1                                   2\n",
    );
    text += &simple(" CA ", "SER", 'A', "52", ' ', 0.0);
    text += &simple(" CA ", "GLY", 'A', "52", 'A', 4.0);
    text += &simple(" CA ", "THR", 'A', "53", ' ', 8.0);
    text += &simple(" CA ", "GLY", 'A', "54", ' ', 12.0);
    text += "TER\n";
    text += &atom(
        "HETATM", "9", " C1 ", ' ', "LIG", 'A', "52", ' ', 20.0, "C", "",
    );
    let t = pdb(&text).topology;
    let ss: Vec<_> = t.residues.iter().map(|r| r.ss).collect();
    use SecondaryStructure::*;
    assert_eq!(ss, [Helix, Helix, Helix, Unknown, Unknown]);
}

#[test]
fn pdb_sheet_strand_read_from_v33_columns() {
    let mut text = String::from("SHEET    1   A 2 THR A   2  CYS A   3  0\n");
    for (i, res) in ["THR", "THR", "CYS", "CYS"].iter().enumerate() {
        text += &simple(
            " CA ",
            res,
            'A',
            &format!("{}", i / 2 + 2),
            ' ',
            i as f32 * 4.0,
        );
    }
    let t = pdb(&text).topology;
    assert!(t
        .residues
        .iter()
        .all(|r| r.ss == SecondaryStructure::Strand));
}

#[test]
fn pdb_ssbond_is_a_disulfide_and_link_to_a_metal_is_metal() {
    let mut text = String::new();
    text += &atom("ATOM", "1", " SG ", ' ', "CYS", 'A', "3", ' ', 0.0, "S", "");
    text += &atom(
        "ATOM", "2", " SG ", ' ', "CYS", 'A', "40", ' ', 50.0, "S", "",
    );
    text += &atom(
        "HETATM", "3", "ZN  ", ' ', "ZN", 'A', "101", ' ', 90.0, "ZN", "2+",
    );
    text += "SSBOND   1 CYS A    3    CYS A   40\n";
    text += "CONECT    1    2\n";
    let t = pdb(&text).topology;
    assert_eq!(
        t.explicit_bonds.len(),
        1,
        "CONECT and SSBOND agree on one bond"
    );
    assert_eq!(t.explicit_bonds[0].kind, ExplicitBondKind::Disulfide);
}

#[test]
fn pdb_link_between_metal_and_atom_has_metal_kind() {
    let mut text = String::new();
    text += &atom("ATOM", "1", " NE2", ' ', "HIS", 'A', "7", ' ', 0.0, "N", "");
    text += &atom(
        "HETATM", "2", "ZN  ", ' ', "ZN", 'A', "101", ' ', 2.0, "ZN", "2+",
    );
    let mut link = vec![b' '; 80];
    link[..6].copy_from_slice(b"LINK  ");
    let put =
        |l: &mut Vec<u8>, at: usize, s: &str| l[at..at + s.len()].copy_from_slice(s.as_bytes());
    put(&mut link, 12, " NE2");
    put(&mut link, 17, "HIS");
    link[21] = b'A';
    put(&mut link, 22, "   7");
    put(&mut link, 42, "ZN  ");
    put(&mut link, 47, "ZN");
    link[51] = b'A';
    put(&mut link, 52, " 101");
    text += std::str::from_utf8(&link).unwrap();
    text += "\n";
    let t = pdb(&text).topology;
    assert_eq!(t.explicit_bonds.len(), 1);
    assert_eq!(t.explicit_bonds[0].kind, ExplicitBondKind::Metal);
}

// ---- mmCIF ----------------------------------------------------------------

const CIF_HEAD: &str = "loop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_alt_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_entity_id\n_atom_site.label_seq_id\n_atom_site.pdbx_PDB_ins_code\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.occupancy\n_atom_site.B_iso_or_equiv\n_atom_site.pdbx_formal_charge\n_atom_site.auth_seq_id\n_atom_site.auth_comp_id\n_atom_site.auth_asym_id\n_atom_site.auth_atom_id\n_atom_site.pdbx_PDB_model_num\n";

#[test]
fn cif_label_and_auth_identifiers_are_both_kept() {
    let text = format!(
        "data_x\n{CIF_HEAD}ATOM 1 C CA . ALA A 1 1 ? 0 0 0 1.0 5.0 ? 27 ALA H CA 1\nHETATM 2 O O . HOH B 2 . ? 5 0 0 1.0 5.0 ? 301 HOH H O 1\n"
    );
    let t = cif(&text).topology;
    assert_eq!(t.chain_name(0), "A");
    assert_eq!(t.residues[0].seq_id, 1);
    assert_eq!(t.residues[0].auth_seq_id, 27);
    let auth = t.names.get(t.chains[0].auth_asym);
    assert_eq!(auth, "H");
    // A non-polymer has no label_seq_id: its author number keys the residue.
    assert_eq!(t.residues[1].seq_id, 301);
}

#[test]
fn cif_column_order_nulls_charge_altloc_and_insertion_code() {
    let text = format!(
        "data_x\n{CIF_HEAD}HETATM 1 ZN ZN A ZN A 1 . A 0 0 0 0.5 5.0 2 101 ZN A ZN 1\nATOM 2 C CA . ALA B 2 3 B 5 0 0 1.0 5.0 ? 3 ALA B CA 1\n"
    );
    let t = cif(&text).topology;
    assert_eq!(t.alt_loc, [b'A', 0]);
    assert_eq!(t.charge, [2, 0]);
    assert_eq!(t.residues[1].ins_code, b'B');
    assert_eq!(t.occupancy[0], 0.5);
}

#[test]
fn cif_quoted_values_semicolon_text_and_primed_atom_names() {
    let text = "data_x\n_entry.id X\n_struct.title\n;A long\ntitle over two lines\n;\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C \"C1'\" DA A 0 0 0\nATOM 2 O 'O5''' DA A 1 0 0\nATOM 3 C 'A B' DA A 2 0 0\n";
    let t = cif(text).topology;
    assert!(t.title.starts_with("A long"));
    assert_eq!(t.atom_name(0), "C1'");
    assert_eq!(t.atom_name(2), "A B");
}

#[test]
fn cif_rows_wrapped_over_several_lines_still_parse() {
    let text = "data_x\nloop_\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n1 C CA ALA A\n  1.0 2.0\n  3.0\n2 N N ALA A 4.0 5.0 6.0\n";
    let s = cif(text);
    assert_eq!(s.topology.atom_count(), 2);
    assert_eq!(s.frame(0).positions()[0], Vec3::new(1.0, 2.0, 3.0));
    assert_eq!(s.frame(0).positions()[1], Vec3::new(4.0, 5.0, 6.0));
}

#[test]
fn cif_single_atom_written_as_key_value_pairs() {
    let text = "data_x\n_atom_site.id 1\n_atom_site.type_symbol ZN\n_atom_site.label_atom_id ZN\n_atom_site.label_comp_id ZN\n_atom_site.label_asym_id A\n_atom_site.Cartn_x 1.5\n_atom_site.Cartn_y 2.5\n_atom_site.Cartn_z 3.5\n";
    let s = cif(text);
    assert_eq!(s.topology.atom_count(), 1);
    assert_eq!(s.frame(0).positions()[0].z, 3.5);
}

#[test]
fn cif_structure_in_a_later_data_block_is_found() {
    let text = format!(
        "data_dictionary\n_entry.id NOPE\nloop_\n_chem_comp.id\nALA\ndata_real\n_entry.id REAL\n{CIF_HEAD}ATOM 1 C CA . ALA A 1 1 ? 0 0 0 1.0 5.0 ? 1 ALA A CA 1\n"
    );
    let t = cif(&text).topology;
    assert_eq!((t.id.as_str(), t.atom_count()), ("REAL", 1));
}

#[test]
fn cif_models_secondary_structure_and_struct_conn_kinds() {
    let mut text = format!("data_x\n{CIF_HEAD}");
    for (model, dx) in [(1, 0.0), (2, 0.5)] {
        text += &format!("ATOM 1 S SG . CYS A 1 3 ? {dx} 0 0 1.0 5.0 ? 3 CYS A SG {model}\n");
        text += &format!(
            "ATOM 2 S SG . CYS A 1 40 ? {} 0 0 1.0 5.0 ? 40 CYS A SG {model}\n",
            50.0 + dx
        );
    }
    text += "loop_\n_struct_conf.conf_type_id\n_struct_conf.beg_label_asym_id\n_struct_conf.beg_label_seq_id\n_struct_conf.end_label_seq_id\nHELX_P A 3 3\n";
    text += "loop_\n_struct_conn.conn_type_id\n_struct_conn.ptnr1_label_asym_id\n_struct_conn.ptnr1_label_seq_id\n_struct_conn.ptnr1_label_atom_id\n_struct_conn.ptnr2_label_asym_id\n_struct_conn.ptnr2_label_seq_id\n_struct_conn.ptnr2_label_atom_id\ndisulf A 3 SG A 40 SG\ncovale_base A 3 SG A 40 SG\nhydrog A 3 SG A 40 SG\n";
    let s = cif(&text);
    assert_eq!(s.frame_count(), 2);
    assert_eq!(s.topology.residues[0].ss, SecondaryStructure::Helix);
    assert_eq!(s.topology.residues[1].ss, SecondaryStructure::Unknown);
    let kinds: Vec<_> = s.topology.explicit_bonds.iter().map(|b| b.kind).collect();
    assert_eq!(
        kinds,
        [ExplicitBondKind::Disulfide, ExplicitBondKind::Covalent]
    );
}

// ---- round trips ----------------------------------------------------------

static NEXT_FILE: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);

fn round_trip(s: &Structure, format: Format) -> Structure {
    let frames: Vec<usize> = (0..s.frame_count()).collect();
    let opts = SaveOptions {
        atoms: None,
        frames: &frames,
    };
    let ext = match format {
        Format::Pdb => "pdb",
        Format::Mmcif => "cif",
        Format::Bcif => unreachable!("BinaryCIF is read-only"),
    };
    let path = std::env::temp_dir().join(format!(
        "vv_audit_{}_{}.{ext}",
        std::process::id(),
        NEXT_FILE.fetch_add(1, std::sync::atomic::Ordering::Relaxed)
    ));
    vv_io::save(s, &path, &opts).unwrap();
    let back = vv_io::load(&path).unwrap();
    std::fs::remove_file(&path).ok();
    back.topology.validate().unwrap();
    back
}

fn assert_same_structure(a: &Structure, b: &Structure) {
    let (ta, tb) = (&a.topology, &b.topology);
    assert_eq!(ta.atom_count(), tb.atom_count());
    assert_eq!(ta.residue_count(), tb.residue_count());
    assert_eq!(ta.chain_count(), tb.chain_count());
    for i in 0..ta.atom_count() {
        assert_eq!(ta.atom_name(i), tb.atom_name(i), "atom {i} name");
        assert_eq!(ta.element[i], tb.element[i], "atom {i} element");
        assert_eq!(ta.alt_loc[i], tb.alt_loc[i], "atom {i} altloc");
        assert_eq!(ta.charge[i], tb.charge[i], "atom {i} charge");
        assert_eq!(ta.flags[i], tb.flags[i], "atom {i} flags");
    }
    for (ra, rb) in ta.residues.iter().zip(&tb.residues) {
        assert_eq!(ta.names.get(ra.comp), tb.names.get(rb.comp));
        assert_eq!(
            (ra.auth_seq_id, ra.ins_code, ra.ss),
            (rb.auth_seq_id, rb.ins_code, rb.ss)
        );
    }
    assert_eq!(a.frame_count(), b.frame_count());
    for f in 0..a.frame_count() {
        for (pa, pb) in a.frame(f).positions().iter().zip(b.frame(f).positions()) {
            assert!(
                (*pa - *pb).abs().max_element() < 1e-3,
                "frame {f}: {pa} vs {pb}"
            );
        }
    }
}

#[test]
fn pdb_fixtures_round_trip_through_the_pdb_writer() {
    for name in ["1CRN.pdb", "1AKE.pdb", "4HHB.pdb"] {
        let original = load(name);
        let back = round_trip(&original, Format::Pdb);
        assert_same_structure(&original, &back);
    }
}

#[test]
fn pdb_writer_round_trips_four_character_names_hybrid36_and_multiple_models() {
    let mut text = String::from("MODEL        1\n");
    for (i, res) in ["ALA2", "ZN", "LIG3"].iter().enumerate() {
        let serial = 99998 + i as i64;
        let s = vv_io::pdb::encode_hybrid36(serial, 5).unwrap();
        text += &atom(
            "ATOM",
            &s,
            " CA ",
            ' ',
            res,
            'A',
            &vv_io::pdb::encode_hybrid36(9999 + i as i64, 4).unwrap(),
            ' ',
            i as f32 * 5.0,
            "C",
            "",
        );
    }
    text += "ENDMDL\nMODEL        2\n";
    text += &atom(
        "ATOM", "99998", " CA ", ' ', "ALA2", 'A', "9999", ' ', 1.0, "C", "",
    );
    text += &atom(
        "ATOM", "99999", " CA ", ' ', "ZN", 'A', "A000", ' ', 6.0, "C", "",
    );
    text += &atom(
        "ATOM", "A0000", " CA ", ' ', "LIG3", 'A', "A001", ' ', 11.0, "C", "",
    );
    text += "ENDMDL\n";
    let original = pdb(&text);
    assert_eq!(original.topology.serial, [99998, 99999, 100000]);
    let back = round_trip(&original, Format::Pdb);
    assert_eq!(back.topology.residue_name(0), "ALA2");
    assert_eq!(back.topology.residue_name(2), "LIG3");
    assert_eq!(back.topology.residues[2].auth_seq_id, 10001);
    assert_same_structure(&original, &back);
}

#[test]
fn pdb_writer_round_trips_secondary_structure_insertion_codes_cell_and_disulfides() {
    let original = load("1CRN.pdb");
    let back = round_trip(&original, Format::Pdb);
    assert_eq!(
        back.topology.annotations.get("cell", "length_a"),
        original.topology.annotations.get("cell", "length_a")
    );
    let kinds = |t: &Topology| {
        t.explicit_bonds
            .iter()
            .filter(|b| b.kind == ExplicitBondKind::Disulfide)
            .count()
    };
    assert_eq!(kinds(&back.topology), kinds(&original.topology));
    assert!(kinds(&original.topology) >= 3);

    let text =
        simple(" CA ", "SER", 'H', "52", ' ', 0.0) + &simple(" CA ", "SER", 'H', "52", 'A', 4.0);
    let ins = round_trip(&pdb(&text), Format::Pdb);
    assert_eq!(ins.topology.residue_count(), 2);
    assert_eq!(ins.topology.residues[1].ins_code, b'A');
}

#[test]
fn pdb_writer_puts_helix_and_sheet_fields_in_v33_columns() {
    let text = "HELIX    1   1 SER A   52  THR A   53  1                                   2\nSHEET    1   A 1 SER A  60  THR A  61  0\n".to_string()
        + &simple(" CA ", "SER", 'A', "52", ' ', 0.0)
        + &simple(" CA ", "THR", 'A', "53", ' ', 4.0)
        + &simple(" CA ", "SER", 'A', "60", ' ', 8.0)
        + &simple(" CA ", "THR", 'A', "61", ' ', 12.0);
    let s = pdb(&text);
    let mut out = Vec::new();
    vv_io::pdb_write::write(&s, None, &[0], &mut out).unwrap();
    let out = String::from_utf8(out).unwrap();
    let line = |rec: &str| {
        out.lines()
            .find(|l| l.starts_with(rec))
            .unwrap()
            .to_string()
    };
    let (h, sh) = (line("HELIX"), line("SHEET"));
    assert_eq!((&h[15..18], &h[19..20], &h[21..25]), ("SER", "A", "  52"));
    assert_eq!((&h[27..30], &h[31..32], &h[33..37]), ("THR", "A", "  53"));
    assert_eq!(&h[38..40], " 1");
    assert_eq!(&h[71..76], "    2");
    assert_eq!(
        (&sh[17..20], &sh[21..22], &sh[22..26]),
        ("SER", "A", "  60")
    );
    assert_eq!(
        (&sh[28..31], &sh[32..33], &sh[33..37]),
        ("THR", "A", "  61")
    );
    assert_eq!(&sh[38..40], " 0");
}

#[test]
fn mmcif_fixtures_round_trip_through_the_mmcif_writer() {
    for name in ["1CRN.cif", "1AKE.cif"] {
        let original = load(name);
        assert_same_structure(&original, &round_trip(&original, Format::Mmcif));
    }
}

// ---- cross-format agreement -----------------------------------------------

type AtomKey = (String, i32, u8, String, String, u8);

fn keyed(s: &Structure) -> HashMap<AtomKey, Vec3> {
    let t = &s.topology;
    let mut out = HashMap::new();
    for (ri, res) in t.residues.iter().enumerate() {
        let chain = t
            .names
            .get(t.chains[res.chain as usize].auth_asym)
            .to_string();
        for a in res.atoms.clone() {
            let a = a as usize;
            let key = (
                chain.clone(),
                res.auth_seq_id,
                res.ins_code,
                t.residue_name(ri).to_string(),
                t.atom_name(a).to_string(),
                t.alt_loc[a],
            );
            let dup = out.insert(key.clone(), s.frame(0).positions()[a]);
            assert!(dup.is_none(), "duplicate atom key {key:?}");
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
            (
                t.names
                    .get(t.chains[r.chain as usize].auth_asym)
                    .to_string(),
                r.auth_seq_id,
                r.ins_code,
            )
        })
        .collect();
    v.sort();
    v
}

#[test]
fn the_same_entry_agrees_between_pdb_and_mmcif() {
    for id in ["1CRN", "1AKE", "4HHB"] {
        let (p, c) = (load(&format!("{id}.pdb")), load(&format!("{id}.cif")));
        assert_eq!(
            p.topology.atom_count(),
            c.topology.atom_count(),
            "{id} atoms"
        );
        assert_eq!(
            p.topology.residue_count(),
            c.topology.residue_count(),
            "{id} residues"
        );
        let (kp, kc) = (keyed(&p), keyed(&c));
        assert_eq!(kp.len(), kc.len());
        for (key, pos) in &kp {
            let other = kc
                .get(key)
                .unwrap_or_else(|| panic!("{id}: {key:?} missing in mmCIF"));
            assert!((*pos - *other).abs().max_element() < 1e-3, "{id} {key:?}");
        }
        for ss in [SecondaryStructure::Helix, SecondaryStructure::Strand] {
            assert_eq!(
                ss_keys(&p.topology, ss),
                ss_keys(&c.topology, ss),
                "{id} {ss:?}"
            );
        }
        assert_eq!(p.topology.element, c.topology.element, "{id} elements");
        assert_eq!(p.topology.charge, c.topology.charge, "{id} charges");
    }
}

// ---- real entries (downloaded, not committed) -------------------------------

fn real(name: &str) -> Option<PathBuf> {
    let p = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../fixtures/real")
        .join(name);
    p.exists().then_some(p)
}

/// Like [`keyed`] but tolerates duplicate keys (last wins).
fn keyed_lossy(s: &Structure) -> HashMap<AtomKey, Vec3> {
    let t = &s.topology;
    let mut out = HashMap::new();
    for (ri, res) in t.residues.iter().enumerate() {
        let chain = t
            .names
            .get(t.chains[res.chain as usize].auth_asym)
            .to_string();
        for a in res.atoms.clone() {
            let a = a as usize;
            let key = (
                chain.clone(),
                res.auth_seq_id,
                res.ins_code,
                t.residue_name(ri).to_string(),
                t.atom_name(a).to_string(),
                t.alt_loc[a],
            );
            out.insert(key, s.frame(0).positions()[a]);
        }
    }
    out
}

/// Compares one entry's `.pdb` and `.cif`; returns the discrepancies.
fn diff_formats(id: &str) -> Vec<String> {
    let (Some(pp), Some(cp)) = (real(&format!("{id}.pdb")), real(&format!("{id}.cif"))) else {
        return vec![format!("{id}: files not downloaded")];
    };
    let (p, c) = (vv_io::load(pp).unwrap(), vv_io::load(cp).unwrap());
    let mut out = Vec::new();
    let (tp, tc) = (&p.topology, &c.topology);
    if tp.atom_count() != tc.atom_count() {
        out.push(format!(
            "{id}: atoms {} vs {}",
            tp.atom_count(),
            tc.atom_count()
        ));
    }
    if p.frame_count() != c.frame_count() {
        out.push(format!(
            "{id}: frames {} vs {}",
            p.frame_count(),
            c.frame_count()
        ));
    }
    let (kp, kc) = (keyed_lossy(&p), keyed_lossy(&c));
    let only_pdb: Vec<_> = kp.keys().filter(|k| !kc.contains_key(*k)).collect();
    let only_cif: Vec<_> = kc.keys().filter(|k| !kp.contains_key(*k)).collect();
    if !only_pdb.is_empty() || !only_cif.is_empty() {
        out.push(format!(
            "{id}: only in pdb {} {:?}; only in cif {} {:?}",
            only_pdb.len(),
            &only_pdb[..only_pdb.len().min(3)],
            only_cif.len(),
            &only_cif[..only_cif.len().min(3)]
        ));
    }
    let far = kp
        .iter()
        .filter(|(k, v)| {
            kc.get(*k)
                .is_some_and(|o| (**v - *o).abs().max_element() >= 1e-3)
        })
        .count();
    if far > 0 {
        out.push(format!("{id}: {far} atoms differ in position"));
    }
    for ss in [SecondaryStructure::Helix, SecondaryStructure::Strand] {
        let (a, b) = (ss_keys(tp, ss), ss_keys(tc, ss));
        if a != b {
            out.push(format!("{id}: {ss:?} residues {} vs {}", a.len(), b.len()));
        }
    }
    out
}

const REAL_ENTRIES: [&str; 7] = ["1IGY", "3NIR", "2K39", "1OKC", "2RH1", "1BNA", "1ZNI"];

/// (first-model atom rows, model count) counted straight from the text.
fn raw_pdb_counts(text: &str) -> (usize, usize) {
    let mut first = 0;
    let mut models = 0;
    for line in text.lines() {
        if line.starts_with("ENDMDL") {
            models += 1;
        }
        if models == 0 && (line.starts_with("ATOM  ") || line.starts_with("HETATM")) {
            first += 1;
        }
    }
    (first, models.max(1))
}

#[test]
#[ignore = "needs fixtures/real downloads (see docs/FORMATS.md)"]
fn real_entries_agree_between_formats_and_with_raw_row_counts() {
    for id in REAL_ENTRIES {
        let problems = diff_formats(id);
        assert!(problems.is_empty(), "{problems:#?}");
        let text = std::fs::read_to_string(real(&format!("{id}.pdb")).unwrap()).unwrap();
        let (atoms, models) = raw_pdb_counts(&text);
        let s = vv_io::load(real(&format!("{id}.pdb")).unwrap()).unwrap();
        assert_eq!(
            s.topology.atom_count(),
            atoms,
            "{id} atoms vs raw ATOM/HETATM rows"
        );
        assert_eq!(s.frame_count(), models, "{id} models vs raw ENDMDL count");
        let cif_text = std::fs::read_to_string(real(&format!("{id}.cif")).unwrap()).unwrap();
        let rows = cif_text
            .lines()
            .filter(|l| l.starts_with("ATOM ") || l.starts_with("HETATM "))
            .count();
        assert_eq!(rows, atoms * models, "{id} mmCIF atom_site rows");
    }
}

#[test]
#[ignore = "needs fixtures/real downloads (see docs/FORMATS.md)"]
fn real_entries_have_the_features_they_were_chosen_for() {
    let load_real = |n: &str| vv_io::load(real(n).unwrap()).unwrap();
    let ig = load_real("1IGY.pdb").topology;
    assert!(
        ig.residues.iter().any(|r| r.ins_code != 0),
        "antibody insertion codes"
    );
    let nir = load_real("3NIR.pdb").topology;
    assert!(nir.alt_loc.iter().any(|&a| a != 0), "alternate locations");
    assert_eq!(load_real("2K39.pdb").frame_count(), 116);
    let dna = load_real("1BNA.pdb").topology;
    assert!((0..dna.residue_count()).any(|r| dna.residue_name(r).trim() == "DC"));
    let okc = load_real("1OKC.pdb").topology;
    assert!(okc.flags.iter().any(|&f| f & flags::HETERO != 0));
}

#[test]
#[ignore = "needs fixtures/real/4V6X.cif (~29 MB)"]
fn real_large_structure_round_trips_through_pdb_with_hybrid36() {
    let original = vv_io::load(real("4V6X.cif").unwrap()).unwrap();
    assert!(original.topology.atom_count() > 99_999);
    let back = round_trip(&original, Format::Pdb);
    let (a, b) = (&original.topology, &back.topology);
    assert_eq!(a.atom_count(), b.atom_count());
    assert!(b.serial.windows(2).all(|w| w[0] < w[1]));
    assert!(b.serial.last().is_some_and(|&n| n > 99_999));
    for (pa, pb) in original
        .frame(0)
        .positions()
        .iter()
        .zip(back.frame(0).positions())
    {
        assert!((*pa - *pb).abs().max_element() < 1e-3);
    }
    let seq = |t: &Topology| t.residues.iter().map(|r| r.auth_seq_id).collect::<Vec<_>>();
    assert_eq!(seq(a), seq(b));
}

#[test]
#[ignore = "needs fixtures/real downloads (see docs/FORMATS.md)"]
fn real_entries_round_trip_through_pdb_and_mmcif_writers() {
    for id in REAL_ENTRIES {
        let original = vv_io::load(real(&format!("{id}.pdb")).unwrap()).unwrap();
        let back = round_trip(&original, Format::Pdb);
        assert_same_structure(&original, &back);
        assert_same_structure(&original, &round_trip(&original, Format::Mmcif));
    }
}
