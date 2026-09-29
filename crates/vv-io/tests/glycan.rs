//! `vv_core::glycan` against a real N-glycosylated structure: ground
//! truth below is what `vv_core::bonds::perceive` +
//! `vv_core::glycan::GlycanPlan` find on the real coordinates, not a
//! read of the file's own `LINK` records.

use std::time::Instant;

use vv_core::glam::Vec3;
use vv_core::{
    bonds, build_glycan_mesh, Attachment, GlycanFrame, GlycanPlan, PolytopeMesh, Shape, SnfgColor,
    Structure,
};

fn fixture_6x3z() -> Structure {
    let path =
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/glycan/6X3Z.pdb");
    vv_io::load(path).unwrap()
}

fn plan_6x3z() -> (Structure, GlycanPlan) {
    let structure = fixture_6x3z();
    let coords = structure.frame(0);
    let bonds = bonds::perceive(&structure.topology, coords.positions());
    let plan = GlycanPlan::build(&structure.topology, &bonds);
    (structure, plan)
}

/// 6X3Z (a GABA-A receptor cryo-EM structure) carries 11 NAG, 7 MAN and
/// 3 BMA residues -- 21 glycan residues total, each named by its own PDB
/// CCD code (the file's own `HETATM` records, not GLYCAM or CHARMM
/// names: the "common/PDB CCD" column of the source script's table).
/// 154 NAG *atoms* / ~14 atoms per NAG residue is exactly 11 residues.
#[test]
fn detects_every_glycan_residue_in_6x3z_by_name() {
    let (_structure, plan) = plan_6x3z();
    assert_eq!(plan.residues.len(), 21, "11 NAG + 7 MAN + 3 BMA");

    let by_shape = |shape: Shape| plan.residues.iter().filter(|r| r.shape == shape).count();
    assert_eq!(by_shape(Shape::Cube), 11, "NAG -> GlcNAc, blue cube");
    assert_eq!(
        by_shape(Shape::Sphere),
        10,
        "7 MAN + 3 BMA -> Man, green sphere"
    );
    assert_eq!(
        plan.residues
            .iter()
            .filter(|r| r.shape == Shape::Hexagon && r.label == "Unknown glycan")
            .count(),
        0,
        "every residue in this fixture is a recognized name"
    );
}

/// Every N-glycosylation site's innermost NAG links to its ASN (7 sites,
/// 7 `ProteinCa` attachments); every other residue links to the sugar
/// one step closer to the protein (`Residue`); none are free-standing
/// (`None`) or capped (`Terminal`) -- this fixture has no O-linked
/// glycans, GLYCAM caps or unlinked monosaccharides.
#[test]
fn every_residue_resolves_a_linkage() {
    let (_structure, plan) = plan_6x3z();
    let protein_ca = plan
        .residues
        .iter()
        .filter(|r| matches!(r.attachment, Attachment::ProteinCa(_)))
        .count();
    let residue_to_residue = plan
        .residues
        .iter()
        .filter(|r| matches!(r.attachment, Attachment::Residue(_)))
        .count();
    let unresolved = plan
        .residues
        .iter()
        .filter(|r| matches!(r.attachment, Attachment::None | Attachment::Terminal(_)))
        .count();
    assert_eq!(protein_ca, 7, "one N-glycosylation site per root NAG");
    assert_eq!(residue_to_residue, 14, "21 total - 7 protein-linked roots");
    assert_eq!(unresolved, 0);
}

/// A specific residue, chosen from the file's own numbering (chain A,
/// residue 404, the site on ASN A 80 per the file's `LINK` record) is a
/// NAG resolved as the blue GlcNAc cube and linked straight to a protein
/// alpha carbon, not to another sugar.
#[test]
fn a_named_root_nag_is_a_blue_cube_linked_to_the_protein() {
    let (structure, plan) = plan_6x3z();
    let root = plan
        .residues
        .iter()
        .find(|r| {
            let rec = &structure.topology.residues[r.residue as usize];
            structure.topology.residue_name(r.residue as usize) == "NAG" && rec.auth_seq_id == 404
        })
        .expect("residue 404 should be a detected NAG");
    assert_eq!(root.shape, Shape::Cube);
    assert_eq!(root.color1, SnfgColor::Blue);
    assert!(matches!(root.attachment, Attachment::ProteinCa(_)));
}

/// `(chain1, seq1, atom1, chain2, seq2, atom2)` per LINK line.
fn parse_link_records(raw: &str) -> Vec<(&str, i32, &str, &str, i32, &str)> {
    raw.lines()
        .filter(|l| l.starts_with("LINK"))
        .map(|line| {
            let tok: Vec<&str> = line.split_whitespace().collect();
            // LINK atom1 res1 chain1 seq1 atom2 res2 chain2 seq2 sym1 sym2 dist
            (
                tok[3],
                tok[4].parse().unwrap(),
                tok[1],
                tok[7],
                tok[8].parse().unwrap(),
                tok[5],
            )
        })
        .collect()
}

/// The 21 `LINK` records 6X3Z's raw PDB text carries (7 ASN ND2-C1,
/// 14 sugar O*-C1), parsed directly from the file here -- independent
/// of `vv_io::pdb`'s own LINK reader, since the file's 303 `CONECT`
/// records also encode some of this connectivity and would make
/// `topology.explicit_bonds` an unreliable ground truth for "exactly
/// these 21". Each one is a perceived bond, and with `explicit_bonds`
/// cleared, `vv_core::bonds`'s geometric glycosidic rule (glycan
/// anomeric carbon to an acceptor oxygen or an N-glycosylation ND2)
/// recovers every one of them from coordinates alone.
#[test]
fn every_link_record_is_a_perceived_bond_and_the_geometric_rule_alone_recovers_them() {
    let path =
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../fixtures/glycan/6X3Z.pdb");
    let raw = std::fs::read_to_string(&path).unwrap();
    let links = parse_link_records(&raw);
    assert_eq!(links.len(), 21, "7 ASN-NAG + 14 sugar-sugar LINK records");

    let structure = fixture_6x3z();
    let t = &structure.topology;
    let coords = structure.frame(0);
    let p = coords.positions();

    let find_atom = |chain: &str, seq: i32, name: &str| -> u32 {
        t.residues
            .iter()
            .find(|r| r.auth_seq_id == seq && t.chain_name(r.chain as usize) == chain)
            .and_then(|r| r.atoms.clone().find(|&a| t.atom_name(a as usize) == name))
            .unwrap_or_else(|| panic!("no atom {name} in {chain}/{seq}"))
    };
    let link_pairs: Vec<[u32; 2]> = links
        .iter()
        .map(|&(chain1, seq1, atom1, chain2, seq2, atom2)| {
            [
                find_atom(chain1, seq1, atom1),
                find_atom(chain2, seq2, atom2),
            ]
        })
        .collect();

    let bonds = bonds::perceive(t, p);
    for &[a, b] in &link_pairs {
        assert!(
            bonds.contains(a, b),
            "LINK {a}-{b} missing from perceived bonds"
        );
    }

    let mut no_links = structure.clone();
    let mut topology = no_links.topology.as_ref().clone();
    topology.explicit_bonds.clear();
    no_links.topology = std::sync::Arc::new(topology);
    let geometric_only = bonds::perceive(&no_links.topology, p);
    for &[a, b] in &link_pairs {
        assert!(
            geometric_only.contains(a, b),
            "geometric glycosidic rule alone should recover {a}-{b}"
        );
    }
}

/// The per-frame path (`GlycanPlan::update_into` + `build_glycan_mesh`)
/// at a realistic few-hundred-residue scale: replicate 6X3Z's 21-residue
/// glycan tree (ring atoms and linkages, offset apart so none collide)
/// enough times to reach ~300 residues, matching a heavily glycosylated
/// structure. Loose bound (10x the task's 1ms target) so a debug-profile
/// CI run doesn't flake; the release-profile number is what matters and
/// is printed for the record.
#[test]
fn per_frame_update_is_fast_at_a_few_hundred_residues() {
    let (structure, plan) = plan_6x3z();
    let base_positions = structure.frame(0).positions().to_vec();
    let repeats = 300 / plan.residues.len().max(1) + 1;

    let mut topology = structure.topology.as_ref().clone();
    let mut positions = Vec::new();
    let atoms_per_copy = topology.atom_count() as u32;
    let residues_per_copy = topology.residues.len() as u32;
    let chains_per_copy = topology.chains.len() as u32;
    for copy in 1..repeats {
        let offset = Vec3::new(copy as f32 * 200.0, 0.0, 0.0);
        positions.extend(base_positions.iter().map(|&p| p + offset));
        for rec in structure.topology.residues.clone() {
            let mut rec = rec;
            rec.atoms = (rec.atoms.start + copy as u32 * atoms_per_copy)
                ..(rec.atoms.end + copy as u32 * atoms_per_copy);
            rec.chain += copy as u32 * chains_per_copy;
            topology.residues.push(rec);
        }
        for rec in structure.topology.chains.clone() {
            let mut rec = rec;
            rec.residues = (rec.residues.start + copy as u32 * residues_per_copy)
                ..(rec.residues.end + copy as u32 * residues_per_copy);
            topology.chains.push(rec);
        }
        topology
            .element
            .extend(structure.topology.element.iter().copied());
        topology
            .name
            .extend(structure.topology.name.iter().copied());
        topology
            .flags
            .extend(structure.topology.flags.iter().copied());
        topology.residue_index.extend(
            structure
                .topology
                .residue_index
                .iter()
                .map(|&r| r + copy as u32 * residues_per_copy),
        );
    }
    positions.splice(0..0, base_positions.iter().copied());

    let bonds = bonds::perceive(&topology, &positions);
    let big_plan = GlycanPlan::build(&topology, &bonds);
    assert!(
        big_plan.residues.len() >= 300,
        "expected >= 300 residues, got {}",
        big_plan.residues.len()
    );

    let mut frame = GlycanFrame::default();
    let mut mesh = PolytopeMesh::default();
    // Warm up (first call sizes the reused buffers).
    big_plan.update_into(&positions, &mut frame);
    build_glycan_mesh(&big_plan, &frame, 4.0, |_| true, &mut mesh);

    const ITERS: usize = 200;
    let start = Instant::now();
    for _ in 0..ITERS {
        big_plan.update_into(&positions, &mut frame);
        build_glycan_mesh(&big_plan, &frame, 4.0, |_| true, &mut mesh);
    }
    let per_frame = start.elapsed() / ITERS as u32;
    eprintln!(
        "glycan per-frame update+mesh: {:?} for {} residues ({:?}/residue)",
        per_frame,
        big_plan.residues.len(),
        per_frame / big_plan.residues.len() as u32
    );
    assert!(
        per_frame.as_micros() < 10_000,
        "per-frame update+mesh took {per_frame:?}, expected well under 10ms even in debug"
    );
}
