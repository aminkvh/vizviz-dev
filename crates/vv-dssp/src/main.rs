//! Standalone secondary-structure assignment: a thin CLI over
//! `vv_core::dssp` (Kabsch & Sander 1983), for anyone who wants the
//! kernel without the rest of vizviz. vizviz itself calls the same
//! function in-process; this binary exists so the algorithm is not
//! locked to being a vizviz user; it is a separate `[[bin]]` rather than
//! vizviz shelling out to itself.

use std::process::ExitCode;

use vv_core::dssp::{assign, DsspCode};

const USAGE: &str = "\
usage: vv-dssp FILE [--frame N]
  FILE      structure to read (.cif/.mmcif/.pdb, optionally .gz)
  --frame N   which model/frame to assign from (default: 0)

Prints one line per protein residue: chain, author residue number,
residue name, and the DSSP code (H alpha helix, G 3-10 helix, I pi
helix, E strand, B isolated bridge, T turn, S bend, - coil). Residues
outside the protein backbone (ligands, water, nucleic acids) are
skipped, not printed as coil.";

fn main() -> ExitCode {
    let mut file = None;
    let mut frame = 0usize;
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "-h" | "--help" => {
                println!("{USAGE}");
                return ExitCode::SUCCESS;
            }
            "--frame" => match args.next().and_then(|v| v.parse().ok()) {
                Some(n) => frame = n,
                None => {
                    eprintln!("--frame needs a number\n\n{USAGE}");
                    return ExitCode::from(2);
                }
            },
            other if other.starts_with('-') => {
                eprintln!("unknown option `{other}`\n\n{USAGE}");
                return ExitCode::from(2);
            }
            path if file.is_none() => file = Some(path.to_owned()),
            _ => {
                eprintln!("only one FILE at a time\n\n{USAGE}");
                return ExitCode::from(2);
            }
        }
    }
    let Some(file) = file else {
        eprintln!("{USAGE}");
        return ExitCode::from(2);
    };

    let structure = match vv_io::load(&file) {
        Ok(s) => s,
        Err(e) => {
            eprintln!("{file}: {e}");
            return ExitCode::FAILURE;
        }
    };
    if frame >= structure.frame_count() {
        eprintln!(
            "{file}: frame {frame} out of range (0..{})",
            structure.frame_count()
        );
        return ExitCode::FAILURE;
    }
    let coords = structure.frame(frame);
    let codes = assign(&structure.topology, coords.positions());

    println!("CHAIN\tSEQID\tRES\tSS");
    for (i, rec) in structure.topology.residues.iter().enumerate() {
        if codes[i] == DsspCode::None {
            continue;
        }
        let chain = structure
            .topology
            .names
            .get(structure.topology.chains[rec.chain as usize].auth_asym);
        println!(
            "{chain}\t{}\t{}\t{}",
            rec.auth_seq_id,
            structure.topology.residue_name(i),
            codes[i].letter()
        );
    }
    ExitCode::SUCCESS
}
