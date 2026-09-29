//! `cargo xtask <command>` — repository automation.

use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::time::Instant;

use vv_io::synth::{protein_like, SynthParams};

const USAGE: &str = "\
usage: cargo xtask <command>

commands:
  fetch [ID...]              download large benchmark structures (default 4V6X 3J3Q)
                             from RCSB into fixtures/large/ as .cif.gz
  synth [--atoms N] [--out PATH]
                             write a synthetic protein-density mmCIF
                             (default 10000000 atoms, fixtures/synth/synth<N>.cif)
  parse FILE...              time loading + bond perception of structure files
  help                       show this message";

const LARGE: [&str; 2] = ["4V6X", "3J3Q"];

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn fetch(ids: &[String]) -> Result<(), String> {
    let dir = repo_root().join("fixtures/large");
    std::fs::create_dir_all(&dir).map_err(|e| e.to_string())?;
    let ids: Vec<String> = if ids.is_empty() {
        LARGE.iter().map(|s| s.to_string()).collect()
    } else {
        ids.to_vec()
    };
    for id in ids {
        let id = id.to_ascii_uppercase();
        let path = dir.join(format!("{id}.cif.gz"));
        if path.exists() {
            println!("{} already present", path.display());
            continue;
        }
        let url = format!("https://files.rcsb.org/download/{id}.cif.gz");
        print!("downloading {url} ... ");
        std::io::stdout().flush().ok();
        let response = ureq::get(&url).call().map_err(|e| format!("{url}: {e}"))?;
        let mut bytes = Vec::new();
        response
            .into_body()
            .into_reader()
            .read_to_end(&mut bytes)
            .map_err(|e| format!("{url}: {e}"))?;
        std::fs::write(&path, &bytes).map_err(|e| e.to_string())?;
        println!("{} bytes -> {}", bytes.len(), path.display());
    }
    Ok(())
}

fn synth(args: &[String]) -> Result<(), String> {
    let mut atoms = 10_000_000usize;
    let mut out: Option<PathBuf> = None;
    let mut it = args.iter();
    while let Some(flag) = it.next() {
        match flag.as_str() {
            "--atoms" => {
                atoms = it
                    .next()
                    .ok_or("--atoms needs a value")?
                    .parse()
                    .map_err(|e| format!("{e}"))?
            }
            "--out" => out = Some(PathBuf::from(it.next().ok_or("--out needs a value")?)),
            other => return Err(format!("unknown option `{other}`")),
        }
    }
    let path = out.unwrap_or_else(|| repo_root().join(format!("fixtures/synth/synth{atoms}.cif")));
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent).map_err(|e| e.to_string())?;
    }
    let t = Instant::now();
    let structure = protein_like(&SynthParams::new(atoms));
    println!(
        "generated {} atoms in {:.2?}",
        structure.atom_count(),
        t.elapsed()
    );
    let t = Instant::now();
    let file = std::fs::File::create(&path).map_err(|e| e.to_string())?;
    let mut writer = std::io::BufWriter::with_capacity(1 << 20, file);
    vv_io::mmcif_write::write(&structure, None, &[0], &mut writer).map_err(|e| e.to_string())?;
    writer.flush().map_err(|e| e.to_string())?;
    let size = std::fs::metadata(&path).map(|m| m.len()).unwrap_or(0);
    println!(
        "wrote {} ({} MB) in {:.2?}",
        path.display(),
        size >> 20,
        t.elapsed()
    );
    Ok(())
}

fn parse(files: &[String]) -> Result<(), String> {
    if files.is_empty() {
        return Err("parse needs at least one file".into());
    }
    for file in files {
        let path = Path::new(file);
        let t = Instant::now();
        let structure = vv_io::load(path).map_err(|e| e.to_string())?;
        let load = t.elapsed();
        let t = Instant::now();
        let bonds = vv_core::bonds::perceive(&structure.topology, structure.frame(0).positions());
        let perceive = t.elapsed();
        let top = &structure.topology;
        println!(
            "{}: {} atoms, {} residues, {} chains, {} frames | load {:.2?} ({:.1} M atoms/s) | bonds {} in {:.2?}",
            path.display(),
            top.atom_count(),
            top.residue_count(),
            top.chain_count(),
            structure.frame_count(),
            load,
            top.atom_count() as f64 / load.as_secs_f64() / 1e6,
            bonds.len(),
            perceive
        );
        if let Err(e) = top.validate() {
            println!("  INVALID topology: {e}");
        }
    }
    Ok(())
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let (command, rest) = match args.split_first() {
        Some((c, r)) => (c.as_str(), r),
        None => ("help", &[][..]),
    };
    let result = match command {
        "fetch" => fetch(rest),
        "synth" => synth(rest),
        "parse" => parse(rest),
        "help" | "--help" | "-h" => {
            println!("{USAGE}");
            Ok(())
        }
        other => Err(format!("unknown command `{other}`\n\n{USAGE}")),
    };
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("{e}");
            ExitCode::FAILURE
        }
    }
}
