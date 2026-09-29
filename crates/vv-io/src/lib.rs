//! File-format readers and writers, and the synthetic structure generator.
//!
//! All parsers are written from the format specifications (clean-room).

use std::io::{Read, Write};
use std::path::Path;
use std::sync::Arc;

use vv_core::fixedbitset::FixedBitSet;
use vv_core::{Element, Structure, StructureError, Topology};

mod atom_site;
pub mod bcif;
pub mod cif;
pub mod dcd;
#[cfg(feature = "fetch")]
pub mod fetch;
pub mod float;
pub mod gro_write;
pub mod mmcif;
mod mmcif_entity;
mod mmcif_entity_write;
pub mod mmcif_write;
mod msgpack;
pub mod netcdf;
pub mod pdb;
mod pdb_names;
mod pdb_seqres;
pub mod pdb_write;
mod polymer_layout;
pub mod pqr_write;
pub mod prmtop;
pub mod psf;
mod ss_range;
pub mod synth;
pub mod trajectory;
pub mod trr;
pub mod write;
pub mod xtc;
pub mod xyz_write;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Format {
    Mmcif,
    /// BinaryCIF (MessagePack-encoded PDBx).
    Bcif,
    Pdb,
}

#[derive(thiserror::Error, Debug)]
pub enum ParseError {
    #[error("cannot read {path}: {source}")]
    Io {
        path: String,
        source: std::io::Error,
    },
    #[error("cannot tell the file format of {0}: supported are mmCIF (.cif, .mmcif, .pdbx), BinaryCIF (.bcif) and PDB (.pdb, .ent), optionally gzipped; MMTF is not")]
    UnknownFormat(String),
    #[error("file contains no atoms")]
    NoAtoms,
    #[error("atom table is missing the required column {0}")]
    MissingColumn(&'static str),
    #[error("missing required section {0}")]
    MissingSection(&'static str),
    #[error("binary CIF: {0}")]
    BinaryCif(String),
    #[error("line {line}: {message}")]
    Malformed { line: usize, message: String },
    #[error(transparent)]
    Structure(#[from] StructureError),
}

/// Best guess at an atom's element from its MD topology mass, the atom
/// name as fallback (`vv_io::psf`, `vv_io::prmtop`).
pub(crate) fn guess_element(mass: f32, name: &str) -> Element {
    let by_mass = Element::from_mass(mass);
    if !by_mass.is_unknown() {
        by_mass
    } else {
        Element::from_atom_name(name.as_bytes())
    }
}

/// An MD topology's partial (force-field) charge, kept only when it is
/// within 0.1 e of an integer -- a monatomic ion's formal charge, not a
/// polarized covalent atom's -- since `Topology::charge` is a small
/// formal-charge column (as in PDB), not a float. Anything else reads as
/// neutral rather than as a misleading rounded formal charge.
pub(crate) fn round_partial_charge(charge: f32) -> i8 {
    let rounded = charge.round();
    if (charge - rounded).abs() < 0.1 {
        rounded.clamp(i8::MIN as f32, i8::MAX as f32) as i8
    } else {
        0
    }
}

/// Truncates (never wraps) a name into the 4-byte space-padded column
/// every other reader uses (`Topology::name`).
pub(crate) fn pad_name4(name: &str) -> [u8; 4] {
    let bytes = name.as_bytes();
    let n = bytes.len().min(4);
    let mut out = [b' '; 4];
    out[..n].copy_from_slice(&bytes[..n]);
    out
}

/// Loads just the topology half of a structure file, for
/// [`trajectory`]/`vv_scene::Command::LoadTrajectory`'s topology argument:
/// PSF and PRMTOP (no coordinates of their own) alongside every format
/// [`load`] reads (whose coordinates are then unused). PSF/PRMTOP are
/// deliberately not part of [`Format`]/[`load`] itself, since neither is
/// a valid single-file structure -- loading one directly (`loadstructure`)
/// still reports "unknown format".
pub fn load_topology(path: impl AsRef<Path>) -> Result<Arc<Topology>, ParseError> {
    let path = path.as_ref();
    let io = |source| ParseError::Io {
        path: path.display().to_string(),
        source,
    };
    let name = path
        .file_name()
        .and_then(|n| n.to_str())
        .map(|n| n.to_ascii_lowercase())
        .unwrap_or_default();
    let name = name.strip_suffix(".gz").unwrap_or(&name);
    let ext = name.rsplit_once('.').map(|(_, e)| e);
    match ext {
        Some("psf") => {
            let text = read_to_string(path).map_err(io)?;
            psf::parse(&text).map(Arc::new)
        }
        Some("prmtop" | "parm7") => {
            let text = read_to_string(path).map_err(io)?;
            prmtop::parse(&text).map(Arc::new)
        }
        _ => load(path).map(|s| s.topology.clone()),
    }
}

/// Reads a file (optionally gzip-compressed) as text. Separate from
/// [`load`]'s zero-copy mmap path: PSF/PRMTOP are text formats parsed with
/// `str` methods, not `load`'s binary-safe columnar parsers, and are
/// small enough next to a coordinate-bearing structure file that an owned
/// copy costs nothing that matters.
fn read_to_string(path: &Path) -> std::io::Result<String> {
    let bytes = std::fs::read(path)?;
    if bytes.starts_with(&[0x1f, 0x8b]) {
        let mut out = String::new();
        flate2::read::GzDecoder::new(&bytes[..]).read_to_string(&mut out)?;
        Ok(out)
    } else {
        String::from_utf8(bytes)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))
    }
}

impl Format {
    /// Guesses the format from the file name, ignoring a `.gz` suffix.
    pub fn from_path(path: &Path) -> Option<Format> {
        let name = path.file_name()?.to_str()?.to_ascii_lowercase();
        let name = name.strip_suffix(".gz").unwrap_or(&name);
        let ext = name.rsplit_once('.')?.1;
        match ext {
            "cif" | "mmcif" | "pdbx" => Some(Format::Mmcif),
            "bcif" => Some(Format::Bcif),
            "pdb" | "ent" => Some(Format::Pdb),
            _ => None,
        }
    }

    /// Guesses the format from the first bytes.
    pub fn sniff(bytes: &[u8]) -> Option<Format> {
        let head = &bytes[..bytes.len().min(4096)];
        if head.starts_with(b"data_") || head.windows(5).any(|w| w == b"loop_") {
            Some(Format::Mmcif)
        } else if head.windows(10).any(|w| w == b"dataBlocks") {
            Some(Format::Bcif)
        } else if head.starts_with(b"HEADER")
            || head.starts_with(b"ATOM")
            || head.starts_with(b"HETATM")
            || head.starts_with(b"REMARK")
            || head.starts_with(b"MODEL")
        {
            Some(Format::Pdb)
        } else {
            None
        }
    }
}

pub fn parse(bytes: &[u8], format: Format) -> Result<Structure, ParseError> {
    match format {
        Format::Mmcif => mmcif::parse(bytes),
        Format::Bcif => bcif::parse(bytes),
        Format::Pdb => pdb::parse(bytes),
    }
}

/// Reads a structure file (optionally gzip-compressed), choosing the
/// format from the extension or, failing that, the content.
pub fn load(path: impl AsRef<Path>) -> Result<Structure, ParseError> {
    let path = path.as_ref();
    let io = |source| ParseError::Io {
        path: path.display().to_string(),
        source,
    };
    let file = std::fs::File::open(path).map_err(io)?;
    let map = unsafe { memmap2::Mmap::map(&file) }.map_err(io)?;
    let bytes: &[u8] = &map;
    let decompressed;
    let bytes = if bytes.starts_with(&[0x1f, 0x8b]) {
        let mut out = Vec::with_capacity(bytes.len() * 4);
        flate2::read::GzDecoder::new(bytes)
            .read_to_end(&mut out)
            .map_err(io)?;
        decompressed = out;
        &decompressed[..]
    } else {
        bytes
    };
    let format = Format::from_path(path)
        .or_else(|| Format::sniff(bytes))
        .ok_or_else(|| ParseError::UnknownFormat(path.display().to_string()))?;
    parse(bytes, format)
}

/// A format `vv_io::save` can write.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SaveFormat {
    Pdb,
    Mmcif,
    Xyz,
    Pqr,
    Gro,
}

impl SaveFormat {
    /// Guesses the format from the file name, ignoring a `.gz` suffix.
    pub fn from_path(path: &Path) -> Option<SaveFormat> {
        let name = path.file_name()?.to_str()?.to_ascii_lowercase();
        let name = name.strip_suffix(".gz").unwrap_or(&name);
        let ext = name.rsplit_once('.')?.1;
        match ext {
            "cif" | "mmcif" | "pdbx" => Some(SaveFormat::Mmcif),
            "pdb" | "ent" => Some(SaveFormat::Pdb),
            "xyz" => Some(SaveFormat::Xyz),
            "pqr" => Some(SaveFormat::Pqr),
            "gro" => Some(SaveFormat::Gro),
            _ => None,
        }
    }
}

#[derive(thiserror::Error, Debug)]
pub enum SaveError {
    #[error("cannot tell the output format of {0}: supported are .pdb/.ent, .cif/.mmcif/.pdbx, .xyz, .pqr, .gro, optionally .gz")]
    UnknownFormat(String),
    #[error("cannot write {path}: {source}")]
    Io {
        path: String,
        source: std::io::Error,
    },
}

/// What to write: a subset of atoms (`None` is every atom) and which
/// coordinate sets, MODEL/ENDMDL- or `pdbx_PDB_model_num`-wrapped when
/// there is more than one (PDB/mmCIF); the other formats just repeat the
/// per-frame block.
pub struct SaveOptions<'a> {
    pub atoms: Option<&'a FixedBitSet>,
    pub frames: &'a [usize],
}

fn write_format(
    structure: &Structure,
    format: SaveFormat,
    opts: &SaveOptions,
    out: &mut impl Write,
) -> std::io::Result<Vec<String>> {
    match format {
        SaveFormat::Pdb => pdb_write::write(structure, opts.atoms, opts.frames, out),
        SaveFormat::Mmcif => mmcif_write::write(structure, opts.atoms, opts.frames, out),
        SaveFormat::Xyz => xyz_write::write(structure, opts.atoms, opts.frames, out),
        SaveFormat::Pqr => pqr_write::write(structure, opts.atoms, opts.frames, out),
        SaveFormat::Gro => gro_write::write(structure, opts.atoms, opts.frames, out),
    }
}

/// Writes a structure file, choosing the format from `path`'s extension
/// (optionally gzip-compressed when it ends in `.gz`). Returns one
/// warning per thing the format couldn't represent exactly (e.g. a PDB
/// chain id that had to be remapped).
pub fn save(
    structure: &Structure,
    path: impl AsRef<Path>,
    opts: &SaveOptions,
) -> Result<Vec<String>, SaveError> {
    let path = path.as_ref();
    let format = SaveFormat::from_path(path)
        .ok_or_else(|| SaveError::UnknownFormat(path.display().to_string()))?;
    let io_err = |source| SaveError::Io {
        path: path.display().to_string(),
        source,
    };
    let file = std::fs::File::create(path).map_err(io_err)?;
    let buffered = std::io::BufWriter::new(file);
    let gz = path
        .file_name()
        .and_then(|n| n.to_str())
        .is_some_and(|n| n.to_ascii_lowercase().ends_with(".gz"));
    if gz {
        let mut enc = flate2::write::GzEncoder::new(buffered, flate2::Compression::default());
        let warnings = write_format(structure, format, opts, &mut enc).map_err(io_err)?;
        enc.finish().map_err(io_err)?;
        Ok(warnings)
    } else {
        let mut buffered = buffered;
        let warnings = write_format(structure, format, opts, &mut buffered).map_err(io_err)?;
        std::io::Write::flush(&mut buffered).map_err(io_err)?;
        Ok(warnings)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn format_from_extension_and_content() {
        assert_eq!(
            Format::from_path(Path::new("x/1abc.cif")),
            Some(Format::Mmcif)
        );
        assert_eq!(
            Format::from_path(Path::new("1abc.CIF.gz")),
            Some(Format::Mmcif)
        );
        assert_eq!(Format::from_path(Path::new("1abc.pdb")), Some(Format::Pdb));
        assert_eq!(
            Format::from_path(Path::new("pdb1abc.ent.gz")),
            Some(Format::Pdb)
        );
        assert_eq!(Format::from_path(Path::new("notes.txt")), None);
        assert_eq!(Format::sniff(b"data_1ABC\n#"), Some(Format::Mmcif));
        assert_eq!(Format::sniff(b"HEADER    X"), Some(Format::Pdb));
        assert_eq!(Format::sniff(b"hello"), None);
    }
}
