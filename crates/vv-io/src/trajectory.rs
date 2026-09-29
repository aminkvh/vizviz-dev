//! Every trajectory format behind one reader: DCD, XTC, TRR and AMBER
//! NetCDF, chosen by extension or, failing that, by the file's first
//! bytes -- plus a PDB or mmCIF structure file, each of its models one
//! frame, for a topology (typically PSF or PRMTOP) that has no
//! trajectory of its own yet, only a coordinate file. All read any frame
//! directly, without decoding the ones before it, except the structure
//! file, whose frames are ordinary in-memory models.

use std::path::Path;

use vv_core::glam::Vec3;
use vv_core::Structure;

use crate::dcd::{DcdError, DcdReader};
use crate::netcdf::{NetCdfError, NetCdfReader};
use crate::trr::{TrrError, TrrReader};
use crate::xtc::{XtcError, XtcReader};

#[derive(thiserror::Error, Debug)]
pub enum TrajectoryError {
    #[error(transparent)]
    Dcd(#[from] DcdError),
    #[error(transparent)]
    Xtc(#[from] XtcError),
    #[error(transparent)]
    Trr(#[from] TrrError),
    #[error(transparent)]
    NetCdf(#[from] NetCdfError),
    #[error(transparent)]
    Structure(#[from] crate::ParseError),
    #[error("cannot read {0}: not a DCD, XTC, TRR, NetCDF, PDB or mmCIF trajectory")]
    UnknownFormat(String),
}

pub enum Trajectory {
    Dcd(DcdReader),
    Xtc(XtcReader),
    Trr(TrrReader),
    NetCdf(NetCdfReader),
    /// A structure file's own models, read as coordinate frames.
    Structure(Structure),
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Kind {
    Dcd,
    Xtc,
    Trr,
    NetCdf,
    Structure,
}

fn kind_by_extension(path: &Path) -> Option<Kind> {
    let ext = path.extension()?.to_str()?.to_ascii_lowercase();
    Some(match ext.as_str() {
        "dcd" => Kind::Dcd,
        "xtc" => Kind::Xtc,
        "trr" => Kind::Trr,
        "nc" | "ncdf" | "netcdf" => Kind::NetCdf,
        _ if crate::Format::from_path(path).is_some() => Kind::Structure,
        _ => return None,
    })
}

fn kind_by_content(head: &[u8]) -> Option<Kind> {
    let be = |at: usize| {
        head.get(at..at + 4)
            .map(|b| i32::from_be_bytes([b[0], b[1], b[2], b[3]]))
    };
    if head.get(4..8) == Some(b"CORD") {
        Some(Kind::Dcd)
    } else if be(0) == Some(1995) {
        Some(Kind::Xtc)
    } else if be(0) == Some(1993) {
        Some(Kind::Trr)
    } else if head.starts_with(b"CDF") || head.starts_with(b"\x89HDF") {
        Some(Kind::NetCdf)
    } else if crate::Format::sniff(head).is_some() {
        Some(Kind::Structure)
    } else {
        None
    }
}

impl Trajectory {
    pub fn open(path: impl AsRef<Path>) -> Result<Self, TrajectoryError> {
        let path = path.as_ref();
        let kind = match kind_by_extension(path) {
            Some(kind) => kind,
            None => {
                let mut head = [0u8; 4096];
                let read = std::fs::File::open(path)
                    .and_then(|mut f| std::io::Read::read(&mut f, &mut head))
                    .unwrap_or(0);
                kind_by_content(&head[..read])
                    .ok_or_else(|| TrajectoryError::UnknownFormat(path.display().to_string()))?
            }
        };
        Ok(match kind {
            Kind::Dcd => Trajectory::Dcd(DcdReader::open(path)?),
            Kind::Xtc => Trajectory::Xtc(XtcReader::open(path)?),
            Kind::Trr => Trajectory::Trr(TrrReader::open(path)?),
            Kind::NetCdf => Trajectory::NetCdf(NetCdfReader::open(path)?),
            Kind::Structure => Trajectory::Structure(crate::load(path)?),
        })
    }

    pub fn atom_count(&self) -> usize {
        match self {
            Trajectory::Dcd(r) => r.atom_count(),
            Trajectory::Xtc(r) => r.atom_count(),
            Trajectory::Trr(r) => r.atom_count(),
            Trajectory::NetCdf(r) => r.atom_count(),
            Trajectory::Structure(s) => s.atom_count(),
        }
    }

    pub fn frame_count(&self) -> usize {
        match self {
            Trajectory::Dcd(r) => r.frame_count(),
            Trajectory::Xtc(r) => r.frame_count(),
            Trajectory::Trr(r) => r.frame_count(),
            Trajectory::NetCdf(r) => r.frame_count(),
            Trajectory::Structure(s) => s.frame_count(),
        }
    }

    /// Frame `index`'s positions, in Å.
    pub fn read_frame(&self, index: usize) -> Result<Vec<Vec3>, TrajectoryError> {
        Ok(match self {
            Trajectory::Dcd(r) => r.read_frame(index)?,
            Trajectory::Xtc(r) => r.read_frame(index)?,
            Trajectory::Trr(r) => r.read_frame(index)?,
            Trajectory::NetCdf(r) => r.read_frame(index)?,
            Trajectory::Structure(s) => s.frame(index).positions().to_vec(),
        })
    }
}

/// A trajectory streams a structure's frames (`Structure::streamed`).
impl vv_core::FrameSource for Trajectory {
    fn frame_count(&self) -> usize {
        Trajectory::frame_count(self)
    }

    fn read(&self, index: usize) -> Result<Vec<Vec3>, String> {
        self.read_frame(index).map_err(|e| e.to_string())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn formats_are_told_apart_by_their_first_bytes() {
        let mut dcd = vec![0u8; 8];
        dcd[4..8].copy_from_slice(b"CORD");
        assert_eq!(kind_by_content(&dcd), Some(Kind::Dcd));
        assert_eq!(kind_by_content(&1995i32.to_be_bytes()), Some(Kind::Xtc));
        assert_eq!(kind_by_content(&1993i32.to_be_bytes()), Some(Kind::Trr));
        assert_eq!(kind_by_content(b"CDF\x02"), Some(Kind::NetCdf));
        // A PDB (or mmCIF) file is a valid trajectory argument too: its
        // models are read as coordinate frames.
        assert_eq!(kind_by_content(b"HEADER  "), Some(Kind::Structure));
        assert_eq!(kind_by_content(b"data_1ABC\n#"), Some(Kind::Structure));
        assert_eq!(kind_by_content(b"not a recognized format"), None);
    }

    #[test]
    fn extension_recognizes_a_structure_file_as_a_trajectory_source() {
        assert_eq!(
            kind_by_extension(Path::new("coords.pdb")),
            Some(Kind::Structure)
        );
        assert_eq!(
            kind_by_extension(Path::new("coords.cif.gz")),
            Some(Kind::Structure)
        );
        assert_eq!(kind_by_extension(Path::new("traj.dcd")), Some(Kind::Dcd));
        assert_eq!(kind_by_extension(Path::new("notes.txt")), None);
    }
}
