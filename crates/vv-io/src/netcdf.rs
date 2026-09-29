//! AMBER NetCDF trajectories (`.nc`, `.ncdf`): the classic NetCDF format
//! (CDF-1, or CDF-2 with 64-bit offsets) with the AMBER convention's
//! `coordinates` variable, frames x atoms x 3, in Å. Every frame is one
//! record of a fixed size, so any frame is read directly. NetCDF-4 files
//! (HDF5 underneath) are not read.
//!
//! Header grammar from the NetCDF classic format specification;
//! `vv-io/tests/traj.rs` checks it against MDAnalysis on
//! `fixtures/traj/1crn.nc`.

use std::io;
use std::path::Path;

use memmap2::Mmap;
use vv_core::glam::Vec3;

const DIMENSION: u32 = 0x0A;
const VARIABLE: u32 = 0x0B;
const ATTRIBUTE: u32 = 0x0C;
const FLOAT: u32 = 5;
const DOUBLE: u32 = 6;
const STREAMING: u32 = u32::MAX;

#[derive(thiserror::Error, Debug)]
pub enum NetCdfError {
    #[error("cannot read {path}: {source}")]
    Io { path: String, source: io::Error },
    #[error("not a classic NetCDF file")]
    NotNetCdf,
    #[error(
        "NetCDF-4 (HDF5) files are not supported; convert with `cpptraj` or `nccopy -k classic`"
    )]
    Hdf5,
    #[error("no `coordinates` variable (frame, atom, spatial) of floats or doubles")]
    NoCoordinates,
    #[error("truncated or malformed NetCDF file: {0}")]
    Malformed(&'static str),
    #[error("frame {index} is out of range (this trajectory has {frame_count})")]
    FrameOutOfRange { index: usize, frame_count: usize },
}

pub struct NetCdfReader {
    mmap: Mmap,
    atom_count: usize,
    frame_count: usize,
    /// Offset of frame 0's coordinates, and the stride between frames.
    begin: usize,
    record_size: usize,
    double: bool,
    scale: f32,
}

struct Header<'a> {
    bytes: &'a [u8],
    at: usize,
}

impl Header<'_> {
    fn u32(&mut self) -> Result<u32, NetCdfError> {
        let b = self
            .bytes
            .get(self.at..self.at + 4)
            .ok_or(NetCdfError::Malformed(
                "header runs past the end of the file",
            ))?;
        self.at += 4;
        Ok(u32::from_be_bytes([b[0], b[1], b[2], b[3]]))
    }

    fn u64(&mut self) -> Result<u64, NetCdfError> {
        Ok((self.u32()? as u64) << 32 | self.u32()? as u64)
    }

    fn name(&mut self) -> Result<&str, NetCdfError> {
        let len = self.u32()? as usize;
        let b = self
            .bytes
            .get(self.at..self.at + len)
            .ok_or(NetCdfError::Malformed("name runs past the end"))?;
        self.at += len.div_ceil(4) * 4;
        std::str::from_utf8(b).map_err(|_| NetCdfError::Malformed("name is not UTF-8"))
    }

    /// A list's tag and length; an absent list is two zeros.
    fn list(&mut self, tag: u32) -> Result<usize, NetCdfError> {
        let t = self.u32()?;
        let n = self.u32()? as usize;
        if t != tag && !(t == 0 && n == 0) {
            return Err(NetCdfError::Malformed("unexpected list tag"));
        }
        Ok(n)
    }

    /// Reads an attribute list; returns the value of a numeric
    /// `scale_factor` in it, if any.
    fn attributes(&mut self) -> Result<Option<f32>, NetCdfError> {
        let mut scale = None;
        for _ in 0..self.list(ATTRIBUTE)? {
            let is_scale = self.name()? == "scale_factor";
            let kind = self.u32()?;
            let n = self.u32()? as usize;
            let size = type_size(kind)?;
            let start = self.at;
            self.at += (n * size).div_ceil(4) * 4;
            if is_scale && n >= 1 {
                let b = &self.bytes[start..start + size];
                scale = match kind {
                    FLOAT => Some(f32::from_be_bytes(b.try_into().expect("4 bytes"))),
                    DOUBLE => Some(f64::from_be_bytes(b.try_into().expect("8 bytes")) as f32),
                    _ => None,
                };
            }
        }
        Ok(scale)
    }
}

fn type_size(kind: u32) -> Result<usize, NetCdfError> {
    Ok(match kind {
        1 | 2 => 1,
        3 => 2,
        4 | FLOAT => 4,
        DOUBLE => 8,
        _ => return Err(NetCdfError::Malformed("unknown value type")),
    })
}

impl NetCdfReader {
    pub fn open(path: impl AsRef<Path>) -> Result<Self, NetCdfError> {
        let path = path.as_ref();
        let io_err = |source| NetCdfError::Io {
            path: path.display().to_string(),
            source,
        };
        let file = std::fs::File::open(path).map_err(io_err)?;
        // SAFETY: read-only mapping, as `dcd`.
        let mmap = unsafe { Mmap::map(&file) }.map_err(io_err)?;
        if mmap.starts_with(b"\x89HDF") {
            return Err(NetCdfError::Hdf5);
        }
        let wide = match mmap.get(..4) {
            Some(b"CDF\x01") => false,
            Some(b"CDF\x02") => true,
            _ => return Err(NetCdfError::NotNetCdf),
        };
        let mut h = Header {
            bytes: &mmap,
            at: 4,
        };
        let numrecs = h.u32()?;
        let mut dims = Vec::new();
        for _ in 0..h.list(DIMENSION)? {
            h.name()?;
            dims.push(h.u32()? as usize);
        }
        h.attributes()?;
        // Record variables' sizes add up to one record.
        let mut record_size = 0;
        let mut record_vars = 0;
        let mut coordinates = None;
        for _ in 0..h.list(VARIABLE)? {
            let name = h.name()?.to_owned();
            let ids: Vec<usize> = (0..h.u32()?)
                .map(|_| Ok(h.u32()? as usize))
                .collect::<Result<_, NetCdfError>>()?;
            let scale = h.attributes()?;
            let kind = h.u32()?;
            let vsize = h.u32()? as usize;
            let begin = if wide {
                h.u64()? as usize
            } else {
                h.u32()? as usize
            };
            let is_record = ids.first().and_then(|&d| dims.get(d)) == Some(&0);
            if is_record {
                record_size += vsize;
                record_vars += 1;
            }
            if name == "coordinates" && is_record && ids.len() == 3 {
                let shape = [dims.get(ids[1]), dims.get(ids[2])];
                if let [Some(&atoms), Some(&3)] = shape {
                    if kind == FLOAT || kind == DOUBLE {
                        coordinates = Some((atoms, begin, kind == DOUBLE, vsize, scale));
                    }
                }
            }
        }
        let (atom_count, begin, double, vsize, scale) =
            coordinates.ok_or(NetCdfError::NoCoordinates)?;
        let frame_bytes = atom_count * 3 * if double { 8 } else { 4 };
        // A lone record variable is not padded to four bytes.
        if record_vars == 1 {
            record_size = frame_bytes;
        } else if vsize == 0 {
            return Err(NetCdfError::Malformed("coordinates have no size"));
        }
        // Frames whose coordinates are wholly in the file (the last record
        // may end with other variables after them).
        let fits = match mmap.len().checked_sub(begin + frame_bytes) {
            Some(rest) => rest / record_size.max(1) + 1,
            None => 0,
        };
        let frame_count = if numrecs == STREAMING {
            fits
        } else {
            (numrecs as usize).min(fits)
        };
        Ok(Self {
            mmap,
            atom_count,
            frame_count,
            begin,
            record_size,
            double,
            scale: scale.unwrap_or(1.0),
        })
    }

    pub fn atom_count(&self) -> usize {
        self.atom_count
    }

    pub fn frame_count(&self) -> usize {
        self.frame_count
    }

    /// Frame `index`'s positions, in Å.
    pub fn read_frame(&self, index: usize) -> Result<Vec<Vec3>, NetCdfError> {
        if index >= self.frame_count {
            return Err(NetCdfError::FrameOutOfRange {
                index,
                frame_count: self.frame_count,
            });
        }
        let bytes = &self.mmap[self.begin + index * self.record_size..];
        let value = |k: usize| -> f32 {
            if self.double {
                f64::from_be_bytes(bytes[k * 8..k * 8 + 8].try_into().expect("8 bytes")) as f32
            } else {
                f32::from_be_bytes(bytes[k * 4..k * 4 + 4].try_into().expect("4 bytes"))
            }
        };
        Ok((0..self.atom_count)
            .map(|i| Vec3::new(value(3 * i), value(3 * i + 1), value(3 * i + 2)) * self.scale)
            .collect())
    }
}
