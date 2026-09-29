//! GROMACS TRR trajectories: uncompressed frames, each a header saying
//! which blocks follow (box, virial, pressure, positions, velocities,
//! forces) and whether reals are single or double precision. Opening
//! scans the headers into an index of frames that carry positions; any of
//! them is then read directly. Positions come back in Å.
//!
//! Layout as xdrfile's `do_trnheader` / `do_htrn` (GROMACS, BSD licence);
//! `vv-io/tests/traj.rs` checks it against MDAnalysis on
//! `fixtures/traj/1crn.trr`.

use std::io;
use std::path::Path;

use memmap2::Mmap;
use vv_core::glam::Vec3;

const MAGIC: i32 = 1993;

#[derive(thiserror::Error, Debug)]
pub enum TrrError {
    #[error("cannot read {path}: {source}")]
    Io { path: String, source: io::Error },
    #[error("not a TRR file (bad magic number)")]
    NotTrr,
    #[error("truncated or malformed TRR file: {0}")]
    Malformed(&'static str),
    #[error("frame {index} is out of range (this trajectory has {frame_count})")]
    FrameOutOfRange { index: usize, frame_count: usize },
}

/// Where one frame's positions are.
#[derive(Clone, Copy)]
struct FrameAt {
    positions: usize,
    double: bool,
}

pub struct TrrReader {
    mmap: Mmap,
    atom_count: usize,
    frames: Vec<FrameAt>,
}

fn be_i32(bytes: &[u8], at: usize) -> Result<i32, TrrError> {
    let b = bytes
        .get(at..at + 4)
        .ok_or(TrrError::Malformed("header runs past the end of the file"))?;
    Ok(i32::from_be_bytes([b[0], b[1], b[2], b[3]]))
}

impl TrrReader {
    pub fn open(path: impl AsRef<Path>) -> Result<Self, TrrError> {
        let path = path.as_ref();
        let io_err = |source| TrrError::Io {
            path: path.display().to_string(),
            source,
        };
        let file = std::fs::File::open(path).map_err(io_err)?;
        // SAFETY: read-only mapping, as `dcd`.
        let mmap = unsafe { Mmap::map(&file) }.map_err(io_err)?;
        let mut frames = Vec::new();
        let mut atom_count = None;
        let mut at = 0;
        while at < mmap.len() {
            if be_i32(&mmap, at)? != MAGIC {
                return Err(if at == 0 {
                    TrrError::NotTrr
                } else {
                    TrrError::Malformed("bad magic number mid-file")
                });
            }
            // Magic, the version string's length plus one, then the string
            // itself (XDR: length, bytes padded to four).
            let len = be_i32(&mmap, at + 8)? as usize;
            let mut p = at + 12 + len.div_ceil(4) * 4;
            let mut ints = [0usize; 13];
            for v in ints.iter_mut() {
                *v = be_i32(&mmap, p)? as usize;
                p += 4;
            }
            // ir, e, box, vir, pres, top, sym, x, v, f, natoms, step, nre.
            let [_, _, box_size, vir, pres, _, _, x, v, f, natoms, _, _] = ints;
            let real = if box_size > 0 {
                box_size / 9
            } else if natoms > 0 && (x | v | f) > 0 {
                x.max(v).max(f) / (natoms * 3)
            } else {
                return Err(TrrError::Malformed("frame has no box or coordinates"));
            };
            if real != 4 && real != 8 {
                return Err(TrrError::Malformed("reals are neither 4 nor 8 bytes"));
            }
            // Time and lambda.
            p += 2 * real;
            let positions = p + box_size + vir + pres;
            if x > 0 {
                if *atom_count.get_or_insert(natoms) != natoms {
                    return Err(TrrError::Malformed("atom count changes between frames"));
                }
                if x != natoms * 3 * real {
                    return Err(TrrError::Malformed("position block has the wrong size"));
                }
                frames.push(FrameAt {
                    positions,
                    double: real == 8,
                });
            }
            at = positions + x + v + f;
            if at > mmap.len() {
                return Err(TrrError::Malformed("last frame is truncated"));
            }
        }
        let atom_count = atom_count.ok_or(TrrError::Malformed("no frame has positions"))?;
        Ok(Self {
            mmap,
            atom_count,
            frames,
        })
    }

    pub fn atom_count(&self) -> usize {
        self.atom_count
    }

    /// Frames with positions.
    pub fn frame_count(&self) -> usize {
        self.frames.len()
    }

    /// Frame `index`'s positions, in Å.
    pub fn read_frame(&self, index: usize) -> Result<Vec<Vec3>, TrrError> {
        let f = *self.frames.get(index).ok_or(TrrError::FrameOutOfRange {
            index,
            frame_count: self.frames.len(),
        })?;
        let bytes = &self.mmap[f.positions..];
        let value = |k: usize| -> f32 {
            if f.double {
                let b = &bytes[k * 8..k * 8 + 8];
                f64::from_be_bytes(b.try_into().expect("8 bytes")) as f32
            } else {
                let b = &bytes[k * 4..k * 4 + 4];
                f32::from_be_bytes(b.try_into().expect("4 bytes"))
            }
        };
        Ok((0..self.atom_count)
            .map(|i| Vec3::new(value(3 * i), value(3 * i + 1), value(3 * i + 2)) * 10.0)
            .collect())
    }
}
