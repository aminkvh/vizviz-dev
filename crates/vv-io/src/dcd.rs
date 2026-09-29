//! Fixed-record-size trajectory reader for the CHARMM/NAMD binary DCD
//! format ("CORD" Fortran-unformatted records with 32-bit markers), as
//! written by CHARMM, NAMD, and MDAnalysis. Every frame occupies the
//! same number of bytes, so any frame can be read directly at
//! `header_len + index * frame_len` without decoding the frames before it
//! -- the property that makes multi-hundred-GB trajectories viable on a
//! laptop. The whole file
//! is memory-mapped (same convention as `load()` for PDB/mmCIF below), so
//! the OS pages in only the bytes a given `read_frame` call touches.
//!
//! Verified against a real DCD file written by MDAnalysis's own
//! `DCDWriter`, not just against our own writer: `vv-io/tests/dcd.rs`
//! checks known per-atom, per-frame coordinates byte-for-byte against
//! `fixtures/small/sample.dcd`.
//!
//! Scope, explicit: reads the common modern layout (`ICNTRL(20) != 0`,
//! i.e. "CHARMM-format" per every DCD reader's own convention) with
//! single-precision coordinates and 32-bit record markers -- what every
//! DCD-writing tool in practice produces today. Deliberately not read:
//! the pre-CHARMM22 layout (`ICNTRL(20) == 0`, a different and essentially
//! extinct header shape), the fixed-atom optimization (`ICNTRL(9)` names a
//! free-atom count smaller than `NATOM`, so only the first frame carries
//! every atom), 64-bit record markers, and the unit-cell record's actual
//! box vectors (its presence is read so frames are sized correctly; the
//! six numbers themselves are skipped, not exposed).

use std::io;
use std::path::Path;

use memmap2::Mmap;
use vv_core::glam::Vec3;

#[derive(thiserror::Error, Debug)]
pub enum DcdError {
    #[error("cannot read {path}: {source}")]
    Io { path: String, source: io::Error },
    #[error("not a DCD file (missing CORD magic)")]
    NotDcd,
    #[error("unsupported DCD variant: pre-CHARMM22 header (ICNTRL(20) == 0)")]
    UnsupportedVariant,
    #[error("unsupported DCD variant: fixed-atom trajectories (a free-atom count in ICNTRL(9) smaller than NATOM) are not read")]
    FixedAtoms,
    #[error("truncated or malformed DCD file: {0}")]
    Malformed(&'static str),
    #[error("frame {index} is out of range (this trajectory has {frame_count})")]
    FrameOutOfRange { index: usize, frame_count: usize },
}

/// An open DCD file: header parsed, ready for `O(1)` random-access reads
/// of any frame.
pub struct DcdReader {
    mmap: Mmap,
    natom: usize,
    frame_count: usize,
    has_unit_cell: bool,
    frames_start: usize,
    frame_len: usize,
}

impl DcdReader {
    pub fn open(path: impl AsRef<Path>) -> Result<Self, DcdError> {
        let path_ref = path.as_ref();
        let io_err = |source| DcdError::Io {
            path: path_ref.display().to_string(),
            source,
        };
        let file = std::fs::File::open(path_ref).map_err(io_err)?;
        let mmap = unsafe { Mmap::map(&file) }.map_err(io_err)?;
        Self::from_mmap(mmap)
    }

    fn from_mmap(mmap: Mmap) -> Result<Self, DcdError> {
        let bytes: &[u8] = &mmap;
        let mut pos = 0usize;

        let header = read_record(bytes, &mut pos)?;
        if header.len() != 84 || &header[0..4] != b"CORD" {
            return Err(DcdError::NotDcd);
        }
        let icntrl: Vec<i32> = header[4..84]
            .chunks_exact(4)
            .map(|c| i32::from_le_bytes(c.try_into().unwrap()))
            .collect();
        let nset = icntrl[0].max(0) as usize;
        let has_unit_cell = icntrl[10] != 0; // ICNTRL(11)
        let charmm_version = icntrl[19]; // ICNTRL(20)
        if charmm_version == 0 {
            return Err(DcdError::UnsupportedVariant);
        }

        let title = read_record(bytes, &mut pos)?;
        let ntitle = i32::from_le_bytes(
            title
                .get(0..4)
                .and_then(|b| b.try_into().ok())
                .ok_or(DcdError::Malformed("title"))?,
        )
        .max(0) as usize;
        if title.len() != 4 + ntitle * 80 {
            return Err(DcdError::Malformed("title length"));
        }

        let natom_rec = read_record(bytes, &mut pos)?;
        let natom = i32::from_le_bytes(
            natom_rec
                .try_into()
                .map_err(|_| DcdError::Malformed("natom"))?,
        )
        .max(0) as usize;

        let nfreat = icntrl[8]; // ICNTRL(9): 0 means "no fixed atoms", not "zero free"
        if nfreat != 0 && (nfreat as usize) < natom {
            return Err(DcdError::FixedAtoms);
        }

        let frames_start = pos;
        let cell_record_len = if has_unit_cell { 4 + 48 + 4 } else { 0 };
        let coord_record_len = 4 + natom * 4 + 4;
        let frame_len = cell_record_len + 3 * coord_record_len;
        let available = bytes.len().saturating_sub(frames_start);
        // Trust the smaller of the declared frame count and what the file
        // actually holds, rather than reading past a truncated file.
        let frame_count = nset.min(available / frame_len);

        Ok(Self {
            mmap,
            natom,
            frame_count,
            has_unit_cell,
            frames_start,
            frame_len,
        })
    }

    pub fn atom_count(&self) -> usize {
        self.natom
    }

    pub fn frame_count(&self) -> usize {
        self.frame_count
    }

    /// Reads frame `index` directly by computed byte offset -- no earlier
    /// frame is read or decoded first.
    pub fn read_frame(&self, index: usize) -> Result<Vec<Vec3>, DcdError> {
        if index >= self.frame_count {
            return Err(DcdError::FrameOutOfRange {
                index,
                frame_count: self.frame_count,
            });
        }
        let bytes: &[u8] = &self.mmap;
        let mut pos = self.frames_start + index * self.frame_len;
        if self.has_unit_cell {
            let cell = read_record(bytes, &mut pos)?;
            if cell.len() != 48 {
                return Err(DcdError::Malformed("unit cell record"));
            }
        }
        let x = read_axis(bytes, &mut pos, self.natom)?;
        let y = read_axis(bytes, &mut pos, self.natom)?;
        let z = read_axis(bytes, &mut pos, self.natom)?;
        Ok((0..self.natom)
            .map(|i| Vec3::new(x[i], y[i], z[i]))
            .collect())
    }
}

/// Reads one Fortran-unformatted record at `*pos`, validates its leading
/// and trailing length markers agree, and advances `*pos` past it.
fn read_record<'a>(bytes: &'a [u8], pos: &mut usize) -> Result<&'a [u8], DcdError> {
    let start = *pos;
    let head = bytes
        .get(start..start + 4)
        .ok_or(DcdError::Malformed("truncated record marker"))?;
    let len = i32::from_le_bytes(head.try_into().unwrap());
    if len < 0 {
        return Err(DcdError::Malformed("negative record length"));
    }
    let len = len as usize;
    let body_start = start + 4;
    let body_end = body_start + len;
    let tail = bytes
        .get(body_end..body_end + 4)
        .ok_or(DcdError::Malformed("truncated record"))?;
    let len2 = i32::from_le_bytes(tail.try_into().unwrap());
    if len2 as usize != len {
        return Err(DcdError::Malformed("mismatched record markers"));
    }
    *pos = body_end + 4;
    Ok(&bytes[body_start..body_end])
}

fn read_axis(bytes: &[u8], pos: &mut usize, natom: usize) -> Result<Vec<f32>, DcdError> {
    let body = read_record(bytes, pos)?;
    if body.len() != natom * 4 {
        return Err(DcdError::Malformed("coordinate record length"));
    }
    Ok(body
        .chunks_exact(4)
        .map(|c| f32::from_le_bytes(c.try_into().unwrap()))
        .collect())
}
