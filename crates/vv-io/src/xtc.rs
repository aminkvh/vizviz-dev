//! GROMACS XTC trajectories: per frame a header, the box, and the
//! positions in nm, lossily compressed ("xdr3dfcoord": coordinates as
//! integers at a stated precision, packed in bits, with runs of small
//! differences between neighbouring atoms). Frames vary in length, so
//! opening scans the headers once into a frame index (each header says
//! how many compressed bytes follow) and any frame is then decoded
//! directly. Positions come back in Å.
//!
//! The decoder follows the xdrfile library's `xdrfile_decompress_coord_
//! float` (GROMACS, BSD licence) step for step; `vv-io/tests/traj.rs`
//! checks it against MDAnalysis on `fixtures/traj/1crn.xtc`.

use std::io;
use std::path::Path;

use memmap2::Mmap;
use vv_core::glam::Vec3;

const MAGIC: i32 = 1995;
/// Frames of at most this many atoms are stored uncompressed.
const UNCOMPRESSED_MAX: usize = 9;
const FIRST_IDX: usize = 9;
const MAGIC_INTS: [u32; 73] = [
    0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 10, 12, 16, 20, 25, 32, 40, 50, 64, 80, 101, 128, 161, 203, 256,
    322, 406, 512, 645, 812, 1024, 1290, 1625, 2048, 2580, 3250, 4096, 5060, 6501, 8192, 10321,
    13003, 16384, 20642, 26007, 32768, 41285, 52015, 65536, 82570, 104031, 131072, 165140, 208063,
    262144, 330280, 416127, 524287, 660561, 832255, 1048576, 1321122, 1664510, 2097152, 2642245,
    3329021, 4194304, 5284491, 6658042, 8388607, 10568983, 13316085, 16777216,
];

#[derive(thiserror::Error, Debug)]
pub enum XtcError {
    #[error("cannot read {path}: {source}")]
    Io { path: String, source: io::Error },
    #[error("not an XTC file (bad magic number)")]
    NotXtc,
    #[error("truncated or malformed XTC file: {0}")]
    Malformed(&'static str),
    #[error("frame {index} is out of range (this trajectory has {frame_count})")]
    FrameOutOfRange { index: usize, frame_count: usize },
}

/// An open XTC file with its frame index.
pub struct XtcReader {
    mmap: Mmap,
    atom_count: usize,
    /// Byte offset of each frame.
    offsets: Vec<usize>,
}

/// Big-endian XDR reads over a byte slice.
struct Xdr<'a> {
    bytes: &'a [u8],
    at: usize,
}

impl Xdr<'_> {
    fn u32(&mut self) -> Result<u32, XtcError> {
        let b = self
            .bytes
            .get(self.at..self.at + 4)
            .ok_or(XtcError::Malformed("frame runs past the end of the file"))?;
        self.at += 4;
        Ok(u32::from_be_bytes([b[0], b[1], b[2], b[3]]))
    }

    fn i32(&mut self) -> Result<i32, XtcError> {
        Ok(self.u32()? as i32)
    }

    fn f32(&mut self) -> Result<f32, XtcError> {
        Ok(f32::from_bits(self.u32()?))
    }

    fn opaque(&mut self, len: usize) -> Result<&[u8], XtcError> {
        let data = self
            .bytes
            .get(self.at..self.at + len)
            .ok_or(XtcError::Malformed("compressed data runs past the end"))?;
        self.at += len.div_ceil(4) * 4;
        Ok(data)
    }
}

/// Frame header: magic, atoms, step, time, 3x3 box, atoms again.
const HEADER: usize = 4 * (4 + 9 + 1);

impl XtcReader {
    pub fn open(path: impl AsRef<Path>) -> Result<Self, XtcError> {
        let path = path.as_ref();
        let io_err = |source| XtcError::Io {
            path: path.display().to_string(),
            source,
        };
        let file = std::fs::File::open(path).map_err(io_err)?;
        // SAFETY: read-only mapping; a file truncated underneath us is a
        // documented hazard of every mmap-based reader (as `dcd`).
        let mmap = unsafe { Mmap::map(&file) }.map_err(io_err)?;
        let mut offsets = Vec::new();
        let mut atom_count = None;
        let mut at = 0;
        while at < mmap.len() {
            let mut x = Xdr { bytes: &mmap, at };
            if x.i32()? != MAGIC {
                return Err(if offsets.is_empty() {
                    XtcError::NotXtc
                } else {
                    XtcError::Malformed("bad magic number mid-file")
                });
            }
            let n = x.i32()? as usize;
            if *atom_count.get_or_insert(n) != n {
                return Err(XtcError::Malformed("atom count changes between frames"));
            }
            x.at = at + HEADER;
            if n <= UNCOMPRESSED_MAX {
                x.at += n * 12;
            } else {
                // precision, minint[3], maxint[3], smallidx, then the bytes.
                x.at += 4 * 8;
                let len = x.u32()? as usize;
                x.at += len.div_ceil(4) * 4;
            }
            if x.at > mmap.len() {
                return Err(XtcError::Malformed("last frame is truncated"));
            }
            offsets.push(at);
            at = x.at;
        }
        let atom_count = atom_count.ok_or(XtcError::Malformed("no frames"))?;
        Ok(Self {
            mmap,
            atom_count,
            offsets,
        })
    }

    pub fn atom_count(&self) -> usize {
        self.atom_count
    }

    pub fn frame_count(&self) -> usize {
        self.offsets.len()
    }

    /// Frame `index`'s positions, in Å.
    pub fn read_frame(&self, index: usize) -> Result<Vec<Vec3>, XtcError> {
        let &at = self.offsets.get(index).ok_or(XtcError::FrameOutOfRange {
            index,
            frame_count: self.offsets.len(),
        })?;
        let mut x = Xdr {
            bytes: &self.mmap,
            at: at + HEADER,
        };
        let n = self.atom_count;
        if n <= UNCOMPRESSED_MAX {
            return (0..n)
                .map(|_| Ok(Vec3::new(x.f32()?, x.f32()?, x.f32()?) * 10.0))
                .collect();
        }
        let precision = x.f32()?;
        let minint = [x.i32()?, x.i32()?, x.i32()?];
        let maxint = [x.i32()?, x.i32()?, x.i32()?];
        let smallidx = x.u32()? as usize;
        let len = x.u32()? as usize;
        let data = x.opaque(len)?;
        let ints = decompress(data, n, minint, maxint, smallidx)?;
        // nm -> Å.
        let scale = 10.0 / precision;
        Ok(ints
            .iter()
            .map(|c| Vec3::new(c[0] as f32, c[1] as f32, c[2] as f32) * scale)
            .collect())
    }
}

/// Reads bits most-significant first, as xdrfile's `receivebits`.
struct Bits<'a> {
    data: &'a [u8],
    at: usize,
    last_bits: u32,
    last_byte: u32,
}

impl Bits<'_> {
    fn next_byte(&mut self) -> Result<u32, XtcError> {
        let b = *self
            .data
            .get(self.at)
            .ok_or(XtcError::Malformed("compressed data ends early"))?;
        self.at += 1;
        Ok(b as u32)
    }

    fn read(&mut self, mut bits: u32) -> Result<u32, XtcError> {
        let mask = if bits >= 32 {
            u32::MAX
        } else {
            (1u32 << bits) - 1
        };
        let mut num = 0u32;
        while bits >= 8 {
            self.last_byte = (self.last_byte << 8) | self.next_byte()?;
            num |= (self.last_byte >> self.last_bits) << (bits - 8);
            bits -= 8;
        }
        if bits > 0 {
            if self.last_bits < bits {
                self.last_bits += 8;
                self.last_byte = (self.last_byte << 8) | self.next_byte()?;
            }
            self.last_bits -= bits;
            num |= (self.last_byte >> self.last_bits) & ((1 << bits) - 1);
        }
        Ok(num & mask)
    }

    /// Three integers packed together in `bits` bits, each below its
    /// `sizes` entry (xdrfile's `receiveints`).
    fn read_ints(&mut self, mut bits: u32, sizes: [u32; 3]) -> Result<[i32; 3], XtcError> {
        let mut bytes = [0u32; 32];
        let mut count = 0;
        while bits > 8 {
            bytes[count] = self.read(8)?;
            count += 1;
            bits -= 8;
        }
        if bits > 0 {
            bytes[count] = self.read(bits)?;
            count += 1;
        }
        let mut out = [0i32; 3];
        for i in (1..3).rev() {
            let mut num: u64 = 0;
            for j in (0..count).rev() {
                num = (num << 8) | bytes[j] as u64;
                let p = num / sizes[i] as u64;
                bytes[j] = p as u32;
                num -= p * sizes[i] as u64;
            }
            out[i] = num as i32;
        }
        out[0] = (bytes[0] | bytes[1] << 8 | bytes[2] << 16 | bytes[3] << 24) as i32;
        Ok(out)
    }
}

/// Bits needed for values below `size` (xdrfile's `sizeofint`).
fn bits_for(size: u32) -> u32 {
    let mut num: u64 = 1;
    let mut bits = 0;
    while size as u64 >= num && bits < 32 {
        bits += 1;
        num <<= 1;
    }
    bits
}

/// Bits needed for three values packed as one number below the product of
/// `sizes` (xdrfile's `sizeofints`).
fn bits_for_all(sizes: [u32; 3]) -> u32 {
    let mut bytes = [0u32; 32];
    bytes[0] = 1;
    let mut count = 1;
    for &size in &sizes {
        let mut tmp: u64 = 0;
        let mut k = 0;
        while k < count {
            tmp += bytes[k] as u64 * size as u64;
            bytes[k] = (tmp & 0xff) as u32;
            tmp >>= 8;
            k += 1;
        }
        while tmp != 0 {
            bytes[k] = (tmp & 0xff) as u32;
            k += 1;
            tmp >>= 8;
        }
        count = k;
    }
    let mut bits = 0;
    let mut num: u32 = 1;
    count -= 1;
    while bytes[count] >= num {
        bits += 1;
        num *= 2;
    }
    bits + count as u32 * 8
}

/// Integer coordinates of `n` atoms from xdr3dfcoord data.
fn decompress(
    data: &[u8],
    n: usize,
    minint: [i32; 3],
    maxint: [i32; 3],
    mut smallidx: usize,
) -> Result<Vec<[i32; 3]>, XtcError> {
    let sizeint = [0, 1, 2].map(|k| (maxint[k] - minint[k] + 1) as u32);
    // Very large ranges are sent one integer at a time.
    let big = (sizeint[0] | sizeint[1] | sizeint[2]) > 0xff_ffff;
    let bitsizeint = sizeint.map(bits_for);
    let bitsize = if big { 0 } else { bits_for_all(sizeint) };
    if !(FIRST_IDX..MAGIC_INTS.len()).contains(&smallidx) {
        return Err(XtcError::Malformed("small-difference index out of range"));
    }
    let mut smaller = MAGIC_INTS[smallidx.saturating_sub(1).max(FIRST_IDX)] as i32 / 2;
    let mut smallnum = MAGIC_INTS[smallidx] as i32 / 2;
    let mut sizesmall = [MAGIC_INTS[smallidx]; 3];

    let mut bits = Bits {
        data,
        at: 0,
        last_bits: 0,
        last_byte: 0,
    };
    let mut out = Vec::with_capacity(n);
    let mut run = 0i32;
    while out.len() < n {
        let mut this = if big {
            [
                bits.read(bitsizeint[0])? as i32,
                bits.read(bitsizeint[1])? as i32,
                bits.read(bitsizeint[2])? as i32,
            ]
        } else {
            bits.read_ints(bitsize, sizeint)?
        };
        for k in 0..3 {
            this[k] += minint[k];
        }
        let mut prev = this;
        let mut is_smaller = 0i32;
        if bits.read(1)? == 1 {
            run = bits.read(5)? as i32;
            is_smaller = run % 3;
            run -= is_smaller;
            is_smaller -= 1;
        }
        if run > 0 {
            let mut k = 0;
            while k < run {
                let d = bits.read_ints(smallidx as u32, sizesmall)?;
                let mut next = [0; 3];
                for c in 0..3 {
                    next[c] = d[c] + prev[c] - smallnum;
                }
                if k == 0 {
                    // The writer swaps the first two atoms of a run (it
                    // compresses water better): the run's first atom comes
                    // out before the big one.
                    std::mem::swap(&mut next, &mut prev);
                    out.push(prev);
                } else {
                    prev = next;
                }
                out.push(next);
                k += 3;
            }
        } else {
            out.push(this);
        }
        if out.len() > n {
            return Err(XtcError::Malformed("more atoms than the header says"));
        }
        smallidx = (smallidx as i32 + is_smaller) as usize;
        if !(FIRST_IDX..MAGIC_INTS.len()).contains(&smallidx) {
            return Err(XtcError::Malformed("small-difference index out of range"));
        }
        if is_smaller < 0 {
            smallnum = smaller;
            smaller = if smallidx > FIRST_IDX {
                MAGIC_INTS[smallidx - 1] as i32 / 2
            } else {
                0
            };
        } else if is_smaller > 0 {
            smaller = smallnum;
            smallnum = MAGIC_INTS[smallidx] as i32 / 2;
        }
        sizesmall = [MAGIC_INTS[smallidx]; 3];
    }
    Ok(out)
}
