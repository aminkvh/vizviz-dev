//! Per-atom value channels: numbers computed outside vizviz (a SASA from
//! FastSASA, a conservation score, an RMSF) attached to a loaded
//! structure by name so the viewer can color by them. This is the return
//! half of "pipe a structure through another package": positions and
//! names go out as NumPy arrays, a column of numbers comes back.
//!
//! A channel is one value per atom, or one per atom per frame for a
//! trajectory (`frames * atoms` values, frame-major). The color range is
//! fixed over all frames so a residue that opens up over time visibly
//! changes color instead of every frame being re-normalized.
//!
//! A session file persists a channel's numbers too, as a `.npy` sidecar
//! next to the session file (`session.rs`), so `color values NAME`
//! survives a save/load round trip without re-running the source tool.

use std::path::Path;
use std::sync::Arc;

use vv_core::Structure;

/// One named column of per-atom (or per-atom-per-frame) numbers.
#[derive(Clone, Debug, PartialEq)]
pub struct ValueChannel {
    /// `frames * atoms` values, frame-major.
    data: Arc<[f32]>,
    frames: usize,
    atoms: usize,
    /// Range over all finite values, for a stable color scale across frames.
    range: (f32, f32),
}

#[derive(Debug, thiserror::Error)]
pub enum ValuesError {
    #[error("cannot read {path}: {source}")]
    Io {
        path: std::path::PathBuf,
        source: std::io::Error,
    },
    #[error("{path}: {reason}")]
    Format {
        path: std::path::PathBuf,
        reason: String,
    },
    #[error(
        "expected one value per atom ({atoms}) or per atom per frame ({frames} x {atoms}), got shape {shape:?}"
    )]
    Shape {
        atoms: usize,
        frames: usize,
        shape: Vec<usize>,
    },
    #[error(
        "expected one value per atom ({atoms}) or per residue ({residues}), or one per frame of those ({frames} x ...), got shape {shape:?}"
    )]
    StructureShape {
        atoms: usize,
        residues: usize,
        frames: usize,
        shape: Vec<usize>,
    },
}

impl ValueChannel {
    /// Builds a channel from `data` of the given `shape` for a structure
    /// with `atoms` atoms and `frames` frames. Accepted shapes: `(atoms,)`,
    /// `(1, atoms)`, `(atoms, 1)`, and `(frames, atoms)`.
    pub fn new(
        data: Vec<f32>,
        shape: &[usize],
        atoms: usize,
        frames: usize,
    ) -> Result<Self, ValuesError> {
        let bad = || ValuesError::Shape {
            atoms,
            frames,
            shape: shape.to_vec(),
        };
        let channel_frames = match shape {
            [n] if *n == atoms => 1,
            [1, n] | [n, 1] if *n == atoms => 1,
            [f, n] if *n == atoms && *f == frames => frames,
            _ => return Err(bad()),
        };
        if data.len() != channel_frames * atoms {
            return Err(bad());
        }
        let range = range_of(&data);
        Ok(Self {
            data: data.into(),
            frames: channel_frames,
            atoms,
            range,
        })
    }

    /// The values for `frame` (clamped; a single-frame channel serves
    /// every frame).
    pub fn frame(&self, frame: usize) -> &[f32] {
        let f = frame.min(self.frames - 1);
        &self.data[f * self.atoms..(f + 1) * self.atoms]
    }

    pub fn frames(&self) -> usize {
        self.frames
    }

    pub fn atoms(&self) -> usize {
        self.atoms
    }

    /// `(min, max)` over all frames; `(0, 0)` when nothing is finite.
    pub fn range(&self) -> (f32, f32) {
        self.range
    }

    /// The whole column, frame-major.
    pub fn data(&self) -> &Arc<[f32]> {
        &self.data
    }

    /// `new` for a structure, also accepting one value per *residue*
    /// (`(residues,)` or `(frames, residues)`), expanded to the residue's
    /// atoms. Per-residue results (a residue SASA, a conservation score)
    /// are what most tools produce, and they read better on a tube.
    pub fn for_structure(
        data: Vec<f32>,
        shape: &[usize],
        structure: &Structure,
    ) -> Result<Self, ValuesError> {
        let atoms = structure.atom_count();
        let residues = structure.topology.residue_count();
        let frames = structure.frame_count();
        let bad = || ValuesError::StructureShape {
            atoms,
            residues,
            frames,
            shape: shape.to_vec(),
        };
        let per_residue = match shape {
            [n] | [1, n] | [n, 1] => *n == residues && residues != atoms,
            [f, n] => *n == residues && residues != atoms && (*f == frames || *f == 1),
            _ => false,
        };
        if !per_residue {
            return Self::new(data, shape, atoms, frames).map_err(|_| bad());
        }
        let rows = data.len() / residues.max(1);
        if data.len() != rows * residues {
            return Err(bad());
        }
        let index = &structure.topology.residue_index;
        let mut expanded = Vec::with_capacity(rows * atoms);
        for row in data.chunks_exact(residues) {
            expanded.extend(index.iter().map(|&r| row[r as usize]));
        }
        Self::new(expanded, &[rows, atoms], atoms, frames).map_err(|_| bad())
    }

    /// Reads `path` (`.npy`, or text with one number per line or per
    /// residue) for `structure`.
    pub fn read(path: &Path, structure: &Structure) -> Result<Self, ValuesError> {
        let (data, shape) = read_values(path)?;
        Self::for_structure(data, &shape, structure)
    }
}

/// Writes `channel` as a float32 `.npy` array: `(atoms,)` for a
/// single-frame channel, `(frames, atoms)` otherwise. Session-file sidecar
/// format; read back by `read_values`/`ValueChannel::new`, which accept
/// exactly these two shapes.
pub fn write_values(path: &Path, channel: &ValueChannel) -> std::io::Result<()> {
    let (frames, atoms) = (channel.frames(), channel.atoms());
    let shape = if frames > 1 {
        format!("({frames}, {atoms})")
    } else {
        format!("({atoms},)")
    };
    let mut data = Vec::with_capacity(channel.data().len() * 4);
    for v in channel.data().iter() {
        data.extend_from_slice(&v.to_le_bytes());
    }
    std::fs::write(path, build_npy("<f4", &shape, &data))
}

/// A minimal version-1.0 `.npy` file: C order, the dict header padded to a
/// 64-byte boundary the way `numpy.save` does.
fn build_npy(descr: &str, shape: &str, data: &[u8]) -> Vec<u8> {
    let dict = format!("{{'descr': '{descr}', 'fortran_order': False, 'shape': {shape}, }}");
    let mut header = dict.into_bytes();
    while (10 + header.len() + 1) % 64 != 0 {
        header.push(b' ');
    }
    header.push(b'\n');
    let mut out = b"\x93NUMPY\x01\x00".to_vec();
    out.extend((header.len() as u16).to_le_bytes());
    out.extend(header);
    out.extend_from_slice(data);
    out
}

fn range_of(values: &[f32]) -> (f32, f32) {
    let (min, max) = values
        .iter()
        .filter(|v| v.is_finite())
        .fold((f32::MAX, f32::MIN), |(lo, hi), &v| (lo.min(v), hi.max(v)));
    if min > max {
        (0.0, 0.0)
    } else {
        (min, max)
    }
}

/// Numbers and their shape from a `.npy` file (any numeric dtype, C
/// order) or a text file (one number per line, `#` comments, blank lines
/// skipped; shape `(n,)`).
pub fn read_values(path: &Path) -> Result<(Vec<f32>, Vec<usize>), ValuesError> {
    let bytes = std::fs::read(path).map_err(|source| ValuesError::Io {
        path: path.to_path_buf(),
        source,
    })?;
    let format = |reason: String| ValuesError::Format {
        path: path.to_path_buf(),
        reason,
    };
    if bytes.starts_with(b"\x93NUMPY") {
        return parse_npy(&bytes).map_err(format);
    }
    let text = std::str::from_utf8(&bytes).map_err(|_| format("not UTF-8 text".into()))?;
    let mut out = Vec::new();
    for (i, line) in text.lines().enumerate() {
        let line = line.split('#').next().unwrap_or("").trim();
        if line.is_empty() {
            continue;
        }
        // Take the last whitespace/comma separated field, so a two-column
        // `index value` file or a one-column CSV both work.
        let field = line
            .rsplit(|c: char| c.is_whitespace() || c == ',')
            .find(|f| !f.is_empty())
            .unwrap_or(line);
        let v: f32 = field
            .parse::<f64>()
            .map_err(|_| format(format!("line {}: `{field}` is not a number", i + 1)))?
            as f32;
        out.push(v);
    }
    let n = out.len();
    Ok((out, vec![n]))
}

/// Minimal `.npy` reader: versions 1.0-3.0, little-endian numeric dtypes,
/// C order. Enough for `np.save(path, values)`; anything else is a clear
/// error rather than a guess.
fn parse_npy(bytes: &[u8]) -> Result<(Vec<f32>, Vec<usize>), String> {
    if bytes.len() < 12 {
        return Err("truncated header".into());
    }
    let major = bytes[6];
    let (header_len, start) = match major {
        1 => {
            let n = u16::from_le_bytes([bytes[8], bytes[9]]) as usize;
            (n, 10)
        }
        2 | 3 => {
            let n = u32::from_le_bytes([bytes[8], bytes[9], bytes[10], bytes[11]]) as usize;
            (n, 12)
        }
        v => return Err(format!("unsupported .npy version {v}")),
    };
    let header = bytes
        .get(start..start + header_len)
        .ok_or("truncated header")?;
    let header = std::str::from_utf8(header).map_err(|_| "header is not UTF-8")?;
    let descr = dict_value(header, "descr").ok_or("header has no descr")?;
    let descr = descr.trim_matches(|c| c == '\'' || c == '"');
    if dict_value(header, "fortran_order").is_some_and(|v| v.trim() == "True") {
        return Err(
            "Fortran-ordered arrays are not supported; save with np.ascontiguousarray".into(),
        );
    }
    let shape_text = dict_value(header, "shape").ok_or("header has no shape")?;
    let shape: Vec<usize> = shape_text
        .trim()
        .trim_start_matches('(')
        .trim_end_matches(')')
        .split(',')
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .map(|s| {
            s.parse::<usize>()
                .map_err(|_| format!("bad shape `{shape_text}`"))
        })
        .collect::<Result<_, _>>()?;
    let count: usize = shape.iter().product();
    let data = &bytes[start + header_len..];

    let (endian, kind, size) = match descr.as_bytes() {
        [e @ (b'<' | b'|' | b'='), k, rest @ ..] => (
            *e,
            *k,
            std::str::from_utf8(rest)
                .ok()
                .and_then(|s| s.parse::<usize>().ok())
                .ok_or_else(|| format!("unsupported dtype {descr}"))?,
        ),
        [b'>', ..] => return Err("big-endian arrays are not supported".into()),
        _ => return Err(format!("unsupported dtype {descr}")),
    };
    let _ = endian;
    if data.len() < count * size {
        return Err(format!(
            "expected {} bytes of data for shape {shape:?}, found {}",
            count * size,
            data.len()
        ));
    }
    let chunks = data[..count * size].chunks_exact(size);
    let values: Vec<f32> = match (kind, size) {
        (b'f', 4) => chunks
            .map(|c| f32::from_le_bytes(c.try_into().unwrap()))
            .collect(),
        (b'f', 8) => chunks
            .map(|c| f64::from_le_bytes(c.try_into().unwrap()) as f32)
            .collect(),
        (b'f', 2) => chunks
            .map(|c| half_to_f32(u16::from_le_bytes([c[0], c[1]])))
            .collect(),
        (b'i', 1) => chunks.map(|c| c[0] as i8 as f32).collect(),
        (b'i', 2) => chunks
            .map(|c| i16::from_le_bytes([c[0], c[1]]) as f32)
            .collect(),
        (b'i', 4) => chunks
            .map(|c| i32::from_le_bytes(c.try_into().unwrap()) as f32)
            .collect(),
        (b'i', 8) => chunks
            .map(|c| i64::from_le_bytes(c.try_into().unwrap()) as f32)
            .collect(),
        (b'u', 1) | (b'b', 1) => chunks.map(|c| c[0] as f32).collect(),
        (b'u', 2) => chunks
            .map(|c| u16::from_le_bytes([c[0], c[1]]) as f32)
            .collect(),
        (b'u', 4) => chunks
            .map(|c| u32::from_le_bytes(c.try_into().unwrap()) as f32)
            .collect(),
        (b'u', 8) => chunks
            .map(|c| u64::from_le_bytes(c.try_into().unwrap()) as f32)
            .collect(),
        _ => {
            return Err(format!(
                "unsupported dtype {descr}; save as float32 or float64"
            ))
        }
    };
    Ok((values, shape))
}

/// The raw text after `'key':` in a Python dict literal, up to the next
/// top-level comma. Good enough for the three keys NumPy writes.
fn dict_value<'a>(header: &'a str, key: &str) -> Option<&'a str> {
    let pos = header
        .find(&format!("'{key}'"))
        .or_else(|| header.find(&format!("\"{key}\"")))?;
    let rest = &header[pos + key.len() + 2..];
    let rest = rest.trim_start().strip_prefix(':')?;
    let mut depth = 0i32;
    let mut in_quote: Option<char> = None;
    for (i, c) in rest.char_indices() {
        match (in_quote, c) {
            (Some(q), c) if c == q => in_quote = None,
            (Some(_), _) => {}
            (None, '\'' | '"') => in_quote = Some(c),
            (None, '(' | '[' | '{') => depth += 1,
            (None, ')' | ']' | '}') => {
                if depth == 0 {
                    return Some(rest[..i].trim());
                }
                depth -= 1;
            }
            (None, ',') if depth == 0 => return Some(rest[..i].trim()),
            _ => {}
        }
    }
    Some(rest.trim())
}

fn half_to_f32(h: u16) -> f32 {
    let sign = if h & 0x8000 != 0 { -1.0 } else { 1.0 };
    let exp = ((h >> 10) & 0x1f) as i32;
    let frac = (h & 0x3ff) as f32;
    match exp {
        // Zero and subnormals: frac * 2^-24.
        0 => sign * frac * 2f32.powi(-24),
        31 => {
            if frac == 0.0 {
                sign * f32::INFINITY
            } else {
                f32::NAN
            }
        }
        _ => sign * (1.0 + frac / 1024.0) * 2f32.powi(exp - 15),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn npy(descr: &str, shape: &str, data: &[u8]) -> Vec<u8> {
        build_npy(descr, shape, data)
    }

    #[test]
    fn reads_float64_and_int_npy() {
        let mut data = Vec::new();
        for v in [1.5f64, -2.0, 3.25] {
            data.extend(v.to_le_bytes());
        }
        let (values, shape) = parse_npy(&npy("<f8", "(3,)", &data)).unwrap();
        assert_eq!(values, [1.5, -2.0, 3.25]);
        assert_eq!(shape, [3]);

        let mut data = Vec::new();
        for v in [7i32, 8, 9, 10] {
            data.extend(v.to_le_bytes());
        }
        let (values, shape) = parse_npy(&npy("<i4", "(2, 2)", &data)).unwrap();
        assert_eq!(values, [7.0, 8.0, 9.0, 10.0]);
        assert_eq!(shape, [2, 2]);
    }

    #[test]
    fn rejects_fortran_order_and_odd_dtypes() {
        let bytes = npy("<f4", "(1,)", &1.0f32.to_le_bytes());
        let text = String::from_utf8(bytes.clone()).unwrap_or_default();
        let _ = text;
        let mut fortran = bytes.clone();
        let pos = fortran.windows(5).position(|w| w == b"False").unwrap();
        fortran.splice(pos..pos + 5, b"True ".iter().copied());
        assert!(parse_npy(&fortran).unwrap_err().contains("Fortran"));
        let err = parse_npy(&npy("|S4", "(1,)", b"abcd")).unwrap_err();
        assert!(err.contains("unsupported dtype"), "{err}");
        let err = parse_npy(&npy("<f8", "(4,)", &[0u8; 8])).unwrap_err();
        assert!(err.contains("expected 32 bytes"), "{err}");
    }

    #[test]
    fn text_files_take_the_last_field_per_line() {
        let dir = std::env::temp_dir().join(format!("vizviz_values_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("v.txt");
        std::fs::write(&path, "# atom sasa\n0 1.5\n1, 2.5\n\n2\t3.5\n").unwrap();
        let (values, shape) = read_values(&path).unwrap();
        assert_eq!(values, [1.5, 2.5, 3.5]);
        assert_eq!(shape, [3]);
        std::fs::write(&path, "1.0\nabc\n").unwrap();
        let err = read_values(&path).unwrap_err().to_string();
        assert!(err.contains("line 2"), "{err}");
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn channel_shapes_and_frames() {
        let c = ValueChannel::new(vec![1.0, 2.0, 3.0], &[3], 3, 1).unwrap();
        assert_eq!(c.frames(), 1);
        assert_eq!(c.frame(5), [1.0, 2.0, 3.0]);
        assert_eq!(c.range(), (1.0, 3.0));
        let c = ValueChannel::new(vec![1.0, 2.0, 3.0], &[1, 3], 3, 4).unwrap();
        assert_eq!(c.frames(), 1);
        let c = ValueChannel::new(vec![0.0, 0.0, 0.0, 9.0, 9.0, 9.0], &[2, 3], 3, 2).unwrap();
        assert_eq!(c.frames(), 2);
        assert_eq!(c.frame(1), [9.0, 9.0, 9.0]);
        assert_eq!(c.range(), (0.0, 9.0));
        assert!(matches!(
            ValueChannel::new(vec![1.0, 2.0], &[2], 3, 1),
            Err(ValuesError::Shape { .. })
        ));
        assert!(matches!(
            ValueChannel::new(vec![0.0; 6], &[2, 3], 3, 5),
            Err(ValuesError::Shape { .. })
        ));
        let c = ValueChannel::new(vec![f32::NAN, f32::INFINITY], &[2], 2, 1).unwrap();
        assert_eq!(c.range(), (0.0, 0.0));
    }

    #[test]
    fn write_values_round_trips_through_read_values() {
        let dir = std::env::temp_dir().join(format!("vizviz_write_values_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();

        let single = ValueChannel::new(vec![1.5, -2.0, 3.25], &[3], 3, 1).unwrap();
        let path = dir.join("single.npy");
        write_values(&path, &single).unwrap();
        let (data, shape) = read_values(&path).unwrap();
        assert_eq!(shape, [3]);
        assert_eq!(data, [1.5, -2.0, 3.25]);

        let per_frame =
            ValueChannel::new(vec![0.0, 0.0, 0.0, 9.0, 9.0, 9.0], &[2, 3], 3, 2).unwrap();
        let path = dir.join("frames.npy");
        write_values(&path, &per_frame).unwrap();
        let (data, shape) = read_values(&path).unwrap();
        assert_eq!(shape, [2, 3]);
        assert_eq!(data, [0.0, 0.0, 0.0, 9.0, 9.0, 9.0]);

        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn per_residue_columns_expand_to_atoms() {
        let path = concat!(env!("CARGO_MANIFEST_DIR"), "/../../fixtures/small/1CRN.cif");
        let s = vv_io::load(std::path::Path::new(path)).unwrap();
        let (atoms, residues) = (s.atom_count(), s.topology.residue_count());
        assert_ne!(atoms, residues);
        let per_res: Vec<f32> = (0..residues).map(|r| r as f32).collect();
        let c = ValueChannel::for_structure(per_res.clone(), &[residues], &s).unwrap();
        assert_eq!(c.atoms(), atoms);
        assert_eq!(c.frame(0)[atoms - 1], (residues - 1) as f32);
        assert_eq!(c.frame(0)[0], 0.0);
        let two = [per_res.clone(), per_res].concat();
        let err = ValueChannel::for_structure(two, &[2, residues], &s).unwrap_err();
        assert!(err.to_string().contains("per residue"), "{err}");
        let per_atom: Vec<f32> = vec![1.0; atoms];
        assert_eq!(
            ValueChannel::for_structure(per_atom, &[atoms], &s)
                .unwrap()
                .frames(),
            1
        );
    }

    #[test]
    fn dict_values_stop_at_top_level_commas() {
        let h = "{'descr': '<f8', 'fortran_order': False, 'shape': (2, 3), }";
        assert_eq!(dict_value(h, "descr"), Some("'<f8'"));
        assert_eq!(dict_value(h, "fortran_order"), Some("False"));
        assert_eq!(dict_value(h, "shape"), Some("(2, 3)"));
        assert_eq!(dict_value("{'shape': (7,), }", "shape"), Some("(7,)"));
    }

    #[test]
    fn half_floats_convert() {
        assert_eq!(half_to_f32(0x3c00), 1.0);
        assert_eq!(half_to_f32(0xc000), -2.0);
        assert_eq!(half_to_f32(0x0000), 0.0);
        assert!(half_to_f32(0x7c00).is_infinite());
        assert!((half_to_f32(0x0001) - 5.960_464_5e-8).abs() < 1e-12);
    }
}
