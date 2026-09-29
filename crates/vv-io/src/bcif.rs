//! BinaryCIF reader, written from the BinaryCIF specification.
//!
//! A file is MessagePack: data blocks of categories of columns, each
//! column a byte array plus the list of encodings that produced it.
//! Decoding undoes the encodings last to first. The decoded columns are
//! re-emitted as CIF text of the first data block with an `_atom_site`
//! table and read by `mmcif::parse`, so both syntaxes share one reader for
//! atoms, secondary structure, connectivity, entities and annotations.

use std::io::Write;

use crate::msgpack::{self, Value};
use crate::ParseError;

/// A decoded column. Strings are `None` where the file has a null index.
enum Cells<'a> {
    Ints(Vec<i64>),
    /// `single`: the values are 32-bit floats and print at that precision.
    Floats {
        values: Vec<f64>,
        single: bool,
    },
    Strings(Vec<Option<&'a str>>),
}

impl Cells<'_> {
    fn len(&self) -> usize {
        match self {
            Cells::Ints(v) => v.len(),
            Cells::Floats { values, .. } => values.len(),
            Cells::Strings(v) => v.len(),
        }
    }

    fn ints(self) -> Result<Vec<i64>, String> {
        match self {
            Cells::Ints(v) => Ok(v),
            _ => Err("an integer encoding was applied to non-integer data".into()),
        }
    }
}

fn bad(message: impl Into<String>) -> ParseError {
    ParseError::BinaryCif(message.into())
}

/// The fixed-width little-endian number types of `ByteArray`.
fn from_bytes<'a>(bytes: &[u8], type_code: i64) -> Result<Cells<'a>, String> {
    fn chunks<const N: usize, T>(bytes: &[u8], f: impl Fn([u8; N]) -> T) -> Result<Vec<T>, String> {
        if bytes.len() % N != 0 {
            return Err(format!("byte array length is not a multiple of {N}"));
        }
        Ok(bytes
            .chunks_exact(N)
            .map(|c| f(c.try_into().expect("exact chunk")))
            .collect())
    }
    Ok(match type_code {
        1 => Cells::Ints(chunks::<1, _>(bytes, |b| i8::from_le_bytes(b) as i64)?),
        2 => Cells::Ints(chunks::<2, _>(bytes, |b| i16::from_le_bytes(b) as i64)?),
        3 => Cells::Ints(chunks::<4, _>(bytes, |b| i32::from_le_bytes(b) as i64)?),
        4 => Cells::Ints(bytes.iter().map(|&b| b as i64).collect()),
        5 => Cells::Ints(chunks::<2, _>(bytes, |b| u16::from_le_bytes(b) as i64)?),
        6 => Cells::Ints(chunks::<4, _>(bytes, |b| u32::from_le_bytes(b) as i64)?),
        32 => Cells::Floats {
            values: chunks::<4, _>(bytes, |b| f32::from_le_bytes(b) as f64)?,
            single: true,
        },
        33 => Cells::Floats {
            values: chunks::<8, _>(bytes, f64::from_le_bytes)?,
            single: false,
        },
        other => return Err(format!("unknown ByteArray type {other}")),
    })
}

fn int_param(enc: &Value<'_>, key: &str) -> Result<i64, String> {
    enc.get(key)
        .and_then(Value::as_i64)
        .ok_or_else(|| format!("encoding lacks `{key}`"))
}

fn float_param(enc: &Value<'_>, key: &str) -> Result<f64, String> {
    enc.get(key)
        .and_then(Value::as_f64)
        .ok_or_else(|| format!("encoding lacks `{key}`"))
}

/// Sums each run of saturated values with the value that ends it. An
/// unsigned stream never holds the (negative) lower bound, so only its
/// upper bound saturates.
fn unpack_integers(packed: &[i64], unsigned: bool, byte_count: i64) -> Result<Vec<i64>, String> {
    let (upper, lower) = match (byte_count, unsigned) {
        (1, true) => (0xff, -0x100),
        (1, false) => (0x7f, -0x80),
        (2, true) => (0xffff, -0x10000),
        (2, false) => (0x7fff, -0x8000),
        _ => return Err(format!("IntegerPacking byteCount {byte_count}")),
    };
    let mut out = Vec::new();
    let mut i = 0;
    while i < packed.len() {
        let mut value = 0i64;
        while i < packed.len() && (packed[i] == upper || packed[i] == lower) {
            value += packed[i];
            i += 1;
        }
        value += *packed.get(i).ok_or("IntegerPacking ends inside a run")?;
        i += 1;
        out.push(value);
    }
    Ok(out)
}

fn run_length(pairs: &[i64], size: usize) -> Result<Vec<i64>, String> {
    let mut out = Vec::with_capacity(size);
    for pair in pairs.chunks(2) {
        let [value, count] = pair else {
            return Err("RunLength data has an odd length".into());
        };
        if *count < 0 || out.len() + *count as usize > size {
            return Err("RunLength expands past its declared size".into());
        }
        out.extend(std::iter::repeat_n(*value, *count as usize));
    }
    Ok(out)
}

fn delta(diffs: Vec<i64>, origin: i64) -> Vec<i64> {
    let mut acc = origin;
    diffs
        .into_iter()
        .map(|d| {
            acc += d;
            acc
        })
        .collect()
}

/// Applies the inverse of one non-`ByteArray`, non-`StringArray` encoding.
fn undo<'a>(cells: Cells<'a>, enc: &Value<'_>, kind: &str) -> Result<Cells<'a>, String> {
    let ints = cells.ints()?;
    Ok(match kind {
        "IntegerPacking" => Cells::Ints(unpack_integers(
            &ints,
            matches!(enc.get("isUnsigned"), Some(Value::Bool(true))),
            int_param(enc, "byteCount")?,
        )?),
        "Delta" => Cells::Ints(delta(ints, int_param(enc, "origin")?)),
        "RunLength" => Cells::Ints(run_length(&ints, int_param(enc, "srcSize")? as usize)?),
        "FixedPoint" => {
            let factor = float_param(enc, "factor")?;
            Cells::Floats {
                values: ints.iter().map(|&i| i as f64 / factor).collect(),
                single: int_param(enc, "srcType")? == 32,
            }
        }
        "IntervalQuantization" => {
            let (min, max) = (float_param(enc, "min")?, float_param(enc, "max")?);
            let steps = float_param(enc, "numSteps")?;
            let scale = (max - min) / (steps - 1.0);
            Cells::Floats {
                values: ints.iter().map(|&i| min + i as f64 * scale).collect(),
                single: int_param(enc, "srcType")? == 32,
            }
        }
        other => return Err(format!("unknown encoding `{other}`")),
    })
}

fn string_array<'a>(indices: &[u8], enc: &Value<'a>) -> Result<Cells<'a>, String> {
    let chars = enc
        .get("stringData")
        .and_then(Value::as_str)
        .ok_or("StringArray lacks `stringData`")?;
    let offsets = decode_bytes(
        enc.get("offsets")
            .and_then(Value::as_bin)
            .ok_or("StringArray lacks `offsets`")?,
        enc.get("offsetEncoding").map_or(&[][..], Value::as_array),
    )?
    .ints()?;
    let indices = decode_bytes(
        indices,
        enc.get("dataEncoding").map_or(&[][..], Value::as_array),
    )?
    .ints()?;
    let piece = |index: i64| -> Result<Option<&'a str>, String> {
        let Ok(i) = usize::try_from(index) else {
            return Ok(None);
        };
        let bounds = |i: usize| usize::try_from(*offsets.get(i)?).ok();
        bounds(i)
            .zip(bounds(i + 1))
            .and_then(|(start, end)| chars.get(start..end))
            .map(Some)
            .ok_or_else(|| "StringArray offset outside its string data".to_string())
    };
    indices
        .into_iter()
        .map(piece)
        .collect::<Result<_, _>>()
        .map(Cells::Strings)
}

/// Undoes `encodings` (last to first) on `bytes`.
fn decode_bytes<'a>(bytes: &'a [u8], encodings: &[Value<'a>]) -> Result<Cells<'a>, String> {
    let mut cells: Option<Cells<'a>> = None;
    for enc in encodings.iter().rev() {
        let kind = enc
            .get("kind")
            .and_then(Value::as_str)
            .ok_or("encoding lacks `kind`")?;
        cells = Some(match (kind, cells.take()) {
            ("StringArray", None) => return string_array(bytes, enc),
            ("ByteArray", None) => from_bytes(bytes, int_param(enc, "type")?)?,
            (_, Some(prev)) => undo(prev, enc, kind)?,
            (_, None) => return Err(format!("`{kind}` needs a ByteArray beneath it")),
        });
    }
    cells.ok_or_else(|| "column has no encoding".to_string())
}

fn decode_data<'a>(data: &Value<'a>) -> Result<Cells<'a>, String> {
    let bytes = data
        .get("data")
        .and_then(Value::as_bin)
        .ok_or("column lacks `data`")?;
    decode_bytes(bytes, data.get("encoding").map_or(&[][..], Value::as_array))
}

struct Column<'a> {
    name: &'a str,
    cells: Cells<'a>,
    /// 0 present, 1 `.`, 2 `?`; `None` when every value is present.
    mask: Option<Vec<i64>>,
}

fn decode_column<'a>(column: &Value<'a>, rows: usize) -> Result<Column<'a>, String> {
    let name = column
        .get("name")
        .and_then(Value::as_str)
        .ok_or("column lacks a name")?;
    let cells = decode_data(column.get("data").ok_or("column lacks `data`")?)
        .map_err(|e| format!("column `{name}`: {e}"))?;
    let mask = match column.get("mask") {
        None | Some(Value::Nil) => None,
        Some(m) => Some(
            decode_data(m)
                .and_then(Cells::ints)
                .map_err(|e| format!("mask of `{name}`: {e}"))?,
        ),
    };
    if cells.len() != rows || mask.as_ref().is_some_and(|m| m.len() != rows) {
        return Err(format!("column `{name}` does not have {rows} rows"));
    }
    Ok(Column { name, cells, mask })
}

/// A CIF value token: bare when the syntax allows, else quoted, else a
/// text field.
fn push_value(out: &mut Vec<u8>, s: &str) {
    let lower = s.to_ascii_lowercase();
    let reserved = ["data_", "loop_", "save_", "global_", "stop_"]
        .iter()
        .any(|p| lower.starts_with(p));
    let bare = !s.is_empty()
        && !reserved
        && !s.bytes().any(|b| b.is_ascii_whitespace())
        && !matches!(
            s.as_bytes()[0],
            b'_' | b'#' | b'$' | b'\'' | b'"' | b';' | b'[' | b']'
        );
    if bare {
        out.extend_from_slice(s.as_bytes());
    } else if !s.contains(['\n', '\r']) && !s.contains('\'') {
        let _ = write!(out, "'{s}'");
    } else if !s.contains(['\n', '\r']) && !s.contains('"') {
        let _ = write!(out, "\"{s}\"");
    } else {
        let _ = write!(out, "\n;{s}\n;\n");
    }
}

fn push_cell(out: &mut Vec<u8>, column: &Column<'_>, row: usize) {
    match column.mask.as_ref().map_or(0, |m| m[row]) {
        0 => {}
        1 => return out.push(b'.'),
        _ => return out.push(b'?'),
    }
    match &column.cells {
        Cells::Ints(v) => {
            let _ = write!(out, "{}", v[row]);
        }
        Cells::Floats { values, single } => {
            let v = values[row];
            let _ = match (v.is_finite(), single) {
                (false, _) => write!(out, "?"),
                (true, true) => write!(out, "{}", v as f32),
                (true, false) => write!(out, "{v}"),
            };
        }
        Cells::Strings(v) => match v[row] {
            Some(s) => push_value(out, s),
            None => out.push(b'?'),
        },
    }
}

fn write_category(out: &mut Vec<u8>, category: &Value<'_>) -> Result<(), String> {
    let name = category
        .get("name")
        .and_then(Value::as_str)
        .ok_or("category lacks a name")?
        .trim_start_matches('_');
    let rows = category
        .get("rowCount")
        .and_then(Value::as_i64)
        .and_then(|n| usize::try_from(n).ok())
        .ok_or("category lacks `rowCount`")?;
    let columns = category
        .get("columns")
        .map_or(&[][..], Value::as_array)
        .iter()
        .map(|c| decode_column(c, rows))
        .collect::<Result<Vec<_>, _>>()
        .map_err(|e| format!("category `{name}`: {e}"))?;
    if columns.is_empty() || rows == 0 {
        return Ok(());
    }
    out.extend_from_slice(b"loop_\n");
    for column in &columns {
        let _ = writeln!(out, "_{name}.{}", column.name);
    }
    for row in 0..rows {
        for (i, column) in columns.iter().enumerate() {
            push_cell(out, column, row);
            out.push(if i + 1 == columns.len() { b'\n' } else { b' ' });
        }
    }
    out.extend_from_slice(b"#\n");
    Ok(())
}

fn has_atom_site(block: &Value<'_>) -> bool {
    block
        .get("categories")
        .map_or(&[][..], Value::as_array)
        .iter()
        .any(|c| {
            c.get("name")
                .and_then(Value::as_str)
                .map(|n| n.trim_start_matches('_'))
                == Some("atom_site")
        })
}

/// The first data block with an `_atom_site` table as CIF text.
pub(crate) fn to_cif(bytes: &[u8]) -> Result<Vec<u8>, ParseError> {
    let file = msgpack::decode(bytes).map_err(bad)?;
    let block = file
        .get("dataBlocks")
        .map_or(&[][..], Value::as_array)
        .iter()
        .find(|b| has_atom_site(b))
        .ok_or(ParseError::NoAtoms)?;
    let header = block
        .get("header")
        .and_then(Value::as_str)
        .unwrap_or("bcif");
    let mut out = Vec::with_capacity(bytes.len() * 8);
    let _ = writeln!(out, "data_{}", header.replace(char::is_whitespace, "_"));
    for category in block.get("categories").map_or(&[][..], Value::as_array) {
        write_category(&mut out, category).map_err(bad)?;
    }
    Ok(out)
}

pub fn parse(bytes: &[u8]) -> Result<vv_core::Structure, ParseError> {
    crate::mmcif::parse(&to_cif(bytes)?)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn unpacked(v: &[i64], unsigned: bool, bytes: i64) -> Vec<i64> {
        unpack_integers(v, unsigned, bytes).unwrap()
    }

    #[test]
    fn integer_packing_sums_saturated_runs() {
        assert_eq!(unpacked(&[1, 2, -3, 127, 1], false, 1), vec![1, 2, -3, 128]);
        assert_eq!(unpacked(&[255, 255, 3, 7], true, 1), vec![513, 7]);
        assert_eq!(unpacked(&[-128, -1], false, 1), vec![-129]);
        assert!(unpack_integers(&[127], false, 1).is_err());
    }

    #[test]
    fn run_length_and_delta_invert_their_examples() {
        assert_eq!(
            run_length(&[1, 3, 2, 1, 3, 2], 6).unwrap(),
            vec![1, 1, 1, 2, 3, 3]
        );
        assert!(run_length(&[1, 9], 6).is_err());
        assert_eq!(delta(vec![0, 3, 2, 1], 1000), vec![1000, 1003, 1005, 1006]);
    }

    #[test]
    fn fixed_point_and_interval_quantization_decode_to_floats() {
        let enc = Value::Map(vec![
            (Value::Str("factor"), Value::Int(100)),
            (Value::Str("srcType"), Value::Int(33)),
        ]);
        let Cells::Floats { values, .. } =
            undo(Cells::Ints(vec![120, 123, 12]), &enc, "FixedPoint").unwrap()
        else {
            panic!("floats expected");
        };
        assert_eq!(values, vec![1.2, 1.23, 0.12]);
        let enc = Value::Map(vec![
            (Value::Str("min"), Value::Int(1)),
            (Value::Str("max"), Value::Int(2)),
            (Value::Str("numSteps"), Value::Int(3)),
            (Value::Str("srcType"), Value::Int(33)),
        ]);
        let Cells::Floats { values, .. } =
            undo(Cells::Ints(vec![0, 1, 2]), &enc, "IntervalQuantization").unwrap()
        else {
            panic!("floats expected");
        };
        assert_eq!(values, vec![1.0, 1.5, 2.0]);
    }

    #[test]
    fn values_are_quoted_when_cif_syntax_needs_it() {
        let render = |s: &str| {
            let mut out = Vec::new();
            push_value(&mut out, s);
            String::from_utf8(out).unwrap()
        };
        assert_eq!(render("CA"), "CA");
        assert_eq!(render("O5'"), "O5'");
        assert_eq!(render("A B"), "'A B'");
        assert_eq!(render("it's here"), "\"it's here\"");
        assert_eq!(render("_x"), "'_x'");
        assert_eq!(render(""), "''");
        assert_eq!(render("a\nb"), "\n;a\nb\n;\n");
    }

    #[test]
    fn a_non_bcif_input_is_an_error_not_a_panic() {
        assert!(parse(b"data_x\n").is_err());
        assert!(parse(&[0x81, 0xa1, b'x', 0x01]).is_err());
    }
}
