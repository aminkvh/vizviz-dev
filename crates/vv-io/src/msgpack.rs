//! A MessagePack reader (msgpack.org specification): just enough to walk
//! a BinaryCIF file. Binary payloads borrow from the input.

const MAX_DEPTH: usize = 32;

#[derive(Debug, Clone, PartialEq)]
pub(crate) enum Value<'a> {
    Nil,
    Bool(bool),
    Int(i64),
    Float(f64),
    Str(&'a str),
    Bin(&'a [u8]),
    Array(Vec<Value<'a>>),
    Map(Vec<(Value<'a>, Value<'a>)>),
}

impl<'a> Value<'a> {
    /// The value under string key `key` of a map.
    pub(crate) fn get(&self, key: &str) -> Option<&Value<'a>> {
        match self {
            Value::Map(entries) => entries
                .iter()
                .find(|(k, _)| matches!(k, Value::Str(s) if *s == key))
                .map(|(_, v)| v),
            _ => None,
        }
    }

    pub(crate) fn as_str(&self) -> Option<&'a str> {
        match self {
            Value::Str(s) => Some(s),
            _ => None,
        }
    }

    pub(crate) fn as_i64(&self) -> Option<i64> {
        match self {
            Value::Int(i) => Some(*i),
            Value::Float(f) if f.fract() == 0.0 => Some(*f as i64),
            _ => None,
        }
    }

    pub(crate) fn as_f64(&self) -> Option<f64> {
        match self {
            Value::Int(i) => Some(*i as f64),
            Value::Float(f) => Some(*f),
            _ => None,
        }
    }

    pub(crate) fn as_bin(&self) -> Option<&'a [u8]> {
        match self {
            Value::Bin(b) => Some(b),
            _ => None,
        }
    }

    pub(crate) fn as_array(&self) -> &[Value<'a>] {
        match self {
            Value::Array(items) => items,
            _ => &[],
        }
    }
}

/// Decodes the single value at the start of `bytes`.
pub(crate) fn decode(bytes: &[u8]) -> Result<Value<'_>, String> {
    Reader { rest: bytes }.value(0)
}

struct Reader<'a> {
    rest: &'a [u8],
}

impl<'a> Reader<'a> {
    fn take(&mut self, n: usize) -> Result<&'a [u8], String> {
        if n > self.rest.len() {
            return Err("MessagePack data ends early".into());
        }
        let (head, tail) = self.rest.split_at(n);
        self.rest = tail;
        Ok(head)
    }

    fn be(&mut self, n: usize) -> Result<u64, String> {
        Ok(self
            .take(n)?
            .iter()
            .fold(0u64, |acc, &b| acc << 8 | b as u64))
    }

    fn len(&mut self, width: usize) -> Result<usize, String> {
        usize::try_from(self.be(width)?).map_err(|_| "MessagePack length too large".to_string())
    }

    fn str(&mut self, n: usize) -> Result<Value<'a>, String> {
        std::str::from_utf8(self.take(n)?)
            .map(Value::Str)
            .map_err(|_| "MessagePack string is not UTF-8".to_string())
    }

    fn array(&mut self, n: usize, depth: usize) -> Result<Value<'a>, String> {
        // Each element takes at least one byte, so a larger claim is corrupt.
        if n > self.rest.len() {
            return Err("MessagePack array longer than the data".into());
        }
        (0..n)
            .map(|_| self.value(depth + 1))
            .collect::<Result<_, _>>()
            .map(Value::Array)
    }

    fn map(&mut self, n: usize, depth: usize) -> Result<Value<'a>, String> {
        if n > self.rest.len() {
            return Err("MessagePack map longer than the data".into());
        }
        (0..n)
            .map(|_| Ok((self.value(depth + 1)?, self.value(depth + 1)?)))
            .collect::<Result<_, String>>()
            .map(Value::Map)
    }

    fn value(&mut self, depth: usize) -> Result<Value<'a>, String> {
        if depth > MAX_DEPTH {
            return Err("MessagePack nests too deeply".into());
        }
        let tag = self.take(1)?[0];
        match tag {
            0x00..=0x7f => Ok(Value::Int(tag as i64)),
            0x80..=0x8f => self.map((tag & 0x0f) as usize, depth),
            0x90..=0x9f => self.array((tag & 0x0f) as usize, depth),
            0xa0..=0xbf => self.str((tag & 0x1f) as usize),
            0xc0 => Ok(Value::Nil),
            0xc2 => Ok(Value::Bool(false)),
            0xc3 => Ok(Value::Bool(true)),
            0xc4..=0xc6 => {
                let n = self.len(1 << (tag - 0xc4))?;
                self.take(n).map(Value::Bin)
            }
            0xca => Ok(Value::Float(f32::from_bits(self.be(4)? as u32) as f64)),
            0xcb => Ok(Value::Float(f64::from_bits(self.be(8)?))),
            0xcc..=0xcf => {
                let v = self.be(1 << (tag - 0xcc))?;
                i64::try_from(v)
                    .map(Value::Int)
                    .map_err(|_| "MessagePack integer exceeds 63 bits".to_string())
            }
            0xd0 => Ok(Value::Int(self.be(1)? as u8 as i8 as i64)),
            0xd1 => Ok(Value::Int(self.be(2)? as u16 as i16 as i64)),
            0xd2 => Ok(Value::Int(self.be(4)? as u32 as i32 as i64)),
            0xd3 => Ok(Value::Int(self.be(8)? as i64)),
            0xd9..=0xdb => {
                let n = self.len(1 << (tag - 0xd9))?;
                self.str(n)
            }
            0xdc | 0xdd => {
                let n = self.len(if tag == 0xdc { 2 } else { 4 })?;
                self.array(n, depth)
            }
            0xde | 0xdf => {
                let n = self.len(if tag == 0xde { 2 } else { 4 })?;
                self.map(n, depth)
            }
            0xe0..=0xff => Ok(Value::Int(tag as i8 as i64)),
            _ => Err(format!("unsupported MessagePack type 0x{tag:02x}")),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn decodes_scalars_strings_and_containers() {
        // {"a": [1, -2, 300, 1.5], "b": nil, "c": bin[9]}
        let bytes = [
            0x83, 0xa1, b'a', 0x94, 0x01, 0xfe, 0xcd, 0x01, 0x2c, 0xcb, 0x3f, 0xf8, 0, 0, 0, 0, 0,
            0, 0xa1, b'b', 0xc0, 0xa1, b'c', 0xc4, 0x01, 0x09,
        ];
        let v = decode(&bytes).unwrap();
        let a = v.get("a").unwrap().as_array();
        assert_eq!(a[1].as_i64(), Some(-2));
        assert_eq!(a[2].as_i64(), Some(300));
        assert_eq!(a[3].as_f64(), Some(1.5));
        assert_eq!(v.get("b"), Some(&Value::Nil));
        assert_eq!(v.get("c").unwrap().as_bin(), Some(&[9u8][..]));
    }

    #[test]
    fn rejects_truncated_and_oversized_input() {
        assert!(decode(&[0xcd, 0x01]).is_err());
        assert!(decode(&[0xdd, 0xff, 0xff, 0xff, 0xff]).is_err());
        assert!(decode(&[0xc1]).is_err());
    }
}
