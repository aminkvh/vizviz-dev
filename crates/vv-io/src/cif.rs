//! CIF 1.1 syntax: a tokenizer and a small generic category reader.
//!
//! Written from the CIF specification (IUCr). The atom-site table of large
//! files is handled by the fast path in `mmcif.rs`; this module serves
//! everything else (headers, secondary structure, connectivity).

use std::collections::HashMap;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TokenKind {
    /// `_category.item`
    Tag,
    /// Any value: bare, quoted (quotes stripped), or a `;` text field.
    Value,
    Loop,
    /// `data_xxx`, `save_xxx`, `global_`, `stop_`
    Block,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Token<'a> {
    pub kind: TokenKind,
    pub text: &'a [u8],
    /// Byte offset of the token's first character.
    pub start: usize,
    /// Byte offset just past the token (past the closing quote/newline).
    pub end: usize,
}

pub struct Tokenizer<'a> {
    src: &'a [u8],
    pos: usize,
}

impl<'a> Tokenizer<'a> {
    pub fn new(src: &'a [u8]) -> Self {
        Self { src, pos: 0 }
    }

    pub fn position(&self) -> usize {
        self.pos
    }

    pub fn seek(&mut self, pos: usize) {
        self.pos = pos.min(self.src.len());
    }

    fn at_line_start(&self, pos: usize) -> bool {
        pos == 0 || self.src[pos - 1] == b'\n'
    }
}

impl<'a> Iterator for Tokenizer<'a> {
    type Item = Token<'a>;

    fn next(&mut self) -> Option<Token<'a>> {
        let src = self.src;
        loop {
            while self.pos < src.len() && src[self.pos].is_ascii_whitespace() {
                self.pos += 1;
            }
            if self.pos >= src.len() {
                return None;
            }
            let c = src[self.pos];
            if c == b'#' {
                while self.pos < src.len() && src[self.pos] != b'\n' {
                    self.pos += 1;
                }
                continue;
            }
            let start = self.pos;
            if c == b';' && self.at_line_start(start) {
                // Text field: from after ';' to the next line starting with ';'.
                let mut end = start + 1;
                loop {
                    match memchr::memchr(b'\n', &src[end..]) {
                        Some(nl) => {
                            end += nl + 1;
                            if end < src.len() && src[end] == b';' {
                                let text = &src[start + 1..end - 1];
                                self.pos = end + 1;
                                return Some(Token {
                                    kind: TokenKind::Value,
                                    text: trim_cr(text),
                                    start,
                                    end: self.pos,
                                });
                            }
                        }
                        None => {
                            self.pos = src.len();
                            return Some(Token {
                                kind: TokenKind::Value,
                                text: &src[start + 1..],
                                start,
                                end: src.len(),
                            });
                        }
                    }
                }
            }
            if c == b'\'' || c == b'"' {
                // Quoted: closes at a matching quote followed by whitespace/EOF.
                let mut end = start + 1;
                while end < src.len() {
                    if src[end] == c && (end + 1 >= src.len() || src[end + 1].is_ascii_whitespace())
                    {
                        self.pos = end + 1;
                        return Some(Token {
                            kind: TokenKind::Value,
                            text: &src[start + 1..end],
                            start,
                            end: self.pos,
                        });
                    }
                    if src[end] == b'\n' {
                        break;
                    }
                    end += 1;
                }
                // Unterminated quote: take the rest of the line as the value.
                self.pos = end;
                return Some(Token {
                    kind: TokenKind::Value,
                    text: &src[start..end],
                    start,
                    end,
                });
            }
            let mut end = start;
            while end < src.len() && !src[end].is_ascii_whitespace() {
                end += 1;
            }
            let text = &src[start..end];
            self.pos = end;
            let kind = if text[0] == b'_' {
                TokenKind::Tag
            } else if text.eq_ignore_ascii_case(b"loop_") {
                TokenKind::Loop
            } else if is_block_keyword(text) {
                TokenKind::Block
            } else {
                TokenKind::Value
            };
            return Some(Token {
                kind,
                text,
                start,
                end,
            });
        }
    }
}

fn is_block_keyword(t: &[u8]) -> bool {
    let lower = |p: &[u8]| t.len() >= p.len() && t[..p.len()].eq_ignore_ascii_case(p);
    lower(b"data_") || lower(b"save_") || lower(b"global_") || lower(b"stop_")
}

fn trim_cr(t: &[u8]) -> &[u8] {
    match t.last() {
        Some(b'\r') => &t[..t.len() - 1],
        _ => t,
    }
}

/// `.` and `?` mean "inapplicable" and "unknown".
pub fn is_null(v: &[u8]) -> bool {
    v == b"." || v == b"?"
}

/// A category's rows as raw byte slices, whether it came from a loop or
/// from single key/value items.
#[derive(Debug, Default, Clone)]
pub struct Category<'a> {
    /// Item names without the `_category.` prefix.
    pub items: Vec<&'a str>,
    pub rows: Vec<Vec<&'a [u8]>>,
}

impl<'a> Category<'a> {
    pub fn column(&self, item: &str) -> Option<usize> {
        self.items.iter().position(|i| i.eq_ignore_ascii_case(item))
    }

    pub fn get(&self, row: usize, item: &str) -> Option<&'a [u8]> {
        let col = self.column(item)?;
        self.rows
            .get(row)?
            .get(col)
            .copied()
            .filter(|v| !is_null(v))
    }

    pub fn get_str(&self, row: usize, item: &str) -> Option<&'a str> {
        std::str::from_utf8(self.get(row, item)?).ok()
    }
}

/// Location of a loop's data rows for a category the caller wants to
/// parse itself (the atom-site fast path).
#[derive(Debug, Clone)]
pub struct LoopLocation<'a> {
    pub items: Vec<&'a str>,
    pub body: std::ops::Range<usize>,
}

/// One data block's categories.
pub struct BlockRead<'a> {
    pub categories: HashMap<String, Category<'a>>,
    pub skipped: HashMap<String, LoopLocation<'a>>,
    /// Data blocks in the whole file up to and including the one read
    /// (fewer than `block + 1` means the block does not exist).
    pub blocks_seen: usize,
}

/// Reads every category of the first data block. See [`read_block`].
pub fn read_categories<'a>(
    src: &'a [u8],
    skip_body: &[&str],
) -> (
    HashMap<String, Category<'a>>,
    HashMap<String, LoopLocation<'a>>,
) {
    let r = read_block(src, skip_body, 0);
    (r.categories, r.skipped)
}

/// Reads every category of data block number `block` (0-based; content
/// before the first `data_` line counts as block 0). Categories listed in
/// `skip_body` are not materialized: their loop header is returned in
/// `skipped` together with the byte range of the rows, and the rows are
/// skipped line-by-line (they must be one row per line, which is true of
/// every PDB-issued file).
pub fn read_block<'a>(src: &'a [u8], skip_body: &[&str], block: usize) -> BlockRead<'a> {
    let mut categories: HashMap<String, Category<'a>> = HashMap::new();
    let mut skipped = HashMap::new();
    let mut tokens = Tokenizer::new(src);
    let mut blocks_seen = 0usize;
    let mut active = block == 0;

    while let Some(tok) = tokens.next() {
        match tok.kind {
            TokenKind::Block => {
                if tok.text.len() >= 5 && tok.text[..5].eq_ignore_ascii_case(b"data_") {
                    blocks_seen += 1;
                    if blocks_seen > block + 1 {
                        blocks_seen = block + 1;
                        break;
                    }
                    active = blocks_seen == block + 1;
                }
            }
            TokenKind::Tag => {
                // Key/value item.
                let (cat, item) = split_tag(tok.text);
                let value = match tokens.next() {
                    Some(v) if v.kind == TokenKind::Value => v.text,
                    Some(other) => {
                        tokens.seek(other.start);
                        b"?"
                    }
                    None => b"?",
                };
                if !active {
                    continue;
                }
                let entry = categories.entry(cat.to_string()).or_default();
                if entry.rows.is_empty() {
                    entry.rows.push(Vec::new());
                }
                entry.items.push(item);
                entry.rows[0].push(value);
            }
            TokenKind::Loop => {
                let mut items: Vec<&str> = Vec::new();
                let mut cat = "";
                let mut body_start = tokens.position();
                loop {
                    let save = tokens.position();
                    match tokens.next() {
                        Some(t) if t.kind == TokenKind::Tag => {
                            let (c, i) = split_tag(t.text);
                            cat = c;
                            items.push(i);
                            body_start = t.end;
                        }
                        Some(t) => {
                            tokens.seek(save);
                            body_start = t.start.min(body_start.max(save));
                            let _ = body_start;
                            body_start = save;
                            break;
                        }
                        None => break,
                    }
                }
                if items.is_empty() {
                    continue;
                }
                let skip = active && skip_body.iter().any(|s| s.eq_ignore_ascii_case(cat));
                if skip {
                    let start = skip_whitespace_and_comments(src, body_start);
                    let end = find_loop_end(src, start);
                    skipped.insert(
                        cat.to_string(),
                        LoopLocation {
                            items,
                            body: start..end,
                        },
                    );
                    tokens.seek(end);
                    continue;
                }
                let mut rows = Vec::new();
                let mut row = Vec::with_capacity(items.len());
                loop {
                    let save = tokens.position();
                    match tokens.next() {
                        Some(t) if t.kind == TokenKind::Value => {
                            row.push(t.text);
                            if row.len() == items.len() {
                                rows.push(std::mem::take(&mut row));
                            }
                        }
                        Some(_) => {
                            tokens.seek(save);
                            break;
                        }
                        None => break,
                    }
                }
                if !active {
                    continue;
                }
                let entry = categories.entry(cat.to_string()).or_default();
                entry.items = items;
                entry.rows = rows;
            }
            TokenKind::Value => {}
        }
    }
    BlockRead {
        categories,
        skipped,
        blocks_seen,
    }
}

fn split_tag(tag: &[u8]) -> (&str, &str) {
    let s = std::str::from_utf8(tag).unwrap_or("_");
    let s = s.strip_prefix('_').unwrap_or(s);
    match s.split_once('.') {
        Some((cat, item)) => (cat, item),
        None => (s, ""),
    }
}

fn skip_whitespace_and_comments(src: &[u8], mut pos: usize) -> usize {
    loop {
        while pos < src.len() && src[pos].is_ascii_whitespace() {
            pos += 1;
        }
        if pos < src.len() && src[pos] == b'#' {
            while pos < src.len() && src[pos] != b'\n' {
                pos += 1;
            }
            continue;
        }
        return pos;
    }
}

/// End of a loop body that has one row per line: the first line that
/// starts a tag, a keyword, a comment, or a text field.
fn find_loop_end(src: &[u8], start: usize) -> usize {
    let mut pos = start;
    while pos < src.len() {
        let line_end = memchr::memchr(b'\n', &src[pos..]).map_or(src.len(), |n| pos + n);
        let line = &src[pos..line_end];
        let first = line.iter().position(|b| !b.is_ascii_whitespace());
        match first {
            None => {}
            Some(i) => {
                let rest = &line[i..];
                let c = rest[0];
                if c == b'_'
                    || c == b'#'
                    || (i == 0 && c == b';')
                    || rest.len() >= 5 && rest[..5].eq_ignore_ascii_case(b"loop_")
                    || is_block_keyword(rest)
                {
                    return pos;
                }
            }
        }
        pos = line_end + 1;
    }
    src.len()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tokenizes_quotes_text_fields_and_comments() {
        let src =
            b"data_x\n_a.b 'hello world' # comment\n_a.c \"it's\"\n_a.d\n;multi\nline\n;\n_a.e 5\n";
        let kinds: Vec<(TokenKind, &[u8])> =
            Tokenizer::new(src).map(|t| (t.kind, t.text)).collect();
        assert_eq!(
            kinds,
            vec![
                (TokenKind::Block, &b"data_x"[..]),
                (TokenKind::Tag, b"_a.b"),
                (TokenKind::Value, b"hello world"),
                (TokenKind::Tag, b"_a.c"),
                (TokenKind::Value, b"it's"),
                (TokenKind::Tag, b"_a.d"),
                (TokenKind::Value, b"multi\nline"),
                (TokenKind::Tag, b"_a.e"),
                (TokenKind::Value, b"5"),
            ]
        );
    }

    #[test]
    fn apostrophe_inside_bare_value_is_literal() {
        let src = b"_a.x C1' _a.y O5'\n";
        let vals: Vec<&[u8]> = Tokenizer::new(src)
            .filter(|t| t.kind == TokenKind::Value)
            .map(|t| t.text)
            .collect();
        assert_eq!(vals, vec![&b"C1'"[..], b"O5'"]);
    }

    #[test]
    fn reads_loops_and_items_and_skips_requested_bodies() {
        let src = b"data_t\n_entry.id  1ABC\nloop_\n_s.a\n_s.b\n1 x\n2 'y z'\n#\nloop_\n_atom_site.id\n_atom_site.x\n1 0.5\n2 1.5\n#\n_struct.title 'T'\n";
        let (cats, skipped) = read_categories(src, &["atom_site"]);
        assert_eq!(cats["entry"].get_str(0, "id"), Some("1ABC"));
        assert_eq!(cats["struct"].get_str(0, "title"), Some("T"));
        let s = &cats["s"];
        assert_eq!(s.items, vec!["a", "b"]);
        assert_eq!(s.rows, vec![vec![&b"1"[..], b"x"], vec![b"2", b"y z"]]);
        let atoms = &skipped["atom_site"];
        assert_eq!(atoms.items, vec!["id", "x"]);
        assert_eq!(&src[atoms.body.clone()], b"1 0.5\n2 1.5\n");
    }

    #[test]
    fn null_values_read_as_none() {
        let src = b"data_t\n_a.b ?\n_a.c .\n_a.d 3\n";
        let (cats, _) = read_categories(src, &[]);
        assert_eq!(cats["a"].get(0, "b"), None);
        assert_eq!(cats["a"].get(0, "c"), None);
        assert_eq!(cats["a"].get(0, "d"), Some(&b"3"[..]));
    }
}
