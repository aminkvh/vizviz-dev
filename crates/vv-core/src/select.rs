//! Selection expressions: `protein and within 4 of resname HEM`.
//!
//! A recursive-descent parser turns the text into an [`Expr`] tree, and
//! evaluation produces one bit per atom. Keywords are case-insensitive;
//! `not` binds tighter than `and`, which binds tighter than `or`. Errors
//! carry a byte span so a UI can underline the offending token.
//!
//! Residue-level tests (`protein`, `chain`, `resid`, ...) decide once per
//! residue or chain and fill whole atom ranges, so they cost O(residues)
//! rather than a string compare per atom.

use std::ops::Range;

use fixedbitset::FixedBitSet;
use glam::Vec3;
use rayon::prelude::*;

use crate::{
    flags, Element, Grid, InternId, ResidueClass, ResidueRec, Roles, SecondaryStructure, Topology,
};

/// Parses and evaluates `expr` in one go.
pub fn select(
    topology: &Topology,
    positions: &[Vec3],
    expr: &str,
) -> Result<FixedBitSet, SelectError> {
    Ok(parse(expr)?.evaluate(topology, positions))
}

#[derive(thiserror::Error, Debug, Clone, PartialEq)]
#[error("{message}")]
pub struct SelectError {
    pub message: String,
    /// Byte range in the source expression.
    pub span: Range<usize>,
}

impl SelectError {
    fn new(message: impl Into<String>, span: Range<usize>) -> Self {
        Self {
            message: message.into(),
            span,
        }
    }
}

// ---------------------------------------------------------------------------
// AST

/// Argument-free keywords.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Class {
    All,
    None,
    /// The residue-class keywords select residues whose
    /// `Topology::residue_class` is that class (`crate::residue_class`).
    Protein,
    Nucleic,
    Lipid,
    Glycan,
    Water,
    Ion,
    /// The residue-role keywords select residues holding that role
    /// (`Topology::residue_roles`); they overlap the class keywords, so a
    /// lipid in a pocket is both `lipid` and `ligand`.
    Ligand,
    Membrane,
    Additive,
    Cofactor,
    /// Water or ion.
    Solvent,
    /// Protein or nucleic.
    Polymer,
    /// HETATM records: a file-level flag, unrelated to the classes above.
    Hetero,
    Backbone,
    Sidechain,
    Hydrogen,
    Helix,
    Strand,
    Coil,
}

/// Integer fields selectable by inclusive ranges.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Field {
    /// Author residue number (what papers cite).
    Resid,
    /// Label (sequence-position) residue number.
    Seqid,
    /// Atom serial from the file.
    Serial,
    /// 0-based atom index.
    Index,
}

/// Float per-atom columns selectable by comparison.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Column {
    Bfactor,
    Occupancy,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Cmp {
    Lt,
    Le,
    Gt,
    Ge,
    Eq,
    Ne,
}

#[derive(Clone, Debug, PartialEq)]
pub enum Expr {
    Class(Class),
    /// Matched exactly against the label *or* auth chain id.
    Chain(Vec<String>),
    /// Segment id, matched exactly.
    Segname(Vec<String>),
    /// Case-insensitive.
    Resname(Vec<String>),
    /// Trimmed, case-insensitive.
    Name(Vec<String>),
    Element(Vec<Element>),
    /// Single alternate-location characters, matched exactly.
    Altloc(Vec<u8>),
    /// Inclusive `(lo, hi)` ranges.
    Range(Field, Vec<(i64, i64)>),
    Compare(Column, Cmp, f32),
    /// Atoms within the distance (inclusive) of any atom of the inner set;
    /// the inner atoms themselves are included.
    Within(f32, Box<Expr>),
    /// Every atom of every residue touched by the inner set.
    Byres(Box<Expr>),
    Not(Box<Expr>),
    And(Box<Expr>, Box<Expr>),
    Or(Box<Expr>, Box<Expr>),
}

// ---------------------------------------------------------------------------
// Lexer

#[derive(Clone, Copy, Debug, PartialEq)]
enum Tok<'a> {
    Word(&'a str),
    LParen,
    RParen,
    Cmp(Cmp),
}

#[derive(Clone, Debug)]
struct Token<'a> {
    tok: Tok<'a>,
    span: Range<usize>,
}

fn is_delimiter(b: u8) -> bool {
    b.is_ascii_whitespace() || matches!(b, b'(' | b')' | b'<' | b'>' | b'=' | b'!')
}

/// Splits on whitespace, parentheses and comparison operators. Everything
/// else is a word; keywords and numbers are recognised by the parser.
fn lex(src: &str) -> Result<Vec<Token<'_>>, SelectError> {
    let bytes = src.as_bytes();
    let mut tokens = Vec::new();
    let mut i = 0;
    while i < bytes.len() {
        let start = i;
        let tok = match bytes[i] {
            b if b.is_ascii_whitespace() => {
                i += 1;
                continue;
            }
            b'(' => {
                i += 1;
                Tok::LParen
            }
            b')' => {
                i += 1;
                Tok::RParen
            }
            b @ (b'<' | b'>' | b'=' | b'!') => {
                let two = bytes.get(i + 1) == Some(&b'=');
                i += if two { 2 } else { 1 };
                let cmp = match (b, two) {
                    (b'<', false) => Cmp::Lt,
                    (b'<', true) => Cmp::Le,
                    (b'>', false) => Cmp::Gt,
                    (b'>', true) => Cmp::Ge,
                    (b'=', _) => Cmp::Eq,
                    (b'!', true) => Cmp::Ne,
                    _ => {
                        return Err(SelectError::new(
                            "unexpected `!`; did you mean `!=`?",
                            start..i,
                        ))
                    }
                };
                Tok::Cmp(cmp)
            }
            _ => {
                while i < bytes.len() && !is_delimiter(bytes[i]) {
                    i += 1;
                }
                Tok::Word(&src[start..i])
            }
        };
        tokens.push(Token {
            tok,
            span: start..i,
        });
    }
    Ok(tokens)
}

// ---------------------------------------------------------------------------
// Parser

#[derive(Clone, Copy, Debug, PartialEq)]
enum Kw {
    And,
    Or,
    Not,
    Of,
    Class(Class),
    Chain,
    Segname,
    Resname,
    Name,
    Element,
    Altloc,
    Range(Field),
    Compare(Column),
    Within,
    Byres,
}

fn keyword(word: &str) -> Option<Kw> {
    Some(match word.to_ascii_lowercase().as_str() {
        "and" => Kw::And,
        "or" => Kw::Or,
        "not" => Kw::Not,
        "of" => Kw::Of,
        "all" => Kw::Class(Class::All),
        "none" => Kw::Class(Class::None),
        "protein" => Kw::Class(Class::Protein),
        "nucleic" => Kw::Class(Class::Nucleic),
        "lipid" => Kw::Class(Class::Lipid),
        "membrane" => Kw::Class(Class::Membrane),
        "additive" => Kw::Class(Class::Additive),
        "cofactor" => Kw::Class(Class::Cofactor),
        "water" => Kw::Class(Class::Water),
        "ion" => Kw::Class(Class::Ion),
        "glycan" => Kw::Class(Class::Glycan),
        "ligand" => Kw::Class(Class::Ligand),
        "solvent" => Kw::Class(Class::Solvent),
        "polymer" => Kw::Class(Class::Polymer),
        "hetero" => Kw::Class(Class::Hetero),
        "backbone" => Kw::Class(Class::Backbone),
        "sidechain" => Kw::Class(Class::Sidechain),
        "hydrogen" => Kw::Class(Class::Hydrogen),
        "helix" => Kw::Class(Class::Helix),
        "strand" => Kw::Class(Class::Strand),
        "coil" => Kw::Class(Class::Coil),
        "chain" => Kw::Chain,
        "segname" | "segid" => Kw::Segname,
        "resname" => Kw::Resname,
        "name" => Kw::Name,
        "element" => Kw::Element,
        "altloc" => Kw::Altloc,
        "resid" => Kw::Range(Field::Resid),
        "seqid" => Kw::Range(Field::Seqid),
        "serial" => Kw::Range(Field::Serial),
        "index" => Kw::Range(Field::Index),
        "bfactor" => Kw::Compare(Column::Bfactor),
        "occupancy" => Kw::Compare(Column::Occupancy),
        "within" => Kw::Within,
        "byres" => Kw::Byres,
        _ => return None,
    })
}

/// Parses an expression without evaluating it.
pub fn parse(src: &str) -> Result<Expr, SelectError> {
    let tokens = lex(src)?;
    if tokens.is_empty() {
        return Err(SelectError::new("empty selection", 0..src.len()));
    }
    let mut parser = Parser {
        src,
        tokens,
        pos: 0,
    };
    let expr = parser.or()?;
    if let Some(t) = parser.peek() {
        return Err(SelectError::new(
            format!(
                "unexpected `{}`; expected `and`, `or` or the end",
                parser.text(t)
            ),
            t.span.start..src.len(),
        ));
    }
    Ok(expr)
}

struct Parser<'a> {
    src: &'a str,
    tokens: Vec<Token<'a>>,
    pos: usize,
}

type Parsed<T> = Result<T, SelectError>;

impl<'a> Parser<'a> {
    fn peek(&self) -> Option<&Token<'a>> {
        self.tokens.get(self.pos)
    }

    fn next(&mut self) -> Option<Token<'a>> {
        let t = self.tokens.get(self.pos).cloned();
        self.pos += t.is_some() as usize;
        t
    }

    fn peek_kw(&self) -> Option<Kw> {
        match self.peek()?.tok {
            Tok::Word(w) => keyword(w),
            _ => None,
        }
    }

    fn eat_kw(&mut self, kw: Kw) -> bool {
        let hit = self.peek_kw() == Some(kw);
        self.pos += hit as usize;
        hit
    }

    fn text(&self, t: &Token<'a>) -> &'a str {
        &self.src[t.span.clone()]
    }

    /// "expected X after `after`", pointing at the token found instead, or
    /// at `after` itself when the input ends there.
    fn expected(&self, after: &Token<'a>, what: &str) -> SelectError {
        let (found, span) = match self.peek() {
            Some(t) => (format!(", found `{}`", self.text(t)), t.span.clone()),
            None => (String::new(), after.span.clone()),
        };
        let message = format!("expected {what} after `{}`{found}", self.text(after));
        SelectError::new(message, span)
    }

    fn or(&mut self) -> Parsed<Expr> {
        let mut lhs = self.and()?;
        while self.eat_kw(Kw::Or) {
            lhs = Expr::Or(Box::new(lhs), Box::new(self.and()?));
        }
        Ok(lhs)
    }

    fn and(&mut self) -> Parsed<Expr> {
        let mut lhs = self.not()?;
        while self.eat_kw(Kw::And) {
            lhs = Expr::And(Box::new(lhs), Box::new(self.not()?));
        }
        Ok(lhs)
    }

    fn not(&mut self) -> Parsed<Expr> {
        if self.eat_kw(Kw::Not) {
            Ok(Expr::Not(Box::new(self.not()?)))
        } else {
            self.primary()
        }
    }

    fn primary(&mut self) -> Parsed<Expr> {
        let Some(tok) = self.next() else {
            // `parse` rejects empty input, so a previous token exists.
            let prev = self.tokens[self.pos - 1].clone();
            return Err(self.expected(&prev, "a selection"));
        };
        let word = match tok.tok {
            Tok::LParen => return self.parenthesised(&tok),
            Tok::RParen => return Err(SelectError::new("unexpected `)`", tok.span)),
            Tok::Cmp(_) => {
                let message = format!("unexpected `{}`", self.text(&tok));
                return Err(SelectError::new(message, tok.span));
            }
            Tok::Word(w) => w,
        };
        match keyword(word) {
            None => Err(SelectError::new(
                format!("unknown keyword `{word}`"),
                tok.span,
            )),
            Some(Kw::And | Kw::Or | Kw::Not | Kw::Of) => Err(SelectError::new(
                format!("unexpected `{word}`; expected a selection"),
                tok.span,
            )),
            Some(Kw::Class(c)) => Ok(Expr::Class(c)),
            Some(Kw::Chain) => Ok(Expr::Chain(self.strings(&tok)?)),
            Some(Kw::Segname) => Ok(Expr::Segname(self.strings(&tok)?)),
            Some(Kw::Resname) => Ok(Expr::Resname(self.strings(&tok)?)),
            Some(Kw::Name) => Ok(Expr::Name(self.strings(&tok)?)),
            Some(Kw::Element) => self.elements(&tok),
            Some(Kw::Altloc) => self.altlocs(&tok),
            Some(Kw::Range(field)) => self.ranges(&tok, field),
            Some(Kw::Compare(column)) => self.compare(&tok, column),
            Some(Kw::Within) => self.within(&tok),
            // The operand binds tighter than `and`/`or`, like `not`.
            Some(Kw::Byres) => Ok(Expr::Byres(Box::new(self.not()?))),
        }
    }

    fn parenthesised(&mut self, open: &Token<'a>) -> Parsed<Expr> {
        let inner = self.or()?;
        match self.next() {
            Some(Token {
                tok: Tok::RParen, ..
            }) => Ok(inner),
            Some(t) => Err(SelectError::new(
                format!("expected `)`, found `{}`", self.text(&t)),
                t.span,
            )),
            None => Err(SelectError::new("unclosed `(`", open.span.clone())),
        }
    }

    /// The value tokens after a keyword: every word up to the next keyword,
    /// operator, parenthesis or the end. At least one is required.
    fn values(&mut self, kw: &Token<'a>) -> Parsed<Vec<Token<'a>>> {
        let mut out = Vec::new();
        while let Some(t) = self.peek() {
            match t.tok {
                Tok::Word(w) if keyword(w).is_none() => out.push(t.clone()),
                _ => break,
            }
            self.pos += 1;
        }
        if out.is_empty() {
            return Err(self.expected(kw, "one or more values"));
        }
        Ok(out)
    }

    fn strings(&mut self, kw: &Token<'a>) -> Parsed<Vec<String>> {
        Ok(self
            .values(kw)?
            .iter()
            .map(|t| self.text(t).to_owned())
            .collect())
    }

    fn elements(&mut self, kw: &Token<'a>) -> Parsed<Expr> {
        let elements = self
            .values(kw)?
            .iter()
            .map(|t| {
                let symbol = self.text(t);
                match Element::from_symbol(symbol.as_bytes()) {
                    Element::UNKNOWN => Err(SelectError::new(
                        format!("unknown element symbol `{symbol}`"),
                        t.span.clone(),
                    )),
                    e => Ok(e),
                }
            })
            .collect::<Parsed<_>>()?;
        Ok(Expr::Element(elements))
    }

    fn altlocs(&mut self, kw: &Token<'a>) -> Parsed<Expr> {
        let chars = self
            .values(kw)?
            .iter()
            .map(|t| match self.text(t).as_bytes() {
                [c] => Ok(*c),
                _ => Err(SelectError::new(
                    "altloc takes single characters, e.g. `altloc A B`",
                    t.span.clone(),
                )),
            })
            .collect::<Parsed<_>>()?;
        Ok(Expr::Altloc(chars))
    }

    fn ranges(&mut self, kw: &Token<'a>, field: Field) -> Parsed<Expr> {
        let ranges = self
            .values(kw)?
            .iter()
            .map(|t| parse_range(self.text(t), t.span.clone()))
            .collect::<Parsed<_>>()?;
        Ok(Expr::Range(field, ranges))
    }

    fn compare(&mut self, kw: &Token<'a>, column: Column) -> Parsed<Expr> {
        let op = self.peek().cloned();
        let Some((cmp, op)) = op.and_then(|t| match t.tok {
            Tok::Cmp(c) => Some((c, t)),
            _ => None,
        }) else {
            return Err(self.expected(kw, "a comparison (<, <=, >, >=, ==, !=)"));
        };
        self.pos += 1;
        let (value, _) = self.number(&op)?;
        Ok(Expr::Compare(column, cmp, value))
    }

    fn within(&mut self, kw: &Token<'a>) -> Parsed<Expr> {
        let (radius, number) = self.number(kw)?;
        if radius < 0.0 {
            return Err(SelectError::new(
                "distance must not be negative",
                number.span,
            ));
        }
        if !self.eat_kw(Kw::Of) {
            return Err(self.expected(&number, "`of`"));
        }
        Ok(Expr::Within(radius, Box::new(self.not()?)))
    }

    /// A finite number token; also returned so callers can point at it.
    fn number(&mut self, after: &Token<'a>) -> Parsed<(f32, Token<'a>)> {
        let Some(t) = self.peek().cloned() else {
            return Err(self.expected(after, "a number"));
        };
        let value = match t.tok {
            Tok::Word(w) => w.parse::<f32>().ok().filter(|v| v.is_finite()),
            _ => None,
        };
        match value {
            Some(v) => {
                self.pos += 1;
                Ok((v, t))
            }
            None => Err(self.expected(after, "a number")),
        }
    }
}

/// `5`, `1-10`, `1:10`, `-3`, `-3--1`. A `-` in first position is a sign;
/// the first later `-` or `:` separates the two ends.
fn parse_range(word: &str, span: Range<usize>) -> Parsed<(i64, i64)> {
    let bad = || {
        SelectError::new(
            format!("bad range `{word}`; use N, N-M or N:M"),
            span.clone(),
        )
    };
    let sep = word
        .char_indices()
        .skip(1)
        .find(|&(_, c)| c == '-' || c == ':')
        .map(|(i, _)| i);
    let (lo, hi) = match sep {
        Some(i) => (&word[..i], &word[i + 1..]),
        None => (word, word),
    };
    let lo: i64 = lo.parse().map_err(|_| bad())?;
    let hi: i64 = hi.parse().map_err(|_| bad())?;
    if lo > hi {
        let message = format!("range `{word}` runs backwards; write {hi}-{lo}");
        return Err(SelectError::new(message, span));
    }
    Ok((lo, hi))
}

// ---------------------------------------------------------------------------
// Evaluation

const PROTEIN_BACKBONE: &[&str] = &["N", "CA", "C", "O", "OXT", "OT1", "OT2"];
const NUCLEIC_BACKBONE: &[&str] = &[
    "P", "OP1", "OP2", "O1P", "O2P", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'",
];

fn in_table(table: &[&str], name: &str) -> bool {
    table.iter().any(|t| t.eq_ignore_ascii_case(name))
}

impl Expr {
    /// One bit per atom. `positions` is only read by `within`, and must
    /// then hold one entry per atom.
    pub fn evaluate(&self, t: &Topology, positions: &[Vec3]) -> FixedBitSet {
        let n = t.atom_count();
        match self {
            Expr::Class(class) => class_mask(t, *class),
            Expr::Chain(ids) => chain_mask(t, ids),
            Expr::Segname(ids) => segname_mask(t, ids),
            Expr::Resname(names) => {
                comp_mask(t, |comp| names.iter().any(|w| w.eq_ignore_ascii_case(comp)))
            }
            Expr::Name(names) => atom_mask(n, |i| {
                let name = t.atom_name(i);
                names.iter().any(|w| w.eq_ignore_ascii_case(name))
            }),
            Expr::Element(elements) => atom_mask(n, |i| elements.contains(&t.element[i])),
            Expr::Altloc(chars) => {
                atom_mask(n, |i| t.alt_loc.get(i).is_some_and(|a| chars.contains(a)))
            }
            Expr::Range(field, ranges) => range_mask(t, *field, ranges),
            Expr::Compare(column, cmp, value) => {
                let column = match column {
                    Column::Bfactor => &t.b_factor,
                    Column::Occupancy => &t.occupancy,
                };
                atom_mask(n, |i| column.get(i).is_some_and(|&x| cmp.test(x, *value)))
            }
            Expr::Within(radius, inner) => {
                within(t, positions, *radius, &inner.evaluate(t, positions))
            }
            Expr::Byres(inner) => byres(t, &inner.evaluate(t, positions)),
            Expr::Not(inner) => {
                let mut m = inner.evaluate(t, positions);
                m.toggle_range(..);
                m
            }
            Expr::And(a, b) => {
                let mut m = a.evaluate(t, positions);
                m.intersect_with(&b.evaluate(t, positions));
                m
            }
            Expr::Or(a, b) => {
                let mut m = a.evaluate(t, positions);
                m.union_with(&b.evaluate(t, positions));
                m
            }
        }
    }
}

impl Cmp {
    fn test(self, x: f32, v: f32) -> bool {
        match self {
            Cmp::Lt => x < v,
            Cmp::Le => x <= v,
            Cmp::Gt => x > v,
            Cmp::Ge => x >= v,
            Cmp::Eq => x == v,
            Cmp::Ne => x != v,
        }
    }
}

fn atom_mask(n: usize, pred: impl Fn(usize) -> bool) -> FixedBitSet {
    let mut m = FixedBitSet::with_capacity(n);
    m.extend((0..n).filter(|&i| pred(i)));
    m
}

/// Whole-residue fill for every residue passing `pred(index, residue)`.
fn residue_mask(t: &Topology, pred: impl Fn(usize, &ResidueRec) -> bool) -> FixedBitSet {
    let mut m = FixedBitSet::with_capacity(t.atom_count());
    for (i, r) in t.residues.iter().enumerate() {
        if pred(i, r) {
            m.insert_range(r.atoms.start as usize..r.atoms.end as usize);
        }
    }
    m
}

/// One flag per interned name, so residue tests compare an id, not a string.
fn name_table(t: &Topology, pred: impl Fn(&str) -> bool) -> Vec<bool> {
    (0..t.names.len())
        .map(|i| pred(t.names.get(InternId(i as u32))))
        .collect()
}

/// Residues whose component name passes `pred`.
fn comp_mask(t: &Topology, pred: impl Fn(&str) -> bool) -> FixedBitSet {
    let table = name_table(t, pred);
    residue_mask(t, |_, r| table[r.comp.0 as usize])
}

fn chain_mask(t: &Topology, ids: &[String]) -> FixedBitSet {
    let wanted: Vec<InternId> = ids.iter().filter_map(|id| t.names.lookup(id)).collect();
    let hit: Vec<bool> = t
        .chains
        .iter()
        .map(|c| wanted.contains(&c.label_asym) || wanted.contains(&c.auth_asym))
        .collect();
    residue_mask(t, |_, r| hit[r.chain as usize])
}

fn segname_mask(t: &Topology, ids: &[String]) -> FixedBitSet {
    let hit: Vec<bool> = (0..t.chain_count())
        .map(|c| ids.iter().any(|id| id == t.segid(c)))
        .collect();
    residue_mask(t, |_, r| hit[r.chain as usize])
}

fn class_mask(t: &Topology, class: Class) -> FixedBitSet {
    let n = t.atom_count();
    match class {
        Class::All => {
            let mut m = FixedBitSet::with_capacity(n);
            m.insert_range(..);
            m
        }
        Class::None => FixedBitSet::with_capacity(n),
        Class::Protein => residue_class_mask(t, &[ResidueClass::Protein]),
        Class::Nucleic => residue_class_mask(t, &[ResidueClass::Nucleic]),
        Class::Lipid => residue_class_mask(t, &[ResidueClass::Lipid]),
        Class::Glycan => residue_class_mask(t, &[ResidueClass::Glycan]),
        Class::Water => residue_class_mask(t, &[ResidueClass::Water]),
        Class::Ion => residue_class_mask(t, &[ResidueClass::Ion]),
        Class::Ligand => role_mask(t, Roles::LIGAND),
        Class::Membrane => role_mask(t, Roles::MEMBRANE),
        Class::Additive => role_mask(t, Roles::ADDITIVE),
        Class::Cofactor => role_mask(t, Roles::COFACTOR),
        Class::Solvent => residue_class_mask(t, &[ResidueClass::Water, ResidueClass::Ion]),
        Class::Polymer => residue_class_mask(t, &[ResidueClass::Protein, ResidueClass::Nucleic]),
        Class::Hetero => atom_mask(n, |i| {
            t.flags.get(i).is_some_and(|f| (f & flags::HETERO) != 0)
        }),
        Class::Backbone => backbone(t),
        Class::Sidechain => {
            let mut m = class_mask(t, Class::Protein);
            m.difference_with(&backbone(t));
            m
        }
        Class::Hydrogen => atom_mask(n, |i| is_hydrogen(t, i)),
        Class::Helix => residue_mask(t, |_, r| r.ss == SecondaryStructure::Helix),
        Class::Strand => residue_mask(t, |_, r| r.ss == SecondaryStructure::Strand),
        Class::Coil => residue_mask(t, |_, r| r.ss == SecondaryStructure::Coil),
    }
}

fn residue_class_mask(t: &Topology, classes: &[ResidueClass]) -> FixedBitSet {
    residue_mask(t, |i, _| classes.contains(&t.residue_class(i)))
}

fn role_mask(t: &Topology, role: Roles) -> FixedBitSet {
    residue_mask(t, |i, _| t.residue_roles(i).contains(role))
}

/// Backbone atoms of protein and nucleic residues, by atom name.
fn backbone(t: &Topology) -> FixedBitSet {
    let mut m = FixedBitSet::with_capacity(t.atom_count());
    for (i, r) in t.residues.iter().enumerate() {
        let names = match t.residue_class(i) {
            ResidueClass::Protein => PROTEIN_BACKBONE,
            ResidueClass::Nucleic => NUCLEIC_BACKBONE,
            _ => continue,
        };
        m.extend(
            r.atoms
                .clone()
                .map(|a| a as usize)
                .filter(|&a| in_table(names, t.atom_name(a))),
        );
    }
    m
}

/// Element H, or an unknown element whose name starts with D: deuterium
/// has no `Element` yet, so neutron structures land here.
fn is_hydrogen(t: &Topology, i: usize) -> bool {
    let e = t.element[i];
    e.is_hydrogen()
        || (e.is_unknown()
            && t.atom_name(i)
                .as_bytes()
                .first()
                .is_some_and(|b| b.eq_ignore_ascii_case(&b'D')))
}

fn range_mask(t: &Topology, field: Field, ranges: &[(i64, i64)]) -> FixedBitSet {
    let n = t.atom_count();
    let hit = |v: i64| ranges.iter().any(|&(lo, hi)| lo <= v && v <= hi);
    match field {
        Field::Resid => residue_mask(t, |_, r| hit(r.auth_seq_id as i64)),
        Field::Seqid => residue_mask(t, |_, r| hit(r.seq_id as i64)),
        Field::Serial => atom_mask(n, |i| t.serial.get(i).is_some_and(|&s| hit(s as i64))),
        Field::Index => {
            let mut m = FixedBitSet::with_capacity(n);
            for &(lo, hi) in ranges {
                let lo = lo.clamp(0, n as i64) as usize;
                let hi = hi.saturating_add(1).clamp(0, n as i64) as usize;
                if lo < hi {
                    m.insert_range(lo..hi);
                }
            }
            m
        }
    }
}

/// Grid over the inner set, then every atom is tested in parallel.
fn within(t: &Topology, positions: &[Vec3], radius: f32, inner: &FixedBitSet) -> FixedBitSet {
    let n = t.atom_count();
    if radius <= 0.0 || inner.is_clear() {
        return inner.clone();
    }
    assert_eq!(positions.len(), n, "`within` needs one position per atom");
    let indices: Vec<u32> = inner.ones().map(|i| i as u32).collect();
    let grid = Grid::build(positions, &indices, radius);
    let hits: Vec<bool> = (0..n)
        .into_par_iter()
        .map(|i| inner.contains(i) || grid.any_within(positions, positions[i], radius))
        .collect();
    atom_mask(n, |i| hits[i])
}

fn byres(t: &Topology, inner: &FixedBitSet) -> FixedBitSet {
    let mut hit = vec![false; t.residue_count()];
    for a in inner.ones() {
        hit[t.residue_index[a] as usize] = true;
    }
    residue_mask(t, |i, _| hit[i])
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{AtomRow, TopologyBuilder};

    struct Atom {
        chain: (&'static str, &'static str),
        comp: &'static str,
        seq: (i32, i32),
        name: &'static str,
        element: &'static str,
        hetero: bool,
        b: f32,
    }

    const fn atom(
        chain: (&'static str, &'static str),
        comp: &'static str,
        seq: (i32, i32),
        name: &'static str,
        element: &'static str,
        hetero: bool,
        b: f32,
    ) -> Atom {
        Atom {
            chain,
            comp,
            seq,
            name,
            element,
            hetero,
            b,
        }
    }

    /// Chain A: ALA 1 (helix) + GLY 2 (strand); chain B: DA 5 (coil);
    /// chain L (auth A): HEM with resid -3; chain W (auth B): two waters.
    /// Atoms sit on the x axis 1.5 apart within a residue, far apart between
    /// residues, so `within` distances are easy to reason about.
    const ATOMS: [Atom; 18] = [
        atom(("A", "A"), "ALA", (1, 1), "N", "N", false, 10.0),
        atom(("A", "A"), "ALA", (1, 1), "CA", "C", false, 20.0),
        atom(("A", "A"), "ALA", (1, 1), "C", "C", false, 30.0),
        atom(("A", "A"), "ALA", (1, 1), "O", "O", false, 40.0),
        atom(("A", "A"), "ALA", (1, 1), "CB", "C", false, 50.0),
        atom(("A", "A"), "ALA", (1, 1), "HB1", "H", false, 5.0),
        atom(("A", "A"), "GLY", (2, 2), "N", "N", false, 15.0),
        atom(("A", "A"), "GLY", (2, 2), "CA", "C", false, 15.0),
        atom(("A", "A"), "GLY", (2, 2), "C", "C", false, 15.0),
        atom(("A", "A"), "GLY", (2, 2), "O", "O", false, 15.0),
        atom(("B", "B"), "DA", (5, 5), "P", "P", false, 15.0),
        atom(("B", "B"), "DA", (5, 5), "OP1", "O", false, 15.0),
        atom(("B", "B"), "DA", (5, 5), "C1'", "C", false, 15.0),
        atom(("B", "B"), "DA", (5, 5), "N9", "N", false, 15.0),
        atom(("L", "A"), "HEM", (1, -3), "FE", "Fe", true, 15.0),
        atom(("L", "A"), "HEM", (1, -3), "NA", "N", true, 15.0),
        atom(("W", "B"), "HOH", (1, 201), "O", "O", true, 15.0),
        atom(("W", "B"), "HOH", (2, 202), "O", "O", true, 15.0),
    ];

    fn fixture() -> (Topology, Vec<Vec3>) {
        let mut b = TopologyBuilder::new();
        let mut x = 0.0;
        let mut last_seq = None;
        for (i, a) in ATOMS.iter().enumerate() {
            if last_seq.replace((a.chain.0, a.seq.0)) != Some((a.chain.0, a.seq.0)) {
                x += 20.0;
            }
            let mut name = [b' '; 4];
            name[..a.name.len()].copy_from_slice(a.name.as_bytes());
            let mut position = Vec3::new(x, 0.0, 0.0);
            if a.name == "CB" {
                position = Vec3::new(x - 1.5 * 3.0, 1.5, 0.0); // 1.5 above CA
            } else if a.name == "HB1" {
                position = Vec3::new(x - 1.5 * 4.0, 2.5, 0.0); // 2.5 above CA
            }
            b.push(&AtomRow {
                element: Element::from_symbol(a.element.as_bytes()),
                name,
                serial: i as u32 + 1,
                alt_loc: if a.name == "CB" { b'A' } else { 0 },
                comp: a.comp,
                asym: a.chain.0,
                auth_asym: a.chain.1,
                seq_id: a.seq.0,
                auth_seq_id: a.seq.1,
                ins_code: 0,
                entity: 1,
                position,
                occupancy: if a.name == "CB" { 0.5 } else { 1.0 },
                b_factor: a.b,
                charge: 0,
                hetero: a.hetero,
            });
            x += 1.5;
        }
        let mut t = b.topology;
        t.residues[0].ss = SecondaryStructure::Helix;
        t.residues[1].ss = SecondaryStructure::Strand;
        t.residues[2].ss = SecondaryStructure::Coil;
        t.validate().unwrap();
        (t, b.positions)
    }

    fn ones(expr: &str) -> Vec<usize> {
        let (t, p) = fixture();
        select(&t, &p, expr)
            .unwrap_or_else(|e| panic!("{expr}: {e}"))
            .ones()
            .collect()
    }

    fn err(expr: &str) -> SelectError {
        parse(expr).expect_err(expr)
    }

    fn range(r: Range<usize>) -> Vec<usize> {
        r.collect()
    }

    #[test]
    fn classes() {
        assert_eq!(ones("all"), range(0..18));
        assert_eq!(ones("none"), vec![]);
        assert_eq!(ones("protein"), range(0..10));
        assert_eq!(ones("nucleic"), range(10..14));
        assert_eq!(ones("water"), vec![16, 17]);
        assert_eq!(ones("hetero"), vec![14, 15, 16, 17]);
        assert_eq!(ones("backbone"), vec![0, 1, 2, 3, 6, 7, 8, 9, 10, 11, 12]);
        assert_eq!(ones("sidechain"), vec![4, 5]);
        assert_eq!(ones("hydrogen"), vec![5]);
        assert_eq!(ones("helix"), range(0..6));
        assert_eq!(ones("strand"), range(6..10));
        assert_eq!(ones("coil"), range(10..14));
        assert_eq!(
            ones("PROTEIN"),
            range(0..10),
            "keywords are case-insensitive"
        );
    }

    /// One atom per residue (`comp`, `name`), same chain, for the naming-
    /// convention tests below: they only care about resname/atom-name
    /// classification, not real bonded geometry.
    fn one_atom_per_residue(entries: &[(&'static str, &'static str)]) -> Topology {
        let mut b = TopologyBuilder::new();
        for (i, &(comp, name)) in entries.iter().enumerate() {
            let mut n = [b' '; 4];
            n[..name.len()].copy_from_slice(name.as_bytes());
            b.push(&AtomRow {
                element: Element::UNKNOWN,
                name: n,
                serial: i as u32 + 1,
                alt_loc: 0,
                comp,
                asym: "A",
                auth_asym: "A",
                seq_id: i as i32 + 1,
                auth_seq_id: i as i32 + 1,
                ins_code: 0,
                entity: 1,
                position: Vec3::new(i as f32 * 20.0, 0.0, 0.0),
                occupancy: 1.0,
                b_factor: 0.0,
                charge: 0,
                hetero: true,
            });
        }
        (*b.finish().unwrap().topology).clone()
    }

    /// The residue names `expr` matches, out of `entries`.
    fn matches(entries: &[(&'static str, &'static str)], expr: &str) -> Vec<&'static str> {
        let t = one_atom_per_residue(entries);
        select(&t, &[], expr)
            .unwrap_or_else(|e| panic!("{expr}: {e}"))
            .ones()
            .map(|i| entries[i].0)
            .collect()
    }

    #[test]
    fn ion_keyword_covers_pdb_charmm_and_amber_names() {
        let entries = [
            ("NA", "NA"),
            ("SOD", "SOD"),
            ("CA", "CA"), // calcium, not an alpha carbon: no other atom in its residue
            ("ZN2+", "ZN2+"),
            ("ALA", "N"),
            ("HOH", "O"),
        ];
        assert_eq!(matches(&entries, "ion"), vec!["NA", "SOD", "CA", "ZN2+"]);
    }

    #[test]
    fn glycan_keyword_delegates_to_the_snfg_table() {
        let entries = [("NAG", "C1"), ("BMA", "C1"), ("ALA", "N")];
        assert_eq!(matches(&entries, "glycan"), vec!["NAG", "BMA"]);
    }

    #[test]
    fn molecule_class_keywords_partition_by_residue_class() {
        let entries = [
            ("ALA", "N"),
            ("DA", "P"),
            ("POPC", "P"),
            ("NAG", "C1"),
            ("TIP3", "OH2"),
            ("SOD", "SOD"),
        ];
        let names = |expr| matches(&entries, expr);
        assert_eq!(names("lipid"), vec!["POPC"]);
        assert_eq!(names("solvent"), vec!["TIP3", "SOD"]);
        assert_eq!(names("polymer"), vec!["ALA", "DA"]);
        assert_eq!(names("polymer or lipid or glycan or solvent").len(), 6);
    }

    #[test]
    fn roles_overlap_the_classes() {
        let entries = [
            ("POPC", "P"), // one lipid: too few for a membrane, so a ligand
            ("NAG", "C1"), // a lone glycan is not attached to anything
            ("GOL", "C1"), // unknown element: class other, no role
            ("ALA", "N"),
        ];
        let names = |expr| matches(&entries, expr);
        assert_eq!(names("ligand"), vec!["POPC", "NAG"]);
        assert_eq!(names("lipid"), vec!["POPC"]);
        assert_eq!(names("membrane"), Vec::<&str>::new());
        assert_eq!(names("ligand and glycan"), vec!["NAG"]);
    }

    #[test]
    fn protein_keyword_covers_charmm_and_amber_protonation_states() {
        let entries = [
            ("HSD", "N"), // CHARMM neutral histidine
            ("HID", "N"), // Amber neutral histidine
            ("CYX", "N"), // Amber disulfide cystine
            ("ASH", "N"), // Amber protonated aspartate
            ("NLN", "N"), // GLYCAM N-linked asparagine
            ("HOH", "O"),
        ];
        assert_eq!(
            matches(&entries, "protein"),
            vec!["HSD", "HID", "CYX", "ASH", "NLN"]
        );
    }

    #[test]
    fn protein_keyword_covers_amber_terminal_residue_templates() {
        let entries = [
            ("NALA", "N"), // Amber's own N-terminal alanine template
            ("CHIS", "N"), // and a C-terminal histidine
            ("NAG", "C1"), // NOT protein: a glycan, not an N-terminal "AG"
            ("NME", "N"),  // NOT protein: an amine cap, "ME" is not a residue code
        ];
        assert_eq!(matches(&entries, "protein"), vec!["NALA", "CHIS"]);
    }

    #[test]
    fn water_keyword_covers_other_force_fields_models() {
        let entries = [
            ("SPC", "OW"),
            ("TIP4", "OW"),
            ("T3P", "O"), // T3P: TIP3P as some MD packages name it
            ("ALA", "N"),
        ];
        assert_eq!(matches(&entries, "water"), vec!["SPC", "TIP4", "T3P"]);
    }

    #[test]
    fn nucleic_keyword_covers_amber_terminal_residue_templates() {
        let entries = [("DA5", "P"), ("U3", "P"), ("ALA", "N")];
        assert_eq!(matches(&entries, "nucleic"), vec!["DA5", "U3"]);
    }

    #[test]
    fn backbone_keyword_covers_charmm_terminal_oxygen_names() {
        let entries = [("ALA", "OT1"), ("ALA", "OT2"), ("ALA", "CB")];
        assert_eq!(matches(&entries, "backbone"), vec!["ALA", "ALA"]);
        assert_eq!(matches(&entries, "sidechain"), vec!["ALA"]);
    }

    #[test]
    fn chain_matches_label_or_auth_id_exactly() {
        assert_eq!(ones("chain A"), vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 14, 15]);
        assert_eq!(ones("chain B"), vec![10, 11, 12, 13, 16, 17]);
        assert_eq!(ones("chain L"), vec![14, 15]);
        assert_eq!(ones("chain W L"), vec![14, 15, 16, 17]);
        assert_eq!(ones("chain a"), vec![], "chain ids are case-sensitive");
    }

    #[test]
    fn segname_needs_a_value_and_matches_nothing_without_segids() {
        assert_eq!(ones("segname PROA"), vec![]);
        assert_eq!(ones("segid PROA"), vec![]);
        assert!(parse("segname").is_err());
    }

    #[test]
    fn names_and_elements() {
        assert_eq!(ones("resname hem"), vec![14, 15]);
        assert_eq!(ones("resname ALA GLY"), range(0..10));
        assert_eq!(ones("name ca"), vec![1, 7]);
        assert_eq!(ones("name C1' n9"), vec![12, 13]);
        assert_eq!(ones("element fe"), vec![14]);
        assert_eq!(ones("element C"), vec![1, 2, 4, 7, 8, 12]);
        assert_eq!(ones("element O H"), vec![3, 5, 9, 11, 16, 17]);
        assert_eq!(ones("altloc A"), vec![4]);
        assert_eq!(ones("altloc B"), vec![]);
    }

    #[test]
    fn ranges_including_negative_numbers() {
        assert_eq!(ones("resid 1-2"), range(0..10));
        assert_eq!(ones("resid 1:2"), range(0..10));
        assert_eq!(ones("resid -3"), vec![14, 15]);
        assert_eq!(ones("resid -3--1"), vec![14, 15]);
        assert_eq!(ones("resid -3-1"), vec![0, 1, 2, 3, 4, 5, 14, 15]);
        assert_eq!(ones("resid 200:202 5"), vec![10, 11, 12, 13, 16, 17]);
        assert_eq!(ones("seqid 1"), vec![0, 1, 2, 3, 4, 5, 14, 15, 16]);
        assert_eq!(ones("serial 1-3 18"), vec![0, 1, 2, 17]);
        assert_eq!(ones("index 0-2"), vec![0, 1, 2]);
        assert_eq!(ones("index 17-100"), vec![17]);
        assert_eq!(ones("index -5-0"), vec![0]);
    }

    #[test]
    fn comparisons() {
        assert_eq!(ones("bfactor < 10"), vec![5]);
        assert_eq!(ones("bfactor <= 10"), vec![0, 5]);
        assert_eq!(ones("bfactor >= 40"), vec![3, 4]);
        assert_eq!(ones("bfactor > 40"), vec![4]);
        assert_eq!(ones("bfactor == 20"), vec![1]);
        assert_eq!(ones("bfactor != 15"), vec![0, 1, 2, 3, 4, 5]);
        assert_eq!(ones("occupancy<1"), vec![4], "operators need no spaces");
        assert_eq!(ones("occupancy = 0.5"), vec![4]);
    }

    #[test]
    fn precedence_and_grouping() {
        // and binds tighter than or.
        assert_eq!(ones("water or protein and name CA"), vec![1, 7, 16, 17]);
        assert_eq!(ones("(water or protein) and name CA"), vec![1, 7]);
        // not binds tighter than and.
        assert_eq!(ones("not protein and water"), vec![16, 17]);
        assert_eq!(ones("not (protein or water)"), vec![10, 11, 12, 13, 14, 15]);
        assert_eq!(ones("not not protein"), range(0..10));
        assert_eq!(ones("protein and not backbone"), vec![4, 5]);
        let expr = parse("all or none and all").unwrap();
        assert!(matches!(expr, Expr::Or(_, ref rhs) if matches!(**rhs, Expr::And(..))));
    }

    #[test]
    fn within_and_byres() {
        // CA atoms at 1.5 from N, C and CB; O is 3.0 away and HB1 2.5.
        assert_eq!(ones("within 2 of name CA"), vec![0, 1, 2, 4, 6, 7, 8]);
        assert_eq!(
            ones("within 1.5 of name CA"),
            vec![0, 1, 2, 4, 6, 7, 8],
            "inclusive"
        );
        assert_eq!(ones("within 1.4 of name CA"), vec![1, 7]);
        assert_eq!(ones("within 0 of name CA"), vec![1, 7]);
        assert_eq!(ones("within 2 of none"), vec![]);
        // The operand binds tighter than `and`.
        assert_eq!(ones("within 2 of name CA and element N"), vec![0, 6]);
        assert_eq!(ones("byres name CA"), range(0..10));
        assert_eq!(ones("byres element fe"), vec![14, 15]);
        assert_eq!(
            ones("byres within 2 of name CB and not name CA"),
            vec![0, 2, 3, 4, 5]
        );
    }

    #[test]
    fn error_messages_and_spans() {
        let e = err("");
        assert_eq!((e.message.as_str(), e.span), ("empty selection", 0..0));
        assert_eq!(err("   ").span, 0..3);

        let e = err("chain A and foo");
        assert_eq!(e.message, "unknown keyword `foo`");
        assert_eq!(e.span, 12..15);

        let e = err("within of protein");
        assert_eq!(e.message, "expected a number after `within`, found `of`");
        assert_eq!(e.span, 7..9);

        let e = err("within 3 protein");
        assert_eq!(e.message, "expected `of` after `3`, found `protein`");
        assert_eq!(e.span, 9..16);

        assert_eq!(err("within -1 of all").span, 7..9);

        let e = err("(protein");
        assert_eq!((e.message.as_str(), e.span), ("unclosed `(`", 0..1));

        let e = err("protein)");
        assert!(e.message.starts_with("unexpected `)`"), "{}", e.message);
        assert_eq!(e.span, 7..8);

        let e = err("protein chain A");
        assert!(e.message.starts_with("unexpected `chain`"), "{}", e.message);
        assert_eq!(e.span, 8..15);

        let e = err("resid 1-");
        assert!(e.message.starts_with("bad range `1-`"), "{}", e.message);
        assert_eq!(e.span, 6..8);
        assert_eq!(err("resid 5-1").span, 6..9);
        assert_eq!(err("resid x").span, 6..7);

        let e = err("bfactor 10");
        assert!(
            e.message.starts_with("expected a comparison"),
            "{}",
            e.message
        );
        assert_eq!(e.span, 8..10);
        assert_eq!(err("bfactor <").span, 8..9);

        let e = err("chain");
        assert_eq!(e.message, "expected one or more values after `chain`");
        assert_eq!(e.span, 0..5);
        assert_eq!(err("chain and protein").span, 6..9);

        let e = err("element Xx");
        assert_eq!(e.message, "unknown element symbol `Xx`");
        assert_eq!(e.span, 8..10);

        assert_eq!(err("altloc AB").span, 7..9);

        let e = err("not");
        assert_eq!(e.message, "expected a selection after `not`");
        assert_eq!(e.span, 0..3);
        assert_eq!(err("protein and").span, 8..11);
        assert_eq!(err("and protein").span, 0..3);
        assert_eq!(err("protein ! water").span, 8..9);
    }

    #[test]
    fn evaluate_without_positions_is_fine_unless_within_is_used() {
        let (t, _) = fixture();
        let m = parse("protein and name CA").unwrap().evaluate(&t, &[]);
        assert_eq!(m.ones().collect::<Vec<_>>(), vec![1, 7]);
    }
}
