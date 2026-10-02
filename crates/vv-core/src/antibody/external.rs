//! Domains numbered by an external program instead of the native profiles:
//! the same annotation shape as [`Domain`](super::Domain), built from
//! labels the program reported and regions from the usual CDR rules.
//!
//! The program is ANARCI (Dunbar & Deane 2016, Bioinformatics 32:298), run
//! by the application; this module only holds and reads its results.

use std::ops::Range;

use super::cdr::CdrDefinition;
use super::numbering::{Label, Scheme};
use super::topology::CdrResidue;
use super::{Annotation, ChainType};

/// A variable domain with labels for one or more schemes.
#[derive(Clone, Debug)]
pub struct ExternalDomain {
    pub chain: ChainType,
    /// Residue range `start..end` in the sequence the program was given.
    pub start: usize,
    pub end: usize,
    labels: Vec<(Scheme, Vec<Label>)>,
}

impl ExternalDomain {
    pub fn new(chain: ChainType, range: Range<usize>) -> Self {
        Self {
            chain,
            start: range.start,
            end: range.end,
            labels: Vec::new(),
        }
    }

    /// Adds the labels of each domain residue in `scheme`; `None` if the
    /// count does not match the domain (see [`Self::set_labels`]).
    pub fn with_labels(mut self, scheme: Scheme, labels: Vec<Label>) -> Option<Self> {
        self.set_labels(scheme, labels).then_some(self)
    }

    /// Replaces the labels in `scheme`; false, changing nothing, if the
    /// count does not match the domain.
    pub fn set_labels(&mut self, scheme: Scheme, labels: Vec<Label>) -> bool {
        let fits = labels.len() == self.end - self.start;
        if fits {
            self.labels.retain(|(s, _)| *s != scheme);
            self.labels.push((scheme, labels));
        }
        fits
    }

    /// The labels for `scheme`: receptor domains only have IMGT, as in the
    /// native numbering.
    fn labels_in(&self, scheme: Scheme) -> Option<&[Label]> {
        let wanted = match self.chain.is_antibody() {
            true => scheme,
            false => Scheme::Imgt,
        };
        let found = self.labels.iter().find(|(s, _)| *s == wanted);
        found.map(|(_, l)| l.as_slice())
    }

    /// Whether [`Self::annotate`] can answer for `scheme` under `definition`.
    pub fn covers(&self, scheme: Scheme, definition: CdrDefinition) -> bool {
        self.labels_in(scheme).is_some() && self.labels_in(definition.native_scheme()).is_some()
    }

    /// Label in `scheme` and region under `definition` for each residue, as
    /// [`Domain::annotate`](super::Domain::annotate); `None` when the
    /// program was not run for one of the two schemes.
    pub fn annotate(&self, scheme: Scheme, definition: CdrDefinition) -> Option<Vec<Annotation>> {
        let shown = self.labels_in(scheme)?;
        let native = self.labels_in(definition.native_scheme())?;
        let regions = definition.regions(self.chain, native);
        Some(
            (self.start..self.end)
                .zip(shown)
                .zip(regions)
                .map(|((index, &label), region)| Annotation {
                    index,
                    label,
                    region,
                })
                .collect(),
        )
    }
}

/// Every CDR residue under `definition`, for chains given as the first
/// residue index of each with its domains.
pub fn cdr_residues_external<'a>(
    chains: impl IntoIterator<Item = (u32, &'a [ExternalDomain])>,
    definition: CdrDefinition,
) -> Vec<CdrResidue> {
    let scheme = definition.native_scheme();
    let mut out = Vec::new();
    for (first, domains) in chains {
        for domain in domains {
            let notes = domain.annotate(scheme, definition).unwrap_or_default();
            out.extend(notes.iter().filter_map(|a| {
                Some(CdrResidue {
                    residue: first + a.index as u32,
                    chain: domain.chain,
                    cdr: a.region.cdr()?,
                })
            }));
        }
    }
    out
}

/// One numbered domain read from ANARCI's text output.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct NumberedDomain {
    /// The FASTA name of the sequence.
    pub name: String,
    pub chain: ChainType,
    pub range: Range<usize>,
    /// One label per residue of `range`, gaps dropped.
    pub labels: Vec<Label>,
}

fn chain_of(letter: &str) -> Option<ChainType> {
    Some(match letter {
        "H" => ChainType::Heavy,
        "K" => ChainType::Kappa,
        "L" => ChainType::Lambda,
        "A" => ChainType::TcrAlpha,
        "B" => ChainType::TcrBeta,
        _ => return None,
    })
}

fn parse_label(number: &str, insertion: Option<&str>) -> Result<Label, String> {
    let number: u16 = number
        .parse()
        .map_err(|_| format!("bad residue number `{number}`"))?;
    match insertion.and_then(|s| s.chars().next()) {
        Some(c) => Ok(Label::with_insertion(number, c)),
        None => Ok(Label::new(number)),
    }
}

/// `|species|type|evalue|score|start|end|` row: chain letter and the
/// inclusive residue range.
fn parse_hit(row: &str) -> Result<(String, Range<usize>), String> {
    let f: Vec<&str> = row.trim_start_matches('#').split('|').collect();
    let num = |i: usize| {
        f.get(i)
            .and_then(|s| s.trim().parse::<usize>().ok())
            .ok_or_else(|| format!("bad hit row `{row}`"))
    };
    let chain = f.get(2).map(|s| s.trim().to_string()).unwrap_or_default();
    Ok((chain, num(5)?..num(6)? + 1))
}

/// Reads the residue lines of one domain: `H 52 Y`, `H 52 A P` (insertion
/// letter), or a `-` where the scheme has a gap. Longer lines are notices
/// ANARCI prints to the same stream; the residue count catches any damage.
fn parse_residues(lines: &[&str]) -> Result<Vec<Label>, String> {
    let mut labels = Vec::new();
    for line in lines {
        let t: Vec<&str> = line.split_whitespace().collect();
        let (insertion, residue) = match t.len() {
            3 => (None, t[2]),
            4 => (Some(t[2]), t[3]),
            5.. => continue,
            _ => return Err(format!("unexpected line `{line}`")),
        };
        if residue != "-" {
            labels.push(parse_label(t[1], insertion)?);
        }
    }
    Ok(labels)
}

struct Block<'a> {
    name: String,
    hit: Option<(String, Range<usize>)>,
    residues: Vec<&'a str>,
}

/// Comment lines ANARCI writes around every domain; any other `# text`
/// line is a sequence's name.
fn is_header(text: &str) -> bool {
    ["ANARCI numbered", "Domain ", "Most significant", "Scheme ="]
        .iter()
        .any(|h| text.starts_with(h))
}

fn blocks(text: &str) -> Result<Vec<Block<'_>>, String> {
    let mut out = Vec::new();
    let mut name = String::new();
    let mut current: Option<Block> = None;
    let mut hit_next = false;
    for line in text.lines() {
        if line.starts_with("//") {
            out.extend(current.take());
        } else if line.starts_with("#|species") {
            hit_next = true;
        } else if line.starts_with("#|") && hit_next {
            hit_next = false;
            if let Some(b) = current.as_mut() {
                b.hit = Some(parse_hit(line)?);
            }
        } else if line.starts_with("# Domain ") {
            out.extend(current.take());
            current = Some(Block {
                name: name.clone(),
                hit: None,
                residues: Vec::new(),
            });
        } else if let Some(n) = line.strip_prefix("# ").filter(|n| !is_header(n)) {
            name = n.trim().to_string();
        } else if !line.starts_with('#') && !line.trim().is_empty() {
            if let Some(b) = current.as_mut() {
                b.residues.push(line);
            }
        }
    }
    out.extend(current);
    Ok(out)
}

/// Domains in ANARCI's text output (its default stdout format). Domains of
/// chain types with no native equivalent (gamma and delta receptors) are
/// left out; a malformed record is an error naming the problem.
pub fn parse_anarci(text: &str) -> Result<Vec<NumberedDomain>, String> {
    let mut out = Vec::new();
    for block in blocks(text)? {
        let (letter, range) = block
            .hit
            .ok_or_else(|| format!("no hit row for `{}`", block.name))?;
        let Some(chain) = chain_of(&letter) else {
            continue;
        };
        let labels = parse_residues(&block.residues)?;
        if labels.len() != range.len() {
            return Err(format!(
                "`{}`: {} numbered residues for a {}-residue domain",
                block.name,
                labels.len(),
                range.len()
            ));
        }
        out.push(NumberedDomain {
            name: block.name,
            chain,
            range,
            labels,
        });
    }
    Ok(out)
}
