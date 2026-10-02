//! Reading the plain-text reply of the scheme authors' numbering program
//! served at bioinf.org.uk/abs/abnum (Abhinandan & Martin 2008, Mol Immunol
//! 45:3832). The application sends one chain per request; this module
//! only reads what comes back. Layout and limits: docs/ANTIBODY.md.

use super::external::{parse_label, NumberedDomain};
use super::{find_domains, ChainType, Label, Scheme};

/// The numberings the program offers.
pub const SCHEMES: [Scheme; 3] = [Scheme::Kabat, Scheme::Chothia, Scheme::Martin];

/// The `scheme` request field for `scheme`; `None` for IMGT and AHo.
pub fn flag(scheme: Scheme) -> Option<&'static str> {
    match scheme {
        Scheme::Kabat => Some("-k"),
        Scheme::Chothia => Some("-c"),
        Scheme::Martin => Some("-m"),
        Scheme::Imgt | Scheme::Aho => None,
    }
}

/// `H52A P`: chain letter, number, optional insertion letter, then the
/// residue.
fn parse_row(line: &str) -> Result<(char, Label, char), String> {
    let bad = || format!("unexpected line `{line}`");
    let mut fields = line.split_whitespace();
    let (Some(label), Some(residue), None) = (fields.next(), fields.next(), fields.next()) else {
        return Err(bad());
    };
    let mut letters = label.chars();
    let chain = letters.next().ok_or_else(bad)?;
    let rest = letters.as_str();
    let digits = rest
        .find(|c: char| !c.is_ascii_digit())
        .unwrap_or(rest.len());
    let insertion = rest[digits..].chars().next().map(|c| c.to_string());
    let label = parse_label(&rest[..digits], insertion.as_deref())?;
    let residue = residue.chars().next().ok_or_else(bad)?;
    Ok((chain, label, residue))
}

/// Light chains come back as `L`; kappa or lambda is decided by the
/// native profiles, kappa when they find nothing.
fn chain_type(letter: char, residues: &str) -> Result<ChainType, String> {
    match letter {
        'H' => Ok(ChainType::Heavy),
        'L' => Ok(find_domains(residues)
            .first()
            .map(|d| d.chain)
            .filter(|c| matches!(c, ChainType::Kappa | ChainType::Lambda))
            .unwrap_or(ChainType::Kappa)),
        other => Err(format!("unknown chain letter `{other}`")),
    }
}

/// The numbered domain in the reply to `query`; `None` when the program
/// found no variable domain (its reply is then a `# Error` comment). It
/// leaves residues outside the domain unnumbered, so the range is located
/// in `query`.
pub fn parse_abnum(query: &str, text: &str) -> Result<Option<NumberedDomain>, String> {
    let mut letter = None;
    let mut labels = Vec::new();
    let mut residues = String::new();
    let rows = text
        .lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty());
    for line in rows {
        let (chain, label, residue) = parse_row(line)?;
        if *letter.get_or_insert(chain) != chain {
            return Err("reply mixes chain letters".into());
        }
        labels.push(label);
        residues.push(residue);
    }
    let Some(letter) = letter else {
        return Ok(None);
    };
    let start = query
        .find(&residues)
        .ok_or("numbered residues are not in the sequence sent")?;
    Ok(Some(NumberedDomain {
        name: String::new(),
        chain: chain_type(letter, &residues)?,
        range: start..start + residues.len(),
        labels,
    }))
}
