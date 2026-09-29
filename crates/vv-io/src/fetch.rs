//! Fetch entries from the RCSB PDB by ID into a per-user cache: the
//! asymmetric unit (what the file deposits) or a biological assembly
//! (the functional oligomer RCSB builds from it, the PDB's
//! "biological unit" download). Density maps (EMDB, 2Fo-Fc) are not
//! fetched yet: nothing can render a volume.
//!
//! Behind the `fetch` feature (on by default) so a build that must not
//! link a TLS stack can leave it out.

use std::path::{Path, PathBuf};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Assembly {
    /// The deposited coordinates.
    AsymmetricUnit,
    /// Biological assembly `n` (1-based, as RCSB numbers them).
    Biological(u32),
}

#[derive(Debug, thiserror::Error)]
pub enum FetchError {
    #[error("`{0}` is not a PDB ID (4 characters like 4HHB, or pdb_0000XXXX)")]
    BadId(String),
    #[error("RCSB has no {what} for {id} (HTTP {status})")]
    NotFound {
        id: String,
        what: String,
        status: u16,
    },
    #[error("download of {url} failed: {source}")]
    Http { url: String, source: ureq::Error },
    #[error("cannot write {path}: {source}")]
    Io {
        path: PathBuf,
        source: std::io::Error,
    },
}

/// Per-user cache directory: `%LOCALAPPDATA%\vizviz\cache` on Windows,
/// `~/Library/Caches/vizviz` on macOS, `$XDG_CACHE_HOME/vizviz` (or
/// `~/.cache/vizviz`) elsewhere. `None` when no home can be found.
pub fn cache_dir() -> Option<PathBuf> {
    let var = |name: &str| std::env::var_os(name).filter(|v| !v.is_empty());
    let base = if cfg!(windows) {
        var("LOCALAPPDATA").map(PathBuf::from)
    } else if cfg!(target_os = "macos") {
        var("HOME").map(|h| PathBuf::from(h).join("Library/Caches"))
    } else {
        var("XDG_CACHE_HOME")
            .map(PathBuf::from)
            .or_else(|| var("HOME").map(|h| PathBuf::from(h).join(".cache")))
    }?;
    Some(
        base.join("vizviz")
            .join(if cfg!(windows) { "cache" } else { "" }),
    )
}

/// Upper-cased, validated ID: four alphanumerics whose first character is
/// a digit, or an extended `pdb_` + 8 characters.
pub fn normalize_id(id: &str) -> Result<String, FetchError> {
    let id = id.trim().to_ascii_uppercase();
    let classic = id.len() == 4
        && id.chars().all(|c| c.is_ascii_alphanumeric())
        && id.chars().next().is_some_and(|c| c.is_ascii_digit());
    let extended = id.len() == 12
        && id.starts_with("PDB_")
        && id[4..].chars().all(|c| c.is_ascii_alphanumeric());
    if classic || extended {
        Ok(id)
    } else {
        Err(FetchError::BadId(id))
    }
}

/// File name in the cache and URL at RCSB for one request.
pub fn locate(id: &str, assembly: Assembly) -> (String, String) {
    let file = match assembly {
        Assembly::AsymmetricUnit => format!("{id}.cif.gz"),
        Assembly::Biological(n) => format!("{id}-assembly{n}.cif.gz"),
    };
    (
        file.clone(),
        format!("https://files.rcsb.org/download/{file}"),
    )
}

/// File name in the cache and URL at RCSB's model server for the
/// asymmetric unit as BinaryCIF.
pub fn locate_bcif(id: &str) -> (String, String) {
    let file = format!("{id}.bcif");
    (file.clone(), format!("https://models.rcsb.org/{file}"))
}

/// The cached path for `id`, downloading it into `cache` first if it is
/// not there yet. Never re-downloads: delete the file to refresh.
pub fn fetch(id: &str, assembly: Assembly, cache: &Path) -> Result<PathBuf, FetchError> {
    let id = normalize_id(id)?;
    let (file, url) = locate(&id, assembly);
    let what = match assembly {
        Assembly::AsymmetricUnit => "entry".to_string(),
        Assembly::Biological(n) => format!("assembly {n}"),
    };
    download(&id, &what, &url, &cache.join(file))
}

/// [`fetch`] for the asymmetric unit as BinaryCIF (smaller and faster
/// to read than the text file).
pub fn fetch_bcif(id: &str, cache: &Path) -> Result<PathBuf, FetchError> {
    let id = normalize_id(id)?;
    let (file, url) = locate_bcif(&id);
    download(&id, "BinaryCIF entry", &url, &cache.join(file))
}

fn download(id: &str, what: &str, url: &str, path: &Path) -> Result<PathBuf, FetchError> {
    if path.exists() {
        return Ok(path.to_path_buf());
    }
    let response = match ureq::get(url).call() {
        Ok(r) => r,
        Err(ureq::Error::StatusCode(status)) => {
            return Err(FetchError::NotFound {
                id: id.to_string(),
                what: what.to_string(),
                status,
            })
        }
        Err(source) => {
            return Err(FetchError::Http {
                url: url.to_string(),
                source,
            })
        }
    };
    let mut bytes = Vec::new();
    std::io::Read::read_to_end(&mut response.into_body().into_reader(), &mut bytes).map_err(
        |source| FetchError::Http {
            url: url.to_string(),
            source: source.into(),
        },
    )?;
    let io = |source| FetchError::Io {
        path: path.to_path_buf(),
        source,
    };
    if let Some(dir) = path.parent() {
        std::fs::create_dir_all(dir).map_err(io)?;
    }
    // Write to a temp name and rename, so a killed download never leaves
    // a truncated file the cache would trust next time.
    let partial = path.with_extension("part");
    std::fs::write(&partial, &bytes).map_err(io)?;
    std::fs::rename(&partial, path).map_err(io)?;
    Ok(path.to_path_buf())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bcif_names_and_urls() {
        let (file, url) = locate_bcif("4HHB");
        assert_eq!(file, "4HHB.bcif");
        assert_eq!(url, "https://models.rcsb.org/4HHB.bcif");
    }

    #[test]
    fn ids_are_validated_and_upper_cased() {
        assert_eq!(normalize_id(" 4hhb ").unwrap(), "4HHB");
        assert_eq!(normalize_id("pdb_00004hhb").unwrap(), "PDB_00004HHB");
        for bad in ["", "hhb", "4hhbx", "ahhb", "4hh-", "pdb_1"] {
            assert!(
                matches!(normalize_id(bad), Err(FetchError::BadId(_))),
                "{bad}"
            );
        }
    }

    #[test]
    fn cache_names_and_urls() {
        let (file, url) = locate("4HHB", Assembly::AsymmetricUnit);
        assert_eq!(file, "4HHB.cif.gz");
        assert_eq!(url, "https://files.rcsb.org/download/4HHB.cif.gz");
        let (file, url) = locate("4HHB", Assembly::Biological(2));
        assert_eq!(file, "4HHB-assembly2.cif.gz");
        assert!(url.ends_with("/4HHB-assembly2.cif.gz"));
    }

    #[test]
    fn a_cached_file_is_returned_without_a_network() {
        let dir = std::env::temp_dir().join(format!("vizviz_fetch_test_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(dir.join("1CRN.cif.gz"), b"placeholder").unwrap();
        let path = fetch("1crn", Assembly::AsymmetricUnit, &dir).unwrap();
        assert_eq!(path, dir.join("1CRN.cif.gz"));
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// Needs the network; run with `--ignored`.
    #[test]
    #[ignore]
    fn downloads_crambin_and_its_assembly() {
        let dir = std::env::temp_dir().join(format!("vizviz_fetch_net_{}", std::process::id()));
        let asym = fetch("1CRN", Assembly::AsymmetricUnit, &dir).unwrap();
        let s = crate::load(&asym).unwrap();
        assert_eq!(s.atom_count(), 327);
        let bcif = fetch_bcif("1CRN", &dir).unwrap();
        assert_eq!(crate::load(&bcif).unwrap().atom_count(), 327);
        let bio = fetch("1CRN", Assembly::Biological(1), &dir).unwrap();
        assert!(crate::load(&bio).unwrap().atom_count() >= 327);
        assert!(matches!(
            fetch("0ZZZ", Assembly::AsymmetricUnit, &dir),
            Err(FetchError::NotFound { status: 404, .. })
        ));
        let _ = std::fs::remove_dir_all(&dir);
    }
}
