//! Raw JSON from two public services, cached on disk like [`crate::fetch`]:
//! the PDBe SIFTS mapping of a PDB entry's chains to UniProt accessions,
//! and a UniProt record's feature table. Parsing is the caller's; a cached
//! file is returned without touching the network, and any failure is an
//! error the caller can ignore to stay offline-safe.

use std::path::{Path, PathBuf};
use std::time::Duration;

use crate::fetch::{cache_dir, normalize_id, FetchError};

const TIMEOUT: Duration = Duration::from_secs(20);

/// The UniProt fields the strip draws as features.
const FEATURE_FIELDS: &str = "accession,id,ft_domain,ft_region,ft_repeat,ft_motif,ft_act_site,\
ft_binding,ft_site,ft_mod_res,ft_carbohyd,ft_disulfid,ft_lipid,ft_crosslnk,ft_signal,\
ft_transit,ft_propep,ft_transmem,ft_intramem,ft_topo_dom,ft_variant";

fn sifts_url(pdb_id: &str) -> String {
    format!(
        "https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{}",
        pdb_id.to_ascii_lowercase()
    )
}

fn features_url(accession: &str) -> String {
    format!("https://rest.uniprot.org/uniprotkb/{accession}.json?fields={FEATURE_FIELDS}")
}

/// A UniProt accession: 6 or 10 upper-case alphanumerics. Rejects anything
/// that would alter the request path.
fn valid_accession(accession: &str) -> bool {
    matches!(accession.len(), 6 | 10)
        && accession
            .bytes()
            .all(|b| b.is_ascii_uppercase() || b.is_ascii_digit())
}

/// The SIFTS mapping of entry `pdb_id`, as JSON text.
pub fn fetch_sifts(pdb_id: &str) -> Result<String, FetchError> {
    let id = normalize_id(pdb_id)?;
    cached_text(&sifts_url(&id), &file_in_cache(&format!("{id}.sifts.json")))
}

/// The feature table of UniProt entry `accession`, as JSON text.
pub fn fetch_features(accession: &str) -> Result<String, FetchError> {
    if !valid_accession(accession) {
        return Err(FetchError::BadId(accession.to_string()));
    }
    cached_text(
        &features_url(accession),
        &file_in_cache(&format!("{accession}.uniprot.json")),
    )
}

fn file_in_cache(name: &str) -> PathBuf {
    cache_dir()
        .unwrap_or_else(std::env::temp_dir)
        .join("uniprot")
        .join(name)
}

fn cached_text(url: &str, path: &Path) -> Result<String, FetchError> {
    if let Ok(text) = std::fs::read_to_string(path) {
        return Ok(text);
    }
    let text = get(url)?;
    let io = |source| FetchError::Io {
        path: path.to_path_buf(),
        source,
    };
    if let Some(dir) = path.parent() {
        std::fs::create_dir_all(dir).map_err(io)?;
    }
    let partial = path.with_extension("part");
    std::fs::write(&partial, &text).map_err(io)?;
    std::fs::rename(&partial, path).map_err(io)?;
    Ok(text)
}

fn get(url: &str) -> Result<String, FetchError> {
    let agent: ureq::Agent = ureq::Agent::config_builder()
        .timeout_global(Some(TIMEOUT))
        .build()
        .into();
    let response = match agent.get(url).header("Accept", "application/json").call() {
        Ok(r) => r,
        Err(ureq::Error::StatusCode(status)) => {
            return Err(FetchError::NotFound {
                id: url.to_string(),
                what: "record".to_string(),
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
    response
        .into_body()
        .read_to_string()
        .map_err(|source| FetchError::Http {
            url: url.to_string(),
            source,
        })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn urls_name_the_service_and_the_entry() {
        assert!(sifts_url("4HHB").ends_with("/mappings/uniprot/4hhb"));
        let url = features_url("P69905");
        assert!(url.starts_with("https://rest.uniprot.org/uniprotkb/P69905.json?fields="));
        assert!(url.contains("ft_binding"));
    }

    #[test]
    fn accessions_are_checked_before_they_reach_a_url() {
        assert!(valid_accession("P69905"));
        assert!(valid_accession("A0A024R161"));
        for bad in ["", "p69905", "P6990", "../etc/x", "P69905?x=1"] {
            assert!(!valid_accession(bad), "{bad}");
        }
    }

    #[test]
    fn a_cached_file_is_returned_without_a_network() {
        let dir = std::env::temp_dir().join(format!("vizviz_up_test_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("X.json");
        std::fs::write(&path, "{\"cached\":true}").unwrap();
        assert_eq!(
            cached_text("http://invalid.invalid/", &path).unwrap(),
            "{\"cached\":true}"
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn an_unreachable_service_is_an_error_not_a_panic() {
        let dir = std::env::temp_dir().join(format!("vizviz_up_off_{}", std::process::id()));
        let path = dir.join("Y.json");
        assert!(cached_text("http://127.0.0.1:9/none", &path).is_err());
        assert!(!path.exists());
    }
}
