//! The public Abnum web service of the scheme authors' group (UCL,
//! http://www.bioinf.org.uk/abs/abnum/): one chain and one scheme per GET
//! request, the plain-text reply cached on disk so a repeat is offline.
//! The site serves plain HTTP only and says its data are not encrypted;
//! callers must have the user's explicit consent before sending anything.
//! Reading the reply: `vv_core::antibody::abnum`.

use std::path::{Path, PathBuf};
use std::sync::{Mutex, PoisonError};
use std::time::{Duration, Instant};

use crate::fetch::{cache_dir, FetchError};

pub const URL: &str = "http://www.bioinf.org.uk/abs/abnum/abnum.cgi";

const TIMEOUT: Duration = Duration::from_secs(20);
/// Least time between two requests to the server.
const SPACING: Duration = Duration::from_millis(500);
const USER_AGENT: &str = concat!("vizviz/", env!("CARGO_PKG_VERSION"), " (antibody numbering)");

/// When the last request finished; held across a request so concurrent
/// callers queue behind it.
static LAST_REQUEST: Mutex<Option<Instant>> = Mutex::new(None);

/// Abnum's reply for `sequence` under the scheme field `flag` (`-k`, `-c`
/// or `-m`), from the cache when present.
pub fn number(sequence: &str, flag: &str) -> Result<String, FetchError> {
    let cache = cache_dir().unwrap_or_else(std::env::temp_dir);
    number_at(URL, &cache, sequence, flag)
}

/// [`number`] against another server and cache directory.
pub fn number_at(
    url: &str,
    cache: &Path,
    sequence: &str,
    flag: &str,
) -> Result<String, FetchError> {
    let valid = !sequence.is_empty()
        && sequence.bytes().all(|b| b.is_ascii_uppercase())
        && matches!(flag, "-k" | "-c" | "-m");
    if !valid {
        return Err(FetchError::BadRequest(format!("{sequence} {flag}")));
    }
    let path = cache_file(cache, sequence, flag);
    if let Some(body) = read_cached(&path, sequence) {
        return Ok(body);
    }
    let body = get(&format!("{url}?plain=1&aaseq={sequence}&scheme={flag}"))?;
    if !body.trim().is_empty() {
        let _ = write_cached(&path, sequence, &body);
    }
    Ok(body)
}

/// FNV-1a, which unlike the std hashers is fixed across releases.
fn fnv(text: &str) -> u64 {
    text.bytes().fold(0xcbf2_9ce4_8422_2325, |h, b| {
        (h ^ u64::from(b)).wrapping_mul(0x0100_0000_01b3)
    })
}

fn cache_file(cache: &Path, sequence: &str, flag: &str) -> PathBuf {
    let name = format!("{}-{:016x}.txt", flag.trim_start_matches('-'), fnv(sequence));
    cache.join("abnum").join(name)
}

/// The cached body, whose first line names the sequence it answers so a
/// hash collision cannot return another chain's numbers.
fn read_cached(path: &Path, sequence: &str) -> Option<String> {
    let text = std::fs::read_to_string(path).ok()?;
    let (head, body) = text.split_once('\n')?;
    (head.strip_prefix("# ")? == sequence).then(|| body.to_string())
}

fn write_cached(path: &Path, sequence: &str, body: &str) -> std::io::Result<()> {
    if let Some(dir) = path.parent() {
        std::fs::create_dir_all(dir)?;
    }
    let partial = path.with_extension("part");
    std::fs::write(&partial, format!("# {sequence}\n{body}"))?;
    std::fs::rename(&partial, path)
}

fn get(url: &str) -> Result<String, FetchError> {
    let mut last = LAST_REQUEST.lock().unwrap_or_else(PoisonError::into_inner);
    if let Some(done) = *last {
        std::thread::sleep(SPACING.saturating_sub(done.elapsed()));
    }
    let reply = request(url);
    *last = Some(Instant::now());
    reply
}

fn request(url: &str) -> Result<String, FetchError> {
    let agent: ureq::Agent = ureq::Agent::config_builder()
        .timeout_global(Some(TIMEOUT))
        .user_agent(USER_AGENT)
        .build()
        .into();
    let http = |source| FetchError::Http {
        url: url.to_string(),
        source,
    };
    let response = match agent.get(url).call() {
        Ok(r) => r,
        Err(ureq::Error::StatusCode(status)) => {
            return Err(FetchError::NotFound {
                id: url.to_string(),
                what: "numbering".to_string(),
                status,
            })
        }
        Err(source) => return Err(http(source)),
    };
    response.into_body().read_to_string().map_err(http)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn scratch(name: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("vizviz-abnum-{name}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        dir
    }

    #[test]
    fn a_cached_reply_is_returned_without_a_request() {
        let dir = scratch("hit");
        let path = cache_file(&dir, "EVQL", "-k");
        write_cached(&path, "EVQL", "H1 E\n").unwrap();
        let got = number_at("http://127.0.0.1:1/never", &dir, "EVQL", "-k").unwrap();
        std::fs::remove_dir_all(&dir).unwrap();
        assert_eq!(got, "H1 E\n");
    }

    #[test]
    fn a_cache_entry_for_another_sequence_is_not_used() {
        let dir = scratch("clash");
        let path = cache_file(&dir, "EVQL", "-k");
        write_cached(&path, "DIQM", "L1 D\n").unwrap();
        assert!(read_cached(&path, "EVQL").is_none());
        std::fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn an_unreachable_server_is_an_error_and_caches_nothing() {
        let dir = scratch("down");
        let got = number_at("http://127.0.0.1:1/abnum.cgi", &dir, "EVQL", "-c");
        assert!(matches!(got, Err(FetchError::Http { .. })));
        assert!(!dir.join("abnum").exists());
    }

    #[test]
    fn only_capitals_and_known_schemes_are_sent() {
        let dir = scratch("bad");
        for (seq, flag) in [("evql", "-k"), ("EV&QL", "-k"), ("EVQL", "-x"), ("", "-k")] {
            let got = number_at("http://127.0.0.1:1/", &dir, seq, flag);
            assert!(matches!(got, Err(FetchError::BadRequest(_))), "{seq} {flag}");
        }
    }
}
