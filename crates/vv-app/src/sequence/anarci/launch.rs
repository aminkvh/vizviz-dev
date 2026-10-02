//! Finding and running the ANARCI executable.

use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};
use std::sync::atomic::{AtomicU32, Ordering};

use vv_core::antibody::Scheme;

use super::Error;

pub const ENV_VAR: &str = "VIZVIZ_ANARCI";

/// Reads the FASTA on stdin into a temp file (ANARCI only takes a file
/// or a literal sequence) and runs `$1` on it; `$1`'s directory joins the
/// PATH so a conda install finds its `hmmscan` and `python`. Free of
/// quotes, which `wsl.exe` would mangle.
const WSL_SCRIPT: &str =
    "exe=$1; shift; case $exe in /*) PATH=${exe%/*}:$PATH;; esac; f=$(mktemp) || exit 1; cat >$f; $exe -i $f $@; st=$?; unlink $f; exit $st";

/// Where ANARCI runs.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Launcher {
    Direct(PathBuf),
    /// An executable inside WSL (a bare name is looked up on its PATH).
    Wsl(String),
}

impl Launcher {
    /// The setting, then `VIZVIZ_ANARCI`, then `ANARCI` on the PATH, then
    /// `ANARCI` inside WSL.
    pub fn locate(configured: Option<&Path>) -> Result<Launcher, Error> {
        let env = std::env::var(ENV_VAR).ok();
        let dirs: Vec<PathBuf> = std::env::var_os("PATH")
            .map(|p| std::env::split_paths(&p).collect())
            .unwrap_or_default();
        Self::locate_with(configured, env.as_deref(), &dirs, cfg!(windows))
    }

    pub fn locate_with(
        configured: Option<&Path>,
        env: Option<&str>,
        dirs: &[PathBuf],
        try_wsl: bool,
    ) -> Result<Launcher, Error> {
        let explicit = configured
            .map(|p| p.to_string_lossy().into_owned())
            .or_else(|| env.map(str::to_string))
            .filter(|s| !s.trim().is_empty());
        if let Some(text) = explicit {
            return Self::explicit(text.trim());
        }
        if let Some(found) = on_path("ANARCI", dirs) {
            return Ok(Launcher::Direct(found));
        }
        if try_wsl && wsl_has("ANARCI") {
            return Ok(Launcher::Wsl("ANARCI".into()));
        }
        Err(Error::NotFound)
    }

    /// A path to use as given. On Windows a path that starts with `/` or
    /// `wsl:` is inside WSL.
    fn explicit(text: &str) -> Result<Launcher, Error> {
        let inside = text.strip_prefix("wsl:");
        let wsl_path = inside.or_else(|| (cfg!(windows) && text.starts_with('/')).then_some(text));
        if let Some(path) = wsl_path {
            return match path.contains(char::is_whitespace) {
                true => Err(Error::Failed(format!("`{path}` has a space in it"))),
                false => Ok(Launcher::Wsl(path.to_string())),
            };
        }
        match Path::new(text).is_file() {
            true => Ok(Launcher::Direct(PathBuf::from(text))),
            false => Err(Error::Failed(format!("`{text}` is not a file"))),
        }
    }

    /// ANARCI's text output for the FASTA `fasta` under `scheme`. Receptor
    /// chains only take IMGT, so the other schemes are restricted to
    /// antibodies.
    pub fn run(&self, fasta: &str, scheme: Scheme) -> Result<String, Error> {
        let args = scheme_args(scheme);
        let output = match self {
            Launcher::Direct(exe) => run_direct(exe, fasta, &args),
            Launcher::Wsl(exe) => run_wsl(exe, fasta, &args),
        };
        finish(output)
    }
}

fn scheme_args(scheme: Scheme) -> Vec<&'static str> {
    let word = match scheme {
        Scheme::Imgt => "imgt",
        Scheme::Kabat => "kabat",
        Scheme::Chothia => "chothia",
        Scheme::Martin => "martin",
        Scheme::Aho => "aho",
    };
    let mut args = vec!["-s", word];
    if scheme != Scheme::Imgt {
        args.extend(["-r", "ig"]);
    }
    args
}

fn on_path(name: &str, dirs: &[PathBuf]) -> Option<PathBuf> {
    let suffixes: &[&str] = if cfg!(windows) {
        &["", ".exe", ".cmd", ".bat"]
    } else {
        &[""]
    };
    dirs.iter()
        .flat_map(|d| suffixes.iter().map(move |s| d.join(format!("{name}{s}"))))
        .find(|p| p.is_file())
}

/// Stops a console window flashing up when the GUI starts a child.
fn quiet(cmd: &mut Command) -> &mut Command {
    #[cfg(windows)]
    {
        use std::os::windows::process::CommandExt;
        cmd.creation_flags(0x0800_0000);
    }
    cmd
}

fn wsl_has(name: &str) -> bool {
    let mut cmd = Command::new("wsl.exe");
    cmd.args(["-e", "sh", "-lc"])
        .arg(format!("command -v {name}"));
    quiet(&mut cmd)
        .stdin(Stdio::null())
        .output()
        .is_ok_and(|o| o.status.success() && !o.stdout.is_empty())
}

/// A FASTA file in the temp directory, removed when dropped.
struct TempFasta(PathBuf);

impl TempFasta {
    fn write(fasta: &str) -> std::io::Result<Self> {
        static COUNTER: AtomicU32 = AtomicU32::new(0);
        let n = COUNTER.fetch_add(1, Ordering::Relaxed);
        let name = format!("vizviz-anarci-{}-{n}.fa", std::process::id());
        let path = std::env::temp_dir().join(name);
        std::fs::write(&path, fasta)?;
        Ok(Self(path))
    }
}

impl Drop for TempFasta {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.0);
    }
}

fn run_direct(exe: &Path, fasta: &str, args: &[&str]) -> std::io::Result<Output> {
    let file = TempFasta::write(fasta)?;
    let mut cmd = Command::new(exe);
    cmd.arg("-i").arg(&file.0).args(args);
    if let (Some(dir), Some(path)) = (
        exe.parent().filter(|d| d.is_absolute()),
        std::env::var_os("PATH"),
    ) {
        let joined = std::env::join_paths(
            std::iter::once(dir.to_path_buf()).chain(std::env::split_paths(&path)),
        );
        if let Ok(joined) = joined {
            cmd.env("PATH", joined);
        }
    }
    quiet(&mut cmd).stdin(Stdio::null()).output()
}

fn run_wsl(exe: &str, fasta: &str, args: &[&str]) -> std::io::Result<Output> {
    let mut cmd = Command::new("wsl.exe");
    cmd.args(["-e", "sh", "-lc", WSL_SCRIPT, "sh", exe])
        .args(args);
    let mut child = quiet(&mut cmd)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()?;
    if let Some(mut stdin) = child.stdin.take() {
        stdin.write_all(fasta.as_bytes())?;
    }
    child.wait_with_output()
}

fn finish(output: std::io::Result<Output>) -> Result<String, Error> {
    let output = output.map_err(|e| Error::Failed(e.to_string()))?;
    if output.status.success() {
        return Ok(String::from_utf8_lossy(&output.stdout).into_owned());
    }
    let stderr = String::from_utf8_lossy(&output.stderr);
    let reason = stderr
        .lines()
        .rev()
        .find(|l| !l.trim().is_empty())
        .map_or_else(|| output.status.to_string(), |l| l.trim().to_string());
    Err(Error::Failed(reason))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nothing_configured_nor_on_the_path_is_not_found() {
        let found = Launcher::locate_with(None, None, &[], false);
        assert!(matches!(found, Err(Error::NotFound)));
        let blank = Launcher::locate_with(None, Some("  "), &[], false);
        assert!(matches!(blank, Err(Error::NotFound)));
    }

    #[test]
    fn the_setting_wins_over_the_environment_and_must_be_a_file() {
        let here = Path::new(env!("CARGO_MANIFEST_DIR")).join("Cargo.toml");
        let found = Launcher::locate_with(Some(&here), Some("/nowhere"), &[], false);
        assert_eq!(found.unwrap(), Launcher::Direct(here));
        let missing = Launcher::locate_with(Some(Path::new("no-such-anarci")), None, &[], false);
        assert!(matches!(missing, Err(Error::Failed(m)) if m.contains("not a file")));
    }

    #[test]
    fn a_wsl_prefix_names_a_path_inside_wsl() {
        let found = Launcher::locate_with(None, Some("wsl:/opt/bin/ANARCI"), &[], false);
        assert_eq!(found.unwrap(), Launcher::Wsl("/opt/bin/ANARCI".into()));
    }

    #[test]
    fn the_path_is_searched_for_the_executable() {
        let dir = std::env::temp_dir().join(format!("vizviz-find-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let name = if cfg!(windows) {
            "ANARCI.exe"
        } else {
            "ANARCI"
        };
        std::fs::write(dir.join(name), "").unwrap();
        let found = Launcher::locate_with(None, None, std::slice::from_ref(&dir), false);
        std::fs::remove_dir_all(&dir).unwrap();
        assert_eq!(found.unwrap(), Launcher::Direct(dir.join(name)));
    }

    #[test]
    fn receptor_chains_keep_imgt_unrestricted() {
        assert_eq!(scheme_args(Scheme::Imgt), ["-s", "imgt"]);
        assert_eq!(scheme_args(Scheme::Kabat), ["-s", "kabat", "-r", "ig"]);
    }
}
