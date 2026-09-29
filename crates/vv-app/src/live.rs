//! `--listen`: other programs drive the running app (docs/LIVE.md).
//!
//! A TCP server on 127.0.0.1 only, on a port the OS picks. The port, a
//! random token, the process id and the command-language version go in a
//! per-user file (`live.json` in the config directory) that clients read;
//! the file goes when the app exits. A client sends one JSON object per
//! line: first `{"hello": TOKEN}`, answered with the versions, then
//! `{"exec": SCRIPT}` as often as it likes, each answered once every line
//! of the script has run -- and everything it started (a background
//! build, a render) has finished -- with the log lines it produced.
//! Requests run one at a time, in order, through the same runner as
//! `--exec`; a failing line ends its request, not the app.

use std::io::{BufRead, BufReader, Write};
use std::net::{TcpListener, TcpStream};
use std::path::PathBuf;
use std::sync::mpsc;

use serde_json::{json, Value};

use crate::gpu_cache::Waker;

/// What one `exec` produced.
pub struct Reply {
    pub ok: bool,
    pub output: Vec<String>,
}

/// A script from a client, waiting for the UI thread.
pub struct Request {
    pub lines: Vec<String>,
    pub reply: mpsc::Sender<Reply>,
}

pub struct LiveServer {
    requests: mpsc::Receiver<Request>,
    info: PathBuf,
    pub port: u16,
}

/// Where a running app says how to reach it.
pub fn info_path() -> Option<PathBuf> {
    crate::layout::config_dir().map(|d| d.join("live.json"))
}

/// 128 random bits, hex: the std hasher's per-process random keys.
fn token() -> String {
    use std::hash::{BuildHasher, Hasher};
    (0..2)
        .map(|k| {
            let mut h = std::collections::hash_map::RandomState::new().build_hasher();
            h.write_u64(k);
            format!("{:016x}", h.finish())
        })
        .collect()
}

impl LiveServer {
    /// Listens on 127.0.0.1 (`port`, or one the OS picks for 0) and writes
    /// the connection file; `waker` wakes the window when a request comes.
    pub fn start(port: u16, waker: Waker) -> std::io::Result<Self> {
        let listener = TcpListener::bind(("127.0.0.1", port))?;
        let port = listener.local_addr()?.port();
        let token = token();
        let info = info_path().ok_or_else(|| std::io::Error::other("no config directory"))?;
        if let Some(dir) = info.parent() {
            std::fs::create_dir_all(dir)?;
        }
        let file = json!({
            "port": port,
            "token": token,
            "pid": std::process::id(),
            "language": vv_scene::LANGUAGE_VERSION,
        });
        write_private(&info, serde_json::to_string_pretty(&file)?.as_bytes())?;
        let (tx, requests) = mpsc::channel();
        std::thread::spawn(move || {
            for stream in listener.incoming().flatten() {
                let (tx, token, waker) = (tx.clone(), token.clone(), waker.clone());
                std::thread::spawn(move || serve(stream, &token, &tx, &waker));
            }
        });
        Ok(Self {
            requests,
            info,
            port,
        })
    }

    /// The next waiting request, if any.
    pub fn poll(&self) -> Option<Request> {
        self.requests.try_recv().ok()
    }
}

/// Writes `bytes` to `path` readable by this user only: whoever has the
/// token can drive the app, and so read and write files as this user.
/// `%APPDATA%` is already per-user on Windows; elsewhere the file itself
/// is made 0600 (a shared machine, an HPC login node).
fn write_private(path: &std::path::Path, bytes: &[u8]) -> std::io::Result<()> {
    let mut options = std::fs::OpenOptions::new();
    options.write(true).create(true).truncate(true);
    #[cfg(unix)]
    std::os::unix::fs::OpenOptionsExt::mode(&mut options, 0o600);
    options.open(path)?.write_all(bytes)
}

impl Drop for LiveServer {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.info);
    }
}

fn send(stream: &mut TcpStream, value: Value) -> bool {
    let mut line = value.to_string();
    line.push('\n');
    stream.write_all(line.as_bytes()).is_ok()
}

/// One client connection: the handshake, then its requests in order.
fn serve(stream: TcpStream, token: &str, requests: &mpsc::Sender<Request>, waker: &Waker) {
    let Ok(mut out) = stream.try_clone() else {
        return;
    };
    let mut greeted = false;
    for line in BufReader::new(stream).lines() {
        let Ok(line) = line else {
            return;
        };
        let message: Value = match serde_json::from_str(&line) {
            Ok(m) => m,
            Err(e) => {
                send(
                    &mut out,
                    json!({"ok": false, "error": format!("not JSON: {e}")}),
                );
                continue;
            }
        };
        if !greeted {
            if message.get("hello").and_then(Value::as_str) != Some(token) {
                send(
                    &mut out,
                    json!({"ok": false, "error": "bad or missing token"}),
                );
                return;
            }
            greeted = true;
            let hello = json!({
                "ok": true,
                "vizviz": env!("CARGO_PKG_VERSION"),
                "language": vv_scene::LANGUAGE_VERSION,
            });
            if !send(&mut out, hello) {
                return;
            }
            continue;
        }
        let Some(script) = message.get("exec").and_then(Value::as_str) else {
            send(
                &mut out,
                json!({"ok": false, "error": "expected {\"exec\": SCRIPT}"}),
            );
            continue;
        };
        let lines = vv_scene::split_script(script);
        if let Some(e) = lines
            .iter()
            .find_map(|l| crate::commands::validate_unattended(l).err())
        {
            send(&mut out, json!({"ok": false, "output": [], "error": e}));
            continue;
        }
        let (reply, answer) = mpsc::channel();
        if requests.send(Request { lines, reply }).is_err() {
            return;
        }
        waker();
        let Ok(done) = answer.recv() else {
            return;
        };
        let body = json!({"ok": done.ok, "output": done.output});
        if !send(&mut out, body) {
            return;
        }
    }
}
