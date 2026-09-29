# Driving a running app

Start the app with `--listen` and other programs can run commands in its
window: a Python script, an AI agent through MCP, anything that can open a
TCP socket. Everything goes through the command language
([COMMANDS.md](COMMANDS.md)), the same one `--exec`, the console and the
ribbon use.

## Python

```python
import vizviz

app = vizviz.launch()  # starts `vizviz --listen`; or vizviz.connect()
app.load("fixtures/small/4HHB.cif")
app.exec("rep cartoon; color chain; lighting full")
app.render("hemoglobin.png", 3840, 2160, samples=64)
app.screenshot("view.png")
```

- `exec(script)` runs command lines and returns the log lines they wrote.
  A failing line raises `vizviz.LiveError` (its `output` has the log) and
  skips the rest of that script; the window stays open.
- A call returns once everything it started is finished: a background
  surface build, a render. Long renders keep the call waiting, so there is
  no timeout by default (`connect(timeout=...)` sets one).
- Paths in `exec` scripts are opened by the app, so relative paths resolve
  in the app's working directory. `load`, `screenshot`, `render`,
  `save_session` and `load_session` make theirs absolute first.
- Commands that would open a dialog (`screenshot`, `open`, `fetch`,
  `savesession`, `loadsession`, `savelayout` with no path) are refused.
- No compiled extension is needed: `vizviz.live` is plain Python.

## MCP

`python -m vizviz.mcp` is a Model Context Protocol server on stdio. It
connects to the window this user started with `--listen` (or starts one)
and offers four tools: `run_commands`, `render_image`, `screenshot`,
`list_commands`. To use it from an MCP client, register the command, for
example:

```json
{"mcpServers": {"vizviz": {"command": "python", "args": ["-m", "vizviz.mcp"]}}}
```

(with `PYTHONPATH` pointing at `crates/vv-py/python` when running from a
checkout).

## The protocol

- The app listens on 127.0.0.1 only, on a port the OS picks, and writes
  `live.json` to its config directory (`%APPDATA%\vizviz`,
  `~/Library/Application Support/vizviz`, `~/.config/vizviz`): the port,
  a random token, the process id and the command-language version. The
  file is removed when the app exits.
- One JSON object per line each way. First `{"hello": TOKEN}`, answered
  `{"ok": true, "vizviz": VERSION, "language": VERSION}`. Then
  `{"exec": SCRIPT}`, answered `{"ok": BOOL, "output": [LINES]}` once
  the script has run.
- Requests run one at a time, in order. `quit` closes the window before
  answering.

`examples/live_check.py` runs all of this end to end.
