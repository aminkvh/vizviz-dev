"""End-to-end check of the live SDK (docs/LIVE.md), with no vizviz running
at the start: an MCP client (through `python -m vizviz.mcp`, which starts
the window) loads, styles and renders; then a Python script connects to
that window and does the same, checks that failures fail only their
request, and closes it.

    PYTHONPATH=crates/vv-py/python py -3 examples/live_check.py OUT_DIR

`vizviz` on PATH is the app it starts. Run it from any directory: every
path it hands the app is absolute.
"""

import contextlib
import json
import os
import subprocess
import sys
from pathlib import Path

import vizviz

repo = Path(__file__).resolve().parent.parent
out = Path(sys.argv[1]).resolve()
out.mkdir(parents=True, exist_ok=True)
with contextlib.suppress(vizviz.LiveError):
    vizviz.connect()
    raise SystemExit("a vizviz is already listening; close it first")

# An MCP client; the server starts the window.
server = subprocess.Popen(
    [sys.executable, "-m", "vizviz.mcp"],
    stdin=subprocess.PIPE,
    stdout=subprocess.PIPE,
    text=True,
    env={**os.environ, "PYTHONPATH": str(repo / "crates/vv-py/python")},
)


def rpc(method, params=None, id_=None):
    message = {"jsonrpc": "2.0", "method": method, "params": params or {}}
    if id_ is not None:
        message["id"] = id_
    server.stdin.write(json.dumps(message) + "\n")
    server.stdin.flush()
    return json.loads(server.stdout.readline()) if id_ is not None else None


def tool(name, arguments, id_):
    return rpc("tools/call", {"name": name, "arguments": arguments}, id_)["result"]


init = rpc(
    "initialize",
    {"protocolVersion": "2025-06-18", "capabilities": {}, "clientInfo": {"name": "check"}},
    1,
)
assert init["result"]["protocolVersion"] == "2025-06-18", init
rpc("notifications/initialized")
tools = {t["name"] for t in rpc("tools/list", id_=2)["result"]["tools"]}
assert {"run_commands", "render_image", "screenshot", "list_commands"} <= tools, tools
script = f"load {repo / 'fixtures/small/1UBQ.cif'}; rep cartoon; color rainbow"
assert not tool("run_commands", {"script": script}, 3)["isError"]
path = out / "mcp_render.png"
args = {"path": str(path), "width": 800, "height": 600, "samples": 16}
rendered = tool("render_image", args, 4)
assert not rendered["isError"] and path.stat().st_size > 10_000, rendered
assert tool("run_commands", {"script": "color nonsense"}, 5)["isError"]
print("mcp:", rendered["content"][0]["text"].splitlines()[-1])
server.stdin.close()
server.wait()

# Python connects to the window the MCP server started.
app = vizviz.connect()
print("connected: vizviz", app.version, "language", app.language)
app.load(repo / "fixtures/small/4HHB.cif")
app.exec("rep cartoon; color chain; addrep sticks hetero and not water; lighting full")
app.screenshot(out / "live_screenshot.png")
app.render(out / "live_render.png", 1600, 1200, samples=32)
for name in ["live_screenshot.png", "live_render.png"]:
    assert (out / name).stat().st_size > 10_000, name

# A failing line fails its request, not the app; dialogs are refused.
for bad in ["rep nonsense", "screenshot"]:
    try:
        app.exec(bad)
    except vizviz.LiveError as e:
        print(f"refused `{bad}`: {e}")
    else:
        raise AssertionError(f"`{bad}` should fail")
assert app.exec("version") == ["> version", "command language 1.0.0"]

with contextlib.suppress(vizviz.LiveError):
    app.exec("quit")  # the window closes before answering
print("live SDK check passed")
