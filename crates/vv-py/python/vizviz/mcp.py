"""An MCP server for a running vizviz window (docs/LIVE.md).

    python -m vizviz.mcp

Speaks the Model Context Protocol over stdio (JSON-RPC 2.0, one message
per line; nothing else goes to stdout) and forwards its tools to the
window through `vizviz.live`: it connects to the vizviz this user started
with `--listen`, or starts one. Plain Python, no compiled extension.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

from .live import LiveApp, LiveError, connect, launch

SUPPORTED = ["2025-06-18", "2025-03-26", "2024-11-05"]


def _text(schema_props: dict, required: list[str]) -> dict:
    return {"type": "object", "properties": schema_props, "required": required}


TOOLS = [
    {
        "name": "run_commands",
        "description": (
            "Run vizviz commands in the live window, one per line or `;`-separated "
            "(e.g. `load /abs/4HHB.cif; rep cartoon; color chain; lighting full`). "
            "Returns the log lines. `list_commands` lists every command. Use "
            "absolute paths: relative ones resolve in the app's directory."
        ),
        "inputSchema": _text({"script": {"type": "string"}}, ["script"]),
    },
    {
        "name": "render_image",
        "description": (
            "Path-trace the current view (real shadows and ambient occlusion) to a "
            "PNG or JPEG at `path`; returns once the file is written."
        ),
        "inputSchema": _text(
            {
                "path": {"type": "string"},
                "width": {"type": "integer", "minimum": 1},
                "height": {"type": "integer", "minimum": 1},
                "samples": {"type": "integer", "minimum": 1, "default": 64},
                "transparent": {"type": "boolean", "default": False},
            },
            ["path"],
        ),
    },
    {
        "name": "screenshot",
        "description": "Save the viewport as it looks now (PNG, JPEG or SVG) at `path`.",
        "inputSchema": _text({"path": {"type": "string"}}, ["path"]),
    },
    {
        "name": "list_commands",
        "description": "Every vizviz command with its usage.",
        "inputSchema": _text({}, []),
    },
]


class Server:
    def __init__(self) -> None:
        self.app: LiveApp | None = None

    def window(self) -> LiveApp:
        if self.app is None:
            try:
                self.app = connect()
            except LiveError:
                self.app = launch()
        return self.app

    def call(self, name: str, args: dict) -> list[str]:
        app = self.window()
        if name == "run_commands":
            return app.exec(args["script"])
        if name == "render_image":
            return app.render(
                args["path"],
                args.get("width"),
                args.get("height"),
                int(args.get("samples", 64)),
                bool(args.get("transparent", False)),
            )
        if name == "screenshot":
            return app.screenshot(args["path"])
        if name == "list_commands":
            return app.exec("help")
        raise LiveError(f"no tool named `{name}`")

    def handle(self, message: dict) -> dict | None:
        method, id_ = message.get("method"), message.get("id")
        if id_ is None:
            return None  # a notification
        if method == "initialize":
            asked = message.get("params", {}).get("protocolVersion")
            result = {
                "protocolVersion": asked if asked in SUPPORTED else SUPPORTED[0],
                "capabilities": {"tools": {}},
                "serverInfo": {"name": "vizviz", "version": "1"},
            }
        elif method == "ping":
            result = {}
        elif method == "tools/list":
            result = {"tools": TOOLS}
        elif method == "tools/call":
            params = message.get("params", {})
            try:
                lines = self.call(params.get("name", ""), params.get("arguments") or {})
                result = {"content": [{"type": "text", "text": "\n".join(lines)}], "isError": False}
            except (LiveError, OSError, KeyError, ValueError) as e:
                output = getattr(e, "output", [])
                text = "\n".join([*output, str(e)] if output else [str(e)])
                result = {"content": [{"type": "text", "text": text}], "isError": True}
        else:
            return {
                "jsonrpc": "2.0",
                "id": id_,
                "error": {"code": -32601, "message": f"no method `{method}`"},
            }
        return {"jsonrpc": "2.0", "id": id_, "result": result}


def main() -> None:
    server = Server()
    for line in sys.stdin:
        if not line.strip():
            continue
        try:
            message = json.loads(line)
        except json.JSONDecodeError as e:
            reply = {"jsonrpc": "2.0", "id": None, "error": {"code": -32700, "message": str(e)}}
        else:
            reply = server.handle(message)
        if reply is not None:
            sys.stdout.write(json.dumps(reply) + "\n")
            sys.stdout.flush()


if __name__ == "__main__":
    print(f"vizviz MCP server ({Path(__file__).name}) on stdio", file=sys.stderr)
    main()
