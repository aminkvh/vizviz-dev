"""The live client and MCP server without a running app (docs/LIVE.md)."""

import os
import subprocess
import sys
from pathlib import Path

import pytest
from vizviz import mcp
from vizviz.live import LiveError

PYTHON_SOURCE = Path(__file__).resolve().parents[1] / "python"


class FakeApp:
    """Stands in for a `LiveApp`: records scripts, fails on `nonsense`."""

    def __init__(self):
        self.scripts = []

    def exec(self, script):
        self.scripts.append(script)
        if "nonsense" in script:
            raise LiveError("error: unknown", [f"> {script}", "error: unknown"])
        return [f"> {script}", "ok"]

    def render(self, path, width, height, samples, transparent):
        return self.exec(f"render {path} {width}x{height} {samples}")

    def screenshot(self, path):
        return self.exec(f"screenshot {path}")


@pytest.fixture
def server(monkeypatch):
    app = FakeApp()
    s = mcp.Server()
    monkeypatch.setattr(s, "window", lambda: app)
    s.app_for_test = app
    return s


def call(server, method, params=None, id_=1):
    return server.handle({"jsonrpc": "2.0", "id": id_, "method": method, "params": params or {}})


def test_initialize_echoes_a_supported_version_and_falls_back_otherwise(server):
    assert call(server, "initialize", {"protocolVersion": "2025-03-26"})["result"][
        "protocolVersion"
    ] == ("2025-03-26")
    fallback = call(server, "initialize", {"protocolVersion": "1999-01-01"})["result"]
    assert fallback["protocolVersion"] == mcp.SUPPORTED[0]
    assert "tools" in fallback["capabilities"]


def test_notifications_get_no_reply(server):
    assert server.handle({"jsonrpc": "2.0", "method": "notifications/initialized"}) is None


def test_tools_are_listed_with_schemas(server):
    tools = call(server, "tools/list")["result"]["tools"]
    names = {t["name"] for t in tools}
    assert names == {"run_commands", "render_image", "screenshot", "list_commands"}
    assert all(t["inputSchema"]["type"] == "object" for t in tools)


def test_unknown_methods_are_json_rpc_errors(server):
    assert call(server, "resources/list")["error"]["code"] == -32601


def test_tool_results_and_tool_errors(server):
    ok = call(
        server, "tools/call", {"name": "run_commands", "arguments": {"script": "rep cartoon"}}
    )
    assert ok["result"]["isError"] is False
    assert "ok" in ok["result"]["content"][0]["text"]
    bad = call(server, "tools/call", {"name": "run_commands", "arguments": {"script": "nonsense"}})
    assert "error" not in bad, "a tool failure is a result, not a protocol error"
    assert bad["result"]["isError"] is True
    assert "error: unknown" in bad["result"]["content"][0]["text"]
    missing = call(server, "tools/call", {"name": "no_such_tool", "arguments": {}})
    assert missing["result"]["isError"] is True


def test_render_passes_its_arguments_through(server):
    call(
        server,
        "tools/call",
        {"name": "render_image", "arguments": {"path": "/tmp/x.png", "width": 64, "height": 48}},
    )
    assert server.app_for_test.scripts[-1] == "render /tmp/x.png 64x48 64"


def test_the_live_client_needs_no_compiled_extension():
    code = (
        "import sys; import vizviz, vizviz.live, vizviz.mcp; "
        "assert 'vizviz._core' not in sys.modules, 'loaded _core'"
    )
    env = {**os.environ, "PYTHONPATH": str(PYTHON_SOURCE)}
    subprocess.run([sys.executable, "-c", code], check=True, env=env)
