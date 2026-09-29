"""Drive a running vizviz window from Python (docs/LIVE.md).

    import vizviz
    app = vizviz.launch()                  # or vizviz.connect() to a running one
    app.load("fixtures/small/4HHB.cif")
    app.exec("rep cartoon; color chain; lighting full")
    app.render("hemoglobin.png", 3840, 2160, samples=64)

The window must have been started with `vizviz --listen` (`launch()` does
that). Plain Python: no compiled extension or NumPy needed.

Paths in `exec` scripts are read by the app, so relative ones resolve in
the app's working directory; `load`, `screenshot`, `render` and the
session helpers make theirs absolute first. A request returns once
everything it started has finished -- a render holds it open until the
image is written -- so there is no timeout by default.
"""

from __future__ import annotations

import json
import os
import shutil
import socket
import subprocess
import sys
import time
from pathlib import Path

__all__ = ["LiveApp", "LiveError", "connect", "info_path", "launch", "LANGUAGE_MAJOR"]

#: The command-language major version this client speaks.
LANGUAGE_MAJOR = 2


class LiveError(RuntimeError):
    """A command failed in the app, or the app could not be reached.

    `output` holds the log lines the request produced, the error last."""

    def __init__(self, message: str, output: list[str] | None = None):
        super().__init__(message)
        self.output = output or []


def info_path() -> Path:
    """Where a listening app writes how to reach it (the app's config
    directory, `live.json`)."""
    if sys.platform == "win32":
        base = Path(os.environ["APPDATA"])
    elif sys.platform == "darwin":
        base = Path.home() / "Library" / "Application Support"
    else:
        base = Path(os.environ.get("XDG_CONFIG_HOME") or Path.home() / ".config")
    return base / "vizviz" / "live.json"


def _alive(pid: int) -> bool:
    if sys.platform == "win32":
        import ctypes

        handle = ctypes.windll.kernel32.OpenProcess(0x1000, False, pid)  # QUERY_LIMITED_INFORMATION
        if not handle:
            return False
        code = ctypes.c_ulong()
        ctypes.windll.kernel32.GetExitCodeProcess(handle, ctypes.byref(code))
        ctypes.windll.kernel32.CloseHandle(handle)
        return code.value == 259  # STILL_ACTIVE
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


class LiveApp:
    """A connection to a running vizviz window."""

    def __init__(self, port: int, token: str, timeout: float | None = None):
        self._sock = socket.create_connection(("127.0.0.1", port))
        self._sock.settimeout(timeout)
        self._file = self._sock.makefile("rwb")
        hello = self._send({"hello": token})
        if not hello.get("ok"):
            raise LiveError(f"vizviz refused the connection: {hello.get('error')}")
        self.version = hello["vizviz"]
        self.language = hello["language"]
        major = int(self.language.split(".")[0])
        if major != LANGUAGE_MAJOR:
            self.close()
            raise LiveError(
                f"the app speaks command language {self.language}; this client speaks "
                f"{LANGUAGE_MAJOR}.x"
            )

    def _send(self, message: dict) -> dict:
        self._file.write(json.dumps(message).encode() + b"\n")
        self._file.flush()
        line = self._file.readline()
        if not line:
            raise LiveError("the app closed the connection")
        return json.loads(line)

    def exec(self, script: str) -> list[str]:
        """Runs a command script (docs/COMMANDS.md; lines or `;`-separated)
        and returns its log lines. Raises `LiveError` when a line fails;
        the rest of that script is skipped."""
        reply = self._send({"exec": script})
        output = reply.get("output", [])
        if not reply.get("ok"):
            error = reply.get("error") or (output[-1] if output else "command failed")
            raise LiveError(error, output)
        return output

    def load(self, path: str | os.PathLike) -> list[str]:
        return self.exec(f"load {Path(path).resolve()}")

    def screenshot(self, path: str | os.PathLike, ssaa: bool = True) -> list[str]:
        return self.exec(f"screenshot {Path(path).resolve()}" + ("" if ssaa else " nossaa"))

    def render(
        self,
        path: str | os.PathLike,
        width: int | None = None,
        height: int | None = None,
        samples: int = 64,
        transparent: bool = False,
    ) -> list[str]:
        """Path-traces the view to `path` (the viewport's size unless
        `width` and `height` are given); returns once it is written."""
        size = f" {width}x{height}" if width and height else ""
        extra = " transparent" if transparent else ""
        return self.exec(f"render {Path(path).resolve()}{size} {samples}{extra}")

    def save_session(self, path: str | os.PathLike) -> list[str]:
        return self.exec(f"savesession {Path(path).resolve()}")

    def load_session(self, path: str | os.PathLike) -> list[str]:
        return self.exec(f"loadsession {Path(path).resolve()}")

    def close(self) -> None:
        self._file.close()
        self._sock.close()

    def __enter__(self) -> LiveApp:
        return self

    def __exit__(self, *exc) -> None:
        self.close()


def connect(timeout: float | None = None) -> LiveApp:
    """Connects to the vizviz window this user started with `--listen`."""
    path = info_path()
    try:
        info = json.loads(path.read_text())
    except FileNotFoundError:
        raise LiveError(
            f"no listening vizviz ({path} is missing): start it with `vizviz --listen`"
        ) from None
    if not _alive(int(info["pid"])):
        raise LiveError(f"the vizviz that wrote {path} (pid {info['pid']}) is gone")
    return LiveApp(int(info["port"]), info["token"], timeout)


def launch(exe: str | os.PathLike | None = None, wait: float = 60.0) -> LiveApp:
    """Starts a vizviz window with `--listen` and connects to it. `exe`
    defaults to `vizviz` on PATH."""
    exe = exe or shutil.which("vizviz")
    if exe is None:
        raise LiveError("no `vizviz` on PATH; pass exe=")
    path = info_path()
    before = path.stat().st_mtime_ns if path.exists() else None
    # Not our stdout: under the MCP server that carries protocol messages only.
    process = subprocess.Popen(
        [str(exe), "--listen"], stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL
    )
    deadline = time.monotonic() + wait
    while time.monotonic() < deadline:
        if process.poll() is not None:
            raise LiveError(f"vizviz exited with code {process.returncode} before listening")
        if path.exists() and path.stat().st_mtime_ns != before:
            try:
                info = json.loads(path.read_text())
            except json.JSONDecodeError:
                pass  # half written
            else:
                if int(info["pid"]) == process.pid:
                    return LiveApp(int(info["port"]), info["token"])
        time.sleep(0.1)
    raise LiveError(f"vizviz did not start listening within {wait} s")
