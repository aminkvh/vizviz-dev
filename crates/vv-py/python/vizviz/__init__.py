"""vizviz: fast, open molecular visualization and analysis.

    import vizviz
    s = vizviz.load("4HHB.cif")
    s.positions           # (n, 3) float32 — a zero-copy, read-only NumPy view
    s.element             # (n,) uint8 atomic numbers
    s.residue_names       # ["VAL", "LEU", ...]
    img = s.render(800, 600, coloring="chain")   # (600, 800, 4) uint8, CPU-rendered
    vizviz.write_png("hemoglobin.png", img)

`vizviz.from_arrays(...)` builds a Structure from your own NumPy columns;
`from_mdanalysis` / `to_mdanalysis` move structures to and from MDAnalysis.
`s.render(..., values=sasa)` colors by any per-atom number. A `Session`
drives the same command bus as the desktop app (load, select, undo/redo,
`set_values`). `vizviz.connect()` / `vizviz.launch()` drive a running
window (docs/LIVE.md). See docs/PYTHON.md in the repository for the full
API.
"""

from .live import LiveApp, LiveError, connect, launch

# The compiled half loads on first use (PEP 562), so the live client and
# the MCP server (`python -m vizviz.mcp`) work with plain Python.
_CORE = {
    "HAS_RENDER",
    "Session",
    "Structure",
    "__version__",
    "element_symbol",
    "fetch",
    "load",
    "load_many",
    "write_png",
}
_INTEROP = {"from_arrays", "from_mdanalysis", "to_mdanalysis"}

__all__ = [
    "HAS_RENDER",
    "LiveApp",
    "LiveError",
    "Session",
    "Structure",
    "__version__",
    "connect",
    "element_symbol",
    "fetch",
    "from_arrays",
    "from_mdanalysis",
    "launch",
    "load",
    "load_many",
    "to_mdanalysis",
    "write_png",
]


def __getattr__(name: str):
    if name in _CORE:
        from . import _core

        if name == "write_png" and not _core.HAS_RENDER:
            raise AttributeError("vizviz was built without the `render` feature")
        return getattr(_core, name)
    if name in _INTEROP:
        from . import interop

        return getattr(interop, name)
    raise AttributeError(f"module 'vizviz' has no attribute {name!r}")
