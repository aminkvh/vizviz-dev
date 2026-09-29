# vizviz (Python package)

The Python API of [vizviz](https://github.com/aminkvh/vizviz-dev): fast,
open molecular visualization and analysis in Rust, with zero-copy NumPy
views of every per-atom column and a headless CPU renderer.

```python
import vizviz

s = vizviz.load("4HHB.cif")
s.positions  # (n, 3) float32, read-only view — no copy
s.element, s.b_factor, s.name  # more per-atom columns
img = s.render(1200, 900, coloring="chain")
vizviz.write_png("4hhb.png", img)
```

Build from a checkout:

```bash
uv sync            # builds the extension with maturin
uv run pytest
```

Full API: `docs/PYTHON.md` in the repository.
