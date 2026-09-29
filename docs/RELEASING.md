# Releasing

A release is a `v*` git tag. `.github/workflows/release.yml` (from
`dist generate`) builds the `vizviz` desktop binary per OS and attaches
archives and shell / PowerShell installers to a GitHub Release, with
CHANGELOG.md's section as the notes. It also releases the `vv-dssp` CLI
the same way; set `dist = false` in its Cargo.toml to stop that.

On every push to `main`, `.github/workflows/docs.yml` publishes the docs
site (`mkdocs.yml`, MkDocs Material) to GitHub Pages.

## One-time setup

1. Create the repository.
2. Set `repository` in the workspace `Cargo.toml`, then run `dist generate`
   and commit `.github/workflows/release.yml` (config is in
   `dist-workspace.toml`; `vv-bench` and `xtask` are excluded).
3. In the repository settings, set Pages > Source to "GitHub Actions".
4. Push `main` and watch `ci.yml` go green on all three OSes before tagging.

## Cutting a release

```bash
# bump `version` in the workspace Cargo.toml (all crates inherit it) and the
# internal crates' `version = ...` pins under [workspace.dependencies]; add
# the release's section to CHANGELOG.md (the GitHub Release notes)
git commit -am "Release vX.Y.Z"
git tag vX.Y.Z
git push && git push --tags
```

Then check that the GitHub Release has five archives and two installers.

The Python package is not published; build its wheel locally:

## Local build

```bash
cargo install --path crates/vv-app --locked --force    # vizviz on PATH
uv run --project crates/vv-py maturin build --release --manifest-path crates/vv-py/Cargo.toml --out dist
pip install dist/vizviz-*-cp310-abi3-*.whl             # into any Python >= 3.10
```

Check that the wheel, installed into a fresh venv, passes the full pytest
suite with `vizviz.__file__` inside that venv's `site-packages`, and that
the installed `vizviz`, run from an unrelated directory, produces a figure
via `--exec` and keeps its layout in the per-user config directory.
