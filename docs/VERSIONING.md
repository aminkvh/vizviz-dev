# Versioning policy

App-level SemVer matters least here. What breaks projects in this domain
is interface drift and file-format drift, so those get explicit contracts.

## App releases

Standard SemVer. Pre-1.0 while the architecture is unstable: `0.MINOR`
bumps may break the plugin API or file schema, but never silently (see
below).

## Command language and Python API

The two interfaces scripts use today, each versioned with SemVer:

- **The command language** ([COMMANDS.md](COMMANDS.md)): command names,
  their arguments and what they print, used by `--exec`, the console,
  sessions, `Session.exec` and live clients ([LIVE.md](LIVE.md)). Its
  version is `vv_scene::LANGUAGE_VERSION` (now 2.0.0); the `version`
  command prints it and a live client checks it when connecting.
- **The Python API** ([PYTHON.md](PYTHON.md)): the `vizviz` package,
  versioned with the package.

Patch releases fix things; minor releases add commands, arguments or
output lines appended at the end; a major release may remove or rename a
command, change what an argument means, or reword output scripts may
parse. Before a major change the old form keeps working for at least one
minor release and says in its output that it is deprecated and what
replaces it. A live client refuses a different major version.

## Plugin API

Versioned independently of the app (`plugin-api: MAJOR.MINOR`).

- Each plugin manifest declares a compatible `api_min`/`api_max` range. The
  host refuses to load a plugin outside that range rather than failing
  unpredictably later.
- Breaking interface changes bump MAJOR; additive changes bump MINOR.
- The app exposes its `plugin-api` version to plugins and to the agent
  tool-call layer at connect time, so an AI agent can detect a mismatch
  instead of guessing.

## Session / project files

Every saved file stamps a `schema_version`.

- Migrations are forward-only. Opening an older file upgrades it and says
  so. Opening a newer file in an older build is a hard error, not a
  best-effort partial load that silently drops data.
- Each version transition gets its own named, tested migration function
  (e.g. `migrate_v3_to_v4`), not one catch-all parser.

## GPU capability

Hardware ray tracing (wgpu's `EXPERIMENTAL_RAY_QUERY`, Vulkan-only) is
optional, not required:

- Core rendering (impostor/LOD) must run on baseline wgpu backends
  (Vulkan/Metal/DX12), no RT extension needed.
- RT and any other optional GPU feature is probed at startup; the app
  falls back cleanly and reports what's active, never a silent slowdown.
- Same treatment for CUDA-based acceleration if added later: optional,
  probed, with a CPU fallback.

## Environment mismatches

See [ENVIRONMENT.md](ENVIRONMENT.md): a separate concern, but a plugin
manifest also pins the environment it expects.
