# Containerized tool plugins

Third integration mechanism alongside native (linked) and subprocess/IPC
plugins (see ENVIRONMENT.md). Containers are for heavy third-party tools
with their own complex dependency stacks — structure prediction, protein
design, docking, MD — where matching their exact CUDA/PyTorch/conda
environment natively isn't realistic or maintainable.

## Why containers, specifically

- Tools like AlphaFold, RFdiffusion, or GROMACS pin exact library versions
  that conflict with each other and with our own toolchain. A container
  sidesteps that instead of us maintaining N incompatible environments.
- Same image runs locally (Docker) and on HPC (Apptainer/Singularity, see
  HPC.md) — one artifact, two runtimes, no separate maintenance.
- Matches how this ecosystem already works in practice — SaaS platforms in
  this space orchestrate exactly this way.

## Plugin contract (draft)

A container plugin is a manifest + image, not custom code per tool:

- **Manifest** declares: image reference, required GPU (yes/no, min VRAM),
  pinned CUDA/driver version (this ecosystem's GPU support breaks per
  release — treat as a per-tool moving target, not a one-time pin), input
  mount (read-only), output mount, and a job description format.
- **I/O contract:** input files mounted read-only into a fixed path,
  output written to a fixed path, job completion signaled by exit code +
  an output manifest file — no custom RPC per tool.
- **Weights/binary-provisioning contract:** a declared mount point/cache
  dir + terms-acceptance flag, separate from the image itself. Covers two
  distinct cases, both real and both first-class manifest fields, not
  afterthoughts:
  - *Gated-but-permissive* (e.g. OpenFold3 weights) — accept terms once,
    then cache.
  - *Never redistributable* — AlphaFold3's weights, and entire binaries
    like **NAMD, AMBER's `pmemd`, ORCA, and commercial docking programs**: license
    terms (registration requirements, EULAs, no-redistribution clauses)
    mean no container anywhere in the ecosystem — not us, not
    BioContainers, not NGC — is legally allowed to ship them prebuilt.
    The platform's own "bring-your-own-licensed-binary" path (user
    supplies their own registered install, we build/mount at install
    time) is mandatory for these, not optional polish.
- Runs through the same plugin/agent tool-call API as everything else (see
  VERSIONING.md) — submitting a containerized docking run looks the same
  to the UI or an AI agent as submitting an HPC job, because under the
  hood it often *is* one (see HPC.md).
- **The manifest doubles as the agent-readable description** (typed
  inputs/outputs, one-line purpose, example invocation) — the same file a
  human reads is what an AI agent reads to decide when to call the tool.
  No separate "agent integration" step per plugin.
- **Run fingerprinting:** a job's inputs + parameters + image digest are
  hashed; an identical request returns the cached result instead of
  re-running a 2-hour job. Opt-out per job.

## Container convention: build our own, on plain OCI

The OCI/Docker format itself is a solved problem — Apptainer/Singularity
converts `docker://` images natively, nothing to invent there. But no
CADD-specific *convention* exists to build the manifest/runtime contract
on: checked BioContainers/bioconda (thin, years stale for most tools —
e.g. Vina 1.1.2 vs current 1.2.7), NVIDIA BioNeMo/NIM (vendor-locked,
partly paywalled, missing tools), and nf-core/proteinfold
(Nextflow-specific, prediction only, excludes design tools). conda-forge,
not BioContainers, is actually the current channel for classical MD
engines — useful to know, not something to build our format on. Decision
stands: **plain OCI images, our own thin manifest** as defined above.

Known Apptainer/HPC gotchas to design around, not discover in production:

- `--nv` exposes **all** host GPUs by default — opposite of Docker's
  opt-in model. Use `CUDA_VISIBLE_DEVICES` to scope it.
- Host `$HOME`/`/tmp` auto-mount into the container by default — a data
  leakage risk if not disabled explicitly.
- Dockerfile `USER` directives are silently ignored.
- MPI ABI mismatches between container and host cause **silent** failures
  (e.g. every rank reporting as rank 0), not a clean error — a multi-node
  MD job's container needs to match the host's MPI ABI, which sometimes
  means building MPI from source inside the image for that specific
  cluster (this is what NERSC documents doing).

## Known risk: license traps hide below the top-level LICENSE

A tool's own license can be clean while an installed sub-dependency isn't
— e.g. BindCraft is MIT/CC-BY but its install script pulls in PyRosetta,
which requires a paid license for commercial use. Dependency license
review (CONTRIBUTING.md) has to check what a container actually installs,
not just the wrapped tool's own repo.
