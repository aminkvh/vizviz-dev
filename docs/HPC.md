# HPC integration strategy

## The constraint

HPC compute nodes usually have no outbound network and no GPU suited for
interactive rendering. Streaming rendered frames back from a cluster isn't
realistic for most academic HPC. The real architecture is:

**Compute remote, render local.**

- Heavy compute (MD runs, docking, ML inference, large-scale trajectory
  analysis) runs as a scheduled job on the cluster.
- The interactive viewer stays on the user's own machine and GPU.
- Results (trajectories, structures, analysis output) sync back to local
  storage for the viewer to load — not streamed.

Job submission, status, and result staging are first-class from the start,
not retrofitted.

## Scheduler support

- Slurm first (most common in structural biology HPC).
- Connector interface designed to add PBS/Torque and LSF later without a
  redesign.
- Job submission goes through the same plugin/agent tool-call API as
  everything else — an AI agent can submit and monitor a job the same way
  a human does through the UI.

## Environment parity

The same `pixi` lockfiles (see ENVIRONMENT.md) make a job submitted from a
laptop reproducible on a cluster node. Cluster execution uses
Apptainer/Singularity images built from the same environment spec.

## Out of scope for v1

Full remote-rendering/pixel-streaming from a cluster GPU is a different,
harder problem than job submission + result sync. Not needed for the model
above; revisit later if a real use case demands it.
