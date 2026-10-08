# BreezeLab case implementation plans

Prepared 2026-10-07 for agents launched by Greg on wpcluster.

**Start here:** `/shared/home/greg/BreezeLab.jl/plans/README.md`.
Repository: https://github.com/PRECISION-GM/BreezeLab.jl

| Agent assignment | Plan | First concrete deliverable |
| --- | --- | --- |
| ENA-LASSO | [ENA-LASSO.md](ENA-LASSO.md) | A separate LASSO example using an authentic ARM member; data acquisition is the next joint task with Greg |
| SEA STARR | [SEA-STARR.md](SEA-STARR.md) | A driver-backed CTRL constructor and short GPU run, with remaining protocol differences explicit |
| TRACER-MIP | [TRACER-MIP.md](TRACER-MIP.md) | ERA5-backed outer-domain run followed by the specified nested coastal CTRL experiment |
| TRACER–DP-SCREAM | [TRACER-DP-SCREAM.md](TRACER-DP-SCREAM.md) | A periodic Breeze case using the actual published experiment's forcing |

Example instruction for an agent:

> Read /shared/home/greg/BreezeLab.jl/plans/README.md and /shared/home/greg/BreezeLab.jl/plans/SEA-STARR.md. Implement that plan in your own checkout. Keep a concise status file in the location below, and proceed to the first meaningful simulation as soon as its prerequisites are met.

These files are plans, not evidence that any new case is implemented or running. No new agents or simulations were launched while writing them.

## Shared working agreements

- Keep **ENA-covert**, **ENA-LASSO**, and **ENA-SCREAM** separate experiment identities, run scripts, input provenance, and output directories. Shared constructors/forcing code are encouraged. Never substitute Covert inputs for a missing LASSO or SCREAM bundle. ENA-SCREAM is a separate future assignment pending its reference experiment; it is not the periodic TRACER assignment here.
- Preserve the existing Covert entry point until a coordinated rename provides compatibility. Do not turn the generic ENA script into another protocol by changing defaults.
- Use your own checkout under `/shared/home/greg/breezelab-work/<case>/`, based on current origin/main; inspect repository instructions first. The existing `/shared/home/greg/BreezeLab.jl` may be stale or used by another agent. Do not reset it. Commit case-specific work on a distinct branch; direct integration to main is authorized, but serialize integration and reconcile the latest main before pushing. Never force-push or overwrite another agent's changes.
- **Do not edit or switch the source checkout at `/shared/home/greg/breezelab-runs/docs-20261007-8dc824d/source`.** Jobs 165–169 share it and were running when this index was written. Job 165 builds and publishes GPU documentation; leave its scripts and publication workflow alone. Do not cancel or duplicate existing simulations.
- Use **H100 or A100**, whichever is available. Slurm is under `/opt/slurm/bin/`. Inspect both allocations and physical GPU processes: some other workloads do not declare GPU resources in Slurm. Follow the existing UUID-pinning launch pattern, accounting for other campaigns. A scheduler allocation is not proof a particular GPU is free. No heavy runs on the login node.
- Use Julia for all figures and analysis. Prefer NumericalEarth data metadata, field time series, regridding, boundaries, and nesting where applicable. Extend public APIs by Julia method extension or package extensions when appropriate; do not copy whole subsystems or use method piracy.
- Model constructors belong in `src/` and return a simulation without running it. Readable scripts in `cases/` should construct, run, and analyze like Breeze examples. Avoid a Python-style main guard for test reuse. Put ingestion/download tooling in `data_wrangling/`, meaningful CPU tests in `test/` under `Pkg.test()`, and Slurm helpers in `execution/`. Add Literate/Documenter examples that can be generated manually on GPU; do not pretend GPU CI exists.
- Data/cache/output directories stay outside Git. Record input IDs, URLs, checksums, units, timestamps, code commits and environment manifests. Do not put credentials into code, logs, Git, or these plans. Proprietary Aeolus LESbrary material may inform implementation but must not be explicitly cited or reproduced in shared documentation.
- Do not contact people by Slack/email without Greg's explicit authorization. Record a precise question only when it actually blocks a required choice. Continue independent implementation work.

## Progress and completion

Each agent owns `status/<case>.md` in this folder (for example `status/SEA-STARR.md`). Record only: commit/branch, implemented behavior, exact launch command, job/output path, checks and results, and any concrete blocker. Do not edit the other agents' status files.

Use a bounded progression: reference-value/CPU tests, one short GPU integration, then the full baseline. Repeat checks only for changes or unresolved failures. A configuration tuple is not a runnable case; a short pilot is not a completed MIP. Clearly label exploratory departures and distinguish simulation completion from scientific agreement.

For cross-case shared changes (exports, dependencies, docs navigation, shared forcing utilities), announce ownership in your status file and coordinate integration rather than making conflicting broad refactors. Keep cases independently runnable.

## Existing cluster evidence

- Current campaign: `/shared/home/greg/breezelab-runs/docs-20261007-8dc824d` (five-run manifest and logs).
- Earlier ENA pilots and staged resources: `/shared/home/greg/breezelab-runs/20261006`.
- Known running source commit: `3330f3cf98a62766bf139ce04276745d00cadd48`; it imports CUDA to register `GPU()` and fixes CLI initial `--dt` handling. New work should inspect current main instead of assuming this remains the latest commit.
- Repository starting references: `docs/production-readiness.md`, `docs/cases/ena.md`, and the case input manifests. The readiness document contains older completed ENA checklist items; treat logs/status as the live evidence.
