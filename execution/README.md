# Cluster launchers

Case physics and settings live in `cases/`; this directory only handles submission.
Run from the repository root with Julia available and the environment instantiated:

```sh
mkdir -p output
# H100 on wpcluster:
sbatch --partition=gpu-prod execution/submit_gpu.sbatch cases/ena_covert.jl
# Or one A100 on a shared node:
sbatch --partition=gpu-p4de-2c execution/submit_gpu.sbatch cases/ena_covert.jl
```

Override wall time, CPU and memory requests for the measured case. The wrapper runs
the Julia script and arguments passed to it. For the general CLI, pass
`cases/cli/run_case.jl --arch gpu --protocol ... --data ...`.

`submit_ena_lasso.sbatch` launches the official LASSO-ENA member with the wpcluster
UUID-pinning procedure (shared GPU lock directory, physical idleness check, per-CPU depot):
`sbatch --partition=<p> --nodelist=<node> execution/submit_ena_lasso.sbatch <GPU UUID> <experiment>`
runs `cases/ena_lasso.jl` into `runs/<experiment>_job<id>/`; extra arguments are passed to
`cases/cli/run_case.jl` instead. `wpcluster_env.sh` sets the same Julia/depot environment
for the login node.

`production_runs.sh` is the inherited ENA parameter-sweep launcher. It submits all
selected members without dependencies; `RUNS=p3_n75` selects just one. Its 72-hour
default is a historical allocation ceiling, not a runtime estimate. Prefer the
single-case command above for a first run and set `TIME` explicitly for a sweep.
The H100/A100 policy applies to new work; no T4 launcher is selected by default.
