# Cluster launchers

Case physics and settings live in `cases/`; this directory only handles submission.
Run from the repository root with Julia available and the environment instantiated:

```sh
mkdir -p output
# H100 on wpcluster:
sbatch --partition=gpu-prod execution/submit_gpu.sbatch cases/eastern_north_atlantic.jl
# Or one A100 on a shared node:
sbatch --partition=gpu-p4de-2c execution/submit_gpu.sbatch cases/eastern_north_atlantic.jl
```

Override wall time, CPU and memory requests for the measured case. The wrapper runs
the Julia script and arguments passed to it. For the general CLI, pass
`cases/cli/run_case.jl --arch gpu --protocol ... --data ...`.

`production_runs.sh` is the inherited ENA parameter-sweep launcher. It submits all
selected members without dependencies; `RUNS=p3_n75` selects just one. Its 72-hour
default is a historical allocation ceiling, not a runtime estimate. Prefer the
single-case command above for a first run and set `TIME` explicitly for a sweep.
The H100/A100 policy applies to new work; no T4 launcher is selected by default.
