# Tests

```sh
julia data_wrangling/fetch_covert_inputs.jl
julia --project -e 'using Pkg; Pkg.test(; allow_reresolve=false)'
```

`Pkg.test()` runs input/forcing/microphysics regressions, protocol checks, offline
ARM downloader checks and ENA CPU execution/output tests for 1M, P3-N75 and P3-aer2.
CI fetches the pinned public Covert inputs first, so the ENA checks run on every push
and PR. Without those inputs, the data-dependent checks are explicitly skipped;
tests never silently download data or require ARM credentials.

On an allocated H100 or A100, the same suite can also exercise the three GPU cases:

```sh
julia --project -e 'using Pkg; Pkg.test(; allow_reresolve=false, test_args=["gpu"])'
```

GitHub-hosted CPU CI does not claim GPU coverage. `ena_execution.jl` checks the
actual case builder, completed time/steps, finite fields, expected diagnostic fields,
consistent saved iterations, and final output time. Temporary outputs are cleaned up.

[diagnostics/](diagnostics/) contains historical investigation scripts: tendency
scans, stage traces, rain-number/runaway probes, mass-budget probes, limiter fuzzing,
and longer staged smokes. They are retained for diagnosis, not counted as passing
regression tests merely because they print results or exit successfully. Their shell
wrappers run from the repository root; override the partition for an available A100.
Promote a probe into `Pkg.test()` when it has a bounded runtime and a justified
assertion; the repeatable case execution checks have already been promoted.
