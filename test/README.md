# Tests

```sh
julia data_wrangling/fetch_covert_inputs.jl
julia --project -e 'using Pkg; Pkg.test(; allow_reresolve=false)'
```

`Pkg.test()` runs input/forcing/microphysics regressions, protocol checks, the LASSO-ENA
bundle parser/validator and adapter checks on a synthetic SAM-style bundle
(`fixtures/lasso_bundle`: analytic forcing profiles, time-interpolated sounding, unit and
tendency conversions, forcing assembly, provenance), ARM station-series readers on
generated NetCDF files, offline ARM downloader checks and ENA CPU execution/output tests
for 1M, P3-N75 and P3-aer2.
CI fetches the pinned public Covert inputs first, so the ENA checks run on every push
and PR. Without those inputs, the data-dependent checks are explicitly skipped;
tests never silently download data or require ARM credentials.

On an allocated H100 or A100, the same suite can also exercise the three GPU cases:

```sh
julia --project -e 'using Pkg; Pkg.test(; allow_reresolve=false, test_args=["gpu"])'
```

GitHub-hosted CPU CI does not claim GPU coverage. `ena_execution.jl` checks the
`ena_covert` constructor (including its initially unadvanced clock), completed time/steps, finite fields, expected diagnostic fields,
consistent saved iterations, and final output time. Temporary outputs are cleaned up.
Tests do not include executable case scripts. The full example and its analysis are
exercised by the [manual Documenter/Literate build](../docs/README.md), with an optional
CPU smoke mode.

[diagnostics/](diagnostics/) contains two investigation scripts: the forcing-free rain-shaft
mass budget (`mass_budget_probe.jl`) and the bounds-preserving WENO limiter fuzz
(`limiter_fuzz_gpu.jl`). They are retained for diagnosis, not counted as passing regression
tests. The earlier P3 probes drove the switch-based ENA builder that each case's own
constructor replaced; they remain in the git history (commit cfb6c98).
