# Cases

| File | Status |
| --- | --- |
| [eastern_north_atlantic.jl](eastern_north_atlantic.jl) | Runnable public Covert ENA case; edit the readable settings, then run with `julia --project cases/eastern_north_atlantic.jl` |
| [sea_starr.jl](sea_starr.jl) | Specification only; model assembly is not implemented and direct execution fails explicitly |
| [tracer_mip.jl](tracer_mip.jl) | Specification only; model assembly is not implemented and direct execution fails explicitly |

The ENA example defines `eastern_north_atlantic(; arch, microphysics, ...)` and runs
only when executed directly. Tests call the same builder with a reduced grid and duration.
For command-line parameter sweeps and official LASSO input handling, use
[cli/run_case.jl](cli/run_case.jl). Slurm launchers live in [execution/](../execution/).
MIP input manifests stay in the corresponding subdirectories; download tools live in
[data_wrangling/](../data_wrangling/).
