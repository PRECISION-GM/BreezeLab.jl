# Cases

| File | Status |
| --- | --- |
| [eastern_north_atlantic.jl](eastern_north_atlantic.jl) | Runnable public Covert ENA case; edit the readable settings, then run with `julia --project cases/eastern_north_atlantic.jl` |
| [ena_lasso.jl](ena_lasso.jl) | Official LASSO-ENA member (`20170718era5d25x100_sbmwrm-aer2-flxsst` by default) from its staged `samin` bundle; fails fast with staging instructions when the bundle is absent. Requires the ARM archive (not yet staged online as of 7 October 2026). |
| [sea_starr.jl](sea_starr.jl) | Runnable SEA STARR CTRL/N100/N030 from the official DEPHY drivers via the exported `sea_starr(; member, ...)` constructor; exploratory CTRL with labelled departures (see [docs/cases/sea_starr.md](../docs/cases/sea_starr.md)) |
| [tracer_dp_scream.jl](tracer_dp_scream.jl) | Runnable periodic LES driven by the DP-SCREAM TRACER IOP forcing (`tracer_dp_scream` constructor); inputs pinned in `tracer_dp_scream/inputs.toml`; see [docs/cases/tracer_dp_scream.md](../docs/cases/tracer_dp_scream.md) |
| [tracer_mip.jl](tracer_mip.jl) | Regional nested coastal control under implementation (branch `tracer-mip`): protocol metadata, grids and aerosol profiles are implemented and tested; see [docs/cases/tracer_mip.md](../docs/cases/tracer_mip.md) |

The ENA example runs from top to bottom: call the exported
`eastern_north_atlantic(; arch, microphysics, ...)` constructor, save provenance,
run the simulation, and plot cloud water and rain. The constructor lives in `src/`
and does not advance the simulation. Tests call it with a reduced grid and duration;
they do not include the example.

[Documenter builds](../docs/README.md) execute the example with Literate and include
its figure in the manual. Documentation is generated manually on an H100 or A100;
an explicit CPU smoke build checks the example and plotting workflow locally.
For command-line parameter sweeps and official LASSO input handling, use
[cli/run_case.jl](cli/run_case.jl). Slurm launchers live in [execution/](../execution/).
MIP input manifests stay in the corresponding subdirectories; download tools live in
[data_wrangling/](../data_wrangling/).
