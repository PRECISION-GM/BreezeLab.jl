# Analysis scripts

Julia scripts that read BreezeLab run output (`*_profiles.jld2`, `*_timeseries.jld2`,
`*_slices.jld2`, `provenance.toml`) and the staged ARM / SAM reference data, and write
figures and TOML summaries under `runs/` (not versioned). Run them on a CPU node
(`sbatch execution/cpu_analysis.sbatch analysis/<script>.jl ...`), never on the login node.

| script | purpose |
| --- | --- |
| `compare_ena_observations.jl` | one run vs ARM MWRRET LWP, VDIS rain, ARSCL/ceilometer cloud boundaries over the run's UTC window |
| `compare_sam_reference.jl` | one LASSO run vs the member's SAM `samstat` statistics (and MWRRET) |
| `ena_covert_analysis.jl` | several covert runs: time-series overlays vs observations, hourly-mean profiles |
| `covert_paper_comparison.jl` | covert runs vs the numbers stated by Covert et al. (2022) over 09–12 UTC, with rain sections |
| `ena_aerosol_audit.jl` | aerosol/droplet-number profiles and conservation of P3-aer runs |
| `check_lasso_bundle.jl` | validate a staged LASSO bundle and build/step tiny CPU cases from it |
| `dump_netcdf_header.jl` | dimensions/variables/attributes of NetCDF files |
| `plot_results.jl`, `animate_*.jl`, `microphysics_snapshot.jl` | inherited plotting/animation tools (separate `analysis/Project.toml`) |

## Pitfall: compute before you plot (Makie state stall)

On the campaign environment (Julia 1.12.7, CairoMakie from the pinned Manifest) a plain
loop such as `ks = findall(≥(threshold), view(data, i, :))` over a few hundred columns,
which takes milliseconds on its own, **stalled for more than an hour** (97 % CPU, no
progress, every Julia thread idle in the scheduler) once `Figure`/`Axis` objects had been
created earlier in the same session — reproduced three times on the CPU node and bisected
to the Makie objects (`runs/scratch/bisect_columns.jl`, 8 October 2026). The scripts here
therefore compute every statistic first and create figures last; keep that order in new
scripts, and if a script hangs after `Figure()` was called, that is the first suspect.
