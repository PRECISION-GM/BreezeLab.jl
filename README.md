# BreezeLab.jl

[Documentation](https://precision-gm.github.io/BreezeLab.jl/)

BreezeLab develops observation-comparable LES cases and differentiable cloud-physics
workflows with [Breeze.jl](https://github.com/NumericalEarth/Breeze.jl) for PRECISION/GENESIS.
The first implementation ports the ENA machinery from
[BreezyLASSO.jl](https://github.com/glwagner/BreezyLASSO.jl); see [NOTICE](NOTICE).

## ENA: one site, two experiment protocols

| Protocol | Inputs and forcing | Current status |
| --- | --- | --- |
| `ena_covert` | Public Covert et al. (2022) bin-repository inputs; prescribed surface fluxes, simple longwave radiation, no wind nudging | Runnable development benchmark. Its domain/duration differ from the published experiment, and its missing vertical grid is reconstructed. |
| `ena_lasso` | Selected ARM `samin` bundle and grid; SST-based bulk fluxes, RRTMGP longwave/shortwave, mean-wind nudging | Breeze implementation of the LASSO forcing pathway with a validating bundle parser (`inspect_lasso_bundle`), the `ena_lasso` constructor and `cases/ena_lasso.jl`; the adopted member `20170718era5d25x100_sbmwrm-aer2-flxsst` is listed by ARM but not yet staged online, so no official run exists. See [the LASSO audit](docs/cases/ena_lasso.md). |

Each protocol has its own constructor. A missing LASSO bundle never falls back to the Covert
case. Every configuration value, including the microphysics, is recorded in each run's provenance.
**Neither pathway is presently a validated reproduction of the original SAM simulations.**
See [ENA protocols and the reproduction checklist](docs/cases/ena.md) for the scientific
differences and the work needed to establish reproduction.

## Quick start

Use Julia 1.12 and Git. The checked-in `Manifest.toml` pins the environment, including
Breeze and Oceananigans revisions required by the imported implementation.

```sh
julia --project -e 'using Pkg; Pkg.instantiate()'
julia --project data_wrangling/fetch_covert_inputs.jl
julia --project -e 'using Pkg; Pkg.test(; allow_reresolve=false)'
```

The fetch script downloads pinned public inputs and writes their revisions/checksums
to `data/MANIFEST.txt`. `Pkg.test()` includes four-second CPU simulations on an
8×8×24 grid for one-moment, P3-N75 and aerosol-coupled P3, and verifies their JLD2
outputs. These check execution and output, not cloud-physics fidelity. First-time
compilation can take several minutes. See [test/](test/) for optional GPU checks.

Run the public Covert configuration on an NVIDIA GPU:

```sh
julia --project cases/ena_covert.jl
```

The default case is 256×256×192 for six simulated hours and is a substantial run.
The example builds the case, runs it, and plots saved cloud-water and rain diagnostics.
Edit the settings in [the readable case file](cases/ena_covert.jl), or use
[the CLI](cases/cli/run_case.jl) for parameter sweeps and protocol/input selection.
For Slurm, adapt the partition/resources in `execution/submit_gpu.sbatch`, create
`output/` before submitting, and pass the case script:

```sh
mkdir -p output
sbatch --partition=gpu-prod execution/submit_gpu.sbatch cases/ena_covert.jl
```

For an official LASSO experiment, obtain the selected bundle from the
[ARM LASSO-ENA Bundle Browser](https://lasso-ena.svcs.arm.gov/latest/bundle_browser.html),
stage it with `julia data_wrangling/stage_lasso_bundle.jl ARCHIVE.tar`, then run
`julia --project cases/ena_lasso.jl`; see [the LASSO instructions](docs/cases/ena.md#lasso-ena)
and [the LASSO audit](docs/cases/ena_lasso.md). A missing bundle stops the script with
staging instructions; the Covert inputs are never substituted.

## Package layout

Each case has one constructor that reads its inputs and builds the grid, model and
`Simulation` in its own body, top to bottom:

| Constructor | Case | Source |
| --- | --- | --- |
| `ena_covert` | Public Covert ENA development benchmark | `src/ena/covert.jl` |
| `ena_lasso` | Official LASSO-ENA member from its staged bundle | `src/ena/lasso.jl` (bundle validator in `src/ena/lasso_bundle.jl`) |
| `tracer_dp_scream` | Periodic TRACER case with the DP-SCREAM IOP forcing | `src/tracer_dp_scream/` |
| `sea_starr` | SEA STARR CTRL / N100 / N030 | `src/sea_starr/` |
| `tracer_mip_outer_simulation`, `tracer_mip_inner_simulation` | Regional TRACER-MIP control | `ext/tracer_mip_regional/` |

Microphysics is passed as a Breeze object, for example
`P3Microphysics(; cloud = CloudDroplets(; number_concentration = 75e6))`, and precision follows
`Oceananigans.defaults.FloatType`, so a script sets `Oceananigans.defaults.FloatType = Float32`
once at the top. The aerosol helpers `lasso_aerosol`, `covert_aerosol` and
`kappa_aerosol_activation` return Breeze aerosol objects.

- `src/components/`: the shared physics components: SAM input readers, vertical grids,
  large-scale forcing operators, surface fluxes, simple longwave radiation, sounding
  initial state, microphysics helpers, diagnostics, output writers and provenance.
- `src/observations/`: ARM station-series and DP-SCREAM output readers.
- `cases/`: readable experiment entry points and MIP specifications; the command-line
  runner `cases/cli/run_case.jl` covers both ENA protocols. See [case status](cases/README.md).
- `data_wrangling/`: input acquisition and manifests; see the
  [NumericalEarth observation-extension design](data_wrangling/README.md).
- `execution/`: Slurm launchers and batch submission helpers.
- `analysis/`: Julia plotting and animation tools inherited from BreezyLASSO.
- `test/`: `Pkg.test()` regressions, including case execution and downloads;
  mass-budget and limiter probes live under `test/diagnostics/`.
- `docs/`: a manually generated Documenter/Literate manual that executes the case
  examples, plus formulations, reproduction requirements, and historical port notes.
  See [documentation build instructions](docs/README.md).

## Next milestones

See the [dated production-readiness assessment and concrete TODOs](docs/production-readiness.md)
for ENA, SEA STARR, and TRACER-MIP. Public MIP inputs/references are checksum-pinned:

```sh
julia data_wrangling/fetch_manifest.jl cases/seastarr/inputs.toml data/seastarr_22241697
julia data_wrangling/fetch_manifest.jl cases/tracer_mip/inputs.toml data/tracer_mip_416423a
```

These stage source files; SEA STARR and TRACER-MIP case adapters are not yet implemented.

The periodic **TRACER–DP-SCREAM** case (a Breeze LES driven by the published DP-SCREAM
TRACER forcing; neither TRACER-MIP nor ENA-SCREAM) is runnable:

```sh
julia data_wrangling/fetch_manifest.jl cases/tracer_dp_scream/inputs.toml data/tracer_dp_scream
julia --project cases/tracer_dp_scream.jl
```

See [docs/cases/tracer_dp_scream.md](docs/cases/tracer_dp_scream.md) for the pinned sources,
what DP-SCREAM actually did, and the intentional physics/resolution differences.

1. Reproduce one official LASSO-ENA ensemble member's forcing and grid, then compare
   Breeze diagnostics to its SAM reference and ARM observations.
2. Add SEASTARR and TRACER-MIP as separately documented experiments using shared
   setup and analysis components.
3. Add a small end-to-end differentiable calibration experiment, including a
   gradient check and independent evaluation. Online training, downscaling, and
   regional LES are planned extensions, not implemented features of this bootstrap.

Large inputs, model output, and ARM archives stay outside version control. Every
production result should retain its inputs/checksums, source revisions, configuration,
and protocol overrides alongside the diagnostics.
