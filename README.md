# BreezeLab.jl

BreezeLab develops observation-comparable LES cases and differentiable cloud-physics
workflows with [Breeze.jl](https://github.com/NumericalEarth/Breeze.jl) for PRECISION/GENESIS.
The first implementation ports the ENA machinery from
[BreezyLASSO.jl](https://github.com/glwagner/BreezyLASSO.jl); see [NOTICE](NOTICE).

## ENA: one site, two experiment protocols

| Protocol | Inputs and forcing | Current status |
| --- | --- | --- |
| `covert_public_bin` | Public Covert et al. (2022) bin-repository inputs; prescribed surface fluxes, simple longwave radiation, no wind nudging | Runnable development benchmark. Its domain/duration differ from the published experiment, and its missing vertical grid is reconstructed. |
| `lasso_ena_official` | Selected ARM `samin` bundle and grid; SST-based bulk fluxes, RRTMGP longwave/shortwave, mean-wind nudging | Breeze implementation of the LASSO forcing pathway; requires official inputs and validation against SAM. |

Protocol selection is explicit. A missing LASSO bundle never falls back to the Covert
case. Grid, physics, and duration overrides are recorded in each run's provenance.
**Neither pathway is presently a validated reproduction of the original SAM simulations.**
See [ENA protocols and the reproduction checklist](docs/cases/ena.md) for the scientific
differences and the work needed to establish reproduction.

## Quick start

Use Julia 1.12 and Git. The checked-in `Manifest.toml` pins the environment, including
Breeze and Oceananigans revisions required by the imported implementation.

```sh
julia --project -e 'using Pkg; Pkg.instantiate()'
julia --project scripts/fetch_inputs.jl
julia --project -e 'using Pkg; Pkg.test(; allow_reresolve=false)'
julia --project scripts/smoke_test.jl
```

The fetch script downloads pinned public inputs and writes their revisions/checksums
to `data/MANIFEST.txt`. The CPU smoke runs four simulated seconds on an 8×8×24 grid
with one-moment microphysics. It writes `output/cpu_smoke/summary.toml`, provenance,
and JLD2 profiles, time series, and slices. This checks execution and output, not
cloud-physics fidelity. First-time compilation can take several minutes.

Run the public Covert configuration on an NVIDIA GPU:

```sh
julia --project scripts/run_case.jl \
  --protocol covert_public_bin --data data/covert2022_bin \
  --arch gpu --float Float32 --output output/covert_public_bin
```

The default case is 256×256×192 for six simulated hours and is a substantial run.
For Slurm, adapt the partition/resources in `scripts/submit_gpu.sbatch`, create
`output/` before submitting, and pass the same protocol/data arguments:

```sh
mkdir -p output
sbatch scripts/submit_gpu.sbatch \
  --protocol covert_public_bin --data data/covert2022_bin \
  --output output/covert_public_bin
```

For an official LASSO experiment, obtain the selected bundle from the
[ARM LASSO-ENA Bundle Browser](https://lasso-ena.svcs.arm.gov/latest/bundle_browser.html),
then follow [the LASSO instructions](docs/cases/ena.md#lasso-ena).

## Package layout

- `src/ena_protocols.jl`: experiment definitions and input validation.
- `src/case_setup.jl`: shared Breeze model assembly, output, and provenance.
- Other `src/` files: reusable SAM readers, forcing operators, grids, surface fluxes,
  initial conditions, and diagnostics.
- `scripts/`: fetch/run/smoke entry points and inherited GPU/P3 diagnostic tools.
- `analysis/`: Julia plotting and animation tools inherited from BreezyLASSO.
- `test/`: synthetic input/physics tests, protocol checks, and public-input tests
  enabled when the Covert data are downloaded.
- `docs/`: case formulations, reproduction requirements, and historical port notes.

The Julia entry point is `ena_simulation(data_dir; protocol, kwargs...)`.
`ena_protocol_settings(data_dir; protocol)` inspects defaults without allocating a
model. The legacy `lasso_ena_simulation(...; preset=...)` remains for imported scripts, with an explicit `preset` required.

## Next milestones

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
