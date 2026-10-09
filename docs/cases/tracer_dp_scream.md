# TRACER – DP-SCREAM-matched periodic case

`cases/tracer_dp_scream.jl` runs a doubly periodic Breeze LES with the forcing of the
published DP-SCREAM TRACER experiment (Oware et al. 2025, JGR Atmospheres,
doi 10.1029/2025JD044113). It is a separate experiment identity from ENA (Covert, LASSO) and
from regional TRACER-MIP; it is not the ENA-SCREAM assignment. Constructor:
`tracer_dp_scream(; ...)` in `src/tracer_dp_scream.jl`; inputs pinned in
`cases/tracer_dp_scream/inputs.toml` (fetch with `data_wrangling/fetch_manifest.jl`).

## Sources (pinned 2026-10-07)

| Item | Pin |
| --- | --- |
| Paper | doi 10.1029/2025JD044113 (JGR Atmos. 130(20), 2025-10-18); accepted manuscript OSTI 3016077 |
| Reference outputs | Zenodo 10.5281/zenodo.15271730 (v2 of concept 15133178; v1 = 15133179), CC-BY-4.0 |
| Run script | E3SM-Project/scmlib `DPxx_SCREAM_SCRIPTS/run_dpxx_scream_TRACER.csh` @ 2dc3f1073a5d03b5f32617e79d68fb20a63cfb43 |
| IOP forcing | `TRACER_iopfile_4scam.nc` from the public E3SM inputdata server (sha256 1a95b3b9…e70e), ARM VARANAL v1 `hou60v1varanaecmwfX12.c1`, doi 10.5439/1860369 |
| Model code of the archived runs | E3SM `git_version` 1d551ea2b0 (2024-03-12), recorded in every archived NetCDF |

## What DP-SCREAM actually did (from the script, the IOP file and the E3SM source at 1d551ea2b0)

- 200 km × 200 km, 20 × 20 spectral elements (3.33 km), 128 hybrid levels to ~3 hPa, physics
  step 100 s, dynamics 8.33 s, RRTMGP every 3 steps (300 s), SHOC + P3, prescribed aerosol.
  The 0.5 km run used 15 s / 1.25 s and started cold on 2022-08-05 00 UTC (Table 1 of the
  paper; the archive's first records show the spin-up). The 3 km August run started on
  2022-08-01 00 UTC; its days 5–15 are the paper's comparison window.
- Forcing: hourly VARANAL profiles on 40 pressure levels. EAMxx forms
  `divT3d = divT + vertdivT` and `divq3d = divq + vertdivq` when the vertical terms exist and
  adds only those to `T` and `qv` (`advance_iop_forcing`); `iop_dosubsidence = false`, so the
  omega profile is unused. Data are read piecewise-constant over each hourly interval.
- `iop_srf_prop = true`: sensible heat flux, evaporation (`lhflx/Lv`) and the radiative
  surface temperature `Tg` come from the file; the momentum stress comes from the land/ocean
  coupler (not reproducible from public inputs).
- Planar HOMME sets `fcor = 0`; `iop_coriolis = false`: no Coriolis force.
- Wind nudging: the current script sets `iop_nudge_uv = true` and the paper states the
  winds are nudged, but the archived code revision predates the DP-EAMxx nudging port
  (E3SM 3f7eee0053, 2024-04-16; the 1d551ea2b0 namelist has no `iop_nudge_*` option). The
  archived runs therefore ran **without** wind nudging; `wind_nudging_timescale = nothing`
  is the default and `10800` reproduces the script's intent as a sensitivity.
- Archived outputs: 30-minute domain means (PRECL, TMQ, TGCLDLWP, TGCLDIWP; TOT_CLOUD_FRAC,
  Z3, T, RELHUM, W_SEC, Q on 128 levels), time axis labelled in local time (CDT, UTC−5).
  No 3-D fields, no winds, no surface fluxes; PRECL has no units attribute (mm/day inferred).

## Breeze implementation

`iop_sam_inputs` maps the IOP window onto SAM-format records, which `tracer_dp_scream` reads
before building its model: initial
sounding from the start record (θ from T, specific humidity as mass fraction), `tls/qls =`
the 3-D tendencies, `wls = 0`, no geostrophic columns, `sfc = (Tg, shflx, lhflx)`, every
record duplicated 1 s before the next so the linear interpolation holds the hourly value.
The large-scale transport therefore enters exactly once (`vertical_advection = nothing`,
`geostrophic = false`, `coriolis = nothing`, `upper_boundary_relaxation = false`).

Intentional differences (recorded in the run label and provenance): Smagorinsky–Lilly LES at
200 m on 51.2 km (160 levels, 50 m near the surface, top 22 km, sponge above ~16.5 km)
versus SHOC at 3.33 km on 200 km; P3 with a prescribed droplet number (200 cm⁻³, assumption)
versus prognostic droplet number with prescribed CCN; constant-coefficient bulk drag for a
0.1 m roughness length versus the coupler stress; IOP-implied albedo 0.15 and emissivity 0.98;
anelastic dynamics; no aerosol radiative effects.

## Checks

`Pkg.test()` covers the IOP reader against reference values, the adapter (3-D transport once,
hold expansion, basis), the archive reader (local time → UTC) and a 4-step CPU case with
RRTMGP, outputs and a checkpoint pickup, all on excerpts in `test/fixtures/`. The GPU script
writes `budget.toml` with column water and (approximate) static-energy closure residuals and
plots water paths, precipitation and time–height cloud fraction against the archived 3 km and
0.5 km runs. Agreement with DP-SCREAM is not observational validation; the ARM products used
in the paper (ARMBECLDRAD, VARANAL precipitation, MRMS) are not staged here.
