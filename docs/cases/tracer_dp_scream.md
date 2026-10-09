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

## Output cadence

Animation fields (xy maps at the slice height, the xz section, LWP and rain maps) are saved every 5 min at 3000 m by default (`TRACER_DP_SCREAM_SLICE_MINUTES`, `TRACER_DP_SCREAM_SLICE_HEIGHT`); profiles are 30-min means and time series are 60 s. The case script writes the slice-file storage estimate (saves, bytes per save, GB) into `provenance.toml`; job 202 (256² × 160, 30-min cadence) measured 1.175 MB per save.

## 200 m baseline (job 202, pre-fix): budgets and comparison

Job 202 (256² × 200 m, 2022-08-01 00 UTC – 08-15 00 UTC, A100, 1.21 days wall, branch commit bea5e8a)
and the 100 m run job 416 (512² × 100 m, H100, main a627610) are **pre-fix** runs: their energy forcing and
prescribed surface energy flux carried the `(cᵖᵛ − cᵖᵈ) T ×` vapor-rate cross term that belongs to the
static-energy formulation only, a spurious heating in the potential-temperature formulation of ≈ +13 W m⁻²
over the 14 days (surface ≈ 9, large-scale ≈ 4; branch `vapor-cross-term-formulation`). They form a matched
pair for the 200 m → 100 m comparison; fixed reruns follow.

**Budgets** (`analysis/tracer_dp_scream_budget.jl`, `src/tracer_dp_scream/budget.jl`). The column-static-energy
budget written by the run (27 % residual) omitted the condensate removed by precipitation (ℒP = 125 of the
146 MJ m⁻² residual) and is replaced by a moist-enthalpy budget: 14-day residual 15 MJ m⁻² (1.8 % of the
sources, ≈ 12 W m⁻²), ≤ 0.9 % on dry days and 3–14 % (≈ 0.1–0.2 ℒP, a gain) on heavy-rain days — the
potential-temperature formulation's sedimentation/phase-change energy coupling in Breeze (open, upstream).
Water closes to 0.6 %.

**Comparison** (`analysis/tracer_dp_scream_compare.jl`; hourly, 1–15 August UTC unless noted; VARANAL =
the MRMS-constrained forcing analysis, domain mean; ARM = AMF1 site):

| | LES 200 m | DP-SCREAM 3 km | DP-SCREAM 0.5 km (5–15 Aug) | VARANAL | ARM site |
| --- | --- | --- | --- | --- | --- |
| precipitation (mm day⁻¹) | 3.58 | 3.22 | 4.23 | 4.12 | 7.57 (gauge) |
| hourly correlation with VARANAL | 0.65 | 0.26 | 0.14 | 1 | 0.43 |
| daily-total correlation with VARANAL | 0.97 | 0.84 | 0.76 | 1 | 0.77 |
| mean daily peak hour (CDT) / mean abs. difference from VARANAL | 15.4 / 1.9 h | 18.2 / 3.3 h | 16.3 / 2.7 h | 15.8 | 22.4 / 4.5 h |
| LWP (g m⁻²) | 40 | 19 | 11 | 64 | 31 (MWR) |
| IWP (g m⁻², 5–15 Aug) | 200 | 200 | 200 | — | — |
| precipitable water (kg m⁻²) | 56.1 | 53.4 | — | 51.6 | 52.7 (MWR) |
| T − sondes 0–2 / 2–6 / 6–12 km (K) | +1.7 / −0.1 / −0.9 | +3.1 / +1.8 / +0.3 | | −0.1 / −0.2 / −0.5 | |
| RH − sondes 0–2 / 2–6 km (%) | −7 / +10 | −8 / −7 | | | |
| cloud fraction <3 / 3–8 / >8 km | 0.020 / 0.067 / 0.111 | 0.018 / 0.024 / 0.100 | | | 0.118 / 0.080 / 0.119 (ARSCL occurrence) |

The convective episodes and their timing follow the forcing in all three models; the LES tracks the
MRMS-constrained precipitation more closely than DP-SCREAM 3 km (whose daily peak is ≈ 2.5 h late) and
has more mid-level (congestus) cloud, closer to ARSCL, while both underrepresent cloud below 3 km relative
to ARSCL (which also counts precipitation at a point). The LES is too moist at 2–6 km (+10 % RH, +4.5 kg m⁻²
PW) and warm in the lowest 2 km (+1.7 K; the surface part of the cross term heats the boundary layer —
re-evaluate with the fixed rerun). Without wind nudging (as the archived DP-SCREAM code) the LES mean wind
stays near the 1 August profile: upper-tropospheric easterlies reach −14 m s⁻¹ at 12 km against ≈ −5 m s⁻¹
in VARANAL (RMS 4.2 m s⁻¹ in u over 0–12 km); the archive has no winds to compare. Agreement with DP-SCREAM
is not observational validation, and the point observations are not domain means.

**Fixed reruns** (main 0bb3406, after the formulation-aware cross-term fix, PR #21), queued 2026-10-09 on
the A100 partition with `execution/tracer_dp_scream_autogpu.sbatch` (`--no-requeue`, 6-h checkpoints,
slices every 5 min at 3 km): job 514 `baseline_fixed_aug01_15` (256² × 200 m) and job 515
`baseline512_fixed_aug01_15` (512² × 100 m, starts after 514). Their budgets are computed with
`analysis/tracer_dp_scream_budget.jl` (the case script on this branch writes the moist-enthalpy budget).
