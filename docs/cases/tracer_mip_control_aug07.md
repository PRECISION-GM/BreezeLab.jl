# TRACER-MIP August 7 2022 outer-domain control: first 24-h attempt (job 233)

Status record of the first full-length run of the regional outer domain, its failure at 13 simulated hours,
the diagnosis, and the comparison of its valid window with the TRACER AMF1 observations at La Porte.
Figures are under `runs/outer_era5_control_aug07_job233/{analysis,diagnostics,animations}/` of the agent
checkout (`/shared/home/greg/breezelab-work/tracer-mip/`). This is not a MIP result.

## What ran

| Item | Value |
| --- | --- |
| Job | Slurm 233, A100-SXM4-80GB (GPU-bdbd7e9e), commit 3e7cec0, started 2026-10-08 04:53 UTC |
| Domain | Grid-1: 750 × 750 × 94, 2 km cells at 29.47 N on a lat-lon grid (22.73–36.22 N, 102.83–87.33 W), ACPC levels, top 22181 m, ETOPO2022 terrain (terrain-following) |
| Boundaries | ERA5 hourly pressure-level parent (NumericalEarth `nested_atmosphere_model`), open lateral BCs + Davies relaxation over 5 cells (rate 1/300 s⁻¹), ρw sponge over the top Lz/4 (≈ 5.5 km) and the lateral zone, rate 1/5 s⁻¹ |
| Physics | Breeze compressible split-explicit; P3 two-moment microphysics with the Tier-1 prescribed two-mode aerosol (fixed, radiatively inactive); RRTMGP all-sky every 60 s; `SlabLand` from ERA5 skin temperature and ERA5-Land soil moisture, sea cells pinned to ERA5 skin temperature; **no turbulence closure**; WENO(5) advection (bounded for water species); fixed Δt = 3 s |
| Outputs | hourly 3-D state (19 fields, 1.4 GB each) and process-rate accumulators; 10-min surface fields and 2-km slices; 10-min 3-D state on the inner 250 km box + 20 km halo (for the offline nest) |
| Throughput | 4.3 s wall per 3 s step (≈ 86 min per simulated hour); 39.3–40.3 GB resident |
| End | NaN in the slab-land temperature at iteration 15700, t = 47100 s (13.08 h, 19:05 UTC); 50 GB of output through 19:00 UTC; no `COMPLETE` marker |

## Failure diagnosis (`analysis/tracer_mip_blowup_diagnosis.jl`; `analysis/blowup_diagnosis.{toml,png}`, `blowup_map.png`)

Hourly 3-D extremes (location of the maximum, whether in the 10-cell lateral rim, in the upper sponge (z > 16.6 km), over sea):

| hour (UTC) | max |w| m/s | where | min T (K) at z | max |u| m/s | where | cells with |w| > 10 (fraction over sea / in rim) | max vertical CFL |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 0 (06) | 1.6 | 379 m, interior | 201 at 16.6 km | 25.6 | 11.8 km, rim (NE) | 0 | 0.09 |
| 1 (07) | 6.2 | 448 m, rim | 201 | 59.4 | 82 m, west rim (i = 2) | 0 | 0.27 |
| 2–4 (08–10) | 8–10 | 21.4 km, SE corner (rim + sponge) | 201 | 59–61 | 30–82 m, west rim | ≤ 1 | 0.27 |
| 5 (11) | 19.3 | 20.2 km, SE corner | 202 | 64.5 | 22.0 km, SE corner | 16 (0.5 in rim) | 0.27 |
| 6 (12) | 32.1 | 10.0 km, interior Gulf | 202 at 17.2 km | 63.9 | 22.0 km, SE corner | 670 (0.99 sea, 0.01 rim) | 0.32 |
| 7 (13) | 37.3 | 13.0 km, interior Gulf | 190.6 at 15.7 km | 56.0 | 22.0 km, SE corner | 1056 (1.00 sea) | 0.37 |
| 8 (14) | 39.2 | 12.4 km, Gulf | 187.8 at 16.6 km | 88.7 | 22.0 km, SE corner | 1355 | 0.39 |
| 10 (16) | 45.1 | 13.3 km, land | 189.1 at 16.6 km | 90.0 | 22.0 km, SE corner | 882 (0.81 sea) | 0.45 |
| 12 (18) | 45.7 | 13.0 km, Gulf | 184.9 at 16.6 km | 105.6 | 22.0 km, SE corner | 2242 (0.86 sea) | 0.46 |
| 13 (19) | 78.3 | 15.1 km, cell NE of Houston | 160.1 at 16.0 km, same cell | 98.2 | 22.0 km, SE corner | 3887 (0.66 sea) | 0.78 |

Reading:

1. **Two independent pathologies.** (a) From the first hours the lateral relaxation zone carries winds of 60 m/s near the
   surface on the inflow (west) boundary, and from hour 5 the **model-top/wall corner** (i = 750, j = 1, z = 22 km, the
   SE corner over the Gulf) holds the domain's wind maximum (88–107 m/s from hour 8) and |w| up to 19 m/s at 20–21 km.
   This is the lid/wall corner mode NumericalEarth's nesting code already documents (the ρw sponge covers it, the
   horizontal momentum there is only Davies-relaxed at 1/300 s⁻¹). (b) From hour 6 (12 UTC, 07 CDT, i.e. **before land
   convection**) hundreds of isolated grid cells over the warm Gulf (SST pinned at ≈ 303 K) develop updrafts of 30–45 m/s
   at 10–13 km, with overshooting tops cooling the layer just below the sponge base to 185–190 K; `blowup_map.png` shows
   these as single-cell "popcorn" storms. The final event is one such cell NE of Houston: w = 78 m/s and T = 160 K at
   15–16 km at 19 UTC, after which the surface coupling received NaNs.
2. **Why the no-closure configuration is at fault.** Without vertical mixing the surface fluxes accumulate in the 50 m
   lowest layer: the surface-wind check (`diagnostics/surface_wind_check.png`) shows the 25 m wind at 0.6 of ERA5's 10 m
   wind with the level above untouched, and the site column saturates (RH 100 %, model "cloud base" at the lowest level
   from 14 UTC). Over the Gulf this builds a thin super-moist layer under a dry profile; grid-scale convection then
   releases it explosively with no horizontal or vertical diffusion acting on the 2 km cells (WENO(5) alone).
3. **Time step.** Horizontal CFL stayed ≤ 0.16 (|u| 107 m/s · 3 s / 2 km); the vertical advective CFL |w|Δt/Δz rose from
   0.27 to 0.46 by hour 12 and 0.78 at hour 13 (Δz = 300 m aloft) — close to, but not beyond, the explicit limit. The
   fixed 3 s step did not cause the blow-up; it merely stopped bounding it (an adaptive CFL step would have slowed down
   instead of producing 78 m/s).
4. **Not implicated**: the terrain (maxima are over sea or at 15 km), the lateral relaxation as the *source* of the deep
   cells (≤ 3 % of |w| > 10 cells in the rim), the microphysics limiter (droplet number reached 4·10¹¹ kg⁻¹ only inside the
   blown-up cells).

## Is the TKE closure alone enough? (job 339, the replacement control)

`TKEBasedTurbulenceClosure` adds vertical eddy diffusion (prognostic TKE, vertically implicit): it removes pathology 2
(surface-layer accumulation: the CPU smokes and the surface-wind reasoning), so the explosive marine convection should
be much weaker. It adds no horizontal diffusion and does nothing for the top/wall corner (pathology 1a), so it may not
be sufficient. Job 339 (started 01:22 UTC 2026-10-09 on the same GPU) is watched every 30 min (`monitor.csv`: 2-km and
inner-box 3-D max |w|, min T with heights); it is cancelled and replaced if the inner-box |w| exceeds 15 m/s or the
minimum temperature falls below 195 K before 14 UTC event time, i.e. the job-233 trajectory.

Prepared stabilised configuration (knobs in `cases/tracer_mip.jl`, defaults unchanged): adaptive time step with a CFL
wizard (`TRACER_MIP_CFL=0.5`, max 3 s), lateral relaxation and ρw sponge over 20 cells (40 km) instead of 5, upper
ρw sponge 7 km deep (base ≈ 15.2 km) instead of Lz/4, all with the TKE closure. A horizontal (2-D) Smagorinsky or a
constant horizontal diffusivity alongside the TKE closure would be the LES-consistent companion, but the pinned Breeze
(5264d3c) accepts a single closure (tuples arrived upstream in #1036); this remains the main structural gap.

## Comparison of the valid window (06–19 UTC) with the TRACER AMF1 observations (`analysis/tracer_mip_site_comparison.jl`)

Observations (ARM Live, `data/arm_tracer_20220807/`, checksummed manifest): houmetM1.b1 (2 m T/RH, 10 m wind,
precipitation), houceilM1.b1 (first cloud base), houarsclkazr1kolliasM1.c0 (cloud-base best estimate), houldM1.b1
(disdrometer), housondewnpnM1.b1 (05:29, 11:30, 17:30, 20:30, 23:29 UTC). Model: the 2 km cell containing La Porte
(29.68 N, 95.07 W), lowest level 25 m vs the 2 m / 10 m measurements. `analysis/site_timeseries.png`, `site_soundings.png`,
`site_comparison.toml`.

| hour (UTC) | T model / obs 2 m (K) | wind model 25 m / obs 10 m (m/s) | direction model / obs (°) | RH model / obs (%) |
| --- | --- | --- | --- | --- |
| 06 | 299.9 / 300.3 | 4.4 / 2.1 | 149 / 161 | 87 / 83 |
| 09 | 299.1 / 299.2 | 3.0 / 0.9 | 120 / 119 | 92 / 88 |
| 12 | 298.4 / 298.2 | 3.6 / 0.1 | 91 / 82 | 97 / 94 |
| 14 | 298.0 / 302.6 | 2.0 / 2.5 | 74 / 125 | 100 / 73 |
| 16 | 300.4 / 304.1 | 1.4 / 3.9 | 30 / 108 | 100 / 67 |
| 18 | 303.8 / 304.9 | 3.3 / 5.5 | 88 / 132 | 96 / 61 |
| 19 | 302.4 / 305.0 | 4.2 / 5.9 | 230 / 125 | 85 / 62 |

- **Night (06–12 UTC)**: temperature within 1 K; the model's lowest level is 2–3.5 m/s where the tower measured a calm
  (< 1 m/s) night — the initial ERA5 wind in the lowest cell relaxes only through surface drag.
- **Sea-breeze onset**: observed onshore (ESE, 105–135°) flow is established by ≈ 13:30–14 UTC (08:30–09 CDT) and
  strengthens to 6–8 m/s by 19–21 UTC. The model keeps a weak (1.4–2 m/s) easterly that veers to NE/N (30–18°) at
  16–17 UTC — offshore at La Porte — and only turns to E/SW at 18–19 UTC: **no sea-breeze front reached the site in the
  model before it failed**, and the inner-box animation shows the coastal gradient but no organised onshore surge.
- **Daytime warming lags by ≈ 2 h and 4–5 K** (298.0 vs 302.6 K at 14 UTC); the model lowest layer saturates (RH 100 %,
  fog-like "cloud base" at 25 m from 14 UTC) where the observed RH drops to 60–70 % — the no-closure surface layer.
- **Clouds**: ceilometer first cloud base 700–900 m through the night and morning; the model column has no cloud until
  the spurious saturated surface layer (no 700–900 m cloud base).
- **Rain**: none observed at La Porte on August 7 (tipping bucket 0.0 mm; disdrometer 0); model 0.02 mm (one 10-min
  shower at 18:30 UTC).
- **Soundings**: see `site_soundings.png` (model columns at 06, 12, 18 UTC vs the 05:29, 11:30, 17:30 UTC sondes).
- **ERA5 10 m wind** (whole domain, `diagnostics/surface_wind_check.png`): model 25 m wind 0.58–0.83 of ERA5's 10 m wind
  after the first hour (see `tracer_mip.md`).

## Land water budget: rain-to-land sign error (found 2026-10-10)

Job 233 (and every outer run before commit 4d6f477 / PR #26) coupled the child's rain to the slab land with the wrong
sign: the rain-to-land shim multiplied Breeze's `bottom_precipitation_flux`, which is already positive downward, by −1,
and NumericalEarth adds the coupler's `Jʳⁿ` to the bucket as precipitation. Rain over land therefore *removed* water
from the bucket (≈ 2 × the accumulated rain relative to a correct run, until the bucket emptied). Sea cells are
unaffected (their bucket is reset to saturation every step).

Size in job 233 (`analysis/tracer_mip_rain_sign_impact.jl` on the run's own correctly signed `surface_rain_flux`;
`analysis/rain_sign_impact.toml`), 06:00–19:00 UTC, 383 680 land cells:

| quantity | value |
| --- | --- |
| land-mean accumulated rain | 0.069 mm (0.0014 mm by 12 UTC, 0.014 mm by 16 UTC, 0.03 mm by 18 UTC) |
| land cells with > 1 mm / > 10 mm | 1.02 % (3903 cells) / 0.15 % |
| largest column accumulation | 120.9 mm (bucket capacity 150 kg m⁻²) |
| implied bucket error, land mean / max | 0.14 kg m⁻² (0.1 % of capacity) / up to emptying the bucket |
| bucket saturation change, rainy (> 1 mm) vs dry land cells | −0.046 vs −0.015 (rainy cells dried faster instead of wetting) |

Assessment: negligible for the domain-mean land state and for the La Porte comparison (0.02 mm of model rain at the
site), but locally wrong in the ~1 % of land cells under the afternoon convection after ≈ 16 UTC, where the surface
became drier (less evaporation, warmer skin) than it should have. It cannot have caused the 13 h failure: the runaway
cells were 99 % over the Gulf from 12 UTC, where land rain was still 0.0014 mm, and sea buckets are pinned.

## Departures from the reference treatments (this run)

Lat-lon grid (not polar stereographic); 94 cells from the 95 ACPC scalar levels; P3 bulk two-moment (not bin);
slab land with bucket hydrology, no vegetation/urban physics, sea pinned to ERA5 skin temperature with land roughness;
hourly ERA5 boundaries; **no turbulence/PBL closure (the defect)**; **rain removed from instead of added to the land
bucket (sign error, see above)**; aerosol activation spectrum scaled by the prescribed
profile; process rates re-evaluated per step rather than P3's internal budget; inner nest not run.
