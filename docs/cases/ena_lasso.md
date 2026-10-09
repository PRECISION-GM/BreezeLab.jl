# ENA-LASSO: the official LASSO-ENA member in Breeze

Companion to [ENA experiment definitions](ena_protocols.md). This page records, against the LASSO
SAM source, what the `lasso_ena_official` adapter does and does not do, which inputs it
requires explicitly, and every physics difference that matters for comparing Breeze with
the SAM member. It is an audit, not a claim of reproduction: no official bundle had been
run when it was written (7 October 2026).

## 1. Member identity and data status

| Item | Value |
| --- | --- |
| Dataset | LASSO-ENA, DOI [10.5439/2572661](https://doi.org/10.5439/2572661); bundle browser https://lasso-ena.svcs.arm.gov/latest/bundle_browser.html |
| SAM source | https://code.arm.gov/lasso/lasso-ena-codes/lasso_sam_sbm , branch `lasso_ena_noice`, commit `12d02446a2147388dc89d828e6e0553106abea0f` (SAM 6.10.8 with LASSO modifications; `Build`: `ADV_MPDATA`, `SGS_TKE`, `RAD_RRTM`, `MICRO_HUJISBM`) |
| Plan candidate | `20170718era5s1n0d25x100_sbmwrm-aer2-flxsst` — **does not exist**. The authenticated ARM Live listing of the `enalasso` datastream (7 October 2026, 2502 files for 2017-07-17..19) holds 18 `samin` archives: `era5d25x100` and `merra2d25x100` × {`sbmwrm`, `morr`} × {`aer1`, `aer2`, `aer3`}`-flxsst`, and `era5s1n0d25x100`/`merra2s1n0d25x100` × `morr` × `aer1..3-flxsst`. The `s1n0` variants were run with Morrison microphysics only. |
| Adopted member | **`20170718era5d25x100_sbmwrm-aer2-flxsst`** (explicit decision, recorded in `BreezeLab.DEFAULT_LASSO_MEMBER`, the status file and every provenance record). Files: `enalasso_samin_…C1.m0.20170718.000000.tar` (inputs), `enalasso_samstat_…C1.m1.20170718.000000.nc` and `enalasso_sam2d_…C1.m1.20170718.000000.nc` (reference outputs). |
| Token meanings (LASSO-ENA *Simulation Ensembles*, tables 8–10) | `era5`/`merra2`: forcing reanalysis sampled over a 5° box with the default wind-only nudging; `era5s1n0`: the same forcing sampled over a 1° box, default nudging (`n0`); `n1`/`n2`: enhanced T/qv nudging above the inversion / above 2 km; `d25x100`: 25 km wide (256 columns) at 100 m, 8 km top; `sbmwrm`: warm-phase HUJI spectral bin; `morr`: Morrison; `aer1/2/3`: low/medium/high aerosol; `flxsst`: fluxes computed online from the reanalysis SST. So the adopted member differs from the plan's candidate only in the forcing sampling region (5° instead of 1°); nudging, microphysics, aerosol and surface treatment are the same. |
| Availability | Listed but **not staged online**: `armlive/saveData` and `armlive/livedata/saveData` return HTTP 404 (empty body) for the samin tar and even the small samstat file, while a control request for a known-online ENA observation file succeeds with the same credentials. Staging requires an ARM order (bundle browser / Data Discovery); nothing was ordered by the agent. |

## 1b. The staged bundle (ARM order 284988, staged 8 October 2026 05:06 UTC)

`data/lasso/20170718era5d25x100_sbmwrm-aer2-flxsst/` (`stage_lasso_bundle.jl`; archive SHA-256
`b7594c2d…`, 53 739 520 B; `samstat` `f5cb2bed…`, `sam2d` `4a4c6682…` under `data/lasso/archives/`).
Bundle README: SAM v6.10.3 plus LASSO modifications, `model_source_git_hash f83adf58`
(the public branch head audited above is `12d0244`; the member was produced by the earlier
revision), `stagesam_ena 8404ad7`. Facts read from the files:

| item | value |
| --- | --- |
| `day0`, `dt`, `nstop` | 199.0 (2017-07-18 00:00 UTC), 1 s, 86400 (24 h) |
| `nrad`, `tauls` | 30 (30 s), 7200 s (`n0` wind nudging, whole column) |
| `snd` | 3 records (199.0, 199.5, 200.0), 1077 **height** levels to 10 755 m (`p = -999.9`, unused) |
| `lsf` | 37 records (198.75–200.25), 431 height levels to 10 762 m, 9 columns (`read_in_geostrophic_wind = .true.`) |
| `sfc` | 37 samples (198.75–200.25): SST 295.15 → 295.10 K, ERA5 H/LE for reference, `TAU = -999.9` |
| `grd` | 260 scalar levels 12.5 … 8087.5 m (25 m to 6012.5 m, then stretched); Breeze centres differ from SAM levels by ≤ 4 m in the stretched part |
| dimensions | `caseid` is a name → `dimensions = (256, 256, 260)` passed explicitly (`lasso_documented_dimensions`) |
| `&SGS_TKE dosmagor = .false.` | prognostic 1.5-order TKE in SAM; Breeze Smagorinsky–Lilly (recorded) |
| `UNIFORM_SFC_FLX = .false.` | per-column `oceflx` in SAM → adapter `surface_flux_law = :sam_oceflx` (Large & Pond drag exactly; Stanton `0.0327√cdn` and Dalton `0.0346√cdn` fitted to Breeze's `a₀ + a₁U + a₂/U` within 0.4 % over 2–15 m s⁻¹; `umin = 1` m s⁻¹; residual: Breeze's log-profile height shift and Li et al. stability vs SAM's two Monin–Obukhov iterations and the stable-branch Stanton `0.018√cdn`) |
| `fcor = 9.19626e-5` | used directly as the f-plane parameter (SAM does so); it equals `2Ω sin φ` with the sidereal Ω, 0.27 % above SAM's own `4π/86400 sin φ` |
| `compute_reffc = .true.` | SAM diagnoses the liquid effective radius from the SBM spectrum; Breeze prescribes 10 μm (recorded departure) |
| `perturb_type = 5`, `timelargescale = 0`, `ug = vg = 0`, `nxco2 = 1`, `doseasons = .false.` | supported as is |
| aerosol | aer2 modes with the SBM cap emulation (`aerosol_supersaturation_cap = 0.003`): 10.8 % of mode 1 and 91.0 % of mode 2 activatable, 2.384×10⁸ kg⁻¹ (≈ 286 cm⁻³) in total; fully soluble chemistry |

`validate_lasso_bundle` reports 0 problems and 8 recorded warnings (CPU job 241;
`analysis/check_lasso_bundle.jl`), and tiny CPU cases built on the full 260-level column step
finitely.

## 1c. Short GPU integration (plan step 6, first part)

Slurm job 243 on A100 `GPU-c764cc9c` (gpu-p4de-2c-2), source `fa26811`, `cases/cli/run_case.jl`
with the documented command and `--hours 1 --profile_interval 600 --slice_interval 1800`:
3600 steps of 1 s in 44.1 min (0.69 s per step, RRTMGP every 30 steps), finite throughout,
cloud fraction 1 from t = 10 min, LWP 53 → 91 → 61 g m⁻², max |w| ≤ 2.9 m s⁻¹, max
`qᶜˡ` 1.6 g kg⁻¹, rain 0.002 mm d⁻¹, T ∈ [252.3, 295.6] K. Aerosol check (`runs/lasso_aer2_1h_job243/checks/`):
`nᵃ + nᶜˡ` equals the capped total everywhere (2.384×10⁸ kg⁻¹), `nᵃ` constant outside
cloud, in-cloud `nᶜˡ`/CF ≈ 2.1–2.7×10⁸ kg⁻¹ (≈ 250–320 cm⁻³), i.e. the activation still
reaches the cap (the Breeze ratchet, now bounded by the SBM's own maximum).

First-hour comparison with the SAM member (`runs/lasso_aer2_1h_job243/checks/sam/`,
`analysis/compare_sam_reference.jl`; SAM `samstat` is a 2-min domain statistic, Breeze a
60-s domain mean): cloud water path 72.8 (SAM) vs 67.3 g m⁻² (Breeze) over the hour, both
rising to a peak of ≈ 115 / 100 g m⁻² at 00:37 UTC and relaxing; shaded cloud fraction 1.00
vs 0.98; surface precipitation ≈ 0 in both; Breeze's rain water path reaches 0.06 g m⁻²
against SAM's 0.002 (earlier drizzle onset in P3). The hour-1 profiles of θ, qᵛ, u and v
coincide with SAM's to within the line width up to 3 km, and both models carry two cloud
layers (≈ 500–800 m and ≈ 1200–1500 m) with the same shapes; Breeze additionally holds a
thin layer at ≈ 250 m that SAM does not, and its lower-layer `qᶜˡ` peaks at 0.13 vs SAM's
0.17 g kg⁻¹ (`sam_profiles.png`). The member starts with near-saturated air below the
stratocumulus, which is why both models report very low cloud bases in hour 1 (SAM's GCSS
`ZCB/ZCT` are undefined/≈ 0 there). MWRRET has no samples before 03:10 UTC.

The full 24-h run (Slurm job 249, `runs/lasso_aer2_24h_job249/`, `--time=22:00:00`,
≈ 17 h expected) was queued for the next physically idle A100 at 06:25 UTC; its
comparison with `samstat`/`sam2d` and the ARM observations over 03:10–24:00 UTC
(`analysis/compare_sam_reference.jl`, `analysis/compare_ena_observations.jl`) is plan step 7.

## 1d. Full 24-h run and plan step 7: Breeze vs the SAM member and ARM observations

Run: Slurm job 249, A100 `GPU-c764cc9c`, source `fa26811` (adapter as in § 1b/1c with the
per-column `sam_oceflx` law, namelist `fcor`, SS-cap emulation, prescribed 10 μm r_eff),
`cases/cli/run_case.jl` with the documented command: 86 400 steps of 1 s in 17.2 h
(0.72 s per step), finite throughout, `runs/lasso_aer2_24h_job249/`. Reference: the member's
`samstat` (2-min domain statistics) and `sam2d` (5-min 2-D fields), SAM `f83adf58`.
Observations: MWRRET `phys_lwp` (from 03:10 UTC), VDIS rain, ARSCL/ceilometer boundaries.
Analysis `analysis/lasso_step7.jl` (CPU node); figures and `step7.toml` in
`runs/lasso_aer2_24h_job249/checks/step7/`, copies in `breezelab-work/campaign-figures/lasso24_*`.
Cloud boundaries use one definition for both models: the *dominant* contiguous cloudy
layer (largest integrated qᶜˡ; SAM `QCL ≥ 0.01 g kg⁻¹`, Breeze cloudy-cell fraction > 0.05).
The GCSS `ZCB/ZCT` series stored in `samstat` average 0.03/0.10 km, inconsistent with its
own `ZINV` (1.5 km), `ZCTMAX` (1.7 km) and `QCL` profiles, so they are not used.

### Window statistics

| 00–24 UTC | SAM | Breeze | ARM |
| --- | --- | --- | --- |
| cloud water path (g m⁻²) | 84.9 | 129.7 | MWRRET 97.6 (03:10–24 UTC, 1-σ ≈ 13) |
| rain water path (g m⁻²) | 4.77 | 0.78 | — |
| surface precipitation (mm d⁻¹) | 0.263 | 0.031 | VDIS 0.012 |
| cloud fraction | 1.00 (CLDSHD) | 0.93 (LWP > 5 g m⁻²) | ARSCL cloudy time fraction 0.84 |
| dominant-layer base / top (m) | 678 / 1227 | 693 / 1216 | ceilometer 699 / ARSCL top 1021 |
| SAM `ZINV` (m) | 1533 | — | — |
| SHF / LHF (W m⁻²) | 1.3 / 31.7 | not saved | — |
| net LW / SW at surface (W m⁻²) | 16.8 / 214.9 | not saved | — |

| 06–12 UTC | SAM | Breeze | ARM |
| --- | --- | --- | --- |
| CWP (g m⁻²) | 64.9 | 149.7 | 168.3 |
| RWP (g m⁻²) | 1.90 | 0.56 | — |
| precipitation (mm d⁻¹) | 0.187 | 0.028 | 0.037 |
| cloud fraction | 1.00 | 0.98 | 1.00 |
| base / top (m) | 526 / 1006 | 308 / 919 | 667 / 1087 |
| SHF / LHF, LWNS / SWNS (W m⁻²) | 3.0 / 45.3, 18.7 / 248.4 | not saved | — |

| 09–12 UTC | SAM | Breeze | ARM |
| --- | --- | --- | --- |
| CWP (g m⁻²) | 30.1 | 81.4 | 205.8 |
| RWP (g m⁻²) | 0.34 | 0.30 | — |
| precipitation (mm d⁻¹) | 0.024 | 0.015 | 0.066 |
| cloud fraction | 0.99 | 0.96 | 1.00 |
| base / top (m) | 736 / 1069 | 397 / 969 | 670 / 1161 |
| SHF / LHF, LWNS / SWNS (W m⁻²) | 3.8 / 52.5, 25.6 / 439.6 | not saved | — |

### Time evolution (`lasso24_timeseries.png`)

Both models carry the same diurnal cycle: a thick nocturnal deck (SAM CWP peak 151 g m⁻² at
03 UTC, Breeze 244 g m⁻² at 06 UTC), midday thinning under the ≈ 750 W m⁻² net surface
shortwave (SAM 15 g m⁻² at 12–13 UTC with CLDSHD 0.93, Breeze 14 g m⁻² at 14 UTC with cloud
fraction 0.50), and recovery after sunset (SAM 110–135 g m⁻², Breeze 276 g m⁻² at 23 UTC).
Breeze is thicker at night by 60–100 %, within the MWRRET envelope at 06–08 UTC and above
it after 21 UTC (observed 0–130 g m⁻² at 19–24 UTC); its midday cloud-fraction dip is deeper
and 2 h later than SAM's. Drizzle differs in timing and amount: SAM has two drizzle episodes,
06–10 UTC (hourly PREC up to 0.61 mm d⁻¹, RWP 5.7 g m⁻²) and 19–24 UTC (0.8–1.6 mm d⁻¹, RWP
19–27 g m⁻²); Breeze's are 0.04 mm d⁻¹ (RWP ≈ 1 g m⁻²) and 0.05–0.20 mm d⁻¹ (RWP up to 8
g m⁻², 0.247 mm d⁻¹ at the end). Over 24 h SAM precipitates 8× more (0.26 vs 0.03 mm d⁻¹);
the disdrometer saw 0.012 mm d⁻¹ (three brief showers). In both models the evening deck
rises to ≈ 1.7–1.9 km (SAM `ZINV` 1.5 km, `ZCTMAX` 1.7 km; Breeze top 1.6–1.9 km) while
ARSCL tops stay near 1.0–1.3 km — a shared forcing/entrainment issue, not an adapter one.

### Profiles (`lasso24_profiles.png`, 03/06/09/12/18/23 UTC)

θ and qᵛ agree throughout: rms differences below 2 km of 0.4–1.0 K and 0.24–0.45 g kg⁻¹,
the inversion at the same height within one level at every hour, the same two-layer cloud
structure (a shallow fog/stratus layer at 300–900 m until 09 UTC below the deck at
1.3–1.7 km, then a single deck at 0.9–1.1 km at midday and 1.3–1.9 km in the evening).
Winds: u within 1 m s⁻¹ in the boundary layer; v differs by 2–3 m s⁻¹ between 1 and 3 km
after 06 UTC (SAM's free-tropospheric v is more negative), the one systematic dynamical
difference. Cloud water: Breeze's deck is 1.5–2× denser in the evening (max qᶜˡ 0.81 vs
0.43 g kg⁻¹ at 23 UTC, 0.49 vs 0.29 at 06 UTC) and its shallow layer wetter (0.4 vs 0.2
g kg⁻¹ at 03–06 UTC). Rain water: SAM's qʳ peaks 4–20× higher in the evening deck
(2.4×10⁻² vs 6.4×10⁻³ g kg⁻¹ at 23 UTC) while Breeze's early-morning fog layer rains more
(1.5×10⁻³ vs 1.9×10⁻⁴ g kg⁻¹ at 06 UTC). Radiative heating: the cloud-top cooling spikes sit
at the same heights; Breeze's are deeper (−122 vs −79 K d⁻¹ at 06 UTC, −90 vs −65 at 23 UTC,
but −21 vs −26 at 12 UTC) with SAM's clear-air LW cooling and daytime SW warming of the
column reproduced. In-cloud droplet number (Breeze): 215–245 cm⁻³, i.e. at the SBM-cap
ceiling; `samstat` carries no droplet number (it is in `sam3dmicro`, not downloaded).

### Plan views (`lasso24_planviews.png`, 06 and 18 UTC)

06 UTC: SAM CLWP 129 ± 36 g m⁻² (every column > 20), Breeze 243 ± 91 g m⁻² — a smooth closed
deck in SAM against a strongly cellular deck in Breeze with 12 % of columns raining above
0.1 mm d⁻¹ (SAM 0 %). 18 UTC: SAM 90 ± 61 g m⁻² with open-cell patches of 200–300 g m⁻²,
Breeze 63 ± 43 g m⁻² with sparser, smaller cells; rain in both confined to a few cells.
So the horizontal organisation differs most at night (Breeze cellular where SAM is
homogeneous), the mean LWP ratio flips sign between 06 and 18 UTC.

### Attribution and verdict

| disagreement | size | attributed to (recorded departure) | status |
| --- | --- | --- | --- |
| nocturnal LWP 60–100 % higher, denser deck, deeper cloud-top cooling, cellular structure | CWP 150 vs 65 g m⁻² (06–12) | prescribed 10 μm r_eff with RRTMGP vs SAM's SBM-diagnosed r_eff (≈ 12–15 μm at 220 cm⁻³; smaller drops → optically thicker cloud → stronger LW cooling → more condensation); the one-layer RRTMGP column; Smagorinsky–Lilly/WENO vs SAM TKE/MPDATA entrainment | accepted physics difference, first candidate for a sensitivity (r_eff 14 μm or diagnosed from P3's droplet number) |
| drizzle 8× lower over 24 h; evening RWP 20× lower | 0.03 vs 0.26 mm d⁻¹ | P3 warm-rain (autoconversion/accretion at Nc ≈ 220 cm⁻³) vs the SBM's explicit collision–coalescence; the SS-cap emulation keeps Breeze's Nc at the ceiling (SAM's actual Nc unknown without sam3dmicro); Breeze's early-morning drizzle comes from the fog layer instead | accepted microphysics difference (P3 vs bin); not an adapter flaw |
| midday cloud-fraction minimum 0.50 vs 0.93, 2 h later | CF 0.96 vs 0.99 (09–12) | thinner deck in SAM breaks up less because its LWP is lower and more uniform; shortwave absorption with the fixed r_eff | consequence of the first row |
| free-tropospheric v 2–3 m s⁻¹ different after 06 UTC | 1–3 km | forcing above the LES: the lsf wind nudging/geostrophic columns are identical, so the difference arises from the sponge/upper boundary (SAM's `upperbound` on the top two levels vs Breeze's) and the RRTMGP column top | to check against the forcing input before a fix is attempted |
| surface fluxes and net radiation not comparable | — | run saved no flux/radiation series | writer gap: add SHF/LHF (bulk BC diagnostics) and surface/top net LW/SW to the time-series writer before the next run |

Verdict: the adapter reproduces the member's thermodynamic structure, diurnal cycle and
cloud geometry within the stated departures; nothing indicates an input, unit, forcing or
coding flaw that would warrant a relaunch. The LWP/drizzle differences are the P3-vs-bin and
effective-radius differences the audit predicted, now quantified (CWP +130 %, precipitation
−88 % over 24 h). Recommended next steps, in order: (1) add the flux/radiation writers;
(2) rerun with `liquid_effective_radius = 14e-6` (SAM's ocean default when `compute_reffc`
is off) and, if Breeze exposes it, a droplet-number-diagnosed r_eff, to isolate the radiation
contribution; (3) download one `sam3dmicro` time for the SBM droplet number.

## 2. What `protocol = :lasso_ena_official` does

`inspect_lasso_bundle` / `validate_lasso_bundle` (`src/lasso_bundle.jl`) read the five
bundle files and either return the protocol defaults for `build_case` or reject the
bundle with every offending setting named. The namelist conventions were taken from
`setparm.f90` (`&PARAMETERS`), `SGS_TKE/sgs.f90` (`&SGS_TKE`) and `params.f90`/`grid.f90`
(defaults for absent keys).

**Required explicitly (never inferred):**

- the member run ID (`member=`; aerosol setting and SAM microphysics come from it, because
  the HUJI-SBM spectrum is compiled into the SAM binary, not read from `prm`);
- the UTC `epoch` (SAM's `day0` has no year; it is checked against `day0` and the member date);
- the domain `dimensions=(Nx, Ny, Nz)` unless `prm` carries `nx_gl/ny_gl/nz_gl` or a numeric
  `NxxNyxNz` case ID (SAM compiles the domain into `domain.f90`; the documented d25x100 domain
  256×256×260 can be requested as `dimensions=:documented` in `ena_lasso` and is recorded as such).

**Read from the bundle:** `dx`, `dy`, `dt` (fixed step: `Δt = max_Δt = dt`), `nstop`
(duration `nstop × dt`), `day0`, `latitude0`/`longitude0` (f-plane), `tauls` (wind-nudging
timescale), `nrad` (radiation every `nrad × dt`), `doupperbound`, `dodamping`,
`read_in_geostrophic_wind` (7- or 9-column `lsf`), `perturb_type` (cases 0 and 5 implemented),
`compute_reffc` (prescribed effective radius: 14 μm when `.false.`, SAM's ocean value), `ug`/`vg`
(translation frame, recorded, not applied), `dosmagor` (recorded). Scalar levels come from
`grd` with SAM's interface rule (`setgrid.f90`: interfaces halfway between levels, top
interface one half-spacing above the top level; missing lines extended with the last spacing).

**Rejected with the reason:** prescribed fluxes (`SFC_FLX_FXD`, `SFC_TAU_FXD`), land or CEM
modes, per-column fluxes (`UNIFORM_SFC_FLX = .false.`), any T/q nudging (`donudging_tq/t/q`,
transient or inversion-following nudging — the `n1`/`n2` members), a wind-nudging height window,
delayed forcing (`timelargescale > 0`), `rad_simple` or missing LW/SW, prescribed radiative
forcing, restarts, subensembles, tracers, smoke, column mode, walls, SCAM input, perpetual or
fixed-sun insolation, homogenized radiation/SST, slab ocean, scaled CO₂, removed trace gases,
a Coriolis parameter detached from the latitude, a constant-`dz` grid, `perturb_type ∉ {0, 5}`,
a sounding that does not reach the model top (SAM would continue with the 1976 standard
atmosphere), `snd` records not bracketing `day0`, `lsf` starting after `day0` and `sfc` not
covering the whole run (SAM extrapolates these linearly; Breeze clamps).

**Conversions:** pressures hPa → Pa, moisture g kg⁻¹ → kg kg⁻¹ *dry mixing ratio*; the
initial sounding is the `snd` pair bracketing `day0` interpolated in time (`setdata.f90`), its
mixing ratio interpolated to the grid *before* conversion to a mass fraction `q = r/(1+r)`;
heights of pressure-level records from SAM's dry hydrostatic recipe (`R = 287`, `cp = 1004`,
`g = 9.81`); `lsf` columns interpolated linearly in height (height records) or pressure
(pressure records) with zero `tls/qls/wls` and held winds above the record top (`forcing.f90`);
`qls` (a dry-basis vapor source) converted to a mass-fraction rate `dq/dt = qᵈ²/(1−qᶜ) qls`
(`= qls/(1+r)²` in clear air); `tls` converted to the energy prognostic so that `dT/dt = tls`;
`wls` applied once, as full-field first-order upwind advection of every prognostic
(`subsidence.f90`, levels 2…Nz−1), never as mean-profile subsidence in addition; winds are
ground-relative (no translation frame; bulk fluxes need the absolute wind). All of this is
checked by `test/ena_lasso.jl` on the synthetic fixture at several heights and times.

## 3. Physics differences that matter (audit)

| Component | SAM member (`lasso_ena_noice` @ 12d0244) | Breeze adapter | Status |
| --- | --- | --- | --- |
| Surface fluxes (`flxsst`) | `surface.f90` LES branch: **one** `oceflx` call per step with the *domain-mean* lowest-level wind (plus `ug,vg`), mean `T`, `q` and the SST; fluxes and stress are horizontally uniform. `oceflx.f90`: Large & Pond neutral drag `cdn = 0.0027/U + 0.000142 + 0.0000764U`, Stanton `0.0327√cdn` (unstable) / `0.018√cdn` (stable), Dalton `0.0346√cdn`, two Monin–Obukhov iterations, `umin = 1 m/s`, `ssq = qsatw(SST, p₁)` with no salinity factor, "potential temperature" taken as `t₀(1) = T + gz/cp`. Latent heat enters only through the vapor flux (T-neutral evaporation by construction of `t`). | Breeze `BulkDrag`/`BulkSensibleHeatFlux`/`BulkVaporFlux` with `PolynomialCoefficient`: **per-column** fluxes from local wind/θ/q, Large & Yeager neutral polynomials (drag identical; scalar coefficients differ), Li et al. (2010) non-iterative stability, `gustiness = 0.1 m/s`, saturated wall at the SST; θ-flux formulation keeps evaporation temperature-neutral. | Deliberate difference: local vs uniform fluxes and a different stability closure. Expect differences in flux variability (cold pools, gusts) and a few W m⁻² in the mean. A SAM-`oceflx` boundary condition on the domain-mean state is the direct fix. |
| Radiation column above the LES | `RAD_RRTM/rad.f90`: RRTMG LW+SW on `nzm+1` layers — the LES plus **one** layer from the model top (≈350 hPa at 8 km) to `≤1e-4 hPa`, with `T` extrapolated from the two top levels, `h₂o` held at the top-level value, trace gases (O₃, CO₂, CH₄, N₂O, O₂, CFCs) path-integrated from the RRTMG `rrtmg_lw.nc` climatology; `LWP = qcl + qpl` (rain included), ice from `qci`; surface emissivity **0.95**; ocean albedo zenith-dependent (`albedo()` from CAM); SW only when the sun is up; `doseasons = .false.` keeps `day0`'s declination. | RRTMGP all-sky on the LES layers only (`domain_nlay = Nz`, TOA at the top face), well-mixed CO₂/CH₄/N₂O numbers and an O₃ profile from `BackgroundAtmosphere`, cloud optics from the model's liquid condensate; emissivity set to 0.95 by the adapter; constant albedo 0.07; actual date. | Both truncate the column crudely but differently; the SAM top layer carries the LES top humidity to TOA. Compare clear-sky LW/SW fluxes at the surface and LES top before attributing cloud differences. Rain in SAM's optical depth and the albedo law are recorded differences. |
| Effective radii | `compute_reffc = .false.` (SAM default): `computeRe_Liquid` → **14 μm over ocean**; `.true.`: diagnosed from the SBM spectrum, clamped 2.5–60 μm. Ice irrelevant (no ice branch). | Prescribed constant; the adapter sets 14 μm when `compute_reffc = .false.` and keeps 10 μm (recorded) when SAM diagnoses it. | Resolved for the default; a diagnosed radius needs a P3-based `effective_radius` model. |
| Wind nudging | `nudging.f90`: domain-mean `u₀/v₀` → `ul0/vl0` (the `lsf` winds, time-interpolated, minus `ug/vg`) with `1/tauls`, whole column by default. | `MeanProfileNudging` of the horizontal mean toward `uls/vls`, same timescale, whole column (a height window is rejected). | Equivalent. |
| T/q nudging | Available (`n1`/`n2` members) but off for the `n0`/plain members. | Not implemented; rejected. | Not needed for the adopted member. |
| Geostrophic forcing / Coriolis | `coriolis.f90`: `f (v − vg0)`, `−f (u − ug0)` with `ug0/vg0` from the `lsf` geostrophic columns (or aliased to `uls/vls` when `read_in_geostrophic_wind = .false.`); `f = 4π/86400 sin φ`. | `FPlane(latitude)` + `TimeVaryingGeostrophicForcing` with the same columns. | Equivalent. |
| Large-scale tendencies and subsidence | `forcing.f90` adds `tls` to `t` and `qls` to vapor each step; `subsidence.f90` advects `u, v, t` and every microphysical field upwind with `wls`, levels 2…Nz−1, once per step after the tendencies. | Same terms as forcings on `θˡⁱ` (energy), vapor and every microphysical prognostic; upwind on Breeze cell centers. | Equivalent up to `−w∂z t` vs `−w∂z θ` (identical to first order in `Tv/T`) and the center/level mismatch below. |
| Grid mapping | Scalar levels `z(k)` from `grd`; interfaces `zi(k+1) = ½(z(k)+z(k+1))`, top `zi(nz) = z(nzm) + ½(z(nzm) − z(nzm−1))`. On the stretched part the levels are **not** the interface midpoints. SAM's `pres`/`rho` come from a virtual-temperature hydrostatic integration of the sounding (`setdata.f90`, `epsv = 0.61`); `pres0` follows the `lsf` surface pressure in time. | Faces from the same interface rule; Breeze fields live at face midpoints (offset recorded, 0 on the uniform 25 m part); Breeze `ReferenceState` from the sounding with Breeze thermodynamics; fixed surface pressure. | Small, recorded; forcing profiles are interpolated to the Breeze centers, so pressure-level records see Breeze's reference pressure. |
| Upper boundary and sponge | `upperbound.f90`: top two levels relaxed to the sounding (`τ = 3600 s`); `damping.f90`: top 30 % geometric Rayleigh damping 1800 → 60 s on `u−u₀`, `v−v₀`, `w`. | Same (`upper_boundary_relaxation_forcings`, `SAMSponge`). | Equivalent. |
| Initial perturbation | `setperturb.f90` case `perturb_type` (namelist; default 0: ±0.02(6−k) K in the lowest five levels; case 5: ±0.1 K, ±0.025 g/kg below 600 m), one uniform random number per cell. | Both cases implemented from the namelist value; a deterministic host RNG. | Equivalent in distribution, not in realization. |
| Aerosol and activation | HUJI-SBM: CCN spectrum `FCCNR_mp` compiled in (the public commit hard-codes a DYCOMS-II RF01 "low" spectrum — 125 cm⁻³ at 0.011 μm, σ 1.2, plus 65 cm⁻³ at 0.06 μm, σ 1.7 — while the LASSO methodology table gives aer1/2/3 as 138/276/552 cm⁻³ at 0.018 μm (σ 1.53) plus 140.5/281/562 cm⁻³ at 0.066 μm (σ 1.78); the production members therefore used a per-member edit of this block). Constant mixing ratio with height (`FCCN0 = FCCNR_mp ρ(k)/ρ(1)`). `diagCCN = .true.` (compile-time): each microphysics call rebuilds the CCN spectrum from the initial one minus the droplet number, removing the **largest** bins first; maximum activation supersaturation 0.3 %; solute density 1.79 g cm⁻³, molecular weight 115, 3 ions. | P3 `AerosolActivation` with the methodology's two lognormal modes for the member's `aerN`, same chemistry, prognostic depleting number, `DiagnosticCCNProjection` once per step restoring the *total* number (mode shape kept). | Material: the SBM's residual CCN are the small, hard-to-activate particles; a number-only reservoir keeps the easy ones, so Breeze activates more readily. Timing (per step vs per call) is a smaller difference. The spectrum used by the actual member must be confirmed from the member's `prm`/metadata when the bundle arrives. |
| Microphysics | HUJI spectral bin, warm, 33 mass bins, sedimentation of the spectrum, explicit collision kernels; rain = drops above a size threshold. | P3 (two-moment cloud/rain with P3 ice categories unused in warm conditions). | Model-form difference by design; compare LWP, cloud fraction, precipitation statistics, not spectra. |
| SGS and advection | `SGS_TKE` with `dosmagor = .true.` by default (**diagnostic Smagorinsky**, not prognostic TKE unless the namelist says so); MPDATA scalar advection. | Smagorinsky–Lilly; WENO (bounds-preserving for water masses). | Recorded; `dosmagor` read from `&SGS_TKE` and written to provenance. |
| Time stepping | Fixed `dt`, third-order Adams–Bashforth; forcing/nudging/damping/surface/Coriolis/pressure/advection/upper boundary/microphysics/radiation in `main.f90` order. | Fixed `dt` (adaptive stepping is an explicit, labelled override), RK3. | Recorded. |

## 4. Observations for 18 July 2017

Staged at `/shared/home/greg/breezelab-runs/20261006/inputs/observations/20170718/` (checksums in
its `manifest.json`): `enamwrret2turnC1.c1` (`phys_lwp` g m⁻² with `qc_phys_lwp`,
`phys_qc_flag`, 1-σ uncertainty; irregular ~1-min samples from 03:10 UTC),
`enaarsclkazrbnd1kolliasC1.c0` (4-s `cloud_layer_base_height`/`cloud_layer_top_height`
for up to 10 layers, `cloud_base_best_estimate`), `enavdisC1.b1` (1-min `rain_rate` mm hr⁻¹
with `qc_rain_rate`). `src/arm_observations.jl` reads them with their QC;
`analysis/compare_ena_observations.jl` compares a run's domain-mean LWP, rain rate, cloud
fraction and cloud base/top with them over the run's own UTC window and labels the run by
its protocol. The LES is a 25.6 km domain mean; the observations are a point/zenith time
series, so the comparison is statistical (means, quantiles, cloudy fractions) and carries
the retrieval uncertainty; matching snapshots is not the criterion.

## 5. Status against the plan's "done" criteria

- Separate run script (`cases/ena_lasso.jl`), constructor (`ena_lasso`), staging/freeze
  tooling (`data_wrangling/stage_lasso_bundle.jl`), launcher (`execution/submit_ena_lasso.sbatch`),
  comparison tooling and this audit: done; CPU tests in `Pkg.test()`.
- Real-bundle validation, the 1-h GPU integration and the full 24-h baseline with its
  SAM/ARM comparison: done (sections 1b–1d). Remaining: flux/radiation writers, the
  r_eff sensitivity, the SBM droplet number from `sam3dmicro`. The Covert benchmark runs are
  not relabelled as LASSO.
