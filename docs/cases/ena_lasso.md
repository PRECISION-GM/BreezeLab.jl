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
- Real-bundle validation and the 1-h GPU integration: done (sections 1b, 1c). The 24-h
  baseline and its SAM/ARM comparison: queued/pending; see the status file. The Covert
  benchmark runs are not relabelled as LASSO.
