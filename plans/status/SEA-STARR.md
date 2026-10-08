# SEA-STARR status

Agent checkout: `/shared/home/greg/breezelab-work/sea-starr/` (branch `sea-starr` from origin/main 3330f3c).
Status updated: 2026-10-08 01:40 UTC.

## Commit / branch
- `sea-starr` @ b54dafd (pushed; merged with origin/main 293081f at 7505dbe, no conflicts). PR #3 open: https://github.com/PRECISION-GM/BreezeLab.jl/pull/3 (no direct push to main, per the coordinator).

## Inputs
- Staged drivers at `/shared/home/greg/breezelab-runs/20261006/inputs/mip_sources/seastarr_22241697/`;
  all six SHA-256 checksums verified against `cases/seastarr/inputs.toml` (22:22 UTC). Covert inputs copied into
  `data/` and verified against `data/MANIFEST.txt` so the existing ENA tests run.

## Implemented behaviour (files in the checkout)
- `docs/cases/sea_starr.md`: plan step 1 — variable/forcing table from the setup PDF and driver metadata, member
  differences (only `na`/`na_nud` differ), initial-state facts, departures.
- `src/dephy_driver.jl`: `read_dephy_driver` (DEPHY v1; top-down levels reversed, `na` mg⁻¹ → kg⁻¹, uniform
  nudging constants checked), `driver_profile_time_series`, `sea_starr_vertical_faces` (10 m to 2500 m, 10 % to 6500 m, Nz = 288).
- `src/sea_starr_aerosol.jl`: `KappaAerosolMode` + `kappa_aerosol_activation` (κ-Köhler activation of the *local*
  prognostic reservoir, a Breeze method specialized on a BreezeLab type — not piracy), `SurfaceAerosolSource`
  (bottom-cell tendency; P3 ignores flux BCs on ρnᵃ), `EvaporationRegeneration` (droplet number removed ∝ P3's own
  diagnosed cloud-evaporation rate and returned to nᵃ), `aerosol_number_columns`. ENA's `DiagnosticCCNProjection` is NOT used.
- `src/sea_starr_forcings.jl`: per-column inversion height (max ∂θₗ/∂z), `InversionMaskUpdater` (domain max + 100 m,
  linear 200 m ramp, refreshed every step), `InversionFollowingNudging` (horizontal-mean nudging × mask).
- `src/sea_starr.jl`: `sea_starr(; member, arch, FT, Nx, Ny, ...)` constructor (returns the simulation unrun): protocol grid,
  anelastic reference state from `ps`/θₗ/qₜ, P3 + κ aerosol, full-field subsidence on every prognostic exactly once,
  time-varying geostrophic wind, bulk SST fluxes (driver `ts_force`, z₀ = 1e-4), RRTMGP LW+SW with driver ozone,
  nudging of θ/qᵛ/nᵃ (1800 s) and u/v (10800 s), equilibrium initial cloud, JLD2 writers (15-min statistics, 1-min
  time series incl. aerosol budget columns, 15-min 2D, hourly 3D), Checkpointer (3 h, cleanup).
- `cases/sea_starr.jl` (construct → provenance → run → Makie figures), `execution/sea_starr.sbatch` (UUID-pinned launcher),
  `test/sea_starr.jl` (grid, reader on a synthetic DEPHY file, interpolation, mask/nudging, κ activation, source,
  regeneration, subsidence-once, staged-driver constructor + 2 s CPU run).

## Explicit departures from the protocol (recorded in `config.departures` / provenance)
1. No aerosol optics in Breeze's RRTMGP (SSA 0.85 not representable) — smoke absorption absent.
2. RRTMGP column ends at the LES top (6.5 km), no atmosphere above, zero downwelling LW at the top.
3. Fixed solar coordinate (driver lat/lon −11.65, −8.34) rather than the moving composite trajectory.
4. Nudging ramp linear over 200 m (Blossey et al. 2013 paywalled from the cluster); mean-profile nudging; θ ← thetal_nud, qᵛ ← qt_nud.
5. Proportional evaporation regeneration; no regeneration from rain evaporation; no interstitial scavenging.
6. Numerical upper sponge over the top 15 % of the domain (not in the protocol).
7. Initial cloud from warm-phase saturation adjustment of (θₗ, qₜ); all aerosol in cloudy cells initially activated.

## Launch commands / jobs
- none yet (short GPU CTRL planned on gpu-p4de-2c-2 GPU-25cad61d-d859-da9a-cbea-5ca3dc568420, verified idle 22:25 UTC).

## Checks and results
- `Pkg.test(; allow_reresolve=false)` on the login node (CPU), merged tree 7505dbe: PASSED (log `logs/pkgtest_2.log`, 01:27 UTC) —
  BreezeLab 238 pass / 1 broken (pre-existing skip), ENA protocols 28, ARM download 16, manifest 4, ENA execution 132, SEA STARR 98.
  (An earlier run on 3670ff9 also passed, `logs/pkgtest_1.log`.) b54dafd changed only `cases/sea_starr.jl` and `analysis/`.
- **Short GPU pilot, job 193** (A100 GPU-25cad61d, 192²×288, 2 simulated hours, commit 89b15f5): the *integration completed*
  (12 933 steps, 1.128 h wall ≈ 0.31 s/step, Δt steady at 0.55 s under the CFL wizard with the 1 s cap; max|w| 3.6 m/s;
  T ∈ [265, 292] K; all prognostics and all written time series finite). The job exited 1 only afterwards, in the script's
  plotting section (`hour` ambiguous between Dates and Oceananigans.Units) — fixed in b54dafd; figures regenerated on the
  login node from the saved output: `runs/ctrl_pilot_2h_job193/sea_starr_{timeseries,profiles}.png`.
  Checks (`analysis/check_sea_starr_run.jl`, log `logs/check_job193.log`):
  - cloud: CWP 101 → 134 g m⁻² (max 134), RWP ≤ 0.6 g m⁻², cloud fraction 0.997–1.0, surface rain ≤ 0.009 mm d⁻¹;
    cloud base 735 → 755 m, cloud top 1035 → 1165 m (15-min means); initial one-cell inversion smears over ~150 m within the
    first 15 min (diagnosed zᵢ dips 1160 → 1040 m, back to 1160 m by 2 h, max column 1230 m) — an initial transient, not yet validated.
  - aerosol: BL nᵃ 1.00 → 1.52×10⁸ kg⁻¹ (consistent with ~6 % entrainment of the 10⁹ kg⁻¹ smoke layer); column totals finite;
    column total change is 16× the surface source because FT nudging/subsidence dominate it (their integrals are not yet output, so
    the budget is not *closed* in the output — gap).
  - radiation: cloud-top LW cooling −1.5 W m⁻³ at 1125 m; LW↑ 420 (sfc) / 335 (top) W m⁻², LW↓ 394 W m⁻² at the surface, 0 at the
    LES top (departure); SW = 0 (night, 21–23 UTC) as expected.
  - surface fluxes: bulk formulae active (SST 293.3 K vs air 291.4 K) but the heat/moisture fluxes are not written to output — gap.
  - memory 7.4 GiB on the A100; output 2 h: 3D 586 MB, checkpoint 2.9 GB (iteration 0), statistics/timeseries/2D < 10 MB.
- **Full 66 h CTRL, job 211** (H100 gpu-prod-1-1 GPU-d640ffb6, verified physically idle with no pending gpu-prod job at 01:34 UTC;
  commit b54dafd): `sbatch --partition=gpu-prod --nodelist=gpu-prod-st-gpu-prod-1-1 --gres=gpu:1 --cpus-per-task=12 --mem=96G --time=48:00:00 execution/sea_starr.sbatch GPU-d640ffb6-f5b1-efc0-dcc6-232466150c94 ctrl_full_66h 66 192`
  → output `runs/ctrl_full_66h_job211/`, log `logs/slurm-211.out`, live `runs/ctrl_full_66h_job211/run.log`.
  Estimate: 66 h / 0.55 s ≈ 4.3×10⁵ steps; at 0.2–0.3 s/step on the H100 → 24–36 h wall (pilot: 0.31 s/step on the A100).
  Declared storage budget ≈ 20 GB (hourly 3D ≈ 0.2 GB/h → 13 GB; latest 2.9 GB checkpoint kept, `cleanup=true`; rest < 1 GB).
  Checkpoints every 3 h; restart with `SEA_STARR_PICKUP=true` (radiation recomputes at its next 60 s interval after pickup).

## Shared-file ownership announcements
- `src/BreezeLab.jl`: appended SEA STARR includes/exports only (no change to ENA exports).
- `test/runtests.jl`: appended `include("sea_starr.jl")`.
- `cases/README.md`: SEA STARR row updated. No `Project.toml` changes.

## Blockers
- Blossey et al. (2013) ramp form unverifiable offline (paywall) — implemented linear ramp, labelled.
- Aerosol optics / atmosphere above the LES top / moving solar geometry need Breeze-side features (RRTMGP aerosol state,
  extended column, time-varying coordinate); recorded as departures, not implemented here.

## Launch log (appended 2026-10-08T00:02:31Z)
- Branch pushed: sea-starr @ 3670ff9 (origin/sea-starr).
- Short GPU pilot submitted (GPU-25cad61d physically idle, no lock, re-verified immediately before submit):
  `sbatch --partition=gpu-p4de-2c --nodelist=gpu-p4de-2c-st-gpu-p4de-2c-2 --gres=gpu:1 --cpus-per-task=12 --mem=96G --time=04:00:00 execution/sea_starr.sbatch GPU-25cad61d-d859-da9a-cbea-5ca3dc568420 ctrl_pilot_2h 2 192`
  → job 190; output `/shared/home/greg/breezelab-work/sea-starr/runs/ctrl_pilot_2h_job190/`; log `logs/slurm-190.out`.
- `Pkg.test()` running on the login node (CPU only), log `logs/pkgtest_1.log`.
- Job 190 FAILED at RRTMGP GPU kernel compile (ozone closure over the driver → non-isbits FunctionField); fixed in 89b15f5 (ozone as a Field).
- Resubmitted (2026-10-08T00:13:27Z, GPU re-verified idle): same sbatch command → job 193; output `runs/ctrl_pilot_2h_job193/`, log `logs/slurm-193.out`.
- `Pkg.test(; allow_reresolve=false)` on the login node (CPU, source 3670ff9 + working tree at start): PASSED —
  BreezeLab 238 pass / 1 broken (pre-existing skip), ENA protocols 28, ARM download 16, manifest 4, ENA execution 132,
  SEA STARR 98 (log `logs/pkgtest_1.log`, 00:51 UTC). A rerun on the final commit is scheduled after the pilot.
- Job 193 (pilot) stepping: at 00:43 UTC t = 30 min simulated, GPU 100 %, 7.4 GiB; ≈0.38 s wall/step, mean Δt ≈ 0.6 s
  (CFL wizard below the 1 s cap once drizzle forms).
- Job 211 confirmed stepping (01:48 UTC): iteration 702, t = 420 s, CWP 118 g m⁻², H100 at 99 %, 7.6 GiB; output files growing
  (`runs/ctrl_full_66h_job211/sea_starr_ctrl_{timeseries,statistics,2d,3d}.jld2`, `run.log`). Early rate ≈ 0.4 s/step including
  kernel warm-up; if it does not improve toward the pilot's 0.31 s/step, 66 h (≈4.3×10⁵ steps) needs ≈ 48 h, i.e. the wall limit —
  in that case restart from the latest 3-hourly checkpoint with `SEA_STARR_PICKUP=true SEA_STARR_OUTPUT=<same dir>`.
  Estimated completion: 2026-10-09 ~02–14 UTC depending on the sustained rate. Not waited for; "simulation completed" ≠ validated.

## Remaining gaps / blockers (concrete)
- Breeze-side features needed for a protocol-complete CTRL: RRTMGP aerosol optics (SSA 0.85), an atmosphere above the LES top in the
  radiation column, and a time-varying solar coordinate.
- Not yet output: surface sensible/latent heat fluxes, nudging/subsidence tendency integrals (needed to close the aerosol column budget),
  the workbook's clear-sky/TOA radiation set and 3D optical thickness; 3D output is a 10-variable subset.
- Blossey et al. (2013) ramp form unverified (paywalled); linear ramp implemented.
- Restart equivalence (checkpoint pickup) not yet exercised on this case.
- N100/N030 members: constructor accepts `member=:N100/:N030` with their own drivers (reader tested); no runs launched.

## Coordinator launch (2026-10-08 ~02:05 UTC)
- **Job 214: N100 full 66 h** on A100 GPU-25cad61d (node gpu-p4de-2c-2 idx 6, verified 0 MiB / no compute apps / lock free): `SEA_STARR_MEMBER=N100 sbatch --export=ALL --partition=gpu-p4de-2c --nodelist=gpu-p4de-2c-st-gpu-p4de-2c-2 --gres=gpu:1 --cpus-per-task=12 --mem=96G --time=96:00:00 execution/sea_starr.sbatch GPU-25cad61d-d859-da9a-cbea-5ca3dc568420 n100_full_66h 66 192` at b54dafd. Output runs/n100_full_66h_job214/. N030 to follow when the next A100 frees (ENA jobs 168/169).
- Note: the CTRL full run (job 211) is on the H100 that job 165 released; Greg asked for A100s — left running (moving it would cost the progress so far and A100 is ≈1.7× slower); flagged to Greg.
- **Job 218: N030 full 66 h** on A100 GPU-898f1f68 (node gpu-p4de-2c-1 idx 3, released by ENA job 168; verified no compute apps, no lock): `SEA_STARR_MEMBER=N030 sbatch --export=ALL --partition=gpu-p4de-2c --nodelist=gpu-p4de-2c-st-gpu-p4de-2c-1 --gres=gpu:1 --cpus-per-task=12 --mem=96G --time=96:00:00 execution/sea_starr.sbatch GPU-898f1f68-a597-3afb-61c7-3245cce93e5d n030_full_66h 66 192` at b54dafd. Output runs/n030_full_66h_job218/.
- Job 214 (N100) confirmed stepping at 02:20 UTC: iter 1014, t = 10 min, Δt 0.58 s, LWP 121 g m⁻², CF 1.0, finite.
- Job 218 (N030) confirmed stepping at 02:40 UTC: iter 1020, t = 10 min, Δt 0.60 s, LWP 118.5 g m⁻², CF 1.0, surface rain 0.13 mm d⁻¹ (vs 0.003 for N100 at the same time), finite.

## Double-resolution runs (Greg request, 2026-10-08 ~04:40 UTC)
- Interpretation: the protocol's standard 384² grid at 50 m (19.2 km domain), i.e. double the grid points of the 192² runs; same vertical grid.
- **Job 230: CTRL 384² full 66 h** on H100 GPU-0dd861f5 (gpu-prod-1-2, released by ENA job 166; verified 0 MiB / no compute apps / no pending gpu-prod jobs): `SEA_STARR_MEMBER=CTRL sbatch --export=ALL --partition=gpu-prod --nodelist=gpu-prod-st-gpu-prod-1-2 --gres=gpu:1 --cpus-per-task=12 --mem=96G --time=120:00:00 execution/sea_starr.sbatch GPU-0dd861f5-113e-1654-d3b8-ff2f7b81fdab ctrl_full_66h_384 66 384`. Expected ≈4× the 192² cost per step (≈0.7 s/step on H100 ⇒ ≈3.5 days). N100/N030 at 384² follow as GPUs free up (no A100 is idle beyond the one reserved for TRACER-MIP).
