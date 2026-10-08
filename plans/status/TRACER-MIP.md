# TRACER-MIP status

Agent checkout: `/shared/home/greg/breezelab-work/tracer-mip/` (branch `tracer-mip` from origin/main 3330f3c).
Last update: 2026-10-08 05:33 UTC (control job 233 stepping: first 10-min write at 05:31 UTC, same rate as the pilot). Pilot completed; 24-h control running (job 233); offline inner nest tested; Pkg.test passes on the merged tree; PR opened.

## Shared-file ownership announcements
- `Project.toml`: added `SpecialFunctions` as a direct dependency (already in the Manifest; no version changed) for the
  aerosol-profile `erf`. Planned: `[weakdeps] NumericalEarth` + `[extensions] BreezeLabNumericalEarthExt` (no effect on the
  main Manifest; NumericalEarth itself is pinned only in a case sub-environment `cases/tracer_mip/Project.toml`, because
  NumericalEarth's weakdep compat forces RRTMGP 1.0.1 -> 0.22.3 and that must not touch the shared environment).
- `src/BreezeLab.jl`: appended a clearly-grouped TRACER-MIP export block and two includes (`src/tracer_mip/*.jl`).
- `test/runtests.jl`: appended `include("tracer_mip.jl")`.

## Capability inventory (plan step 1) — done
- The ERA5/Atmospheres/NestedModels/Lands modules and the three examples named in the plan live in
  **NumericalEarth.jl** (github.com/NumericalEarth/NumericalEarth.jl), not in Breeze.jl. Pinned revision inspected and
  resolved: `d07eb240a5679b7a6913e4d979e8378e76786ad7` (main, 2026-10-06, v0.8.1). Local read-only clone:
  `/shared/home/greg/breezelab-work/tracer-mip-external/NumericalEarth.jl`.
- Compatibility with the project pins: NumericalEarth compat `Breeze = "0.11"` (pin 5264d3c is 0.11.3) and
  `Oceananigans = "0.111-0.113"` (pin 847e125 is 0.113.3) — resolves. Every `Breeze` symbol its Breeze extension imports
  exists at 5264d3c except `Breeze.AtmosphereModels.surface_precipitation_flux` (renamed upstream to
  `bottom_precipitation_flux` before the pin; Julia only warns and the extension falls back to a zero field, so the child's
  rain never reaches the slab-land bucket — a coupling gap to patch by method extension).
- Side effect: adding NumericalEarth to the BreezeLab environment downgrades RRTMGP 1.0.1 -> 0.22.3 (NumericalEarth weakdep
  compat "0.21.9, 0.22"; Breeze accepts both). Decision: isolate in `cases/tracer_mip/Project.toml` (sub-environment).
- What exists upstream: `nested_atmosphere_model(child_grid, ERA5HourlyPressureLevels(); dates, dir, terrain=ETOPO2022(), ...)`
  builds an ERA5 `PrescribedAtmosphere` parent on its native pressure-level grid, a Breeze compressible child
  (`LatitudeLongitudeGrid`, terrain-following, open lateral BCs + Davies relaxation from parent-derived prognostics,
  ERA5 initialization + adiabatic balance); `SlabLand` (prognostic skin T + bucket/variably-saturated hydrology);
  `AtmosphereLandModel(Simulation(nest), land; radiation=RadiativeTransferModel)` with Monin-Obukhov land fluxes and
  RRTMGP surface coupling; `PrescribedOcean` (SST) for atmosphere-ocean fluxes.
- Precisely missing for the protocol: (a) Breeze-parent -> Breeze-child nesting (parent must be a `PrescribedAtmosphere`;
  a live prognostic parent is noted as future work) — one-way 2 km -> 500 m therefore has to be offline: save the outer
  prognostic state on the inner region and build a `PrescribedAtmosphere` from it; (b) ocean+land under one atmosphere is
  "at most one surface type per cell" (no tile fractions) — a coastal mask is needed; (c) rain-to-land coupling name
  mismatch above; (d) no two-moment/P3 defaults in the nested path (1-moment default; P3 must be passed explicitly);
  (e) no process-rate accumulators in Breeze P3 output; (f) Tier-1 height-dependent prescribed aerosol: Breeze P3
  `AerosolActivation` modes are vertically uniform (prognostic reservoir caps activation but does not scale the modal
  activation spectrum).

## Protocol re-read (plan step 1/2) — done
`cases/tracer_mip/protocol.toml` transcribes Table 1-3 + the aerosol section of the pinned Roadmap.md and the 2025-03 PDF
(identical numbers): 06 UTC starts on 2022-06-17 and 2022-08-07, 24 h; 750²×2 km / 500²×500 m one-way nests centered
29.4719, -95.0792; 95 ACPC scalar levels (first below ground), top ~22 km; Δt 3/1.5 s; outputs 60 min / 10 min / 2 min
(17-01 UTC); radiation every 60 s; Tier 1 fixed aerosols, radiatively inactive; two lognormal modes per case
(Aug 7: 1425 cm⁻³ Dm 49 nm σ 1.8 and 263 cm⁻³ Dm 175 nm σ 1.4; Jun 17: 3443/30 nm/1.5 and 531/136 nm/1.5), ρ_sfc 1.159,
κ 0.26, 50 mg⁻¹ total floor, floor only above 6 km.

## Implemented so far (branch `tracer-mip`, commit 387ba65, pushed 23:14 UTC)
- `src/tracer_mip/protocol.jl`: protocol loader/validator, `acpc_vertical_faces` (94 cells, faces at scalar-level midpoints,
  top 22181 m), `tracer_mip_horizontal_extent`/`tracer_mip_grid` (lat-lon grid with 2 km / 500 m cells at the center
  latitude — a documented "similar model option" to the roadmap's polar-stereographic grid).
- `src/tracer_mip/aerosol_profiles.jl`: exact port of the notebook's two-mode profiles (shape fit, 1.159 conversion,
  50 mg⁻¹ floor split by modal fraction, floor-only above 6 km), sensitivity multiplier, P3 `AerosolMode`s with β_act = κ.
- `test/tracer_mip.jl`: notebook-table reproduction (both cases), floor/ratio invariants, extents, faces, P3 modes.
- `data_wrangling/fetch_era5_tracer_mip.jl`: staging script (pressure levels over the padded outer box, single levels,
  ERA5-Land) with a checksum manifest.
- `src/tracer_mip/prescribed_aerosol.jl` (`PrescribedAerosolProfile`): Tier-1 aerosol for Breeze P3 by method extension on
  P3's aerosol hooks — droplet number prognostic, no reservoir advected/depleted, a static `nᵃ(z AGL)` field, and the
  activation spectrum scaled by n(z)/n(0) (exact for the protocol's shared-shape modes). Test: `test/tracer_mip_aerosol_model.jl`.
- `src/tracer_mip/process_diagnostics.jl` (`ProcessRateAccumulators`): the 12 Table-2 process rates re-evaluated from
  P3's `compute_p3_process_rates` each step and accumulated as kg/kg (or #/mg, K) per output interval; reset callback
  after each write; restart loses only the in-flight interval. Mapping of P3 rates documented in the file header.
- `ext/BreezeLabNumericalEarthExt.jl` (+ `ext/tracer_mip_regional/`): `tracer_mip_outer_simulation(arch; case, parent=:era5|:synthetic, ...)`
  = NumericalEarth `nested_atmosphere_model` (ERA5 pressure-level parent, open BCs + Davies relaxation, ETOPO terrain) with
  P3 + Tier-1 aerosol, BreezeLab bounded/positive scalar advection, SmagorinskyLilly, `SlabLand` initialized from ERA5 skin T
  and ERA5-Land soil moisture, sea cells (terrain ≤ 0) pinned each step to ERA5 skin temperature with a saturated bucket,
  RRTMGP every 60 s with land/sea albedo, fixed Δt = 3 s, process accumulators. The `:synthetic` parent (analytic
  hydrostatic state, labelled `SyntheticTestParent`, `exploratory_synthetic_boundaries = true`) exists for CPU software tests
  and exploratory throughput runs only. Test: `test/tracer_mip_regional.jl` (runs only with NumericalEarth loadable).
- `cases/tracer_mip.jl`: readable outer-domain runner (hourly 3-D state/water/process output, surface fields, provenance,
  snapshot figure); `execution/tracer_mip_outer.sbatch`: UUID-pinned launcher per COORDINATION.md.
- `cases/tracer_mip/Project.toml`: case sub-environment (NumericalEarth d07eb240, CopernicusClimateDataStore 0.2.2, RRTMGP 0.22.3).
- Docs: `docs/cases/tracer_mip.md` (+ `docs/make.jl` page, `cases/README.md` row).
All of the above is written; CPU verification is in progress (see Checks). None of it has been run on a GPU.

## Inner 500 m nest (plan step 4) — implemented and CPU-tested on the pilot's saved state
`ext/tracer_mip_regional/inner_domain.jl`: `outer_run_parent(run_dir)` turns `inner_region/inner_region_state.jld2`
(T, p, qᵛ, qᶜˡ, qʳ, qⁱ, u, v, w, ρ every 10 min on the inner extent + 10 outer cells) into a `PrescribedAtmosphere` on a
lat-lon grid whose vertical is NumericalEarth's `PressureLevelVerticalDiscretization` carrying the outer terrain-following
cell heights (static) and the terrain as surface geopotential; `tracer_mip_inner_simulation(arch; outer_run_dir, …)`
nests the 500 m child in it with the same physics as the outer (P3 Tier-1, RRTMGP, slab land + sea pinning, Δt = 1.5 s),
land initialized from the child's lowest level (ERA5 at the inner box not staged yet — recorded in `config.land_init`).
Test `test/tracer_mip_inner.jl` against the pilot's 7-time inner-region output (145×145×94 window): parent built,
heights monotonic to > 20 km, child initialized within 5 K of the parent at the surface, 6 s coupled run finite —
**10/10 pass (job 227)**. Not yet: inner 10-min / 2-min tracking writers and a GPU inner run. (Earlier design text:)
Offline one-way nesting: the outer run saves T, p, qᵛ, qᶜˡ, qʳ, qⁱ, u, v, w, ρ every 10 min on the inner extent + 10 outer
cells (`runs/<exp>/inner_region/inner_region_state.jld2`, README alongside). A parent `PrescribedAtmosphere` is then built
on a `LatitudeLongitudeGrid` whose vertical is NumericalEarth's `PressureLevelVerticalDiscretization` with a static
per-cell geopotential g·z(i,j,k) taken from the outer terrain-following heights (so `surface_elevation(parent)` and the
terrain blend work as for ERA5), and `nested_atmosphere_model(parent, inner_grid; terrain = ETOPO2022(), …)` nests the
500 m child exactly as the outer is nested in ERA5. Still to do: write the static outer height field with the saved state,
the parent-from-JLD2 constructor, the 10-min/2-min inner writers, the `bottom_precipitation_flux` rain-to-land patch.

## Storage estimate for the full protocol output (Float32, 94 levels, ~22 3-D fields)
Outer 750² hourly: ≈4.6 GB/write × 25 ≈ 115 GB. Inner 500² 10-min: ≈2.1 GB × 145 ≈ 300 GB; 2-min tracking 17-01 UTC
(240 writes) ≈ 500 GB unless restricted to a tracking subset. Total ≈ 0.9 TB per case; /shared has 2.1 TB free (7 % used).

## ERA5 access — staged (no blocker)
Credentials present (`~/.cdsapirc`, `~/.config/era5cli/cds_key.txt`, env in `~/.bashrc`; never logged). Backend:
CopernicusClimateDataStore.jl 0.2.2 (pure-Julia CDS API v2) via NumericalEarth's extension, run on the CPU partition.
Staged under `/shared/home/greg/breezelab-work/tracer-mip/data/era5/` with checksum manifests `MANIFEST_era5_aug07.toml`
(job 199) and `MANIFEST_era5_jun17.toml` (job 201): per case 333 pressure-level files (9 variables × 37 levels × 37 hourly
times 00 UTC event day → 12 UTC next day, box 103.33–86.83 W / 22.23–36.72 N = outer box + 0.5°), 444 single-level files
(12 variables, unpadded outer box), 9 ERA5-Land files at 06 UTC (skin/soil T, soil moisture 4 layers). The single-level
geopotential at the padded box (NumericalEarth clips sub-surface levels with it) was fetched by the smoke run and is being
added to the staging script's top-up (jobs submitted 01:20 UTC). Upstream gap fixed locally: `nc_varnames` has no
ERA5-Land method in NumericalEarth d07eb240 (the staging script defines it; candidate upstream PR).

## Checks and results (all on the CPU partition cpu-c6a24xlarge; login node is 4 cores shared by four agents)
- `Pkg.test(; allow_reresolve=false)` on the package environment: **passes** (job 203, commit 20eb05e; again job 234 on the
  merged tree 3e7cec0 that the control runs from): BreezeLab 237 + ENA
  protocols 28 + downloads 16+4 + ENA case execution 132 + TRACER-MIP protocol 27 + aerosol profiles 93 + Tier-1 aerosol
  in P3 24 + process accumulators 13. Log: `logs/slurm-203.out`.
- Regional machinery tests (`test/tracer_mip_regional.jl`, case sub-environment, synthetic parent, 10×8×24 outer grid,
  20 s coupled run with slab land, sea pinning, open boundaries, P3 Tier-1 aerosol, accumulators): **29/29 pass**
  (job 196, commit 57750ed). Log: `logs/slurm-196.out`.
- ERA5-driven CPU smokes of the real outer domain (32×32×94 reduction of the protocol extent, 36 s = 12 × 3 s steps):
  **with no explicit closure the real-ERA5 outer domain time-steps stably**, with RRTMGP every 60 s (job 212) and without
  radiation (job 213): p ∈ [41, 1018] hPa, T ∈ [201, 308] K, ρ ∈ [0.068, 1.17], droplets activate to the prescribed
  surface aerosol number (nᶜˡ ≤ 1.23e9 kg⁻¹ = 1229 mg⁻¹). Every run with `SmagorinskyLilly` (jobs 204, 207, 209 flat
  terrain, 210 Δt = 1 s) died in the first step with a negative pressure (`DomainError` in P3's `(p/pˢᵗ)^κ`): the
  isotropic Smagorinsky length on the 50 m × 47 km (test) / 50 m × 2 km (protocol) grid is the culprit, so the
  constructor default is now `closure = nothing` (numerical diffusion only, as NumericalEarth's own ERA5 example) and
  `TKEBasedTurbulenceClosure` (vertical, implicit) also steps stably on the same smoke (job 215, 36 s, same envelope as no closure) and is the
  candidate PBL scheme for the control (knob `TRACER_MIP_CLOSURE=tke`).
  Job 208 died with a bus error while four smokes wrote the same ETOPO field cache concurrently (artifact, not a code bug).
  Outputs of the smokes: `runs/era5_cpu_smoke_job<id>/` (state/process/surface/slices JLD2, inner_region, provenance.toml).
- Fixes found by these runs: `convert_profile` signature, reduced grids below the halo, `fill_halo_regions!` import,
  synthetic-parent spacing, `parent_condensates = nothing` for condensate-free parents, the `surface_precipitation_flux`
  rain-to-land shim (NumericalEarth probes `applicable(f, model, μ)` then calls `f(model)`), ERA5-Land `nc_varnames`,
  `ApparentSolarPosition` keyword `epoch` in the shared ENA path (same fix as branch tracer-dp-scream).

## GPU runs
- **Job 216** (submitted 02:07 UTC): real-ERA5 outer-domain pilot, full 750×750×94 grid, 1 simulated hour, Δt = 3 s,
  RRTMGP every 60 s, P3 + Tier-1 aerosol, process accumulators, no explicit closure, hourly 3-D output + 10-min
  surface/slices + 10-min inner-region state. GPU GPU-bdbd7e9e-5b00-25f4-cbf3-a2954645c086 (node gpu-p4de-2c-st-gpu-p4de-2c-1,
  idx 2) verified physically idle (0 MiB, 0 %, no compute apps, no lock) at 02:07:37 UTC immediately before submission.
  Command (from the checkout at commit of 02:07 UTC):
  `sbatch --partition=gpu-p4de-2c --nodelist=gpu-p4de-2c-st-gpu-p4de-2c-1 --gres=gpu:1 --cpus-per-task=12 --mem=96G --time=04:00:00 --export=ALL,EXPERIMENT=outer_era5_pilot_750,TRACER_MIP_NX=750,TRACER_MIP_NY=750,TRACER_MIP_STOP_HOURS=1,TRACER_MIP_PARENT=era5,TRACER_MIP_CLOSURE=none,TRACER_MIP_PROGRESS_SECONDS=60 execution/tracer_mip_outer.sbatch GPU-bdbd7e9e-5b00-25f4-cbf3-a2954645c086`
  Output: `/shared/home/greg/breezelab-work/tracer-mip/runs/outer_era5_pilot_750_job216/` (COMMIT, COMMAND, Manifests,
  gpu_memory.csv sampled every 30 s, JLD2 outputs, provenance.toml, COMPLETE on exit 0). Log: `logs/slurm-216.out`.
  Result: **failed at load** ("failed to find source of parent package: RRTMGP"): the per-CPU depot of node 1
  (CPU tag `Intel_R__Xeon_R__Platinum_8275CL_CPU___3_00GHz_`, not AMD as COORDINATION.md assumed) held no package sources,
  which Pkg had installed into the c6a (AMD) depot. Fixed by `execution/depot_layout.sh`: `packages/`, `artifacts/`,
  `clones/`, `registries/`, `scratchspaces/` of every per-CPU depot are symlinks into `depot/shared/` (sourced by all
  launchers). No GPU compute was used.
- **Job 217** (submitted 02:10 UTC, same command, GPU re-verified idle at 02:10:15 UTC; started 02:18 UTC on
  GPU-bdbd7e9e, A100-SXM4-80GB): **running the full 750×750×94 real-ERA5 outer domain**. Model built and t = 0 output
  written at 02:35 UTC (≈13 min compile/construction for the node's CPU tag); first 10-min output (200 steps of 3 s) at
  02:54:25 UTC → 1136 s wall for the first 200 steps including first-step kernel compilation. **GPU memory: 40.3 GB
  resident** (of 81.9 GB; 27.6 GB after construction, before the first step), GPU utilization 100 %. Output so far:
  `outer_state.jld2` 1.23 GB per hourly 3-D write (19 fields), `inner_region_state.jld2` 45 MB per 10-min write.
  **Steady-state throughput: 4.3 s wall per 3 s step** (20 steps per 1.43 min from iteration 20 to 160; the first 20 steps
  took 7.8 min with kernel compilation) → ≈86 min per simulated hour, ≈35 h for the 24-h control on one A100 with
  Float32, P3 + Tier-1 aerosol, RRTMGP every 60 s, accumulators, no closure. The 500² inner nest (0.44× the cells) will
  need ≈18 GB and, run offline, ≈15 h per 24 h at the same per-cell cost with Δt = 1.5 s (2× the steps, 0.44× cells).
  **Completed (COMPLETE marker, exit 0) at 04:10 UTC: 1200 steps = 1 simulated hour in 1.586 h wall** (1.33 h of
  stepping + 13 min build + 10 min first-step compile). State stayed bounded and finite throughout: T ∈ [200.8, 305.2] K,
  p ∈ [40.7, 1019] hPa, ρ ∈ [0.067, 1.175]; max|u| rose from 25.6 to 64 m/s within 8 min then relaxed to ≈59 m/s;
  max|w| ≈ 6 m/s with one 12 m/s excursion at 29 min; cloud water ≤ 1.5 g/kg, droplet number up to 2.7e9 kg⁻¹.
  The wind maximum's location (relaxation zone / lid sponge / interior) and hydrostatic residual are being diagnosed
  by `analysis/tracer_mip_outer_pilot.jl` (job submitted 04:11; writes `pilot_summary.toml` and figures into the run dir).
  Output: hourly 3-D state 1.38 GB/write, inner-region state (145×145×94 window + halo) 45 MB per 10 min, surface and
  2-km slices every 10 min, process-rate accumulators per hourly interval. Note: harmless JLD2Writer warnings about
  `thermodynamic_constants` (the writer inspects the coupled model, not the child).
- **Job 219** (submitted 03:06 UTC with `--dependency=afterok:217`, started 04:10 UTC on the same GPU after the pilot
  exited; the in-job lock/idleness check passed): **full August 7 control, 750×750×94, 24 simulated hours**, same
  configuration as the pilot, `--time=48:00:00`, progress every 10 simulated minutes. Expected ≈35 h wall (finish
  ≈15:00 UTC 2026-10-09). Output: `runs/outer_era5_control_aug07_job219/`. It is an *outer-domain* control; the inner
  500 m nest is run offline from its saved inner-region state afterwards.
  **Cancelled at 04:33 UTC after 23 min** (`scancel 219`): the pilot diagnostics showed the hourly process-rate
  accumulators saved as zeros — Oceananigans runs callbacks before output writers within a step, so my reset callback
  zeroed them before the writer; fixed by deferring the reset to the next accumulation (commit of 04:32, with a
  writer-order test, job 228 passes; accumulation once per completed step). Also the inner-region velocities are now
  saved at cell centers. **Relaunched as job 233** at 04:53 UTC from commit 3e7cec0 (merged `tracer-mip-inner`), GPU
  re-verified idle (0 MiB, 0 %, no compute apps, no lock) immediately before submission; same command with
  `EXPERIMENT=outer_era5_control_aug07`; output `runs/outer_era5_control_aug07_job233/`; expected ≈35 h wall.
- Pilot diagnostics (`runs/outer_era5_pilot_750_job217/pilot_summary.toml`, `pilot_overview.png`, `pilot_land_response.png`,
  job 223): at t = 1 h the wind maximum (59.4 m/s) sits at cell i = 2, 76 m height — inside the 5-cell lateral relaxation
  zone of the west boundary (rim max 60 m/s vs **interior max 36.5 m/s**, the jet); max|w| 6.2 m/s at 440 m near the south
  boundary; hydrostatic residual |∂p/∂z + ρg|/(ρg) median 0.37 % (p99 9 %; the diagnostic uses reference heights, so
  terrain columns inflate it); land skin temperature cools 298.9 → 297.8 K over the nocturnal hour, sea cells pinned at
  302.75 K (sea fraction 31.8 %); surface rain starts (max 1.4e-6 kg m⁻² s⁻¹); cloud water ≤ 1.7 g/kg in 0.9 % of columns.
  Boundary-zone wind amplification at the inflow boundary is flagged for the control's evaluation (relaxation
  rate/width, surface drag in the rim).

## Integration
Branch `tracer-mip` contains `origin/main` (293081f). Pull request https://github.com/PRECISION-GM/BreezeLab.jl/pull/4 opened 05:20 UTC against `main`, per the
coordinator's instruction not to push to main; it lists implemented behaviour, test evidence, run job IDs/paths and the
explicit protocol departures. Shared files touched: `Project.toml` (SpecialFunctions dep, NumericalEarth weakdep/extension
declaration), `src/BreezeLab.jl` (one grouped export block + includes + two stub declarations), `src/case_setup.jl`
(`ApparentSolarPosition` keyword, same as tracer-dp-scream), `test/runtests.jl` (three includes), `docs/make.jl`, `cases/README.md`.

## Next
1. Monitor job 233 to completion (≈15:00 UTC 2026-10-09); run `analysis/tracer_mip_outer_pilot.jl` on it (sea-breeze
   development, boundary-zone winds, surface energy/water budget, storage: ≈36 GB).
2. Inner 500 m GPU run from its saved state: `tracer_mip_inner_simulation`, add the 10-min full and 2-min (17–01 UTC)
   tracking writers, stage ERA5 at the inner box for land initialization, size the run (≈18 GB, ≈15 h).
3. June 17 outer control (ERA5 staged) and the aerosol variants once the LOW factor is resolved.
4. Evaluate `TKEBasedTurbulenceClosure` in a GPU pilot; investigate the inflow-rim wind amplification.

## Blockers
None concrete at this time.
