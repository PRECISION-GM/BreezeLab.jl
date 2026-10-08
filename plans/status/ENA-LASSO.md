# ENA-LASSO status

Agent: ENA-LASSO. Checkout: `/shared/home/greg/breezelab-work/ena-lasso` (branch `ena-lasso` from origin/main 3330f3c).
Last update: 2026-10-08 ~09:15 UTC. PR: https://github.com/PRECISION-GM/BreezeLab.jl/pull/2 (branch merged with origin/main 293081f; no direct push to main per coordinator).

## Blocker (external, exact) — updated 22:35 UTC
- Plan candidate `20170718era5s1n0d25x100_sbmwrm-aer2-flxsst` DOES NOT EXIST at ARM. Authenticated ARM Live listing
  (`livedata/query?ds=enalasso&start=2017-07-17&end=2017-07-19&wt=json`, HTTP 200, 2502 files; `ds=enalassoC1` → 0 files)
  shows 18 samin archives: `s1n0` members exist only with Morrison (`morr`); the spectral-bin warm members are
  `era5d25x100_sbmwrm-aer{1,2,3}-flxsst` and `merra2d25x100_sbmwrm-aer{1,2,3}-flxsst`. Listing saved (no token) at
  /tmp/claude-1237001154/-shared-home-greg-BreezeLab-jl/34bc4fe9-74de-48df-b2f3-df7170261300/scratchpad/arm_query/.
- Member adopted EXPLICITLY per coordinator: `20170718era5d25x100_sbmwrm-aer2-flxsst`
  (archive `enalasso_samin_20170718era5d25x100_sbmwrm-aer2-flxsstC1.m0.20170718.000000.tar`, reference outputs
  `enalasso_samstat_...C1.m1.20170718.000000.nc`, `enalasso_sam2d_...C1.m1.20170718.000000.nc`).
- Download attempts (one per route, 22:31–22:33 UTC, credentials from ~/.bashrc, never logged):
  `data_wrangling/fetch_arm_inputs.jl` → "ARM returned HTTP 404"; direct `armlive/saveData` and `armlive/livedata/saveData`
  → HTTP 404 with EMPTY body for both the samin tar and the (small) samstat .nc. A control request for a known-online
  observation file (`enavdisC1.b1.20170718.000000.cdf`) on the same route with the same credentials succeeded
  (see below), so credentials and route work: the LASSO files are listed but NOT staged online (bundle browser:
  "the data typically resides on tape"). Staging needs an ARM order (Data Discovery / bundle browser YAML order) — a
  human action for Greg; no ordering campaign was started by this agent.
- Therefore no LASSO GPU run is possible yet; no job submitted. ENA-covert jobs 165–169 untouched.
- Re-checked for the coordinator (01:10 UTC): still no .tar under data/ — the download did NOT succeed (HTTP 404 on both
  routes, see above); the next step is an ARM staging order by Greg, then the launch command below.

## Implemented behaviour (branch `ena-lasso`, see `git log`)
- `src/lasso_bundle.jl`: `parse_lasso_member` (run-ID tokens), `read_sam_namelist_groups`, `lasso_scalar_levels`
  (setgrid.f90 extension rule), `inspect_lasso_bundle`/`validate_lasso_bundle`: required explicit member/epoch/dimensions,
  ~40 namelist switches checked against SAM defaults with the reason, time coverage of snd/lsf/sfc (setdata/forcing/
  setforcing.f90 semantics), sounding-vs-model-top check, lsf layout vs read_in_geostrophic_wind, epoch vs day0 vs member
  date; settings derived: fixed dt, nstop*dt, nrad*dt, tauls, SAM emissivity 0.95, 14 μm ocean effective radius
  (compute_reffc=.false.), perturb_type 0/5, translation frame recorded not applied. Provenance gets `bundle` (member,
  checksums, SAM revision, grid/time facts, warnings) and `staging` (bundle.toml).
- `src/ena_lasso.jl`: `ena_lasso(bundle_dir; member, epoch, dimensions=:documented, ...)`, `lasso_bundle_directory`,
  `lasso_bundle_available`, `lasso_bundle_missing_message` (fail-fast with archive name/DOI/staging command).
- Fidelity fixes in shared code: initial sounding time-interpolated to day0 (setdata.f90), file mixing ratio interpolated
  before the mass-fraction conversion, setperturb.f90 case 0 added. ENA-covert path unchanged in behaviour (its first
  sounding is at day0).
- `src/arm_observations.jl`: MWRRET LWP, VDIS rain, ARSCL cloud-boundary readers with QC, window statistics,
  BreezeLab run-output readers; `analysis/compare_ena_observations.jl` (labels Covert runs as the benchmark they are).
- `cases/ena_lasso.jl`, `data_wrangling/stage_lasso_bundle.jl` (freeze: sha256, member tokens, DOI, SAM revision),
  `data_wrangling/fetch_arm_inputs.jl` now also accepts samstat/sam2d names, `execution/submit_ena_lasso.sbatch`
  (COORDINATION.md UUID-lock procedure), `execution/wpcluster_env.sh`, `docs/cases/ena_lasso.md` (audit).
- Deviation note: the ARM discovery query and single download attempt were made as soon as the coordinator confirmed
  credentials (22:28–22:33 UTC), before the code work finished, because the result (member does not exist; files not
  staged) determined the member identity used throughout the code and docs.

## Checks
- Targeted session (login node, 2026-10-08 00:40 UTC, commit after a78588c): ENA protocols 13/13, member identity 20/20,
  namelist groups 7/7, grd rule 4/4, initial sounding 7/7, perturbation cases 4/4, bundle inspection/validation 61/61,
  adapter units/winds/forcing assembly 67/67, RRTMGP CPU construction+update 5/5 (covers the ApparentSolarPosition
  keyword fix flagged by the coordinator), ena_lasso wrapper 14/14, ARM readers 12/12. Logs: logs/lasso_tests.log.
- Full `Pkg.test()` with the final code (logs/test3.log, exit 0, 2026-10-08 ~01:05 UTC): BreezeLab 238 passed/1 broken
  (pre-existing), ENA protocols 13, member identity 20, namelist groups 7, grd 4, initial sounding 7, perturbation 4,
  bundle validation 61, adapter 67, RRTMGP CPU 5, ena_lasso wrapper 14, ARM readers 12, ARM download safety 16,
  manifest 4, ENA case execution 132 — all passed.
- Comparison tooling exercised on the completed Covert-public-bin benchmark job167 (1M, 06–12 UTC; NOT LASSO):
  /shared/home/greg/breezelab-work/ena-lasso/runs/covert_benchmark_job167_comparison/{comparison.toml,comparison.png}.
  Window 06–12 UTC: MWRRET phys_lwp mean 168 ± 13 (1σ) g m⁻² (median 158) vs LES domain-mean cloud LWP 100 g m⁻²;
  VDIS rain 0.0015 mm hr⁻¹ vs LES 0.003 mm hr⁻¹; ARSCL lowest-layer top 1087 m vs LES cloud top 1218 m; ARSCL lowest-layer
  base median 160 m (drizzle/first radar gate, not the cloud base) vs LES cloud base 732 m — the base comparison needs the
  ceilometer `cloud_base_best_estimate` (reader returns it; the script's next revision should use it). Benchmark evidence
  for the tooling only; no LASSO conclusion.
- No GPU job submitted (no bundle). Jobs 152–170 untouched.

## Exact launch command once the bundle is staged
```sh
cd /shared/home/greg/breezelab-work/ena-lasso && source execution/wpcluster_env.sh
set -a; source ~/.bashrc; set +a   # ARM_USERNAME/ARM_TOKEN (never logged)
$JULIA --startup-file=no data_wrangling/fetch_arm_inputs.jl \
    enalasso_samin_20170718era5d25x100_sbmwrm-aer2-flxsstC1.m0.20170718.000000.tar data/lasso/archives
$JULIA --startup-file=no data_wrangling/stage_lasso_bundle.jl \
    data/lasso/archives/enalasso_samin_20170718era5d25x100_sbmwrm-aer2-flxsstC1.m0.20170718.000000.tar data/lasso
$JULIA --startup-file=no --project -e 'using BreezeLab; b = inspect_lasso_bundle(lasso_bundle_directory(BreezeLab.DEFAULT_LASSO_MEMBER); member=BreezeLab.DEFAULT_LASSO_MEMBER, dimensions=(256,256,260)); foreach(println, b.problems); foreach(println, b.warnings); println(b.time)'
# short integration (1 h) then full nstop*dt on a physically idle A100/H100 (re-verify with ssh nvidia-smi first):
sbatch --partition=gpu-p4de-2c --nodelist=gpu-p4de-2c-st-gpu-p4de-2c-2 --time=02:00:00 \
    execution/submit_ena_lasso.sbatch GPU-25cad61d-d859-da9a-cbea-5ca3dc568420 lasso_aer2_1h \
    --protocol lasso_ena_official --member 20170718era5d25x100_sbmwrm-aer2-flxsst \
    --data data/lasso/20170718era5d25x100_sbmwrm-aer2-flxsst --epoch 2017-07-18T00:00:00 --dimensions 256,256,260 \
    --arch gpu --float Float32 --hours 1
sbatch --partition=gpu-p4de-2c --nodelist=<idle node> --time=24:00:00 execution/submit_ena_lasso.sbatch <idle GPU UUID> lasso_aer2_full
```
(The epoch 2017-07-18T00:00:00 assumes day0 = 199.0; the inspect step prints the namelist day0 and the adapter rejects a
mismatch. Reference outputs samstat/sam2d: `ARM_MAX_BYTES=4000000000 $JULIA data_wrangling/fetch_arm_inputs.jl <name> data/lasso/reference`.)


## Aerosol audit of the ENA-covert P3-aer2 runs (Greg's request, 2026-10-08 ~02:30 UTC)
- Verdict (docs/cases/ena_aerosol_audit.md, evidence runs/ena_aerosol_audit/): fields are not corrupted (nᵃ+nᶜˡ held exactly at
  the initial 4.618e8 kg⁻¹ = 557 cm⁻³ by the diagCCN projection; Δt 0.5 vs 0.25 agree to <3 %), but in-cloud nᶜˡ ≈ 330–400×10⁶ kg⁻¹
  (≈ 350–430 cm⁻³, 60–75 % of ALL aerosol incl. the Aitken mode) vs the observed Nc = 75 cm⁻³ of the Covert case and vs ≤ 270 cm⁻³
  that the SBM's ss_max = 0.3 % allows from aer2. Causes: LASSO aer2 (557 cm⁻³) applied to the 75 cm⁻³ Covert case; Breeze
  activation has no supersaturation cap and ratchets one way to N_act(S_grid); condensate-free start supersaturated by ~6 %
  activates ≈ 430 cm⁻³ in the first seconds; Breeze default mass_fraction_soluble = 0.9 (SBM: fully soluble).
- Fix on branch (ef6a037): lasso_aerosol_modes mass_fraction_soluble = 1 and maximum_supersaturation (SBM-cap emulation via
  Breeze's activated_number; aer2 → ≈ 8 % of mode 1 + 89 % of mode 2 activatable), build_case aerosol_supersaturation_cap,
  CLI --aerosol_ss_cap, recorded in provenance. Pkg.test on the CPU partition: job 235 (result below).
- Relaunch of p3_aer2 Δt 0.5 with the fix: BLOCKED at 02:40 UTC — every A100 (16) and H100 (2) carries a compute process
  (ssh nvidia-smi --query-compute-apps); will be submitted with execution/submit_ena_lasso.sbatch to the first physically
  idle A100 (output runs/covert_p3_aer2_sscap_job<id>/), never touching the campaign directory.
- Five-run ENA-covert analysis vs ARM obs: analysis/ena_covert_analysis.jl, CPU job 236 → runs/ena_covert_analysis/,
  report docs/cases/ena_covert_analysis.md (pending the job).


## Bundle staged (coordinator, ARM order 284988, 2026-10-08 05:06 UTC) — LASSO work resumed
- data/lasso/20170718era5d25x100_sbmwrm-aer2-flxsst/{snd,lsf,sfc,prm,grd,bundle.toml,extra/}; archive sha256 b7594c2d…
  (53 739 520 B); samstat f5cb2bed… (94 426 591 B) and sam2d 4a4c6682… (991 972 348 B) in data/lasso/archives/.
  Bundle README: SAM v6.10.3 + LASSO mods, model_source_git_hash f83adf58 (not the public HEAD 12d0244 audited), stagesam 8404ad7.
- Member facts from prm/files: day0 = 199.0 (epoch 2017-07-18T00:00 UTC), dt = 1 s, nstop = 86400 (24 h), nrad = 30 (30 s),
  tauls = 7200 s, perturb_type = 5, read_in_geostrophic_wind = .true. (lsf 9 columns, 37 records 198.75–200.25, 431 height levels
  to 10.76 km), snd 3 records (199.0/199.5/200.0; 1077 HEIGHT levels to 10.76 km, p = -999.9 unused), sfc 37 samples (SST, ERA5
  H/LE for reference, TAU = -999.9), grd 260 levels 12.5…8087.5 m, caseid non-numeric → dimensions=(256,256,260) explicit.
  &SGS_TKE dosmagor = .false. (prognostic TKE), UNIFORM_SFC_FLX = .false. (per-column oceflx), fcor = 9.19626e-5, compute_reffc = .true.
- Decisions implemented (commit 162f402/b458839): (1) uniform_sfc_flx no longer rejected; .false. → adapter's per-column bulk
  fluxes with the new `surface_flux_law = :sam_oceflx` (Large & Pond drag exactly, Stanton 0.0327√cdn and Dalton 0.0346√cdn fitted
  to Breeze's (a₀+a₁U+a₂/U) form within 2 %, umin = 1 m/s, no gustiness; departures recorded in provenance: Breeze log-profile
  height shift + Li et al. stability vs SAM's two MO iterations and stable-branch Stanton 0.018√cdn); .true. recorded as a
  departure (uniform vs per-column). (2) namelist fcor accepted as the f-plane parameter (`coriolis_parameter`, FPlane(f=fcor)),
  recorded with the latitude-derived value as a warning. Also: bundle README metadata into provenance; LASSO default
  aerosol_supersaturation_cap = 0.003 (SBM ss_max emulation) and fully soluble SBM chemistry. Tests added (test/ena_lasso.jl,
  "SAM oceflx surface law" testset).
- Real-bundle validation: analysis/check_lasso_bundle.jl on the CPU partition (job 241; 239 failed on a precompile defect, fixed);
  Pkg.test job 242. Results recorded below when done.
- ENA-covert five-run analysis delivered: docs/cases/ena_covert_analysis.md, runs/ena_covert_analysis/ (CPU job 238). Verdict:
  1M and P3-N75 consistent with the case and within the obs envelope (LWP 160–169 vs 168 ± 13 g m⁻²; drizzle 0.002–0.003 vs
  0.0015 mm hr⁻¹; cloud top ~150 m high); P3-aer2 members unphysical in droplet number (≈ 400 cm⁻³ vs 75 observed); Δt 0.25 vs
  0.5 within chaotic spread; no numerical flaw.
- Real-bundle validation PASSED (CPU job 241): 0 problems, 8 warnings (fcor, 4 m level offset, TKE vs Smagorinsky, per-column
  oceflx, SS cap, compute_reffc, doseasons, diagCCN); tiny CPU cases (1M and capped P3-aer2 on the full 260-level column with the
  1077-level sounding / 431-level lsf) step finitely; capped aer2 modes 25.0e6 + 213.4e6 kg⁻¹ (10.8 % / 91.0 % activatable).
- LASSO 1-h short integration: Slurm job 243, GPU-c764cc9c (A100, gpu-p4de-2c-st-gpu-p4de-2c-2), source fa26811, output
  runs/lasso_aer2_1h_job243/, command = documented launch command with --hours 1 --profile_interval 600 --slice_interval 1800.
  Note: the bundle's fcor = 9.19626e-5 equals 2Ω sin φ with the SIDEREAL Ω (Oceananigans' FPlane value); the validator now accepts
  either convention silently (commit after fa26811; the 1-h job ran the 1–20 m/s oceflx fit, 2 % vs 0.4 % — immaterial).
- ENA-covert p3_aer2 Δt 0.5 relaunch with the aerosol fix (--aerosol_ss_cap 0.003, fully soluble): Slurm job 244 queued on
  gpu-p4de-2c (auto GPU selection at start), output runs/covert_p3_aer2_sscap_job244/.
- Pkg.test with the final adapter code: CPU job 248 (246 had one test-side expectation wrong; fixed) — result below.
- Job 243 timing: stepping from ~05:40 UTC, t = 30 min reached 06:02 → ≈ 0.7 s per 1-s step on the A100 (256²×260, P3-aer2,
  RRTMGP every 30 steps); the 24-h run needs ≈ 17 h wall → will be submitted with --time=22:00:00 right after the short run passes.
- SHORT INTEGRATION PASSED (job 243, 06:23 UTC): 3600 steps of 1 s in 44.1 min wall (0.69 s/step), finite throughout, cloud
  fraction 1.0 from t = 10 min, LWP 53 → 91 → 61 g m⁻² over the hour, max|w| ≤ 2.9 m s⁻¹, max qᶜˡ 1.6 g kg⁻¹, rain 0.002 mm d⁻¹,
  T ∈ [252.3, 295.6] K. Output runs/lasso_aer2_1h_job243/ (provenance.toml with bundle/readme/surface-law/fcor/SS-cap records).
  Acceptance checks (profiles of nᵃ/nᶜˡ/radiative flux divergence, SAM samstat comparison over the first hour): CPU job submitted
  (runs/lasso_short_checks_job243/).
- FULL 24-h RUN: Slurm job 249 (`blab-lasso-24h`, gpu-p4de-2c, auto GPU selection), queued since 06:25 UTC behind the other
  campaigns' jobs; output runs/lasso_aer2_24h_job249/; --time=22:00:00 (≈ 17 h expected once started).
- First-hour comparison with the SAM member (runs/lasso_aer2_1h_job243/checks/sam/, `compare_sam_reference.jl`): CWP SAM 72.8 vs
  Breeze 67.3 g m⁻² (−8 %), shaded cloud fraction 1.00 vs 0.98, PREC ≈ 0 in both, RWP 0.0004 vs 0.024 g m⁻²; SAM's GCSS ZCB/ZCT
  are ≈ 0 in hour 1 (undefined early), Breeze base 150 m / top 1312 m (the member starts with near-saturated air below the
  stratocumulus — both models carry very low cloud base in hour 1). MWRRET has no samples before 03:10 UTC.
- Pkg.test with the final adapter code (CPU job 248): BreezeLab 246/1 broken (pre-existing), protocols 13, member identity 20,
  namelist groups 7, grd 4, initial sounding 7, perturbation 4, bundle validation 69, adapter 71, RRTMGP CPU 5, SAM oceflx 19,
  ena_lasso wrapper 14, ARM readers 12, ARM download 16, manifest 4, ENA execution 132 — all passed.
- Plan step 7 tooling ready: analysis/compare_sam_reference.jl (Breeze vs samstat CWP/RWP/PREC/CLDSHD/ZCB/ZCT/SHF/LHF/LWNS/SWNS and
  θ/qᵛ/qᶜˡ/u/v profiles, plus MWRRET); samstat = 720 samples (2-min), sam2d = 288 (5-min) 2-D fields on day 199.0–200.0.


## In-flight jobs at hand-back (2026-10-08 08:05 UTC) and how to finish them
- Job 244 `blab-covert-aer2fix` (A100 GPU-c764cc9c, started 06:23, ≈ 5.3 h): ENA-covert p3_aer2 Δt 0.5 with --aerosol_ss_cap 0.003
  (provenance records maximum_supersaturation 0.003, mass_fraction_soluble 1, n₁ 2.48e7, n₂ 2.12e8 kg⁻¹). When COMPLETE:
  `sbatch execution/cpu_analysis.sbatch analysis/ena_aerosol_audit.jl runs/ena_aerosol_audit_fix runs/covert_p3_aer2_sscap_job244`
  and `... analysis/compare_ena_observations.jl runs/covert_p3_aer2_sscap_job244 /shared/home/greg/breezelab-runs/20261006/inputs/observations/20170718`
  → expect in-cloud nᶜˡ ≤ ≈ 290 cm⁻³ (was ≈ 400) and drizzle closer to the N75 members; append to docs/cases/ena_aerosol_audit.md.
- Job 249 `blab-lasso-24h` (pending, auto GPU, 22 h wall): the full LASSO member. When COMPLETE (runs/lasso_aer2_24h_job249/):
  `sbatch --wrap="bash execution/lasso_short_checks.sh runs/lasso_aer2_24h_job249" ...` (aerosol/radiation checks + samstat comparison)
  and `analysis/compare_ena_observations.jl runs/lasso_aer2_24h_job249 <obs dir>`; then the plan-step-7 report in docs/cases/ena_lasso.md.
  If the launcher exits 1 with "No physically idle GPU", resubmit the same sbatch line (it only runs on a GPU with no compute process).
- PR https://github.com/PRECISION-GM/BreezeLab.jl/pull/2 is current (branch ena-lasso at 4c089b8, merged with origin/main 293081f).


## Job 244 analysis (capped aerosol) and the Covert-consistent relaunch (2026-10-08 09:15 UTC)
- Job 244 (p3_aer2 + ss cap 0.003 + full solubility, Δt 0.5, 6 h, COMPLETE): in-cloud nᶜˡ 279 cm⁻³ (was 394 uncapped), i.e. at
  the SBM-activatable aer2 total (≈ 286); nᵃ+nᶜˡ conserved at 2.369e8 kg⁻¹; LWP mean 157 / 193 g m⁻² at 12 UTC (obs 168 ± 13);
  rain 2e-5 mm hr⁻¹ (obs 0.0015; N75 members 0.002–0.003) — drizzle still suppressed; base/top 795/1265 m vs obs ≈650/1090.
  Figures: runs/ena_aerosol_audit_sscap/, runs/ena_covert_analysis_sscap/; addendum in docs/cases/ena_covert_analysis.md.
- Cause (decided): the member choice — LASSO aer2 is not the Covert case aerosol. The Covert bulk paper prescribes the observed
  Nc = 75 cm⁻³; the public bin repo carries no aerosol spec (SBM spectrum compiled in; no bin ENA paper in Mechem's 2024 vita).
- Fix (commit 814dead, test fix 10ec1e3): `:p3_covert_n75` = LASSO mode shapes + SBM chemistry scaled so the SBM-capped
  activatable number is 75 cm⁻³ at the surface density (≈ 146 cm⁻³ total, factor ≈ 0.26), diagnostic-CCN projection; labelled
  as this package's configuration. eastern_north_atlantic(; microphysics=:p3_covert_n75) and the CLI accept it. Pkg.test (CPU job 264,
  09:40 UTC) with the new scheme: all pass — BreezeLab 290/1 pre-existing broken, ENA case execution 176 (now four schemes), LASSO
  testsets unchanged (69/71/5/19/14).
- Relaunch: Slurm job 263 `blab-covert-n75a` (gpu-p4de-2c, auto GPU, 6 h Δt 0.5, --aerosol_ss_cap 0.003), output
  runs/covert_p3_covert_n75_job263/. When COMPLETE: `sbatch execution/cpu_analysis.sbatch analysis/ena_aerosol_audit.jl
  runs/ena_aerosol_audit_n75 runs/covert_p3_covert_n75_job263` + the six-run overlay (analysis/ena_covert_analysis.jl with the
  new run); expect in-cloud nᶜˡ ≤ 75 cm⁻³ and drizzle like the N75 members.
- LASSO 24-h job 249 running on GPU-c764cc9c since ~08:40 UTC (undisturbed).

## Shared-file ownership announcements
- Will touch: `src/BreezeLab.jl` exports (new LASSO constructor/wrapper names), `cases/README.md`, `docs/cases/ena.md` cross-link,
  `data_wrangling/README.md` (new staging script entry), `test/runtests.jl` include list. No Project.toml deps added so far.

## Data order (coordinator, 2026-10-08 ~03:30 UTC)
- ARM LASSO-ENA order **284988** placed for member `20170718era5d25x100_sbmwrm-aer2-flxsst`: samin tar (53.74 MB), samstat .nc (94.43 MB), sam2d .nc (991.97 MB); FTP delivery, email link to follow. Coordinator polls the `armlive/saveData` route for HTTP 200 and will run the staging sequence above, then the short A100 integration and the full configured duration.
- 05:25 UTC (coordinator): order 284988 delivered (anonymous FTP, per-datastream subdirs, `.v0` suffix stripped). Files in data/lasso/archives/ with SHA256SUMS: samin 53 739 520 B b7594c2d…, samstat 94 426 591 B f5cb2bed…, sam2d 991 972 348 B 4a4c6682… (sizes match the order). README: SAM v6.10.3 + LASSO mods, model_source_git_hash f83adf58 (lasso_ena_noice), stagesam_ena 8404ad7. Staged with `stage_lasso_bundle.jl` → data/lasso/20170718era5d25x100_sbmwrm-aer2-flxsst/. `inspect_lasso_bundle`: day0 199.0, dt 1 s, nstop 86400, nrad 30, tauls 7200 s; **rejections to resolve before any run:** `uniform_sfc_flx = false` (per-column SAM fluxes, adapter implements uniform branch) and `fcor = 9.19626e-5` detached from latitude. Handed to the ENA-LASSO agent with the launch; no GPU idle at 05:30 UTC.
