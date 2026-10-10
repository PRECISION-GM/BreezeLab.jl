# Moving BreezeLab from Breeze 5264d3c (0.11.3 line) to Breeze 0.12.0

Status: October 2026, branch `breeze-0.12`. This note maps every Breeze change between the old pin and
the 0.12.0 release onto BreezeLab code, records the environment split forced by NumericalEarth, and
reports a CPU before/after physics check. It is the reference for deciding which reruns use 0.12.

## What changes in the environments

| environment | before | after |
| --- | --- | --- |
| package (`Project.toml`, `Manifest.toml`), used by ENA Covert, ENA LASSO, TRACER–DP-SCREAM, SEA STARR and `docs/` | Breeze `5264d3c` via `[sources]` (0.11.3 line: PR 959 head + four ENA fixes), compat `0.11` | registered Breeze **0.12.0** (tag `aded2b9`, commit `740c693`), compat `0.11.3, 0.12` (see below) |
| Oceananigans (both environments) | `847e125` via `[sources]` (0.113.3 + Float32 WENO-Z cap + limiter 0/0 fix) | unchanged; Breeze 0.12 needs Oceananigans ≥ 0.113.2, which `847e125` (0.113.3) satisfies |
| CloudMicrophysics | 0.41.0 | 0.41.0 (Breeze 0.12 allows 0.41–0.43; the bump is left for a separate change so that no physics difference here comes from CloudMicrophysics) |
| `cases/tracer_mip` sub-environment (NumericalEarth) | Breeze `5264d3c`, NumericalEarth `d07eb240` | **unchanged** (see *NumericalEarth*) |

Only the Breeze entry of `Manifest.toml` changes (plus the package's own entry, which the old manifest
recorded without its weak dependency).

### NumericalEarth blocks 0.12 for TRACER-MIP

NumericalEarth declares Breeze as a weak dependency. Every registered NumericalEarth version from 0.8.0
to 0.8.2 (`WeakCompat.toml` in General, `["0.8 - 0"] Breeze = "0.11"`), NumericalEarth `main`
(`02eaa19a`, 2026-10-09, `Breeze = "0.11"`) and the pinned `d07eb240` all require **Breeze 0.11**. No
NumericalEarth version allows Breeze 0.12. The `cases/tracer_mip` environment therefore stays on Breeze
`5264d3c`. To keep it resolvable, the package compat is `Breeze = "0.11.3, 0.12"`; BreezeLab's source
needs no version branches (no API change in 0.12 reaches the code BreezeLab uses, see the table), and the
full test suite runs in both environments (results below).

NumericalEarth PR 753 (open, branch `glw/bottom-precipitation-flux`) resolves `surface_precipitation_flux`
or `bottom_precipitation_flux`, whichever the loaded Breeze defines; a NumericalEarth release with a
`Breeze = "0.11, 0.12"` weak compat would let TRACER-MIP follow. Until then:

- the rain-to-land shim in `ext/BreezeLabNumericalEarthExt.jl` stays: NumericalEarth (all versions,
  `main` included) still imports `surface_precipitation_flux`, which Breeze removed in PR 959 (in the
  pin and in 0.12);
- TRACER-MIP cannot use PR 1036 (tuples of closures: horizontal Smagorinsky + TKE) — a follow-up once
  NumericalEarth allows 0.12.

### Three ENA fixes that 0.12.0 does not contain

The old pin was not a release: it was the head of PR 959 plus five commits on the Breeze branch
`glw/lasso-ena-pr959`. One of them reached Breeze `main` in another form, three did not:

| pin commit | fix | in 0.12.0? |
| --- | --- | --- |
| `16c8e0e`, `304bef8` | P3: `cloud_number_per_cloud_mass` returns 0 below `minimum_mass_mixing_ratio` (subnormal cloud mass ⇒ `Nᶜˡ/(ρqᶜˡ)` = Inf ⇒ `Inf × 0` = NaN in `ρnᶜˡ`; P3-aer2 went NaN after 8–9 steps) | **no** (0.12.0 still divides with `safe_divide`) |
| `14da778` | P3: prescribed-Nᶜˡ path re-diagnoses the number against the residual cloud mass before homogeneous freezing (ENA 6 km P3-N75 went NaN on its third step) | **no** |
| `b633386` | `AtmosphereModel` moves the microphysics to the grid architecture *before* the boundary conditions capture it (θ-formulation energy-flux conditions carry the microphysics into halo kernels; P3's host lookup tables made the GPU launch fail with a non-bitstype argument) | **no** (0.12.0 still calls `on_architecture(arch, microphysics)` after materializing the boundary conditions) |
| `5264d3c` | microphysical terminal velocities in the advective time-step limit | yes, rewritten in PR 997 (`maximum_vertical_transport_speed`) |
| `29bcc41` | comment style | — |

Consequences for registered 0.12.0: GPU runs with P3 and the potential-temperature formulation with an
energy-flux boundary condition (every BreezeLab P3 case) are expected to fail at the first halo fill, and
P3 with a prognostic aerosol reservoir (LASSO default `aer2`, SEA STARR) loses the NaN guard. The five
pin-only commits other than `5264d3c` apply cleanly onto `v0.12.0` (`git cherry-pick 16c8e0e9 304bef8c 14da7781 b6333861 29bcc419`;
the regression tests they add come with them). **Before GPU reruns move to 0.12, that branch has to be
published (e.g. `glw/lasso-ena-0.12.0` on NumericalEarth/Breeze.jl) and pinned in `[sources]`**, or the
fixes have to be released upstream. This change could not publish the branch itself.

## Commit map: `5264d3c..v0.12.0` (25 commits)

The pin and 0.12.0 diverge at `d0afaf4` (PR 1027): the pin carries PR 959 up to `cc004bc`, while 0.12.0
carries the squash-merged PR 959 and 24 other commits.

| Breeze PR (commit) | change | BreezeLab code it reaches | effect here |
| --- | --- | --- | --- |
| #959 (`0e3d153`) unify sedimentation transport, condensate thermal coupling | the pin already contains the branch up to `cc004bc` (shared sedimentation interface, enthalpy transport, implicit remainder after the solve, `bottom_precipitation_flux` rename). New relative to the pin: `49c9381` (θ content reads the diagnosed T instead of re-inverting θ: roundoff), `3c69e30` (sedimentation tendency in its own kernel: performance, GPU register pressure), `dbb55e8` (`condensate_phase` → `condensate_liquid_fraction`), tests and a benchmark | none of the renamed internals are used by BreezeLab | roundoff only; faster thermodynamic tendency kernel on GPU |
| #997 (`ebdf116`) moisture partitioning, unified moist-state conversions | saturation adjustment keeps prognostic precipitation in the partition (fixes Breeze issue 1013: 1M cloud liquid was diagnosed as adjusted liquid *minus* rain); `wall_potential_temperature` uses the near-wall moist composition (bulk sensible heat flux in the θ formulation was biased by a dry Exner function); `set!` converts `qᵗ` for every scheme; `filter_timescale = Inf` samples instead of freezing filtered surface fields; sedimentation speeds in the CFL; `adjust_thermodynamic_state` gains an optional precipitation argument | ENA one-moment members (`cloud_formation = SaturationAdjustment`); `bulk_surface_flux_boundary_conditions` (SEA STARR `surface = :bulk_sst`, LASSO `:sam_oceflx`/`:breeze`) and DP-SCREAM's `BulkDrag` (identical over the 30-minute check); `initial_conditions.jl` calls `adjust_thermodynamic_state(𝒰, SaturationAdjustment, constants)` (still valid: the new argument defaults to no precipitation). BreezeLab never passes `filtered_velocities`, so the filter change does not apply | physics change for 1M (cloud/rain partition); bulk sensible heat flux change < 0.1 W m⁻² over the ocean (see *Attribution*) |
| #1031 (`f33703a`) WP1M/MP1M draw precipitation from `ρqᵉ` (Breeze issue 1012) | autoconversion, accretion, rain evaporation and the negative-rain guard changed `ρqʳ` without the opposite change in `ρqᵉ`: water was created (or destroyed) | ENA Covert/LASSO one-moment members | physics change (water conservation) |
| #1006 (`925af9d`) boundary conditions on microphysical prognostics | P3, 1M and 2M fields are rebuilt with the boundary conditions given for them (they were silently dropped) | BreezeLab gives none today. SEA STARR's `SurfaceAerosolSource` (bottom-cell forcing, written because P3 ignored a flux condition on `ρnᵃ`) can become `FluxBoundaryCondition(7e5)` at the bottom of `ρnᵃ`: the bottom-cell tendency `F/Δz₁` is the same in exact arithmetic | none until adopted (follow-up, below) |
| #1036 (`cdacc35`) tuples of turbulence closures | `closure = (c₁, c₂)` | TRACER-MIP outer domain (TKE only today) | not usable while TRACER-MIP stays on 0.11 (follow-up) |
| #1047 (`fa81035`), #1050 (`033489e`), #872 (`0469b01`) | two-moment aerosol activation reads `air_pressure` (unchanged on the anelastic core); density-closed states allowed; reference-pressure saturation adjustment restricted to states that have a reference pressure | BreezeLab uses P3's own activation (and `KappaAerosolMode` extends P3's `aerosol_activation_rate`/`activated_number`, whose signatures are unchanged); no 2M scheme | none |
| #1060 (`7545908`) refresh the anelastic pressure solver when ρᵣ changes | `set_to_mean!` / `set!(…; compute_reference_state=true)` / `compute_reference_state!` now refresh the tridiagonal coefficients | no anelastic BreezeLab case rewrites the reference state after construction; TRACER-MIP calls `set!(nest; …, compute_reference_state=true)` on its NumericalEarth child, which is compressible | none |
| #1048 (`fe0951f`) non-precipitating `BulkMicrophysics` steps | method fixes | not used | none |
| #1045 (`3286991`) split-explicit tendencies and `TendencyCallsite` callbacks | `AcousticRungeKutta3` builds the slow tendencies once (about 33 % less time per step in the PR's CPU case) and `TendencyCallsite` callbacks reach the stage | compressible split-explicit only: the TRACER-MIP child, which stays on 0.11; BreezeLab registers no `TendencyCallsite` callbacks | none now; a speed-up for TRACER-MIP once it can move |
| #998 (`bd3d6e2`) RRTMGP column batching | opt-in `column_batches = n` (default `nothing`: unchanged path) | `RadiativeTransferModel` in SEA STARR, DP-SCREAM, LASSO | none by default; an option to cut RRTMGP workspace memory on large GPU domains |
| #1017 (`4c0f6a3`) | `clamp_negative_numbers!` as tuple recursion (Reactant) | — | none |
| #1040 (`720e477`) | CloudMicrophysics compat 0.41–0.43 | — | none (kept at 0.41.0) |
| #1032, #1033, #1037, #1038, #1039, #1041, #1043, #1044, #1046, #1058, #1061 | CI, benchmarks, examples, Reactant, tests, version | — | none |

### API breaks

One, from PR 998: `RadiativeTransferModel` gained a `column_batches` field before `schedule`. SEA STARR's
`extended_column_radiation` builds the struct positionally (to append the driver's upper atmosphere to
the RRTMGP columns) and failed with a `MethodError`; it now passes `nothing` (one unbatched solve, which
keeps the solver's own extended state) when the field exists, so it also still builds on Breeze 0.11 in
`cases/tracer_mip`. Provenance needed a change of its own: it read the Breeze revision from a
`[sources]` pin, so it now records a registered dependency as `registered v0.12.0 (git-tree-sha1 …)` and
reads only a Manifest entry whose version is the one loaded (in `cases/tracer_mip` that is the
sub-environment's `5264d3c`). Everything else BreezeLab imports or extends keeps its call form:
`materialize_atmosphere_model_forcing`, `compute_forcing!`, `_update_radiation!`, P3's
`activated_number`/`aerosol_activation_rate`/`compute_p3_process_rates`/`p3_process_properties`/
`p3_adiabatic_temperature_tendency`, `grid_microphysical_state`, `diagnose_thermodynamic_state`,
`maybe_adjust_thermodynamic_state` (0.12 adds a six-argument form that falls back to the four-argument
one BreezeLab calls), `grid_moisture_fractions`, `adjust_thermodynamic_state` (new optional argument),
`moisture_specific_name`, `prognostic_field_names`, `bottom_precipitation_flux`, the RRTMGP extension
helpers, and the forcing and boundary-condition constructors.

## Follow-ups deliberately not in the version bump

- **ENA fixes branch** (above): required before GPU P3 runs on 0.12.
- **SEA STARR surface aerosol source as a boundary condition** (#1006): replace `SurfaceAerosolSource`
  under `nᵃ` with `boundary_conditions = (; ρnᵃ = FieldBoundaryConditions(bottom = FluxBoundaryCondition(7e5)))`
  merged with the bulk conditions, and drop the corresponding entry from the run's `departures`. Same
  tendency in exact arithmetic; it changes the provenance and the code path, so it belongs to the SEA
  STARR owner and a CPU regression against the current forcing.
- **TRACER-MIP composite closure** (#1036): a tuple of a horizontal Smagorinsky closure and `TKEBasedTurbulenceClosure()`
  once NumericalEarth allows Breeze 0.12.
- **CloudMicrophysics 0.42/0.43**: a separate bump with its own regression check.
- **RRTMGP column batching** (#998) for the 512² DP-SCREAM and 384² SEA STARR grids if GPU memory is tight.

## Before/after physics check (CPU)

`analysis/breeze_upgrade_physics_check.jl` builds each case with identical settings in the old
environment (`main` at `afb84a9`: Breeze `5264d3c`) and the new one (registered Breeze 0.12.0), steps it
with a fixed Δt (adaptive stepping removed so both take the same steps) and records domain means;
`analysis/breeze_upgrade_physics_compare.jl` tabulates the two. Float32 as in production, 16 × 16
columns at the protocol Δx (100 m for ENA, 50 m SEA STARR, 200 m DP-SCREAM), the full protocol
vertical grids, CPU partition (AMD EPYC 7R13). Cases: ENA Covert with the one-moment member and P3-N75,
ENA LASSO with its protocol default (P3 with the capped `aer2` prognostic aerosol, RRTMGP), TRACER–DP-SCREAM
(P3 N = 200 cm⁻³, RRTMGP), SEA STARR CTRL (P3 with κ-Köhler reservoir, bulk SST fluxes, extended-column
RRTMGP). Two closed boxes isolate the microphysics: a resting, horizontally uniform cloudy layer
(LWP ≈ 450 g m⁻²) with no forcing, surface flux or radiation, once with the one-moment scheme and once
with P3 (the P3 box also starts with rain in its lowest 200 m).

### 50 steps

| case | quantity | old (5264d3c) | new (0.12.0) | new − old |
| --- | --- | ---: | ---: | ---: |
| box, one-moment | LWP / RWP [g m⁻²] | 437.291 / 9.0832 | 437.264 / 9.0784 | −0.027 / −0.005 |
| | change of prognostic water ∫Σρq dz [g m⁻²] | **+9.09** | **+0.001** | −9.09 |
| | change of the diagnosed water ∫ρ(qᵛ+qᶜˡ+qʳ) dz [g m⁻²] | 0.000 | 0.000 | 0 |
| box, P3 | all quantities | | | bitwise identical |
| | rain through the bottom [mm d⁻¹]: `bottom_precipitation_flux` / `surface_rain_flux` | +32.05 / +31.98 | +32.05 / +31.98 | 0 |
| Covert one-moment | LWP / RWP [g m⁻²] | 102.313 / 0.1537 | 102.316 / 0.1535 | +0.004 / −0.0002 |
| Covert P3-N75 | LWP [g m⁻²] | 114.9243 | 114.9242 | −9e-5 |
| LASSO default (P3 aer2) | LWP [g m⁻²], Nᶜˡ in cloud [cm⁻³] | 42.3318, 247.79 | 42.3314, 247.79 | −4e-4, +5e-5 |
| DP-SCREAM | all quantities (IWP 0.73 g m⁻², no liquid yet) | | | identical to printed precision |
| SEA STARR CTRL | LWP [g m⁻²], Nᶜˡ [cm⁻³] | 111.8770, 118.122 | 111.8769, 118.123 | −4e-5, +1e-4 |

Domain-mean T, θ and qᵛ agree to ≤ 1e-4 K and ≤ 2e-6 g kg⁻¹ in every case.

### 30 minutes (1800 steps; DP-SCREAM 900 steps of 2 s)

| case | LWP [g m⁻²] old → new | RWP [g m⁻²] old → new | other |
| --- | --- | --- | --- |
| Covert one-moment | 91.11 → 90.08 (−1.1 %) | 0.988 → 0.825 (−17 %) | ∫ρθ dz change −3.47 → −5.14 kg K m⁻² (rain evaporation now cools); prognostic water +1.0 g m⁻² higher in old |
| Covert P3-N75 | 110.66 → 110.42 (−0.2 %) | 1.070 → 1.059 | surface rain 0.0137 → 0.0102 mm d⁻¹ (16² domain, chaotic); water budget equal to 2e-5 kg m⁻² |
| LASSO default (P3 aer2) | 82.59 → 82.42 (−0.2 %) | 0.0114 → 0.0114 | Nᶜˡ 220.6 → 220.6 cm⁻³; no NaN in either |
| DP-SCREAM | 0 → 0 | 0 → 0 | IWP 15.61 → 15.61 g m⁻², T, θ, qᵛ, ∫ρθ dz identical to printed precision |
| SEA STARR CTRL | 114.56 → 114.88 (+0.3 %) | 0.581 → 0.578 | Nᶜˡ 98.9 → 99.7 cm⁻³; surface rain 6.8e-4 → 9.9e-4 mm d⁻¹; no NaN in either |

### Attribution

- **One-moment members: Breeze #1031 and #997 (real physics change).** In the old code the
  saturation-adjustment one-moment scheme gave `ρqᵉ` no tendency while autoconversion, accretion and
  rain evaporation changed `ρqʳ` (issue 1012), and the adjustment diagnosed cloud as adjusted liquid
  minus rain (issue 1013). The two errors cancel in every *diagnosed* field — the box's diagnosed water
  is constant and its LWP and RWP agree with 0.12 to 0.01 % — but the prognostic water grows by exactly
  the rain produced (9.1 g m⁻² in 50 s in the box), and evaporating rain neither moistens nor cools the
  sub-cloud layer. In the 30-minute Covert run this shows as 17 % less rain water, 1 % less LWP and a
  stronger evaporative cooling of the column on 0.12. Any one-moment run from the old pin carries this
  water source; its size scales with the rain production.
- **P3 members: roundoff only.** The P3 box is bitwise identical, DP-SCREAM is identical to printed
  precision for 30 minutes (no liquid yet), and the P3 cases with liquid differ by 1e-7–1e-5
  (relative) after 50 steps (#959's `49c9381`, which reads the diagnosed temperature in the
  θ condensate content instead of re-inverting θ, and floating-point reordering in #997's moist
  conversions). The turbulence on a 16² domain amplifies that seed to ±0.3 % in LWP, ±0.7 % in Nᶜˡ and
  tens of percent in the tiny surface rain after 30 minutes — the chaotic noise floor, not a model
  change. Water budgets agree to ≤ 6e-5 kg m⁻².
- **Bulk surface fluxes (#997 wall θ):** the moist instead of dry Exner function at the wall changes
  θ₀ by about θ₀ · Δκ · ln(p₀/pˢᵗ) ≈ 290 K × 7e-4 × 0.015 ≈ 0.003 K over the ocean (qᵛ ≈ 0.01,
  p₀ ≈ 1015 hPa), i.e. < 0.1 W m⁻² of sensible heat flux. Not detectable in SEA STARR or LASSO above
  the noise floor.
- **Missing ENA fixes:** neither the LASSO `aer2` nor the SEA STARR run went NaN in 30 minutes on
  the CPU without the subnormal-cloud guard; the guard protects against a state that occurs in the
  full-size GPU runs, and the GPU energy-flux/P3 failure cannot show on the CPU. Neither result
  removes the need for the fix branch.

### Sign of the rain-to-land shim (TRACER-MIP, not changed here)

Breeze's `bottom_precipitation_flux` is positive **downward** (docstring, and the P3 box: +32.05 mm d⁻¹
against BreezeLab's `surface_rain_flux` +31.98 mm d⁻¹ for the same falling rain, on 5264d3c and on
0.12.0 alike). NumericalEarth reads `Jʳⁿ` as positive downward too. The shim in
`ext/BreezeLabNumericalEarthExt.jl` multiplies by −1 on the assumption that Breeze's flux is positive
upward, so rain reaching the TRACER-MIP slab land enters with the wrong sign. The regional test only
checks `Jʳⁿ ≥ 0` without rain. This is a TRACER-MIP physics fix (drop the `-1 *`, add a raining-column
test) to make separately; it is not part of this version bump.

## Recommendations for the queued corrected reruns

All queued corrected ENA, LASSO, DP-SCREAM and SEA STARR runs use P3, for which 0.12.0 changes nothing
beyond roundoff, while registered 0.12.0 lacks the GPU and NaN fixes the pin carries:

- 532–535 (Covert P3-N75 base, inversion 5 m, no closure, RRTMGP LW), 536–537 (LASSO aer2 24 h,
  r_eff 14 μm): **keep on the old pin** (their `ena-fixed` checkout). On registered 0.12.0 they would
  fail at the first halo fill on the GPU (missing `b633386`), and there is no physics gain for P3.
- 514–515 (DP-SCREAM 200 m / 100 m): **keep on the old pin**; identical physics over 30 minutes.
- 491–493 (SEA STARR resumes): must continue with the code they started from (coordination rule).
- 501 (TRACER-MIP outer domain): cannot move (NumericalEarth requires Breeze 0.11).
- Move to 0.12 only after the ENA-fix branch on top of 0.12.0 is published and pinned; first
  candidates are new runs that use the **one-moment** scheme (the five-run ENA campaign's 1M member):
  they gain water conservation and rain-evaporation thermodynamics (#1031, #997).

## Tests

All on the CPU partition (`cpu-c6a24xlarge`, AMD EPYC 7R13), Julia 1.12.7, own depot first in
`JULIA_DEPOT_PATH` (new compile caches for 0.12):

| environment | command | result |
| --- | --- | --- |
| package, Breeze 0.12.0 (this branch) | `Pkg.test()` | job 571 (commit `9ba0209`): 1183 passed, 1 broken (the existing `@test_broken`), 0 failed, 1 h 19 min |
| package, old pin (main `afb84a9`, baseline) | `Pkg.test()` | job 551: 1183 passed, 1 broken (the same `@test_broken`), 1 h 22 min |
| `cases/tracer_mip`, Breeze 5264d3c + NumericalEarth d07eb240 + this branch's source | `julia --project=cases/tracer_mip test/runtests.jl` (includes the 29 regional machinery tests) | job 572 (commit `9ba0209`): all testsets pass, including the 29 regional tests; the sub-environment Manifest is unchanged by `Pkg.resolve()` |
| `docs` (derived by `docs/setup.jl`) | instantiate, precompile, load BreezeLab/Breeze/Documenter/Literate | job 564: loads with Breeze 0.12.0 |

The first test run on 0.12.0 failed only in the provenance test, which expected a `[sources]` pin for
Breeze; the SEA STARR physics check then found the `RadiativeTransferModel` field added by PR 998. Both
are fixed on this branch. Physics-check output: `/shared/home/greg/breezelab-work/breeze012-runs/physics`
(50 steps) and `physics_long` (30 minutes), each with the environment's `Manifest.toml`.
