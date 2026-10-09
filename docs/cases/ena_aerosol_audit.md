# Aerosol configuration audit of the ENA-covert P3-aer2 runs

Runs audited (Covert-public-bin development benchmark, source commit `3330f3c`, A100,
06–12 UTC 18 July 2017, 256×256×192 at 35 m): `p3_aer2_dt0.5_job168` and
`p3_aer2_dt0.25_job166` under `/shared/home/greg/breezelab-runs/docs-20261007-8dc824d/production/`.
Evidence: hourly-mean profiles of `nᵃ`, `nᶜˡ`, `qᶜˡ`, cloudy fraction and the slices,
read on a CPU node by `analysis/ena_aerosol_audit.jl`; figures and the numeric summary in
`runs/ena_aerosol_audit/` (`aerosol_profiles.png`, `slices_*.png`, `aerosol_audit.toml`).
Code paths: `src/case_setup.jl` (`lasso_aerosol_modes`, `build_microphysics`, the
`p3_initialization = :condensate_free` branch, the `DiagnosticCCNProjection` callback),
`src/large_scale_forcings.jl` (`DiagnosticCCNProjection`, `AerosolReplenishment`), Breeze
`5264d3c` `PredictedParticleProperties/aerosol_activation.jl` and
`cloud_droplet_activation_rates.jl`, and the LASSO SAM `MICRO_HUJISBM` (`12d0244`).

## Verdict

The aerosol fields are **not corrupted** (no NaN, collapse or blow-up; `nᵃ + nᶜˡ` is held
exactly at the initial total by the diagnostic-CCN projection; the two time steps agree to
three digits), but they are **not physically meaningful for this case**:

1. **Wrong aerosol for the case.** The member applies the LASSO-ENA "aer2" prescription
   (276 + 281 = 557 cm⁻³ at the surface) to the Covert et al. (2022) case, whose control
   droplet number is the *observed* `Nc = 75 cm⁻³` ("based on airborne in situ measurements
   of the 18 July 2017 case"; sensitivities 50 and 100 cm⁻³). The P3-aer2 member of the
   Covert benchmark is therefore a LASSO-aerosol sensitivity, not the Covert case; the
   case-consistent P3 member is `p3_n75`.
2. **Over-activation.** In cloud, `nᶜˡ ≈ 3.3–3.6×10⁸ kg⁻¹` (≈ 350–400 cm⁻³, 60–75 % of all
   aerosol, including part of the Aitken mode) throughout the six hours, 5× the observed
   75 cm⁻³, and above the ≈ 270 cm⁻³ the SBM itself can activate from aer2 (its
   `ss_max = 0.3 %` cap). Breeze's Morrison–Grabowski activation has no supersaturation
   cap and is a one-way relaxation (`ncnuc = max(0, N_act(S_grid) − nᶜˡ)/τ_act`,
   `τ_act = 1 s`) toward the equilibrium count at the *grid-cell* supersaturation, so every
   overshoot ratchets `nᶜˡ` up. The largest overshoot is the **initial burst**: with
   `p3_initialization = :condensate_free` the cloud layer starts with ≈ 0.55 g kg⁻¹ of
   excess vapor (S ≈ 6 %), which activates ≈ 390×10⁶ kg⁻¹ (≈ 430 cm⁻³) in the first
   seconds (first-hour mean peak, both runs). Later, the cloud-base/top supersaturation
   spikes keep `nᶜˡ` near 330–360×10⁶ kg⁻¹.
3. **Reservoir semantics.** `nᵃ` is initialized everywhere (free troposphere included) at
   the constant mixing ratio `4.618×10⁸ kg⁻¹` (SAM does the same: `FCCN0 = FCCNR·ρ(k)/ρ(1)`),
   stays exactly constant outside cloud, and is depleted to 140–200×10⁶ kg⁻¹ in the cloud
   layer. Because the projection rewrites `ρnᵃ ← ρ n₀ − ρnᶜˡ − ρnʳ` every step, the
   reservoir is a pure diagnostic of the droplet number: activation's own depletion and the
   advection of `ρnᵃ` have no lasting effect (no double counting, no drift), and CCN are
   "regenerated" wherever droplets evaporate or coalesce, as in the SBM's `diagCCN`. This
   is the intended emulation, but it also means the reservoir cannot limit the ratchet.

So the symptom Greg saw ("aerosol fields wrong") is real and expected from the
configuration: a 557 cm⁻³ LASSO aerosol applied to a 75 cm⁻³ case, activated by a scheme
without the SBM's supersaturation cap from a supersaturated condensate-free start.

## Evidence

`runs/ena_aerosol_audit/aerosol_audit.toml` (both runs; hour k = mean over hour k−1…k):

| hour | in-cloud `nᵃ` (10⁶ kg⁻¹) | in-cloud `nᶜˡ`/CF (10⁶ kg⁻¹) | `nᶜˡ + nᵃ` in cloud | `nᵃ` min / max | cloud base / top (m) | max `qᶜˡ` (g kg⁻¹) |
| --- | --- | --- | --- | --- | --- | --- |
| 1 (job168) | — (no cloud) | — | — | 461.8 / 461.8 | — | 0 |
| 2 | 144 | 402 | 462.3 | 71.7 / 461.8 | 765 / 1175 | 0.55 |
| 3 | 264 | 351 | 461.8 | 142 / 461.8 | 695 / 1225 | 0.56 |
| 4 | 226 | 355 | 461.8 | 135 / 461.8 | 745 / 1245 | 0.61 |
| 5 | 195 | 357 | 461.8 | 136 / 461.8 | 775 / 1255 | 0.64 |
| 6 | 186 | 335 | 461.8 | 138 / 461.8 | 795 / 1265 | 0.63 |
| 7 | 181 | 327 | 461.8 | 135 / 461.8 | 805 / 1275 | 0.64 |

job166 (Δt = 0.25 s) differs by < 3 % in every entry. Initial total `n₀ = 4.618×10⁸ kg⁻¹`
= 557 cm⁻³ / 1.206 kg m⁻³ (surface reference density), recorded in provenance as `n₁ + n₂`.
Time series: mean LWP 155 g m⁻² (job168) / 165 g m⁻² (job166), final 190 / 210 g m⁻²,
cloud fraction 1, rain 2–3×10⁻⁴ mm day⁻¹ (drizzle suppressed, consistent with 350–400 cm⁻³
droplets; the MWRRET 06–12 UTC mean is 168 ± 13 g m⁻² and the one-moment run gives 100 g m⁻²).

Units check: LASSO quotes cm⁻³; `lasso_aerosol_modes` converts with the surface reference
density to kg⁻¹ (`1e6 N/ρ₀`), which is SAM's constant-mixing-ratio initialization, and
Breeze's `AerosolMode.number_mixing_ratio` is per kg of air and `set!` writes
`ρnᵃ = ρ Σ modes` everywhere. Modal parameters (r = 0.018 / 0.066 μm, σ = 1.53 / 1.78) match
the LASSO methodology table 6. Chemistry: density 1790 kg m⁻³, M = 0.115 kg mol⁻¹,
ν = 3 match `micro_prm.f90` (`RO_SOLUTE = 1.79`, `mwaero = 115`, `ions = 3`).

## Specific misconfigurations (campaign source `3330f3c`)

| # | Where | What | Effect |
| --- | --- | --- | --- |
| 1 | `src/case_setup.jl` L60–72 `lasso_aerosol_modes`: `AerosolMode(...; chemistry..., kwargs...)` leaves Breeze's default `mass_fraction_soluble = 0.9` | The SBM solute is fully soluble (`micro_prm.f90` has no insoluble fraction) | `β_act` 10 % low → critical supersaturations 5 % high; minor |
| 2 | `src/case_setup.jl` L85–90 `build_microphysics` + Breeze `aerosol_activation.jl` L185–202 (`activation_supersaturation_threshold = 1e-6`, no maximum) | No counterpart of the SBM `ss_max = 0.003` activation cap; the Aitken mode (median s_c ≈ 0.73 %) activates in any cell with S ≳ 0.5 % | In-cloud `nᶜˡ` 350–400 cm⁻³ instead of ≤ 270 cm⁻³ |
| 3 | `src/case_setup.jl` L443–448 (`p3_initialization = :condensate_free`, the default) | Cloud layer starts ≈ 0.55 g kg⁻¹ supersaturated (S ≈ 6 %) and activates ≈ 430 cm⁻³ at once; the one-way activation then never releases them | First-hour `nᶜˡ` peak 390×10⁶ kg⁻¹; sets the level for the run |
| 4 | `src/ena_protocols.jl` / `execution/production_runs.sh` member choice `p3_aer2` for the Covert benchmark (`eastern_north_atlantic(; microphysics=:p3_aer2)`) | LASSO aer2 (557 cm⁻³) is not the Covert case aerosol (observed `Nc = 75 cm⁻³`) | The member answers a LASSO-aerosol question, not the Covert one |
| 5 | `src/large_scale_forcings.jl` L569–607 `DiagnosticCCNProjection` once per step (callback, L481–485 of `case_setup.jl`) | Correct emulation of `diagCCN` on the total number, but the SBM removes the *largest* CCN bins first, so its residual CCN are the hardest to activate; the number-only reservoir keeps the lognormal shape | Secondary; makes re-activation easier than in the SBM |

Not a misconfiguration: `CloudDroplets(number_concentration = 75e6)` passed alongside
`AerosolActivation` only sets the Liu–Daum shape parameter and the prescribed-number
fallback, which the prognostic path does not use (`process_rates.jl` L1055–1059);
`AerosolReplenishment` was not active (`aerosol_replenishment = diagnostic_ccn`).

## Fix on this branch

- `lasso_aerosol(; reference_density, mass_fraction_soluble = 1, maximum_supersaturation = nothing)`
  (at the time of the audit `lasso_aerosol_modes`):
  explicit full solubility (fixes #1) and an optional SBM-cap emulation that scales each
  mode's number to the fraction activatable below `maximum_supersaturation`
  (`activated_fraction`, evaluated with Breeze's own `activated_number`): for aer2 at
  0.3 % that is ≈ 8 % of mode 1 and ≈ 89 % of mode 2, i.e. ≈ 270 cm⁻³ activatable in total
  (fixes #2 as an explicit, recorded approximation: the lognormal shape is kept, so activation
  at low supersaturation is slightly underestimated). The cap is a keyword of the aerosol
  helper, `--aerosol_ss_cap` in the CLI, and the resulting mode numbers are recorded in
  provenance. `ena_lasso` applies it by default.
- #3 is left as the protocol default (the SBM also starts condensate-free), but with the
  cap the initial burst cannot exceed the SBM's activatable total; the alternative
  `initialization = :equilibrium, initial_droplet_number = 75e6` remains available.
- #4 is a labelling matter: the Covert benchmark's case-consistent P3 member is `p3_n75`;
  `p3_aer2` is documented here as a LASSO-aerosol sensitivity.
- A proper fix of #2 belongs in Breeze (`AerosolActivation(; maximum_supersaturation)`
  applied inside `activated_number`), which would also bound the in-cloud ratchet.

Relaunch (plan: the affected configuration `p3_aer2`, Δt = 0.5 s, with the cap and full
solubility, on a physically idle A100 under `runs/`): see `plans/status/ENA-LASSO.md`
for the job ID and output path once submitted.
