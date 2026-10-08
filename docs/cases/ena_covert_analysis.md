# ENA-covert campaign (docs-20261007-8dc824d): five completed runs vs ARM observations

Covert-public-bin development benchmark, 18 July 2017 06–12 UTC, 256×256×192 at 35 m
(8.96 km), reconstructed 192-level grid, prescribed surface fluxes, simple longwave
radiation, no wind nudging; source `3330f3c`; A100/H100. **Not LASSO.** Runs:

| label | job | microphysics | Δt (s) | output |
| --- | --- | --- | --- | --- |
| one_moment_dt0.5 | 167 | one-moment (saturation adjustment) | 0.5 | `production/one_moment_dt0.5_job167` |
| p3_n75_dt0.5_docs165 | 165 | P3, prescribed Nᶜˡ = 75 cm⁻³ | 0.5 | `source/output/documentation/ena` |
| p3_n75_dt0.25 | 169 | P3, Nᶜˡ = 75 cm⁻³ | 0.25 | `production/p3_n75_dt0.25_job169` |
| p3_aer2_dt0.5 | 168 | P3, LASSO aer2 aerosol, diagnostic CCN | 0.5 | `production/p3_aer2_dt0.5_job168` |
| p3_aer2_dt0.25 | 166 | P3, LASSO aer2 aerosol, diagnostic CCN | 0.25 | `production/p3_aer2_dt0.25_job166` |

Tooling: `analysis/ena_covert_analysis.jl` (CPU node) calls `analysis/compare_ena_observations.jl`
per run and draws the overlays; outputs in `runs/ena_covert_analysis/` (`timeseries_all_runs.png`,
`profiles_all_runs.png`, `summary.toml`, `comparison_<run>/`). Observations:
MWRRET `phys_lwp` (good retrievals, 1-σ ≈ 13 g m⁻²), VDIS `rain_rate`, ARSCL lowest-layer top
and ceilometer/MPL `cloud_base_best_estimate`, all 06–12 UTC.

## Results (06–12 UTC means unless stated)

| quantity | obs | one_moment | n75 Δt 0.5 | n75 Δt 0.25 | aer2 Δt 0.5 | aer2 Δt 0.25 |
| --- | --- | --- | --- | --- | --- | --- |
| LWP (g m⁻²), mean | 168 ± 13 (median 158, q10–q90 56–292) | 100 | 160 | 169 | 155 | 165 |
| LWP at 12 UTC | ≈ 120 (noisy) | 113 | 195 | 208 | 190 | 210 |
| cloud fraction | 1 (zenith) | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |
| rain (mm hr⁻¹), mean | 0.0015 (intermittent, max 0.39) | 0.0029 | 0.0023 | 0.0031 | 1×10⁻⁵ | 1×10⁻⁵ |
| cloud top (m), 12 UTC | ≈ 1090 (ARSCL, 1000–1250) | 1265 | 1245 | 1245 | 1275 | 1275 |
| cloud base (m), 12 UTC | ≈ 600–750 (ceilometer) | 735 | 765 | 745 | 805 | 775 |
| in-cloud Nᶜˡ (cm⁻³) | 75 (airborne, Covert et al. 2022) | — | 75 (prescribed) | 75 | 394 | 425 |

- **LWP.** The P3 members reproduce the observed six-hour mean (160–169 vs 168 g m⁻²) but
  rise monotonically to ≈ 200 g m⁻² while the observations fluctuate between 0 and 590 g m⁻²
  about a flat mean; the one-moment control is ≈ 40 % too thin. The domain mean cannot show
  the zenith variability; the comparison is a mean/quantile one.
- **Drizzle.** The 75 cm⁻³ members drizzle at 0.002–0.003 mm hr⁻¹, the same order as the
  disdrometer mean (0.0015, three showers). The aer2 members suppress drizzle by 300×
  (`nᶜˡ` ≈ 400 cm⁻³, see below), which is unphysical for this day.
- **Cloud boundaries.** All runs deepen the layer: the inversion rises from 1130 to
  ≈ 1250–1275 m while ARSCL tops stay near 1000–1250 m (mean 1087 m), and the LES base
  (735–805 m) sits 50–150 m above the ceilometer base (600–750 m). The first-hour profiles
  show the initial cloud at 800–1140 m; the deepening is a model response to the large-scale
  forcing/surface fluxes of the public configuration, not a Δt effect.
- **Δt sensitivity.** Δt = 0.25 vs 0.5 s changes the six-hour LWP mean by 6 % (n75) and
  6 % (aer2), the final LWP by 7–10 %, cloud base by ±30 m, rain by 30 %: within the
  chaotic spread, with no sign of a time-step instability (the earlier P3 runaways are gone).
- **Aerosol (the flaw).** `nᵃ` is uniform at 557 cm⁻³-equivalent (4.62×10⁸ kg⁻¹) outside
  cloud and depleted to ≈ 220 cm⁻³ in cloud while `nᶜˡ` ≈ 390–430 cm⁻³ (70–75 % of all
  aerosol, Aitken mode included), 5× the observed droplet number. Diagnosis and fix in
  [the aerosol audit](ena_aerosol_audit.md): LASSO aer2 (557 cm⁻³) applied to a 75 cm⁻³
  case, no supersaturation cap in Breeze's activation (the SBM has `ss_max = 0.3 %`), a
  condensate-free start supersaturated by ≈ 6 %, and Breeze's default soluble fraction.
  The relaunch of `p3_aer2` Δt = 0.5 s with the SBM-cap emulation and full solubility
  (`--aerosol_ss_cap 0.003`) is recorded in the status file (job ID, `runs/` path) when a
  GPU is physically idle.
- **Other departures to keep in mind** (not flaws of these runs): the public configuration
  is the 8.96 km / 06–12 UTC sensitivity domain of the paper, not its 30.24 km / 06–15 UTC
  control; the vertical grid is a 10-m reconstruction where the paper used 5 m near the
  surface and inversion; surface fluxes and `rad_simple` longwave follow the public namelist.

## Verdict

The one-moment and P3-N75 runs are physically consistent with the case definition and
within the observational envelope for LWP, drizzle and cloud fraction, with a ≈ 150 m
high bias in cloud top/base that the public forcing produces in all members. The P3-aer2
runs are a LASSO-aerosol sensitivity whose droplet number is unphysical for 18 July 2017;
they should not be read as the Covert case. No numerical flaw (NaN, runaway, Δt
dependence) was found in any of the five.

## Addendum (8 October 2026): the capped-aerosol relaunch, job 244

`runs/covert_p3_aer2_sscap_job244/` — `p3_aer2`, Δt = 0.5 s, with the aerosol-audit fix
(`--aerosol_ss_cap 0.003`, fully soluble SBM chemistry; provenance: `maximum_supersaturation
= 0.003`, `n₁ + n₂ = 2.37×10⁸ kg⁻¹` ≈ 286 cm⁻³ activatable), on A100 `GPU-c764cc9c`,
43 200 steps, finite. Analysis: `runs/ena_aerosol_audit_sscap/` and
`runs/ena_covert_analysis_sscap/` (six-run overlay including this run).

| quantity (06–12 UTC) | obs | p3_aer2 (168, uncapped) | p3_aer2 + cap (244) | p3_n75 (165) |
| --- | --- | --- | --- | --- |
| in-cloud nᶜˡ, last hour (cm⁻³) | 75 | 394 | **279** | 75 (prescribed) |
| in-cloud nᵃ, last hour (cm⁻³) | — | 219 | 49 | — |
| nᵃ + nᶜˡ (10⁸ kg⁻¹, conserved) | — | 4.618 | 2.369 | — |
| LWP mean / at 12 UTC (g m⁻²) | 168 ± 13 / ≈ 120 | 155 / 190 | 157 / 193 | 160 / 195 |
| rain mean / max (mm hr⁻¹) | 0.0015 / 0.39 | 1×10⁻⁵ / 2×10⁻⁵ | 2×10⁻⁵ / 1×10⁻⁴ | 0.0023 / 0.0037 |
| cloud base / top at 12 UTC (m) | ≈ 650 / ≈ 1090 | 805 / 1275 | 795 / 1265 | 765 / 1245 |

The cap works as designed: the in-cloud droplet number falls from 394 to 279 cm⁻³, i.e. to
the SBM's activatable aer2 total (≈ 286 cm⁻³; the diagnostic-CCN projection keeps
`nᵃ + nᶜˡ` at that total, so nearly the whole reservoir is activated in cloud), and the
SAM-side expectation for aer2 is met. But for the **Covert case** the droplet number is
still 3.7× the observed 75 cm⁻³ and drizzle is still suppressed (2×10⁻⁵ vs 0.002–0.003
mm hr⁻¹ for the 75 cm⁻³ members and 0.0015 observed); LWP and cloud boundaries are
indistinguishable from the other P3 members. **Remaining cause: the member choice** —
LASSO aer2 (557 cm⁻³ total, ≈ 286 cm⁻³ activatable) is not the aerosol of the Covert case.

What the Covert case prescribes: the bulk paper (Covert, Mechem & Zhang 2022, ACP 22, 1159)
uses the observed `Nc = 75 cm⁻³` ("based on airborne in situ measurements of the 18 July
2017 case"; sensitivities 50 and 100 cm⁻³); the public bin repository carries only
`snd/lsf/sfc/prm` with no aerosol specification (the HUJI-SBM spectrum is compiled into
SAM, and no bin-microphysics ENA paper by these authors exists in Mechem's 2024
publication list). A Covert-consistent aerosol therefore has to be *built*: the new
`:p3_covert_n75` member (`covert_aerosol_modes`) keeps the LASSO mode shapes and SBM
chemistry and scales their number so that the SBM-capped activatable total is exactly
75 cm⁻³ at the surface density (≈ 146 cm⁻³ total aerosol, scale factor ≈ 0.26), with the
diagnostic-CCN projection; it is labelled as this package's configuration, not as a
Covert or ARM prescription. Relaunch: Slurm job 263 (`runs/covert_p3_covert_n75_job263/`,
6 h, Δt = 0.5 s, next physically idle A100). Expected: in-cloud nᶜˡ ≤ 75 cm⁻³ and drizzle
comparable with the N75 members; the result will be appended here.
