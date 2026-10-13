# SEA STARR CTRL / N100 / N030: aerosol–precipitation response (66 h)

All three members ran the full 66 h with the same code (b54dafd): the original radiation configuration (column ending at
6.5 km, fixed solar coordinate) and the original evaporation regeneration (see the artifact below). CTRL is job 211. N100 and N030
are two segments each: jobs 214/218 to their 60-h checkpoints (the A100 nodes went down at t ≈ 62 h on 2026-10-09), then jobs
491/492 resumed from those checkpoints into new directories (`execution/sea_starr_resume.sbatch`). The segments are stitched at the
restore time (`analysis/sea_starr_segments.jl`). The CTRL rerun with the radiation fix (jobs 428 → 493) is not part of this comparison.

Restart continuity: the restored state reproduced the original 1-min records exactly (N100 ⟨LWP⟩ 85.40 g m⁻²/CF 0.994; N030 58.06/0.817).
The first progress line after the restart matched the original's line at the checkpoint except max|w| (5.69 vs 5.75 and 3.88 vs 3.87 m s⁻¹).
The resumed N100 tracked the original to 0.1 g m⁻² over the 2 h the original had continued past its checkpoint.

Products: `analysis/sea_starr_members.jl` → `runs/members_analysis/{members_timeseries,members_state_space}.png`,
`members_summary.{md,toml}` (copies `campaign-figures/seastarr_members_*`). Animations (`campaign-figures/animations/`):
`seastarr_planview_lwp_rain_ctrl_n100_n030.mp4` (5.0 MB, 265 frames, 15 min, 0–66 h, the three members side by side),
`seastarr_ctrl_xz_cloud_rain_aerosol.mp4` (0.48 MB), `seastarr_ctrl_time_height_aerosol_droplets.mp4` (0.25 MB).

## Numbers (15-min means; SEVIRI = composite of 40 trajectories, screened retrievals for LWP and N_d)

| | CTRL | N100 | N030 | SEVIRI |
| --- | --- | --- | --- | --- |
| LWP day 1 / 2 / 3 (g m⁻²) | 113 / 82 / 58 | 123 / 106 / 77 | 112 / 82 / 73 | 70 / 44 / 29 |
| cloud fraction day 1 / 2 / 3 | 0.99 / 0.92 / 0.93 | 1.00 / 0.98 / 0.98 | 0.97 / 0.83 / 0.82 | 0.86 / 0.66 / 0.41 |
| surface rain day 1 / 2 / 3 (mm d⁻¹) | 0.0004 / 0.00004 / 0.000005 | 0.005 / 0.026 / 0.009 | 0.43 / 0.68 / 0.62 | — |
| max rain (mm d⁻¹, 15-min) | 0.006 | 0.22 | 2.0 | — |
| rain fraction of liquid, day 1 / 2 / 3 | 0.001 / 0.0003 / 0.0002 | 0.007 / 0.011 / 0.007 | 0.11 / 0.17 / 0.18 | — |
| LWC-weighted in-cloud N_d day 1 / 2 / 3 (cm⁻³) | 261 / 475 / 501* | 87 / 79 / 108* | 18 / 13 / 13 | 136 / 186 / 108 |
| aerosol below zᵢ day 1 / 2 / 3 (mg⁻¹) | 275 / 590 / 673* | 89 / 102 / 124* | 21 / 19 / 22 | — |
| zᵢ rise (m h⁻¹) / zᵢ at 66 h (m) | 15.8 / 2205 | 15.8 / 2200 | 8.1 / 1690 | cloud tops 1.2 → 1.5–1.9 km |
| mean entrainment velocity (mm s⁻¹) | 8.4 | 8.3 | 5.6 | — |
| first sustained CF < 0.9 (≥ 2 h) | 40.8 h | none | 28.5 h | — |
| CF / LWP at 66 h | 0.29 / 7.8 | 0.72 / 16 | 0.76 / 36 | 0.23 / 16 |

\* after 54.5 h (CTRL) and 59.2 h (N100) biased high by the regeneration artifact (below).

## Interpretation and comparison with the intercomparison (Diamond et al., EGUsphere 2026-5350)

- **CTRL** (smoke ≈ 1000 mg⁻¹ entrained into the boundary layer, BL aerosol 100 → 670 mg⁻¹): deepening–warming transition with
  negligible drizzle, as in all five LES of the preprint. N_d of 260–500 cm⁻³ exceeds SEVIRI (100–200), the same high bias the preprint
  reports for every model ("several hundred mg⁻¹", maxima 423–716). `sea_starr_analysis.md` has the full CTRL verdict.
- **N100** (100 mg⁻¹ free troposphere): Breeze stays on the deepening–warming branch. Drizzle never exceeds 0.22 mm d⁻¹ (rain fraction
  ≈ 1 %), the deck stays nearly overcast through day 3 (CF ≥ 0.88 until the final afternoon) and the inversion deepens like CTRL. The
  final-afternoon thinning reaches CF 0.72 at 66 h, less than CTRL's 0.29, because N100 holds more liquid (LWP 77 vs 58 g m⁻² on day 3).
  This matches UW-SAM/MIMICA/DALES in the preprint; NOAA-SAM and DHARMA instead break up into open cells during the second night.
  N_d ≈ 80–110 cm⁻³ is close to the observed range ("produces Nc values closer to observed for UW-SAM and MIMICA").
- **N030** (30 mg⁻¹): drizzle persists from the first hours (0.4–0.7 mm d⁻¹ day means, peaks 2 mm d⁻¹, 11–18 % of the liquid in
  rain). N_d is 10–18 cm⁻³ (preprint: UW-SAM 14–30, DALES 10–20, NOAA-SAM/DHARMA "several mg⁻¹"). Cloud cover is lost gradually
  (CF 0.97 → 0.83 → 0.82; minimum 0.69 at 63 h) and the boundary layer deepens at half the CTRL rate (8 vs 16 m h⁻¹), ending 500 m
  shallower. This is the preprint's UW-SAM N030: "a prolonged breakup over the first day and a half … deepen[s] at a markedly slower
  pace". The boundary-layer aerosol stays near 20 mg⁻¹ because the surface source (70 cm⁻² s⁻¹) plus entrainment balance collection.
  The drizzle–depletion feedback that drives NOAA-SAM/DHARMA to open cells therefore does not run away.
- **Drizzle state space** (`members_state_space.png`, preprint Fig. 7): N030 sits beyond the adiabatic rₑ = 16 µm line with rain
  fractions 0.1–0.25; N100 straddles 12–14 µm at ≈ 1 %; CTRL lies well below 12 µm with no rain. This is the clear separation at
  rₑ ≈ 12–16 µm that NOAA-SAM, UW-SAM and DHARMA show. Breeze's warm rain (P3 with Khairoutdinov–Kogan autoconversion/accretion) sits in the
  middle of the ensemble's rain-efficiency spread: weaker than NOAA-SAM/DHARMA, stronger than MIMICA/DALES.
- **Against SEVIRI:** all three members are too cloudy on days 2–3. CTRL ends closest (CF 0.29 vs 0.23–0.29 observed in the last hours);
  N100 and N030 keep 0.72–0.76. SEVIRI screening (cloud fraction > 90 %, SZA < 70°) biases its LWP/N_d toward overcast daytime scenes,
  so the observational comparison is qualitative, as the preprint stresses.

**Numerical artifact.** CTRL (54.5 h) and N100 (59.2 h, in segment 1) each contain one event in which `EvaporationRegeneration` injected
aerosol from nothing: +7.7 % and +15 % of the domain aerosol in one minute, concentrated in a few columns. The cause is a lagged evaporation
rate divided by a collapsed cloud-water mass. N_d and boundary-layer aerosol after those times are biased high. Cloud fraction, LWP and
drizzle are not visibly affected. The fix is in PR #31 (bounded fractional rate). N030 and the resumed segments have no such event.

Not established: agreement with any individual model (no model output is distributed); the effect of the radiation fix on the members
(only the CTRL rerun uses it); results free of the regeneration artifact (needs reruns with PR #31).
