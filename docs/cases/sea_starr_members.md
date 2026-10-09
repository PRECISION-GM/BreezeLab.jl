# SEA STARR CTRL / N100 / N030: aerosol–precipitation response (INTERIM, to 62 h)

**Status: interim.** CTRL (job 211) is complete (66 h). N100 (job 214) and N030 (job 218) reached t ≈ 62 h before both
A100 nodes went down for replacement on 2026-10-09 11:00 UTC; they are resumed from their 60-h checkpoints into new segment
directories (jobs 491, 492; `execution/sea_starr_resume.sbatch`). This page is regenerated from the stitched segments by the
post-processing job (`execution/sea_starr_postprocess.sbatch`, dependency on the resume jobs) and will be finalised then.
All three runs use the same code (b54dafd): the original radiation configuration (column ending at 6.5 km, fixed solar
coordinate) — the CTRL rerun with the radiation fix (job 428 → 493) is not part of this comparison.

Products: `analysis/sea_starr_members.jl` → `runs/members_analysis_interim/{members_timeseries,members_state_space}.png`,
`members_summary.{md,toml}`; animations: `campaign-figures/animations/seastarr_planview_lwp_rain_ctrl_n100_n030.mp4`.

## Numbers (15-min means; SEVIRI = composite of 40 trajectories, screened retrievals for LWP and N_d)

| | CTRL | N100 | N030 | SEVIRI |
| --- | --- | --- | --- | --- |
| LWP day 1 / 2 / 3 (g m⁻²) | 113 / 82 / 58 | 123 / 106 / 91 | 112 / 82 / 83 | 70 / 44 / 29 |
| cloud fraction day 1 / 2 / 3 | 0.99 / 0.92 / 0.93 | 1.00 / 0.98 / 0.99 | 0.97 / 0.83 / 0.85 | 0.86 / 0.66 / 0.41 |
| surface rain day 1 / 2 / 3 (mm d⁻¹) | 0.0004 / 0.00004 / 0.000005 | 0.005 / 0.026 / 0.011 | 0.43 / 0.68 / 0.74 | — |
| max rain (mm d⁻¹, 15-min) | 0.006 | 0.22 | 2.0 | — |
| rain fraction of liquid, day 2 | 0.0003 | 0.011 | 0.17 | — |
| LWC-weighted in-cloud N_d day 1 / 2 / 3 (cm⁻³) | 261 / 475 / 501 | 87 / 79 / 104 | 18 / 13 / 13 | 136 / 186 / 108 |
| aerosol below zᵢ day 1 / 2 / 3 (mg⁻¹) | 275 / 590 / 673 | 89 / 102 / 116 | 21 / 19 / 22 | — |
| zᵢ rise (m h⁻¹) / final zᵢ (m) | 15.8 / 2205 | 17.1 / 2220 (62 h) | 9.2 / 1730 (62 h) | cloud tops 1.2 → 1.5–1.9 km |
| mean entrainment velocity (mm s⁻¹) | 8.4 | 8.6 | 5.9 | — |
| first sustained CF < 0.9 (≥ 2 h) | 40.8 h | none (min CF 0.88 at 44 h) | 28.5 h | — |
| final CF / LWP | 0.29 / 7.8 (66 h) | 1.00 / 43 (62 h) | 0.72 / 40 (62 h) | 0.23 / 16 (66 h) |

## Interpretation and comparison with the intercomparison (Diamond et al., EGUsphere 2026-5350)

- **CTRL** (smoke ≈ 1000 mg⁻¹ entrained into the boundary layer, BL aerosol 100 → 670 mg⁻¹): deepening–warming transition with
  negligible drizzle, as in all five LES of the preprint; N_d 260–500 cm⁻³ is above SEVIRI (100–200), the same high bias the preprint
  reports for every model ("several hundred mg⁻¹", maxima 423–716). See `sea_starr_analysis.md` for the full CTRL verdict.
- **N100** (100 mg⁻¹ free troposphere): Breeze stays on the deepening–warming branch — drizzle never exceeds 0.22 mm d⁻¹ (rain
  fraction ≈ 1 %), the cloud deck persists (CF ≥ 0.88) and the inversion deepens like CTRL. This is the UW-SAM / MIMICA / DALES
  behaviour in the preprint; NOAA-SAM and DHARMA instead break up into open cells during the second night. N_d ≈ 80–100 cm⁻³ is close
  to the observed range ("produces Nc values closer to observed for UW-SAM and MIMICA").
- **N030** (30 mg⁻¹): persistent drizzle from the first hours (0.4–0.7 mm d⁻¹ day means, peaks 2 mm d⁻¹, 10–20 % of the liquid in
  rain), N_d 10–18 cm⁻³ (preprint: UW-SAM 14–30, DALES 10–20, NOAA-SAM/DHARMA "several mg⁻¹"), a gradual loss of cloud cover
  (CF 0.97 → 0.83–0.85 → 0.72) and markedly slower deepening (9 m h⁻¹ vs 16–17) — the preprint's UW-SAM N030 ("a prolonged breakup
  over the first day and a half … deepen[s] at a markedly slower pace"). The boundary-layer aerosol stays near 20 mg⁻¹ (surface source
  70 cm⁻² s⁻¹ plus entrainment balance collection), so the positive drizzle–depletion feedback of NOAA-SAM/DHARMA does not run away.
- **Drizzle state space** (`members_state_space.png`, preprint Fig. 7): N030 sits beyond the adiabatic rₑ = 16 µm line with rain fractions
  0.1–0.25; N100 straddles 12–14 µm with ≈ 1 %; CTRL lies well below 12 µm with no rain — the clear separation at rₑ ≈ 12–16 µm that
  NOAA-SAM, UW-SAM and DHARMA show. Breeze's warm rain (P3 with Khairoutdinov–Kogan autoconversion/accretion) is therefore in the
  middle of the ensemble's rain-efficiency spread: weaker than NOAA-SAM/DHARMA, stronger than MIMICA/DALES.
- **Against SEVIRI** all three members are too cloudy on days 2–3 (CTRL is closest at the end); the composite's earlier breakup is
  better matched by none of the Breeze members. As the preprint stresses, SEVIRI screening (cloud fraction > 90 %, SZA < 70°) biases
  its LWP/N_d to overcast, daytime scenes; the observational comparison is qualitative.

What is not established: agreement with any individual model (no model output is distributed); the N100/N030 final 4 h; the effect
of the radiation fix on the members (only the CTRL rerun uses it).
