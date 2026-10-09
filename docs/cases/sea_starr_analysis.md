# SEA STARR CTRL (job 211): full-run analysis, comparison with references, and verdict on the breakup

Run: `/shared/home/greg/breezelab-work/sea-starr/runs/ctrl_full_66h_job211/` (commit b54dafd, H100, 192² × 288,
50 m / 10 m, Δt ≈ 0.56 s under the 1 s cap, 431 310 steps, 22.1 h wall, `COMPLETE`). Analysis:
`analysis/sea_starr_ctrl_analysis.jl` (CPU partition, job 378), figures and `summary.{md,toml}` under
`runs/ctrl_full_66h_job211/analysis/`; animations under `runs/ctrl_full_66h_job211/animations/` and
`/shared/home/greg/breezelab-work/campaign-figures/animations/`. Local solar time (LST) uses the hourly mean
longitude of the 40 GEOS-5 trajectories (−1.2 °E at t = 0 → −14.9 °E at 66 h; LST = UTC + lon/15 h).

Figures: `ctrl_timeseries_vs_seviri.png` (15-min statistics vs SEVIRI), `ctrl_time_height.png` (θₗ, qₜ, qᶜˡ,
nᵃ, nᶜˡ, radiative heating, nudging mask, cloud fraction), `ctrl_profiles_vs_targets.png` (model vs driver
nudging targets).

## 1. Protocol diagnostics over the 66 h

| quantity | day 1 (0–24 h) | day 2 (24–48 h) | day 3 (48–66 h) | start → end |
| --- | --- | --- | --- | --- |
| LWP cloud + rain (g m⁻²), model | 113 | 82 | 58 | 101 → 7.8 (max 171 at 11 h) |
| LWP, SEVIRI composite (all / screened) | 59 / 70 | 54 / 44 | 26 / 29 | 55 → 16 |
| cloud fraction (LWP > 5 g m⁻²), model | 0.99 | 0.92 | 0.93 | 1.00 → 0.29 (min 0.29 at 66 h; 0.74 at 42 h) |
| cloud fraction, SEVIRI | 0.86 | 0.66 | 0.41 | 0.92 → 0.23 (last 6 h mean 0.29) |
| surface rain (mm d⁻¹) | 0.0004 | 0.00004 | 0.000005 | max 0.006 (first hour) |
| ⟨zᵢ⟩ (m) | 1160 → 1546 | → 1950 | → 2207 | rise 17.4 m h⁻¹ |
| cloud top / base (m), model | 1035 / 735 → | 1565 / 1000 → | 2195 / 2095 | |
| cloud top, SEVIRI screened (km) | 1.2–1.3 | 1.4–1.7 | 1.5–1.9 | |
| N_d at cloud top (cm⁻³), model | 100–200 | 231 mean | 385 at 66 h | column max 300–600 |
| N_d, SEVIRI screened (cm⁻³) | 130–140 | 186 mean | 100–150 | |
| nᵃ below zᵢ (mg⁻¹) | 118 → 427 | → 652 | → 731 | smoke entrained from the FT |
| nᵃ at max zᵢ + 500 m (mg⁻¹) | 995 → 1300 | → 400 | → 364 | follows the driver's na_nud |

Aerosol budget: the column total (interstitial + droplets + rain) goes 2.63 → 2.58 ×10¹² m⁻² while the surface
source adds 1.66 ×10¹¹ m⁻²; the balance is the free-tropospheric nudging/subsidence of the smoke layer
(whose target first rises to ≈1480 mg⁻¹ and then falls), not scavenging — surface rain is < 0.01 mm d⁻¹ throughout.
Droplets hold 1–6 % of the column total.

Radiation (15-min means; signs positive-up in the files): LW↓ surface 394–409 W m⁻², LW↑ surface 436, LW↑ at the
LES top 329–337, LW↓ at the LES top **0 (configuration departure, see §3)**; SW↓ at the LES top peaks at 1233
(TOA value, no attenuation above 6.5 km — also a departure), SW↓ surface up to 887 W m⁻². Net LW flux jump across
the cloud top ≈ 100 W m⁻² at night, 25–70 W m⁻² by day; the strongest radiative cooling rate is −265 K d⁻¹ in the
10 m cloud-top cell.

Entrainment: w_e = dzᵢ/dt − w_s(zᵢ) with the driver's subsidence at zᵢ (−1.7 → −4.3 mm s⁻¹, mean −3.5 mm s⁻¹):
mean 8.4 mm s⁻¹ (8.3 / 8.5 / 8.2 by day), nocturnal maxima ≈ 11–12 mm s⁻¹, daytime minima 2–4 mm s⁻¹ (max 16
during the first hour's adjustment). For comparison, DYCOMS-II RF01 nocturnal Sc entrains at ≈ 4 mm s⁻¹ and
Lagrangian Sc-to-Cu transitions (ASTEX, CGILS S12 transitions) at ≈ 5–15 mm s⁻¹; Diamond et al. (2022) report
≈ 1 km of MBL deepening over 3 days for this trajectory (ours: 1.05 km).

## 2. References: what exists and what does not

- **SEA STARR intercomparison preprint** (Diamond et al., EGUsphere 2026-5350, discussion started 6 Oct 2026; five
  LES: NOAA-SAM, UW-SAM, MIMICA, DHARMA, DALES — this is the project's own result for exactly this driver): in CTRL
  "all models show qualitatively similar behavior: mainly overcast with limited breakup during the first day (6–18 hours),
  full recovery over the second night (18–30 hours), greater breakup during the second day (30–42 hours), at least
  partial recovery at night (42–54 hours), and extensive breakup over the final day near Ascension Island. This pattern
  fits the SEVIRI observations qualitatively well, albeit with a generally stronger nocturnal recovery before Day 3."
  "All LES models except DALES deepen more rapidly than the clouds tracked by SEVIRI"; the models "trend on the high
  end" of the LASIC radiosonde inversion heights, DHARMA "deepens well beyond". N_d of "several hundred mg⁻¹
  (maximum values of 716 … 423 …)" versus observed 100–200 mg⁻¹. All models are "deepening-warming" transitions in
  CTRL; drizzle-depletion transitions appear only in N100/N030 and only in some models. Their radiation uses "the
  simulated properties below 6.5 km and a constant atmospheric profile above 6.5 km taken from the GEOS-5
  composite", the solar zenith angle follows the trajectory, and a >10 W m⁻² inter-code difference in LW↓ at 6.5 km
  changed MBL deepening by "several hundred meter[s]". Lowering the smoke SSA 0.85 → 0.80 "amplifies the diurnal
  cycle … but does not fundamentally change the pacing of the SCT" and "prevent[s] ~100–200 m of MBL growth".
  The preprint shows curves (Figs. 3–5), not tabulated numbers; no model output is distributed with the Zenodo record.
- **SEVIRI composite along the 40 trajectories** (`SEA_STARR_Raw_Trajectories.nc`, the observational constraint the
  protocol supplies): cloud fraction, LWP, N_d, cloud-top height (table in `summary.md`). Day means above.
- **Diamond et al. (2022, ACP 22, 12113)**: the single-trajectory SAM predecessor of this case (same dates, 50 m,
  12 km, 324 levels to 7 km): "almost completely overcast throughout" for 3 days, ≈ 1 km of deepening, day-3 LWP
  ≈ 175–200 g m⁻², N_c strongly enhanced by smoke, precipitation suppressed.
- Setup PDF and Zenodo record: no expected-evolution description, no reference output (the record description is empty).
- Not available: the preprint's model output, ORACLES/CLARIFY/LASIC in-situ data (not staged here), the Blossey (2013) ramp.

## 3. Is the breakup physical or a flaw?

**Verdict: the breakup is physical.** The CTRL follows the deepening–warming transition that all five intercomparison
LES simulate for this driver and that the SEVIRI composite shows: daytime thinning on day 1 (LWP 168 → 26 g m⁻² between
09 and 15 UTC, ≈ 08:40–14:40 LST, CF stays ≥ 0.98), full nocturnal recovery (149 g m⁻² at 00 UTC), a deeper day-2
collapse (13 g m⁻², CF 0.74 at 42 h ≈ 14:20 LST), partial recovery (CF 1.0, LWP 77 at 54 h) and the final afternoon
breakup to CF 0.29 and 8 g m⁻² at 66 h (14:00 LST near Ascension). The SEVIRI composite breaks up *earlier and more*
(CF 0.66 on day 2, 0.41 on day 3, 0.23 at 63–66 h) with roughly half the LWP; the model's end state (CF 0.29) matches
the observed last-6-h mean (0.29). Relative to the observations and to the preprint's ensemble, Breeze is therefore on
the cloudy/slow side of the transition, not an outlier that breaks up too early. Candidate causes, checked:

1. **Missing smoke absorption** (Breeze RRTMGP has no aerosol optics; SSA effectively 1). Unfixable here. Its sign is
   *more* deepening and a *weaker* diurnal amplitude (above-cloud heating stabilizes the inversion; in-BL heating thins
   cloud). Quantified by the preprint's SSA sensitivity: ≈ 1 K d⁻¹ peak absorption difference, 100–200 m of deepening.
   Note the FT thermodynamic state is nudged to the GEOS composite (which contains the real smoke heating), so only the
   direct in-LES absorption is missing.
2. **Inversion rise 17 m h⁻¹ / entrainment 8.4 mm s⁻¹.** Checked that the subsidence applied is the driver's `wa`
   (one `LargeScaleVerticalAdvection` per prognostic, max |Δw| = 0 vs the driver, −4.1 mm s⁻¹ at 2 km, warming tendency
   +2.6×10⁻⁵ K s⁻¹ on the stable FT — correct sign), that the nudging acts on u, v, θ, qᵛ, nᵃ only above max zᵢ + 100 m
   (mask base − max zᵢ = 79 m on average, the first cell above) and holds the FT to the targets (|Δθ| ≤ 0.09 K,
   |Δqᵛ| ≤ 0.04 g kg⁻¹, |Δnᵃ| ≤ 31 mg⁻¹ in the 300 m–4 km layer above the inversion). The deepening (1.05 km in 66 h)
   equals Diamond et al. (2022)'s ≈ 1 km and lies with the SAMs/DHARMA on the high side of the preprint's envelope;
   SEVIRI cloud tops reach 1.5–1.9 km by day 3 versus our 2.2 km. The entrainment rate is high but within the range of
   Lagrangian transitions and is consistent with the cloud-top LW cooling of ≈ 100 W m⁻² at night.
3. **Radiation column ending at 6.5 km** (LW↓ = 0 at the top, TOA SW at the top, fixed solar coordinate). Confirmed as
   a departure from what every intercomparison model does. Quantified with the new `extended_column_radiation`
   (driver's time-mean atmosphere, 48 layers to 78 km): LW↓ at 6.5 km becomes 113 W m⁻² (was 0), but only +2.5 W m⁻²
   more reaches 1.4 km and +4.8 W m⁻² reaches 2.5 km (the FT water vapour absorbs the window); SW↓ at the LES top falls
   7 % (1124 vs 1209 W m⁻² at cos θ ≈ 0.88) and at the surface 5 %. The trajectory solar position changes cos θ by
   < 0.01. These are a few per cent of the cloud-top energy budget — a real bias toward slightly stronger cloud-top
   cooling and SW absorption, i.e. toward faster deepening, but not the cause of the breakup.
4. **Initial inversion smearing**: the one-cell initial inversion spreads over ≈ 150 m in the first 15 min (zᵢ 1160 →
   1040 m, back to 1160 m by 2 h) — a spin-up transient inside the protocol's 3 h spin-up.
5. Not a cause: drizzle (≤ 0.006 mm d⁻¹; the model stays on the entrainment branch, consistent with N_d of 150–400 cm⁻³).

**Fix applied and relaunch.** The radiation configuration is corrected on the branch (`radiation = :rrtmgp_extended`,
`solar = :trajectory` are now the defaults of `sea_starr`): the RRTMGP column carries the driver's upper atmosphere above
the LES top and the zenith angle follows the composite trajectory. CPU tests cover the upper layers, the top-of-domain
fluxes and the zenith angle. The corrected CTRL is queued as `ctrl_r2_66h` (see the status file for the job ID). The
remaining departures are the absent smoke absorption (needs aerosol optics in Breeze), the linear nudging ramp, the
proportional evaporation regeneration and the upper sponge; the N_d excess (≈ 2× SEVIRI on day 3) is shared with the
intercomparison models.

## 4. Numbers behind the verdict (from `summary.toml`)

`w_entrainment_mean_mm_s = 8.36`, `w_subsidence_at_zi_mean_mm_s = −3.54`, `zi_mean_rise_m_per_h = 17.4`,
`subsidence_max_abs_diff_vs_driver_m_s = 0`, `subsidence_sign_ok = true`, `nudging_base_minus_zimax_mean_m = 79`,
`cloud_top_lw_cooling_night_W_m2 = 99`, `lwd_top_mean_W_m2 = 0`, `swd_top_max_W_m2 = 1233`,
`day3_cf_mean = 0.93` vs `day3_seviri_cf_mean = 0.41`, `cf_final = 0.29` vs `seviri_cf_final_6h_mean = 0.29`,
`day3_lwp_mean_g_m2 = 57.5` vs `day3_seviri_lwp_good_mean_g_m2 = 29.3`, `nd_top_day2_cm3 = 231` vs
`seviri_nd_good_day2_cm3 = 186`.
