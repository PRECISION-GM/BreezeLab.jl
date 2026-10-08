# TRACER-MIP August 7 2022 outer-domain control — animations

Produced by `analysis/animate_tracer_mip.jl` from the 10-minute output of the 2 km outer domain
(750 × 750 × 94, real ERA5 boundaries, P3 Tier-1 aerosol, RRTMGP, slab land with sea cells pinned to
ERA5 skin temperature; BreezeLab job 233 on an A100). This is an outer-domain control run, not a
completed MIP submission; nothing has been compared with observations.

| File | Content | Fields and units | Cadence |
| --- | --- | --- | --- |
| `tracer_mip_aug07_plan_view_domain.mp4` / `_montage.png` | whole 1500 km domain, coastline (sea mask = 0.5 contour), site (star), inner box (dashed) | left: lowest-model-level (25 m) air temperature (K) with lowest-level wind arrows every 30 cells (60 km); right: surface rain rate (mm h⁻¹) with cloud water at 2 km (0.05 and 0.5 g kg⁻¹ grey contours) | 10 min |
| `tracer_mip_aug07_plan_view_inner_box.mp4` / `_montage.png` | the 250 km inner (500 m nest) box ± 20 km | left: lowest-level temperature (K) and wind every 8 cells (16 km); right: liquid water path (g m⁻²) with surface rain-rate contours 1 / 5 / 20 mm h⁻¹ | 10 min |
| `tracer_mip_aug07_section_site.mp4` / `_montage.png` | x–z section along 29.47 N through the site, ± 145 km, 0–16 km | top: cloud + rain + ice water mixing ratio (g kg⁻¹); bottom: vertical velocity (m s⁻¹); terrain line | 10 min |

Titles give UTC and local (CDT = UTC − 5 h) time. Montages show 06, 12, 18 and 00 UTC (or the latest
available frame while the run is in progress, labelled). MP4: H.264, 1000 px wide, 8 frames/s,
compression raised until each file is under 12 MB. `animations.toml` records frame counts and sizes.
Section heights are the reference-column heights at the site (terrain here is < 50 m).
