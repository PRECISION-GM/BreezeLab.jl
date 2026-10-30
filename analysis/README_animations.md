# TRACER-MIP August 7 2022 outer-domain control — animations

Produced by `analysis/animate_tracer_mip.jl` from the fast output (every `TRACER_MIP_FAST_MINUTES`: 2.5 min in job 501,
1 min by default from PR #28 on) of the 2 km outer domain (750 × 750 × 94, real ERA5 boundaries, P3 Tier-1 aerosol,
RRTMGP, slab land with sea cells pinned to ERA5 skin temperature). `PROVENANCE.txt` names the run and the commit,
and `animations.toml` records the output cadence, the frames and the playback of each movie. This is an outer-domain
control run, not a completed MIP submission. The comparison with ERA5 and the AMF1 site is in the `control_comparison`
and `site_comparison` files next to these movies.

| File | Content | Fields and units | Cadence |
| --- | --- | --- | --- |
| `tracer_mip_aug07_plan_view_domain.mp4` / `_montage.png` | whole 1500 km domain, coastline (sea mask = 0.5 contour), site (star), inner box (dashed) | left: lowest-model-level (25 m) air temperature (K) with lowest-level wind arrows every 30 cells (60 km); right: surface rain rate (mm h⁻¹) with cloud water at 2 km (0.05 and 0.5 g kg⁻¹ grey contours) | every output time |
| `tracer_mip_aug07_plan_view_inner_box.mp4` / `_montage.png` | the 250 km inner (500 m nest) box ± 20 km | left: lowest-level temperature (K) and wind every 8 cells (16 km); right: liquid water path (g m⁻²) with surface rain-rate contours 1 / 5 / 20 mm h⁻¹ | every output time |
| `tracer_mip_aug07_plan_view_inner_box_17-01UTC.mp4` | the inner box during the sea-breeze/convection window, 17–01 UTC (12–20 CDT) | as the inner-box movie | every output time, slower playback |
| `tracer_mip_aug07_section_site.mp4` / `_montage.png` | x–z section along 29.47 N through the site, ± 145 km, 0–16 km | top: cloud + rain + ice water mixing ratio (g kg⁻¹); bottom: vertical velocity (m s⁻¹); terrain line | every output time |

Titles give UTC and local (CDT = UTC − 5 h) time. Montages show 06, 12, 18 and 00 UTC (or the latest
available frame while the run is in progress, labelled). MP4: H.264, 1000 px wide.
Playback is set in simulated minutes per second of video: ≤ 25 min/s for the plan views, ≤ 12.5 min/s for the section and
the 17–01 UTC window. The frame rate follows from the output spacing (10 fps for 2.5-min output, 12 fps for 1-min output),
and compression is raised until each file is under 12 MB. `animations.toml` records frame counts and sizes.
Section heights are the reference-column heights at the site (terrain here is < 50 m).
