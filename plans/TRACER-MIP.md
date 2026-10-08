# TRACER-MIP — regional nested coastal controls

Read [README.md](README.md) first. Status file: `status/TRACER-MIP.md`.

## Objective and scope

Implement the actual regional TRACER-MIP Tier 1 control, initially August 7, 2022, followed by June 17. This is separate from the periodic TRACER–DP-SCREAM experiment. Do not present a doubly periodic ENA-like run as regional TRACER-MIP.

Official roadmap: https://arm-synergy.github.io/tracer-mip/Roadmap.html . Pinned source: https://github.com/ARM-Synergy/tracer-mip/tree/416423ac1c5ea7673a07dacf14622af74ebd3b3e . `cases/tracer_mip/inputs.toml` pins the roadmap and aerosol notebook, **not ERA5 initialization/boundaries**. `cases/tracer_mip.jl` currently contains a configuration tuple and deliberately errors when executed.

Current protocol summary to verify against the pinned files: two 24-hour cases starting 06 UTC on June 17 and August 7, 2022; 750² outer grid at 2 km and 500² inner grid at 500 m, protocol vertical levels (95 noted in earlier inspection), outer/inner timesteps 3/1.5 s. Tier 1 prescribes aerosols while prognosing droplet number; aerosol radiation is disabled. Re-read exact grid extents, levels and aerosol definitions before implementation.

## Implementation sequence

1. Inventory existing NumericalEarth capabilities before writing adapters. Inspect `examples/breeze_downscaling_era5.jl`, `examples/era5_forced_slab_land.jl`, `examples/breeze_over_slab_land.jl`, `src/DataWrangling/ERA5/`, `src/Atmospheres/`, `src/NestedModels/`, and `src/Lands/`. Pin a compatible dependency revision and identify the precise missing coupling/API pieces. Do not assume nesting or land exchange is absent, or assume existing examples already satisfy this protocol.
2. Retrieve/cache the required ERA5 atmospheric, surface and land data with NumericalEarth metadata/download APIs where possible; use existing authorized credentials without exposing them. Include temporal coverage needed for boundary interpolation. Derive terrain/land-sea mask and land initial state from the specified sources. Record provenance, units, time basis, interpolation and missing variables. Report a concrete access blocker only after checking current setup.
3. Implement the outer domain first: initialization, time-dependent open boundaries, land/sea surface exchange and coastal geometry. Run a short integration checking hydrostatic balance, boundary transients, surface energy/water exchange and plausible sea-breeze development. This is an implementation milestone, not MIP completion.
4. Add one-way 2 km to 500 m nesting using NumericalEarth. Test staggered-grid interpolation, reference-state/moisture conventions, boundary timing, terrain alignment and parent/child synchronization. Confirm the intended forcing treatment; avoid applying the same large-scale forcing both through boundaries and an unintended interior tendency.
5. Port the exact prescribed two-mode aerosol height profiles from the reference notebook. Test number/mass units, modal parameters, floors and height interpolation. Apply the required prescription cadence while allowing cloud droplet number to evolve. Do not import SEA STARR's absorbing interactive-aerosol choices. Start with CTRL; the low-aerosol factor (0.3 versus 1/3) and floor scaling remain questions for later variants, not reasons to delay CTRL.
6. Implement required separate microphysical-process diagnostics, accumulated with correct units and restart behavior. Follow prescribed outer/inner output cadence and two-minute cell-tracking output during 17–01 UTC. Determine total storage before full runs and retain required fields rather than substituting only domain averages.
7. Add source constructors and a readable `cases/tracer_mip.jl` that constructs, runs and analyzes the chosen control. CPU unit tests should exercise transformations and boundary transfer; GPU integration is manual. Once the coupled short run passes, execute the full August control, then June, then aerosol variants once definitions are resolved.

## Done means

Reproducible data staging, runnable regional/nested control, exact case metadata and dependency environment, tests and a completed 24-hour baseline with required outputs and Julia diagnostic figures. Summarize material departures from the reference land, radiation and microphysics treatments. Report observed memory/throughput and choose an H100/A100 allocation based on them; do not assume this large pair fits one GPU.

The public roadmap includes an old 2025 submission deadline. Scientific case implementation can proceed independently, but current submission arrangements and acceptance of model departures need confirmation before any formal submission. Do not contact organizers without Greg's authorization.
