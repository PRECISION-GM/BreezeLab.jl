# Production readiness — 7 October 2026

ENA is ready for a longer **development production pilot**. Neither MIP has an
executable case adapter in this repository yet. Passing ENA tests does not establish
SEA STARR or TRACER-MIP readiness. This is a dated assessment, not a live job monitor.

| Experiment | Evidence / present status | Next acceptance milestone |
| --- | --- | --- |
| ENA public Covert configuration | 12 GPU checks passed: three microphysics choices, 4 s and 600 s on 32², 60 s on 256² on T4/H100; 192 vertical levels | Restart equivalence, one-hour full-domain pilot, precipitation/timestep and budget checks, then six-hour baselines |
| Official LASSO-ENA | Selected `samin` archive unavailable through tested order/API routes | Obtain original inputs and grid, inspect metadata, then run and evaluate against matching SAM member |
| SEA STARR | Three official DEPHY drivers, setup PDF, output workbook and trajectories downloaded and checksum-verified | Implement and test DEPHY forcing, moving inversion nudging, interactive aerosol lifecycle and absorbing radiation; execute CTRL pilot |
| TRACER-MIP | Official roadmap and aerosol notebook pinned; no executable regional case | Implement one-way nested coastal control, prescribed two-mode aerosol and process-rate exports; execute August control pilot |

## ENA: finish the transition to longer runs

Measured on H100 at 256×256×192, 35 m, Δt=0.5 s: median stepping times are
0.049 s (1M), 0.165 s (P3-N75), and 0.262 s (P3-aer2). Six-hour stepping estimates
are 0.58, 1.97 and 3.14 hours. Initial wall-time requests of 1.5, 3.5 and 5 hours
allow overhead and evolving microphysics cost; re-estimate from the one-hour pilot.
Peak H100 memory was 3.6–7.2 GiB. These are short-run measurements, not guarantees.

- [ ] Verify a stopped/restarted run against an uninterrupted run, including aerosol,
  radiation/forcing state and output continuity; record any state not restored.
- [ ] Run a one-hour simulated-time pilot on the full public domain, initially 1M and
  P3-aer2. Track water/number budgets and sedimentation-inclusive CFL as rain develops.
  Short P3 tests reached CFL≈0.55, above the bounds-preserving limiter's 5/18 guarantee;
  observed positivity alone does not resolve this. Check timestep sensitivity.
- [ ] Run the six-hour public baseline members sequentially with checkpoints; evaluate
  against staged ARM MWR LWP, VDIS rain and appropriately qualified ARSCL boundaries.

This is the public 8.96 km / six-hour setup with a reconstructed vertical grid,
not the 30.24 km / nine-hour Covert paper configuration. See [ENA](cases/ena.md).
Evidence: wpcluster campaign `breezelab-runs/20261006`, jobs 94/98/106/107,
`RESULTS.md` and per-run `summary.toml`; source `e43676c`.

## SEA STARR: first MIP implementation priority

Source: [official protocol and drivers](https://zenodo.org/records/22241697).
Input inventory: [pinned manifest](../cases/seastarr/inputs.toml).

- [ ] Implement DEPHY ingestion and reference-value tests: θl/total-water initialization,
  number-per-mass conversion, time interpolation, subsidence applied once, trajectory SST,
  and LES grid built from the specification rather than driver midlevels.
- [ ] Implement inversion-following free-tropospheric nudging and test its ramp/timescales.
- [ ] Verify aerosol activation, surface source, scavenging and evaporation/regeneration
  budgets. **Do not reuse ENA's `diagnostic_ccn` reset for interactive SEA STARR aerosol.**
- [ ] Validate absorbing-aerosol shortwave/longwave coupling and upper radiation column.
  Resolve spectral optical properties and trajectory solar geometry before production.
- [ ] Implement workbook output definitions, units, 15-minute statistics/2D and hourly 3D
  output; validate restarts and diagnostic accumulators on a short GPU CTRL pilot.
- [ ] Benchmark the allowed 192² domain, then choose it or standard 384² at 50 m.
  Run CTRL for 66 h (including 3 h spinup), then N100/N030. Compare aerosol/drizzle
  budgets and cloud evolution to reference outputs; do not substitute an ENA aerosol setup.

Specification/manuscript nudging differences and incomplete aerosol optics remain
protocol questions. Source files being available does not mean these choices are resolved.

## TRACER-MIP: regional workflow, separate from the periodic LES adapter

Source: [official roadmap](https://arm-synergy.github.io/tracer-mip/Roadmap.html),
revision [416423a](https://github.com/ARM-Synergy/tracer-mip/tree/416423ac1c5ea7673a07dacf14622af74ebd3b3e).
[Pinned references](../cases/tracer_mip/inputs.toml) do not include ERA5 forcing data.

- [ ] Stage ERA5 initialization/boundaries and surface/terrain data for the two dates.
  Build the 2 km → 500 m one-way nested workflow using NumericalEarth; verify boundary
  state transfer, coastal land/sea exchange and sea-breeze evolution.
- [ ] Port the reference two-mode aerosol height profiles with exact unit/floor tests.
  Tier 1 prescribes aerosol each step but retains prognostic droplet number; aerosol
  radiation is disabled. Do not copy SEA STARR's interactive absorbing aerosol choices.
- [ ] Export separately integrated microphysical processes, correct units and restart
  accumulators. Preserve inner-domain two-minute output during the cell-tracking window.
- [ ] Run August 7, 2022 CTRL as a 24-hour pilot, then June 17 CTRL and aerosol variants.
  Resolve 0.3 versus 1/3 low-aerosol factor and background-floor scaling before variants.
- [ ] Confirm current submission arrangements and accepted model departures with project
  coordinators; the public roadmap still gives a 2025 deadline. No contact sent by this work.

Priorities: advance ENA checkpoint/pilot work while implementing SEA STARR; prepare
TRACER inputs and aerosol/profile tests alongside it. MIP production follows case-specific
physics, restart, budget and output acceptance, not merely successful GPU execution.
