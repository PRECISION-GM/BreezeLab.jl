# ENA-LASSO — reproduce an official LASSO member in Breeze

Read [README.md](README.md) first. Status file: `status/ENA-LASSO.md`.

## Objective and current boundary

Implement a distinct `cases/ena_lasso.jl` example and run an authentic LASSO-ENA forcing experiment, comparing Breeze statistics with the corresponding SAM member and ARM observations. Retain ENA-covert separately. Matching forcing does not make Breeze's P3 scheme identical to SAM's bin microphysics.

Existing code already supports `ena_simulation(...; protocol=:lasso_ena_official, epoch=...)`. This is a partial adapter with synthetic-fixture tests, not an established full LASSO reproduction. Read `src/ena_protocols.jl`, `sam_input_files.jl`, `case_setup.jl`, `surface_fluxes.jl`, `large_scale_forcings.jl`, and `docs/cases/ena.md` before extending it.

**Data acquisition is Greg's next joint task with the coordinating assistant. Do not start another browser/API retrieval campaign merely because this plan exists.** Proceed with code inspection and supported-input validation; pick up the acquired bundle when available.

## Inputs needed

- Official `samin` archive containing `snd`, `lsf`, `sfc`, `prm`, and `grd`, plus corresponding SAM output/reference diagnostics.
- Candidate member: `20170718era5s1n0d25x100_sbmwrm-aer2-flxsst`. Another obtainable official warm-cloud/SST-flux member is acceptable if its identity is explicit; do not silently change members.
- Candidate filename: `enalasso_samin_20170718era5s1n0d25x100_sbmwrm-aer2-flxsstC1.m0.20170718.000000.tar`. On October 6 this exact file returned 404 from tested ARM API routes, despite working credentials for another dataset. Browser ordering also failed. This could be availability/staging or an incorrect archive selection; do not infer bad credentials or invent inputs.
- ARM dataset: https://doi.org/10.5439/2572661 ; browser: https://lasso-ena.arm.gov/ ; SAM source: https://code.arm.gov/lasso/lasso-ena-codes/lasso_sam_sbm . Existing downloader: `data_wrangling/fetch_arm_inputs.jl`.

## Implementation sequence

1. Freeze the obtained member: archive checksum, member ID/DOI, date/UTC epoch, original SAM revision and outputs. Inspect actual metadata before adopting any candidate start time. Stage under a case-specific external data directory.
2. Parse the original grid and namelist, including horizontal dimensions, vertical staggering, duration, timestep, aerosol choices, surface treatment, nudging and radiation cadence. SAM `day0` alone does not give a year. Require explicit missing dimensions rather than infer them from Covert. Reject unsupported configurations with an actionable message.
3. Test interpolation, units and tendency conversion against original file values at several heights/times. Distinguish dry-air mixing ratios from mass fractions; avoid double-counting subsidence or large-scale advection; verify ground-relative surface winds and pressure/height conversion.
4. Audit the existing differences that materially affect comparison: SST bulk exchange versus SAM stability/gustiness treatment; radiation atmosphere above the LES top and prescribed effective radii; nudging; aerosol replenishment timing; grid mapping; P3 versus SAM bin microphysics. Match the imposed experiment where feasible and record deliberate model-physics differences. Do not let a successful parser stand in for this audit.
5. Add a source constructor/wrapper and readable `cases/ena_lasso.jl`, retaining shared ENA utilities. Example should construct, save provenance, run, then plot. Add targeted `Pkg.test()` coverage for the real adapter using redistributable small fixtures; keep private/large ARM files out of CI and Git.
6. Run one short H100/A100 integration with the actual bundle, then its full configured duration. Check finite state, water/number bounds, forcing and radiation/surface-flux behavior, and restart/output continuity. Use a separate output directory from ENA-covert.
7. Compare common-time-window LWP, rain, cloud fraction/boundaries, thermodynamic/wind profiles and surface fluxes against the matching SAM member and available ARM observations. Use NumericalEarth metadata adapters where practical. Report reference uncertainty and unresolved physical differences instead of expecting chaotic snapshots to match.

## Done means

- Separate ENA-LASSO run script, constructor, documented input provenance and exact reproducible launch command.
- Targeted CPU tests pass; a real-bundle GPU baseline has completed, with figures and a concise comparison report.
- If the bundle is still unavailable, report that exact external blocker and completed independent code work. Do not relabel a Covert simulation as LASSO or mark the case complete.
