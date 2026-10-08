# SEA STARR — implement and run the CTRL LES

Read [README.md](README.md) first. Status file: `status/SEA-STARR.md`.

## Objective and available inputs

Turn `cases/sea_starr.jl` from a specification into a readable, executable Breeze example. First deliver CTRL, then N100/N030 using their own drivers. No ARM account is required for these public inputs.

Official protocol/drivers: https://zenodo.org/records/22241697 . The pinned checksums and URLs are in `cases/seastarr/inputs.toml`. Three DEPHY driver NetCDFs, setup PDF, standard-variable workbook and raw trajectories have already been downloaded and verified by this project. Locate the staged files in the campaign data directories; if absent from your checkout, use `data_wrangling/fetch_manifest.jl` with the pinned manifest rather than redownloading arbitrary latest files. Verify against the manifest.

Nominal case in the current specification: 2017-08-15 21 UTC to 2017-08-18 15 UTC, 66 hours including 3 hours spinup; 50 m horizontal spacing, standard 384² or protocol-permitted 192²; 10 m lower vertical spacing to 2500 m, stretching above to 6500 m. Reconcile every setting with the official setup and driver before freezing CTRL.

## Implementation sequence

1. Read the setup PDF and NetCDF metadata together. Make a compact table of variable names, dimensions, units, initialization/forcing role, interpolation and target Breeze quantity. Identify whether each thermodynamic tendency already includes vertical transport or radiation. Record exact ambiguous nudging/optical-property questions; do not stall unrelated code.
2. Implement a DEPHY reader and `sea_starr` constructor in `src/`. Initialize liquid-water potential temperature/total water consistently with Breeze thermodynamics, pressure/reference density and moisture basis. Convert number per mass/volume explicitly. Build the LES grid from the protocol, not from driver midlevel coordinates.
3. Wire time-varying trajectory SST, surface exchange, prescribed large-scale tendencies/subsidence and wind nudging. Implement inversion-following free-tropospheric thermodynamic nudging, including inversion diagnosis, ramp and timescales from the authoritative specification. Write reference-value tests for interpolation, units, nudging-mask behavior and applying subsidence exactly once.
4. Implement the interactive aerosol budget: activation, specified surface source, scavenging and evaporation/regeneration. Check number and water transfers process by process. **Do not reuse ENA's diagnostic CCN reset**, which would erase the intended aerosol evolution. Inspect current Breeze/CloudMicrophysics capabilities first and add only missing mechanisms.
5. Connect absorbing-aerosol radiation, trajectory-dependent solar geometry and the required atmosphere above the LES domain. Confirm spectral properties with the supplied protocol/data. If optics or conflicting nudging definitions remain unspecified, explicitly isolate that blocker; any provisional run must be labeled exploratory, not a protocol-complete CTRL.
6. Implement standard-variable definitions/units from the workbook, including required profiles, budgets, 15-minute statistics/2D and hourly 3D outputs as specified. Ensure restarts preserve evolving aerosol, forcing position and accumulated diagnostics. Add small CPU tests to `Pkg.test()` and a manually generated Documenter/Literate example.
7. Run a short GPU CTRL on the allowed 192² grid, checking cloud development, budgets, radiation/surface fluxes, finite state and timestep behavior. Then run the full 66-hour CTRL, with checkpoints and a declared output/storage budget. Move to N100/N030 only after the baseline setup is checked; use their actual driver files.

## First deliverable and completion

First concrete deliverable: a driver-backed constructor, tests of imposed forcing, and a short GPU integration with explicit physics coverage. This should precede broad refactoring or a general case-framework redesign.

Complete baseline deliverable: `cases/sea_starr.jl` runs top to bottom; all required CTRL mechanisms are implemented or any agreed departure is explicit; full integration and required diagnostics are present; Julia figures show LWP, precipitation, aerosol/droplet evolution and budget closure. Include commit, command, job IDs, input hashes and output paths. Do not claim agreement with another model unless its reference outputs have actually been obtained and compared.
