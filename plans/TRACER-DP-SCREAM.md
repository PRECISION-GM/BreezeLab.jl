# TRACER–DP-SCREAM — matched periodic comparison

Read [README.md](README.md) first. Status file: `status/TRACER-DP-SCREAM.md`.

## Objective and evidence

Implement `cases/tracer_dp_scream.jl` as a separate periodic Breeze experiment using the initialization, forcing, surface treatment and comparison windows of the published DP-SCREAM TRACER study. This provides a direct comparison path without making regional TRACER-MIP nesting a prerequisite. It is neither TRACER-MIP nor ENA-SCREAM.

Sources:

- Oware et al., *Evaluating E3SM Global Storm-Resolving Model Simulations of Deep Convection: Insights From DP-SCREAM During TRACER*: https://doi.org/10.1029/2025JD044113 . Read the paper and supplements, not just the abstract/archive description.
- Output archive inspected: https://zenodo.org/records/15133179 . This is version 1 and advertises a newer version; resolve the current record, compare contents and pin the chosen version before using it.
- Case catalog: https://github.com/E3SM-Project/scmlib/wiki/E3SM-Intensive-Observation-Period-Case-Library . It identifies `DPxx_SCREAM_scripts/run_dpxx_scream_TRACER.csh` and a companion SCM script. Inspect the actual current repository and pin a revision; a catalog entry alone is not verification of input availability.

The inspected archive contains model **outputs**, including 3 km and 0.5 km runs for August 5–15, 2022; a longer 3 km run and SCM comparator; SHOC mixing-length/buoyancy-flux sensitivities; and August 6 initialization variants labeled in local time. These are not yet verified as a complete forcing bundle or full 3D state archive.

## Implementation sequence

1. Obtain and inspect the publication, run script, namelist and linked IOP inputs. Identify exact DP-SCREAM version/commit (do not conflate Fortran v0 and GPU v1), domain dimensions, vertical levels, initialization/forcing time basis, spinup, surface fluxes, nudging, radiation, aerosols and transport choices. Follow the script's actual input locations. Convert local-time labels only after checking NetCDF time units/calendar and the paper's conventions.
2. Download the chosen public reference output files and inspect variables, dimensions, averaging intervals and units. Determine which comparisons they actually support. Do not infer full 3D fields from small archive file sizes or promise diagnostics absent from the archive.
3. Freeze a baseline from the published August 5–15 comparison. Preserve the actual initialization and continuous forcing history; selecting a one-day demonstration may require a full spinup or restart and must be labeled separately. If reference forcing cannot be retrieved, record exact filename/location and what was attempted. Do not reconstruct it from model outputs or silently replace it with ERA5.
4. Implement an IOP forcing reader and periodic model constructor in BreezeLab, reusing existing forcing/thermodynamic utilities. Check vertical pressure/height coordinates, moisture basis, imposed tendencies, surface fluxes, wind frame and nudging; apply large-scale transport exactly once. Write reference-value tests against small redistributable input excerpts.
5. Define two distinct scientific comparisons: a Breeze LES driven by the same experiment, and, if desired, a like-resolution 500 m/3 km run with an appropriate closure. A LES closure at 3 km is not automatically comparable to SHOC. Record intentional physics and resolution differences and use spatial/time filtering when comparing fine LES to coarse DP-SCREAM output. Start with one defensible baseline rather than immediately reproducing all sensitivities.
6. Build a readable `cases/tracer_dp_scream.jl` around a source constructor; include provenance and Julia analysis. Run one short GPU integration, check water/energy budgets and finite state, then the selected full baseline on available H100/A100 resources. Save checkpoints so long periods can resume without losing forcing/diagnostic state.
7. Compare available cloud fraction/time-height structure, LWP/IWP, precipitation, thermodynamic profiles, fog occurrence and convective timing against DP-SCREAM and ARM observations where available. Use identical windows and explicit definitions. The objective is understanding differences; agreement with DP-SCREAM alone is not observational validation.

## Done means

A distinct executable periodic TRACER case with authentic pinned inputs, a verified reference-output reader, targeted CPU tests, completed baseline output and comparison plots, plus exact launch command/commit/job paths. Document remaining discrepancies and data gaps. This agent does not own regional nesting, TRACER-MIP submission, or the ENA-SCREAM reference experiment.
