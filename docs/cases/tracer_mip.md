# TRACER-MIP: regional nested coastal Tier 1 control

This page records the verified protocol of the TRACER model intercomparison (Houston, 17 June and 7 August 2022), the
NumericalEarth capabilities BreezeLab builds on, what is implemented, and every material departure from the reference
treatments. It is a working record, not a claim of a completed MIP submission.

Pinned sources (`cases/tracer_mip/inputs.toml`, checksum-verified): `Roadmap.md`, the March-2025 roadmap PDF and the
aerosol-profile notebook of [ARM-Synergy/tracer-mip @ 416423ac](https://github.com/ARM-Synergy/tracer-mip/tree/416423ac1c5ea7673a07dacf14622af74ebd3b3e).
The numbers below were re-read from those files; `cases/tracer_mip/protocol.toml` holds the transcription that the code
reads (`tracer_mip_protocol()` validates it).

## Protocol (Roadmap Tables 1-3)

| Item | Protocol value |
| --- | --- |
| Cases | 24 h from 06 UTC 2022-06-17 and 06 UTC 2022-08-07 |
| Initial and boundary data | 3-hourly ERA5 0.25° pressure-level reanalysis |
| Nests | 2, one-way only, common center 29.4719 N, 95.0792 W; polar stereographic "or similar" |
| Grid-1 | 750 × 750 at 2000 m, Δt = 3 s, full output every 60 min |
| Grid-2 | 500 × 500 at 500 m (the innermost ACPC-MIP nest), Δt = 1.5 s, full output every 10 min and every 2 min from 17 UTC to 01 UTC for cell tracking |
| Vertical | 95 ACPC scalar levels (m AGL, first level -24 m), top ≈ 22 km / 50 hPa |
| Physics | Coriolis on, no convection scheme, interactive land surface (urban physics if available), two-moment or bin microphysics with prognostic droplet number, best PBL/diffusion every step, LW+SW radiation every 60 s, radiatively inactive aerosols |
| Aerosol (Tier 1) | fixed at each step; two lognormal modes: Aug 7: 1425 cm⁻³ (Dm 49 nm, σ 1.8) + 263 cm⁻³ (175 nm, 1.4); Jun 17: 3443 cm⁻³ (30 nm, 1.5) + 531 cm⁻³ (136 nm, 1.5); converted with ρ = 1.159 kg m⁻³ to mg⁻¹; lidar-derived shape, total floor 50 mg⁻¹, floor only above 6 km; κ = 0.26 |
| Variants | ×3 and "3× lower (coefficient 0.3)" surface numbers, same shape; the roadmap also writes ⅓ — unresolved for LOW |
| Required outputs | 3-D state and water species with numbers, 2-D surface/radiation fields, and separate process rates (condensation, evaporation, deposition, sublimation, melting, freezing, droplet and ice nucleation, riming of cloud and rain, autoconversion + accretion) as kg/kg per output interval |

## NumericalEarth capability inventory (verified 2026-10-07)

The ERA5, `Atmospheres`, `NestedModels` and `Lands` modules and the examples `breeze_downscaling_era5.jl`,
`era5_forced_slab_land.jl` and `breeze_over_slab_land.jl` live in
[NumericalEarth.jl](https://github.com/NumericalEarth/NumericalEarth.jl), not in Breeze. Pinned revision
`d07eb240a5679b7a6913e4d979e8378e76786ad7` (v0.8.1, 2026-10-06): its compat (`Breeze = "0.11"`,
`Oceananigans = "0.111-0.113"`) accepts the project pins Breeze 5264d3c (0.11.3) and Oceananigans 847e125 (0.113.3).

Available: an ERA5 pressure-level `PrescribedAtmosphere` parent on the native 0.25° grid with time-varying geopotential
heights; `nested_atmosphere_model(child_grid, ERA5HourlyPressureLevels(); dates, dir, terrain = ETOPO2022(), …)` that builds a
Breeze compressible split-explicit child on a terrain-following `LatitudeLongitudeGrid` with open lateral boundary conditions,
Davies relaxation, a ρw lid/wall sponge, ERA5 initialization and adiabatic balancing; `SlabLand` (prognostic skin temperature,
bucket or variably-saturated hydrology) coupled through `AtmosphereLandModel` with Monin-Obukhov fluxes and RRTMGP surface
radiation; `PrescribedOcean` for SST-driven atmosphere-ocean fluxes.

Missing or mismatched for this protocol:

1. Nesting a Breeze child in a *Breeze* parent is not implemented (the parent must be a `PrescribedAtmosphere`); one-way
   2 km → 500 m nesting is therefore done offline: the outer run writes its prognostic state over the inner region and a
   `PrescribedAtmosphere` is built from that output.
2. Land and ocean fluxes under one Breeze atmosphere are "at most one surface type per cell"; the coastal domain needs a mask
   (BreezeLab pins sea cells of the slab to ERA5 SST with a saturated bucket and sea roughness).
3. The extension imports `Breeze.AtmosphereModels.surface_precipitation_flux`, renamed `bottom_precipitation_flux` before the
   pin; Julia warns and falls back to a zero rain flux into the land bucket.
4. P3 `AerosolActivation` modes are vertically uniform; Tier-1 height-dependent prescribed aerosol needs a scaled activation.
5. No accumulated process-rate diagnostics; BreezeLab adds them.
6. NumericalEarth's weak-dependency compat forces RRTMGP 1.0.1 → 0.22.3 when added to the BreezeLab environment, so it is
   pinned only in the case sub-environment `cases/tracer_mip/Project.toml`.

## BreezeLab implementation choices and departures

- **Horizontal grid**: Oceananigans `LatitudeLongitudeGrid` instead of polar stereographic. Extents are chosen so that cells are
  2 km × 2 km (outer) and 500 m × 500 m (inner) at the center latitude: outer 22.727-36.217 N, 102.827-87.332 W; inner
  28.348-30.596 N, 96.370-93.788 W (inside the roadmap's corner envelope). Zonal spacing varies from 2.12 km (south) to
  1.85 km (north) across the outer domain.
- **Vertical grid**: the 95 scalar levels are RAMS `zt` heights; interfaces are their midpoints, the bottom face is 0 m and the
  top face extrapolates to 22181 m: 94 cells, Δz = 50 m at the surface, 300 m above ≈ 5.5 km (`acpc_vertical_faces`).
- **Aerosol**: exact port of the notebook (`tracer_mip_aerosol_profile`), tested against its printed tables. The bulk κ = 0.26 is
  imposed by setting the P3 solute parameter β_act = κ (van 't Hoff factor derived from the chosen density/molar mass).
- **Microphysics**: Breeze P3 (two-moment ice with predicted properties, prognostic droplet number). Not a bin scheme.
- **Land**: `SlabLand` (skin temperature + bucket) is a slab, not a multi-layer soil/vegetation/urban model; initialized from
  ERA5(-Land) skin temperature and soil moisture. No urban physics.
- **Radiation**: RRTMGP all-sky every 60 s; aerosols are radiatively inactive (protocol).
- **Boundaries**: hourly ERA5 (the protocol minimum is 3-hourly), interpolated linearly in time; Davies relaxation over the
  outermost cells; no interior large-scale tendency is applied (forcing enters only through the boundaries and the surface).
- **Turbulence/PBL**: `TKEBasedTurbulenceClosure` (vertical eddy diffusivity with prognostic TKE, vertically implicit) is
  the control's closure. The first 24-h control (job 233) ran without a closure and shows a surface-layer defect
  (`runs/outer_era5_control_aug07_job233/diagnostics/surface_wind_check.{toml,png}`): the lowest-level (25 m) wind,
  initialised at ERA5's 1000 hPa value (domain mean 7.2 m/s vs ERA5 10 m 4.9 m/s — the ERA5 pressure levels carry no
  surface layer), falls within one hour to 2.6–3.0 m/s, 0.58–0.83 of ERA5's 10 m wind for the rest of the run, while the
  second level (76 m) keeps the free-stream value (site profiles: 3.1 vs 8.0 m/s at 3 h, 1.9 vs 4.4 m/s at 9 h, with
  ERA5 5.3 and 3.8 m/s at ≈160 m). The surface stress is confined to the 50 m lowest layer because nothing mixes
  momentum vertically; the Monin-Obukhov fluxes then see an unphysically weak wind. Judged a defect of the no-closure
  configuration; the control is rerun with the TKE closure. `SmagorinskyLilly` cannot be used on 50 m × 2 km cells: its
  isotropic filter width (≈585 m) makes the explicit vertical diffusion unstable (νΔt/Δz² > 1/2) in the first step. A
  horizontal-only Smagorinsky with the 1-D TKE closure would be the LES-consistent mesoscale choice, but the pinned
  Breeze cannot combine closures (post-pin #1036); horizontal mixing is the WENO scheme's.
- **Rain to land**: NumericalEarth imports the renamed `surface_precipitation_flux`; the extension defines it at run
  time in `Breeze.AtmosphereModels` as the sign-flipped `bottom_precipitation_flux` so the child's rain reaches the bucket.
- **Nesting**: one-way and offline — the outer run saves its state on the inner region + halo every 10 min, and
  `outer_run_parent` rebuilds a `PrescribedAtmosphere` on a static 3-D-height grid for `tracer_mip_inner_simulation`.

## Runs so far (all on wpcluster; see the status file for job IDs)

| Run | Configuration | Result |
| --- | --- | --- |
| CPU smokes (jobs 204–215) | real ERA5, 32×32×94, 36 s | stable without closure and with TKE closure; Smagorinsky fails in step 1 |
| GPU pilot (job 217) | real ERA5, 750×750×94, 1 h, A100-80GB | completed; 40.3 GB; 4.3 s wall per 3 s step (≈86 min per simulated hour) |
| GPU control (job 233) | same, 24 h, no closure | running; surface-layer wind defect diagnosed (see Turbulence/PBL) |
| GPU control (job 339) | same, 24 h, TKE closure | queued after job 233 (replacement August 7 control) |
| Inner-nest CPU test (job 227) | parent from the pilot's saved state, 8×8×24 child, 6 s | passes |

## Status

See `plans/status/TRACER-MIP.md` in the coordination checkout for the live status, job IDs and results. Implemented and
tested on CPU: protocol metadata, grids, vertical faces, aerosol profiles and P3 modes, the Tier-1 prescribed aerosol in P3,
process-rate accumulators, ERA5 staging, the outer-domain constructor with coastal exchange, and the offline inner nest.
The 24-h August 7 outer control is running; the inner 500 m run, its 10-min/2-min writers and the June 17 case are next.
A pilot is not a completed MIP; nothing here has been compared against observations.
