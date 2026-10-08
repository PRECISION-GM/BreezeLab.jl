# ENA experiment definitions

ENA identifies the Eastern North Atlantic site. A case additionally specifies a date,
forcing dataset, grid, physics, and experiment protocol. BreezeLab keeps these choices
separate from its shared model assembly. Changing an option produces a recorded
sensitivity experiment; it does not silently change the protocol's definition.

## Covert et al. (2022)

Sources: [paper](https://acp.copernicus.org/articles/22/1159/2022/),
[public bin inputs](https://github.com/dmechem/ENA_variability_LES_bin_paper), and
[public bulk inputs](https://github.com/dmechem/ENA_variability_LES_bulk_paper).
The fetch script pins each repository revision and records SHA-256 checksums.

The supported `covert_public_bin` protocol follows the runnable public bin input
configuration: 18 July 2017, 06–12 UTC, 256×256×192, 35 m horizontal spacing
(8.96 km square). It uses time-dependent prescribed sensible/latent heat and stress,
simple longwave radiation, full-field large-scale vertical advection, and no
mean-wind nudging. The default Breeze microphysics is P3 with N75 initialization.

The paper describes 864×864×192 (30.24 km square), 06–15 UTC experiments. The public
configuration is therefore a **development benchmark**, not the complete paper setup.
It omits the original vertical-grid file: BreezeLab reconstructs 192 cells, with
10 m spacing to 1.5 km and stretching to 20 km. Recovering and checking the actual
paper grid and experiment variants is required before claiming paper reproduction.
The bulk repository is also downloaded for comparison, but is not a separately
validated protocol in this package.

## LASSO-ENA

See [the ENA-LASSO audit](ena_lasso.md) for the adopted member, the bundle validation
rules and the SAM-vs-Breeze physics differences recorded against the SAM source.

Sources: [ARM bundle browser](https://lasso-ena.svcs.arm.gov/latest/bundle_browser.html),
[LASSO-ENA dataset](https://doi.org/10.5439/2572661), and
[LASSO SAM source](https://code.arm.gov/lasso/lasso-ena-codes/lasso_sam_sbm)
(`lasso_ena_noice`, reference revision `12d02446a2147388dc89d828e6e0553106abea0f`).

Obtain a `samin` bundle through an ARM account and stage it with
`julia data_wrangling/stage_lasso_bundle.jl ARCHIVE.tar`, which extracts `snd`, `lsf`,
`sfc`, `prm` and `grd` into `data/lasso/<run ID>/` and freezes the archive checksum, member
tokens, DOI and SAM revision in `bundle.toml`. The adopted member is
`20170718era5d25x100_sbmwrm-aer2-flxsst` (the plan's `…era5s1n0…` spectral-bin candidate
does not exist at ARM; see [the audit](ena_lasso.md)). Do not substitute the public
Covert inputs.

The current adapter targets warm-cloud, ocean, SST-flux, wind-nudged LASSO members.
[`inspect_lasso_bundle`](ena_lasso.md) reads the grid from `grd`, horizontal spacing and
timestep from `prm`, duration from `nstop*dt`, radiation cadence from `nrad*dt`, the
nudging timescale, `perturb_type` and `compute_reffc`, and rejects every unsupported
namelist setting with its reason. The member run ID, the UTC epoch (SAM's `day0` does not
identify the year) and, when the namelist carries neither `nx_gl,ny_gl,nz_gl` nor a numeric
`Nx×Ny×Nz` case ID, the dimensions must be given explicitly; they are recorded in
provenance and checked against the bundle. Use `Nx`, `Ny`, and `z_faces` overrides for
grid sensitivities.

With the staged bundle, the readable script runs the member end to end:

```sh
julia --project cases/ena_lasso.jl        # member from ENA_LASSO_MEMBER, bundle from ENA_LASSO_BUNDLE
```

or, from the CLI, **after confirming the bundle's date, inputs and grid**:

```sh
julia --project cases/cli/run_case.jl \
  --protocol lasso_ena_official \
  --member 20170718era5d25x100_sbmwrm-aer2-flxsst \
  --data data/lasso/20170718era5d25x100_sbmwrm-aer2-flxsst \
  --epoch 2017-07-18T00:00:00 --dimensions 256,256,260 \
  --arch gpu --float Float32 --output output/ena_lasso
```

(the epoch above is the day-of-year of `day0` in the member's year; the adapter rejects a
mismatch). The API equivalent is:

```julia
using BreezeLab, Oceananigans, Dates
case = ena_lasso(; member="20170718era5d25x100_sbmwrm-aer2-flxsst",
                   epoch=DateTime(2017, 7, 18, 0), arch=GPU())   # dimensions=:documented → 256×256×260
write_provenance("provenance.toml", case)
run!(case.simulation)
```

The adapter defaults to SST-based bulk surface exchange, RRTMGP longwave/shortwave,
mean-wind nudging on the `tauls` timescale, full-field vertical advection, P3-aer2,
and diagnostic CCN replenishment. These defaults describe a Breeze counterpart of
the LASSO forcing protocol; the symbol `lasso_ena_official` identifies the input
protocol, not a certification of equivalent results.

## Shared forcing formulation

For a transported variable φ, the imposed large-scale vertical-advection contribution is

```math
(∂φ/∂t)_{LS,v} = -w_{LS}(z,t) ∂φ/∂z.
```

It acts on the resolved field. Mean-wind nudging contributes a horizontally uniform
acceleration, independently for each horizontal component:

```math
(∂u/∂t)_{nudge} = [u_{LS}(z,t) - ⟨u⟩(z,t)] / τ_{LS}.
```

SAM moisture inputs are mixing ratios r per dry-air mass. In clear air the Breeze
specific-humidity conversion is `q = r/(1+r)`, with tendency
`dq/dt = (dr/dt)/(1+r)^2`. In cloudy air, holding condensate mass fraction qᶜ fixed,
the vapor tendency is `(qᵈ)^2/(1-qᶜ) * drᵛ/dt`, where qᵈ is dry-air mass fraction.
Thermodynamic forcing is converted to the chosen Breeze
prognostic variable while preserving the specified physical temperature tendency.
See `large_scale_forcings.jl` and its tests for the implemented conversions.

## Known differences to resolve or quantify

- Breeze P3/one-moment microphysics differs from SAM's HUJI bin microphysics.
- The default subgrid closure/advection are Smagorinsky–Lilly/WENO; they differ from
  the SAM reference's TKE closure and transport.
- RRTMGP columns currently stop at the LES top and use prescribed effective radii.
  Matching the reference requires checking the atmospheric column above the LES,
  radiation state, optical assumptions, and heating/flux profiles.
- Bulk surface exchange differs from SAM's iterative stability/gustiness treatment.
  The current bulk-flux path uses a zero translation velocity to retain
  ground-relative surface winds; the SAM translation is recorded in the case label.
- The diagnostic CCN reservoir projection is applied once per model timestep,
  rather than at every microphysics call as in the reference algorithm.
- Converting SAM scalar levels to Breeze finite-volume cells needs comparison with
  the selected bundle's staggered grid, especially near the surface and model top.

## LASSO reproduction acceptance checklist

1. **Freeze one ensemble member.** Archive the bundle ID, DOI, input archive checksum,
   SAM revision, and reference outputs. Confirm date, spin-up/evaluation windows,
   `day0`, dimensions, scalar/face heights, aerosols, SST, nudging, and radiation flags.
2. **Match imposed tendencies.** Compare interpolated sounding/forcing and surface
   time series at selected times. Check units/signs, wind frame, pressure/height
   conversion, column mass/water budgets, and horizontal/vertical advection separately.
3. **Check each physics component.** Compare surface H/LE/stress and radiation
   fluxes/heating against SAM before the full coupled run. Resolve the differences
   above or label them as controlled model differences with sensitivity tests.
4. **Run a short GPU test, then the complete member.** Inspect finite values,
   timestep stability, negative water/number concentrations, budgets, and output
   completeness. Preserve both Float32 and Float64 short-run checks where practical.
5. **Evaluate statistically.** Compare LWP, cloud fraction and boundaries,
   precipitation, mean thermodynamic/wind profiles, and surface fluxes over common
   sampling windows against SAM and ARM observations. Agree on tolerances and
   uncertainty estimates before declaring reproduction; matching chaotic snapshots
   is not the criterion.

The automated protocol test uses a synthetic SAM-style fixture and a tiny sensitivity
with one-moment microphysics, radiation disabled, and no aerosol replenishment. The public-input CPU smoke also
uses a reduced grid/duration. Neither establishes full LASSO or Covert scientific
fidelity, and no official ARM bundle is checked into this repository.
