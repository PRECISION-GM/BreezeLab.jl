# Case API

Each case has one constructor. It reads the case inputs and builds the grid, the
`AtmosphereModel` and the `Simulation` in its own body, then returns the case without
advancing it. A caller can inspect or modify the simulation, attach additional
diagnostics, and then call `run!`.

Microphysics is passed as a Breeze object, and precision follows
`Oceananigans.defaults.FloatType`.

```@docs
ena_covert
ena_lasso
tracer_dp_scream
sea_starr
```

## Constructor inputs

```@docs
lasso_aerosol
covert_aerosol
first_level_reference_density
kappa_aerosol_activation
inspect_lasso_bundle
validate_lasso_bundle
parse_lasso_member
write_provenance
```
