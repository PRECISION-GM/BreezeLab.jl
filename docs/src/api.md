# Case API

Constructors return a configured case without advancing it. A caller can inspect or
modify the simulation, attach additional diagnostics, and then call `run!`.

```@docs
eastern_north_atlantic
ena_simulation
ena_protocol_settings
ena_lasso
inspect_lasso_bundle
validate_lasso_bundle
parse_lasso_member
sea_starr
tracer_dp_scream
```
