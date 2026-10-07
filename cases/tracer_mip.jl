# TRACER-MIP Tier 1 case specification — MODEL ASSEMBLY NOT IMPLEMENTED YET.
# References: data_wrangling/fetch_manifest.jl cases/tracer_mip/inputs.toml data/tracer_mip_416423a
# This requires regional nesting and a land surface, not a periodic ENA configuration.
using Dates

tracer_mip = (
    start_time = DateTime(2022, 8, 7, 6),    # second control case: June 17, 2022, 06 UTC
    duration_hours = 24,
    latitude = 29.4719, longitude = -95.0792,
    outer = (Nx=750, Ny=750, Δx=2000.0, Δy=2000.0, Δt=3.0),
    inner = (Nx=500, Ny=500, Δx=500.0, Δy=500.0, Δt=1.5),
    aerosol_tier = 1,                       # prescribe aerosols each step; evolve droplets
    aerosol_radiation = false,
    aerosol_multiplier = 1.0,               # start with CTRL; low-case factor is unresolved
    radiation_interval = 60.0,
    outer_output_interval = 3600.0,
    inner_output_interval = 600.0,
    tracking_output_interval = 120.0,       # 17–01 UTC cell-tracking window
)

# Assembly still needs ERA5 forcing, one-way outer → inner nesting, terrain/land
# exchange, the prescribed two-mode aerosol profiles, and microphysical-process exports.
# Use the protocol's vertical levels and notebook profiles; do not infer them from ENA.
if abspath(PROGRAM_FILE) == @__FILE__
    error("TRACER-MIP is a specification only: its regional model adapter is not implemented yet.")
end
