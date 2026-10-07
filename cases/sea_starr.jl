# SEA STARR case specification — MODEL ASSEMBLY NOT IMPLEMENTED YET.
# Official inputs: data_wrangling/fetch_manifest.jl cases/seastarr/inputs.toml data/seastarr_22241697
# This file records the intended experiment, not a runnable replacement using ENA physics.
using Dates

sea_starr = (
    member = :CTRL,                         # :CTRL, :N100, :N030; each has its own driver
    start_time = DateTime(2017, 8, 15, 21),
    stop_time = DateTime(2017, 8, 18, 15),    # 66 hours including 3 hours spinup
    Nx = 384, Ny = 384,                     # 192 × 192 is also permitted by the protocol
    Δx = 50.0, Δy = 50.0,
    Δz_lower = 10.0, uniform_top = 2500.0,
    vertical_stretch = 1.1, domain_top = 6500.0,
    Δt = 1.0,
    thermodynamic_nudging_timescale = 1800.0,
    wind_nudging_timescale = 10800.0,
    aerosol_surface_source = 7e5,            # number m⁻² s⁻¹
)

# Assembly still needs DEPHY ingestion, moving-inversion nudging, trajectory SST,
# absorbing-aerosol radiation and aerosol source/scavenging/regeneration budgets.
# ENA's diagnostic_ccn reset is inappropriate for this interactive aerosol experiment.
# See docs/production-readiness.md for the protocol questions and output requirements.
if abspath(PROGRAM_FILE) == @__FILE__
    error("SEA STARR is a specification only: its model adapter is not implemented yet.")
end
