"""
    BreezeLab

Observation-comparable LES experiments built on Breeze.jl. Each case has one constructor that
reads its inputs and builds the grid, model and `Simulation` in its own body:

- [`ena_covert`](@ref): the public Covert et al. (2022) ENA development benchmark, 18 July 2017
- [`ena_lasso`](@ref): an official LASSO-ENA member from its staged SAM bundle
- [`tracer_dp_scream`](@ref): the doubly periodic TRACER case driven by the DP-SCREAM IOP forcing
- [`sea_starr`](@ref): the SEA STARR CTRL / N100 / N030 members from their DEPHY drivers
- `tracer_mip_outer_simulation`, `tracer_mip_inner_simulation`: the regional TRACER-MIP
  control (NumericalEarth extension)

Microphysics is passed as a Breeze object, and precision follows
`Oceananigans.defaults.FloatType`. The shared physics components live in `src/components`.
"""
module BreezeLab

export
    # case constructors
    ena_covert, ena_lasso, tracer_dp_scream, sea_starr, tracer_mip_outer_simulation, tracer_mip_inner_simulation,
    # constructor inputs: vertical grids, perturbation, sponge, aerosol
    covert_public_bin_vertical_faces, covert_inversion_refined_vertical_faces, ena_vertical_faces,
    lasso_ena_vertical_faces, uniform_then_stretched_faces,
    tracer_dp_scream_vertical_faces, sea_starr_vertical_faces, acpc_vertical_faces,
    InitialPerturbation, SAMSponge, neutral_drag_coefficient,
    lasso_aerosol, covert_aerosol, first_level_reference_density, activated_fraction,
    kappa_aerosol_activation, tracer_mip_aerosol_profile, PrescribedAerosolProfile,
    # case inputs and their validation
    read_sam_sounding, read_sam_large_scale_forcing, read_sam_surface_forcing, read_sam_namelist, read_sam_grd,
    LassoMember, parse_lasso_member, inspect_lasso_bundle, validate_lasso_bundle, lasso_bundle_directory,
    read_iop_forcing, tracer_dp_scream_settings, read_dephy_driver, tracer_mip_protocol, tracer_mip_case_window,
    # observations and run output
    ARMSeries, read_arm_lwp, read_arm_rain_rate, read_arm_cloud_boundaries, window_statistics,
    cloud_fraction_from_boundaries, breezelab_timeseries, breezelab_cloud_boundaries,
    read_dp_scream_output, dp_scream_window,
    # diagnostics
    cloud_liquid, rain_mass_fraction, liquid_water_path, ice_water_path, precipitable_water, total_condensate,
    cloud_fraction, cloud_fraction_profile, total_cloud_fraction_profile,
    surface_rain_flux, surface_ice_flux, cloud_boundaries, ProgressMessenger,
    ProcessRateAccumulators, process_rate_output_callback, PROCESS_RATE_NAMES, PROCESS_RATE_UNITS,
    # provenance
    write_provenance

using Dates: Dates, DateTime
using TOML: TOML
using Printf: @sprintf
using Oceananigans
using Oceananigans.Units
using Oceananigans.Fields: interior
using Oceananigans.Grids: znodes
using CUDA # Register Oceananigans' GPU() constructor for the default case architecture.
using Breeze
using Breeze.Microphysics.PredictedParticleProperties: PredictedParticlePropertiesMicrophysics, CloudDroplets,
                                                       AerosolActivation, AerosolMode, has_prognostic_aerosol,
                                                       activated_number
using CloudMicrophysics: CloudMicrophysics   # loads Breeze's one-moment extension
using RRTMGP: RRTMGP                         # loads Breeze's RRTMGP extension
using NCDatasets: NCDatasets
using ClimaComms: ClimaComms

"""
    package_path(parts...)

A path inside the BreezeLab checkout (e.g. `package_path("data", "covert2022_bin")`).
"""
package_path(parts...) = joinpath(dirname(@__DIR__), parts...)

# Shared building blocks: SAM inputs, grids, forcing operators, surface fluxes, radiation,
# initial state, microphysics helpers, diagnostics, output and provenance.
include("components/sam_input_files.jl")
include("components/vertical_grid.jl")
include("components/forcing_profiles.jl")
include("components/large_scale_forcings.jl")
include("components/simple_longwave_radiation.jl")
include("components/surface_fluxes.jl")
include("components/initial_conditions.jl")
include("components/microphysics.jl")
include("components/advection.jl")
include("components/diagnostics.jl")
include("components/energy_budget_diagnostics.jl")
include("components/output.jl")
include("components/provenance.jl")
include("components/checkpointing.jl")

# Observations and reference model output.
include("observations/arm_observations.jl")
include("observations/dp_scream_output.jl")

# Cases: one constructor each.
include("ena/covert.jl")
include("ena/lasso_bundle.jl")
include("ena/lasso.jl")
include("tracer_dp_scream/iop_forcing.jl")
include("tracer_dp_scream/tracer_dp_scream.jl")
include("sea_starr/dephy_driver.jl")
include("sea_starr/aerosol.jl")
include("sea_starr/forcings.jl")
include("sea_starr/radiation.jl")
include("sea_starr/sea_starr.jl")
include("tracer_mip/protocol.jl")
include("tracer_mip/aerosol_profiles.jl")
include("tracer_mip/prescribed_aerosol.jl")
include("tracer_mip/process_diagnostics.jl")

# Regional (ERA5-nested, land/sea-coupled) TRACER-MIP constructors: implemented in the NumericalEarth extension.
function tracer_mip_outer_simulation end
function tracer_mip_inner_simulation end

end # module
