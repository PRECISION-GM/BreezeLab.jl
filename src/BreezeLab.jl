"""
    BreezeLab

Observation-comparable LES experiments and workflows built on Breeze.jl.
The initial ENA implementation ports the LASSO and Covert forcing pathways: SAM input-file readers, the large-scale forcing operators
(time-varying geostrophic wind, domain-mean wind nudging, full-field upwind vertical
advection, thermodynamically consistent horizontal advective tendencies, SAM sponge),
simple/RRTMGP radiation, prescribed or bulk surface fluxes, and the 18 July 2017 case
driver with one-moment and P3 microphysics.
"""
module BreezeLab

export
    # SAM input files
    SAMSounding, SAMLargeScaleForcing, SAMSurfaceForcing,
    read_sam_sounding, read_sam_large_scale_forcing, read_sam_surface_forcing, read_sam_namelist,
    sam_hydrostatic_heights, record_heights, interpolate_profile, day_to_seconds, file_sha256, mass_fraction_from_mixing_ratio,
    # grids and profiles
    lasso_ena_cell_centers, lasso_ena_vertical_faces, uniform_vertical_faces, uniform_then_stretched_faces,
    LargeScaleForcingProfiles, profile_time_series, surface_time_series, sam_interpolate_column,
    # forcings
    TimeVaryingGeostrophicForcing, time_varying_geostrophic_forcings,
    MeanProfileNudging, LargeScaleVerticalAdvection, large_scale_thermodynamic_forcings,
    LargeScaleEnergyForcing, LargeScaleMoistureForcing, SAMSponge, sam_sponge_rates,
    upper_boundary_relaxation_forcings, UpperBoundaryEnergyRelaxation, UpperBoundaryMoistureRelaxation, SoundingTargetProfiles,
    SimpleLongwaveRadiation,
    # surface
    prescribed_surface_flux_boundary_conditions, bulk_surface_flux_boundary_conditions, sam_oceflx_neutral_polynomials,
    prescribed_heat_flux_boundary_conditions,
    PrescribedStressUpdater, prescribed_stress_updater,
    SeaSurfaceTemperatureUpdater,
    # initial state
    SoundingProfiles, saturation_partition, InitialPerturbation, initial_state_columns, perturbation_array,
    # case
    eastern_north_atlantic, ena_simulation, ena_protocol_settings, lasso_ena_simulation, build_case, lasso_aerosol_modes, covert_aerosol_modes, activated_fraction, write_provenance, read_sam_grd, faces_from_centers, epoch_from_day_of_year,
    covert_public_bin_vertical_faces, AerosolReplenishment, DiagnosticCCNProjection, initial_sounding, perturbation_amplitudes,
    # LASSO-ENA bundle and case
    LassoMember, parse_lasso_member, lasso_variant_tokens, lasso_samin_filename, lasso_reference_filenames,
    read_sam_namelist_groups, lasso_scalar_levels, inspect_lasso_bundle, validate_lasso_bundle, lasso_bundle_record,
    ena_lasso, lasso_bundle_directory, lasso_bundle_available, lasso_bundle_missing_message, lasso_documented_dimensions,
    # ARM observations and run output readers
    ARMSeries, read_arm_lwp, read_arm_rain_rate, read_arm_cloud_boundaries, window_statistics,
    cloud_fraction_from_boundaries, breezelab_timeseries, breezelab_cloud_boundaries,
    read_case_inputs,
    # TRACER–DP-SCREAM (IOP forcing, archived reference outputs, periodic case)
    IOPForcing, read_iop_forcing, iop_datetimes, iop_index, iop_sam_inputs, iop_surface_albedo, iop_column_integral,
    iop_potential_temperature, fractional_day_of_year,
    read_dp_scream_output, dp_scream_window, parse_cf_time_units,
    tracer_dp_scream, tracer_dp_scream_settings, tracer_dp_scream_vertical_faces, neutral_drag_coefficient,
    # diagnostics
    cloud_liquid, rain_mass_fraction, liquid_water_path, ice_water_path, precipitable_water, total_condensate,
    cloud_fraction, cloud_fraction_profile, total_cloud_fraction_profile,
    surface_rain_flux, surface_ice_flux, cloud_boundaries, ProgressMessenger,
    # SEA STARR
    DEPHYDriver, read_dephy_driver, aerosol_number_per_kg, driver_profile_time_series, driver_initial_profile,
    inversion_height, sea_starr_vertical_faces,
    KappaAerosolMode, kappa_aerosol_activation, kappa_critical_supersaturation,
    SurfaceAerosolSource, EvaporationRegeneration, EvaporationRateUpdater, cloud_evaporation_rate_field, aerosol_number_columns,
    InversionFollowingNudging, InversionMaskUpdater, inversion_height_field, nudging_mask_weights,
    sea_starr, sea_starr_driver_path,
    # TRACER-MIP (regional nested coastal control; the regional constructor lives in the NumericalEarth extension)
    tracer_mip_protocol, tracer_mip_case_window, acpc_vertical_faces, tracer_mip_horizontal_extent, tracer_mip_grid,
    TracerMIPAerosolProfile, tracer_mip_aerosol_profile, aerosol_shape, surface_number_mixing_ratios,
    mode_number_mixing_ratios, total_number_mixing_ratio, aerosol_profile_table, tracer_mip_p3_aerosol_modes,
    PrescribedAerosolProfile,
    ProcessRateAccumulators, accumulate_process_rates!, reset_process_accumulators!, process_rate_output_callback,
    PROCESS_RATE_NAMES, PROCESS_RATE_UNITS, tracer_mip_outer_simulation, tracer_mip_inner_simulation

using Oceananigans
using Oceananigans.Units
using CUDA # Register Oceananigans' GPU() constructor for the default case architecture.
using Breeze

include("sam_input_files.jl")
include("vertical_grid.jl")
include("forcing_profiles.jl")
include("large_scale_forcings.jl")
include("simple_longwave_radiation.jl")
include("surface_fluxes.jl")
include("initial_conditions.jl")
include("diagnostics.jl")
include("checkpointing.jl")
include("case_setup.jl")
include("ena_protocols.jl")
include("eastern_north_atlantic.jl")
include("lasso_bundle.jl")
include("ena_lasso.jl")
include("arm_observations.jl")
include("iop_forcing.jl")
include("dp_scream_output.jl")
include("tracer_dp_scream.jl")
include("dephy_driver.jl")
include("sea_starr_aerosol.jl")
include("sea_starr_forcings.jl")
include("sea_starr.jl")
include("tracer_mip/protocol.jl")
include("tracer_mip/aerosol_profiles.jl")
include("tracer_mip/prescribed_aerosol.jl")
include("tracer_mip/process_diagnostics.jl")

# Regional (ERA5-nested, land/sea-coupled) constructor: implemented in the NumericalEarth extension.
function tracer_mip_outer_simulation end
function tracer_mip_inner_simulation end

end # module
