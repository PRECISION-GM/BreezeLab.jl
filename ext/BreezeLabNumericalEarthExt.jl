module BreezeLabNumericalEarthExt

# Regional (nested, ERA5-driven, land/sea-coupled) TRACER-MIP constructors. Active only when
# NumericalEarth is loaded (the case sub-environment `cases/tracer_mip/Project.toml`); the package
# environment itself does not depend on NumericalEarth. Everything here extends NumericalEarth's
# public constructors with BreezeLab-owned types and callbacks; no NumericalEarth or Breeze method
# is redefined.

using BreezeLab
using BreezeLab: tracer_mip_protocol, tracer_mip_case_window, tracer_mip_grid, tracer_mip_aerosol_profile,
                 PrescribedAerosolProfile, ProcessRateAccumulators, scalar_advection_schemes
using NumericalEarth
using NumericalEarth: BoundingBox, Metadata, Metadatum, SlabLand, AtmosphereLandModel, ETOPO2022,
                      NestedModel, nested_atmosphere_model
using NumericalEarth.DataWrangling: default_horizontal_padding
using NumericalEarth.DataWrangling.ERA5: ERA5HourlyPressureLevels, ERA5HourlySingleLevel, ERA5HourlyLand
using NumericalEarth.Atmospheres: PrescribedAtmosphere
using Oceananigans
using Oceananigans.Units
using Oceananigans.Units: Time
using Oceananigans.Fields: interior, interpolate!
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Grids: znode, znodes
using Oceananigans.Simulations: Callback
using Oceananigans.Utils: launch!
using KernelAbstractions: @kernel, @index
using Breeze
using Breeze.Microphysics.PredictedParticleProperties: CloudDroplets
using Dates: DateTime, Hour, Second
using Printf: @sprintf
using Statistics: mean

# Rain-to-land shim. NumericalEarth d07eb240's Breeze extension imports
# `Breeze.AtmosphereModels.surface_precipitation_flux`, renamed `bottom_precipitation_flux` before the
# pinned Breeze, so the binding is undefined and `applicable(surface_precipitation_flux, …)` throws at
# coupling time. Supply it: Breeze's diagnostic is the bottom-face advective flux (positive upward),
# NumericalEarth's `Jʳⁿ` is positive downward (kg m⁻² s⁻¹), hence the sign flip. Remove once upstream
# follows the rename. The child's rain then reaches the slab bucket as intended.
surface_precipitation_flux_shim(model, microphysics) =
    Field(-1 * Breeze.AtmosphereModels.bottom_precipitation_flux(model, microphysics).operand)

function __init__()
    # Define the name where NumericalEarth's import points (`Breeze.AtmosphereModels`), at run time, so
    # it works regardless of extension load order and is never baked into a precompile image.
    AM = Breeze.AtmosphereModels
    if !isdefined(AM, :surface_precipitation_flux)
        # NumericalEarth probes `applicable(f, model, microphysics)` and then calls `f(model)`: provide both forms.
        Core.eval(AM, :(surface_precipitation_flux(model, microphysics) = $(surface_precipitation_flux_shim)(model, microphysics)))
        Core.eval(AM, :(surface_precipitation_flux(model) = $(surface_precipitation_flux_shim)(model, model.microphysics)))
    end
    return nothing
end

include("tracer_mip_regional/coastal_surface.jl")
include("tracer_mip_regional/synthetic_parent.jl")
include("tracer_mip_regional/outer_domain.jl")
include("tracer_mip_regional/inner_domain.jl")

end # module
