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

include("tracer_mip_regional/coastal_surface.jl")
include("tracer_mip_regional/synthetic_parent.jl")
include("tracer_mip_regional/outer_domain.jl")

end # module
