#####
##### Surface turbulent fluxes, radiative fluxes at the surface and the LES top, and the
##### in-cloud droplet number: the series SAM writes as SHF/LHF, LWNS/SWNS, LWNT/SWNT and
##### the in-cloud Nc statistic, needed to attribute differences against the LASSO runs.
#####

using Oceananigans: Field, Average, Center, Face, KernelFunctionOperation
using Oceananigans.BoundaryConditions: getbc
using Oceananigans.Models: boundary_condition_args

#####
##### Surface sensible and latent heat fluxes from the bottom boundary conditions
#####

# The fluxes are evaluated from the model's own materialized bottom boundary conditions with
# the arguments Breeze passes them, so they are exactly the fluxes the tendencies receive
# (bulk SST formulas for LASSO flxsst, the prescribed sfc series for Covert).
@inline function surface_latent_heat_flux_kernel(i, j, k, grid, bcq, clock, fields, dynamics_fields, ℒ)
    return ℒ * getbc(bcq, i, j, grid, clock, fields, dynamics_fields)
end

# Breeze converts an enthalpy flux 𝒬ᵀ = ρ cᵖᵐ w'T' to the θ flux as Jᶿ = 𝒬ᵀ / (cᵖᵐ Π); the
# inverse with vapor-only mixture properties at the first cell gives the sensible heat flux.
@inline function surface_sensible_heat_flux_kernel(i, j, k, grid, bcθ, clock, fields, dynamics_fields,
                                                   cᵖᵈ, cᵖᵛ, Rᵈ, Rᵛ, p₀)
    Jᶿ = getbc(bcθ, i, j, grid, clock, fields, dynamics_fields)
    qᵛ = @inbounds fields.qᵛ[i, j, 1]
    p = @inbounds dynamics_fields.p[i, j, 1]
    cᵖᵐ = (1 - qᵛ) * cᵖᵈ + qᵛ * cᵖᵛ
    Rᵐ = (1 - qᵛ) * Rᵈ + qᵛ * Rᵛ
    Π = (p / p₀)^(Rᵐ / cᵖᵐ)
    return cᵖᵐ * Π * Jᶿ
end

function moisture_prognostic_name(fields)
    for name in (:ρqᵛ, :ρqᵗ, :ρqᵉ)
        haskey(fields, name) && return name
    end
    throw(ArgumentError("surface_heat_fluxes needs a ρqᵛ, ρqᵗ or ρqᵉ prognostic"))
end

"""
    surface_heat_fluxes(model)

2D fields `(; sensible, latent)` of the upward surface sensible and latent heat fluxes
[W m⁻²], evaluated from the model's bottom boundary conditions (the fluxes the model
actually applies): the latent flux is ℒ E with Breeze's reference latent heat and the
sensible flux is cᵖᵐ Π Jᶿ at the first cell of the liquid-ice potential-temperature model.
For prescribed fluxes these are the sfc file's H and LE.
"""
function surface_heat_fluxes(model)
    grid = model.grid
    FT = eltype(grid)
    prognostic = Oceananigans.prognostic_fields(model)
    haskey(prognostic, :ρθ) || throw(ArgumentError("surface_heat_fluxes needs the liquid-ice potential-temperature prognostic ρθ"))
    bcθ = prognostic.ρθ.boundary_conditions.bottom
    bcq = prognostic[moisture_prognostic_name(prognostic)].boundary_conditions.bottom
    clock, fields, dynamics_fields = boundary_condition_args(model)
    constants = model.thermodynamic_constants
    ℒ = FT(constants.liquid.reference_latent_heat)
    latent = KernelFunctionOperation{Center, Center, Nothing}(surface_latent_heat_flux_kernel, grid,
                                                              bcq, clock, fields, dynamics_fields, ℒ)
    cᵖᵈ = FT(constants.dry_air.heat_capacity)
    cᵖᵛ = FT(constants.vapor.heat_capacity)
    Rᵈ = FT(constants.molar_gas_constant / constants.dry_air.molar_mass)
    Rᵛ = FT(constants.molar_gas_constant / constants.vapor.molar_mass)
    p₀ = FT(model.dynamics.reference_state.standard_pressure)
    sensible = KernelFunctionOperation{Center, Center, Nothing}(surface_sensible_heat_flux_kernel, grid, bcθ,
                                                                clock, fields, dynamics_fields, cᵖᵈ, cᵖᵛ, Rᵈ, Rᵛ, p₀)
    return (; sensible = Field(sensible), latent = Field(latent))
end

#####
##### Radiative fluxes at the surface and at the LES top
#####

@inline face_level_sum(i, j, k, grid, a, b, level, sign) = @inbounds sign * (a[i, j, level] + b[i, j, level])
@inline face_level_value(i, j, k, grid, a, level) = @inbounds a[i, j, level]

"""
    radiative_boundary_fluxes(radiation, grid)

2D fields of the net radiative fluxes [W m⁻²] at the surface (bottom face) and at the LES top
(face Nz+1; SAM's LWNT/SWNT are at the top of its extended radiation column instead), in SAM's
sign conventions: longwave net positive upward (LWNS, LWNT), shortwave net positive downward
(SWNS, SWNT), plus the downwelling shortwave at the LES top. RRTMGP models give all five;
[`SimpleLongwaveRadiation`](@ref) gives the two longwave nets; `nothing` gives an empty tuple.
"""
radiative_boundary_fluxes(::Nothing, grid) = NamedTuple()

function radiative_boundary_fluxes(radiation::SimpleLongwaveRadiation, grid)
    Nz = size(grid, 3)
    F = radiation.flux
    net(level) = Field(KernelFunctionOperation{Center, Center, Nothing}(face_level_value, grid, F, level))
    return (; surface_net_longwave = net(1), top_net_longwave = net(Nz + 1))
end

function radiative_boundary_fluxes(radiation, grid)
    hasproperty(radiation, :upwelling_longwave_flux) || return NamedTuple()
    Nz = size(grid, 3)
    FT = eltype(grid)
    # Breeze stores upwelling fluxes positive and downwelling fluxes negative (net = up + dn, positive up)
    net(up, dn, level, sign) = Field(KernelFunctionOperation{Center, Center, Nothing}(face_level_sum, grid, up, dn,
                                                                                      level, FT(sign)))
    lw_up, lw_dn = radiation.upwelling_longwave_flux, radiation.downwelling_longwave_flux
    sw_up, sw_dn = radiation.upwelling_shortwave_flux, radiation.downwelling_shortwave_flux
    return (; surface_net_longwave = net(lw_up, lw_dn, 1, 1),
              surface_net_shortwave = net(sw_up, sw_dn, 1, -1),
              top_net_longwave = net(lw_up, lw_dn, Nz + 1, 1),
              top_net_shortwave = net(sw_up, sw_dn, Nz + 1, -1),
              top_downwelling_shortwave = Field(KernelFunctionOperation{Center, Center, Nothing}(face_level_sum, grid,
                                                                       sw_dn, sw_dn, Nz + 1, FT(-1 / 2))))
end

#####
##### In-cloud droplet number
#####

@inline cloudy_number_density(i, j, k, grid, ρ, n, q, threshold) =
    @inbounds ifelse(q[i, j, k] > threshold, ρ[i, j, k] * n[i, j, k], zero(eltype(grid)))

@inline safe_ratio(i, j, k, grid, a, b) = @inbounds ifelse(b[i, j, k] > 0, a[i, j, k] / max(b[i, j, k], eps(eltype(grid))),
                                                           zero(eltype(grid)))

"""
    in_cloud_droplet_number(model; threshold=1e-5)

Statistics of the cloud droplet number concentration ρᵣ nᶜˡ [m⁻³] over cloudy cells
(`qᶜˡ > threshold`, 0.01 g kg⁻¹ by default: the Covert et al. 2022 and `cloud_fraction_profile`
threshold), for schemes with a prognostic `nᶜˡ`; `nothing` otherwise. Returns
`(; mean, cloudy_fraction, profile)`: the volume-weighted in-cloud domain mean (0 when there is
no cloud), the cloudy volume fraction, and the horizontal mean of ρᵣ nᶜˡ over cloudy cells per
level (divide by the cloud-fraction profile for the in-cloud mean at each level).
"""
function in_cloud_droplet_number(model; threshold=1e-5)
    μ = model.microphysical_fields
    haskey(μ, :nᶜˡ) || return nothing
    grid = model.grid
    FT = eltype(grid)
    ρ = reference_density(model)
    cloudy_n = KernelFunctionOperation{Center, Center, Center}(cloudy_number_density, grid, ρ, μ.nᶜˡ,
                                                               cloud_liquid(model), FT(threshold))
    cloudy = KernelFunctionOperation{Center, Center, Center}(cell_indicator, grid, cloud_liquid(model), FT(threshold))
    numerator = Field(Average(cloudy_n))
    cloudy_fraction = Field(Average(cloudy))
    mean = Field(KernelFunctionOperation{Nothing, Nothing, Nothing}(safe_ratio, grid, numerator, cloudy_fraction))
    return (; mean, cloudy_fraction, profile = Field(Average(cloudy_n, dims=(1, 2))))
end
