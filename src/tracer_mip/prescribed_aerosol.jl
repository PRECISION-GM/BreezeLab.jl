#####
##### Tier-1 prescribed aerosol for Breeze P3: height-dependent, fixed at every step.
#####
##### Breeze's `AerosolActivation` describes vertically uniform modes. The TRACER-MIP control
##### prescribes the two modes' number mixing ratios as a common shape times the surface values
##### (`TracerMIPAerosolProfile`), so the Morrison-Grabowski activation spectrum of the whole
##### population at height z is the surface spectrum scaled by f(z) = n(z)/n(0). This type
##### carries a surface `AerosolActivation` plus a static `nᵃ` field holding the prescribed total
##### number mixing ratio at every cell (AGL), and extends P3's aerosol hooks so that:
#####   * droplet number ρnᶜˡ is prognostic (activation, collection, ... as with `AerosolActivation`),
#####   * no reservoir ρnᵃ is advected or depleted ("fixed aerosols at each model step"),
#####   * the activation target is f(z) · N_act(S, T) — capped by the local prescribed number.
#####

using KernelAbstractions: @kernel, @index
using Oceananigans: CenterField, Center, Face
using Oceananigans.Grids: znode
using Oceananigans.Utils: launch!
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Breeze.Microphysics.PredictedParticleProperties: PredictedParticleProperties, AerosolActivation, AerosolMode,
                                                       PredictedParticlePropertiesMicrophysics, total_activated_number
import Breeze.Microphysics.PredictedParticleProperties: cloud_prognostic_names, aerosol_prognostic_names,
    cloud_number_correction_pairs, cloud_number_correction_fields, aerosol_correction_fields,
    cloud_number_fields, aerosol_fields, grid_cloud_droplet_number, grid_aerosol_number,
    write_cloud_number_diagnostic!, write_aerosol_diagnostic!, add_cloud_number_tendency!, add_aerosol_tendency!,
    p3_aerosol_tendency_fields, compute_cloud_droplet_activation, initial_aerosol_number, has_prognostic_aerosol

const P3Scheme = PredictedParticlePropertiesMicrophysics

struct PrescribedAerosolProfile{FT, A, P}
    activation :: A          # AerosolActivation at the surface number mixing ratios (fixed reservoir)
    profile :: P             # TracerMIPAerosolProfile (isbits; evaluated per cell on the device)
    surface_number :: FT     # Σ modes n(0) [kg⁻¹]
end

"""
    PrescribedAerosolProfile(profile::TracerMIPAerosolProfile; FT = Float32, kappa = 0.26, activation_kw...)

Tier-1 aerosol for `P3Microphysics(FT; cloud, aerosol = PrescribedAerosolProfile(...))`: the two
protocol modes at their surface numbers (`tracer_mip_p3_aerosol_modes`) wrapped with the height
profile. `activation_kw` go to `AerosolActivation` (e.g. `activation_timescale`).
"""
function PrescribedAerosolProfile(profile::TracerMIPAerosolProfile; FT = Float32, kappa = 0.26,
                                  thermodynamic_constants = ThermodynamicConstants(FT), activation_kw...)
    modes = tracer_mip_p3_aerosol_modes(profile; FT, kappa, thermodynamic_constants)
    activation = AerosolActivation(modes...; prognostic = false, thermodynamic_constants, activation_kw...)
    surface_number = FT(sum(surface_number_mixing_ratios(profile)))
    profile_FT = convert_profile(FT, profile)
    return PrescribedAerosolProfile{FT, typeof(activation), typeof(profile_FT)}(activation, profile_FT, surface_number)
end

convert_profile(FT, p::TracerMIPAerosolProfile{FT}) = p
function convert_profile(FT, p::TracerMIPAerosolProfile)
    modes = map(m -> TracerMIPAerosolMode{FT}(m.number_cm3, m.median_diameter, m.sigma), p.modes)
    return TracerMIPAerosolProfile{FT, typeof(modes)}(modes, p.reference_density, p.minimum_total_number,
                                                      p.floor_height, FT.(p.shape), p.multiplier)
end

Base.summary(a::PrescribedAerosolProfile) =
    string("PrescribedAerosolProfile(", length(a.activation.modes), " modes, TRACER-MIP height profile, fixed)")
Base.show(io::IO, a::PrescribedAerosolProfile) = print(io, summary(a))

# Prognostic set: droplet number is predicted; there is no aerosol reservoir prognostic.
@inline cloud_prognostic_names(::PrescribedAerosolProfile) = (:ρqᶜˡ, :ρnᶜˡ)
@inline aerosol_prognostic_names(::PrescribedAerosolProfile) = ()
@inline cloud_number_correction_pairs(::PrescribedAerosolProfile, μ) = tuple((μ.ρnᶜˡ, μ.ρqᶜˡ))
@inline cloud_number_correction_fields(::PrescribedAerosolProfile, μ) = tuple(μ.ρnᶜˡ)
@inline aerosol_correction_fields(::PrescribedAerosolProfile, μ) = ()
@inline has_prognostic_aerosol(::PrescribedAerosolProfile) = false
@inline initial_aerosol_number(a::PrescribedAerosolProfile) = a.surface_number
@inline p3_aerosol_tendency_fields(G, ::PrescribedAerosolProfile) = (; G.ρnᶜˡ)

@inline cloud_number_fields(::PrescribedAerosolProfile, grid) = (; ρnᶜˡ = CenterField(grid), nᶜˡ = CenterField(grid))

# The prescribed number mixing ratio is written once, per cell, in meters above ground (the
# bottom face of the column: the terrain height on a terrain-following grid, 0 otherwise).
@kernel function _fill_prescribed_aerosol!(nᵃ, grid, profile)
    i, j, k = @index(Global, NTuple)
    z  = znode(i, j, k, grid, Center(), Center(), Center())
    z₀ = znode(i, j, 1, grid, Center(), Center(), Face())
    @inbounds nᵃ[i, j, k] = total_number_mixing_ratio(profile, z - z₀)
end

function aerosol_fields(a::PrescribedAerosolProfile, grid)
    nᵃ = CenterField(grid)
    launch!(grid.architecture, grid, :xyz, _fill_prescribed_aerosol!, nᵃ, grid, a.profile)
    fill_halo_regions!(nᵃ)
    return (; nᵃ)
end

@inline grid_cloud_droplet_number(p3::P3Scheme, ::PrescribedAerosolProfile, μ, i, j, k, ρ) = @inbounds μ.ρnᶜˡ[i, j, k] / ρ
@inline grid_aerosol_number(::PrescribedAerosolProfile, μ, i, j, k, ρ) = @inbounds μ.nᵃ[i, j, k]

@inline function write_cloud_number_diagnostic!(μ, i, j, k, grid, ::PrescribedAerosolProfile, ℳ)
    @inbounds μ.nᶜˡ[i, j, k] = ℳ.nᶜˡ
    return nothing
end
@inline write_aerosol_diagnostic!(μ, i, j, k, grid, ::PrescribedAerosolProfile, ℳ) = nothing

@inline function add_cloud_number_tendency!(G, i, j, k, grid, ::PrescribedAerosolProfile, result)
    @inbounds G.ρnᶜˡ[i, j, k] += result.tendency_ρnᶜˡ
    return nothing
end
@inline add_aerosol_tendency!(G, i, j, k, grid, ::PrescribedAerosolProfile, result) = nothing

"""
Activation with the surface spectrum scaled by the local prescribed number: the equilibrium
droplet count is `f(z) N_act(S, T)` with `f(z) = nᵃ(z) / nᵃ(0)` (both modes share the shape, so
the per-mode scaling and the total scaling coincide); the droplet number relaxes toward it over
the activation timescale, as in Breeze's `aerosol_activation_rate`.
"""
@inline function compute_cloud_droplet_activation(a::PrescribedAerosolProfile, p3, qᶜˡ, nᶜˡ, nᵃ,
                                                  qᵛ, qᵛ⁺ˡ, T, ρ, constants)
    act = a.activation
    FT = typeof(T)
    S = (qᵛ - qᵛ⁺ˡ) / max(qᵛ⁺ˡ, act.minimum_saturation_mass_fraction)
    scale = max(zero(FT), nᵃ) / a.surface_number
    N_target = scale * total_activated_number(act, T, S)
    ncnuc = max(zero(FT), N_target - nᶜˡ) / act.activation_timescale
    seed_mass = 4 * FT(π) / 3 * act.liquid_water_density * act.activated_droplet_radius^3
    qcnuc = ncnuc * seed_mass
    is_supersaturated = S > act.activation_supersaturation_threshold
    return (; mass = ifelse(is_supersaturated, qcnuc, zero(FT)),
              number = ifelse(is_supersaturated, ncnuc, zero(FT)))
end
