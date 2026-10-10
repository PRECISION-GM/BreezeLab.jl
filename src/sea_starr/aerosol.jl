#####
##### SEA STARR interactive aerosol: κ-Köhler activation from the local prognostic reservoir,
##### a prescribed surface number source, and aerosol regeneration on droplet evaporation.
#####
##### What Breeze's P3 provides (written against 0.11.3 / 5264d3c; unchanged in 0.12.0 except as noted), and what
##### is added here:
#####
#####   * Breeze `AerosolActivation(...; prognostic=true)` carries the unactivated reservoir
#####     `ρnᵃ` and removes one aerosol per activated droplet. Its activated number is the
#####     *constant* mode `number_mixing_ratio × f(S)` only *capped* by the local reservoir,
#####     so a height-varying or depleted `nᵃ` would not scale activation. The method
#####     `aerosol_activation_rate` is extended below for a BreezeLab-owned mode type
#####     ([`KappaAerosolMode`](@ref)) so that the activated target is `f(S) (nᶜˡ + nᵃ)` with
#####     the κ-Köhler critical supersaturation. (A Breeze function specialized on a
#####     BreezeLab type: an extension, not piracy.)
#####   * P3 has no surface aerosol flux. Before Breeze 0.12 a `FluxBoundaryCondition` under `ρnᵃ`
#####     was accepted but ignored (P3 materialized its fields with default no-flux conditions), so
#####     the source is applied as a bottom-cell forcing ([`SurfaceAerosolSource`](@ref)). Breeze
#####     0.12 (PR 1006) applies such a condition; the forcing is kept until the switch is made and
#####     checked separately (docs/breeze_0.12_upgrade.md): `F / Δz₁` either way.
#####   * P3 never returns aerosol: droplets that evaporate carry their aerosol away
#####     (`ρnᶜˡ` is zeroed when `ρqᶜˡ` is exhausted) and the reservoir only decreases.
#####     [`EvaporationRegeneration`](@ref) removes droplet number in proportion to the
#####     cloud-water evaporation rate that P3 itself diagnoses, and returns the same number
#####     to `nᵃ`. Coalescence (autoconversion, accretion, self-collection) reduces droplet and
#####     rain number without regeneration — that is the collision scavenging of aerosol — and
#####     rain reaching the surface removes the aerosol it carries. Rain evaporation does not
#####     regenerate aerosol here (nʳ ~ 10⁴ m⁻³ ≪ nᵃ ~ 10⁸ m⁻³; labelled approximation).
#####
##### Nothing here reuses ENA's `DiagnosticCCNProjection`, which would reset the reservoir
##### to its initial value every step.
#####

using Adapt: Adapt, adapt
using Oceananigans: Field, KernelFunctionOperation, compute!
using Oceananigans.Grids: Center
using Oceananigans.Operators: Δzᶜᶜᶜ
using Oceananigans.Architectures: on_architecture, architecture
using Breeze.AtmosphereModels: AtmosphereModels
using Breeze.Microphysics.PredictedParticleProperties: PredictedParticleProperties, AerosolActivation,
                                                        aerosol_activation_rate, activated_number, sum_aerosol_number
using Breeze.Microphysics.PredictedParticleProperties: compute_p3_process_rates, p3_process_properties, p3_adiabatic_temperature_tendency

const PPP = PredictedParticleProperties

#####
##### κ-Köhler lognormal mode with activation from the local reservoir
#####

"""
    KappaAerosolMode(FT=Oceananigans.defaults.FloatType; number_mixing_ratio, mean_radius, geometric_std, kappa)

One lognormal aerosol mode described by its geometric mean radius, geometric standard
deviation and hygroscopicity parameter κ (Petters & Kreidenweis 2007). `number_mixing_ratio`
[kg⁻¹] is only the value P3 seeds the reservoir with before `set!(model; nᵃ=...)`; activation
uses the local prognostic reservoir instead. SEA STARR: diameter 185 nm (`mean_radius =
92.5e-9`), σ = 1.5, κ = 0.2.
"""
struct KappaAerosolMode{FT}
    number_mixing_ratio :: FT
    mean_radius :: FT
    geometric_std :: FT
    kappa :: FT
end

KappaAerosolMode(FT=Oceananigans.defaults.FloatType; number_mixing_ratio, mean_radius=92.5e-9, geometric_std=1.5, kappa=0.2) =
    KappaAerosolMode{FT}(FT(number_mixing_ratio), FT(mean_radius), FT(geometric_std), FT(kappa))

Base.summary(::KappaAerosolMode) = "KappaAerosolMode"
Base.show(io::IO, m::KappaAerosolMode) =
    print(io, "KappaAerosolMode(rm=", m.mean_radius, " m, σg=", m.geometric_std, ", κ=", m.kappa, ", seed Na=", m.number_mixing_ratio, " kg⁻¹)")

const KappaAerosolActivation{FT, P} = AerosolActivation{FT, P, <:Tuple{Vararg{KappaAerosolMode{FT}}}}

"""
    kappa_aerosol_activation(FT=Oceananigans.defaults.FloatType; number_mixing_ratio, mean_radius=92.5e-9, geometric_std=1.5, kappa=0.2,
                             thermodynamic_constants=ThermodynamicConstants(FT), kwargs...)

A prognostic-reservoir `AerosolActivation` built on one [`KappaAerosolMode`](@ref). The
remaining keywords are those of Breeze's `AerosolActivation` (activation timescale, surface
tension fit, activated droplet radius, supersaturation threshold and floors).
"""
function kappa_aerosol_activation(FT=Oceananigans.defaults.FloatType; number_mixing_ratio, mean_radius=92.5e-9, geometric_std=1.5, kappa=0.2,
                                  thermodynamic_constants=ThermodynamicConstants(FT),
                                  activation_timescale=1,
                                  surface_tension_reference=0.0761,
                                  surface_tension_temperature_derivative=-1.55e-4,
                                  lognormal_activation_factor=4.242,
                                  activated_droplet_radius=1e-6,
                                  activation_supersaturation_threshold=1e-6,
                                  minimum_supersaturation=1e-20,
                                  minimum_saturation_mass_fraction=1e-20)
    mode = KappaAerosolMode(FT; number_mixing_ratio, mean_radius, geometric_std, kappa)
    modes = (mode,)
    c = thermodynamic_constants
    return AerosolActivation{FT, true, typeof(modes)}(
        modes, FT(c.vapor.molar_mass), FT(c.molar_gas_constant), FT(activation_timescale), FT(c.liquid.density),
        FT(surface_tension_reference), FT(surface_tension_temperature_derivative), FT(c.energy_reference_temperature),
        FT(lognormal_activation_factor), FT(activated_droplet_radius), FT(activation_supersaturation_threshold),
        FT(minimum_supersaturation), FT(minimum_saturation_mass_fraction))
end

"""
    kappa_critical_supersaturation(mode::KappaAerosolMode, aerosol, T)

κ-Köhler critical supersaturation of the geometric-mean particle,
`s_m = √(4A³ / (27 κ r_m³))` with the Kelvin parameter `A = 2 M_w σ(T) / (ρ_w R T)` and the
linear surface-tension fit of `aerosol`. For σ_v, M_w, ρ_w of water and κ = β_act this is
identical to Morrison & Grabowski's `2 β⁻¹ᐟ² (A / 3 r_m)³ᐟ²`.
"""
@inline function kappa_critical_supersaturation(mode::KappaAerosolMode, aerosol, T)
    σ = aerosol.surface_tension_reference +
        aerosol.surface_tension_temperature_derivative * (T - aerosol.surface_tension_reference_temperature)
    A = 2 * aerosol.molecular_weight_water * σ / (aerosol.liquid_water_density * aerosol.universal_gas_constant * T)
    x = A / (3 * mode.mean_radius)
    return 2 / sqrt(mode.kappa) * sqrt(x)^3
end

"""
    activated_fraction(mode::KappaAerosolMode, aerosol, T, S)

Fraction of the lognormal mode activated at supersaturation `S` (M&G2007 erf form):
`½ [1 − erf(2 ln(s_m / S) / (4.242 ln σ_g))]`.
"""
@inline function activated_fraction(mode::KappaAerosolMode, aerosol, T, S)
    FT = typeof(S)
    s_m = kappa_critical_supersaturation(mode, aerosol, T)
    S_safe = max(S, aerosol.minimum_supersaturation)
    u = 2 * log(s_m / S_safe) / (aerosol.lognormal_activation_factor * log(mode.geometric_std))
    return FT(0.5) * (1 - PPP.erf(u))
end

# Breeze's per-mode activated number (used by its fixed-reservoir paths and diagnostics):
# the seed number times the κ-Köhler fraction.
@inline PPP.activated_number(mode::KappaAerosolMode, aerosol::AerosolActivation, T, S) =
    mode.number_mixing_ratio * activated_fraction(mode, aerosol, T, S)

# Prognostic reservoir: the equilibrium target is the activated fraction of the particles
# actually present in the cell (interstitial reservoir + those already in droplets).
@inline function PPP.aerosol_activation_rate(aerosol::KappaAerosolActivation{FT, true}, nᶜˡ, nᵃ, qᵛ, qᵛ⁺ˡ, T) where FT
    S = (qᵛ - qᵛ⁺ˡ) / max(qᵛ⁺ˡ, aerosol.minimum_saturation_mass_fraction)
    nᵃ_available = max(zero(nᵃ), nᵃ)
    total = nᶜˡ + nᵃ_available
    f = zero(S)
    for mode in aerosol.modes
        f += activated_fraction(mode, aerosol, T, S)
    end
    f = min(f, one(f))
    N_target = min(f * total, total)
    ncnuc = max(zero(N_target), N_target - nᶜˡ) / aerosol.activation_timescale
    seed_mass = 4 * FT(π) / 3 * aerosol.liquid_water_density * aerosol.activated_droplet_radius^3
    qcnuc = ncnuc * seed_mass
    is_supersaturated = S > aerosol.activation_supersaturation_threshold
    ncnuc = ifelse(is_supersaturated, ncnuc, zero(FT))
    qcnuc = ifelse(is_supersaturated, qcnuc, zero(FT))
    return (; ncnuc, qcnuc)
end

#####
##### Surface aerosol number source
#####

"""
    SurfaceAerosolSource(flux)

Specific tendency of `nᵃ` representing a prescribed surface number flux `flux` [m⁻² s⁻¹]
(SEA STARR: 70 cm⁻² s⁻¹ = 7×10⁵ m⁻² s⁻¹) deposited in the bottom cell: `flux / (ρ₁ Δz₁)`
[kg⁻¹ s⁻¹] at `k = 1`, zero elsewhere. Supply under the `nᵃ` key (Breeze multiplies by ρ).
"""
struct SurfaceAerosolSource{FT, D}
    flux :: FT
    density :: D
end

SurfaceAerosolSource(flux) = SurfaceAerosolSource(flux, nothing)

Adapt.adapt_structure(to, f::SurfaceAerosolSource) = SurfaceAerosolSource(f.flux, adapt(to, f.density))
Base.summary(f::SurfaceAerosolSource) = string("SurfaceAerosolSource(", f.flux, " m⁻² s⁻¹)")
Base.show(io::IO, f::SurfaceAerosolSource) = print(io, summary(f))

@inline function (f::SurfaceAerosolSource)(i, j, k, grid, clock, fields)
    ρ₁ = @inbounds f.density[i, j, 1]
    source = f.flux / (ρ₁ * Δzᶜᶜᶜ(i, j, 1, grid))
    return ifelse(k == 1, source, zero(source))
end

function AtmosphereModels.materialize_atmosphere_model_forcing(f::SurfaceAerosolSource,
                                                               field, name, model_field_names, context::NamedTuple)
    name === :nᵃ || throw(ArgumentError("SurfaceAerosolSource must be supplied under the `nᵃ` key, got $name"))
    FT = eltype(field.grid)
    return SurfaceAerosolSource(convert(FT, f.flux), context.total_density)
end

#####
##### Regeneration on cloud-droplet evaporation
#####

# Cloud-water evaporation rate [kg/kg/s] (≥ 0) that P3 would apply in this cell, from one
# evaluation of its process rates with the same per-cell setup as `_p3_add_tendencies_kernel!`.
@inline function p3_cloud_evaporation_rate(i, j, k, grid, μ, formulation, dynamics, constants, p3, ρ_field, velocities)
    @inbounds begin
        ρ = ρ_field[i, j, k]
        qᵛᵉ = μ.qᵛ[i, j, k]
        ℳ = AtmosphereModels.grid_microphysical_state(i, j, k, grid, p3, μ, ρ, nothing, velocities)
        q = AtmosphereModels.moisture_fractions(p3, ℳ, qᵛᵉ)
        𝒰₀ = AtmosphereModels.diagnose_thermodynamic_state(i, j, k, grid, formulation, dynamics, q)
        𝒰 = AtmosphereModels.maybe_adjust_thermodynamic_state(𝒰₀, p3, qᵛᵉ, constants)
        temperature_tendency = p3_adiabatic_temperature_tendency(ℳ, 𝒰, constants)
        vapor_tendency = zero(temperature_tendency)
        surface_temperature = μ.surface_temperature[i, j, 1]
    end
    properties = p3_process_properties(p3, ρ, ℳ)
    rates = compute_p3_process_rates(p3, ρ, ℳ, 𝒰, constants, properties,
                                     surface_temperature, temperature_tendency, vapor_tendency)
    return max(zero(ρ), -rates.condensation)
end

"""
    cloud_evaporation_rate_field(model)

`Field` of P3's cloud-water evaporation rate [kg/kg/s] (the negative part of its
condensation rate), evaluated from the current model state when `compute!`d.
"""
function cloud_evaporation_rate_field(model)
    grid = model.grid
    p3 = model.microphysics
    is_p3(p3) || throw(ArgumentError("cloud_evaporation_rate_field needs P3 microphysics"))
    ρ = AtmosphereModels.total_density(model.dynamics)
    op = KernelFunctionOperation{Center, Center, Center}(p3_cloud_evaporation_rate, grid,
                                                        model.microphysical_fields, model.formulation, model.dynamics,
                                                        model.thermodynamic_constants, p3, ρ, model.velocities)
    return Field(op)
end

"""
    EvaporationRegeneration(evaporation_rate; minimum_mass_mixing_ratio=1e-14)

Forcing pair returning aerosol to the reservoir as droplets evaporate. With `E` the
cloud-water evaporation rate [kg/kg/s] (a `Field` refreshed by
[`EvaporationRateUpdater`](@ref) once per time step) the droplet number is removed in
proportion to the evaporated mass and the same number is regenerated:

    R = nᶜˡ E / max(qᶜˡ, q_min)      F_nᶜˡ = -R,   F_nᵃ = +R

Supply `regeneration.nᶜˡ` under `nᶜˡ` and `regeneration.nᵃ` under `nᵃ` (specific keys).
**Labelled departure:** P3 itself leaves `nᶜˡ` unchanged while cloud water evaporates and
discards the number when the mass is exhausted; the proportional rule is the usual two-moment
bulk convention and makes the aerosol number budget closed under evaporation.
"""
struct EvaporationRegeneration{S, E, FT}
    sign :: S
    evaporation_rate :: E
    minimum_mass_mixing_ratio :: FT
end

EvaporationRegeneration(evaporation_rate; minimum_mass_mixing_ratio=1e-14) =
    (; nᶜˡ = EvaporationRegeneration(Val(-1), evaporation_rate, minimum_mass_mixing_ratio),
       nᵃ = EvaporationRegeneration(Val(+1), evaporation_rate, minimum_mass_mixing_ratio))

Adapt.adapt_structure(to, f::EvaporationRegeneration) =
    EvaporationRegeneration(f.sign, adapt(to, f.evaporation_rate), f.minimum_mass_mixing_ratio)
Base.summary(f::EvaporationRegeneration{Val{s}}) where s = string("EvaporationRegeneration(", s > 0 ? "+" : "-", "nᶜˡ E/qᶜˡ)")
Base.show(io::IO, f::EvaporationRegeneration) = print(io, summary(f))

@inline function (f::EvaporationRegeneration{Val{s}})(i, j, k, grid, clock, fields) where s
    @inbounds begin
        E = f.evaporation_rate[i, j, k]
        nᶜˡ = fields.nᶜˡ[i, j, k]
        qᶜˡ = fields.qᶜˡ[i, j, k]
    end
    R = max(zero(nᶜˡ), nᶜˡ) * E / max(qᶜˡ, f.minimum_mass_mixing_ratio)
    return s * R
end

function AtmosphereModels.materialize_atmosphere_model_forcing(f::EvaporationRegeneration{Val{s}},
                                                               field, name, model_field_names, context::NamedTuple) where s
    expected = s > 0 ? :nᵃ : :nᶜˡ
    name === expected || throw(ArgumentError("this EvaporationRegeneration must be supplied under `$expected`, got $name"))
    FT = eltype(field.grid)
    return EvaporationRegeneration(f.sign, f.evaporation_rate, convert(FT, f.minimum_mass_mixing_ratio))
end

"""
    EvaporationRateUpdater(rate_field, evaporation_rate)

Callback (`IterationInterval(1)`) that `compute!`s `rate_field` (from
[`cloud_evaporation_rate_field`](@ref)) and copies it into the plain `CenterField`
`evaporation_rate` held by the [`EvaporationRegeneration`](@ref) forcings before each time
step. The regeneration therefore uses the rate of the state at the start of the step
(one-step lag; the same approximation as ENA's once-per-step projection).
"""
struct EvaporationRateUpdater{R, E}
    rate_field :: R
    evaporation_rate :: E
end

function (u::EvaporationRateUpdater)(simulation)
    compute!(u.rate_field)
    parent(u.evaporation_rate) .= parent(u.rate_field)
    return nothing
end

#####
##### Aerosol number budget diagnostics
#####

"""
    aerosol_number_columns(model)

Vertically integrated number densities [m⁻²] of the interstitial aerosol (`∫ρnᵃ dz`),
cloud droplets and rain drops, plus their sum (the total aerosol-containing particle
number), as 2D `Field`s for the budget check
`d/dt ∫(nᵃ + nᶜˡ + nʳ) = surface source − coalescence − surface rain number flux + nudging + subsidence`.
"""
function aerosol_number_columns(model)
    μ = model.microphysical_fields
    haskey(μ, :ρnᵃ) || throw(ArgumentError("aerosol_number_columns needs a prognostic aerosol reservoir"))
    nᵃ = Field(Integral(μ.ρnᵃ, dims=3))
    nᶜˡ = Field(Integral(μ.ρnᶜˡ, dims=3))
    nʳ = Field(Integral(μ.ρnʳ, dims=3))
    total = Field(Integral(μ.ρnᵃ + μ.ρnᶜˡ + μ.ρnʳ, dims=3))
    return (; nᵃ, nᶜˡ, nʳ, total)
end
