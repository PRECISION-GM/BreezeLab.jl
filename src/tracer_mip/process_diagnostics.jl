#####
##### TRACER-MIP microphysical process-rate diagnostics for Breeze P3.
#####
##### Table 2 of the roadmap asks for separate process rates accumulated between output writes
##### ("kg/kg per output time interval": the sum of rate × Δt over the interval). Breeze's P3
##### adds only net tendencies to the prognostics, so the accumulators re-evaluate
##### `compute_p3_process_rates` from the current state once per time step (the same state and
##### driver the tendency kernel uses) and add rate × Δt. Mapping of P3 rates to the requested
##### quantities (all kg/kg per interval unless noted; number rates in #/mg per interval):
#####   liquid_condensation  = max(condensation, 0) + ccn_activation_mass + rain_condensation + coating_condensation
#####   liquid_evaporation   = max(-condensation, 0) + rain_evaporation + coating_evaporation
#####   ice_deposition       = max(deposition, 0) + nucleation_mass      (deposition nucleation is vapor → ice)
#####   ice_sublimation      = max(-deposition, 0)
#####   melting              = partial_melting + complete_melting        (all melting sources)
#####   freezing             = cloud + rain immersion freezing + cloud + rain homogeneous freezing
#####                          (in-situ freezing of cloud and rain; riming excluded; refreezing of the
#####                          liquid coating on ice, which is not a droplet/raindrop freezing, excluded)
#####   droplet_nucleation   = ccn_activation_number                     [#/mg per interval]
#####   ice_nucleation       = deposition nucleation + immersion + homogeneous freezing numbers; rime
#####                          splintering (secondary) excluded          [#/mg per interval]
#####   cloud_riming         = cloud_riming + wet_growth_cloud + cloud_warm_collection
#####   rain_riming          = rain_riming + wet_growth_rain + rain_warm_collection
#####   autoconversion_accretion = autoconversion + accretion
#####   latent_heating       = net (ℒˡ·liquid net + ℒⁱ·ice net + ℒᶠ·(freezing − melting)) / cᵖᵐ, accumulated
#####                          as K per interval; divide by the interval for K/s.
##### The accumulators are reset by `reset_process_accumulators!` (call it after each output write
##### through the writer's schedule, see `process_rate_output_callback`). On restart the partial
##### interval that was in flight is lost unless the accumulator fields are checkpointed.
#####

using KernelAbstractions: @kernel, @index
using Oceananigans: CenterField, Center
using Oceananigans.Utils: launch!
using Oceananigans.Fields: interior
using Oceananigans.Simulations: Callback
using Breeze: ThermodynamicConstants
using Breeze.AtmosphereModels: AtmosphereModels
using Breeze.Thermodynamics: liquid_latent_heat, ice_latent_heat, mixture_heat_capacity, temperature
using Breeze.Microphysics.PredictedParticleProperties: P3MicrophysicalState, compute_p3_process_rates,
                                                       p3_process_properties, p3_adiabatic_temperature_tendency

const PROCESS_RATE_NAMES = (:liquid_condensation, :liquid_evaporation, :ice_deposition, :ice_sublimation,
                            :melting, :freezing, :droplet_nucleation, :ice_nucleation,
                            :cloud_riming, :rain_riming, :autoconversion_accretion, :latent_heating)

const PROCESS_RATE_UNITS = (liquid_condensation = "kg/kg per output interval", liquid_evaporation = "kg/kg per output interval",
                            ice_deposition = "kg/kg per output interval", ice_sublimation = "kg/kg per output interval",
                            melting = "kg/kg per output interval", freezing = "kg/kg per output interval",
                            droplet_nucleation = "#/mg per output interval", ice_nucleation = "#/mg per output interval",
                            cloud_riming = "kg/kg per output interval", rain_riming = "kg/kg per output interval",
                            autoconversion_accretion = "kg/kg per output interval", latent_heating = "K per output interval")

struct ProcessRateAccumulators{F, M, C}
    fields :: F           # NamedTuple of CenterFields keyed by PROCESS_RATE_NAMES
    model :: M
    constants :: C
    pending_reset :: Base.RefValue{Bool}   # set by the output-schedule callback, applied before the next accumulation
end

"""
    ProcessRateAccumulators(model)

Allocate the twelve TRACER-MIP process-rate accumulators for a Breeze `AtmosphereModel` using P3.
"""
function ProcessRateAccumulators(model)
    fields = NamedTuple{PROCESS_RATE_NAMES}(ntuple(_ -> CenterField(model.grid), length(PROCESS_RATE_NAMES)))
    return ProcessRateAccumulators(fields, model, model.thermodynamic_constants, Ref(false))
end

@inline function mip_process_rates(p3, ρ, ℳ, 𝒰, constants, surface_temperature)
    properties = p3_process_properties(p3, ρ, ℳ)
    temperature_tendency = p3_adiabatic_temperature_tendency(ℳ, 𝒰, constants)
    r = compute_p3_process_rates(p3, ρ, ℳ, 𝒰, constants, properties, surface_temperature,
                                 temperature_tendency, zero(temperature_tendency))
    FT = typeof(ρ)
    cond = r.condensation
    dep = r.deposition
    liquid_condensation = max(cond, zero(FT)) + r.ccn_activation_mass + r.rain_condensation + r.coating_condensation
    liquid_evaporation = max(-cond, zero(FT)) + r.rain_evaporation + r.coating_evaporation
    ice_deposition = max(dep, zero(FT)) + r.nucleation_mass
    ice_sublimation = max(-dep, zero(FT))
    melting = r.partial_melting + r.complete_melting
    freezing = r.cloud_freezing_mass + r.rain_freezing_mass + r.cloud_homogeneous_mass + r.rain_homogeneous_mass
    droplet_nucleation = r.ccn_activation_number
    ice_nucleation = r.nucleation_number + r.cloud_freezing_number + r.rain_freezing_number +
                     r.cloud_homogeneous_number + r.rain_homogeneous_number
    cloud_riming = r.cloud_riming + r.wet_growth_cloud + r.cloud_warm_collection
    rain_riming = r.rain_riming + r.wet_growth_rain + r.rain_warm_collection
    autoconversion_accretion = r.autoconversion + r.accretion
    T = temperature(𝒰, constants)
    ℒˡ = liquid_latent_heat(T, constants)
    ℒⁱ = ice_latent_heat(T, constants)
    ℒᶠ = ℒⁱ - ℒˡ
    cᵖᵐ = mixture_heat_capacity(𝒰.moisture_mass_fractions, constants)
    riming = cloud_riming + rain_riming + r.refreezing
    latent_heating = (ℒˡ * (liquid_condensation - liquid_evaporation) + ℒⁱ * (ice_deposition - ice_sublimation) +
                      ℒᶠ * (freezing + riming - melting)) / cᵖᵐ
    return (; liquid_condensation, liquid_evaporation, ice_deposition, ice_sublimation, melting, freezing,
              droplet_nucleation, ice_nucleation, cloud_riming, rain_riming, autoconversion_accretion, latent_heating)
end

@kernel function _accumulate_process_rates!(acc, Δt, μ, formulation, dynamics, grid, constants, p3, ρ_field, velocities)
    i, j, k = @index(Global, NTuple)
    @inbounds begin
        ρ = ρ_field[i, j, k]
        qᵛᵉ = μ.qᵛ[i, j, k]
        ℳ = AtmosphereModels.grid_microphysical_state(i, j, k, grid, p3, μ, ρ, nothing, velocities)
        q = AtmosphereModels.moisture_fractions(p3, ℳ, qᵛᵉ)
        𝒰₀ = AtmosphereModels.diagnose_thermodynamic_state(i, j, k, grid, formulation, dynamics, q)
        𝒰 = AtmosphereModels.maybe_adjust_thermodynamic_state(𝒰₀, p3, qᵛᵉ, constants)
        surface_temperature = μ.surface_temperature[i, j, 1]
    end
    rates = mip_process_rates(p3, ρ, ℳ, 𝒰, constants, surface_temperature)
    FT = typeof(ρ)
    per_mg = FT(1e-6)
    @inbounds begin
        acc.liquid_condensation[i, j, k]      += Δt * rates.liquid_condensation
        acc.liquid_evaporation[i, j, k]       += Δt * rates.liquid_evaporation
        acc.ice_deposition[i, j, k]           += Δt * rates.ice_deposition
        acc.ice_sublimation[i, j, k]          += Δt * rates.ice_sublimation
        acc.melting[i, j, k]                  += Δt * rates.melting
        acc.freezing[i, j, k]                 += Δt * rates.freezing
        acc.droplet_nucleation[i, j, k]       += Δt * rates.droplet_nucleation * per_mg
        acc.ice_nucleation[i, j, k]           += Δt * rates.ice_nucleation * per_mg
        acc.cloud_riming[i, j, k]             += Δt * rates.cloud_riming
        acc.rain_riming[i, j, k]              += Δt * rates.rain_riming
        acc.autoconversion_accretion[i, j, k] += Δt * rates.autoconversion_accretion
        acc.latent_heating[i, j, k]           += Δt * rates.latent_heating
    end
end

"""
    accumulate_process_rates!(acc::ProcessRateAccumulators, Δt)

Add `rate × Δt` of every process to the accumulators from the model's current state.
"""
function accumulate_process_rates!(acc::ProcessRateAccumulators, Δt)
    model = acc.model
    grid = model.grid
    μ = model.microphysical_fields
    ρ_field = AtmosphereModels.total_density(model.dynamics)
    launch!(grid.architecture, grid, :xyz, _accumulate_process_rates!,
            acc.fields, convert(eltype(grid), Δt), μ, model.formulation, model.dynamics, grid,
            acc.constants, model.microphysics, ρ_field, model.velocities)
    return nothing
end

# Callback form: accumulate over the step that just completed (the simulation's Δt). A reset requested
# by `process_rate_output_callback` on the previous iteration is applied first, so the writer that ran
# after the callbacks on that iteration saw the complete interval.
function (acc::ProcessRateAccumulators)(simulation)
    if acc.pending_reset[]
        reset_process_accumulators!(acc)
        acc.pending_reset[] = false
    end
    return accumulate_process_rates!(acc, simulation.Δt)
end

"""
    reset_process_accumulators!(acc)

Zero every accumulator (to be called right after each output write).
"""
function reset_process_accumulators!(acc::ProcessRateAccumulators)
    for f in acc.fields
        fill!(f, 0)
    end
    return nothing
end

"""
    process_rate_output_callback(acc, schedule)

A callback on the output `schedule` that *requests* a reset of the accumulators: Oceananigans runs
callbacks before output writers within a step, so the accumulators are zeroed only at the start of
the next accumulation, after the writer on the same schedule has saved the completed interval.
"""
process_rate_output_callback(acc::ProcessRateAccumulators, schedule) = Callback(sim -> (acc.pending_reset[] = true), schedule)
