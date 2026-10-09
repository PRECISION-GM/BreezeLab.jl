# Tier-1 prescribed aerosol inside a Breeze P3 model (CPU, small grid): fields, prognostic set,
# height scaling of the activation rate, and a few time steps of a supersaturated column.
using BreezeLab: tracer_mip_aerosol_profile, total_number_mixing_ratio, surface_number_mixing_ratios, PrescribedAerosolProfile
using Breeze.Microphysics.PredictedParticleProperties: CloudDroplets, AerosolActivation, compute_cloud_droplet_activation,
                                                       has_prognostic_aerosol, aerosol_activation_rate
using Breeze.Thermodynamics: saturation_specific_humidity, PlanarLiquidSurface
using Oceananigans.TimeSteppers: time_step!

@testset "Tier-1 prescribed aerosol in P3" begin
    FT = Float64
    profile = tracer_mip_aerosol_profile(:aug07; FT)
    aerosol = with_float_type(() -> PrescribedAerosolProfile(profile), FT)
    @test !has_prognostic_aerosol(aerosol)
    @test aerosol.surface_number ≈ sum(surface_number_mixing_ratios(profile))
    @test length(aerosol.activation.modes) == 2

    grid = test_grid(; Nz=40, Lz=8000)
    constants = ThermodynamicConstants(FT)
    reference_state = ReferenceState(grid, constants; base_pressure=101300, potential_temperature=z -> 295 + 0.003z)
    dynamics = AnelasticDynamics(reference_state)
    p3 = P3Microphysics(FT; cloud=CloudDroplets(FT; number_concentration=100e6), aerosol)
    names = Breeze.AtmosphereModels.prognostic_field_names(p3)
    @test :ρnᶜˡ ∈ names && :ρnᵃ ∉ names
    scalar_advection = BreezeLab.scalar_advection_schemes(5, p3, :qᵛ)
    model = AtmosphereModel(grid; dynamics, microphysics=p3, thermodynamic_constants=constants,
                            momentum_advection=WENO(order=5), scalar_advection)
    μ = model.microphysical_fields
    @test haskey(μ, :nᵃ) && haskey(μ, :ρnᶜˡ) && haskey(μ, :nᶜˡ) && !haskey(μ, :ρnᵃ)

    @testset "prescribed number field follows the profile in height AGL" begin
        zc = Array(znodes(grid, Center()))
        nᵃ = Array(interior(μ.nᵃ, 1, 1, :))
        @test all(nᵃ .≈ [total_number_mixing_ratio(profile, z) for z in zc])
        @test nᵃ[1] < aerosol.surface_number                       # first cell center is above the surface
        @test nᵃ[end] ≈ 50e6                                        # floor aloft
        @test interior(μ.nᵃ)[3, 5, 7] == interior(μ.nᵃ)[1, 1, 7]    # horizontally homogeneous
    end

    @testset "activation scales with the local prescribed number" begin
        T = 290.0; ρ = 1.1
        qᵛ⁺ˡ = saturation_specific_humidity(T, ρ, constants, PlanarLiquidSurface())
        qᵛ = 1.004 * qᵛ⁺ˡ                                             # 0.4 % supersaturation
        surface = compute_cloud_droplet_activation(aerosol, p3, 0.0, 0.0, aerosol.surface_number, qᵛ, qᵛ⁺ˡ, T, ρ, constants)
        reference = aerosol_activation_rate(aerosol.activation, 0.0, qᵛ, qᵛ⁺ˡ, T)
        @test surface.number ≈ reference.ncnuc && surface.mass ≈ reference.qcnuc   # surface: Breeze's own rate
        @test surface.number > 0
        for z in (1000.0, 2500.0, 7000.0)
            n = total_number_mixing_ratio(profile, z)
            aloft = compute_cloud_droplet_activation(aerosol, p3, 0.0, 0.0, n, qᵛ, qᵛ⁺ˡ, T, ρ, constants)
            @test aloft.number ≈ surface.number * n / aerosol.surface_number
            @test aloft.number ≤ n / aerosol.activation.activation_timescale   # never more droplets than aerosols
        end
        # Existing droplets reduce the nucleation; subsaturation switches it off
        partial = compute_cloud_droplet_activation(aerosol, p3, 1e-4, 0.5 * surface.number * aerosol.activation.activation_timescale,
                                                   aerosol.surface_number, qᵛ, qᵛ⁺ˡ, T, ρ, constants)
        @test 0 < partial.number < surface.number
        dry = compute_cloud_droplet_activation(aerosol, p3, 0.0, 0.0, aerosol.surface_number, 0.9 * qᵛ⁺ˡ, qᵛ⁺ˡ, T, ρ, constants)
        @test dry.number == 0 && dry.mass == 0
    end

    @testset "time stepping a supersaturated layer activates droplets, nᵃ stays fixed" begin
        zc = Array(znodes(grid, Center()))
        ρᵣ = Array(interior(reference_state.density, 1, 1, :))
        T_profile = [295.0 - 0.0065z for z in zc]
        qᵛ₀ = [(1000 ≤ z ≤ 2000 ? 1.003 : 0.8) * saturation_specific_humidity(T_profile[k], ρᵣ[k], constants, PlanarLiquidSurface())
               for (k, z) in enumerate(zc)]
        col(v) = repeat(reshape(v, 1, 1, length(v)), size(grid, 1), size(grid, 2), 1)
        set!(model; T=col(T_profile), qᵛ=col(qᵛ₀))
        nᵃ_before = copy(interior(μ.nᵃ))
        for _ in 1:5
            time_step!(model, 1.0)
        end
        @test all(isfinite, interior(model.temperature))
        @test interior(μ.nᵃ) == nᵃ_before                            # fixed at each model step
        nᶜˡ = Array(interior(μ.nᶜˡ, 1, 1, :))
        k_cloud = findall(z -> 1000 ≤ z ≤ 2000, zc)
        @test maximum(nᶜˡ[k_cloud]) > 0
        @test maximum(nᶜˡ[k_cloud]) ≤ maximum(Array(interior(μ.nᵃ, 1, 1, :))[k_cloud]) * (1 + 1e-6)
        k_clear = findall(z -> z > 4000, zc)
        @test all(nᶜˡ[k_clear] .≤ 1e-6 * maximum(nᶜˡ[k_cloud]))
    end
end

@testset "TRACER-MIP process-rate accumulators (P3)" begin
    using BreezeLab: ProcessRateAccumulators, accumulate_process_rates!, reset_process_accumulators!,
                     process_rate_output_callback, PROCESS_RATE_NAMES, PROCESS_RATE_UNITS
    using Breeze.Microphysics.PredictedParticleProperties: CloudDroplets
    using Breeze.Thermodynamics: saturation_specific_humidity, PlanarLiquidSurface
    FT = Float64
    grid = test_grid(; Nz=40, Lz=8000)
    constants = ThermodynamicConstants(FT)
    reference_state = ReferenceState(grid, constants; base_pressure=101300, potential_temperature=z -> 295 + 0.003z)
    dynamics = AnelasticDynamics(reference_state)
    aerosol = with_float_type(() -> PrescribedAerosolProfile(tracer_mip_aerosol_profile(:aug07; FT)), FT)
    p3 = P3Microphysics(FT; cloud=CloudDroplets(FT; number_concentration=100e6), aerosol)
    model = AtmosphereModel(grid; dynamics, microphysics=p3, thermodynamic_constants=constants,
                            momentum_advection=WENO(order=5), scalar_advection=BreezeLab.scalar_advection_schemes(5, p3, :qᵛ))
    zc = Array(znodes(grid, Center()))
    ρᵣ = Array(interior(reference_state.density, 1, 1, :))
    T_profile = [295.0 - 0.0065z for z in zc]
    qᵛ₀ = [(1000 ≤ z ≤ 2000 ? 1.003 : 0.8) * saturation_specific_humidity(T_profile[k], ρᵣ[k], constants, PlanarLiquidSurface())
           for (k, z) in enumerate(zc)]
    col(v) = repeat(reshape(v, 1, 1, length(v)), size(grid, 1), size(grid, 2), 1)
    set!(model; T=col(T_profile), qᵛ=col(qᵛ₀))

    acc = ProcessRateAccumulators(model)
    @test keys(acc.fields) == PROCESS_RATE_NAMES && keys(PROCESS_RATE_UNITS) == PROCESS_RATE_NAMES
    @test all(f -> all(iszero, interior(f)), acc.fields)
    simulation = Simulation(model; Δt=1.0, stop_iteration=3)
    add_callback!(simulation, acc, IterationInterval(1))
    run!(simulation)
    k_cloud = findall(z -> 1000 ≤ z ≤ 2000, zc)
    k_clear = findall(z -> z > 4000, zc)
    cond = Array(interior(acc.fields.liquid_condensation, 1, 1, :))
    nuc = Array(interior(acc.fields.droplet_nucleation, 1, 1, :))
    heat = Array(interior(acc.fields.latent_heating, 1, 1, :))
    @test all(f -> all(isfinite, interior(f)), acc.fields)
    @test maximum(cond[k_cloud]) > 0 && all(cond[k_clear] .== 0)        # condensation only in the supersaturated layer
    @test maximum(nuc[k_cloud]) > 0                                       # droplets nucleated (#/mg per interval)
    @test maximum(heat[k_cloud]) > 0                                      # condensational heating, K per interval
    k_warm = findall(z -> 295.0 - 0.0065z > 273.15, zc)                    # ice processes only below 0 °C
    @test all(Array(interior(acc.fields.ice_deposition, 1, 1, :))[k_warm] .== 0)
    @test all(Array(interior(acc.fields.melting, 1, 1, :))[k_warm] .== 0)
    @test all(interior(acc.fields.ice_deposition) .>= 0) && all(interior(acc.fields.ice_sublimation) .>= 0)
    @test all(f -> all(≥(0), interior(f)), (acc.fields.liquid_condensation, acc.fields.liquid_evaporation,
                                            acc.fields.cloud_riming, acc.fields.autoconversion_accretion))
    # Explicit single accumulation adds exactly one rate × Δt on top of the current state
    before = copy(interior(acc.fields.liquid_condensation))
    accumulate_process_rates!(acc, 2.0)
    @test all(interior(acc.fields.liquid_condensation) .≥ before)
    reset_process_accumulators!(acc)
    @test all(f -> all(iszero, interior(f)), acc.fields)
    cb = process_rate_output_callback(acc, IterationInterval(2))
    @test cb isa Oceananigans.Simulations.Callback

    # Writer and reset on the same schedule: the saved interval must be the completed, non-zero one
    model2 = AtmosphereModel(grid; dynamics, microphysics=p3, thermodynamic_constants=constants,
                             momentum_advection=WENO(order=5), scalar_advection=BreezeLab.scalar_advection_schemes(5, p3, :qᵛ))
    set!(model2; T=col(T_profile), qᵛ=col(qᵛ₀))
    acc2 = ProcessRateAccumulators(model2)
    sim2 = Simulation(model2; Δt=1.0, stop_iteration=4)
    add_callback!(sim2, acc2, IterationInterval(1))
    tmp = mktempdir()
    sim2.output_writers[:acc] = JLD2Writer(model2, acc2.fields; schedule=IterationInterval(2), filename=joinpath(tmp, "acc.jld2"), overwrite_files=true)
    add_callback!(sim2, process_rate_output_callback(acc2, IterationInterval(2)))
    run!(sim2)
    saved = FieldTimeSeries(joinpath(tmp, "acc.jld2"), "liquid_condensation")
    @test length(saved.times) == 3                                   # iterations 0, 2, 4
    @test maximum(interior(saved[1])) == 0                           # nothing accumulated before the first step
    @test maximum(interior(saved[2])) > 0 && maximum(interior(saved[3])) > 0   # each hourly-like interval saved before its reset
    rm(tmp; recursive=true)
end
