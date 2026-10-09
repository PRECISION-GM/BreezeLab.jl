# SEA STARR reference-value tests: DEPHY reader, units, interpolation, nudging mask,
# κ-Köhler activation, surface source and subsidence applied exactly once.
using Test
using BreezeLab
using Breeze
using Oceananigans
using Oceananigans.Units
using Oceananigans.Units: Time
using Oceananigans.Fields: interior
using Statistics
using Dates: DateTime
using NCDatasets
using Breeze.Microphysics.PredictedParticleProperties: aerosol_activation_rate, AerosolActivation

const SEASTARR_DIRS = (joinpath(@__DIR__, "..", "data", "seastarr_22241697"),
                       "/shared/home/greg/breezelab-runs/20261006/inputs/mip_sources/seastarr_22241697")
const SEASTARR_DIR = something(findfirst(d -> isfile(joinpath(d, "SEA_STARR_CTRL_SCM_driver.nc")), SEASTARR_DIRS), 0) == 0 ?
                     nothing : SEASTARR_DIRS[findfirst(d -> isfile(joinpath(d, "SEA_STARR_CTRL_SCM_driver.nc")), SEASTARR_DIRS)]
const HAVE_SEASTARR = !isnothing(SEASTARR_DIR)

# A tiny synthetic DEPHY file (top-down levels, mg⁻¹ aerosol) for reader tests that must not
# depend on the staged inputs.
function write_synthetic_dephy(path)
    zh = [3000.0, 2000.0, 1000.0, 500.0, 100.0, 5.0]          # top-down, as the SEA STARR files
    nt = 3
    NCDataset(path, "c") do ds
        ds.attrib["format_version"] = "DEPHY SCM format version 1"
        ds.attrib["case"] = "SYNTHETIC/CTRL"
        ds.attrib["end_date"] = "2017-08-15 23:00:00"
        defDim(ds, "time", nt); defDim(ds, "t0", 1); defDim(ds, "zh", length(zh)); defDim(ds, "lat", 1); defDim(ds, "lon", 1)
        t = defVar(ds, "time", Float64, ("time",); attrib=Dict("units" => "seconds since 2017-08-15 21:00:00", "calendar" => "gregorian"))
        t[:] = [0.0, 3600.0, 7200.0]
        t0 = defVar(ds, "t0", Float64, ("t0",); attrib=Dict("units" => "seconds since 2017-08-15 21:00:00", "calendar" => "gregorian"))
        t0[:] = [0.0]
        defVar(ds, "zh", Float64, ("zh",); attrib=Dict("units" => "m"))[:] = zh
        defVar(ds, "lat", Float64, ("lat",))[:] = [-11.0]
        defVar(ds, "lon", Float64, ("lon",))[:] = [-8.0]
        defVar(ds, "ps", Float64, ("t0",))[:] = [101000.0]
        defVar(ds, "ts", Float64, ("t0",))[:] = [293.0]
        defVar(ds, "z0", Float64, ("t0",))[:] = [1e-4]
        prof(name, values; units="") = (v = defVar(ds, name, Float64, ("zh", "t0"); attrib=Dict("units" => units)); v[:, 1] = values)
        prof("thetal", [310.0, 305.0, 300.0, 290.0, 290.0, 290.0])
        prof("ta", [285.0, 287.0, 289.0, 287.0, 290.0, 291.0])
        prof("qt", [0.002, 0.002, 0.003, 0.009, 0.009, 0.009]; units="1")
        prof("na", [1000.0, 900.0, 800.0, 100.0, 100.0, 100.0]; units="mg-1")
        prof("ua", fill(-5.0, 6)); prof("va", fill(3.0, 6)); prof("pa", [70000.0, 80000.0, 90000.0, 95000.0, 99800.0, 100900.0])
        force(name, values; units="") = (v = defVar(ds, name, Float64, ("zh", "time"); attrib=Dict("units" => units)); v[:, :] = values)
        w = repeat([-0.004, -0.004, -0.002, -0.001, -0.0002, -0.00001], 1, nt)
        force("wa", w); force("wap", -w .* 1.0); force("ug", fill(-5.0, 6, nt)); force("vg", fill(3.0, 6, nt))
        force("thetal_nud", repeat([310.0, 305.0, 300.0, 290.0, 290.0, 290.0], 1, nt) .+ [0.0 1.0 2.0])
        force("ta_nud", fill(288.0, 6, nt)); force("qt_nud", fill(0.002, 6, nt); units="1")
        force("na_nud", fill(900.0, 6, nt); units="mg-1"); force("ua_nud", fill(-5.0, 6, nt)); force("va_nud", fill(3.0, 6, nt))
        force("o3", fill(3e-8, 6, nt)); force("pa_force", repeat([70000.0, 80000.0, 90000.0, 95000.0, 99800.0, 100900.0], 1, nt))
        for n in ("ta", "thetal", "qt", "na"); force("nudging_constant_" * n, fill(1 / 1800, 6, nt)); end
        for n in ("ua", "va"); force("nudging_constant_" * n, fill(1 / 10800, 6, nt)); end
        defVar(ds, "ts_force", Float64, ("time",))[:] = [293.0, 293.5, 294.0]
        defVar(ds, "ps_force", Float64, ("time",))[:] = [101000.0, 100990.0, 100980.0]
    end
    return path
end

@testset "SEA STARR" begin
    @testset "protocol vertical grid" begin
        f = sea_starr_vertical_faces()
        @test f[1] == 0 && f[end] == 6500
        @test length(f) == 289                                   # Nz = 288
        @test all(diff(f[1:251]) .≈ 10)                          # 250 cells of 10 m to 2500 m
        @test f[251] == 2500
        d = diff(f[251:end])
        @test all(isapprox.(d[2:end-1] ./ d[1:end-2], 1.1; rtol=1e-6))   # 10 % growth
        @test d[1] ≈ 11                                          # first stretched cell
        @test 360 < d[end] < 375                                  # last cell clipped to 6500 m
        @test issorted(f)
    end

    @testset "DEPHY reader (synthetic)" begin
        path = write_synthetic_dephy(joinpath(mktempdir(), "synthetic_driver.nc"))
        d = read_dephy_driver(path)
        @test d.case == "SYNTHETIC/CTRL"
        @test d.start_time == DateTime(2017, 8, 15, 21) && d.end_time == DateTime(2017, 8, 15, 23)
        @test d.times == [0.0, 3600.0, 7200.0]
        @test issorted(d.z) && d.z[1] == 5.0 && d.z[end] == 3000.0      # reversed to ascending
        @test d.initial.thetal[1] == 290.0 && d.initial.thetal[end] == 310.0
        @test d.initial.na[1] == 1e8 && d.initial.na[end] == 1e9           # mg⁻¹ → kg⁻¹
        @test aerosol_number_per_kg(100.0) == 1e8
        @test d.forcing.w[end, 1] == -0.004 && d.forcing.w[1, 1] == -0.00001
        @test d.forcing.thetal_nud[1, 3] == 292.0                            # time axis kept
        @test d.nudging_rates.thetal == 1 / 1800 && d.nudging_rates.u == 1 / 10800
        @test d.sst == [293.0, 293.5, 294.0]
        @test d.roughness_length == 1e-4 && d.surface_pressure == 101000.0
        @test inversion_height(d.z, d.initial.thetal) == 750.0                # between 500 and 1000 m
        # interpolation onto a model grid: linear in height, constant extrapolation, linear in time
        grid = RectilinearGrid(CPU(), Float64; size=(8, 8, 10), x=(0, 800), y=(0, 800), z=(0, 4000), halo=(5, 5, 5), topology=(Periodic, Periodic, Bounded))
        zc = Array(znodes(grid, Center()))                                   # 200, 600, 1000, ...
        w = driver_profile_time_series(grid, d, :w, zc)
        @test w[1, 1, 1, Time(0.0)] ≈ -0.0002 + (-0.001 + 0.0002) * (200 - 100) / 400
        @test w[1, 1, 3, Time(0.0)] ≈ -0.002                                  # exactly on a level
        @test w[1, 1, 10, Time(0.0)] ≈ -0.004                                 # above the top level: held
        θn = driver_profile_time_series(grid, d, :thetal_nud, zc)
        @test θn[1, 1, 3, Time(1800.0)] ≈ 300.5                               # half-way between records
        @test driver_initial_profile(d, :qt, 750.0) ≈ 0.006
    end

    @testset "nudging mask and masked mean-profile nudging" begin
        z = collect(50.0:100.0:2950.0)
        w = nudging_mask_weights(z, 1300.0; ramp_depth=200)
        @test all(w[z .≤ 1300] .== 0)
        @test all(w[z .≥ 1500] .== 1)
        @test w[findfirst(==(1350.0), z)] ≈ 0.25 && w[findfirst(==(1450.0), z)] ≈ 0.75
        grid = RectilinearGrid(CPU(), Float64; size=(8, 8, 30), x=(0, 800), y=(0, 800), z=(0, 3000), halo=(5, 5, 5), topology=(Periodic, Periodic, Bounded))
        zc = Array(znodes(grid, Center()))
        constants = ThermodynamicConstants(Float64)
        reference_state = ReferenceState(grid, constants; base_pressure=101000, potential_temperature=z -> z < 1000 ? 290.0 : 300.0 + 0.003z)
        dynamics = AnelasticDynamics(reference_state)
        microphysics = SaturationAdjustment(Float64; equilibrium=WarmPhaseEquilibrium())
        times = [0.0, 7200.0]
        target = profile_time_series(grid, times, [fill(305.0, 30), fill(305.0, 30)])
        mask = Field{Nothing, Nothing, Center}(grid)
        nudge = InversionFollowingNudging(target, mask; timescale=1800)
        model = AtmosphereModel(grid; dynamics, microphysics, thermodynamic_constants=constants, forcing=(; θ=nudge))
        # inversion (10 K jump) between 950 and 1050 m in one column shifted up by 200 m in another
        set!(model; θ=(x, y, z) -> (z < (x < 100 ? 1200 : 1000) ? 290.0 : 300.0 + 0.003z), qᵗ=1e-3)
        updater = InversionMaskUpdater(model; mask, offset=100, ramp_depth=200)
        updater(model)
        zi = interior(updater.inversion_height)
        @test zi[1, 1, 1] == 1200.0 && zi[8, 1, 1] == 1000.0              # per-column maximum gradient
        @test updater.inversion_height_max[] == 1200.0 && updater.mask_base[] == 1300.0
        m = vec(interior(mask))
        @test all(m[zc .< 1300] .== 0) && all(m[zc .> 1500] .== 1) && 0 < m[findfirst(==(1350.0), zc)] < 1
        Oceananigans.TimeSteppers.update_state!(model)
        f = model.forcing.ρθ.forcing
        @test f isa InversionFollowingNudging
        fields = Oceananigans.fields(model)
        θ̄ = vec(mean(interior(fields.θ), dims=(1, 2)))
        k_in = findfirst(==(2050.0), zc)
        @test f(1, 1, k_in, grid, model.clock, fields) ≈ -(θ̄[k_in] - 305.0) / 1800
        @test f(2, 3, k_in, grid, model.clock, fields) ≈ f(1, 1, k_in, grid, model.clock, fields)   # mean-based, uniform
        k_below = findfirst(==(650.0), zc)
        @test f(1, 1, k_below, grid, model.clock, fields) == 0                                       # masked below the inversion
        k_ramp = findfirst(==(1350.0), zc)
        @test f(1, 1, k_ramp, grid, model.clock, fields) ≈ -0.25 * (θ̄[k_ramp] - 305.0) / 1800
    end

    @testset "κ-Köhler activation from the local reservoir" begin
        FT = Float64
        constants = ThermodynamicConstants(FT)
        aerosol = kappa_aerosol_activation(FT; number_mixing_ratio=1e9, thermodynamic_constants=constants)
        @test aerosol isa AerosolActivation{FT, true}
        @test Breeze.Microphysics.PredictedParticleProperties.has_prognostic_aerosol(aerosol)
        mode = aerosol.modes[1]
        @test mode.mean_radius == 92.5e-9 && mode.geometric_std == 1.5 && mode.kappa == 0.2
        T = 285.0
        s_m = kappa_critical_supersaturation(mode, aerosol, T)
        # κ-Köhler: s_c = √(4A³ / (27 κ D³/8)) with A = 4 σ M_w / (R T ρ_w)
        σ = 0.0761 - 1.55e-4 * (T - constants.energy_reference_temperature)
        A = 4 * σ * constants.vapor.molar_mass / (constants.molar_gas_constant * T * constants.liquid.density)
        D = 185e-9
        @test s_m ≈ sqrt(4 * A^3 / (27 * 0.2 * D^3)) rtol=1e-10
        @test 1e-3 < s_m < 2e-3                                                 # ≈ 0.12 %
        @test activated_fraction(mode, aerosol, T, s_m) ≈ 0.5
        @test activated_fraction(mode, aerosol, T, 10 * s_m) > 0.99
        @test activated_fraction(mode, aerosol, T, s_m / 10) < 0.01
        # the activation target scales with the particles present, not with the seed number
        qᵛ⁺ˡ = 8e-3
        S = 3e-3
        qᵛ = qᵛ⁺ˡ * (1 + S)
        f = activated_fraction(mode, aerosol, T, S)
        for (nᶜˡ, nᵃ) in ((0.0, 1e8), (0.0, 1e9), (2e7, 1e8))
            r = aerosol_activation_rate(aerosol, nᶜˡ, nᵃ, qᵛ, qᵛ⁺ˡ, T)
            @test r.ncnuc ≈ max(0, f * (nᶜˡ + nᵃ) - nᶜˡ) / aerosol.activation_timescale
            @test r.qcnuc ≈ r.ncnuc * 4π / 3 * constants.liquid.density * (1e-6)^3
        end
        @test aerosol_activation_rate(aerosol, 0.0, 1e8, qᵛ⁺ˡ * (1 - 1e-3), qᵛ⁺ˡ, T).ncnuc == 0      # subsaturated
        @test aerosol_activation_rate(aerosol, 1.9e8, 1e7, qᵛ, qᵛ⁺ˡ, T).ncnuc == 0                   # already above the target f (nᶜˡ + nᵃ)
        @test aerosol_activation_rate(aerosol, 0.0, 0.0, qᵛ, qᵛ⁺ˡ, T).ncnuc == 0                     # empty reservoir
    end

    @testset "surface source, regeneration and subsidence in a P3 model" begin
        FT = Float64
        grid = RectilinearGrid(CPU(), FT; size=(8, 8, 20), x=(0, 800), y=(0, 800), z=(0, 2000), halo=(5, 5, 5), topology=(Periodic, Periodic, Bounded))
        zc = Array(znodes(grid, Center()))
        constants = ThermodynamicConstants(FT)
        reference_state = ReferenceState(grid, constants; base_pressure=101000, potential_temperature=290)
        dynamics = AnelasticDynamics(reference_state)
        ρᵣ = Array(interior(reference_state.density, 1, 1, :))
        aerosol = kappa_aerosol_activation(FT; number_mixing_ratio=1e8, thermodynamic_constants=constants)
        p3 = P3Microphysics(FT; cloud=Breeze.Microphysics.PredictedParticleProperties.CloudDroplets(FT; number_concentration=100e6), aerosol)
        times = [0.0, 3600.0]
        w = profile_time_series(grid, times, [fill(-0.01, 20), fill(-0.01, 20)])
        vadv = LargeScaleVerticalAdvection(w)
        evaporation_rate = CenterField(grid)
        regen = EvaporationRegeneration(evaporation_rate)
        source = SurfaceAerosolSource(7e5)
        forcing = (; θ=vadv, qᵛ=vadv, nᵃ=(vadv, source, regen.nᵃ), nᶜˡ=(vadv, regen.nᶜˡ), qᶜˡ=vadv)
        model = AtmosphereModel(grid; dynamics, microphysics=p3, thermodynamic_constants=constants, forcing)
        # cloud layer between 800 and 1200 m with droplets, interstitial aerosol everywhere
        set!(model; T=(x, y, z) -> 288 - 0.006z, qᵛ=(x, y, z) -> 800 < z < 1200 ? 9.5e-3 : 8e-3,
                    qᶜˡ=(x, y, z) -> 800 < z < 1200 ? 3e-4 : 0.0, nᶜˡ=(x, y, z) -> 800 < z < 1200 ? 8e7 : 0.0,
                    nᵃ=(x, y, z) -> 800 < z < 1200 ? 2e7 : 1e8)
        Oceananigans.TimeSteppers.update_state!(model; compute_tendencies=false)
        fields = Oceananigans.fields(model)
        # subsidence is applied exactly once per field: each materialized forcing holds one LargeScaleVerticalAdvection
        count_vadv(f::Breeze.Forcings.SpecificForcing) = count_vadv(f.forcing)
        count_vadv(f::LargeScaleVerticalAdvection) = 1
        count_vadv(f::Oceananigans.Forcings.MultipleForcings) = sum(count_vadv, f.forcings)
        count_vadv(f) = 0
        for name in (:ρθ, :ρqᵛ, :ρnᵃ, :ρnᶜˡ, :ρqᶜˡ)
            @test count_vadv(model.forcing[name]) == 1
        end
        @test !any(f -> f isa Breeze.Forcings.SubsidenceForcing, values(model.forcing))
        # the surface source: flux / (ρ₁ Δz₁) in the bottom cell only, ×ρ by Breeze
        Δz₁ = 100.0
        unwrap(f::Breeze.Forcings.SpecificForcing) = unwrap(f.forcing)
        unwrap(f) = f
        members(F) = F isa Oceananigans.Forcings.MultipleForcings ? map(unwrap, F.forcings) :
                     (G = unwrap(F); G isa Oceananigans.Forcings.MultipleForcings ? map(unwrap, G.forcings) : (G,))
        parts = members(model.forcing.ρnᵃ)
        src = parts[findfirst(f -> f isa SurfaceAerosolSource, parts)]
        @test src(1, 1, 1, grid, model.clock, fields) ≈ 7e5 / (ρᵣ[1] * Δz₁)
        @test src(1, 1, 2, grid, model.clock, fields) == 0
        # regeneration: with E = 1e-6 kg/kg/s in the cloud layer, R = nᶜˡ E / qᶜˡ leaves nᶜˡ and enters nᵃ
        set!(evaporation_rate, (x, y, z) -> 800 < z < 1200 ? 1e-6 : 0.0)
        rn = parts[findfirst(f -> f isa EvaporationRegeneration, parts)]
        k_cloud = findfirst(z -> 800 < z < 1200, zc)
        R = 8e7 * 1e-6 / 3e-4
        @test rn(1, 1, k_cloud, grid, model.clock, fields) ≈ R
        @test rn(1, 1, 1, grid, model.clock, fields) == 0
        cparts = members(model.forcing.ρnᶜˡ)
        rc = cparts[findfirst(f -> f isa EvaporationRegeneration, cparts)]
        @test rc(1, 1, k_cloud, grid, model.clock, fields) ≈ -R
        # the whole aerosol-number forcing is conservative: source + regeneration sum to the flux
        total = sum(rn(1, 1, k, grid, model.clock, fields) * ρᵣ[k] * Δz₁ + rc(1, 1, k, grid, model.clock, fields) * ρᵣ[k] * Δz₁ for k in 1:20)
        @test abs(total) < 1e-6 * 7e5
        # the diagnosed evaporation rate is zero in saturated/clear cells and finite
        E = cloud_evaporation_rate_field(model)
        compute!(E)
        @test all(isfinite, interior(E)) && all(interior(E) .≥ 0)
        # ten steps in a clear, subsaturated column: the column aerosol number grows by the surface
        # flux and nothing else (in cloud, P3's cloud self-collection is a genuine coalescence sink)
        set!(model; T=(x, y, z) -> 288 - 0.006z, qᵛ=3e-3, qᶜˡ=0, nᶜˡ=0, nᵃ=1e8)   # subsaturated everywhere
        set!(evaporation_rate, 0)
        cols = aerosol_number_columns(model)
        compute!(cols.total); n0 = mean(interior(cols.total))
        for _ in 1:10
            time_step!(model, 1.0)
        end
        compute!(cols.total); n1 = mean(interior(cols.total))
        @test n1 - n0 ≈ 7e5 * 10 rtol=1e-3
        compute!(cols.nᶜˡ); @test mean(interior(cols.nᶜˡ)) == 0
    end

    if HAVE_SEASTARR
        @testset "staged drivers and constructor" begin
            d = read_dephy_driver(joinpath(SEASTARR_DIR, "SEA_STARR_CTRL_SCM_driver.nc"))
            @test d.case == "SEA_STARR/CTRL"
            @test d.start_time == DateTime(2017, 8, 15, 21) && d.end_time == DateTime(2017, 8, 18, 15)
            @test length(d.times) == 67 && d.times[end] == 66 * 3600
            @test length(d.z) == 336 && d.z[1] == 5.0 && all(diff(d.z[1:250]) .≈ 10)
            @test isapprox(d.latitude, -11.649; atol=1e-3) && isapprox(d.longitude, -8.336; atol=1e-3)
            @test d.initial.na[1] == 1e8                                   # 100 mg⁻¹ in the boundary layer
            @test 9.9e8 < maximum(d.initial.na) < 1.2e9                     # smoke layer
            @test inversion_height(d.z, d.initial.thetal; z_max=3000) == 1160.0
            @test d.nudging_rates.thetal == 1 / 1800 && d.nudging_rates.u == 1 / 10800
            @test isapprox(d.sst[1], 293.33; atol=0.01) && isapprox(d.sst[end], 298.06; atol=0.01)
            @test all(d.forcing.w[:, 1] .== d.forcing.w[:, end])            # time-constant subsidence profile
            @test d.forcing.ug == d.forcing.u_nud
            n100 = read_dephy_driver(joinpath(SEASTARR_DIR, "SEA_STARR_N100_SCM_driver.nc"))
            @test all(n100.initial.na .== 1e8) && n100.initial.thetal == d.initial.thetal
            @test_throws ArgumentError sea_starr(; member=:N100, driver_path=joinpath(SEASTARR_DIR, "SEA_STARR_CTRL_SCM_driver.nc"), arch=CPU(), write_output=false)

            case = sea_starr(; member=:CTRL, data_dir=SEASTARR_DIR, arch=CPU(), FT=Float32, Nx=8, Ny=8,
                               z_faces=collect(range(0, 3000, length=31)), stop_time=2.0, write_output=true,
                               output_dir=mktempdir(), checkpoint_interval=nothing, progress_interval=100,
                               statistics_interval=2.0, timeseries_interval=1.0, fields_2d_interval=2.0, fields_3d_interval=2.0)
            @test case.model.clock.time == 0 && case.model.clock.iteration == 0
            @test case.config.member == "CTRL" && case.config.aerosol_surface_flux == 7e5
            @test case.config.thermodynamic_nudging_timescale == 1800 && case.config.wind_nudging_timescale == 10800
            @test case.config.roughness_length == 1e-4
            @test length(case.config.departures) ≥ 5
            @test !any(occursin("zero downwelling LW", dep) for dep in case.config.departures)
            μ = case.model.microphysical_fields
            @test haskey(μ, :ρnᵃ) && haskey(μ, :ρnᶜˡ)
            nᵃ = Array(interior(μ.nᵃ, 1, 1, :)); nᶜˡ = Array(interior(μ.nᶜˡ, 1, 1, :)); qᶜˡ = Array(interior(μ.qᶜˡ, 1, 1, :))
            # equilibrium start: cloud where (θₗ, qₜ) is saturated, droplets = aerosol there, reservoir elsewhere
            @test any(qᶜˡ .> 0)
            @test all(nᶜˡ[qᶜˡ .> 0] .≈ Float32(1e8)) && all(nᵃ[qᶜˡ .> 0] .≈ 0)
            @test nᵃ[1] ≈ Float32(1e8) && nᶜˡ[1] == 0
            @test maximum(nᵃ) > 9e8                                         # the smoke layer above the inversion
            @test case.mask_updater.inversion_height_max[] ≈ 1150 atol=100   # ≈ driver inversion (1160 m)
            @test Set(keys(case.simulation.output_writers)) == Set((:statistics, :timeseries, :fields_2d, :fields_3d))
            # radiation: the driver's upper atmosphere above the LES top and the trajectory solar position
            @test case.config.radiation == "rrtmgp_extended" && case.config.radiation_layers_above > 10 && case.config.radiation_column_top > 70000
            up = upper_atmosphere_layers(d, 3000.0)
            @test up.faces[1] == 3000.0 && issorted(up.z) && all(diff(up.p_faces) .< 0) && all(up.T .> 150) && all(up.q .≥ 0)
            r = case.model.radiation
            @test r.solar_position isa Breeze.AtmosphereModels.FixedCosineZenith
            cb = [c.func for c in values(case.simulation.callbacks) if c.func isa TrajectorySolarPosition]
            @test length(cb) == 1
            cosz, lon, lat = trajectory_cos_zenith(cb[1], 15 * 3600.0)        # 2017-08-16 12 UTC, lon ≈ −4.4°E → ≈ 11:40 LST
            @test 0.8 < cosz < 0.95 && -5.5 < lon < -3.5 && -15 < lat < -13
            @test trajectory_cos_zenith(cb[1], 0.0)[1] == 0                   # 21 UTC: night
            Nz = case.model.grid.Nz
            case.model.clock.time = 15 * 3600.0; cb[1](case.model)
            Breeze.AtmosphereModels.update_radiation!(r, case.model)
            lwd_top = -mean(interior(r.downwelling_longwave_flux, :, :, Nz + 1))
            swd_top = -mean(interior(r.downwelling_shortwave_flux, :, :, Nz + 1))
            @test 50 < lwd_top < 250                                           # downwelling LW from the atmosphere above the LES top
            @test 0.85 * 1361 * cosz < swd_top < 1361 * cosz                   # SW attenuated above the LES top (O₃, Rayleigh)
            @test all(isfinite, interior(r.flux_divergence))
            plain = sea_starr(; member=:CTRL, data_dir=SEASTARR_DIR, arch=CPU(), FT=Float32, Nx=8, Ny=8,
                                z_faces=collect(range(0, 3000, length=31)), stop_time=2.0, write_output=false,
                                radiation=:rrtmgp, solar=:fixed, checkpoint_interval=nothing, progress_interval=100)
            @test plain.config.radiation_layers_above == 0
            Breeze.AtmosphereModels.update_radiation!(plain.model.radiation, plain.model)
            @test -mean(interior(plain.model.radiation.downwelling_longwave_flux, :, :, Nz + 1)) == 0
            case.model.clock.time = 0.0; cb[1](case.model)
            run!(case.simulation)
            @test case.model.clock.time ≈ 2.0
            @test all(f -> all(isfinite, Array(interior(f))), values(Oceananigans.prognostic_fields(case.model)))
            path = joinpath(case.config.output_dir, "provenance.toml")
            write_provenance(path, case)
            @test isfile(path)
        end
    else
        @test_skip HAVE_SEASTARR
        @info "SEA STARR drivers not staged; skipping the staged-driver tests (julia data_wrangling/fetch_manifest.jl cases/seastarr/inputs.toml data/seastarr_22241697)"
    end
end
