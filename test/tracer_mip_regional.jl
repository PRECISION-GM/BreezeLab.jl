# Regional nesting / open-boundary / land-sea exchange software tests for the TRACER-MIP outer domain.
# They need NumericalEarth (the case sub-environment `cases/tracer_mip`), so `runtests.jl` includes this
# file only when NumericalEarth is loadable. The parent is the SYNTHETIC analytic state: these tests
# check the machinery (interpolation, hydrostatic initialization, boundary transfer, coupling, stepping),
# not any MIP result.
using NumericalEarth
using NumericalEarth: PrescribedAtmosphere, NestedModel, AtmosphereLandModel
using BreezeLab: tracer_mip_outer_simulation, tracer_mip_grid, tracer_mip_protocol
using Oceananigans.Units: Time
using Oceananigans.Fields: interior
using Statistics: mean, maximum, minimum

ext = Base.get_extension(BreezeLab, :BreezeLabNumericalEarthExt)

@testset "TRACER-MIP regional machinery (synthetic parent, software test only)" begin
    @test !isnothing(ext)
    arch = CPU()
    FT = Float32
    Nx, Ny = 10, 8
    z_faces = collect(range(0, 12_000, length = 25))          # reduced vertical for speed; protocol faces are tested elsewhere
    protocol = tracer_mip_protocol()
    Oceananigans.defaults.FloatType = FT   # restored in the `finally` below
    try
    grid = tracer_mip_grid(arch, :outer; Nx, Ny, z_faces)

    @testset "synthetic parent and sea mask helpers" begin
        parent = ext.synthetic_parent_atmosphere(grid; times = [0.0, 3600.0, 7200.0])
        @test parent.source isa ext.SyntheticTestParent
        T = parent.temperature[1]
        p = parent.pressure[1]
        @test minimum(interior(T)) > 190 && maximum(interior(T)) ≤ 300
        @test all(diff(Array(interior(p, 1, 1, :))) .< 0)       # pressure decreases with height
        elevation = ext.synthetic_coastal_elevation(grid)
        @test minimum(interior(elevation)) == 0 && maximum(interior(elevation)) > 0
        flat_mask = ext.sea_mask(grid)
        @test all(==(1), interior(flat_mask))                   # flat grid: all sea
        α = ext.surface_albedo_field(flat_mask)
        @test all(≈(0.06f0), interior(α))
    end

    @testset "outer-domain constructor (synthetic boundaries, no radiation)" begin
        case = tracer_mip_outer_simulation(arch; parent = :synthetic, Nx, Ny, z_faces,
                                           stop_time = 20.0, radiation = false, closure = nothing,
                                           terrain = ext.synthetic_coastal_elevation(grid))
        @test case.config.exploratory_synthetic_boundaries == true
        @test case.model isa NumericalEarth.EarthSystemModel
        @test case.nest isa NestedModel
        child = case.child
        @test size(child.grid) == (Nx, Ny, 24)
        # coastal mask: sea east of the coastline, land west
        mask = Array(interior(case.sea_mask))
        @test any(==(1), mask) && any(==(0), mask)
        # land pinned at sea: temperature equals the idealized SST there, bucket saturated
        T_land = Array(interior(case.land.temperature)); 𝒮 = Array(interior(case.land.saturation))
        sea = mask .== 1
        @test all(T_land[sea] .≈ case.pinning.sea_surface_temperature)
        @test all(𝒮[sea] .≈ 1)
        @test all(𝒮[.!sea] .< 1)
        # child initialized from the parent: temperature close to the parent profile, finite density
        ρ = child.dynamics.total_density
        @test all(isfinite, interior(ρ)) && minimum(interior(ρ)) > 0
        Tc = Array(interior(child.temperature))
        @test 190 < minimum(Tc) && maximum(Tc) < 305
        @test abs(mean(Tc[:, :, 1]) - 300) < 3                  # near-surface temperature from the parent
        # open-boundary transfer: west boundary condition is parent-derived
        bcs = child.momentum.ρu.boundary_conditions
        @test !(bcs.west isa Oceananigans.BoundaryConditions.DefaultBoundaryCondition)
        @test !(child.dynamics isa Breeze.AnelasticDynamics)    # compressible, as in the nested path
        # Tier-1 aerosol is in place
        μ = child.microphysical_fields
        @test haskey(μ, :nᵃ) && haskey(μ, :ρnᶜˡ) && !haskey(μ, :ρnᵃ)
        # a few coupled steps: finite state, boundary values track the parent, no NaNs
        run!(case.simulation)
        @test case.simulation.model.clock.time ≈ 20
        @test all(isfinite, interior(child.temperature)) && all(isfinite, interior(ρ))
        @test all(isfinite, interior(case.land.temperature))
        u = Array(interior(child.velocities.u))
        @test abs(mean(u) - 3) < 1.5                            # the parent's uniform 3 m/s zonal wind persists
        # rain-to-land shim: the coupler's rain diagnostic is Breeze's (sign-flipped) bottom precipitation flux
        Jʳⁿ = case.model.interfaces.exchanger.atmosphere.state.Jʳⁿ
        @test Jʳⁿ isa Field && !isnothing(Jʳⁿ.operand)
        @test all(interior(Jʳⁿ) .>= 0)
        acc = case.accumulators
        @test all(f -> all(isfinite, interior(f)), acc.fields)
        T_land = Array(interior(case.land.temperature))
        @test all(T_land[sea] .≈ case.pinning.sea_surface_temperature)   # still pinned after stepping
    end
    finally
        Oceananigans.defaults.FloatType = Float64
    end
end
