# Offline one-way nest: parent from a saved outer inner-region state (CPU smoke job 213 output), inner child.
using NumericalEarth
using NumericalEarth: PrescribedAtmosphere
using BreezeLab: tracer_mip_inner_simulation, tracer_mip_grid, tracer_mip_protocol
using Oceananigans.Fields: interior
using Statistics: mean
ext = Base.get_extension(BreezeLab, :BreezeLabNumericalEarthExt)
outer_dir = get(ENV, "TRACER_MIP_OUTER_RUN", "/shared/home/greg/breezelab-work/tracer-mip/runs/era5_cpu_smoke_job213")

@testset "offline inner nest from saved outer state" begin
    parent = ext.outer_run_parent(outer_dir; FT = Float32)
    @test parent.source isa ext.OuterRunParent
    @test NumericalEarth.surface_elevation(parent.grid) !== nothing
    T = parent.temperature[1]
    @test 190 < minimum(interior(T)) && maximum(interior(T)) < 320
    @test all(diff(Array(interior(parent.pressure[1], 1, 1, :))) .< 0)
    z = [NumericalEarth.Grids.znode(1, 1, k, parent.grid, Center(), Center(), Center()) for k in 1:size(parent.grid, 3)]
    @test issorted(z) && z[end] > 20_000
    case = tracer_mip_inner_simulation(CPU(); outer_run_dir = outer_dir, FT = Float32, Nx = 8, Ny = 8,
                                       z_faces = collect(range(0, 12_000, length = 25)), stop_time = 6.0,
                                       radiation = nothing, terrain = nothing, process_rates = false)
    @test case.config.nesting[1:7] == "one-way"
    child = case.child
    @test all(isfinite, interior(child.temperature))
    @test abs(mean(Array(interior(child.temperature))[:, :, 1]) - mean(Array(interior(parent.temperature[1]))[:, :, 1])) < 5
    run!(case.simulation)
    @test case.simulation.model.clock.time ≈ 6
    @test all(isfinite, interior(child.temperature)) && minimum(interior(child.dynamics.total_density)) > 0
end
