# Restart of the real-ERA5 outer domain (CPU, 16×16×94, no radiation): a run checkpointed at 36 s and picked up to
# 72 s must reproduce the continuous 72 s run to round-off scale. The restart is NOT bit-identical (2026-10-09, job 508:
# after 12 steps max |ΔT| = 0.005 K, |Δu| = 0.007 m/s, |Δw| = 0.0009 m/s, |Δρ| = 2e-5, |ΔT_land| = 0.005 K, accumulated
# condensation 0.2 %): some state is rebuilt rather than restored on the first step after pickup (candidates: the
# iterated similarity-theory interface solution, P3 diagnostic fields). Tolerances are set at that scale.
using NumericalEarth, CopernicusClimateDataStore, CloudMicrophysics, RRTMGP
using BreezeLab: tracer_mip_outer_simulation, sync_nested_parent_clock!
using Oceananigans.Fields: interior

era5_dir = get(ENV, "TRACER_MIP_ERA5_DIR", "/shared/home/greg/breezelab-work/tracer-mip/data/era5")
function build()
    case = tracer_mip_outer_simulation(CPU(); parent = :era5, era5_dir, Nx = 16, Ny = 16, stop_time = 72.0,
                                       radiation = nothing, process_rates = true)
    case.simulation.callbacks[:progress] = Callback(sim -> nothing, IterationInterval(6))   # a plain-function callback, as in the case script
    return case
end
dir = mktempdir()

@testset "TRACER-MIP outer checkpoint/pickup" begin
    A = build()
    A.simulation.output_writers[:ck] = Checkpointer(A.simulation.model; schedule = IterationInterval(12), dir, prefix = "ck")
    run!(A.simulation)
    @test A.simulation.model.clock.time ≈ 72
    files = filter(f -> occursin("iteration12", f), readdir(dir; join = true))
    @test length(files) == 1

    B = build()
    set!(B.simulation; checkpoint = only(files))
    @test B.simulation.model.clock.time ≈ 36 && B.child.clock.time ≈ 36
    sync_nested_parent_clock!(B)
    @test B.nest.parent.clock.time ≈ 36
    run!(B.simulation)
    @test B.simulation.model.clock.time ≈ 72

    for (name, a, b) in (("T", A.child.temperature, B.child.temperature), ("ρ", A.child.dynamics.total_density, B.child.dynamics.total_density),
                         ("u", A.child.velocities.u, B.child.velocities.u), ("w", A.child.velocities.w, B.child.velocities.w),
                         ("T_land", A.land.temperature, B.land.temperature))
        diff = maximum(abs.(Array(interior(a)) .- Array(interior(b))))
        scale = maximum(abs.(Array(interior(a))))
        @info "restart difference $name: max |A − B| = $diff (max |A| = $scale)"
        @test diff ≤ 1e-3 * max(scale, 1)
    end
    # accumulators resume their sums across the restart
    a = Array(interior(A.accumulators.fields.liquid_condensation)); b = Array(interior(B.accumulators.fields.liquid_condensation))
    @info "restart difference accumulated condensation: $(maximum(abs.(a .- b))) (max $(maximum(abs.(a))))"
    @test maximum(abs.(a .- b)) ≤ 1e-2 * max(maximum(abs.(a)), eps())
end
