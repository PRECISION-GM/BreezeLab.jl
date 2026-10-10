# Rain over a land column must ADD water to the slab bucket (the coupler's Jʳⁿ is positive downward).
# Synthetic parent (software test only), CPU, two identical outer domains, one with a rain shaft over a land cell.
using NumericalEarth
using BreezeLab: tracer_mip_outer_simulation, tracer_mip_grid
using Oceananigans.Fields: interior
using Oceananigans.TimeSteppers: update_state!

ext = Base.get_extension(BreezeLab, :BreezeLabNumericalEarthExt)

@testset "TRACER-MIP rain reaches the slab land with the right sign" begin
    previous_FT = Oceananigans.defaults.FloatType
    Oceananigans.defaults.FloatType = Float32   # restored in the `finally` below
    try
        Nx, Ny = 10, 8
        z_faces = collect(range(0, 12_000, length = 25))
        grid = tracer_mip_grid(CPU(), :outer; Nx, Ny, z_faces)
        build() = tracer_mip_outer_simulation(CPU(); parent = :synthetic, Nx, Ny, z_faces, stop_time = 30.0,
                                              radiation = false, closure = nothing, process_rates = false,
                                              terrain = ext.synthetic_coastal_elevation(grid))
        wet = build(); dry = build()
        mask = Array(interior(wet.sea_mask))[:, :, 1]
        i, j = 3, Ny ÷ 2                       # west of the synthetic coastline: land
        @test mask[i, j] == 0

        # A rain shaft over (i, j): 2 g/kg of 1 mm drops in the lowest three cells (1.5 km)
        μ = wet.child.microphysical_fields
        ρ = Array(interior(wet.child.dynamics.total_density))
        qʳ = 2e-3; nʳ = qʳ / (4 / 3 * π * 1000 * (0.5e-3)^3)
        for k in 1:3
            interior(μ.ρqʳ)[i, j, k] = ρ[i, j, k] * qʳ
            interior(μ.ρnʳ)[i, j, k] = ρ[i, j, k] * nʳ
        end
        update_state!(wet.nest); update_state!(wet.model)

        Jʳⁿ = Array(interior(wet.model.interfaces.exchanger.atmosphere.state.Jʳⁿ))[:, :, 1]
        @info "coupler rain flux at the shaft: $(Jʳⁿ[i, j]) kg m⁻² s⁻¹ (expected ≈ ρ qʳ wʳ ≈ 0.01)"
        @test Jʳⁿ[i, j] > 1e-3                 # positive downward
        @test all(Jʳⁿ .>= 0)

        M₀ = Array(interior(wet.land.water_storage))[i, j, 1]
        run!(wet.simulation); run!(dry.simulation)
        ΔM_wet = Array(interior(wet.land.water_storage))[i, j, 1] - M₀
        ΔM_dry = Array(interior(dry.land.water_storage))[i, j, 1] - M₀
        @info "bucket change over 30 s at the shaft: with rain $(ΔM_wet), without $(ΔM_dry) kg m⁻²"
        @test ΔM_wet - ΔM_dry > 0.1            # the shaft's rain (≈ 0.3 kg m⁻² in 30 s) was added, not removed
    finally
        Oceananigans.defaults.FloatType = previous_FT
    end
end
