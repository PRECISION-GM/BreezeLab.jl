# Checks of the Covert sonde-nudging + RRTMGP (LW+SW) path without time stepping:
#   1. nudging tendencies (θ, moisture, u, v) evaluated with the materialized forcing kernels at several times;
#   2. RRTMGP shortwave over a day (TOA and surface downwelling SW, cos zenith from the apparent solar position).
#   julia --project=<tree> check_forcing.jl <data dir>
using BreezeLab, Breeze, Oceananigans, Oceananigans.Units, CUDA, Statistics, Printf
using Oceananigans.AbstractOperations: KernelFunctionOperation
using Breeze.Microphysics.PredictedParticleProperties: CloudDroplets

data = ARGS[1]
Oceananigans.defaults.FloatType = Float32
r = COVERT_RESOLUTIONS.r128
case = ena_covert(; arch = GPU(), data_dir = data, microphysics = P3Microphysics(; cloud = CloudDroplets(; number_concentration = 75e6)),
                    Nx = 32, Ny = 32, Nz = r.Nz, Lx = r.Lx / 4, Ly = r.Ly / 4, Δt = 1.0, max_Δt = 5.0, cfl = 0.8,
                    closure = nothing, scalar_advection_order = 9, radiation = :rrtmgp,
                    nudging_tq = 3600, nudging_uv = 3600, write_output = false, stop_time = 48hours)
model = case.model
grid = model.grid
println(case.config.radiation); println(case.config.tq_nudging); println(case.config.wind_nudging)

nudgings = Pair{Symbol, Any}[]
function walk!(name, x)
    if x isa BreezeLab.MeanProfileNudging
        push!(nudgings, name => x)
    elseif x isa Tuple
        foreach(y -> walk!(name, y), x)
    elseif hasfield(typeof(x), :forcings)
        walk!(name, x.forcings)
    elseif hasfield(typeof(x), :forcing)
        walk!(name, x.forcing)
    end
end
for (name, f) in pairs(model.forcing)
    walk!(name, f)
end
println("nudgings found: ", [first(p) for p in nudgings])

z = Array(znodes(grid, Center()))
levels = [findmin(abs.(z .- h))[2] for h in (300, 1000, 1300, 1500, 1800, 2500, 4000, 8000, 15000, 19500)]
scale = Dict(:ρθ => 86400, :ρqᵛ => 86400e3, :ρu => 86400, :ρv => 86400, :ρE => 86400)
unit = Dict(:ρθ => "K/day", :ρqᵛ => "g/kg/day", :ρu => "m/s/day", :ρv => "m/s/day")

for t in (0.0, 3hours, 6hours, 12hours, 24hours, 47hours, 50hours)
    model.clock.time = t
    println(@sprintf("\n=== t = %.1f h", t / 3600))
    for (name, f) in nudgings
        Oceananigans.Forcings  # (namespace touch)
        BreezeLab.AtmosphereModels.compute_forcing!(f)
        tendency = Field(KernelFunctionOperation{Center, Center, Center}(f, grid, model.clock, nothing))
        compute!(tendency)
        column = Array(interior(tendency, 1, 1, :))
        mean_profile = Array(interior(f.averaged_field, 1, 1, :))
        target = [CUDA.@allowscalar(f.target[1, 1, k, Oceananigans.Units.Time(t)]) for k in 1:length(z)]
        s = get(scale, name, 1)
        println("  ", name, " [", get(unit, name, "?"), "]  (z: tendency | mean → target)")
        for k in levels
            println(@sprintf("    z=%7.0f  %+9.3f | %10.5g → %10.5g", z[k], s * column[k], mean_profile[k], target[k]))
        end
        # tendency must equal −𝒲 (mean − target)/τ (and be identical in every column)
        all_columns = Array(interior(tendency))
        spread = maximum(abs.(all_columns .- all_columns[1:1, 1:1, :]))
        println(@sprintf("    column spread %.2e, max |tendency| %.3g %s", spread * s, maximum(abs.(column)) * s, get(unit, name, "")))
    end
end

println("\n=== RRTMGP shortwave over 24 h (domain means), UTC = 06 + t")
rad = model.radiation
Nz = size(grid, 3)
for h in 0:2:24
    model.clock.time = h * 3600.0
    Breeze.AtmosphereModels._update_radiation!(rad, model)
    sw_dn = Array(interior(rad.downwelling_shortwave_flux))
    sw_up = Array(interior(rad.upwelling_shortwave_flux))
    lw_up = Array(interior(rad.upwelling_longwave_flux))
    lw_dn = Array(interior(rad.downwelling_longwave_flux))
    heating = Array(interior(rad.flux_divergence))
    println(@sprintf("  UTC %05.2f  SW↓ top %7.1f  SW↓ sfc %7.1f  SW↑ sfc %6.1f  LW↑ sfc %6.1f  LW↓ sfc %6.1f  W/m²  (column flux divergence %.1f W/m²)",
                     mod(6 + h, 24), -mean(sw_dn[:, :, Nz + 1]), -mean(sw_dn[:, :, 1]), mean(sw_up[:, :, 1]),
                     mean(lw_up[:, :, 1]), -mean(lw_dn[:, :, 1]), sum(mean(heating, dims = (1, 2))[1, 1, :] .* diff(Array(znodes(grid, Face()))))))
end
