# # TRACER-MIP outer domain (Grid-1) — August 7, 2022 Tier 1 control
#
# Builds the regional 2 km outer domain of the TRACER-MIP control with `tracer_mip_outer_simulation`
# (ERA5 pressure-level parent, open boundaries with Davies relaxation, terrain-following grid on ETOPO2022,
# P3 with the prescribed Tier-1 two-mode aerosol, RRTMGP every 60 s, slab land with sea cells pinned to
# the ERA5 skin temperature), writes the protocol's hourly outer-domain output plus provenance, runs, and
# plots a surface snapshot. Run it from the case sub-environment:
#
#     julia --project=cases/tracer_mip cases/tracer_mip.jl
#
# The inner 500 m nest is a separate, one-way (offline) step that consumes this run's output; it is not
# part of this script. Reduce `Nx`, `Ny` and `stop_time` for a smoke test; the full 750² domain needs a
# large GPU (see docs/cases/tracer_mip.md and the status file for measured memory/throughput).

using BreezeLab
using NumericalEarth
using CopernicusClimateDataStore
using CloudMicrophysics
using RRTMGP
using CUDA
using Oceananigans
using Oceananigans.Units
using Breeze
using CairoMakie
using Printf
using Dates: DateTime, now, UTC
using TOML: TOML

# ## Settings (edit here)

# Pilots override the protocol size/duration through the environment (the Slurm launcher sets them):
# TRACER_MIP_NX/NY (cells), TRACER_MIP_STOP_HOURS, TRACER_MIP_PARENT = era5 | synthetic (synthetic is
# an exploratory software test of the machinery, never a MIP control), TRACER_MIP_ARCH = gpu | cpu.
arch = get(ENV, "TRACER_MIP_ARCH", "gpu") == "cpu" ? CPU() : GPU()
case = :aug07
Nx = parse(Int, get(ENV, "TRACER_MIP_NX", "750"))
Ny = parse(Int, get(ENV, "TRACER_MIP_NY", "750"))
stop_time = parse(Float64, get(ENV, "TRACER_MIP_STOP_HOURS", "24")) * 3600
parent = Symbol(get(ENV, "TRACER_MIP_PARENT", "era5"))
# Pilot/diagnostic knobs (protocol values by default): radiation rrtmgp|none, terrain etopo|flat, Δt in seconds.
radiation = get(ENV, "TRACER_MIP_RADIATION", "rrtmgp") == "none" ? nothing : :rrtmgp
terrain_choice = get(ENV, "TRACER_MIP_TERRAIN", "etopo")
Δt = parse(Float64, get(ENV, "TRACER_MIP_DT", "3"))
closure_choice = get(ENV, "TRACER_MIP_CLOSURE", "tke")   # tke (default, PBL) | none | smagorinsky
era5_dir = get(ENV, "TRACER_MIP_ERA5_DIR", joinpath(pkgdir(BreezeLab), "data", "era5"))
output_dir = get(ENV, "TRACER_MIP_OUTPUT_DIR", joinpath(pkgdir(BreezeLab), "output", "tracer_mip_outer_$(case)" * (parent === :era5 ? "" : "_EXPLORATORY_$(parent)")))
output_interval = 1hour             # Grid-1: 60-min full output
slice_interval = 10minutes          # lightweight monitoring slices
mkpath(output_dir)

# ## Build

extra = terrain_choice == "flat" ? (; terrain = nothing) : NamedTuple()
closure_choice == "none" && (extra = merge(extra, (; closure = nothing)))
closure_choice == "smagorinsky" && (extra = merge(extra, (; closure = SmagorinskyLilly(Float32))))
run_case = tracer_mip_outer_simulation(arch; case, parent, era5_dir, Nx, Ny, stop_time, radiation, Δt, extra...)
parent === :era5 || @warn "EXPLORATORY run with synthetic boundaries: software test of the machinery, not a TRACER-MIP control"
simulation = run_case.simulation
child = run_case.child
land = run_case.land
acc = run_case.accumulators

# ## Output: the roadmap's 3-D state/water fields hourly, 2-D surface fields, process accumulators

u, v, w = child.velocities
μ = child.microphysical_fields
state_3d = (; ρ = child.dynamics.total_density, T = child.temperature, p = Breeze.AtmosphereModels.dynamics_pressure(child.dynamics),
              u, v, w, qᵛ = μ.qᵛ, qᶜˡ = μ.qᶜˡ, nᶜˡ = μ.nᶜˡ, qʳ = μ.qʳ, nʳ = μ.nʳ, qⁱ = μ.qⁱ, nⁱ = μ.nⁱ, qᶠ = μ.qᶠ, bᶠ = μ.bᶠ, qʷⁱ = μ.qʷⁱ, nᵃ = μ.nᵃ)
process_3d = isnothing(acc) ? NamedTuple() : acc.fields
surface_2d = (; T_land = land.temperature, saturation = land.saturation, sea = run_case.sea_mask,
                rain = surface_rain_flux(child))
radiation_2d = isnothing(run_case.radiation) ? NamedTuple() :
               (; SW_down = run_case.radiation.downwelling_shortwave_flux, LW_down = run_case.radiation.downwelling_longwave_flux)

simulation.output_writers[:state] = JLD2Writer(child, state_3d; schedule = TimeInterval(output_interval),
                                               filename = joinpath(output_dir, "outer_state.jld2"), overwrite_files = true)
simulation.output_writers[:process] = JLD2Writer(child, process_3d; schedule = TimeInterval(output_interval),
                                                 filename = joinpath(output_dir, "outer_process_rates.jld2"), overwrite_files = true)
isnothing(acc) || add_callback!(simulation, process_rate_output_callback(acc, TimeInterval(output_interval)))
simulation.output_writers[:surface] = JLD2Writer(child, merge(surface_2d, radiation_2d); schedule = TimeInterval(slice_interval),
                                                 filename = joinpath(output_dir, "outer_surface.jld2"), overwrite_files = true)
# Inner-region prognostic state for the offline one-way 500 m nest: the raw thermodynamic/wind state
# (what a `PrescribedAtmosphere` parent needs) over the inner extent plus a halo of `inner_halo`
# outer cells (relaxation zone + interpolation stencil), every `inner_interval` (≥ 10 min).
inner_interval = 10minutes
inner_halo = 10
inner_extent = tracer_mip_horizontal_extent(tracer_mip_protocol(), :inner)
λc = Array(Oceananigans.Grids.λnodes(child.grid, Center()))
φc = Array(Oceananigans.Grids.φnodes(child.grid, Center()))
i_range = max(1, searchsortedfirst(λc, inner_extent.longitude[1]) - inner_halo):min(Nx, searchsortedlast(λc, inner_extent.longitude[2]) + inner_halo)
j_range = max(1, searchsortedfirst(φc, inner_extent.latitude[1]) - inner_halo):min(Ny, searchsortedlast(φc, inner_extent.latitude[2]) + inner_halo)
inner_dir = joinpath(output_dir, "inner_region"); mkpath(inner_dir)
# velocities moved to cell centers so the saved window has the same (i, j) extent for every field
inner_state = (; T = child.temperature, p = Breeze.AtmosphereModels.dynamics_pressure(child.dynamics), qᵛ = μ.qᵛ,
                 qᶜˡ = μ.qᶜˡ, qʳ = μ.qʳ, qⁱ = μ.qⁱ, u = @at((Center, Center, Center), u), v = @at((Center, Center, Center), v),
                 w = @at((Center, Center, Center), w), ρ = child.dynamics.total_density)
simulation.output_writers[:inner_region] = JLD2Writer(child, inner_state; schedule = TimeInterval(inner_interval),
                                                      indices = (i_range, j_range, :),
                                                      filename = joinpath(inner_dir, "inner_region_state.jld2"), overwrite_files = true)
open(joinpath(inner_dir, "README.txt"), "w") do io
    println(io, "Outer-domain state saved on the inner nest region every $(prettytime(inner_interval)) for the offline one-way nest.")
    println(io, "Outer cell index ranges: i = $(i_range), j = $(j_range) (inner extent plus $(inner_halo) outer cells); all 94 levels.")
    println(io, "Fields: T (K), p (Pa), qᵛ qᶜˡ qʳ qⁱ (kg/kg), u v w (m/s, cell centers), ρ (kg/m³) on the outer terrain-following lat-lon grid.")
end

k_aloft = searchsortedfirst(Array(znodes(child.grid, Center())), 2000)
slices = (; w_2km = view(w, :, :, k_aloft), qᶜˡ_2km = view(μ.qᶜˡ, :, :, k_aloft), T_sfc = view(child.temperature, :, :, 1),
            u_sfc = view(u, :, :, 1), v_sfc = view(v, :, :, 1))
simulation.output_writers[:slices] = JLD2Writer(child, slices; schedule = TimeInterval(slice_interval),
                                                filename = joinpath(output_dir, "outer_slices.jld2"), overwrite_files = true)

wall = Ref(time_ns())
function progress(sim)
    ρ = child.dynamics.total_density
    p = Breeze.AtmosphereModels.dynamics_pressure(child.dynamics)
    @info @sprintf("iter %6d  t = %s  Δt = %s  wall %s | max|u| %.1f max|w| %.2f  T ∈ [%.1f, %.1f]  ρ ∈ [%.3f, %.3f]  p ∈ [%.0f, %.0f]  max qᶜˡ %.2e  max nᶜˡ %.2e",
                   iteration(sim), prettytime(sim), prettytime(sim.Δt), prettytime(1e-9 * (time_ns() - wall[])),
                   maximum(abs, u), maximum(abs, w), minimum(child.temperature), maximum(child.temperature),
                   minimum(ρ), maximum(ρ), minimum(p), maximum(p), maximum(μ.qᶜˡ), maximum(μ.nᶜˡ))
end
add_callback!(simulation, progress, TimeInterval(parse(Float64, get(ENV, "TRACER_MIP_PROGRESS_SECONDS", "300"))))
Oceananigans.Diagnostics.erroring_NaNChecker!(simulation)

# ## Provenance

provenance = Dict("generated_utc" => string(now(UTC)), "config" => Dict(string(k) => BreezeLab.toml_value(v) for (k, v) in pairs(run_case.config)),
                  "era5_manifest" => joinpath(era5_dir, "MANIFEST_era5_$(case).toml"),
                  "software" => Dict("BreezeLab" => BreezeLab.package_source_description(BreezeLab),
                                     "Breeze" => string(Base.pkgversion(Breeze)), "Oceananigans" => string(Base.pkgversion(Oceananigans)),
                                     "NumericalEarth" => string(Base.pkgversion(NumericalEarth)), "julia" => string(VERSION)))
open(joinpath(output_dir, "provenance.toml"), "w") do io
    TOML.print(io, provenance)
end

# ## Run

run!(simulation)

# ## A surface snapshot

T_series = FieldTimeSeries(joinpath(output_dir, "outer_slices.jld2"), "T_sfc")
w_series = FieldTimeSeries(joinpath(output_dir, "outer_slices.jld2"), "w_2km")
fig = Figure(size = (1100, 480))
ax1 = Axis(fig[1, 1]; title = "surface T (K), t = $(prettytime(T_series.times[end]))", xlabel = "longitude", ylabel = "latitude")
ax2 = Axis(fig[1, 2]; title = "w at 2 km (m/s)", xlabel = "longitude")
hm1 = heatmap!(ax1, T_series[end]; colormap = :thermal)
hm2 = heatmap!(ax2, w_series[end]; colormap = :balance, colorrange = (-3, 3))
Colorbar(fig[1, 3], hm2)
save(joinpath(output_dir, "outer_surface_snapshot.png"), fig)
fig
