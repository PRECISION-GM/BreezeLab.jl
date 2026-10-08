# Diagnostics for an outer-domain pilot run directory (CPU; reads the JLD2 output, writes figures + a TOML summary).
#   julia --project=cases/tracer_mip analysis/tracer_mip_outer_pilot.jl runs/<experiment>_job<id>
using Oceananigans, Oceananigans.Units
using Oceananigans.Fields: interior
using CairoMakie, Printf, Statistics
using TOML: TOML
using Dates: now, UTC

dir = abspath(get(ARGS, 1, "."))
state = FieldTimeSeries(joinpath(dir, "outer_state.jld2"), "u"; backend = OnDisk())
times = state.times
grid = state.grid
Nx, Ny, Nz = size(grid)
zc = Array(znodes(grid, Center()))
summary = Dict{String, Any}("directory" => dir, "generated_utc" => string(now(UTC)), "times_s" => collect(times), "grid" => [Nx, Ny, Nz])

function load(name, n)
    fts = FieldTimeSeries(joinpath(dir, "outer_state.jld2"), name; backend = OnDisk())
    return Array(interior(fts[n]))
end

for (label, n) in (("initial", 1), ("final", length(times)))
    u = load("u", n); v = load("v", n); w = load("w", n); T = load("T", n); p = load("p", n); ρ = load("ρ", n)
    qᶜˡ = load("qᶜˡ", n); qʳ = load("qʳ", n); nᶜˡ = load("nᶜˡ", n)
    speed = sqrt.(u[1:Nx, :, :] .^ 2 .+ v[:, 1:Ny, :] .^ 2)
    imax = argmax(abs.(u)); jmax = argmax(abs.(w))
    # boundary-zone (outermost 10 cells) versus interior wind extremes
    rim = falses(Nx, Ny); rim[1:10, :] .= true; rim[end-9:end, :] .= true; rim[:, 1:10] .= true; rim[:, end-9:end] .= true
    rim3 = repeat(rim, 1, 1, Nz)
    # hydrostatic residual: |∂p/∂z + ρ g| / (ρ g), domain median over the lower 10 km
    dpdz = (p[:, :, 2:end] .- p[:, :, 1:end-1]) ./ reshape(diff(zc), 1, 1, :)
    ρm = 0.5 .* (ρ[:, :, 2:end] .+ ρ[:, :, 1:end-1])
    kk = findall(z -> z < 10_000, 0.5 .* (zc[1:end-1] .+ zc[2:end]))
    residual = abs.(dpdz[:, :, kk] .+ 9.81 .* ρm[:, :, kk]) ./ (9.81 .* ρm[:, :, kk])
    summary[label] = Dict(
        "time_s" => times[n],
        "max_abs_u" => maximum(abs, u), "max_abs_u_index" => collect(Tuple(imax)), "max_abs_u_height_m" => zc[min(imax[3], Nz)],
        "max_abs_w" => maximum(abs, w), "max_abs_w_index" => collect(Tuple(jmax)), "max_abs_w_height_m" => zc[min(jmax[3], Nz)],
        "max_speed_rim" => maximum(speed[rim3[1:Nx, 1:Ny, :]]), "max_speed_interior" => maximum(speed[.!rim3[1:Nx, 1:Ny, :]]),
        "T_range" => [minimum(T), maximum(T)], "p_range" => [minimum(p), maximum(p)], "ρ_range" => [minimum(ρ), maximum(ρ)],
        "hydrostatic_residual_median" => median(residual), "hydrostatic_residual_p99" => quantile(vec(residual), 0.99),
        "max_qcl" => maximum(qᶜˡ), "max_qr" => maximum(qʳ), "max_ncl" => maximum(nᶜˡ),
        "cloud_fraction_columns" => mean(maximum(qᶜˡ; dims = 3) .> 1e-5), "finite" => all(isfinite, T) && all(isfinite, p))
end

surface = FieldTimeSeries(joinpath(dir, "outer_surface.jld2"), "T_land"); sea = FieldTimeSeries(joinpath(dir, "outer_surface.jld2"), "sea")
rain = FieldTimeSeries(joinpath(dir, "outer_surface.jld2"), "rain")
mask = Array(interior(sea[1]))[:, :, 1] .== 1
Tl = [Array(interior(surface[n]))[:, :, 1] for n in 1:length(surface.times)]
summary["land"] = Dict("times_s" => collect(surface.times), "sea_fraction" => mean(mask),
                       "land_T_mean" => [mean(T[.!mask]) for T in Tl], "sea_T_mean" => [mean(T[mask]) for T in Tl],
                       "max_rain_flux" => [maximum(Array(interior(rain[n]))) for n in 1:length(rain.times)])

procs = FieldTimeSeries(joinpath(dir, "outer_process_rates.jld2"), "liquid_condensation"; backend = OnDisk())
summary["process_rates"] = Dict(name => maximum(Array(interior(FieldTimeSeries(joinpath(dir, "outer_process_rates.jld2"), name; backend = OnDisk())[end])))
                                for name in ("liquid_condensation", "liquid_evaporation", "droplet_nucleation", "autoconversion_accretion", "latent_heating", "ice_deposition", "freezing"))

open(joinpath(dir, "pilot_summary.toml"), "w") do io
    TOML.print(io, summary)
end

# Figures: surface wind speed and temperature, w at 2 km, cloud water at 2 km (final time), land temperature change.
slices = joinpath(dir, "outer_slices.jld2")
u_s = FieldTimeSeries(slices, "u_sfc"); v_s = FieldTimeSeries(slices, "v_sfc"); T_s = FieldTimeSeries(slices, "T_sfc")
w2 = FieldTimeSeries(slices, "w_2km"); q2 = FieldTimeSeries(slices, "qᶜˡ_2km")
n = length(u_s.times)
fig = Figure(size = (1400, 900))
ax = Axis(fig[1, 1]; title = @sprintf("surface T (K), t = %s", prettytime(T_s.times[n])), xlabel = "lon", ylabel = "lat"); hm = heatmap!(ax, T_s[n]; colormap = :thermal); Colorbar(fig[1, 2], hm)
ax = Axis(fig[1, 3]; title = "surface |U| (m/s)"); hm = heatmap!(ax, Field(sqrt(u_s[n]^2 + v_s[n]^2)); colormap = :speed); Colorbar(fig[1, 4], hm)
ax = Axis(fig[2, 1]; title = "w at 2 km (m/s)"); hm = heatmap!(ax, w2[n]; colormap = :balance, colorrange = (-3, 3)); Colorbar(fig[2, 2], hm)
ax = Axis(fig[2, 3]; title = "cloud water at 2 km (g/kg)"); hm = heatmap!(ax, Field(1e3 * q2[n]); colormap = :dense); Colorbar(fig[2, 4], hm)
save(joinpath(dir, "pilot_overview.png"), fig)
fig2 = Figure(size = (700, 500))
ax = Axis(fig2[1, 1]; title = "land skin temperature change (K)"); hm = heatmap!(ax, Field(surface[end] - surface[1]); colormap = :balance, colorrange = (-5, 5)); Colorbar(fig2[1, 2], hm)
save(joinpath(dir, "pilot_land_response.png"), fig2)
println("wrote ", joinpath(dir, "pilot_summary.toml"))
