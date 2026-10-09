# Where and when the August 7 no-closure control (job 233) blew up. Reads the hourly 3-D state, the 10-min
# slices/surface fields and the 10-min inner-region 3-D state; writes blowup_diagnosis.toml and figures.
#   julia --project=cases/tracer_mip analysis/tracer_mip_blowup_diagnosis.jl <run dir> [<out dir>]
using Oceananigans, Breeze, CairoMakie, Statistics, Printf
using Oceananigans.Fields: interior
using Oceananigans.Grids: λnodes, φnodes, znode, Center, Face
using TOML: TOML

run_dir = abspath(ARGS[1]); out_dir = abspath(get(ARGS, 2, joinpath(run_dir, "analysis"))); mkpath(out_dir)
Δt = 3.0; rim_cells = 10
state(name) = FieldTimeSeries(joinpath(run_dir, "outer_state.jld2"), name; backend = OnDisk())
slice(name) = FieldTimeSeries(joinpath(run_dir, "outer_slices.jld2"), name; backend = OnDisk())
inner(name) = FieldTimeSeries(joinpath(run_dir, "inner_region", "inner_region_state.jld2"), name; backend = OnDisk())
sea = dropdims(Array(interior(FieldTimeSeries(joinpath(run_dir, "outer_surface.jld2"), "sea"; backend = OnDisk())[1])); dims = 3) .> 0.5
w_fts = state("w"); T_fts = state("T"); u_fts = state("u")
grid = w_fts.grid; Nx, Ny, Nz = size(grid)
λ = Array(λnodes(grid, Center())); φ = Array(φnodes(grid, Center()))
zc = [znode(Nx ÷ 2, Ny ÷ 2, k, grid, Center(), Center(), Center()) for k in 1:Nz]
zf = [znode(Nx ÷ 2, Ny ÷ 2, k, grid, Center(), Center(), Face()) for k in 1:Nz+1]
dz = diff(zf); Lz = zf[end]
in_rim(i, j) = i ≤ rim_cells || i > Nx - rim_cells || j ≤ rim_cells || j > Ny - rim_cells
in_sponge(k) = zc[k] > 0.75Lz
describe(i, j, k) = Dict("i" => i, "j" => j, "k" => k, "lon" => λ[i], "lat" => φ[j], "z_m" => zc[k],
                         "rim" => in_rim(i, j), "sponge" => in_sponge(k), "sea" => sea[min(i, Nx), min(j, Ny)])

hours = Float64.(collect(w_fts.times)) ./ 3600
hourly = Dict{String, Any}[]
per_level_maxw = zeros(length(hours), Nz)
for n in eachindex(hours)
    w = Array(interior(w_fts[n])); T = Array(interior(T_fts[n])); u = Array(interior(u_fts[n]))
    aw = abs.(w); iw = argmax(aw); iT = argmin(T); iu = argmax(abs.(u))
    for k in 1:Nz; per_level_maxw[n, k] = maximum(@view aw[:, :, min(k, size(aw, 3))]); end
    big = findall(>(10), aw)   # cells with |w| > 10 m/s
    cfl_v = maximum(maximum(@view aw[:, :, k]) / dz[min(k, Nz)] * Δt for k in 1:size(aw, 3))
    push!(hourly, Dict("hour" => hours[n], "max_abs_w" => maximum(aw), "max_w_at" => describe(Tuple(iw)...),
                       "min_T" => minimum(T), "min_T_at" => describe(Tuple(iT)...), "max_T" => maximum(T),
                       "max_abs_u" => maximum(abs.(u)), "max_u_at" => describe(min(Tuple(iu)[1], Nx), Tuple(iu)[2], Tuple(iu)[3]),
                       "n_cells_w_gt_10" => length(big),
                       "frac_w_gt_10_in_rim" => isempty(big) ? 0.0 : mean(in_rim(c[1], c[2]) for c in big),
                       "frac_w_gt_10_in_sponge" => isempty(big) ? 0.0 : mean(in_sponge(min(c[3], Nz)) for c in big),
                       "frac_w_gt_10_over_sea" => isempty(big) ? 0.0 : mean(sea[c[1], c[2]] for c in big),
                       "max_vertical_advective_cfl" => cfl_v, "max_horizontal_cfl" => maximum(abs.(u)) * Δt / 2000))
    @info @sprintf("h=%4.1f max|w|=%6.2f at z=%6.0f m (rim=%s sponge=%s sea=%s) minT=%6.1f at z=%6.0f  max|u|=%5.1f  CFLv=%5.2f  n(|w|>10)=%d",
                   hours[n], maximum(aw), zc[min(Tuple(iw)[3], Nz)], in_rim(Tuple(iw)[1], Tuple(iw)[2]), in_sponge(min(Tuple(iw)[3], Nz)), sea[Tuple(iw)[1], Tuple(iw)[2]],
                   minimum(T), zc[Tuple(iT)[3]], maximum(abs.(u)), cfl_v, length(big))
end

# 10-min: 2-km slices and inner-region 3-D
w2 = slice("w_2km"); Ts = slice("T_sfc"); wi = inner("w"); Ti = inner("T")
t10 = Float64.(collect(w2.times)) ./ 3600
i_range, j_range, _ = wi.indices
ten = Dict{String, Vector{Float64}}("hour" => t10, "max_abs_w_2km" => Float64[], "min_T_sfc" => Float64[], "max_T_sfc" => Float64[],
                                     "max_abs_w_2km_rim" => Float64[], "max_abs_w_2km_interior" => Float64[],
                                     "inner_max_abs_w" => Float64[], "inner_max_w_height" => Float64[], "inner_min_T" => Float64[])
rim2 = [in_rim(i, j) for i in 1:Nx, j in 1:Ny]
for n in eachindex(t10)
    w = dropdims(Array(interior(w2[n])); dims = 3); T = dropdims(Array(interior(Ts[n])); dims = 3)
    push!(ten["max_abs_w_2km"], maximum(abs.(w))); push!(ten["max_abs_w_2km_rim"], maximum(abs.(w[rim2]))); push!(ten["max_abs_w_2km_interior"], maximum(abs.(w[.!rim2])))
    push!(ten["min_T_sfc"], minimum(T)); push!(ten["max_T_sfc"], maximum(T))
    if n ≤ length(wi.times)
        wI = abs.(Array(interior(wi[n]))); iw = argmax(wI)
        push!(ten["inner_max_abs_w"], maximum(wI)); push!(ten["inner_max_w_height"], zc[min(Tuple(iw)[3], Nz)]); push!(ten["inner_min_T"], minimum(Array(interior(Ti[n]))))
    end
end

report = Dict("run_dir" => run_dir, "dt" => Δt, "Lz" => Lz, "sponge_base_m" => 0.75Lz, "rim_cells" => rim_cells, "hourly" => hourly, "ten_minute" => ten,
              "per_level_max_abs_w" => Dict("z_m" => zc, "hours" => hours, "max_abs_w" => [per_level_maxw[n, :] for n in eachindex(hours)]))
open(joinpath(out_dir, "blowup_diagnosis.toml"), "w") do io; TOML.print(io, report); end

fig = Figure(size = (1200, 800))
ax1 = Axis(fig[1, 1]; xlabel = "hours after 06 UTC", ylabel = "max |w| (m/s)", title = "extremes vs time", yscale = log10)
lines!(ax1, t10, max.(ten["max_abs_w_2km"], 1e-2); label = "2 km slice, whole domain"); lines!(ax1, t10, max.(ten["max_abs_w_2km_rim"], 1e-2); label = "2 km slice, rim"; linestyle = :dash)
lines!(ax1, t10[1:length(ten["inner_max_abs_w"])], max.(ten["inner_max_abs_w"], 1e-2); label = "inner box 3-D")
scatterlines!(ax1, hours, [h["max_abs_w"] for h in hourly]; label = "3-D hourly", marker = :circle); axislegend(ax1; position = :lt)
ax2 = Axis(fig[1, 2]; xlabel = "hours after 06 UTC", ylabel = "T (K)", title = "min T (3-D hourly) and surface T range")
scatterlines!(ax2, hours, [h["min_T"] for h in hourly]; label = "min T 3-D"); lines!(ax2, t10, ten["min_T_sfc"]; label = "min T_sfc"); lines!(ax2, t10, ten["max_T_sfc"]; label = "max T_sfc"); axislegend(ax2; position = :lb)
ax3 = Axis(fig[2, 1]; xlabel = "max |w| per level (m/s)", ylabel = "z (km)", title = "hourly profiles of max |w|", xscale = log10)
for n in eachindex(hours); lines!(ax3, max.(per_level_maxw[n, :], 1e-2), zc ./ 1e3; color = n, colorrange = (1, length(hours)), colormap = :viridis); end
hlines!(ax3, [0.75Lz / 1e3]; color = :red, linestyle = :dash); Colorbar(fig[2, 2], limits = (hours[1], hours[end]), colormap = :viridis, label = "hour")
ax4 = Axis(fig[2, 3]; xlabel = "hours", ylabel = "CFL", title = "max vertical advective CFL |w|Δt/Δz (hourly)")
scatterlines!(ax4, hours, [h["max_vertical_advective_cfl"] for h in hourly]); hlines!(ax4, [1.0]; color = :red)
save(joinpath(out_dir, "blowup_diagnosis.png"), fig; px_per_unit = 1)
# map of where |w| > 10 at the last two hours
fig2 = Figure(size = (1000, 450))
for (p, n) in enumerate(max(1, length(hours) - 1):length(hours))
    w = Array(interior(w_fts[n])); colmax = dropdims(maximum(abs.(w); dims = 3); dims = 3)
    ax = Axis(fig2[1, p]; title = @sprintf("column max |w| (m/s), hour %.0f", hours[n]), aspect = DataAspect())
    hm = heatmap!(ax, λ, φ, colmax; colormap = :magma, colorrange = (0, 20)); contour!(ax, λ, φ, Float32.(sea); levels = [0.5], color = :cyan); Colorbar(fig2[1, p + 2], hm)
end
save(joinpath(out_dir, "blowup_map.png"), fig2; px_per_unit = 1)
@info "written" out_dir
