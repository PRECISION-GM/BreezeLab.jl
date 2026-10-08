# Lowest-level wind of a TRACER-MIP outer run versus ERA5 10 m and pressure-level winds.
#   julia --project=cases/tracer_mip analysis/tracer_mip_surface_wind_check.jl <run or snapshot dir> <era5 dir> [<out dir>]
# Writes surface_wind_check.toml and surface_wind_check.png. Reads the run's 10-min outer_slices (u_sfc, v_sfc at the
# lowest model level, 25 m) and inner_region_state (3-D u, v on the inner box), and ERA5 hourly 10 m winds (outer box)
# and pressure-level u, v, geopotential (padded box). ERA5 is the model's own initial/boundary data, so at 06 UTC the
# two must agree up to the 10 m → 25 m log-layer difference; later divergence measures the model's surface-layer physics.
using Oceananigans, Breeze, NCDatasets, CairoMakie, Statistics, Printf
using Oceananigans.Fields: interior
using Oceananigans.Grids: λnodes, φnodes, znode, Center, Face
using Dates: DateTime, Hour, format
using TOML: TOML

run_dir = abspath(ARGS[1]); era5_dir = abspath(ARGS[2]); out_dir = abspath(get(ARGS, 3, joinpath(run_dir, "diagnostics")))
mkpath(out_dir)
case_start = DateTime(2022, 8, 7, 6)
site = (λ = -95.0792, φ = 29.4719)
outer_box = "_-102.8_-87.3_22.7_36.2"; parent_box = "_-103.3_-86.8_22.2_36.7"
stamp(d) = format(d, "yyyy-mm-ddTHH")

u_s = FieldTimeSeries(joinpath(run_dir, "outer_slices.jld2"), "u_sfc"; backend = OnDisk())
v_s = FieldTimeSeries(joinpath(run_dir, "outer_slices.jld2"), "v_sfc"; backend = OnDisk())
grid = u_s.grid; λ = Array(λnodes(grid, Center())); φ = Array(φnodes(grid, Center()))
times = Float64.(collect(u_s.times)); nt = length(times)
inner = Dict(n => FieldTimeSeries(joinpath(run_dir, "inner_region", "inner_region_state.jld2"), n; backend = OnDisk()) for n in ("u", "v"))
i_range, j_range, _ = inner["u"].indices
λi, φi = λ[i_range], φ[j_range]
i_site = argmin(abs.(λi .- site.λ)); j_site = argmin(abs.(φi .- site.φ))
I_site = argmin(abs.(λ .- site.λ)); J_site = argmin(abs.(φ .- site.φ))
Nz = size(grid, 3)
z_site = [znode(i_range[i_site], j_range[j_site], k, grid, Center(), Center(), Center()) for k in 1:Nz]
z_sfc = znode(i_range[i_site], j_range[j_site], 1, grid, Center(), Center(), Face())

centered(a, dim) = dim == 1 ? 0.5 .* (a[1:end-1, :] .+ a[2:end, :]) : 0.5 .* (a[:, 1:end-1] .+ a[:, 2:end])
function speed2d(n)   # u_sfc / v_sfc are face-located; average to cell centers
    u = dropdims(Array(interior(u_s[n])); dims = 3); v = dropdims(Array(interior(v_s[n])); dims = 3)
    size(u, 1) == length(λ) + 1 && (u = centered(u, 1)); size(v, 2) == length(φ) + 1 && (v = centered(v, 2))
    return sqrt.(u .^ 2 .+ v .^ 2)
end
rim = falses(length(λ), length(φ)); rim[1:10, :] .= true; rim[end-9:end, :] .= true; rim[:, 1:10] .= true; rim[:, end-9:end] .= true
model_mean = Float64[]; model_site = Float64[]; model_interior_mean = Float64[]
for n in 1:nt
    s = speed2d(n); push!(model_mean, mean(s)); push!(model_interior_mean, mean(s[.!rim])); push!(model_site, s[I_site, J_site])
end

# ERA5 10 m winds, hourly, over the (unpadded) outer box; nearest-cell at the site.
function ncvar(path, candidates)
    ds = NCDataset(path)
    name = first(filter(c -> haskey(ds, c), candidates))
    lon = Array(ds["longitude"][:]); lat = Array(ds["latitude"][:])
    data = Array(ds[name][:])
    close(ds)
    return lon, lat, dropdims(data; dims = Tuple(findall(==(1), size(data))))
end
hours = 0:floor(Int, times[end] / 3600)
era_mean = Float64[]; era_site = Float64[]; era_hours_s = Float64[]
for h in hours
    d = case_start + Hour(h)
    lon, lat, u10 = ncvar(joinpath(era5_dir, "10m_u_component_of_wind_ERA5HourlySingleLevel_$(stamp(d))$(outer_box).nc"), ["u10"])
    _, _, v10 = ncvar(joinpath(era5_dir, "10m_v_component_of_wind_ERA5HourlySingleLevel_$(stamp(d))$(outer_box).nc"), ["v10"])
    s = sqrt.(u10 .^ 2 .+ v10 .^ 2)
    push!(era_mean, mean(s)); push!(era_hours_s, 3600h)
    push!(era_site, s[argmin(abs.(lon .- site.λ)), argmin(abs.(lat .- site.φ))])
end

# Profiles at the site through the lowest 1.5 km: model (inner-region 3-D state) vs ERA5 pressure levels
profile_hours = filter(h -> 3600h ≤ times[end], [0, 3, 6, 9, 12, 18])
model_profiles = Dict{Int, Vector{Float64}}(); era_profiles = Dict{Int, Tuple{Vector{Float64}, Vector{Float64}}}()
for h in profile_hours
    n = findfirst(t -> abs(t - 3600h) < 1, times); isnothing(n) && continue
    u = Array(interior(inner["u"][n]))[i_site, j_site, :]; v = Array(interior(inner["v"][n]))[i_site, j_site, :]
    model_profiles[h] = sqrt.(u .^ 2 .+ v .^ 2)
    d = case_start + Hour(h)
    lon, lat, up = ncvar(joinpath(era5_dir, "u_component_of_wind_ERA5HourlyPressureLevels_$(stamp(d))$(parent_box).nc"), ["u"])
    _, _, vp = ncvar(joinpath(era5_dir, "v_component_of_wind_ERA5HourlyPressureLevels_$(stamp(d))$(parent_box).nc"), ["v"])
    _, _, zp = ncvar(joinpath(era5_dir, "geopotential_ERA5HourlyPressureLevels_$(stamp(d))$(parent_box).nc"), ["z"])
    i = argmin(abs.(lon .- site.λ)); j = argmin(abs.(lat .- site.φ))
    zs = zp[i, j, :] ./ 9.80665; sp = sqrt.(up[i, j, :] .^ 2 .+ vp[i, j, :] .^ 2)
    keep = findall(z -> -50 < z < 3000, zs)
    era_profiles[h] = (zs[keep], sp[keep])
end

k_1km = findall(z -> z - z_sfc ≤ 1500, z_site)
report = Dict{String, Any}("run_dir" => run_dir, "times_h" => times ./ 3600, "model_domain_mean_lowest_level_speed" => model_mean,
    "model_domain_mean_lowest_level_speed_excluding_rim" => model_interior_mean, "model_site_lowest_level_speed" => model_site,
    "era5_hours" => era_hours_s ./ 3600, "era5_domain_mean_10m_speed" => era_mean, "era5_site_10m_speed" => era_site,
    "lowest_level_height_agl_m" => z_site[1] - z_sfc, "site_terrain_m" => z_sfc,
    "profiles" => Dict(string(h) => Dict("model_z_agl" => z_site[k_1km] .- z_sfc, "model_speed" => model_profiles[h][k_1km],
                                         "era5_z" => era_profiles[h][1], "era5_speed" => era_profiles[h][2]) for h in keys(model_profiles)))
ratio_hours = [h for h in hours if 3600h ≤ times[end]]
report["ratio_model_to_era5_domain_mean_by_hour"] = [model_mean[findfirst(t -> abs(t - 3600h) < 1, times)] / era_mean[h + 1] for h in ratio_hours]
open(joinpath(out_dir, "surface_wind_check.toml"), "w") do io; TOML.print(io, report); end

fig = Figure(size = (1100, 450))
ax = Axis(fig[1, 1]; xlabel = "hours after 06 UTC", ylabel = "wind speed (m/s)", title = "lowest model level (25 m) vs ERA5 10 m")
lines!(ax, times ./ 3600, model_mean; label = "model, domain mean"); lines!(ax, times ./ 3600, model_interior_mean; label = "model, excluding 10-cell rim", linestyle = :dash)
lines!(ax, times ./ 3600, model_site; label = "model, site")
scatterlines!(ax, era_hours_s ./ 3600, era_mean; label = "ERA5 10 m, domain mean", marker = :circle)
scatterlines!(ax, era_hours_s ./ 3600, era_site; label = "ERA5 10 m, site", marker = :utriangle)
axislegend(ax; position = :rt, fontsize = 9)
ax2 = Axis(fig[1, 2]; xlabel = "wind speed (m/s)", ylabel = "height AGL (m)", title = "site profiles (solid model, dashed ERA5)", limits = (nothing, (0, 1500)))
for (k, h) in enumerate(sort(collect(keys(model_profiles))))
    c = Makie.wong_colors()[mod1(k, 7)]
    lines!(ax2, model_profiles[h][k_1km], z_site[k_1km] .- z_sfc; color = c, label = "$(h) h")
    scatterlines!(ax2, era_profiles[h][2], era_profiles[h][1]; color = c, linestyle = :dash, marker = :circle, markersize = 6)
end
axislegend(ax2; position = :rb, fontsize = 9)
save(joinpath(out_dir, "surface_wind_check.png"), fig; px_per_unit = 1)
@info "written" out_dir report["ratio_model_to_era5_domain_mean_by_hour"]
