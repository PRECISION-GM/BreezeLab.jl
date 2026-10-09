# Model column at the TRACER AMF1 site (La Porte, 29.670 N 95.059 W) versus ARM observations and the soundings.
#   julia --project=cases/tracer_mip analysis/tracer_mip_site_comparison.jl <run dir> <arm dir> [<out dir>]
# Model: 10-min inner-region 3-D state (T, p, qᵛ, qᶜˡ, qʳ, u, v at cell centers), 10-min surface rain flux.
# Obs: houmetM1.b1 (2 m T/RH, 10 m wind, tipping-bucket + PWD precipitation), houceilM1.b1 (first cloud base),
# houarsclkazr1kolliasM1.c0 (cloud base best estimate), houldM1.b1 (disdrometer rain rate), housondewnpnM1.b1.
using Oceananigans, Breeze, NCDatasets, CairoMakie, Statistics, Printf
using Oceananigans.Fields: interior
using Oceananigans.Grids: λnodes, φnodes, znode, Center, Face
using Dates
using TOML: TOML

run_dir = abspath(ARGS[1]); arm_dir = abspath(ARGS[2]); out_dir = abspath(get(ARGS, 3, joinpath(run_dir, "analysis"))); mkpath(out_dir)
t0 = DateTime(2022, 8, 7, 6); site = (λ = -95.059, φ = 29.670)   # AMF1 La Porte
inner(name) = FieldTimeSeries(joinpath(run_dir, "inner_region", "inner_region_state.jld2"), name; backend = OnDisk())
F = Dict(n => inner(n) for n in ("T", "p", "qᵛ", "qᶜˡ", "qʳ", "u", "v"))
rain_fts = FieldTimeSeries(joinpath(run_dir, "outer_surface.jld2"), "rain"; backend = OnDisk())
grid = F["T"].grid; λ = Array(λnodes(grid, Center())); φ = Array(φnodes(grid, Center()))
i_range, j_range, _ = F["T"].indices
is = argmin(abs.(λ[i_range] .- site.λ)); js = argmin(abs.(φ[j_range] .- site.φ)); I = i_range[is]; J = j_range[js]
Nz = size(grid, 3); zc = [znode(I, J, k, grid, Center(), Center(), Center()) for k in 1:Nz]; z0 = znode(I, J, 1, grid, Center(), Center(), Face())
times = Float64.(collect(F["T"].times)); nt = min(length(times), length(rain_fts.times))
mtime = [t0 + Second(round(Int, t)) for t in times[1:nt]]
col(name, n) = Array(interior(F[name][n]))[is, js, :]
T1 = Float64[]; U1 = Float64[]; D1 = Float64[]; RH1 = Float64[]; CB = Float64[]; RR = Float64[]; profiles = Dict{Int, Dict{String, Vector{Float64}}}()
# Saturation vapour pressure (Bolton) for RH from qᵛ, p, T
es(T) = 611.2 * exp(17.67 * (T - 273.15) / (T - 29.65))
for n in 1:nt
    T = col("T", n); p = col("p", n); qv = col("qᵛ", n); qc = col("qᶜˡ", n); u = col("u", n); v = col("v", n)
    push!(T1, T[1]); push!(U1, hypot(u[1], v[1])); push!(D1, mod(270 - atand(v[1], u[1]), 360))
    e = qv[1] * p[1] / (0.622 + 0.378 * qv[1]); push!(RH1, 100 * e / es(T[1]))
    k = findfirst(>(1e-5), qc); push!(CB, isnothing(k) ? NaN : zc[k] - z0)
    push!(RR, 3600 * Array(interior(rain_fts[n]))[I, J, 1])
    h = round(Int, times[n] / 3600)
    if abs(times[n] - 3600h) < 1 && h in (0, 6, 12)
        e = qv .* p ./ (0.622 .+ 0.378 .* qv); rh = 100 .* e ./ es.(T)
        profiles[h] = Dict("z" => zc .- z0, "T" => T, "rh" => rh, "wspd" => hypot.(u, v), "p" => p)
    end
end

# ARM observations (2022-08-07 and 08-08 files concatenated); times decoded by NCDatasets
function arm(prefix, vars)
    out = Dict{String, Vector}(v => Float64[] for v in vars); out["time"] = DateTime[]
    for f in sort(filter(x -> startswith(x, prefix), readdir(arm_dir)))
        ds = NCDataset(joinpath(arm_dir, f))
        append!(out["time"], DateTime.(ds["time"][:]))
        for v in vars; append!(out[v], Float64.(coalesce.(ds[v][:], NaN))); end
        close(ds)
    end
    return out
end
met = arm("houmetM1.b1", ["temp_mean", "rh_mean", "wspd_vec_mean", "wdir_vec_mean", "atmos_pressure", "tbrg_precip_total", "pwd_precip_rate_mean_1min"])
ceil = arm("houceilM1.b1", ["first_cbh"]); ld = arm("houldM1.b1", ["precip_rate"]); arscl = arm("houarsclkazr1kolliasM1.c0", ["cloud_base_best_estimate"])
window(d) = (d["time"] .>= t0) .& (d["time"] .<= t0 + Hour(24))
sondes = Dict{String, Dict{String, Vector{Float64}}}()
for f in sort(filter(x -> startswith(x, "housondewnpnM1.b1.20220807"), readdir(arm_dir)))
    ds = NCDataset(joinpath(arm_dir, f)); alt = Float64.(ds["alt"][:]); keep = findall(a -> a < 20000, alt)
    sondes[f[22:27]] = Dict("z" => alt[keep] .- alt[1], "T" => Float64.(ds["tdry"][keep]) .+ 273.15, "rh" => Float64.(ds["rh"][keep]), "wspd" => Float64.(ds["wspd"][keep]), "p" => 100 .* Float64.(ds["pres"][keep]))
    close(ds)
end

# Sea-breeze onset: first time after 12 UTC with sustained (≥ 30 min) onshore wind direction (90°–200° at La Porte)
onshore(d) = 90 ≤ d ≤ 200
function onset(times, dirs, after)
    idx = findall(t -> t ≥ after, times)
    for (q, i) in enumerate(idx)
        win = idx[q:min(q + 29, end)]; length(win) < 10 && break
        all(onshore.(dirs[win])) && return times[i]
    end
    return nothing
end
m10 = round.(Int, times[1:nt] ./ 600)
obs_onset = onset(met["time"][window(met)], met["wdir_vec_mean"][window(met)], t0 + Hour(6))
model_onset = let idx = findall(t -> t ≥ t0 + Hour(6), mtime), found = nothing
    for q in eachindex(idx); i = idx[q]; win = idx[q:min(q + 2, end)]; (length(win) == 3 && all(onshore.(D1[win]))) && (found = mtime[i]; break); end; found
end
obs_rain_mm = sum(filter(!isnan, met["tbrg_precip_total"][window(met)])); model_rain_mm = sum(RR) / 6
sel(d, v) = (m = window(d); (d["time"][m], d[v][m]))
hourly_mean(times, vals, h) = (m = [t0 + Hour(h) ≤ t < t0 + Hour(h + 1) for t in times]; v = filter(!isnan, Float64.(vals[m])); isempty(v) ? NaN : mean(v))
hours = 0:floor(Int, times[nt] / 3600)
summary = Dict("site" => "AMF1 La Porte 29.670N 95.059W; model cell ($I, $J) at $(λ[I]), $(φ[J]); lowest level $(round(zc[1]-z0)) m vs met 2 m T / 10 m wind",
    "model_last_time" => string(mtime[end]), "obs_seabreeze_onset" => string(obs_onset), "model_seabreeze_onset" => string(model_onset),
    "obs_rain_mm_24h" => obs_rain_mm, "model_rain_mm_to_last_time" => model_rain_mm,
    "hourly" => Dict(string(h) => Dict("model_T" => hourly_mean(mtime, T1, h), "obs_T2m" => hourly_mean(sel(met, "temp_mean")..., h) + 273.15,
                                         "model_wspd" => hourly_mean(mtime, U1, h), "obs_wspd10m" => hourly_mean(sel(met, "wspd_vec_mean")..., h),
                                         "model_wdir" => hourly_mean(mtime, D1, h), "obs_wdir" => hourly_mean(sel(met, "wdir_vec_mean")..., h),
                                         "model_rh" => hourly_mean(mtime, RH1, h), "obs_rh" => hourly_mean(sel(met, "rh_mean")..., h),
                                         "model_cloud_base_m" => hourly_mean(mtime, CB, h), "obs_ceil_first_cbh_m" => hourly_mean(sel(ceil, "first_cbh")..., h),
                                         "obs_arscl_cbh_m" => hourly_mean(sel(arscl, "cloud_base_best_estimate")..., h)) for h in hours),
    "sonde_files" => collect(keys(sondes)))
open(joinpath(out_dir, "site_comparison.toml"), "w") do io; TOML.print(io, summary); end

hrs(t) = Dates.value.(t .- t0) ./ 3.6e6
fig = Figure(size = (1200, 1000))
ax = Axis(fig[1, 1]; ylabel = "T (K)", title = "La Porte: model lowest level (25 m) vs met 2 m"); lines!(ax, hrs(mtime), T1; label = "model"); lines!(ax, hrs(sel(met, "temp_mean")[1]), sel(met, "temp_mean")[2] .+ 273.15; label = "obs 2 m"); axislegend(ax; position = :rb)
ax = Axis(fig[1, 2]; ylabel = "RH (%)"); lines!(ax, hrs(mtime), RH1); lines!(ax, hrs(sel(met, "rh_mean")...))
ax = Axis(fig[2, 1]; ylabel = "wind speed (m/s)"); lines!(ax, hrs(mtime), U1; label = "model 25 m"); lines!(ax, hrs(sel(met, "wspd_vec_mean")...); label = "obs 10 m"); axislegend(ax)
ax = Axis(fig[2, 2]; ylabel = "wind direction (°)"); scatter!(ax, hrs(mtime), D1; markersize = 4, label = "model"); scatter!(ax, hrs(sel(met, "wdir_vec_mean")...); markersize = 2, label = "obs"); hspan!(ax, 90, 200; color = (:green, 0.1)); axislegend(ax)
ax = Axis(fig[3, 1]; ylabel = "cloud base (m)", xlabel = "hours after 06 UTC", limits = (nothing, (0, 6000))); scatter!(ax, hrs(sel(ceil, "first_cbh")...); markersize = 2, label = "ceilometer"); scatter!(ax, hrs(sel(arscl, "cloud_base_best_estimate")...); markersize = 2, label = "ARSCL"); scatter!(ax, hrs(mtime), CB; markersize = 6, color = :red, label = "model (qᶜˡ > 1e-5)"); axislegend(ax)
ax = Axis(fig[3, 2]; ylabel = "rain rate (mm/h)", xlabel = "hours after 06 UTC"); lines!(ax, hrs(mtime), RR; label = "model"); lines!(ax, hrs(sel(ld, "precip_rate")...); label = "disdrometer"); lines!(ax, hrs(sel(met, "pwd_precip_rate_mean_1min")...); label = "PWD"); axislegend(ax)
save(joinpath(out_dir, "site_timeseries.png"), fig; px_per_unit = 1)

fig2 = Figure(size = (1200, 500))
pairs = [(0, "052900"), (6, "113000"), (12, "173000")]
for (c, (h, sf)) in enumerate(pairs)
    haskey(profiles, h) && haskey(sondes, sf) || continue
    ax = Axis(fig2[1, c]; title = "model $(h) h vs sonde $(sf[1:2]):$(sf[3:4]) UTC", xlabel = "T (K) / wind (m/s)", ylabel = "z AGL (m)", limits = (nothing, (0, 15000)))
    lines!(ax, profiles[h]["T"], profiles[h]["z"]; color = :red, label = "model T"); lines!(ax, sondes[sf]["T"], sondes[sf]["z"]; color = :red, linestyle = :dash, label = "sonde T")
    lines!(ax, 10 .* profiles[h]["wspd"], profiles[h]["z"]; color = :blue, label = "model 10×wind"); lines!(ax, 10 .* sondes[sf]["wspd"], sondes[sf]["z"]; color = :blue, linestyle = :dash, label = "sonde 10×wind")
    c == 1 && axislegend(ax; position = :rt, fontsize = 9)
end
save(joinpath(out_dir, "site_soundings.png"), fig2; px_per_unit = 1)
@info "written" out_dir obs_onset model_onset obs_rain_mm model_rain_mm
