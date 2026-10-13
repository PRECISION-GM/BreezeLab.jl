# TRACER-MIP August 7 outer controls side by side: model runs vs ERA5 vs the TRACER AMF1 (La Porte) observations.
#   julia --project=cases/tracer_mip analysis/tracer_mip_control_comparison.jl <out dir> <era5 dir> <arm dir> <label=run dir> [<label=run dir> ...]
# Per run (10-min or 2.5-min records): La Porte column from the inner-region 3-D state (lowest level ≈ 25 m: T, wind,
# RH, onshore component, cloud base, column cloud water) and surface rain; domain-mean lowest-level wind from the slices;
# a Hovmöller of the onshore lowest-level wind along the coast-normal line through the site; 2-km max|w| (stability).
# Sea-breeze onset (all sources, same rule): first time ≥ 12 UTC with the onshore component toward 315° (NW, the coast
# normal at La Porte) ≥ 1.5 m s⁻¹ sustained for 30 min. Writes control_comparison.toml and figures into <out dir>.
using Oceananigans, Breeze, NCDatasets, CairoMakie, Statistics, Printf
using Oceananigans.Fields: interior
using Oceananigans.Grids: λnodes, φnodes, znode, Center, Face
using Dates
using TOML: TOML

out_dir = abspath(ARGS[1]); era5_dir = abspath(ARGS[2]); arm_dir = abspath(ARGS[3]); mkpath(out_dir)
runs = [(label = split(a, "=")[1], dir = abspath(split(a, "=")[2])) for a in ARGS[4:end]]
t0 = DateTime(2022, 8, 7, 6); t_end = t0 + Hour(24)
site = (λ = -95.059, φ = 29.670)                      # AMF1 La Porte
normal = (-sind(45), cosd(45))                        # unit vector toward 315° (onshore at La Porte)
onshore(u, v) = normal[1] * u + normal[2] * v
es(T) = 611.2 * exp(17.67 * (T - 273.15) / (T - 29.65))
hrs(t) = Dates.value.(t .- t0) ./ 3.6e6

function onset(times, comp; after = t0 + Hour(6), threshold = 1.5, window = Minute(30))
    for (i, t) in enumerate(times)
        t < after && continue
        w = findall(s -> t ≤ s ≤ t + window, times)
        times[w[end]] - t < window - Minute(3) && return nothing     # not enough data left
        all(comp[w] .≥ threshold) && return t
    end
    return nothing
end

#####
##### Observations (ARM) and ERA5
#####

function arm(prefix, vars)
    out = Dict{String, Vector}(v => Float64[] for v in vars); out["time"] = DateTime[]
    for f in sort(filter(x -> startswith(x, prefix), readdir(arm_dir)))
        NCDataset(joinpath(arm_dir, f)) do ds
            append!(out["time"], DateTime.(ds["time"][:]))
            for v in vars; append!(out[v], Float64.(coalesce.(ds[v][:], NaN))); end
        end
    end
    m = (out["time"] .>= t0) .& (out["time"] .<= t_end)
    return Dict(k => v[m] for (k, v) in out)
end
met = arm("houmetM1.b1", ["temp_mean", "rh_mean", "wspd_vec_mean", "wdir_vec_mean", "tbrg_precip_total"])
met_u = -met["wspd_vec_mean"] .* sind.(met["wdir_vec_mean"]); met_v = -met["wspd_vec_mean"] .* cosd.(met["wdir_vec_mean"])
met_c = onshore.(met_u, met_v)
ceil = arm("houceilM1.b1", ["first_cbh"])
obs_onset = onset(met["time"], met_c)

function era5_series(name, i_site = nothing)
    out = Float64[]; site_out = Float64[]; times = DateTime[]
    for h in 0:24
        d = t0 + Hour(h); stamp = Dates.format(d, "yyyy-mm-ddTHH")
        f = joinpath(era5_dir, "$(name)_ERA5HourlySingleLevel_$(stamp)_-102.8_-87.3_22.7_36.2.nc")
        isfile(f) || continue
        NCDataset(f) do ds
            key = first(filter(k -> !(k in ("longitude", "latitude", "valid_time", "time", "number", "expver")), keys(ds)))
            a = dropdims(Array(ds[key]); dims = Tuple(findall(==(1), size(Array(ds[key])))))
            lon = Array(ds["longitude"][:]); lat = Array(ds["latitude"][:])
            push!(out, NaN); push!(site_out, a[argmin(abs.(lon .- site.λ)), argmin(abs.(lat .- site.φ))])
            out[end] = mean(a); push!(times, d)
        end
    end
    return times, out, site_out
end
e_t, e_u10_mean, e_u10_site = era5_series("10m_u_component_of_wind")
_, e_v10_mean, e_v10_site = era5_series("10m_v_component_of_wind")
_, e_t2_mean, e_t2_site = era5_series("2m_temperature")
era5_speed_site = hypot.(e_u10_site, e_v10_site); era5_c_site = onshore.(e_u10_site, e_v10_site)
# Domain-mean 10 m wind SPEED (mean of |U|, comparable to the model's domain-mean lowest-level speed)
era5_mean_speed = map(e_t) do d
    stamp = Dates.format(d, "yyyy-mm-ddTHH"); box = "_-102.8_-87.3_22.7_36.2.nc"
    u = NCDataset(ds -> Array(ds["u10"]), joinpath(era5_dir, "10m_u_component_of_wind_ERA5HourlySingleLevel_" * stamp * box))
    v = NCDataset(ds -> Array(ds["v10"]), joinpath(era5_dir, "10m_v_component_of_wind_ERA5HourlySingleLevel_" * stamp * box))
    mean(hypot.(u, v))
end
era5_onset = onset(e_t, era5_c_site; window = Minute(60))

#####
##### Model runs
#####

struct RunSeries
    label; time; T; U; dir; c; rh; cbh; lwp; rain; dom_U; dom_time; w2km; w2km_interior; hov; hov_s; hov_time; onset; last
end

function load_run(label, dir)
    inner(n) = FieldTimeSeries(joinpath(dir, "inner_region", "inner_region_state.jld2"), n; backend = OnDisk())
    F = Dict(n => inner(n) for n in ("T", "p", "qᵛ", "qᶜˡ", "qʳ", "u", "v", "ρ"))
    grid = F["T"].grid; λ = Array(λnodes(grid, Center())); φ = Array(φnodes(grid, Center()))
    i_range, j_range, _ = F["T"].indices
    is = argmin(abs.(λ[i_range] .- site.λ)); js = argmin(abs.(φ[j_range] .- site.φ)); I = i_range[is]; J = j_range[js]
    Nz = size(grid, 3); zc = [znode(I, J, k, grid, Center(), Center(), Center()) for k in 1:Nz]
    zf = [znode(I, J, k, grid, Center(), Center(), Face()) for k in 1:Nz+1]; z0 = zf[1]; dz = diff(zf)
    surf = FieldTimeSeries(joinpath(dir, "outer_surface.jld2"), "rain"; backend = OnDisk())
    nt = min(length(F["T"].times), length(surf.times))
    times = [t0 + Millisecond(round(Int, 1000 * t)) for t in Float64.(collect(F["T"].times))[1:nt]]
    col(name, n) = Array(interior(F[name][n]))[is, js, :]
    T1 = Float64[]; U1 = Float64[]; D1 = Float64[]; C1 = Float64[]; RH1 = Float64[]; CB = Float64[]; LWP = Float64[]; RR = Float64[]
    for n in 1:nt
        T = col("T", n); p = col("p", n); qv = col("qᵛ", n); qc = col("qᶜˡ", n); u = col("u", n); v = col("v", n); ρ = col("ρ", n)
        push!(T1, T[1]); push!(U1, hypot(u[1], v[1])); push!(D1, mod(270 - atand(v[1], u[1]), 360)); push!(C1, onshore(u[1], v[1]))
        e = qv[1] * p[1] / (0.622 + 0.378 * qv[1]); push!(RH1, 100 * e / es(T[1]))
        k = findfirst(>(1e-5), qc); push!(CB, isnothing(k) ? NaN : zc[k] - z0)
        push!(LWP, 1e3 * sum(ρ .* qc .* dz)); push!(RR, 3600 * Array(interior(surf[n]))[I, J, 1])
    end
    # Slices: domain-mean lowest-level wind, 2-km max|w|, onshore Hovmöller along the coast normal through the site
    sl(n) = FieldTimeSeries(joinpath(dir, "outer_slices.jld2"), n; backend = OnDisk())
    us = sl("u_sfc"); vs = sl("v_sfc"); w2 = sl("w_2km")
    Nx, Ny, _ = size(us.grid); λo = Array(λnodes(us.grid, Center())); φo = Array(φnodes(us.grid, Center()))
    rim = [i ≤ 10 || i > Nx - 10 || j ≤ 10 || j > Ny - 10 for i in 1:Nx, j in 1:Ny]
    km_per_deg = π * 6371 / 180
    s_points = collect(-150:4:150)                                       # km along the normal (positive = inland, NW)
    hov_ij = [(argmin(abs.(λo .- (site.λ + s * normal[1] / (km_per_deg * cosd(site.φ))))), argmin(abs.(φo .- (site.φ + s * normal[2] / km_per_deg)))) for s in s_points]
    ns = length(us.times); stimes = [t0 + Millisecond(round(Int, 1000 * t)) for t in Float64.(collect(us.times))]
    domU = Float64[]; wmax = Float64[]; wint = Float64[]; hov = zeros(ns, length(s_points))
    center(a, d) = d == 1 ? 0.5 .* (a[1:end-1, :] .+ a[2:end, :]) : 0.5 .* (a[:, 1:end-1] .+ a[:, 2:end])
    for n in 1:ns
        u = dropdims(Array(interior(us[n])); dims = 3); v = dropdims(Array(interior(vs[n])); dims = 3)
        size(u, 1) == Nx + 1 && (u = center(u, 1)); size(v, 2) == Ny + 1 && (v = center(v, 2))
        push!(domU, mean(hypot.(u, v)))
        w = abs.(dropdims(Array(interior(w2[n])); dims = 3)); push!(wmax, maximum(w)); push!(wint, maximum(w[.!rim]))
        hov[n, :] = [onshore(u[i, j], v[i, j]) for (i, j) in hov_ij]
    end
    model_onset = onset(times, C1)
    return RunSeries(label, times, T1, U1, D1, C1, RH1, CB, LWP, RR, domU, stimes, wmax, wint, hov, s_points, stimes, model_onset, times[end])
end
series = [load_run(r.label, r.dir) for r in runs]

#####
##### Statistics against observations (hourly means, common hours)
#####

hourly(times, vals, h) = (m = [t0 + Hour(h) ≤ t < t0 + Hour(h + 1) for t in times]; v = filter(!isnan, Float64.(vals[m])); isempty(v) ? NaN : mean(v))
function skill(s::RunSeries)
    hs = [h for h in 0:23 if t0 + Hour(h + 1) ≤ s.last]
    dT = [hourly(s.time, s.T, h) - (hourly(met["time"], met["temp_mean"], h) + 273.15) for h in hs]
    dU = [hourly(s.time, s.U, h) - hourly(met["time"], met["wspd_vec_mean"], h) for h in hs]
    dRH = [hourly(s.time, s.rh, h) - hourly(met["time"], met["rh_mean"], h) for h in hs]
    ok(x) = filter(!isnan, x)
    return Dict("hours_compared" => length(hs), "T_bias_K" => mean(ok(dT)), "T_rmse_K" => sqrt(mean(ok(dT) .^ 2)),
                "wind_bias_m_s" => mean(ok(dU)), "wind_rmse_m_s" => sqrt(mean(ok(dU) .^ 2)),
                "RH_bias_pct" => mean(ok(dRH)), "RH_rmse_pct" => sqrt(mean(ok(dRH) .^ 2)),
                "daytime_T_max_K" => maximum(s.T[[t0 + Hour(4) ≤ t ≤ t0 + Hour(16) for t in s.time]]; init = NaN),
                "seabreeze_onset" => string(s.onset), "last_time" => string(s.last),
                "rain_at_site_mm" => sum(s.rain[1:end-1] .* diff(hrs(s.time))),
                "max_w_2km" => maximum(s.w2km), "max_w_2km_interior" => maximum(s.w2km_interior),
                "domain_mean_lowest_wind_hourly" => [hourly(s.dom_time, s.dom_U, h) for h in hs],
                "era5_domain_mean_10m_wind_hourly" => [era5_mean_speed[k] for k in eachindex(e_t) if hrs([e_t[k]])[1] in hs])
end
obs_daytime_Tmax = maximum(met["temp_mean"][[t0 + Hour(4) ≤ t ≤ t0 + Hour(16) for t in met["time"]]]) + 273.15
report = Dict("observations" => Dict("seabreeze_onset" => string(obs_onset), "daytime_T2m_max_K" => obs_daytime_Tmax,
                                     "rain_mm" => sum(filter(!isnan, met["tbrg_precip_total"]))),
              "era5" => Dict("seabreeze_onset_hourly" => string(era5_onset)),
              "runs" => Dict(s.label => skill(s) for s in series),
              "onset_rule" => "first time ≥ 12 UTC with onshore component toward 315° ≥ 1.5 m/s sustained 30 min (ERA5 hourly: 60 min)")
open(joinpath(out_dir, "control_comparison.toml"), "w") do io; TOML.print(io, report); end

#####
##### Figures
#####

colors = Makie.wong_colors()
fig = Figure(size = (1300, 1100))
axT = Axis(fig[1, 1]; ylabel = "T (K)", title = "La Porte: model lowest level (25 m), AMF1 2 m / 10 m, ERA5")
axU = Axis(fig[1, 2]; ylabel = "wind speed (m/s)")
axC = Axis(fig[2, 1]; ylabel = "onshore wind toward 315° (m/s)")
axR = Axis(fig[2, 2]; ylabel = "RH (%)")
axB = Axis(fig[3, 1]; ylabel = "cloud base (m)", xlabel = "hours after 06 UTC", limits = (nothing, (0, 4000)))
axP = Axis(fig[3, 2]; ylabel = "rain rate at site (mm/h)", xlabel = "hours after 06 UTC")
lines!(axT, hrs(met["time"]), met["temp_mean"] .+ 273.15; color = :black, label = "AMF1")
lines!(axU, hrs(met["time"]), met["wspd_vec_mean"]; color = :black, label = "AMF1")
lines!(axC, hrs(met["time"]), met_c; color = :black, label = "AMF1")
lines!(axR, hrs(met["time"]), met["rh_mean"]; color = :black)
scatter!(axB, hrs(ceil["time"]), ceil["first_cbh"]; color = :black, markersize = 2, label = "ceilometer")
scatterlines!(axT, hrs(e_t), e_t2_site; color = :grey, marker = :diamond, label = "ERA5 2 m")
scatterlines!(axU, hrs(e_t), era5_speed_site; color = :grey, marker = :diamond, label = "ERA5 10 m")
scatterlines!(axC, hrs(e_t), era5_c_site; color = :grey, marker = :diamond, label = "ERA5 10 m")
for (k, s) in enumerate(series)
    c = colors[k]
    lines!(axT, hrs(s.time), s.T; color = c, label = s.label); lines!(axU, hrs(s.time), s.U; color = c, label = s.label)
    lines!(axC, hrs(s.time), s.c; color = c, label = s.label); lines!(axR, hrs(s.time), s.rh; color = c)
    scatter!(axB, hrs(s.time), s.cbh; color = c, markersize = 4, label = s.label); lines!(axP, hrs(s.time), s.rain; color = c, label = s.label)
    isnothing(s.onset) || vlines!(axC, hrs([s.onset]); color = c, linestyle = :dash)
end
isnothing(obs_onset) || vlines!(axC, hrs([obs_onset]); color = :black, linestyle = :dash)
hlines!(axC, [1.5]; color = :grey, linestyle = :dot)
for ax in (axT, axU, axC); axislegend(ax; position = :lt, labelsize = 9); end
axislegend(axB; position = :rt, labelsize = 9)
save(joinpath(out_dir, "site_comparison_runs.png"), fig; px_per_unit = 1)

fig2 = Figure(size = (450 * length(series) + 100, 500))
for (k, s) in enumerate(series)
    ax = Axis(fig2[1, k]; title = "$(s.label): onshore lowest-level wind (m/s)", xlabel = "distance along 315° from La Porte (km; + inland)",
              ylabel = k == 1 ? "hours after 06 UTC" : "", limits = (nothing, (0, 24)))
    hm = heatmap!(ax, s.hov_s, hrs(s.hov_time), permutedims(s.hov); colormap = :balance, colorrange = (-8, 8))
    vlines!(ax, [0]; color = :black, linestyle = :dot)
    isnothing(obs_onset) || hlines!(ax, hrs([obs_onset]); color = :black, linestyle = :dash)
    k == length(series) && Colorbar(fig2[1, k + 1], hm)
end
save(joinpath(out_dir, "seabreeze_hovmoller.png"), fig2; px_per_unit = 1)

fig3 = Figure(size = (1100, 420))
ax1 = Axis(fig3[1, 1]; xlabel = "hours after 06 UTC", ylabel = "max |w| at 2 km (m/s)", yscale = log10, title = "stability: 2-km max |w|")
ax2 = Axis(fig3[1, 2]; xlabel = "hours after 06 UTC", ylabel = "domain-mean lowest-level wind (m/s)", title = "domain-mean wind vs ERA5 10 m")
for (k, s) in enumerate(series)
    lines!(ax1, hrs(s.dom_time), max.(s.w2km, 1e-2); color = colors[k], label = s.label)
    lines!(ax1, hrs(s.dom_time), max.(s.w2km_interior, 1e-2); color = colors[k], linestyle = :dash)
    lines!(ax2, hrs(s.dom_time), s.dom_U; color = colors[k], label = s.label)
end
scatterlines!(ax2, hrs(e_t), era5_mean_speed; color = :grey, marker = :diamond, label = "ERA5 10 m, domain-mean speed")
axislegend(ax1; position = :lt); axislegend(ax2; position = :rt, labelsize = 9)
save(joinpath(out_dir, "stability_and_domain_wind.png"), fig3; px_per_unit = 1)
@info "written" out_dir report["observations"] report["era5"] Dict(k => (v["seabreeze_onset"], v["T_rmse_K"], v["wind_rmse_m_s"]) for (k, v) in report["runs"])
