# Compare one BreezeLab ENA run with the staged ARM ENA observations of its day over the
# common UTC window (plan step 7, observation side). Usage from the repository root:
#
#   julia --project analysis/compare_ena_observations.jl RUN_DIR OBS_DIR [OUTPUT_DIR]
#
# RUN_DIR holds provenance.toml and the *_timeseries.jld2 / *_profiles.jld2 writers;
# OBS_DIR holds the ARM files of that day (enamwrret2turnC1.c1 LWP, enavdisC1.b1 rain,
# enaarsclkazrbnd1kolliasC1.c0 cloud boundaries). The run's protocol label is carried into
# every output: a Covert-public-bin run is reported as the development benchmark it is,
# never as a LASSO result. Writes comparison.toml and comparison.png to OUTPUT_DIR
# (default: RUN_DIR/comparison).
using BreezeLab
using CairoMakie
using Dates
using TOML
using Statistics

length(ARGS) ≥ 2 || error("Usage: julia --project analysis/compare_ena_observations.jl RUN_DIR OBS_DIR [OUTPUT_DIR]")
run_dir, obs_dir = ARGS[1], ARGS[2]
output_dir = length(ARGS) ≥ 3 ? ARGS[3] : joinpath(run_dir, "comparison")
mkpath(output_dir)

find_file(dir, prefix) = (files = filter(f -> startswith(f, prefix), readdir(dir));
                          isempty(files) ? nothing : joinpath(dir, first(sort(files))))
todict(nt) = Dict{String, Any}(string(k) => v for (k, v) in pairs(nt))   # TOML needs string keys

series = breezelab_timeseries(run_dir)
bounds = breezelab_cloud_boundaries(run_dir)
provenance = TOML.parsefile(joinpath(run_dir, "provenance.toml"))
protocol = series.protocol
benchmark = protocol == "covert_public_bin" ?
    "Covert-public-bin development benchmark (public inputs; NOT an official LASSO-ENA member)" :
    protocol == "lasso_ena_official" ? "LASSO-ENA member $(get(provenance["config"], "member", get(provenance, "protocol_member", "undeclared"))) (Breeze adapter; not a SAM result)" :
    "protocol $protocol"

t0 = first(series.time)
t1 = last(series.time)
window = (; start = string(t0), stop = string(t1))
println("run: ", run_dir, "\n  ", benchmark, "\n  window ", t0, " – ", t1, " UTC (run output cadence ", round(Int, series.seconds[2] - series.seconds[1]), " s)")

les_stats(values) = window_statistics(series.time, values, BitVector(isfinite.(values)), t0, t1)
result = Dict{String, Any}(
    "run" => Dict("directory" => abspath(run_dir), "protocol" => protocol, "benchmark" => benchmark, "label" => series.label,
                  "epoch" => string(series.epoch), "window" => todict(window)),
    "les" => Dict("lwp_g_m2" => todict((les_stats(series.lwp))),
                  "rwp_g_m2" => todict((les_stats(series.rwp))),
                  "cloud_fraction_lwp_gt_5" => todict((les_stats(series.cloud_fraction))),
                  "rain_rate_mm_hr" => todict((les_stats(series.rain_rate))),
                  "cloud_base_m" => todict((window_statistics(bounds.time, bounds.base, BitVector(isfinite.(bounds.base)), t0, t1 + Second(1)))),
                  "cloud_top_m" => todict((window_statistics(bounds.time, bounds.top, BitVector(isfinite.(bounds.top)), t0, t1 + Second(1)))),
                  "definitions" => Dict("lwp" => "domain-mean cloud liquid water path (rain excluded)",
                                        "cloud_fraction" => "fraction of columns with cloud LWP > 5 g m⁻²",
                                        "rain_rate" => "domain-mean surface rain mass flux, 1 kg m⁻² = 1 mm",
                                        "cloud_boundaries" => "lowest/highest level with cloudy-cell fraction > $(bounds.threshold) in each profile-averaging interval")),
    "observations" => Dict{String, Any}())

fig = Figure(size = (1000, 1000))
hours(t) = [Dates.value(x - t0) / 3.6e6 for x in t]
ax_lwp = Axis(fig[1, 1], ylabel = "LWP (g m⁻²)", title = "Liquid water path")
ax_rain = Axis(fig[2, 1], ylabel = "rain rate (mm hr⁻¹)", title = "Surface rain rate")
ax_cloud = Axis(fig[3, 1], xlabel = "hours since $(t0) UTC", ylabel = "height (m)", title = "Cloud base and top")
lines!(ax_lwp, hours(series.time), series.lwp; label = "LES domain mean", color = :black)
lines!(ax_rain, hours(series.time), series.rain_rate; label = "LES domain mean", color = :black)
lines!(ax_cloud, hours(bounds.time), bounds.base; label = "LES base", color = :black)
lines!(ax_cloud, hours(bounds.time), bounds.top; label = "LES top", color = :gray40)

mwr = find_file(obs_dir, "enamwrret2turn")
if !isnothing(mwr)
    lwp = read_arm_lwp(mwr)
    stats = window_statistics(lwp, t0, t1)
    unc = window_statistics(lwp.time, lwp.uncertainty, lwp.good, t0, t1)
    result["observations"]["mwrret_lwp_g_m2"] = merge(todict((stats)), Dict("datastream" => lwp.datastream, "variable" => lwp.variable,
                                                                                "mean_1sigma_uncertainty" => unc.mean,
                                                                                "note" => "zenith point retrieval at the site; the LES value is a domain mean"))
    sel = lwp.good .& [t0 ≤ t ≤ t1 for t in lwp.time]
    scatter!(ax_lwp, hours(lwp.time[sel]), lwp.value[sel]; markersize = 3, color = (:dodgerblue, 0.6), label = "MWRRET phys_lwp (good)")
    println("  MWRRET LWP: obs mean $(round(stats.mean; digits=1)) ± $(round(unc.mean; digits=1)) g m⁻² (n_good = $(stats.n_good)/$(stats.n)); LES mean $(round(les_stats(series.lwp).mean; digits=1)) g m⁻²")
end
vdis = find_file(obs_dir, "enavdis")
if !isnothing(vdis)
    rain = read_arm_rain_rate(vdis)
    stats = window_statistics(rain, t0, t1)
    result["observations"]["vdis_rain_rate_mm_hr"] = merge(todict((stats)), Dict("datastream" => rain.datastream,
                                                                                     "note" => "one-minute disdrometer rain rate at the site"))
    sel = rain.good .& [t0 ≤ t ≤ t1 for t in rain.time]
    lines!(ax_rain, hours(rain.time[sel]), rain.value[sel]; color = (:dodgerblue, 0.8), label = "VDIS rain_rate (good)")
    println("  VDIS rain: obs mean $(round(stats.mean; digits=3)) mm hr⁻¹ (n_good = $(stats.n_good)/$(stats.n)); LES mean $(round(les_stats(series.rain_rate).mean; digits=3)) mm hr⁻¹")
end
arscl = find_file(obs_dir, "enaarsclkazrbnd")
if !isnothing(arscl)
    cb = read_arm_cloud_boundaries(arscl)
    cf = cloud_fraction_from_boundaries(cb, t0, t1; max_base = 3000)
    result["observations"]["arscl_cloud"] = Dict("datastream" => cb.datastream, "n" => cf.n, "n_cloudy" => cf.n_cloudy,
                                                 "time_fraction_cloudy_base_below_3km" => cf.cloud_fraction,
                                                 "base_m" => todict((cf.base)), "top_m" => todict((cf.top)),
                                                 "note" => "lowest hydrometeor layer (radar+lidar; drizzle below cloud base lowers the radar base); the LES cloud fraction is an areal LWP-based fraction")
    sel = cb.cloudy .& [t0 ≤ t ≤ t1 for t in cb.time]
    scatter!(ax_cloud, hours(cb.time[sel]), cb.base[sel]; markersize = 2, color = (:dodgerblue, 0.5), label = "ARSCL base")
    scatter!(ax_cloud, hours(cb.time[sel]), cb.top[sel]; markersize = 2, color = (:orange, 0.5), label = "ARSCL top")
    println("  ARSCL: cloudy time fraction $(round(cf.cloud_fraction; digits=2)), base $(round(cf.base.mean; digits=0)) m, top $(round(cf.top.mean; digits=0)) m; LES base $(round(result["les"]["cloud_base_m"]["mean"]; digits=0)) m, top $(round(result["les"]["cloud_top_m"]["mean"]; digits=0)) m, LWP cloud fraction $(round(les_stats(series.cloud_fraction).mean; digits=2))")
end
axislegend(ax_lwp, position = :lt); axislegend(ax_rain, position = :lt); axislegend(ax_cloud, position = :lt)
Label(fig[0, 1], benchmark, fontsize = 13, tellwidth = false)
save(joinpath(output_dir, "comparison.png"), fig; px_per_unit = 2)
open(io -> TOML.print(io, result), joinpath(output_dir, "comparison.toml"), "w")
println("wrote ", joinpath(output_dir, "comparison.toml"), " and comparison.png")
