# Five-run analysis of the completed ENA-covert campaign against the staged ARM observations:
# per-run comparison (compare_ena_observations.jl), overlays of LWP/cloud fraction/rain,
# hourly-mean profiles (θ, qᵗ, qᶜˡ, nᶜˡ, nᵃ) and the Δt sensitivity pairs.
#   julia --project analysis/ena_covert_analysis.jl OUTPUT_DIR OBS_DIR LABEL=RUN_DIR [LABEL=RUN_DIR ...]
using BreezeLab, CairoMakie, JLD2, Oceananigans, Oceananigans.Units, Statistics, TOML, Dates, Printf
using Oceananigans.Grids: Center, znodes
using Oceananigans.Fields: interior

output_dir, obs_dir = ARGS[1], ARGS[2]
runs = [(split(a, "=", limit=2)[1], split(a, "=", limit=2)[2]) for a in ARGS[3:end]]
mkpath(output_dir)
julia = Base.julia_cmd()
for (label, dir) in runs
    out = joinpath(output_dir, "comparison_$label")
    isfile(joinpath(out, "comparison.toml")) && continue
    run(`$julia --startup-file=no --project=$(dirname(dirname(pathof(BreezeLab)))) $(joinpath(dirname(dirname(pathof(BreezeLab))), "analysis", "compare_ena_observations.jl")) $dir $obs_dir $out`)
end

lwp_obs = read_arm_lwp(only(filter(f -> startswith(basename(f), "enamwrret2turn"), readdir(obs_dir; join=true))))
rain_obs = read_arm_rain_rate(only(filter(f -> startswith(basename(f), "enavdis"), readdir(obs_dir; join=true))))
cb_obs = read_arm_cloud_boundaries(only(filter(f -> startswith(basename(f), "enaarsclkazrbnd"), readdir(obs_dir; join=true))))
colors = Makie.wong_colors()
series = Dict(label => breezelab_timeseries(dir) for (label, dir) in runs)
t0 = minimum(first(s.time) for s in values(series)); t1 = maximum(last(s.time) for s in values(series))
hours(t) = [Dates.value(x - t0) / 3.6e6 for x in t]

fig = Figure(size=(1300, 1000), fontsize=13)
ax1 = Axis(fig[1, 1], ylabel="LWP (g m⁻²)", title="Domain-mean cloud LWP vs MWRRET (zenith)")
ax2 = Axis(fig[1, 2], ylabel="cloud fraction", title="Cloud fraction (LWP > 5 g m⁻²)")
ax3 = Axis(fig[2, 1], xlabel="hours since $(t0) UTC", ylabel="rain (mm hr⁻¹)", title="Surface rain rate vs VDIS")
ax4 = Axis(fig[2, 2], xlabel="hours since $(t0) UTC", ylabel="height (m)", title="Cloud base/top vs ARSCL top and ceilometer base")
sel = lwp_obs.good .& [t0 ≤ t ≤ t1 for t in lwp_obs.time]
scatter!(ax1, hours(lwp_obs.time[sel]), lwp_obs.value[sel]; markersize=2, color=(:gray50, 0.5), label="MWRRET phys_lwp")
# 10-min running median of the observations for readability
selr = rain_obs.good .& [t0 ≤ t ≤ t1 for t in rain_obs.time]
lines!(ax3, hours(rain_obs.time[selr]), rain_obs.value[selr]; color=(:gray50, 0.8), label="VDIS")
selc = cb_obs.cloudy .& [t0 ≤ t ≤ t1 for t in cb_obs.time]
scatter!(ax4, hours(cb_obs.time[selc]), cb_obs.top[selc]; markersize=1.5, color=(:gray50, 0.4), label="ARSCL top")
selb = [t0 ≤ t ≤ t1 && isfinite(b) for (t, b) in zip(cb_obs.time, cb_obs.base_best_estimate)]
scatter!(ax4, hours(cb_obs.time[selb]), cb_obs.base_best_estimate[selb]; markersize=1.5, color=(:black, 0.4), label="ceilometer/MPL base")
summary = Dict{String, Any}()
for (i, (label, dir)) in enumerate(runs)
    s = series[label]; c = colors[mod1(i, 7)]
    lines!(ax1, hours(s.time), s.lwp; color=c, label)
    lines!(ax2, hours(s.time), s.cloud_fraction; color=c, label)
    lines!(ax3, hours(s.time), s.rain_rate; color=c, label)
    b = breezelab_cloud_boundaries(dir)
    lines!(ax4, hours(b.time), b.base; color=c, label); lines!(ax4, hours(b.time), b.top; color=c, linestyle=:dash)
    summary[label] = Dict("lwp_mean" => mean(s.lwp), "lwp_final" => s.lwp[end], "cf_mean" => mean(s.cloud_fraction),
                          "rain_mean_mm_hr" => mean(s.rain_rate), "rain_max_mm_hr" => maximum(s.rain_rate),
                          "base_final" => b.base[end], "top_final" => b.top[end], "protocol" => s.protocol)
end
axislegend(ax1, position=:lt, labelsize=10); axislegend(ax3, position=:lt, labelsize=10); axislegend(ax4, position=:lb, labelsize=9)
Label(fig[0, :], "ENA-covert public-bin development benchmark runs (NOT LASSO) vs ARM ENA observations, 18 July 2017", fontsize=14)
save(joinpath(output_dir, "timeseries_all_runs.png"), fig; px_per_unit=2)

# profiles at the last hourly mean (and the first) for θ, qᵗ, qᶜˡ, nᶜˡ, nᵃ
function prof(dir, name)
    file = only(filter(f -> endswith(f, "_profiles.jld2"), readdir(dir; join=true)))
    present = jldopen(f -> haskey(f["timeseries"], name), file)
    return present ? FieldTimeSeries(file, name) : nothing
end
col(f) = vec(Array(interior(f)))
fp = Figure(size=(1600, 900), fontsize=13)
axes = [Axis(fp[1, k], ylabel="z (m)", limits=(nothing, (0, 2500)), title=t) for (k, t) in enumerate(("θ (K)", "qᵗ = qᵛ+qᶜˡ (g kg⁻¹)", "qᶜˡ (g kg⁻¹)", "cloudy fraction", "nᶜˡ (cm⁻³, ρ₀-scaled)", "nᵃ (cm⁻³, ρ₀-scaled)"))]
for (i, (label, dir)) in enumerate(runs)
    c = colors[mod1(i, 7)]
    θ = prof(dir, "θ"); qv = prof(dir, "qᵛ"); qc = prof(dir, "qᶜˡ"); cf = prof(dir, "cloud_fraction"); nc = prof(dir, "nᶜˡ"); na = prof(dir, "nᵃ")
    z = collect(znodes(θ.grid, Center()))
    ρ0 = get(TOML.parsefile(joinpath(dir, "provenance.toml"))["config"], "reference_density", 1.206)   # only P3-aer runs record it; 1.206 is the ENA surface value
    for (n, ls) in ((1, :dot), (length(θ.times), :solid))
        lines!(axes[1], col(θ[n]), z; color=c, linestyle=ls, label=(ls == :solid ? label : nothing))
        lines!(axes[2], 1e3 .* (col(qv[n]) .+ col(qc[n])), z; color=c, linestyle=ls)
        lines!(axes[3], 1e3 .* col(qc[n]), z; color=c, linestyle=ls)
        lines!(axes[4], col(cf[n]), z; color=c, linestyle=ls)
        isnothing(nc) || lines!(axes[5], 1e-6 .* ρ0 .* col(nc[n]), z; color=c, linestyle=ls)
        isnothing(na) || lines!(axes[6], 1e-6 .* ρ0 .* col(na[n]), z; color=c, linestyle=ls)
    end
    k = length(θ.times); cfk = col(cf[k]); cloudy = cfk .> 0.05
    summary[label]["final_hour"] = Dict("cloud_base_m" => any(cloudy) ? z[findfirst(cloudy)] : NaN, "cloud_top_m" => any(cloudy) ? z[findlast(cloudy)] : NaN,
                                        "qc_max_gkg" => 1e3 * maximum(col(qc[k])),
                                        "nc_incloud_cm3" => isnothing(nc) ? NaN : 1e-6 * ρ0 * mean((col(nc[k]) ./ max.(cfk, 1e-6))[cloudy]),
                                        "na_incloud_cm3" => isnothing(na) ? NaN : 1e-6 * ρ0 * mean(col(na[k])[cloudy]),
                                        "na_min_cm3" => isnothing(na) ? NaN : 1e-6 * ρ0 * minimum(col(na[k])),
                                        "theta_surface" => col(θ[k])[1], "qt_surface_gkg" => 1e3 * (col(qv[k])[1] + col(qc[k])[1]))
end
axislegend(axes[1], position=:rb, labelsize=10)
Label(fp[0, :], "Hourly-mean profiles: first hour (dotted) and last hour (solid); ENA-covert benchmark runs", fontsize=14)
save(joinpath(output_dir, "profiles_all_runs.png"), fp; px_per_unit=2)
open(io -> TOML.print(io, summary), joinpath(output_dir, "summary.toml"), "w")
for (label, s) in summary
    @printf("%-22s LWP mean %.1f final %.1f g/m²; CF %.2f; rain mean %.4f max %.4f mm/hr; base/top %.0f/%.0f m; nc in-cloud %.0f cm⁻³; na in-cloud %.0f cm⁻³ (min %.0f)\n",
            label, s["lwp_mean"], s["lwp_final"], s["cf_mean"], s["rain_mean_mm_hr"], s["rain_max_mm_hr"], s["final_hour"]["cloud_base_m"], s["final_hour"]["cloud_top_m"],
            s["final_hour"]["nc_incloud_cm3"], s["final_hour"]["na_incloud_cm3"], s["final_hour"]["na_min_cm3"])
end
println("wrote ", output_dir)
