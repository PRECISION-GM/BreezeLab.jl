# Compare a BreezeLab LASSO-ENA run with the member's SAM reference statistics (samstat) and
# the ARM observations over the common UTC window (plan step 7). Usage:
#   julia --project analysis/compare_sam_reference.jl RUN_DIR SAMSTAT_NC [OBS_DIR] [OUTPUT_DIR]
# SAM's samstat carries domain statistics every nstat steps (time in fractional day of year);
# BreezeLab's run carries its UTC epoch in provenance. Writes sam_comparison.{toml,png}.
using BreezeLab, NCDatasets, CairoMakie, Dates, TOML, Statistics
using Oceananigans
using Oceananigans.Grids: Center, znodes
using Oceananigans.Fields: interior

stage(msg) = (println(stderr, "[", Dates.format(Dates.now(), "HH:MM:SS"), "] ", msg); flush(stderr))
stage("loaded")
run_dir, samstat = ARGS[1], ARGS[2]
obs_dir = length(ARGS) ≥ 3 ? ARGS[3] : nothing
output_dir = length(ARGS) ≥ 4 ? ARGS[4] : joinpath(run_dir, "sam_comparison")
mkpath(output_dir)

stage("reading run output")
series = breezelab_timeseries(run_dir)
bounds = breezelab_cloud_boundaries(run_dir)
prov = TOML.parsefile(joinpath(run_dir, "provenance.toml"))
member = get(prov, "protocol_member", "undeclared")
year = Dates.year(series.epoch)
t0, t1 = first(series.time), last(series.time)

# SAM time axis: fractional day of year, 1.0 = 1 January 00 UTC (as day0)
stage("opening samstat")
sam = NCDataset(samstat)
sam_days = Float64.(sam["time"][:])
sam_time = [DateTime(year, 1, 1) + Millisecond(round(Int, (d - 1) * 86400e3)) for d in sam_days]
sam_sim = string(get(sam.attrib, "sim_name", "?"))
sam_hash = string(get(sam.attrib, "model_source_git_hash", "?"))
member == sam_sim || @warn "run member $member differs from samstat sim_name $sam_sim"
inwin = [t0 ≤ t ≤ t1 for t in sam_time]
sv(name) = Float64[ismissing(v) ? NaN : Float64(v) for v in sam[name][:]]
todict(nt) = Dict{String, Any}(string(k) => v for (k, v) in pairs(nt))   # TOML needs string keys
function stats(x)
    v = filter(isfinite, x)
    isempty(v) && return (; mean = NaN, min = NaN, max = NaN, n = 0)
    return (; mean = mean(v), min = minimum(v), max = maximum(v), n = length(v))
end

result = Dict{String, Any}(
    "run" => Dict("directory" => abspath(run_dir), "member" => member, "label" => series.label, "window" => Dict("start" => string(t0), "stop" => string(t1))),
    "sam" => Dict("file" => abspath(samstat), "sim_name" => sam_sim, "model_source_git_hash" => sam_hash, "samples_in_window" => count(inwin)))

quantities = [("CWP", "lwp", "cloud water path (g m⁻²)", 1.0),
         ("RWP", "rwp", "rain water path (g m⁻²)", 1.0),
         ("PREC", "rain_rate", "surface precipitation (mm hr⁻¹)", 1 / 24),   # SAM mm/day → mm/hr
         ("CLDSHD", "cloud_fraction", "shaded cloud fraction", 1.0)]
stage("time series figure")
fig = Figure(size = (1300, 1300), fontsize = 13)
hours(t) = Float64[Dates.value(x - t0) / 3.6e6 for x in t]
# Makie cannot set limits from an all-NaN series: draw only the finite points, skip empty series.
function plot_finite!(ax, x, y; kwargs...)
    keep = BitVector(isfinite.(y))
    any(keep) && lines!(ax, x[keep], y[keep]; kwargs...)
    return nothing
end
for (k, (sname, bname, title, scale)) in enumerate(quantities)
    stage("panel $sname: axis")
    ax = Axis(fig[(k - 1) ÷ 2 + 1, (k - 1) % 2 + 1], title = title, xlabel = "hours since $(t0) UTC")
    stage("panel $sname: read")
    s = sv(sname) .* scale
    stage("panel $sname: lines (SAM)")
    plot_finite!(ax, hours(sam_time[inwin]), s[inwin]; color = :gray30, label = "SAM $sname")
    b = getproperty(series, Symbol(bname))
    stage("panel $sname: lines (Breeze)")
    plot_finite!(ax, hours(series.time), b; color = :dodgerblue, label = "Breeze $bname")
    stage("panel $sname: legend")
    isempty(ax.scene.plots) || axislegend(ax, position = :lt, labelsize = 10)
    stage("panel $sname: stats")
    result[sname] = Dict{String, Any}("sam" => todict(stats(s[inwin])), "breeze" => todict(stats(b)))
end
stage("cloud boundaries panel")
ax5 = Axis(fig[3, 1], title = "cloud base/top (m): SAM GCSS ZCB/ZCT vs Breeze profile boundaries", xlabel = "hours since $(t0) UTC")
zcb = 1e3 .* sv("ZCB"); zct = 1e3 .* sv("ZCT")
plot_finite!(ax5, hours(sam_time[inwin]), zcb[inwin]; color = :gray30, label = "SAM ZCB")
plot_finite!(ax5, hours(sam_time[inwin]), zct[inwin]; color = :gray30, linestyle = :dash, label = "SAM ZCT")
plot_finite!(ax5, hours(bounds.time), bounds.base; color = :dodgerblue, label = "Breeze base")
plot_finite!(ax5, hours(bounds.time), bounds.top; color = :dodgerblue, linestyle = :dash, label = "Breeze top")
isempty(ax5.scene.plots) || axislegend(ax5, position = :lt, labelsize = 10)
result["ZCB_m"] = Dict("sam" => todict((stats(zcb[inwin]))), "breeze" => todict((stats(bounds.base))))
result["ZCT_m"] = Dict("sam" => todict((stats(zct[inwin]))), "breeze" => todict((stats(bounds.top))))
ax6 = Axis(fig[3, 2], title = "SAM surface fluxes and radiation (W m⁻²; Breeze run output has no flux series yet)", xlabel = "hours since $(t0) UTC")
for (name, c) in (("SHF", :orange), ("LHF", :red), ("LWNS", :purple), ("SWNS", :gold))
    v = sv(name); plot_finite!(ax6, hours(sam_time[inwin]), v[inwin]; color = c, label = "SAM $name")
    result[name] = Dict("sam" => todict((stats(v[inwin]))))
end
isempty(ax6.scene.plots) || axislegend(ax6, position = :lt, labelsize = 10)
if !isnothing(obs_dir)
    mwr = filter(f -> startswith(basename(f), "enamwrret2turn"), readdir(obs_dir; join = true))
    if !isempty(mwr)
        lwp = read_arm_lwp(first(mwr)); sel = lwp.good .& [t0 ≤ t ≤ t1 for t in lwp.time]
        any(sel) && scatter!(content(fig[1, 1]), hours(lwp.time[sel]), lwp.value[sel]; markersize = 2, color = (:black, 0.3), label = "MWRRET")
        result["CWP"]["mwrret"] = todict(window_statistics(lwp, t0, t1))
        println("MWRRET good samples in window: ", count(sel))
    end
end
Label(fig[0, :], "Breeze vs SAM reference $sam_sim (SAM $sam_hash); Breeze is an adapter of the forcing, not the SAM model", fontsize = 13, tellwidth = false)
stage("saving time series figure")
save(joinpath(output_dir, "sam_comparison.png"), fig; px_per_unit = 2)
stage("profiles")

# profiles at the last common hour: θ, qᵛ, qᶜˡ, u, v
file = only(filter(f -> endswith(f, "_profiles.jld2"), readdir(run_dir; join = true)))
pf = Dict(n => FieldTimeSeries(file, n) for n in ("θ", "qᵛ", "qᶜˡ", "u", "v"))
zb = collect(znodes(pf["θ"].grid, Center()))
nb = length(pf["θ"].times)
tb_end = series.epoch + Millisecond(round(Int, 1e3 * pf["θ"].times[nb]))
js = argmin(abs.(Dates.value.(sam_time .- tb_end)))
zs = Float64.(sam["z"][:])
col(f, n) = vec(Array(interior(f[n])))
fp = Figure(size = (1500, 600), fontsize = 13)
for (k, (sname, bname, scale, title)) in enumerate((("THETA", "θ", 1.0, "θ (K)"), ("QV", "qᵛ", 1e3, "qᵛ (g kg⁻¹)"), ("QCL", "qᶜˡ", 1e3, "qᶜˡ (g kg⁻¹)"), ("U", "u", 1.0, "u (m s⁻¹)"), ("V", "v", 1.0, "v (m s⁻¹)")))
    ax = Axis(fp[1, k], title = title, ylabel = "z (m)", limits = (nothing, (0, 3000)))
    s = Float64[ismissing(v) ? NaN : Float64(v) for v in sam[sname][:, js]]
    keep = BitVector(isfinite.(s)); any(keep) && lines!(ax, s[keep], zs[keep]; color = :gray30, label = "SAM $(sam_time[js])")
    lines!(ax, scale .* col(pf[bname], nb), zb; color = :dodgerblue, label = "Breeze $(tb_end)")
    k == 1 && !isempty(ax.scene.plots) && axislegend(ax, position = :rb, labelsize = 9)
end
Label(fp[0, :], "Profiles at the last common hour (SAM instantaneous statistics sample vs Breeze hourly mean)", fontsize = 13, tellwidth = false)
stage("saving profiles figure")
save(joinpath(output_dir, "sam_profiles.png"), fp; px_per_unit = 2)
close(sam)
open(io -> TOML.print(io, result), joinpath(output_dir, "sam_comparison.toml"), "w")
for (k, v) in result
    v isa Dict && haskey(v, "sam") && haskey(v, "breeze") && println(rpad(k, 8), " SAM mean ", round(v["sam"]["mean"]; sigdigits = 4), "  Breeze mean ", round(v["breeze"]["mean"]; sigdigits = 4))
end
println("wrote ", output_dir)
