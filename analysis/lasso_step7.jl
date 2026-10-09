# Plan step 7 for a full LASSO-ENA member run: Breeze vs the SAM member's samstat/sam2d
# reference outputs and the ARM observations, over the full run and the 06–12 / 09–12 UTC
# windows. Usage (CPU node):
#   julia --project analysis/lasso_step7.jl RUN_DIR SAMSTAT_NC SAM2D_NC OBS_DIR OUTPUT_DIR [PLAN_HOURS=6,18]
# All statistics are computed before any figure is created (see analysis/README.md).
using BreezeLab, NCDatasets, JLD2, CairoMakie, Dates, TOML, Statistics, Printf
using Oceananigans
using Oceananigans.Grids: Center, Face, znodes
using Oceananigans.Fields: interior

stage(msg) = (println(stderr, "[", Dates.format(Dates.now(), "HH:MM:SS"), "] ", msg); flush(stderr))
run_dir, samstat_file, sam2d_file, obs_dir, output_dir = ARGS[1:5]
plan_hours = length(ARGS) ≥ 6 ? parse.(Float64, split(ARGS[6], ",")) : [6.0, 18.0]
mkpath(output_dir)
todict(nt) = Dict{String, Any}(string(k) => v for (k, v) in pairs(nt))
fin(x) = Float64[ismissing(v) ? NaN : Float64(v) for v in x]
function stats(x)
    v = filter(isfinite, x)
    isempty(v) && return (; mean = NaN, min = NaN, max = NaN, n = 0)
    return (; mean = mean(v), min = minimum(v), max = maximum(v), n = length(v))
end

# ---------------------------------------------------------------- Breeze run
stage("run output")
series = breezelab_timeseries(run_dir)
# Cloud boundaries of the *stratocumulus* layer: the highest contiguous cloudy layer (cloudy-cell
# fraction > threshold) of each profile, so that the surface fog / drizzle-moistened layer below
# a clear gap does not set the base (SAM's GCSS ZCB/ZCT are column-wise means with their own rules).
function layer_boundaries(run_dir; threshold = 0.05)
    file = only(filter(f -> endswith(f, "_profiles.jld2"), readdir(run_dir; join = true)))
    cf = FieldTimeSeries(file, "cloud_fraction"); z = collect(znodes(cf.grid, Center()))
    base = Float64[]; top = Float64[]; fog_top = Float64[]
    for n in eachindex(cf.times)
        c = vec(Array(interior(cf[n]))) .> threshold
        ks = findall(c)
        if isempty(ks)
            push!(base, NaN); push!(top, NaN); push!(fog_top, NaN); continue
        end
        # split into contiguous runs; take the highest run as the stratocumulus layer
        runs = Vector{UnitRange{Int}}(); start = ks[1]; prev = ks[1]
        for k in ks[2:end]
            k == prev + 1 || (push!(runs, start:prev); start = k)
            prev = k
        end
        push!(runs, start:prev)
        main = last(runs)
        push!(base, z[first(main)]); push!(top, z[last(main)])
        push!(fog_top, length(runs) > 1 ? z[last(runs[1])] : NaN)
    end
    return (; seconds = collect(Float64, cf.times), base, top, fog_top, threshold)
end
bounds = layer_boundaries(run_dir)
bounds_all = breezelab_cloud_boundaries(run_dir)   # lowest cloudy level of any layer (fog included)
prov = TOML.parsefile(joinpath(run_dir, "provenance.toml"))
member = get(prov, "protocol_member", "undeclared")
epoch = series.epoch
year = Dates.year(epoch)
hours_b = series.seconds ./ 3600
pfile = only(filter(f -> endswith(f, "_profiles.jld2"), readdir(run_dir; join = true)))
pnames = jldopen(f -> [k for k in keys(f["timeseries"]) if k != "t"], pfile)
P = Dict(n => FieldTimeSeries(pfile, n) for n in ("θ", "qᵛ", "qᶜˡ", "qʳ", "u", "v", "T", "cloud_fraction", "nᶜˡ", "radiative_flux_divergence") if n in pnames)
zb = collect(znodes(P["θ"].grid, Center()))
tp = collect(Float64, P["θ"].times)                   # seconds, profile (10-min) averaging intervals
col(f, n) = vec(Array(interior(f[n])))
# 2-D slices (lwp, rain) read raw from JLD2 (halo 5), with their times
sfile = only(filter(f -> endswith(f, "_slices.jld2"), readdir(run_dir; join = true)))
slice_times, slice_lwp, slice_rain = jldopen(sfile) do f
    ts = f["timeseries"]
    its = sort(parse.(Int, collect(keys(ts["t"]))))
    t = [Float64(ts["t/$i"]) for i in its]
    halo = 5
    lw = [Array{Float64}(ts["lwp/$i"][halo+1:end-halo, halo+1:end-halo, 1]) for i in its]
    rn = [Array{Float64}(ts["rain/$i"][halo+1:end-halo, halo+1:end-halo, 1]) for i in its]
    t, lw, rn
end

# ---------------------------------------------------------------- SAM samstat
stage("samstat")
sam = NCDataset(samstat_file)
sam_days = fin(sam["time"][:])
sam_hours = (sam_days .- sam_days[1]) .* 24          # hours since 00 UTC of day0
zs = fin(sam["z"][:])
sv(n) = fin(sam[n][:])
samp(n, j) = fin(sam[n][:, j])
sam_rho = samp("RHO", 1)
sam_attr = Dict(k => string(v) for (k, v) in sam.attrib)
# SAM layer boundaries from the QCL profiles with the paper/Breeze definition (highest contiguous layer
# with qc ≥ 0.01 g kg⁻¹); the GCSS ZCB/ZCT series stored in samstat (units "km") average 0.03/0.10 km,
# inconsistent with ZINV (1.5 km), ZCTMAX (1.7 km) and the QCL profiles, so they are quoted raw only.
function sam_layer_boundaries(sam, zs; threshold = 0.01)
    qcl = sam["QCL"][:, :]
    base = Float64[]; top = Float64[]; fog_top = Float64[]
    for j in axes(qcl, 2)
        c = [!ismissing(v) && v ≥ threshold for v in qcl[:, j]]
        ks = findall(c)
        if isempty(ks)
            push!(base, NaN); push!(top, NaN); push!(fog_top, NaN); continue
        end
        runs = Vector{UnitRange{Int}}(); start = ks[1]; prev = ks[1]
        for k in ks[2:end]
            k == prev + 1 || (push!(runs, start:prev); start = k)
            prev = k
        end
        push!(runs, start:prev)
        push!(base, zs[first(last(runs))]); push!(top, zs[last(last(runs))])
        push!(fog_top, length(runs) > 1 ? zs[last(runs[1])] : NaN)
    end
    return (; base, top, fog_top)
end
sam_layers = sam_layer_boundaries(sam, zs)

# ---------------------------------------------------------------- windows
windows = Dict("00-24" => (0.0, 24.0), "06-12" => (6.0, 12.0), "09-12" => (9.0, 12.0))
inwin(h, w) = w[1] ≤ h ≤ w[2]
breeze_window(h, x, w) = stats(x[[inwin(t, w) for t in h]])
sam_window(n, w, scale = 1.0) = stats((sv(n) .* scale)[[inwin(t, w) for t in sam_hours]])
result = Dict{String, Any}("run" => Dict("directory" => abspath(run_dir), "member" => member, "epoch" => string(epoch), "label" => series.label),
                           "sam" => Dict("samstat" => abspath(samstat_file), "sam2d" => abspath(sam2d_file), "sim_name" => get(sam_attr, "sim_name", "?"),
                                         "model_source_git_hash" => get(sam_attr, "model_source_git_hash", "?")),
                           "windows" => Dict{String, Any}())
ts_quantities = [("CWP", series.lwp, "cloud water path (g m⁻²)", 1.0), ("RWP", series.rwp, "rain water path (g m⁻²)", 1.0),
                 ("PREC", series.rain_rate .* 24, "surface precipitation (mm day⁻¹)", 1.0), ("CLDSHD", series.cloud_fraction, "cloud fraction", 1.0)]
for (wname, w) in windows
    d = Dict{String, Any}()
    for (sname, bvalues, _, scale) in ts_quantities
        d[sname] = Dict("sam" => todict(sam_window(sname, w, scale)), "breeze" => todict(breeze_window(hours_b, bvalues, w)))
    end
    hb = bounds.seconds ./ 3600
    d["ZCB_m"] = Dict("sam" => todict(breeze_window(sam_hours, sam_layers.base, w)), "breeze" => todict(breeze_window(hb, bounds.base, w)),
                      "sam_gcss_zcb_raw_km" => todict(sam_window("ZCB", w)),
                      "breeze_lowest_cloudy_level" => todict(breeze_window(hb, bounds_all.base, w)),
                      "breeze_fog_layer_top" => todict(breeze_window(hb, bounds.fog_top, w)), "sam_fog_layer_top" => todict(breeze_window(sam_hours, sam_layers.fog_top, w)),
                      "definition" => "highest contiguous layer with qc ≥ 0.01 g kg⁻¹ (SAM: QCL profile; Breeze: cloudy-cell fraction > 0.05)")
    d["ZCT_m"] = Dict("sam" => todict(breeze_window(sam_hours, sam_layers.top, w)), "breeze" => todict(breeze_window(hb, bounds.top, w)),
                      "sam_gcss_zct_raw_km" => todict(sam_window("ZCT", w)), "sam_zinv_m" => todict(sam_window("ZINV", w, 1e3)), "sam_zctmax_m" => todict(sam_window("ZCTMAX", w, 1e3)))
    for n in ("SHF", "LHF", "LWNS", "SWNS", "LWNT", "SWNT")
        d[n] = Dict("sam" => todict(sam_window(n, w)), "breeze" => "not saved by the run")
    end
    result["windows"][wname] = d
end

# ---------------------------------------------------------------- ARM observations over the windows
stage("observations")
obs = Dict{String, Any}()
mwr = filter(f -> startswith(basename(f), "enamwrret2turn"), readdir(obs_dir; join = true))
vdis = filter(f -> startswith(basename(f), "enavdis"), readdir(obs_dir; join = true))
arscl = filter(f -> startswith(basename(f), "enaarsclkazrbnd"), readdir(obs_dir; join = true))
lwp_obs = isempty(mwr) ? nothing : read_arm_lwp(first(mwr))
rain_obs = isempty(vdis) ? nothing : read_arm_rain_rate(first(vdis))
cb_obs = isempty(arscl) ? nothing : read_arm_cloud_boundaries(first(arscl))
for (wname, w) in windows
    t0 = epoch + Second(round(Int, 3600w[1])); t1 = epoch + Second(round(Int, 3600w[2]))
    d = Dict{String, Any}()
    isnothing(lwp_obs) || (d["mwrret_lwp_g_m2"] = todict(window_statistics(lwp_obs, t0, t1)))
    isnothing(rain_obs) || (d["vdis_rain_mm_day"] = todict(map(x -> x isa Number && isfinite(x) ? 24x : x, window_statistics(rain_obs, t0, t1))))
    if !isnothing(cb_obs)
        cf = cloud_fraction_from_boundaries(cb_obs, t0, t1; max_base = 3000)
        d["arscl"] = Dict("time_fraction_cloudy" => cf.cloud_fraction, "lowest_layer_top_m" => todict(cf.top),
                          "ceilometer_base_m" => todict(window_statistics(cb_obs.time, cb_obs.base_best_estimate, BitVector(isfinite.(cb_obs.base_best_estimate)), t0, t1)))
    end
    obs[wname] = d
end
result["observations"] = obs

# ---------------------------------------------------------------- profiles at selected hours
stage("profiles")
profile_hours = [3.0, 6.0, 9.0, 12.0, 18.0, 23.0]
nearest(ts, h) = argmin(abs.(ts .- h))
profiles = Dict{String, Any}()
rows = Any[]
for h in profile_hours
    nb = nearest(tp ./ 3600, h); js = nearest(sam_hours, h)
    bθ = col(P["θ"], nb); bqv = col(P["qᵛ"], nb); bqc = col(P["qᶜˡ"], nb); bqr = col(P["qʳ"], nb); bu = col(P["u"], nb); bv = col(P["v"], nb)
    bcf = col(P["cloud_fraction"], nb)
    sθ = samp("THETA", js); sqv = samp("QV", js); sqc = samp("QCL", js); sqr = samp("QPL", js); su = samp("U", js); svv = samp("V", js)
    # radiative heating (K day⁻¹): Breeze flux divergence [W m⁻³] / (ρ cp); ρ from SAM's RHO interpolated to the Breeze levels
    heat_b = haskey(P, "radiative_flux_divergence") ? col(P["radiative_flux_divergence"], nb) ./ (interpolate_profile(zs, sam_rho, zb) .* 1004) .* 86400 : fill(NaN, length(zb))
    heat_s = samp("RADQR", js); heat_lw = samp("RADQRLW", js); heat_sw = samp("RADQRSW", js)
    # in-cloud droplet number: Breeze nᶜˡ/CF (kg⁻¹ → cm⁻³ with ρ); SAM samstat carries no droplet number
    ρb = interpolate_profile(zs, sam_rho, zb)
    cloudy = bcf .> 0.05
    nc_b = haskey(P, "nᶜˡ") && any(cloudy) ? 1e-6 * mean((ρb .* col(P["nᶜˡ"], nb) ./ max.(bcf, 1e-6))[cloudy]) : NaN
    push!(rows, (; h, zb, zs, bθ, bqv, bqc, bqr, bu, bv, sθ, sqv, sqc, sqr, su, svv, heat_b, heat_s, heat_lw, heat_sw,
                   tb = epoch + Second(round(Int, tp[nb])), tsam = epoch + Millisecond(round(Int, 3.6e6 * sam_hours[js]))))
    profiles["h$(Int(h))"] = Dict("breeze_time" => string(epoch + Second(round(Int, tp[nb]))), "sam_time" => string(epoch + Millisecond(round(Int, 3.6e6 * sam_hours[js]))),
                                  "breeze_max_qc_gkg" => 1e3 * maximum(bqc), "sam_max_qc_gkg" => maximum(sqc),
                                  "breeze_max_qr_gkg" => 1e3 * maximum(bqr), "sam_max_qr_gkg" => maximum(sqr),
                                  "breeze_in_cloud_nc_cm3" => nc_b, "sam_nc" => "not in samstat (sam3dmicro not downloaded)",
                                  "breeze_min_heating_K_day" => minimum(filter(isfinite, heat_b)), "sam_min_heating_K_day" => minimum(filter(isfinite, heat_s)),
                                  "theta_rms_below_2km_K" => sqrt(mean((bθ[zb .< 2000] .- interpolate_profile(zs, sθ, zb[zb .< 2000])) .^ 2)),
                                  "qv_rms_below_2km_gkg" => sqrt(mean((1e3 .* bqv[zb .< 2000] .- interpolate_profile(zs, sqv, zb[zb .< 2000])) .^ 2)))
end
result["profiles"] = profiles

# ---------------------------------------------------------------- plan views (sam2d vs Breeze 2-D slices)
stage("sam2d")
s2 = NCDataset(sam2d_file)
s2_hours = (fin(s2["time"][:]) .- sam_days[1]) .* 24
x2 = fin(s2["x"][:]) ./ 1e3; y2 = fin(s2["y"][:]) ./ 1e3
plans = Any[]
plan_stats = Dict{String, Any}()
for h in plan_hours
    j2 = nearest(s2_hours, h); jb = nearest(slice_times ./ 3600, h)
    clwp = 1e3 .* Float64.(s2["CLWP"][:, :, j2])      # mm → g m⁻²
    prec = Float64.(s2["Prec"][:, :, j2])               # mm day⁻¹
    zc = 1e3 .* Float64.(s2["ZC"][:, :, j2])            # km → m
    blwp = 1e3 .* slice_lwp[jb]; brain = 86400 .* slice_rain[jb]
    push!(plans, (; h, x2, y2, clwp, prec, zc, blwp, brain, tsam = s2_hours[j2], tb = slice_times[jb] / 3600))
    plan_stats["h$(Int(h))"] = Dict("sam_time_h" => s2_hours[j2], "breeze_time_h" => slice_times[jb] / 3600,
                                   "sam_clwp" => todict(stats(vec(clwp))), "breeze_lwp" => todict(stats(vec(blwp))),
                                   "sam_clwp_std" => std(vec(clwp)), "breeze_lwp_std" => std(vec(blwp)),
                                   "sam_fraction_lwp_gt_20" => mean(vec(clwp) .> 20), "breeze_fraction_lwp_gt_20" => mean(vec(blwp) .> 20),
                                   "sam_prec_mm_day" => todict(stats(vec(prec))), "breeze_prec_mm_day" => todict(stats(vec(brain))),
                                   "sam_fraction_raining_gt_0.1" => mean(vec(prec) .> 0.1), "breeze_fraction_raining_gt_0.1" => mean(vec(brain) .> 0.1),
                                   "sam_cloud_top_m" => todict(stats(vec(zc))))
end
result["plan_views"] = plan_stats
close(s2)

# SAM series for the figures (read before plotting)
sam_series = Dict(n => sv(n) for n in ("CWP", "RWP", "PREC", "CLDSHD", "ZCB", "ZCT", "ZINV", "SHF", "LHF", "LWNS", "SWNS"))
close(sam)
open(io -> TOML.print(io, result), joinpath(output_dir, "step7.toml"), "w")

# ---------------------------------------------------------------- figures (after all statistics)
stage("figures: time series")
hb = bounds.seconds ./ 3600
fig = Figure(size = (1600, 1500), fontsize = 13)
function shade!(ax)
    vspan!(ax, [6.0], [12.0]; color = (:gray, 0.12)); vspan!(ax, [9.0], [12.0]; color = (:gray, 0.12))
end
ax1 = Axis(fig[1, 1], title = "cloud water path (g m⁻²)"); shade!(ax1)
lines!(ax1, sam_hours, sam_series["CWP"]; color = :gray30, label = "SAM CWP"); lines!(ax1, hours_b, series.lwp; color = :dodgerblue, label = "Breeze")
if !isnothing(lwp_obs)
    sel = lwp_obs.good
    scatter!(ax1, [Dates.value(t - epoch) / 3.6e6 for t in lwp_obs.time[sel]], lwp_obs.value[sel]; markersize = 1.5, color = (:black, 0.25), label = "MWRRET")
end
axislegend(ax1, position = :lt, labelsize = 9)
ax2 = Axis(fig[1, 2], title = "rain water path (g m⁻²)"); shade!(ax2)
lines!(ax2, sam_hours, sam_series["RWP"]; color = :gray30, label = "SAM"); lines!(ax2, hours_b, series.rwp; color = :dodgerblue, label = "Breeze"); axislegend(ax2, position = :lt, labelsize = 9)
ax3 = Axis(fig[2, 1], title = "surface precipitation (mm day⁻¹)"); shade!(ax3)
lines!(ax3, sam_hours, sam_series["PREC"]; color = :gray30, label = "SAM"); lines!(ax3, hours_b, 24 .* series.rain_rate; color = :dodgerblue, label = "Breeze")
if !isnothing(rain_obs)
    sel = rain_obs.good
    lines!(ax3, [Dates.value(t - epoch) / 3.6e6 for t in rain_obs.time[sel]], 24 .* rain_obs.value[sel]; color = (:black, 0.5), label = "VDIS")
end
axislegend(ax3, position = :lt, labelsize = 9)
ax4 = Axis(fig[2, 2], title = "cloud fraction"); shade!(ax4)
lines!(ax4, sam_hours, sam_series["CLDSHD"]; color = :gray30, label = "SAM CLDSHD"); lines!(ax4, hours_b, series.cloud_fraction; color = :dodgerblue, label = "Breeze (LWP > 5 g m⁻²)"); axislegend(ax4, position = :lb, labelsize = 9)
ax5 = Axis(fig[3, 1], title = "stratocumulus base / top (m): highest contiguous cloudy layer", xlabel = "hours since $(epoch) UTC", limits = (nothing, (0, 2500))); shade!(ax5)
keep = isfinite.(sam_layers.base); lines!(ax5, sam_hours[keep], sam_layers.base[keep]; color = :gray30, label = "SAM base (QCL ≥ 0.01 g/kg)"); keep = isfinite.(sam_layers.top); lines!(ax5, sam_hours[keep], sam_layers.top[keep]; color = :gray30, linestyle = :dash, label = "SAM top")
lines!(ax5, sam_hours, 1e3 .* sam_series["ZINV"]; color = :gray60, linestyle = :dot, label = "SAM ZINV")
keep = isfinite.(bounds.base); lines!(ax5, hb[keep], bounds.base[keep]; color = :dodgerblue, label = "Breeze base"); keep = isfinite.(bounds.top); lines!(ax5, hb[keep], bounds.top[keep]; color = :dodgerblue, linestyle = :dash, label = "Breeze top")
if !isnothing(cb_obs)
    sel = cb_obs.cloudy
    scatter!(ax5, [Dates.value(t - epoch) / 3.6e6 for t in cb_obs.time[sel]], cb_obs.top[sel]; markersize = 1, color = (:orange, 0.3), label = "ARSCL top")
    selb = isfinite.(cb_obs.base_best_estimate)
    scatter!(ax5, [Dates.value(t - epoch) / 3.6e6 for t in cb_obs.time[selb]], cb_obs.base_best_estimate[selb]; markersize = 1, color = (:black, 0.3), label = "ceilometer base")
end
axislegend(ax5, position = :lt, labelsize = 8)
ax6 = Axis(fig[3, 2], title = "SAM surface fluxes and net radiation (W m⁻²); Breeze run saved none", xlabel = "hours since $(epoch) UTC"); shade!(ax6)
for (n, c) in (("SHF", :orange), ("LHF", :red), ("LWNS", :purple), ("SWNS", :gold))
    lines!(ax6, sam_hours, sam_series[n]; color = c, label = "SAM $n")
end
axislegend(ax6, position = :lt, labelsize = 9)
Label(fig[0, :], "Breeze vs SAM member $member (SAM $(get(sam_attr, "model_source_git_hash", "?"))) and ARM observations, 18 July 2017; shaded: 06–12 and 09–12 UTC windows", fontsize = 13, tellwidth = false)
save(joinpath(output_dir, "step7_timeseries.png"), fig; px_per_unit = 2)

stage("figures: profiles")
fp = Figure(size = (1800, 300 * length(rows) + 60), fontsize = 12)
for (i, r) in enumerate(rows)
    panels = (("θ (K)", r.bθ, r.sθ, 1.0, (285, 320)), ("qᵛ (g kg⁻¹)", r.bqv, r.sqv, 1e3, (0, 17)), ("qᶜˡ (g kg⁻¹)", r.bqc, r.sqc, 1e3, nothing),
              ("qʳ (g kg⁻¹)", r.bqr, r.sqr, 1e3, nothing), ("u (m s⁻¹)", r.bu, r.su, 1.0, nothing), ("v (m s⁻¹)", r.bv, r.svv, 1.0, nothing),
              ("radiative heating (K day⁻¹)", r.heat_b, r.heat_s, 1.0, nothing))
    for (k, (title, b, s, scale, xl)) in enumerate(panels)
        ax = Axis(fp[i, k], title = i == 1 ? title : "", ylabel = k == 1 ? "z (m), $(Int(r.h)) UTC" : "", limits = (xl, (0, 3000)))
        keep = isfinite.(s); lines!(ax, s[keep], r.zs[keep]; color = :gray30)
        keep = isfinite.(b); lines!(ax, scale .* b[keep], r.zb[keep]; color = :dodgerblue)
        if k == 7
            keep = isfinite.(r.heat_lw); lines!(ax, r.heat_lw[keep], r.zs[keep]; color = :gray60, linestyle = :dash)
        end
    end
end
Label(fp[0, :], "Profiles at selected hours: SAM samstat sample (grey; dashed = LW only in the heating panel) vs Breeze 10-min mean (blue)", fontsize = 13, tellwidth = false)
save(joinpath(output_dir, "step7_profiles.png"), fp; px_per_unit = 2)

stage("figures: plan views")
fv = Figure(size = (1800, 560 * length(plans) + 60), fontsize = 12)
for (i, p) in enumerate(plans)
    lwmax = max(quantile(vec(p.clwp), 0.99), quantile(vec(p.blwp), 0.99), 1.0)
    ax1 = Axis(fv[i, 1], title = @sprintf("SAM CLWP (g m⁻²), %.2f UTC", p.tsam), xlabel = "x (km)", ylabel = "y (km)", aspect = DataAspect())
    hm = heatmap!(ax1, p.x2, p.y2, p.clwp; colormap = :Blues, colorrange = (0, lwmax)); Colorbar(fv[i, 2], hm)
    ax2 = Axis(fv[i, 3], title = @sprintf("Breeze LWP (g m⁻²), %.2f UTC", p.tb), xlabel = "x (km)", aspect = DataAspect())
    nx, ny = size(p.blwp); xb = range(0.05, 25.6; length = nx); yb = range(0.05, 25.6; length = ny)
    hm2 = heatmap!(ax2, xb, yb, p.blwp; colormap = :Blues, colorrange = (0, lwmax)); Colorbar(fv[i, 4], hm2)
    pmax = max(quantile(vec(p.prec), 0.99), quantile(vec(p.brain), 0.99), 0.1)
    ax3 = Axis(fv[i, 5], title = "SAM Prec (mm day⁻¹)", xlabel = "x (km)", aspect = DataAspect())
    hm3 = heatmap!(ax3, p.x2, p.y2, p.prec; colormap = :viridis, colorrange = (0, pmax)); Colorbar(fv[i, 6], hm3)
    ax4 = Axis(fv[i, 7], title = "Breeze surface rain (mm day⁻¹)", xlabel = "x (km)", aspect = DataAspect())
    hm4 = heatmap!(ax4, xb, yb, p.brain; colormap = :viridis, colorrange = (0, pmax)); Colorbar(fv[i, 8], hm4)
end
Label(fv[0, :], "Plan views from sam2d (5-min instantaneous) and Breeze 2-D slices (30-min instantaneous): cloud liquid water path and surface rain", fontsize = 13, tellwidth = false)
save(joinpath(output_dir, "step7_planviews.png"), fv; px_per_unit = 2)

# ---------------------------------------------------------------- console summary
for (wname, d) in sort(collect(result["windows"]); by = first)
    @printf("%s UTC: CWP SAM %.1f Breeze %.1f | RWP %.3f / %.3f | PREC %.3f / %.3f mm/d | CLDSHD %.2f / %.2f | ZCB %.0f / %.0f m | ZCT %.0f / %.0f m | SAM SHF %.1f LHF %.1f LWNS %.1f SWNS %.1f\n",
            wname, d["CWP"]["sam"]["mean"], d["CWP"]["breeze"]["mean"], d["RWP"]["sam"]["mean"], d["RWP"]["breeze"]["mean"], d["PREC"]["sam"]["mean"], d["PREC"]["breeze"]["mean"],
            d["CLDSHD"]["sam"]["mean"], d["CLDSHD"]["breeze"]["mean"], d["ZCB_m"]["sam"]["mean"], d["ZCB_m"]["breeze"]["mean"], d["ZCT_m"]["sam"]["mean"], d["ZCT_m"]["breeze"]["mean"],
            d["SHF"]["sam"]["mean"], d["LHF"]["sam"]["mean"], d["LWNS"]["sam"]["mean"], d["SWNS"]["sam"]["mean"])
    o = obs[wname]
    haskey(o, "mwrret_lwp_g_m2") && @printf("   obs: MWRRET LWP %.1f (n %d) | VDIS %.3f mm/d | ARSCL top %.0f m, ceilometer base %.0f m, cloudy %.2f\n",
            o["mwrret_lwp_g_m2"]["mean"], o["mwrret_lwp_g_m2"]["n_good"], get(get(o, "vdis_rain_mm_day", Dict()), "mean", NaN),
            get(get(get(o, "arscl", Dict()), "lowest_layer_top_m", Dict()), "mean", NaN), get(get(get(o, "arscl", Dict()), "ceilometer_base_m", Dict()), "mean", NaN), get(get(o, "arscl", Dict()), "time_fraction_cloudy", NaN))
end
for (k, v) in sort(collect(profiles); by = first)
    @printf("profile %s: Breeze/SAM %s / %s  max qc %.2f / %.2f g/kg  max qr %.2e / %.2e  Nc(Breeze) %.0f cm⁻³  min heating %.1f / %.1f K/d  θ rms %.2f K  qv rms %.2f g/kg\n",
            k, v["breeze_time"], v["sam_time"], v["breeze_max_qc_gkg"], v["sam_max_qc_gkg"], v["breeze_max_qr_gkg"], v["sam_max_qr_gkg"], v["breeze_in_cloud_nc_cm3"],
            v["breeze_min_heating_K_day"], v["sam_min_heating_K_day"], v["theta_rms_below_2km_K"], v["qv_rms_below_2km_gkg"])
end
for (k, v) in sort(collect(plan_stats); by = first)
    @printf("plan %s: LWP mean SAM %.1f (std %.1f, >20: %.2f) Breeze %.1f (std %.1f, >20: %.2f); rain mean %.3f / %.3f mm/d, raining >0.1: %.2f / %.2f; SAM cloud top mean %.0f m\n",
            k, v["sam_clwp"]["mean"], v["sam_clwp_std"], v["sam_fraction_lwp_gt_20"], v["breeze_lwp"]["mean"], v["breeze_lwp_std"], v["breeze_fraction_lwp_gt_20"],
            v["sam_prec_mm_day"]["mean"], v["breeze_prec_mm_day"]["mean"], v["sam_fraction_raining_gt_0.1"], v["breeze_fraction_raining_gt_0.1"], v["sam_cloud_top_m"]["mean"])
end
println("wrote ", output_dir)
