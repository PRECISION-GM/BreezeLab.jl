# Compare a TRACER–DP-SCREAM Breeze run with the archived DP-SCREAM 3 km (August run) and 0.5 km
# (5–15 August cold start) domain means, the ARM VARANAL forcing analysis (its precipitation, LWP, PW and
# T/q/u/v state; domain means over the 150-km analysis domain) and ARM TRACER observations at the AMF1 site
# (rain gauge, MWR LWP/PWV, ARSCL hydrometeor layers, radiosondes), on identical UTC windows.
#
#   julia --project analysis/tracer_dp_scream_compare.jl <run dir> <data dir: IOP + DP-SCREAM files> <ARM dir> [<out dir>]
#
# Conventions (all printed into comparison.toml):
# - UTC throughout; CDT = UTC − 5 for diurnal composites and daily timing (CDT days 1–13 August).
# - Hourly series are means over (t − 1 h, t]; DP-SCREAM 30-min means and VARANAL hourly fluxes are
#   end-of-interval stamped (EAMxx averaged history; VARANAL "average over the hour leading up to").
# - Precipitation in mm day⁻¹: LES surface rain + ice flux; DP-SCREAM PRECL (unit inferred); VARANAL Prec
#   (MRMS-constrained domain mean); ARM tipping-bucket gauge (point).
# - Cloud fraction: LES fraction of cells with qᶜˡ + qⁱ > 10⁻⁵ kg kg⁻¹ (domain mean); DP-SCREAM TOT_CLOUD_FRAC
#   (model cloud fraction, domain mean); ARSCL hourly frequency of hydrometeor-layer occurrence at the site
#   (radar + lidar; includes precipitation; point) — different definitions, compared for structure and timing.
# - Profiles are sampled at the radiosonde launch times (nearest record) and averaged; RH is over liquid
#   (Bolton) for LES, VARANAL and computed consistently; the sondes report RH over liquid; DP-SCREAM RELHUM as archived.
using BreezeLab, Oceananigans, NCDatasets, CairoMakie, Statistics, Dates, TOML, Printf
using BreezeLab: saved_series, iop_datetimes

run_dir, data_dir, arm_dir = ARGS[1:3]
out_dir = length(ARGS) ≥ 4 ? ARGS[4] : joinpath(run_dir, "analysis")
mkpath(out_dir)
start = DateTime(2022, 8, 1); stop = DateTime(2022, 8, 15)
hours = collect(start + Hour(1):Hour(1):stop)                 # end-of-hour stamps
cdt(t) = t - Hour(5)
nanmean(x) = (y = filter(isfinite, x); isempty(y) ? NaN : mean(y))

# Mean of `values` over (h − Δ, h] for each stamp in `grid`
function binned(times, values, grid; Δ = Hour(1))
    order = sortperm(times); t = times[order]; v = values[order]
    map(grid) do h
        a = searchsortedfirst(t, h - Δ + Millisecond(1)); b = searchsortedlast(t, h)
        b ≥ a ? nanmean(v[a:b]) : NaN
    end
end
function binned_columns(times, M, grid; Δ = Hour(1))
    reduce(hcat, [begin
        sel = findall(t -> h - Δ < t ≤ h, times)
        isempty(sel) ? fill(NaN, size(M, 1)) : [nanmean(M[k, sel]) for k in axes(M, 1)]
    end for h in grid])
end
interp(z, v, zq) = BreezeLab.interpolate_profile(z, v, zq)
zgrid = collect(0.25:0.25:16.0) .* 1000

#####
##### Breeze LES
#####

config = TOML.parsefile(joinpath(run_dir, "provenance.toml"))["config"]
Oceananigans.defaults.FloatType = config["FT"] == "Float32" ? Float32 : Float64
case = tracer_dp_scream(; iop_path = joinpath(data_dir, "TRACER_iopfile_4scam.nc"), start, stop, arch = CPU(),
                        Nx = 8, Ny = 8, Δx = Float64(config["Lx"]) / 8, radiation = false, write_output = false)
grid = case.grid; Nz = size(grid, 3); Hz = grid.Hz
zc = Float64.(Array(znodes(grid, Center())))
pᵣ = Float64.(Array(interior(case.model.dynamics.reference_state.pressure, 1, 1, :)))
ts, s = saved_series(joinpath(run_dir, "tracer_timeseries.jld2"), (:lwp, :iwp, :rwp, :rain_flux, :ice_flux, :precipitable_water))
tl = start .+ Millisecond.(round.(Int, 1000 .* ts))
les = (; precip = binned(tl, 86400 .* (s.rain_flux .+ something(s.ice_flux, zero(s.rain_flux))), hours),
         lwp = binned(tl, 1e3 .* s.lwp, hours), iwp = binned(tl, 1e3 .* something(s.iwp, zero(s.lwp)), hours),
         pw = binned(tl, s.precipitable_water, hours))
tp, P = saved_series(joinpath(run_dir, "tracer_profiles.jld2"), (:T, :qᵛ, :total_cloud_fraction, :u, :v); interior_range = (Hz + 1):(Hz + Nz))
tpl = start .+ Millisecond.(round.(Int, 1000 .* tp))
esat(T) = 611.2 * exp(17.67 * (T - 273.15) / (T - 29.65))           # Bolton (1980), over liquid
rh(T, q, p) = 100 * (q * p / (0.622 + 0.378 * q)) / esat(T)
les_rh = [rh(P.T[k, n], P.qᵛ[k, n], pᵣ[k]) for k in 1:Nz, n in eachindex(tp)]
les_cf_hourly = binned_columns(tpl, P.total_cloud_fraction, hours)

#####
##### DP-SCREAM archive
#####

dp3 = read_dp_scream_output(joinpath(data_dir, "DP-SCREAM 3km_August.nc"))
dp05 = read_dp_scream_output(joinpath(data_dir, "DP-SCREAM 0.5km_August_5_to_August_15_run.nc"))
dpseries(o, x; scale = 1) = binned(o.time, scale .* x, hours)
dp = Dict(:dp3 => dp3, :dp05 => dp05)
dps = Dict(k => (; precip = dpseries(o, o.PRECL), lwp = dpseries(o, o.TGCLDLWP; scale = 1e3), iwp = dpseries(o, o.TGCLDIWP; scale = 1e3),
                   pw = haskey(o, :TMQ) ? dpseries(o, o.TMQ) : fill(NaN, length(hours))) for (k, o) in dp)
relhum_scale = maximum(filter(isfinite, dp3.RELHUM)) > 2 ? 1.0 : 100.0     # % or fraction
function dp_cf_on_grid(o)
    z = vec(mean(o.Z3, dims = 2)); order = sortperm(z)
    M = reduce(hcat, [interp(z[order], o.TOT_CLOUD_FRAC[order, n], zgrid) for n in eachindex(o.time)])
    binned_columns(o.time, M, hours)
end
dp_cf = Dict(k => dp_cf_on_grid(o) for (k, o) in dp)

#####
##### VARANAL (the forcing analysis)
#####

iop = read_iop_forcing(joinpath(data_dir, "TRACER_iopfile_4scam.nc"); profile_variables = (:T, :q, :u, :v),
                       surface_variables = (:Ps, :Prec, :LWP, :prew))
ti = iop_datetimes(iop)
va = (; precip = binned(ti, 86400 .* iop.surface.Prec, hours), lwp = binned(ti, 1e4 .* iop.surface.LWP, hours),
        pw = binned(ti, 10 .* iop.surface.prew, hours))
function varanal_heights(n)
    p = Float64.(iop.levels); ps = Float64(iop.surface.Ps[n])
    keep = findall(≤(ps), p); order = keep[sortperm(p[keep]; rev = true)]
    pk = vcat(ps, p[order]); Tv = (Float64.(iop.profiles.T[order, n]) .* (1 .+ 0.608 .* Float64.(iop.profiles.q[order, n])))
    Tv = vcat(Tv[1], Tv); z = zeros(length(pk))
    for k in 2:length(pk)
        z[k] = z[k - 1] + 287.04 / 9.81 * (Tv[k] + Tv[k - 1]) / 2 * log(pk[k - 1] / pk[k])
    end
    return order, z[2:end], p[order]
end

#####
##### ARM TRACER (AMF1 site)
#####

armfiles(stream) = sort(filter(f -> startswith(f, stream), readdir(arm_dir)))
function arm_series(stream, names)
    times = DateTime[]; vals = Dict(n => Float64[] for n in names)
    for f in armfiles(stream)
        NCDataset(joinpath(arm_dir, f)) do ds
            base = DateTime(1970) + Second(round(Int, ds["base_time"][]))
            t = base .+ Millisecond.(round.(Int, 1000 .* Float64.(nomissing(ds["time_offset"][:], NaN))))
            append!(times, t)
            for n in names
                x = Float64.(nomissing(ds[n][:], NaN)); x[x .< -9000] .= NaN
                append!(vals[n], x)
            end
        end
    end
    return times, vals
end
tm, met = arm_series("houmetM1.b1", ["tbrg_precip_total"])
gauge = binned(tm, met["tbrg_precip_total"], hours) .* 60 .* 24           # mm per minute → mm day⁻¹
tw, mwr = arm_series("houmwrret1liljclouM1.c2", ["be_lwp", "be_pwv"])
mwr["be_lwp"][(mwr["be_lwp"] .< -50) .| (mwr["be_lwp"] .> 3000)] .= NaN
arm = (; precip = gauge, lwp = binned(tw, mwr["be_lwp"], hours), pw = binned(tw, 10 .* mwr["be_pwv"], hours))
function arscl_occurrence()
    occ = zeros(length(zgrid), length(hours)); cnt = zeros(length(hours))
    for f in armfiles("houarsclkazr1kolliasM1.c0")
        NCDataset(joinpath(arm_dir, f)) do ds
            base = DateTime(1970) + Second(round(Int, ds["base_time"][]))
            t = base .+ Millisecond.(round.(Int, 1000 .* Float64.(ds["time_offset"][:])))
            B = Float64.(nomissing(ds["cloud_layer_base_height"][:, :], NaN)); Tt = Float64.(nomissing(ds["cloud_layer_top_height"][:, :], NaN))
            for n in eachindex(t)
                h = cld(Dates.value(t[n] - start), 3_600_000)          # hourly bin (h − 1 h, h]
                1 ≤ h ≤ length(hours) || continue
                cnt[h] += 1
                for l in axes(B, 1)
                    b = B[l, n]; top = Tt[l, n]
                    (isfinite(b) && isfinite(top) && b > -1000) || continue
                    occ[:, h] .+= (zgrid .≥ b) .& (zgrid .≤ top)
                end
            end
        end
    end
    return occ ./ max.(cnt', 1), cnt
end
arscl_cf, arscl_count = arscl_occurrence()
arscl_cf[:, arscl_count .== 0] .= NaN

sondes = map(armfiles("housondewnpnM1.b1")) do f
    NCDataset(joinpath(arm_dir, f)) do ds
        base = DateTime(1970) + Second(round(Int, ds["base_time"][]))
        z = Float64.(nomissing(ds["alt"][:], NaN)) .- 8.0
        T = Float64.(nomissing(ds["tdry"][:], NaN)) .+ 273.15; R = Float64.(nomissing(ds["rh"][:], NaN))
        good = isfinite.(z) .& isfinite.(T) .& isfinite.(R) .& (z .> 0)
        z, T, R = z[good], T[good], R[good]
        binmean(v) = [(sel = findall(x -> zz - 125 < x ≤ zz + 125, z); isempty(sel) ? NaN : mean(v[sel])) for zz in zgrid]
        (; launch = base, T = binmean(T), RH = binmean(R))
    end
end
filter!(s -> start ≤ s.launch < stop, sondes)

#####
##### Profiles at the sonde launch times
#####

nearest(times, t) = argmin(abs.(Dates.value.(times .- t)))
function sampled_profiles()
    out = Dict{Symbol, Vector{Vector{Float64}}}(k => Vector{Float64}[] for k in (:les_T, :les_RH, :les_u, :les_v, :dp3_T, :dp3_RH, :va_T, :va_RH, :va_u, :va_v, :obs_T, :obs_RH))
    for snd in sondes
        n = nearest(tpl, snd.launch)
        push!(out[:les_T], interp(zc, P.T[:, n], zgrid)); push!(out[:les_RH], interp(zc, les_rh[:, n], zgrid))
        push!(out[:les_u], interp(zc, P.u[:, n], zgrid)); push!(out[:les_v], interp(zc, P.v[:, n], zgrid))
        m = nearest(dp3.time, snd.launch); z = dp3.Z3[:, m]; o = sortperm(z)
        push!(out[:dp3_T], interp(z[o], dp3.T[o, m], zgrid)); push!(out[:dp3_RH], interp(z[o], relhum_scale .* dp3.RELHUM[o, m], zgrid))
        i = nearest(ti, snd.launch); order, zv, pv = varanal_heights(i)
        Tv = Float64.(iop.profiles.T[order, i]); qv = Float64.(iop.profiles.q[order, i])
        push!(out[:va_T], interp(zv, Tv, zgrid)); push!(out[:va_RH], interp(zv, rh.(Tv, qv, pv), zgrid))
        push!(out[:va_u], interp(zv, Float64.(iop.profiles.u[order, i]), zgrid)); push!(out[:va_v], interp(zv, Float64.(iop.profiles.v[order, i]), zgrid))
        push!(out[:obs_T], snd.T); push!(out[:obs_RH], snd.RH)
    end
    Dict(k => [nanmean([p[j] for p in v]) for j in eachindex(zgrid)] for (k, v) in out)
end
prof = sampled_profiles()

#####
##### Metrics
#####

window5 = findall(h -> h > DateTime(2022, 8, 5), hours)
valid(a, b) = findall(i -> isfinite(a[i]) && isfinite(b[i]), eachindex(a))
corr(a, b) = (v = valid(a, b); length(v) > 2 ? cor(a[v], b[v]) : NaN)
rmse(a, b) = (v = valid(a, b); sqrt(mean((a[v] .- b[v]) .^ 2)))
layer(prof_key, ref_key, z1, z2) = (sel = findall(z -> z1 ≤ z < z2, zgrid); nanmean(prof[prof_key][sel] .- prof[ref_key][sel]))

# CDT days 1–13 August (05 UTC to 05 UTC)
cdt_days = [DateTime(2022, 8, d, 5) for d in 1:13]
function daily_timing(x)
    map(cdt_days) do d0
        sel = findall(h -> d0 < h ≤ d0 + Day(1), hours)
        y = x[sel]
        all(isnan, y) && return (; total_mm = NaN, peak_hour_cdt = NaN)
        k = argmax(replace(y, NaN => -Inf))
        (; total_mm = sum(filter(isfinite, y)) / 24, peak_hour_cdt = Float64(Dates.hour(cdt(hours[sel][k])) + 0.5))
    end
end
circmean(h) = (v = filter(isfinite, h); isempty(v) ? NaN : mod(angle(mean(cis.(2π .* v ./ 24))) * 24 / (2π), 24))
circdiff(a, b) = (d = mod(a - b + 12, 24) - 12; abs(d))
precips = Dict("Breeze LES 200 m" => les.precip, "DP-SCREAM 3 km" => dps[:dp3].precip, "DP-SCREAM 0.5 km" => dps[:dp05].precip,
               "VARANAL (MRMS-constrained, domain)" => va.precip, "ARM gauge (site)" => arm.precip)
timing = Dict(k => daily_timing(v) for (k, v) in precips)
diurnal(x) = [nanmean(x[findall(h -> Dates.hour(cdt(h - Minute(30))) == H && h > cdt_days[1], hours)]) for H in 0:23]

summary = Dict{String, Any}()
summary["window"] = Dict("start" => string(start), "stop" => string(stop), "hours" => length(hours), "sondes" => length(sondes),
                         "arscl_hours" => count(>(0), arscl_count))
means(x; sel = eachindex(hours)) = nanmean(x[sel])
for (name, X) in (("les", les), ("dp3", dps[:dp3]), ("dp05", dps[:dp05]), ("varanal", va), ("arm", arm))
    d = Dict{String, Any}()
    for f in (:precip, :lwp, :iwp, :pw)
        haskey(X, f) || continue
        d["$(f)_mean_aug01_15"] = means(X[f]); d["$(f)_mean_aug05_15"] = means(X[f]; sel = window5)
    end
    summary[name] = d
end
summary["precip_hourly_correlation_vs_varanal"] = Dict("les" => corr(les.precip, va.precip), "dp3" => corr(dps[:dp3].precip, va.precip),
    "dp05_aug05_15" => corr(dps[:dp05].precip[window5], va.precip[window5]), "les_aug05_15" => corr(les.precip[window5], va.precip[window5]),
    "arm_gauge" => corr(arm.precip, va.precip), "les_vs_dp3" => corr(les.precip, dps[:dp3].precip))
summary["precip_hourly_rmse_vs_varanal_mm_day"] = Dict("les" => rmse(les.precip, va.precip), "dp3" => rmse(dps[:dp3].precip, va.precip))
daily_totals(k) = [t.total_mm for t in timing[k]]
summary["daily_total_correlation_vs_varanal"] = Dict(k => corr(daily_totals(k), daily_totals("VARANAL (MRMS-constrained, domain)")) for k in keys(timing))
summary["daily_peak_hour_cdt"] = Dict(k => [t.peak_hour_cdt for t in v] for (k, v) in timing)
summary["daily_total_mm"] = Dict(k => daily_totals(k) for k in keys(timing))
summary["mean_peak_hour_cdt_circular"] = Dict(k => circmean([t.peak_hour_cdt for t in v]) for (k, v) in timing)
ref_peaks = [t.peak_hour_cdt for t in timing["VARANAL (MRMS-constrained, domain)"]]
summary["mean_abs_peak_hour_difference_vs_varanal_h"] = Dict(k => nanmean([circdiff(a, b) for (a, b) in zip([t.peak_hour_cdt for t in v], ref_peaks)]) for (k, v) in timing)
for (ref, label) in ((:va, "varanal"), (:obs, "sondes"))
    summary["T_bias_K_vs_$label"] = Dict("les" => Dict("0-2km" => layer(:les_T, Symbol(ref, "_T"), 0, 2000), "2-6km" => layer(:les_T, Symbol(ref, "_T"), 2000, 6000), "6-12km" => layer(:les_T, Symbol(ref, "_T"), 6000, 12000)),
                                         "dp3" => Dict("0-2km" => layer(:dp3_T, Symbol(ref, "_T"), 0, 2000), "2-6km" => layer(:dp3_T, Symbol(ref, "_T"), 2000, 6000), "6-12km" => layer(:dp3_T, Symbol(ref, "_T"), 6000, 12000)))
    summary["RH_bias_pct_vs_$label"] = Dict("les" => Dict("0-2km" => layer(:les_RH, Symbol(ref, "_RH"), 0, 2000), "2-6km" => layer(:les_RH, Symbol(ref, "_RH"), 2000, 6000)),
                                            "dp3" => Dict("0-2km" => layer(:dp3_RH, Symbol(ref, "_RH"), 0, 2000), "2-6km" => layer(:dp3_RH, Symbol(ref, "_RH"), 2000, 6000)))
end
summary["T_bias_K_varanal_vs_sondes"] = Dict("0-2km" => layer(:va_T, :obs_T, 0, 2000), "2-6km" => layer(:va_T, :obs_T, 2000, 6000), "6-12km" => layer(:va_T, :obs_T, 6000, 12000))
wind_sel = findall(z -> z ≤ 12000, zgrid)
summary["mean_wind_les_minus_varanal_rms_0_12km"] = Dict("u" => sqrt(nanmean((prof[:les_u][wind_sel] .- prof[:va_u][wind_sel]) .^ 2)),
                                                           "v" => sqrt(nanmean((prof[:les_v][wind_sel] .- prof[:va_v][wind_sel]) .^ 2)))
cf_profile(M) = [nanmean(M[k, :]) for k in axes(M, 1)]
les_cf_grid = reduce(hcat, [interp(zc, les_cf_hourly[:, n], zgrid) for n in eachindex(hours)])
for (name, M) in (("les", les_cf_grid), ("dp3", dp_cf[:dp3]), ("arscl", arscl_cf))
    p = cf_profile(M); k = argmax(replace(p, NaN => -Inf))
    low = findall(z -> z < 3000, zgrid); mid = findall(z -> 3000 ≤ z < 8000, zgrid); high = findall(z -> z ≥ 8000, zgrid)
    summary["cloud_fraction_profile_$name"] = Dict("max" => p[k], "height_of_max_km" => zgrid[k] / 1000,
        "mean_below_3km" => nanmean(p[low]), "mean_3_8km" => nanmean(p[mid]), "mean_above_8km" => nanmean(p[high]))
end
open(io -> TOML.print(io, summary), joinpath(out_dir, "comparison.toml"), "w")

#####
##### Figures
#####

th = [Dates.value(h - start) / 3.6e6 for h in hours]
colors = Dict("les" => :black, "dp3" => :dodgerblue, "dp05" => :purple, "varanal" => :darkorange, "arm" => :gray50)
fig = Figure(size = (1100, 1100), fontsize = 12)
for (row, (field, label)) in enumerate(((:precip, "precipitation (mm day⁻¹)"), (:lwp, "LWP (g m⁻²)"), (:iwp, "IWP (g m⁻²)"), (:pw, "precipitable water (kg m⁻²)")))
    ax = Axis(fig[row, 1], ylabel = label, xlabel = row == 4 ? "hours since 2022-08-01 00 UTC" : "")
    lines!(ax, th, les[field], color = colors["les"], label = "Breeze LES 200 m (51 km)")
    lines!(ax, th, dps[:dp3][field], color = colors["dp3"], label = "DP-SCREAM 3 km (200 km)")
    lines!(ax, th, dps[:dp05][field], color = colors["dp05"], label = "DP-SCREAM 0.5 km")
    haskey(va, field) && lines!(ax, th, va[field], color = colors["varanal"], label = "VARANAL (domain)")
    haskey(arm, field) && lines!(ax, th, arm[field], color = (colors["arm"], 0.7), label = "ARM site")
    row == 1 && axislegend(ax, position = :lt, framevisible = false, nbanks = 2)
end
save(joinpath(out_dir, "timeseries_comparison.png"), fig)

fig = Figure(size = (1100, 1100), fontsize = 12)
for (row, (M, title)) in enumerate(((les_cf_grid, "Breeze LES 200 m: total cloud fraction (qᶜˡ + qⁱ > 10⁻⁵)"),
                                    (dp_cf[:dp3], "DP-SCREAM 3 km: TOT_CLOUD_FRAC"),
                                    (dp_cf[:dp05], "DP-SCREAM 0.5 km: TOT_CLOUD_FRAC"),
                                    (arscl_cf, "ARSCL (site): hourly hydrometeor occurrence")))
    ax = Axis(fig[row, 1], title = title, ylabel = "z (km)", xlabel = row == 4 ? "hours since 2022-08-01 00 UTC" : "")
    hm = heatmap!(ax, th, zgrid ./ 1e3, permutedims(M), colormap = :Blues, colorrange = (0, 1))
    row == 1 && Colorbar(fig[1:4, 2], hm)
end
save(joinpath(out_dir, "cloud_fraction_time_height_comparison.png"), fig)

fig = Figure(size = (1100, 800), fontsize = 12)
ax = Axis(fig[1, 1], xlabel = "hour (CDT)", ylabel = "precipitation (mm day⁻¹)", title = "Diurnal composite, CDT days 1–13 Aug")
for (k, x) in precips
    lines!(ax, 0.5:1:23.5, diurnal(x), label = k)
end
axislegend(ax, position = :lt, framevisible = false)
ax2 = Axis(fig[1, 2], xlabel = "hour (CDT)", ylabel = "LWP (g m⁻²)")
lines!(ax2, 0.5:1:23.5, diurnal(les.lwp), color = :black); lines!(ax2, 0.5:1:23.5, diurnal(dps[:dp3].lwp), color = :dodgerblue)
lines!(ax2, 0.5:1:23.5, diurnal(va.lwp), color = :darkorange); lines!(ax2, 0.5:1:23.5, diurnal(arm.lwp), color = :gray50)
for (col, (M, title)) in enumerate(((les_cf_grid, "LES"), (dp_cf[:dp3], "DP-SCREAM 3 km"), (arscl_cf, "ARSCL")))
    ax = Axis(fig[2, col > 2 ? 2 : 1][1, col > 2 ? 1 : col], title = "cloud fraction diurnal composite: $title", xlabel = "hour (CDT)", ylabel = "z (km)")
    D = reduce(hcat, [[nanmean(M[k, findall(h -> Dates.hour(cdt(h - Minute(30))) == H && h > cdt_days[1], hours)]) for k in axes(M, 1)] for H in 0:23])
    heatmap!(ax, 0.5:1:23.5, zgrid ./ 1e3, permutedims(D), colormap = :Blues, colorrange = (0, 0.6))
end
save(joinpath(out_dir, "diurnal_composites.png"), fig)

fig = Figure(size = (1100, 650), fontsize = 12)
axT = Axis(fig[1, 1], xlabel = "T − sondes (K)", ylabel = "z (km)", title = "mean at sonde launch times ($(length(sondes)) sondes)")
axR = Axis(fig[1, 2], xlabel = "RH (%)", ylabel = "z (km)")
axU = Axis(fig[1, 3], xlabel = "u (m s⁻¹)", ylabel = "z (km)")
for (key, label, c) in ((:les, "Breeze LES", :black), (:dp3, "DP-SCREAM 3 km", :dodgerblue), (:va, "VARANAL", :darkorange))
    lines!(axT, prof[Symbol(key, "_T")] .- prof[:obs_T], zgrid ./ 1e3, color = c, label = label)
    lines!(axR, prof[Symbol(key, "_RH")], zgrid ./ 1e3, color = c, label = label)
    haskey(prof, Symbol(key, "_u")) && lines!(axU, prof[Symbol(key, "_u")], zgrid ./ 1e3, color = c, label = label)
end
lines!(axR, prof[:obs_RH], zgrid ./ 1e3, color = :gray50, label = "sondes"); vlines!(axT, [0], color = :gray50)
axislegend(axR, position = :rt, framevisible = false)
save(joinpath(out_dir, "profiles_comparison.png"), fig)

@info "wrote $(joinpath(out_dir, "comparison.toml")) and four figures"
println(read(joinpath(out_dir, "comparison.toml"), String))
