# Compare the ENA-covert runs with the numbers Covert, Mechem & Zhang (2022, ACP 22, 1159)
# report for their 09:00–12:00 UTC analysis window of 18 July 2017. Usage:
#   julia --project analysis/covert_paper_comparison.jl OUTPUT_DIR SFC_FILE LABEL=RUN_DIR [LABEL=RUN_DIR ...]
# Every statistic is computed over the same 09–12 UTC window (model hours 3–6) from the
# saved hourly-mean profiles (θ, qᵛ, qᶜˡ, qʳ, w²) and the hourly x–z slices, with the
# paper's cloud definition (qᶜˡ ≥ 0.01 g kg⁻¹); the surface fluxes of the covert protocol
# are prescribed inputs (SFC_FLX_FXD), so they are read from the `sfc` file. Writes
# paper_comparison.toml, mean_profiles.png, rain_sections.png and nc_sensitivity.png.
# All statistics are computed before any Makie object is created: with Makie axes alive in the
# session, the column loop (`findall(≥(t), view(...))`) stalled for >1 h (reproduced in isolation
# on the CPU node, 8 Oct 2026); the same loop takes milliseconds when it runs first.
using BreezeLab, CairoMakie, JLD2, Oceananigans, Statistics, TOML, Dates, Printf
using BreezeLab: interpolate_profile
using Oceananigans.Grids: Center, Face, znodes
using Oceananigans.Fields: interior

stage(msg) = (println(stderr, "[", Dates.format(Dates.now(), "HH:MM:SS"), "] ", msg); flush(stderr))
stage("loaded")
output_dir, sfc_file = ARGS[1], ARGS[2]
runs = [(split(a, "=", limit=2)[1], split(a, "=", limit=2)[2]) for a in ARGS[3:end]]
mkpath(output_dir)

const PAPER = (; cloud_base = 821.0, cloud_top = 1109.0, inversion = 1132.5, inversion_hPa = 895.0,
                 peak_lwc_gkg = 0.5, shf = 11.8, lhf = 105.8, window = "09:00–12:00 UTC",
                 base_shift_per_25 = 10.0, top_shift_per_25 = 25.0,
                 mixed_layer_theta_l = 292.2, mixed_layer_qt = 11.2,
                 cloud_threshold_gkg = 0.01)   # "points having qc of 0.01 g kg⁻¹ or greater"
window = (3 * 3600.0, 6 * 3600.0)             # seconds since 06 UTC
const QC_THRESHOLD = 1e-3 * PAPER.cloud_threshold_gkg   # kg/kg

# prescribed surface fluxes (sfc file: day, SST, H, LE, tau) averaged over the window
sfc = read_sam_surface_forcing(sfc_file)
t_sfc = day_to_seconds.(sfc.day, 199.25)
window_mean(ts, values) = begin
    tt = range(window[1], window[2]; length = 361)
    mean(interpolate_profile(ts, values, t) for t in tt)
end
shf_in, lhf_in = window_mean(t_sfc, sfc.sensible_heat_flux), window_mean(t_sfc, sfc.latent_heat_flux)
stage("fluxes done")

col(f, n) = vec(Array(interior(f[n])))
in_window(t) = window[1] < t ≤ window[2] + 1    # hourly means ending at 4, 5, 6 h
function window_profiles(dir)
    file = only(filter(f -> endswith(f, "_profiles.jld2"), readdir(dir; join = true)))
    names = ["θ", "qᵛ", "qᶜˡ", "qʳ", "cloud_fraction", "w²", "T", "nᶜˡ", "nᵃ"]
    have(n) = jldopen(f -> haskey(f["timeseries"], n), file)
    fts = Dict(n => FieldTimeSeries(file, n) for n in names if have(n))
    θ = fts["θ"]; z = collect(znodes(θ.grid, Center())); zf = collect(znodes(θ.grid, Face()))
    idx = [n for n in eachindex(θ.times) if in_window(θ.times[n])]
    isempty(idx) && error("$dir has no hourly profile inside 09–12 UTC")
    mean_profile(name) = haskey(fts, name) ? mean(col(fts[name], n) for n in idx) : nothing
    qr = fts["qʳ"]
    section = (; hours = qr.times ./ 3600, qʳ = hcat([col(qr, n) for n in eachindex(qr.times)]...))   # (z, t) hourly means, full run
    return (; z, zf, hours = θ.times[idx] ./ 3600, θ = mean_profile("θ"), qᵛ = mean_profile("qᵛ"), qᶜˡ = mean_profile("qᶜˡ"),
              qʳ = mean_profile("qʳ"), cf = mean_profile("cloud_fraction"), w² = mean_profile("w²"), T = mean_profile("T"),
              nᶜˡ = (haskey(fts, "nᶜˡ") ? mean_profile("nᶜˡ") : nothing), section)
end
# column-wise ("ARSCL-like") cloud base/top from the hourly x–z slices: per column the lowest/highest
# cell with qᶜˡ ≥ threshold, averaged over cloudy columns and the window's slice times. The slice
# arrays are read raw from JLD2 (halo of 5 stripped) with the heights of the profile grid `z`.
function column_boundaries(dir, z; threshold = QC_THRESHOLD, halo = 5)
    file = only(filter(f -> endswith(f, "_slices.jld2"), readdir(dir; join = true)))
    nz = length(z)
    # read the window's slices first (no counters inside the closure), then count
    slices = jldopen(file) do f
        ts = f["timeseries"]
        iterations = sort(parse.(Int, filter(!=("serialized"), collect(keys(ts["qᶜˡ_xz"])))))
        inside = [it for it in iterations if in_window(ts["t/$it"])]
        [Array{Float64}(ts["qᶜˡ_xz/$it"][halo+1:end-halo, 1, halo+1:halo+nz]) for it in inside]
    end
    bases = Float64[]; tops = Float64[]; total = 0
    for data in slices, i in axes(data, 1)
        total += 1
        ks = findall(≥(threshold), view(data, i, :))
        isempty(ks) && continue
        push!(bases, z[first(ks)]); push!(tops, z[last(ks)])
    end
    return (; base = isempty(bases) ? NaN : mean(bases), top = isempty(tops) ? NaN : mean(tops),
              column_cloud_fraction = total == 0 ? NaN : length(bases) / total, n_columns = total)
end
inversion_height(θ, z) = (dθ = diff(θ) ./ diff(z); k = argmax(dθ); (z[k] + z[k+1]) / 2)
function nc_estimate(dir, p)
    prov = TOML.parsefile(joinpath(dir, "provenance.toml"))["config"]
    if isnothing(p.nᶜˡ)
        return (; label = "prescribed", value_cm3 = get(prov, "droplet_number", NaN) / 1e6)
    end
    ρ0 = get(prov, "reference_density", 1.206)
    cloudy = p.cf .> 0.05
    return (; label = "in-cloud mean", value_cm3 = 1e-6 * ρ0 * mean((p.nᶜˡ ./ max.(p.cf, 1e-6))[cloudy]))
end

results = Dict{String, Any}("paper" => Dict(string(k) => v for (k, v) in pairs(PAPER)),
                            "prescribed_fluxes_W_m2" => Dict("sensible" => shf_in, "latent" => lhf_in, "source" => sfc_file,
                                                             "note" => "SFC_FLX_FXD inputs averaged over 09–12 UTC; the runs save no flux output"))
colors = Makie.wong_colors()
sections = Any[]
profiles_to_plot = Any[]
summary_rows = String[]
for (i, (label, dir)) in enumerate(runs)
    # a run is usable when its saved time series reaches the end of the window (the docs run has no COMPLETE marker)
    tsfile = filter(f -> endswith(f, "_timeseries.jld2"), readdir(dir; join = true))
    (length(tsfile) == 1 && last(FieldTimeSeries(only(tsfile), "lwp").times) ≥ window[2] - 1) ||
        (println("skipping $label: $dir has no time series reaching $(window[2]) s"); continue)
    stage("$label: profiles")
    p = window_profiles(dir); c = colors[mod1(i, 7)]
    cloudy = p.qᶜˡ .≥ QC_THRESHOLD                       # the paper's definition applied to the mean profile
    base = any(cloudy) ? p.z[findfirst(cloudy)] : NaN; top = any(cloudy) ? p.z[findlast(cloudy)] : NaN
    stage("$label: columns")
    cb = column_boundaries(dir, p.z)
    stage("$label: inversion")
    zi = inversion_height(p.θ, p.z)
    kpeak = argmax(p.qᶜˡ)
    nc = nc_estimate(dir, p)
    qt = 1e3 .* (p.qᵛ .+ p.qᶜˡ)
    θl = p.θ   # the saved θ is Breeze's liquid-ice potential temperature (θˡⁱ), SAM's θl analogue
    # rain: window-mean profile, its maximum, and the rain water path series
    stage("$label: time series")
    tsf = only(filter(f -> endswith(f, "_timeseries.jld2"), readdir(dir; join = true)))
    rwp = FieldTimeSeries(tsf, "rwp"); rain = FieldTimeSeries(tsf, "rain_flux")
    rwp_series = [1e3 * Array(interior(rwp[n]))[1] for n in eachindex(rwp.times)]
    rain_series = [3600 * Array(interior(rain[n]))[1] for n in eachindex(rain.times)]
    win = [in_window(t) for t in rwp.times]
    kr = argmax(p.qʳ)
    r = Dict("run" => dir, "cloud_base_m_profile" => base, "cloud_top_m_profile" => top,
             "cloud_base_m_columns" => cb.base, "cloud_top_m_columns" => cb.top, "column_cloud_fraction" => cb.column_cloud_fraction,
             "inversion_m" => zi, "peak_lwc_gkg" => 1e3 * p.qᶜˡ[kpeak], "peak_lwc_height_m" => p.z[kpeak],
             "nc_cm3" => nc.value_cm3, "nc_definition" => nc.label, "hours_used" => collect(p.hours),
             "qt_below_cloud_gkg" => mean(qt[p.z .< 600]), "theta_l_below_cloud_K" => mean(θl[p.z .< 600]),
             "max_w_variance_m2s2" => maximum(p.w²), "max_w_variance_height_m" => p.zf[argmax(p.w²)],
             "max_qr_gkg" => 1e3 * p.qʳ[kr], "max_qr_height_m" => p.z[kr], "qr_surface_gkg" => 1e3 * p.qʳ[1],
             "rwp_g_m2_window_mean" => mean(rwp_series[win]), "rwp_g_m2_max" => maximum(rwp_series),
             "rain_mm_hr_window_mean" => mean(rain_series[win]), "rain_mm_hr_max" => maximum(rain_series))
    results[label] = r
    push!(profiles_to_plot, (; label, color = c, p.z, p.zf, qc = 1e3 .* p.qᶜˡ, qt, θl, qr = max.(1e3 .* p.qʳ, 1e-7), w2 = p.w²))
    push!(sections, (; label, color = c, p.section, z = p.z, rwp_hours = rwp.times ./ 3600, rwp = rwp_series))
    push!(summary_rows, @sprintf("%-22s base %4.0f/%4.0f top %4.0f/%4.0f (profile/columns) zi %4.0f peakLWC %.2f g/kg @%4.0f m  Nc %5.0f cm⁻³ (%s)  w²max %.3f  qr max %.2e g/kg @%4.0f m  RWP %.2f g/m²  rain %.4f mm/hr",
                                 label, base, cb.base, top, cb.top, zi, 1e3 * p.qᶜˡ[kpeak], p.z[kpeak], nc.value_cm3, nc.label[1:min(end, 9)], maximum(p.w²),
                                 1e3 * p.qʳ[kr], p.z[kr], mean(rwp_series[win]), mean(rain_series[win])))
end
stage("profile figure")
fig = Figure(size = (1800, 700), fontsize = 13)
titles = ("qᶜˡ (g kg⁻¹)", "qᵗ = qᵛ + qᶜˡ (g kg⁻¹)", "θˡ (K)", "qʳ (g kg⁻¹, log)", "w variance (m² s⁻²)")
xlims = (nothing, (0, 13), (288, 312), nothing, nothing)
axes = [Axis(fig[1, k], title = t, ylabel = k == 1 ? "z (m)" : "", limits = (xlims[k], (0, 1600)),
             xscale = (k == 4 ? log10 : identity)) for (k, t) in enumerate(titles)]
for ax in axes
    hlines!(ax, [PAPER.cloud_base, PAPER.cloud_top]; color = :black, linestyle = :dot)
    hlines!(ax, [PAPER.inversion]; color = :black, linestyle = :dash)
end
vlines!(axes[1], [PAPER.peak_lwc_gkg]; color = :black, linestyle = :dot)
vlines!(axes[2], [PAPER.mixed_layer_qt]; color = :black, linestyle = :dot); vlines!(axes[3], [PAPER.mixed_layer_theta_l]; color = :black, linestyle = :dot)
for q in profiles_to_plot
    lines!(axes[1], q.qc, q.z; color = q.color, label = q.label)
    lines!(axes[2], q.qt, q.z; color = q.color); lines!(axes[3], q.θl, q.z; color = q.color)
    lines!(axes[4], q.qr, q.z; color = q.color)
    lines!(axes[5], q.w2, q.zf; color = q.color)
end
axislegend(axes[1], position = :rt, labelsize = 9)
Label(fig[0, :], "ENA-covert runs, 09–12 UTC means, vs Covert et al. (2022) stated values (dotted: cloud base 821 m / top 1109 m, peak LWC 0.5 g kg⁻¹, mixed-layer qᵗ 11.2 g kg⁻¹ and θˡ 292.2 K; dashed: inversion 1132.5 m)", fontsize = 13, tellwidth = false)
save(joinpath(output_dir, "mean_profiles.png"), fig; px_per_unit = 2)

stage("rain sections")
# rain mass fraction time–height sections (hourly-mean profiles, the run's saved cadence) with the RWP series, one row per member
nsec = length(sections)
fs = Figure(size = (1800, 320 * nsec + 80), fontsize = 13)
qr_floor, qr_ceiling = 1e-7, max(1e-6, maximum(maximum(1e3 .* s.section.qʳ) for s in sections))
for (i, s) in enumerate(sections)
    ax = Axis(fs[i, 1], title = "$(s.label): rain mass fraction qʳ (g kg⁻¹), hourly means", ylabel = "z (m)", xlabel = i == nsec ? "hours since 06 UTC" : "", limits = (nothing, (0, 1600)))
    hm = heatmap!(ax, s.section.hours, s.z, max.(1e3 .* s.section.qʳ', qr_floor); colormap = :viridis, colorscale = log10, colorrange = (qr_floor, qr_ceiling))
    i == 1 && Colorbar(fs[1:nsec, 2], hm; label = "qʳ (g kg⁻¹)")
    ax2 = Axis(fs[i, 3], title = "rain water path (g m⁻²)", xlabel = i == nsec ? "hours since 06 UTC" : "")
    lines!(ax2, s.rwp_hours, s.rwp; color = s.color)
    vspan!(ax2, [window[1] / 3600], [window[2] / 3600]; color = (:gray, 0.15))
end
colsize!(fs.layout, 1, Relative(0.62))
Label(fs[0, :], "Rain in the ENA-covert runs (shared log scale); the paper's Fig. 6 shows instantaneous x–z sections of qc and qr, with drizzle described in the text (see docs)", fontsize = 13, tellwidth = false)
save(joinpath(output_dir, "rain_sections.png"), fs; px_per_unit = 2)

stage("nc sensitivity")
# Nc sensitivity across the P3 members (prescribed or diagnosed in-cloud Nc)
pts = [(results[l]["nc_cm3"], results[l]["cloud_base_m_profile"], results[l]["cloud_top_m_profile"], l) for (l, _) in runs
       if haskey(results, l) && isfinite(results[l]["nc_cm3"]) && !startswith(l, "one_moment")]   # P3 members only
if length(pts) ≥ 2
    fn = Figure(size = (900, 500), fontsize = 13)
    ax = Axis(fn[1, 1], xlabel = "droplet number Nc (cm⁻³)", ylabel = "height (m)", title = "Cloud base/top (mean qᶜˡ ≥ 0.01 g kg⁻¹) vs Nc, 09–12 UTC; paper: base +10 m, top +25 m per +25 cm⁻³")
    scatter!(ax, [p[1] for p in pts], [p[2] for p in pts]; color = :dodgerblue, label = "base")
    scatter!(ax, [p[1] for p in pts], [p[3] for p in pts]; color = :orange, label = "top")
    for p in pts
        text!(ax, p[1], p[3]; text = p[4], fontsize = 9, offset = (3, 3))
    end
    nc = [p[1] for p in pts]; x = range(minimum(nc), maximum(nc); length = 2)
    lines!(ax, x, PAPER.cloud_base .+ PAPER.base_shift_per_25 / 25 .* (x .- 75); color = :dodgerblue, linestyle = :dot, label = "paper base slope through 821 m")
    lines!(ax, x, PAPER.cloud_top .+ PAPER.top_shift_per_25 / 25 .* (x .- 75); color = :orange, linestyle = :dot, label = "paper top slope through 1109 m")
    axislegend(ax, position = :lt, labelsize = 9)
    save(joinpath(output_dir, "nc_sensitivity.png"), fn; px_per_unit = 2)
    sorted = sort(pts; by = first)
    slope(y) = (A = hcat(ones(length(nc)), nc); c = A \ y; 25 * c[2])      # least-squares m per 25 cm⁻³
    results["nc_sensitivity"] = Dict("members" => [p[4] for p in sorted], "nc_cm3" => [p[1] for p in sorted],
                                     "base_m" => [p[2] for p in sorted], "top_m" => [p[3] for p in sorted],
                                     "base_shift_per_25_cm3" => slope([p[2] for p in pts]),
                                     "top_shift_per_25_cm3" => slope([p[3] for p in pts]), "fit" => "least squares over the P3 members")
end
open(io -> TOML.print(io, results), joinpath(output_dir, "paper_comparison.toml"), "w")
println("prescribed fluxes 09–12 UTC: H = ", round(shf_in; digits = 1), " W/m², LE = ", round(lhf_in; digits = 1), " W/m² (paper: 11.8 / 105.8)")
foreach(println, summary_rows)
haskey(results, "nc_sensitivity") && println("Nc sensitivity (least squares, P3 members): base ", round(results["nc_sensitivity"]["base_shift_per_25_cm3"]; digits = 1), " m, top ", round(results["nc_sensitivity"]["top_shift_per_25_cm3"]; digits = 1), " m per 25 cm⁻³ (paper 10 / 25)")
println("wrote ", output_dir)
