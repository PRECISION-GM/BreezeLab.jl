# Resolution pair of TRACER–DP-SCREAM runs (e.g. 200 m job 202 vs 100 m job 416, same physics and forcing):
# domain-mean series and timing from each run's analysis/comparison.toml (analysis/tracer_dp_scream_compare.jl),
# plus convective-core statistics from the x–z sections at y = Ly/2 sampled at a common 30-min cadence.
#
#   julia --project analysis/tracer_dp_scream_resolution_pair.jl <coarse run dir> <fine run dir> <out dir> [labels "200 m,100 m"]
#
# Core: a cloudy updraft cell of the section, w > 2 m s⁻¹ (w averaged to cell centres) with qᶜˡ > 10⁻⁵ kg kg⁻¹.
# Statistics over all sampled sections of 1–15 August: the core fraction per height, the 99.9th percentile of w
# per height, the widths of contiguous cores along x at 5 km (periodic), and the daily maximum of w.
using JLD2, TOML, Statistics, Printf, CairoMakie

coarse, fine, out_dir = ARGS[1:3]
labels = length(ARGS) ≥ 4 ? split(ARGS[4], ",") : ["coarse", "fine"]
mkpath(out_dir)
cadence = 1800.0

function sections(run_dir)
    config = TOML.parsefile(joinpath(run_dir, "provenance.toml"))["config"]
    Nx, Nz = config["Nx"], config["Nz"]
    Δx = config["Lx"] / Nx
    file = joinpath(run_dir, "tracer_slices.jld2")
    jldopen(file) do f
        ts = f["timeseries"]
        iterations = sort(parse.(Int, collect(keys(ts["t"]))))
        times = [ts["t/$i"] for i in iterations]
        keep = findall(t -> isapprox(mod(t + 1, cadence), 1; atol = 2), times)   # multiples of 30 min
        H = (size(ts["qᶜˡ_xz/$(iterations[1])"], 1) - Nx) ÷ 2
        Hz = (size(ts["qᶜˡ_xz/$(iterations[1])"], 3) - Nz) ÷ 2
        grid_faces = nothing
        W = Array{Float32}(undef, Nx, Nz, length(keep)); Q = similar(W)
        for (m, n) in enumerate(keep)
            i = iterations[n]
            w = ts["w_xz/$i"][H+1:H+Nx, 1, Hz+1:Hz+Nz+1]
            W[:, :, m] .= (w[:, 1:Nz] .+ w[:, 2:Nz+1]) ./ 2
            Q[:, :, m] .= ts["qᶜˡ_xz/$i"][H+1:H+Nx, 1, Hz+1:Hz+Nz]
        end
        return (; times = times[keep], W, Q, Δx, Nx, Nz, config)
    end
end

zfaces(Nz) = begin
    # the case grid: 50 m to 2 km, stretched to 22 km (tracer_dp_scream_vertical_faces)
    N_uniform = 40; faces = collect(0.0:50.0:2000.0)
    target = (22000 - 2000) / 50; N = Nz - N_uniform
    lo, hi = 1.0, 2.0
    for _ in 1:200
        r = (lo + hi) / 2; sum(r^k for k in 1:N) < target ? (lo = r) : (hi = r)
    end
    r = (lo + hi) / 2
    for k in 1:N
        push!(faces, faces[end] + 50 * r^k)
    end
    faces[end] = 22000.0; faces
end

function core_statistics(s)
    zf = zfaces(s.Nz); zc = (zf[1:end-1] .+ zf[2:end]) ./ 2
    core = (s.W .> 2) .& (s.Q .> 1e-5)
    fraction = [mean(core[:, k, :]) for k in 1:s.Nz]
    w999 = [quantile(vec(s.W[:, k, :]), 0.999) for k in 1:s.Nz]
    k5 = argmin(abs.(zc .- 5000))
    widths = Float64[]
    for m in axes(core, 3)
        row = core[:, k5, m]
        all(row) && (push!(widths, s.Nx * s.Δx); continue)
        any(row) || continue
        start = findfirst(!, row)                     # walk the periodic row from a non-core cell
        run = 0
        for j in 0:s.Nx-1
            c = row[mod1(start + j, s.Nx)]
            if c
                run += 1
            elseif run > 0
                push!(widths, run * s.Δx); run = 0
            end
        end
        run > 0 && push!(widths, run * s.Δx)
    end
    days = floor.(Int, s.times ./ 86400)
    daily_wmax = [maximum(s.W[:, :, findall(==(d), days)]; init = -Inf) for d in 0:13]
    return (; zc, fraction, w999, widths, daily_wmax, k5_height = zc[k5],
              cores_per_section_at_5km = length(widths) / size(core, 3))
end

runs = [sections(coarse), sections(fine)]
stats = [core_statistics(r) for r in runs]
comparisons = [TOML.parsefile(joinpath(d, "analysis", "comparison.toml")) for d in (coarse, fine)]

summary = Dict{String, Any}()
for (label, r, st, c) in zip(labels, runs, stats, comparisons)
    sel(z1, z2) = findall(z -> z1 ≤ z < z2, st.zc)
    summary[label] = Dict(
        "Δx_m" => r.Δx, "sections" => length(r.times),
        "precip_mean_mm_day" => c["les"]["precip_mean_aug01_15"], "lwp_mean_g_m2" => c["les"]["lwp_mean_aug01_15"],
        "iwp_mean_g_m2" => c["les"]["iwp_mean_aug01_15"], "pw_mean_kg_m2" => c["les"]["pw_mean_aug01_15"],
        "precip_hourly_corr_vs_varanal" => c["precip_hourly_correlation_vs_varanal"]["les"],
        "daily_total_corr_vs_varanal" => c["daily_total_correlation_vs_varanal"]["Breeze LES"],
        "mean_peak_hour_cdt" => c["mean_peak_hour_cdt_circular"]["Breeze LES"],
        "mean_abs_peak_hour_diff_vs_varanal_h" => c["mean_abs_peak_hour_difference_vs_varanal_h"]["Breeze LES"],
        "daily_total_mm" => c["daily_total_mm"]["Breeze LES"],
        "daily_peak_hour_cdt" => c["daily_peak_hour_cdt"]["Breeze LES"],
        "cloud_fraction_below_3km" => c["cloud_fraction_profile_les"]["mean_below_3km"],
        "cloud_fraction_3_8km" => c["cloud_fraction_profile_les"]["mean_3_8km"],
        "cloud_fraction_above_8km" => c["cloud_fraction_profile_les"]["mean_above_8km"],
        "T_bias_vs_sondes_K" => c["T_bias_K_vs_sondes"]["les"], "RH_bias_vs_sondes_pct" => c["RH_bias_pct_vs_sondes"]["les"],
        "core_fraction_1_5km" => mean(st.fraction[sel(1000, 5000)]), "core_fraction_5_10km" => mean(st.fraction[sel(5000, 10000)]),
        "w999_max_m_s" => maximum(st.w999), "height_of_w999_max_km" => st.zc[argmax(st.w999)] / 1000,
        "core_width_mean_m_at_5km" => isempty(st.widths) ? NaN : mean(st.widths),
        "core_width_median_m_at_5km" => isempty(st.widths) ? NaN : median(st.widths),
        "cores_per_section_at_5km" => st.cores_per_section_at_5km,
        "daily_max_w_m_s" => st.daily_wmax)
end
summary["reference"] = Dict("dp3_precip_mean_mm_day" => comparisons[1]["dp3"]["precip_mean_aug01_15"],
                            "varanal_precip_mean_mm_day" => comparisons[1]["varanal"]["precip_mean_aug01_15"],
                            "dp3_mean_peak_hour_cdt" => comparisons[1]["mean_peak_hour_cdt_circular"]["DP-SCREAM 3 km"],
                            "varanal_mean_peak_hour_cdt" => comparisons[1]["mean_peak_hour_cdt_circular"]["VARANAL (MRMS-constrained, domain)"],
                            "dp3_cloud_fraction_3_8km" => comparisons[1]["cloud_fraction_profile_dp3"]["mean_3_8km"],
                            "arscl_cloud_fraction_3_8km" => comparisons[1]["cloud_fraction_profile_arscl"]["mean_3_8km"])
open(io -> TOML.print(io, summary), joinpath(out_dir, "resolution_pair.toml"), "w")

fig = Figure(size = (1100, 520), fontsize = 12)
ax1 = Axis(fig[1, 1], xlabel = "cloudy-updraft core fraction (w > 2 m s⁻¹, qᶜˡ > 10⁻⁵)", ylabel = "z (km)", xscale = log10)
ax2 = Axis(fig[1, 2], xlabel = "99.9th percentile of w (m s⁻¹)", ylabel = "z (km)")
ax3 = Axis(fig[1, 3], xlabel = "core width at 5 km (km)", ylabel = "cores")
for (label, st, c) in zip(labels, stats, (:black, :crimson))
    lines!(ax1, max.(st.fraction, 1e-6), st.zc ./ 1e3, color = c, label = label)
    lines!(ax2, st.w999, st.zc ./ 1e3, color = c, label = label)
    isempty(st.widths) || hist!(ax3, st.widths ./ 1e3, bins = 0:0.2:6, color = (c, 0.4), label = label)
end
ylims!(ax1, 0, 16); ylims!(ax2, 0, 16); axislegend(ax2, position = :rb, framevisible = false)
save(joinpath(out_dir, "core_statistics.png"), fig)
println(read(joinpath(out_dir, "resolution_pair.toml"), String))
