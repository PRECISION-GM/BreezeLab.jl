# SEA STARR aerosol–precipitation response across CTRL, N100 and N030: domain-mean evolution, drizzle state space
# and transition type, against the SEVIRI composite along the 40 trajectories and the statements of the SEA STARR
# intercomparison preprint (Diamond et al., EGUsphere 2026-5350).
#
#   julia --project analysis/sea_starr_members.jl --ctrl <dirs> --n100 <dirs> --n030 <dirs> [--out <dir>]
#
# Each run is a comma-separated list of restart-segment directories (segment 1 first; see sea_starr_segments.jl).
# Partial runs are fine: every quantity is reported up to the last record available.

using BreezeLab
using BreezeLab: interpolate_profile, driver_initial_profile
using Oceananigans
using Oceananigans.Fields: interior
using CairoMakie
using NCDatasets
using Statistics
using Printf
using Dates
using TOML

include(joinpath(@__DIR__, "sea_starr_segments.jl"))

function parse_args(args)
    opts = Dict{String, String}()
    i = 1
    while i ≤ length(args)
        startswith(args[i], "--") && i < length(args) || error("unexpected argument $(args[i])")
        opts[args[i][3:end]] = args[i + 1]
        i += 2
    end
    return opts
end
opts = parse_args(ARGS)
specs = [(:CTRL, get(opts, "ctrl", "runs/ctrl_full_66h_job211")),
         (:N100, get(opts, "n100", "runs/n100_full_66h_job214")),
         (:N030, get(opts, "n030", "runs/n030_full_66h_job218"))]
out_dir = get(opts, "out", "runs/members_analysis")
inputs = get(opts, "inputs", "/shared/home/greg/breezelab-runs/20261006/inputs/mip_sources/seastarr_22241697")
mkpath(out_dir)
commit = try strip(read(`git -C $(@__DIR__) rev-parse --short HEAD`, String)) catch; "unknown" end
CairoMakie.activate!(type = "png", px_per_unit = 1.5)

#####
##### SEVIRI composite and local solar time
#####

tr = NCDataset(joinpath(inputs, "SEA_STARR_Raw_Trajectories.nc"))
traj_hours = Float64.(tr["t"][:])
get_tr(n) = (a = Float64.(coalesce.(tr[n][:, :], NaN)); a[a .<= -900] .= NaN; a)
sev = Dict(n => get_tr(n) for n in ("cf", "lwp", "lwp_good", "nd_good", "height_good", "lon"))
close(tr)
sev_mean(n) = [(v = filter(isfinite, sev[n][i, :]); isempty(v) ? NaN : mean(v)) for i in axes(sev[n], 1)]
mean_lon = sev_mean("lon")
local_solar_hour(s) = mod(s / 3600 + 21 + interpolate_profile(traj_hours .* 3600, mean_lon, s) / 15, 24)
driver = read_dephy_driver(joinpath(inputs, "SEA_STARR_CTRL_SCM_driver.nc"))
w_sub(ζ) = interpolate_profile(driver.z, view(driver.forcing.w, :, 1), ζ)   # identical in all three drivers

#####
##### Per-member diagnostics
#####

value(f) = [f[n][1, 1, 1] for n in eachindex(f.times)]
column(f, n) = vec(Array(interior(f[n], 1, 1, :)))
smooth(x, n) = [mean(x[max(1, i - n):min(length(x), i + n)]) for i in eachindex(x)]

# Reference density profile (anelastic, as in the case) for number concentrations per volume
ρᵣ_profile = let
    grid = RectilinearGrid(CPU(), Float64; size = (8, 8, 288), x = (0, 400), y = (0, 400), z = sea_starr_vertical_faces(),
                           halo = (5, 5, 5), topology = (Periodic, Periodic, Bounded))
    rs = BreezeLab.Breeze.ReferenceState(grid, BreezeLab.Breeze.ThermodynamicConstants(Float64); base_pressure = driver.surface_pressure,
                                         potential_temperature = ζ -> driver_initial_profile(driver, :thetal, ζ),
                                         vapor_mass_fraction = ζ -> driver_initial_profile(driver, :qt, ζ))
    vec(Array(interior(rs.density, 1, 1, :)))
end

function member_diagnostics(name, spec)
    dirs = segment_dirs(spec)
    prefix = "sea_starr_" * lowercase(string(name))
    ts(n) = stitched_series(dirs, prefix * "_timeseries.jld2", n)
    st(n) = stitched_series(dirs, prefix * "_statistics.jld2", n)
    cwp_s = ts("cwp"); t1 = cwp_s.times
    cwp = 1e3 .* value(cwp_s); rwp = 1e3 .* value(ts("rwp")); cf = value(ts("cloud_fraction"))
    rain = 86400 .* value(ts("rain_flux")); zi = value(ts("zi")); nt_col = value(ts("n_total_column")); na_col = value(ts("nᵃ_column"))
    qc_s = st("qᶜˡ"); t15 = qc_s.times
    z = Array(znodes(qc_s.grid, Center())); Δz = diff(Array(znodes(qc_s.grid, Face())))
    n15 = length(t15)
    QC = hcat([column(qc_s, n) for n in 1:n15]...); QR = hcat([column(st("qʳ"), n) for n in 1:n15]...)
    NC = hcat([column(st("nᶜˡ"), n) for n in 1:n15]...); NA = hcat([column(st("nᵃ"), n) for n in 1:n15]...)
    CFP = hcat([column(st("cloud_fraction"), n) for n in 1:n15]...)
    ρ = ρᵣ_profile[1:length(z)]
    bin(x) = [mean(x[(t1 .> t - 900) .& (t1 .≤ t)]) for t in t15]
    lwp15 = bin(cwp) .+ bin(rwp); cf15 = bin(cf); rain15 = bin(rain); zi15 = bin(zi); rwp15 = bin(rwp)
    thr = 1e-5
    top15 = [(k = findlast(>(thr), QC[:, n]); isnothing(k) ? NaN : z[k]) for n in 1:n15]
    base15 = [(k = findfirst(>(thr), QC[:, n]); isnothing(k) ? NaN : z[k]) for n in 1:n15]
    # LWC-weighted in-cloud droplet number [cm⁻³] (preprint Fig. 5/7 definition, here from horizontal means)
    nd15 = map(1:n15) do n
        w = QC[:, n] .* ρ .* Δz
        sum(w) > 0 || return NaN
        sum(1e-6 .* NC[:, n] .* ρ ./ max.(CFP[:, n], 0.05) .* w) / sum(w)
    end
    rain_fraction15 = [(c = sum(ρ .* QC[:, n] .* Δz); r = sum(ρ .* QR[:, n] .* Δz); (c + r) > 0 ? r / (c + r) : NaN) for n in 1:n15]
    incloud_lwp15 = lwp15 ./ max.(cf15, 0.05)
    kbl(n) = max(2, something(findlast(≤(zi15[n]), z), 2))
    na_bl15 = [1e-6 * mean(NA[1:kbl(n), n]) for n in 1:n15]
    zi_s = smooth(zi, 30); dzidt = similar(zi_s)
    for i in eachindex(zi_s)
        a = max(1, i - 90); b = min(length(zi_s), i + 90)
        dzidt[i] = (zi_s[b] - zi_s[a]) / (t1[b] - t1[a])
    end
    we = dzidt .- w_sub.(zi_s)
    return (; name, dirs, t1, cwp, rwp, cf, rain, zi, nt_col, na_col, t15, lwp15, cf15, rain15, rwp15, zi15, top15, base15,
              nd15, rain_fraction15, incloud_lwp15, na_bl15, we, end_hours = t1[end] / 3600)
end

members = [member_diagnostics(name, spec) for (name, spec) in specs]

#####
##### Numbers
#####

numbers = Dict{String, Any}("commit" => commit, "generated" => string(now(UTC)))
day_ranges = (("day1", 0, 24), ("day2", 24, 48), ("day3", 48, 66.01))
for m in members
    k = string(m.name)
    numbers["$(k)_segments"] = segment_label(m.dirs)
    numbers["$(k)_end_hours"] = m.end_hours
    th = m.t15 ./ 3600
    for (lab, a, b) in day_ranges
        s = (th .≥ a) .& (th .< b)
        any(s) || continue
        numbers["$(k)_$(lab)_lwp"] = mean(m.lwp15[s]); numbers["$(k)_$(lab)_cf"] = mean(m.cf15[s])
        numbers["$(k)_$(lab)_rain_mm_d"] = mean(m.rain15[s]); numbers["$(k)_$(lab)_rwp"] = mean(m.rwp15[s])
        numbers["$(k)_$(lab)_nd_cm3"] = mean(filter(isfinite, m.nd15[s]))
        numbers["$(k)_$(lab)_na_bl_mg"] = mean(m.na_bl15[s])
        numbers["$(k)_$(lab)_rain_fraction"] = mean(filter(isfinite, m.rain_fraction15[s]))
    end
    numbers["$(k)_rain_max_mm_d"] = maximum(m.rain15); numbers["$(k)_rain_max_hour"] = th[argmax(m.rain15)]
    numbers["$(k)_cf_min"] = minimum(m.cf15); numbers["$(k)_cf_min_hour"] = th[argmin(m.cf15)]
    numbers["$(k)_cf_final"] = m.cf15[end]; numbers["$(k)_lwp_final"] = m.lwp15[end]
    numbers["$(k)_zi_final_m"] = m.zi15[end]; numbers["$(k)_zi_rise_m_per_h"] = (m.zi15[end] - m.zi15[1]) / (th[end] - th[1])
    numbers["$(k)_w_entrainment_mean_mm_s"] = 1e3 * mean(m.we)
    numbers["$(k)_nd_min_cm3"] = minimum(filter(isfinite, m.nd15)); numbers["$(k)_na_bl_min_mg"] = minimum(m.na_bl15)
    numbers["$(k)_aerosol_column_change_m2"] = m.nt_col[end] - m.nt_col[1]
    # first time the 15-min cloud fraction stays below 0.9 for ≥ 2 h
    below = m.cf15 .< 0.9
    first_break = nothing
    for n in eachindex(below)
        stop = findfirst(≥(m.t15[n] + 7200), m.t15)
        isnothing(stop) && break
        all(below[n:stop]) && (first_break = th[n]; break)
    end
    numbers["$(k)_first_sustained_cf_below_0p9_hour"] = something(first_break, "none")
    # transition type: drizzle-depletion if sustained surface rain ≥ 0.5 mm/d coincides with BL aerosol depletion
    # below 50 % of its maximum, otherwise deepening–warming (preprint terminology)
    drizzle = any(smooth(m.rain15, 4) .≥ 0.5)
    depleted = minimum(m.na_bl15[th .≥ 6]) < 0.5 * maximum(m.na_bl15)
    numbers["$(k)_transition_type"] = drizzle && depleted ? "drizzle-depletion" : drizzle ? "drizzling, aerosol not depleted (gradual)" : "deepening-warming"
end
for (lab, a, b) in day_ranges
    sv = (traj_hours .≥ a) .& (traj_hours .< b)
    numbers["SEVIRI_$(lab)_cf"] = mean(filter(isfinite, sev["cf"][sv, :]))
    numbers["SEVIRI_$(lab)_lwp_screened"] = mean(filter(isfinite, sev["lwp_good"][sv, :]))
    numbers["SEVIRI_$(lab)_nd_screened_cm3"] = mean(filter(isfinite, sev["nd_good"][sv, :]))
end

#####
##### Figures
#####

colors = Dict(:CTRL => Makie.wong_colors()[1], :N100 => Makie.wong_colors()[2], :N030 => Makie.wong_colors()[3])
hours_ticks(ax) = (ax.xticks = (collect(0:6:66), [@sprintf("%d\n%02d LST", h, round(Int, local_solar_hour(h * 3600))) for h in 0:6:66]))
fig = Figure(size = (1150, 1450), fontsize = 12)
Label(fig[0, 1:2], "SEA STARR CTRL / N100 / N030 (15-min means; segments stitched at the restore time) vs SEVIRI composite", font = :bold, tellwidth = false)
panels = ((1, 1, "LWP cloud + rain (g m⁻²)", m -> m.lwp15, "lwp"), (1, 2, "Cloud fraction (LWP > 5 g m⁻²)", m -> m.cf15, "cf"),
          (2, 1, "Surface rain (mm day⁻¹)", m -> m.rain15, nothing), (2, 2, "Rain fraction of liquid water", m -> m.rain_fraction15, nothing),
          (3, 1, "Mean inversion height (m)", m -> m.zi15, nothing), (3, 2, "LWC-weighted in-cloud N_d (cm⁻³)", m -> m.nd15, "nd_good"),
          (4, 1, "Aerosol below zᵢ (mg⁻¹)", m -> m.na_bl15, nothing), (4, 2, "In-cloud LWP (g m⁻²)", m -> m.incloud_lwp15, "lwp_good"))
for (r, c, label, getter, sevkey) in panels
    ax = Axis(fig[r, c], ylabel = label); hours_ticks(ax)
    for m in members
        lines!(ax, m.t15 ./ 3600, getter(m), color = colors[m.name], label = string(m.name))
    end
    isnothing(sevkey) || scatter!(ax, traj_hours, sev_mean(sevkey), color = :black, markersize = 5, label = "SEVIRI")
    (r, c) == (1, 1) && axislegend(ax, position = :lt, labelsize = 10)
    label == "Surface rain (mm day⁻¹)" && (ax.yscale = Makie.pseudolog10)
end
save(joinpath(out_dir, "members_timeseries.png"), fig)

# Drizzle state space (preprint Fig. 7): hourly markers, colour = rain fraction; adiabatic cloud-top rₑ isolines
cw = 2.0e-6; ρw = 1000.0; kdisp = 0.8     # adiabatic condensate gradient [kg m⁻⁴], water density, dispersion factor
lwp_for_re(N, re) = (4π / 3 * ρw * kdisp * N * re^3)^2 / (2cw) * 1e3     # g m⁻² for N in m⁻³, rₑ in m
fig2 = Figure(size = (1200, 430), fontsize = 12)
Label(fig2[0, 1:4], "Fraction of liquid in rain vs LWC-weighted N_d and in-cloud LWP (hourly); lines: adiabatic cloud-top rₑ = 12, 14, 16 µm", font = :bold, tellwidth = false)
Nd_axis = 10 .^ range(log10(5), log10(800), length = 100)
for (i, m) in enumerate(members)
    ax = Axis(fig2[1, i], title = string(m.name), xscale = log10, yscale = log10, xlabel = "N_d (cm⁻³)", ylabel = i == 1 ? "in-cloud LWP (g m⁻²)" : "")
    for (re, ls) in ((12e-6, :solid), (14e-6, :dash), (16e-6, :dot))
        lines!(ax, Nd_axis, lwp_for_re.(Nd_axis .* 1e6, re), color = :black, linestyle = ls)
    end
    hourly = 1:4:length(m.t15)
    ok = [n for n in hourly if isfinite(m.nd15[n]) && m.nd15[n] > 0 && m.incloud_lwp15[n] > 1]
    sc = scatter!(ax, m.nd15[ok], m.incloud_lwp15[ok], color = m.rain_fraction15[ok], colormap = :viridis, colorrange = (0, 0.3),
                  marker = (:circle, :rect, :utriangle)[i], markersize = 8)
    xlims!(ax, 5, 800); ylims!(ax, 5, 1000)
    i == 3 && Colorbar(fig2[1, 4], sc, label = "rain / (cloud + rain)")
end
save(joinpath(out_dir, "members_state_space.png"), fig2)

#####
##### Summary
#####

open(joinpath(out_dir, "members_summary.toml"), "w") do io
    TOML.print(io, Dict(k => (v isa Tuple || v isa Vector ? string(v) : v) for (k, v) in numbers))
end
open(joinpath(out_dir, "members_summary.md"), "w") do io
    println(io, "# SEA STARR CTRL / N100 / N030 summary (commit ", commit, ", ", Dates.format(now(UTC), "yyyy-mm-dd HH:MM"), " UTC)\n")
    for m in members
        println(io, "- ", m.name, ": ", segment_label(m.dirs), @sprintf(" (to t = %.2f h)", m.end_hours))
    end
    println(io, "\n| quantity | CTRL | N100 | N030 | SEVIRI |\n| --- | --- | --- | --- | --- |")
    fmt(v) = v isa AbstractFloat ? @sprintf("%.3g", v) : string(v)
    rows = ["day1_lwp", "day2_lwp", "day3_lwp", "day1_cf", "day2_cf", "day3_cf", "day1_rain_mm_d", "day2_rain_mm_d", "day3_rain_mm_d",
            "day1_rwp", "day2_rwp", "day3_rwp", "day1_rain_fraction", "day2_rain_fraction", "day3_rain_fraction",
            "day1_nd_cm3", "day2_nd_cm3", "day3_nd_cm3", "day1_na_bl_mg", "day2_na_bl_mg", "day3_na_bl_mg",
            "rain_max_mm_d", "rain_max_hour", "cf_min", "cf_min_hour", "cf_final", "lwp_final", "zi_final_m", "zi_rise_m_per_h",
            "w_entrainment_mean_mm_s", "nd_min_cm3", "na_bl_min_mg", "first_sustained_cf_below_0p9_hour", "transition_type", "end_hours"]
    sevmap = Dict("day1_lwp" => "SEVIRI_day1_lwp_screened", "day2_lwp" => "SEVIRI_day2_lwp_screened", "day3_lwp" => "SEVIRI_day3_lwp_screened",
                  "day1_cf" => "SEVIRI_day1_cf", "day2_cf" => "SEVIRI_day2_cf", "day3_cf" => "SEVIRI_day3_cf",
                  "day1_nd_cm3" => "SEVIRI_day1_nd_screened_cm3", "day2_nd_cm3" => "SEVIRI_day2_nd_screened_cm3", "day3_nd_cm3" => "SEVIRI_day3_nd_screened_cm3")
    for r in rows
        vals = [fmt(get(numbers, "$(m.name)_$(r)", "—")) for m in members]
        s = haskey(sevmap, r) ? fmt(get(numbers, sevmap[r], "—")) : ""
        println(io, "| ", r, " | ", join(vals, " | "), " | ", s, " |")
    end
end
println("members analysis written to ", out_dir)
for k in sort(collect(keys(numbers)))
    println(k, " = ", numbers[k])
end
