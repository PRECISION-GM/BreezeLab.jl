# Full analysis of a completed SEA STARR run: workbook-style diagnostics over the whole run,
# comparison with the SEVIRI retrievals along the composite trajectories, and the checks behind
# the breakup verdict (entrainment vs subsidence, subsidence as applied vs the driver, nudging
# mask, radiation at the surface and the LES top, cloud-top radiative cooling).
#
#   julia --project analysis/sea_starr_ctrl_analysis.jl --run runs/ctrl_full_66h_job211 [--member CTRL] [--out <run>/analysis]
#
# Writes figures (PNG) and `summary.md` / `summary.toml` with the numbers under --out.

using BreezeLab
using BreezeLab: interpolate_profile
using Breeze
using Oceananigans
using Oceananigans.Units
using Oceananigans.Fields: interior
using Oceananigans.Units: Time
using CairoMakie
using NCDatasets
using JLD2
using Statistics
using Printf
using Dates
using TOML

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
run_spec = get(opts, "run", "runs/ctrl_full_66h_job211")   # comma-separated segment directories, segment 1 first
include(joinpath(@__DIR__, "sea_starr_segments.jl"))
run_dirs = segment_dirs(run_spec)
run_dir = first(run_dirs)
member = Symbol(get(opts, "member", "CTRL"))
out_dir = get(opts, "out", joinpath(run_dir, "analysis"))
inputs = get(opts, "inputs", "/shared/home/greg/breezelab-runs/20261006/inputs/mip_sources/seastarr_22241697")
mkpath(out_dir)
prefix = joinpath(run_dir, "sea_starr_" * lowercase(string(member)))
base_prefix = "sea_starr_" * lowercase(string(member))
epoch = DateTime(2017, 8, 15, 21)
CairoMakie.activate!(type = "png", px_per_unit = 1.5)
commit = try strip(read(`git rev-parse --short HEAD`, String)) catch; "unknown" end
numbers = Dict{String, Any}()
note(k, v) = (numbers[k] = v; v)

#####
##### Driver, trajectories (local solar time, SEVIRI composite)
#####

driver = read_dephy_driver(joinpath(inputs, "SEA_STARR_$(member)_SCM_driver.nc"))
tr = NCDataset(joinpath(inputs, "SEA_STARR_Raw_Trajectories.nc"))
trajectory_hours = Float64.(tr["t"][:])
get_tr(n) = (a = Float64.(coalesce.(tr[n][:, :], NaN)); a[a .<= -900] .= NaN; a)
lon_tr = get_tr("lon"); lat_tr = get_tr("lat")
mean_lon = [mean(filter(isfinite, lon_tr[n, :])) for n in axes(lon_tr, 1)]
mean_lat = [mean(filter(isfinite, lat_tr[n, :])) for n in axes(lat_tr, 1)]
sev = Dict(n => get_tr(n) for n in ("cf", "cf_warm", "lwp", "lwp_good", "nd", "nd_good", "height", "height_good", "sza", "cer_good"))
close(tr)
longitude_at(s) = interpolate_profile(trajectory_hours .* 3600, mean_lon, s)
local_solar_hour(s) = (h = (s / 3600 + 21 + longitude_at(s) / 15); mod(h, 24))
sev_mean(n) = [(v = filter(isfinite, sev[n][i, :]); isempty(v) ? NaN : mean(v)) for i in axes(sev[n], 1)]
sev_count(n) = [count(isfinite, sev[n][i, :]) for i in axes(sev[n], 1)]

#####
##### Model output
#####

ts_file = prefix * "_timeseries.jld2"
series(name) = stitched_series(run_dirs, base_prefix * "_timeseries.jld2", name)
value(f) = [f[n][1, 1, 1] for n in eachindex(f.times)]
cwp_s = series("cwp"); t1 = cwp_s.times
cwp = 1e3 .* value(cwp_s); rwp = 1e3 .* value(series("rwp")); cf = value(series("cloud_fraction"))
rain = 86400 .* value(series("rain_flux")); zi = value(series("zi")); zimax = value(series("zi_max"))
na_col = value(series("nᵃ_column")); nc_col = value(series("nᶜˡ_column")); nr_col = value(series("nʳ_column")); nt_col = value(series("n_total_column"))
th = t1 ./ 3600

st_file = prefix * "_statistics.jld2"
prof(name) = stitched_series(run_dirs, base_prefix * "_statistics.jld2", name)
θp = prof("θ"); t15 = θp.times; th15 = t15 ./ 3600
z = Array(znodes(θp.grid, Center())); zf = Array(znodes(θp.grid, Face())); Nz = length(z)
column(fts, n) = vec(Array(interior(fts[n], 1, 1, :)))
matrix(name) = (f = prof(name); hcat([column(f, n) for n in eachindex(f.times)]...))   # (z, t)
Θ = matrix("θ"); T = matrix("T"); QV = matrix("qᵛ"); QC = matrix("qᶜˡ"); QR = matrix("qʳ")
NA = matrix("nᵃ"); NC = matrix("nᶜˡ"); NR = matrix("nʳ"); CFP = matrix("cloud_fraction"); MASK = matrix("nudging_mask")
RFD = matrix("radiative_flux_divergence")
LWU = matrix("lw_up"); LWD = -matrix("lw_down"); SWU = matrix("sw_up"); SWD = -matrix("sw_down")   # Breeze stores downwelling as negative; make all positive magnitudes
U = matrix("u"); V = matrix("v"); W2 = matrix("w²")
# reference density profile from the driver surface pressure and θ/qᵗ (same construction as the case)
ρᵣ = let
    grid = RectilinearGrid(CPU(), Float64; size=(8, 8, Nz), x=(0, 400), y=(0, 400), z=zf, halo=(5, 5, 5), topology=(Periodic, Periodic, Bounded))
    constants = ThermodynamicConstants(Float64)
    rs = ReferenceState(grid, constants; base_pressure = driver.surface_pressure,
                        potential_temperature = ζ -> driver_initial_profile(driver, :thetal, ζ),
                        vapor_mass_fraction = ζ -> driver_initial_profile(driver, :qt, ζ))
    vec(Array(interior(rs.density, 1, 1, :)))
end

n15 = length(t15)
note("run_hours", th[end]); note("statistics_records", n15); note("timeseries_records", length(t1))

#####
##### 15-minute workbook-style statistics
#####

bin15(x) = [mean(x[(t1 .> t - 900) .& (t1 .≤ t)]) for t in t15]   # mean over the 15 min ending at each statistics time
cwp15 = bin15(cwp); rwp15 = bin15(rwp); lwp15 = cwp15 .+ rwp15; cf15 = bin15(cf); rain15 = bin15(rain); zi15 = bin15(zi); zimax15 = bin15(zimax)
cloud_threshold = 1e-5                                 # 0.01 g kg⁻¹ (workbook `cl` threshold)
base15 = [(k = findfirst(>(cloud_threshold), QC[:, n]); isnothing(k) ? NaN : z[k]) for n in 1:n15]
top15 = [(k = findlast(>(cloud_threshold), QC[:, n]); isnothing(k) ? NaN : z[k]) for n in 1:n15]
# cloud-top droplet number [cm⁻³]: nᶜˡ ρ at the highest level with ⟨qᶜˡ⟩ > threshold, in-cloud by dividing by the layer cloud fraction
nd_top15 = [(k = findlast(>(cloud_threshold), QC[:, n]); isnothing(k) ? NaN : 1e-6 * NC[k, n] * ρᵣ[k] / max(CFP[k, n], 0.05)) for n in 1:n15]
nd_max15 = [1e-6 * maximum(NC[:, n] .* ρᵣ) for n in 1:n15]
kbl(n) = max(2, findlast(≤(zi15[n]), z))
na_bl15 = [1e-6 * mean(NA[1:kbl(n), n]) for n in 1:n15]                     # mg⁻¹ (mean below ⟨zᵢ⟩)
nc_bl15 = [1e-6 * mean(NC[1:kbl(n), n]) for n in 1:n15]
na_ft15 = [1e-6 * NA[findfirst(≥(zimax15[n] + 500), z), n] for n in 1:n15]    # mg⁻¹ 500 m above the max inversion
mask_base15 = [(k = findfirst(>(0), MASK[:, n]); isnothing(k) ? NaN : z[k]) for n in 1:n15]
tke15 = [sum((W2[1:end-1, n] .+ W2[2:end, n]) ./ 2 .* diff(zf)) for n in 1:n15]   # ∫⟨w²⟩dz (m³ s⁻²); w² lives on faces

# radiation at the surface (k = 1) and the LES top (k = Nz + 1), W m⁻²
lwd_sfc = LWD[1, :]; lwu_sfc = LWU[1, :]; swd_sfc = SWD[1, :]; swu_sfc = SWU[1, :]
lwd_top = LWD[end, :]; lwu_top = LWU[end, :]; swd_top = SWD[end, :]; swu_top = SWU[end, :]
# cloud-top LW cooling: net LW flux (up − down) jump between 200 m below and 100 m above the mean cloud top
net_lw(k, n) = LWU[k, n] - LWD[k, n]
lw_cooling15 = [(isnan(top15[n]) ? NaN : (ka = min(Nz + 1, findfirst(≥(top15[n] + 100), zf)); kb = max(1, findlast(≤(top15[n] - 200), zf)); net_lw(ka, n) - net_lw(kb, n))) for n in 1:n15]
# strongest radiative cooling rate (K day⁻¹) in the column and its height
cool15 = [minimum(RFD[:, n]) for n in 1:n15]; cool_z15 = [z[argmin(RFD[:, n])] for n in 1:n15]
cp = 1005.0
cool_Kday15 = cool15 ./ (ρᵣ[[argmin(RFD[:, n]) for n in 1:n15]] .* cp) .* 86400

#####
##### Entrainment: w_e = dzᵢ/dt − w_s(zᵢ)
#####

# hourly-smoothed mean inversion height from the 1-min series, derivative over ±1.5 h
smooth(x, n) = [mean(x[max(1, i - n):min(length(x), i + n)]) for i in eachindex(x)]
zi_s = smooth(zi, 30)
dzidt = similar(zi_s)
Δ = 90
for i in eachindex(zi_s)
    a = max(1, i - Δ); b = min(length(zi_s), i + Δ)
    dzidt[i] = (zi_s[b] - zi_s[a]) / (t1[b] - t1[a])
end
w_sub(ζ) = interpolate_profile(driver.z, view(driver.forcing.w, :, 1), ζ)   # the driver's wa (identical at all times)
ws_zi = w_sub.(zi_s)
we = dzidt .- ws_zi
hourly = [findfirst(≥(h * 3600), t1) for h in 0:floor(Int, th[end])]
hourly = filter(!isnothing, hourly)
note("zi_initial_m", zi[1]); note("zi_final_m", zi_s[end]); note("zi_mean_rise_m_per_h", (zi_s[end] - zi_s[1]) / (th[end] - th[1]))
note("w_subsidence_at_zi_mean_mm_s", 1e3 * mean(ws_zi)); note("w_entrainment_mean_mm_s", 1e3 * mean(we)); note("w_entrainment_max_mm_s", 1e3 * maximum(we))
note("w_entrainment_day1_mm_s", 1e3 * mean(we[th .< 24])); note("w_entrainment_day2_mm_s", 1e3 * mean(we[(th .≥ 24) .& (th .< 48)])); note("w_entrainment_day3_mm_s", 1e3 * mean(we[th .≥ 48]))

#####
##### Subsidence as applied: rebuild the forcing on a small CPU case and compare with the driver
#####

subsidence_check = let
    Oceananigans.defaults.FloatType = Float64
    case = sea_starr(; member, arch = CPU(), Nx = 8, Ny = 8, driver_path = driver.path,
                       write_output = false, radiation = false, stop_time = 1.0, progress_interval = 1e9, checkpoint_interval = nothing)
    model = case.model
    grid = model.grid
    zc = Array(znodes(grid, Center()))
    # the materialized vertical-advection forcing of θ and its velocity time series
    unwrap(f::Breeze.Forcings.SpecificForcing) = unwrap(f.forcing)
    unwrap(f::Oceananigans.Forcings.MultipleForcings) = reduce(vcat, [collect(unwrap(g)) for g in f.forcings])
    unwrap(f) = [f]
    parts = unwrap(model.forcing.ρθ)
    vadv = parts[findfirst(p -> p isa LargeScaleVerticalAdvection, parts)]
    nvadv = count(p -> p isa LargeScaleVerticalAdvection, parts)
    w_applied = [vadv.velocity[1, 1, k, Time(3600.0)] for k in 1:length(zc)]
    w_driver = interpolate_profile(driver.z, view(driver.forcing.w, :, 2), zc)
    maxdiff = maximum(abs.(w_applied .- w_driver))
    # sign: a subsidence (w < 0) tendency on a θ profile increasing with height must be positive (warming)
    Oceananigans.TimeSteppers.update_state!(model; compute_tendencies = false)
    fields = Oceananigans.fields(model)
    k = findfirst(≥(2000), zc)
    Fθ = vadv(1, 1, k, grid, model.clock, fields)
    ∂zθ = (fields.θ[1, 1, k + 1] - fields.θ[1, 1, k]) / (zc[k + 1] - zc[k])
    (; nvadv, maxdiff, w_applied_2000 = w_applied[k], Fθ, ∂zθ, sign_ok = Fθ > 0 && w_applied[k] < 0 && ∂zθ > 0,
       nudged = [n for n in keys(model.forcing) if any(p -> p isa InversionFollowingNudging, unwrap(model.forcing[n]))])
end
note("subsidence_forcings_per_field", subsidence_check.nvadv); note("subsidence_max_abs_diff_vs_driver_m_s", subsidence_check.maxdiff)
note("subsidence_w_at_2000m_m_s", subsidence_check.w_applied_2000); note("subsidence_theta_tendency_at_2000m_K_s", subsidence_check.Fθ)
note("subsidence_sign_ok", subsidence_check.sign_ok); note("nudged_fields", string(subsidence_check.nudged))

#####
##### Model free troposphere vs the driver nudging targets
#####

ft_check = map((12, 24, 36, 48, 60)) do h
    n = findfirst(≥(h * 3600), t15); isnothing(n) && return nothing
    m = findfirst(≥(h * 3600), driver.times)
    θ_nud = interpolate_profile(driver.z, view(driver.forcing.thetal_nud, :, m), z)
    q_nud = interpolate_profile(driver.z, view(driver.forcing.qt_nud, :, m), z)
    na_nud = interpolate_profile(driver.z, view(driver.forcing.na_nud, :, m), z)
    sel = (z .≥ zimax15[n] + 300) .& (z .≤ 4000)
    (; h, dθ = mean(Θ[sel, n] .- θ_nud[sel]), dq = 1e3 * mean(QV[sel, n] .- q_nud[sel]), dna = 1e-6 * mean(NA[sel, n] .- na_nud[sel]),
       dθ_max = maximum(abs.(Θ[sel, n] .- θ_nud[sel])))
end
ft_check = filter(!isnothing, ft_check)
note("ft_theta_minus_target_K", [(c.h, round(c.dθ; digits = 3)) for c in ft_check])
note("ft_qv_minus_target_g_kg", [(c.h, round(c.dq; digits = 3)) for c in ft_check])
note("ft_na_minus_target_mg", [(c.h, round(c.dna; digits = 1)) for c in ft_check])

#####
##### Numbers
#####

idx(h) = findfirst(≥(h * 3600), t15)
for (lab, a, b) in (("day1", 0, 24), ("day2", 24, 48), ("day3", 48, 66.01))
    s = (th15 .≥ a) .& (th15 .< b)
    note("$(lab)_lwp_mean_g_m2", mean(lwp15[s])); note("$(lab)_cf_mean", mean(cf15[s])); note("$(lab)_rain_mean_mm_d", mean(rain15[s]))
    sv = (trajectory_hours .≥ a) .& (trajectory_hours .< b)
    note("$(lab)_seviri_cf_mean", mean(filter(isfinite, sev["cf"][sv, :]))); note("$(lab)_seviri_lwp_good_mean_g_m2", mean(filter(isfinite, sev["lwp_good"][sv, :])))
    note("$(lab)_seviri_lwp_all_mean_g_m2", mean(filter(isfinite, sev["lwp"][sv, :])))
end
note("lwp_max_g_m2", maximum(lwp15)); note("lwp_max_hour", th15[argmax(lwp15)]); note("lwp_min_g_m2", minimum(lwp15)); note("lwp_min_hour", th15[argmin(lwp15)])
note("cf_min", minimum(cf15)); note("cf_min_hour", th15[argmin(cf15)]); note("cf_final", cf15[end]); note("lwp_final_g_m2", lwp15[end])
note("rain_max_mm_d", maximum(rain15)); note("rwp_max_g_m2", maximum(rwp15))
note("cloud_base_initial_m", base15[1]); note("cloud_top_initial_m", top15[1]); note("cloud_top_final_m", top15[end]); note("cloud_base_final_m", base15[end])
note("nd_top_initial_cm3", nd_top15[1]); note("nd_top_day2_cm3", mean(filter(isfinite, nd_top15[(th15 .≥ 24) .& (th15 .< 48)]))); note("nd_top_final_cm3", nd_top15[end])
note("seviri_nd_good_day2_cm3", mean(filter(isfinite, sev["nd_good"][(trajectory_hours .≥ 24) .& (trajectory_hours .< 48), :])))
note("na_bl_initial_mg", na_bl15[1]); note("na_bl_final_mg", na_bl15[end]); note("na_ft_initial_mg", na_ft15[1]); note("na_ft_final_mg", na_ft15[end])
note("aerosol_column_initial_m2", nt_col[1]); note("aerosol_column_final_m2", nt_col[end]); note("surface_source_total_m2", 7e5 * t1[end])
note("lwd_top_mean_W_m2", mean(lwd_top)); note("lwu_top_mean_W_m2", mean(lwu_top)); note("swd_top_max_W_m2", maximum(swd_top))
note("lwd_sfc_mean_W_m2", mean(lwd_sfc)); note("lwu_sfc_mean_W_m2", mean(lwu_sfc)); note("swd_sfc_max_W_m2", maximum(swd_sfc))
note("cloud_top_lw_cooling_night_W_m2", mean(filter(isfinite, lw_cooling15[swd_top .< 1]))); note("cloud_top_lw_cooling_min_W_m2", minimum(filter(isfinite, lw_cooling15)))
note("max_radiative_cooling_K_day", minimum(cool_Kday15)); note("max_radiative_cooling_height_m", cool_z15[argmin(cool_Kday15)])
note("nudging_base_minus_zimax_mean_m", mean(filter(isfinite, mask_base15 .- zimax15)))
note("mean_lon_start_end", (mean_lon[1], mean_lon[end])); note("mean_lat_start_end", (mean_lat[1], mean_lat[end]))
note("sst_start_end_K", (driver.sst[1], driver.sst[end]))
note("seviri_cf_final_6h_mean", mean(filter(isfinite, sev["cf"][61:67, :]))); note("seviri_cloud_top_good_day3_km", mean(filter(isfinite, sev["height_good"][49:67, :])))

#####
##### Figures
#####

lst_ticks(ax) = begin
    hours = 0:6:66
    ax.xticks = (collect(hours), [@sprintf("%d\n%02d LST", h, round(Int, local_solar_hour(h * 3600))) for h in hours])
end
night = [(a, b) for (a, b) in zip(th15[1:end-1], th15[2:end]) if swd_top[findfirst(==(a), th15)] < 1]
shade_night!(ax) = for (a, b) in night; vspan!(ax, a, b, color = (:gray, 0.12)); end

fig = Figure(size = (1100, 1500), fontsize = 12)
Label(fig[0, 1:2], "SEA STARR $(member) ($(segment_label(run_dirs))): 15-min statistics vs SEVIRI along the composite trajectories (grey: night at the LES top)", font = :bold, tellwidth = false)
ax = Axis(fig[1, 1], ylabel = "Water path (g m⁻²)"); shade_night!(ax)
lines!(ax, th15, lwp15, label = "model LWP (cloud + rain)"); lines!(ax, th15, rwp15, label = "model RWP")
scatter!(ax, trajectory_hours, sev_mean("lwp"), color = :black, markersize = 6, label = "SEVIRI LWP (all)")
scatter!(ax, trajectory_hours, sev_mean("lwp_good"), color = :red, marker = :diamond, markersize = 7, label = "SEVIRI LWP (screened)")
axislegend(ax, position = :rt, labelsize = 10); lst_ticks(ax)
ax = Axis(fig[1, 2], ylabel = "Cloud fraction"); shade_night!(ax)
lines!(ax, th15, cf15, label = "model (LWP > 5 g m⁻²)"); scatter!(ax, trajectory_hours, sev_mean("cf"), color = :black, markersize = 6, label = "SEVIRI cf")
axislegend(ax, position = :lb, labelsize = 10); lst_ticks(ax)
ax = Axis(fig[2, 1], ylabel = "Surface rain (mm day⁻¹)"); shade_night!(ax); lines!(ax, th15, rain15); lst_ticks(ax)
ax = Axis(fig[2, 2], ylabel = "Height (m)"); shade_night!(ax)
lines!(ax, th15, zi15, label = "⟨zᵢ⟩"); lines!(ax, th15, zimax15, label = "max zᵢ"); lines!(ax, th15, top15, label = "cloud top"); lines!(ax, th15, base15, label = "cloud base")
lines!(ax, th15, mask_base15, label = "nudging base", linestyle = :dash, color = :gray)
scatter!(ax, trajectory_hours, 1e3 .* sev_mean("height_good"), color = :black, markersize = 6, label = "SEVIRI cloud top (screened)")
axislegend(ax, position = :lt, labelsize = 10); lst_ticks(ax)
ax = Axis(fig[3, 1], ylabel = "Droplet number (cm⁻³)"); shade_night!(ax)
lines!(ax, th15, nd_top15, label = "model, cloud top (in-cloud)"); lines!(ax, th15, nd_max15, label = "model, column max", linestyle = :dot)
scatter!(ax, trajectory_hours, sev_mean("nd_good"), color = :black, markersize = 6, label = "SEVIRI Nd (screened)")
axislegend(ax, position = :lt, labelsize = 10); lst_ticks(ax); ylims!(ax, 0, 600)
ax = Axis(fig[3, 2], ylabel = "Aerosol number (mg⁻¹)"); shade_night!(ax)
lines!(ax, th15, na_bl15, label = "⟨nᵃ⟩ below zᵢ"); lines!(ax, th15, nc_bl15, label = "⟨nᶜˡ⟩ below zᵢ"); lines!(ax, th15, na_ft15, label = "nᵃ at max zᵢ + 500 m")
axislegend(ax, position = :rt, labelsize = 10); lst_ticks(ax)
ax = Axis(fig[4, 1], ylabel = "Column number (10¹² m⁻²)"); shade_night!(ax)
lines!(ax, th, 1e-12 .* na_col, label = "nᵃ"); lines!(ax, th, 1e-12 .* nc_col, label = "nᶜˡ"); lines!(ax, th, 1e-12 .* nr_col .* 1e3, label = "nʳ × 10³"); lines!(ax, th, 1e-12 .* nt_col, label = "total", linestyle = :dash)
axislegend(ax, position = :lt, labelsize = 10); lst_ticks(ax)
ax = Axis(fig[4, 2], ylabel = "Radiation (W m⁻²)"); shade_night!(ax)
lines!(ax, th15, swd_top, label = "SW↓ LES top"); lines!(ax, th15, swd_sfc, label = "SW↓ surface"); lines!(ax, th15, lwd_sfc, label = "LW↓ surface"); lines!(ax, th15, lwu_top, label = "LW↑ LES top"); lines!(ax, th15, lwd_top, label = "LW↓ LES top (= 0: departure)")
axislegend(ax, position = :rt, labelsize = 9); lst_ticks(ax)
ax = Axis(fig[5, 1], ylabel = "Entrainment (mm s⁻¹)", xlabel = "hours since 2017-08-15 21 UTC / local solar time"); shade_night!(ax)
lines!(ax, th, 1e3 .* dzidt, label = "dzᵢ/dt"); lines!(ax, th, 1e3 .* ws_zi, label = "w_s(zᵢ) (driver)"); lines!(ax, th, 1e3 .* we, label = "w_e = dzᵢ/dt − w_s")
axislegend(ax, position = :rt, labelsize = 10); lst_ticks(ax)
ax = Axis(fig[5, 2], ylabel = "Cloud-top LW cooling (W m⁻²)", xlabel = "hours / local solar time"); shade_night!(ax)
lines!(ax, th15, lw_cooling15, label = "ΔF_LW across cloud top (−200 m, +100 m)"); lines!(ax, th15, cool_Kday15 ./ 10, label = "min heating rate / 10 (K day⁻¹)", linestyle = :dot)
axislegend(ax, position = :rb, labelsize = 10); lst_ticks(ax)
save(joinpath(out_dir, "ctrl_timeseries_vs_seviri.png"), fig)

# time–height
fig2 = Figure(size = (1200, 1400), fontsize = 12)
Label(fig2[0, 1:4], "SEA STARR $(member): horizontal means (15-min) below 3.5 km, with ⟨zᵢ⟩ (white) and the nudging base (dashed)", font = :bold, tellwidth = false)
kz = findlast(≤(3500), z)
panels = (("θₗ (K)", Θ, :thermal, (288, 318)), ("qᵗ = qᵛ + qᶜˡ (g kg⁻¹)", 1e3 .* (QV .+ QC), :viridis, (0, 12)),
          ("qᶜˡ (g kg⁻¹)", 1e3 .* QC, :Blues, (0, 0.5)), ("nᵃ (mg⁻¹)", 1e-6 .* NA, :magma, (0, 1200)),
          ("nᶜˡ (mg⁻¹)", 1e-6 .* NC, :viridis, (0, 300)), ("radiative heating (K day⁻¹)", RFD ./ (ρᵣ .* cp) .* 86400, :RdBu, (-40, 40)),
          ("nudging mask", MASK, :grays, (0, 1)), ("layer cloud fraction", CFP, :Blues, (0, 1)))
for (i, (ttl, A, cmap, rng)) in enumerate(panels)
    r, c = fldmod1(i, 2)
    ax = Axis(fig2[r, 2c - 1], title = ttl, ylabel = c == 1 ? "z (m)" : "", xlabel = r == 4 ? "hours since 2017-08-15 21 UTC" : "")
    hm = heatmap!(ax, th15, z[1:kz], permutedims(A[1:kz, :]), colormap = cmap, colorrange = rng)
    lines!(ax, th15, zi15, color = :white, linewidth = 1); lines!(ax, th15, mask_base15, color = :white, linestyle = :dash, linewidth = 1)
    Colorbar(fig2[r, 2c], hm, width = 8)
end
save(joinpath(out_dir, "ctrl_time_height.png"), fig2)

# profiles vs nudging targets
fig3 = Figure(size = (1200, 450), fontsize = 12)
Label(fig3[0, 1:3], "SEA STARR $(member): mean profiles (solid) and driver nudging targets (dashed) — nudging acts only above max zᵢ + 100 m", font = :bold, tellwidth = false)
axθ = Axis(fig3[1, 1], xlabel = "θₗ (K)", ylabel = "z (m)"); axq = Axis(fig3[1, 2], xlabel = "qᵛ (g kg⁻¹)"); axn = Axis(fig3[1, 3], xlabel = "nᵃ (mg⁻¹)")
cols = Makie.wong_colors()
for (i, h) in enumerate((0, 24, 48, 66))
    n = idx(min(h, th15[end])); isnothing(n) && continue
    m = findfirst(≥(h * 3600), driver.times); isnothing(m) && (m = length(driver.times))
    lines!(axθ, Θ[1:kz, n], z[1:kz], color = cols[i], label = "t = $h h"); lines!(axθ, interpolate_profile(driver.z, view(driver.forcing.thetal_nud, :, m), z[1:kz]), z[1:kz], color = cols[i], linestyle = :dash)
    lines!(axq, 1e3 .* QV[1:kz, n], z[1:kz], color = cols[i]); lines!(axq, 1e3 .* interpolate_profile(driver.z, view(driver.forcing.qt_nud, :, m), z[1:kz]), z[1:kz], color = cols[i], linestyle = :dash)
    lines!(axn, 1e-6 .* NA[1:kz, n], z[1:kz], color = cols[i]); lines!(axn, 1e-6 .* interpolate_profile(driver.z, view(driver.forcing.na_nud, :, m), z[1:kz]), z[1:kz], color = cols[i], linestyle = :dash)
end
axislegend(axθ, position = :rb); xlims!(axθ, 286, 320)
save(joinpath(out_dir, "ctrl_profiles_vs_targets.png"), fig3)

#####
##### Summary files
#####

open(joinpath(out_dir, "summary.toml"), "w") do io
    TOML.print(io, Dict("run" => run_spec, "member" => string(member), "commit" => commit, "generated" => string(now(UTC)),
                        "numbers" => Dict(k => (v isa Tuple || v isa Vector ? string(v) : v) for (k, v) in numbers)))
end
open(joinpath(out_dir, "summary.md"), "w") do io
    println(io, "# SEA STARR $(member) analysis summary (", segment_label(run_dirs), ", commit ", commit, ", ", Dates.format(now(UTC), "yyyy-mm-dd HH:MM"), " UTC)\n")
    println(io, "| quantity | value |\n| --- | --- |")
    for k in sort(collect(keys(numbers)))
        v = numbers[k]
        println(io, "| ", k, " | ", v isa AbstractFloat ? @sprintf("%.4g", v) : string(v), " |")
    end
    println(io, "\n## Hourly table (15-min statistics at the hour)\n")
    println(io, "| t (h) | UTC | LST | LWP | CF | rain | ⟨zᵢ⟩ | cloud top | Nd top | nᵃ BL | w_e (mm/s) | LW↓ sfc | SW↓ top | cloud-top ΔF_LW |")
    println(io, "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |")
    for h in 0:3:floor(Int, th15[end])
        n = idx(h); i1 = findfirst(≥(h * 3600), t1)
        @printf(io, "| %d | %s | %05.2f | %.0f | %.2f | %.3f | %.0f | %.0f | %.0f | %.0f | %.2f | %.0f | %.0f | %.1f |\n", h, Dates.format(epoch + Hour(h), "dd HH"), local_solar_hour(h * 3600),
                lwp15[n], cf15[n], rain15[n], zi15[n], top15[n], nd_top15[n], na_bl15[n], 1e3 * we[i1], lwd_sfc[n], swd_top[n], lw_cooling15[n])
    end
    println(io, "\n## SEVIRI composite along the 40 trajectories (mean over trajectories)\n")
    println(io, "| t (h) | cf | LWP all (g m⁻²) | LWP screened | Nd screened (cm⁻³) | cloud top screened (km) | n(cf) | n(LWP screened) |\n| --- | --- | --- | --- | --- | --- | --- | --- |")
    for i in 1:3:length(trajectory_hours)
        @printf(io, "| %d | %.2f | %.1f | %.1f | %.1f | %.2f | %d | %d |\n", trajectory_hours[i], sev_mean("cf")[i], sev_mean("lwp")[i], sev_mean("lwp_good")[i], sev_mean("nd_good")[i], sev_mean("height_good")[i], sev_count("cf")[i], sev_count("lwp_good")[i])
    end
end
println("analysis written to ", out_dir)
for k in sort(collect(keys(numbers)))
    println(k, " = ", numbers[k])
end
