# Aerosol/droplet audit of completed ENA P3-aer2 runs: profiles of nᵃ, nᶜˡ, qᶜˡ, cloud
# fraction (hourly means), the time series, and cloud-slab slices. Writes figures and a TOML
# summary to OUTPUT_DIR. Usage:
#   julia --project analysis/ena_aerosol_audit.jl OUTPUT_DIR RUN_DIR [RUN_DIR ...]
using CairoMakie, JLD2, Oceananigans, Oceananigans.Units, Statistics, TOML, Printf
using Oceananigans.Grids: Center, Face, znodes
using Oceananigans.Fields: interior

output_dir = ARGS[1]; runs = ARGS[2:end]
mkpath(output_dir)
summary = Dict{String, Any}()

function profiles(dir)
    file = only(filter(f -> endswith(f, "_profiles.jld2"), readdir(dir; join=true)))
    names = ["nᵃ", "nᶜˡ", "qᶜˡ", "qʳ", "qᵛ", "T", "cloud_fraction", "w²"]
    d = Dict(n => FieldTimeSeries(file, n) for n in names)
    return d
end
column(f) = vec(Array(interior(f)))
fig_all = Figure(size=(1500, 1100), fontsize=13)
colors = Makie.wong_colors()
for (r, dir) in enumerate(runs)
    label = basename(rstrip(dir, '/'))
    prov = TOML.parsefile(joinpath(dir, "provenance.toml"))
    cfg = prov["config"]
    ρ_surface = cfg["reference_density"]
    n_total = cfg["n₁"] + cfg["n₂"]                       # kg⁻¹, from the provenance conversion
    N_total = cfg["N₁"] + cfg["N₂"]                       # cm⁻³ at the surface
    p = profiles(dir)
    t = p["nᵃ"].times
    z = collect(znodes(p["nᵃ"].grid, Center()))
    nz = length(z)
    # reference density profile from the saved grid is not stored; use the surface value for cm⁻³ scaling (noted)
    s = Dict{String, Any}("label" => label, "n_initial_kg" => n_total, "N_initial_cm3_surface" => N_total,
                          "reference_density_surface" => ρ_surface, "hours" => collect(t ./ 3600))
    axn = Axis(fig_all[1, r], xlabel="nᵃ (10⁶ kg⁻¹)", ylabel="z (m)", title="$label\naerosol reservoir nᵃ (hourly means)", limits=(nothing, (0, 3000)))
    axc = Axis(fig_all[2, r], xlabel="nᶜˡ (10⁶ kg⁻¹)", ylabel="z (m)", title="droplet number nᶜˡ", limits=(nothing, (0, 3000)))
    axq = Axis(fig_all[3, r], xlabel="qᶜˡ (g kg⁻¹)  /  cloud fraction", ylabel="z (m)", title="cloud liquid (solid), cloudy fraction (dashed)", limits=(nothing, (0, 3000)))
    vlines!(axn, [n_total / 1e6]; color=:black, linestyle=:dot, label="initial total")
    for n in eachindex(t)
        c = get(colors, mod1(n, length(colors)), :gray)
        na = column(p["nᵃ"][n]); nc = column(p["nᶜˡ"][n]); qc = column(p["qᶜˡ"][n]); cf = column(p["cloud_fraction"][n])
        lines!(axn, na ./ 1e6, z; color=c, label=@sprintf("%.0f–%.0f h", (n == 1 ? 0 : t[n-1]) / 3600, t[n] / 3600))
        lines!(axc, nc ./ 1e6, z; color=c)
        lines!(axq, 1e3 .* qc, z; color=c)
        lines!(axq, cf, z; color=c, linestyle=:dash)
        cloudy = cf .> 0.05
        incloud_nc = any(cloudy) ? mean(nc[cloudy] ./ max.(cf[cloudy], 1e-6)) : NaN   # mean over cloudy cells ≈ ⟨nᶜˡ⟩/CF
        s["hour_$(n)"] = Dict("t_end_s" => t[n],
                               "na_min" => minimum(na), "na_max" => maximum(na), "na_surface" => na[1], "na_top" => na[end],
                               "na_mean_below_3km" => mean(na[z .< 3000]),
                               "na_incloud_mean" => any(cloudy) ? mean(na[cloudy]) : NaN,
                               "nc_max" => maximum(nc), "nc_incloud_per_cloudy_cell" => incloud_nc,
                               "nc_plus_na_incloud" => any(cloudy) ? mean((na .+ nc)[cloudy]) : NaN,
                               "qc_max_gkg" => 1e3 * maximum(qc), "cf_max" => maximum(cf),
                               "cloud_base_m" => any(cloudy) ? z[findfirst(cloudy)] : NaN,
                               "cloud_top_m" => any(cloudy) ? z[findlast(cloudy)] : NaN)
    end
    axislegend(axn, position=:rb, labelsize=10)
    # time series
    tsf = only(filter(f -> endswith(f, "_timeseries.jld2"), readdir(dir; join=true)))
    lwp = FieldTimeSeries(tsf, "lwp"); cf = FieldTimeSeries(tsf, "cloud_fraction"); rain = FieldTimeSeries(tsf, "rain_flux")
    tt = lwp.times
    s["timeseries"] = Dict("lwp_g_m2_mean" => mean(1e3 * Array(interior(lwp[n]))[1] for n in eachindex(tt)),
                           "lwp_g_m2_final" => 1e3 * Array(interior(lwp[end]))[1],
                           "cloud_fraction_final" => Array(interior(cf[end]))[1],
                           "rain_mm_day_mean" => mean(86400 * Array(interior(rain[n]))[1] for n in eachindex(tt)))
    # slices: cloud liquid at 900 m and the xz section at the end
    slf = only(filter(f -> endswith(f, "_slices.jld2"), readdir(dir; join=true)))
    qxz = FieldTimeSeries(slf, "qᶜˡ_xz"); qxy = FieldTimeSeries(slf, "qᶜˡ_xy")
    fs = Figure(size=(1200, 500))
    ax1 = Axis(fs[1, 1], title="$label qᶜˡ (g kg⁻¹) at 900 m, t = $(round(qxy.times[end]/3600; digits=1)) h", xlabel="x (km)", ylabel="y (km)")
    ax2 = Axis(fs[1, 2], title="qᶜˡ x–z section", xlabel="x (km)", ylabel="z (m)", limits=(nothing, (0, 2500)))
    g = qxy.grid
    x = collect(Oceananigans.Grids.xnodes(g, Center())) ./ 1e3; y = collect(Oceananigans.Grids.ynodes(g, Center())) ./ 1e3
    hm = heatmap!(ax1, x, y, 1e3 .* Array(interior(qxy[end]))[:, :, 1]; colormap=:Blues); Colorbar(fs[1, 0], hm)
    zz = collect(znodes(qxz.grid, Center()))
    hm2 = heatmap!(ax2, x, zz, 1e3 .* Array(interior(qxz[end]))[:, 1, :]; colormap=:Blues); Colorbar(fs[1, 3], hm2)
    save(joinpath(output_dir, "slices_$label.png"), fs; px_per_unit=2)
    summary[label] = s
    println(label, ": initial nᵃ = ", n_total, " kg⁻¹ (", N_total, " cm⁻³ at ρ₀ = ", ρ_surface, ")")
    for n in eachindex(t)
        h = s["hour_$(n)"]
        @printf("  hour %d: nᵃ min/max %.3e/%.3e, in-cloud nᵃ %.3e, in-cloud nᶜˡ/CF %.3e, nᶜˡ+nᵃ in cloud %.3e kg⁻¹, qc max %.3f g/kg, CF max %.2f, base/top %.0f/%.0f m\n",
                n, h["na_min"], h["na_max"], h["na_incloud_mean"], h["nc_incloud_per_cloudy_cell"], h["nc_plus_na_incloud"], h["qc_max_gkg"], h["cf_max"], h["cloud_base_m"], h["cloud_top_m"])
    end
end
save(joinpath(output_dir, "aerosol_profiles.png"), fig_all; px_per_unit=2)
open(io -> TOML.print(io, summary), joinpath(output_dir, "aerosol_audit.toml"), "w")
println("wrote ", output_dir)
