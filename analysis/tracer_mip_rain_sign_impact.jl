# Size of the rain-to-land sign error in an outer run made before 2026-10-10 (rain removed from the slab bucket
# instead of added). Uses the run's own correctly signed surface rain flux (BreezeLab `surface_rain_flux`, saved as
# `rain` in outer_surface.jld2) over land cells (sea mask 0); the bucket error is ≈ 2 × the accumulated rain where the
# bucket stayed between empty and full.
#   julia --project=cases/tracer_mip analysis/tracer_mip_rain_sign_impact.jl <run dir> [<out dir>]
using Oceananigans, Statistics, Printf
using Oceananigans.Fields: interior
using TOML: TOML

run_dir = abspath(ARGS[1]); out_dir = abspath(get(ARGS, 2, joinpath(run_dir, "analysis"))); mkpath(out_dir)
file = joinpath(run_dir, "outer_surface.jld2")
rain = FieldTimeSeries(file, "rain"; backend = OnDisk())
sat = FieldTimeSeries(file, "saturation"; backend = OnDisk())
sea = dropdims(Array(interior(FieldTimeSeries(file, "sea"; backend = OnDisk())[1])); dims = 3) .> 0.5
land = .!sea
times = Float64.(collect(rain.times)); nt = length(times)
Δts = vcat(diff(times), times[end] > times[1] ? times[end] - times[end-1] : 0.0)
acc = zeros(size(sea))                 # kg m⁻² = mm
land_mean_rate = Float64[]; land_max_rate = Float64[]
for n in 1:nt
    r = dropdims(Array(interior(rain[n])); dims = 3)
    n < nt && (acc .+= max.(r, 0) .* Δts[n])
    push!(land_mean_rate, 3600 * mean(r[land])); push!(land_max_rate, 3600 * maximum(r[land]))
end
a = acc[land]
s₀ = dropdims(Array(interior(sat[1])); dims = 3)[land]; s₁ = dropdims(Array(interior(sat[nt])); dims = 3)[land]
rainy = a .> 1
M_max = 150.0                          # BucketHydrology default capacity [kg m⁻²]
report = Dict("run_dir" => run_dir, "hours" => (times[end] - times[1]) / 3600, "n_land_cells" => count(land),
    "land_mean_accumulated_rain_mm" => mean(a), "land_max_accumulated_rain_mm" => maximum(a),
    "land_fraction_rain_gt_1mm" => mean(a .> 1), "land_fraction_rain_gt_10mm" => mean(a .> 10),
    "land_p99_accumulated_rain_mm" => quantile(a, 0.99),
    "implied_bucket_error_land_mean_kg_m2" => 2 * mean(a), "implied_bucket_error_max_kg_m2" => 2 * maximum(a),
    "implied_saturation_error_land_mean" => 2 * mean(a) / M_max, "implied_saturation_error_max" => 2 * maximum(a) / M_max,
    "saturation_land_mean_start" => mean(s₀), "saturation_land_mean_end" => mean(s₁),
    "saturation_change_rainy_cells_mean" => any(rainy) ? mean(s₁[rainy] .- s₀[rainy]) : NaN,
    "saturation_change_dry_cells_mean" => mean(s₁[.!rainy] .- s₀[.!rainy]),
    "n_rainy_land_cells_gt_1mm" => count(rainy),
    "land_mean_rain_rate_mm_h_by_record" => land_mean_rate, "land_max_rain_rate_mm_h_by_record" => land_max_rate,
    "record_hours" => (times .- times[1]) ./ 3600)
open(joinpath(out_dir, "rain_sign_impact.toml"), "w") do io; TOML.print(io, report); end
@info @sprintf("land: mean accumulated rain %.3f mm, max %.1f mm, %.2f%% of cells > 1 mm, %.3f%% > 10 mm; implied bucket error mean %.3f, max %.1f kg/m² (saturation %.4f / %.3f); saturation mean %.3f → %.3f (rainy cells Δ %.3f, dry cells Δ %.3f)",
               mean(a), maximum(a), 100mean(a .> 1), 100mean(a .> 10), 2mean(a), 2maximum(a), 2mean(a) / M_max, 2maximum(a) / M_max,
               mean(s₀), mean(s₁), report["saturation_change_rainy_cells_mean"], report["saturation_change_dry_cells_mean"])
