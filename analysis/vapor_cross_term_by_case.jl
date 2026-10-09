# Size of the "(cᵖᵛ − cᵖᵈ) T × vapor-rate" heating that the energy forcing and the prescribed surface
# energy flux applied in the liquid-ice potential-temperature formulation before the
# formulation-aware fix (where it is spurious: cᵖᵐ dT = F there, and a vapor source at fixed θ leaves T
# unchanged). Computed from each constructor's own inputs on a tiny CPU grid (the vertical grid and
# forcing are the cases' own), using the initial temperature/moisture columns, as window means:
#
#   large scale:  (cᵖᵛ − cᵖᵈ) ∫ ρᵣ T r dz,  r = vapor rate of the qls source (mixing-ratio basis: (1 − qᵛ)² R)
#   surface:      (cᵖᵛ − cᵖᵈ) SST LE / ℒˡᵣ  (prescribed-flux cases with temperature_neutral_evaporation)
#
#   julia --project analysis/vapor_cross_term_by_case.jl <covert data dir> <LASSO bundle dir> <TRACER IOP file> [out.toml]
using BreezeLab, Oceananigans, Breeze, TOML, Printf, Dates
using BreezeLab: day_to_seconds, lasso_documented_dimensions, epoch_from_day_of_year, parse_lasso_member, DEFAULT_LASSO_MEMBER
covert_dir, lasso_dir, iop_path = ARGS[1:3]
out = length(ARGS) ≥ 4 ? ARGS[4] : "vapor_cross_term_by_case.toml"
Oceananigans.defaults.FloatType = Float64
constants = ThermodynamicConstants(Float64)
cᵖᵈ = constants.dry_air.heat_capacity; cᵖᵛ = constants.vapor.heat_capacity; ℒ = constants.liquid.reference_latent_heat

function large_scale_cross(case, stop_time; basis)
    grid = case.grid
    ρ = Array(interior(case.model.dynamics.reference_state.density, 1, 1, :))
    Δz = diff(Array(znodes(grid, Face())))
    T = case.columns.T; qᵛ = case.columns.qᵛ
    factor = basis === :mixing_ratio ? (1 .- qᵛ) .^ 2 : ones(length(qᵛ))
    qls = case.forcing_profiles.qls
    ts = range(0, stop_time, length = 241)
    vals = map(ts) do t
        n = clamp(searchsortedlast(qls.times, t), 1, length(qls.times))
        r = Array(interior(qls[n], 1, 1, :)) .* factor
        (cᵖᵛ - cᵖᵈ) * sum(ρ .* T .* r .* Δz)
    end
    return sum(vals) / length(vals), maximum(abs, vals)
end

function surface_cross(sfc, day0, stop_time)
    times = day_to_seconds.(sfc.day, day0)
    ts = range(0, stop_time, length = 241)
    held(v, t) = v[clamp(searchsortedlast(times, t), 1, length(v))]
    interp(v, t) = BreezeLab.interpolate_profile(times, v, t)
    vals = [(cᵖᵛ - cᵖᵈ) * interp(sfc.sst, t) * interp(sfc.latent_heat_flux, t) / ℒ for t in ts]
    return sum(vals) / length(vals), maximum(abs, vals)
end

results = Dict{String, Any}()

covert = ena_covert(; arch = CPU(), data_dir = covert_dir, Nx = 8, Ny = 8, Lx = 280, Ly = 280, write_output = false)
st = covert.config.stop_time
ls, ls_max = large_scale_cross(covert, st; basis = :mixing_ratio)
sf, sf_max = surface_cross(covert.sfc, covert.config.day0, st)
results["ena_covert"] = Dict("window_hours" => st / 3600, "large_scale_W_m2" => ls, "large_scale_max_W_m2" => ls_max,
                             "surface_W_m2" => sf, "surface_max_W_m2" => sf_max, "total_W_m2" => ls + sf,
                             "paths" => "LargeScaleEnergyForcing (mixing-ratio basis) + prescribed H/LE with temperature_neutral_evaporation + top-2-level relaxation cross term")

member = DEFAULT_LASSO_MEMBER
inspection = inspect_lasso_bundle(lasso_dir; member, dimensions = lasso_documented_dimensions(member))
epoch = epoch_from_day_of_year(inspection.time.day0; year = Dates.year(parse_lasso_member(member).date))
lasso = ena_lasso(; member, epoch, bundle_dir = lasso_dir, arch = CPU(), Nx = 8, Ny = 8, Lx = 800, Ly = 800,
                  radiation = false, write_output = false)
st = lasso.config.stop_time
ls, ls_max = large_scale_cross(lasso, st; basis = :mixing_ratio)
results["ena_lasso"] = Dict("window_hours" => st / 3600, "large_scale_W_m2" => ls, "large_scale_max_W_m2" => ls_max,
                            "surface_W_m2" => 0.0, "total_W_m2" => ls,
                            "paths" => "LargeScaleEnergyForcing (mixing-ratio basis) + top-2-level relaxation cross term; Breeze bulk SST fluxes (unaffected)")

tracer = tracer_dp_scream(; iop_path, arch = CPU(), Nx = 8, Ny = 8, Δx = 6400.0, radiation = false, write_output = false)
st = tracer.config.stop_time
ls, ls_max = large_scale_cross(tracer, st; basis = :mass_fraction)
sf, sf_max = surface_cross(tracer.sfc, tracer.config.day0, st)
results["tracer_dp_scream"] = Dict("window_hours" => st / 3600, "large_scale_W_m2" => ls, "large_scale_max_W_m2" => ls_max,
                                   "surface_W_m2" => sf, "surface_max_W_m2" => sf_max, "total_W_m2" => ls + sf,
                                   "paths" => "LargeScaleEnergyForcing (mass-fraction basis) + prescribed IOP H/LE with temperature_neutral_evaporation",
                                   "note" => "initial temperature column used for T; the run-based value is in the enthalpy budget")
results["sea_starr"] = Dict("total_W_m2" => 0.0, "paths" => "θ nudging, full-field vertical advection, Breeze bulk fluxes: no cross term")
results["tracer_mip"] = Dict("total_W_m2" => 0.0, "paths" => "ext/tracer_mip_regional: no LargeScaleEnergyForcing or prescribed-flux cross term")
open(io -> TOML.print(io, results), out, "w")
for (k, v) in sort(collect(results); by = first)
    @printf("%-18s total %7.2f W/m²  (large scale %s, surface %s)\n", k, v["total_W_m2"],
            haskey(v, "large_scale_W_m2") ? @sprintf("%.2f", v["large_scale_W_m2"]) : "—",
            haskey(v, "surface_W_m2") ? @sprintf("%.2f", v["surface_W_m2"]) : "—")
end
