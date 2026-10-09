# Column water and moist-enthalpy budgets of a TRACER–DP-SCREAM run (see src/tracer_dp_scream/budget.jl),
# totals and per day, written to <out>/budget_enthalpy.toml and printed.
#   julia --project analysis/tracer_dp_scream_budget.jl <run dir> <IOP file> [<out dir> = <run dir>] [applied_vapor_cross_term = true]
using BreezeLab, Oceananigans, TOML, Printf, Dates
run_dir, iop_path = ARGS[1], ARGS[2]
out_dir = length(ARGS) ≥ 3 ? ARGS[3] : run_dir
cross = length(ARGS) ≥ 4 ? parse(Bool, ARGS[4]) : true
mkpath(out_dir)
prov = TOML.parsefile(joinpath(run_dir, "provenance.toml"))
config = prov["config"]
FT = config["FT"] == "Float32" ? Float32 : Float64
Oceananigans.defaults.FloatType = FT
start = DateTime(get(config, "start", "2022-08-01T00:00:00")); stop = DateTime(get(config, "stop", "2022-08-15T00:00:00"))
# Same vertical grid, window and forcing as the run; a small horizontal grid suffices for ρᵣ and the forcing profiles.
case = tracer_dp_scream(; iop_path, start, stop, arch = CPU(), Nx = 8, Ny = 8, Δx = Float64(config["Lx"]) / 8,
                        radiation = false, write_output = false)
@assert size(case.grid, 3) == config["Nz"]
budget = tracer_dp_scream_budget(case, joinpath(run_dir, "tracer_timeseries.jld2"), joinpath(run_dir, "tracer_profiles.jld2");
                                 applied_vapor_cross_term = cross)
toml(x::NamedTuple) = Dict(string(k) => toml(v) for (k, v) in pairs(x))
toml(x::AbstractVector) = [toml(v) for v in x]
toml(x::Real) = Float64(x)
toml(x) = x
open(joinpath(out_dir, "budget_enthalpy.toml"), "w") do io
    TOML.print(io, Dict("run" => abspath(run_dir), "total" => toml(budget.total), "daily" => toml(budget.periods),
                        "applied_vapor_cross_term" => cross, "constants" => toml(budget.constants),
                        "method" => "moist enthalpy h = cpm(T−Tr) + ℒˡ qv + (ℒˡ−ℒⁱ) qi from the 60-s series; see src/tracer_dp_scream/budget.jl"))
end
@printf("%-4s %10s %10s %10s %10s %10s %10s | %9s %9s %9s | %10s %10s %8s\n", "day", "ΔH MJ", "sfc", "rad", "LS", "prec", "resid", "rel", "W resid", "W rel", "S resid", "ℒP", "Xsfc+Xls")
for (d, b) in enumerate(vcat(budget.periods, [budget.total]))
    e = b.energy; w = b.water
    @printf("%-4s %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f | %9.4f %9.3f %9.4f | %10.2f %10.2f %8.2f\n",
            d ≤ length(budget.periods) ? string(d - 1) : "all", e.ΔH / 1e6, e.surface / 1e6, e.radiation / 1e6, e.large_scale / 1e6,
            e.precipitation / 1e6, e.residual / 1e6, e.relative, w.residual, w.relative, b.legacy_static_energy_residual / 1e6,
            b.latent_heat_of_precipitation / 1e6, (b.vapor_cross_term_surface + b.vapor_cross_term_large_scale) / 1e6)
end
