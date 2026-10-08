# Acceptance checks of a short LASSO integration on a CPU node: aerosol/droplet profiles, radiation and
# state extrema, and the SAM samstat comparison. Usage: sbatch ... --wrap="bash execution/lasso_short_checks.sh [RUN_DIR]"
set -e
cd /shared/home/greg/breezelab-work/ena-lasso
source execution/wpcluster_env.sh; export JULIA_NUM_THREADS=4
R=${1:-runs/lasso_aer2_1h_job243}
"$JULIA" --startup-file=no --project=. analysis/ena_aerosol_audit.jl $R/checks $R
"$JULIA" --startup-file=no --project=. analysis/compare_sam_reference.jl $R data/lasso/archives/enalasso_samstat_20170718era5d25x100_sbmwrm-aer2-flxsstC1.m1.20170718.000000.nc /shared/home/greg/breezelab-runs/20261006/inputs/observations/20170718 $R/checks/sam
"$JULIA" --startup-file=no --project=. -e ' "$R"
using Oceananigans, JLD2, Statistics, TOML
using Oceananigans.Grids: Center, znodes
using Oceananigans.Fields: interior
R = ARGS[1]
f = joinpath(R, "lasso_ena_profiles.jld2")
for name in ("radiative_flux_divergence", "T", "u", "v", "qᵛ", "nᵃ", "nᶜˡ", "qʳ")
    fts = FieldTimeSeries(f, name); z = collect(znodes(fts.grid, Center()))
    for n in (1, length(fts.times))
        c = vec(Array(interior(fts[n])))
        println(rpad(name, 26), " t=", fts.times[n], " min ", minimum(c), " max ", maximum(c), " z(max) ", z[argmax(c)], " surface ", c[1], " top ", c[end])
    end
end
p = TOML.parsefile(joinpath(R, "provenance.toml"))
println("config: ", (; p["config"]["surface_flux_law"], p["config"]["coriolis_parameter"], p["config"]["aerosol_supersaturation_cap"], p["config"]["radiation_interval"], p["config"]["liquid_effective_radius"], p["config"]["surface_emissivity"], p["config"]["n₁"], p["config"]["n₂"]))
println("bundle readme hash: ", get(get(p["bundle"], "readme", Dict()), "model_source_git_hash", "?"))' "$R"
