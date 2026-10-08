# Validate a staged LASSO-ENA bundle and build/step a tiny CPU case from it (no GPU needed):
#   julia --project analysis/check_lasso_bundle.jl [BUNDLE_DIR] [MEMBER]
using BreezeLab, Oceananigans, Dates, TOML
member = length(ARGS) ≥ 2 ? ARGS[2] : BreezeLab.DEFAULT_LASSO_MEMBER
dir = length(ARGS) ≥ 1 ? ARGS[1] : lasso_bundle_directory(member)
m = parse_lasso_member(member)
dims = lasso_documented_dimensions(m)
b = inspect_lasso_bundle(dir; member, dimensions=dims)
epoch = epoch_from_day_of_year(b.time.day0; year=year(m.date))
println("bundle: ", dir, "\nmember: ", member, "\nepoch: ", epoch, "\ndimensions: ", dims)
println("time: ", b.time); println("switches: ", b.switches); println("readme: ", b.readme)
println("grid: Nz=", length(b.grid.levels), " first=", b.grid.first_level, " top=", b.grid.top_level, " uniform Δz=", b.grid.uniform_spacing, " max level offset=", b.grid.max_level_offset)
println("PROBLEMS (", length(b.problems), "):"); foreach(p -> println("  - ", p), b.problems)
println("WARNINGS (", length(b.warnings), "):"); foreach(w -> println("  - ", w), b.warnings)
isempty(b.problems) || error("bundle rejected")
println("settings: ", (; b.settings.Nx, b.settings.Ny, b.settings.Lx, b.settings.Δt, b.settings.stop_time, b.settings.radiation_interval,
                        b.settings.wind_nudging_timescale, b.settings.microphysics, b.settings.surface_flux_law, b.settings.coriolis_parameter,
                        b.settings.aerosol_supersaturation_cap, b.settings.liquid_effective_radius, b.settings.surface_emissivity))
# tiny CPU construction on the full grd column: exercises the 1077-level sounding and 431-level lsf interpolation
case = ena_lasso(dir; member, epoch, dimensions=dims, arch=CPU(), FT=Float64, Nx=8, Ny=8, Lx=800, Ly=800,
                 microphysics=:one_moment, radiation=nothing, aerosol_replenishment=nothing, stop_time=3.0, write_output=false, progress_interval=100)
println("label: ", case.config.label)
println("columns: z1=", case.columns.z[1], " T1=", case.columns.T[1], " qt1=", case.columns.qᵗ[1], " u1=", case.columns.u[1], " ztop=", case.columns.z[end], " Ttop=", case.columns.T[end], " qtop=", case.columns.qᵗ[end])
run!(case.simulation)
println("finite after 3 steps: ", all(f -> all(isfinite, interior(f)), values(Oceananigans.prognostic_fields(case.model))))
# P3-aer2 construction with the LASSO defaults (radiation off, 1 step) to check the capped aerosol path
case2 = ena_lasso(dir; member, epoch, dimensions=dims, arch=CPU(), FT=Float32, Nx=8, Ny=8, Lx=800, Ly=800, radiation=nothing,
                  stop_time=1.0, write_output=false, progress_interval=100)
println("aer2 modes: ", case2.model.microphysics.aerosol)
println("aerosol record: ", (; case2.config.N₁, case2.config.N₂, case2.config.n₁, case2.config.n₂, case2.config.activatable_fraction₁, case2.config.activatable_fraction₂, case2.config.maximum_supersaturation))
run!(case2.simulation)
println("finite aer2 after 1 step: ", all(f -> all(isfinite, interior(f)), values(Oceananigans.prognostic_fields(case2.model))))
mkpath("runs/lasso_bundle_check"); write_provenance("runs/lasso_bundle_check/provenance_cpu_check.toml", case2)
println("OK")
