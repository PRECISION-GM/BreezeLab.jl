# Fast public-input CPU smoke test. Reduced grid and duration are recorded as overrides.
using BreezeLab, Breeze, Oceananigans, Statistics, TOML

data = joinpath(@__DIR__, "..", "data", "covert2022_bin")
output = isempty(ARGS) ? joinpath(@__DIR__, "..", "output", "cpu_smoke") : only(ARGS)
case = ena_simulation(data; protocol=:covert_public_bin, arch=CPU(), FT=Float64,
                      Nx=8, Ny=8, Lx=280, Ly=280, z_faces=collect(range(0, 6000, length=25)),
                      microphysics=:one_moment, stop_time=4.0, output_dir=output,
                      timeseries_interval=1.0, profile_interval=4.0, slice_interval=4.0)
mkpath(output)
write_provenance(joinpath(output, "provenance.toml"), case)
run!(case.simulation)
finite = all(f -> all(isfinite, interior(f)), values(Oceananigans.prognostic_fields(case.model)))
summary = Dict("finite" => finite, "time_seconds" => case.model.clock.time,
               "iterations" => case.model.clock.iteration,
               "mean_lwp_kg_m2" => mean(Array(interior(liquid_water_path(case.model; species=:cloud)))),
               "purpose" => "CPU execution check; not a cloud-physics validation")
open(joinpath(output, "summary.toml"), "w") do io
    TOML.print(io, summary)
end
finite || error("non-finite prognostic field in CPU smoke test")
@info "CPU smoke test passed" output summary
