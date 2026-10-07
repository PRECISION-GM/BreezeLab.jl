# Eastern North Atlantic: the public Covert et al. (2022) configuration.
# Run: julia --project cases/eastern_north_atlantic.jl
# Inputs: julia data_wrangling/fetch_covert_inputs.jl
#
# 18 July 2017, 06–12 UTC; 8.96 km square at 35 m, 192 reconstructed vertical levels.
# Prescribed surface fluxes, simple longwave radiation, no mean-wind nudging,
# fixed 0.5 s timestep. This is the public development case, not official LASSO
# or the larger, nine-hour published experiment. See docs/cases/ena.md.
using BreezeLab, Oceananigans, Oceananigans.Units

function eastern_north_atlantic(;
    arch = GPU(),
    FT = Float32,
    microphysics = :p3_n75,              # :one_moment, :p3_n75, or :p3_aer2
    Nx = 256,
    Ny = 256,
    stop_time = 6hours,
    data_dir = joinpath(@__DIR__, "..", "data", "covert2022_bin"),
    output_dir = joinpath(@__DIR__, "..", "output", "ena_$(microphysics)"),
    kwargs...)

    # The aerosol-coupled ENA member uses a prescribed reservoir. This is not
    # the evolving aerosol lifecycle required by SEA STARR.
    aerosol_replenishment = microphysics === :p3_aer2 ? :diagnostic_ccn : nothing

    return ena_simulation(data_dir;
                          protocol = :covert_public_bin,
                          arch, FT, microphysics, aerosol_replenishment,
                          Nx, Ny, Lx = 35Nx, Ly = 35Ny,
                          stop_time, output_dir, kwargs...)
end

if abspath(PROGRAM_FILE) == @__FILE__
    case = eastern_north_atlantic()
    mkpath(case.config.output_dir)
    write_provenance(joinpath(case.config.output_dir, "provenance.toml"), case)
    run!(case.simulation)
end
