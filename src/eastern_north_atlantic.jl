"""
    eastern_north_atlantic(; arch=GPU(), FT=Float32, microphysics=:p3_n75,
                            Nx=256, Ny=256, stop_time=6hours,
                            data_dir, output_dir, kwargs...)

Build the public Covert ENA development case for 18 July 2017, starting at 06 UTC.
Return a named tuple containing `simulation`, `model`, `config`, and input metadata,
without advancing the simulation. Call `run!(case.simulation)` to run it.

Defaults use 35 m horizontal spacing, 192 reconstructed vertical levels, prescribed
surface fluxes, simple longwave radiation, no mean-wind nudging, and a fixed 0.5 s
step. Choose `:one_moment`, `:p3_n75`, `:p3_covert_n75` (prognostic aerosol activating
the case's observed 75 cm⁻³) or `:p3_aer2` (LASSO aer2 aerosol sensitivity); the aerosol
members use diagnostic CCN replenishment. Extra keywords are passed to [`ena_simulation`](@ref).

Inputs default to `data/covert2022_bin` in the package checkout; fetch them with
`julia data_wrangling/fetch_covert_inputs.jl`. Outputs default to
`output/ena_<microphysics>` in the working directory. Input acquisition and
provenance writing are explicit caller actions.

This is the public development configuration, not the larger published experiment
or official LASSO. Use [`ena_simulation`](@ref) with `protocol=:lasso_ena_official`
and an official input bundle for LASSO.
"""
function eastern_north_atlantic(;
    arch = GPU(),
    FT = Float32,
    microphysics = :p3_n75,              # :one_moment, :p3_n75, or :p3_aer2
    Nx = 256,
    Ny = 256,
    stop_time = 6hours,
    data_dir = joinpath(@__DIR__, "..", "data", "covert2022_bin"),
    output_dir = joinpath(pwd(), "output", "ena_$(microphysics)"),
    kwargs...)

    # The aerosol-coupled ENA members use a prescribed (diagnostic-CCN) reservoir. This is not
    # the evolving aerosol lifecycle required by SEA STARR. `:p3_aer2` is the LASSO aer2
    # sensitivity; `:p3_covert_n75` is the case-consistent aerosol (75 cm⁻³ activatable).
    aerosol_replenishment = microphysics ∈ (:p3_aer2, :p3_covert_n75) ? :diagnostic_ccn : nothing

    return ena_simulation(data_dir;
                          protocol = :covert_public_bin,
                          arch, FT, microphysics, aerosol_replenishment,
                          Nx, Ny, Lx = 35Nx, Ly = 35Ny,
                          stop_time, output_dir, kwargs...)
end
