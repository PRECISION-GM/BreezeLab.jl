# # Eastern North Atlantic
#
# This example builds, runs, and analyzes the public Covert et al. (2022)
# development case for 18 July 2017, 06–12 UTC. The domain is 8.96 km square
# with 35 m horizontal spacing and 192 reconstructed vertical levels.
# It uses prescribed surface fluxes, simple longwave radiation, and no
# mean-wind nudging. This configuration differs from both official LASSO
# and the larger, nine-hour published experiment.
#
# Fetch the inputs with `julia data_wrangling/fetch_covert_inputs.jl`, then run
# `julia --project cases/eastern_north_atlantic.jl` on an NVIDIA GPU.

using BreezeLab
using Oceananigans, Oceananigans.Units
using CairoMakie

# ## Build the simulation
#
# The package constructor returns a case without advancing its clock. Here we
# choose fixed-number P3 microphysics (`:p3_n75`); `:one_moment` and `:p3_aer2`
# are also available. The aerosol-coupled member uses a prescribed CCN reservoir.

arch = GPU()
Nx, Ny = 256, 256
z_faces = covert_public_bin_vertical_faces()
stop_time = 6hours
timeseries_interval = 60seconds
profile_interval = 1hour
slice_interval = 1hour

output_dir = joinpath(pkgdir(BreezeLab), "output", "eastern_north_atlantic")
case = eastern_north_atlantic(; arch, Nx, Ny, z_faces, stop_time,
                              microphysics = :p3_n75,
                              output_dir, output_prefix = "ena",
                              timeseries_interval, profile_interval, slice_interval)
simulation = case.simulation;

# ## Run
#
# Save the configuration, input checksums, and protocol overrides alongside the
# output. This six-hour simulation is a substantial GPU run; output writers
# replace files with the same prefix when the example is rerun.

mkpath(output_dir)
write_provenance(joinpath(output_dir, "provenance.toml"), case)
run!(simulation)

# ## Analyze cloud water and rain
#
# Read the saved horizontal means. Cloud liquid water path is stored in kg m⁻²;
# multiply by 1000 for g m⁻². The downward surface rain mass flux is in
# kg m⁻² s⁻¹; multiplying by seconds per day gives mm day⁻¹ of liquid water.
# These are model diagnostics; comparison with ARM requires matching its sampling
# and retrieval definitions.

filename = joinpath(output_dir, "ena_timeseries.jld2")
lwp = FieldTimeSeries(filename, "lwp")
rain = FieldTimeSeries(filename, "rain_flux")
times = lwp.times ./ hour
cloud_water_path = [1e3 * lwp[n][1, 1, 1] for n in eachindex(times)]
rain_rate = [86400 * rain[n][1, 1, 1] for n in eachindex(times)]

fig = Figure(size = (900, 650))
ax_lwp = Axis(fig[1, 1], xlabel = "Hours since 06 UTC, 18 July 2017",
              ylabel = "Cloud liquid water path (g m⁻²)")
ax_rain = Axis(fig[2, 1], xlabel = "Hours since 06 UTC, 18 July 2017",
               ylabel = "Surface rain rate (mm day⁻¹)")
lines!(ax_lwp, times, cloud_water_path)
lines!(ax_rain, times, rain_rate)
save(joinpath(output_dir, "cloud_water_and_rain.png"), fig)
#md cp(joinpath(output_dir, "cloud_water_and_rain.png"), "eastern_north_atlantic.png"; force=true);
fig #src

#md # ![Cloud water path and surface rain rate](eastern_north_atlantic.png)
