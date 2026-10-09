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
# `julia --project cases/ena_covert.jl` on an NVIDIA GPU.

using BreezeLab
using Breeze
using Breeze.Microphysics.PredictedParticleProperties: CloudDroplets
using Oceananigans, Oceananigans.Units
using CairoMakie

# ## Build the simulation
#
# The package constructor returns a case without advancing its clock. Precision
# is set once through Oceananigans' default, and the microphysics is a Breeze
# object: here P3 with the observed 75 cm⁻³ droplet number. An aerosol-coupled
# member would pass `P3Microphysics(; cloud, aerosol = covert_aerosol(; reference_density))`,
# with `reference_density = first_level_reference_density(data_dir, z_faces)`.

Oceananigans.defaults.FloatType = Float32
arch = GPU()
Nx, Ny = 256, 256
z_faces = covert_public_bin_vertical_faces()
stop_time = 6hours
timeseries_interval = 60seconds
profile_interval = 1hour
slice_interval = 1hour

microphysics = P3Microphysics(; cloud = CloudDroplets(; number_concentration = 75e6))

output_dir = joinpath(pkgdir(BreezeLab), "output", "ena_covert")
case = ena_covert(; arch, Nx, Ny, z_faces, stop_time, microphysics,
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
#md cp(joinpath(output_dir, "cloud_water_and_rain.png"), "ena_covert.png"; force=true);
fig #src

#md # ![Cloud water path and surface rain rate](ena_covert.png)
