# # Eastern North Atlantic — official LASSO-ENA member
#
# This example builds, runs, and analyzes one official LASSO-ENA ensemble member from
# its ARM `samin` input bundle (DOI 10.5439/2572661). It is a different experiment
# from the Covert development benchmark in `eastern_north_atlantic.jl`: the grid,
# duration, surface fluxes (bulk from SST), radiation (RRTMGP longwave and shortwave),
# wind nudging and aerosol setting follow the member's namelist and files. The member
# identity is explicit; a missing bundle stops the script with staging instructions and
# nothing is substituted.
#
# Stage the bundle with `julia data_wrangling/stage_lasso_bundle.jl ARCHIVE.tar`, then
# run `julia --project cases/ena_lasso.jl` on an NVIDIA GPU (H100 or A100). Set
# `ENA_LASSO_MEMBER`, `ENA_LASSO_BUNDLE` and `ENA_LASSO_EPOCH` to override the member,
# the bundle directory and the UTC start.

using BreezeLab
using Oceananigans, Oceananigans.Units
using CairoMakie
using Dates

# ## Identify the member and its staged bundle
#
# `20170718era5d25x100_sbmwrm-aer2-flxsst` is the spectral-bin, warm-cloud, medium-aerosol,
# SST-flux member of 18 July 2017 forced with ERA5 (5° sampling, wind nudging only). The
# plan's original candidate `…era5s1n0d25x100_sbmwrm-…` does not exist at ARM (see
# `docs/cases/ena_lasso.md`).

member = get(ENV, "ENA_LASSO_MEMBER", "20170718era5d25x100_sbmwrm-aer2-flxsst")
bundle_dir = lasso_bundle_directory(member)
lasso_bundle_available(bundle_dir) || error(lasso_bundle_missing_message(member, bundle_dir))

# ## Inspect the bundle before building anything
#
# The namelist `day0` and the member date fix the UTC start; SAM's `day0` alone does not
# carry the year, so the epoch is formed explicitly here and then verified by the adapter.
# The documented d25x100 domain (256 × 256 × 260) is requested explicitly and checked
# against the bundle's `grd`.

inspection = inspect_lasso_bundle(bundle_dir; member, dimensions=lasso_documented_dimensions(member))
epoch = haskey(ENV, "ENA_LASSO_EPOCH") ? DateTime(ENV["ENA_LASSO_EPOCH"]) :
        epoch_from_day_of_year(inspection.time.day0; year=year(parse_lasso_member(member).date))
foreach(w -> println("note: ", w), inspection.warnings)
isempty(inspection.problems) || error("bundle cannot be reproduced:\n  - " * join(inspection.problems, "\n  - "))
println("member $member starts $epoch UTC for $(inspection.time.stop_time / 3600) h at Δt = $(inspection.time.dt) s")

# ## Build the simulation
#
# Everything below the member line is the protocol default read from the bundle; only
# the architecture, precision and output cadence are chosen here (and recorded).

arch = GPU()
FT = Float32
output_dir = joinpath(pkgdir(BreezeLab), "output", "ena_lasso", member)
case = ena_lasso(bundle_dir; member, epoch, arch, FT, output_dir, output_prefix = "ena_lasso",
                 timeseries_interval = 60seconds, profile_interval = 10minutes, slice_interval = 30minutes)
simulation = case.simulation;

# ## Run
#
# Provenance holds the input checksums, the archive checksum and download date from
# `bundle.toml`, the member identity, the SAM reference revision and every recorded
# protocol difference.

mkpath(output_dir)
write_provenance(joinpath(output_dir, "provenance.toml"), case)
run!(simulation)

# ## Analyze cloud water, rain and cloud boundaries
#
# Model diagnostics only; the comparison against the SAM member and ARM observations
# over the common window is `analysis/compare_ena_observations.jl`.

series = breezelab_timeseries(output_dir)
bounds = breezelab_cloud_boundaries(output_dir)
hours = series.seconds ./ 3600

fig = Figure(size = (900, 900))
ax_lwp = Axis(fig[1, 1], xlabel = "Hours since $(epoch) UTC", ylabel = "Cloud liquid water path (g m⁻²)")
ax_rain = Axis(fig[2, 1], xlabel = "Hours since $(epoch) UTC", ylabel = "Surface rain rate (mm hr⁻¹)")
ax_cloud = Axis(fig[3, 1], xlabel = "Hours since $(epoch) UTC", ylabel = "Cloud base / top (m)")
lines!(ax_lwp, hours, series.lwp)
lines!(ax_rain, hours, series.rain_rate)
lines!(ax_cloud, bounds.seconds ./ 3600, bounds.base, label = "base")
lines!(ax_cloud, bounds.seconds ./ 3600, bounds.top, label = "top")
axislegend(ax_cloud, position = :rb)
Label(fig[0, 1], "LASSO-ENA member $member (Breeze; not a SAM result)", fontsize = 14)
save(joinpath(output_dir, "cloud_water_rain_boundaries.png"), fig)
fig #src
