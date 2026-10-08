# Stage the ERA5 inputs of one TRACER-MIP case with NumericalEarth's download API.
#
#   set -a; source ~/.bashrc; set +a     # or ~/.cdsapirc; credentials are never printed
#   julia --project=cases/tracer_mip data_wrangling/fetch_era5_tracer_mip.jl aug07 data/era5
#
# What is staged (per-datetime NetCDF files named by NumericalEarth's `metadata_filename`):
#   * pressure levels (37 levels, hourly, run window ± margin) over the OUTER grid's bounding box
#     padded by ERA5's default 1/2°: T, u, v, q, cloud liquid/ice, rain, snow water contents and
#     geopotential — exactly the set `PrescribedAtmosphere(region, dates, ERA5HourlyPressureLevels())`
#     loads, so the regional constructor finds every file cached;
#   * single levels over the unpadded outer box, hourly: MSLP, surface pressure, skin and sea-surface
#     temperature, 2 m T / dewpoint, 10 m winds, total precipitation, downwelling SW/LW, orography;
#   * ERA5-Land at the start hour: skin/soil temperatures and soil moisture (land initial state).
# The regions are built from `tracer_mip_grid` itself so they match the constructor's by construction.
# A TOML manifest with sizes, SHA-256 checksums, request parameters and software revisions is written
# next to the files. Credentials come from ~/.cdsapirc or CDSAPI_URL/CDSAPI_KEY; nothing is echoed.

using BreezeLab: tracer_mip_protocol, tracer_mip_case_window, tracer_mip_grid, file_sha256
using NumericalEarth
using NumericalEarth.DataWrangling: Metadata, MetadataSet, BoundingBox, default_horizontal_padding, metadata_path
using NumericalEarth.DataWrangling.ERA5: ERA5HourlyPressureLevels, ERA5HourlySingleLevel, ERA5HourlyLand, ERA5
# Upstream gap (NumericalEarth d07eb240): the CDS splitter looks up NetCDF short names with
# `nc_varnames(dataset)`, which has no ERA5-Land method, so soil variables raise a KeyError. Supply the
# land mapping that NumericalEarth itself defines; remove once upstream adds it.
ERA5.nc_varnames(::ERA5.ERA5LandDataset) = ERA5.ERA5Land_netcdf_variable_names   # more specific than the generic ERA5Dataset method
using CopernicusClimateDataStore   # activates the ERA5 download extension
using Oceananigans: CPU
using Downloads: Downloads
using Dates: DateTime, Hour, now, UTC
using TOML: TOML
using Printf: @sprintf

case = Symbol(get(ARGS, 1, "aug07"))
dir = abspath(get(ARGS, 2, joinpath(@__DIR__, "..", "data", "era5")))
margin = Hour(parse(Int, get(ARGS, 3, "6")))   # hours before/after the case window for boundary interpolation

have_credentials = isfile(joinpath(homedir(), ".cdsapirc")) ||
                   (haskey(ENV, "CDSAPI_URL") && haskey(ENV, "CDSAPI_KEY")) ||
                   isfile(joinpath(homedir(), ".config", "era5cli", "cds_key.txt"))
have_credentials || error("No CDS credentials: expected ~/.cdsapirc (url:/key:), or CDSAPI_URL and CDSAPI_KEY in the environment")

protocol = tracer_mip_protocol()
start, stop = tracer_mip_case_window(protocol, case)
dates = collect((start - margin):Hour(1):(stop + margin))

# Regions derived from the actual outer grid (Float32, CPU, reduced horizontal size: only the extent matters).
grid = tracer_mip_grid(CPU(), :outer; protocol, Nx = 8, Ny = 8, terrain_following = false)   # ≥ halo cells; only the extent matters
pressure_dataset = ERA5HourlyPressureLevels()
parent_region = BoundingBox(grid; padding = default_horizontal_padding(pressure_dataset))
child_region = BoundingBox(grid)
land_dataset = ERA5HourlyLand()
single_dataset = ERA5HourlySingleLevel()

mkpath(dir)
@info "Staging ERA5 for TRACER-MIP $case" dir start stop margin parent_region child_region

pressure_variables = (:temperature, :eastward_velocity, :northward_velocity, :specific_humidity,
                      :specific_cloud_liquid_water_content, :specific_rain_water_content,
                      :specific_cloud_ice_water_content, :specific_snow_water_content, :geopotential)
single_variables = (:mean_sea_level_pressure, :surface_pressure, :skin_temperature, :sea_surface_temperature,
                    :temperature, :dewpoint_temperature, :eastward_velocity, :northward_velocity,
                    :total_precipitation, :downwelling_shortwave_radiation, :downwelling_longwave_radiation, :topography)
land_variables = (:skin_temperature, :soil_temperature_level_1, :soil_temperature_level_2, :soil_temperature_level_3,
                  :soil_temperature_level_4, :volumetric_soil_water_layer_1, :volumetric_soil_water_layer_2,
                  :volumetric_soil_water_layer_3, :volumetric_soil_water_layer_4)

requests = (
    (label = "pressure_levels", dataset = pressure_dataset, names = pressure_variables, region = parent_region, dates = dates),
    (label = "single_levels",   dataset = single_dataset,   names = single_variables,   region = child_region,  dates = dates),
    (label = "land",            dataset = land_dataset,     names = land_variables,     region = child_region,  dates = [start]),
)

staged = Dict{String, Any}()
for r in requests
    t₀ = time()
    mset = MetadataSet(r.names...; dataset = r.dataset, dates = r.dates, region = r.region, dir)
    Downloads.download(mset)
    files = String[]
    for name in r.names, d in r.dates
        push!(files, metadata_path(Metadata(name; dataset = r.dataset, dates = d, region = r.region, dir)))
    end
    missing_files = filter(!isfile, files)
    isempty(missing_files) || error("$(r.label): $(length(missing_files)) expected files are missing, e.g. $(first(missing_files))")
    bytes = sum(filesize, files)
    @info @sprintf("%s: %d files, %.1f MB, %.0f s", r.label, length(files), bytes / 1e6, time() - t₀)
    staged[r.label] = Dict("dataset" => string(nameof(typeof(r.dataset))),
                           "variables" => [string(n) for n in r.names],
                           "region_longitude" => collect(r.region.longitude), "region_latitude" => collect(r.region.latitude),
                           "first_date" => string(first(r.dates)), "last_date" => string(last(r.dates)), "n_dates" => length(r.dates),
                           "n_files" => length(files), "bytes" => bytes,
                           "files" => Dict(basename(f) => Dict("bytes" => filesize(f), "sha256" => file_sha256(f)) for f in files))
end

manifest = Dict("case" => string(case), "generated_utc" => string(now(UTC)), "directory" => dir,
                "cds_dataset_pressure" => "reanalysis-era5-pressure-levels", "cds_dataset_single" => "reanalysis-era5-single-levels",
                "cds_dataset_land" => "reanalysis-era5-land", "time_basis" => "UTC, hourly instantaneous (accumulations are CDS hourly accumulations)",
                "pressure_levels_hPa" => Int.(pressure_dataset.pressure_levels ./ 100),
                "case_window" => Dict("start" => string(start), "stop" => string(stop), "margin_hours" => margin.value),
                "software" => Dict("NumericalEarth" => string(pkgversion(NumericalEarth)),
                                   "CopernicusClimateDataStore" => string(pkgversion(CopernicusClimateDataStore)),
                                   "julia" => string(VERSION)),
                "requests" => staged)
manifest_path = joinpath(dir, "MANIFEST_era5_$(case).toml")
open(manifest_path, "w") do io
    TOML.print(io, manifest)
end
@info "Wrote $manifest_path"
