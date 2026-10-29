#####
##### Standard output writers of the SAM-forced cases (ENA, TRACER–DP-SCREAM).
#####

"""
    add_output_writers!(simulation; output_dir, output_prefix, profile_interval, timeseries_interval,
                        slice_interval, slice_height, extra_slice_heights = (;), y_averaged_sections = false,
                        section_fields = (;), energy_budget_series = false, in_cloud_threshold = 1e-5)

Hourly-averaged profiles, domain-mean time series and the `_slices.jld2` snapshots, which always
hold `qᶜˡ_xz`, `qʳ_xz`, `w_xz`, `θ_xz` (at y = Ly/2), `qᶜˡ_xy`, `w_xy` (at `slice_height`), and the
`lwp` and surface `rain` maps. Optional visualization slices:

- `extra_slice_heights`: more xy slices, e.g. `(; w = (500, 750), qᶜˡ = (1100,))` writes
  `w_xy_500m`, `w_xy_750m` and `qᶜˡ_xy_1100m` (the nearest cell center at or above each height).
- `y_averaged_sections = true`: y-averaged xz sections `qᶜˡ_xz_yavg`, `qʳ_xz_yavg`, `w_xz_yavg`.
- `section_fields`: a `NamedTuple` of further 3D fields or operations (process rates, for example),
  each written like `qᶜˡ`: `name_xz`, `name_xz_yavg` (with `y_averaged_sections`) and the xy slices
  listed under its name in `extra_slice_heights`.
"""
function add_output_writers!(simulation; output_dir, output_prefix, profile_interval,
                             timeseries_interval, slice_interval, slice_height,
                             extra_slice_heights = (;),
                             y_averaged_sections = false,
                             section_fields = (;),
                             energy_budget_series = false,
                             in_cloud_threshold = 1e-5)
    model = simulation.model
    grid = model.grid
    u, v, w = model.velocities
    μ = model.microphysical_fields
    qᶜˡ = cloud_liquid(model)
    qʳ = rain_mass_fraction(model)
    qᵛ = μ.qᵛ
    θ = liquid_ice_potential_temperature(model)
    T = model.temperature
    # Breeze's static-energy *diagnostic* s = cᵖᵐ T + gz − ℒˡqˡ − ℒⁱqⁱ: its column integral, with the water
    # paths, reconstructs the column moist enthalpy of the budgets (tracer_dp_scream_budget)
    s = Breeze.AtmosphereModels.Diagnostics.StaticEnergy(model)

    profile_fields = (; u, v, w² = w^2, uw = u * w, vw = v * w, θ, T, s, qᵛ, qᶜˡ, qʳ,
                        cloud_fraction = cloud_fraction_profile(model),
                        total_cloud_fraction = total_cloud_fraction_profile(model))
    haskey(μ, :nᶜˡ) && (profile_fields = merge(profile_fields, (; nᶜˡ = μ.nᶜˡ)))
    haskey(μ, :nᵃ) && (profile_fields = merge(profile_fields, (; nᵃ = μ.nᵃ)))
    haskey(μ, :qⁱ) && (profile_fields = merge(profile_fields, (; qⁱ = μ.qⁱ)))
    if !isnothing(model.radiation)
        profile_fields = merge(profile_fields, (; radiative_flux_divergence = model.radiation.flux_divergence))
    end
    ρ = reference_density(model)
    # Fields that are already horizontal means (cloud fraction) are written as they are:
    # Oceananigans main refuses to average over dimensions a field no longer has.
    profiles = NamedTuple(name => (horizontally_reduced(f) ? f : Average(f, dims=(1, 2)))
                          for (name, f) in pairs(profile_fields))

    in_cloud = energy_budget_series ? in_cloud_droplet_number(model; threshold = in_cloud_threshold) : nothing
    isnothing(in_cloud) || (profiles = merge(profiles, (; cloudy_droplet_number = in_cloud.profile)))

    simulation.output_writers[:profiles] =
        JLD2Writer(model, profiles; filename = joinpath(output_dir, output_prefix * "_profiles.jld2"),
                   schedule = AveragedTimeInterval(profile_interval), overwrite_files = true)

    lwp = liquid_water_path(model; species=:cloud)
    rwp = liquid_water_path(model; species=:rain)
    rain = surface_rain_flux(model)
    timeseries = (; lwp = Average(lwp, dims=(1, 2)),
                    rwp = Average(rwp, dims=(1, 2)),
                    cloud_fraction = cloud_fraction(model),
                    rain_flux = Average(rain, dims=(1, 2)),
                    # water and energy budget terms: precipitable water, column static energy,
                    # and the column-integrated radiative flux convergence (W m⁻², positive warming)
                    precipitable_water = Average(precipitable_water(model), dims=(1, 2)),
                    column_static_energy = Average(Field(Integral(ρ * s, dims=3)), dims=(1, 2)))
    if haskey(μ, :qⁱ)
        timeseries = merge(timeseries, (; iwp = Average(ice_water_path(model), dims=(1, 2)),
                                          ice_flux = Average(surface_ice_flux(model), dims=(1, 2))))
    end
    if !isnothing(model.radiation)
        timeseries = merge(timeseries, (; column_radiative_heating = Average(Field(Integral(model.radiation.flux_divergence, dims=3)), dims=(1, 2))))
    end

    if energy_budget_series
        # SAM's SHF/LHF, LWNS/SWNS, LWNT/SWNT (here at the LES top) and the in-cloud droplet number
        fluxes = surface_heat_fluxes(model)
        timeseries = merge(timeseries, (; surface_sensible_heat_flux = Average(fluxes.sensible, dims=(1, 2)),
                                          surface_latent_heat_flux = Average(fluxes.latent, dims=(1, 2))))
        radiative = radiative_boundary_fluxes(model.radiation, grid)
        timeseries = merge(timeseries, NamedTuple(name => Average(f, dims=(1, 2)) for (name, f) in pairs(radiative)))
        isnothing(in_cloud) || (timeseries = merge(timeseries, (; in_cloud_droplet_number = in_cloud.mean,
                                                                  cloudy_volume_fraction = in_cloud.cloudy_fraction)))
    end

    simulation.output_writers[:timeseries] =
        JLD2Writer(model, timeseries; filename = joinpath(output_dir, output_prefix * "_timeseries.jld2"),
                   schedule = TimeInterval(timeseries_interval), overwrite_files = true)

    z = Array(znodes(grid, Center()))
    level(height) = min(searchsortedfirst(z, height), length(z))
    k = level(slice_height)
    j = max(1, size(grid, 2) ÷ 2)
    slices = (; qᶜˡ_xz = view(qᶜˡ, :, j, :), qʳ_xz = view(qʳ, :, j, :), w_xz = view(w, :, j, :),
                θ_xz = section_slice(θ, (:, j, :)),
                qᶜˡ_xy = view(qᶜˡ, :, :, k), w_xy = view(w, :, :, k), lwp, rain)

    # Visualization extras: more xy heights, y-averaged sections, and the `section_fields` hook
    sectioned = merge((; qᶜˡ, qʳ, w, θ), section_fields)
    for (name, f) in pairs(section_fields)
        slices = merge(slices, NamedTuple{(Symbol(name, :_xz),)}((section_slice(f, (:, j, :)),)))
    end
    if y_averaged_sections
        for name in (:qᶜˡ, :qʳ, :w, keys(section_fields)...)
            slices = merge(slices, NamedTuple{(Symbol(name, :_xz_yavg),)}((Average(sectioned[name], dims=2),)))
        end
    end
    for (name, heights) in pairs(extra_slice_heights)
        haskey(sectioned, name) ||
            throw(ArgumentError("extra_slice_heights names $name; available: $(join(keys(sectioned), ", "))"))
        for height in heights
            label = Symbol(name, :_xy_, round(Int, height), :m)
            slices = merge(slices, NamedTuple{(label,)}((section_slice(sectioned[name], (:, :, level(height))),)))
        end
    end

    simulation.output_writers[:slices] =
        JLD2Writer(model, slices; filename = joinpath(output_dir, output_prefix * "_slices.jld2"),
                   schedule = TimeInterval(slice_interval), overwrite_files = true)
    return nothing
end

# A slice of a Field is a view; any other operation (a diagnostic KernelFunctionOperation, say)
# becomes a Field computed on the slice only.
section_slice(f::Oceananigans.Fields.Field, indices) = view(f, indices...)
section_slice(f, indices) = Field(f; indices)

horizontally_reduced(f) = (loc = Oceananigans.Fields.location(f); loc[1] === Nothing && loc[2] === Nothing)

