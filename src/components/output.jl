#####
##### Standard output writers of the SAM-forced cases (ENA, TRACER–DP-SCREAM).
#####

function add_output_writers!(simulation; output_dir, output_prefix, profile_interval,
                             timeseries_interval, slice_interval, slice_height,
                             energy_budget_series = false,
                             temperature_neutral_evaporation = false,
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
    s = Breeze.AtmosphereModels.Diagnostics.StaticEnergy(model)   # a diagnostic in either formulation

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
        fluxes = surface_heat_fluxes(model; temperature_neutral_evaporation)
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
    k = searchsortedfirst(z, slice_height)
    j = max(1, size(grid, 2) ÷ 2)
    slices = (; qᶜˡ_xz = view(qᶜˡ, :, j, :), qʳ_xz = view(qʳ, :, j, :), w_xz = view(w, :, j, :),
                qᶜˡ_xy = view(qᶜˡ, :, :, k), w_xy = view(w, :, :, k), lwp, rain)

    simulation.output_writers[:slices] =
        JLD2Writer(model, slices; filename = joinpath(output_dir, output_prefix * "_slices.jld2"),
                   schedule = TimeInterval(slice_interval), overwrite_files = true)
    return nothing
end

horizontally_reduced(f) = (loc = Oceananigans.Fields.location(f); loc[1] === Nothing && loc[2] === Nothing)

