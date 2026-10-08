#####
##### SEA STARR CTRL / N100 / N030: driver-backed Breeze case constructor.
#####
##### Protocol (SEA_STARR_Model_Setup_Specifications.pdf + DEPHY drivers, Zenodo 22241697):
#####   50 m horizontal spacing, 384² (or 192²) columns; 10 m vertical spacing to 2500 m then
#####   10 % stretching to 6500 m; Δt = 1 s; 2017-08-15 21 UTC → 2017-08-18 15 UTC (66 h);
#####   prescribed subsidence (full-field vertical advection), geostrophic wind, bulk surface
#####   fluxes from the trajectory SST (z₀ = 10⁻⁴ m), inversion-following nudging of θₗ, qₜ,
#####   Nₐ (τ = 1800 s) and u, v (τ = 10800 s) above max(zᵢ) + 100 m with a 200 m ramp,
#####   interactive radiation, one-mode κ-Köhler aerosol (d = 185 nm, σ = 1.5, κ = 0.2) with a
#####   70 cm⁻² s⁻¹ surface source.
#####
##### See docs/cases/sea_starr.md for the input inventory and the list of labelled
##### departures (no aerosol optics, radiation column ending at the LES top, fixed solar
##### coordinate, linear nudging ramp, proportional evaporation regeneration).
#####

using Dates: Dates, DateTime
using Oceananigans
using Oceananigans.Units
using Oceananigans.Fields: interior
using Oceananigans.Grids: znodes
using Oceananigans.OutputWriters: Checkpointer
using Breeze
using Breeze.Microphysics.PredictedParticleProperties: CloudDroplets

const SEA_STARR_MEMBERS = (:CTRL, :N100, :N030)
const SEA_STARR_DRIVER_FILES = Dict(:CTRL => "SEA_STARR_CTRL_SCM_driver.nc",
                                    :N100 => "SEA_STARR_N100_SCM_driver.nc",
                                    :N030 => "SEA_STARR_N030_SCM_driver.nc")

"""
    sea_starr_driver_path(member; data_dir)

Path of the DEPHY driver of `member` (`:CTRL`, `:N100`, `:N030`) in `data_dir`.
"""
function sea_starr_driver_path(member; data_dir)
    member ∈ SEA_STARR_MEMBERS || throw(ArgumentError("SEA STARR member must be one of $(SEA_STARR_MEMBERS), got $member"))
    return joinpath(data_dir, SEA_STARR_DRIVER_FILES[member])
end

"""
    sea_starr(; member=:CTRL, arch=GPU(), FT=Float32, data_dir, kwargs...)

Build the SEA STARR experiment `member` from its DEPHY driver and return
`(; simulation, model, grid, driver, config, inputs, ...)` without advancing it.

Keyword arguments (protocol defaults):

- `data_dir`: directory holding the three `SEA_STARR_*_SCM_driver.nc` files
  (default `data/seastarr_22241697` in the package; stage them with
  `julia data_wrangling/fetch_manifest.jl cases/seastarr/inputs.toml data/seastarr_22241697`)
- `Nx = 192, Ny = 192`, `Δx = Δy = 50`, `z_faces = sea_starr_vertical_faces()` (Nz = 288)
- `stop_time = 66hours`, `Δt = 1`, `max_Δt = 1` (fixed 1 s step unless the CFL wizard must shorten it), `cfl = 0.7`
- `initialization = :equilibrium`: partition (θₗ, qₜ) into (T, qᵛ, qᶜˡ) with Breeze's warm-phase
  saturation adjustment and activate the fraction `initial_activated_fraction` (default 1) of the
  aerosol in cloudy cells; `:condensate_free` starts from the supersaturated vapor and lets P3 condense
- `aerosol_surface_flux = 7e5` m⁻² s⁻¹ (70 cm⁻² s⁻¹); `regeneration = true`
- `mean_diameter = 185e-9`, `geometric_std = 1.5`, `kappa = 0.2`
- `radiation = :rrtmgp` (all-sky LW+SW, `radiation_interval = 60` s) or `nothing`; `ozone = :driver`
- `liquid_effective_radius = 10e-6`, `surface_albedo = 0.07`, `surface_emissivity = 0.98`
- `nudging = true`, `nudging_offset = 100`, `nudging_ramp = 200`, `inversion_search_top = 4000` m
- `vertical_advection = :full_field` or `nothing`; `geostrophic = true`; `surface = :bulk_sst` or `nothing`
- `sponge = SAMSponge(damping_depth_fraction=0.15)`: numerical damping of the top 15 % (not in the protocol)
- `closure = :smagorinsky_lilly`, `advection_order = 5`, `perturbation = InitialPerturbation()`
- output: `output_dir`, `output_prefix`, `statistics_interval = 15minutes`, `timeseries_interval = 60`,
  `fields_2d_interval = 15minutes`, `fields_3d_interval = 1hour`, `checkpoint_interval = 3hours` (or `nothing`),
  `write_output = true`, `progress_interval = 10minutes`
"""
function sea_starr(; member = :CTRL,
                     arch = GPU(),
                     FT = Float32,
                     data_dir = joinpath(dirname(dirname(pathof(BreezeLab))), "data", "seastarr_22241697"),
                     driver_path = sea_starr_driver_path(member; data_dir),
                     Nx = 192, Ny = 192, Δx = 50.0, Δy = 50.0,
                     z_faces = sea_starr_vertical_faces(),
                     stop_time = 66hours,
                     Δt = 1.0, max_Δt = 1.0, cfl = 0.7,
                     initialization = :equilibrium,
                     initial_activated_fraction = 1.0,
                     aerosol_surface_flux = 7e5,
                     regeneration = true,
                     mean_diameter = 185e-9, geometric_std = 1.5, kappa = 0.2,
                     activation_timescale = 1.0,
                     radiation = :rrtmgp,
                     radiation_interval = 60,
                     ozone = :driver,
                     liquid_effective_radius = 10e-6,
                     ice_effective_radius = 30e-6,
                     surface_albedo = 0.07,
                     surface_emissivity = 0.98,
                     CO₂ = 405e-6, CH₄ = 1.85e-6, N₂O = 330e-9,
                     nudging = true,
                     nudging_offset = 100.0,
                     nudging_ramp = 200.0,
                     inversion_search_top = 4000.0,
                     thermodynamic_nudging_timescale = nothing,  # nothing → driver (1800 s)
                     wind_nudging_timescale = nothing,           # nothing → driver (10800 s)
                     vertical_advection = :full_field,
                     geostrophic = true,
                     surface = :bulk_sst,
                     roughness_length = nothing,                 # nothing → driver z0
                     gustiness = 0.1,
                     sponge = SAMSponge(damping_depth_fraction=0.15),
                     closure = :smagorinsky_lilly,
                     advection_order = 5,
                     perturbation = InitialPerturbation(),
                     label = "SEA STARR $(member) (Breeze; exploratory CTRL: see config.departures)",
                     output_dir = joinpath(pwd(), "output", "sea_starr_$(lowercase(string(member)))"),
                     output_prefix = "sea_starr_$(lowercase(string(member)))",
                     statistics_interval = 15minutes,
                     timeseries_interval = 60,
                     fields_2d_interval = 15minutes,
                     fields_3d_interval = 1hour,
                     checkpoint_interval = 3hours,
                     write_output = true,
                     progress_interval = 10minutes)

    Oceananigans.defaults.FloatType = FT
    closure = closure === :smagorinsky_lilly ? SmagorinskyLilly(FT) : closure
    initialization ∈ (:equilibrium, :condensate_free) ||
        throw(ArgumentError("initialization must be :equilibrium or :condensate_free, got $initialization"))

    #####
    ##### Driver
    #####

    driver = read_dephy_driver(driver_path)
    expected_case = "SEA_STARR/$(member)"
    driver.case == expected_case ||
        throw(ArgumentError("driver $(driver_path) is for case $(driver.case), not $expected_case; no other member's driver is substituted"))
    epoch = driver.start_time
    τ_thermo = something(thermodynamic_nudging_timescale, 1 / driver.nudging_rates.thetal)
    τ_wind = something(wind_nudging_timescale, 1 / driver.nudging_rates.u)
    z₀ = something(roughness_length, driver.roughness_length)
    latitude = driver.latitude
    longitude = driver.longitude

    #####
    ##### Grid, reference state, dynamics
    #####

    Nz = length(z_faces) - 1
    Lx = Nx * Δx
    Ly = Ny * Δy
    grid = RectilinearGrid(arch, FT; size=(Nx, Ny, Nz), x=(0, Lx), y=(0, Ly), z=z_faces,
                           halo=(5, 5, 5), topology=(Periodic, Periodic, Bounded))
    constants = ThermodynamicConstants(FT)
    constants64 = ThermodynamicConstants(Float64)

    reference_state = ReferenceState(grid, constants;
                                     base_pressure = driver.surface_pressure,
                                     potential_temperature = z -> driver_initial_profile(driver, :thetal, z),
                                     vapor_mass_fraction = z -> driver_initial_profile(driver, :qt, z))
    dynamics = AnelasticDynamics(reference_state)
    coriolis = FPlane(FT; latitude)

    z_centers = Array(znodes(grid, Center()))
    ρᵣ = Array(interior(reference_state.density, 1, 1, :))
    pᵣ = Array(interior(reference_state.pressure, 1, 1, :))

    #####
    ##### Microphysics: P3 with the κ-Köhler mode and a prognostic reservoir
    #####

    na₀ = driver_initial_profile(driver, :na, z_centers)        # kg⁻¹, height-varying
    aerosol = kappa_aerosol_activation(FT; number_mixing_ratio = maximum(na₀),
                                       mean_radius = mean_diameter / 2, geometric_std, kappa,
                                       thermodynamic_constants = constants, activation_timescale)
    cloud = CloudDroplets(FT; number_concentration = 100e6)      # only the construction-time DSD shape with aerosol present
    microphysics_model = P3Microphysics(FT; cloud, aerosol)
    moisture_name = Breeze.AtmosphereModels.moisture_specific_name(microphysics_model)   # :qᵛ
    momentum_advection = WENO(order=advection_order)
    scalar_advection = scalar_advection_schemes(advection_order, microphysics_model, moisture_name; energy_name=:ρθ)

    #####
    ##### Large-scale forcing
    #####

    w_ls = driver_profile_time_series(grid, driver, :w, z_centers)
    ug = driver_profile_time_series(grid, driver, :ug, z_centers)
    vg = driver_profile_time_series(grid, driver, :vg, z_centers)
    θ_nud = driver_profile_time_series(grid, driver, :thetal_nud, z_centers)
    q_nud = driver_profile_time_series(grid, driver, :qt_nud, z_centers)
    na_nud = driver_profile_time_series(grid, driver, :na_nud, z_centers)
    u_nud = driver_profile_time_series(grid, driver, :u_nud, z_centers)
    v_nud = driver_profile_time_series(grid, driver, :v_nud, z_centers)

    geostrophic_forcing = geostrophic ? time_varying_geostrophic_forcings(ug, vg) : (; u=nothing, v=nothing)
    vadv = vertical_advection === :full_field ? LargeScaleVerticalAdvection(w_ls) :
           isnothing(vertical_advection) ? nothing :
           throw(ArgumentError("vertical_advection must be :full_field or nothing (SubsidenceForcing mis-weights microphysical fields in the pinned Breeze)"))

    # The nudging mask is shared by all nudged variables and updated once per step.
    mask = Field{Nothing, Nothing, Center}(grid)
    nudge(target, τ) = nudging ? InversionFollowingNudging(target, mask; timescale=τ) : nothing

    # Aerosol: surface source and evaporation regeneration (plus nudging and subsidence below).
    aerosol_source = aerosol_surface_flux > 0 ? SurfaceAerosolSource(aerosol_surface_flux) : nothing
    evaporation_rate = nothing
    regen = (; nᶜˡ=nothing, nᵃ=nothing)
    if regeneration
        evaporation_rate = CenterField(grid)   # filled by EvaporationRateUpdater from the model state
        regen = EvaporationRegeneration(evaporation_rate)
    end

    compact(args...) = Tuple(a for a in args if !isnothing(a))
    forcing = Dict{Symbol, Tuple}()
    forcing[:u] = compact(geostrophic_forcing.u, nudge(u_nud, τ_wind), vadv, sponge)
    forcing[:v] = compact(geostrophic_forcing.v, nudge(v_nud, τ_wind), vadv, sponge)
    forcing[:w] = compact(sponge)
    forcing[:θ] = compact(nudge(θ_nud, τ_thermo), vadv)
    forcing[moisture_name] = compact(nudge(q_nud, τ_thermo), vadv)
    for ρname in Breeze.AtmosphereModels.prognostic_field_names(microphysics_model)
        name = Symbol(string(ρname)[nextind(string(ρname), 1):end])
        name === moisture_name && continue
        forcing[name] = compact(vadv)
    end
    forcing[:nᵃ] = (forcing[:nᵃ]..., compact(nudge(na_nud, τ_thermo), aerosol_source, regen.nᵃ)...)
    forcing[:nᶜˡ] = (forcing[:nᶜˡ]..., compact(regen.nᶜˡ)...)
    forcing = NamedTuple(name => value for (name, value) in forcing if !isempty(value))

    #####
    ##### Surface
    #####

    Tₛ = Field{Center, Center, Nothing}(grid)
    set!(Tₛ, FT(driver.sst[1]))
    sst_updater = SeaSurfaceTemperatureUpdater(Tₛ, driver.times, FT.(driver.sst))
    boundary_conditions = if surface === :bulk_sst
        bulk_surface_flux_boundary_conditions(grid, Tₛ; moisture_name, roughness_length=z₀, gustiness)
    elseif isnothing(surface)
        NamedTuple()
    else
        throw(ArgumentError("surface must be :bulk_sst or nothing, got $surface"))
    end

    #####
    ##### Radiation
    #####

    o3_profile = ozone === :driver ? (z -> interpolate_profile(driver.z, view(driver.forcing.o3, :, 1), z)) :
                 ozone isa Number || ozone isa Function ? ozone :
                 throw(ArgumentError("ozone must be :driver, a number or a function of z"))
    background_atmosphere = BackgroundAtmosphere(; CO₂, CH₄, N₂O, O₃ = o3_profile)
    radiation_model = if radiation === :rrtmgp
        RadiativeTransferModel(grid, AllSkyOptics(), constants;
                               surface_temperature = Tₛ,
                               surface_albedo, surface_emissivity, background_atmosphere,
                               solar_position = ApparentSolarPosition(; coordinate=(longitude, latitude), epoch),
                               schedule = TimeInterval(radiation_interval),
                               liquid_effective_radius = ConstantRadiusParticles(liquid_effective_radius),
                               ice_effective_radius = ConstantRadiusParticles(ice_effective_radius))
    elseif isnothing(radiation)
        nothing
    else
        throw(ArgumentError("radiation must be :rrtmgp or nothing, got $radiation"))
    end

    #####
    ##### Model
    #####

    model = AtmosphereModel(grid; formulation = :LiquidIcePotentialTemperature, dynamics, coriolis, closure,
                            microphysics = microphysics_model, radiation = radiation_model,
                            momentum_advection, scalar_advection, forcing, boundary_conditions,
                            thermodynamic_constants = constants)

    #####
    ##### Initial condition
    #####

    θₗ₀ = driver_initial_profile(driver, :thetal, z_centers)
    qₜ₀ = driver_initial_profile(driver, :qt, z_centers)
    u₀ = driver_initial_profile(driver, :u, z_centers)
    v₀ = driver_initial_profile(driver, :v, z_centers)
    ϵ = perturbation_array(Nx, Ny, z_centers, perturbation)
    column(values) = reshape(values, 1, 1, Nz)
    δT = perturbation.amplitude_T
    δq = perturbation.amplitude_q
    T = zeros(Nz); qᵛ = zeros(Nz); qᶜˡ = zeros(Nz)
    for k in 1:Nz
        T[k], qᵛ[k], qᶜˡ[k] = saturation_partition(θₗ₀[k], qₜ₀[k], pᵣ[k]; constants=constants64)
    end
    cloudy = qᶜˡ .> 0
    if initialization === :equilibrium
        nᶜˡ₀ = initial_activated_fraction .* na₀ .* cloudy
        nᵃ₀ = na₀ .- nᶜˡ₀
        set!(model; T = column(T) .+ δT .* ϵ,
                    qᵛ = max.(0, column(qᵛ) .+ δq .* ϵ),
                    qᶜˡ = repeat(column(qᶜˡ), Nx, Ny, 1),
                    nᶜˡ = repeat(column(nᶜˡ₀), Nx, Ny, 1),
                    nᵃ = repeat(column(nᵃ₀), Nx, Ny, 1),
                    u = repeat(column(u₀), Nx, Ny, 1), v = repeat(column(v₀), Nx, Ny, 1))
    else
        Tᵈ = [Breeze.Thermodynamics.temperature(Breeze.Thermodynamics.LiquidIcePotentialTemperatureState{Float64}(θₗ₀[k], Breeze.Thermodynamics.MoistureMassFractions(qₜ₀[k]), 1e5, pᵣ[k]), constants64) for k in 1:Nz]
        set!(model; T = column(Tᵈ) .+ δT .* ϵ,
                    qᵛ = max.(0, column(qₜ₀) .+ δq .* ϵ),
                    nᵃ = repeat(column(na₀), Nx, Ny, 1),
                    u = repeat(column(u₀), Nx, Ny, 1), v = repeat(column(v₀), Nx, Ny, 1))
    end

    #####
    ##### Simulation and callbacks
    #####

    simulation = Simulation(model; Δt, stop_time)
    conjure_time_step_wizard!(simulation; cfl, max_Δt)
    Oceananigans.Diagnostics.erroring_NaNChecker!(simulation)
    add_callback!(simulation, sst_updater, IterationInterval(1))
    mask_updater = InversionMaskUpdater(model; mask, offset=nudging_offset, ramp_depth=nudging_ramp, z_max=inversion_search_top)
    mask_updater(simulation)                      # mask for the first step
    nudging && add_callback!(simulation, mask_updater, IterationInterval(1))
    if regeneration
        rate_field = cloud_evaporation_rate_field(model)
        rate_updater = EvaporationRateUpdater(rate_field, evaporation_rate)
        rate_updater(simulation)
        add_callback!(simulation, rate_updater, IterationInterval(1))
    end
    add_callback!(simulation, ProgressMessenger(model), TimeInterval(progress_interval))

    if write_output
        mkpath(output_dir)
        add_sea_starr_output_writers!(simulation; output_dir, output_prefix, statistics_interval, timeseries_interval,
                                      fields_2d_interval, fields_3d_interval, mask_updater)
        if !isnothing(checkpoint_interval)
            simulation.output_writers[:checkpointer] =
                Checkpointer(model; schedule=TimeInterval(checkpoint_interval), dir=output_dir,
                             prefix=output_prefix * "_checkpoint", cleanup=true)
        end
    end

    departures = (
        "no aerosol optical properties in Breeze RRTMGP (protocol: single-scatter albedo 0.85 at 550 nm); aerosol absorption is absent",
        "RRTMGP column ends at the LES top ($(z_faces[end]) m); no atmosphere above it and zero downwelling LW at the top",
        "fixed solar coordinate (driver lat/lon) instead of the moving trajectory",
        "nudging ramp linear over $(nudging_ramp) m (Blossey et al. 2013 form not verified); nudging acts on the horizontal mean; θ nudged to thetal_nud, qᵛ to qt_nud",
        "κ-Köhler activation of the local reservoir replaces Breeze's constant-mode activation (BreezeLab extension)",
        "aerosol regeneration: droplet number removed ∝ evaporated cloud mass and returned to nᵃ (P3 default discards it); no regeneration from rain evaporation; no interstitial scavenging",
        "surface aerosol source applied as a bottom-cell tendency (P3 ignores flux boundary conditions on ρnᵃ)",
        "numerical upper sponge (SAMSponge over the top $(isnothing(sponge) ? 0 : sponge.parameters.damping_depth_fraction * 100) %) not in the protocol",
        "initial cloud from warm-phase saturation adjustment of (θₗ, qₜ) with activated fraction $(initial_activated_fraction) (initialization = $(initialization))",
        "surface fluxes: Breeze bulk formulae (Large & Yeager polynomials) with z₀ = $(z₀) m; protocol only states ts and z0",
    )

    config = (; label, member=string(member), arch=string(typeof(arch)), FT=string(FT), Nx, Ny, Nz, Lx, Ly, Δx, Δy,
                z_top=z_faces[end], epoch=string(epoch), end_time=string(driver.end_time), latitude, longitude,
                microphysics="P3 + KappaAerosolMode (prognostic reservoir)", mean_diameter, geometric_std, kappa,
                activation_timescale, aerosol_surface_flux, regeneration, initialization=string(initialization),
                initial_activated_fraction, initial_aerosol_max=maximum(na₀), initial_aerosol_surface=na₀[1],
                radiation=string(radiation), radiation_interval, ozone=string(ozone), liquid_effective_radius, ice_effective_radius,
                surface_albedo, surface_emissivity, CO₂, CH₄, N₂O,
                surface=string(surface), roughness_length=z₀, gustiness,
                nudging, nudging_offset, nudging_ramp, inversion_search_top,
                thermodynamic_nudging_timescale=τ_thermo, wind_nudging_timescale=τ_wind,
                vertical_advection=string(vertical_advection), geostrophic,
                sponge=isnothing(sponge) ? "nothing" : summary(sponge), closure=isnothing(closure) ? "nothing" : summary(closure),
                advection_order, stop_time, Δt_initial=Δt, max_Δt, cfl,
                time_stepping = max_Δt == Δt ? "fixed Δt = $Δt s (protocol dt = 1 s; the CFL wizard may only shorten it)" : "adaptive (initial Δt = $Δt s, cfl = $cfl, max_Δt = $max_Δt s)",
                perturbation=string(perturbation),
                write_output, output_dir=abspath(output_dir), output_prefix, statistics_interval, timeseries_interval,
                fields_2d_interval, fields_3d_interval, checkpoint_interval=something(checkpoint_interval, 0), progress_interval,
                driver_case=driver.case, driver_version=string(get(driver.attributes, "version", "")),
                departures=departures)

    inputs = (; driver = driver_path)

    return (; simulation, model, grid, driver, config, inputs, surface_temperature=Tₛ, mask_updater,
              evaporation_rate, protocol=:sea_starr, preset=member, protocol_dimensions=(Nx, Ny, Nz),
              protocol_overrides=NamedTuple(), initial_columns=(; z=z_centers, θₗ=θₗ₀, qₜ=qₜ₀, T, qᵛ, qᶜˡ, nᵃ=na₀, u=u₀, v=v₀, ρ=ρᵣ, p=pᵣ))
end

#####
##### Output
#####

function add_sea_starr_output_writers!(simulation; output_dir, output_prefix, statistics_interval, timeseries_interval,
                                       fields_2d_interval, fields_3d_interval, mask_updater)
    model = simulation.model
    grid = model.grid
    Nz = size(grid, 3)
    u, v, w = model.velocities
    μ = model.microphysical_fields
    qᶜˡ = cloud_liquid(model)
    qʳ = rain_mass_fraction(model)
    θ = liquid_ice_potential_temperature(model)
    T = model.temperature
    mean(f) = horizontally_reduced(f) ? f : Average(f, dims=(1, 2))

    profile_fields = (; u, v, w² = w^2, uw = u * w, vw = v * w, θ, T, qᵛ = μ.qᵛ, qᶜˡ, qʳ,
                        nᵃ = μ.nᵃ, nᶜˡ = μ.nᶜˡ, nʳ = μ.nʳ,
                        cloud_fraction = cloud_fraction_profile(model; threshold=1e-5),
                        nudging_mask = mask_updater.mask)
    if !isnothing(model.radiation)
        r = model.radiation
        profile_fields = merge(profile_fields, (; radiative_flux_divergence = r.flux_divergence,
                                                  lw_up = r.upwelling_longwave_flux, lw_down = r.downwelling_longwave_flux,
                                                  sw_up = r.upwelling_shortwave_flux, sw_down = r.downwelling_shortwave_flux))
    end
    profiles = NamedTuple(name => mean(f) for (name, f) in pairs(profile_fields))
    simulation.output_writers[:statistics] =
        JLD2Writer(model, profiles; filename = joinpath(output_dir, output_prefix * "_statistics.jld2"),
                   schedule = AveragedTimeInterval(statistics_interval), overwrite_files = true)

    cwp = liquid_water_path(model; species=:cloud)
    rwp = liquid_water_path(model; species=:rain)
    rain = surface_rain_flux(model)
    zi = mask_updater.inversion_height
    columns = aerosol_number_columns(model)
    timeseries = (; cwp = Average(cwp, dims=(1, 2)), rwp = Average(rwp, dims=(1, 2)),
                    cloud_fraction = cloud_fraction(model), rain_flux = Average(rain, dims=(1, 2)),
                    zi = Average(zi, dims=(1, 2)),
                    zi_max = Field(Oceananigans.Fields.Reduction(maximum!, zi, dims=(1, 2))),
                    nᵃ_column = Average(columns.nᵃ, dims=(1, 2)), nᶜˡ_column = Average(columns.nᶜˡ, dims=(1, 2)),
                    nʳ_column = Average(columns.nʳ, dims=(1, 2)), n_total_column = Average(columns.total, dims=(1, 2)))
    simulation.output_writers[:timeseries] =
        JLD2Writer(model, timeseries; filename = joinpath(output_dir, output_prefix * "_timeseries.jld2"),
                   schedule = TimeInterval(timeseries_interval), overwrite_files = true)

    fields_2d = (; cwp, rwp, rain, zi, nᵃ_column = columns.nᵃ, n_total_column = columns.total)
    simulation.output_writers[:fields_2d] =
        JLD2Writer(model, fields_2d; filename = joinpath(output_dir, output_prefix * "_2d.jld2"),
                   schedule = TimeInterval(fields_2d_interval), overwrite_files = true)

    fields_3d = (; u, v, w, T, qᵛ = μ.qᵛ, qᶜˡ, qʳ, nᵃ = μ.nᵃ, nᶜˡ = μ.nᶜˡ, nʳ = μ.nʳ)
    simulation.output_writers[:fields_3d] =
        JLD2Writer(model, fields_3d; filename = joinpath(output_dir, output_prefix * "_3d.jld2"),
                   schedule = TimeInterval(fields_3d_interval), overwrite_files = true)
    return nothing
end
