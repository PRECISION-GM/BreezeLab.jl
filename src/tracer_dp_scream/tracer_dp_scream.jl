#####
##### TRACER–DP-SCREAM: a doubly periodic Breeze case driven by the forcing of the published
##### DP-SCREAM TRACER experiment (Oware et al. 2025, doi 10.1029/2025JD044113; outputs
##### 10.5281/zenodo.15271730; run script scmlib DPxx_SCREAM_SCRIPTS/run_dpxx_scream_TRACER.csh).
#####
##### Protocol facts taken from the script, the EAMxx source at the archived commit
##### (git_version 1d551ea2b0, 2024-03-12) and the IOP file:
#####   * forcing file TRACER_iopfile_4scam.nc (ARM VARANAL ERA5-constrained, hourly, 40 levels,
#####     150-km domain around 29.75 N, 95.45 W); DP-SCREAM 200 km × 200 km at 3.33 km, 128 levels,
#####     dt = 100 s, RRTMGP every 300 s, SHOC + P3, land/ocean coupler stress;
#####   * large-scale transport = the file's 3-D advective tendencies (divT + vertdivT,
#####     divq + vertdivq) applied once; no omega subsidence (`iop_dosubsidence = false`);
#####   * `iop_srf_prop = true`: sensible/latent heat fluxes and the radiative surface
#####     temperature (Tg) prescribed from the file;
#####   * planar HOMME has `fcor = 0` and `iop_coriolis = false`: no Coriolis force;
#####   * `iop_nudge_uv = true` in the script, but the archived E3SM revision predates the
#####     DP-EAMxx wind-nudging port (E3SM 3f7eee0053, 2024-04-16): the archived runs applied
#####     **no** wind nudging. `wind_nudging_timescale = nothing` reproduces that; 10800 s is
#####     the script's intended (EAMxx default) relaxation.
#####
##### Deliberate differences of this Breeze case (recorded in the label and provenance):
#####   LES closure (Smagorinsky–Lilly) at Δx = 200 m on a 51.2 km domain instead of SHOC at
#####   3.33 km on 200 km; 160 levels to 22 km (DP-SCREAM: 128 hybrid levels to ~40 km) with a
#####   Rayleigh sponge above ~16.5 km; P3 with a prescribed droplet number (SCREAM P3: prognostic
#####   droplet number with prescribed CCN); constant bulk drag with a land roughness length
#####   instead of the ELM/ocean coupler stress; anelastic dynamics; no aerosol radiative effects.
#####

using Dates: Dates, DateTime
using Oceananigans.Units: minutes

"""
    tracer_dp_scream_vertical_faces(; Δz=50, uniform_top=2000, top=22000, Nz=160)

Cell interfaces: uniform `Δz` to `uniform_top`, then a constant-ratio stretch to `top`
(50 m → ≈ 400 m at 22 km with the defaults). Deep-convective LES grid; the DP-SCREAM
128-level hybrid grid is not reproduced.
"""
tracer_dp_scream_vertical_faces(; Δz=50.0, uniform_top=2000.0, top=22000.0, Nz=160) =
    uniform_then_stretched_faces(; Nz, top, Δz, uniform_top)

const TRACER_DP_SCREAM_DATA = package_path("data", "tracer_dp_scream")

"""
    neutral_drag_coefficient(z₁, roughness_length; κ=0.4)

`(κ / ln(z₁ / z₀))²`, the neutral log-law drag coefficient referenced to the first cell
center height `z₁`.
"""
neutral_drag_coefficient(z₁, roughness_length; κ=0.4) = (κ / log(z₁ / roughness_length))^2

"""
    tracer_dp_scream(; iop_path, start=DateTime(2022, 8, 1), stop=DateTime(2022, 8, 15), kwargs...)

Build the periodic TRACER case driven by the DP-SCREAM IOP forcing and return
`(; simulation, model, grid, config, inputs, ...)` without advancing it. Precision follows
`Oceananigans.defaults.FloatType`.

Defaults (the DP-SCREAM protocol facts are listed at the top of this file):

- `start`/`stop`: UTC `DateTime`s on IOP records. The DP-SCREAM August run started at
  2022-08-01 00 UTC and the paper's 5–15 August window is its days 5–15, so the default
  reproduces that initialization and continuous forcing history (14 days); the archived
  0.5-km run started cold at 2022-08-05 00 UTC. `stop_time` defaults to the window length.
- `arch = GPU()`, `Nx = Ny = 256`, `Δx = Δy = 200` m, `z_faces = tracer_dp_scream_vertical_faces()`
- `microphysics = P3Microphysics(; cloud = CloudDroplets(; number_concentration = 200e6))`:
  any Breeze microphysics; the default prescribes a continental 200 cm⁻³ droplet number
- `large_scale_transport = :iop_3d` (`divT + vertdivT`, `divq + vertdivq` applied once, no
  subsidence) or `:iop_horizontal`; `hold_forcing = true` (piecewise-constant hourly forcing)
- `wind_nudging_timescale = nothing` (the archived runs applied none), `coriolis = nothing`
- `radiation = true`: RRTMGP every `radiation_interval = 300` s with the IOP `Tg` as surface
  temperature, `surface_albedo` from the IOP shortwave fluxes (≈ 0.15), `surface_emissivity = 0.98`
- surface: the IOP sensible and latent heat fluxes, and a neutral log-law bulk drag for
  `roughness_length = 0.1` m referenced to the first cell center (or `drag_coefficient`),
  `gustiness = 1` m/s
- `sponge = SAMSponge(damping_depth_fraction = 0.25)`, `closure = SmagorinskyLilly()`,
  `advection_order = 5`, `Δt = 2`, `max_Δt = 5`, `cfl = 0.7`, `perturbation = InitialPerturbation()`
- outputs: profiles every 30 min, time series every 60 s, slices (animation fields) every 5 min at 3000 m,
  `checkpoint_interval = nothing` (set e.g. `6hours` for resumable runs)
"""
function tracer_dp_scream(;
    iop_path = joinpath(TRACER_DP_SCREAM_DATA, "TRACER_iopfile_4scam.nc"),
    start = DateTime(2022, 8, 1),
    stop = DateTime(2022, 8, 15),
    arch = GPU(),
    Nx = 256, Ny = 256, Δx = 200.0, Δy = Δx,
    z_faces = tracer_dp_scream_vertical_faces(),
    microphysics = P3Microphysics(; cloud = CloudDroplets(; number_concentration = 200e6)),
    large_scale_transport = :iop_3d,
    hold_forcing = true,
    wind_nudging_timescale = nothing,
    coriolis = nothing,
    radiation = true,
    radiation_interval = 300.0,
    surface_albedo = nothing,
    surface_emissivity = 0.98,
    liquid_effective_radius = 10e-6,
    ice_effective_radius = 30e-6,
    background_atmosphere = BackgroundAtmosphere(CO₂ = 405e-6, CH₄ = 1.85e-6, N₂O = 330e-9),
    roughness_length = 0.1,
    drag_coefficient = nothing,
    gustiness = 1.0,
    sponge = SAMSponge(damping_depth_fraction = 0.25),
    closure = SmagorinskyLilly(),
    advection_order = 5,
    perturbation = InitialPerturbation(),
    initialization = :condensate_free,
    initial_droplet_number = nothing,
    Δt = 2.0, max_Δt = 5.0, cfl = 0.7,
    stop_time = nothing,
    write_output = true,
    output_dir = joinpath(pwd(), "output", "tracer_dp_scream"),
    output_prefix = "tracer",
    profile_interval = 30minutes,
    timeseries_interval = 60.0,
    slice_interval = 5minutes,
    slice_height = 3000.0,
    progress_interval = 10minutes,
    checkpoint_interval = nothing,
    checkpoint_cleanup = true)

    #####
    ##### Inputs: the IOP forcing file, as SAM-format records over the window
    #####

    isfile(iop_path) || throw(ArgumentError("IOP forcing file not found: $iop_path (fetch it with `julia data_wrangling/fetch_manifest.jl cases/tracer_dp_scream/inputs.toml data/tracer_dp_scream`)"))
    iop = read_iop_forcing(iop_path)
    inputs = iop_sam_inputs(iop; start, stop, transport = large_scale_transport, hold = hold_forcing)
    (; soundings, lsf, sfc, day0, epoch) = inputs
    albedo = isnothing(surface_albedo) ? iop_surface_albedo(iop; start, stop) : surface_albedo   # night windows have none
    duration = Dates.value(stop - start) / 1000
    stop_time = something(stop_time, duration)
    stop_time ≤ duration || throw(ArgumentError("stop_time = $stop_time s exceeds the forcing window ($duration s)"))
    latitude = Float64(iop.latitude)
    longitude = Float64(iop.longitude)
    moisture_basis = :mass_fraction           # the IOP q is specific humidity
    sounding = initial_sounding(soundings, day0)
    profiles = SoundingProfiles(sounding; moisture_basis)

    #####
    ##### Grid, reference state, dynamics
    #####

    Nz = length(z_faces) - 1
    Lx, Ly = Nx * Δx, Ny * Δy
    grid = RectilinearGrid(arch; size = (Nx, Ny, Nz), x = (0, Lx), y = (0, Ly), z = z_faces,
                           halo = (5, 5, 5), topology = (Periodic, Periodic, Bounded))
    FT = eltype(grid)
    constants = ThermodynamicConstants()
    reference_state = sounding_reference_state(grid, profiles, constants)
    dynamics = AnelasticDynamics(reference_state)

    z_centers = Array(znodes(grid, Center()))
    pᵣ = Array(interior(reference_state.pressure, 1, 1, :))
    ρᵣ = Array(interior(reference_state.density, 1, 1, :))

    #####
    ##### Microphysics and advection
    #####

    check_precision(microphysics, FT)
    moisture_name = Breeze.AtmosphereModels.moisture_specific_name(microphysics)
    momentum_advection = WENO(order = advection_order)
    scalar_advection = scalar_advection_schemes(advection_order, microphysics, moisture_name; energy_name = :ρθ)

    #####
    ##### Large-scale forcing: the IOP transport tendencies, applied once (no subsidence)
    #####

    forcing_profiles = LargeScaleForcingProfiles(grid, lsf, z_centers, pᵣ; day0)
    thermodynamic = large_scale_thermodynamic_forcings(forcing_profiles.tls, forcing_profiles.qls;
                                                       microphysics, thermodynamic_constants = constants,
                                                       moisture_name, moisture_basis)
    nudge(target) = isnothing(wind_nudging_timescale) ? nothing : MeanProfileNudging(target; timescale = wind_nudging_timescale)
    forcing = compact_forcing(Dict{Symbol, Tuple}(
        :u => (nudge(forcing_profiles.uls), sponge),
        :v => (nudge(forcing_profiles.vls), sponge),
        :w => (sponge,),
        :E => (thermodynamic.E,),
        moisture_name => (thermodynamic[moisture_name],)))

    #####
    ##### Surface: IOP sensible and latent heat fluxes (iop_srf_prop), bulk drag, IOP Tg
    #####

    surface_series = surface_time_series(grid, sfc, day0)
    Tₛ = Field{Center, Center, Nothing}(grid)
    set!(Tₛ, FT(sfc.sst[1]))
    sst_updater = SeaSurfaceTemperatureUpdater(Tₛ, surface_series.times, FT.(sfc.sst))
    Cᴰ = something(drag_coefficient, neutral_drag_coefficient(z_faces[2] / 2, roughness_length))
    heat = prescribed_heat_flux_boundary_conditions(grid, sfc, day0; thermodynamic_constants = constants, moisture_name)
    drag = BulkDrag(; coefficient = FT(Cᴰ), gustiness = FT(gustiness), surface_temperature = Tₛ)
    boundary_conditions = merge((; ρu = FieldBoundaryConditions(bottom = drag),
                                   ρv = FieldBoundaryConditions(bottom = drag)), heat.bcs)

    #####
    ##### Radiation: RRTMGP LW+SW with the IOP ground temperature
    #####

    radiation_model = if radiation
        RadiativeTransferModel(grid, AllSkyOptics(), constants;
                               surface_temperature = Tₛ,
                               surface_albedo = albedo, surface_emissivity, background_atmosphere,
                               solar_position = ApparentSolarPosition(; coordinate = (longitude, latitude), epoch),
                               schedule = TimeInterval(radiation_interval),
                               liquid_effective_radius = ConstantRadiusParticles(liquid_effective_radius),
                               ice_effective_radius = ConstantRadiusParticles(ice_effective_radius))
    else
        nothing
    end

    #####
    ##### Model and initial state
    #####

    model = AtmosphereModel(grid; formulation = :LiquidIcePotentialTemperature, dynamics, coriolis, closure,
                            microphysics, radiation = radiation_model, momentum_advection, scalar_advection,
                            forcing, boundary_conditions, thermodynamic_constants = constants)

    columns = initial_state_columns(profiles, z_centers, pᵣ; constants = ThermodynamicConstants(Float64))
    set_sounding_initial_state!(model, columns; perturbation, moisture_basis, initialization, initial_droplet_number)

    #####
    ##### Simulation
    #####

    simulation = Simulation(model; Δt, stop_time)
    conjure_time_step_wizard!(simulation; cfl, max_Δt)
    Oceananigans.Diagnostics.erroring_NaNChecker!(simulation)
    add_callback!(simulation, sst_updater, IterationInterval(1))
    add_callback!(simulation, ProgressMessenger(model), TimeInterval(progress_interval))

    if write_output
        mkpath(output_dir)
        add_output_writers!(simulation; output_dir, output_prefix, profile_interval,
                            timeseries_interval, slice_interval, slice_height)
        if !isnothing(checkpoint_interval)
            simulation.output_writers[:checkpointer] =
                Checkpointer(model; schedule = TimeInterval(checkpoint_interval), dir = output_dir,
                             prefix = output_prefix * "_checkpoint", overwrite_files = true, cleanup = checkpoint_cleanup)
        end
    end

    label = string("TRACER DP-SCREAM-matched periodic case: IOP ", get(iop.attributes, "datastream", "?"),
                   " from ", start, " to ", stop, " UTC; transport ", large_scale_transport, " applied once, no subsidence, ",
                   isnothing(coriolis) ? "no Coriolis" : "Coriolis $(summary(coriolis))", ", wind nudging ",
                   isnothing(wind_nudging_timescale) ? "off (archived E3SM 1d551ea2b0)" : "τ = $wind_nudging_timescale s",
                   "; IOP H/LE + Tg, bulk drag Cᴰ = ", round(Cᴰ; sigdigits=3), " (z₀ = $roughness_length m); ",
                   "LES Smagorinsky at Δx = $Δx m on $(Nx * Δx / 1000) km (DP-SCREAM: SHOC at 3.33 km on 200 km, 128 levels, ",
                   "prognostic droplet number) — intentional physics/resolution differences, not a DP-SCREAM reproduction")

    config = (; protocol = "tracer_dp_scream", label, arch = string(typeof(arch)), FT = string(FT),
                Nx, Ny, Nz, Lx, Ly, Δx, Δy, z_top = z_faces[end], day0, epoch = string(epoch), latitude, longitude,
                start = string(start), stop = string(stop), window = inputs.window, forcing_records = length(inputs.lsf),
                large_scale_transport, hold_forcing, moisture_basis = string(moisture_basis),
                iop_datastream = get(iop.attributes, "datastream", ""), iop_doi = get(iop.attributes, "doi", ""),
                dp_scream_git_version = "1d551ea2b0", scmlib_revision = "2dc3f1073a5d03b5f32617e79d68fb20a63cfb43",
                zenodo_record = "10.5281/zenodo.15271730",
                coriolis = isnothing(coriolis) ? "nothing" : summary(coriolis),
                wind_nudging_timescale = something(wind_nudging_timescale, 0),
                microphysics_record(microphysics; reference_density = ρᵣ[1])...,
                initialization = string(initialization),
                radiation = radiation ? "RRTMGP all-sky LW+SW" : "none", radiation_interval,
                surface_albedo = albedo, surface_albedo_source = isnothing(surface_albedo) ? "IOP srfswup/srfswdn" : "caller",
                surface_emissivity, liquid_effective_radius, ice_effective_radius,
                background_CO₂ = background_atmosphere.CO₂, background_CH₄ = background_atmosphere.CH₄,
                background_N₂O = background_atmosphere.N₂O, background_O₃ = string(background_atmosphere.O₃),
                surface = "IOP sensible and latent heat fluxes, constant-coefficient bulk drag, IOP Tg",
                drag_coefficient = Cᴰ, roughness_length, gustiness,
                sponge = isnothing(sponge) ? "nothing" : summary(sponge),
                closure = isnothing(closure) ? "nothing" : summary(closure), advection_order,
                stop_time, Δt, max_Δt, cfl, perturbation = string(perturbation),
                write_output, output_dir = abspath(output_dir), output_prefix, profile_interval, timeseries_interval,
                slice_interval, slice_height, progress_interval,
                checkpoint_interval = something(checkpoint_interval, 0), checkpoint_cleanup)

    return (; simulation, model, grid, config, inputs = (; iop = iop_path),
              iop, start, stop, window = inputs.window, iop_inputs = inputs,
              soundings, sounding, lsf, sfc, profiles, columns, forcing_profiles, forcing,
              surface_series, surface_temperature = Tₛ)
end

"""
    tracer_dp_scream_settings(iop_path; start, stop)

Inspect the forcing window without allocating a model: records, IOP-implied albedo, and
the initial surface state.
"""
function tracer_dp_scream_settings(iop_path = joinpath(TRACER_DP_SCREAM_DATA, "TRACER_iopfile_4scam.nc");
                                   start = DateTime(2022, 8, 1), stop = DateTime(2022, 8, 15))
    iop = read_iop_forcing(iop_path)
    i0 = iop_index(iop, start); i1 = iop_index(iop, stop)
    return (; iop, start, stop, window = (i0, i1), records = i1 - i0 + 1,
              duration_seconds = Dates.value(stop - start) / 1000,
              surface_albedo = iop_surface_albedo(iop; start, stop),
              initial_surface_pressure = iop.surface.Ps[i0], initial_ground_temperature = iop.surface.Tg[i0],
              latitude = iop.latitude, longitude = iop.longitude,
              datastream = get(iop.attributes, "datastream", ""), doi = get(iop.attributes, "doi", ""))
end
