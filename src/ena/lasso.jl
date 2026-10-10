#####
##### ENA LASSO: an official LASSO-ENA member from its staged SAM `samin` bundle.
#####
##### Protocol (lasso-ena.svcs.arm.gov modeling_methodology + the lasso_ena_noice SAM):
#####   anelastic, doubly periodic 25.6 km × 25.6 km at 100 m, 260 levels from the bundle's
#####   grd, Coriolis from the namelist fcor, time-height large-scale forcing (tls, qls, wls,
#####   uls/vls, ug/vg), domain-mean wind nudging with τ = tauls, upper-boundary relaxation and
#####   sponge per the namelist, bulk surface fluxes from the SST (flxsst), RRTMG(P) LW+SW every
#####   nrad × dt, P3 with the member's aerosol and the SBM diagCCN reservoir rule.
##### Every namelist setting the constructor cannot represent is rejected by
##### `validate_lasso_bundle` (lasso_bundle.jl) before anything is allocated.
#####

"""
The member adopted after the authenticated ARM listing of 7 October 2026: the plan's
candidate `20170718era5s1n0d25x100_sbmwrm-aer2-flxsst` does not exist (the `s1n0` members
were run with Morrison microphysics only); the closest obtainable spectral-bin warm-cloud
SST-flux member is this one. Changing it is an explicit decision, never a default.
"""
const DEFAULT_LASSO_MEMBER = "20170718era5d25x100_sbmwrm-aer2-flxsst"

"""
    lasso_bundle_directory(member)

Where the staged bundle of `member` is expected: `\$ENA_LASSO_BUNDLE` when set, otherwise
`data/lasso/<member>` in the package checkout (populated by
`data_wrangling/stage_lasso_bundle.jl`).
"""
function lasso_bundle_directory(member)
    id = member isa LassoMember ? member.id : parse_lasso_member(member).id
    return get(ENV, "ENA_LASSO_BUNDLE", package_path("data", "lasso", id))
end

"""
    lasso_bundle_available(directory)

`true` when `snd`, `lsf`, `sfc`, `prm` and `grd` are all present in `directory`.
"""
lasso_bundle_available(directory) = all(isfile(joinpath(directory, f)) for f in LASSO_BUNDLE_FILES)

"""
    lasso_bundle_missing_message(member, directory)

The actionable message printed when a member's bundle is not staged: which files are
missing, the exact ARM archive name, the DOI and the staging command. It never suggests
substituting other inputs.
"""
function lasso_bundle_missing_message(member, directory)
    m = member isa LassoMember ? member : parse_lasso_member(member)
    missing_files = [f for f in LASSO_BUNDLE_FILES if !isfile(joinpath(directory, f))]
    archive = lasso_samin_filename(m)
    return string("LASSO-ENA bundle for member $(m.id) is not staged: $directory lacks ",
                  join(missing_files, ", "), ".\n",
                  "  Obtain the ARM archive $archive (DOI $LASSO_ENA_DOI; it must be staged online by ARM — ",
                  "the ARM Live `saveData` route returned HTTP 404 for it on 2026-10-07 although the listing shows it), then run\n",
                  "    julia data_wrangling/stage_lasso_bundle.jl data/lasso/archives/$archive\n",
                  "  (set ENA_LASSO_BUNDLE to use another directory). The Covert public inputs are a different protocol and are never substituted.")
end

"""
    lasso_documented_dimensions(member)

The documented LASSO-ENA domain for a member's `d<width>x<spacing>` token (modeling
methodology, domain table): `d25x100` is 25.6 km with 256 columns at 100 m and 260 levels
to 8 km. Returns `(Nx, Ny, Nz)`, used only when the caller asks for the documented domain
explicitly; the bundle's `grd` must agree.
"""
function lasso_documented_dimensions(member)
    m = member isa LassoMember ? member : parse_lasso_member(member)
    m.domain == "d25x100" && return (256, 256, 260)
    throw(ArgumentError("no documented dimensions for domain token $(m.domain); pass dimensions=(Nx, Ny, Nz) explicitly"))
end

"""
    ena_lasso(; member = DEFAULT_LASSO_MEMBER, epoch::DateTime, bundle_dir = lasso_bundle_directory(member),
                dimensions = :documented, arch = GPU(), microphysics = nothing, kwargs...)

Build the official LASSO-ENA member `member` from its staged `samin` bundle and return
`(; simulation, model, grid, config, inputs, ...)` without advancing it. Precision follows
`Oceananigans.defaults.FloatType`.

The bundle (`snd`, `lsf`, `sfc`, `prm`, `grd`) is validated by [`validate_lasso_bundle`](@ref):
every namelist setting the constructor cannot represent is rejected with the reason. A missing
bundle raises an `ArgumentError` with `lasso_bundle_missing_message` before the
architecture is touched, so the staging instructions also appear on machines without a GPU.

- `epoch::DateTime` is the UTC start; it must agree with the namelist `day0` and the member
  date (checked). The year is never inferred.
- `dimensions = :documented` uses `lasso_documented_dimensions` for the member's domain
  token; pass `(Nx, Ny, Nz)` to be explicit, or `nothing` when the namelist carries them.
- `microphysics = nothing` builds the protocol P3: the member's LASSO aerosol
  ([`lasso_aerosol`](@ref)) capped at the SBM `ss_max = 0.3 %`, converted with the
  first-level reference density, and a 75 cm⁻³ initial droplet number. Pass any Breeze
  microphysics object instead for a sensitivity experiment.
- `wind_nudging = true` nudges the domain-mean winds to `uls`, `vls` (namelist `donudging_uv`);
  `false` turns it off for a sensitivity experiment.
- `diagnostic_ccn = nothing` applies the SBM diagCCN reservoir rule (`DiagnosticCCNProjection`)
  whenever the microphysics carries a prognostic aerosol reservoir.

Keywords whose default is `nothing` take the bundle's value (any value the caller replaces is
listed in `config.overrides` and appended to the label): `Nx`, `Ny`, `Lx`, `Ly` (`dx`, `dy`),
`z_faces` (`grd`), `stop_time` (`nstop × dt`), `Δt` (`dt`, fixed: `max_Δt = Δt`),
`radiation_interval` (`nrad × dt`), `wind_nudging_timescale` (`tauls`),
`liquid_effective_radius` (14 μm, SAM's RRTMG ocean value, unless `compute_reffc`) and
`perturbation` (the namelist `perturb_type`). The others are the protocol's:
`surface_flux_law = :sam_oceflx`, `surface_emissivity = 0.95` (SAM RRTM), `surface_albedo = 0.07`,
`ice_effective_radius = 30e-6`, `background_atmosphere`, `radiation = true` (RRTMGP LW+SW;
`false` only for software tests), `closure = SmagorinskyLilly()`, `advection_order = 5`,
`cfl = 0.7`, `initialization = :condensate_free`, and the output keywords of [`ena_covert`](@ref).
"""
function ena_lasso(; member = DEFAULT_LASSO_MEMBER,
                     epoch,
                     bundle_dir = nothing,
                     dimensions = :documented,
                     arch = nothing,
                     microphysics = nothing,
                     diagnostic_ccn = nothing,
                     Nx = nothing,
                     Ny = nothing,
                     Lx = nothing,
                     Ly = nothing,
                     z_faces = nothing,
                     stop_time = nothing,
                     Δt = nothing,
                     max_Δt = Δt,
                     cfl = 0.7,
                     radiation = true,
                     radiation_interval = nothing,
                     liquid_effective_radius = nothing,
                     ice_effective_radius = 30e-6,
                     surface_albedo = 0.07,
                     surface_emissivity = 0.95,
                     background_atmosphere = BackgroundAtmosphere(CO₂ = 405e-6, CH₄ = 1.85e-6, N₂O = 330e-9),
                     surface_flux_law = :sam_oceflx,
                     wind_nudging = true,
                     wind_nudging_timescale = nothing,
                     closure = SmagorinskyLilly(),
                     advection_order = 5,
                     perturbation = nothing,
                     initialization = :condensate_free,
                     initial_droplet_number = nothing,
                     write_output = true,
                     output_dir = nothing,
                     output_prefix = "ena_lasso",
                     profile_interval = 1hour,
                     timeseries_interval = 60,
                     slice_interval = 10minutes,
                     slice_height = 900,
                     progress_interval = 10minutes,
                     checkpoint_interval = nothing,
                     checkpoint_cleanup = true)

    #####
    ##### Inputs: the validated samin bundle
    #####

    isnothing(member) && throw(ArgumentError("ena_lasso needs member = \"<run ID>\" naming the samin bundle; the identity is never inferred"))
    member = member isa LassoMember ? member : parse_lasso_member(member)
    directory = something(bundle_dir, lasso_bundle_directory(member))
    lasso_bundle_available(directory) || throw(ArgumentError(lasso_bundle_missing_message(member, directory)))
    epoch isa DateTime || throw(ArgumentError("ena_lasso needs epoch = DateTime(...) in UTC from the ARM bundle"))
    isnothing(arch) && (arch = GPU())
    requested_dimensions = dimensions === :documented ? lasso_documented_dimensions(member) : dimensions
    bundle = validate_lasso_bundle(directory; member, dimensions = requested_dimensions, epoch)

    paths = (; snd = joinpath(directory, "snd"), lsf = joinpath(directory, "lsf"), sfc = joinpath(directory, "sfc"),
               prm = joinpath(directory, "prm"), grd = joinpath(directory, "grd"))
    soundings = read_sam_sounding(paths.snd)
    lsf = read_sam_large_scale_forcing(paths.lsf)
    sfc = read_sam_surface_forcing(paths.sfc)
    staged = joinpath(directory, "bundle.toml")
    staging = isfile(staged) ? TOML.parsefile(staged) :
              "no bundle.toml: the directory was not staged with data_wrangling/stage_lasso_bundle.jl (archive checksum/download provenance unavailable)"

    overrides = [name for (name, value) in pairs((; Nx, Ny, Lx, Ly, z_faces, stop_time, Δt, radiation_interval,
                                                  wind_nudging_timescale, liquid_effective_radius, perturbation, microphysics))
                 if !isnothing(value)]
    radiation || push!(overrides, :radiation)
    wind_nudging || push!(overrides, :wind_nudging)
    (; domain, switches, location) = bundle
    timing = bundle.time
    Nx = something(Nx, domain.Nx)
    Ny = something(Ny, domain.Ny)
    Lx = something(Lx, domain.Nx * domain.dx)
    Ly = something(Ly, domain.Ny * domain.dy)
    z_faces = something(z_faces, bundle.grid.faces)
    Δt = something(Δt, timing.dt)
    max_Δt = something(max_Δt, Δt)
    stop_time = something(stop_time, Float64(timing.stop_time))
    radiation_interval = something(radiation_interval, Float64(timing.radiation_interval))
    wind_nudging_timescale = something(wind_nudging_timescale, timing.tauls)
    # cam_rad_parameterizations: SAM's RRTMG uses rliqocean = 14 μm unless compute_reffc
    liquid_effective_radius = something(liquid_effective_radius, switches.compute_reffc ? 10e-6 : 14e-6)
    perturbation = something(perturbation, InitialPerturbation(; sam_perturb_type = switches.perturb_type))
    day0 = timing.day0
    (; latitude, longitude) = location
    moisture_basis = :mixing_ratio            # SAM files carry mixing ratios per kg of dry air
    output_dir = something(output_dir, joinpath(pwd(), "output", "ena_lasso", member.id))

    sounding = initial_sounding(soundings, day0)
    profiles = SoundingProfiles(sounding; moisture_basis)

    #####
    ##### Grid, reference state, dynamics
    #####

    Nz = length(z_faces) - 1
    grid = RectilinearGrid(arch; size = (Nx, Ny, Nz), x = (0, Lx), y = (0, Ly), z = z_faces,
                           halo = (5, 5, 5), topology = (Periodic, Periodic, Bounded))
    FT = eltype(grid)
    constants = ThermodynamicConstants()
    reference_state = sounding_reference_state(grid, profiles, constants)
    dynamics = AnelasticDynamics(reference_state)
    # SAM uses the namelist fcor directly when it is given (setgrid.f90: only fcor = -999 is
    # replaced by 4π/86400 sin φ); the LASSO bundles carry it.
    coriolis = isnothing(bundle.coriolis_parameter) ? FPlane(; latitude) : FPlane(; f = bundle.coriolis_parameter)

    z_centers = Array(znodes(grid, Center()))
    ρᵣ = Array(interior(reference_state.density, 1, 1, :))
    pᵣ = Array(interior(reference_state.pressure, 1, 1, :))

    #####
    ##### Microphysics and advection
    #####

    protocol_microphysics = isnothing(microphysics)
    if protocol_microphysics
        aerosol = lasso_aerosol(; setting = member.aerosol, reference_density = ρᵣ[1],
                                maximum_supersaturation = 0.003)       # HUJI-SBM ss_max (microphysics.f90)
        microphysics = P3Microphysics(; cloud = CloudDroplets(; number_concentration = 75e6), aerosol)
    end
    check_precision(microphysics, FT)
    diagnostic_ccn = something(diagnostic_ccn, has_aerosol_reservoir(microphysics))
    diagnostic_ccn && !has_aerosol_reservoir(microphysics) &&
        throw(ArgumentError("diagnostic_ccn = true needs P3 with a prognostic aerosol reservoir"))
    moisture_name = Breeze.AtmosphereModels.moisture_specific_name(microphysics)
    momentum_advection = WENO(order = advection_order)
    scalar_advection = scalar_advection_schemes(advection_order, microphysics, moisture_name; energy_name = :ρθ)

    #####
    ##### Large-scale forcing (SAM forcing.f90, nudging.f90, subsidence.f90, upperbound.f90, damping.f90)
    #####

    # The namelist translation frame (ug, vg) is not applied: the bulk fluxes need the
    # ground-relative wind, so winds and targets stay ground-relative (the same physical state).
    forcing_profiles = LargeScaleForcingProfiles(grid, lsf, z_centers, pᵣ; day0)
    geostrophic = time_varying_geostrophic_forcings(forcing_profiles.ug, forcing_profiles.vg)
    nudging_u = wind_nudging ? MeanProfileNudging(forcing_profiles.uls; timescale = wind_nudging_timescale) : nothing
    nudging_v = wind_nudging ? MeanProfileNudging(forcing_profiles.vls; timescale = wind_nudging_timescale) : nothing
    thermodynamic = large_scale_thermodynamic_forcings(forcing_profiles.tls, forcing_profiles.qls;
                                                       microphysics, thermodynamic_constants = constants,
                                                       moisture_name, moisture_basis)
    subsidence = LargeScaleVerticalAdvection(forcing_profiles.wls)
    upper = if switches.doupperbound
        targets = SoundingTargetProfiles(grid, soundings, z_centers, pᵣ; day0, moisture_basis)
        upper_boundary_relaxation_forcings(targets.T, targets.q; microphysics,
                                           thermodynamic_constants = constants, moisture_name)
    else
        NamedTuple{(:s, moisture_name)}((nothing, nothing))
    end
    sponge = switches.dodamping ? SAMSponge() : nothing

    forcing = Dict{Symbol, Tuple}(
        :u => (geostrophic.u, nudging_u, subsidence, sponge),
        :v => (geostrophic.v, nudging_v, subsidence, sponge),
        :w => (sponge,),
        :E => (thermodynamic.E, upper.E),            # tls and the top relaxation, as energy tendencies
        :θ => (subsidence,))
    for name in specific_prognostic_names(microphysics)
        forcing[name] = (subsidence,)                # subsidence.f90 advects every microphysical field
    end
    forcing[moisture_name] = (thermodynamic[moisture_name], subsidence, upper[moisture_name])
    forcing = compact_forcing(forcing)

    #####
    ##### Surface: bulk fluxes from the sfc SST series (flxsst)
    #####

    surface_series = surface_time_series(grid, sfc, day0)
    Tₛ = Field{Center, Center, Nothing}(grid)
    set!(Tₛ, FT(sfc.sst[1]))
    sst_updater = SeaSurfaceTemperatureUpdater(Tₛ, surface_series.times, FT.(sfc.sst))
    boundary_conditions, surface_record = bulk_surface_flux_boundary_conditions(grid, Tₛ; moisture_name,
                                                                               law = surface_flux_law)

    #####
    ##### Radiation: RRTMGP LW+SW every nrad × dt
    #####

    radiation_model = if radiation
        RadiativeTransferModel(grid, AllSkyOptics(), constants;
                               surface_temperature = Tₛ,
                               surface_albedo, surface_emissivity, background_atmosphere,
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
    if diagnostic_ccn
        n_initial = sum(mode.number_mixing_ratio for mode in microphysics.aerosol.modes)
        add_callback!(simulation, DiagnosticCCNProjection(model, n_initial), IterationInterval(1))
    end
    add_callback!(simulation, ProgressMessenger(model), TimeInterval(progress_interval))

    if write_output
        mkpath(output_dir)
        add_output_writers!(simulation; output_dir, output_prefix, profile_interval,
                            timeseries_interval, slice_interval, slice_height,
                            energy_budget_series = true)
        if !isnothing(checkpoint_interval)
            simulation.output_writers[:checkpointer] =
                Checkpointer(model; schedule = TimeInterval(checkpoint_interval), dir = output_dir,
                             prefix = output_prefix * "_checkpoint", overwrite_files = true, cleanup = checkpoint_cleanup)
        end
    end

    label = isempty(overrides) ? bundle.label : string(bundle.label, " [overrides: ", join(overrides, ", "), "]")
    config = (; protocol = "lasso_ena_official", label, overrides = string.(overrides), member = member.id,
                dimension_request = string(dimensions), arch = string(typeof(arch)), FT = string(FT),
                Nx, Ny, Nz, Lx, Ly, z_top = z_faces[end], day0, epoch = string(epoch), latitude, longitude,
                coriolis_parameter = Float64(coriolis.f),
                coriolis_parameter_source = isnothing(bundle.coriolis_parameter) ? "2Ω sin(latitude)" : "namelist fcor",
                microphysics_record(microphysics; reference_density = ρᵣ[1])...,
                protocol_microphysics,
                aerosol_supersaturation_cap = protocol_microphysics ? 0.003 : "set by the caller's microphysics",
                diagnostic_ccn, initialization = string(initialization),
                initial_droplet_number = something(initial_droplet_number, 0),
                radiation = radiation ? "RRTMGP all-sky LW+SW" : "none (software test)", radiation_interval,
                liquid_effective_radius, ice_effective_radius, surface_albedo, surface_emissivity,
                background_CO₂ = background_atmosphere.CO₂, background_CH₄ = background_atmosphere.CH₄,
                background_N₂O = background_atmosphere.N₂O, background_O₃ = string(background_atmosphere.O₃),
                surface = "bulk fluxes from the SST series (flxsst)", surface_record...,
                wind_nudging, wind_nudging_timescale, upper_boundary_relaxation = switches.doupperbound,
                sponge = isnothing(sponge) ? "nothing" : summary(sponge),
                closure = isnothing(closure) ? "nothing" : summary(closure), advection_order,
                stop_time, Δt, max_Δt, cfl,
                time_stepping = max_Δt == Δt ? "fixed Δt = $Δt s" :
                                "adaptive (initial Δt = $Δt s, cfl = $cfl, max_Δt = $max_Δt s); SAM used fixed dt",
                perturbation = string(perturbation),
                write_output, output_dir = abspath(output_dir), output_prefix, profile_interval, timeseries_interval,
                slice_interval, slice_height, progress_interval,
                checkpoint_interval = something(checkpoint_interval, 0), checkpoint_cleanup,
                bundle = lasso_bundle_record(bundle), staging)

    return (; simulation, model, grid, config, inputs = paths, bundle,
              soundings, sounding, lsf, sfc, profiles, columns, forcing_profiles, forcing,
              surface_series, surface_temperature = Tₛ)
end
