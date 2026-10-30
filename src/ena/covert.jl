#####
##### ENA Covert: the public Covert et al. (2022) development benchmark, 18 July 2017 06 UTC.
#####
##### From the public files of the Covert bin-paper repository (snd, lsf, sfc, prm):
#####   caseid 256x256x192 at dx = dy = 35 m (8.96 km), day0 = 199.25, nstop × dt = 6 h,
#####   prescribed H, LE and τ (SFC_FLX_FXD, SFC_TAU_FXD), SAM rad_simple longwave only,
#####   no wind nudging, upper-boundary relaxation (doupperbound) and sponge (dodamping),
#####   SAM's translating frame (namelist ug, vg).
##### Not the 864²×192, 30.24-km, 06-15 UTC configuration of the published paper, and not
##### an official LASSO-ENA reproduction (see `ena_lasso`).
#####

"""
    covert_public_bin_vertical_faces(; Nz=192, top=20000, uniform_top=1500,
                                       Δz=uniform_top / (Nz - round(Int, 7Nz / 32)))

A labelled **reconstruction** of the 192-level, 20-km-top vertical grid of the Covert et al.
(2022) SAM runs (the public repository advertises but does not ship its `grd`/`domain.f90`):
uniform `Δz` to `uniform_top`, then a constant-ratio stretch to `top`.

With the default `Nz = 192` this is 150 cells of 10 m to 1.5 km and 42 stretched cells to 20 km.
Other `Nz` keep the design (the boundary layer and inversion uniformly resolved, the same 7/32
of the levels stretched above `uniform_top`): `Nz = 128` gives 100 cells of 15 m and 28 stretched
cells, `Nz = 256` gives 200 cells of 7.5 m and 56 stretched cells.
"""
function covert_public_bin_vertical_faces(; Nz=192, top=20000.0, uniform_top=1500.0,
                                            Δz=uniform_top / (Nz - round(Int, 7Nz / 32)))
    faces = uniform_then_stretched_faces(; Nz, top, Δz, uniform_top)
    faces[end] = top
    length(faces) == Nz + 1 || error("vertical grid construction produced $(length(faces) - 1) cells, expected $Nz")
    return faces
end

"""
    covert_inversion_refined_vertical_faces(; Δz=10.0, Δz_fine=5.0, fine_bottom=800.0, fine_top=1400.0,
                                             uniform_top=1500.0, top=20000.0, N_stretch=42)

The [`covert_public_bin_vertical_faces`](@ref) reconstruction with the inversion layer
refined to `Δz_fine` between `fine_bottom` and `fine_top` (5 m over 800–1400 m, the spacing
Covert et al. (2022) use "near the surface and inversion layer"), `Δz` elsewhere below
`uniform_top`, and the same `N_stretch` geometrically stretched cells to `top` as the
reference grid. With the defaults this gives 80 + 120 + 10 + 42 = 252 cells (the reference
has 192). A labelled sensitivity grid for the cloud-top/inversion bias test, not the
paper's grid (which is 5 m near the surface as well and 192 levels in total).
"""
function covert_inversion_refined_vertical_faces(; Δz=10.0, Δz_fine=5.0, fine_bottom=800.0, fine_top=1400.0,
                                                   uniform_top=1500.0, top=20000.0, N_stretch=42)
    0 < fine_bottom < fine_top ≤ uniform_top < top || throw(ArgumentError("need 0 < fine_bottom < fine_top ≤ uniform_top < top"))
    faces = collect(0.0:Δz:fine_bottom)
    append!(faces, (fine_bottom + Δz_fine):Δz_fine:fine_top)
    fine_top < uniform_top && append!(faces, (fine_top + Δz):Δz:uniform_top)
    r = geometric_growth_ratio((top - uniform_top) / Δz, N_stretch)
    for n in 1:N_stretch
        push!(faces, faces[end] + Δz * r^n)
    end
    faces[end] = top
    issorted(faces) && all(>(0), diff(faces)) || error("refined vertical grid is not strictly increasing")
    return faces
end

"""
    ena_vertical_faces(name)

Named vertical grids of the ENA cases: `:covert` ([`covert_public_bin_vertical_faces`](@ref)),
`:covert_inversion_5m` ([`covert_inversion_refined_vertical_faces`](@ref)), `:covert_r128_inversion_5m`
(the `r128` levels, 15 m, with 5 m over 990–1395 m: 182 cells), `:covert_r256_inversion_5m` (the `r256`
levels, 10 m, with 5 m over 990–1400 m: 233 cells) and `:lasso`
([`lasso_ena_vertical_faces`](@ref)). The name is what run scripts record in provenance.
"""
function ena_vertical_faces(name)
    name = Symbol(name)
    name === :covert && return covert_public_bin_vertical_faces()
    name === :covert_inversion_5m && return covert_inversion_refined_vertical_faces()
    name === :covert_r128_inversion_5m &&
        return covert_inversion_refined_vertical_faces(; Δz = 15.0, fine_bottom = 990.0, fine_top = 1395.0, N_stretch = 28)
    name === :covert_r256_inversion_5m &&
        return covert_inversion_refined_vertical_faces(; Δz = 10.0, fine_bottom = 990.0, fine_top = 1400.0, N_stretch = 42)
    name === :lasso && return lasso_ena_vertical_faces()
    throw(ArgumentError("unknown vertical grid $name (covert, covert_inversion_5m, covert_r128_inversion_5m, covert_r256_inversion_5m, lasso)"))
end

"""
    COVERT_RESOLUTIONS

The resolution ladder of the Covert case at the fixed 8.96 km domain: `r128` (128² × 128, 70 m,
Δz = 15 m), `r256` (256² × 192, 35 m, Δz = 10 m: the public `caseid` grid) and `r512` (512² × 256,
17.5 m, Δz = 7.5 m), each with its time step (the protocol's 0.5 s at 35 m, scaled with Δx).
`Lx` and `Ly` are given explicitly: the default `Lx = Nx dx` would keep the namelist's 35 m
spacing and change the domain instead.
Pass `COVERT_RESOLUTIONS.r512...` to [`ena_covert`](@ref).
"""
const COVERT_RESOLUTIONS = (r128 = (Nx = 128, Ny = 128, Nz = 128, Lx = 8960.0, Ly = 8960.0, Δt = 1.0),
                            r256 = (Nx = 256, Ny = 256, Nz = 192, Lx = 8960.0, Ly = 8960.0, Δt = 0.5),
                            r512 = (Nx = 512, Ny = 512, Nz = 256, Lx = 8960.0, Ly = 8960.0, Δt = 0.25))

# The namelist switches of the Covert protocol: the constructor implements exactly these.
const COVERT_NAMELIST = (sfc_flx_fxd = true, sfc_tau_fxd = true, doradsimple = true, dolongwave = true,
                         doshortwave = false, donudging_uv = false, dolargescale = true, dosfcforcing = true,
                         doupperbound = true, dodamping = true, docoriolis = true)

function check_covert_namelist(namelist, path)
    required = ("caseid", "dx", "dy", "dt", "nstop", "day0", "latitude0", "longitude0", "ug", "vg",
                string.(keys(COVERT_NAMELIST))...)
    missing_keys = [k for k in required if !haskey(namelist, k)]
    isempty(missing_keys) ||
        throw(ArgumentError("$path is not a Covert-protocol namelist: it lacks $(missing_keys)"))
    for (key, value) in pairs(COVERT_NAMELIST)
        namelist[string(key)] == value ||
            throw(ArgumentError("$path is not a Covert-protocol namelist: $key = $(namelist[string(key)]), the protocol requires $value"))
    end
    return nothing
end

"""
    ena_covert(; arch = GPU(),
                 data_dir = <package>/data/covert2022_bin,
                 microphysics = P3Microphysics(; cloud = CloudDroplets(; number_concentration = 75e6)),
                 kwargs...)

Build the public Covert ENA development case for 18 July 2017, starting at 06 UTC, and return
`(; simulation, model, grid, config, inputs, ...)` without advancing it. Call
`run!(case.simulation)` to run it. Precision follows `Oceananigans.defaults.FloatType`.

The inputs are the SAM `snd`, `lsf`, `sfc` and `prm` files in `data_dir` (fetch them with
`julia data_wrangling/fetch_covert_inputs.jl`). The namelist must carry the Covert protocol
switches (prescribed fluxes and stress, `rad_simple` longwave, no wind nudging, upper-boundary
relaxation and sponge); anything else is rejected.

Keyword arguments (the namelist value is used where the default is `nothing`; any such value the
caller replaces is listed in `config.overrides` and appended to the label):

- `microphysics`: any Breeze microphysics. The default is P3 with the observed 75 cm⁻³
  droplet number. For an aerosol-coupled member pass P3 with `aerosol = covert_aerosol(...)`
  or `lasso_aerosol(...)`, converting with [`first_level_reference_density`](@ref).
- `diagnostic_ccn`: apply the SBM-style `DiagnosticCCNProjection` to the aerosol reservoir;
  on by default when `microphysics` carries one.
- `Nx`, `Ny`, `Nz` (namelist `caseid`), `Lx = Nx dx`, `Ly = Ny dy`,
  `z_faces` (`covert_public_bin_vertical_faces(; Nz)`; give `Nz` or `z_faces`, not both)
- `stop_time` (`nstop × dt`), `Δt` (`dt`), `max_Δt = Δt` (fixed step), `cfl = 0.7`
- `translation_velocity` (namelist `ug`, `vg`): SAM's translating frame
- `sponge = SAMSponge()`, `closure = SmagorinskyLilly()`, `advection_order = 5` (WENO order of momentum
  and, by default, `scalar_advection_order` of the scalars)
- `reference_extension = nothing`: `(z₁, z₂)` extends the anelastic reference profile's `z₁` values up
  to `z₂` (see `sounding_reference_state`; a sensitivity test, recorded as an override)
- `radiation = :simple_longwave` (the protocol's SAM `rad_simple`, every step) or `:rrtmgp_longwave`
  (a sensitivity, recorded as an override: RRTMGP all-sky optics with the solar constant set to 0 so
  only longwave acts, every 60 s, over the sfc SST with emissivity 0.95, 10 μm / 30 μm liquid/ice
  effective radii and 405 ppm CO₂; the column ends at the LES top)
- `perturbation = InitialPerturbation()` (`setperturb.f90` case 5)
- `initialization = :condensate_free` or `:equilibrium` (P3 only), `initial_droplet_number`
- output: `write_output`, `output_dir`, `output_prefix`, `profile_interval = 1hour`,
  `timeseries_interval = 60`, `slice_interval = 1minute` (the movie standard: ≤ 15 simulated minutes per
  second of video at 10–15 fps), `slice_height = 900` (the `qᶜˡ_xy` and `w_xy` level),
  `extra_slice_heights = (; w = (500, 750), qᶜˡ = (1100,))` (mid-boundary layer, near cloud base and
  near cloud top), `y_averaged_sections = true`, `section_fields = model -> (;)` (a function of the model
  returning further 3D fields or operations, such as process rates, written as xz sections, y-averaged
  sections and the xy slices named in `extra_slice_heights`; see `add_output_writers!`),
  `process_rates` (on for P3: the `process_rate_outputs` in `<prefix>_process_rates.jld2` every
  `process_rate_interval = slice_interval`, xy slices at `process_rate_heights`),
  `overwrite_output = true` (`false` appends to existing output files, for a run picked up from a checkpoint),
  `progress_interval = 10minutes`, `checkpoint_interval = nothing`, `checkpoint_cleanup = true`
"""
function ena_covert(; arch = GPU(),
                      data_dir = package_path("data", "covert2022_bin"),
                      microphysics = P3Microphysics(; cloud = CloudDroplets(; number_concentration = 75e6)),
                      diagnostic_ccn = has_aerosol_reservoir(microphysics),
                      Nx = nothing,
                      Ny = nothing,
                      Nz = nothing,
                      Lx = nothing,
                      Ly = nothing,
                      z_faces = nothing,
                      stop_time = nothing,
                      Δt = nothing,
                      max_Δt = Δt,
                      cfl = 0.7,
                      translation_velocity = nothing,
                      sponge = SAMSponge(),
                      closure = SmagorinskyLilly(),
                      radiation = :simple_longwave,
                      advection_order = 5,
                      scalar_advection_order = advection_order,
                      reference_extension = nothing,
                      perturbation = InitialPerturbation(),
                      initialization = :condensate_free,
                      initial_droplet_number = nothing,
                      label = "Covert-public-bin development benchmark (not an official LASSO-ENA reproduction)",
                      write_output = true,
                      output_dir = joinpath(pwd(), "output", "ena_covert"),
                      output_prefix = "ena_covert",
                      profile_interval = 1hour,
                      timeseries_interval = 60,
                      slice_interval = 1minute,
                      slice_height = 900,
                      extra_slice_heights = (; w = (500, 750), qᶜˡ = (1100,)),
                      y_averaged_sections = true,
                      section_fields = model -> (;),
                      process_rates = microphysics isa PredictedParticlePropertiesMicrophysics,
                      process_rate_interval = slice_interval,
                      process_rate_heights = PROCESS_RATE_SLICE_HEIGHTS,
                      overwrite_output = true,
                      progress_interval = 10minutes,
                      checkpoint_interval = nothing,
                      checkpoint_cleanup = true)

    #####
    ##### Inputs: SAM snd, lsf, sfc and the Covert namelist
    #####

    paths = (; snd = joinpath(data_dir, "snd"), lsf = joinpath(data_dir, "lsf"),
               sfc = joinpath(data_dir, "sfc"), prm = joinpath(data_dir, "prm"))
    missing_files = [name for (name, path) in pairs(paths) if !isfile(path)]
    isempty(missing_files) ||
        throw(ArgumentError("the Covert protocol needs $(join(missing_files, ", ")) in $data_dir (julia data_wrangling/fetch_covert_inputs.jl); no other inputs are substituted"))
    namelist = read_sam_namelist(paths.prm)
    check_covert_namelist(namelist, paths.prm)
    soundings = read_sam_sounding(paths.snd)
    lsf = read_sam_large_scale_forcing(paths.lsf)
    sfc = read_sam_surface_forcing(paths.sfc)

    overrides = [name for (name, value) in pairs((; Nx, Ny, Nz, Lx, Ly, z_faces, stop_time, Δt, translation_velocity))
                 if !isnothing(value)]
    radiation === :simple_longwave || push!(overrides, :radiation)
    isnothing(reference_extension) || push!(overrides, :reference_extension)
    nx, ny, nz = parse_caseid(namelist["caseid"])
    isnothing(Nz) || isnothing(z_faces) || throw(ArgumentError("give Nz or z_faces, not both"))
    z_faces = isnothing(z_faces) ? covert_public_bin_vertical_faces(; Nz = something(Nz, nz)) : z_faces
    Nx = something(Nx, nx)
    Ny = something(Ny, ny)
    Lx = something(Lx, Nx * namelist["dx"])
    Ly = something(Ly, Ny * namelist["dy"])
    Δt = something(Δt, Float64(namelist["dt"]))
    max_Δt = something(max_Δt, Δt)
    stop_time = something(stop_time, Float64(namelist["nstop"] * namelist["dt"]))
    translation_velocity = something(translation_velocity, (Float64(namelist["ug"]), Float64(namelist["vg"])))
    day0 = Float64(namelist["day0"])
    epoch = epoch_from_day_of_year(day0)
    latitude = Float64(namelist["latitude0"])
    longitude = Float64(namelist["longitude0"])
    moisture_basis = :mixing_ratio            # SAM files carry mixing ratios per kg of dry air

    # setdata.f90 interpolates the bracketing snd records to day0
    sounding = initial_sounding(soundings, day0)
    profiles = SoundingProfiles(sounding; moisture_basis)

    #####
    ##### Grid, reference state, dynamics
    #####

    Nz = length(z_faces) - 1
    # WENO orders above 7 need a sixth horizontal halo point
    halo = max(advection_order, scalar_advection_order) > 7 ? (6, 6, 5) : (5, 5, 5)
    grid = RectilinearGrid(arch; size = (Nx, Ny, Nz), x = (0, Lx), y = (0, Ly), z = z_faces,
                           halo, topology = (Periodic, Periodic, Bounded))
    FT = eltype(grid)
    constants = ThermodynamicConstants()
    reference_state = sounding_reference_state(grid, profiles, constants; extension = reference_extension)
    dynamics = AnelasticDynamics(reference_state)
    coriolis = FPlane(; latitude)

    z_centers = Array(znodes(grid, Center()))
    ρᵣ = Array(interior(reference_state.density, 1, 1, :))
    pᵣ = Array(interior(reference_state.pressure, 1, 1, :))

    #####
    ##### Microphysics and advection
    #####

    check_precision(microphysics, FT)
    diagnostic_ccn && !has_aerosol_reservoir(microphysics) &&
        throw(ArgumentError("diagnostic_ccn = true needs P3 with a prognostic aerosol reservoir"))
    moisture_name = Breeze.AtmosphereModels.moisture_specific_name(microphysics)
    momentum_advection = WENO(order = advection_order)
    scalar_advection = scalar_advection_schemes(scalar_advection_order, microphysics, moisture_name; energy_name = :ρθ)

    #####
    ##### Large-scale forcing (SAM forcing.f90, subsidence.f90, upperbound.f90, damping.f90)
    #####

    forcing_profiles = LargeScaleForcingProfiles(grid, lsf, z_centers, pᵣ; day0)

    # SAM solves for the winds relative to the translating frame (namelist ug, vg): the pair is
    # subtracted from the initial winds and the geostrophic profiles, and added back for the
    # surface stress. Vertical advection, the sponge and Coriolis are frame-invariant.
    uᶠ, vᶠ = FT.(translation_velocity)
    geostrophic = time_varying_geostrophic_forcings(shifted_profile_time_series(forcing_profiles.ug, uᶠ),
                                                    shifted_profile_time_series(forcing_profiles.vg, vᶠ))
    thermodynamic = large_scale_thermodynamic_forcings(forcing_profiles.tls, forcing_profiles.qls;
                                                       microphysics, thermodynamic_constants = constants,
                                                       moisture_name, moisture_basis)
    subsidence = LargeScaleVerticalAdvection(forcing_profiles.wls)
    targets = SoundingTargetProfiles(grid, soundings, z_centers, pᵣ; day0, moisture_basis)
    upper = upper_boundary_relaxation_forcings(targets.T, targets.q; microphysics,
                                               thermodynamic_constants = constants, moisture_name)

    forcing = Dict{Symbol, Tuple}(
        :u => (geostrophic.u, subsidence, sponge),
        :v => (geostrophic.v, subsidence, sponge),
        :w => (sponge,),
        :E => (thermodynamic.E, upper.E),            # tls and the top relaxation, as energy tendencies
        :θ => (subsidence,))
    for name in specific_prognostic_names(microphysics)
        forcing[name] = (subsidence,)                # subsidence.f90 advects every microphysical field
    end
    forcing[moisture_name] = (thermodynamic[moisture_name], subsidence, upper[moisture_name])
    forcing = compact_forcing(forcing)

    #####
    ##### Surface: prescribed H, LE and wind-aligned stress τ (SFC_FLX_FXD, SFC_TAU_FXD)
    #####

    boundary_conditions, stress = prescribed_surface_flux_boundary_conditions(grid, sfc, day0;
                                                                             thermodynamic_constants = constants,
                                                                             surface_density = ρᵣ[1], moisture_name,
                                                                             frame_velocity = (uᶠ, vᶠ))

    #####
    ##### Radiation: SAM rad_simple longwave, every step
    #####

    radiation in (:simple_longwave, :rrtmgp_longwave) ||
        throw(ArgumentError("radiation must be :simple_longwave (the protocol) or :rrtmgp_longwave, got $radiation"))
    radiation_scheme = radiation
    sst_updater = nothing
    radiation = if radiation_scheme === :simple_longwave
        SimpleLongwaveRadiation(grid; schedule = IterationInterval(1))
    else
        surface_series = surface_time_series(grid, sfc, day0)
        Tₛ = Field{Center, Center, Nothing}(grid)
        set!(Tₛ, FT(sfc.sst[1]))
        sst_updater = SeaSurfaceTemperatureUpdater(Tₛ, surface_series.times, FT.(sfc.sst))
        RadiativeTransferModel(grid, AllSkyOptics(), constants;
                               surface_temperature = Tₛ, surface_emissivity = 0.95, surface_albedo = 0.07,
                               solar_constant = 0,
                               background_atmosphere = BackgroundAtmosphere(CO₂ = 405e-6, CH₄ = 1.85e-6, N₂O = 330e-9),
                               solar_position = ApparentSolarPosition(; coordinate = (longitude, latitude), epoch),
                               schedule = TimeInterval(60),
                               liquid_effective_radius = ConstantRadiusParticles(10e-6),
                               ice_effective_radius = ConstantRadiusParticles(30e-6))
    end

    #####
    ##### Model and initial state
    #####

    model = AtmosphereModel(grid; formulation = :LiquidIcePotentialTemperature, dynamics, coriolis, closure,
                            microphysics, radiation, momentum_advection, scalar_advection, forcing,
                            boundary_conditions, thermodynamic_constants = constants)

    columns = initial_state_columns(profiles, z_centers, pᵣ; constants = ThermodynamicConstants(Float64))
    set_sounding_initial_state!(model, columns; perturbation, moisture_basis, initialization,
                                frame_velocity = (uᶠ, vᶠ), initial_droplet_number)

    #####
    ##### Simulation
    #####

    simulation = Simulation(model; Δt, stop_time)
    conjure_time_step_wizard!(simulation; cfl, max_Δt)
    Oceananigans.Diagnostics.erroring_NaNChecker!(simulation)
    stress_updater = prescribed_stress_updater(stress, model.velocities)
    stress_updater(simulation)                            # the stress of the initial wind
    add_callback!(simulation, stress_updater, IterationInterval(1))
    isnothing(sst_updater) || add_callback!(simulation, sst_updater, IterationInterval(1))
    if diagnostic_ccn
        n_initial = sum(mode.number_mixing_ratio for mode in microphysics.aerosol.modes)
        add_callback!(simulation, DiagnosticCCNProjection(model, n_initial), IterationInterval(1))
    end
    add_callback!(simulation, ProgressMessenger(model), TimeInterval(progress_interval))

    if write_output
        mkpath(output_dir)
        add_output_writers!(simulation; output_dir, output_prefix, profile_interval,
                            timeseries_interval, slice_interval, slice_height,
                            extra_slice_heights, y_averaged_sections,
                            section_fields = section_fields(model),
                            energy_budget_series = true, overwrite_files = overwrite_output)
        if process_rates
            simulation.output_writers[:process_rates] =
                JLD2Writer(model, process_rate_outputs(model; heights = process_rate_heights);
                           filename = joinpath(output_dir, output_prefix * "_process_rates.jld2"),
                           schedule = TimeInterval(process_rate_interval), overwrite_files = overwrite_output)
        end
        if !isnothing(checkpoint_interval)
            simulation.output_writers[:checkpointer] =
                Checkpointer(model; schedule = TimeInterval(checkpoint_interval), dir = output_dir,
                             prefix = output_prefix * "_checkpoint", overwrite_files = true, cleanup = checkpoint_cleanup)
        end
    end

    label = isempty(overrides) ? label : string(label, " [overrides: ", join(overrides, ", "), "]")
    config = (; protocol = "covert_public_bin", label, overrides = string.(overrides), arch = string(typeof(arch)), FT = string(FT),
                Nx, Ny, Nz, Lx, Ly, z_top = z_faces[end], day0, epoch = string(epoch), latitude, longitude,
                coriolis_parameter = Float64(coriolis.f),
                translation_velocity_u = Float64(uᶠ), translation_velocity_v = Float64(vᶠ),
                microphysics_record(microphysics; reference_density = ρᵣ[1])...,
                diagnostic_ccn, initialization = string(initialization),
                initial_droplet_number = something(initial_droplet_number, 0),
                radiation = radiation_scheme === :simple_longwave ? "SAM rad_simple longwave (every step)" :
                            "RRTMGP all-sky longwave only (solar constant 0), every 60 s, SST surface, emissivity 0.95, r_eff 10/30 μm",
                surface = "prescribed H, LE and wind-aligned τ (SFC_FLX_FXD, SFC_TAU_FXD)",
                wind_nudging = "none (donudging_uv = .false.)",
                sponge = isnothing(sponge) ? "nothing" : summary(sponge),
                closure = isnothing(closure) ? "nothing" : string(closure), advection_order, scalar_advection_order,
                reference_extension = isnothing(reference_extension) ? "none" : collect(reference_extension),
                stop_time, Δt, max_Δt, cfl,
                time_stepping = max_Δt == Δt ? "fixed Δt = $Δt s" :
                                "adaptive (initial Δt = $Δt s, cfl = $cfl, max_Δt = $max_Δt s); SAM used fixed dt",
                perturbation = string(perturbation),
                write_output, output_dir = abspath(output_dir), output_prefix, profile_interval, timeseries_interval,
                slice_interval, slice_height, extra_slice_heights,
                y_averaged_sections, process_rates, process_rate_interval, overwrite_output, progress_interval,
                checkpoint_interval = something(checkpoint_interval, 0), checkpoint_cleanup)

    return (; simulation, model, grid, config, inputs = paths,
              namelist, soundings, sounding, lsf, sfc, profiles, columns, forcing_profiles, forcing)
end
