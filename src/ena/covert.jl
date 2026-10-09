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
    covert_public_bin_vertical_faces(; Nz=192, top=20000, Δz=10, uniform_top=1500)

A labelled **reconstruction** of the 192-level, 20-km-top vertical grid of the Covert et al.
(2022) SAM runs (the public repository advertises but does not ship its `grd`/`domain.f90`):
uniform `Δz` to `uniform_top`, then a constant-ratio stretch to `top`.
"""
function covert_public_bin_vertical_faces(; Nz=192, top=20000.0, Δz=10.0, uniform_top=1500.0)
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
`:covert_inversion_5m` ([`covert_inversion_refined_vertical_faces`](@ref)) and `:lasso`
([`lasso_ena_vertical_faces`](@ref)). The name is what run scripts record in provenance.
"""
function ena_vertical_faces(name)
    name = Symbol(name)
    name === :covert && return covert_public_bin_vertical_faces()
    name === :covert_inversion_5m && return covert_inversion_refined_vertical_faces()
    name === :lasso && return lasso_ena_vertical_faces()
    throw(ArgumentError("unknown vertical grid $name (covert, covert_inversion_5m, lasso)"))
end

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
- `Nx`, `Ny` (namelist `caseid`), `Lx = Nx dx`, `Ly = Ny dy`,
  `z_faces` (`covert_public_bin_vertical_faces` with the `caseid` level count)
- `stop_time` (`nstop × dt`), `Δt` (`dt`), `max_Δt = Δt` (fixed step), `cfl = 0.7`
- `translation_velocity` (namelist `ug`, `vg`): SAM's translating frame
- `sponge = SAMSponge()`, `closure = SmagorinskyLilly()`, `advection_order = 5`
- `perturbation = InitialPerturbation()` (`setperturb.f90` case 5)
- `initialization = :condensate_free` or `:equilibrium` (P3 only), `initial_droplet_number`
- output: `write_output`, `output_dir`, `output_prefix`, `profile_interval = 1hour`,
  `timeseries_interval = 60`, `slice_interval = 10minutes`, `slice_height = 900`,
  `progress_interval = 10minutes`, `checkpoint_interval = nothing`, `checkpoint_cleanup = true`
"""
function ena_covert(; arch = GPU(),
                      data_dir = package_path("data", "covert2022_bin"),
                      microphysics = P3Microphysics(; cloud = CloudDroplets(; number_concentration = 75e6)),
                      diagnostic_ccn = has_aerosol_reservoir(microphysics),
                      Nx = nothing,
                      Ny = nothing,
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
                      advection_order = 5,
                      perturbation = InitialPerturbation(),
                      initialization = :condensate_free,
                      initial_droplet_number = nothing,
                      label = "Covert-public-bin development benchmark (not an official LASSO-ENA reproduction)",
                      write_output = true,
                      output_dir = joinpath(pwd(), "output", "ena_covert"),
                      output_prefix = "ena_covert",
                      profile_interval = 1hour,
                      timeseries_interval = 60,
                      slice_interval = 10minutes,
                      slice_height = 900,
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

    overrides = [name for (name, value) in pairs((; Nx, Ny, Lx, Ly, z_faces, stop_time, Δt, translation_velocity))
                 if !isnothing(value)]
    nx, ny, nz = parse_caseid(namelist["caseid"])
    z_faces = isnothing(z_faces) ? covert_public_bin_vertical_faces(; Nz = nz) : z_faces
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
    grid = RectilinearGrid(arch; size = (Nx, Ny, Nz), x = (0, Lx), y = (0, Ly), z = z_faces,
                           halo = (5, 5, 5), topology = (Periodic, Periodic, Bounded))
    FT = eltype(grid)
    constants = ThermodynamicConstants()
    reference_state = sounding_reference_state(grid, profiles, constants)
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
    scalar_advection = scalar_advection_schemes(advection_order, microphysics, moisture_name; energy_name = :ρθ)

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
        :E => (thermodynamic.s, upper.s),            # tls and the top relaxation, as energy tendencies
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

    radiation = SimpleLongwaveRadiation(grid; schedule = IterationInterval(1))

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
    if diagnostic_ccn
        n_initial = sum(mode.number_mixing_ratio for mode in microphysics.aerosol.modes)
        add_callback!(simulation, DiagnosticCCNProjection(model, n_initial), IterationInterval(1))
    end
    add_callback!(simulation, ProgressMessenger(model), TimeInterval(progress_interval))

    if write_output
        mkpath(output_dir)
        add_output_writers!(simulation; output_dir, output_prefix, profile_interval,
                            timeseries_interval, slice_interval, slice_height,
                            energy_budget_series = true,
                            temperature_neutral_evaporation = true)
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
                radiation = "SAM rad_simple longwave (every step)",
                surface = "prescribed H, LE and wind-aligned τ (SFC_FLX_FXD, SFC_TAU_FXD)",
                wind_nudging = "none (donudging_uv = .false.)",
                sponge = isnothing(sponge) ? "nothing" : summary(sponge),
                closure = isnothing(closure) ? "nothing" : summary(closure), advection_order,
                stop_time, Δt, max_Δt, cfl,
                time_stepping = max_Δt == Δt ? "fixed Δt = $Δt s" :
                                "adaptive (initial Δt = $Δt s, cfl = $cfl, max_Δt = $max_Δt s); SAM used fixed dt",
                perturbation = string(perturbation),
                write_output, output_dir = abspath(output_dir), output_prefix, profile_interval, timeseries_interval,
                slice_interval, slice_height, progress_interval,
                checkpoint_interval = something(checkpoint_interval, 0), checkpoint_cleanup)

    return (; simulation, model, grid, config, inputs = paths,
              namelist, soundings, sounding, lsf, sfc, profiles, columns, forcing_profiles, forcing)
end
