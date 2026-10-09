#####
##### TRACER-MIP outer domain (Grid-1): ERA5-driven Breeze child with coastal slab-land exchange.
#####

using NumericalEarth.DataWrangling: NearestNeighborInpainting
using Breeze: materialize_terrain!, moisture_specific_name
using Oceananigans.TimeSteppers: update_state!

"""
    initialize_child_from_parent!(nest; balancer = true)

Initialize a `NestedModel` child from its exchanger's parent-derived prognostics (the state that also
drives the lateral boundaries): interpolate ρᵈ, ρθ, ρqᵛᵉ, ρu, ρv to the child, set the state with a
recomputed reference, remove the terrain-induced contravariant vertical momentum, and run Breeze's
adiabatic balancer. Mirrors NumericalEarth's dataset-path initialization for a hand-built parent.
"""
function initialize_child_from_parent!(nest; balancer = true)
    child = nest.child
    grid = child.grid
    prognostic = nest.exchanger.prognostic
    t₀ = first(prognostic.ρᵈ.times)
    to_child(fts, loc = (Center, Center, Center)) = (f = Field{loc...}(grid); interpolate!(f, fts[Time(t₀)]); f)
    ρᵈ   = to_child(prognostic.ρᵈ)
    ρθ   = to_child(prognostic.ρθ)
    ρqᵛᵉ = to_child(prognostic.ρqᵛᵉ)
    ρu   = to_child(prognostic.ρu, (Face, Center, Center))
    ρv   = to_child(prognostic.ρv, (Center, Face, Center))
    ρ   = Field(ρᵈ + ρqᵛᵉ)
    qᵛᵉ = Field(ρqᵛᵉ / ρ)
    θˡⁱ = Field(ρθ / ρᵈ)
    moisture = NamedTuple{(moisture_specific_name(child.microphysics),)}((qᵛᵉ,))
    set!(nest; ρ, ρu, ρv, θˡⁱ, moisture..., compute_reference_state = true)
    update_state!(nest)
    if !isnothing(child.dynamics.contravariant_vertical_momentum)
        interior(child.momentum.ρw) .-= interior(child.dynamics.contravariant_vertical_momentum)
        update_state!(nest)
    end
    set!(nest; balancer)
    return nest
end

# Hydrometeor fields of a hand-built parent, or `nothing` when it carries none (synthetic parents).
function parent_condensates_of(parent)
    μ = parent.microphysical_variables
    names = (:qᶜˡ, :qʳ, :qᶜⁱ, :qˢ)
    any(n -> haskey(μ, n), names) || return nothing
    return NamedTuple{names}(map(n -> get(μ, n, nothing), names))
end

"""
    tracer_mip_outer_simulation(arch; case = :aug07, parent = :era5, era5_dir, stop_time = 24hours, kwargs...)

Build the TRACER-MIP Grid-1 (outer, 2 km) control as an Oceananigans `Simulation` of a NumericalEarth
`AtmosphereLandModel` wrapping a Breeze compressible child nested in the ERA5 pressure-level parent.
Returns a named tuple with `simulation`, `model`, `nest`, `land`, `grid`, `accumulators`, `sea_mask`,
`config` and `parent`. Nothing is run.

- `parent = :era5` uses `ERA5HourlyPressureLevels()` files in `era5_dir` (stage them with
  `data_wrangling/fetch_era5_tracer_mip.jl`); `parent = :synthetic` (software testing only) builds
  `synthetic_parent_atmosphere` and labels the configuration `exploratory_synthetic_boundaries = true`.
- `Nx`, `Ny` default to the protocol's 750²; reduce them for tests (extent stays the protocol's).
- `microphysics = nothing` builds the Tier-1 protocol P3: prognostic droplet number activated from the
  case's prescribed two-mode profile (`PrescribedAerosolProfile(tracer_mip_aerosol_profile(case))`) with a
  100 cm⁻³ construction-time droplet number. For the aerosol sensitivities pass the Breeze object, e.g.
  `P3Microphysics(; cloud, aerosol = PrescribedAerosolProfile(tracer_mip_aerosol_profile(case; multiplier = 3)))`
  (HIGH). Precision follows `Oceananigans.defaults.FloatType`.
- Land: `SlabLand` initialized from ERA5 skin temperature and ERA5-Land soil moisture (or constants for
  the synthetic parent); sea cells (terrain ≤ 0 m) pinned to the ERA5 skin temperature each step.
- Radiation: RRTMGP all-sky every `radiation_interval` (protocol 60 s), constant land/sea albedo.
- `Δt` fixed at the protocol's 3 s. `closure` defaults to Breeze's `TKEBasedTurbulenceClosure()` (vertical eddy
  diffusivity with prognostic TKE, vertically implicit): without a closure (job 233) the surface stress is deposited only
  in the 50 m lowest layer, which decelerates to ≈0.6 of ERA5's 10 m wind while the level above stays at the free-stream
  value (a 2–3× jump across one cell). `SmagorinskyLilly` is unusable here: its isotropic filter width (2000·2000·50)^(1/3)
  ≈ 585 m gives ν ≈ (0.16·585)²|S| ≳ 10³ m²/s near the surface, so the explicit vertical diffusion violates νΔt/Δz² ≲ 1/2
  on 50 m cells and the first step blows up (negative pressure). A horizontal-only Smagorinsky would be the consistent
  mesoscale companion to a 1-D PBL closure, but the pinned Breeze cannot combine two closures (tuples of closures are
  post-pin, #1036), so horizontal mixing remains the advection scheme's.
"""
function BreezeLab.tracer_mip_outer_simulation(arch;
        case = :aug07,
        parent = :era5,
        era5_dir = joinpath(pkgdir(BreezeLab), "data", "era5"),
        protocol = tracer_mip_protocol(),
        Nx = protocol["grids"]["outer"]["Nx"],
        Ny = protocol["grids"]["outer"]["Ny"],
        z_faces = BreezeLab.acpc_vertical_faces(protocol["vertical"]["scalar_levels_m_agl"]),
        stop_time = Hour(protocol["cases"][String(case)]["duration_hours"]).value * 3600,
        Δt = protocol["grids"]["outer"]["timestep_s"],
        microphysics = nothing,
        radiation = true,
        radiation_interval = protocol["forcing"]["radiation_interval_s"],
        closure = TKEBasedTurbulenceClosure(),
        terrain = parent === :era5 ? ETOPO2022() : nothing,
        relaxation_width = 5,
        relaxation_rate = 1/300,
        advection_order = 5,
        land_albedo = 0.17, sea_albedo = 0.06, surface_emissivity = 0.98,
        soil_porosity = 0.45,
        synthetic_kw = NamedTuple(),
        balancer = true,
        process_rates = true,
        label = "tracer_mip_outer")

    start, stop = tracer_mip_case_window(protocol, case)
    grid = tracer_mip_grid(arch, :outer; protocol, Nx, Ny, z_faces)
    FT = eltype(grid)
    constants = ThermodynamicConstants()

    # Tier-1 microphysics: P3 with droplet number activated from the prescribed profile.
    protocol_microphysics = isnothing(microphysics)
    if protocol_microphysics
        aerosol = PrescribedAerosolProfile(tracer_mip_aerosol_profile(case; protocol); thermodynamic_constants = constants)
        microphysics = P3Microphysics(; cloud = CloudDroplets(; number_concentration = 100e6), aerosol)
    end
    BreezeLab.check_precision(microphysics, FT)
    scalar_advection = scalar_advection_schemes(advection_order, microphysics, :qᵛ; energy_name = :ρθ)
    momentum_advection = WENO(order = advection_order)

    nest_kw = (; relaxation_rate, relaxation_width, terrain, terrain_blend_width = relaxation_width,
                 microphysics = microphysics, momentum_advection, scalar_advection, closure, thermodynamic_constants = constants)

    dates = (start, stop)
    exploratory = parent !== :era5
    nest, parent_atmosphere = if parent === :era5
        nest = nested_atmosphere_model(grid, ERA5HourlyPressureLevels(); dates, dir = era5_dir, balancer, nest_kw...)
        nest, nest.parent
    elseif parent === :synthetic
        times = collect(0.0:3600.0:(stop_time + 3600.0))
        parent_atmosphere = synthetic_parent_atmosphere(grid; times, synthetic_kw...)
        if terrain isa Field
            # mirror NumericalEarth: materialize the child's terrain before the model is built
            materialize_terrain!(grid, terrain)
        end
        nest = nested_atmosphere_model(parent_atmosphere, grid; parent_condensates = parent_condensates_of(parent_atmosphere),
                                       base_pressure = FT(get(synthetic_kw, :surface_pressure, 101300.0)), nest_kw...)
        initialize_child_from_parent!(nest; balancer)
        nest, parent_atmosphere
    elseif parent isa PrescribedAtmosphere
        nest = nested_atmosphere_model(parent, grid; parent_condensates = parent_condensates_of(parent), nest_kw...)
        initialize_child_from_parent!(nest; balancer)
        nest, parent
    else
        throw(ArgumentError("parent must be :era5, :synthetic or a PrescribedAtmosphere"))
    end
    child = nest.child

    #####
    ##### Land, sea mask, pinning
    #####

    extent = BreezeLab.tracer_mip_horizontal_extent(protocol, :outer; FT)
    land_grid = LatitudeLongitudeGrid(arch, FT; size = (Nx, Ny), halo = (5, 5),
                                      longitude = extent.longitude, latitude = extent.latitude,
                                      topology = (Bounded, Bounded, Flat))
    land = SlabLand(land_grid)
    mask = sea_mask(grid)
    land_region = BoundingBox(land_grid)
    inpainting = NearestNeighborInpainting(200)
    if parent === :era5
        single = ERA5HourlySingleLevel()
        set!(land.temperature, Metadatum(:skin_temperature; dataset = single, date = start, region = land_region, dir = era5_dir))
        swvl = Field{Center, Center, Nothing}(land_grid)
        set!(swvl, Metadatum(:volumetric_soil_water_layer_1; dataset = ERA5HourlyLand(), date = start, region = land_region, dir = era5_dir);
             inpainting)
        Mmax = default_saturated_storage(land)
        M₀ = Field{Center, Center, Nothing}(land_grid)
        interior(M₀) .= Mmax .* clamp.(interior(swvl) ./ soil_porosity, 0, 1)   # NaN (sea) cells are pinned below
        interior(M₀)[isnan.(interior(M₀))] .= Mmax
        set!(land; M = M₀)
        skt = FieldTimeSeries(Metadata(:skin_temperature; dataset = single, dates, region = land_region, dir = era5_dir), arch;
                              time_indices_in_memory = 8)
        sst = skt
    else
        T_sfc = FT(get(synthetic_kw, :surface_temperature, 300.0))
        set!(land; T = T_sfc, M = 0.5 * default_saturated_storage(land))
        sst = T_sfc - 1   # idealized: a sea 1 K cooler than the initial land
    end
    pinning = SeaSurfacePinning(land, mask, sst)
    pin_sea_surface!(pinning, 0.0)

    #####
    ##### Radiation and coupled model
    #####

    albedo = surface_albedo_field(mask; sea = sea_albedo, land = land_albedo)
    rtm = if radiation
        RadiativeTransferModel(grid, AllSkyOptics(), constants;
                               solar_position = ApparentSolarPosition(epoch = start),
                               surface_albedo = albedo, surface_emissivity,
                               schedule = TimeInterval(radiation_interval))
    else
        nothing
    end

    atmosphere = Simulation(nest; Δt)
    model = AtmosphereLandModel(atmosphere, land; radiation = rtm)
    simulation = Simulation(model; Δt, stop_time)
    add_callback!(simulation, pinning, IterationInterval(1))

    accumulators = process_rates ? ProcessRateAccumulators(child) : nothing
    isnothing(accumulators) || add_callback!(simulation, accumulators, IterationInterval(1))

    config = (; label, case = string(case), parent = parent === :era5 ? "ERA5HourlyPressureLevels" : summary(parent_atmosphere.source),
                exploratory_synthetic_boundaries = exploratory, start = string(start), stop = string(stop),
                arch = string(typeof(arch)), FT = string(FT), Nx, Ny, Nz = length(z_faces) - 1,
                longitude = collect(extent.longitude), latitude = collect(extent.latitude),
                Δt, stop_time, protocol_microphysics,
                microphysics = protocol_microphysics ? "P3 two-moment, Tier-1 PrescribedAerosolProfile (fixed, height-dependent), radiatively inactive" : summary(microphysics),
                aerosol = isnothing(microphysics.aerosol) ? "none" : summary(microphysics.aerosol),
                radiation = radiation ? "RRTMGP all-sky" : "none", radiation_interval,
                closure = summary(closure), relaxation_width, relaxation_rate, advection_order,
                terrain = isnothing(terrain) ? "none (flat)" : summary(terrain), land_albedo, sea_albedo, surface_emissivity,
                land = "SlabLand (skin temperature + bucket); sea cells pinned to prescribed skin temperature with saturated bucket",
                process_rates, era5_dir = parent === :era5 ? abspath(era5_dir) : "none")

    return (; simulation, model, nest, child, parent = parent_atmosphere, land, grid, land_grid, sea_mask = mask,
              pinning, accumulators, radiation = rtm, config)
end
