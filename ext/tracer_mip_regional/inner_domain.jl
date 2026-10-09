#####
##### TRACER-MIP inner domain (Grid-2, 500 m): offline one-way nest driven by the saved outer state.
#####
##### NumericalEarth nests a Breeze child only in a `PrescribedAtmosphere`. The outer run therefore
##### saves its raw state (T, p, qᵛ, qᶜˡ, qʳ, qⁱ, u, v, w, ρ) on the inner region + halo every 10 min
##### (`inner_region/inner_region_state.jld2`), and `outer_run_parent` turns that file into a
##### `PrescribedAtmosphere` whose grid carries the outer terrain-following heights as a static 3-D
##### `PressureLevelVerticalDiscretization` — the same machinery the ERA5 parent uses, so lateral
##### boundaries, Davies relaxation, terrain blending and initialization apply unchanged.
#####

using NumericalEarth.Grids: PressureLevelVerticalDiscretization
using Oceananigans.OutputReaders: FieldTimeSeries, OnDisk, InMemory
using Oceananigans.Grids: λnodes, φnodes, znode, Center, Face

struct OuterRunParent
    directory :: String
end
Base.summary(p::OuterRunParent) = string("OuterRunParent(", p.directory, ")")

const INNER_REGION_FILE = joinpath("inner_region", "inner_region_state.jld2")

# Move a saved x-face (or y-face) window onto the cell centers of the (ni, nj) sub-region. The writer's
# index window selects faces i₁…i₂ (one fewer than the centers need), so the last center takes the
# one-sided face value; a window that already holds ni + 1 faces is averaged exactly.
function centered_data(fts, n, (ni, nj))
    data = Array(interior(fts[n]))
    LX, LY, _ = Oceananigans.Fields.location(fts)
    if LX === Face
        data = size(data, 1) == ni + 1 ? 0.5 .* (data[1:ni, :, :] .+ data[2:ni+1, :, :]) :
               cat(0.5 .* (data[1:end-1, :, :] .+ data[2:end, :, :]), data[end:end, :, :]; dims = 1)
    elseif LY === Face
        data = size(data, 2) == nj + 1 ? 0.5 .* (data[:, 1:nj, :] .+ data[:, 2:nj+1, :]) :
               cat(0.5 .* (data[:, 1:end-1, :] .+ data[:, 2:end, :]), data[:, end:end, :]; dims = 2)
    end
    size(data)[1:2] == (ni, nj) || throw(DimensionMismatch("saved $(fts.name) window $(size(data)) does not match the ($ni, $nj) sub-region"))
    return data
end

"""
    outer_run_parent(run_dir; arch = CPU(), gravitational_acceleration = 9.80665)

Build a `PrescribedAtmosphere` from an outer run's `inner_region/inner_region_state.jld2`: all saved times, on a
bounded `LatitudeLongitudeGrid` spanning the saved cells whose vertical coordinate is the outer grid's per-column
cell-center height (static), with the terrain height as surface geopotential so `surface_elevation(parent)` works.
Velocities are moved to cell centers; `qⁱ` (total ice) is passed as the parent's `qᶜⁱ`; no snow category.
"""
function outer_run_parent(run_dir; arch = CPU(), gravitational_acceleration = 9.80665)
    FT = Oceananigans.defaults.FloatType
    file = joinpath(run_dir, INNER_REGION_FILE)
    isfile(file) || throw(ArgumentError("no saved inner-region state at $file"))
    T_fts = FieldTimeSeries(file, "T"; backend = OnDisk())
    outer = T_fts.grid                       # the outer terrain-following grid (needs Breeze loaded)
    i_range, j_range, _ = T_fts.indices
    ni, nj = length(i_range), length(j_range)
    Nz = size(outer, 3)
    times = FT.(T_fts.times)

    λf = Array(λnodes(outer, Face())); φf = Array(φnodes(outer, Face()))
    longitude = (λf[first(i_range)], λf[last(i_range) + 1])
    latitude  = (φf[first(j_range)], φf[last(j_range) + 1])

    # Heights of the outer cells on the sub-region (static), as a geopotential on a provisional grid.
    z_static = collect(FT, range(0, 1, length = Nz + 1))
    provisional = LatitudeLongitudeGrid(CPU(), FT; size = (ni, nj, Nz), halo = (5, 5, 5), longitude, latitude,
                                        z = z_static, topology = (Bounded, Bounded, Bounded))
    Φ = CenterField(provisional)
    Φ_sfc = Field{Center, Center, Nothing}(provisional)
    g = FT(gravitational_acceleration)
    for (jj, j) in enumerate(j_range), (ii, i) in enumerate(i_range)
        Φ_sfc[ii, jj, 1] = g * znode(i, j, 1, outer, Center(), Center(), Face())
        for k in 1:Nz
            Φ[ii, jj, k] = g * znode(i, j, k, outer, Center(), Center(), Center())
        end
    end
    fill_halo_regions!(Φ); fill_halo_regions!(Φ_sfc)
    z = PressureLevelVerticalDiscretization(Φ; gravitational_acceleration = g, surface_geopotential = Φ_sfc)
    grid = LatitudeLongitudeGrid(arch, FT; size = (ni, nj, Nz), halo = (5, 5, 5), longitude, latitude, z,
                                 topology = (Bounded, Bounded, Bounded))

    function series(name)
        src = FieldTimeSeries(file, name; backend = OnDisk())
        fts = FieldTimeSeries{Center, Center, Center}(grid, times)
        for n in eachindex(times)
            interior(fts[n]) .= Oceananigans.on_architecture(arch, FT.(centered_data(src, n, (ni, nj))))
        end
        fill_halo_regions!(fts)
        return fts
    end

    u = series("u"); v = series("v"); T = series("T"); qᵛ = series("qᵛ"); p = series("p")
    qᶜˡ = series("qᶜˡ"); qʳ = series("qʳ"); qⁱ = series("qⁱ")
    return PrescribedAtmosphere(grid, times; source = OuterRunParent(abspath(run_dir)),
                                velocities = (; u, v), temperature = T, specific_humidity = qᵛ, pressure = p,
                                microphysical_variables = (; qᶜˡ, qʳ, qᶜⁱ = qⁱ))
end

"""
    tracer_mip_inner_simulation(arch; outer_run_dir, case = :aug07, ...)

Grid-2 (500 m) control nested one-way in the saved outer state (`outer_run_parent`). Same physics
choices as `tracer_mip_outer_simulation` (P3 + Tier-1 aerosol, RRTMGP, slab land with sea pinning, no closure by
default), Δt = 1.5 s, inner extent and 500 m cells from the protocol. Land is initialized from ERA5 at the inner box
(stage it with the fetch script's `inner` group) or, with `land_init = :outer_skin`, from the parent's lowest-level
temperature as a fallback that is recorded in `config`.
"""
function BreezeLab.tracer_mip_inner_simulation(arch; outer_run_dir,
        case = :aug07, protocol = tracer_mip_protocol(),
        Nx = protocol["grids"]["inner"]["Nx"], Ny = protocol["grids"]["inner"]["Ny"],
        z_faces = BreezeLab.acpc_vertical_faces(protocol["vertical"]["scalar_levels_m_agl"]),
        stop_time = nothing, Δt = protocol["grids"]["inner"]["timestep_s"],
        microphysics = nothing,
        radiation = true, radiation_interval = protocol["forcing"]["radiation_interval_s"],
        closure = nothing, terrain = ETOPO2022(), relaxation_width = 5, relaxation_rate = 1/300, advection_order = 5,
        land_albedo = 0.17, sea_albedo = 0.06, surface_emissivity = 0.98,
        land_init = :outer_skin, sea_surface_temperature = nothing, balancer = true, process_rates = true,
        label = "tracer_mip_inner")

    start, stop = tracer_mip_case_window(protocol, case)
    parent = outer_run_parent(outer_run_dir; arch)
    stop_time = something(stop_time, Float64(last(parent.temperature.times)))
    grid = tracer_mip_grid(arch, :inner; protocol, Nx, Ny, z_faces)
    FT = eltype(grid)
    constants = ThermodynamicConstants()
    protocol_microphysics = isnothing(microphysics)
    if protocol_microphysics
        aerosol = PrescribedAerosolProfile(tracer_mip_aerosol_profile(case; protocol); thermodynamic_constants = constants)
        microphysics = P3Microphysics(; cloud = CloudDroplets(; number_concentration = 100e6), aerosol)
    end
    BreezeLab.check_precision(microphysics, FT)
    scalar_advection = scalar_advection_schemes(advection_order, microphysics, :qᵛ; energy_name = :ρθ)
    nest = nested_atmosphere_model(parent, grid; parent_condensates = parent_condensates_of(parent),
                                   base_pressure = FT(mean(interior(parent.pressure[1], :, :, 1))),
                                   relaxation_rate, relaxation_width, terrain, terrain_blend_width = relaxation_width,
                                   microphysics, momentum_advection = WENO(order = advection_order), scalar_advection,
                                   closure, thermodynamic_constants = constants)
    initialize_child_from_parent!(nest; balancer)
    child = nest.child

    extent = BreezeLab.tracer_mip_horizontal_extent(protocol, :inner; FT)
    land_grid = LatitudeLongitudeGrid(arch, FT; size = (Nx, Ny), halo = (5, 5), longitude = extent.longitude,
                                      latitude = extent.latitude, topology = (Bounded, Bounded, Flat))
    land = SlabLand(land_grid)
    mask = sea_mask(grid)
    # Lowest-level air temperature of the initialized child (the parent's state interpolated onto the inner
    # grid) as the slab's initial skin temperature and the default sea surface temperature (static).
    T₁ = Field{Center, Center, Nothing}(land_grid)
    interior(T₁) .= interior(child.temperature, :, :, 1)
    fill_halo_regions!(T₁)
    if land_init === :outer_skin
        set!(land; T = T₁, M = 0.5 * default_saturated_storage(land))
    else
        throw(ArgumentError("land_init must be :outer_skin (ERA5 inner-box initialization not staged yet)"))
    end
    sst = something(sea_surface_temperature, T₁)
    pinning = SeaSurfacePinning(land, mask, sst)
    pin_sea_surface!(pinning, 0.0)

    albedo = surface_albedo_field(mask; sea = sea_albedo, land = land_albedo)
    rtm = radiation ? RadiativeTransferModel(grid, AllSkyOptics(), constants;
                                       solar_position = ApparentSolarPosition(epoch = start), surface_albedo = albedo,
                                       surface_emissivity, schedule = TimeInterval(radiation_interval)) : nothing
    atmosphere = Simulation(nest; Δt)
    model = AtmosphereLandModel(atmosphere, land; radiation = rtm)
    simulation = Simulation(model; Δt, stop_time)
    add_callback!(simulation, pinning, IterationInterval(1))
    accumulators = process_rates ? ProcessRateAccumulators(child) : nothing
    isnothing(accumulators) || add_callback!(simulation, accumulators, IterationInterval(1))
    config = (; label, case = string(case), parent = summary(parent.source), outer_run_dir = abspath(outer_run_dir),
                start = string(start), stop = string(stop), arch = string(typeof(arch)), FT = string(FT), Nx, Ny,
                Nz = length(z_faces) - 1, Δt, stop_time, protocol_microphysics,
                microphysics = protocol_microphysics ? "P3 two-moment, Tier-1 PrescribedAerosolProfile (fixed, height-dependent), radiatively inactive" : summary(microphysics),
                aerosol = isnothing(microphysics.aerosol) ? "none" : summary(microphysics.aerosol),
                radiation = radiation ? "RRTMGP all-sky" : "none", radiation_interval, closure = isnothing(closure) ? "nothing" : summary(closure),
                relaxation_width, relaxation_rate, advection_order, terrain = isnothing(terrain) ? "none (flat)" : summary(terrain),
                land_init = string(land_init), land_albedo, sea_albedo, surface_emissivity,
                nesting = "one-way, offline: parent = outer state saved every 10 min on the inner region + halo")
    return (; simulation, model, nest, child, parent, land, grid, land_grid, sea_mask = mask, pinning, accumulators,
              radiation = rtm, config)
end
