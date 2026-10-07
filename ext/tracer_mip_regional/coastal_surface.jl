#####
##### Coastal land/sea exchange for a Breeze atmosphere over `SlabLand`.
#####
##### NumericalEarth computes atmosphere-land fluxes on every cell of the land grid ("at most one
##### surface type per cell"). For the Houston domain (half Gulf of Mexico) BreezeLab keeps a single
##### `SlabLand` under the whole atmosphere and, at every step, pins the sea cells to a prescribed sea
##### surface temperature with a saturated bucket, so Monin-Obukhov fluxes over water see the SST and
##### evaporate at the atmospheric-demand limit. Sea cells are those whose terrain is at or below sea
##### level. Departures: sea roughness is the land closure's; no wave-state dependence.
#####

"""
    sea_mask(grid; threshold = 0)

`Field{Center, Center, Nothing}` equal to 1 where the bottom face of the (terrain-following) grid is at or
below `threshold` meters, 0 elsewhere. On a flat grid without terrain every cell is sea.
"""
function sea_mask(grid; threshold = 0)
    mask = Field{Center, Center, Nothing}(grid)
    launch!(grid.architecture, grid, :xy, _sea_mask!, mask, grid, convert(eltype(grid), threshold))
    fill_halo_regions!(mask)
    return mask
end

@kernel function _sea_mask!(mask, grid, threshold)
    i, j = @index(Global, NTuple)
    z₀ = znode(i, j, 1, grid, Center(), Center(), Face())
    @inbounds mask[i, j, 1] = ifelse(z₀ <= threshold, one(z₀), zero(z₀))
end

"""
    surface_albedo_field(mask; sea = 0.06, land = 0.17)

Broadband surface albedo from the sea mask (constant values; the CGLS albedo product is not used).
"""
function surface_albedo_field(mask; sea = 0.06, land = 0.17)
    α = Field{Center, Center, Nothing}(mask.grid)
    set!(α, mask * sea + (1 - mask) * land)
    fill_halo_regions!(α)
    return α
end

struct SeaSurfacePinning{M, S, L, W, FT}
    mask :: M                   # 1 at sea
    sea_surface_temperature :: S # FieldTimeSeries (any grid), 2-D Field, or Number [K]
    land :: L
    work :: W                   # 2-D scratch on the land grid
    saturated_storage :: FT     # bucket content imposed at sea [kg m⁻²]
end

"""
    SeaSurfacePinning(land, mask, sst; saturated_storage = maximum bucket capacity)

Callback (register with `IterationInterval(1)`) that overwrites the slab's temperature with `sst` and its
water storage with `saturated_storage` wherever `mask == 1`, then refreshes the slab diagnostics.
"""
function SeaSurfacePinning(land, mask, sst; saturated_storage = default_saturated_storage(land))
    work = Field{Center, Center, Nothing}(land.grid)
    FT = eltype(land.grid)
    return SeaSurfacePinning(mask, sst, land, work, convert(FT, saturated_storage))
end

default_saturated_storage(land) = land.hydrology isa NumericalEarth.Lands.BucketHydrology ?
    maximum_storage_value(land.hydrology.maximum_water_storage) : 150

maximum_storage_value(x::Number) = x
maximum_storage_value(x) = maximum(x)   # per-cell Field capacity

sea_surface_now!(work, sst::Number, t) = (fill!(interior(work), sst); work)
sea_surface_now!(work, sst::Field, t) = sst
function sea_surface_now!(work, sst::FieldTimeSeries, t)
    interpolate!(work, sst[Time(t)])
    return work
end

@kernel function _pin_sea!(T, M, mask, sst, saturated)
    i, j = @index(Global, NTuple)
    @inbounds begin
        at_sea = mask[i, j, 1] > 0.5
        T[i, j, 1] = ifelse(at_sea, sst[i, j, 1], T[i, j, 1])
        M[i, j, 1] = ifelse(at_sea, saturated, M[i, j, 1])
    end
end

function (p::SeaSurfacePinning)(simulation)
    pin_sea_surface!(p, simulation.model.clock.time)
    return nothing
end

function pin_sea_surface!(p::SeaSurfacePinning, t)
    land = p.land
    sst = sea_surface_now!(p.work, p.sea_surface_temperature, t)
    launch!(land.grid.architecture, land.grid, :xy, _pin_sea!, land.temperature, land.water_storage, p.mask, sst, p.saturated_storage)
    Oceananigans.TimeSteppers.update_state!(land)
    return nothing
end
