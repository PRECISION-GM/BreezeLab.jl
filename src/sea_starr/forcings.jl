#####
##### SEA STARR inversion-following free-tropospheric nudging.
#####
##### Protocol (setup PDF / driver attributes): nudging only above the *domain maximum*
##### inversion height (per-column height of the maximum ∂θₗ/∂z, maximum over columns)
##### plus 100 m, ramped up between that level and 200 m higher following Blossey et al.
##### (2013); θ, qₜ, Nₐ with τ = 1800 s and u, v with τ = 10800 s. The driver's
##### `nudging_constant_*` are uniform, so the mask is diagnosed here every time step.
#####
##### The nudging relaxes the *horizontal mean* profile (as `MeanProfileNudging` does for the
##### LASSO winds), so resolved free-tropospheric eddies are not damped. Labelled choice: the
##### Blossey et al. ramp is implemented as linear in height (the paper could not be
##### retrieved; the driver only states the 200 m gap).
#####

using Adapt: Adapt, adapt
using KernelAbstractions: @kernel, @index
using Oceananigans: Field, Average, KernelFunctionOperation, compute!, set!
using Oceananigans.Grids: Center, Face, znode, znodes
using Oceananigans.Fields: interior
using Oceananigans.Units: Time
using Oceananigans.Utils: prettysummary
using Breeze.AtmosphereModels: AtmosphereModels

#####
##### Per-column inversion height
#####

# Height of the maximum vertical gradient of ϕ in column (i, j): the midpoint between the
# two cell centers bracketing the largest (ϕ[k+1] - ϕ[k]) / Δz, searched below z_max.
@inline function column_inversion_height(i, j, k, grid, ϕ, z_max)
    Nz = size(grid, 3)
    FT = eltype(grid)
    best = FT(-Inf)
    z_best = zero(FT)
    for kk in 1:Nz-1
        z_lo = znode(i, j, kk, grid, Center(), Center(), Center())
        z_hi = znode(i, j, kk + 1, grid, Center(), Center(), Center())
        @inbounds g = (ϕ[i, j, kk + 1] - ϕ[i, j, kk]) / (z_hi - z_lo)
        take = (g > best) & (z_hi ≤ z_max)
        best = ifelse(take, g, best)
        z_best = ifelse(take, (z_lo + z_hi) / 2, z_best)
    end
    return z_best
end

"""
    inversion_height_field(model; z_max=Inf)

2D `Field` of the per-column inversion height [m]: the height of the maximum vertical
gradient of the liquid-ice potential temperature `θ` (cell-center differences), restricted
to levels below `z_max`. `compute!` it to refresh.
"""
function inversion_height_field(model; z_max=Inf)
    grid = model.grid
    θ = Oceananigans.fields(model).θ
    op = KernelFunctionOperation{Center, Center, Nothing}(column_inversion_height, grid, θ, convert(eltype(grid), z_max))
    return Field(op)
end

#####
##### Nudging mask
#####

"""
    nudging_mask_weights(z, z_base; ramp_depth=200)

Weights in [0, 1] at heights `z` for nudging that starts at `z_base` and reaches full
strength at `z_base + ramp_depth` (linear ramp).
"""
nudging_mask_weights(z, z_base; ramp_depth=200) = [clamp((ζ - z_base) / ramp_depth, 0, 1) for ζ in z]

"""
    InversionMaskUpdater(model; mask=Field{Nothing, Nothing, Center}(model.grid), offset=100, ramp_depth=200, z_max=Inf)

Callback (`IterationInterval(1)`) maintaining the shared nudging mask
`Field{Nothing, Nothing, Center}`: it computes the per-column inversion height from the
current `θ`, takes the maximum over the domain, and sets the mask to the linear ramp from
`max(zᵢ) + offset` to `max(zᵢ) + offset + ramp_depth`. The latest values are kept in
`inversion_height_max[]` and `mask_base[]` for diagnostics; `inversion_height` is the
2D field (workbook `zi`).
"""
struct InversionMaskUpdater{M, Z, FT}
    mask :: M
    inversion_height :: Z
    z_centers :: Vector{FT}
    offset :: FT
    ramp_depth :: FT
    inversion_height_max :: Base.RefValue{FT}
    mask_base :: Base.RefValue{FT}
end

function InversionMaskUpdater(model; mask=Field{Nothing, Nothing, Center}(model.grid), offset=100, ramp_depth=200, z_max=Inf)
    grid = model.grid
    FT = eltype(grid)
    zi = inversion_height_field(model; z_max)
    z = FT.(Array(znodes(grid, Center())))
    updater = InversionMaskUpdater(mask, zi, z, FT(offset), FT(ramp_depth), Ref(FT(NaN)), Ref(FT(NaN)))
    return updater
end

function (u::InversionMaskUpdater)(simulation_or_model)
    compute!(u.inversion_height)
    zmax = maximum(interior(u.inversion_height))
    u.inversion_height_max[] = zmax
    base = zmax + u.offset
    u.mask_base[] = base
    w = nudging_mask_weights(u.z_centers, base; ramp_depth=u.ramp_depth)
    set!(u.mask, reshape(w, 1, 1, length(w)))
    return nothing
end

#####
##### Masked mean-profile nudging
#####

"""
    InversionFollowingNudging(target, mask; timescale)

Relax the horizontal mean of a prognostic toward the profile time series `target(z, t)`
with rate `mask(z) / timescale`, where `mask` is the shared
`Field{Nothing, Nothing, Center}` maintained by [`InversionMaskUpdater`](@ref):

    Fϕ(z, t) = -mask(z) (⟨ϕ⟩(z, t) - target(z, t)) / timescale

Supply under a *specific* key (`θ`, `qᵛ`, `nᵃ`, `u`, `v`). For microphysical names the
horizontal mean of the specific quantity is formed as ⟨ρϕ⟩ / ρᵣ (anelastic reference density).
"""
struct InversionFollowingNudging{T, A, M, R}
    target :: T
    averaged_field :: A
    mask :: M
    rate :: R
end

InversionFollowingNudging(target, mask; timescale) = InversionFollowingNudging(target, nothing, mask, 1 / timescale)

Adapt.adapt_structure(to, f::InversionFollowingNudging) =
    InversionFollowingNudging(adapt(to, f.target), adapt(to, f.averaged_field), adapt(to, f.mask), adapt(to, f.rate))

Base.summary(f::InversionFollowingNudging) = string("InversionFollowingNudging(timescale=", prettysummary(1 / f.rate), " s, masked above the inversion)")
Base.show(io::IO, f::InversionFollowingNudging) = print(io, summary(f))

@inline function (f::InversionFollowingNudging)(i, j, k, grid, clock, fields)
    @inbounds begin
        ϕ̄ = f.averaged_field[1, 1, k]
        m = f.mask[1, 1, k]
    end
    ϕᵗ = profile_value(f.target, k, clock.time)
    return - f.rate * m * (ϕ̄ - ϕᵗ)
end

function AtmosphereModels.materialize_atmosphere_model_forcing(f::InversionFollowingNudging,
                                                               field, name, model_field_names, context::NamedTuple)
    startswith(string(name), "ρ") &&
        throw(ArgumentError("InversionFollowingNudging returns a specific tendency; supply it under the specific key (e.g. `θ`, `nᵃ`)"))
    averaged_field = if name ∈ keys(context.specific_fields)
        Field(Average(context.specific_fields[name], dims=(1, 2)))
    else
        # Breeze hands the ρ-weighted prognostic for the other specific names.
        Field(Average(field / context.total_density, dims=(1, 2)))
    end
    FT = eltype(field.grid)
    return InversionFollowingNudging(f.target, averaged_field, f.mask, convert(FT, f.rate))
end

function AtmosphereModels.compute_forcing!(f::InversionFollowingNudging)
    compute!(f.averaged_field)
    return nothing
end
