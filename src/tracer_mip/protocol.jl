#####
##### TRACER-MIP protocol metadata: grids, vertical levels and case dates.
#####
##### The numbers come from `cases/tracer_mip/protocol.toml`, a transcription of the pinned
##### roadmap (ARM-Synergy/tracer-mip @ 416423ac). Functions here turn them into Breeze
##### grids; the choices that are not protocol values (lat-lon map, face derivation) are
##### documented in the TOML `[breezelab]` table and in docs/cases/tracer_mip.md.
#####

using TOML: TOML
using Dates: DateTime, Hour

const TRACER_MIP_PROTOCOL_PATH = normpath(joinpath(@__DIR__, "..", "..", "cases", "tracer_mip", "protocol.toml"))

"""
    tracer_mip_protocol(path = TRACER_MIP_PROTOCOL_PATH)

Parse the transcribed TRACER-MIP protocol (a `Dict`). Validates the invariants used
downstream: two 24-hour cases starting 06 UTC, 750²/500² nests at 2000/500 m, 95 quoted
scalar levels, and two aerosol modes per case.
"""
function tracer_mip_protocol(path = TRACER_MIP_PROTOCOL_PATH)
    protocol = TOML.parsefile(path)
    for case in ("jun17", "aug07")
        c = protocol["cases"][case]
        c["duration_hours"] == 24 || error("TRACER-MIP cases are 24 h, got $(c["duration_hours"]) for $case")
        Dates.hour(c["start_utc"]) == 6 || error("TRACER-MIP cases start at 06 UTC, got $(c["start_utc"]) for $case")
        length(protocol["aerosol"][case]) == 2 || error("expected two aerosol modes for $case")
    end
    g = protocol["grids"]
    (g["outer"]["Nx"], g["outer"]["Ny"], g["outer"]["spacing_m"]) == (750, 750, 2000.0) || error("outer grid is 750² at 2 km")
    (g["inner"]["Nx"], g["inner"]["Ny"], g["inner"]["spacing_m"]) == (500, 500, 500.0) || error("inner grid is 500² at 500 m")
    length(protocol["vertical"]["scalar_levels_m_agl"]) == protocol["vertical"]["n_levels_quoted"] == 95 ||
        error("expected the 95 ACPC scalar levels")
    return protocol
end

"""
    tracer_mip_case_window(protocol, case)

`(start, stop)` UTC `DateTime`s of a case (`:aug07` or `:jun17`).
"""
function tracer_mip_case_window(protocol, case)
    c = protocol["cases"][String(case)]
    start = DateTime(c["start_utc"])
    return start, start + Hour(c["duration_hours"])
end

"""
    acpc_vertical_faces(scalar_levels = protocol levels)

Cell interfaces (m AGL) for the 95 ACPC/TRACER-MIP scalar levels. The scalar levels are RAMS
`zt` heights whose first entry (-24 m) lies below ground; interfaces are their midpoints
(`zm`), the bottom face is 0 m and the top face extrapolates the last spacing. 94 cells,
Δz = 50 m at the surface, 300 m aloft, top at 22181 m.
"""
function acpc_vertical_faces(scalar_levels = tracer_mip_protocol()["vertical"]["scalar_levels_m_agl"])
    zt = Float64.(scalar_levels)
    issorted(zt) || throw(ArgumentError("scalar levels must increase with height"))
    zt[1] < 0 < zt[2] || throw(ArgumentError("the first scalar level must be below ground and the second above"))
    faces = [(zt[k] + zt[k + 1]) / 2 for k in 1:length(zt) - 1]
    push!(faces, zt[end] + (zt[end] - zt[end - 1]) / 2)
    faces[1] == 0 || (faces[1] = 0.0)   # -24/+24 midpoint is exactly the surface
    return faces
end

const EARTH_RADIUS = 6371e3
meters_per_degree(R = EARTH_RADIUS) = π * R / 180

"""
    tracer_mip_horizontal_extent(protocol, nest; FT=Float64)

Longitude and latitude bounds `(; longitude, latitude)` of `nest` (`:outer` or `:inner`) on a
latitude-longitude grid whose cells are `spacing_m` square at the shared center. This is
BreezeLab's "similar model option" to the roadmap's polar-stereographic grid: a `Nx × Ny`
`LatitudeLongitudeGrid` spanning `Ny Δ` meridionally and `Nx Δ / cos φ₀` zonally.
"""
function tracer_mip_horizontal_extent(protocol, nest; FT = Float64)
    g = protocol["grids"]
    n = g[String(nest)]
    φ₀, λ₀ = g["center_latitude"], g["center_longitude"]
    Lφ = n["Ny"] * n["spacing_m"] / meters_per_degree()
    Lλ = n["Nx"] * n["spacing_m"] / (meters_per_degree() * cosd(φ₀))
    longitude = (FT(λ₀ - Lλ / 2), FT(λ₀ + Lλ / 2))
    latitude  = (FT(φ₀ - Lφ / 2), FT(φ₀ + Lφ / 2))
    return (; longitude, latitude)
end

"""
    tracer_mip_grid(arch, nest; protocol, FT=Float32, Nx, Ny, z_faces, terrain_following=true, halo=(5,5,5))

The bounded `LatitudeLongitudeGrid` of `nest` (`:outer` 750² at 2 km, `:inner` 500² at 500 m) over
the protocol extent with the ACPC vertical faces. `Nx`/`Ny` may be reduced for tests (the
spacing then coarsens; the extent is the protocol's). With `terrain_following = true` the
vertical is wrapped in Breeze's `TerrainFollowingVerticalDiscretization` so terrain can be
materialized later by the regional constructor.
"""
function tracer_mip_grid(arch, nest; protocol = tracer_mip_protocol(), FT = Float32,
                         Nx = protocol["grids"][String(nest)]["Nx"],
                         Ny = protocol["grids"][String(nest)]["Ny"],
                         z_faces = acpc_vertical_faces(protocol["vertical"]["scalar_levels_m_agl"]),
                         terrain_following = true, halo = (5, 5, 5))
    extent = tracer_mip_horizontal_extent(protocol, nest; FT)
    Nz = length(z_faces) - 1
    z = terrain_following ? TerrainFollowingVerticalDiscretization(FT.(z_faces)) : FT.(z_faces)
    return LatitudeLongitudeGrid(arch, FT; size = (Nx, Ny, Nz), halo,
                                 longitude = extent.longitude, latitude = extent.latitude, z,
                                 topology = (Bounded, Bounded, Bounded))
end
