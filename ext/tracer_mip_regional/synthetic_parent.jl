#####
##### SYNTHETIC parent atmosphere — SOFTWARE TESTING ONLY.
#####
##### An analytic, horizontally uniform, hydrostatic `PrescribedAtmosphere` on a coarse lat-lon grid
##### covering the child's padded box, so that the open-boundary / nesting / land-exchange machinery
##### can be exercised on CPU (and in exploratory GPU throughput runs) without reanalysis data. It is
##### never a stand-in for ERA5 and must not be presented as a MIP control. Every object built from it
##### carries `source = SyntheticTestParent()`.
#####

struct SyntheticTestParent end
Base.summary(::SyntheticTestParent) = "SyntheticTestParent (idealized, software testing only)"

"""
    synthetic_parent_atmosphere(child_grid; times = 0:3600:7200, padding, spacing = 0.25, Nz = 48, z_top = 24e3,
                                surface_pressure = 101300, surface_temperature = 300, lapse_rate = 6.5e-3,
                                tropopause = 15e3, surface_specific_humidity = 0.016, humidity_scale_height = 2500,
                                u = 3, v = 2)

Hydrostatic analytic state: T(z) = T₀ − Γ·min(z, z_trop); p from the hydrostatic integral of that profile with the
dry gas constant; qᵛ(z) = q₀ exp(−z/H_q); uniform winds (u, v); no condensate. Times are seconds.
"""
function synthetic_parent_atmosphere(child_grid; times = collect(0.0:3600.0:7200.0),
                                     padding = default_horizontal_padding(ERA5HourlyPressureLevels()),
                                     spacing = 0.25, Nz = 48, z_top = 24e3,
                                     surface_pressure = 101300.0, surface_temperature = 300.0, lapse_rate = 6.5e-3,
                                     tropopause = 15e3, surface_specific_humidity = 0.016, humidity_scale_height = 2500.0,
                                     u = 3.0, v = 2.0)
    arch = child_grid.architecture
    FT = eltype(child_grid)
    box = BoundingBox(child_grid; padding)
    # ERA5-like spacing: parent nodes must bracket the child domain, so the half-cell must stay below `padding`
    Nx = max(5, ceil(Int, (box.longitude[2] - box.longitude[1]) / spacing))
    Ny = max(5, ceil(Int, (box.latitude[2] - box.latitude[1]) / spacing))
    parent_grid = LatitudeLongitudeGrid(arch, FT; size = (Nx, Ny, Nz), halo = (5, 5, 5),
                                        longitude = box.longitude, latitude = box.latitude, z = (0, z_top),
                                        topology = (Bounded, Bounded, Bounded))
    g, Rᵈ = 9.80665, 287.05
    T_of(z) = surface_temperature - lapse_rate * min(z, tropopause)
    T_trop = surface_temperature - lapse_rate * tropopause
    p_trop = surface_pressure * (T_trop / surface_temperature)^(g / (Rᵈ * lapse_rate))
    p_of(z) = z ≤ tropopause ? surface_pressure * (T_of(z) / surface_temperature)^(g / (Rᵈ * lapse_rate)) :
                               p_trop * exp(-g * (z - tropopause) / (Rᵈ * T_trop))
    q_of(z) = surface_specific_humidity * exp(-z / humidity_scale_height)

    parent = PrescribedAtmosphere(parent_grid, FT.(times); source = SyntheticTestParent())
    set!(parent.temperature,       (λ, φ, z, t) -> T_of(z))
    set!(parent.specific_humidity, (λ, φ, z, t) -> q_of(z))
    set!(parent.pressure,          (λ, φ, z, t) -> p_of(z))
    set!(parent.velocities.u,      (λ, φ, z, t) -> u)
    set!(parent.velocities.v,      (λ, φ, z, t) -> v)
    return parent
end

"""
    synthetic_coastal_elevation(grid; coast_longitude, slope = 1e-3, maximum = 60)

An analytic elevation `Field` (m): sea (0) east/south of a straight coastline rising linearly inland. Testing only.
"""
function synthetic_coastal_elevation(grid; coast_longitude = nothing, slope = 1e-3, maximum = 60.0)
    λ₁, λ₂ = grid.λᶠᵃᵃ[1], grid.λᶠᵃᵃ[grid.Nx + 1]
    λc = something(coast_longitude, (λ₁ + λ₂) / 2)
    elevation = Field{Center, Center, Nothing}(grid)
    m_per_deg = π * 6371e3 / 180
    set!(elevation, (λ, φ) -> clamp(slope * (λc - λ) * m_per_deg, 0, maximum))
    fill_halo_regions!(elevation)
    return elevation
end
