#####
##### SEA STARR radiation: the composite atmosphere above the LES top inside the RRTMGP column,
##### and the solar position following the composite trajectory.
#####
##### The intercomparison models compute radiation "for the simulated properties below 6.5 km and
##### a constant atmospheric profile above 6.5 km taken from the GEOS-5 composite" and vary the
##### solar zenith angle "with the latitude and longitude of the Lagrangian trajectory"
##### (Diamond et al. 2026, egusphere-2026-5350, Sect. 2.2). The pinned Breeze RRTMGP column ends
##### at the LES top with zero incident longwave and a TOA solar flux, and takes a fixed
##### coordinate. `extended_column_radiation` builds Breeze's all-sky RadiativeTransferModel
##### with extra RRTMGP layers above the LES top filled once from the driver's time-mean upper
##### atmosphere (p, T, qᵛ, O₃ up to 76 km); Breeze's own update kernels keep writing the LES
##### layers, and the fluxes copied back to the LES faces now include the downwelling radiation
##### from above. `TrajectorySolarPosition` sets the solver's cos(zenith) from the trajectory's
##### mean (lon, lat) every step.
#####

using Dates: DateTime, Millisecond
using ClimaComms: ClimaComms
using RRTMGP: RRTMGP, AllSkyRadiation, RRTMGPSolver, RRTMGPGridParams, lookup_tables
using RRTMGP.AtmosphericStates: AtmosphericState, CloudState, MaxRandomOverlap
using RRTMGP.BCs: LwBCs, SwBCs
using Oceananigans.Fields: ZFaceField, CenterField, ConstantField
using Oceananigans.Architectures: architecture
using Breeze.AtmosphereModels: RadiativeTransferModel, SurfaceRadiation, FixedCosineZenith, materialize_background_atmosphere
using Breeze.CelestialMechanics: cos_solar_zenith_angle

function rrtmgp_extension()
    ext = Base.get_extension(Breeze, :BreezeRRTMGPExt)
    isnothing(ext) && error("Breeze's RRTMGP extension is not loaded (using RRTMGP, NCDatasets)")
    return ext
end

const AVOGADRO = 6.02214076e23

"""
    upper_atmosphere_layers(driver, z_top; constants)

Time-mean driver column above `z_top` as RRTMGP layers: `(; z, p, T, q, o3, p_faces, T_faces)`
with layer centers at the driver mid-levels above `z_top`, faces midway between them (the
lowest face at `z_top` itself), pressures interpolated in log p and temperatures linearly in z.
"""
function upper_atmosphere_layers(driver::DEPHYDriver, z_top)
    sel = findall(>(z_top), driver.z)
    isempty(sel) && throw(ArgumentError("the driver has no levels above $z_top m"))
    zc = driver.z[sel]
    p̄ = vec(mean(driver.forcing.p, dims=2)); T̄ = vec(mean(driver.forcing.T_nud, dims=2))
    q̄ = vec(mean(driver.forcing.qt_nud, dims=2)); ō = vec(mean(driver.forcing.o3, dims=2))
    faces = vcat(z_top, [(zc[k] + zc[k+1]) / 2 for k in 1:length(zc)-1], zc[end] + (zc[end] - zc[end-1]) / 2)
    logp(ζ) = interpolate_profile(driver.z, log.(p̄), ζ)
    p_faces = exp.(logp.(faces))
    T_faces = interpolate_profile(driver.z, T̄, faces)
    return (; z = zc, p = p̄[sel], T = T̄[sel], q = max.(q̄[sel], 0), o3 = ō[sel], faces, p_faces, T_faces)
end

"""
    extended_column_radiation(grid, constants; driver, surface_temperature, kwargs...)

Breeze's all-sky `RadiativeTransferModel` with the driver's time-mean atmosphere above the LES
top appended to every RRTMGP column (see the file header). Keyword arguments as for
`RadiativeTransferModel(grid, AllSkyOptics(), constants; ...)`.
"""
function extended_column_radiation(grid, constants; driver,
                                   surface_temperature,
                                   surface_albedo = 0.07,
                                   surface_emissivity = 0.98,
                                   background_atmosphere = BackgroundAtmosphere(),
                                   solar_position = FixedCosineZenith(0),
                                   solar_constant = 1361,
                                   schedule = IterationInterval(1),
                                   liquid_effective_radius = ConstantRadiusParticles(10e-6),
                                   ice_effective_radius = ConstantRadiusParticles(30e-6),
                                   ice_roughness = 2)
    ext = rrtmgp_extension()
    FT = eltype(grid)
    parameters = ext.RRTMGPParameters(constants)
    solar_position = ext.maybe_infer_solar_position(solar_position, grid)
    background_atmosphere = materialize_background_atmosphere(background_atmosphere, grid)
    surface_albedo = ext.materialize_surface_property(surface_albedo, grid, solar_position)
    surface_emissivity = ext.materialize_surface_property(surface_emissivity, grid, solar_position)

    arch = architecture(grid)
    Nx, Ny, Nz = size(grid)
    Nc = Nx * Ny
    z_top = Array(znodes(grid, Face()))[Nz + 1]
    upper = upper_atmosphere_layers(driver, z_top)
    Nₑ = length(upper.z)
    nlay = Nz + Nₑ

    context = ext.rrtmgp_context(arch)
    ArrayType = ClimaComms.array_type(context.device)
    grid_params = RRTMGPGridParams(FT; context, domain_nlay = nlay, ncol = Nc)
    radiation_method = AllSkyRadiation(false, false)
    luts = lookup_tables(grid_params, radiation_method)
    Nband_lw = luts.nbnd_lw; Nband_sw = luts.nbnd_sw; Ngas = luts.ngas_sw

    # Host arrays for the full column, upper layers filled from the driver, then moved to the device.
    g = constants.gravitational_acceleration
    mᵈ = constants.dry_air.molar_mass
    mᵛ = constants.vapor.molar_mass
    layerdata = zeros(FT, 4, nlay, Nc)
    pᶠ = zeros(FT, nlay + 1, Nc)
    Tᶠ = zeros(FT, nlay + 1, Nc)
    vmr_h2o = zeros(FT, nlay, Nc)
    vmr_o3 = zeros(FT, nlay, Nc)
    for m in 1:Nₑ
        k = Nz + m
        Δp = upper.p_faces[m] - upper.p_faces[m + 1]
        Δp > 0 || throw(ArgumentError("upper-atmosphere pressure is not decreasing with height at layer $m"))
        q = upper.q[m]
        column_dry = (Δp / g) * (1 - q) / mᵈ * AVOGADRO / 1e4
        Tlayer = (upper.T_faces[m] + upper.T_faces[m + 1]) / 2
        layerdata[1, k, :] .= column_dry
        layerdata[2, k, :] .= upper.p[m]
        layerdata[3, k, :] .= clamp(Tlayer, 160, 355)
        layerdata[4, k, :] .= 0
        pᶠ[k + 1, :] .= upper.p_faces[m + 1]
        Tᶠ[k + 1, :] .= clamp(upper.T_faces[m + 1], 160, 355)
        vmr_h2o[k, :] .= q / (1 - q) * (mᵈ / mᵛ)
        vmr_o3[k, :] .= upper.o3[m]
    end
    rrtmgp_λ = ArrayType{FT}(undef, Nc); rrtmgp_φ = ArrayType{FT}(undef, Nc)
    ext.set_longitude!(rrtmgp_λ, solar_position, grid); ext.set_latitude!(rrtmgp_φ, solar_position, grid)
    rrtmgp_layerdata = ArrayType(layerdata); rrtmgp_pᶠ = ArrayType(pᶠ); rrtmgp_Tᶠ = ArrayType(Tᶠ)
    rrtmgp_T₀ = ArrayType{FT}(undef, Nc)
    vmr = ext.initialize_global_mean_vmr(Ngas, nlay, Nc, FT, ArrayType)
    ext.set_global_mean_gases!(vmr, luts.idx_gases_sw, background_atmosphere)
    copyto!(vmr.vmr_h2o, ArrayType(vmr_h2o)); copyto!(vmr.vmr_o3, ArrayType(vmr_o3))

    zero_layers() = (a = ArrayType{FT}(undef, nlay, Nc); fill!(a, zero(FT)); a)
    false_layers() = (a = ArrayType{Bool}(undef, nlay, Nc); fill!(a, false); a)
    cloud_state = CloudState(zero_layers(), zero_layers(), zero_layers(), zero_layers(), zero_layers(),
                             false_layers(), false_layers(), MaxRandomOverlap(), ice_roughness)
    atmospheric_state = AtmosphericState(rrtmgp_λ, rrtmgp_φ, rrtmgp_layerdata, rrtmgp_pᶠ, rrtmgp_Tᶠ, rrtmgp_T₀, vmr, cloud_state, nothing)

    cos_zenith = ArrayType{FT}(undef, Nc); ext.initialize_cos_zenith!(cos_zenith, solar_position)
    rrtmgp_ℐ₀ = ArrayType{FT}(undef, Nc); rrtmgp_ℐ₀ .= convert(FT, solar_constant)
    rrtmgp_ε₀ = ArrayType{FT}(undef, Nband_lw, Nc); rrtmgp_αb₀ = ArrayType{FT}(undef, Nband_sw, Nc); rrtmgp_αw₀ = ArrayType{FT}(undef, Nband_sw, Nc)
    surface_emissivity = ext.constant_field_property(surface_emissivity, FT)
    surface_albedo = ext.constant_field_property(surface_albedo, FT)
    if surface_temperature isa Number
        surface_temperature = ConstantField(convert(FT, surface_temperature))
        rrtmgp_T₀ .= surface_temperature.constant
    end
    lw_bcs = LwBCs(rrtmgp_ε₀, nothing)
    sw_bcs = SwBCs(cos_zenith, rrtmgp_ℐ₀, rrtmgp_αb₀, nothing, rrtmgp_αw₀)
    solver = RRTMGPSolver(grid_params, radiation_method, parameters, lw_bcs, sw_bcs, atmospheric_state)

    surface_radiation = SurfaceRadiation(surface_temperature, surface_emissivity, surface_albedo, surface_albedo)
    ext.update_rrtmgp_surface_boundary_conditions!(solver, surface_radiation, grid)
    radius(r) = r isa ConstantRadiusParticles ? ConstantRadiusParticles(convert(FT, r.radius)) : r
    # Breeze 0.12 (PR 998) added `column_batches` before `schedule`; `nothing` is its one unbatched solve,
    # which keeps the solver's own (extended) state. Breeze 0.11 (cases/tracer_mip) has no such field.
    column_batches = hasfield(RadiativeTransferModel, :column_batches) ? (nothing,) : ()
    rtm = RadiativeTransferModel(convert(FT, solar_constant), solar_position, surface_radiation, background_atmosphere,
                                 atmospheric_state, solver, nothing,
                                 ZFaceField(grid), ZFaceField(grid), ZFaceField(grid), ZFaceField(grid), CenterField(grid),
                                 radius(liquid_effective_radius), radius(ice_effective_radius), column_batches..., schedule)
    return rtm, (; layers_above = Nₑ, z_top, column_top = upper.faces[end], p_top = upper.p_faces[end])
end

#####
##### Solar position along the trajectory
#####

"""
    TrajectorySolarPosition(radiation, epoch, times, longitude, latitude)

Callback (`IterationInterval(1)`) writing `max(cos θ_z, 0)` for the composite trajectory's
mean position at the current time into the RRTMGP shortwave boundary condition, using Breeze's
`cos_solar_zenith_angle` (Spencer 1971). The radiation model must have been built with
`solar_position = FixedCosineZenith(...)` so that Breeze's own update is a no-op.
"""
struct TrajectorySolarPosition{R, T, L}
    radiation :: R
    epoch :: DateTime
    times :: T
    longitude :: L
    latitude :: L
end

function TrajectorySolarPosition(radiation, epoch, times, longitude, latitude)
    radiation.solar_position isa FixedCosineZenith ||
        throw(ArgumentError("TrajectorySolarPosition needs a radiation model built with solar_position = FixedCosineZenith(...)"))
    return TrajectorySolarPosition{typeof(radiation), typeof(times), typeof(longitude)}(radiation, epoch, times, longitude, latitude)
end

function trajectory_cos_zenith(u::TrajectorySolarPosition, t)
    datetime = u.epoch + Millisecond(round(Int, 1000t))
    lon = interpolate_profile(u.times, u.longitude, t)
    lat = interpolate_profile(u.times, u.latitude, t)
    return max(cos_solar_zenith_angle(datetime, lon, lat), 0), lon, lat
end

function (u::TrajectorySolarPosition)(simulation_or_model)
    model = simulation_or_model isa Simulation ? simulation_or_model.model : simulation_or_model
    cosθ, _, _ = trajectory_cos_zenith(u, model.clock.time)
    u.radiation.longwave_solver.sws.bcs.cos_zenith .= convert(eltype(u.radiation.solar_constant), cosθ)
    return nothing
end

"""
    composite_trajectory_path(path)

Hourly mean (longitude, latitude) of the GEOS-5 trajectories in `SEA_STARR_Raw_Trajectories.nc`
(missing values ignored): `(; times [s], longitude, latitude)`.
"""
function composite_trajectory_path(path)
    ds = NCDataset(path)
    lon = Float64.(coalesce.(ds["lon"][:, :], NaN)); lat = Float64.(coalesce.(ds["lat"][:, :], NaN))
    hours = Float64.(ds["t"][:])
    close(ds)
    finite_mean(v) = mean(filter(isfinite, v))
    return (; times = hours .* 3600, longitude = [finite_mean(lon[n, :]) for n in axes(lon, 1)], latitude = [finite_mean(lat[n, :]) for n in axes(lat, 1)])
end
