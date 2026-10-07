#####
##### TRACER–DP-SCREAM: a doubly periodic Breeze case driven by the forcing of the published
##### DP-SCREAM TRACER experiment (Oware et al. 2025, doi 10.1029/2025JD044113; outputs
##### 10.5281/zenodo.15271730; run script scmlib DPxx_SCREAM_SCRIPTS/run_dpxx_scream_TRACER.csh).
#####
##### Protocol facts taken from the script, the EAMxx source at the archived commit
##### (git_version 1d551ea2b0, 2024-03-12) and the IOP file:
#####   * forcing file TRACER_iopfile_4scam.nc (ARM VARANAL ERA5-constrained, hourly, 40 levels,
#####     150-km domain around 29.75 N, 95.45 W); DP-SCREAM 200 km × 200 km at 3.33 km, 128 levels,
#####     dt = 100 s, RRTMGP every 300 s, SHOC + P3, land/ocean coupler stress;
#####   * large-scale transport = the file's 3-D advective tendencies (divT + vertdivT,
#####     divq + vertdivq) applied once; no omega subsidence (`iop_dosubsidence = false`);
#####   * `iop_srf_prop = true`: sensible/latent heat fluxes and the radiative surface
#####     temperature (Tg) prescribed from the file;
#####   * planar HOMME has `fcor = 0` and `iop_coriolis = false`: no Coriolis force;
#####   * `iop_nudge_uv = true` in the script, but the archived E3SM revision predates the
#####     DP-EAMxx wind-nudging port (E3SM 3f7eee0053, 2024-04-16): the archived runs applied
#####     **no** wind nudging. `wind_nudging_timescale = nothing` reproduces that; 10800 s is
#####     the script's intended (EAMxx default) relaxation.
#####
##### Deliberate differences of this Breeze case (recorded in the label and provenance):
#####   LES closure (Smagorinsky–Lilly) at Δx = 200 m on a 51.2 km domain instead of SHOC at
#####   3.33 km on 200 km; 160 levels to 22 km (DP-SCREAM: 128 hybrid levels to ~40 km) with a
#####   Rayleigh sponge above ~16.5 km; P3 with a prescribed droplet number (SCREAM P3: prognostic
#####   droplet number with prescribed CCN); constant bulk drag with a land roughness length
#####   instead of the ELM/ocean coupler stress; anelastic dynamics; no aerosol radiative effects.
#####

using Dates: Dates, DateTime
using Oceananigans.Units: minutes

"""
    tracer_dp_scream_vertical_faces(; Δz=50, uniform_top=2000, top=22000, Nz=160)

Cell interfaces: uniform `Δz` to `uniform_top`, then a constant-ratio stretch to `top`
(50 m → ≈ 400 m at 22 km with the defaults). Deep-convective LES grid; the DP-SCREAM
128-level hybrid grid is not reproduced.
"""
tracer_dp_scream_vertical_faces(; Δz=50.0, uniform_top=2000.0, top=22000.0, Nz=160) =
    uniform_then_stretched_faces(; Nz, top, Δz, uniform_top)

const TRACER_DP_SCREAM_DATA = joinpath(@__DIR__, "..", "data", "tracer_dp_scream")

"""
    neutral_drag_coefficient(z₁, roughness_length; κ=0.4)

`(κ / ln(z₁ / z₀))²`, the neutral log-law drag coefficient referenced to the first cell
center height `z₁`.
"""
neutral_drag_coefficient(z₁, roughness_length; κ=0.4) = (κ / log(z₁ / roughness_length))^2

"""
    tracer_dp_scream(; iop_path, start=DateTime(2022, 8, 1), stop=DateTime(2022, 8, 15), kwargs...)

Build the periodic TRACER case driven by the DP-SCREAM IOP forcing without running it.
Returns the `build_case` bundle plus `iop`, `start`, `stop`, `window` and protocol records.

Defaults (all overridable, overrides are recorded in the provenance):

- `start`/`stop`: UTC `DateTime`s on IOP records. The DP-SCREAM August run started at
  2022-08-01 00 UTC and the paper's 5–15 August window is its days 5–15, so the default
  reproduces that initialization and continuous forcing history (14 days); the archived
  0.5-km run started cold at 2022-08-05 00 UTC.
- `arch = GPU()`, `FT = Float32`, `Nx = Ny = 256`, `Δx = Δy = 200` m, `z_faces = tracer_dp_scream_vertical_faces()`
- `microphysics = :p3_n75` with `droplet_number = 200e6` m⁻³ (continental assumption, see above)
- `large_scale_transport = :iop_3d`, `hold_forcing = true` (piecewise-constant hourly forcing)
- `wind_nudging_timescale = nothing` (archived runs), `coriolis = nothing`
- `radiation = :rrtmgp` every `radiation_interval = 300` s with the IOP `Tg` as surface
  temperature, `surface_albedo` from the IOP shortwave fluxes (≈ 0.15), `surface_emissivity = 0.98`
- `surface = :prescribed_heat_fluxes_bulk_drag`: IOP sensible/latent heat fluxes, neutral
  log-law drag for `roughness_length = 0.1` m referenced to the first cell center, `gustiness = 1` m/s
- `sponge = SAMSponge(damping_depth_fraction = 0.25)`, `Δt = 2`, `max_Δt = 5`, `cfl = 0.7`
- outputs: profiles every 30 min, time series every 60 s, slices every 30 min at 1500 m,
  `checkpoint_interval = nothing` (set e.g. `6hours` for resumable runs)
"""
function tracer_dp_scream(;
    iop_path = joinpath(TRACER_DP_SCREAM_DATA, "TRACER_iopfile_4scam.nc"),
    start = DateTime(2022, 8, 1),
    stop = DateTime(2022, 8, 15),
    arch = GPU(),
    FT = Float32,
    Nx = 256, Ny = 256, Δx = 200.0, Δy = Δx,
    z_faces = tracer_dp_scream_vertical_faces(),
    microphysics = :p3_n75,
    droplet_number = 200e6,
    large_scale_transport = :iop_3d,
    hold_forcing = true,
    wind_nudging_timescale = nothing,
    coriolis = nothing,
    radiation = :rrtmgp,
    radiation_interval = 300.0,
    surface_albedo = nothing,
    surface_emissivity = 0.98,
    roughness_length = 0.1,
    drag_coefficient = nothing,
    gustiness = 1.0,
    sponge = SAMSponge(damping_depth_fraction=0.25),
    Δt = 2.0, max_Δt = 5.0, cfl = 0.7,
    stop_time = nothing,
    output_dir = joinpath(pwd(), "output", "tracer_dp_scream"),
    output_prefix = "tracer",
    profile_interval = 30minutes,
    timeseries_interval = 60.0,
    slice_interval = 30minutes,
    slice_height = 1500.0,
    checkpoint_interval = nothing,
    kwargs...)

    isfile(iop_path) || throw(ArgumentError("IOP forcing file not found: $iop_path (fetch it with `julia data_wrangling/fetch_manifest.jl cases/tracer_dp_scream/inputs.toml data/tracer_dp_scream`)"))
    iop = read_iop_forcing(iop_path)
    inputs = iop_sam_inputs(iop; start, stop, transport=large_scale_transport, hold=hold_forcing)
    albedo = isnothing(surface_albedo) ? iop_surface_albedo(iop; start, stop) : surface_albedo
    z₁ = z_faces[2] / 2
    Cᴰ = isnothing(drag_coefficient) ? neutral_drag_coefficient(z₁, roughness_length) : drag_coefficient
    duration = Dates.value(stop - start) / 1000
    stop_time = isnothing(stop_time) ? duration : stop_time
    stop_time ≤ duration || throw(ArgumentError("stop_time = $stop_time s exceeds the forcing window ($duration s)"))

    label = string("TRACER DP-SCREAM-matched periodic case: IOP ", get(iop.attributes, "datastream", "?"),
                   " from ", start, " to ", stop, " UTC; transport ", large_scale_transport, " applied once, no subsidence, ",
                   "no Coriolis, wind nudging ", isnothing(wind_nudging_timescale) ? "off (archived E3SM 1d551ea2b0)" : "τ = $wind_nudging_timescale s",
                   "; IOP H/LE + Tg, bulk drag Cᴰ = ", round(Cᴰ; sigdigits=3), " (z₀ = $roughness_length m); ",
                   "LES Smagorinsky at Δx = $Δx m on $(Nx * Δx / 1000) km (DP-SCREAM: SHOC at 3.33 km on 200 km, 128 levels, ",
                   "prognostic droplet number) — intentional physics/resolution differences, not a DP-SCREAM reproduction")

    case = build_case(inputs; arch, FT, Nx, Ny, Lx=Nx * Δx, Ly=Ny * Δy, z_faces,
                      day0=inputs.day0, epoch=inputs.epoch, moisture_basis=:mass_fraction,
                      latitude=Float64(iop.latitude), longitude=Float64(iop.longitude),
                      microphysics, droplet_number, radiation, radiation_interval,
                      surface_albedo=albedo, surface_emissivity,
                      surface=:prescribed_heat_fluxes_bulk_drag, drag_coefficient=Cᴰ, gustiness,
                      wind_nudging_timescale, translation_velocity=(0.0, 0.0),
                      vertical_advection=nothing, geostrophic=false, coriolis, thermodynamic_tendencies=true,
                      upper_boundary_relaxation=false, sponge,
                      stop_time, Δt, max_Δt, cfl, label, output_dir, output_prefix,
                      profile_interval, timeseries_interval, slice_interval, slice_height,
                      checkpoint_interval, kwargs...)

    protocol = (; large_scale_transport, hold_forcing, roughness_length, surface_albedo_source = isnothing(surface_albedo) ? "IOP srfswup/srfswdn" : "override",
                  iop_datastream = get(iop.attributes, "datastream", ""), iop_doi = get(iop.attributes, "doi", ""),
                  dp_scream_git_version = "1d551ea2b0", scmlib_revision = "2dc3f1073a5d03b5f32617e79d68fb20a63cfb43",
                  zenodo_record = "10.5281/zenodo.15271730", start = string(start), stop = string(stop),
                  window = inputs.window, forcing_records = length(inputs.lsf))
    return merge(case, (; preset = :tracer_dp_scream, protocol = :tracer_dp_scream, protocol_dimensions = nothing,
                          protocol_overrides = merge(protocol, NamedTuple(kwargs)), iop, start, stop, window = inputs.window,
                          iop_inputs = inputs))
end

"""
    tracer_dp_scream_settings(iop_path; start, stop)

Inspect the forcing window without allocating a model: records, IOP-implied albedo, and
the initial surface state.
"""
function tracer_dp_scream_settings(iop_path = joinpath(TRACER_DP_SCREAM_DATA, "TRACER_iopfile_4scam.nc");
                                   start = DateTime(2022, 8, 1), stop = DateTime(2022, 8, 15))
    iop = read_iop_forcing(iop_path)
    i0 = iop_index(iop, start); i1 = iop_index(iop, stop)
    return (; iop, start, stop, window = (i0, i1), records = i1 - i0 + 1,
              duration_seconds = Dates.value(stop - start) / 1000,
              surface_albedo = iop_surface_albedo(iop; start, stop),
              initial_surface_pressure = iop.surface.Ps[i0], initial_ground_temperature = iop.surface.Tg[i0],
              latitude = iop.latitude, longitude = iop.longitude,
              datastream = get(iop.attributes, "datastream", ""), doi = get(iop.attributes, "doi", ""))
end
