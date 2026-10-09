#####
##### ARM VARANAL / E3SM "iopfile_4scam" forcing: reader, and adapter to the SAM-format
##### records that `tracer_dp_scream` consumes.
#####
##### The TRACER file (`TRACER_iopfile_4scam.nc`, datastream hou60v1varanaecmwfX12.c1,
##### doi 10.5439/1860369) is hourly, on 40 pressure levels (50–1025 hPa, ascending in the
##### file, i.e. model top first), with ERA5-based, precipitation-constrained large-scale
##### tendencies. DP-SCREAM/EAMxx (`intensive_observation_period.cpp`, E3SM 1d551ea2b0)
##### reads `T, q, u, v, Ps, Tg, shflx, lhflx, divT, vertdivT, divq, vertdivq, omega`:
#####
#####   * with `vertdivT`/`vertdivq` present it forms `divT3d = divT + vertdivT` and
#####     `divq3d = divq + vertdivq` and adds *only* those to T and qv (`advance_iop_forcing`);
#####     omega-based subsidence is a separate option (`iop_dosubsidence`, false for TRACER);
#####   * data are read at the IOP record whose interval contains the model time
#####     (piecewise constant in time, `get_iop_file_time_idx`);
#####   * `iop_srf_prop = true` overwrites the surface sensible heat flux, evaporation
#####     (lhflx / Lv) and the radiative surface temperature (Tg) with the file values.
#####
##### `iop_sam_inputs` maps this onto the SAM record structs so the forcing operators of
##### `large_scale_forcings.jl` and the surface machinery of `surface_fluxes.jl` are reused
##### unchanged: tls = divT3d, qls = divq3d (mass-fraction basis), uls/vls = u/v (nudging
##### targets, if nudging is enabled), wls = 0 (no separate vertical transport), sfc = (Tg,
##### shflx, lhflx). The piecewise-constant time dependence is reproduced by duplicating each
##### record just before the next one (`hold_gap` seconds earlier).
#####

using NCDatasets: NCDataset, nomissing
using Dates: Dates, DateTime
using Breeze: ThermodynamicConstants
using Breeze.Thermodynamics: dry_air_gas_constant

"""
    IOPForcing

The contents of an E3SM/SCAM IOP forcing file. `times` are seconds since `base_time`
(the file's `bdate`), `levels` are the file's pressure levels [Pa] in file order, `profiles`
holds `nlev × ntime` matrices and `surface` holds `ntime` vectors, both keyed by the file
variable names. Fill values become `NaN`.
"""
struct IOPForcing{FT}
    path :: String
    base_time :: DateTime
    times :: Vector{Float64}
    levels :: Vector{FT}
    latitude :: FT
    longitude :: FT
    profiles :: NamedTuple
    surface :: NamedTuple
    attributes :: Dict{String, String}
end

Base.summary(iop::IOPForcing) =
    string("IOPForcing(", basename(iop.path), ": ", length(iop.times), " records from ", iop.base_time,
           ", ", length(iop.levels), " levels, ", join(string.(keys(iop.profiles)), " "), ")")
Base.show(io::IO, iop::IOPForcing) = print(io, summary(iop))

const IOP_PROFILE_VARIABLES = (:T, :q, :u, :v, :omega, :divT, :vertdivT, :divq, :vertdivq)
const IOP_SURFACE_VARIABLES = (:Ps, :Tg, :Tsair, :shflx, :lhflx, :Prec, :srfswdn, :srfswup, :srflwdn, :srflwup, :windsrf)

"""
    read_iop_forcing(path; FT=Float64, profile_variables=IOP_PROFILE_VARIABLES,
                     surface_variables=IOP_SURFACE_VARIABLES)

Read an IOP forcing NetCDF file (`bdate`, `tsec`, `lev`, `lat`, `lon`, profile variables
dimensioned `(time, lev, lat, lon)` and surface variables `(time, lat, lon)`).
"""
function read_iop_forcing(path; FT=Float64, profile_variables=IOP_PROFILE_VARIABLES,
                          surface_variables=IOP_SURFACE_VARIABLES)
    NCDataset(path) do ds
        bdate = Int(ds["bdate"][])
        base_time = DateTime(bdate ÷ 10000, (bdate ÷ 100) % 100, bdate % 100)
        times = Float64.(nomissing(vec(Array(ds["tsec"])), NaN))
        issorted(times) && allunique(times) || error("IOP times must be strictly increasing in $path")
        levels = FT.(nomissing(vec(Array(ds["lev"])), NaN))
        latitude = FT(ds["lat"][1])
        longitude = FT(ds["lon"][1])
        read_profile(name) = FT.(nomissing(ds[String(name)][1, 1, :, :], NaN))
        read_surface(name) = FT.(vec(nomissing(ds[String(name)][1, 1, :], NaN)))
        profiles = NamedTuple{profile_variables}(map(read_profile, profile_variables))
        surface = NamedTuple{surface_variables}(map(read_surface, surface_variables))
        attributes = Dict{String, String}(string(k) => string(v) for (k, v) in ds.attrib)
        return IOPForcing(abspath(path), base_time, times, levels, latitude, longitude, profiles, surface, attributes)
    end
end

"""
    iop_datetimes(iop)

UTC `DateTime` of every record (`base_time + tsec`).
"""
iop_datetimes(iop::IOPForcing) = [iop.base_time + Dates.Millisecond(round(Int, 1000t)) for t in iop.times]

"""
    iop_index(iop, t::DateTime)

Index of the record at exactly `t` (within one second); IOP windows start and end on
records, as EAMxx requires the model time to lie inside the file's time span.
"""
function iop_index(iop::IOPForcing, t::DateTime)
    seconds = Dates.value(t - iop.base_time) / 1000
    i = findfirst(τ -> abs(τ - seconds) ≤ 1, iop.times)
    isnothing(i) && throw(ArgumentError("$t is not a record of $(basename(iop.path)) (records every $(iop.times[2] - iop.times[1]) s from $(iop.base_time))"))
    return i
end

"""
    fractional_day_of_year(t::DateTime)

SAM `day` convention: 1.0 at 00 UTC on 1 January of `t`'s year.
"""
fractional_day_of_year(t::DateTime) = 1 + Dates.value(t - DateTime(Dates.year(t), 1, 1)) / 86_400_000

"""
    iop_potential_temperature(T, p; constants=ThermodynamicConstants(Float64), standard_pressure=1e5)

Dry potential temperature `T (p₀/p)^(Rᵈ/cᵖᵈ)` with Breeze's constants, the inverse of the
conversion `initial_state_columns` applies when it recovers `T` from `θ` at the reference pressure.
"""
function iop_potential_temperature(T, p; constants=ThermodynamicConstants(Float64), standard_pressure=1e5)
    κ = dry_air_gas_constant(constants) / constants.dry_air.heat_capacity
    return T * (standard_pressure / p)^κ
end

# Duplicate every record (except the last) `gap` seconds before the next record so that
# linear time interpolation reproduces EAMxx's piecewise-constant reading of the IOP file.
function hold_records(records::Vector{SAMLargeScaleForcing{FT}}, gap) where FT
    held = SAMLargeScaleForcing{FT}[]
    for (n, r) in enumerate(records)
        push!(held, r)
        n == length(records) && break
        day′ = records[n + 1].day - gap / 86400
        day′ > r.day || throw(ArgumentError("hold gap $gap s is not smaller than the record spacing"))
        push!(held, SAMLargeScaleForcing(FT(day′), r.surface_pressure, r.z, r.p, r.tls, r.qls, r.uls, r.vls, r.wls, r.ug, r.vg, r.has_geostrophic_columns))
    end
    return held
end

function hold_surface(sfc::SAMSurfaceForcing{FT}, gap) where FT
    n = length(sfc.day)
    idx = Int[]; day = FT[]
    for i in 1:n
        push!(idx, i); push!(day, sfc.day[i])
        i == n && break
        push!(idx, i); push!(day, sfc.day[i + 1] - gap / 86400)
    end
    return SAMSurfaceForcing(day, sfc.sst[idx], sfc.sensible_heat_flux[idx], sfc.latent_heat_flux[idx], sfc.kinematic_stress[idx])
end

"""
    iop_sam_inputs(iop; start, stop, transport=:iop_3d, hold=true, hold_gap=1.0,
                   constants=ThermodynamicConstants(Float64))

The IOP records from `start` to `stop` (UTC `DateTime`s on records) as the input bundle
of the SAM-forced cases: one `SAMSounding` at `start` (pressure levels, `θ` from `T`, `q`
kept as the file's specific humidity ⇒ use `moisture_basis = :mass_fraction`), the
`SAMLargeScaleForcing` records with

- `transport = :iop_3d` (DP-SCREAM with `vertdivT`/`vertdivq` present): `tls = divT + vertdivT`,
  `qls = divq + vertdivq`, `wls = 0`;
- `transport = :iop_horizontal` (sensitivity): `tls = divT`, `qls = divq`, `wls = 0`;

`uls/vls = u/v` (nudging targets), `ug/vg` aliased to them (no geostrophic columns), and the
surface series `sfc = (day, Tg, shflx, lhflx, 0)`. With `hold = true` every record is
duplicated `hold_gap` seconds before the next one so the linearly interpolating forcing
operators reproduce EAMxx's piecewise-constant hourly forcing. Returns `(; soundings, lsf,
sfc, namelist, paths, day0, epoch, start, stop, window, transport, hold)`.
"""
function iop_sam_inputs(iop::IOPForcing{FT}; start::DateTime, stop::DateTime, transport=:iop_3d,
                        hold=true, hold_gap=1.0, constants=ThermodynamicConstants(Float64)) where FT
    transport ∈ (:iop_3d, :iop_horizontal) || throw(ArgumentError("transport must be :iop_3d or :iop_horizontal, got $transport"))
    i0 = iop_index(iop, start)
    i1 = iop_index(iop, stop)
    i1 > i0 || throw(ArgumentError("stop must be after start"))
    dt = iop_datetimes(iop)
    order = sortperm(iop.levels; rev=true)           # SAM order: highest pressure (bottom) first
    p = Float64.(iop.levels[order])
    P = iop.profiles
    missing_profiles = [name for name in (:T, :q, :u, :v, :divT, :divq) if !haskey(P, name)]
    isempty(missing_profiles) || throw(ArgumentError("IOP file lacks $(missing_profiles)"))
    if transport === :iop_3d
        haskey(P, :vertdivT) && haskey(P, :vertdivq) ||
            throw(ArgumentError("transport = :iop_3d needs vertdivT and vertdivq (EAMxx builds divT3d/divq3d from them)"))
    end
    column(name, i) = Float64.(P[name][order, i])
    nan = fill(NaN, length(p))
    day(i) = fractional_day_of_year(dt[i])
    S = iop.surface
    T₀ = column(:T, i0)
    θ₀ = [iop_potential_temperature(T₀[k], p[k]; constants) for k in eachindex(p)]
    sounding = SAMSounding(day(i0), Float64(S.Ps[i0]), copy(nan), p, θ₀, column(:q, i0), column(:u, i0), column(:v, i0))
    lsf = map(i0:i1) do i
        tls = transport === :iop_3d ? column(:divT, i) .+ column(:vertdivT, i) : column(:divT, i)
        qls = transport === :iop_3d ? column(:divq, i) .+ column(:vertdivq, i) : column(:divq, i)
        u = column(:u, i); v = column(:v, i)
        SAMLargeScaleForcing(day(i), Float64(S.Ps[i]), copy(nan), p, tls, qls, u, v, zeros(length(p)), copy(u), copy(v), false)
    end
    days = [day(i) for i in i0:i1]
    sfc = SAMSurfaceForcing(days, Float64.(S.Tg[i0:i1]), Float64.(S.shflx[i0:i1]), Float64.(S.lhflx[i0:i1]), zeros(length(days)))
    if hold
        lsf = hold_records(lsf, hold_gap)
        sfc = hold_surface(sfc, hold_gap)
    end
    all(r -> all(isfinite, r.tls) && all(isfinite, r.qls) && all(isfinite, r.uls) && all(isfinite, r.vls), lsf) ||
        error("non-finite forcing values in the selected IOP window")
    all(isfinite, θ₀) && all(isfinite, sounding.q) || error("non-finite initial profile at $start")
    namelist = Dict{String, NamelistValue}()
    return (; soundings = [sounding], lsf, sfc, namelist, paths = (; iop = iop.path),
              day0 = day(i0), epoch = dt[i0], start, stop, window = (i0, i1), transport, hold)
end

"""
    iop_surface_albedo(iop; start, stop, minimum_downwelling=10)

Energy-weighted surface shortwave albedo over the window, `Σ srfswup / Σ srfswdn` for records
with downwelling SW above `minimum_downwelling` W/m². The TRACER file labels `srfswup` as
"net" but its values are the upwelling flux (ratio ≈ 0.15, versus 0.85 for the net reading);
the same holds for `srflwup` (≈ σTg⁴).
"""
function iop_surface_albedo(iop::IOPForcing; start::DateTime, stop::DateTime, minimum_downwelling=10)
    i0 = iop_index(iop, start); i1 = iop_index(iop, stop)
    dn = iop.surface.srfswdn[i0:i1]; up = iop.surface.srfswup[i0:i1]
    lit = findall(d -> isfinite(d) && d > minimum_downwelling, dn)
    isempty(lit) && error("no daylit records between $start and $stop")
    return sum(up[lit]) / sum(dn[lit])
end

"""
    iop_column_integral(iop, name, i; g=9.81)

`∫ ϕ dp / g` of profile variable `name` at record `i` over the levels at or above the surface
pressure `Ps` (kg m⁻² × the variable's unit, e.g. kg m⁻² s⁻¹ for a mixing-ratio tendency),
using the trapezoidal rule on the file's pressure levels with the surface value taken from the
lowest level above the surface, as the file's note prescribes.
"""
function iop_column_integral(iop::IOPForcing, name::Symbol, i::Integer; g=9.81)
    ps = Float64(iop.surface.Ps[i])
    order = sortperm(iop.levels)                 # ascending pressure (top → bottom)
    p = Float64.(iop.levels[order]); ϕ = Float64.(iop.profiles[name][order, i])
    inside = findall(≤(ps), p)
    isempty(inside) && return 0.0
    pp = vcat(p[inside], ps); ϕϕ = vcat(ϕ[inside], ϕ[inside[end]])
    return sum((pp[k+1] - pp[k]) * (ϕϕ[k] + ϕϕ[k+1]) / 2 for k in 1:length(pp)-1) / g
end
