#####
##### DEPHY SCM-format driver reader (format version 1) for the SEA STARR drivers.
#####
##### The SEA STARR drivers (`SEA_STARR_{CTRL,N100,N030}_SCM_driver.nc`) store 336 mid-levels
##### `zh` TOP-DOWN (76 km → 5 m), 67 hourly forcing records, one initial record (`t0`), a
##### single trajectory position (`lat`, `lon`) and the aerosol number in "mg-1" (per milligram
##### of air). The reader reverses the levels to ascending order, converts the aerosol number
##### to kg⁻¹ and keeps every field as a plain host array so that the case constructor can
##### interpolate onto the LES grid. See docs/cases/sea_starr.md for the variable inventory.
#####

using Dates: Dates, DateTime
using NCDatasets: NCDataset, dimnames
using Statistics: mean

"""
    DEPHYDriver

In-memory copy of a DEPHY SCM driver file, levels ascending. Fields:

- `path`, `case`, `attributes`: file provenance and the global attributes
- `start_time`, `end_time`: UTC period of the forcing (`time` units and `end_date`)
- `latitude`, `longitude`: the single trajectory position of the file
- `z`: mid-level heights [m], ascending
- `times`: forcing times [s since `start_time`]
- `surface_pressure` [Pa], `surface_temperature` [K], `roughness_length` [m]: `ps`, `ts`, `z0`
- `initial`: profiles at `t0` — `thetal` [K], `qt` [kg/kg, mass fraction], `na` [kg⁻¹],
  `u`, `v` [m/s], `p` [Pa], `T` [K]
- `forcing`: `(Nz, Nt)` matrices — `w` [m/s], `omega` [Pa/s], `ug`, `vg` [m/s], `thetal_nud`
  [K], `T_nud` [K], `qt_nud` [kg/kg], `na_nud` [kg⁻¹], `u_nud`, `v_nud` [m/s], `o3` [mol/mol], `p` [Pa]
- `nudging_rates`: the (uniform) `nudging_constant_*` values [s⁻¹] for `thetal`, `qt`, `na`, `u`, `v`
- `sst` [K], `ps` [Pa]: `ts_force`, `ps_force` time series
"""
struct DEPHYDriver{FT}
    path :: String
    case :: String
    attributes :: Dict{String, Any}
    start_time :: DateTime
    end_time :: DateTime
    latitude :: FT
    longitude :: FT
    z :: Vector{FT}
    times :: Vector{FT}
    surface_pressure :: FT
    surface_temperature :: FT
    roughness_length :: FT
    initial :: NamedTuple
    forcing :: NamedTuple
    nudging_rates :: NamedTuple
    sst :: Vector{FT}
    ps :: Vector{FT}
end

Base.summary(d::DEPHYDriver) =
    string("DEPHYDriver(", d.case, ": ", length(d.z), " levels ", first(d.z), "–", last(d.z), " m, ",
           length(d.times), " records ", d.start_time, " → ", d.end_time, ")")
Base.show(io::IO, d::DEPHYDriver) = print(io, summary(d))

"""
    aerosol_number_per_kg(na_per_mg)

Convert the driver's aerosol number per milligram of air (`units = "mg-1"`) to a number
mixing ratio per kilogram: `nᵃ = 10⁶ na`.
"""
aerosol_number_per_kg(na_per_mg) = 1e6 * na_per_mg

# Read a variable, replace `missing` by NaN, and convert to FT.
function read_values(ds, name, FT)
    v = Array(ds[name])          # full array with its dimensions (`[:]` flattens 2D variables)
    return FT.(coalesce.(v, NaN))
end

"""
    read_dephy_driver(path; FT=Float64)

Read a DEPHY SCM driver file into a [`DEPHYDriver`](@ref) with ascending levels. Checks
that the file declares `format_version = "DEPHY SCM format version 1"`, that the `time`
axis is monotone, that no profile contains missing values, and that every
`nudging_constant_*` field is uniform (the SEA STARR drivers carry the inversion mask in
their global attributes rather than in these fields).
"""
function read_dephy_driver(path; FT=Float64)
    isfile(path) || throw(ArgumentError("DEPHY driver not found: $path"))
    ds = NCDataset(path)
    try
        attributes = Dict{String, Any}(k => v for (k, v) in ds.attrib)
        get(attributes, "format_version", "") == "DEPHY SCM format version 1" ||
            throw(ArgumentError("$path is not a DEPHY SCM format version 1 driver (format_version = $(get(attributes, "format_version", "missing")))"))
        case = string(get(attributes, "case", basename(path)))

        time_values = ds["time"][:]               # DateTime (decoded by NCDatasets)
        t0 = ds["t0"][1]
        all(t isa DateTime for t in time_values) || throw(ArgumentError("$path: time axis did not decode to DateTime"))
        start_time = DateTime(t0)
        times = FT[Dates.value(DateTime(t) - start_time) / 1000 for t in time_values]
        issorted(times) && allunique(times) || throw(ArgumentError("$path: forcing times are not strictly increasing"))
        end_time = haskey(attributes, "end_date") ? DateTime(replace(string(attributes["end_date"]), " " => "T")) : DateTime(last(time_values))

        zh = read_values(ds, "zh", FT)
        order = sortperm(zh)                        # the files are top-down; make ascending
        z = zh[order]
        issorted(z) && allunique(z) || throw(ArgumentError("$path: zh levels are not strictly monotone"))

        profile(name) = (v = read_values(ds, name, FT); vec(v)[order])
        matrix(name) = (m = read_values(ds, name, FT); m[order, :])
        check_finite(name, v) = all(isfinite, v) || throw(ArgumentError("$path: $name contains missing/non-finite values"))

        latitude = read_values(ds, "lat", FT)[1]
        longitude = read_values(ds, "lon", FT)[1]
        surface_pressure = read_values(ds, "ps", FT)[1]
        surface_temperature = read_values(ds, "ts", FT)[1]
        roughness_length = read_values(ds, "z0", FT)[1]

        na_units = get(ds["na"].attrib, "units", "")
        na_units == "mg-1" || throw(ArgumentError("$path: expected aerosol number units mg-1, found $(repr(na_units))"))
        qt_units = get(ds["qt"].attrib, "units", "")
        qt_units == "1" || throw(ArgumentError("$path: expected qt units \"1\" (mass fraction), found $(repr(qt_units))"))

        initial = (; thetal = profile("thetal"), qt = profile("qt"), na = aerosol_number_per_kg.(profile("na")),
                     u = profile("ua"), v = profile("va"), p = profile("pa"), T = profile("ta"))
        for (name, v) in pairs(initial)
            check_finite(name, v)
        end

        forcing = (; w = matrix("wa"), omega = matrix("wap"), ug = matrix("ug"), vg = matrix("vg"),
                     thetal_nud = matrix("thetal_nud"), T_nud = matrix("ta_nud"), qt_nud = matrix("qt_nud"),
                     na_nud = aerosol_number_per_kg.(matrix("na_nud")),
                     u_nud = matrix("ua_nud"), v_nud = matrix("va_nud"), o3 = matrix("o3"), p = matrix("pa_force"))
        for (name, m) in pairs(forcing)
            check_finite(name, m)
            size(m) == (length(z), length(times)) || throw(ArgumentError("$path: $name has size $(size(m)), expected $((length(z), length(times)))"))
        end

        function uniform_rate(name)
            m = read_values(ds, name, FT)
            check_finite(name, m)
            r = m[1]
            all(x -> x == r, m) || throw(ArgumentError("$path: $name is not uniform; the inversion-following mask must then come from the file"))
            return r
        end
        nudging_rates = (; thetal = uniform_rate("nudging_constant_thetal"), T = uniform_rate("nudging_constant_ta"),
                           qt = uniform_rate("nudging_constant_qt"), na = uniform_rate("nudging_constant_na"),
                           u = uniform_rate("nudging_constant_ua"), v = uniform_rate("nudging_constant_va"))

        sst = read_values(ds, "ts_force", FT); check_finite("ts_force", sst)
        ps = read_values(ds, "ps_force", FT); check_finite("ps_force", ps)
        length(sst) == length(times) || throw(ArgumentError("$path: ts_force length mismatch"))

        return DEPHYDriver{FT}(abspath(path), case, attributes, start_time, end_time, latitude, longitude, z, times,
                               surface_pressure, surface_temperature, roughness_length,
                               initial, forcing, nudging_rates, sst, ps)
    finally
        close(ds)
    end
end

"""
    driver_profile_time_series(grid, driver, name, z_centers)

The forcing matrix `driver.forcing[name]` interpolated linearly in height onto the model
cell centers `z_centers` (constant extrapolation beyond the driver's 5 m–76 km range) as a
`FieldTimeSeries{Nothing, Nothing, Center}` at `driver.times`; kernels interpolate linearly
in time through `fts[1, 1, k, Time(t)]`.
"""
function driver_profile_time_series(grid, driver::DEPHYDriver, name::Symbol, z_centers)
    m = driver.forcing[name]
    columns = [interpolate_profile(driver.z, view(m, :, n), z_centers) for n in eachindex(driver.times)]
    return profile_time_series(grid, driver.times, columns)
end

"""
    driver_initial_profile(driver, name, z)

The initial profile `driver.initial[name]` interpolated linearly in height at `z`.
"""
driver_initial_profile(driver::DEPHYDriver, name::Symbol, z) = interpolate_profile(driver.z, driver.initial[name], z)

"""
    inversion_height(z, θ; z_max=Inf)

Height of the maximum vertical gradient of the profile `θ(z)` (midpoint of the two levels),
searching levels below `z_max`. Used for the driver inventory and tests; the online
diagnosis on the model grid is in `sea_starr_forcings.jl`.
"""
function inversion_height(z, θ; z_max=Inf)
    keep = findall(≤(z_max), z)
    length(keep) ≥ 2 || throw(ArgumentError("need at least two levels below z_max"))
    zk = z[keep]; θk = θ[keep]
    g = diff(θk) ./ diff(zk)
    k = argmax(g)
    return (zk[k] + zk[k+1]) / 2
end

#####
##### Protocol vertical grid
#####

"""
    sea_starr_vertical_faces(; Δz=10, uniform_top=2500, stretch=1.1, top=6500)

Cell interfaces of the SEA STARR protocol grid: `Δz` spacing from the surface to
`uniform_top`, then each cell `stretch` times thicker than the one below "until `top` is
reached". The last interface is clipped to exactly `top` (with the defaults 38 stretched cells
would end at 6504.4 m; the clipped last cell is 370 m instead of 374 m). Returns Nz + 1 faces
with Nz = 288 for the defaults.
"""
function sea_starr_vertical_faces(; Δz=10.0, uniform_top=2500.0, stretch=1.1, top=6500.0)
    N_uniform = round(Int, uniform_top / Δz)
    N_uniform * Δz ≈ uniform_top || throw(ArgumentError("uniform_top must be a multiple of Δz"))
    faces = collect(range(0.0, uniform_top; length=N_uniform + 1))
    h = Δz
    while faces[end] < top - 1e-6
        h *= stretch
        push!(faces, faces[end] + h)
    end
    faces[end] = top
    return faces
end
