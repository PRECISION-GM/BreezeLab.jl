#####
##### Reader for the DP-SCREAM TRACER reference outputs archived on Zenodo
##### (concept record 10.5281/zenodo.15133178; pinned version 15271730, CC-BY-4.0).
#####
##### Each NetCDF file holds domain-mean (horizontally averaged) EAMxx history at 30-minute
##### intervals on the 128 hybrid levels: PRECL, TMQ, TGCLDLWP, TGCLDIWP (time) and
##### TOT_CLOUD_FRAC, Z3, T, RELHUM, W_SEC, Q (time, lev). The archived `time` axes are
##### labelled in **local time** (CDT = UTC−5): the August run (case AUGUST_decr) starts at
##### "2022-07-31 19:00", i.e. 2022-08-01 00 UTC, the first IOP record of the month, and the
##### 3-month run at "2022-06-30 19:00" = 2022-07-01 00 UTC, the IOP file's base time. The
##### reader converts the labels back to UTC with `utc_offset_hours` (default 5). PRECL carries
##### no units attribute; its window mean (≈ 4.4 over 5–15 August) matches the IOP
##### precipitation in mm/day, which is the unit assumed here.
#####

using NCDatasets: NCDataset, nomissing
using Dates: Dates, DateTime

const DP_SCREAM_VARIABLES = ("PRECL", "TMQ", "TGCLDLWP", "TGCLDIWP", "TOT_CLOUD_FRAC", "Z3", "T", "RELHUM", "W_SEC", "Q")

const CF_TIME_FACTORS_MS = Dict("nanoseconds" => 1e-6, "microseconds" => 1e-3, "milliseconds" => 1.0,
                                "seconds" => 1e3, "minutes" => 6e4, "hours" => 3.6e6, "days" => 8.64e7)

"""
    parse_cf_time_units(units)

`(unit_ms, origin::DateTime)` from a CF `"<unit> since <date>"` string (the archive uses
"nanoseconds since 2022-08-04T19:00:00" and "minutes since 2022-08-04 19:00:00").
"""
function parse_cf_time_units(units::AbstractString)
    m = match(r"^\s*(\w+)\s+since\s+(.+?)\s*$", units)
    isnothing(m) && throw(ArgumentError("cannot parse time units \"$units\""))
    unit = lowercase(m.captures[1])
    haskey(CF_TIME_FACTORS_MS, unit) || throw(ArgumentError("unsupported time unit \"$unit\""))
    origin = DateTime(replace(m.captures[2], " " => "T"))
    return CF_TIME_FACTORS_MS[unit], origin
end

"""
    read_dp_scream_output(path; utc_offset_hours=5, variables=DP_SCREAM_VARIABLES)

Read one archived DP-SCREAM file. Returns a `NamedTuple` with `time` (UTC `DateTime`s),
`local_time` (the file's labels), `lev` (hybrid level midpoints, hPa, top first as in the
file), the available `variables` as `Vector`s (time series) or `Matrix`es (`lev × time`),
the global attributes, and the `case`/`git_version` strings recorded by EAMxx.
"""
function read_dp_scream_output(path; utc_offset_hours=5, variables=DP_SCREAM_VARIABLES)
    NCDataset(path) do ds
        raw = ds["time"].var[:]
        factor, origin = parse_cf_time_units(ds["time"].attrib["units"])
        local_time = [origin + Dates.Millisecond(round(Int, Float64(t) * factor)) for t in raw]
        time = local_time .+ Dates.Hour(utc_offset_hours)
        lev = Float64.(nomissing(ds["lev"][:], NaN))
        found = Symbol[]; values = Any[]
        for name in variables
            haskey(ds, name) || continue
            push!(found, Symbol(name))
            push!(values, Float64.(nomissing(ds[name][:], NaN)))
        end
        data = NamedTuple{Tuple(found)}(Tuple(values))
        attributes = Dict{String, String}(string(k) => string(v) for (k, v) in ds.attrib)
        return merge((; path = abspath(path), time, local_time, utc_offset_hours, lev,
                        case = get(attributes, "case", ""), git_version = get(attributes, "git_version", ""),
                        attributes), data)
    end
end

"""
    dp_scream_window(output, start, stop)

Indices of the records with `start ≤ time ≤ stop` (UTC).
"""
dp_scream_window(output, start::DateTime, stop::DateTime) = findall(t -> start ≤ t ≤ stop, output.time)
