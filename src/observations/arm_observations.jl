#####
##### ARM ENA station observations and BreezeLab run output as comparable time series.
#####
##### Readers keep the native UTC samples and the product QC; nothing is regridded to the
##### LES. Statistics are taken over an explicit UTC window so that a run and the
##### observations are compared on the same interval. The NumericalEarth metadata adapter
##### design in data_wrangling/README.md is the intended long-term home of the discovery/
##### download side; these readers are the station-series part of that design.
#####

using NCDatasets: NCDataset
using Dates: DateTime, Second, Millisecond
using Statistics: mean, median, std, quantile
using TOML: TOML
using Oceananigans: FieldTimeSeries
using Oceananigans.Grids: Center, znodes
using Oceananigans.Fields: interior

"""
    ARMSeries

One station time series: `datastream`, `variable`, `units`, UTC `time`, `value` (`NaN`
where the file has a fill value), optional 1-σ `uncertainty`, and `good` (the product's
QC passed and the value is finite).
"""
struct ARMSeries{V, U}
    datastream :: String
    variable :: String
    units :: String
    time :: Vector{DateTime}
    value :: V
    uncertainty :: U
    good :: BitVector
end

Base.length(s::ARMSeries) = length(s.time)
Base.summary(s::ARMSeries) = string("ARMSeries(", s.datastream, " ", s.variable, " [", s.units, "], ", length(s), " samples, ",
                                    count(s.good), " good)")
Base.show(io::IO, s::ARMSeries) = print(io, summary(s))

to_float(x) = Float64[ismissing(v) ? NaN : Float64(v) for v in x]
to_datetime(x) = DateTime[DateTime(t) for t in x]
attribute(ds, name, default="") = haskey(ds.attrib, name) ? string(ds.attrib[name]) : default
variable_attribute(v, name, default="") = haskey(v.attrib, name) ? string(v.attrib[name]) : default
# ARM `qc_*` fields are bit-packed test results; zero means every test passed.
qc_passed(ds, name, n) = haskey(ds, name) ? BitVector(Bool[!ismissing(v) && v == 0 for v in ds[name][:]]) : trues(n)

"""
    read_arm_lwp(path; variable="phys_lwp")

MWRRET liquid water path [g m⁻²] (`enamwrret2turnC1.c1`): `phys_lwp` with its
`qc_phys_lwp` and `phys_qc_flag` (0 = retrieval good) and `phys_lwp_uncertainty`. Use
`variable="stat_lwp"` for the statistical retrieval.
"""
function read_arm_lwp(path; variable="phys_lwp")
    NCDataset(path) do ds
        time = to_datetime(ds["time"][:])
        value = to_float(ds[variable][:])
        n = length(time)
        good = qc_passed(ds, "qc_" * variable, n) .& isfinite.(value)
        variable == "phys_lwp" && (good .&= qc_passed(ds, "phys_qc_flag", n))
        uncertainty = haskey(ds, variable * "_uncertainty") ? to_float(ds[variable * "_uncertainty"][:]) : nothing
        return ARMSeries(attribute(ds, "datastream", basename(path)), variable, variable_attribute(ds[variable], "units"),
                         time, value, uncertainty, good)
    end
end

"""
    read_arm_rain_rate(path)

Video disdrometer rain rate [mm hr⁻¹] (`enavdisC1.b1`, one-minute samples) with
`qc_rain_rate`; negative or fill values are not good.
"""
function read_arm_rain_rate(path)
    NCDataset(path) do ds
        time = to_datetime(ds["time"][:])
        value = to_float(ds["rain_rate"][:])
        good = qc_passed(ds, "qc_rain_rate", length(time)) .& isfinite.(value) .& (value .>= 0)
        return ARMSeries(attribute(ds, "datastream", basename(path)), "rain_rate", variable_attribute(ds["rain_rate"], "units"),
                         time, value, nothing, good)
    end
end

"""
    read_arm_cloud_boundaries(path)

KAZR ARSCL cloud boundaries (`enaarsclkazrbnd1kolliasC1.c0`, 4-s samples): the base and
top [m] of the lowest hydrometeor layer (`cloud_layer_base_height[1, :]`,
`cloud_layer_top_height[1, :]`; `NaN` for clear sky or fill), the ceilometer/MPL
`cloud_base_best_estimate`, the number of layers, and `cloudy` (a valid lowest layer).
"""
function read_arm_cloud_boundaries(path)
    NCDataset(path) do ds
        time = to_datetime(ds["time"][:])
        bases = ds["cloud_layer_base_height"][:, :]
        tops = ds["cloud_layer_top_height"][:, :]
        n = length(time)
        valid(x) = !ismissing(x) && x >= 0
        base = Float64[valid(bases[1, i]) ? bases[1, i] : NaN for i in 1:n]
        top = Float64[valid(tops[1, i]) ? tops[1, i] : NaN for i in 1:n]
        n_layers = Int[count(l -> valid(bases[l, i]), axes(bases, 1)) for i in 1:n]
        best = haskey(ds, "cloud_base_best_estimate") ? to_float(ds["cloud_base_best_estimate"][:]) : fill(NaN, n)
        best = Float64[isfinite(b) && b >= 0 ? b : NaN for b in best]
        return (; datastream = attribute(ds, "datastream", basename(path)), time, base, top, n_layers,
                  base_best_estimate = best, cloudy = BitVector(isfinite.(base)))
    end
end

"""
    window_statistics(time, value, good, t0, t1)
    window_statistics(series::ARMSeries, t0, t1)

Statistics of the good samples with `t0 ≤ time < t1`: `n` (samples in the window),
`n_good`, `fraction_good`, `mean`, `median`, `std`, `q10`, `q90`, `min`, `max`
(`NaN` when no good sample falls in the window).
"""
function window_statistics(time, value, good, t0::DateTime, t1::DateTime)
    inside = BitVector(Bool[t0 ≤ t < t1 for t in time])
    selected = inside .& good
    n = count(inside)
    x = value[selected]
    isempty(x) && return (; n, n_good = 0, fraction_good = 0.0, mean = NaN, median = NaN, std = NaN,
                             q10 = NaN, q90 = NaN, min = NaN, max = NaN)
    return (; n, n_good = length(x), fraction_good = length(x) / n,
              mean = mean(x), median = median(x), std = length(x) > 1 ? std(x) : 0.0,
              q10 = quantile(x, 0.1), q90 = quantile(x, 0.9), min = minimum(x), max = maximum(x))
end

window_statistics(s::ARMSeries, t0::DateTime, t1::DateTime) = window_statistics(s.time, s.value, s.good, t0, t1)

"""
    cloud_fraction_from_boundaries(bounds, t0, t1; max_base=3000)

From [`read_arm_cloud_boundaries`](@ref): the fraction of samples in `[t0, t1)` whose
lowest layer base lies at or below `max_base` [m] (a time-fraction cloud cover at the
zenith, the observational counterpart of the LES area cloud fraction), with the base and
top statistics of those samples.
"""
function cloud_fraction_from_boundaries(bounds, t0::DateTime, t1::DateTime; max_base=3000)
    inside = BitVector(Bool[t0 ≤ t < t1 for t in bounds.time])
    low = BitVector(Bool[isfinite(b) && b ≤ max_base for b in bounds.base])
    n = count(inside)
    n_low = count(inside .& low)
    return (; n, n_cloudy = n_low, cloud_fraction = n == 0 ? NaN : n_low / n,
              base = window_statistics(bounds.time, bounds.base, low, t0, t1),
              top = window_statistics(bounds.time, bounds.top, low .& BitVector(isfinite.(bounds.top)), t0, t1))
end

#####
##### BreezeLab run output
#####

run_file(run_dir, suffix) = only(filter(f -> endswith(f, suffix), readdir(run_dir; join=true)))

function run_provenance(run_dir)
    path = joinpath(run_dir, "provenance.toml")
    isfile(path) || throw(ArgumentError("$run_dir has no provenance.toml; the epoch and protocol label come from it"))
    return TOML.parsefile(path)
end

"""
    breezelab_timeseries(run_dir)

The domain-mean time series a BreezeLab ENA run wrote (`*_timeseries.jld2`) as UTC series:
`lwp`, `rwp` [g m⁻²], `cloud_fraction`, `rain_rate` [mm hr⁻¹] (from the surface rain mass
flux, 1 kg m⁻² = 1 mm), with `epoch`, `seconds`, `time`, and the provenance `protocol`
and `label`. The epoch comes from the run's `provenance.toml`, never from a file name.
"""
function breezelab_timeseries(run_dir)
    provenance = run_provenance(run_dir)
    epoch = DateTime(provenance["config"]["epoch"])
    file = run_file(run_dir, "_timeseries.jld2")
    lwp = FieldTimeSeries(file, "lwp")
    rwp = FieldTimeSeries(file, "rwp")
    cf = FieldTimeSeries(file, "cloud_fraction")
    rain = FieldTimeSeries(file, "rain_flux")
    seconds = collect(Float64, lwp.times)
    scalar(fts, n) = Float64(Array(interior(fts[n]))[1])
    time = DateTime[epoch + Millisecond(round(Int, 1000s)) for s in seconds]
    return (; protocol = get(provenance, "protocol", "unknown"), label = get(provenance["config"], "label", ""),
              epoch, seconds, time,
              lwp = [1e3 * scalar(lwp, n) for n in eachindex(seconds)],
              rwp = [1e3 * scalar(rwp, n) for n in eachindex(seconds)],
              cloud_fraction = [scalar(cf, n) for n in eachindex(seconds)],
              rain_rate = [3600 * scalar(rain, n) for n in eachindex(seconds)])
end

"""
    breezelab_cloud_boundaries(run_dir; threshold=0.05)

Cloud base and top [m] of a run from its saved horizontal-mean cloud-fraction profiles
(`*_profiles.jld2`, averaged over each output interval): the lowest and highest cell
centers whose cloudy-cell fraction exceeds `threshold` (`NaN` when none does), at the end
of each averaging interval.
"""
function breezelab_cloud_boundaries(run_dir; threshold=0.05)
    provenance = run_provenance(run_dir)
    epoch = DateTime(provenance["config"]["epoch"])
    file = run_file(run_dir, "_profiles.jld2")
    cf = FieldTimeSeries(file, "cloud_fraction")
    z = collect(znodes(cf.grid, Center()))
    seconds = collect(Float64, cf.times)
    base = Float64[]; top = Float64[]; maximum_fraction = Float64[]
    for n in eachindex(seconds)
        profile = vec(Array(interior(cf[n])))
        cloudy = findall(>(threshold), profile)
        push!(maximum_fraction, maximum(profile))
        push!(base, isempty(cloudy) ? NaN : z[first(cloudy)])
        push!(top, isempty(cloudy) ? NaN : z[last(cloudy)])
    end
    time = DateTime[epoch + Millisecond(round(Int, 1000s)) for s in seconds]
    return (; epoch, seconds, time, base, top, maximum_fraction, threshold)
end
