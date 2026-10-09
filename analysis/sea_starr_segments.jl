# Read SEA STARR output written in several segments (one directory per restart segment) as one time series.
#
# A run interrupted and resumed from a checkpoint has its output in `runs/<run>/` (segment 1) and
# `runs/<run>_seg2_job<id>/` (segment 2, …). Segment k+1 restarts from a checkpoint at time t_k, while segment k
# may have continued past t_k before it stopped. `StitchedSeries` keeps segment k's records with t ≤ t_k and segment
# k+1's records with t > t_k, where t_k is `resume_time_s` from segment k+1's `provenance_resume.toml` (written by
# execution/sea_starr_resume.jl) or, failing that, its first record time. Indexing returns the underlying `Field`,
# so code written for a `FieldTimeSeries` (`s[n]`, `s.times`, `s.grid`, `interior(s[n], …)`) works unchanged.

using Oceananigans
using Oceananigans.OutputReaders: OnDisk, InMemory
using TOML

struct StitchedSeries{G, S}
    times :: Vector{Float64}
    grid :: G
    sources :: S                  # Vector of (FieldTimeSeries, index) per record
    segments :: Vector{String}
end

Base.getindex(s::StitchedSeries, n::Int) = (source = s.sources[n]; source[1][source[2]])
Base.length(s::StitchedSeries) = length(s.times)
Base.lastindex(s::StitchedSeries) = length(s.times)
Base.eachindex(s::StitchedSeries) = eachindex(s.times)

"""
    segment_dirs(spec)

Split a comma-separated list of run directories (segment 1 first).
"""
segment_dirs(spec::AbstractString) = String.(filter(!isempty, strip.(split(spec, ','))))

"""
    restart_time(dir)

Model time [s] from which segment `dir` was resumed (`provenance_resume.toml`), or `nothing` for a first segment.
"""
function restart_time(dir)
    path = joinpath(dir, "provenance_resume.toml")
    isfile(path) || return nothing
    extra = get(TOML.parsefile(path), "extra", Dict())
    t = get(extra, "resume_time_s", nothing)
    return t isa Number ? Float64(t) : nothing
end

"""
    stitched_series(dirs, filename, name; backend = InMemory())

`StitchedSeries` of output `name` from `<dir>/<filename>` over the segment directories `dirs` (missing files in
later segments are skipped, e.g. while a segment has not written yet).
"""
function stitched_series(dirs, filename, name; backend = InMemory())
    series = [FieldTimeSeries(joinpath(d, filename), name; backend) for d in dirs if isfile(joinpath(d, filename))]
    present = [d for d in dirs if isfile(joinpath(d, filename))]
    isempty(series) && error("no $filename in $(dirs)")
    times = Float64[]; sources = Tuple{Any, Int}[]
    for (k, fts) in enumerate(series)
        t_end = k < length(series) ? something(restart_time(present[k + 1]), Float64(first(series[k + 1].times))) : Inf
        t_start = isempty(times) ? -Inf : last(times)
        for (i, t) in enumerate(fts.times)
            (t > t_start && t ≤ t_end) || continue
            push!(times, t); push!(sources, (fts, i))
        end
    end
    return StitchedSeries(times, series[1].grid, sources, present)
end

segment_label(dirs) = join(basename.(rstrip.(dirs, '/')), " + ")
