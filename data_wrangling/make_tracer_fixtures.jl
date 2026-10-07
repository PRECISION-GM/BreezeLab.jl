# Write the small redistributable excerpts used by test/tracer_dp_scream.jl:
#   test/fixtures/tracer_iop_excerpt.nc     — 4 hourly records (2022-08-05 00–03 UTC) of the ARM VARANAL
#                                              TRACER IOP file (public ARM data, doi 10.5439/1860369)
#   test/fixtures/dp_scream_3km_excerpt.nc  — the first 8 half-hourly records of the archived DP-SCREAM
#                                              3 km August 5–15 output (Zenodo 15271730, CC-BY-4.0)
#   test/fixtures/dp_scream_ns_excerpt.nc   — 3 records of the 3 km August run (nanosecond time units)
# Usage: julia --project data_wrangling/make_tracer_fixtures.jl [data/tracer_dp_scream] [test/fixtures]
using NCDatasets

data_dir = length(ARGS) ≥ 1 ? ARGS[1] : joinpath(@__DIR__, "..", "data", "tracer_dp_scream")
fixture_dir = length(ARGS) ≥ 2 ? ARGS[2] : joinpath(@__DIR__, "..", "test", "fixtures")
mkpath(fixture_dir)

function excerpt(src, dst, time_range; variables, attributes)
    NCDataset(src) do ds
        NCDataset(dst, "c") do out
            for (k, v) in ds.attrib
                k in attributes && (out.attrib[k] = v)
            end
            out.attrib["excerpt_of"] = basename(src)
            out.attrib["excerpt_records"] = string(first(time_range), ":", last(time_range))
            for (name, dim) in ds.dim
                defDim(out, name, name == "time" ? length(time_range) : dim)
            end
            for name in variables
                v = ds[name]
                dims = dimnames(v)
                raw = Array(v.var)  # undecoded values keep the file's own time encoding
                sel = "time" in dims ? selectdim(raw, findfirst(==("time"), dims), time_range) : raw
                attrib = Dict{String, Any}(k => a for (k, a) in v.attrib)
                newvar = defVar(out, name, eltype(raw), dims; attrib)
                newvar.var[ntuple(_ -> Colon(), ndims(sel))...] = collect(sel)
            end
        end
    end
    println("wrote ", dst)
end

iop = joinpath(data_dir, "TRACER_iopfile_4scam.nc")
# 2022-08-05 00 UTC is tsec = 35 days → record 841 (1-based)
excerpt(iop, joinpath(fixture_dir, "tracer_iop_excerpt.nc"), 841:844;
        variables = ["bdate", "tsec", "lev", "lat", "lon", "T", "q", "u", "v", "omega", "divT", "vertdivT", "divq", "vertdivq",
                     "Ps", "Tg", "Tsair", "shflx", "lhflx", "Prec", "srfswdn", "srfswup", "srflwdn", "srflwup", "windsrf"],
        attributes = ["datastream", "doi", "description", "references", "averaging_interval", "note", "Conventions"])

excerpt(joinpath(data_dir, "DP-SCREAM 3km_August_5_to_August_15_run.nc"), joinpath(fixture_dir, "dp_scream_3km_excerpt.nc"), 1:8;
        variables = ["time", "lev", "PRECL", "TMQ", "TGCLDLWP", "TGCLDIWP", "TOT_CLOUD_FRAC", "Z3", "T", "RELHUM", "W_SEC"],
        attributes = ["case", "git_version", "source_id", "time_period_freq", "Conventions", "source"])

excerpt(joinpath(data_dir, "DP-SCREAM 3km_August.nc"), joinpath(fixture_dir, "dp_scream_ns_excerpt.nc"), 1:3;
        variables = ["time", "lev", "PRECL", "TGCLDLWP", "T"],
        attributes = ["case", "git_version", "time_period_freq"])
