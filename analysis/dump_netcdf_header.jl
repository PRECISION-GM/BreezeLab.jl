# Print dimensions, variables (dims, type, units, long_name) and selected global attributes of NetCDF files.
#   julia --project analysis/dump_netcdf_header.jl FILE [FILE ...]
using NCDatasets
for path in ARGS
    println("=========== ", path)
    NCDataset(path) do ds
        println("dims: ", Dict(k => length(v) for (k, v) in ds.dim))
        for a in sort(collect(keys(ds.attrib)))
            v = string(ds.attrib[a]); println("GLOBAL ", a, " = ", v[1:min(end, 200)])
        end
        for (n, v) in ds
            at = Dict(a => v.attrib[a] for a in ("units", "long_name") if haskey(v.attrib, a))
            println(n, " ", dimnames(v), " ", eltype(v), " ", at)
        end
        if haskey(ds, "time")
            t = ds["time"][:]; println("time[1:3] = ", t[1:min(3, end)], " ... time[end] = ", t[end], " (", length(t), ")")
        end
    end
end
