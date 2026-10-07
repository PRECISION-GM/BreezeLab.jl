# Run one LASSO-ENA / Covert case from the command line.
#
#   julia --project cases/cli/run_case.jl --data data/covert2022_bin --protocol covert_public_bin \
#         --microphysics p3_n75 --arch gpu --Nx 256 --Ny 256 --hours 6 --output output/p3_n75
#
# Every run writes <output>/provenance.toml (inputs + checksums + configuration).

using BreezeLab, Breeze, Oceananigans, Oceananigans.Units, CUDA, Dates, Random

function parse_args(args)
    opts = Dict{String, String}()
    i = 1
    while i ≤ length(args)
        key = args[i]
        startswith(key, "--") || error("unexpected argument $key")
        i < length(args) || error("missing value for $key")
        haskey(opts, key[3:end]) && error("duplicate option $key")
        opts[key[3:end]] = args[i+1]
        i += 2
    end
    return opts
end

opts = parse_args(ARGS)
getopt(k, default) = get(opts, k, default)
allowed = Set(["data", "protocol", "preset", "arch", "float", "seed", "epoch", "member", "dimensions", "microphysics", "Nx", "Ny", "Lx", "Ly", "hours", "radiation", "surface", "nudging", "vertical_advection", "p3_initialization", "aerosol_replenishment", "cfl", "dt", "max_dt", "lasso_grid", "moment_advection", "formulation", "profile_interval", "closure", "slice_interval", "output"])
unknown = setdiff(Set(keys(opts)), allowed)
isempty(unknown) || error("unknown options: $(join(sort!(collect(unknown)), ", "))")

haskey(opts, "protocol") && haskey(opts, "preset") && error("use --protocol, or the legacy --preset, not both")
preset = Symbol(get(opts, "protocol", get(opts, "preset", "")))
preset in (:covert_public_bin, :lasso_ena_official) || error("select --protocol covert_public_bin or --protocol lasso_ena_official")
haskey(opts, "data") || error("--data must point to the selected protocol's input directory")
data = opts["data"]
getopt("arch", "cpu") in ("cpu", "gpu") || error("--arch must be cpu or gpu")
getopt("float", "Float32") in ("Float32", "Float64") || error("--float must be Float32 or Float64")
arch = lowercase(getopt("arch", "cpu")) == "gpu" ? GPU() : CPU()
FT = getopt("float", "Float32") == "Float64" ? Float64 : Float32
seed = parse(Int, getopt("seed", "1234"))

kw = Dict{Symbol, Any}()
if haskey(opts, "dimensions")
    dims = Tuple(parse.(Int, split(opts["dimensions"], ',')))
    length(dims) == 3 || error("--dimensions must be Nx,Ny,Nz from the bundle domain configuration")
    kw[:dimensions] = dims
end
haskey(opts, "epoch") && (kw[:epoch] = DateTime(opts["epoch"]))
haskey(opts, "member") && (kw[:member] = opts["member"])      # LASSO run ID, e.g. 20170718era5d25x100_sbmwrm-aer2-flxsst
preset === :lasso_ena_official && !haskey(opts, "member") && error("--member <LASSO run ID> is required with --protocol lasso_ena_official")
haskey(opts, "microphysics") && (kw[:microphysics] = Symbol(opts["microphysics"]))   # otherwise the preset decides
haskey(opts, "Nx") && (kw[:Nx] = parse(Int, opts["Nx"]))
haskey(opts, "Ny") && (kw[:Ny] = parse(Int, opts["Ny"]))
haskey(opts, "Lx") && (kw[:Lx] = parse(Float64, opts["Lx"]))
haskey(opts, "Ly") && (kw[:Ly] = parse(Float64, opts["Ly"]))
haskey(opts, "hours") && (kw[:stop_time] = parse(Float64, opts["hours"]) * 3600)
haskey(opts, "radiation") && (kw[:radiation] = opts["radiation"] == "nothing" ? nothing : Symbol(opts["radiation"]))
haskey(opts, "surface") && (kw[:surface] = Symbol(opts["surface"]))
haskey(opts, "nudging") && (kw[:wind_nudging_timescale] = opts["nudging"] == "nothing" ? nothing : parse(Float64, opts["nudging"]))
haskey(opts, "vertical_advection") && (kw[:vertical_advection] = opts["vertical_advection"] == "nothing" ? nothing : Symbol(opts["vertical_advection"]))
haskey(opts, "p3_initialization") && (kw[:p3_initialization] = Symbol(opts["p3_initialization"]))
haskey(opts, "aerosol_replenishment") && (kw[:aerosol_replenishment] = opts["aerosol_replenishment"] == "nothing" ? nothing :
                                          opts["aerosol_replenishment"] == "diagnostic_ccn" ? :diagnostic_ccn : parse(Float64, opts["aerosol_replenishment"]))
haskey(opts, "cfl") && (kw[:cfl] = parse(Float64, opts["cfl"]))
haskey(opts, "max_dt") && (kw[:max_Δt] = parse(Float64, opts["max_dt"]))
haskey(opts, "dt") && (kw[:Δt] = parse(Float64, opts["dt"]))
haskey(opts, "lasso_grid") && opts["lasso_grid"] == "true" && (kw[:z_faces] = lasso_ena_vertical_faces())
haskey(opts, "moment_advection") && (kw[:moment_advection] = Symbol(opts["moment_advection"]))
haskey(opts, "formulation") && (kw[:formulation] = Symbol(opts["formulation"]))   # LiquidIcePotentialTemperature (default) or StaticEnergy
haskey(opts, "profile_interval") && (kw[:profile_interval] = parse(Float64, opts["profile_interval"]))
if haskey(opts, "closure")
    opts["closure"] in ("none", "smagorinsky_lilly") ||
        throw(ArgumentError("closure must be none or smagorinsky_lilly"))
    kw[:closure] = opts["closure"] == "none" ? nothing : :smagorinsky_lilly
end
haskey(opts, "slice_interval") && (kw[:slice_interval] = parse(Float64, opts["slice_interval"]))

microphysics_label = get(opts, "microphysics", preset === :lasso_ena_official ? "p3_aer2" : "p3_n75")
output_dir = getopt("output", "output/$(preset)_$(microphysics_label)")
@info "Building $preset on $(typeof(arch)) with $FT" opts
case = ena_simulation(data; protocol=preset, arch, FT, output_dir,
                            perturbation = InitialPerturbation(seed=seed), kw...)
mkpath(output_dir)
provenance = write_provenance(joinpath(output_dir, "provenance.toml"), case;
                              extra = (; command = join(ARGS, " "), hostname = gethostname(),
                                         gpu = CUDA.functional() ? CUDA.name(CUDA.device()) : "none"))
@info "Provenance written to $provenance"
@info "Label: $(case.config.label)"
println(case.model)

wall = time()
run!(case.simulation)
@info "Finished in $(round((time() - wall) / 60, digits=1)) minutes; output in $output_dir"
