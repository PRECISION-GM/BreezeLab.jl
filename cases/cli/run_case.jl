# Run one ENA case (Covert or LASSO) from the command line.
#
#   julia --project cases/cli/run_case.jl --protocol covert_public_bin --data data/covert2022_bin \
#         --microphysics p3_n75 --arch gpu --Nx 256 --Ny 256 --hours 6 --output output/p3_n75
#   julia --project cases/cli/run_case.jl --protocol lasso_ena_official --data data/lasso/<member> \
#         --member <member> --epoch 2017-07-18T00:00:00 --arch gpu --output output/lasso
#
# This script is the only place where microphysics names are mapped to Breeze objects (a
# shell cannot pass a Julia object); the constructors take the objects. Every run writes
# <output>/provenance.toml (inputs + checksums + configuration).

using BreezeLab, Breeze, Oceananigans, Oceananigans.Units, CUDA, Dates
using Breeze.Microphysics.PredictedParticleProperties: CloudDroplets

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
allowed = Set(["data", "protocol", "arch", "float", "seed", "epoch", "member", "dimensions", "microphysics", "aerosol_ss_cap",
               "Nx", "Ny", "Lx", "Ly", "hours", "cfl", "dt", "max_dt", "lasso_grid", "grid", "initialization", "closure",
               "profile_interval", "slice_interval", "radiation", "nudging", "surface_flux_law", "output"])
unknown = setdiff(Set(keys(opts)), allowed)
isempty(unknown) || error("unknown options: $(join(sort!(collect(unknown)), ", "))")

protocol = Symbol(getopt("protocol", ""))
protocol in (:covert_public_bin, :lasso_ena_official) || error("select --protocol covert_public_bin or --protocol lasso_ena_official")
haskey(opts, "data") || error("--data must point to the selected protocol's input directory")
data = opts["data"]
getopt("arch", "cpu") in ("cpu", "gpu") || error("--arch must be cpu or gpu")
getopt("float", "Float32") in ("Float32", "Float64") || error("--float must be Float32 or Float64")
Oceananigans.defaults.FloatType = getopt("float", "Float32") == "Float64" ? Float64 : Float32
arch = getopt("arch", "cpu") == "gpu" ? GPU() : CPU()
# --grid covert|covert_inversion_5m|lasso: a named vertical grid (recorded in provenance as `vertical_grid`);
# --lasso_grid true is the older spelling of --grid lasso. Without either the protocol's own grid is used.
vertical_grid = get(opts, "grid", getopt("lasso_grid", "false") == "true" ? "lasso" : "protocol")
z_faces = vertical_grid == "protocol" ? nothing : ena_vertical_faces(vertical_grid)

# --microphysics: the members of the ENA campaigns as Breeze objects. The aerosol members convert
# LASSO's per-cm³ numbers with the case's first-level reference density.
function microphysics_member(name; data, z_faces, cap)
    cloud = CloudDroplets(; number_concentration = 75e6)
    name == "one_moment" && return Base.get_extension(Breeze, :BreezeCloudMicrophysicsExt).OneMomentCloudMicrophysics(;
                                       cloud_formation = SaturationAdjustment(; equilibrium = WarmPhaseEquilibrium()))
    name == "p3_n75" && return P3Microphysics(; cloud)
    reference_density = first_level_reference_density(data, z_faces)
    if name in ("p3_aer1", "p3_aer2", "p3_aer3")
        aerosol = lasso_aerosol(; setting = Symbol(name[4:end]), reference_density, maximum_supersaturation = cap)
        return P3Microphysics(; cloud, aerosol)
    elseif name == "p3_covert_n75"
        return P3Microphysics(; cloud, aerosol = covert_aerosol(; reference_density, maximum_supersaturation = something(cap, 0.003)))
    end
    error("--microphysics must be one_moment, p3_n75, p3_aer1, p3_aer2, p3_aer3 or p3_covert_n75")
end

kw = Dict{Symbol, Any}()
haskey(opts, "Nx") && (kw[:Nx] = parse(Int, opts["Nx"]))
haskey(opts, "Ny") && (kw[:Ny] = parse(Int, opts["Ny"]))
haskey(opts, "Lx") && (kw[:Lx] = parse(Float64, opts["Lx"]))
haskey(opts, "Ly") && (kw[:Ly] = parse(Float64, opts["Ly"]))
haskey(opts, "hours") && (kw[:stop_time] = parse(Float64, opts["hours"]) * 3600)
haskey(opts, "cfl") && (kw[:cfl] = parse(Float64, opts["cfl"]))
haskey(opts, "dt") && (kw[:Δt] = parse(Float64, opts["dt"]))
haskey(opts, "max_dt") && (kw[:max_Δt] = parse(Float64, opts["max_dt"]))
haskey(opts, "initialization") && (kw[:initialization] = Symbol(opts["initialization"]))
haskey(opts, "profile_interval") && (kw[:profile_interval] = parse(Float64, opts["profile_interval"]))
haskey(opts, "slice_interval") && (kw[:slice_interval] = parse(Float64, opts["slice_interval"]))
if haskey(opts, "closure")
    opts["closure"] in ("none", "smagorinsky_lilly") || throw(ArgumentError("--closure must be none or smagorinsky_lilly"))
    kw[:closure] = opts["closure"] == "none" ? nothing : SmagorinskyLilly()
end
isnothing(z_faces) || (kw[:z_faces] = z_faces)
# Output files keep the campaign prefix that the analysis scripts read (`lasso_ena_profiles.jld2`, ...).
output_prefix = "lasso_ena"
cap = getopt("aerosol_ss_cap", "nothing") == "nothing" ? nothing : parse(Float64, opts["aerosol_ss_cap"])

if protocol === :covert_public_bin
    for option in ("member", "epoch", "dimensions", "radiation", "nudging", "surface_flux_law")
        haskey(opts, option) && error("--$option applies to --protocol lasso_ena_official only")
    end
    member_name = getopt("microphysics", "p3_n75")
    faces = something(z_faces, covert_public_bin_vertical_faces())
    microphysics = microphysics_member(member_name; data, z_faces = faces, cap)
    output_dir = getopt("output", "output/covert_public_bin_$(getopt("microphysics", "p3_n75"))")
    @info "Building the Covert case on $(typeof(arch)) with $(Oceananigans.defaults.FloatType)" opts
    seed = parse(Int, getopt("seed", "1234"))
    case = ena_covert(; arch, data_dir = data, microphysics, output_dir, output_prefix,
                        perturbation = InitialPerturbation(; seed), kw...)
else
    haskey(opts, "member") || error("--member <LASSO run ID> is required with --protocol lasso_ena_official")
    haskey(opts, "epoch") || error("--epoch <UTC DateTime> is required with --protocol lasso_ena_official")
    haskey(opts, "seed") && error("--seed applies to the Covert protocol; LASSO uses the namelist perturbation")
    member = opts["member"]
    epoch = DateTime(opts["epoch"])
    dimensions = haskey(opts, "dimensions") ? Tuple(parse.(Int, split(opts["dimensions"], ','))) :
                 BreezeLab.lasso_documented_dimensions(member)
    length(dimensions) == 3 || error("--dimensions must be Nx,Ny,Nz from the bundle domain configuration")
    if haskey(opts, "radiation")
        opts["radiation"] in ("rrtmgp", "none", "nothing") || error("--radiation must be rrtmgp or none")
        kw[:radiation] = opts["radiation"] == "rrtmgp"
    end
    if haskey(opts, "nudging")      # a timescale in seconds, or none/nothing to switch the wind nudging off
        opts["nudging"] in ("none", "nothing") ? (kw[:wind_nudging] = false) :
                                                 (kw[:wind_nudging_timescale] = parse(Float64, opts["nudging"]))
    end
    haskey(opts, "surface_flux_law") && (kw[:surface_flux_law] = Symbol(opts["surface_flux_law"]))
    # Without --microphysics the member's protocol P3 is built by ena_lasso (capped at 0.3 %). A
    # LASSO aerosol member from the CLI is capped at the SBM ss_max as well unless --aerosol_ss_cap says otherwise.
    member_name = get(opts, "microphysics", haskey(opts, "aerosol_ss_cap") ? "p3_" * string(parse_lasso_member(member).aerosol) : "protocol")
    cap = haskey(opts, "aerosol_ss_cap") ? cap : 0.003
    if member_name != "protocol"
        faces = isnothing(z_faces) ? validate_lasso_bundle(data; member, dimensions, epoch).grid.faces : z_faces
        kw[:microphysics] = microphysics_member(member_name; data, z_faces = faces, cap)
    end
    output_dir = getopt("output", "output/lasso_ena_official_$(getopt("microphysics", "protocol"))")
    @info "Building LASSO member $member on $(typeof(arch)) with $(Oceananigans.defaults.FloatType)" opts
    case = ena_lasso(; member, epoch, dimensions, bundle_dir = data, arch, output_dir, output_prefix, kw...)
end

mkpath(output_dir)
provenance = write_provenance(joinpath(output_dir, "provenance.toml"), case;
                              extra = (; command = join(ARGS, " "), hostname = gethostname(), vertical_grid,
                                         microphysics_member = member_name, aerosol_supersaturation_cap = something(cap, 0.0),
                                         gpu = CUDA.functional() ? CUDA.name(CUDA.device()) : "none"))
@info "Provenance written to $provenance"
@info "Label: $(case.config.label)"
println(case.model)

wall = time()
run!(case.simulation)
@info "Finished in $(round((time() - wall) / 60, digits=1)) minutes; output in $output_dir"
