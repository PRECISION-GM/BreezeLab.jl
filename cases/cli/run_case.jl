# Run one ENA case (Covert or LASSO) from the command line.
#
#   julia --project cases/cli/run_case.jl --protocol covert_public_bin --data data/covert2022_bin \
#         --microphysics p3_n75 --arch gpu --Nx 256 --Ny 256 --hours 6 --output output/p3_n75
#   julia --project cases/cli/run_case.jl --protocol covert_public_bin --data data/covert2022_bin \
#         --resolution r512 --arch gpu --hours 6 --checkpoint_interval 1800 --output output/r512
#   (add --pickup true to resume from the latest checkpoint in --output, appending to its output files)
#   julia --project cases/cli/run_case.jl --protocol lasso_ena_official --data data/lasso/<member> \
#         --member <member> --epoch 2017-07-18T00:00:00 --arch gpu --output output/lasso
#
# This script is the only place where microphysics names are mapped to Breeze objects (a
# shell cannot pass a Julia object); the constructors take the objects. Every run writes
# <output>/provenance.toml (inputs + checksums + configuration).

using BreezeLab, Breeze, Oceananigans, Oceananigans.Units, CUDA, Dates
using Breeze.Microphysics.PredictedParticleProperties: CloudDroplets, ProcessRate

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
               "resolution", "Nx", "Ny", "Nz", "Lx", "Ly", "hours", "checkpoint_interval", "pickup",
               "process_rates", "process_rate_interval", "prandtl", "scalar_advection_order", "reference_extension", "lw_F0", "dynamic_averaging", "advection_order", "microphysics_advection_order", "cfl", "dt", "max_dt", "lasso_grid", "grid", "initialization", "closure",
               "profile_interval", "slice_interval", "liquid_reff", "radiation", "nudging", "surface_flux_law", "output",
               "nudging_tq", "nudging_z", "nudging_uv", "p3_scale"])
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

# --p3_scale name=factor,...: multiply P3 warm-rain process-rate parameters by the given factors
# (autoconversion_coefficient, accretion_coefficient, rain_self_collection_coefficient,
# initial_rain_drop_mass; e.g. initial_rain_drop_mass=4.096 is a 40 µm rather than 25 µm seed drop).
# Every other process-rate parameter keeps its default.
const P3_SCALABLE = ("autoconversion_coefficient", "accretion_coefficient",
                     "rain_self_collection_coefficient", "initial_rain_drop_mass")

function p3_process_rates(spec)
    isnothing(spec) && return nothing
    defaults = ProcessRate()
    scaled = Dict{Symbol, Any}()
    for entry in split(spec, ',')
        name, factor = split(entry, '=')
        name in P3_SCALABLE || error("--p3_scale: $name is not one of $(join(P3_SCALABLE, ", "))")
        scaled[Symbol(name)] = getproperty(defaults, Symbol(name)) * parse(Float64, factor)
    end
    return ProcessRate(; scaled...)
end

# --microphysics: the members of the ENA campaigns as Breeze objects. The aerosol members convert
# LASSO's per-cm³ numbers with the case's first-level reference density.
function microphysics_member(name; data, z_faces, cap, process_rates = nothing)
    cloud = CloudDroplets(; number_concentration = 75e6)
    !isnothing(process_rates) && !startswith(name, "p3_") && error("--p3_scale needs P3 microphysics")
    name == "one_moment" && return Base.get_extension(Breeze, :BreezeCloudMicrophysicsExt).OneMomentCloudMicrophysics(;
                                       cloud_formation = SaturationAdjustment(; equilibrium = WarmPhaseEquilibrium()))
    name == "p3_n75" && return P3Microphysics(; cloud, process_rates)
    # p3_n<N>: P3 with a prescribed droplet number of N cm⁻³ (p3_n50, p3_n100, ...), as p3_n75
    m = match(r"^p3_n(\d+)$", name)
    isnothing(m) || return P3Microphysics(; cloud = CloudDroplets(; number_concentration = parse(Int, m.captures[1]) * 1e6), process_rates)
    reference_density = first_level_reference_density(data, z_faces)
    if name in ("p3_aer1", "p3_aer2", "p3_aer3")
        aerosol = lasso_aerosol(; setting = Symbol(name[4:end]), reference_density, maximum_supersaturation = cap)
        return P3Microphysics(; cloud, aerosol, process_rates)
    elseif name == "p3_covert_n75"
        return P3Microphysics(; cloud, aerosol = covert_aerosol(; reference_density, maximum_supersaturation = something(cap, 0.003)), process_rates)
    end
    error("--microphysics must be one_moment, p3_n75 (or p3_n<N> for N cm⁻³), p3_aer1, p3_aer2, p3_aer3 or p3_covert_n75")
end

kw = Dict{Symbol, Any}()
# --resolution r128|r256|r512 (Covert only): Nx, Ny, Nz, Lx, Ly and Δt of COVERT_RESOLUTIONS; explicit options override
if haskey(opts, "resolution")
    protocol === :covert_public_bin || error("--resolution applies to --protocol covert_public_bin only")
    haskey(COVERT_RESOLUTIONS, Symbol(opts["resolution"])) ||
        error("--resolution must be one of $(join(keys(COVERT_RESOLUTIONS), ", "))")
    merge!(kw, Dict(pairs(COVERT_RESOLUTIONS[Symbol(opts["resolution"])])))
end
haskey(opts, "Nx") && (kw[:Nx] = parse(Int, opts["Nx"]))
haskey(opts, "Ny") && (kw[:Ny] = parse(Int, opts["Ny"]))
haskey(opts, "Nz") && (kw[:Nz] = parse(Int, opts["Nz"]))
haskey(opts, "checkpoint_interval") && (kw[:checkpoint_interval] = parse(Float64, opts["checkpoint_interval"]))
# --pickup true: resume from the latest checkpoint (needs --checkpoint_interval) and append to the output files
getopt("pickup", "false") in ("true", "false") || error("--pickup must be true or false")
pickup = getopt("pickup", "false") == "true"
pickup && !haskey(kw, :checkpoint_interval) && error("--pickup needs --checkpoint_interval")
pickup && (kw[:overwrite_output] = false)
# --process_rates true|false (Covert only; default: on for P3) and --process_rate_interval in seconds
haskey(opts, "process_rates") && (kw[:process_rates] = parse(Bool, opts["process_rates"]))
haskey(opts, "process_rate_interval") && (kw[:process_rate_interval] = parse(Float64, opts["process_rate_interval"]))
haskey(opts, "Lx") && (kw[:Lx] = parse(Float64, opts["Lx"]))
haskey(opts, "Ly") && (kw[:Ly] = parse(Float64, opts["Ly"]))
haskey(opts, "hours") && (kw[:stop_time] = parse(Float64, opts["hours"]) * 3600)
haskey(opts, "cfl") && (kw[:cfl] = parse(Float64, opts["cfl"]))
haskey(opts, "dt") && (kw[:Δt] = parse(Float64, opts["dt"]))
haskey(opts, "max_dt") && (kw[:max_Δt] = parse(Float64, opts["max_dt"]))
haskey(opts, "initialization") && (kw[:initialization] = Symbol(opts["initialization"]))
haskey(opts, "profile_interval") && (kw[:profile_interval] = parse(Float64, opts["profile_interval"]))
haskey(opts, "slice_interval") && (kw[:slice_interval] = parse(Float64, opts["slice_interval"]))
if haskey(opts, "closure") || haskey(opts, "prandtl")
    closure_name = getopt("closure", "smagorinsky_lilly")
    closure_name in ("none", "smagorinsky_lilly", "moist_smagorinsky_lilly", "dynamic_smagorinsky") ||
        throw(ArgumentError("--closure must be none, smagorinsky_lilly, moist_smagorinsky_lilly or dynamic_smagorinsky"))
    closure_name == "dynamic_smagorinsky" && haskey(opts, "prandtl") && error("--prandtl applies to the Smagorinsky–Lilly closures")
    closure_name == "none" && haskey(opts, "prandtl") && error("--prandtl needs a Smagorinsky–Lilly --closure")
    # --prandtl: the turbulent Prandtl number of Smagorinsky–Lilly (Oceananigans' default 1; SAM uses 1/3).
    # moist_smagorinsky_lilly: Breeze's MoistSmagorinskyLilly, the Lilly stability reduction with the saturated N² in cloud
    Pr = parse(Float64, getopt("prandtl", "1"))
    # dynamic_smagorinsky: Oceananigans' DynamicSmagorinsky, Lagrangian averaging by default;
    # --dynamic_averaging xy averages the coefficient over horizontal planes instead
    averaging = getopt("dynamic_averaging", "lagrangian")
    averaging in ("lagrangian", "xy") || error("--dynamic_averaging must be lagrangian or xy")
    kw[:closure] = closure_name == "none" ? nothing :
                   closure_name == "smagorinsky_lilly" ? SmagorinskyLilly(; Pr) :
                   closure_name == "moist_smagorinsky_lilly" ? MoistSmagorinskyLilly(; Pr) :
                   averaging == "xy" ? DynamicSmagorinsky(; averaging = (1, 2)) : DynamicSmagorinsky()
end
# --lw_F0: the cloud-top longwave flux parameter F₀ of SAM rad_simple [W m⁻²] (Covert default: LASSO's 113; SAM/DYCOMS 70)
haskey(opts, "lw_F0") && (kw[:simple_longwave] = (; F₀ = parse(Float64, opts["lw_F0"])))
# --advection_order: WENO order of momentum (default 5; the scalars follow it unless --scalar_advection_order is given)
haskey(opts, "advection_order") && (kw[:advection_order] = parse(Int, opts["advection_order"]))
# --microphysics_advection_order: WENO order of the microphysical prognostics (default: the scalar order)
haskey(opts, "microphysics_advection_order") && (kw[:microphysics_advection_order] = parse(Int, opts["microphysics_advection_order"]))
# --scalar_advection_order: WENO order of the scalars (default: that of momentum, 5)
haskey(opts, "scalar_advection_order") && (kw[:scalar_advection_order] = parse(Int, opts["scalar_advection_order"]))
# --reference_extension z₁,z₂: hold the anelastic reference profile at its z₁ values up to z₂ (Covert only)
haskey(opts, "reference_extension") && (kw[:reference_extension] = Tuple(parse.(Float64, split(opts["reference_extension"], ','))))
if !isnothing(z_faces)
    haskey(opts, "Nz") && error("--grid and --Nz are exclusive")
    pop!(kw, :Nz, nothing)          # --grid replaces the vertical grid of --resolution
    kw[:z_faces] = z_faces
end
# Output files keep the campaign prefix that the analysis scripts read (`lasso_ena_profiles.jld2`, ...).
output_prefix = "lasso_ena"
cap = getopt("aerosol_ss_cap", "nothing") == "nothing" ? nothing : parse(Float64, opts["aerosol_ss_cap"])

if protocol === :covert_public_bin
    for option in ("member", "epoch", "dimensions", "nudging", "surface_flux_law", "liquid_reff")
        haskey(opts, option) && error("--$option applies to --protocol lasso_ena_official only")
    end
    member_name = getopt("microphysics", "p3_n75")
    faces = something(z_faces, covert_public_bin_vertical_faces(; Nz = get(kw, :Nz, 192)))
    microphysics = microphysics_member(member_name; data, z_faces = faces, cap,
                                      process_rates = p3_process_rates(get(opts, "p3_scale", nothing)))
    output_dir = getopt("output", "output/covert_public_bin_$(getopt("microphysics", "p3_n75"))")
    @info "Building the Covert case on $(typeof(arch)) with $(Oceananigans.defaults.FloatType)" opts
    # --radiation simple (protocol rad_simple, default) | rrtmgp_longwave (RRTMGP LW-only sensitivity)
    #             | rrtmgp (RRTMGP LW+SW with the diurnal cycle of the apparent solar position)
    if haskey(opts, "radiation")
        opts["radiation"] in ("simple", "rrtmgp_longwave", "rrtmgp") ||
            error("--radiation must be simple, rrtmgp_longwave or rrtmgp with the Covert protocol")
        kw[:radiation] = opts["radiation"] == "simple" ? :simple_longwave : Symbol(opts["radiation"])
    end
    # --nudging_tq τ: domain-mean θ and moisture nudged toward the snd records with timescale τ [s],
    #   weighted by a sin² ramp between the heights of --nudging_z z₁,z₂ (default 1200,1800 m: free troposphere only)
    # --nudging_uv τ: domain-mean u, v nudged toward the snd winds over the full column, timescale τ [s]
    haskey(opts, "nudging_tq") && (kw[:nudging_tq] = parse(Float64, opts["nudging_tq"]))
    haskey(opts, "nudging_uv") && (kw[:nudging_uv] = parse(Float64, opts["nudging_uv"]))
    if haskey(opts, "nudging_z")
        haskey(opts, "nudging_tq") || error("--nudging_z needs --nudging_tq")
        heights = Tuple(parse.(Float64, split(opts["nudging_z"], ',')))
        length(heights) == 2 || error("--nudging_z must be z₁,z₂")
        kw[:nudging_heights] = heights
    end
    seed = parse(Int, getopt("seed", "1234"))
    case = ena_covert(; arch, data_dir = data, microphysics, output_dir, output_prefix,
                        perturbation = InitialPerturbation(; seed), kw...)
else
    for option in ("process_rates", "process_rate_interval", "pickup", "scalar_advection_order",
                   "microphysics_advection_order", "reference_extension", "lw_F0", "nudging_tq", "nudging_z", "nudging_uv", "p3_scale")
        haskey(opts, option) && error("--$option applies to --protocol covert_public_bin only")
    end
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
    # --liquid_reff: RRTMGP liquid effective radius [m] in place of the bundle rule (compute_reffc ? 10 : 14 μm);
    # recorded in config.overrides and the label
    haskey(opts, "liquid_reff") && (kw[:liquid_effective_radius] = parse(Float64, opts["liquid_reff"]))
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
provenance_file = pickup ? "provenance_pickup_$(Dates.format(Dates.now(), "yyyymmddTHHMMSS")).toml" : "provenance.toml"
provenance = write_provenance(joinpath(output_dir, provenance_file), case;
                              extra = (; command = join(ARGS, " "), hostname = gethostname(), vertical_grid,
                                         microphysics_member = member_name, aerosol_supersaturation_cap = something(cap, 0.0),
                                         gpu = CUDA.functional() ? CUDA.name(CUDA.device()) : "none"))
@info "Provenance written to $provenance"
@info "Label: $(case.config.label)"
println(case.model)

wall = time()
run!(case.simulation; pickup)
@info "Finished in $(round((time() - wall) / 60, digits=1)) minutes; output in $output_dir"
