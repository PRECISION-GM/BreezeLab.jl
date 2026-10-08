# ENA protocol definitions; shared model assembly lives in case_setup.jl.


"""
    read_sam_grd(path)

Read a SAM `grd` file (one scalar-level height per line, optionally more lines than the
model uses) and return the cell-center heights.
"""
function read_sam_grd(path)
    values = Float64[]
    for line in eachline(path)
        t = split(line)
        isempty(t) && continue
        v = tryparse(Float64, replace(t[1], r"[dD]" => "e"))
        isnothing(v) || push!(values, v)
    end
    isempty(values) && error("no heights found in $path")
    return values
end

"""
    faces_from_centers(centers)

Cell interfaces halfway between consecutive centers (bottom face at 0), for a SAM `grd`.
"""
function faces_from_centers(centers)
    faces = zeros(length(centers) + 1)
    for k in 1:length(centers)-1
        faces[k+1] = (centers[k] + centers[k+1]) / 2
    end
    faces[end] = centers[end] + (centers[end] - faces[end-1])
    return faces
end

"""
    covert_public_bin_vertical_faces(; Nz=192, top=20000, Δz=10, uniform_top=1500)

A labelled **reconstruction** of the 192-level, 20-km-top vertical grid of the Covert et al.
(2022) SAM runs (the public repository advertises but does not ship its `grd`/`domain.f90`):
uniform `Δz = 10 m` from the surface to `uniform_top` (150 cells, covering the surface layer
and the inversion), then the remaining 42 cells stretched with a constant growth ratio so
that the top face lands exactly at `top`. Exactly `Nz + 1` faces are returned; pass the
archive's `grd` (via `faces_from_centers(read_sam_grd(path))`) once it is available.
"""
function covert_public_bin_vertical_faces(; Nz=192, top=20000.0, Δz=10.0, uniform_top=1500.0)
    N_uniform = round(Int, uniform_top / Δz)
    N_stretch = Nz - N_uniform
    N_stretch ≥ 1 || throw(ArgumentError("Nz = $Nz leaves no stretched cells above $uniform_top m"))
    faces = collect(0.0:Δz:uniform_top)
    r = geometric_growth_ratio((top - uniform_top) / Δz, N_stretch)
    for n in 1:N_stretch
        push!(faces, faces[end] + Δz * r^n)
    end
    faces[end] = top
    length(faces) == Nz + 1 || error("vertical grid construction produced $(length(faces) - 1) cells, expected $Nz")
    return faces
end

"""
    covert_inversion_refined_vertical_faces(; Δz=10.0, Δz_fine=5.0, fine_bottom=800.0, fine_top=1400.0,
                                             uniform_top=1500.0, top=20000.0, N_stretch=42)

The [`covert_public_bin_vertical_faces`](@ref) reconstruction with the inversion layer
refined to `Δz_fine` between `fine_bottom` and `fine_top` (5 m over 800–1400 m, the spacing
Covert et al. (2022) use "near the surface and inversion layer"), `Δz` elsewhere below
`uniform_top`, and the same `N_stretch` geometrically stretched cells to `top` as the
reference grid. With the defaults this gives 80 + 120 + 10 + 42 = 252 cells (the reference
has 192). A labelled sensitivity grid for the cloud-top/inversion bias test, not the
paper's grid (which is 5 m near the surface as well and 192 levels in total).
"""
function covert_inversion_refined_vertical_faces(; Δz=10.0, Δz_fine=5.0, fine_bottom=800.0, fine_top=1400.0,
                                                   uniform_top=1500.0, top=20000.0, N_stretch=42)
    0 < fine_bottom < fine_top ≤ uniform_top < top || throw(ArgumentError("need 0 < fine_bottom < fine_top ≤ uniform_top < top"))
    faces = collect(0.0:Δz:fine_bottom)
    append!(faces, (fine_bottom + Δz_fine):Δz_fine:fine_top)
    fine_top < uniform_top && append!(faces, (fine_top + Δz):Δz:uniform_top)
    all(f -> isinteger(round(f; digits=9)) || true, faces)
    r = geometric_growth_ratio((top - uniform_top) / Δz, N_stretch)
    for n in 1:N_stretch
        push!(faces, faces[end] + Δz * r^n)
    end
    faces[end] = top
    issorted(faces) && all(>(0), diff(faces)) || error("refined vertical grid is not strictly increasing")
    return faces
end

"""
    ena_vertical_faces(name)

Named vertical grids of the ENA cases: `:covert` ([`covert_public_bin_vertical_faces`](@ref)),
`:covert_inversion_5m` ([`covert_inversion_refined_vertical_faces`](@ref)) and `:lasso`
([`lasso_ena_vertical_faces`](@ref)). The name is what run scripts record in provenance.
"""
function ena_vertical_faces(name)
    name = Symbol(name)
    name === :covert && return covert_public_bin_vertical_faces()
    name === :covert_inversion_5m && return covert_inversion_refined_vertical_faces()
    name === :lasso && return lasso_ena_vertical_faces()
    throw(ArgumentError("unknown vertical grid $name (covert, covert_inversion_5m, lasso)"))
end

function require_namelist!(namelist, required, preset)
    missing_keys = [k for k in required if !haskey(namelist, k)]
    isempty(missing_keys) || throw(ArgumentError("preset $preset requires namelist parameters $(missing_keys) (found $(sort(collect(keys(namelist)))))"))
    return nothing
end

function require_namelist_value!(namelist, key, value, preset)
    namelist[key] == value || throw(ArgumentError("preset $preset requires $key = $value in the namelist, found $(namelist[key])"))
    return nothing
end

parse_caseid(caseid) = parse.(Int, split(lowercase(String(caseid)), "x"))

"""
    ena_protocol_settings(data_dir; protocol, dimensions=nothing, member=nothing, epoch=nothing)

Read and validate protocol inputs without allocating a model. Return namelist-derived
defaults. If a LASSO bundle omits dimensions from its namelist, supply
`dimensions=(Nx, Ny, Nz)` from the original domain configuration; `member` (the run ID)
and `epoch` (UTC) are checked against it (see [`inspect_lasso_bundle`](@ref)). Neither
path is a numerically validated reproduction of the original SAM runs.

- `:covert_public_bin` — the runnable *development benchmark* from the public files of
  the Covert et al. (2022) bin-paper repository: `caseid = 256x256x192` at `dx = dy = 35 m`
  (8.96 km), `day0 = 199.25` (18 July 2017 06 UTC), `nstop × dt = 6 h`, prescribed
  H/LE/τ (`SFC_FLX_FXD`, `SFC_TAU_FXD`), SAM `rad_simple` longwave only, no wind nudging,
  `doupperbound`, `dodamping`. Not the 864²×192, 30.24-km, 06-15 UTC configuration of the
  published paper, and **not** an official LASSO-ENA reproduction.
- `:lasso_ena_official` — the official protocol from a `samin` bundle directory
  (`snd`, `lsf`, `sfc`, `prm`, `grd`): grid from `grd` and the namelist, bulk fluxes from
  SST, RRTMGP LW+SW with SAM's 0.95 surface emissivity and 14 μm ocean effective radius,
  `tauls` wind nudging (n0), P3 with the member's aerosol setting and the SBM `diagCCN`
  aerosol-reservoir projection, `nrad × dt` radiation cadence, fixed `dt`, duration
  `nstop × dt`, and the namelist's `perturb_type`. Unsupported namelist settings (T/q
  nudging, delayed forcing, restarts, …) are rejected with the reason.
"""
function ena_protocol_settings(data_dir; protocol, dimensions=nothing, member=nothing, epoch=nothing)
    preset = protocol
    preset in (:covert_public_bin, :lasso_ena_official) || throw(ArgumentError("unknown ENA protocol $protocol"))
    preset === :lasso_ena_official || isnothing(member) ||
        throw(ArgumentError("member identifies a LASSO-ENA samin bundle; the $protocol protocol has no member"))
    required = preset === :lasso_ena_official ? ("snd", "lsf", "sfc", "prm", "grd") : ("snd", "lsf", "sfc", "prm")
    missing_files = filter(name -> !isfile(joinpath(data_dir, name)), required)
    isempty(missing_files) || throw(ArgumentError("ENA protocol $protocol requires $(join(missing_files, ", ")) in $data_dir; no alternate inputs will be substituted"))
    prm_path = joinpath(data_dir, "prm")
    isfile(prm_path) || throw(ArgumentError("preset $preset needs the namelist $(prm_path)"))
    namelist = read_sam_namelist(prm_path)

    if preset === :covert_public_bin
        isnothing(dimensions) || throw(ArgumentError("Covert dimensions are specified by its namelist; use Nx, Ny, z_faces overrides for sensitivities"))
        require_namelist!(namelist, ["caseid", "dx", "dy", "dt", "nstop", "day0", "latitude0", "longitude0",
                                     "sfc_flx_fxd", "sfc_tau_fxd", "doradsimple", "dolongwave", "doshortwave",
                                     "donudging_uv", "dolargescale", "dosfcforcing", "doupperbound", "dodamping"], preset)
        require_namelist_value!(namelist, "sfc_flx_fxd", true, preset)
        require_namelist_value!(namelist, "sfc_tau_fxd", true, preset)
        require_namelist_value!(namelist, "doradsimple", true, preset)
        require_namelist_value!(namelist, "dolongwave", true, preset)
        require_namelist_value!(namelist, "doshortwave", false, preset)
        require_namelist_value!(namelist, "donudging_uv", false, preset)
        require_namelist_value!(namelist, "dolargescale", true, preset)
        require_namelist_value!(namelist, "dosfcforcing", true, preset)
        nx, ny, nz = parse_caseid(namelist["caseid"])
        defaults = (; label = "Covert-public-bin development benchmark (not an official LASSO-ENA reproduction)",
                      Nx = nx, Ny = ny, Lx = nx * namelist["dx"], Ly = ny * namelist["dy"],
                      z_faces = covert_public_bin_vertical_faces(; Nz=nz),
                      day0 = Float64(namelist["day0"]),
                      latitude = Float64(namelist["latitude0"]), longitude = Float64(namelist["longitude0"]),
                      stop_time = Float64(namelist["nstop"] * namelist["dt"]),
                      Δt = Float64(namelist["dt"]),
                      max_Δt = Float64(namelist["dt"]),
                      translation_velocity = (Float64(get(namelist, "ug", 0.0)), Float64(get(namelist, "vg", 0.0))),
                      microphysics = :p3_n75,
                      radiation = :simple,
                      surface = :prescribed_fluxes,
                      wind_nudging_timescale = nothing,
                      vertical_advection = :full_field,
                      upper_boundary_relaxation = Bool(namelist["doupperbound"]),
                      sponge = Bool(namelist["dodamping"]) ? SAMSponge() : nothing,
                      aerosol_replenishment = nothing)
        nz == length(defaults.z_faces) - 1 ||
            throw(ArgumentError("Covert-public-bin preset: namelist caseid has $nz levels but the reconstructed vertical grid has $(length(defaults.z_faces) - 1)"))
    elseif preset === :lasso_ena_official
        # Every layout/namelist/grid/time-coverage check lives in lasso_bundle.jl; unsupported
        # settings are reported together with the reason, and nothing is inferred.
        bundle = validate_lasso_bundle(data_dir; member, dimensions, epoch)
        defaults = bundle.settings
    else
        throw(ArgumentError("unknown preset $preset (:covert_public_bin or :lasso_ena_official)"))
    end

    return defaults
end

"""
    ena_simulation(data_dir; protocol, dimensions=nothing, kwargs...)

Build an ENA experiment with an explicit protocol (`:covert_public_bin` or
`:lasso_ena_official`). Overrides create a labelled sensitivity experiment and are
recorded separately from the protocol defaults. For LASSO, supply the run ID `member`,
`epoch::DateTime` (UTC; the year is never inferred from the day of year) and, when the
namelist carries no dimensions, `dimensions=(Nx, Ny, Nz)`; the bundle is validated by
[`validate_lasso_bundle`](@ref) and its record is kept under `case.bundle` for provenance.
This is a Breeze implementation of the forcing protocol, not a claim of SAM equivalence.
"""
function ena_simulation(data_dir; protocol, dimensions=nothing, member=nothing, kwargs...)
    overrides = Dict{Symbol, Any}(kwargs)
    bundle = nothing
    if protocol === :lasso_ena_official
        epoch = get(overrides, :epoch, nothing)
        epoch isa Dates.DateTime || throw(ArgumentError("LASSO requires epoch=DateTime(...) from the ARM bundle in UTC"))
        isnothing(member) &&
            throw(ArgumentError("LASSO requires member=\"<run ID>\" (e.g. 20170718era5d25x100_sbmwrm-aer2-flxsst) naming the samin bundle; the identity is never inferred from a directory"))
        bundle = validate_lasso_bundle(data_dir; member, dimensions, epoch)
        defaults = bundle.settings
        haskey(overrides, :day0) && overrides[:day0] != defaults.day0 &&
            throw(ArgumentError("day0 is fixed by the LASSO namelist ($(defaults.day0)); it cannot be overridden"))
    else
        defaults = ena_protocol_settings(data_dir; protocol, dimensions, member)
    end
    for (name, value) in overrides
        if haskey(defaults, name) && defaults[name] != value
            @warn "protocol $protocol: overriding $name = $(defaults[name]) with $value (recorded in provenance)"
        end
    end
    settings = merge(defaults, NamedTuple(overrides))
    settings = merge(settings, (; label = string(settings.label, isempty(overrides) ? "" : " [overrides: $(join(sort!(string.(keys(overrides))), ", "))]")))
    case = build_case(data_dir; settings...)
    member_id = isnothing(bundle) ? "none" : bundle.member.id
    return merge(case, (; preset=protocol, protocol, protocol_dimensions=dimensions, protocol_member=member_id,
                          protocol_overrides=NamedTuple(overrides),
                          bundle = isnothing(bundle) ? "none" : lasso_bundle_record(bundle)))
end

# Compatibility with the imported BreezyLASSO scripts. New callers should select
# the protocol explicitly with ena_simulation.
lasso_ena_simulation(data_dir; preset, kwargs...) =
    ena_simulation(data_dir; protocol=preset, kwargs...)
