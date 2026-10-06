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
    ena_protocol_settings(data_dir; protocol, dimensions=nothing)

Read and validate protocol inputs without allocating a model. Return namelist-derived
defaults. If a bundle omits dimensions from its namelist, supply
`dimensions=(Nx, Ny, Nz)` from the original domain configuration. Neither path is a numerically validated reproduction of the original SAM runs.

- `:covert_public_bin` — the runnable *development benchmark* from the public files of
  the Covert et al. (2022) bin-paper repository: `caseid = 256x256x192` at `dx = dy = 35 m`
  (8.96 km), `day0 = 199.25` (18 July 2017 06 UTC), `nstop × dt = 6 h`, prescribed
  H/LE/τ (`SFC_FLX_FXD`, `SFC_TAU_FXD`), SAM `rad_simple` longwave only, no wind nudging,
  `doupperbound`, `dodamping`. Not the 864²×192, 30.24-km, 06-15 UTC configuration of the
  published paper, and **not** an official LASSO-ENA reproduction.
- `:lasso_ena_official` — the official protocol from a `samin` bundle directory
  (`snd`, `lsf`, `sfc`, `prm`, `grd`): grid from `grd` and the namelist, bulk fluxes from
  SST, RRTMGP LW+SW, `tauls` wind nudging (n0), P3-aer2 with the SBM `diagCCN`
  aerosol-reservoir projection, duration `nstop × dt`.
"""
function ena_protocol_settings(data_dir; protocol, dimensions=nothing)
    preset = protocol
    preset in (:covert_public_bin, :lasso_ena_official) || throw(ArgumentError("unknown ENA protocol $protocol"))
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
        grd_path = joinpath(data_dir, "grd")
        isfile(grd_path) || throw(ArgumentError("preset $preset needs the vertical grid file $(grd_path) from the samin bundle"))
        require_namelist!(namelist, ["dx", "dy", "dt", "nstop", "day0", "latitude0", "longitude0",
                                     "sfc_flx_fxd", "ocean", "dolargescale", "dosfcforcing", "donudging_uv", "tauls",
                                     "doupperbound", "dodamping", "read_in_geostrophic_wind", "nrad",
                                     "doradsimple", "dolongwave", "doshortwave"], preset)
        require_namelist_value!(namelist, "sfc_flx_fxd", false, preset)
        require_namelist_value!(namelist, "ocean", true, preset)
        require_namelist_value!(namelist, "dolargescale", true, preset)
        require_namelist_value!(namelist, "dosfcforcing", true, preset)
        require_namelist_value!(namelist, "donudging_uv", true, preset)
        require_namelist_value!(namelist, "doradsimple", false, preset)
        require_namelist_value!(namelist, "dolongwave", true, preset)
        require_namelist_value!(namelist, "doshortwave", true, preset)
        centers = read_sam_grd(grd_path)
        lsf_columns = read_sam_large_scale_forcing(joinpath(data_dir, "lsf"))[1].has_geostrophic_columns
        Bool(namelist["read_in_geostrophic_wind"]) == lsf_columns ||
            throw(ArgumentError("namelist read_in_geostrophic_wind = $(namelist["read_in_geostrophic_wind"]) but the lsf file " *
                                (lsf_columns ? "has" : "lacks") * " the ug/vg columns (with .false. the LASSO SAM aliases ug/vg to uls/vls, which the reader reproduces)"))
        sam_ug = Float64(get(namelist, "ug", 0.0))
        sam_vg = Float64(get(namelist, "vg", 0.0))
        if all(k -> haskey(namelist, k), ("nx_gl", "ny_gl", "nz_gl"))
            nx, ny, nz = Int.((namelist["nx_gl"], namelist["ny_gl"], namelist["nz_gl"]))
        elseif haskey(namelist, "caseid") && occursin(r"^\d+x\d+x\d+$", lowercase(namelist["caseid"]))
            nx, ny, nz = parse_caseid(namelist["caseid"])
        elseif dimensions isa NTuple{3, Integer}
            nx, ny, nz = dimensions
        else
            throw(ArgumentError("LASSO bundle must specify nx_gl/ny_gl/nz_gl or a numeric Nx×Ny×Nz caseid; otherwise supply dimensions=(Nx, Ny, Nz) from its domain configuration"))
        end
        isnothing(dimensions) || dimensions == (nx, ny, nz) ||
            throw(ArgumentError("supplied dimensions disagree with the namelist; use model overrides for a sensitivity experiment"))
        nx > 0 && ny > 0 && 2 ≤ nz ≤ length(centers) ||
            throw(ArgumentError("invalid LASSO dimensions ($nx, $ny, $nz) for $(length(centers)) grid levels"))
        all(isfinite, centers) && all(>(0), centers) && all(>(0), diff(centers)) ||
            throw(ArgumentError("LASSO grid centers must be finite, positive and strictly increasing"))
        defaults = (; label = string("LASSO-ENA official protocol (samin bundle; RRTMGP columns end at the LES top: Breeze analog of RRTMG, not padded to TOA; ",
                                     "namelist translation frame (ug, vg) = ($sam_ug, $sam_vg) NOT applied because bulk fluxes need the ground-relative wind)"),
                      Nx = nx, Ny = ny, Lx = nx * namelist["dx"], Ly = ny * namelist["dy"],
                      z_faces = faces_from_centers(centers[1:nz]),
                      day0 = Float64(namelist["day0"]),
                      latitude = Float64(namelist["latitude0"]), longitude = Float64(namelist["longitude0"]),
                      stop_time = Float64(namelist["nstop"] * namelist["dt"]),
                      Δt = Float64(namelist["dt"]),
                      max_Δt = Float64(namelist["dt"]),
                      translation_velocity = (0.0, 0.0), # bulk fluxes use the model wind; a nonzero frame is refused
                      microphysics = :p3_aer2,
                      radiation = :rrtmgp,
                      radiation_interval = Float64(namelist["nrad"] * namelist["dt"]),
                      surface = :bulk_sst,
                      wind_nudging_timescale = Float64(namelist["tauls"]),
                      vertical_advection = :full_field,
                      upper_boundary_relaxation = Bool(namelist["doupperbound"]),
                      sponge = Bool(namelist["dodamping"]) ? SAMSponge() : nothing,
                      aerosol_replenishment = :diagnostic_ccn)
    else
        throw(ArgumentError("unknown preset $preset (:covert_public_bin or :lasso_ena_official)"))
    end

    return defaults
end

"""
    ena_simulation(data_dir; protocol, dimensions=nothing, kwargs...)

Build an ENA experiment with an explicit protocol (`:covert_public_bin` or
`:lasso_ena_official`). Overrides create a labelled sensitivity experiment and are
recorded separately from the protocol defaults. For LASSO, supply `epoch::DateTime`
from the selected ARM bundle (UTC); the year is never inferred from the day of year.
This is a Breeze implementation of the forcing protocol, not a claim of SAM equivalence.
"""
function ena_simulation(data_dir; protocol, dimensions=nothing, kwargs...)
    defaults = ena_protocol_settings(data_dir; protocol, dimensions)
    overrides = Dict{Symbol, Any}(kwargs)
    if protocol === :lasso_ena_official
        epoch = get(overrides, :epoch, nothing)
        epoch isa Dates.DateTime || throw(ArgumentError("LASSO requires epoch=DateTime(...) from the ARM bundle in UTC"))
        day0 = get(overrides, :day0, defaults.day0)
        expected = epoch_from_day_of_year(day0; year=Dates.year(epoch))
        abs(Dates.value(epoch - expected)) ≤ 1000 ||
            throw(ArgumentError("epoch $epoch disagrees with namelist day0=$day0 (expected $expected)"))
    end
    for (name, value) in overrides
        if haskey(defaults, name) && defaults[name] != value
            @warn "protocol $protocol: overriding $name = $(defaults[name]) with $value (recorded in provenance)"
        end
    end
    settings = merge(defaults, NamedTuple(overrides))
    settings = merge(settings, (; label = string(settings.label, isempty(overrides) ? "" : " [overrides: $(join(sort!(string.(keys(overrides))), ", "))]")))
    case = build_case(data_dir; settings...)
    return merge(case, (; preset=protocol, protocol, protocol_dimensions=dimensions, protocol_overrides=NamedTuple(overrides)))
end

# Compatibility with the imported BreezyLASSO scripts. New callers should select
# the protocol explicitly with ena_simulation.
lasso_ena_simulation(data_dir; preset, kwargs...) =
    ena_simulation(data_dir; protocol=preset, kwargs...)
