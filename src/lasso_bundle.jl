#####
##### LASSO-ENA `samin` bundle: member identity, layout, namelist, grid and time-coverage
##### validation (ENA-LASSO plan steps 1–3).
#####
##### Conventions were verified against the LASSO SAM source (lasso_sam_sbm, branch
##### lasso_ena_noice, commit 12d0244):
#####   * setparm.f90      — the PARAMETERS namelist and its parameters
#####   * SGS_TKE/sgs.f90  — the SGS_TKE namelist (`dosmagor`, default .true.)
#####   * setgrid.f90      — `grd` holds scalar-level heights; interfaces are midpoints,
#####                        missing lines are extrapolated with the last spacing
#####   * setforcing.f90 / setdata.f90 / forcing.f90 — snd/lsf/sfc layouts and their time
#####                        handling (sounding interpolated to `day0`; lsf/snd padded with
#####                        a constant record after the last time; sfc extrapolated linearly
#####                        outside its samples and requiring at least two samples)
#####   * main.f90         — `day = day0 + nstep*dt/86400`, `nstop` steps
#####   * MICRO_HUJISBM/micro_prm.f90 — `diagCCN` is a compile-time parameter (.true.)
#####
##### Nothing here substitutes inputs: a bundle that cannot be represented in Breeze is
##### rejected with the list of offending settings.
#####

using Dates: Dates, DateTime, Date

const LASSO_ENA_DOI = "10.5439/2572661"
const LASSO_ENA_BUNDLE_BROWSER = "https://lasso-ena.svcs.arm.gov/latest/bundle_browser.html"
const LASSO_SAM_REFERENCE = (; url = "https://code.arm.gov/lasso/lasso-ena-codes/lasso_sam_sbm.git",
                               branch = "lasso_ena_noice",
                               commit = "12d02446a2147388dc89d828e6e0553106abea0f")
const LASSO_BUNDLE_FILES = ("snd", "lsf", "sfc", "prm", "grd")

#####
##### Member identity
#####

"""
    LassoMember

The identity of one LASSO-ENA ensemble member parsed from its run ID, e.g.
`20170718era5d25x100_sbmwrm-aer2-flxsst`. Tokens follow the ARM Live listing of the
`enalasso` datastream (October 2026): date, forcing reanalysis (`era5`/`merra2`), optional
variant tokens between the forcing and the domain (`s1n0` on some members; `n0` is the
documented wind-only nudging option, `s1` is recorded verbatim because the public
methodology page does not define it), domain token `d<width>x<spacing>` (`d25x100`:
25.6 km with 256 columns at 100 m), microphysics (`sbmwrm` HUJI spectral bin, warm-cloud
branch; `morr` Morrison), aerosol setting (`aer1`–`aer3`) and surface flux method
(`flxsst`: fluxes computed online from the reanalysis SST).
"""
struct LassoMember
    id :: String
    date :: Date
    forcing :: Symbol
    variant :: String
    domain :: String
    nominal_width_km :: Int
    spacing_m :: Int
    microphysics :: Symbol
    aerosol :: Symbol
    surface :: Symbol
end

const LASSO_MEMBER_PATTERN = r"^(\d{8})(era5|merra2)((?:[a-z]\d+)*)(d(\d+)x(\d+))_(sbmwrm|morr)-(aer[123])-(flx[a-z]+)$"

"""
    parse_lasso_member(id)

Parse a LASSO-ENA run ID into a [`LassoMember`](@ref). The ID is never guessed: an ID that
does not follow the listing's pattern is an error.
"""
function parse_lasso_member(id::AbstractString)
    m = match(LASSO_MEMBER_PATTERN, String(id))
    isnothing(m) && throw(ArgumentError("LASSO member ID $(repr(String(id))) does not follow <yyyymmdd><era5|merra2>[variant]d<km>x<m>_<sbmwrm|morr>-aer<1-3>-flx<...>; " *
                                        "use the run ID exactly as the ARM listing / bundle browser spells it"))
    date = Date(m.captures[1], "yyyymmdd")
    return LassoMember(String(id), date, Symbol(m.captures[2]), String(m.captures[3]), String(m.captures[4]),
                       parse(Int, m.captures[5]), parse(Int, m.captures[6]),
                       Symbol(m.captures[7]), Symbol(m.captures[8]), Symbol(m.captures[9]))
end

"""
    lasso_variant_tokens(member)

The variant tokens of a member ID (`["s1", "n0"]` for `…era5s1n0d25x100…`, empty for the
plain `…era5d25x100…` members).
"""
lasso_variant_tokens(member::LassoMember) = [String(t.match) for t in eachmatch(r"[a-z]\d+", member.variant)]

lasso_date_string(member::LassoMember) = Dates.format(member.date, "yyyymmdd")

"""
    lasso_samin_filename(member)

The ARM archive name of the member's `samin` input bundle.
"""
lasso_samin_filename(member::LassoMember) = "enalasso_samin_$(member.id)C1.m0.$(lasso_date_string(member)).000000.tar"

"""
    lasso_reference_filenames(member)

The per-member reference outputs: `samstat` (domain statistics) and `sam2d` (2-D fields).
"""
lasso_reference_filenames(member::LassoMember) =
    (; samstat = "enalasso_samstat_$(member.id)C1.m1.$(lasso_date_string(member)).000000.nc",
       sam2d = "enalasso_sam2d_$(member.id)C1.m1.$(lasso_date_string(member)).000000.nc")

lasso_member_record(member::LassoMember) =
    (; id = member.id, date = string(member.date), forcing = string(member.forcing), variant = member.variant,
       variant_tokens = lasso_variant_tokens(member), domain = member.domain,
       nominal_width_km = member.nominal_width_km, spacing_m = member.spacing_m,
       sam_microphysics = string(member.microphysics), aerosol = string(member.aerosol), surface = string(member.surface),
       samin = lasso_samin_filename(member), lasso_reference_filenames(member)..., doi = LASSO_ENA_DOI)

#####
##### Fortran namelist groups
#####

"""
    read_sam_namelist_groups(path)

Parse a SAM `prm` file into `Dict{String, Dict{String, NamelistValue}}` keyed by lower-case
group name (`parameters`, `sgs_tke`, `micro_hujisbm`, …). Values use the same conversions as
[`read_sam_namelist`](@ref); a key repeated inside one group keeps its last value.
"""
function read_sam_namelist_groups(path)
    groups = Dict{String, Dict{String, NamelistValue}}()
    current = nothing
    body = String[]
    flush!(group, lines) = begin
        isnothing(group) && return
        params = get!(groups, group, Dict{String, NamelistValue}())
        for m in eachmatch(r"([A-Za-z_][A-Za-z0-9_]*)\s*=\s*('[^']*'|\"[^\"]*\"|[^,\s]+)", join(lines, " "))
            params[lowercase(m.captures[1])] = parse_namelist_value(m.captures[2])
        end
    end
    for raw in split(read(path, String), '\n')
        line = strip(first(split(raw, '!')))
        isempty(line) && continue
        if startswith(line, "&") || startswith(line, "\$")
            flush!(current, body)
            name = lowercase(String(strip(line[2:end])))
            name in ("end",) && (current = nothing; empty!(body); continue)
            current = name
            empty!(body)
        elseif line == "/" || lowercase(line) in ("\$end", "&end")
            flush!(current, body)
            current = nothing
            empty!(body)
        else
            push!(body, line)
        end
    end
    flush!(current, body)
    return groups
end

#####
##### Grid (setgrid.f90)
#####

"""
    lasso_scalar_levels(grd_values, nz)

The `nz` scalar-level heights SAM uses from a `grd` file: the first `nz` lines, or, when the
file has fewer lines, the listed levels extended with the last spacing
(`setgrid.f90`: `z(k) = z(k-1) + (z(k-1) - z(k-2))`). Returns `(levels, extrapolated)`.
"""
function lasso_scalar_levels(grd_values, nz)
    n = length(grd_values)
    n ≥ 2 || throw(ArgumentError("grd needs at least two levels, found $n"))
    if n ≥ nz
        return grd_values[1:nz], 0
    end
    levels = copy(grd_values)
    while length(levels) < nz
        push!(levels, levels[end] + (levels[end] - levels[end-1]))
    end
    return levels, nz - n
end

#####
##### Time coverage helpers
#####

day_range(records) = (first(records).day, last(records).day)

#####
##### Inspection
#####

"""
    inspect_lasso_bundle(directory; member=nothing, dimensions=nothing, epoch=nothing)

Read a LASSO-ENA `samin` bundle (`snd`, `lsf`, `sfc`, `prm`, `grd` in `directory`) and
return a record describing what it contains, what the Breeze adapter will use, and every
setting the adapter cannot represent. Nothing is substituted or inferred:

- `member` (run ID or [`LassoMember`](@ref)) fixes the member identity; it is required by
  [`ena_lasso`](@ref) and recorded in provenance. The aerosol setting and the SAM
  microphysics come from it (the namelist does not carry them: the HUJI-SBM spectrum is
  compiled in).
- `dimensions = (Nx, Ny, Nz)` is required when the namelist carries neither `nx_gl/ny_gl/nz_gl`
  nor a numeric `NxxNyxNz` case ID (SAM fixes them at compile time in `domain.f90`).
- `epoch::DateTime` (UTC) is checked against `day0` and the member date when given; SAM's
  `day0` never identifies the year.

The returned named tuple has `problems` (fatal, each an actionable sentence), `warnings`
(deliberate, recorded differences), `settings` (keyword defaults for `build_case`),
`namelist`, `groups`, `grid`, `time`, `switches`, `checksums` and `member`. Use
[`validate_lasso_bundle`](@ref) to turn `problems` into an error.
"""
function inspect_lasso_bundle(directory; member=nothing, dimensions=nothing, epoch=nothing)
    problems = String[]
    warnings = String[]
    member = member isa AbstractString ? parse_lasso_member(member) : member
    member isa Union{Nothing, LassoMember} || throw(ArgumentError("member must be a run ID string or a LassoMember"))

    missing_files = filter(name -> !isfile(joinpath(directory, name)), LASSO_BUNDLE_FILES)
    isempty(missing_files) ||
        throw(ArgumentError("LASSO-ENA bundle directory $directory lacks $(join(missing_files, ", ")); " *
                            "stage the samin archive with data_wrangling/stage_lasso_bundle.jl; no alternate inputs are substituted"))
    files = Dict(name => joinpath(directory, name) for name in LASSO_BUNDLE_FILES)
    checksums = Dict(name => file_sha256(path) for (name, path) in files)

    groups = read_sam_namelist_groups(files["prm"])
    haskey(groups, "parameters") ||
        push!(problems, "prm has no &PARAMETERS namelist group (found $(sort(collect(keys(groups))))); SAM's setparm.f90 reads exactly that group")
    nml = get(groups, "parameters", Dict{String, NamelistValue}())
    sgs = get(groups, "sgs_tke", Dict{String, NamelistValue}())

    getnml(key, default) = get(nml, key, default)
    isset(key) = haskey(nml, key)

    # --- required parameters ---------------------------------------------------------
    required = ("dx", "dy", "dt", "nstop", "day0", "latitude0", "longitude0", "tauls", "nrad")
    missing_keys = [k for k in required if !isset(k)]
    isempty(missing_keys) ||
        push!(problems, "prm lacks required parameters $(missing_keys); the adapter never fills them in")

    # --- switch table: (key, required value, SAM default) ------------------------------
    # A key absent from the namelist takes SAM's default (params.f90 / grid.f90).
    expectations = (
        ("les", true, false, "the adapter implements the LES branch of surface.f90 (uniform fluxes from the domain-mean state)"),
        ("ocean", true, false, "bulk fluxes from SST need OCEAN = .true. (oceflx.f90)"),
        ("land", false, false, "land surfaces (landflx.f90) are not implemented"),
        ("cem", false, false, "CEM mode computes local per-column fluxes (surface.f90)"),
        ("sfc_flx_fxd", false, false, "the LASSO flxsst protocol computes H/LE online; prescribed fluxes are the Covert protocol"),
        ("sfc_tau_fxd", false, false, "stress must come from oceflx (SFC_TAU_FXD = .true. would read tau from sfc)"),
        ("uniform_sfc_flx", true, true, "per-column fluxes (UNIFORM_SFC_FLX = .false.) are not the uniform-flux LES branch the adapter implements"),
        ("dolargescale", true, false, "the lsf forcing must be active"),
        ("dosfcforcing", true, false, "the sfc SST series must be active"),
        ("donudging_uv", true, false, "wind nudging to uls/vls is part of the protocol (nudging.f90)"),
        ("donudging_tq", false, false, "temperature/moisture nudging (nudging options n1/n2) is not implemented in BreezeLab"),
        ("donudging_t", false, false, "temperature nudging is not implemented in BreezeLab"),
        ("donudging_q", false, false, "moisture nudging is not implemented in BreezeLab"),
        ("donudging_transient", false, false, "transient nudging is not implemented in BreezeLab"),
        ("dovariable_tauz", false, false, "inversion-following nudging heights are not implemented in BreezeLab"),
        ("doradsimple", false, false, "the protocol uses RRTMG; rad_simple is the Covert control"),
        ("dolongwave", true, false, "longwave radiation must be on"),
        ("doshortwave", true, false, "shortwave radiation must be on"),
        ("doradforcing", false, false, "a prescribed `rad` tendency file is not supported"),
        ("docloud", true, false, "cloud microphysics must be on"),
        ("doprecip", true, false, "precipitation must be on"),
        ("dosurface", true, false, "surface fluxes must be on"),
        ("docoriolis", true, false, "Coriolis must be on (the geostrophic forcing enters through it)"),
        ("dofplane", true, true, "a latitude-varying Coriolis parameter is not supported"),
        ("docoriolisz", false, false, "the vertical Coriolis component is not supported"),
        ("doconstdz", false, false, "the protocol grid comes from grd; a constant dz namelist grid is not the official grid"),
        ("doperpetual", false, false, "perpetual insolation is not supported"),
        ("dosolarconstant", false, false, "a fixed solar constant/zenith angle is not supported"),
        ("doradhomo", false, false, "horizontally homogenized radiation is not supported"),
        ("dossthomo", false, false, "homogenized SST is not applicable to the uniform SST series"),
        ("dodynamicocean", false, false, "slab/dynamic ocean is not supported"),
        ("doensemble", false, false, "subensemble perturbation files (tqpert) are not supported"),
        ("dotracers", false, false, "passive tracers are not supported"),
        ("dosmoke", false, false, "the smoke-cloud case is not supported"),
        ("docolumn", false, false, "single-column mode is not supported"),
        ("dowallx", false, false, "x walls are not supported (doubly periodic domain)"),
        ("dowally", false, false, "y walls are not supported (doubly periodic domain)"),
        ("doscamiopdata", false, false, "SCAM IOP netCDF forcing is not supported; the text snd/lsf/sfc files are"),
        ("notracegases", false, false, "removing trace gases is not supported"),
    )
    for (key, wanted, default, why) in expectations
        value = getnml(key, default)
        value isa Bool || (push!(problems, "prm $key = $value is not a logical"); continue)
        value == wanted ||
            push!(problems, "prm $key = $value" * (isset(key) ? "" : " (SAM default)") * " is unsupported: $why")
    end

    # --- numeric settings with constraints ---------------------------------------------
    timelargescale = Float64(getnml("timelargescale", 0.0))
    timelargescale == 0 ||
        push!(problems, "prm timelargescale = $timelargescale s delays the large-scale tendencies and subsidence (forcing.f90); a delayed start is not implemented, so this member cannot be reproduced")
    nxco2 = Float64(getnml("nxco2", 1.0))
    nxco2 == 1 || push!(problems, "prm nxco2 = $nxco2 scales the RRTMG CO₂ profile; only the unscaled profile is represented")
    nrestart = Int(getnml("nrestart", 0))
    nrestart == 0 || push!(problems, "prm nrestart = $nrestart is a restarted SAM run; the adapter starts from the sounding")
    perturb_type = Int(getnml("perturb_type", 0))
    perturb_type in (0, 5) ||
        push!(problems, "prm perturb_type = $perturb_type: only setperturb.f90 cases 0 (±0.02 K × (6-k) in the five lowest levels) and 5 (±0.1 K, ±0.025 g/kg below 600 m) are implemented")
    nudging_z1 = Float64(getnml("nudging_uv_z1", -1.0))
    nudging_z2 = Float64(getnml("nudging_uv_z2", 1.0e6))
    fcor = Float64(getnml("fcor", -999.0))
    dosmagor = Bool(get(sgs, "dosmagor", true))
    compute_reffc = Bool(getnml("compute_reffc", false))
    compute_reffi = Bool(getnml("compute_reffi", false))
    doseasons = Bool(getnml("doseasons", false))
    ug = Float64(getnml("ug", 0.0))
    vg = Float64(getnml("vg", 0.0))

    # --- dimensions --------------------------------------------------------------------
    nx = ny = nz = nothing
    if all(isset, ("nx_gl", "ny_gl", "nz_gl"))
        nx, ny, nz = Int.((nml["nx_gl"], nml["ny_gl"], nml["nz_gl"]))
        dimension_source = "namelist nx_gl/ny_gl/nz_gl"
    elseif isset("caseid") && nml["caseid"] isa String && occursin(r"^\d+x\d+x\d+$", lowercase(nml["caseid"]))
        nx, ny, nz = parse_caseid(nml["caseid"])
        dimension_source = "numeric caseid $(nml["caseid"])"
    elseif dimensions isa NTuple{3, Integer}
        nx, ny, nz = dimensions
        dimension_source = "explicit dimensions keyword (SAM domain.f90 is compiled in)"
    else
        push!(problems, "prm carries neither nx_gl/ny_gl/nz_gl nor a numeric NxxNyxNz caseid (SAM fixes the domain at compile time in domain.f90): " *
                        "supply dimensions=(Nx, Ny, Nz) from the member's domain configuration (d25x100 members: 256, 256, 260)")
        dimension_source = "missing"
    end
    if !isnothing(dimensions) && !isnothing(nx) && dimensions != (nx, ny, nz)
        push!(problems, "supplied dimensions $dimensions disagree with the $dimension_source ($nx, $ny, $nz); use the Nx/Ny/z_faces overrides for a sensitivity experiment instead")
    end
    if !isnothing(member) && !isnothing(nx)
        expected_nx = round(Int, 1000 * member.nominal_width_km / member.spacing_m)
        # d25x100 means 25.6 km / 100 m = 256 columns; the token truncates the width.
        nx ∈ (expected_nx, expected_nx + 6) ||
            push!(warnings, "member domain token $(member.domain) suggests ≈$expected_nx columns but the configuration has Nx = $nx")
    end

    # --- files ---------------------------------------------------------------------------
    soundings = read_sam_sounding(files["snd"])
    lsf = read_sam_large_scale_forcing(files["lsf"])
    sfc = read_sam_surface_forcing(files["sfc"])
    grd_values = read_sam_grd(files["grd"])

    dx = Float64(getnml("dx", NaN)); dy = Float64(getnml("dy", NaN)); dt = Float64(getnml("dt", NaN))
    nstop = Int(getnml("nstop", 0)); day0 = Float64(getnml("day0", NaN)); nrad = Int(getnml("nrad", 1))
    tauls = Float64(getnml("tauls", NaN))
    latitude = Float64(getnml("latitude0", NaN)); longitude = Float64(getnml("longitude0", NaN))
    if !isnothing(member) && member.spacing_m > 0 && isfinite(dx) && !(dx ≈ member.spacing_m && dy ≈ member.spacing_m)
        push!(warnings, "member token $(member.domain) states $(member.spacing_m) m spacing but prm has dx = $dx, dy = $dy")
    end
    dt > 0 || push!(problems, "prm dt = $dt must be positive")
    nstop > 0 || push!(problems, "prm nstop = $nstop must be positive")
    nrad ≥ 1 || push!(problems, "prm nrad = $nrad must be at least 1")
    isfinite(tauls) && tauls > 0 || push!(problems, "prm tauls = $tauls must be a positive nudging timescale")
    stop_time = nstop * dt
    day_end = day0 + stop_time / 86400

    if fcor != -999.0
        expected_f = 4π / 86400 * sind(latitude)
        isapprox(fcor, expected_f; rtol=1e-3) ||
            push!(problems, "prm fcor = $fcor differs from 4π/86400 sin(latitude0) = $expected_f; a Coriolis parameter detached from the latitude is not supported")
    end

    # --- vertical grid (setgrid.f90) ------------------------------------------------------
    grid = nothing
    if !isnothing(nz)
        if nz < 2
            push!(problems, "Nz = $nz is too small")
        else
            all(isfinite, grd_values) && all(>(0), grd_values) && issorted(grd_values; lt=(<=)) ||
                push!(problems, "grd scalar levels must be finite, positive and strictly increasing")
            levels, extrapolated = try
                lasso_scalar_levels(grd_values, nz)
            catch err
                push!(problems, sprint(showerror, err)); (grd_values, 0)
            end
            extrapolated > 0 &&
                push!(warnings, "grd lists $(length(grd_values)) levels for Nz = $nz; SAM extends the last $extrapolated with the final spacing (setgrid.f90), reproduced here")
            length(grd_values) > nz + 1 &&
                push!(warnings, "grd lists $(length(grd_values)) levels; SAM reads only the first Nz + 1 = $(nz + 1) (the extra one sets nothing the adapter uses)")
            faces = faces_from_centers(levels)
            midpoints = (faces[1:end-1] .+ faces[2:end]) ./ 2
            offset = maximum(abs, midpoints .- levels)
            grid = (; levels, faces, midpoints, max_level_offset = offset, extrapolated,
                      first_level = levels[1], top_level = levels[end], top_face = faces[end],
                      uniform_spacing = levels[2] - levels[1])
            offset > 1e-6 &&
                push!(warnings, "SAM scalar levels sit up to $(round(offset; digits=3)) m away from the Breeze cell centers (midpoints of the SAM interfaces) in the stretched layers; forcing profiles are interpolated to the Breeze centers")
        end
    end

    # --- time coverage (setdata.f90 / forcing.f90) ---------------------------------------
    snd_days = day_range(soundings)
    lsf_days = day_range(lsf)
    sfc_days = (first(sfc.day), last(sfc.day))
    issorted([r.day for r in soundings]) || push!(problems, "snd records are not in increasing time order")
    issorted([r.day for r in lsf]) || push!(problems, "lsf records are not in increasing time order")
    issorted(sfc.day) || push!(problems, "sfc samples are not in increasing time order")
    tol = 1e-6
    if length(soundings) == 1
        push!(warnings, "snd holds a single record; SAM's setdata.f90 needs two records bracketing day0 and would abort, the adapter holds the profile constant")
    elseif !(snd_days[1] - tol ≤ day0 ≤ snd_days[2] + tol)
        push!(problems, "day0 = $day0 lies outside the snd records $(snd_days); SAM aborts ('day is beyond the sounding time range') and the adapter will not extrapolate")
    end
    lsf_days[1] - tol ≤ day0 ||
        push!(problems, "the first lsf record ($(lsf_days[1])) is later than day0 = $day0; SAM would extrapolate the forcing backwards in time while Breeze clamps it, so the files must start at or before day0")
    lsf_days[2] + tol < day_end &&
        push!(warnings, "lsf ends at day $(lsf_days[2]) before the run end $(round(day_end; digits=4)); SAM pads a constant record after the last time (setforcing.f90) and Breeze clamps, which agree")
    length(sfc.day) ≥ 2 || push!(problems, "sfc needs at least two samples (setforcing.f90 aborts otherwise)")
    (sfc_days[1] - tol ≤ day0 && sfc_days[2] + tol ≥ day_end) ||
        push!(problems, "sfc samples $(sfc_days) do not cover the run $(round(day0; digits=4))–$(round(day_end; digits=4)); SAM extrapolates the SST series linearly outside its samples while Breeze clamps, so the file must cover the run")

    # Geostrophic columns vs namelist (forcing.f90 aliases ug0/vg0 to uls/vls without them)
    read_geostrophic = Bool(getnml("read_in_geostrophic_wind", false))
    read_geostrophic == first(lsf).has_geostrophic_columns ||
        push!(problems, "prm read_in_geostrophic_wind = $read_geostrophic but the lsf file " *
                        (first(lsf).has_geostrophic_columns ? "has" : "lacks") * " the ug/vg columns (setforcing.f90 reads 9 or 7 columns accordingly)")

    # --- vertical coverage ---------------------------------------------------------------
    if !isnothing(grid)
        sounding_top = minimum(maximum(record_heights(r)) for r in soundings)
        sounding_top ≥ grid.top_level ||
            push!(problems, "the sounding reaches $(round(sounding_top; digits=1)) m but the model's top scalar level is $(grid.top_level) m; SAM continues with the 1976 standard atmosphere above the sounding (setdata.f90), which the adapter does not implement")
        lsf_top = minimum(maximum(record_heights(r, soundings[1])) for r in lsf)
        lsf_top < grid.top_level &&
            push!(warnings, "the lsf records end at $(round(lsf_top; digits=1)) m below the model top $(grid.top_level) m; above them SAM zeroes tls/qls/wls and holds the winds (forcing.f90), as the adapter does")
        nudging_z1 ≤ grid.first_level && nudging_z2 ≥ grid.top_level ||
            push!(problems, "prm nudging_uv_z1/z2 = ($nudging_z1, $nudging_z2) m restrict the wind nudging to a height window; only whole-column nudging is implemented")
    end

    # --- epoch / member date ---------------------------------------------------------------
    if !isnothing(epoch)
        epoch isa DateTime || push!(problems, "epoch must be a DateTime in UTC, got $(typeof(epoch))")
    end
    if !isnothing(epoch) && epoch isa DateTime && isfinite(day0)
        expected = epoch_from_day_of_year(day0; year=Dates.year(epoch))
        abs(Dates.value(epoch - expected)) ≤ 1000 ||
            push!(problems, "epoch $epoch disagrees with namelist day0 = $day0 (which is $expected in $(Dates.year(epoch)))")
        if !isnothing(member) && Date(expected) != member.date
            push!(problems, "member date $(member.date) and namelist day0 = $day0 ($expected) name different days")
        end
    end

    # --- deliberate differences recorded as warnings -------------------------------------
    push!(warnings, "SAM SGS: " * (dosmagor ? "SGS_TKE with dosmagor = .true. (diagnostic Smagorinsky)" : "SGS_TKE prognostic 1.5-order TKE") *
                    "; Breeze uses Smagorinsky–Lilly")
    push!(warnings, compute_reffc ?
          "compute_reffc = .true.: SAM diagnoses the liquid effective radius from the SBM spectrum (2.5–60 μm); Breeze prescribes a constant radius" :
          "compute_reffc = .false.: SAM's RRTMG uses computeRe_Liquid = 14 μm over ocean; the adapter prescribes 14 μm")
    compute_reffi && push!(warnings, "compute_reffi = .true. has no effect in the warm-cloud (noice) branch")
    doseasons || push!(warnings, "doseasons = .false.: SAM keeps the solar declination of day0's date (rad.f90 dayForSW); Breeze follows the calendar, a negligible difference over one day")
    (ug != 0 || vg != 0) &&
        push!(warnings, "namelist translation frame (ug, vg) = ($ug, $vg) m/s: SAM integrates winds relative to it and adds it back for the surface fluxes; the adapter integrates ground-relative winds (bulk fluxes need the absolute wind) and leaves the nudging/geostrophic targets unshifted, which is the same physical state")
    push!(warnings, "HUJI-SBM diagCCN is a compile-time parameter (.true. in micro_prm.f90): the adapter's DiagnosticCCNProjection reproduces the reservoir rule once per step")

    # --- settings for build_case -----------------------------------------------------------
    aerosol = isnothing(member) ? nothing : member.aerosol
    microphysics = isnothing(aerosol) ? :p3_aer2 : Symbol("p3_", aerosol)
    settings = isnothing(grid) || !isempty(problems) ? nothing :
        (; label = lasso_label(member, dimension_source, ug, vg),
           Nx = nx, Ny = ny, Lx = nx * dx, Ly = ny * dy,
           z_faces = grid.faces,
           day0, latitude, longitude,
           stop_time = Float64(stop_time), Δt = dt, max_Δt = dt,
           translation_velocity = (0.0, 0.0),
           microphysics,
           radiation = :rrtmgp,
           radiation_interval = Float64(nrad * dt),
           surface = :bulk_sst,
           surface_emissivity = 0.95,                       # RAD_RRTM/rad.f90: surfaceEmissivity = 0.95
           liquid_effective_radius = compute_reffc ? 10e-6 : 14e-6,   # cam_rad_parameterizations: rliqocean = 14 μm
           wind_nudging_timescale = tauls,
           vertical_advection = :full_field,
           upper_boundary_relaxation = Bool(getnml("doupperbound", false)),
           sponge = Bool(getnml("dodamping", false)) ? SAMSponge() : nothing,
           aerosol_replenishment = :diagnostic_ccn,
           perturbation = InitialPerturbation(; sam_perturb_type = perturb_type))

    switches = (; dosmagor, compute_reffc, compute_reffi, doseasons, perturb_type, timelargescale, nxco2,
                  read_in_geostrophic_wind = read_geostrophic, nudging_uv_z1 = nudging_z1, nudging_uv_z2 = nudging_z2,
                  doupperbound = Bool(getnml("doupperbound", false)), dodamping = Bool(getnml("dodamping", false)),
                  sam_translation = (ug, vg), fcor)
    time = (; day0, dt, nstop, stop_time, day_end, nrad, radiation_interval = nrad * dt, tauls,
              snd_days, lsf_days, sfc_days, snd_records = length(soundings), lsf_records = length(lsf), sfc_samples = length(sfc.day))

    return (; directory = abspath(directory), member, files, checksums, namelist = nml, groups, grid, time, switches,
              dimension_source, problems, warnings, settings)
end

function lasso_label(member, dimension_source, ug, vg)
    id = isnothing(member) ? "member not declared" : member.id
    return string("LASSO-ENA official protocol, member $id (samin bundle; dimensions from $dimension_source; ",
                  "RRTMGP columns end at the LES top where SAM's RRTMG adds one layer to TOA; ",
                  "namelist translation frame (ug, vg) = ($ug, $vg) not applied because bulk fluxes need the ground-relative wind)")
end

"""
    validate_lasso_bundle(directory; kwargs...)

[`inspect_lasso_bundle`](@ref) and throw an `ArgumentError` listing every unsupported or
missing setting when the bundle cannot be reproduced. Returns the inspection record.
"""
function validate_lasso_bundle(directory; kwargs...)
    bundle = inspect_lasso_bundle(directory; kwargs...)
    isempty(bundle.problems) ||
        throw(ArgumentError("LASSO-ENA bundle $(bundle.directory) cannot be reproduced as configured:\n  - " *
                            join(bundle.problems, "\n  - ")))
    return bundle
end

"""
    lasso_bundle_record(bundle)

TOML-ready provenance of an inspected bundle: member identity, file checksums, the SAM
reference revision, grid and time facts, switches and the recorded warnings.
"""
function lasso_bundle_record(bundle)
    grid = bundle.grid
    return (; directory = bundle.directory,
              member = isnothing(bundle.member) ? "undeclared" : lasso_member_record(bundle.member),
              doi = LASSO_ENA_DOI, bundle_browser = LASSO_ENA_BUNDLE_BROWSER,
              sam_reference = LASSO_SAM_REFERENCE,
              checksums = NamedTuple{Tuple(Symbol.(sort(collect(keys(bundle.checksums)))))}(Tuple(bundle.checksums[k] for k in sort(collect(keys(bundle.checksums))))),
              dimension_source = bundle.dimension_source,
              grid = isnothing(grid) ? "invalid" :
                     (; Nz = length(grid.levels), first_level = grid.first_level, top_level = grid.top_level, top_face = grid.top_face,
                        uniform_spacing = grid.uniform_spacing, extrapolated_levels = grid.extrapolated, max_level_offset = grid.max_level_offset),
              time = bundle.time, switches = bundle.switches, warnings = bundle.warnings)
end
