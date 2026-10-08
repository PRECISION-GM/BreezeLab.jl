#####
##### The official LASSO-ENA member constructor: a thin, explicit wrapper over
##### `ena_simulation(...; protocol=:lasso_ena_official)` that names the member, locates its
##### staged bundle, fails fast with staging instructions when it is absent, and carries the
##### staging provenance into the case record.
#####

using TOML: TOML

"""
The member adopted after the authenticated ARM listing of 7 October 2026: the plan's
candidate `20170718era5s1n0d25x100_sbmwrm-aer2-flxsst` does not exist (the `s1n0` members
were run with Morrison microphysics only); the closest obtainable spectral-bin warm-cloud
SST-flux member is this one. Changing it is an explicit decision, never a default.
"""
const DEFAULT_LASSO_MEMBER = "20170718era5d25x100_sbmwrm-aer2-flxsst"

"""
    lasso_bundle_directory(member)

Where the staged bundle of `member` is expected: `\$ENA_LASSO_BUNDLE` when set, otherwise
`data/lasso/<member>` in the package checkout (populated by
`data_wrangling/stage_lasso_bundle.jl`).
"""
function lasso_bundle_directory(member)
    id = member isa LassoMember ? member.id : parse_lasso_member(member).id
    return get(ENV, "ENA_LASSO_BUNDLE", joinpath(dirname(@__DIR__), "data", "lasso", id))
end

"""
    lasso_bundle_available(directory)

`true` when `snd`, `lsf`, `sfc`, `prm` and `grd` are all present in `directory`.
"""
lasso_bundle_available(directory) = all(isfile(joinpath(directory, f)) for f in LASSO_BUNDLE_FILES)

"""
    lasso_bundle_missing_message(member, directory)

The actionable message printed when a member's bundle is not staged: which files are
missing, the exact ARM archive name, the DOI and the staging command. It never suggests
substituting other inputs.
"""
function lasso_bundle_missing_message(member, directory)
    m = member isa LassoMember ? member : parse_lasso_member(member)
    missing_files = [f for f in LASSO_BUNDLE_FILES if !isfile(joinpath(directory, f))]
    archive = lasso_samin_filename(m)
    return string("LASSO-ENA bundle for member $(m.id) is not staged: $directory lacks ",
                  join(missing_files, ", "), ".\n",
                  "  Obtain the ARM archive $archive (DOI $LASSO_ENA_DOI; it must be staged online by ARM — ",
                  "the ARM Live `saveData` route returned HTTP 404 for it on 2026-10-07 although the listing shows it), then run\n",
                  "    julia data_wrangling/stage_lasso_bundle.jl data/lasso/archives/$archive\n",
                  "  (set ENA_LASSO_BUNDLE to use another directory). The Covert public inputs are a different protocol and are never substituted.")
end

"""
    lasso_documented_dimensions(member)

The documented LASSO-ENA domain for a member's `d<width>x<spacing>` token (modeling
methodology, domain table): `d25x100` is 25.6 km with 256 columns at 100 m and 260 levels
to 8 km. Returns `(Nx, Ny, Nz)`, used only when the caller asks for the documented domain
explicitly; the bundle's `grd` must agree.
"""
function lasso_documented_dimensions(member)
    m = member isa LassoMember ? member : parse_lasso_member(member)
    m.domain == "d25x100" && return (256, 256, 260)
    throw(ArgumentError("no documented dimensions for domain token $(m.domain); pass dimensions=(Nx, Ny, Nz) explicitly"))
end

"""
    ena_lasso(bundle_dir=lasso_bundle_directory(member); member=DEFAULT_LASSO_MEMBER, epoch,
              dimensions=:documented, arch=nothing, FT=Float32, output_dir, output_prefix="ena_lasso", kwargs...)

Build the official LASSO-ENA member `member` from its staged `samin` bundle without running
it. Returns the [`ena_simulation`](@ref) case extended with `staging` (the `bundle.toml`
written by `data_wrangling/stage_lasso_bundle.jl`, or a note when absent) so that
[`write_provenance`](@ref) records the archive checksum, DOI and download date next to the
input-file checksums.

- `epoch::DateTime` is the UTC start; it must agree with the namelist `day0` and the
  member date (checked).
- `dimensions = :documented` uses [`lasso_documented_dimensions`](@ref) for the member's
  domain token (recorded in provenance as such); pass a tuple to be fully explicit, or
  `nothing` when the namelist carries `nx_gl/ny_gl/nz_gl`.
- A missing bundle raises an `ArgumentError` with [`lasso_bundle_missing_message`](@ref),
  before the architecture is touched (`arch = nothing` becomes `GPU()` only afterwards, so the
  staging message also appears on machines without a GPU).

Everything else (`Nx`, `Ny`, `z_faces`, `microphysics`, `stop_time`, …) is an explicit,
provenance-recorded override of the protocol defaults, exactly as in `ena_simulation`.
"""
function ena_lasso(bundle_dir=nothing;
                   member = DEFAULT_LASSO_MEMBER,
                   epoch,
                   dimensions = :documented,
                   arch = nothing,
                   FT = Float32,
                   output_dir = nothing,
                   output_prefix = "ena_lasso",
                   kwargs...)
    m = member isa LassoMember ? member : parse_lasso_member(member)
    directory = isnothing(bundle_dir) ? lasso_bundle_directory(m) : bundle_dir
    lasso_bundle_available(directory) || throw(ArgumentError(lasso_bundle_missing_message(m, directory)))
    arch = isnothing(arch) ? GPU() : arch
    dims = dimensions === :documented ? lasso_documented_dimensions(m) : dimensions
    output_dir = isnothing(output_dir) ? joinpath(pwd(), "output", "ena_lasso", m.id) : output_dir
    case = ena_simulation(directory; protocol=:lasso_ena_official, member=m, epoch, dimensions=dims,
                          arch, FT, output_dir, output_prefix, kwargs...)
    staged = joinpath(directory, "bundle.toml")
    staging = isfile(staged) ? TOML.parsefile(staged) :
              "no bundle.toml: the directory was not staged with data_wrangling/stage_lasso_bundle.jl (archive checksum/download provenance unavailable)"
    return merge(case, (; staging, dimension_request = string(dimensions)))
end
