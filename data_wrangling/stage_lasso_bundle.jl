#!/usr/bin/env julia
# Freeze one LASSO-ENA `samin` archive into a case input directory (plan step 1):
#
#   julia data_wrangling/stage_lasso_bundle.jl data/lasso/archives/enalasso_samin_<member>C1.m0.<date>.000000.tar [data/lasso]
#
# Standard-library only (SHA, Tar, TOML, Dates): no model initialization. The script
#   * accepts only an exact ARM samin archive name and derives the member ID from it,
#   * records the archive's SHA-256 and size, the `.tar.toml` download sidecar written by
#     fetch_arm_inputs.jl (retrieval time, service, DOI) when present,
#   * extracts into a temporary directory, locates `snd`, `lsf`, `sfc`, `prm` and `grd`
#     wherever the archive places them (exactly one of each), copies them to
#     <output>/<member>/ and every other archive member to <output>/<member>/extra/,
#   * writes <output>/<member>/bundle.toml with per-file checksums, the member tokens,
#     the DOI and the LASSO SAM reference revision, and refuses to overwrite a staged member.
# The namelist is not interpreted here; `BreezeLab.inspect_lasso_bundle` does that and
# reports the UTC epoch, grid, duration, and any unsupported setting.
module LassoStaging

using Dates, SHA, Tar, TOML

const DOI = "10.5439/2572661"
const SAM_REFERENCE = Dict("url" => "https://code.arm.gov/lasso/lasso-ena-codes/lasso_sam_sbm.git",
                           "branch" => "lasso_ena_noice", "commit" => "12d02446a2147388dc89d828e6e0553106abea0f")
const ARCHIVE_PATTERN = r"^enalasso_samin_((\d{8})(era5|merra2)((?:[a-z]\d+)*)(d\d+x\d+)_(sbmwrm|morr)-(aer[123])-(flx[a-z]+))C1\.m0\.(\d{8})\.(\d{6})\.tar$"
const REQUIRED = ("snd", "lsf", "sfc", "prm", "grd")

sha256_hex(path) = bytes2hex(open(sha256, path))

"""Member tokens of an exact samin archive name (error for any other name)."""
function member_from_archive(filename)
    m = match(ARCHIVE_PATTERN, filename)
    isnothing(m) && error("Not an exact LASSO-ENA samin archive name: $filename")
    m.captures[2] == m.captures[9] || error("Archive date $(m.captures[9]) differs from the member date $(m.captures[2])")
    return Dict("id" => m.captures[1], "date" => m.captures[2], "forcing" => m.captures[3], "variant" => m.captures[4],
                "domain" => m.captures[5], "sam_microphysics" => m.captures[6], "aerosol" => m.captures[7],
                "surface" => m.captures[8], "archive" => filename, "doi" => DOI)
end

function stage(archive, output_root)
    isfile(archive) || error("Archive not found: $archive")
    filename = basename(archive)
    member = member_from_archive(filename)
    destination = joinpath(abspath(output_root), member["id"])
    ispath(destination) && error("$(destination) exists; staged members are never overwritten (remove it deliberately first)")
    sidecar = archive * ".toml"
    download = isfile(sidecar) ? TOML.parsefile(sidecar) : Dict{String, Any}("note" => "no fetch_arm_inputs.jl sidecar found")
    digest = sha256_hex(archive)
    if haskey(download, "sha256") && download["sha256"] != digest
        error("Archive checksum $digest differs from the download record $(download["sha256"])")
    end
    headers = Tar.list(archive; strict=true)
    any(h -> h.type == :file, headers) || error("Archive contains no regular files")

    mktempdir(abspath(output_root)) do temporary
        extracted = joinpath(temporary, "extracted")
        Tar.extract(archive, extracted)
        found = Dict{String, Vector{String}}(name => String[] for name in REQUIRED)
        others = String[]
        for (root, _, files) in walkdir(extracted), file in files
            path = joinpath(root, file)
            rel = relpath(path, extracted)
            if file in REQUIRED
                push!(found[file], rel)
            else
                push!(others, rel)
            end
        end
        for name in REQUIRED
            length(found[name]) == 1 || error("Archive holds $(length(found[name])) files named $name (need exactly one): $(found[name])")
        end
        staged = joinpath(temporary, "staged")
        mkpath(joinpath(staged, "extra"))
        checksums = Dict{String, Any}()
        for name in REQUIRED
            cp(joinpath(extracted, found[name][1]), joinpath(staged, name))
            checksums[name] = sha256_hex(joinpath(staged, name))
        end
        for rel in others
            target = joinpath(staged, "extra", rel)
            mkpath(dirname(target))
            cp(joinpath(extracted, rel), target)
        end
        record = Dict{String, Any}(
            "member" => member,
            "archive" => Dict("filename" => filename, "sha256" => digest, "bytes" => filesize(archive),
                              "source_path" => abspath(archive), "download" => download),
            "staged_utc" => string(now(UTC)),
            "files" => Dict(name => found[name][1] for name in REQUIRED),
            "checksums" => checksums,
            "extra_files" => others,
            "sam_reference" => SAM_REFERENCE,
            "note" => "Epoch, grid, duration and unsupported settings are reported by BreezeLab.inspect_lasso_bundle; never substitute Covert inputs for this member.")
        open(io -> TOML.print(io, record), joinpath(staged, "bundle.toml"), "w")
        mkpath(dirname(destination))
        mv(staged, destination)
    end
    return destination
end

function main(args)
    if args == ["--help"] || isempty(args)
        println("Usage: julia data_wrangling/stage_lasso_bundle.jl ARCHIVE.tar [OUTPUT_ROOT=data/lasso]")
        return
    end
    length(args) ≤ 2 || error("Expected ARCHIVE.tar [OUTPUT_ROOT]; use --help")
    output_root = length(args) == 2 ? args[2] : joinpath(@__DIR__, "..", "data", "lasso")
    println("Staged ", stage(args[1], output_root))
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    try
        LassoStaging.main(ARGS)
    catch err
        println(stderr, sprint(showerror, err))
        exit(1)
    end
end
