#####
##### Provenance records.
#####

"""
    write_provenance(path, case; extra=NamedTuple())

Write a TOML provenance record of a case built by any BreezeLab constructor: the input file
paths and SHA-256 checksums (`case.inputs`), Breeze, Oceananigans and BreezeLab sources
(pinned revision and checkout state), the Julia version, the Manifest checksum, and the
configuration record `case.config` (values normalized to TOML scalars), whose `protocol`
entry is also written at the top level. Round-trips through `TOML.parsefile`.
"""
function write_provenance(path, case; extra=NamedTuple())
    inputs = Dict{String, Any}()
    for (name, file) in pairs(case.inputs)
        isnothing(file) && continue
        inputs[string(name)] = abspath(file)
        inputs[string(name, "_sha256")] = file_sha256(file)
    end
    manifest = package_path("Manifest.toml")
    software = Dict{String, Any}(
        "Breeze" => string(Base.pkgversion(Breeze)),
        "Breeze_source" => breeze_source_description(),
        "Oceananigans" => string(Base.pkgversion(Oceananigans)),
        "Oceananigans_source" => oceananigans_source_description(),
        "BreezeLab_source" => package_source_description(BreezeLab),
        "julia" => string(VERSION))
    isfile(manifest) && (software["BreezeLab_manifest_sha256"] = file_sha256(manifest))
    record = Dict{String, Any}(
        "generated" => string(Dates.now()),
        "protocol" => string(get(case.config, :protocol, "custom")),
        "inputs" => inputs,
        "software" => software,
        "config" => Dict{String, Any}(string(k) => toml_value(v) for (k, v) in pairs(case.config)),
        "extra" => Dict{String, Any}(string(k) => toml_value(v) for (k, v) in pairs(extra)))
    open(path, "w") do io
        TOML.print(io, record)
    end
    return path
end

toml_value(x::Union{Bool, Integer, AbstractFloat, AbstractString}) = x
toml_value(x::Symbol) = string(x)
toml_value(::Nothing) = "nothing"
toml_value(x::Tuple) = [toml_value(v) for v in x]
toml_value(x::AbstractVector) = [toml_value(v) for v in x]
toml_value(x::NamedTuple) = Dict{String, Any}(string(k) => toml_value(v) for (k, v) in pairs(x))
toml_value(x::AbstractDict) = Dict{String, Any}(string(k) => toml_value(v) for (k, v) in x)
toml_value(x) = string(x)

# A pinned dependency's git revision (from the active Manifest) and the dirty state of a
# `dev`ed checkout, if that is how the package is being loaded.
breeze_source_description() = package_source_description(Breeze)
oceananigans_source_description() = package_source_description(Oceananigans)

function package_source_description(package::Module)
    dir = Base.pkgdir(package)
    name = string(nameof(package))
    # BreezeLab's own Manifest carries the `[sources]` pins; the active project's Manifest is
    # the fallback (inside `Pkg.test` the active project is a sandbox).
    manifests = (package_path("Manifest.toml"),
                 Base.active_project() === nothing ? "" : joinpath(dirname(Base.active_project()), "Manifest.toml"))
    rev = "unknown"
    block_pattern = Regex("(?s)\\[\\[deps\\.$name\\]\\]\\n(.*?)(?=\\n\\[\\[|\\z)")
    for manifest in manifests
        isfile(manifest) || continue
        block = match(block_pattern, read(manifest, String))
        isnothing(block) && continue
        repo_rev = match(r"repo-rev = \"([^\"]+)\"", block.captures[1])
        repo_url = match(r"repo-url = \"([^\"]+)\"", block.captures[1])
        (isnothing(repo_rev) || isnothing(repo_url)) && continue
        rev = string(repo_url.captures[1], "@", repo_rev.captures[1])
        break
    end
    dirty = ""
    if isdir(joinpath(dir, ".git"))
        status = try
            read(`git -C $dir status --porcelain`, String)
        catch
            ""
        end
        head = try
            strip(read(`git -C $dir rev-parse HEAD`, String))
        catch
            "unknown"
        end
        dirty = string(" (checkout ", head, isempty(strip(status)) ? ", clean)" : ", DIRTY)")
    end
    return string(rev, " at ", dir, dirty)
end
