# Download public case references/inputs from a reviewed, checksum-pinned manifest.
# Usage: julia scripts/fetch_manifest.jl cases/seastarr/inputs.toml data/seastarr_22241697
using Downloads, SHA, TOML

function fetch_manifest(manifest_path, output)
    manifest = TOML.parsefile(manifest_path)
    mkpath(output)
    for entry in manifest["files"]
        name = entry["filename"]
        (basename(name) == name && !occursin('\\', name) && name ∉ (".", "..")) ||
            error("Manifest filenames must be plain basenames")
        startswith(entry["url"], "https://") || error("Expected an HTTPS source")
        destination = joinpath(output, name)
        valid(path) = filesize(path) == entry["bytes"] &&
                      bytes2hex(open(sha256, path)) == entry["sha256"]
        if ispath(destination)
            isfile(destination) && valid(destination) || error("Existing file fails checksum: $name")
        else
            mktempdir(output) do temporary
                path = joinpath(temporary, name)
                Downloads.download(entry["url"], path; timeout=180)
                valid(path) || error("Downloaded file fails checksum: $name")
                mv(path, destination)
            end
        end
        println("Verified ", name)
    end
    println("Verified ", length(manifest["files"]), " files; source: ", manifest["source"])
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 2 || error("Usage: julia scripts/fetch_manifest.jl MANIFEST OUTPUT_DIR")
    fetch_manifest(ARGS...)
end
