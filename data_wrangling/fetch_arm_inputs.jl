#!/usr/bin/env julia
# Standalone standard-library tool: no Breeze/CUDA initialization required.
module ARMInputs

using Dates, Downloads, SHA, TOML, Tar

const SERVICE = "https://adc.arm.gov/armlive/saveData"
# Size cap for one file; `ARM_MAX_BYTES` raises it for the larger per-member reference
# outputs (samstat/sam2d NetCDF). sam3d/sam3dmicro/samrest files are never accepted.
const MAX_BYTES = parse(Int, get(ENV, "ARM_MAX_BYTES", string(200 * 1024^2)))
# Accepted LASSO-ENA filenames, exactly as the ARM Live listing spells them: the samin
# input bundle (tar, level m0) and the two per-member reference outputs (samstat domain
# statistics and sam2d two-dimensional fields, NetCDF, level m1).
const ACCEPTED_FILENAME = r"^enalasso_(samin_[A-Za-z0-9_-]+C1\.m0\.\d{8}\.\d{6}\.tar|(samstat|sam2d)_[A-Za-z0-9_-]+C1\.m1\.\d{8}\.\d{6}\.nc)$"

urlencode(s) = join(isletter(Char(b)) && b < 128 || isdigit(Char(b)) || b in codeunits("-._~") ?
                    string(Char(b)) : "%" * uppercase(string(b; base=16, pad=2)) for b in codeunits(s))

function request_archive(url, output)
    downloader = Downloads.Downloader()
    # The URL contains the API token. Never forward it through a redirect.
    downloader.easy_hook = (easy, info) ->
        Downloads.Curl.setopt(easy, Downloads.Curl.CURLOPT_FOLLOWLOCATION, false)
    progress = (total, now, upload_total, uploaded) ->
        (max(total, now) > MAX_BYTES && error("Archive exceeds download size limit"))
    response = Downloads.request(url; output, downloader, progress, timeout=180,
                                 verbose=false, throw=false)
    response isa Downloads.Response || error("ARM connection failed")
    return response.status
end

"""Download one exact ARM LASSO-ENA file (samin archive or samstat/sam2d reference output); never extract or overwrite it."""
function fetch_archive(filename, output_dir;
                       username=get(ENV, "ARM_USERNAME", ""), token=get(ENV, "ARM_TOKEN", ""),
                       request=request_archive)
    occursin(ACCEPTED_FILENAME, filename) ||
        error("Expected an exact LASSO-ENA samin .tar (or samstat/sam2d .nc) filename from the ARM listing")
    (isempty(username) || isempty(token)) && error("Set ARM_USERNAME and ARM_TOKEN outside the repository")
    mkpath(output_dir)
    destination = abspath(joinpath(output_dir, filename))
    metadata = destination * ".toml"
    (ispath(destination) || ispath(metadata)) && error("Destination already exists; refusing to overwrite")
    url = SERVICE * "?user=" * urlencode(username * ":" * token) * "&file=" * urlencode(filename)

    mktempdir(output_dir) do temporary
        archive = joinpath(temporary, "input.tar")
        status = try
            request(url, archive)
        catch
            # Downloads errors can contain the authenticated URL. Do not print them.
            error("ARM transfer failed; request details suppressed to protect credentials")
        end
        status == 404 && error("ARM returned HTTP 404 for this archive; no input was installed")
        status == 200 || error("ARM returned HTTP $status; no input was installed")
        0 < filesize(archive) <= MAX_BYTES || error("Empty or oversized ARM response")
        headers = try
            Tar.list(archive; strict=true)
        catch
            error("ARM response is not a valid tar archive; no input was installed")
        end
        any(h -> h.type == :file, headers) || error("ARM archive contains no regular files")
        record = Dict("filename" => filename, "service" => SERVICE,
                      "downloaded_utc" => string(now(UTC)), "doi" => "10.5439/2572661",
                      "bytes" => filesize(archive),
                      "sha256" => bytes2hex(open(sha256, archive)))
        staged_metadata = joinpath(temporary, "provenance.toml")
        open(io -> TOML.print(io, record), staged_metadata, "w")
        mv(archive, destination)
        mv(staged_metadata, metadata)
    end
    return destination
end

function main(args)
    if args == ["--help"]
        println("Usage: julia data_wrangling/fetch_arm_inputs.jl FILENAME OUTPUT_DIR")
        println("Load ARM_USERNAME and ARM_TOKEN through your local secret configuration, not command history or Git.")
        println("Downloads one samin tar plus checksum/provenance; does not extract it.")
        return
    end
    length(args) == 2 || error("Expected FILENAME OUTPUT_DIR; use --help")
    println("Saved ", fetch_archive(args[1], args[2]))
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    try
        ARMInputs.main(ARGS)
    catch err
        println(stderr, sprint(showerror, err))
        exit(1)
    end
end
