include(joinpath(@__DIR__, "..", "data_wrangling", "fetch_manifest.jl"))

@testset "Pinned input manifest integrity (offline)" begin
    mktempdir() do root
        output = mkdir(joinpath(root, "data"))
        fixture = joinpath(output, "input.txt")
        write(fixture, "pinned contents")
        entry = Dict("filename" => "input.txt", "bytes" => filesize(fixture),
                     "sha256" => bytes2hex(open(sha256, fixture)),
                     "url" => "https://example.invalid/input.txt")
        manifest = Dict("source" => "offline fixture", "files" => [entry])
        path = joinpath(root, "inputs.toml")
        open(io -> TOML.print(io, manifest), path, "w")
        @test fetch_manifest(path, output) === nothing
        write(fixture, "corruption")
        @test_throws ErrorException fetch_manifest(path, output)
        @test read(fixture, String) == "corruption"  # don't silently replace existing data
        entry["filename"] = "../outside.txt"
        open(io -> TOML.print(io, manifest), path, "w")
        @test_throws ErrorException fetch_manifest(path, output)
    end
end
