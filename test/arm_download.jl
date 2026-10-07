# Run independently of the simulation environment: julia test/arm_download.jl
using Test, Tar, TOML
include(joinpath(@__DIR__, "..", "scripts", "fetch_arm_inputs.jl"))

@testset "ARM archive download safety and provenance" begin
    name = "enalasso_samin_20170718era5s1n0d25x100_sbmwrm-aer2-flxsstC1.m0.20170718.000000.tar"
    credentials = (username="test-user", token="test-secret")
    mktempdir() do root
        source = mkdir(joinpath(root, "source"))
        write(joinpath(source, "snd"), "synthetic fixture, not scientific inputs\n")
        fixture = Tar.create(source, joinpath(root, "fixture.tar"))
        output = joinpath(root, "downloads")
        success = (url, file) -> (cp(fixture, file); 200)
        dest = ARMInputs.fetch_archive(name, output; credentials..., request=success)
        metadata = TOML.parsefile(dest * ".toml")
        @test read(dest) == read(fixture)
        @test metadata["bytes"] == filesize(fixture)
        @test metadata["sha256"] == bytes2hex(open(ARMInputs.sha256, fixture))
        @test !occursin(credentials.token, read(dest * ".toml", String))
        @test_throws ErrorException ARMInputs.fetch_archive(name, output; credentials..., request=success)
        @test_throws ErrorException ARMInputs.fetch_archive("../" * name, output; credentials..., request=success)
        for failure in ((url, file) -> 404,
                        (url, file) -> (write(file, "<html>Unavailable</html>"); 200),
                        (url, file) -> error(url))
            bad_output = mktempdir(root)
            err = try
                ARMInputs.fetch_archive(name, bad_output; credentials..., request=failure)
                nothing
            catch e
                e
            end
            @test err isa ErrorException
            @test !occursin(credentials.token, sprint(showerror, err))
            @test isempty(readdir(bad_output))
        end
    end
    @test ARMInputs.urlencode("a:b+c/@") == "a%3Ab%2Bc%2F%40"
end
