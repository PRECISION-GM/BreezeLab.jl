# Exercise the readable case entry point, including its saved diagnostic outputs.
using JLD2
include(joinpath(@__DIR__, "..", "cases", "eastern_north_atlantic.jl"))

function test_ena_execution(arch, FT, microphysics)
    mktempdir() do output
        case = eastern_north_atlantic(; arch, FT, microphysics, output_dir=output,
                                      Nx=8, Ny=8, z_faces=collect(range(0, 6000, length=25)),
                                      stop_time=4.0, timeseries_interval=1.0,
                                      profile_interval=4.0, slice_interval=4.0)
        write_provenance(joinpath(output, "provenance.toml"), case)
        run!(case.simulation)
        @test case.model.clock.time ≈ 4.0
        @test case.model.clock.iteration == 8
        @test all(f -> all(isfinite, Array(interior(f))), values(Oceananigans.prognostic_fields(case.model)))

        expected = Dict("profiles" => ("u", "v", "T", "qᵛ", "qᶜˡ", "qʳ"),
                        "timeseries" => ("lwp", "rwp", "cloud_fraction", "rain_flux"),
                        "slices" => ("qᶜˡ_xz", "qʳ_xz", "w_xz", "lwp", "rain"))
        for (kind, fields) in expected
            path = joinpath(output, "lasso_ena_$(kind).jld2")
            @test isfile(path)
            jldopen(path) do file
                series = file["timeseries"]
                @test all(name -> haskey(series, name), fields)
                iterations = sort(parse.(Int, collect(keys(series["t"]))))
                @test series["t/$(last(iterations))"] ≈ 4.0
                for name in fields
                    # Oceananigans stores field metadata beside the numbered snapshots.
                    snapshots = filter(!=("serialized"), collect(keys(series[name])))
                    @test sort(parse.(Int, snapshots)) == iterations
                    @test all(i -> all(isfinite, series["$name/$i"]), iterations)
                end
            end
        end
    end
end

@testset "ENA case execution and output" begin
    if HAVE_COVERT
        @testset "$scheme on CPU" for scheme in (:one_moment, :p3_n75, :p3_aer2)
            test_ena_execution(CPU(), scheme === :one_moment ? Float64 : Float32, scheme)
        end
        if "gpu" in ARGS
            @test CUDA.functional()
            @testset "$scheme on GPU" for scheme in (:one_moment, :p3_n75, :p3_aer2)
                test_ena_execution(GPU(), Float32, scheme)
            end
        end
    else
        @test_skip HAVE_COVERT
        @info "Fetch data_wrangling/fetch_covert_inputs.jl to enable ENA execution checks"
    end
end
