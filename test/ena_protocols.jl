@testset "ENA protocols remain distinct" begin
    mktempdir() do dir
        for (source, target) in (("snd_pressure", "snd"), ("lsf_7col", "lsf"), ("sfc", "sfc"))
            cp(joinpath(FIXTURES, source), joinpath(dir, target))
        end
        @test_throws UndefKeywordError ena_simulation(dir)
        @test_throws UndefKeywordError lasso_ena_simulation(dir)
        prm = read(joinpath(FIXTURES, "prm"), String)
        write(joinpath(dir, "prm"), replace(prm, "16x16x24" => "16x16x192"))
        covert = ena_protocol_settings(dir; protocol=:covert_public_bin)
        @test covert.surface == :prescribed_fluxes
        @test covert.radiation == :simple
        @test isnothing(covert.wind_nudging_timescale)
        @test covert.stop_time == 21600
        # Covert files are never accepted as a LASSO bundle: no grd, and Covert switches.
        @test_throws ArgumentError ena_protocol_settings(dir; protocol=:lasso_ena_official)
        @test_throws ArgumentError ena_protocol_settings(dir; protocol=:unknown)
        @test_throws ArgumentError ena_protocol_settings(dir; protocol=:covert_public_bin, member="20170718era5d25x100_sbmwrm-aer2-flxsst")
        write(joinpath(dir, "grd"), join(125:250:5875, '\n'))
        err = try
            ena_protocol_settings(dir; protocol=:lasso_ena_official, dimensions=(16, 16, 24))
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("sfc_flx_fxd", sprint(showerror, err))     # the prescribed-flux Covert switch is named
        @test occursin("donudging_uv", sprint(showerror, err))
    end
    # ... and a LASSO bundle is never accepted as the Covert protocol.
    @test_throws ArgumentError ena_protocol_settings(LASSO_FIXTURE; protocol=:covert_public_bin)
end
