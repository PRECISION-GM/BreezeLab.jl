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
        @test_throws ArgumentError ena_protocol_settings(dir; protocol=:lasso_ena_official)
        @test_throws ArgumentError ena_protocol_settings(dir; protocol=:unknown)

        # Synthetic input fixture only; not an official ARM bundle or a physics validation.
        lasso_prm = replace(prm, "SFC_FLX_FXD = .true." => "SFC_FLX_FXD = .false.",
                           "donudging_uv = .false." => "donudging_uv = .true.",
                           "doradsimple = .true." => "doradsimple = .false.",
                           "doshortwave = .false." => "doshortwave = .true.",
                           "nstop = 43200" => "nstop = 172800",
                           "nrad = 1" => "read_in_geostrophic_wind = .false., nrad = 120")
        write(joinpath(dir, "prm"), lasso_prm)
        write(joinpath(dir, "grd"), join(125:250:5875, '\n'))
        lasso = ena_protocol_settings(dir; protocol=:lasso_ena_official)
        @test lasso.surface == :bulk_sst && lasso.radiation == :rrtmgp
        @test lasso.wind_nudging_timescale == 10800
        @test lasso.radiation_interval == 60 && lasso.stop_time == 86400
        @test lasso.Nx == lasso.Ny == 16
        @test length(lasso.z_faces) == 25
        @test lasso.aerosol_replenishment == :diagnostic_ccn
        @test_throws ArgumentError ena_protocol_settings(dir; protocol=:covert_public_bin)
        @test_throws ArgumentError ena_simulation(dir; protocol=:lasso_ena_official)
        @test_throws ArgumentError ena_simulation(dir; protocol=:lasso_ena_official, epoch=DateTime(2017, 7, 19))

        # Exercise the official-input route with a tiny, explicitly labelled sensitivity.
        case = ena_simulation(dir; protocol=:lasso_ena_official, epoch=DateTime(2017, 7, 18, 6),
                              Nx=8, Ny=8, microphysics=:one_moment, radiation=nothing, aerosol_replenishment=nothing,
                              stop_time=1.0, write_output=false)
        run!(case.simulation)
        @test all(f -> all(isfinite, interior(f)), values(Oceananigans.prognostic_fields(case.model)))
        path = write_provenance(joinpath(dir, "provenance.toml"), case)
        record = TOML.parsefile(path)
        @test record["protocol"] == "lasso_ena_official"
        @test record["protocol_overrides"]["radiation"] == "nothing"
        @test haskey(record["inputs"], "grd_sha256")
        @test haskey(record["software"], "BreezeLab_source")

        # An official run must never silently choose a 256-square domain.
        write(joinpath(dir, "prm"), replace(lasso_prm, "16x16x24" => "named-case"))
        @test_throws ArgumentError ena_protocol_settings(dir; protocol=:lasso_ena_official)
        explicit = ena_protocol_settings(dir; protocol=:lasso_ena_official, dimensions=(16, 16, 24))
        @test explicit.Nx == 16 && length(explicit.z_faces) == 25
        @test_throws ArgumentError ena_protocol_settings(dir; protocol=:lasso_ena_official, dimensions=(16, 16))
        write(joinpath(dir, "prm"), lasso_prm)
        @test_throws ArgumentError ena_protocol_settings(dir; protocol=:lasso_ena_official, dimensions=(32, 32, 24))
        write(joinpath(dir, "prm"), replace(lasso_prm, "doshortwave = .true." => "doshortwave = .false."))
        @test_throws ArgumentError ena_protocol_settings(dir; protocol=:lasso_ena_official)
        write(joinpath(dir, "prm"), lasso_prm)
        write(joinpath(dir, "grd"), "125\n125\n250\n")
        @test_throws ArgumentError ena_protocol_settings(dir; protocol=:lasso_ena_official)
    end
end
