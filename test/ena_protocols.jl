@testset "ENA protocols remain distinct" begin
    mktempdir() do dir
        for (source, target) in (("snd_pressure", "snd"), ("lsf_7col", "lsf"), ("sfc", "sfc"), ("prm", "prm"))
            cp(joinpath(FIXTURES, source), joinpath(dir, target))
        end
        # the fixture namelist carries the Covert switches ...
        namelist = read_sam_namelist(joinpath(dir, "prm"))
        @test isnothing(BreezeLab.check_covert_namelist(namelist, "prm"))
        case = with_float_type(Float64) do
            ena_covert(; arch=CPU(), data_dir=dir, microphysics=one_moment_microphysics(),
                         z_faces=collect(range(0, 6000, length=25)), write_output=false, progress_interval=100)
        end
        @test case.config.Nx == 16 && case.config.Lx == 16 * 35 && case.config.stop_time == 21600
        @test case.config.translation_velocity_u == 5.0
        # ... and any other switch is refused, with its name
        write(joinpath(dir, "prm"), replace(read(joinpath(FIXTURES, "prm"), String), "donudging_uv = .false." => "donudging_uv = .true."))
        err = try; ena_covert(; arch=CPU(), data_dir=dir); nothing; catch e; e; end
        @test err isa ArgumentError && occursin("donudging_uv", sprint(showerror, err))
        # Covert files are never accepted as a LASSO bundle: no grd, and Covert switches.
        @test_throws ArgumentError ena_lasso(; bundle_dir=dir, epoch=DateTime(2017, 7, 18, 6), arch=CPU())
        write(joinpath(dir, "prm"), read(joinpath(FIXTURES, "prm"), String))
        write(joinpath(dir, "grd"), join(125:250:5875, '\n'))
        err = try; validate_lasso_bundle(dir; dimensions=(16, 16, 24)); nothing; catch e; e; end
        @test err isa ArgumentError
        @test occursin("sfc_flx_fxd", sprint(showerror, err))     # the prescribed-flux Covert switch is named
        @test occursin("donudging_uv", sprint(showerror, err))
    end
    # ... and a LASSO bundle is never accepted as the Covert protocol.
    @test_throws ArgumentError ena_covert(; arch=CPU(), data_dir=LASSO_FIXTURE)
end
