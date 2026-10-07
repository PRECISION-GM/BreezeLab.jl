# TRACER–DP-SCREAM: IOP forcing reader, SAM-record adapter, reference-output reader and the
# periodic case constructor, against small redistributable excerpts (data_wrangling/make_tracer_fixtures.jl).
using Dates: Dates, DateTime, Hour
using JLD2
using Oceananigans.Units: Time

const IOP_FIXTURE = joinpath(FIXTURES, "tracer_iop_excerpt.nc")
const DPS_FIXTURE = joinpath(FIXTURES, "dp_scream_3km_excerpt.nc")
const DPS_NS_FIXTURE = joinpath(FIXTURES, "dp_scream_ns_excerpt.nc")

@testset "TRACER–DP-SCREAM" begin
    @testset "IOP forcing reader (reference values, 2022-08-05 00 UTC)" begin
        iop = read_iop_forcing(IOP_FIXTURE)
        @test iop.base_time == DateTime(2022, 7, 1)
        @test length(iop.times) == 4 && iop.times[1] == 35 * 86400 && iop.times[2] - iop.times[1] == 3600
        @test iop_datetimes(iop)[1] == DateTime(2022, 8, 5) && iop_datetimes(iop)[4] == DateTime(2022, 8, 5, 3)
        @test iop_index(iop, DateTime(2022, 8, 5, 2)) == 3
        @test_throws ArgumentError iop_index(iop, DateTime(2022, 8, 5, 2, 30))
        @test iop.latitude == 29.75 && iop.longitude == -95.45
        @test length(iop.levels) == 40 && iop.levels[1] == 5000 && iop.levels[end] == 102500   # Pa, top first
        @test iop.attributes["datastream"] == "hou60v1varanaecmwfX12.c1"
        k1000 = findfirst(==(100000), iop.levels); k500 = findfirst(==(50000), iop.levels)
        @test isapprox(iop.profiles.T[k1000, 1], 304.93; atol=0.01)
        @test isapprox(iop.profiles.q[k1000, 1], 15.822e-3; rtol=1e-3)                    # specific humidity, kg/kg
        @test isapprox(iop.profiles.u[k1000, 1], -1.19; atol=0.01) && isapprox(iop.profiles.v[k1000, 1], 5.19; atol=0.01)
        @test isapprox(iop.profiles.T[k500, 1], 268.12; atol=0.01)
        @test isapprox(iop.profiles.divT[k500, 1], 0.857e-5; rtol=2e-3)
        @test isapprox(iop.profiles.vertdivT[k500, 1], 11.485e-5; rtol=2e-3)
        @test isapprox(iop.profiles.divq[k500, 1], -1.316e-8; rtol=2e-3)
        @test isapprox(iop.profiles.vertdivq[k500, 1], -3.328e-8; rtol=2e-3)
        @test isapprox(iop.surface.Ps[1], 100867.34; atol=0.01)
        @test isapprox(iop.surface.Tg[1], 305.01; atol=0.01) && isapprox(iop.surface.Tsair[1], 305.99; atol=0.01)
        @test isapprox(iop.surface.shflx[1], 13.429; atol=1e-3) && isapprox(iop.surface.lhflx[1], 31.688; atol=1e-3)
        @test all(isfinite, iop.profiles.divT) && all(isfinite, iop.surface.Tg)
        # file-level sanity used by the adapter: below-surface levels repeat the lowest level above the surface
        @test iop.profiles.T[end, 1] == iop.profiles.T[end - 1, 1]
        # the "net" shortwave labels hold upwelling fluxes: albedo ≈ 0.15, not 0.85
        @test iop_surface_albedo(iop; start=DateTime(2022, 8, 5), stop=DateTime(2022, 8, 5, 3)) ≈ 0 atol=1e-12 skip=true
        @test fractional_day_of_year(DateTime(2022, 8, 5)) == 217.0
        @test fractional_day_of_year(DateTime(2022, 8, 1, 6)) == 213.25
        @test iop_potential_temperature(300.0, 1e5) == 300.0
        κ = 287.0 / 1005.0
        @test isapprox(iop_potential_temperature(268.12, 50000.0), 268.12 * 2^κ; rtol=1e-3)
    end

    @testset "SAM-record adapter: 3-D transport once, hold expansion, mass-fraction basis" begin
        iop = read_iop_forcing(IOP_FIXTURE)
        inputs = iop_sam_inputs(iop; start=DateTime(2022, 8, 5), stop=DateTime(2022, 8, 5, 3))
        @test inputs.day0 == 217.0 && inputs.epoch == DateTime(2022, 8, 5)
        @test inputs.window == (1, 4) && inputs.paths.iop == abspath(IOP_FIXTURE)
        snd = inputs.soundings[1]
        @test length(inputs.soundings) == 1
        @test all(isnan, snd.z) && issorted(snd.p; rev=true) && snd.p[1] == 102500 && snd.surface_pressure ≈ 100867.34 atol=0.01
        k = findfirst(==(50000), snd.p)
        @test isapprox(snd.θ[k], iop_potential_temperature(268.12, 50000.0); rtol=1e-3)
        @test snd.q[k] == Float64(iop.profiles.q[findfirst(==(50000), iop.levels), 1])      # no unit conversion: specific humidity
        # tls = divT + vertdivT, qls = divq + vertdivq at the same level; wls = 0; nudging targets = u, v
        r1 = inputs.lsf[1]
        @test isapprox(r1.tls[k], (0.857 + 11.485) * 1e-5; rtol=2e-3)
        @test isapprox(r1.qls[k], (-1.316 - 3.328) * 1e-8; rtol=2e-3)
        @test all(iszero, r1.wls) && r1.uls == r1.ug && r1.vls == r1.vg && !r1.has_geostrophic_columns
        @test isapprox(r1.uls[findfirst(==(100000), snd.p)], -1.19; atol=0.01)
        # hold expansion: 4 records → 7, strictly increasing, each copy 1 s before the next record
        @test length(inputs.lsf) == 7
        days = [r.day for r in inputs.lsf]
        @test issorted(days) && allunique(days)
        @test inputs.lsf[2].tls == inputs.lsf[1].tls && inputs.lsf[2].day ≈ inputs.lsf[3].day - 1 / 86400
        @test inputs.lsf[3].day == 217.0 + 1 / 24
        @test length(inputs.sfc.day) == 7 && inputs.sfc.sst[1] ≈ 305.01 atol=0.01 && inputs.sfc.sensible_heat_flux[2] == inputs.sfc.sensible_heat_flux[1]
        @test all(iszero, inputs.sfc.kinematic_stress)
        linear = iop_sam_inputs(iop; start=DateTime(2022, 8, 5), stop=DateTime(2022, 8, 5, 3), hold=false)
        @test length(linear.lsf) == 4 && length(linear.sfc.day) == 4
        horizontal = iop_sam_inputs(iop; start=DateTime(2022, 8, 5), stop=DateTime(2022, 8, 5, 1), transport=:iop_horizontal)
        @test isapprox(horizontal.lsf[1].tls[k], 0.857e-5; rtol=2e-3)
        @test_throws ArgumentError iop_sam_inputs(iop; start=DateTime(2022, 8, 5, 1), stop=DateTime(2022, 8, 5))
        @test_throws ArgumentError iop_sam_inputs(iop; start=DateTime(2022, 8, 5), stop=DateTime(2022, 8, 5, 1), transport=:omega)
        # the forcing profiles interpolated onto a model grid: linear in pressure, zero above the 50-hPa top
        grid = test_grid(; Nz=24, Lz=24000)
        zc = Array(znodes(grid, Center()))
        pᵣ = 100867.34 .* exp.(-zc ./ 7500)
        profiles = LargeScaleForcingProfiles(grid, inputs.lsf, zc, pᵣ; day0=inputs.day0)
        @test profiles.times[1] == 0 && profiles.times[3] == 3600 && profiles.times[2] == 3599
        @test profiles.tls[1, 1, grid.Nz, Time(0.0)] == 0 && profiles.qls[1, 1, grid.Nz, Time(0.0)] == 0  # above 50 hPa
        @test profiles.tls[1, 1, 1, Time(1800.0)] == profiles.tls[1, 1, 1, Time(0.0)]                     # held within the hour
        @test profiles.tls[1, 1, 1, Time(3600.0)] != profiles.tls[1, 1, 1, Time(0.0)]
        # column integral of the vapor source over the IOP levels (kg m⁻² s⁻¹)
        S = iop_column_integral(iop, :divq, 1) + iop_column_integral(iop, :vertdivq, 1)
        @test isfinite(S) && abs(S) < 1e-3
    end

    @testset "DP-SCREAM reference-output reader (local-time labels → UTC)" begin
        out = read_dp_scream_output(DPS_FIXTURE)
        @test out.git_version == "1d551ea2b0" && out.case == "scream_dp_TRACER_AUGUST_decr"
        @test out.local_time[1] == DateTime(2022, 8, 4, 19) && out.time[1] == DateTime(2022, 8, 5)
        @test out.time[2] - out.time[1] == Dates.Minute(30) && length(out.time) == 8
        @test length(out.lev) == 128 && out.lev[1] < out.lev[end]                      # hPa, top first
        @test size(out.T) == (128, 8) && size(out.TOT_CLOUD_FRAC) == (128, 8) && length(out.PRECL) == 8
        @test isapprox(out.PRECL[2], 0.2099; rtol=1e-3)                                # mm/day (inferred)
        @test isapprox(out.T[end, 2], 312.25; atol=0.01)                               # lowest level, record 2
        @test isapprox(out.TGCLDLWP[2], 0.0; atol=1e-4) && out.TGCLDIWP[2] > 0
        @test all(0 .≤ filter(isfinite, out.TOT_CLOUD_FRAC) .≤ 1)
        @test dp_scream_window(out, DateTime(2022, 8, 5), DateTime(2022, 8, 5, 1)) == 1:3
        ns = read_dp_scream_output(DPS_NS_FIXTURE)
        @test ns.local_time[1] == DateTime(2022, 7, 31, 19) && ns.time[1] == DateTime(2022, 8, 1)
        @test parse_cf_time_units("minutes since 2022-08-04 19:00:00") == (6e4, DateTime(2022, 8, 4, 19))
        @test_throws ArgumentError parse_cf_time_units("fortnights since 2022-08-04")
    end

    @testset "Periodic case constructor builds, steps, writes and checkpoints (CPU)" begin
        z_faces = tracer_dp_scream_vertical_faces()
        @test length(z_faces) == 161 && z_faces[1] == 0 && z_faces[end] == 22000 && issorted(z_faces)
        @test all(diff(z_faces[1:41]) .≈ 50) && diff(z_faces)[end] < 500
        @test isapprox(neutral_drag_coefficient(25.0, 0.1), (0.4 / log(250))^2)
        mktempdir() do output
            case = tracer_dp_scream(; iop_path=IOP_FIXTURE, start=DateTime(2022, 8, 5), stop=DateTime(2022, 8, 5, 3),
                                    arch=CPU(), FT=Float32, Nx=4, Ny=4, Δx=500.0,
                                    z_faces=collect(range(0, 22000, length=25)), microphysics=:p3_n75,
                                    stop_time=4.0, Δt=1.0, max_Δt=1.0, output_dir=output,
                                    timeseries_interval=1.0, profile_interval=4.0, slice_interval=4.0,
                                    checkpoint_interval=2.0, progress_interval=100)
            @test case.preset === :tracer_dp_scream
            @test case.model.clock.time == 0 && case.model.clock.iteration == 0
            @test case.config.day0 == 217.0 && case.config.epoch == string(DateTime(2022, 8, 5))
            @test case.config.moisture_basis == "mass_fraction"
            @test case.config.coriolis == "nothing" && !case.config.geostrophic
            @test case.config.vertical_advection == "nothing" && !case.config.upper_boundary_relaxation
            @test case.config.wind_nudging_timescale == 0
            @test case.config.surface == "prescribed_heat_fluxes_bulk_drag"
            @test isapprox(case.config.drag_coefficient, neutral_drag_coefficient(22000 / 24 / 2, 0.1))
            @test case.config.radiation == "rrtmgp" && case.config.radiation_interval == 300.0
            @test 0.1 < case.config.surface_albedo < 0.2
            @test isapprox(case.config.latitude, 29.75) && isapprox(case.config.longitude, -95.45)
            @test case.config.checkpoint_interval == 2.0
            @test occursin("no Coriolis", case.config.label) && occursin("not a DP-SCREAM reproduction", case.config.label)
            @test case.protocol_overrides.dp_scream_git_version == "1d551ea2b0"
            # large-scale transport enters exactly once: thermodynamic tendencies on the energy and
            # vapor keys, no LargeScaleVerticalAdvection, no geostrophic/Coriolis/nudging momentum forcing
            forcing = case.model.forcing
            @test !haskey(forcing, :ρu) && !haskey(forcing, :ρv)
            @test haskey(forcing, :ρθ) && haskey(forcing, :ρqᵛ)
            @test isnothing(case.model.coriolis)
            @test !any(f -> inner(f) isa LargeScaleVerticalAdvection, values(forcing))
            # surface: prescribed energy/vapor fluxes and bulk drag
            bcs = case.model.momentum.ρu.boundary_conditions.bottom
            @test bcs.condition isa Breeze.BoundaryConditions.BulkDragFunction
            @test case.surface_temperature[1, 1, 1] ≈ 305.01 atol=0.01
            # initial state from the IOP sounding: lowest level near the 1000-hPa record temperature
            T = Array(interior(case.model.temperature))
            @test all(isfinite, T) && 300 < T[1, 1, 1] < 306 && T[1, 1, end] < 230
            write_provenance(joinpath(output, "provenance.toml"), case)
            record = TOML.parsefile(joinpath(output, "provenance.toml"))
            @test record["protocol"] == "tracer_dp_scream" && haskey(record["inputs"], "iop_sha256")
            @test record["protocol_overrides"]["large_scale_transport"] == "iop_3d"
            run!(case.simulation)
            @test case.model.clock.time ≈ 4.0 && case.model.clock.iteration == 4
            @test all(f -> all(isfinite, Array(interior(f))), values(Oceananigans.prognostic_fields(case.model)))
            ts = joinpath(output, "tracer_timeseries.jld2")
            @test isfile(ts)
            jldopen(ts) do file
                series = file["timeseries"]
                for name in ("lwp", "iwp", "precipitable_water", "column_static_energy", "column_radiative_heating", "rain_flux", "ice_flux")
                    @test haskey(series, name)
                end
                iterations = sort(parse.(Int, collect(keys(series["t"]))))
                @test series["t/$(last(iterations))"] ≈ 4.0
                @test series["precipitable_water/$(last(iterations))"][1, 1, 1] > 10   # kg m⁻², moist Houston column
            end
            jldopen(joinpath(output, "tracer_profiles.jld2")) do file
                @test haskey(file["timeseries"], "total_cloud_fraction")
            end
            checkpoints = filter(f -> startswith(f, "tracer_checkpoint") && endswith(f, ".jld2"), readdir(output))
            @test length(checkpoints) == 1                                                  # cleanup keeps the latest
            # pickup: a fresh case continues from the checkpoint to a later stop time
            resumed = tracer_dp_scream(; iop_path=IOP_FIXTURE, start=DateTime(2022, 8, 5), stop=DateTime(2022, 8, 5, 3),
                                       arch=CPU(), FT=Float32, Nx=4, Ny=4, Δx=500.0,
                                       z_faces=collect(range(0, 22000, length=25)), microphysics=:p3_n75,
                                       stop_time=6.0, Δt=1.0, max_Δt=1.0, output_dir=output,
                                       timeseries_interval=1.0, profile_interval=4.0, slice_interval=4.0,
                                       checkpoint_interval=2.0, progress_interval=100)
            run!(resumed.simulation; pickup=true)
            @test resumed.model.clock.time ≈ 6.0 && resumed.model.clock.iteration == 6
            @test all(f -> all(isfinite, Array(interior(f))), values(Oceananigans.prognostic_fields(resumed.model)))
        end
        @test_throws ArgumentError tracer_dp_scream(; iop_path=IOP_FIXTURE, start=DateTime(2022, 8, 5), stop=DateTime(2022, 8, 5, 3),
                                                     arch=CPU(), Nx=4, Ny=4, z_faces=collect(range(0, 22000, length=25)),
                                                     stop_time=1e6, write_output=false)
        @test_throws ArgumentError tracer_dp_scream(; iop_path="/nonexistent/iop.nc", arch=CPU())
        settings = tracer_dp_scream_settings(IOP_FIXTURE; start=DateTime(2022, 8, 5), stop=DateTime(2022, 8, 5, 3))
        @test settings.records == 4 && settings.duration_seconds == 10800 && settings.doi == "10.5439/1860369"
    end
end
