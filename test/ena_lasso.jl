# The official LASSO-ENA adapter exercised on a synthetic SAM-style bundle fixture
# (test/fixtures/lasso_bundle: analytic linear forcing profiles, two soundings, three SST
# samples, a 30-level grd and a LASSO-style namelist). None of this is ARM data; the tests
# establish the parser, validation, unit/tendency conversions and the forcing assembly,
# not physical fidelity.
using Dates: DateTime, Date, Second
using Oceananigans.Units: Time
using Breeze.Thermodynamics: MoistureMassFractions

const LASSO_MEMBER = "20170718era5d25x100_sbmwrm-aer2-flxsst"
const LASSO_DIMS = (16, 16, 30)
const LASSO_EPOCH = DateTime(2017, 7, 18, 0)

lasso_copy(dir) = (foreach(f -> cp(joinpath(LASSO_FIXTURE, f), joinpath(dir, f)), BreezeLab.LASSO_BUNDLE_FILES); dir)
problem_text(dir; kwargs...) = join(inspect_lasso_bundle(dir; kwargs...).problems, "\n")
rewrite_prm(dir, pairs...) = write(joinpath(dir, "prm"),
                                   foldl((s, p) -> replace(s, p), pairs; init=read(joinpath(LASSO_FIXTURE, "prm"), String)))
caught(f) = try; f(); nothing; catch e; e; end

@testset "LASSO-ENA member identity" begin
    m = parse_lasso_member(LASSO_MEMBER)
    @test m.date == Date(2017, 7, 18) && m.forcing == :era5 && m.variant == "" && m.domain == "d25x100"
    @test m.nominal_width_km == 25 && m.spacing_m == 100
    @test m.microphysics == :sbmwrm && m.aerosol == :aer2 && m.surface == :flxsst
    @test lasso_samin_filename(m) == "enalasso_samin_20170718era5d25x100_sbmwrm-aer2-flxsstC1.m0.20170718.000000.tar"
    @test lasso_reference_filenames(m).samstat == "enalasso_samstat_20170718era5d25x100_sbmwrm-aer2-flxsstC1.m1.20170718.000000.nc"
    @test lasso_reference_filenames(m).sam2d == "enalasso_sam2d_20170718era5d25x100_sbmwrm-aer2-flxsstC1.m1.20170718.000000.nc"
    s1n0 = parse_lasso_member("20170718era5s1n0d25x100_morr-aer1-flxsst")
    @test lasso_variant_tokens(s1n0) == ["s1", "n0"] && s1n0.microphysics == :morr && s1n0.aerosol == :aer1
    @test isempty(lasso_variant_tokens(m))
    @test parse_lasso_member("20170718merra2d25x100_sbmwrm-aer3-flxsst").forcing == :merra2
    @test_throws ArgumentError parse_lasso_member("covert2022_bin")
    @test_throws ArgumentError parse_lasso_member("20170718era5d25x100_sbmwrm-aer4-flxsst")
    @test_throws ArgumentError parse_lasso_member("20170718era5d25x100_p3-aer2-flxsst")
    @test lasso_documented_dimensions(m) == (256, 256, 260)
    @test_throws ArgumentError lasso_documented_dimensions("20170718era5d102x100_sbmwrm-aer2-flxsst")
    @test BreezeLab.DEFAULT_LASSO_MEMBER == LASSO_MEMBER
    # the staging/download scripts and the package agree on the accepted ARM file names
    isdefined(Main, :ARMInputs) || include(joinpath(@__DIR__, "..", "data_wrangling", "fetch_arm_inputs.jl"))
    @test occursin(ARMInputs.ACCEPTED_FILENAME, lasso_samin_filename(m))
    @test occursin(ARMInputs.ACCEPTED_FILENAME, lasso_reference_filenames(m).samstat)
    @test occursin(ARMInputs.ACCEPTED_FILENAME, lasso_reference_filenames(m).sam2d)
    @test !occursin(ARMInputs.ACCEPTED_FILENAME, "enalasso_sam3d_20170718era5d25x100_sbmwrm-aer2-flxsstC1.m1.20170718.000000.nc")
    @test !occursin(ARMInputs.ACCEPTED_FILENAME, "enalasso_samrest_20170718era5d25x100_sbmwrm-aer2-flxsstC1.m1.20170718.000000.tar")
end

@testset "Namelist groups (setparm.f90, sgs.f90)" begin
    groups = read_sam_namelist_groups(joinpath(LASSO_FIXTURE, "prm"))
    @test Set(keys(groups)) == Set(["parameters", "sgs_tke", "micro_hujisbm"])
    p = groups["parameters"]
    @test p["caseid"] == "era5d25x100_sbmwrm-aer2-flxsst"
    @test p["sfc_flx_fxd"] === false && p["read_in_geostrophic_wind"] === true && p["ocean"] === true
    @test p["nrad"] == 30 && p["dt"] == 2.0 && p["nstop"] == 1800 && p["day0"] == 199.0 && p["tauls"] == 7200.0
    @test groups["sgs_tke"]["dosmagor"] === true
    @test isempty(groups["micro_hujisbm"])
    # the flat parser (all groups merged) still sees the same values
    flat = read_sam_namelist(joinpath(LASSO_FIXTURE, "prm"))
    @test flat["dosmagor"] === true && flat["latitude0"] == 39.0916
end

@testset "SAM grd levels and interfaces (setgrid.f90)" begin
    levels, extra = lasso_scalar_levels([12.5, 37.5, 62.5], 5)
    @test levels == [12.5, 37.5, 62.5, 87.5, 112.5] && extra == 2
    levels, extra = lasso_scalar_levels(collect(12.5:25:1000), 10)
    @test length(levels) == 10 && extra == 0
    @test_throws ArgumentError lasso_scalar_levels([12.5], 3)
    faces = faces_from_centers([12.5, 37.5, 62.5, 112.5])      # uneven last spacing
    @test faces == [0.0, 25.0, 50.0, 87.5, 137.5]              # zi(nz) = z(nzm) + (z(nzm) - zi(nzm))
end

@testset "Time-interpolated initial sounding (setdata.f90)" begin
    snd = read_sam_sounding(joinpath(LASSO_FIXTURE, "snd"))
    @test initial_sounding(snd, 199.0) === snd[1]
    @test initial_sounding(snd, 200.0) === snd[2]
    half = initial_sounding(snd, 199.5)
    @test half.day == 199.5 && half.θ[1] ≈ 292.7 && half.q[1] ≈ 11.1e-3 && half.p[1] == 104000
    @test all(isnan, half.z) && half.u == snd[1].u && half.surface_pressure == 101930
    @test_throws ErrorException initial_sounding(snd, 198.9)
    @test_throws ErrorException initial_sounding(snd, 200.1)
    @test initial_sounding(snd[1:1], 150.0) === snd[1]
end

@testset "Perturbation cases (setperturb.f90)" begin
    zc = collect(12.5:25:1000)
    δT, δq = perturbation_amplitudes(InitialPerturbation(), zc)
    @test all(δT[zc .≤ 600] .== 0.1) && all(δT[zc .> 600] .== 0) && all(δq[zc .≤ 600] .== 0.025e-3)
    δT0, δq0 = perturbation_amplitudes(InitialPerturbation(sam_perturb_type=0), zc)
    @test δT0[1:5] ≈ [0.1, 0.08, 0.06, 0.04, 0.02] && all(δT0[6:end] .== 0) && all(δq0 .== 0)
    ϵ = perturbation_array(3, 2, zc, InitialPerturbation(sam_perturb_type=0, seed=3))
    @test all(ϵ[:, :, 6:end] .== 0) && any(ϵ[:, :, 1:5] .!= 0)
    @test_throws ArgumentError perturbation_amplitudes(InitialPerturbation(sam_perturb_type=2), zc)
end

@testset "LASSO bundle inspection and validation" begin
    bundle = inspect_lasso_bundle(LASSO_FIXTURE; member=LASSO_MEMBER, dimensions=LASSO_DIMS, epoch=LASSO_EPOCH)
    @test isempty(bundle.problems)
    s = bundle.settings
    @test s.Nx == 16 && s.Ny == 16 && s.Lx == 1600.0 && length(s.z_faces) == 31
    @test s.Δt == 2.0 && s.max_Δt == 2.0 && s.stop_time == 3600.0
    @test s.radiation_interval == 60.0 && s.wind_nudging_timescale == 7200.0
    @test s.surface == :bulk_sst && s.radiation == :rrtmgp && s.vertical_advection == :full_field
    @test s.surface_emissivity == 0.95 && s.liquid_effective_radius == 14e-6
    @test s.surface_flux_law == :sam_oceflx && s.aerosol_supersaturation_cap == 0.003 && isnothing(s.coriolis_parameter)
    @test isempty(bundle.readme)
    @test s.microphysics == :p3_aer2 && s.aerosol_replenishment == :diagnostic_ccn
    @test s.perturbation.sam_perturb_type == 5 && s.translation_velocity == (0.0, 0.0)
    @test s.upper_boundary_relaxation && s.sponge isa SAMSponge
    @test s.latitude == 39.0916 && s.longitude == -28.0257 && s.day0 == 199.0
    @test bundle.grid.levels[1] == 12.5 && bundle.grid.top_level == 987.5 && bundle.grid.extrapolated == 0
    @test bundle.grid.faces[end] == 987.5 + (987.5 - 887.5) / 2
    @test bundle.grid.max_level_offset > 0                      # stretched levels are not face midpoints
    @test bundle.time.day_end ≈ 199.0 + 3600 / 86400 && bundle.time.radiation_interval == 60
    @test bundle.dimension_source == "explicit dimensions keyword (SAM domain.f90 is compiled in)"
    @test any(occursin("dosmagor = .true.", w) for w in bundle.warnings)
    @test any(occursin("14 μm", w) for w in bundle.warnings)
    @test any(occursin("diagCCN", w) for w in bundle.warnings)
    record = lasso_bundle_record(bundle)
    @test record.member.id == LASSO_MEMBER && record.sam_reference.commit == "12d02446a2147388dc89d828e6e0553106abea0f"
    @test length(record.checksums.prm) == 64 && record.doi == "10.5439/2572661"
    @test validate_lasso_bundle(LASSO_FIXTURE; member=LASSO_MEMBER, dimensions=LASSO_DIMS).settings.Ny == 16
    @test inspect_lasso_bundle(LASSO_FIXTURE; dimensions=LASSO_DIMS).settings.microphysics == :p3_aer2   # undeclared member keeps the production default
    @test_throws ArgumentError inspect_lasso_bundle(LASSO_FIXTURE; member=42)

    mktempdir() do dir
        lasso_copy(dir)
        # dimensions are required (SAM compiles them in) and must agree when the namelist has them
        @test occursin("dimensions=(Nx, Ny, Nz)", problem_text(dir))
        rewrite_prm(dir, "nstop = 1800" => "nstop = 1800, nx_gl = 16, ny_gl = 16, nz_gl = 30")
        @test isempty(inspect_lasso_bundle(dir; dimensions=LASSO_DIMS).problems)
        @test inspect_lasso_bundle(dir).dimension_source == "namelist nx_gl/ny_gl/nz_gl"
        @test occursin("disagree", problem_text(dir; dimensions=(32, 32, 30)))
        rewrite_prm(dir, "caseid ='era5d25x100_sbmwrm-aer2-flxsst'" => "caseid ='16x16x30'")
        @test inspect_lasso_bundle(dir).dimension_source == "numeric caseid 16x16x30"
        # every unsupported switch is named with the reason (absent keys take SAM's defaults)
        for (edit, needle) in (("donudging_tq = .false." => "donudging_tq = .true.", "donudging_tq = true"),
                               ("timelargescale = 0." => "timelargescale = 3600.", "timelargescale"),
                               ("SFC_TAU_FXD = .false." => "SFC_TAU_FXD = .true.", "sfc_tau_fxd"),
                               ("SFC_FLX_FXD = .false." => "SFC_FLX_FXD = .true.", "sfc_flx_fxd"),
                               ("doshortwave = .true." => "doshortwave = .false.", "doshortwave"),
                               ("doradsimple = .false." => "doradsimple = .true.", "doradsimple"),
                               ("perturb_type = 5" => "perturb_type = 2", "perturb_type = 2"),
                               ("READ_IN_GEOSTROPHIC_WIND = .true." => "READ_IN_GEOSTROPHIC_WIND = .false.", "ug/vg columns"),
                               ("nrad = 30" => "nrad = 30, nudging_uv_z1 = 500.", "nudging_uv_z1"),
                               ("nrad = 30" => "nrad = 30, nxco2 = 2", "nxco2"),
                               ("nrestart = 0" => "nrestart = 1", "nrestart"),
                               ("nrad = 30" => "nrad = 30, doperpetual = .true.", "doperpetual"),
                               ("nrad = 30" => "nrad = 30, LAND = .true.", "land"),
                               ("LES = .true." => "LES = .false.", "les"),
                               ("dosgs = .true., dodamping = .true." => "dosgs = .true., docolumn = .true., dodamping = .true.", "docolumn"))
            rewrite_prm(dir, edit)
            text = problem_text(dir; dimensions=LASSO_DIMS)
            @test occursin(needle, text)
        end
        # the official members compute per-column fluxes and carry fcor: both accepted, both recorded
        rewrite_prm(dir, "UNIFORM_SFC_FLX = .true." => "UNIFORM_SFC_FLX = .false.", "nrad = 30" => "nrad = 30, fcor = 9.19626e-05")
        bpc = inspect_lasso_bundle(dir; dimensions=LASSO_DIMS)
        @test isempty(bpc.problems) && bpc.settings.coriolis_parameter == 9.19626e-5 && !bpc.switches.uniform_sfc_flx
        @test any(occursin("per column", w) for w in bpc.warnings)
        @test !any(occursin("fcor", w) for w in bpc.warnings)                # the bundle value is 2Ω sin φ (sidereal): no deviation warning
        rewrite_prm(dir, "nrad = 30" => "nrad = 30, fcor = 9.17e-05")       # SAM's own 4π/86400 sin φ: accepted silently too
        @test !any(occursin("fcor", w) for w in inspect_lasso_bundle(dir; dimensions=LASSO_DIMS).warnings)
        rewrite_prm(dir, "nrad = 30" => "nrad = 30, fcor = 1.0e-04")        # detached from the latitude: accepted, warned
        bfc = inspect_lasso_bundle(dir; dimensions=LASSO_DIMS)
        @test isempty(bfc.problems) && bfc.settings.coriolis_parameter == 1e-4 && any(occursin("fcor = 0.0001 is used", w) for w in bfc.warnings)
        rewrite_prm(dir, "nrad = 30" => "nrad = 30, fcor = -999.")          # SAM's sentinel: derive from latitude
        @test isnothing(inspect_lasso_bundle(dir; dimensions=LASSO_DIMS).settings.coriolis_parameter)
        mkpath(joinpath(dir, "extra"))
        write(joinpath(dir, "extra", "README_LASSO-ENA.txt"), "sim_name: x\nmodel_type: SAM v6.10.3 plus LASSO modifications\nmodel_source_git_hash: f83adf58\n")
        @test inspect_lasso_bundle(dir; dimensions=LASSO_DIMS).readme["model_source_git_hash"] == "f83adf58"
        rm(joinpath(dir, "extra"); recursive=true)
        rewrite_prm(dir, "perturb_type = 5," => "")          # absent → SAM case 0, supported
        b0 = inspect_lasso_bundle(dir; dimensions=LASSO_DIMS)
        @test isempty(b0.problems) && b0.settings.perturbation.sam_perturb_type == 0
        rewrite_prm(dir, "compute_reffc = .false." => "compute_reffc = .true.")
        bre = inspect_lasso_bundle(dir; dimensions=LASSO_DIMS)
        @test bre.settings.liquid_effective_radius == 10e-6 && any(occursin("SBM spectrum", w) for w in bre.warnings)
        rewrite_prm(dir, "&SGS_TKE\n dosmagor = .true." => "&SGS_TKE\n dosmagor = .false.")
        @test any(occursin("1.5-order TKE", w) for w in inspect_lasso_bundle(dir; dimensions=LASSO_DIMS).warnings)
        rewrite_prm(dir, "nrad = 30, ug = 0.0, vg = 0.0" => "nrad = 30, ug = 5.0, vg = -8.0")
        bt = inspect_lasso_bundle(dir; dimensions=LASSO_DIMS)
        @test isempty(bt.problems) && bt.switches.sam_translation == (5.0, -8.0) && any(occursin("translation frame", w) for w in bt.warnings)
        # time coverage (setdata.f90 / forcing.f90 / setforcing.f90)
        rewrite_prm(dir, "day0 = 199.0" => "day0 = 198.5")
        text = problem_text(dir; dimensions=LASSO_DIMS)
        @test occursin("snd records", text) && occursin("first lsf record", text) && occursin("sfc samples", text)
        rewrite_prm(dir, "nstop = 1800" => "nstop = 86400")      # two days: sfc must cover the run
        @test occursin("sfc samples", problem_text(dir; dimensions=LASSO_DIMS))
        rewrite_prm(dir, "day0 = 199.0" => "day0 = 199.5")       # inside the records: fine
        @test isempty(inspect_lasso_bundle(dir; dimensions=LASSO_DIMS).problems)
        rewrite_prm(dir)
        # the epoch is checked against day0 and the member date, never inferred
        @test occursin("epoch", problem_text(dir; dimensions=LASSO_DIMS, epoch=DateTime(2017, 7, 18, 6)))
        @test isempty(inspect_lasso_bundle(dir; dimensions=LASSO_DIMS, epoch=DateTime(2018, 7, 18, 0)).problems)   # the year is the caller's
        @test occursin("different days", problem_text(dir; dimensions=LASSO_DIMS, member="20170719era5d25x100_sbmwrm-aer2-flxsst", epoch=LASSO_EPOCH))
        # grid: the sounding must cover the model top; short grd files are extended as SAM does
        write(joinpath(dir, "grd"), join(12.5:25:1000, '\n'))
        @test occursin("standard atmosphere", problem_text(dir; dimensions=(16, 16, 200)))
        b = inspect_lasso_bundle(dir; dimensions=(16, 16, 60))
        @test b.grid.extrapolated == 20 && b.grid.levels[end] == 12.5 + 25 * 59 && isempty(b.problems)
        @test any(occursin("extends", w) for w in b.warnings)
        write(joinpath(dir, "grd"), "125\n125\n250\n")
        @test occursin("strictly increasing", problem_text(dir; dimensions=(16, 16, 3)))
        write(joinpath(dir, "grd"), join(12.5:25:1000, '\n'))
        # a 7-column lsf with read_in_geostrophic_wind = .true. is a layout mismatch
        cp(joinpath(FIXTURES, "lsf_7col"), joinpath(dir, "lsf"); force=true)
        @test occursin("lacks the ug/vg columns", problem_text(dir; dimensions=LASSO_DIMS))
        # missing files are named; nothing is substituted
        rm(joinpath(dir, "grd"))
        @test_throws ArgumentError inspect_lasso_bundle(dir)
        @test occursin("lacks grd", sprint(showerror, caught(() -> validate_lasso_bundle(dir))))
    end
end

@testset "LASSO adapter: units, winds, forcing assembly" begin
    common = (; protocol=:lasso_ena_official, member=LASSO_MEMBER, epoch=LASSO_EPOCH, dimensions=LASSO_DIMS,
                FT=Float64, Nx=8, Ny=8, Lx=800, Ly=800, microphysics=:one_moment, radiation=nothing,
                aerosol_replenishment=nothing, write_output=false, progress_interval=100)
    @test_throws ArgumentError ena_simulation(LASSO_FIXTURE; common..., epoch=nothing)
    @test_throws ArgumentError ena_simulation(LASSO_FIXTURE; common..., member=nothing)
    @test_throws ArgumentError ena_simulation(LASSO_FIXTURE; common..., epoch=DateTime(2017, 7, 18, 6))
    @test_throws ArgumentError ena_simulation(LASSO_FIXTURE; common..., day0=199.5)
    @test_throws ArgumentError ena_simulation(LASSO_FIXTURE; common..., dimensions=nothing)

    case = ena_simulation(LASSO_FIXTURE; common..., stop_time=4.0)
    grid = case.grid
    zc = Array(znodes(grid, Center()))
    @test case.config.Δt_initial == 2.0 && case.config.max_Δt == 2.0 && case.simulation.stop_time == 4.0
    @test case.config.surface == "bulk_sst" && case.config.wind_nudging_timescale == 7200
    @test startswith(case.config.surface_flux_law, "sam_oceflx") && case.config.coriolis_parameter_source == "4π/86400 sin(latitude)"
    @test case.config.coriolis_parameter ≈ FPlane(Float64; latitude=39.0916).f rtol=1e-6    # 2Ω sin φ, sidereal Ω
    @test case.config.minimum_wind_speed == 1.0 && case.config.gustiness == 0.0 && case.config.fit_error < 0.01
    @test case.config.surface_emissivity == 0.95 && case.config.liquid_effective_radius == 14e-6
    @test case.protocol_member == LASSO_MEMBER && occursin("member $LASSO_MEMBER", case.config.label)
    @test case.bundle.member.aerosol == "aer2"
    @test case.config.epoch == string(LASSO_EPOCH) && case.config.day0 == 199.0
    # the vertical grid is the grd file's: 30 levels, uniform 25 m below 600 m
    @test grid.Nz == 30 && zc[1] == 12.5 && zc[24] == 587.5
    @test isapprox(zc[end], 987.5; atol=1e-9)                  # top face at 987.5 + 50, so the top center is the SAM level
    # initial sounding: day0 is the first record; moisture is the dry mixing ratio converted to a mass fraction
    @test case.sounding === case.soundings[1]
    r₁ = interpolate_profile(record_heights(case.soundings[1]), case.soundings[1].q, zc[1])     # file mixing ratio at z(1)
    @test case.columns.qᵗ[1] ≈ r₁ / (1 + r₁) rtol=1e-12
    @test 0 < case.columns.qᵗ[1] < r₁
    # pressure → height conversion (setdata.f90 with dry T, R = 287, cp = 1004, g = 9.81)
    snd = case.soundings[1]
    T₁ = snd.θ[1] * (snd.p[1] / 1e5)^(287 / 1004); T₂ = snd.θ[2] * (snd.p[2] / 1e5)^(287 / 1004)
    z₁ = 287 / 9.81 * T₁ * log(snd.surface_pressure / snd.p[1])
    @test record_heights(snd)[1] ≈ z₁ && z₁ < 0
    @test record_heights(snd)[2] ≈ z₁ + 287 / (2 * 9.81) * (T₁ + T₂) * log(snd.p[1] / snd.p[2])
    # winds are ground-relative: the sounding wind is the initial wind, the nudging target is uls itself
    @test !case.config.translation_frame_applied && case.config.translation_velocity_u == 0
    @test Array(interior(case.model.velocities.u, 1, 1, :)) ≈ case.columns.u
    # forcing profiles at several heights and times reproduce the analytic fixture (linear in z, scale 1 → 2 over the day)
    fp = case.forcing_profiles
    for (k, t) in ((1, 0.0), (10, 43200.0), (24, 86400.0))
        z = zc[k]; s = 1 + t / 86400
        @test fp.tls[1, 1, k, Time(t)] ≈ s * (-1e-5 - 2e-9 * z) rtol=1e-9
        @test fp.qls[1, 1, k, Time(t)] ≈ s * (-1e-8 + 1e-12 * z) rtol=1e-9
        @test fp.wls[1, 1, k, Time(t)] ≈ -0.002 * z / 1000 rtol=1e-9
        @test fp.uls[1, 1, k, Time(t)] ≈ 5 + 0.002 * z rtol=1e-9
        @test fp.vls[1, 1, k, Time(t)] ≈ -3
        @test fp.ug[1, 1, k, Time(t)] ≈ 6 + 0.001 * z rtol=1e-9      # separate geostrophic columns
        @test fp.vg[1, 1, k, Time(t)] ≈ -4
    end
    @test fp.has_geostrophic_columns
    # the SST series feeds the shared surface temperature field (sfc: 294.9 K at day0, 294.95 K at +12 h)
    @test case.surface_series.sst[1, 1, 1, Time(0.0)] ≈ 294.9 && case.surface_series.sst[1, 1, 1, Time(43200.0)] ≈ 294.95
    # dry mixing-ratio source → mass-fraction rate: dq/dt = R / (1 + r)² in clear air
    r = 0.0111; R = -1e-8
    q = MoistureMassFractions(r / (1 + r), 0.0, 0.0)
    @test BreezeLab.moisture_rates(Val(:mixing_ratio), q, R).vapor ≈ R / (1 + r)^2 rtol=1e-12
    @test BreezeLab.moisture_rates(Val(:mass_fraction), q, R).vapor == R
    # forcing assembly: each term once, no mean-profile subsidence, subsidence on the model's own θ
    f = case.forcing
    count_of(T, terms) = count(x -> x isa T, terms)
    @test count_of(LargeScaleVerticalAdvection, f.θ) == 1 && count_of(LargeScaleVerticalAdvection, f.E) == 0
    @test count_of(LargeScaleEnergyForcing, f.E) == 1 && count_of(UpperBoundaryEnergyRelaxation, f.E) == 1
    # the one-moment moisture prognostic is the equilibrium total water qᵉ (P3: qᵛ)
    @test count_of(LargeScaleMoistureForcing, f.qᵉ) == 1 && count_of(LargeScaleVerticalAdvection, f.qᵉ) == 1 &&
          count_of(UpperBoundaryMoistureRelaxation, f.qᵉ) == 1
    @test count_of(TimeVaryingGeostrophicForcing, f.u) == 1 && count_of(MeanProfileNudging, f.u) == 1 &&
          count_of(LargeScaleVerticalAdvection, f.u) == 1 && count_of(SAMSponge, f.u) == 1
    @test count_of(SAMSponge, f.w) == 1 && !any(x -> x isa Breeze.SubsidenceForcing, (f.u..., f.θ..., f.qᵉ...))
    @test count_of(LargeScaleVerticalAdvection, f.qʳ) == 1    # every microphysical field is advected (subsidence.f90)
    # run and record provenance
    run!(case.simulation)
    @test all(x -> all(isfinite, interior(x)), values(Oceananigans.prognostic_fields(case.model)))
    mktempdir() do dir
        record = TOML.parsefile(write_provenance(joinpath(dir, "provenance.toml"), case))
        @test record["protocol"] == "lasso_ena_official" && record["protocol_member"] == LASSO_MEMBER
        @test record["bundle"]["member"]["samin"] == lasso_samin_filename(parse_lasso_member(LASSO_MEMBER))
        @test record["bundle"]["sam_reference"]["commit"] == "12d02446a2147388dc89d828e6e0553106abea0f"
        @test haskey(record["bundle"]["checksums"], "grd") && record["inputs"]["grd_sha256"] == record["bundle"]["checksums"]["grd"]
        @test record["bundle"]["time"]["day_end"] ≈ 199.0 + 3600 / 86400
        @test record["protocol_overrides"]["stop_time"] == 4.0
        @test any(occursin("diagCCN", w) for w in record["bundle"]["warnings"])
    end

    # a namelist fcor becomes the f-plane parameter
    mktempdir() do dir
        lasso_copy(dir)
        rewrite_prm(dir, "nrad = 30" => "nrad = 30, fcor = 9.19626e-05")
        withf = ena_simulation(dir; common..., stop_time=2.0)
        @test withf.model.coriolis.f ≈ 9.19626e-5 && withf.config.coriolis_parameter_source == "namelist fcor"
    end
    # a namelist translation frame is recorded but not applied with bulk fluxes (ground-relative winds)
    mktempdir() do dir
        lasso_copy(dir)
        rewrite_prm(dir, "nrad = 30, ug = 0.0, vg = 0.0" => "nrad = 30, ug = 5.0, vg = -8.0")
        moving = ena_simulation(dir; common..., stop_time=2.0)
        @test moving.config.sam_translation_u == 5.0 && moving.config.sam_translation_v == -8.0
        @test !moving.config.translation_frame_applied
        @test Array(interior(moving.model.velocities.u, 1, 1, :)) ≈ moving.columns.u
        @test moving.forcing_profiles.uls[1, 1, 1, Time(0.0)] ≈ 5 + 0.002 * zc[1]
        @test_throws ArgumentError ena_simulation(dir; common..., translation_velocity=(5.0, -8.0))
        # perturbation case 0 of a namelist without perturb_type reaches the initial state
        rewrite_prm(dir, "perturb_type = 5," => "")
        case0 = ena_simulation(dir; common..., stop_time=2.0)
        @test case0.config.perturbation == string(InitialPerturbation(sam_perturb_type=0))
        T = Array(interior(case0.model.temperature))
        @test maximum(abs, T[:, :, 1] .- case0.columns.T[1]) ≤ 0.1 + 1e-3 && maximum(abs, T[:, :, 1] .- case0.columns.T[1]) > 1e-3
        @test all(isapprox.(T[:, :, 7], case0.columns.T[7]; rtol=1e-6))      # no perturbation above level 5
    end
end

@testset "LASSO protocol defaults build with RRTMGP on the CPU" begin
    # The production radiation path (RRTMGP LW+SW, SST surface temperature, SAM emissivity,
    # 14 μm effective radius, nrad*dt schedule) is constructed and updated once on a tiny grid;
    # the Covert CPU tests use :simple radiation and never reached this branch.
    case = ena_simulation(LASSO_FIXTURE; protocol=:lasso_ena_official, member=LASSO_MEMBER, epoch=LASSO_EPOCH,
                          dimensions=LASSO_DIMS, FT=Float32, Nx=8, Ny=8, Lx=800, Ly=800, microphysics=:one_moment,
                          aerosol_replenishment=nothing, write_output=false, progress_interval=100)
    @test case.config.radiation == "rrtmgp" && case.config.radiation_interval == 60.0
    @test case.config.surface_emissivity == 0.95 && case.config.liquid_effective_radius == 14e-6
    @test !isnothing(case.model.radiation)
    @test case.model.radiation.solar_position.epoch == LASSO_EPOCH
    Oceananigans.TimeSteppers.update_state!(case.model)
    @test all(isfinite, interior(case.model.radiation.flux_divergence))
end

@testset "SAM oceflx surface law" begin
    laws = sam_oceflx_neutral_polynomials()
    cdn(U) = 0.0027 / U + 0.000142 + 0.0000764 * U
    @test laws.drag == (0.000142, 0.0000764, 0.0027) && laws.fit_error < 0.01
    ev(p, U) = p[1] + p[2] * U + p[3] / U
    for U in (2.0, 3.0, 7.0, 12.0, 15.0)
        @test ev(laws.sensible, U) ≈ 0.0327 * sqrt(cdn(U)) rtol=0.01
        @test ev(laws.latent, U) ≈ 0.0346 * sqrt(cdn(U)) rtol=0.01
    end
    @test ev(laws.sensible, 1.0) ≈ 0.0327 * sqrt(cdn(1.0)) rtol=0.1     # outside the fitted range
    # SAM's scalar coefficients are ~10 % above Breeze's Large & Yeager defaults at 7 m/s
    @test ev(laws.sensible, 7.0) > ev((1.28e-4, 6.8e-5, 2.43e-3), 7.0) * 1.05
    grid = test_grid(; Nz=8, Lz=400)
    Tₛ = Field{Center, Center, Nothing}(grid); set!(Tₛ, 295.0)
    bcs, record = bulk_surface_flux_boundary_conditions(grid, Tₛ; moisture_name=:qᵛ, law=:sam_oceflx)
    @test record.minimum_wind_speed == 1.0 && record.gustiness == 0.0 && startswith(record.surface_flux_law, "sam_oceflx")
    @test bcs.ρE.bottom.condition.coefficient.polynomial == laws.sensible
    @test bcs.ρqᵛ.bottom.condition.coefficient.polynomial == laws.latent
    @test bcs.ρu.bottom.condition.coefficient.polynomial == laws.drag
    bcs_b, record_b = bulk_surface_flux_boundary_conditions(grid, Tₛ; moisture_name=:qᵛ)
    @test record_b.gustiness == 0.1 && isnothing(bcs_b.ρE.bottom.condition.coefficient.polynomial)
    @test_throws ArgumentError bulk_surface_flux_boundary_conditions(grid, Tₛ; moisture_name=:qᵛ, law=:other)
end

@testset "ena_lasso constructor wrapper" begin
    err = caught(() -> ena_lasso(joinpath(LASSO_FIXTURE, "absent"); member=LASSO_MEMBER, epoch=LASSO_EPOCH))
    @test err isa ArgumentError
    message = sprint(showerror, err)
    @test occursin("enalasso_samin_20170718era5d25x100_sbmwrm-aer2-flxsstC1.m0.20170718.000000.tar", message)
    @test occursin("stage_lasso_bundle.jl", message) && occursin("never substituted", message)
    @test lasso_bundle_available(LASSO_FIXTURE) && !lasso_bundle_available(FIXTURES)
    withenv("ENA_LASSO_BUNDLE" => nothing) do
        @test endswith(lasso_bundle_directory(LASSO_MEMBER), joinpath("data", "lasso", LASSO_MEMBER))
    end
    withenv("ENA_LASSO_BUNDLE" => LASSO_FIXTURE) do
        @test lasso_bundle_directory(LASSO_MEMBER) == LASSO_FIXTURE
    end
    # the documented d25x100 domain disagrees with the 16×16×30 fixture (grd has 30 levels)
    @test_throws ArgumentError ena_lasso(LASSO_FIXTURE; member=LASSO_MEMBER, epoch=LASSO_EPOCH, arch=CPU())
    mktempdir() do staged
        lasso_copy(staged)
        write(joinpath(staged, "bundle.toml"), "member = \"$LASSO_MEMBER\"\narchive_sha256 = \"test\"\n")
        output = joinpath(staged, "output")
        case = ena_lasso(staged; member=LASSO_MEMBER, epoch=LASSO_EPOCH, dimensions=LASSO_DIMS, arch=CPU(), FT=Float64,
                         Nx=8, Ny=8, Lx=800, Ly=800, microphysics=:one_moment, radiation=nothing, aerosol_replenishment=nothing,
                         stop_time=2.0, timeseries_interval=1.0, profile_interval=2.0, slice_interval=2.0,
                         output_dir=output, progress_interval=100)
        @test case.staging["archive_sha256"] == "test" && case.dimension_request == "(16, 16, 30)"
        @test case.model.clock.iteration == 0
        record = TOML.parsefile(write_provenance(joinpath(output, "provenance.toml"), case))
        @test record["staging"]["member"] == LASSO_MEMBER
        run!(case.simulation)
        series = breezelab_timeseries(output)
        @test series.protocol == "lasso_ena_official" && series.epoch == LASSO_EPOCH
        @test length(series.time) == 3 && series.time[end] == LASSO_EPOCH + Second(2)
        @test all(isfinite, series.lwp) && all(isfinite, series.rain_rate)
        bounds = breezelab_cloud_boundaries(output)
        @test length(bounds.time) ≥ 1 && bounds.time[end] == LASSO_EPOCH + Second(2)
    end
end
