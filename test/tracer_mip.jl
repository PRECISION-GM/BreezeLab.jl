# TRACER-MIP protocol metadata, grids and the prescribed aerosol profiles. CPU only; no data.
using BreezeLab: tracer_mip_protocol, tracer_mip_case_window, acpc_vertical_faces, tracer_mip_horizontal_extent,
                 tracer_mip_grid, tracer_mip_aerosol_profile, aerosol_shape, surface_number_mixing_ratios,
                 mode_number_mixing_ratios, total_number_mixing_ratio, aerosol_profile_table, tracer_mip_p3_aerosol_modes
using Breeze.Microphysics.PredictedParticleProperties: AerosolMode, AerosolActivation
using Breeze.TerrainFollowingDiscretization: TerrainFollowingGrid

@testset "TRACER-MIP protocol" begin
    protocol = tracer_mip_protocol()
    start, stop = tracer_mip_case_window(protocol, :aug07)
    @test start == DateTime(2022, 8, 7, 6) && stop == DateTime(2022, 8, 8, 6)
    @test tracer_mip_case_window(protocol, :jun17)[1] == DateTime(2022, 6, 17, 6)
    @test protocol["grids"]["outer"]["timestep_s"] == 3.0 && protocol["grids"]["inner"]["timestep_s"] == 1.5
    @test protocol["grids"]["inner"]["tracking_output_interval_s"] == 120.0
    @test protocol["forcing"]["radiation_interval_s"] == 60.0

    @testset "ACPC vertical faces" begin
        faces = acpc_vertical_faces()
        @test length(faces) == 95                       # 94 cells from 95 scalar levels
        @test faces[1] == 0 && faces[2] == 50 && faces[3] == 103
        @test faces[end] == 22181 && faces[end-1] == 21881
        @test issorted(faces)
        Δz = diff(faces)
        @test minimum(Δz) == 50 && maximum(Δz) ≈ 300.5
        @test all(Δz[end-10:end] .== 300)               # uniform 300 m in the upper troposphere
        centers = (faces[1:end-1] .+ faces[2:end]) ./ 2
        zt = protocol["vertical"]["scalar_levels_m_agl"][2:end]
        @test maximum(abs.(centers .- zt)) ≤ 5          # cell centers track the scalar levels (≤ 3.5 m on the stretched part)
        @test_throws ArgumentError acpc_vertical_faces([-24, -10, 10])
    end

    @testset "horizontal extents" begin
        outer = tracer_mip_horizontal_extent(protocol, :outer)
        inner = tracer_mip_horizontal_extent(protocol, :inner)
        φ₀, λ₀ = 29.4719, -95.0792
        @test sum(outer.latitude) / 2 ≈ φ₀ && sum(outer.longitude) / 2 ≈ λ₀
        @test sum(inner.latitude) / 2 ≈ φ₀ && sum(inner.longitude) / 2 ≈ λ₀
        m_per_deg = π * 6371e3 / 180
        @test (outer.latitude[2] - outer.latitude[1]) * m_per_deg ≈ 1500e3
        @test (outer.longitude[2] - outer.longitude[1]) * m_per_deg * cosd(φ₀) ≈ 1500e3
        @test (inner.latitude[2] - inner.latitude[1]) * m_per_deg ≈ 250e3
        # the lat-lon box lies within the roadmap's polar-stereographic corner envelope
        @test 22.552 ≤ outer.latitude[1] && outer.latitude[2] ≤ 35.941 + 0.3
        @test -103.376 ≤ outer.longitude[1] && outer.longitude[2] ≤ -86.782
        @test 28.348 - 0.01 ≤ inner.latitude[1] && inner.latitude[2] ≤ 30.602 + 0.01   # 250 km box vs. the stereographic corners (±1 km)
        @test -96.378 ≤ inner.longitude[1] && inner.longitude[2] ≤ -93.759
    end

    @testset "grid construction (reduced size)" begin
        grid = tracer_mip_grid(CPU(), :outer; Nx = 8, Ny = 6, FT = Float64)
        @test size(grid) == (8, 6, 94)
        @test grid isa TerrainFollowingGrid
        @test grid.Lz ≈ 22181
        flat = tracer_mip_grid(CPU(), :inner; Nx = 5, Ny = 5, terrain_following = false, FT = Float32)
        @test size(flat) == (5, 5, 94) && eltype(flat) == Float32
        @test Oceananigans.Grids.topology(flat) == (Bounded, Bounded, Bounded)
    end
end

@testset "TRACER-MIP prescribed aerosol profiles (notebook port)" begin
    aug = tracer_mip_aerosol_profile(:aug07)
    jun = tracer_mip_aerosol_profile(:jun17)

    @testset "units and surface values" begin
        n₀ = surface_number_mixing_ratios(aug)
        @test n₀[1] ≈ 1425e6 / 1.159 && n₀[2] ≈ 263e6 / 1.159     # kg⁻¹
        @test collect(1e-6 .* n₀) ≈ [1229.5, 226.9] atol = 0.06       # mg⁻¹, Table 3 / notebook
        @test collect(1e-6 .* surface_number_mixing_ratios(jun)) ≈ [2970.7, 458.2] atol = 0.06
        @test aerosol_shape(aug, 0.0) == 1
        @test aug.modes[1].median_diameter ≈ 49e-9 && aug.modes[2].sigma == 1.4
        @test jun.modes[1].sigma == 1.5 && jun.modes[2].median_diameter ≈ 136e-9
    end

    # Values printed by the reference notebook (mg⁻¹) at its 66-level sampling (km).
    notebook_aug = [(0.000, 1229.5, 226.9), (0.338, 1054.1, 194.5), (1.015, 810.9, 149.7), (1.692, 559.1, 103.2),
                    (2.031, 354.7, 65.5), (2.369, 189.6, 35.0), (3.046, 82.8, 15.3), (3.723, 51.3, 9.5),
                    (4.062, 42.2, 7.8), (5.754, 42.2, 7.8), (6.092, 42.2, 7.8), (22.000, 42.2, 7.8)]
    notebook_jun = [(0.000, 2970.7, 458.2), (0.677, 2218.8, 342.2), (1.354, 1710.9, 263.9), (2.031, 857.1, 132.2),
                    (2.708, 273.6, 42.2), (4.400, 79.3, 12.2), (5.077, 52.6, 8.1), (5.415, 43.6, 6.7),
                    (5.754, 43.3, 6.7), (10.154, 43.3, 6.7)]
    heights_km = range(0, 22, length = 66)
    @testset "notebook table reproduction" begin
        for (profile, rows) in ((aug, notebook_aug), (jun, notebook_jun))
            for (zkm, n1, n2) in rows
                z = heights_km[argmin(abs.(heights_km .- zkm))] * 1000      # the notebook's own sample height
                n = 1e-6 .* mode_number_mixing_ratios(profile, z)
                @test n[1] ≈ n1 atol = 0.15
                @test n[2] ≈ n2 atol = 0.15
            end
        end
    end

    @testset "floor, modal ratio and monotonicity" begin
        for profile in (aug, jun)
            n₀ = surface_number_mixing_ratios(profile)
            r = n₀ ./ sum(n₀)
            for z in (0.0, 500.0, 2000.0, 4500.0, 6001.0, 12000.0)
                n = mode_number_mixing_ratios(profile, z)
                @test collect(n ./ sum(n)) ≈ collect(r)                 # modal ratio preserved everywhere
                @test total_number_mixing_ratio(profile, z) ≥ 50e6 - 1e-6  # 50 mg⁻¹ total floor
            end
            @test total_number_mixing_ratio(profile, 6001.0) ≈ 50e6
            @test total_number_mixing_ratio(profile, 20000.0) ≈ 50e6
            zs = 0:100:6000
            totals = [total_number_mixing_ratio(profile, z) for z in zs]
            @test issorted(totals; rev = true)
            table = aerosol_profile_table(profile, [0.0, 1000.0])
            @test table[1].total ≈ 1e-6 * sum(n₀) && length(table[2].modes) == 2
        end
    end

    @testset "sensitivity multipliers" begin
        high = tracer_mip_aerosol_profile(:aug07; multiplier = 3)
        @test collect(surface_number_mixing_ratios(high)) ≈ 3 .* collect(surface_number_mixing_ratios(aug))
        @test total_number_mixing_ratio(high, 10000.0) ≈ 50e6                  # floor not scaled by default
        scaled = tracer_mip_aerosol_profile(:aug07; multiplier = 3, scale_floor = true)
        @test total_number_mixing_ratio(scaled, 10000.0) ≈ 150e6
        low = tracer_mip_aerosol_profile(:jun17; multiplier = 0.3)
        @test 1e-6 * surface_number_mixing_ratios(low)[1] ≈ 0.3 * 2970.7 atol = 0.05
    end

    @testset "P3 modes with κ = 0.26" begin
        modes = tracer_mip_p3_aerosol_modes(aug; FT = Float64)
        @test length(modes) == 2 && all(m -> m isa AerosolMode{Float64}, modes)
        @test modes[1].number_mixing_ratio ≈ 1425e6 / 1.159
        @test modes[1].mean_radius ≈ 24.5e-9 && modes[2].mean_radius ≈ 87.5e-9
        @test modes[1].geometric_std == 1.8 && modes[2].geometric_std == 1.4
        @test all(m -> m.solute_activity ≈ 0.26, modes)                     # β_act reproduces the bulk κ
        activation = AerosolActivation(modes...; prognostic = true)
        @test length(activation.modes) == 2
        modes32 = tracer_mip_p3_aerosol_modes(jun; FT = Float32)
        @test eltype(modes32) == AerosolMode{Float32} && modes32[1].solute_activity ≈ 0.26f0
    end
end
