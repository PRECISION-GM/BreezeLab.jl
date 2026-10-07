# Station-series readers and window statistics on synthetic NetCDF files that mimic the
# ARM ENA products' variable names, units and QC conventions (generated here; no ARM data).
using NCDatasets
using Dates: DateTime, Second

@testset "ARM observation readers and window statistics" begin
    mktempdir() do dir
        t0 = DateTime(2017, 7, 18, 6)
        t1 = DateTime(2017, 7, 18, 12)
        units = "seconds since 2017-07-18 00:00:00 0:00"

        # MWRRET2turn: phys_lwp with bit-packed qc and the retrieval flag
        mwr = joinpath(dir, "enamwrret2turnC1.c1.20170718.031053.nc")
        NCDataset(mwr, "c"; attrib=Dict("datastream" => "enamwrret2turnC1.c1")) do ds
            defDim(ds, "time", 6)
            defVar(ds, "time", Float64[21600, 21660, 21720, 25200, 28800, 46800], ("time",); attrib=Dict("units" => units))
            defVar(ds, "phys_lwp", Float32[50, 60, 999, 999, 70, 80], ("time",); attrib=Dict("units" => "g/m^2"))
            defVar(ds, "qc_phys_lwp", Int32[0, 0, 4, 0, 0, 0], ("time",))
            defVar(ds, "phys_qc_flag", Int16[0, 0, 0, 1, 0, 0], ("time",))
            defVar(ds, "phys_lwp_uncertainty", Float32[5, 5, -9999, 5, 5, 5], ("time",); attrib=Dict("missing_value" => -9999f0))
            defVar(ds, "stat_lwp", Float32[55, 65, -9999, 75, 85, 95], ("time",); attrib=Dict("units" => "g/m^2", "missing_value" => -9999f0))
        end
        lwp = read_arm_lwp(mwr)
        @test lwp.datastream == "enamwrret2turnC1.c1" && lwp.units == "g/m^2" && length(lwp) == 6
        @test lwp.good == BitVector([1, 1, 0, 0, 1, 1]) && isnan(lwp.uncertainty[3])
        stats = window_statistics(lwp, t0, t1)
        @test stats.n == 5 && stats.n_good == 3 && stats.fraction_good == 0.6
        @test stats.mean ≈ 60 && stats.median == 60 && stats.min == 50 && stats.max == 70
        @test isnan(window_statistics(lwp, DateTime(2017, 7, 19), DateTime(2017, 7, 20)).mean)
        stat = read_arm_lwp(mwr; variable="stat_lwp")
        @test stat.good == BitVector([1, 1, 0, 1, 1, 1]) && stat.uncertainty === nothing && isnan(stat.value[3])

        # ARSCL KAZR cloud boundaries: layers × time, -1 for clear sky
        arscl = joinpath(dir, "enaarsclkazrbnd1kolliasC1.c0.20170718.000000.nc")
        NCDataset(arscl, "c"; attrib=Dict("datastream" => "enaarsclkazrbnd1kolliasC1.c0")) do ds
            defDim(ds, "time", 5); defDim(ds, "layer", 2)
            defVar(ds, "time", Float64[21600, 21604, 21608, 21612, 43200], ("time",); attrib=Dict("units" => units))
            base = Float32[800 -1 900 4000 850; 2500 -1 -1 -1 -1]
            top = Float32[1200 -1 1300 6000 1250; 3000 -1 -1 -1 -1]
            defVar(ds, "cloud_layer_base_height", base, ("layer", "time"); attrib=Dict("units" => "m", "missing_value" => -9999f0))
            defVar(ds, "cloud_layer_top_height", top, ("layer", "time"); attrib=Dict("units" => "m", "missing_value" => -9999f0))
            defVar(ds, "cloud_base_best_estimate", Float32[790, -1, -2, 3900, 840], ("time",); attrib=Dict("units" => "m"))
        end
        bounds = read_arm_cloud_boundaries(arscl)
        @test bounds.base[1] == 800 && isnan(bounds.base[2]) && bounds.top[3] == 1300 && bounds.n_layers == [2, 0, 1, 1, 1]
        @test bounds.cloudy == BitVector([1, 0, 1, 1, 1]) && isnan(bounds.base_best_estimate[3])
        cf = cloud_fraction_from_boundaries(bounds, t0, t1; max_base=3000)
        @test cf.n == 4 && cf.n_cloudy == 2 && cf.cloud_fraction == 0.5          # the 4 km layer is above max_base
        @test cf.base.mean ≈ 850 && cf.top.mean ≈ 1250

        # VDIS rain rate with qc
        vdis = joinpath(dir, "enavdisC1.b1.20170718.000000.cdf")
        NCDataset(vdis, "c"; attrib=Dict("datastream" => "enavdisC1.b1")) do ds
            defDim(ds, "time", 4)
            defVar(ds, "time", Float64[21600, 21660, 21720, 21780], ("time",); attrib=Dict("units" => units))
            defVar(ds, "rain_rate", Float32[0.0, 0.4, -9999, 0.2], ("time",); attrib=Dict("units" => "mm/hr", "missing_value" => -9999f0))
            defVar(ds, "qc_rain_rate", Int32[0, 0, 0, 2], ("time",))
        end
        rain = read_arm_rain_rate(vdis)
        @test rain.units == "mm/hr" && rain.good == BitVector([1, 1, 0, 0]) && isnan(rain.value[3])
        @test window_statistics(rain, t0, t1).mean ≈ 0.2
    end
end
