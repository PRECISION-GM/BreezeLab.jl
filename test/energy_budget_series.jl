# Surface turbulent fluxes, net radiative fluxes at the surface and LES top, and the in-cloud
# droplet number written by the ENA (LASSO and Covert) output writers.
using BreezeLab: surface_heat_fluxes, radiative_boundary_fluxes, in_cloud_droplet_number, interpolate_time_series
using Oceananigans.Fields: interior
using Oceananigans: compute!
using JLD2: jldopen

domain_mean(f) = (compute!(f); sum(interior(f)) / length(interior(f)))

const ENERGY_SERIES = (:surface_sensible_heat_flux, :surface_latent_heat_flux, :surface_net_longwave,
                       :top_net_longwave, :in_cloud_droplet_number, :cloudy_volume_fraction)

@testset "LASSO energy-budget and droplet-number series" begin
    output = mktempdir()
    case = with_float_type(Float32) do
        ena_lasso(; member=LASSO_MEMBER, epoch=LASSO_EPOCH, bundle_dir=LASSO_FIXTURE, dimensions=LASSO_DIMS, arch=CPU(),
                    Nx=8, Ny=8, Lx=800, Ly=800, stop_time=2.0, timeseries_interval=1.0, profile_interval=2.0,
                    slice_interval=2.0, output_dir=output, progress_interval=100)
    end
    model = case.model
    names = keys(case.simulation.output_writers[:timeseries].outputs)
    @test all(∈(names), ENERGY_SERIES)
    @test all(∈(names), (:surface_net_shortwave, :top_net_shortwave, :top_downwelling_shortwave))
    @test :cloudy_droplet_number ∈ keys(case.simulation.output_writers[:profiles].outputs)
    run!(case.simulation)

    # radiation: SAM sign conventions from Breeze's (up > 0, down < 0) face fluxes
    radiation = model.radiation
    Nz = size(model.grid, 3)
    r = radiative_boundary_fluxes(radiation, model.grid)
    lw_up, lw_dn = Array(interior(radiation.upwelling_longwave_flux)), Array(interior(radiation.downwelling_longwave_flux))
    sw_up, sw_dn = Array(interior(radiation.upwelling_shortwave_flux)), Array(interior(radiation.downwelling_shortwave_flux))
    @test domain_mean(r.surface_net_longwave) ≈ mean(lw_up[:, :, 1] .+ lw_dn[:, :, 1]) rtol=1e-5
    @test domain_mean(r.top_net_longwave) ≈ mean(lw_up[:, :, Nz+1] .+ lw_dn[:, :, Nz+1]) rtol=1e-5
    @test domain_mean(r.surface_net_shortwave) ≈ -mean(sw_up[:, :, 1] .+ sw_dn[:, :, 1]) atol=1e-3
    @test domain_mean(r.top_downwelling_shortwave) ≈ -mean(sw_dn[:, :, Nz+1]) atol=1e-3
    @test domain_mean(r.surface_net_longwave) > 0 && domain_mean(r.top_net_longwave) > domain_mean(r.surface_net_longwave)

    # bulk surface fluxes from the SST: finite, the latent flux upward over the warm ocean,
    # and the sensible flux is cᵖᵐ Π Jᶿ of the ρθ boundary flux
    fluxes = surface_heat_fluxes(model)
    shf, lhf = domain_mean(fluxes.sensible), domain_mean(fluxes.latent)
    @test isfinite(shf) && isfinite(lhf) && lhf > 0
    @test abs(shf) < 200 && lhf < 500

    # in-cloud droplet number: volume-weighted mean of ρᵣ nᶜˡ over cells with qᶜˡ > 1e-5
    μ = model.microphysical_fields
    stats = in_cloud_droplet_number(model)
    q = zeros(Float32, size(μ.qᶜˡ)); n = fill(5f7, size(μ.nᶜˡ))
    q[1:4, :, 10:12] .= 2f-4; n[1:4, :, 10:12] .= 8f7; q[5, 1, 3] = 5f-6       # one sub-threshold cell
    set!(μ.qᶜˡ, q); set!(μ.nᶜˡ, n)
    ρ = Array(interior(model.dynamics.reference_state.density, 1, 1, :))
    Δz = diff(collect(case.bundle.grid.faces))
    cloudy = q .> 1f-5
    weights = cloudy .* reshape(Δz, 1, 1, :)
    expected = sum(weights .* reshape(ρ, 1, 1, :) .* n) / sum(weights)
    @test domain_mean(stats.mean) ≈ expected rtol=1e-4
    @test domain_mean(stats.cloudy_fraction) ≈ sum(weights) / (8 * 8 * sum(Δz)) rtol=1e-4
    compute!(stats.profile)
    @test Array(interior(stats.profile))[1, 1, 11] ≈ 0.5 * ρ[11] * 8f7 rtol=1e-4
    set!(μ.qᶜˡ, 0)
    @test domain_mean(stats.mean) == 0                       # no cloud → 0, not NaN

    # the series reach the time-series file
    jldopen(joinpath(output, "ena_lasso_timeseries.jld2")) do file
        @test all(name -> haskey(file, "timeseries/$name"), ENERGY_SERIES)
        @test length(keys(file["timeseries/surface_latent_heat_flux"])) ≥ 3
    end
end

if HAVE_COVERT
    @testset "Covert energy-budget series recover the prescribed sfc fluxes" begin
        z_faces = collect(range(0, 6000, length=25))
        case = with_float_type(Float64) do
            ρ₁ = first_level_reference_density(COVERT_DIR, z_faces)
            ena_covert(; arch=CPU(), data_dir=COVERT_DIR, Nx=8, Ny=8, z_faces,
                         microphysics=p3_microphysics(; aerosol=lasso_aerosol(; reference_density=ρ₁)),
                         stop_time=1.0, write_output=true, output_dir=mktempdir(), progress_interval=100)
        end
        names = keys(case.simulation.output_writers[:timeseries].outputs)
        @test all(∈(names), ENERGY_SERIES) && :surface_net_shortwave ∉ names   # rad_simple: longwave only
        sfc = case.sfc
        times = [day_to_seconds(d, case.config.day0) for d in sfc.day]
        H = interpolate_time_series(times, sfc.sensible_heat_flux, 0.0)
        LE = interpolate_time_series(times, sfc.latent_heat_flux, 0.0)
        # the prescribed energy flux is the sfc file's H and the latent flux LE
        fluxes = surface_heat_fluxes(case.model)
        @test domain_mean(fluxes.latent) ≈ LE rtol=1e-6
        @test domain_mean(fluxes.sensible) ≈ H rtol=1e-3
        r = radiative_boundary_fluxes(case.model.radiation, case.model.grid)
        F = Array(interior(case.model.radiation.flux))
        @test domain_mean(r.surface_net_longwave) ≈ mean(F[:, :, 1]) && domain_mean(r.top_net_longwave) ≈ mean(F[:, :, end])

        # the RRTMGP longwave-only sensitivity: SST surface, no shortwave, recorded as an override
        @test_throws ArgumentError ena_covert(; arch=CPU(), data_dir=COVERT_DIR, Nx=8, Ny=8, z_faces, radiation=:rrtmgp,
                                               write_output=false)
        rrtmgp = with_float_type(Float32) do
            ena_covert(; arch=CPU(), data_dir=COVERT_DIR, Nx=8, Ny=8, z_faces, radiation=:rrtmgp_longwave,
                         stop_time=2.0, write_output=true, output_dir=mktempdir(), progress_interval=100)
        end
        @test rrtmgp.model.radiation isa Breeze.AtmosphereModels.RadiativeTransferModel
        @test "radiation" ∈ rrtmgp.config.overrides && occursin("longwave only", rrtmgp.config.radiation)
        @test :surface_net_shortwave ∈ keys(rrtmgp.simulation.output_writers[:timeseries].outputs)
        run!(rrtmgp.simulation)
        rr = radiative_boundary_fluxes(rrtmgp.model.radiation, rrtmgp.model.grid)
        @test domain_mean(rr.surface_net_shortwave) == 0 && domain_mean(rr.top_downwelling_shortwave) == 0
        @test domain_mean(rr.top_net_longwave) > domain_mean(rr.surface_net_longwave) > 0
    end
end
