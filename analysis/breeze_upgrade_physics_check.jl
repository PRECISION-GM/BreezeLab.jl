# Before/after physics check for a Breeze upgrade (docs/breeze_0.12_upgrade.md).
#
# Builds each case at a small size with fixed settings, steps it a fixed number of times on the
# CPU, and writes domain-mean diagnostics every `BREEZE_CHECK_EVERY` steps to
# `<BREEZE_CHECK_OUTPUT>/<case>.csv`. Run the same file under the old and the new environment and
# compare the CSVs with `analysis/breeze_upgrade_physics_compare.jl`:
#
#   julia --project=<old checkout> analysis/breeze_upgrade_physics_check.jl
#   julia --project=<new checkout> analysis/breeze_upgrade_physics_check.jl
#
# Environment: BREEZE_CHECK_OUTPUT (required), BREEZE_CHECK_CASES (comma list, default all),
# BREEZE_CHECK_STEPS (default 50), BREEZE_CHECK_EVERY (default 5), BREEZE_CHECK_FLOAT (Float32),
# BREEZE_CHECK_NX (default 16), and the input locations COVERT_DIR, LASSO_DIR, SEA_STARR_DIR, TRACER_IOP.
#
# Besides the cases, two closed "box" experiments isolate the microphysics: a resting, horizontally
# uniform cloudy layer with no forcing, no surface fluxes and no radiation, stepped with the
# one-moment (saturation adjustment) and the P3 scheme. Their total water and temperature drifts
# separate the scheme-level changes from the case-level ones; the P3 box also starts with rain in
# the lowest levels so that the sign of `bottom_precipitation_flux` can be checked against
# BreezeLab's `surface_rain_flux` (both documented as positive downward).

using BreezeLab, Breeze, Oceananigans, Oceananigans.Units, Dates, Printf, Statistics
using Oceananigans.Fields: interior
using Oceananigans: compute!
using Breeze.Microphysics.PredictedParticleProperties: CloudDroplets
using CloudMicrophysics: CloudMicrophysics

const OUT = ENV["BREEZE_CHECK_OUTPUT"]
const STEPS = parse(Int, get(ENV, "BREEZE_CHECK_STEPS", "50"))
const EVERY = parse(Int, get(ENV, "BREEZE_CHECK_EVERY", "5"))
const NXY = parse(Int, get(ENV, "BREEZE_CHECK_NX", "16"))
const CASES = split(get(ENV, "BREEZE_CHECK_CASES", "box_one_moment,box_p3,covert_one_moment,covert_p3_n75,lasso_default,dp_scream,sea_starr_ctrl"), ',')
Oceananigans.defaults.FloatType = get(ENV, "BREEZE_CHECK_FLOAT", "Float32") == "Float64" ? Float64 : Float32
mkpath(OUT)

const COVERT_DIR = get(ENV, "COVERT_DIR", joinpath(pkgdir(BreezeLab), "data", "covert2022_bin"))
const LASSO_DIR = get(ENV, "LASSO_DIR", "/shared/home/greg/breezelab-work/ena-lasso/data/lasso/20170718era5d25x100_sbmwrm-aer2-flxsst")
const SEA_STARR_DIR = get(ENV, "SEA_STARR_DIR", "/shared/home/greg/breezelab-runs/20261006/inputs/mip_sources/seastarr_22241697")
const TRACER_IOP = get(ENV, "TRACER_IOP", "/shared/home/greg/breezelab-work/tracer-dp-scream/data/tracer_dp_scream/TRACER_iopfile_4scam.nc")

one_moment() = Base.get_extension(Breeze, :BreezeCloudMicrophysicsExt).OneMomentCloudMicrophysics(;
                   cloud_formation = SaturationAdjustment(; equilibrium = WarmPhaseEquilibrium()))
p3(N = 75e6; aerosol = nothing) = P3Microphysics(; cloud = CloudDroplets(; number_concentration = N), aerosol)

#####
##### Diagnostics (version-agnostic: only fields both Breeze versions provide)
#####

field_or_nothing(nt, name) = haskey(nt, name) ? nt[name] : nothing
cpu(f) = Array(interior(f))

function diagnostics(model)
    grid = model.grid
    μ = model.microphysical_fields
    ρᵣ = Array(interior(model.dynamics.reference_state.density))           # (1, 1, Nz)
    Δz = reshape([Oceananigans.Operators.Δzᶜᶜᶜ(1, 1, k, grid) for k in 1:size(grid, 3)], 1, 1, :)
    colsum(q) = dropdims(sum(ρᵣ .* q .* Δz, dims = 3), dims = 3)            # kg m⁻²
    qᵛ = cpu(μ.qᵛ)
    qᶜˡ = cpu(μ.qᶜˡ)
    qʳ = cpu(μ.qʳ)
    qⁱ = isnothing(field_or_nothing(μ, :qⁱ)) ? zero(qᶜˡ) : cpu(μ.qⁱ)
    T = cpu(model.temperature)
    f = Oceananigans.fields(model)
    θ = haskey(f, :θ) ? cpu(f.θ) : haskey(f, :ρθ) ? cpu(f.ρθ) ./ ρᵣ : fill(NaN, size(T))
    cloudy = qᶜˡ .> 1e-5
    nᶜˡ_cm3 = if haskey(μ, :ρnᶜˡ)
        n = cpu(μ.ρnᶜˡ) ./ 1e6                                               # cm⁻³
        any(cloudy) ? mean(n[cloudy]) : 0.0
    else
        NaN
    end
    nʳ = haskey(μ, :ρnʳ) ? mean(cpu(μ.ρnʳ)) : NaN                           # m⁻³
    rain = mean(cpu(surface_rain_flux(model))) * 86400                     # mm day⁻¹ (BreezeLab, positive down)
    bottom = try
        flux = Breeze.AtmosphereModels.bottom_precipitation_flux(model)
        compute!(flux)
        mean(cpu(flux)) * 86400
    catch
        NaN
    end
    water = mean(colsum(qᵛ .+ qᶜˡ .+ qʳ .+ qⁱ))                              # kg m⁻², from the diagnosed partition
    # Prognostic water: every prognostic moisture density (ρqᵉ or ρqᵛ, and the condensates), not
    # the diagnosed partition. P3's rime mass ρqᶠ is part of ρqⁱ and is not counted twice.
    prognostic = Oceananigans.prognostic_fields(model)
    water_names = [name for name in keys(prognostic) if startswith(string(name), "ρq") && name !== :ρqᶠ]
    Δz₃ = Δz .* ones(size(T))
    water_prognostic = sum(sum(cpu(prognostic[name]) .* Δz₃) for name in water_names) / (size(T, 1) * size(T, 2))
    ρθ_column = haskey(prognostic, :ρθ) ? sum(cpu(prognostic.ρθ) .* Δz₃) / (size(T, 1) * size(T, 2)) : NaN
    return (; t = model.clock.time, iteration = model.clock.iteration,
              LWP = 1e3 * mean(colsum(qᶜˡ)), RWP = 1e3 * mean(colsum(qʳ)), IWP = 1e3 * mean(colsum(qⁱ)),
              water, water_prognostic, ρθ_column, mean_T = mean(T), mean_θ = mean(θ), mean_qv = 1e3 * mean(qᵛ),
              cloud_fraction = mean(dropdims(sum(ρᵣ .* qᶜˡ .* Δz, dims = 3), dims = 3) .> 5e-3),
              Nc_cloudy = nᶜˡ_cm3, mean_nr = nʳ, max_qr = 1e3 * maximum(qʳ),
              rain_mm_day = rain, bottom_flux_mm_day = bottom,
              max_w = maximum(abs, cpu(model.velocities.w)))
end

function write_row(io, d; header = false)
    header && println(io, join(string.(keys(d)), ","))
    println(io, join((x isa AbstractFloat ? @sprintf("%.10g", x) : string(x) for x in values(d)), ","))
end

function step_and_record(name, model, Δt)
    path = joinpath(OUT, name * ".csv")
    open(path, "w") do io
        write_row(io, diagnostics(model); header = true)
        wall = time()
        for n in 1:STEPS
            time_step!(model, Δt)
            if n % EVERY == 0 || n == STEPS
                d = diagnostics(model)
                write_row(io, d)
                flush(io)
                all(isfinite, (d.LWP, d.mean_T, d.water, d.water_prognostic)) || (@warn "$name: non-finite at step $n"; break)
            end
        end
        @info @sprintf("%s: %d steps of %.2f s in %.1f s wall", name, STEPS, Δt, time() - wall)
    end
    return path
end

# Steps a case's model with the case's own callbacks (nudging masks, SST, solar position, ...),
# by running its Simulation for exactly STEPS fixed steps.
function run_case(name, case, Δt)
    simulation = case.simulation
    simulation.Δt = Δt
    simulation.stop_iteration = 0
    simulation.stop_time = Inf
    # Remove the adaptive time-step wizard so that both environments take identical steps.
    for key in collect(keys(simulation.callbacks))
        cb = simulation.callbacks[key]
        occursin("TimeStepWizard", string(typeof(cb.func))) && pop!(simulation.callbacks, key)
    end
    path = joinpath(OUT, name * ".csv")
    open(path, "w") do io
        write_row(io, diagnostics(case.model); header = true)
        wall = time()
        for n in EVERY:EVERY:STEPS
            simulation.stop_iteration = n
            simulation.running = true
            run!(simulation)
            d = diagnostics(case.model)
            write_row(io, d)
            flush(io)
            all(isfinite, (d.LWP, d.mean_T, d.water, d.water_prognostic)) || (@warn "$name: non-finite at step $n"; break)
        end
        @info @sprintf("%s: %d steps of %.2f s in %.1f s wall", name, STEPS, Δt, time() - wall)
    end
    return path
end

#####
##### Closed box: resting, horizontally uniform cloudy layer, nothing but microphysics and sedimentation
#####

function box_model(microphysics)
    grid = RectilinearGrid(CPU(); size = (2, 2, 40), x = (0, 200), y = (0, 200), z = (0, 2000),
                           topology = (Periodic, Periodic, Bounded))
    constants = ThermodynamicConstants()
    reference_state = ReferenceState(grid, constants; base_pressure = 101500, potential_temperature = 288)
    dynamics = AnelasticDynamics(reference_state)
    model = AtmosphereModel(grid; dynamics, microphysics, formulation = :LiquidIcePotentialTemperature,
                            thermodynamic_constants = constants)
    return model
end

θ_box(z) = 288 + (z > 1200 ? 8 + 0.006 * (z - 1200) : 0.0)
qᵗ_box(z) = z > 1200 ? 0.004 : 0.0085                    # cloudy from the ≈ 400 m condensation level to the 1200 m inversion

function box_one_moment()
    model = box_model(one_moment())
    set!(model; θ = (x, y, z) -> θ_box(z), qᵗ = (x, y, z) -> qᵗ_box(z))
    return model
end

function box_p3()
    model = box_model(p3())
    pᵣ = Array(interior(model.dynamics.reference_state.pressure, 1, 1, :))
    z = Array(znodes(model.grid, Center()))
    T = similar(z); qᵛ = similar(z); qᶜˡ = similar(z)
    for k in eachindex(z)
        T[k], qᵛ[k], qᶜˡ[k] = BreezeLab.saturation_partition(θ_box(z[k]), qᵗ_box(z[k]), pᵣ[k];
                                                             constants = ThermodynamicConstants(Float64))
    end
    column(v) = repeat(reshape(v, 1, 1, :), 2, 2, 1)
    qʳ = [zk < 200 ? 2e-4 : 0.0 for zk in z]               # rain in the lowest levels: bottom-flux sign check
    nʳ = [zk < 200 ? 2e4 : 0.0 for zk in z]                # ≈ 0.6 mm drops
    set!(model; T = column(T), qᵛ = column(qᵛ), qᶜˡ = column(qᶜˡ), qʳ = column(qʳ), nʳ = column(nʳ))
    return model
end

#####
##### Cases
#####

function build(name)
    common = (; arch = CPU(), write_output = false)
    if name == "box_one_moment"
        return box_one_moment(), 1.0
    elseif name == "box_p3"
        return box_p3(), 1.0
    elseif name == "covert_one_moment" || name == "covert_p3_n75"
        microphysics = name == "covert_one_moment" ? one_moment() : p3(75e6)
        case = ena_covert(; common..., data_dir = COVERT_DIR, microphysics, Nx = NXY, Ny = NXY, Lx = 100NXY, Ly = 100NXY,
                            Δt = 1.0, max_Δt = 1.0, stop_time = 1e6, progress_interval = 1e9)
        return case, 1.0
    elseif name == "lasso_default"
        case = ena_lasso(; common..., bundle_dir = LASSO_DIR, member = "20170718era5d25x100_sbmwrm-aer2-flxsst",
                           epoch = DateTime(2017, 7, 18), dimensions = (256, 256, 260),
                           Nx = NXY, Ny = NXY, Lx = 100NXY, Ly = 100NXY, stop_time = 1e6, progress_interval = 1e9)
        return case, 1.0
    elseif name == "dp_scream"
        case = tracer_dp_scream(; common..., iop_path = TRACER_IOP, start = DateTime(2022, 8, 1), stop = DateTime(2022, 8, 15),
                                  Nx = NXY, Ny = NXY, Δt = 2.0, max_Δt = 2.0, stop_time = 1e6, progress_interval = 1e9)
        return case, 2.0
    elseif name == "sea_starr_ctrl"
        case = sea_starr(; common..., member = :CTRL, data_dir = SEA_STARR_DIR, Nx = NXY, Ny = NXY, stop_time = 1e6,
                           checkpoint_interval = nothing, progress_interval = 1e9)
        return case, 1.0
    end
    error("unknown case $name")
end

for name in CASES
    @info "Building $name with Breeze $(pkgversion(Breeze)) ($(Oceananigans.defaults.FloatType))"
    try
        built, Δt = build(name)
        if built isa NamedTuple
            run_case(name, built, Δt)
        else
            step_and_record(name, built, Δt)
        end
    catch err
        @error "$name failed" exception = (err, catch_backtrace())
        open(io -> println(io, sprint(showerror, err)), joinpath(OUT, name * ".error"), "w")
    end
end
