# # TRACER – DP-SCREAM-matched periodic case
#
# This example builds, runs and analyzes a doubly periodic Breeze LES driven by the
# forcing of the published DP-SCREAM TRACER experiment (Oware et al. 2025,
# doi 10.1029/2025JD044113): the ARM VARANAL IOP file `TRACER_iopfile_4scam.nc` used by
# the E3SM scmlib script `run_dpxx_scream_TRACER.csh`. The archived DP-SCREAM August run
# started at 2022-08-01 00 UTC; the paper's 5–15 August window is its days 5–15, so the
# baseline reproduces that initialization and continuous forcing history. The LES
# (Smagorinsky–Lilly at 200 m on 51.2 km, 160 levels to 22 km, P3 with prescribed
# droplet number, prescribed IOP heat fluxes with bulk drag, RRTMGP, no Coriolis, no
# wind nudging — see `tracer_dp_scream` for the evidence) is an intentional
# physics/resolution departure from SHOC at 3.33 km on 200 km, not a DP-SCREAM reproduction.
#
# Fetch the pinned inputs and reference outputs with
# `julia data_wrangling/fetch_manifest.jl cases/tracer_dp_scream/inputs.toml data/tracer_dp_scream`,
# then run `julia --project cases/tracer_dp_scream.jl` on an NVIDIA GPU. The Slurm wrapper
# `execution/tracer_dp_scream.sbatch` sets the environment variables read below.

using BreezeLab
using BreezeLab: day_to_seconds, interpolate_profile
using Breeze
using Breeze.Microphysics.PredictedParticleProperties: CloudDroplets
using Oceananigans, Oceananigans.Units
using Dates: Dates, DateTime, Second, Hour
using CairoMakie
using Statistics: mean
using TOML: TOML
using Printf: @sprintf

# ## Settings
#
# The default run covers the full window (14 days). A short integration sets
# `TRACER_DP_SCREAM_STOP_TIME` (seconds); `TRACER_DP_SCREAM_PICKUP_FROM` names a previous
# output directory whose latest checkpoint is continued (new output files, same forcing).

data_dir = joinpath(pkgdir(BreezeLab), "data", "tracer_dp_scream")
start = DateTime(2022, 8, 1)
stop = DateTime(2022, 8, 15)
full_duration = Dates.value(stop - start) / 1000
stop_time = parse(Float64, get(ENV, "TRACER_DP_SCREAM_STOP_TIME", string(full_duration)))
checkpoint_interval = parse(Float64, get(ENV, "TRACER_DP_SCREAM_CHECKPOINT", "21600"))
output_dir = get(ENV, "TRACER_DP_SCREAM_OUTPUT", joinpath(pkgdir(BreezeLab), "output", "tracer_dp_scream"))
pickup_from = get(ENV, "TRACER_DP_SCREAM_PICKUP_FROM", "")
Oceananigans.defaults.FloatType = Float32
arch = GPU()
# Horizontal grid: the baseline uses 256² cells of 200 m (51.2 km). `TRACER_DP_SCREAM_NX`
# and `TRACER_DP_SCREAM_DX` override both so that a 512² × 100 m run covers the same domain.
Nx = Ny = parse(Int, get(ENV, "TRACER_DP_SCREAM_NX", "256"))
Δx = parse(Float64, get(ENV, "TRACER_DP_SCREAM_DX", "200"))
# Fields used for animation (the xy maps at the slice height, the xz section, LWP and rain
# maps) are saved at least four times finer than the 30-min profile means: every
# `TRACER_DP_SCREAM_SLICE_MINUTES` (default 5) at `TRACER_DP_SCREAM_SLICE_HEIGHT` m (default
# 3000, near the afternoon cloud base of this case; 1500 m sat below it in the 200 m baseline).
# Profiles (30 min) and time series (60 s) are unchanged.
slice_minutes = parse(Float64, get(ENV, "TRACER_DP_SCREAM_SLICE_MINUTES", "5"))
slice_height = parse(Float64, get(ENV, "TRACER_DP_SCREAM_SLICE_HEIGHT", "3000"))
slice_interval = slice_minutes * 60

# ## Build the simulation

case = tracer_dp_scream(; iop_path = joinpath(data_dir, "TRACER_iopfile_4scam.nc"),
                        start, stop, stop_time, arch, Nx, Ny, Δx,
                        microphysics = P3Microphysics(; cloud = CloudDroplets(; number_concentration = 200e6)),
                        output_dir, output_prefix = "tracer", slice_interval, slice_height,
                        checkpoint_interval = checkpoint_interval > 0 ? checkpoint_interval : nothing)
simulation = case.simulation;

# Storage estimate for the slice file: three xz sections (Nx × Nz) and four xy maps (Nx × Ny)
# in Float32 per save; job 202 (256² × 160, 30-min cadence) measured 0.76 of that interior
# size per save (1.175 MB), which calibrates the estimate recorded in the provenance.
Nz = case.config.Nz
slice_saves = floor(Int, stop_time / slice_interval) + 1
slice_bytes_per_save = 0.76 * 4 * (3 * Nx * Nz + 4 * Nx * Ny)
slice_file_estimate_GB = slice_saves * slice_bytes_per_save / 1e9
@info @sprintf("slice file estimate: %d saves every %.0f min at %.0f m ≈ %.1f GB (%.2f MB per save)",
               slice_saves, slice_minutes, slice_height, slice_file_estimate_GB, slice_bytes_per_save / 1e6)

# ## Run
#
# Provenance (inputs with checksums, pinned sources, configuration, overrides) is written
# next to the output before the run starts. A pickup copies the latest checkpoint of the
# previous segment into this output directory; Oceananigans restores the prognostic state
# and clock, the time-step wizard re-adapts Δt within a few steps, and the surface
# temperature callback refreshes Tg after the first step.

mkpath(output_dir)
write_provenance(joinpath(output_dir, "provenance.toml"), case;
                 extra = (; stop_time, checkpoint_interval, pickup_from,
                            slice_minutes, slice_height, slice_saves, slice_bytes_per_save, slice_file_estimate_GB,
                            slice_estimate_basis = "0.76 × Float32 interior bytes of 3 xz + 4 xy fields per save, calibrated on job 202",
                            dp_scream_reference = "Zenodo 10.5281/zenodo.15271730 (3 km August run, case scream_dp_TRACER_AUGUST_decr; 0.5 km 5–15 August run)",
                            scmlib_script = "DPxx_SCREAM_SCRIPTS/run_dpxx_scream_TRACER.csh @ 2dc3f1073a5d03b5f32617e79d68fb20a63cfb43"))
pickup = false
if !isempty(pickup_from)
    checkpoints = filter(f -> startswith(f, "tracer_checkpoint") && endswith(f, ".jld2"), readdir(pickup_from))
    isempty(checkpoints) && error("no tracer_checkpoint*.jld2 in $pickup_from")
    latest = checkpoints[argmax([parse(Int, match(r"iteration(\d+)", f).captures[1]) for f in checkpoints])]
    cp(joinpath(pickup_from, latest), joinpath(output_dir, latest); force=true)
    @info "Picking up from $latest"
    pickup = true
end
run!(simulation; pickup)

# ## Analysis: time series against the archived DP-SCREAM runs
#
# The archived 3 km August run (which contains the paper's 5–15 August period) and the
# cold-started 0.5 km run are read with their local-time labels converted to UTC. The LES
# liquid/ice water paths are kg m⁻² like `TGCLDLWP`/`TGCLDIWP`; the LES surface rain and
# ice fluxes (kg m⁻² s⁻¹) are converted to mm day⁻¹ like the archive's `PRECL`
# (whose unit is inferred, see `read_dp_scream_output`). Horizontal means of a 51.2 km LES
# are compared with 200 km domain means: no spatial filtering is needed for domain means,
# but their sampling of mesoscale organization differs.

ts_file = joinpath(output_dir, "tracer_timeseries.jld2")
read_series(name) = (fts = FieldTimeSeries(ts_file, name); (fts.times, [fts[n][1, 1, 1] for n in eachindex(fts.times)]))
t, lwp = read_series("lwp")
_, iwp = read_series("iwp")
_, rwp = read_series("rwp")
_, pw = read_series("precipitable_water")
_, rain = read_series("rain_flux")
_, ice = read_series("ice_flux")
utc = start .+ Second.(round.(Int, t))
hours = t ./ hour

references = Dict{String, Any}()
for (key, file) in (("DP-SCREAM 3 km (August run)", "DP-SCREAM 3km_August.nc"),
                    ("DP-SCREAM 0.5 km (5–15 Aug cold start)", "DP-SCREAM 0.5km_August_5_to_August_15_run.nc"))
    path = joinpath(data_dir, file)
    isfile(path) && (references[key] = read_dp_scream_output(path))
end

fig = Figure(size = (1000, 900))
ax1 = Axis(fig[1, 1], ylabel = "Liquid water path (g m⁻²)")
ax2 = Axis(fig[2, 1], ylabel = "Ice water path (g m⁻²)")
ax3 = Axis(fig[3, 1], ylabel = "Surface precipitation (mm day⁻¹)", xlabel = "Hours since $(start) UTC")
lines!(ax1, hours, 1e3 .* lwp, label = "Breeze LES 200 m")
lines!(ax2, hours, 1e3 .* iwp, label = "Breeze LES 200 m")
lines!(ax3, hours, 86400 .* (rain .+ ice), label = "Breeze LES 200 m (rain + ice)")
for (key, ref) in references
    sel = dp_scream_window(ref, utc[1], utc[end])
    isempty(sel) && continue
    th = [Dates.value(ref.time[n] - start) / 3.6e6 for n in sel]
    lines!(ax1, th, 1e3 .* ref.TGCLDLWP[sel], label = key)
    lines!(ax2, th, 1e3 .* ref.TGCLDIWP[sel], label = key)
    lines!(ax3, th, ref.PRECL[sel], label = key)
end
axislegend(ax1, position = :lt)
save(joinpath(output_dir, "water_paths_and_precipitation.png"), fig)

# ## Time–height cloud fraction
#
# The LES profile writer stores 30-minute means of the total (liquid + ice) cloudy-cell
# fraction; the archive stores `TOT_CLOUD_FRAC` on hybrid levels with their heights `Z3`.

pf = joinpath(output_dir, "tracer_profiles.jld2")
cf = FieldTimeSeries(pf, "total_cloud_fraction")
z = Array(znodes(cf.grid, Center()))
cf_matrix = [cf[n][1, 1, k] for k in eachindex(z), n in eachindex(cf.times)]
fig2 = Figure(size = (1100, 800))
ax = Axis(fig2[1, 1], ylabel = "Height (km)", title = "Breeze LES total cloud fraction")
hm = heatmap!(ax, cf.times ./ hour, z ./ 1e3, permutedims(cf_matrix), colormap = :Blues, colorrange = (0, 1))
Colorbar(fig2[1, 2], hm)
if haskey(references, "DP-SCREAM 3 km (August run)")
    ref = references["DP-SCREAM 3 km (August run)"]
    sel = dp_scream_window(ref, utc[1], utc[end])
    th = [Dates.value(ref.time[n] - start) / 3.6e6 for n in sel]
    zref = vec(mean(ref.Z3[:, sel], dims=2)) ./ 1e3
    ascending = sortperm(zref)
    ax2 = Axis(fig2[2, 1], ylabel = "Height (km)", xlabel = "Hours since $(start) UTC", title = "DP-SCREAM 3 km TOT_CLOUD_FRAC")
    hm2 = heatmap!(ax2, th, zref[ascending], permutedims(ref.TOT_CLOUD_FRAC[ascending, sel]), colormap = :Blues, colorrange = (0, 1))
    ylims!(ax2, 0, 20)
    Colorbar(fig2[2, 2], hm2)
end
ylims!(ax, 0, 20)
save(joinpath(output_dir, "cloud_fraction_time_height.png"), fig2)

# ## Water and energy budget checks
#
# Column water: d(PW + LWP + RWP + IWP)/dt = E − P + S_q, with E the prescribed surface vapor
# flux, P the surface rain + ice flux and S_q = ∫ρ qls dz the imposed large-scale vapor
# source (from the model's own forcing profiles). Column static energy:
# dS/dt ≈ F_E + ∫(−∂F_rad/∂z) dz + ∫ρ cᵖ tls dz, with F_E the prescribed energy flux;
# the enthalpy carried out by precipitation is not included, so this check is approximate.
# Residuals are reported relative to the integrated source magnitudes.

_, S = read_series("column_static_energy")
_, R = read_series("column_radiative_heating")
ρᵣ = Array(interior(case.model.dynamics.reference_state.density, 1, 1, :))
Δz = diff(Array(znodes(case.grid, Face())))
# the forcing profiles are held piecewise constant, so the record at or before `tt` is exact
function column_source(fts, tt)
    times = Array(fts.times)        # the FieldTimeSeries clock lives on the GPU; no scalar indexing
    n = max(1, searchsortedlast(times, tt))
    column = Array(interior(fts[n], 1, 1, :))
    return sum(ρᵣ[k] * column[k] * Δz[k] for k in eachindex(Δz))
end
constants = case.model.thermodynamic_constants
ℒ = constants.liquid.reference_latent_heat
sfc = case.sfc
surface_value(values, tt) = interpolate_profile(day_to_seconds.(sfc.day, case.config.day0), values, tt)
cᵖᵈ = constants.dry_air.heat_capacity; cᵖᵛ = constants.vapor.heat_capacity
trapz(f) = sum((f[n] + f[n+1]) / 2 * (t[n+1] - t[n]) for n in 1:length(t)-1)
E = [surface_value(sfc.latent_heat_flux, tt) / ℒ for tt in t]
Sq = [column_source(case.forcing_profiles.qls, tt) for tt in t]
P = rain .+ ice
W = pw .+ lwp .+ rwp .+ iwp
water_residual = (W[end] - W[1]) - trapz(E .- P .+ Sq)
FE = [surface_value(sfc.sensible_heat_flux, tt) + (cᵖᵛ - cᵖᵈ) * surface_value(sfc.sst, tt) * E[n] for (n, tt) in enumerate(t)]
St = [cᵖᵈ * column_source(case.forcing_profiles.tls, tt) for tt in t]
energy_residual = (S[end] - S[1]) - trapz(FE .+ R .+ St)
budget = Dict("water" => Dict("delta_total_water_kg_m2" => W[end] - W[1], "integrated_evaporation" => trapz(E),
                              "integrated_precipitation" => trapz(P), "integrated_large_scale_source" => trapz(Sq),
                              "residual_kg_m2" => water_residual,
                              "residual_relative_to_sources" => abs(water_residual) / max(trapz(abs.(E) .+ abs.(P) .+ abs.(Sq)), 1e-12)),
              "energy" => Dict("delta_column_static_energy_J_m2" => S[end] - S[1], "integrated_surface_energy_flux" => trapz(FE),
                               "integrated_radiative_heating" => trapz(R), "integrated_large_scale_heating" => trapz(St),
                               "residual_J_m2" => energy_residual,
                               "residual_relative_to_sources" => abs(energy_residual) / max(trapz(abs.(FE) .+ abs.(R) .+ abs.(St)), 1e-12),
                               "note" => "precipitation enthalpy export not included"),
              "final_time_seconds" => t[end], "all_finite" => all(isfinite, W) && all(isfinite, S))
open(joinpath(output_dir, "budget.toml"), "w") do io
    TOML.print(io, budget)
end
@info "Water budget residual $(water_residual) kg m⁻² (relative $(budget["water"]["residual_relative_to_sources"])); energy residual $(energy_residual) J m⁻² (relative $(budget["energy"]["residual_relative_to_sources"]))"
fig #src
