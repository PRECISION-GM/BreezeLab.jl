# Resume a SEA STARR run from a checkpoint into a NEW output directory (one directory per segment).
#
# Run with the project of the commit the run STARTED from (a git worktree at that commit), e.g.
#   julia --project=<worktree> execution/sea_starr_resume.jl
# Environment (as cases/sea_starr.jl, plus the resume settings):
#   SEA_STARR_MEMBER, SEA_STARR_NX, SEA_STARR_HOURS, SEA_STARR_DATA, SEA_STARR_CHECKPOINT_HOURS,
#   SEA_STARR_TWOD_MINUTES (only if the original run set it), SEA_STARR_OUTPUT (the new segment directory),
#   SEA_STARR_RESUME_CHECKPOINT  checkpoint file to restore,
#   SEA_STARR_RESUME_FROM        the original run directory (read only: run.log for Δt, time series for continuity),
#   SEA_STARR_ARCH = gpu | cpu,  SEA_STARR_RESUME_STEPS (optional; stop after this many steps — CPU pre-check).
#
# Why not `run!(simulation; pickup = true)`: the pinned Oceananigans checkpoints `prognostic_state(simulation)`,
# whose generic fallback serializes every callback's `func` whole (BreezeLab's mask, evaporation-rate, SST and
# solar updaters, including P3 tables and, for the radiation-fix runs, the RRTMGP model); no restore method exists
# for those objects. This driver restores only what is state: the model (clock, prognostic fields, closure fields,
# time stepper) and the output-writer schedules/averages. It then
#   * sets Δt to the last value logged before the checkpoint (`Simulation.Δt` is not checkpointed; the CFL wizard
#     takes over at its next actuation),
#   * runs `update_state!` and the stateless callbacks (SST, inversion mask, evaporation rate, trajectory solar
#     position) so the first restarted step sees forcing consistent with the restored state rather than the
#     constructor's t = 0 state (the ENA campaign found that a bare pickup used build-time callback state),
#   * radiation: the radiative-transfer schedule is not checkpointed; its TimeInterval fires on the first
#     `update_state!` after the restore, so fluxes are recomputed from the restored state before the first step.

using BreezeLab
using Breeze
using Oceananigans
using Oceananigans.Units
using Oceananigans.Fields: interior
using Oceananigans.OutputWriters: load_nested_data, reconcile_restored_output_schedule!
using Oceananigans: restore_prognostic_state!
using Oceananigans.TimeSteppers: update_state!
using CUDA
using JLD2
using Dates
using Printf
using Statistics

env(k, default = nothing) = (v = get(ENV, k, ""); isempty(v) ? default : v)

member = Symbol(env("SEA_STARR_MEMBER", "CTRL"))
arch = lowercase(env("SEA_STARR_ARCH", "gpu")) == "cpu" ? CPU() : GPU()
Nx = Ny = parse(Int, env("SEA_STARR_NX", "192"))
stop_time = parse(Float64, env("SEA_STARR_HOURS", "66")) * hours
data_dir = env("SEA_STARR_DATA")
output_dir = env("SEA_STARR_OUTPUT")
checkpoint_interval = parse(Float64, env("SEA_STARR_CHECKPOINT_HOURS", "3")) * hours
checkpoint = env("SEA_STARR_RESUME_CHECKPOINT")
from_dir = env("SEA_STARR_RESUME_FROM")
steps = env("SEA_STARR_RESUME_STEPS")
prefix = "sea_starr_" * lowercase(string(member))
any(isnothing, (data_dir, output_dir, checkpoint, from_dir)) &&
    error("set SEA_STARR_DATA, SEA_STARR_OUTPUT, SEA_STARR_RESUME_CHECKPOINT and SEA_STARR_RESUME_FROM")
isfile(checkpoint) || error("checkpoint not found: $checkpoint")

# Never let an overwriting writer open the original directory (or any directory that already has output).
mkpath(output_dir)
realpath(output_dir) == realpath(from_dir) && error("SEA_STARR_OUTPUT is the original run directory $from_dir")
existing = filter(f -> endswith(f, ".jld2") && !occursin("_checkpoint_iteration", f), readdir(output_dir))
isempty(existing) || error("segment directory $output_dir already contains output $(existing); refusing to overwrite")

kwargs = (; member, arch, Nx, Ny, stop_time, data_dir, output_dir, checkpoint_interval, output_prefix = prefix)
twod = env("SEA_STARR_TWOD_MINUTES")
isnothing(twod) || (kwargs = merge(kwargs, (; fields_2d_interval = parse(Float64, twod) * minutes)))
case = sea_starr(; kwargs...)
simulation = case.simulation
model = simulation.model

#####
##### Restore model state and output-writer schedules
#####

checkpoint_iteration = parse(Int, match(r"_iteration(\d+)\.jld2$", checkpoint).captures[1])
model_state, writers_state = jldopen(checkpoint, "r") do file
    load_nested_data(file["simulation/model"]), load_nested_data(file["simulation/output_writers"])
end
restore_prognostic_state!(model, model_state)
restore_prognostic_state!(simulation.output_writers, writers_state)
for writer in values(simulation.output_writers)
    reconcile_restored_output_schedule!(writer, model)
end
model.clock.iteration == checkpoint_iteration ||
    error("restored iteration $(model.clock.iteration) ≠ checkpoint iteration $checkpoint_iteration")
t₀ = model.clock.time
@info @sprintf("restored iteration %d, t = %.1f s (%.3f h) from %s", model.clock.iteration, t₀, t₀ / 3600, checkpoint)

# Δt: last value logged at or before the checkpoint iteration in the original run.
function logged_time_step(log, iteration)
    Δt = nothing
    for line in eachline(log)
        m = match(r"iter\s+(\d+), t = .*?, Δt = ([0-9.]+) (ms|μs|seconds|second|s)\b", line)
        isnothing(m) && continue
        parse(Int, m.captures[1]) ≤ iteration || break
        value = parse(Float64, m.captures[2])
        Δt = m.captures[3] == "ms" ? value / 1e3 : m.captures[3] == "μs" ? value / 1e6 : value
    end
    return Δt
end
Δt_restored = something(tryparse(Float64, something(env("SEA_STARR_RESUME_DT"), "")),
                        logged_time_step(joinpath(from_dir, "run.log"), checkpoint_iteration),
                        simulation.Δt)
simulation.Δt = Δt_restored

#####
##### Refresh the stateless callbacks against the restored state
#####

refreshable = [BreezeLab.SeaSurfaceTemperatureUpdater, BreezeLab.InversionMaskUpdater, BreezeLab.EvaporationRateUpdater]
isdefined(BreezeLab, :TrajectorySolarPosition) && push!(refreshable, BreezeLab.TrajectorySolarPosition)
update_state!(model)
refreshed = String[]
for (name, callback) in pairs(simulation.callbacks)
    if any(T -> callback.func isa T, refreshable)
        callback.func(simulation)
        push!(refreshed, string(name, ":", nameof(typeof(callback.func))))
    end
end
update_state!(model)

#####
##### Continuity check against the original run's 1-min time series
#####

mean_lwp = Field(Average(liquid_water_path(model; species = :cloud), dims = (1, 2)))
cf = cloud_fraction(model)
compute!(mean_lwp); compute!(cf)
lwp_restored = 1e3 * sum(interior(mean_lwp)); cf_restored = sum(interior(cf))
original = jldopen(joinpath(from_dir, prefix * "_timeseries.jld2"), "r") do f
    its = sort(parse.(Int, keys(f["timeseries/t"])))
    ts = [f["timeseries/t/$i"] for i in its]
    n = argmin(abs.(ts .- t₀))
    (; t = ts[n], cwp = 1e3 * f["timeseries/cwp/$(its[n])"][1], cf = f["timeseries/cloud_fraction/$(its[n])"][1])
end
@info @sprintf("continuity at t = %.1f s: restored ⟨LWP⟩ %.2f g m⁻², CF %.3f | original record t = %.1f s: ⟨LWP⟩ %.2f g m⁻², CF %.3f | Δt = %.4f s | refreshed %s",
               t₀, lwp_restored, cf_restored, original.t, original.cwp, original.cf, Δt_restored, join(refreshed, ", "))

#####
##### Provenance and run
#####

driver_sha = BreezeLab.file_sha256(@__FILE__)
code_commit = try strip(read(`git -C $(pkgdir(BreezeLab)) rev-parse HEAD`, String)) catch; "unknown" end
write_provenance(joinpath(output_dir, "provenance_resume.toml"), case;
                 extra = (; hostname = gethostname(), started = string(Dates.now()), resume_from_run = abspath(from_dir),
                            resume_checkpoint = abspath(checkpoint), resume_checkpoint_bytes = filesize(checkpoint),
                            resume_iteration = checkpoint_iteration, resume_time_s = t₀, resume_Δt_s = Δt_restored,
                            refreshed_callbacks = join(refreshed, ", "), code_commit, resume_driver_sha256 = driver_sha,
                            continuity_restored_lwp_g_m2 = lwp_restored, continuity_restored_cf = cf_restored,
                            continuity_original_t_s = original.t, continuity_original_lwp_g_m2 = original.cwp,
                            continuity_original_cf = original.cf))

if !isnothing(steps)
    simulation.stop_iteration = model.clock.iteration + parse(Int, steps)
end
wall = time()
run!(simulation)
compute!(mean_lwp); compute!(cf)
@info @sprintf("segment stopped at iteration %d, t = %.1f s (%.3f h) after %.1f min wall; ⟨LWP⟩ %.2f g m⁻², CF %.3f, all prognostic fields finite: %s",
               model.clock.iteration, model.clock.time, model.clock.time / 3600, (time() - wall) / 60,
               1e3 * sum(interior(mean_lwp)), sum(interior(cf)),
               all(f -> all(isfinite, Array(interior(f))), values(Oceananigans.prognostic_fields(model))))
