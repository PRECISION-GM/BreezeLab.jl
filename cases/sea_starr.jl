# # SEA STARR CTRL
#
# This example builds, runs and analyzes the SEA STARR CTRL member (Diamond et al.,
# Zenodo record 22241697): a 66-hour Lagrangian stratocumulus-to-cumulus transition in
# the southeast Atlantic under a biomass-burning smoke layer, driven by the official DEPHY
# driver. The protocol grid is 50 m horizontal (384², or the permitted 192² used here),
# 10 m vertical spacing to 2500 m then 10 % stretching to 6500 m, and a 1 s time step.
# Stage the inputs with
# `julia data_wrangling/fetch_manifest.jl cases/seastarr/inputs.toml data/seastarr_22241697`,
# then run `julia --project cases/sea_starr.jl` on an NVIDIA GPU. Settings can be
# overridden from the environment (`SEA_STARR_HOURS`, `SEA_STARR_NX`, `SEA_STARR_OUTPUT`,
# `SEA_STARR_DATA`, `SEA_STARR_MEMBER`, `SEA_STARR_CHECKPOINT_HOURS`, `SEA_STARR_TWOD_MINUTES`, `SEA_STARR_PICKUP`) so
# the same script serves the short pilot, the full integration and a restart from the latest
# checkpoint in the output directory.
#
# The physics that the pinned Breeze cannot represent is listed in
# `case.config.departures` and written to the provenance file: no aerosol optical
# properties, a radiation column ending at the LES top, a fixed solar coordinate, and the
# labelled nudging-ramp, activation and regeneration choices. A run with these departures
# is an exploratory CTRL, not a protocol-complete one.

using BreezeLab
using Oceananigans, Oceananigans.Units
using Oceananigans.Fields: interior
using CairoMakie
using Statistics
import Dates

# ## Build the simulation
#
# The exported constructor reads the driver, builds the protocol grid, the P3 microphysics
# with the κ-Köhler aerosol mode and prognostic reservoir, the subsidence/geostrophic
# forcing, the inversion-following nudging, the bulk surface fluxes from the trajectory SST,
# RRTMGP radiation and the output writers, and returns the case without advancing it.

member = Symbol(get(ENV, "SEA_STARR_MEMBER", "CTRL"))
arch = GPU()
Nx = Ny = parse(Int, get(ENV, "SEA_STARR_NX", "192"))
stop_time = parse(Float64, get(ENV, "SEA_STARR_HOURS", "66")) * hours
data_dir = get(ENV, "SEA_STARR_DATA", joinpath(pkgdir(BreezeLab), "data", "seastarr_22241697"))
output_dir = get(ENV, "SEA_STARR_OUTPUT", joinpath(pkgdir(BreezeLab), "output", "sea_starr_$(lowercase(string(member)))"))
checkpoint_interval = parse(Float64, get(ENV, "SEA_STARR_CHECKPOINT_HOURS", "3")) * hours
fields_2d_interval = parse(Float64, get(ENV, "SEA_STARR_TWOD_MINUTES", "3")) * minutes   # animation cadence (plan views, lowest level, x–z slices)

case = sea_starr(; member, arch, Nx, Ny, stop_time, data_dir, output_dir, checkpoint_interval, fields_2d_interval,
                   output_prefix = "sea_starr_$(lowercase(string(member)))")
simulation = case.simulation;

# ## Run
#
# Provenance (driver path and checksum, configuration, departures, software revisions) is
# written next to the output before the run starts; `case.config.output_storage_estimate_GB`
# records the storage budget of the writers (at 192² × 288 and the 3-min animation cadence:
# ≈ 14 GB hourly 3D, ≈ 3.4 GB 2-D/slices, 2.9 GB for the retained checkpoint).

pickup = get(ENV, "SEA_STARR_PICKUP", "false") == "true"
mkpath(output_dir)
write_provenance(joinpath(output_dir, pickup ? "provenance_restart_$(Dates.format(Dates.now(), "yyyymmddHHMM")).toml" : "provenance.toml"), case;
                 extra = (; hostname = gethostname(), started = string(Dates.now()), pickup))
# A restart resumes the prognostic state from the latest checkpoint; the forcing position
# follows the restored clock and the radiation recomputes at its next scheduled interval.
run!(simulation; pickup)

# ## Analyze
#
# Read the saved domain means: cloud and rain water paths, surface precipitation,
# inversion height, and the aerosol number budget (interstitial, in droplets, in rain, and
# their column total, which changes only through the surface source, nudging, subsidence,
# collisions and rain removal).

prefix = joinpath(output_dir, "sea_starr_$(lowercase(string(member)))")
ts = prefix * "_timeseries.jld2"
series(name) = FieldTimeSeries(ts, name)
cwp = series("cwp"); rwp = series("rwp"); rain = series("rain_flux"); zi = series("zi"); zimax = series("zi_max")
na = series("nᵃ_column"); nc = series("nᶜˡ_column"); nr = series("nʳ_column"); nt = series("n_total_column")
t = cwp.times ./ 3600
value(f) = [f[n][1, 1, 1] for n in eachindex(f.times)]

fig = Figure(size = (1000, 1100))
ax1 = Axis(fig[1, 1], xlabel = "Hours since 2017-08-15 21 UTC", ylabel = "Water path (g m⁻²)")
lines!(ax1, t, 1e3 .* value(cwp), label = "cloud")
lines!(ax1, t, 1e3 .* value(rwp), label = "rain")
axislegend(ax1, position = :lt)
ax2 = Axis(fig[2, 1], xlabel = "Hours", ylabel = "Surface rain (mm day⁻¹)")
lines!(ax2, t, 86400 .* value(rain))
ax3 = Axis(fig[3, 1], xlabel = "Hours", ylabel = "Inversion height (m)")
lines!(ax3, t, value(zi), label = "mean")
lines!(ax3, t, value(zimax), label = "max (nudging base − 100 m)")
axislegend(ax3, position = :lt)
ax4 = Axis(fig[4, 1], xlabel = "Hours", ylabel = "Column number (10¹⁰ m⁻²)")
lines!(ax4, t, 1e-10 .* value(na), label = "aerosol nᵃ")
lines!(ax4, t, 1e-10 .* value(nc), label = "droplets nᶜˡ")
lines!(ax4, t, 1e-10 .* value(nr), label = "rain nʳ")
lines!(ax4, t, 1e-10 .* value(nt), label = "total", linestyle = :dash)
axislegend(ax4, position = :rt)
save(joinpath(output_dir, "sea_starr_timeseries.png"), fig)

# Profiles of the 15-minute statistics at the start and end of the run: θₗ, moisture, cloud
# water, droplet and aerosol number, and the nudging mask.

st = prefix * "_statistics.jld2"
prof(name) = FieldTimeSeries(st, name)
θp = prof("θ"); qv = prof("qᵛ"); ql = prof("qᶜˡ"); ncp = prof("nᶜˡ"); nap = prof("nᵃ"); mask = prof("nudging_mask")
z = znodes(θp.grid, Center())
n_last = length(θp.times)
column(f, n) = vec(Array(interior(f[n], 1, 1, :)))

fig2 = Figure(size = (1200, 500))
axθ = Axis(fig2[1, 1], xlabel = "θₗ (K)", ylabel = "z (m)")
lines!(axθ, column(θp, 1), z, label = "start"); lines!(axθ, column(θp, n_last), z, label = "end"); axislegend(axθ, position = :rb)
axq = Axis(fig2[1, 2], xlabel = "qᵛ, 10 qᶜˡ (g kg⁻¹)")
lines!(axq, 1e3 .* column(qv, n_last), z, label = "qᵛ"); lines!(axq, 1e4 .* column(ql, n_last), z, label = "10 qᶜˡ"); axislegend(axq, position = :rt)
axn = Axis(fig2[1, 3], xlabel = "number (mg⁻¹)")
lines!(axn, 1e-6 .* column(nap, 1), z, label = "nᵃ start"); lines!(axn, 1e-6 .* column(nap, n_last), z, label = "nᵃ end")
lines!(axn, 1e-6 .* column(ncp, n_last), z, label = "nᶜˡ end"); axislegend(axn, position = :rt)
axm = Axis(fig2[1, 4], xlabel = "nudging mask")
lines!(axm, column(mask, n_last), z)
for ax in (axθ, axq, axn, axm)
    ylims!(ax, 0, 3000)
end
save(joinpath(output_dir, "sea_starr_profiles.png"), fig2)
#md cp(joinpath(output_dir, "sea_starr_timeseries.png"), "sea_starr_timeseries.png"; force=true);
fig #src

#md # ![SEA STARR time series](sea_starr_timeseries.png)
