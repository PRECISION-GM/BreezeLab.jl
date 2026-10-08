# Post-run checks for a SEA STARR run directory: finite state, cloud development, aerosol
# budget, radiation/surface behaviour and time-step behaviour, from the JLD2 outputs and
# the Slurm log. Usage: julia --project analysis/check_sea_starr_run.jl <run_dir> [slurm_log]
using Oceananigans
using Oceananigans.Units
using Oceananigans.Fields: interior
using Statistics
using Printf
using JLD2

run_dir = ARGS[1]
slurm_log = length(ARGS) ≥ 2 ? ARGS[2] : nothing
prefix = joinpath(run_dir, "sea_starr_ctrl")

println("== time series (", prefix, "_timeseries.jld2)")
ts = prefix * "_timeseries.jld2"
series(name) = FieldTimeSeries(ts, name)
value(f) = [f[n][1, 1, 1] for n in eachindex(f.times)]
names = ("cwp", "rwp", "cloud_fraction", "rain_flux", "zi", "zi_max", "nᵃ_column", "nᶜˡ_column", "nʳ_column", "n_total_column")
data = Dict(name => value(series(name)) for name in names)
t = series("cwp").times
finite = all(all(isfinite, v) for v in values(data))
println("all time series finite: ", finite)
@printf("t = %.2f .. %.2f h, %d samples\n", t[1] / 3600, t[end] / 3600, length(t))
for (name, scale, unit) in (("cwp", 1e3, "g m⁻²"), ("rwp", 1e3, "g m⁻²"), ("cloud_fraction", 1, "-"), ("rain_flux", 86400, "mm d⁻¹"),
                            ("zi", 1, "m"), ("zi_max", 1, "m"))
    v = data[name] .* scale
    @printf("%-15s start %10.3f  min %10.3f  max %10.3f  end %10.3f %s\n", name, v[1], minimum(v), maximum(v), v[end], unit)
end
println("-- aerosol number budget (column, m⁻²)")
for name in ("nᵃ_column", "nᶜˡ_column", "nʳ_column", "n_total_column")
    v = data[name]
    @printf("%-15s start %.4e  end %.4e  change %+.3e\n", name, v[1], v[end], v[end] - v[1])
end
Δtotal = data["n_total_column"][end] - data["n_total_column"][1]
source = 7e5 * (t[end] - t[1])
@printf("surface source over the run: %.3e m⁻²; total change / source = %.3f (nudging, subsidence, coalescence and rain removal account for the rest)\n", source, Δtotal / source)

println("== statistics profiles (", prefix, "_statistics.jld2)")
st = prefix * "_statistics.jld2"
prof(name) = FieldTimeSeries(st, name)
column(f, n) = vec(Array(interior(f[n], 1, 1, :)))
θ = prof("θ"); qc = prof("qᶜˡ"); na = prof("nᵃ"); nc = prof("nᶜˡ"); mask = prof("nudging_mask"); cf = prof("cloud_fraction")
z = Array(znodes(θ.grid, Center()))
nlast = length(θ.times)
println("profile records: ", nlast, " (", θ.times[1] / 60, " .. ", θ.times[end] / 60, " min)")
for n in (1, nlast)
    qcol = column(qc, n); ccol = column(cf, n)
    cloudy = findall(>(1e-5), qcol)
    base = isempty(cloudy) ? NaN : z[first(cloudy)]; top = isempty(cloudy) ? NaN : z[last(cloudy)]
    kmax = isempty(cloudy) ? 0 : cloudy[argmax(qcol[cloudy])]
    @printf("record %d (%.0f min): cloud base %.0f m, top %.0f m, max qᶜˡ %.2e at %.0f m, max layer cloud fraction %.2f; nudging base %.0f m\n",
            n, θ.times[n] / 60, base, top, isempty(cloudy) ? 0 : qcol[kmax], isempty(cloudy) ? NaN : z[kmax], maximum(ccol),
            let m = column(mask, n); i = findfirst(>(0), m); isnothing(i) ? NaN : z[i] end)
    nacol = column(na, n); nccol = column(nc, n)
    @printf("   nᵃ: surface %.3e, 500 m %.3e, 2000 m %.3e kg⁻¹; nᶜˡ max %.3e kg⁻¹ at %.0f m\n",
            nacol[1], nacol[argmin(abs.(z .- 500))], nacol[argmin(abs.(z .- 2000))], maximum(nccol), z[argmax(nccol)])
end
try
    rfd = prof("radiative_flux_divergence")
    h = column(rfd, nlast)
    @printf("radiative flux divergence at the end: min %.3e, max %.3e W m⁻³ (min at %.0f m)\n", minimum(h), maximum(h), z[argmin(h)])
    lwu = prof("lw_up"); lwd = prof("lw_down"); swd = prof("sw_down"); swu = prof("sw_up")
    zf = Array(znodes(lwu.grid, Face()))
    for name in ("lw_up", "lw_down", "sw_down", "sw_up")
        f = prof(name); col = vec(Array(interior(f[nlast], 1, 1, :)))
        @printf("   %-8s surface %8.2f  top %8.2f W m⁻² (stored sign convention: positive up)\n", name, col[1], col[end])
    end
catch e
    println("radiation profiles not available: ", e)
end

if !isnothing(slurm_log) && isfile(slurm_log)
    println("== Slurm log (", slurm_log, ")")
    lines = filter(l -> occursin("iter", l) && occursin("Δt", l), readlines(slurm_log))
    println("progress messages: ", length(lines))
    isempty(lines) || (println(lines[1]); println(lines[end]))
    Δts = [parse(Float64, m.captures[1]) * (m.captures[2] == "ms" ? 1e-3 : 1.0) for l in lines for m in (match(r"Δt = ([0-9.]+) (ms|second|s)", l),) if !isnothing(m)]
    isempty(Δts) || @printf("Δt: min %.3f max %.3f s over %d messages\n", minimum(Δts), maximum(Δts), length(Δts))
    walls = [m.captures[1] for l in lines for m in (match(r"wall = ([^|]+)\|", l),) if !isnothing(m)]
    isempty(walls) || println("wall at last message: ", strip(walls[end]))
    Tr = [(parse(Float64, m.captures[1]), parse(Float64, m.captures[2])) for l in lines for m in (match(r"T ∈ \[([0-9.]+), ([0-9.]+)\]", l),) if !isnothing(m)]
    isempty(Tr) || @printf("T range over run: [%.1f, %.1f] K\n", minimum(first.(Tr)), maximum(last.(Tr)))
end
println("done")
