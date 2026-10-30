# Animations of a TRACER–DP-SCREAM-matched periodic run (cases/tracer_dp_scream.jl output):
#
#   1. convection.mp4   — the saved slices: cloud liquid and w at the slice height (xy), the cloud-liquid /
#                          rain / w section at y = Ly/2 (xz), with the domain-mean LWP, IWP and surface
#                          precipitation traces below and a cursor at the frame time (UTC and CDT = UTC−5).
#   2. time_height.mp4  — the LES time–height total (liquid + ice) cloud-fraction section building up over
#                          the run, alongside the archived DP-SCREAM 3 km (August run) and 0.5 km (5–15 August
#                          cold start) TOT_CLOUD_FRAC domain means for the same UTC window, and the
#                          precipitation traces of all three.
#   3. convection_montage.png — four frames of (1) at the times of the four largest domain-mean LWP maxima of
#                          distinct days (fallback for viewers without video), and time_height.png (last frame of 2).
#
#   julia --project analysis/animate_tracer_dp_scream.jl <run dir> [<output dir> = <run dir>/animations]
#
# Playback standard: cloud-scale (episode) movies play at most 15 simulated minutes per second of video at
# 10–15 fps (1-min slices at 12 fps = 12 min/s); multi-day overviews play about 1 simulated hour per second.
# The frame stride and frame rate follow from the saved frame spacing and the target speed:
#   stride = ⌈speed / (max_fps · spacing)⌉,  fps = clamp(speed / (stride · spacing), 2, max_fps)
# with speed = ANIMATION_PLAYBACK_MINUTES_PER_SECOND (default 12 for windows up to 2 days, 60 otherwise).
# Frames coarser than the target allows (e.g. 30-min slices for an episode) play at 2 fps and a warning is printed.
#
# Environment: ANIMATION_PLAYBACK_MINUTES_PER_SECOND, ANIMATION_MAX_FPS (15), ANIMATION_FPS and
# ANIMATION_FRAME_STRIDE (explicit overrides of the computed values), ANIMATION_WINDOW_HOURS ("a:b", hours since
# the start), ANIMATION_TAG (file-name suffix), ANIMATION_PARTS (convection,montage,time_height),
# ANIMATION_RUN_LABEL, ANIMATION_PX_PER_UNIT (2), ANIMATION_WIDTH (1000 points), ANIMATION_COMPRESSION (CRF, 28),
# TIME_HEIGHT_STRIDE (profile records per time–height frame; default from the playback rule), DP_SCREAM_DATA
# (archive directory, default data/tracer_dp_scream of the package), TRACER_START (2022-08-01T00:00:00).
# The run's JLD2 files are read as they are; snapshot copies are recommended while the run is still writing.

using BreezeLab
using CairoMakie, JLD2, Oceananigans, Oceananigans.Units, Printf, Statistics
using Dates: Dates, DateTime, Second, Hour, Day
using TOML: TOML

run_dir = ARGS[1]
out_dir = length(ARGS) ≥ 2 ? ARGS[2] : joinpath(run_dir, "animations")
mkpath(out_dir)
max_fps = parse(Int, get(ENV, "ANIMATION_MAX_FPS", "15"))
min_fps = 2
"""
    playback_plan(spacing, span; stride=nothing, fps=nothing)

Frame stride and frame rate for frames `spacing` seconds apart over a window `span` seconds long, following the
playback standard (see the header). Returns `(; stride, fps, speed)` with `speed` the actual simulated seconds
per second of video.
"""
function playback_plan(spacing, span; stride = nothing, fps = nothing)
    target = 60 * parse(Float64, get(ENV, "ANIMATION_PLAYBACK_MINUTES_PER_SECOND", span ≤ 2 * 86400 ? "12" : "60"))
    s = something(stride, max(1, ceil(Int, target / (max_fps * spacing))))
    f = something(fps, clamp(round(Int, target / (s * spacing)), min_fps, max_fps))
    speed = f * s * spacing
    speed > 1.25 * target && @warn @sprintf("frames every %.0f min cannot meet %.0f sim-min/s; playing at %.0f sim-min/s", spacing / 60, target / 60, speed / 60)
    return (; stride = s, fps = f, speed)
end
env_int(name) = haskey(ENV, name) ? parse(Int, ENV[name]) : nothing
width = parse(Int, get(ENV, "ANIMATION_WIDTH", "1000"))
compression = parse(Int, get(ENV, "ANIMATION_COMPRESSION", "28"))
data_dir = get(ENV, "DP_SCREAM_DATA", joinpath(pkgdir(BreezeLab), "data", "tracer_dp_scream"))
start = DateTime(get(ENV, "TRACER_START", "2022-08-01T00:00:00"))
# A window (hours since start, "a:b"), a frame stride, a file tag and the parts to render let one run give
# both a full-length movie and a short smooth one from fine-cadence slices without reading every frame twice.
window = let w = get(ENV, "ANIMATION_WINDOW_HOURS", ""); isempty(w) ? (-Inf, Inf) : Tuple(parse.(Float64, split(w, ":"))) end
tag = get(ENV, "ANIMATION_TAG", "")
parts = split(get(ENV, "ANIMATION_PARTS", "convection,montage,time_height"), ",")
run_label = get(ENV, "ANIMATION_RUN_LABEL", "")
px_per_unit = parse(Float64, get(ENV, "ANIMATION_PX_PER_UNIT", "2"))   # 2: 2000-px frames for the 1000-point figure; 1 for web versions
tagged(name) = isempty(tag) ? name : replace(name, r"(\.[a-z0-9]+)$" => SubstitutionString("_" * tag * "\\1"))
utc_offset = -5                                    # Houston, CDT

prefix = only(unique(first.(split.(filter(f -> endswith(f, "_slices.jld2"), readdir(run_dir)), "_slices"))))
slice_file = joinpath(run_dir, prefix * "_slices.jld2")
series_file = joinpath(run_dir, prefix * "_timeseries.jld2")
profile_file = joinpath(run_dir, prefix * "_profiles.jld2")

stamp(t) = (utc = start + Second(round(Int, t)); local_t = utc + Hour(utc_offset);
            @sprintf("%s UTC (%s CDT), day %.2f", Dates.format(utc, "yyyy-mm-dd HH:MM"), Dates.format(local_t, "HH:MM"), t / 86400))
slab(f, dims) = dropdims(Array(interior(f)); dims)

#####
##### Time series (60 s) and profiles (30-min means)
#####

series(name) = (fts = FieldTimeSeries(series_file, name); (fts.times, [fts[n][1, 1, 1] for n in eachindex(fts.times)]))
t_s, lwp_s = series("lwp")
_, iwp_s = series("iwp")
_, rain_s = series("rain_flux")
_, ice_s = series("ice_flux")
precip_s = 86400 .* (rain_s .+ ice_s)                      # mm day⁻¹
hours_s = t_s ./ 3600

#####
##### 1. Slice animation
#####

# Slices are read frame by frame from disk (fine-cadence runs hold thousands of 512² frames)
slices(name) = FieldTimeSeries(slice_file, name; backend = OnDisk())
qxy = slices("qᶜˡ_xy")
lxy = slices("lwp")
wxy = slices("w_xy")
qxz = slices("qᶜˡ_xz")
rxz = slices("qʳ_xz")
wxz = slices("w_xz")
times = qxy.times
Nt = length(times)
window_frames = [n for n in 1:Nt if window[1] ≤ times[n] / 3600 ≤ window[2]]
plan = playback_plan(times[2] - times[1], times[last(window_frames)] - times[first(window_frames)];
                     stride = env_int("ANIMATION_FRAME_STRIDE"), fps = env_int("ANIMATION_FPS"))
frames = window_frames[1:plan.stride:end]
framerate = plan.fps
grid = qxy.grid
x = collect(xnodes(grid, Center())) ./ 1e3
y = collect(ynodes(grid, Center())) ./ 1e3
z = collect(znodes(grid, Center())) ./ 1e3
zf = collect(znodes(grid, Face())) ./ 1e3
prov = joinpath(run_dir, "provenance.toml")
config = isfile(prov) ? get(TOML.parsefile(prov), "config", Dict()) : Dict()
slice_height = Float64(get(config, "slice_height", 1500.0))
Δx_run = haskey(config, "Lx") && haskey(config, "Nx") ? config["Lx"] / config["Nx"] : x[2] - x[1]
model_label = string("Breeze LES ", round(Int, haskey(config, "Lx") ? Δx_run : 1000Δx_run), " m", isempty(run_label) ? "" : " ($run_label)")

q_xy(n) = 1e3 .* slab(qxy[n], 3)                 # g kg⁻¹ at the slice height
l_xy(n) = 1e3 .* slab(lxy[n], 3)                 # g m⁻², cloud liquid water path
w_xy(n) = slab(wxy[n], 3)                       # m s⁻¹
q_xz(n) = 1e3 .* slab(qxz[n], 2)                 # g kg⁻¹, (x, z)
r_xz(n) = 1e3 .* slab(rxz[n], 2)                 # g kg⁻¹, (x, z)
w_xz(n) = slab(wxz[n], 2)                       # m s⁻¹, (x, z faces)

sample = 1:max(1, Nt ÷ 40):Nt
qmax = max(0.2, quantile(vcat([vec(q_xz(n)) for n in sample]...), 0.999))
rmax = max(0.05, quantile(vcat([vec(r_xz(n)) for n in sample]...), 0.999))
wmax = 6.0
lmax = max(100.0, quantile(vcat([vec(l_xy(n)) for n in sample]...), 0.995))

function slice_figure()
    n = Observable(1)
    fig = Figure(size = (width, round(Int, 1.05width)), fontsize = 11)
    Label(fig[0, 1:2], @lift(string(model_label, ", TRACER–DP-SCREAM forcing — ", stamp(times[$n]))), fontsize = 15, tellwidth = false)
    ax1 = Axis(fig[1, 1], title = "cloud liquid water path (g m⁻²)", xlabel = "x (km)", ylabel = "y (km)", aspect = DataAspect())
    hm1 = heatmap!(ax1, x, y, @lift(l_xy($n)), colormap = :Blues, colorrange = (0, lmax))
    Colorbar(fig[1, 1][1, 2], hm1, width = 8)
    ax2 = Axis(fig[1, 2], title = @sprintf("w (m s⁻¹) at %.0f m", slice_height), xlabel = "x (km)", ylabel = "y (km)", aspect = DataAspect())
    hm2 = heatmap!(ax2, x, y, @lift(w_xy($n)), colormap = :balance, colorrange = (-wmax, wmax))
    Colorbar(fig[1, 2][1, 2], hm2, width = 8)
    ax3 = Axis(fig[2, 1:2], title = "section at y = Ly/2: qᶜˡ (blue shading, g kg⁻¹), rain qʳ (grey contours at 0.01, 0.1, 0.5 g kg⁻¹), w = ±2 m s⁻¹ (red/purple)",
               xlabel = "x (km)", ylabel = "z (km)")
    hm3 = heatmap!(ax3, x, z, @lift(q_xz($n)), colormap = :Blues, colorrange = (0, qmax))
    contour!(ax3, x, z, @lift(r_xz($n)), levels = [0.01, 0.1, 0.5], color = (:gray20, 0.9), linewidth = 0.8)
    contour!(ax3, x, zf, @lift(w_xz($n)), levels = [-2.0], color = (:purple, 0.8), linewidth = 0.8)
    contour!(ax3, x, zf, @lift(w_xz($n)), levels = [2.0], color = (:red, 0.8), linewidth = 0.8)
    ylims!(ax3, 0, 18)
    Colorbar(fig[2, 3], hm3, width = 8, label = "qᶜˡ (g kg⁻¹)")
    ax4 = Axis(fig[3, 1:2], xlabel = "hours since $(Dates.format(start, "yyyy-mm-dd HH:MM")) UTC", ylabel = "water path (g m⁻²)")
    lines!(ax4, hours_s, 1e3 .* lwp_s, color = :steelblue, label = "LWP")
    lines!(ax4, hours_s, 1e3 .* iwp_s, color = :slategray, label = "IWP")
    axislegend(ax4, position = :lt, framevisible = false)
    ax5 = Axis(fig[3, 1:2], ylabel = "precipitation (mm day⁻¹)", yaxisposition = :right)
    hidespines!(ax5); hidexdecorations!(ax5)
    lines!(ax5, hours_s, precip_s, color = :darkorange, label = "surface rain + ice")
    linkxaxes!(ax4, ax5)
    vlines!(ax4, @lift([times[$n] / 3600]), color = :black, linewidth = 1.5)
    xlims!(ax4, max(0, times[first(frames)] / 3600 - 1), min(max(times[end], t_s[end]), times[last(frames)] + 3600) / 3600)
    rowsize!(fig.layout, 1, Relative(0.42)); rowsize!(fig.layout, 2, Relative(0.3)); rowsize!(fig.layout, 3, Relative(0.2))
    return fig, n
end

if "convection" in parts
    fig1, n1 = slice_figure()
    n1[] = first(frames)
    convection_path = joinpath(out_dir, tagged("convection.mp4"))
    record(fig1, convection_path, frames; framerate, compression, px_per_unit) do i
        n1[] = i
    end
    @info @sprintf("wrote %s (%.1f MB, %d frames every %.0f min at %d fps = %.1f sim-min per second)", convection_path,
                   filesize(convection_path) / 1e6, length(frames), plan.stride * (times[2] - times[1]) / 60, framerate, plan.speed / 60)
end

picks = Int[]
if "montage" in parts
    # Montage: the four strongest cloud-liquid frames on distinct days
    lwp_at(n) = (k = searchsortedlast(t_s, times[n]); k ≥ 1 ? lwp_s[k] : 0.0)
    ranked = sort(1:Nt; by = n -> -lwp_at(n))
    picks = Int[]
    for n in ranked
        any(abs(times[n] - times[m]) < 86400 for m in picks) && continue
        push!(picks, n)
        length(picks) == 4 && break
    end
    sort!(picks)
    figm = Figure(size = (2width, round(Int, 2.1width)), fontsize = 11)
    for (k, n) in enumerate(picks)
        sub = figm[(k - 1) ÷ 2 + 1, (k - 1) % 2 + 1] = GridLayout()
        Label(sub[0, 1:2], stamp(times[n]), fontsize = 14, tellwidth = false)
        a1 = Axis(sub[1, 1], title = "cloud liquid water path (g m⁻²)", aspect = DataAspect(), xlabel = "x (km)", ylabel = "y (km)")
        heatmap!(a1, x, y, l_xy(n), colormap = :Blues, colorrange = (0, lmax))
        a2 = Axis(sub[1, 2], title = @sprintf("w (m s⁻¹) at %.0f m", slice_height), aspect = DataAspect(), xlabel = "x (km)")
        heatmap!(a2, x, y, w_xy(n), colormap = :balance, colorrange = (-wmax, wmax))
        a3 = Axis(sub[2, 1:2], title = "qᶜˡ shading, qʳ contours, w = ±2 m s⁻¹ at y = Ly/2", xlabel = "x (km)", ylabel = "z (km)")
        heatmap!(a3, x, z, q_xz(n), colormap = :Blues, colorrange = (0, qmax))
        contour!(a3, x, z, r_xz(n), levels = [0.01, 0.1, 0.5], color = (:gray20, 0.9), linewidth = 0.8)
        contour!(a3, x, zf, w_xz(n), levels = [-2.0], color = (:purple, 0.8), linewidth = 0.8)
        contour!(a3, x, zf, w_xz(n), levels = [2.0], color = (:red, 0.8), linewidth = 0.8)
        ylims!(a3, 0, 18)
    end
    montage_path = joinpath(out_dir, "convection_montage.png")
    save(montage_path, figm; px_per_unit = 1)
    @info "wrote $montage_path (frames at $(join([stamp(times[n]) for n in picks], "; ")))"
end

if "time_height" in parts
#####
##### 2. Time–height cloud fraction vs the archived DP-SCREAM domain means
#####

cf = FieldTimeSeries(profile_file, "total_cloud_fraction")
tp = cf.times
Np = length(tp)
cf_matrix = [cf[n][1, 1, k] for k in eachindex(z), n in 1:Np]      # (z, t)
stop_utc = start + Second(round(Int, max(tp[end], times[end])))
window_hours = Dates.value(stop_utc - start) / 3.6e6

references = Dict{String, Any}()
for (key, file) in (("DP-SCREAM 3 km (August run, SHOC, 200 km domain mean)", "DP-SCREAM 3km_August.nc"),
                    ("DP-SCREAM 0.5 km (5–15 August cold start, 200 km domain mean)", "DP-SCREAM 0.5km_August_5_to_August_15_run.nc"))
    path = joinpath(data_dir, file)
    if isfile(path)
        ref = read_dp_scream_output(path)
        sel = dp_scream_window(ref, start, start + Day(14))
        th = [Dates.value(ref.time[n] - start) / 3.6e6 for n in sel]
        zref = vec(mean(ref.Z3[:, sel], dims = 2)) ./ 1e3
        order = sortperm(zref)
        references[key] = (; th, zref = zref[order], cf = ref.TOT_CLOUD_FRAC[order, sel], precip = ref.PRECL[sel])
    else
        @warn "reference file not found: $path"
    end
end

plan_th = playback_plan(tp[2] - tp[1], tp[end] - tp[1]; stride = env_int("TIME_HEIGHT_STRIDE"))
stride = plan_th.stride
frames2 = collect(stride:stride:Np)
isempty(frames2) && push!(frames2, Np)
last(frames2) == Np || push!(frames2, Np)
m = Observable(frames2[1])
masked_les = @lift begin
    M = fill(NaN, size(cf_matrix)); M[:, 1:$m] .= cf_matrix[:, 1:$m]; permutedims(M)
end
masked_ref(ref) = @lift begin
    tcur = tp[$m] / 3600
    M = copy(ref.cf); M[:, ref.th .> tcur] .= NaN; permutedims(M)
end

fig2 = Figure(size = (width, round(Int, 1.1width)), fontsize = 11)
Label(fig2[0, 1:2], @lift(string("Total cloud fraction, time–height — through ", stamp(tp[$m]))), fontsize = 15, tellwidth = false)
axA = Axis(fig2[1, 1], title = "Breeze LES 200 m (51.2 km domain mean, cells with qᶜˡ + qⁱ > 10⁻⁵ kg kg⁻¹)", ylabel = "z (km)")
hmA = heatmap!(axA, tp ./ 3600, z, masked_les, colormap = :Blues, colorrange = (0, 1))
ylims!(axA, 0, 18); xlims!(axA, 0, window_hours)
Colorbar(fig2[1, 2], hmA, width = 8)
row = 2
for (key, ref) in sort(collect(references); by = first)
    global row
    ax = Axis(fig2[row, 1], title = key, ylabel = "z (km)")
    heatmap!(ax, ref.th, ref.zref, masked_ref(ref), colormap = :Blues, colorrange = (0, 1))
    ylims!(ax, 0, 18); xlims!(ax, 0, window_hours)
    row += 1
end
axP = Axis(fig2[row, 1], xlabel = "hours since $(Dates.format(start, "yyyy-mm-dd HH:MM")) UTC", ylabel = "precipitation (mm day⁻¹)")
lines!(axP, hours_s, precip_s, color = :darkorange, label = "LES rain + ice")
for (key, ref) in sort(collect(references); by = first)
    lines!(axP, ref.th, ref.precip, label = first(split(key, " (")), linewidth = 1)
end
axislegend(axP, position = :lt, framevisible = false)
vlines!(axP, @lift([tp[$m] / 3600]), color = :black, linewidth = 1.5)
xlims!(axP, 0, window_hours)
time_height_path = joinpath(out_dir, "time_height.mp4")
record(fig2, time_height_path, frames2; framerate = plan_th.fps, compression) do i
    m[] = i
end
m[] = Np
save(joinpath(out_dir, "time_height.png"), fig2; px_per_unit = 1)
@info @sprintf("wrote %s (%.1f MB, %d frames at %d fps = %.1f sim-min per second)", time_height_path,
               filesize(time_height_path) / 1e6, length(frames2), plan_th.fps, plan_th.speed / 60)

#####
##### README
#####

commit = try strip(read(`git -C $(pkgdir(BreezeLab)) rev-parse HEAD`, String)) catch; "unknown" end
open(joinpath(out_dir, "README.md"), "w") do io
    println(io, "# TRACER–DP-SCREAM baseline animations\n")
    println(io, "Run: `$(abspath(run_dir))` ($(model_label), $(get(config, "Nx", "?"))² columns, 160 levels to 22 km; forcing = DP-SCREAM TRACER IOP file; start $(start) UTC; CDT = UTC−5). Rendered by `analysis/animate_tracer_dp_scream.jl` at BreezeLab commit $commit on $(Dates.now(Dates.UTC)) UTC from outputs reaching t = $(round(times[end] / 86400, digits=2)) days.\n")
    println(io, "| File | Content | Cadence |")
    println(io, "| --- | --- | --- |")
    println(io, "| `convection.mp4` | xy: cloud liquid water path (g m⁻², colour range 0–$(round(Int, lmax))) and vertical velocity w (m s⁻¹, ±$wmax) at $(round(Int, slice_height)) m (the saved qᶜˡ slice at that height is mostly below the ~3 km cloud base and is not shown); xz at y = Ly/2: qᶜˡ shading, rain qʳ contours (0.01/0.1/0.5 g kg⁻¹), w = ±2 m s⁻¹ contours; below: domain-mean LWP and IWP (g m⁻²) and surface rain + ice flux (mm day⁻¹) with a cursor | frames every $(round(Int, plan.stride * (times[2] - times[1]) / 60)) min at $framerate fps ($(round(plan.speed / 60, digits=1)) simulated min per second), $(length(frames)) frames |")
    println(io, "| `time_height.mp4` / `time_height.png` | total (liquid + ice) cloud fraction, LES 30-min horizontal means vs archived DP-SCREAM 3 km (August run) and 0.5 km (5–15 Aug cold start) TOT_CLOUD_FRAC on their mean Z3 heights, same UTC window; bottom: precipitation (LES rain + ice; DP-SCREAM PRECL, mm day⁻¹ inferred) | one frame per $(stride) profile records ($(stride * round(Int, (tp[2] - tp[1]) / 60)) min) at $(plan_th.fps) fps ($(round(plan_th.speed / 60, digits=1)) simulated min per second) |")
    println(io, "| `convection_montage.png` | four `convection.mp4` frames at the largest domain-mean LWP of four distinct days: $(join([stamp(times[n]) for n in picks], "; ")) | — |")
    println(io, "\nUnits: mass fractions in g kg⁻¹, water paths in g m⁻², fluxes converted to mm day⁻¹ of liquid water, heights in km above the surface. The DP-SCREAM archive time labels are local (CDT) and were shifted +5 h to UTC. The LES is an intentional physics/resolution departure from DP-SCREAM (see docs/cases/tracer_dp_scream.md); agreement is not implied.")
end

end
@info "done"
