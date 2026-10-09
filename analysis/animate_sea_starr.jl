# SEA STARR animations from the saved JLD2 output of the 66 h runs (CairoMakie `record`,
# H.264 MP4, ≤ 1000 px wide), plus a four-frame PNG montage of each:
#
#   1. plan views of the liquid water path (cloud + rain) and the surface rain rate,
#      CTRL / N100 / N030 side by side at the same times (15 min cadence; UTC and the
#      trajectory's local solar time in the title),
#   2. an x–z cross-section through the CTRL domain: cloud water, rain water and the
#      interstitial aerosol number (hourly 3D output),
#   3. a time–height view of the horizontal-mean aerosol and droplet number (15 min statistics)
#      with the mean inversion height, revealed as time advances.
#
# The script only reads the run directories (the GPU jobs keep writing to them). It is
# re-runnable as the runs progress: it uses the records present when it starts and, for the
# three-member comparison, the times common to all three.
#
#   julia --project analysis/animate_sea_starr.jl \
#       --ctrl runs/ctrl_full_66h_job211 --n100 runs/n100_full_66h_job214 --n030 runs/n030_full_66h_job218 \
#       --copy-to /shared/home/greg/breezelab-work/campaign-figures/animations [--max-frames N] [--framerate 12]
#
# Outputs go to <run>/animations/ (the comparison under the CTRL run) and are copied to
# --copy-to with a README listing fields, units, cadence and the source commit.

using Oceananigans
using Oceananigans.Units
using Oceananigans.Fields: interior
using Oceananigans.OutputReaders: OnDisk
using CairoMakie
using NCDatasets
using Statistics
using Printf
using Dates
using JLD2: jldopen

#####
##### Arguments
#####

function parse_args(args)
    opts = Dict{String, String}()
    i = 1
    while i ≤ length(args)
        startswith(args[i], "--") && i < length(args) || error("unexpected argument $(args[i])")
        opts[args[i][3:end]] = args[i + 1]
        i += 2
    end
    return opts
end

opts = parse_args(ARGS)
ctrl_dir = get(opts, "ctrl", "runs/ctrl_full_66h_job211")
n100_dir = get(opts, "n100", "runs/n100_full_66h_job214")
n030_dir = get(opts, "n030", "runs/n030_full_66h_job218")
copy_to = get(opts, "copy-to", "/shared/home/greg/breezelab-work/campaign-figures/animations")
max_frames = parse(Int, get(opts, "max-frames", "100000"))
framerate = parse(Int, get(opts, "framerate", "12"))
trajectories = get(opts, "trajectories", "/shared/home/greg/breezelab-runs/20261006/inputs/mip_sources/seastarr_22241697/SEA_STARR_Raw_Trajectories.nc")
compression = parse(Int, get(opts, "compression", "24"))
members = [(:CTRL, ctrl_dir, "sea_starr_ctrl"), (:N100, n100_dir, "sea_starr_n100"), (:N030, n030_dir, "sea_starr_n030")]
commit = try strip(read(`git rev-parse --short HEAD`, String)) catch; "unknown" end

epoch = DateTime(2017, 8, 15, 21)     # driver start (seconds since 2017-08-15 21:00:00 UTC)
width = 1000                          # px, the maximum requested width
CairoMakie.activate!(type = "png", px_per_unit = 1)

#####
##### Local solar time along the composite trajectory
#####

# Hourly mean longitude of the 40 GEOS-5 trajectories of the composite (missing values
# ignored); the driver itself carries a single mid-period position (−8.34° E).
function trajectory_longitude(path)
    isfile(path) || return (t -> -8.336)
    ds = NCDataset(path)
    lon = Float64.(coalesce.(ds["lon"][:, :], NaN))
    t = Float64.(ds["t"][:]) .* 3600
    close(ds)
    mean_lon = [mean(filter(isfinite, lon[n, :])) for n in axes(lon, 1)]
    return s -> (s ≤ t[1] ? mean_lon[1] : s ≥ t[end] ? mean_lon[end] :
                 (k = searchsortedlast(t, s); w = (s - t[k]) / (t[k+1] - t[k]); (1 - w) * mean_lon[k] + w * mean_lon[k+1]))
end
longitude_at = trajectory_longitude(trajectories)

function time_label(t)
    utc = epoch + Second(round(Int, t))
    lon = longitude_at(t)
    lst = utc + Second(round(Int, lon / 15 * 3600))
    return @sprintf("%s UTC  |  local solar time %s (lon %.1f°E)  |  t = %.2f h", Dates.format(utc, "yyyy-mm-dd HH:MM"), Dates.format(lst, "HH:MM"), lon, t / 3600)
end

#####
##### Helpers
#####

mkpath_for(path) = (mkpath(dirname(path)); path)

function montage(frames, path; frame_titles=nothing)
    # frames: Vector of RGB images (from `colorbuffer`), laid out 2 × 2 at half size
    n = length(frames)
    h, w = size(frames[1])
    fig = Figure(size = (width, round(Int, width * h / w)), figure_padding = 0)
    for (k, img) in enumerate(frames)
        r, c = fldmod1(k, 2)
        ax = Axis(fig[r, c], aspect = DataAspect())
        hidedecorations!(ax); hidespines!(ax)
        image!(ax, rotr90(img))
    end
    rowgap!(fig.layout, 2); colgap!(fig.layout, 2)
    save(path, fig)
    return path
end

function record_with_montage(fig, update!, frames, mp4_path, png_path; montage_frames)
    # H.264 "high" profile with 4:2:0 pixels plays in browsers (Makie's default is high422)
    record(fig, mp4_path, frames; framerate, compression, profile = "high", pixel_format = "yuv420p") do n
        update!(n)
    end
    # `colorbuffer` returned the last recorded frame after `record`; saving re-renders the scene
    images = map(montage_frames) do n
        update!(n)
        tmp = tempname() * ".png"
        save(tmp, fig)
        img = CairoMakie.Makie.FileIO.load(tmp)
        rm(tmp)
        img
    end
    montage(images, png_path)
    return nothing
end

montage_indices(frames) = unique(round.(Int, range(first(frames), last(frames); length = 4)))

fmt_mb(path) = @sprintf("%.1f MB", filesize(path) / 1e6)

#####
##### 1. Plan views: LWP and surface rain, CTRL / N100 / N030
#####

function plan_view_animation()
    series = map(members) do (name, dir, prefix)
        file = joinpath(dir, prefix * "_2d.jld2")
        (; name, dir,
           cwp = FieldTimeSeries(file, "cwp"; backend = OnDisk()),
           rwp = FieldTimeSeries(file, "rwp"; backend = OnDisk()),
           rain = FieldTimeSeries(file, "rain"; backend = OnDisk()))
    end
    # times common to the three runs (they are at different stages of the 66 h integration)
    common = reduce(intersect, [round.(s.cwp.times) for s in series])
    isempty(common) && error("no common output times across the members")
    index(s, t) = findfirst(==(t), round.(s.cwp.times))
    frames = 1:min(length(common), max_frames)
    grid = series[1].cwp.grid
    x = xnodes(grid, Center()) ./ 1e3
    y = ynodes(grid, Center()) ./ 1e3
    # fixed colour ranges from a subsample of frames, shared by the three members
    sample = frames[1:max(1, length(frames) ÷ 24):end]
    lwp_sample = Float32[]; rain_sample = Float32[]
    for n in sample, s in series
        i = index(s, common[n])
        append!(lwp_sample, vec(1e3 .* (interior(s.cwp[i], :, :, 1) .+ interior(s.rwp[i], :, :, 1))))
        append!(rain_sample, vec(86400 .* interior(s.rain[i], :, :, 1)))
    end
    lwp_max = max(50f0, quantile(lwp_sample, 0.995))
    rain_max = max(0.1f0, quantile(rain_sample, 0.995))

    fig = Figure(size = (width, 790), fontsize = 13)
    title = Observable(time_label(common[1]))
    Label(fig[0, 1:3], title, fontsize = 15, font = :bold, tellwidth = false)
    Label(fig[1, 0], "Liquid water path (cloud + rain, g m⁻²)", rotation = π/2, tellheight = false)
    Label(fig[2, 0], "Surface rain rate (mm day⁻¹)", rotation = π/2, tellheight = false)
    lwp_obs = [Observable(zeros(Float32, length(x), length(y))) for _ in series]
    rain_obs = [Observable(zeros(Float32, length(x), length(y))) for _ in series]
    member_titles = [Observable(string(s.name)) for s in series]
    for (k, s) in enumerate(series)
        ax1 = Axis(fig[1, k], title = member_titles[k], titlesize = 12, aspect = DataAspect(), xlabel = "", ylabel = k == 1 ? "y (km)" : "")
        hm1 = heatmap!(ax1, x, y, lwp_obs[k], colormap = :viridis, colorrange = (0, lwp_max))
        ax2 = Axis(fig[2, k], aspect = DataAspect(), xlabel = "x (km)", ylabel = k == 1 ? "y (km)" : "")
        hm2 = heatmap!(ax2, x, y, rain_obs[k], colormap = :dense, colorrange = (0, rain_max))
        k > 1 && (hideydecorations!(ax1, grid = false); hideydecorations!(ax2, grid = false))
        hidexdecorations!(ax1, grid = false)
        k == length(series) && (Colorbar(fig[1, 4], hm1, width = 10); Colorbar(fig[2, 4], hm2, width = 10))
    end
    function update!(n)
        t = common[n]
        for (k, s) in enumerate(series)
            i = index(s, t)
            lwp = 1e3 .* (interior(s.cwp[i], :, :, 1) .+ interior(s.rwp[i], :, :, 1))
            rr = 86400 .* interior(s.rain[i], :, :, 1)
            lwp_obs[k][] = lwp
            rain_obs[k][] = rr
            member_titles[k][] = @sprintf("%s\n⟨LWP⟩ %.0f g m⁻², ⟨rain⟩ %.3f mm d⁻¹", s.name, mean(lwp), mean(rr))
        end
        title[] = time_label(t)
    end
    mp4 = mkpath_for(joinpath(ctrl_dir, "animations", "seastarr_planview_lwp_rain_ctrl_n100_n030.mp4"))
    png = joinpath(ctrl_dir, "animations", "seastarr_planview_lwp_rain_ctrl_n100_n030_montage.png")
    record_with_montage(fig, update!, frames, mp4, png; montage_frames = montage_indices(frames))
    @info "plan views: $(length(frames)) frames (15 min), t = $(common[first(frames)]/3600) – $(common[last(frames)]/3600) h; LWP range 0–$(round(lwp_max)) g m⁻², rain 0–$(round(rain_max; digits=2)) mm d⁻¹" mp4 fmt_mb(mp4) png
    return (mp4, png, length(frames), common[last(frames)], lwp_max, rain_max)
end

#####
##### 2. CTRL cross-section: cloud water, rain water, aerosol number
#####

function cross_section_animation()
    # Prefer the x–z slices saved at the animation cadence (runs after the 3-min 2-D output was
    # introduced); fall back to the hourly 3D fields for older runs.
    file2d = joinpath(ctrl_dir, "sea_starr_ctrl_2d.jld2")
    has_slices = jldopen(f -> haskey(f["timeseries"], "qᶜˡ_xz"), file2d)
    file = has_slices ? file2d : joinpath(ctrl_dir, "sea_starr_ctrl_3d.jld2")
    suffix = has_slices ? "_xz" : ""
    qc = FieldTimeSeries(file, "qᶜˡ" * suffix; backend = OnDisk())
    qr = FieldTimeSeries(file, "qʳ" * suffix; backend = OnDisk())
    na = FieldTimeSeries(file, "nᵃ" * suffix; backend = OnDisk())
    nc = FieldTimeSeries(file, "nᶜˡ" * suffix; backend = OnDisk())
    grid = qc.grid
    j = size(grid, 2) ÷ 2
    x = xnodes(grid, Center()) ./ 1e3
    z = znodes(grid, Center())
    kmax = findlast(≤(3000), z)
    zc = z[1:kmax]
    frames = 1:min(length(qc.times), max_frames)
    slice(fts, n) = has_slices ? Array(interior(fts[n], :, 1, 1:kmax)) : Array(interior(fts[n], :, j, 1:kmax))
    cadence = length(qc.times) > 1 ? (qc.times[2] - qc.times[1]) / 60 : NaN

    fig = Figure(size = (width, 900), fontsize = 13)
    title = Observable(time_label(qc.times[1]))
    Label(fig[0, 1:2], title, fontsize = 15, font = :bold, tellwidth = false)
    qc_obs = Observable(zeros(Float32, length(x), kmax)); qr_obs = Observable(zeros(Float32, length(x), kmax))
    na_obs = Observable(zeros(Float32, length(x), kmax)); nc_obs = Observable(zeros(Float32, length(x), kmax))
    ax1 = Axis(fig[1, 1], title = "Cloud water qᶜˡ (g kg⁻¹)", ylabel = "z (m)")
    hm1 = heatmap!(ax1, x, zc, qc_obs, colormap = :Blues, colorrange = (0, 1.2))
    Colorbar(fig[1, 2], hm1, width = 10)
    ax2 = Axis(fig[2, 1], title = "Rain water qʳ (g kg⁻¹, log scale)", ylabel = "z (m)")
    hm2 = heatmap!(ax2, x, zc, qr_obs, colormap = :Purples, colorrange = (1e-4, 1e-1), colorscale = log10, lowclip = :transparent)
    Colorbar(fig[2, 2], hm2, width = 10)
    ax3 = Axis(fig[3, 1], title = "Interstitial aerosol nᵃ (mg⁻¹; contours: droplet number nᶜˡ at 50, 100, 200 mg⁻¹)", xlabel = "x (km)", ylabel = "z (m)")
    hm3 = heatmap!(ax3, x, zc, na_obs, colormap = :magma, colorrange = (0, 1200))
    contour!(ax3, x, zc, nc_obs, levels = [50, 100, 200], color = :cyan, linewidth = 0.8)
    Colorbar(fig[3, 2], hm3, width = 10)
    hidexdecorations!(ax1, grid = false); hidexdecorations!(ax2, grid = false)
    function update!(n)
        qc_obs[] = 1e3 .* slice(qc, n)
        qr_obs[] = max.(1e3 .* slice(qr, n), 1f-6)
        na_obs[] = 1e-6 .* slice(na, n)
        nc_obs[] = 1e-6 .* slice(nc, n)
        title[] = string("CTRL x–z section at y = ", round(ynodes(grid, Center())[j] / 1e3; digits = 2), " km   —   ", time_label(qc.times[n]))
    end
    mp4 = mkpath_for(joinpath(ctrl_dir, "animations", "seastarr_ctrl_xz_cloud_rain_aerosol.mp4"))
    png = joinpath(ctrl_dir, "animations", "seastarr_ctrl_xz_cloud_rain_aerosol_montage.png")
    record_with_montage(fig, update!, frames, mp4, png; montage_frames = montage_indices(frames))
    @info "cross-section: $(length(frames)) frames ($(cadence) min cadence, $(has_slices ? "saved slices" : "hourly 3D")) to t = $(qc.times[last(frames)]/3600) h" mp4 fmt_mb(mp4) png
    return (mp4, png, length(frames), qc.times[last(frames)])
end

#####
##### 3. Time–height: horizontal-mean aerosol and droplet number
#####

function time_height_animation()
    file = joinpath(ctrl_dir, "sea_starr_ctrl_statistics.jld2")
    na = FieldTimeSeries(file, "nᵃ"); nc = FieldTimeSeries(file, "nᶜˡ"); qc = FieldTimeSeries(file, "qᶜˡ")
    zi = FieldTimeSeries(joinpath(ctrl_dir, "sea_starr_ctrl_timeseries.jld2"), "zi")
    z = znodes(na.grid, Center())
    kmax = findlast(≤(3000), z)
    zc = z[1:kmax]
    t = na.times ./ 3600
    frames = 1:min(length(t), max_frames)
    column(fts, n) = vec(Array(interior(fts[n], 1, 1, 1:kmax)))
    # (t, z) matrices, as Makie's heatmap(x = t, y = z, data) expects
    NA = Matrix{Float32}(undef, length(frames), kmax); NC = similar(NA); QC = similar(NA)
    for n in frames
        NA[n, :] .= 1e-6 .* column(na, n)     # mg⁻¹
        NC[n, :] .= 1e-6 .* column(nc, n)
        QC[n, :] .= 1e3 .* column(qc, n)      # g kg⁻¹
    end
    zi_t = zi.times ./ 3600
    zi_v = [zi[n][1, 1, 1] for n in eachindex(zi.times)]

    fig = Figure(size = (width, 820), fontsize = 13)
    title = Observable(time_label(na.times[1]))
    Label(fig[0, 1:2], title, fontsize = 15, font = :bold, tellwidth = false)
    na_obs = Observable(fill(NaN32, size(NA))); nc_obs = Observable(fill(NaN32, size(NC)))
    zi_obs = Observable(Point2f[])
    ax1 = Axis(fig[1, 1], title = "Horizontal-mean interstitial aerosol nᵃ (mg⁻¹), with ⟨zᵢ⟩ (white)", ylabel = "z (m)")
    hm1 = heatmap!(ax1, t[frames], zc, na_obs, colormap = :magma, colorrange = (0, 1200), nan_color = :transparent)
    lines!(ax1, zi_obs, color = :white, linewidth = 1.5)
    Colorbar(fig[1, 2], hm1, width = 10)
    ax2 = Axis(fig[2, 1], title = "Horizontal-mean droplet number nᶜˡ (mg⁻¹); contour: ⟨qᶜˡ⟩ = 0.01 g kg⁻¹", xlabel = "hours since 2017-08-15 21 UTC", ylabel = "z (m)")
    hm2 = heatmap!(ax2, t[frames], zc, nc_obs, colormap = :viridis, colorrange = (0, 300), nan_color = :transparent)
    qc_obs = Observable(fill(0f0, size(QC)))
    contour!(ax2, t[frames], zc, qc_obs, levels = [0.01], color = :white, linewidth = 1)
    Colorbar(fig[2, 2], hm2, width = 10)
    xlims!(ax1, t[first(frames)], max(t[last(frames)], t[first(frames)] + 1)); xlims!(ax2, t[first(frames)], max(t[last(frames)], t[first(frames)] + 1))
    hidexdecorations!(ax1, grid = false)
    function update!(n)
        A = fill(NaN32, size(NA)); A[1:n, :] .= NA[1:n, :]; na_obs[] = A
        B = fill(NaN32, size(NC)); B[1:n, :] .= NC[1:n, :]; nc_obs[] = B
        C = fill(0f0, size(QC)); C[1:n, :] .= QC[1:n, :]; qc_obs[] = C
        m = zi_t .≤ t[n]
        zi_obs[] = Point2f.(zi_t[m], zi_v[m])
        title[] = time_label(na.times[n])
    end
    mp4 = mkpath_for(joinpath(ctrl_dir, "animations", "seastarr_ctrl_time_height_aerosol_droplets.mp4"))
    png = joinpath(ctrl_dir, "animations", "seastarr_ctrl_time_height_aerosol_droplets_montage.png")
    record_with_montage(fig, update!, frames, mp4, png; montage_frames = montage_indices(frames))
    @info "time–height: $(length(frames)) frames (15 min) to t = $(t[last(frames)]) h" mp4 fmt_mb(mp4) png
    return (mp4, png, length(frames), na.times[last(frames)])
end

#####
##### Run, copy, README
#####

results = Dict{String, Any}()
results["planview"] = plan_view_animation()
results["xz"] = cross_section_animation()
results["time_height"] = time_height_animation()

mkpath(copy_to)
copied = String[]
for (_, r) in results, path in r[1:2]
    cp(path, joinpath(copy_to, basename(path)); force = true)
    push!(copied, basename(path))
end
open(joinpath(copy_to, "README.md"), "w") do io
    println(io, "# SEA STARR animations\n")
    println(io, "Rendered ", Dates.format(now(UTC), "yyyy-mm-dd HH:MM"), " UTC by `analysis/animate_sea_starr.jl` (BreezeLab commit ", commit, ") from the ",
            "JLD2 output of the 66 h runs on wpcluster: CTRL `", ctrl_dir, "`, N100 `", n100_dir, "`, N030 `", n030_dir, "` (runs still in progress ",
            "when rendered; the frames stop at the last record available). MP4 = H.264, ", width, " px wide, ", framerate, " frames per second.\n")
    r = results["planview"]
    println(io, "- `", basename(r[1]), "` / `", basename(r[2]), "`: plan views of the liquid water path (cloud + rain water, g m⁻², top row, shared range 0–",
            round(Int, r[5]), ") and the surface rain rate (mm day⁻¹, bottom row, shared range 0–", round(r[6]; digits = 2),
            ") for CTRL, N100 and N030 at the same times; 15 min cadence, ", r[3], " frames to t = ", round(r[4] / 3600; digits = 2),
            " h. Title: UTC and the local solar time at the composite trajectory's mean longitude (GEOS-5 raw trajectories).")
    r = results["xz"]
    println(io, "- `", basename(r[1]), "` / `", basename(r[2]), "`: CTRL x–z cross-section through the domain centre (y = 4.8 km) below 3 km: cloud water qᶜˡ (g kg⁻¹), ",
            "rain water qʳ (g kg⁻¹, log scale 10⁻⁴–10⁻¹) and interstitial aerosol number nᵃ (mg⁻¹ = 10⁶ kg⁻¹) with droplet-number contours (mg⁻¹); hourly 3D output, ",
            r[3], " frames to t = ", round(r[4] / 3600; digits = 2), " h.")
    r = results["time_height"]
    println(io, "- `", basename(r[1]), "` / `", basename(r[2]), "`: CTRL time–height of the horizontal-mean interstitial aerosol nᵃ and droplet number nᶜˡ (mg⁻¹, 15 min averages) ",
            "below 3 km, revealed as time advances, with the domain-mean inversion height (white) and the ⟨qᶜˡ⟩ = 0.01 g kg⁻¹ contour; ", r[3], " frames to t = ",
            round(r[4] / 3600; digits = 2), " h.")
    println(io, "\nTime axis: seconds since the driver start 2017-08-15 21:00 UTC. The runs carry the labelled protocol departures listed in docs/cases/sea_starr.md ",
            "(no aerosol optics, radiation column ending at 6.5 km, fixed solar coordinate, linear nudging ramp): these are exploratory CTRL/N100/N030 integrations, not validated results.")
end
println("copied to $copy_to: ", join(copied, ", "))
for (k, r) in results
    println(k, ": ", r[1], " (", fmt_mb(r[1]), "), ", r[2], " (", fmt_mb(r[2]), ")")
end
