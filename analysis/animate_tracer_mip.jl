# Animations of a TRACER-MIP outer-domain run (CairoMakie `record`, H.264 MP4 ≤ 12 MB, ≤ 1000 px wide)
# plus 4-frame PNG montages (06, 12, 18, 00 UTC). Reads the run's 10-min output only:
#   outer_slices.jld2   T_sfc, u_sfc, v_sfc (lowest level), w_2km, qᶜˡ_2km            — whole 1500 km domain
#   outer_surface.jld2  rain (surface rain flux, kg m⁻² s⁻¹), sea (1 at sea), T_land   — whole domain
#   inner_region/inner_region_state.jld2  T p qᵛ qᶜˡ qʳ qⁱ u v w ρ (3-D, inner box + halo)
#
#   julia --project=cases/tracer_mip analysis/animate_tracer_mip.jl <run_dir> [<out_dir>]
#
# Environment: ANIMATION_FPS (8), ANIMATION_COMPRESSION (H.264 CRF-like, 28; raised automatically until each
# MP4 is under ANIMATION_MAX_MB = 12), CASE_START ("2022-08-07T06:00:00"), LOCAL_UTC_OFFSET_HOURS (-5, CDT).
# Run on a CPU node from a *snapshot copy* of a running job's files, never on the files being written.

using CairoMakie, Oceananigans, Oceananigans.Units, Breeze, Printf, Statistics
using Oceananigans.Fields: interior
using Oceananigans.Grids: λnodes, φnodes, znode, Center, Face
using Dates: DateTime, Second, Hour, format
using TOML: TOML

run_dir = abspath(get(ARGS, 1, "."))
out_dir = abspath(get(ARGS, 2, joinpath(run_dir, "animations")))
mkpath(out_dir)
framerate = parse(Int, get(ENV, "ANIMATION_FPS", "8"))
compression₀ = parse(Int, get(ENV, "ANIMATION_COMPRESSION", "28"))
max_mb = parse(Float64, get(ENV, "ANIMATION_MAX_MB", "12"))
case_start = DateTime(get(ENV, "CASE_START", "2022-08-07T06:00:00"))
local_offset = Hour(parse(Int, get(ENV, "LOCAL_UTC_OFFSET_HOURS", "-5")))
site = (λ = -95.0792, φ = 29.4719)              # domain centre (KHGX NEXRAD), the roadmap's grid centre
inner_half_width_deg = (λ = 1.2912, φ = 1.1242)  # the 250 km inner box (protocol extents) in degrees
CairoMakie.activate!(type = "png")

stamp(t) = (utc = case_start + Second(round(Int, t)); @sprintf("%s UTC (%s CDT)", format(utc, "yyyy-mm-dd HH:MM"), format(utc + local_offset, "HH:MM")))

#####
##### Load
#####

slices = Dict(name => FieldTimeSeries(joinpath(run_dir, "outer_slices.jld2"), name; backend = OnDisk())
              for name in ("T_sfc", "u_sfc", "v_sfc", "w_2km", "qᶜˡ_2km"))
surface = Dict(name => FieldTimeSeries(joinpath(run_dir, "outer_surface.jld2"), name; backend = OnDisk())
               for name in ("rain", "sea"))
inner_file = joinpath(run_dir, "inner_region", "inner_region_state.jld2")
inner = Dict(name => FieldTimeSeries(inner_file, name; backend = OnDisk())
             for name in ("T", "qᶜˡ", "qʳ", "qⁱ", "u", "v", "w", "ρ"))

grid = slices["T_sfc"].grid
Nx, Ny, Nz = size(grid)
λ = Array(λnodes(grid, Center())); φ = Array(φnodes(grid, Center()))
times = slices["T_sfc"].times
n_common = minimum(length(fts.times) for fts in (values(slices)..., values(surface)..., values(inner)...))
frames = 1:n_common
@info "Animating $(n_common) frames: $(stamp(times[1])) → $(stamp(times[n_common]))"

field2d(fts, n) = dropdims(Array(interior(fts[n])); dims = 3)
sea = field2d(surface["sea"], 1)
land_mask = sea .< 0.5

# Inner window (indices saved with the fields) and its column heights (static terrain-following grid)
i_range, j_range, _ = inner["T"].indices
ni, nj = length(i_range), length(j_range)
λi, φi = λ[i_range], φ[j_range]
zf = Array{Float64}(undef, ni, nj, Nz + 1)
for (jj, j) in enumerate(j_range), (ii, i) in enumerate(i_range), k in 1:Nz+1
    zf[ii, jj, k] = znode(i, j, k, grid, Center(), Center(), Face())
end
dz = diff(zf; dims = 3)
terrain = zf[:, :, 1]
i_site = argmin(abs.(λi .- site.λ)); j_site = argmin(abs.(φi .- site.φ))
z_site = [znode(i_range[i_site], j_range[j_site], k, grid, Center(), Center(), Center()) for k in 1:Nz]

inner3d(name, n) = Array(interior(inner[name][n]))
function inner_fields(n)
    ρ = inner3d("ρ", n); qᶜˡ = inner3d("qᶜˡ", n); qʳ = inner3d("qʳ", n); qⁱ = inner3d("qⁱ", n)
    lwp = dropdims(sum(ρ .* qᶜˡ .* dz; dims = 3); dims = 3)
    rwp = dropdims(sum(ρ .* qʳ .* dz; dims = 3); dims = 3)
    T₁ = inner3d("T", n)[:, :, 1]; u₁ = inner3d("u", n)[:, :, 1]; v₁ = inner3d("v", n)[:, :, 1]
    condensate = 1e3 .* (qᶜˡ[:, j_site, :] .+ qʳ[:, j_site, :] .+ qⁱ[:, j_site, :])
    w = inner3d("w", n)[:, j_site, :]
    rain = 3600 .* field2d(surface["rain"], n)[i_range, j_range]     # mm/h (1 kg m⁻² = 1 mm)
    return (; T₁, u₁, v₁, lwp = 1e3 .* lwp, rwp = 1e3 .* rwp, condensate, w, rain)
end

# Makie ≥ 0.24 renamed `arrows` → `arrows2d`; support both.
const new_arrows = isdefined(CairoMakie.Makie, :arrows2d!)
arrowfn! = new_arrows ? CairoMakie.Makie.arrows2d! : CairoMakie.Makie.arrows!
arrow_style = new_arrows ? (shaftwidth = 1, tipwidth = 4, tiplength = 5) : (linewidth = 0.8, arrowsize = 6)
function wind_arrows!(ax, xs, ys, u::Observable, v::Observable; step, lengthscale)
    xi = 1:step:length(xs); yi = 1:step:length(ys)
    pts = [Point2f(x, y) for x in xs[xi], y in ys[yi]]
    dirs = lift(u, v) do uu, vv
        [Vec2f(uu[i, j], vv[i, j]) for i in xi, j in yi]
    end
    arrowfn!(ax, vec(pts), lift(vec, dirs); lengthscale, color = :black, arrow_style...)
end

function finish_record(fig, path, draw!, frames)
    compression = compression₀
    while true
        record(fig, path, frames; framerate, compression, px_per_unit = 1) do n
            draw!(n)
        end
        mb = filesize(path) / 1e6
        @info @sprintf("%s: %.1f MB (compression %d)", basename(path), mb, compression)
        (mb ≤ max_mb || compression ≥ 45) && return mb
        compression += 4
    end
end

montage_hours = [0, 6, 12, 18]        # 06, 12, 18, 00 UTC
function montage_frames()
    picks = Int[]; labels = String[]
    for h in montage_hours
        n = findfirst(t -> abs(t - 3600h) < 1, times[frames])
        if isnothing(n)
            push!(picks, n_common); push!(labels, "latest available")
        else
            push!(picks, n); push!(labels, "")
        end
    end
    return picks, labels
end

#####
##### 1. Plan view, whole domain
#####

function plan_view(; zoom = false)
    xs, ys = zoom ? (λi, φi) : (λ, φ)
    fig = Figure(size = (1000, 540), fontsize = 12)
    title = Observable("")
    Label(fig[0, 1:4], title; fontsize = 15, tellwidth = false)
    ax1 = Axis(fig[1, 1]; xlabel = "longitude", ylabel = "latitude", aspect = DataAspect(),
               title = zoom ? "lowest-level T (K), wind" : "lowest-level T (K), 10 m-level wind")
    ax2 = Axis(fig[1, 3]; xlabel = "longitude", aspect = DataAspect(),
               title = zoom ? "LWP (g m⁻²), rain rate contours 1/5/20 mm h⁻¹" : "rain rate (mm h⁻¹); cloud water at 2 km (white)")
    T = Observable(zeros(Float32, length(xs), length(ys))); u = Observable(zeros(Float32, length(xs), length(ys))); v = Observable(zeros(Float32, length(xs), length(ys)))
    B = Observable(zeros(Float32, length(xs), length(ys))); C = Observable(zeros(Float32, length(xs), length(ys)))
    hm1 = heatmap!(ax1, xs, ys, T; colormap = :thermal, colorrange = zoom ? (294, 310) : (285, 312))
    Colorbar(fig[1, 2], hm1)
    if zoom
        hm2 = heatmap!(ax2, xs, ys, B; colormap = :dense, colorrange = (0, 800))
        contour!(ax2, xs, ys, C; levels = [1, 5, 20], color = [:orange, :red, :darkred], linewidth = 1)
    else
        hm2 = heatmap!(ax2, xs, ys, B; colormap = Reverse(:Blues), colorrange = (0, 20), lowclip = :white)
        contour!(ax2, xs, ys, C; levels = [0.05, 0.5], color = :grey40, linewidth = 0.6)
    end
    Colorbar(fig[1, 4], hm2)
    coast = zoom ? sea[i_range, j_range] : sea
    for ax in (ax1, ax2)
        contour!(ax, xs, ys, coast; levels = [0.5], color = :black, linewidth = 1)
        scatter!(ax, [site.λ], [site.φ]; marker = :star5, color = :white, strokecolor = :black, strokewidth = 1, markersize = 14)
        if !zoom   # inner box outline
            λb = site.λ .+ inner_half_width_deg.λ .* [-1, 1, 1, -1, -1]; φb = site.φ .+ inner_half_width_deg.φ .* [-1, -1, 1, 1, -1]
            lines!(ax, λb, φb; color = :black, linestyle = :dash, linewidth = 1)
        end
    end
    wind_arrows!(ax1, xs, ys, u, v; step = zoom ? 8 : 30, lengthscale = zoom ? 0.012 : 0.03)
    function draw!(n)
        if zoom
            f = inner_fields(n)
            T[] = f.T₁; u[] = f.u₁; v[] = f.v₁; B[] = f.lwp; C[] = f.rain
        else
            T[] = field2d(slices["T_sfc"], n); u[] = field2d(slices["u_sfc"], n); v[] = field2d(slices["v_sfc"], n)
            B[] = 3600 .* field2d(surface["rain"], n); C[] = 1e3 .* field2d(slices["qᶜˡ_2km"], n)
        end
        title[] = (zoom ? "TRACER-MIP Aug 7 2022 outer control, inner 250 km box — " : "TRACER-MIP Aug 7 2022 outer control (2 km, 1500 km) — ") * stamp(times[n])
    end
    return fig, draw!
end

#####
##### 3. Vertical section through the site (inner-region 3-D state, j row at 29.47 N)
#####

function section_view()
    fig = Figure(size = (1000, 640), fontsize = 12)
    title = Observable("")
    Label(fig[0, 1:2], title; fontsize = 15, tellwidth = false)
    ax1 = Axis(fig[1, 1]; ylabel = "z (km)", title = "cloud + rain + ice water (g kg⁻¹)", limits = (nothing, (0, 16)))
    ax2 = Axis(fig[2, 1]; xlabel = "longitude (section at 29.47 N through the site)", ylabel = "z (km)", title = "w (m s⁻¹)", limits = (nothing, (0, 16)))
    zkm = z_site ./ 1e3
    Q = Observable(zeros(Float32, ni, Nz)); W = Observable(zeros(Float32, ni, Nz))
    hm1 = heatmap!(ax1, λi, zkm, Q; colormap = :dense, colorrange = (0, 3), lowclip = :white)
    hm2 = heatmap!(ax2, λi, zkm, W; colormap = :balance, colorrange = (-4, 4))
    Colorbar(fig[1, 2], hm1); Colorbar(fig[2, 2], hm2)
    for ax in (ax1, ax2)
        lines!(ax, λi, terrain[:, j_site] ./ 1e3; color = :black)
        vlines!(ax, [site.λ]; color = :grey30, linestyle = :dot)
    end
    function draw!(n)
        f = inner_fields(n)
        Q[] = f.condensate; W[] = f.w
        title[] = "TRACER-MIP Aug 7 2022 outer control — x–z section through the site — " * stamp(times[n])
    end
    return fig, draw!
end

#####
##### Render
#####

products = [("plan_view_domain", () -> plan_view(; zoom = false)),
            ("plan_view_inner_box", () -> plan_view(; zoom = true)),
            ("section_site", section_view)]
report = Dict{String, Any}("run_dir" => run_dir, "n_frames" => n_common, "first" => stamp(times[1]), "last" => stamp(times[n_common]),
                           "cadence_s" => length(times) > 1 ? times[2] - times[1] : 0.0, "framerate" => framerate)
picks, labels = montage_frames()
for (name, make) in products
    fig, draw! = make()
    mp4 = joinpath(out_dir, "tracer_mip_aug07_$(name).mp4")
    mb = finish_record(fig, mp4, draw!, frames)
    report[name * "_mp4_MB"] = round(mb; digits = 2)
    # 4-frame montage: one full figure per pick, rendered to PNG then composed
    panels = map(enumerate(picks)) do (k, n)
        draw!(n)
        isempty(labels[k]) || (fig.content[1].text[] = fig.content[1].text[] * " — " * labels[k])
        colorbuffer(fig)
    end
    montage = Figure(size = (2 * size(panels[1], 2) ÷ 2, 2 * size(panels[1], 1) ÷ 2))
    for (k, img) in enumerate(panels)
        ax = Axis(montage[(k - 1) ÷ 2 + 1, (k - 1) % 2 + 1]; aspect = DataAspect()); hidedecorations!(ax); hidespines!(ax)
        image!(ax, rotr90(img))
    end
    png = joinpath(out_dir, "tracer_mip_aug07_$(name)_montage.png")
    save(png, montage; px_per_unit = 1)
    report[name * "_montage_MB"] = round(filesize(png) / 1e6; digits = 2)
end
open(joinpath(out_dir, "animations.toml"), "w") do io
    TOML.print(io, report)
end
@info "done" report
