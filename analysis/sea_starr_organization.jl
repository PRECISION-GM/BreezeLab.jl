# Mesoscale organization of SEA STARR cloud fields from the saved plan-view cloud water path (2-D output):
# cloud cover, cloud objects and clear holes (connected regions on the periodic domain), and the spectral
# (integral) length scale of the LWP field — e.g. to compare the 9.6 km (192²) and 19.2 km (384²) CTRL domains.
#
#   julia --project analysis/sea_starr_organization.jl --runs "CTRL 192²|CTRL|<dirs>;CTRL 384²|CTRL|<dirs>" [--out <dir>] [--every 1]
#
# A column is cloudy when its cloud water path exceeds 5 g m⁻² (the model's `cloud_fraction` definition).
# Objects use 4-connectivity with periodic wrap-around. Equivalent diameter = √(4A/π). The integral length scale is
# ℓ = 2π Σ P(k) / Σ k P(k) over the radially averaged power spectrum of the LWP anomaly (Jonker et al. 1999 style);
# it is bounded by the domain size, which is itself part of the domain-size sensitivity being measured.

using Oceananigans
using Oceananigans.Fields: interior
using Oceananigans.OutputReaders: OnDisk
using CairoMakie
using Statistics
using Printf
using Dates
using TOML

const FFTW = Oceananigans.Solvers.FFTW
include(joinpath(@__DIR__, "sea_starr_segments.jl"))

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
runs = map(split(get(opts, "runs", "CTRL 192²|CTRL|runs/ctrl_full_66h_job211;CTRL 384²|CTRL|runs/ctrl_full_66h_384_job230"), ';')) do entry
    label, member, dirs = split(entry, '|')
    (String(label), "sea_starr_" * lowercase(member), segment_dirs(dirs))
end
out_dir = get(opts, "out", "runs/organization_analysis")
every = parse(Float64, get(opts, "every", "1")) * 3600
threshold = parse(Float64, get(opts, "threshold", "5")) / 1e3     # kg m⁻²
mkpath(out_dir)
commit = try strip(read(`git -C $(@__DIR__) rev-parse --short HEAD`, String)) catch; "unknown" end
CairoMakie.activate!(type = "png", px_per_unit = 1.5)

# Connected components of `mask` with 4-connectivity on a doubly periodic grid; returns component sizes (cells).
function component_sizes(mask::AbstractMatrix{Bool})
    Nx, Ny = size(mask)
    seen = falses(Nx, Ny)
    sizes = Int[]
    stack = Tuple{Int, Int}[]
    for j in 1:Ny, i in 1:Nx
        (mask[i, j] && !seen[i, j]) || continue
        seen[i, j] = true; push!(stack, (i, j)); n = 0
        while !isempty(stack)
            a, b = pop!(stack); n += 1
            for (da, db) in ((1, 0), (-1, 0), (0, 1), (0, -1))
                ia = mod1(a + da, Nx); jb = mod1(b + db, Ny)
                if mask[ia, jb] && !seen[ia, jb]
                    seen[ia, jb] = true; push!(stack, (ia, jb))
                end
            end
        end
        push!(sizes, n)
    end
    return sizes
end

# Integral length scale [m] from the radially averaged power spectrum of the anomaly of `field` (Δ: grid spacing).
function integral_length_scale(field, Δ)
    Nx, Ny = size(field)
    a = field .- mean(field)
    P = abs2.(FFTW.fft(a))
    kx = 2π .* FFTW.fftfreq(Nx, 1 / Δ); ky = 2π .* FFTW.fftfreq(Ny, 1 / Δ)
    dk = 2π / (Nx * Δ)
    nbins = Nx ÷ 2
    Pk = zeros(nbins); kk = [(b - 0.5) * dk for b in 1:nbins]
    for j in 1:Ny, i in 1:Nx
        k = sqrt(kx[i]^2 + ky[j]^2)
        b = floor(Int, k / dk) + 1
        (1 ≤ b ≤ nbins && k > 0) && (Pk[b] += P[i, j])
    end
    return sum(Pk) > 0 ? 2π * sum(Pk) / sum(kk .* Pk) : NaN
end

results = map(runs) do (label, prefix, dirs)
    cwp = stitched_series(dirs, prefix * "_2d.jld2", "cwp"; backend = OnDisk())
    grid = cwp.grid
    Nx, Ny = size(grid, 1), size(grid, 2)
    Δ = grid.Lx / Nx
    area_km2 = grid.Lx * grid.Ly / 1e6
    picks = [n for n in eachindex(cwp.times) if abs(cwp.times[n] / every - round(cwp.times[n] / every)) < 1e-6]
    rows = map(picks) do n
        L = Array(interior(cwp[n], :, :, 1))
        cloudy = L .> threshold
        objects = component_sizes(cloudy); holes = component_sizes(.!cloudy)
        d(sizes) = isempty(sizes) ? NaN : mean(sqrt.(4 .* sizes .* Δ^2 ./ π)) / 1e3
        (; t = cwp.times[n], cover = mean(cloudy), lwp = 1e3 * mean(L), lwp_std = 1e3 * std(L),
           n_objects = length(objects) / area_km2 * 100, object_diameter_km = d(objects),
           largest_object_fraction = isempty(objects) ? 0.0 : maximum(objects) / (Nx * Ny),
           n_holes = length(holes) / area_km2 * 100, hole_diameter_km = d(holes),
           largest_hole_fraction = isempty(holes) ? 0.0 : maximum(holes) / (Nx * Ny),
           length_scale_km = integral_length_scale(L, Δ) / 1e3)
    end
    @info "$label: $(length(rows)) frames, domain $(grid.Lx / 1e3) km, Δ = $Δ m"
    (; label, dirs, domain_km = grid.Lx / 1e3, rows)
end

#####
##### Figure and summary
#####

fig = Figure(size = (1150, 1000), fontsize = 12)
Label(fig[0, 1:2], "Mesoscale organization from the cloud water path (cloudy: CWP > $(round(Int, 1e3threshold)) g m⁻²)", font = :bold, tellwidth = false)
panels = ((1, 1, "cloud cover", r -> r.cover), (1, 2, "integral length scale ℓ (km)", r -> r.length_scale_km),
          (2, 1, "cloud objects per 100 km²", r -> r.n_objects), (2, 2, "mean object equivalent diameter (km)", r -> r.object_diameter_km),
          (3, 1, "clear holes per 100 km²", r -> r.n_holes), (3, 2, "mean clear-hole equivalent diameter (km)", r -> r.hole_diameter_km),
          (4, 1, "largest cloud object (fraction of domain)", r -> r.largest_object_fraction), (4, 2, "CWP standard deviation (g m⁻²)", r -> r.lwp_std))
for (r, c, label, getter) in panels
    ax = Axis(fig[r, c], ylabel = label, xlabel = r == 4 ? "hours since 2017-08-15 21 UTC" : "")
    for (k, res) in enumerate(results)
        lines!(ax, [x.t / 3600 for x in res.rows], [getter(x) for x in res.rows], color = Makie.wong_colors()[k], label = res.label)
    end
    (r, c) == (1, 1) && axislegend(ax, position = :lb)
end
save(joinpath(out_dir, "organization.png"), fig)

day_ranges = (("day1", 0, 24), ("day2", 24, 48), ("day3", 48, 66.01))
numbers = Dict{String, Any}("commit" => commit, "generated" => string(now(UTC)), "threshold_kg_m2" => threshold)
for res in results
    k = replace(res.label, r"[^A-Za-z0-9]+" => "_")
    numbers["$(k)_domain_km"] = res.domain_km
    for (lab, a, b) in day_ranges, (q, f) in (("cover", r -> r.cover), ("length_scale_km", r -> r.length_scale_km),
                                             ("object_diameter_km", r -> r.object_diameter_km), ("hole_diameter_km", r -> r.hole_diameter_km),
                                             ("n_objects_per_100km2", r -> r.n_objects), ("n_holes_per_100km2", r -> r.n_holes))
        sel = [f(x) for x in res.rows if a ≤ x.t / 3600 < b]
        sel = filter(isfinite, sel)
        isempty(sel) || (numbers["$(k)_$(lab)_$(q)"] = mean(sel))
    end
    last = res.rows[end]
    numbers["$(k)_final_cover"] = last.cover; numbers["$(k)_final_length_scale_km"] = last.length_scale_km
    numbers["$(k)_final_object_diameter_km"] = last.object_diameter_km; numbers["$(k)_final_largest_hole_fraction"] = last.largest_hole_fraction
end
open(joinpath(out_dir, "organization_summary.toml"), "w") do io
    TOML.print(io, numbers)
end
open(joinpath(out_dir, "organization_summary.md"), "w") do io
    println(io, "# Mesoscale organization (commit ", commit, ", ", Dates.format(now(UTC), "yyyy-mm-dd HH:MM"), " UTC)\n")
    labels = [replace(r.label, r"[^A-Za-z0-9]+" => "_") for r in results]
    println(io, "| quantity | ", join([r.label for r in results], " | "), " |\n|", repeat(" --- |", length(results) + 1))
    for q in ("cover", "length_scale_km", "object_diameter_km", "hole_diameter_km", "n_objects_per_100km2", "n_holes_per_100km2"), (lab, _, _) in day_ranges
        vals = [(v = get(numbers, "$(k)_$(lab)_$(q)", NaN); @sprintf("%.3g", v)) for k in labels]
        println(io, "| ", lab, " ", q, " | ", join(vals, " | "), " |")
    end
    for q in ("final_cover", "final_length_scale_km", "final_object_diameter_km", "final_largest_hole_fraction")
        println(io, "| ", q, " | ", join([@sprintf("%.3g", numbers["$(k)_$(q)"]) for k in labels], " | "), " |")
    end
end
println("organization analysis written to ", out_dir)
for k in sort(collect(keys(numbers)))
    println(k, " = ", numbers[k])
end
