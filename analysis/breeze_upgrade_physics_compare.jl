# Compare the CSVs written by `analysis/breeze_upgrade_physics_check.jl` under two environments.
#
#   julia analysis/breeze_upgrade_physics_compare.jl <old dir> <new dir>
#
# Prints, per case, the final-record values (old, new, new − old) and the drift of the prognostic
# water and of ∫ρθ dz over the run, as a Markdown table.

using Printf

function read_csv(path)
    lines = readlines(path)
    header = Symbol.(split(lines[1], ','))
    rows = [parse.(Float64, split(line, ',')) for line in lines[2:end]]
    return header, rows
end

const COLUMNS = (:LWP, :RWP, :IWP, :mean_T, :mean_θ, :mean_qv, :Nc_cloudy, :rain_mm_day, :bottom_flux_mm_day, :max_w)
const UNITS = Dict(:LWP => "g m⁻²", :RWP => "g m⁻²", :IWP => "g m⁻²", :mean_T => "K", :mean_θ => "K", :mean_qv => "g kg⁻¹",
                   :Nc_cloudy => "cm⁻³", :rain_mm_day => "mm d⁻¹", :bottom_flux_mm_day => "mm d⁻¹", :max_w => "m s⁻¹")

old_dir, new_dir = ARGS[1], ARGS[2]
cases = sort!([splitext(f)[1] for f in readdir(new_dir) if endswith(f, ".csv")] ∪
              [splitext(f)[1] for f in readdir(old_dir) if endswith(f, ".csv")])

fmt(x) = isnan(x) ? "—" : abs(x) ≥ 1e4 || (abs(x) < 1e-3 && x != 0) ? @sprintf("%.3e", x) : @sprintf("%.4f", x)

for case in cases
    paths = joinpath.((old_dir, new_dir), case * ".csv")
    println("\n### ", case, "\n")
    if !all(isfile, paths)
        for (label, dir) in (("old", old_dir), ("new", new_dir))
            err = joinpath(dir, case * ".error")
            isfile(err) && println("$label failed: ", first(first(readlines(err)), 300))
        end
        continue
    end
    (h₀, r₀), (h₁, r₁) = read_csv.(paths)
    col(h, r, name) = (i = findfirst(==(name), h); isnothing(i) ? NaN : r[i])
    last₀, last₁ = r₀[end], r₁[end]
    it₀, it₁ = col(h₀, last₀, :iteration), col(h₁, last₁, :iteration)
    println("final record: old iteration ", Int(it₀), " (t = ", col(h₀, last₀, :t), " s), new iteration ", Int(it₁),
            " (t = ", col(h₁, last₁, :t), " s)\n")
    println("| quantity | old | new | new − old |")
    println("| --- | ---: | ---: | ---: |")
    for name in COLUMNS
        a, b = col(h₀, last₀, name), col(h₁, last₁, name)
        isnan(a) && isnan(b) && continue
        println("| ", name, " [", UNITS[name], "] | ", fmt(a), " | ", fmt(b), " | ", fmt(b - a), " |")
    end
    for (name, unit) in ((:water_prognostic, "kg m⁻²"), (:water, "kg m⁻²"), (:ρθ_column, "kg K m⁻²"))
        d₀ = col(h₀, last₀, name) - col(h₀, r₀[1], name)
        d₁ = col(h₁, last₁, name) - col(h₁, r₁[1], name)
        println("| drift of ", name, " [", unit, "] | ", fmt(d₀), " | ", fmt(d₁), " | ", fmt(d₁ - d₀), " |")
    end
end
