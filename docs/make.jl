using BreezeLab
using CairoMakie
using Documenter
using Literate

all(==("--smoke"), ARGS) || error("Usage: julia --project=docs docs/make.jl [--smoke]")
smoke = "--smoke" in ARGS
CairoMakie.activate!(type="png")

# Execute the same top-level script shown in the manual. The smoke build changes
# explicit settings before Literate runs, so its displayed code matches its results.
function prepare_example(content)
    run_name = smoke ? "ena_smoke" : "ena"
    content = replace(content,
        "output_dir = joinpath(pkgdir(BreezeLab), \"output\", \"eastern_north_atlantic\")" =>
        "output_dir = joinpath(pkgdir(BreezeLab), \"output\", \"documentation\", $(repr(run_name)))")
    if smoke
        replacements = [
            "arch = GPU()" => "arch = CPU()",
            "Nx, Ny = 256, 256" => "Nx, Ny = 8, 8",
            "z_faces = covert_public_bin_vertical_faces()" =>
                "z_faces = collect(range(0, 6000, length=25))",
            "stop_time = 6hours" => "stop_time = 4seconds",
            "timeseries_interval = 60seconds" => "timeseries_interval = 1second",
            "profile_interval = 1hour" => "profile_interval = 4seconds",
            "slice_interval = 1hour" => "slice_interval = 4seconds",
        ]
        for (before, after) in replacements
            occursin(before, content) || error("Smoke override no longer matches the example: $before")
            content = replace(content, before => after)
        end
        content = replace(content, "# # Eastern North Atlantic" => """
# # Eastern North Atlantic
#
# !!! warning "CPU documentation smoke build"
#     This page executes an 8×8×24, four-second CPU variant. It checks the script,
#     output reading, plotting, and documentation, not the six-hour cloud evolution
#     described below. The code shown includes the reduced settings.
"""; count=1)
    end
    return content * """

# ## Execution environment

using InteractiveUtils: versioninfo
versioninfo()
import Pkg
Pkg.status()
"""
end

generated = joinpath(@__DIR__, "src", "generated")
mkpath(generated)

# Static case pages written alongside the implementations (audits, inventories, analyses).
static_pages = [
    "ena.md" => "ena_protocols.md",
    "ena_lasso.md" => "ena_lasso.md",
    "ena_covert_analysis.md" => "ena_covert_analysis.md",
    "ena_aerosol_audit.md" => "ena_aerosol_audit.md",
    "sea_starr.md" => "sea_starr.md",
    "tracer_dp_scream.md" => "tracer_dp_scream.md",
    "tracer_mip.md" => "tracer_mip.md",
    "campaign_results.md" => "campaign_results.md",
]
for (source, target) in static_pages
    cp(joinpath(@__DIR__, "cases", source), joinpath(generated, target); force=true)
end

# Campaign figures (Julia/Makie output of the October 2026 runs) are not kept in Git.
# Point `BREEZELAB_CAMPAIGN_FIGURES` at a directory of PNGs to include them; otherwise the
# results page keeps its text and the figure links are left unresolved.
figures_dir = get(ENV, "BREEZELAB_CAMPAIGN_FIGURES", "/shared/home/greg/breezelab-work/campaign-figures")
if isdir(figures_dir)
    target = joinpath(generated, "figures")
    mkpath(target)
    for name in filter(endswith(".png"), readdir(figures_dir))
        cp(joinpath(figures_dir, name), joinpath(target, name); force=true)
    end
    @info "Copied $(count(endswith(".png"), readdir(figures_dir))) campaign figures from $figures_dir"
else
    @warn "Campaign figure directory not found; the results page will have broken figure links" figures_dir
end

# The ENA Covert example is executed (GPU in the full build, a CPU smoke variant otherwise).
Literate.markdown(joinpath(dirname(@__DIR__), "cases", "eastern_north_atlantic.jl"), generated;
                 flavor=Literate.DocumenterFlavor(), execute=true,
                 preprocess=prepare_example)

# The other case scripts are rendered without execution: each needs staged inputs
# (ARM LASSO bundle, DEPHY drivers, the E3SM IOP file, ERA5) and hours on a GPU.
# Their runs are launched with the Slurm scripts under execution/.
for script in ("ena_lasso.jl", "sea_starr.jl", "tracer_dp_scream.jl", "tracer_mip.jl")
    Literate.markdown(joinpath(dirname(@__DIR__), "cases", script), generated;
                     name=replace(script, ".jl" => "_script"),
                     flavor=Literate.DocumenterFlavor(), execute=false,
                     preprocess=content -> content * """

# !!! note "Rendered, not executed"
#     This page shows the case script as committed. It is executed on the cluster with the
#     corresponding `execution/*.sbatch` launcher; results are summarized on the
#     [campaign results](campaign_results.md) page.
""")
end

makedocs(; root=@__DIR__,
         modules=[BreezeLab],
         sitename=smoke ? "BreezeLab — CPU smoke build" : "BreezeLab",
         format=Documenter.HTML(prettyurls=false, edit_link="main"),
         checkdocs=:none, # The initial manual documents selected case APIs.
         pages=["Home" => "index.md",
                "Campaign results" => "generated/campaign_results.md",
                "Eastern North Atlantic" => [
                    "ENA example (executed)" => "generated/eastern_north_atlantic.md",
                    "ENA protocols" => "generated/ena_protocols.md",
                    "ENA-covert analysis" => "generated/ena_covert_analysis.md",
                    "ENA aerosol audit" => "generated/ena_aerosol_audit.md",
                    "ENA-LASSO member" => "generated/ena_lasso.md",
                    "ENA-LASSO case script" => "generated/ena_lasso_script.md",
                ],
                "SEA STARR" => [
                    "Driver inventory and mapping" => "generated/sea_starr.md",
                    "Case script" => "generated/sea_starr_script.md",
                ],
                "TRACER–DP-SCREAM" => [
                    "Case description" => "generated/tracer_dp_scream.md",
                    "Case script" => "generated/tracer_dp_scream_script.md",
                ],
                "TRACER-MIP" => [
                    "Case description" => "generated/tracer_mip.md",
                    "Case script" => "generated/tracer_mip_script.md",
                ],
                "Case API" => "api.md"])
