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
cp(joinpath(@__DIR__, "cases", "ena.md"), joinpath(generated, "ena_protocols.md"); force=true)
cp(joinpath(@__DIR__, "cases", "ena_lasso.md"), joinpath(generated, "ena_lasso.md"); force=true)
# cases/ena_lasso.jl is not executed here: it needs the staged ARM bundle and a GPU.
cp(joinpath(@__DIR__, "cases", "tracer_mip.md"), joinpath(generated, "tracer_mip.md"); force=true)

# Add further runnable case scripts here as their implementations land. SEA STARR
# and TRACER are currently specifications, so they are not executable examples yet.
Literate.markdown(joinpath(dirname(@__DIR__), "cases", "eastern_north_atlantic.jl"), generated;
                 flavor=Literate.DocumenterFlavor(), execute=true,
                 preprocess=prepare_example)

makedocs(; root=@__DIR__,
         modules=[BreezeLab],
         sitename=smoke ? "BreezeLab — CPU smoke build" : "BreezeLab",
         format=Documenter.HTML(prettyurls=false, edit_link="main"),
         checkdocs=:none, # The initial manual documents selected case APIs.
         pages=["Home" => "index.md",
                "ENA example" => "generated/eastern_north_atlantic.md",
                "ENA protocols" => "generated/ena_protocols.md",
                "ENA-LASSO audit" => "generated/ena_lasso.md",
                "TRACER-MIP" => "generated/tracer_mip.md",
                "Case API" => "api.md"])
