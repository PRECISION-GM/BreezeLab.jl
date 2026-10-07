# Derive a docs environment from the pinned simulation environment.
# Re-running this resets docs/Manifest.toml to those pins before adding docs tools.
using Pkg

root = dirname(@__DIR__)
cp(joinpath(root, "Manifest.toml"), joinpath(@__DIR__, "Manifest.toml"); force=true)
Pkg.activate(@__DIR__)
cd(@__DIR__) do
    Pkg.develop(PackageSpec(path=".."); preserve=Pkg.PRESERVE_ALL)
end
Pkg.instantiate()
