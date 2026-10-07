# Manual documentation builds

The instructions and resource requirements are in [the manual home page](src/index.md#build-this-manual).
`make.jl` executes the ENA case with Literate (including analysis and a figure), then
builds local HTML with Documenter. No GPU CI or deployment is configured.

From the repository root, after fetching the Covert inputs:

```sh
julia docs/setup.jl
julia --project=docs docs/make.jl          # full GPU example, H100 or A100
julia --project=docs docs/make.jl --smoke  # optional four-second CPU workflow check
```

The full and smoke commands are alternative builds; each executes the example.
Open `docs/build/index.html` after the chosen build finishes.
