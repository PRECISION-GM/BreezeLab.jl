# Download pinned public Covert et al. (2022) inputs and record provenance.
using SHA

const REPOS = [
    (name = "covert2022_bulk", url = "https://github.com/dmechem/ENA_variability_LES_bulk_paper.git",
     rev = "a403498ff44045b531af5dbfe2ee461a9c878dd9",
     files = ["snd" => "snd", "lsf_time_varying" => "lsf", "sfc_time_varying" => "sfc"]),
    (name = "covert2022_bin", url = "https://github.com/dmechem/ENA_variability_LES_bin_paper.git",
     rev = "07958e6e53543763a2e2d8faf1ea09433c1086c0",
     files = ["snd" => "snd", "lsf" => "lsf", "sfc" => "sfc", "prm" => "prm"]),
]

data_dir = joinpath(@__DIR__, "..", "data")
mkpath(data_dir)
# Assemble the manifest in memory; only publish it after all downloads succeed.
manifest = IOBuffer()
for repo in REPOS
    mktempdir() do tmp
        run(`git clone -q $(repo.url) $tmp`)
        run(`git -C $tmp checkout -q --detach $(repo.rev)`)
        commit = strip(read(`git -C $tmp rev-parse HEAD`, String))
        commit == repo.rev || error("input revision mismatch for $(repo.name)")
        target = joinpath(data_dir, repo.name)
        mkpath(target)
        println(manifest, "\n[", repo.name, "]\nurl = \"", repo.url, "\"\ncommit = \"", commit, "\"")
        for (src, dst) in repo.files
            cp(joinpath(tmp, src), joinpath(target, dst); force=true)
            digest = bytes2hex(open(sha256, joinpath(target, dst)))
            println(manifest, dst, " = \"", src, "\"  # sha256 ", digest)
        end
    end
end
write(joinpath(data_dir, "MANIFEST.txt"), take!(manifest))
println("inputs written to ", abspath(data_dir))
