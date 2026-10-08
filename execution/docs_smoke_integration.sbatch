#!/bin/bash
#SBATCH --job-name=blab-docs-smoke
#SBATCH --partition=cpu-c6a24xlarge
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=16 --mem=64G --time=03:00:00
#SBATCH --output=/shared/home/greg/breezelab-work/integration/logs/slurm-%j.out
set -euo pipefail
W=/shared/home/greg/breezelab-work/integration; C=/shared/home/greg/breezelab-runs/20261006; D=/shared/home/greg/breezelab-runs/docs-20261007-8dc824d
CPU_TAG=$(awk -F ': ' '/model name/ {print $2; exit}' /proc/cpuinfo | tr -c 'A-Za-z0-9' '_')
mkdir -p "$W/depot/$CPU_TAG"
export JULIA_DEPOT_PATH="$W/depot/$CPU_TAG:/shared/home/greg/breezelab-work/tracer-mip/depot/$CPU_TAG:/shared/home/greg/breezelab-work/tracer-dp-scream/depot/$CPU_TAG:$D/depot-compute/$CPU_TAG:$C/depot:$HOME/.julia:"
export JULIA_PKG_PRECOMPILE_AUTO=0 JULIA_NUM_THREADS=8 JULIA_NUM_PRECOMPILE_TASKS=8
cd $W; echo "commit=$(git rev-parse HEAD) node=$(hostname) started=$(date -u +%FT%TZ)"
[ -d data/covert2022_bin ] || "$C/julia-1.12.7/bin/julia" --startup-file=no --project=. data_wrangling/fetch_covert_inputs.jl 2>&1 | tail -3
"$C/julia-1.12.7/bin/julia" --startup-file=no docs/setup.jl 2>&1 | tail -3
time "$C/julia-1.12.7/bin/julia" --startup-file=no --project=docs docs/make.jl --smoke
ls docs/build | head; ls docs/build/generated | head -30
echo "finished=$(date -u +%FT%TZ)"
