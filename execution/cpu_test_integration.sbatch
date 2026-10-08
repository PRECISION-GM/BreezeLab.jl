#!/bin/bash
#SBATCH --job-name=blab-integ-test
#SBATCH --partition=cpu-c6a24xlarge
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=16 --mem=64G --time=03:00:00
#SBATCH --output=/shared/home/greg/breezelab-work/integration/logs/slurm-%j.out
set -euo pipefail
W=/shared/home/greg/breezelab-work/integration; C=/shared/home/greg/breezelab-runs/20261006; D=/shared/home/greg/breezelab-runs/docs-20261007-8dc824d
CPU_TAG=$(awk -F ': ' '/model name/ {print $2; exit}' /proc/cpuinfo | tr -c 'A-Za-z0-9' '_')
mkdir -p "$W/depot/$CPU_TAG"
export JULIA_DEPOT_PATH="$W/depot/$CPU_TAG:/shared/home/greg/breezelab-work/tracer-mip/depot/$CPU_TAG:/shared/home/greg/breezelab-work/tracer-dp-scream/depot/$CPU_TAG:/shared/home/greg/breezelab-work/sea-starr/depot/$CPU_TAG:/shared/home/greg/breezelab-work/ena-lasso/depot/$CPU_TAG:$D/depot-compute/$CPU_TAG:$C/depot:$HOME/.julia:"
export JULIA_PKG_PRECOMPILE_AUTO=0 JULIA_NUM_THREADS=8 JULIA_NUM_PRECOMPILE_TASKS=8
cd $W; echo "commit=$(git rev-parse HEAD) node=$(hostname) started=$(date -u +%FT%TZ)"
"$C/julia-1.12.7/bin/julia" --startup-file=no --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()' 2>&1 | tail -3
time "$C/julia-1.12.7/bin/julia" --startup-file=no --project=. -e 'using Pkg; Pkg.test(; allow_reresolve=false)'
echo "finished=$(date -u +%FT%TZ)"
