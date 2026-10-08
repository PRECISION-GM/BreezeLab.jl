#!/bin/bash
# CPU-partition job: animations of a TRACER-MIP outer run (analysis/animate_tracer_mip.jl) from a snapshot copy.
#   sbatch --partition=cpu-c6a24xlarge --export=ALL,RUN=<snapshot dir>,OUT=<animations dir> execution/tracer_mip_cpu_animate.sbatch
#SBATCH --job-name=blab-tmip-anim
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=06:00:00
#SBATCH --output=/shared/home/greg/breezelab-work/tracer-mip/logs/slurm-%j.out
set -euo pipefail
CHECKOUT=${CHECKOUT:-/shared/home/greg/breezelab-work/tracer-mip-anim}
W=${DEPOT_ROOT:-/shared/home/greg/breezelab-work/tracer-mip}      # compiled caches / shared sources
C=/shared/home/greg/breezelab-runs/20261006
D=/shared/home/greg/breezelab-runs/docs-20261007-8dc824d
CPU_TAG=$(awk -F ': ' '/model name/ {print $2; exit}' /proc/cpuinfo | tr -c 'A-Za-z0-9' '_')
source "$W/execution/depot_layout.sh"
export JULIA_DEPOT_PATH="$W/depot/$CPU_TAG:$D/depot-compute/$CPU_TAG:$D/depot-login:$C/depot:$HOME/.julia:"
export JULIA_PKG_PRECOMPILE_AUTO=0 JULIA_NUM_THREADS=8
cd "$CHECKOUT"
echo "== $(date -u) node=$(hostname) commit=$(git rev-parse HEAD) run=$RUN out=$OUT"
"$C/julia-1.12.7/bin/julia" --startup-file=no --project=cases/tracer_mip -e 'using Pkg; Pkg.precompile("BreezeLab")' 2>&1 | tail -1
"$C/julia-1.12.7/bin/julia" --startup-file=no --project=cases/tracer_mip analysis/animate_tracer_mip.jl "$RUN" "$OUT"
echo "== $(date -u) done"
