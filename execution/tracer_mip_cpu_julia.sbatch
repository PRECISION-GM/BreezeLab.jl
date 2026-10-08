#!/bin/bash
# Generic CPU-partition runner for a Julia script in the TRACER-MIP case environment.
#   sbatch --partition=cpu-c6a24xlarge --export=ALL,SCRIPT=analysis/x.jl,SCRIPT_ARGS="a b" execution/tracer_mip_cpu_julia.sbatch
#SBATCH --job-name=blab-tmip-julia
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=04:00:00
#SBATCH --output=/shared/home/greg/breezelab-work/tracer-mip/logs/slurm-%j.out
set -euo pipefail
CHECKOUT=${CHECKOUT:-$(cd "$(dirname "$0")/.." && pwd)}
W=${DEPOT_ROOT:-/shared/home/greg/breezelab-work/tracer-mip}
C=/shared/home/greg/breezelab-runs/20261006
D=/shared/home/greg/breezelab-runs/docs-20261007-8dc824d
CPU_TAG=$(awk -F ': ' '/model name/ {print $2; exit}' /proc/cpuinfo | tr -c 'A-Za-z0-9' '_')
source "$W/execution/depot_layout.sh"
export JULIA_DEPOT_PATH="$W/depot/$CPU_TAG:$D/depot-compute/$CPU_TAG:$D/depot-login:$C/depot:$HOME/.julia:"
export JULIA_PKG_PRECOMPILE_AUTO=0 JULIA_NUM_THREADS=8
cd "$CHECKOUT"
echo "== $(date -u) node=$(hostname) commit=$(git rev-parse HEAD) script=$SCRIPT args=${SCRIPT_ARGS:-}"
"$C/julia-1.12.7/bin/julia" --startup-file=no --project=cases/tracer_mip "$SCRIPT" ${SCRIPT_ARGS:-}
echo "== $(date -u) done"
