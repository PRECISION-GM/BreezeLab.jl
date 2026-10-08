#!/bin/bash
# CPU-partition job: stage ERA5 for one case (credentials from ~/.bashrc; never printed).
#   sbatch --partition=cpu-c6a24xlarge --export=ALL,CASE=aug07 execution/tracer_mip_cpu_fetch.sbatch
#SBATCH --job-name=blab-tmip-era5
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=/shared/home/greg/breezelab-work/tracer-mip/logs/slurm-%j.out
set -euo pipefail
W=/shared/home/greg/breezelab-work/tracer-mip
C=/shared/home/greg/breezelab-runs/20261006
D=/shared/home/greg/breezelab-runs/docs-20261007-8dc824d
CPU_TAG=$(awk -F ': ' '/model name/ {print $2; exit}' /proc/cpuinfo | tr -c 'A-Za-z0-9' '_')
source "$W/execution/depot_layout.sh"
export JULIA_DEPOT_PATH="$W/depot/$CPU_TAG:$D/depot-compute/$CPU_TAG:$D/depot-login:$C/depot:$HOME/.julia:"
export JULIA_PKG_PRECOMPILE_AUTO=0 JULIA_NUM_THREADS=4
set -a; source "$HOME/.bashrc" 2>/dev/null || true; set +a
cd "$W"
CASE=${CASE:-aug07}
echo "== $(date -u) node=$(hostname) commit=$(git rev-parse HEAD) case=$CASE"
"$C/julia-1.12.7/bin/julia" --startup-file=no --project=cases/tracer_mip data_wrangling/fetch_era5_tracer_mip.jl "$CASE" "$W/data/era5" 6 2>&1 | grep --line-buffered -v -i "key:\|token"
echo "== $(date -u) done"
