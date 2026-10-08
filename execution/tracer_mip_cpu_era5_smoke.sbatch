#!/bin/bash
# CPU-partition smoke of the real ERA5-driven outer domain at reduced size (16x16x94, 30 s):
# exercises the staged ERA5 parent, ETOPO terrain (downloaded into the shared depot scratch),
# land initialization from ERA5/ERA5-Land, sea pinning and RRTMGP, before any GPU run.
#SBATCH --job-name=blab-tmip-era5smoke
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=/shared/home/greg/breezelab-work/tracer-mip/logs/slurm-%j.out
set -euo pipefail
W=/shared/home/greg/breezelab-work/tracer-mip
C=/shared/home/greg/breezelab-runs/20261006
D=/shared/home/greg/breezelab-runs/docs-20261007-8dc824d
CPU_TAG=$(awk -F ': ' '/model name/ {print $2; exit}' /proc/cpuinfo | tr -c 'A-Za-z0-9' '_')
source "$W/execution/depot_layout.sh"
export JULIA_DEPOT_PATH="$W/depot/$CPU_TAG:$D/depot-compute/$CPU_TAG:$D/depot-login:$C/depot:$HOME/.julia:"
export JULIA_PKG_PRECOMPILE_AUTO=0 JULIA_NUM_PRECOMPILE_TASKS=8 JULIA_NUM_THREADS=8
set -a; source "$HOME/.bashrc" 2>/dev/null || true; set +a
cd "$W"
echo "== $(date -u) node=$(hostname) commit=$(git rev-parse HEAD)"
export TRACER_MIP_ARCH=cpu TRACER_MIP_NX=${TRACER_MIP_NX:-16} TRACER_MIP_NY=${TRACER_MIP_NY:-16} TRACER_MIP_STOP_HOURS=${TRACER_MIP_STOP_HOURS:-0.01} TRACER_MIP_PARENT=era5
export TRACER_MIP_ERA5_DIR="$W/data/era5" TRACER_MIP_OUTPUT_DIR="$W/runs/era5_cpu_smoke_job${SLURM_JOB_ID}"
mkdir -p "$TRACER_MIP_OUTPUT_DIR"
"$C/julia-1.12.7/bin/julia" --startup-file=no --project=cases/tracer_mip -e 'using Pkg; Pkg.precompile("BreezeLab")' 2>&1 | tail -2
"$C/julia-1.12.7/bin/julia" --startup-file=no --project=cases/tracer_mip cases/tracer_mip.jl 2>&1 | grep --line-buffered -v -i "key:\|token"
touch "$TRACER_MIP_OUTPUT_DIR/COMPLETE"
echo "== $(date -u) done"
