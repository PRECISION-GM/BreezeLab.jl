#!/bin/bash
#SBATCH --job-name=blab-seastarr
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=12
#SBATCH --mem=96G
#SBATCH --output=/shared/home/greg/breezelab-work/sea-starr/logs/slurm-%j.out
# SEA STARR launcher following the wpcluster UUID-pinning procedure (COORDINATION.md):
#   sbatch --partition=<p> --nodelist=<node> --time=<wall> execution/sea_starr.sbatch <GPU_UUID> <experiment> [hours] [Nx]
# The job takes the shared lock for the GPU, refuses a GPU with running compute processes,
# pins CUDA_VISIBLE_DEVICES to it, accepts only H100/A100, and runs cases/sea_starr.jl with
# the per-CPU compiled-cache overlay. Outputs go to runs/<experiment>_job<id>/; COMPLETE is
# written only on exit 0.
set -euo pipefail
W=/shared/home/greg/breezelab-work/sea-starr
C=/shared/home/greg/breezelab-runs/20261006
GPU_UUID=${1:?GPU UUID}
EXPERIMENT=${2:?experiment name}
HOURS=${3:-66}
NX=${4:-192}
LOCKDIR=/shared/home/greg/breezelab-work/gpu-locks
mkdir -p "$LOCKDIR" "$W/logs" "$W/runs"
exec 9>"$LOCKDIR/$GPU_UUID.lock"
flock -n 9 || { echo "GPU $GPU_UUID is locked by another job"; exit 1; }
[ -z "$(nvidia-smi -i "$GPU_UUID" --query-compute-apps=pid --format=csv,noheader)" ] || { echo 'Selected GPU is busy'; exit 1; }
export CUDA_VISIBLE_DEVICES="$GPU_UUID"
GPU_NAME=$(nvidia-smi -i "$CUDA_VISIBLE_DEVICES" --query-gpu=name --format=csv,noheader)
case "$GPU_NAME" in *H100*|*A100*) ;; *) echo "Unexpected GPU: $GPU_NAME"; exit 1 ;; esac

RUN="$W/runs/${EXPERIMENT}_job${SLURM_JOB_ID}"
mkdir -p "$RUN"
trap 'code=$?; if (( code != 0 )); then printf "exit=%s job=%s\n" "$code" "$SLURM_JOB_ID" > "$RUN/FAILED"; fi' EXIT

CPU_TAG=$(awk -F ': ' '/model name/ {print $2; exit}' /proc/cpuinfo | tr -c 'A-Za-z0-9' '_')
mkdir -p "$W/depot/$CPU_TAG"
export JULIA_DEPOT_PATH="$W/depot/$CPU_TAG:/shared/home/greg/breezelab-runs/docs-20261007-8dc824d/depot-compute/$CPU_TAG:/shared/home/greg/breezelab-runs/docs-20261007-8dc824d/depot-login:$C/depot:$HOME/.julia:"
export JULIA_PKG_PRECOMPILE_AUTO=0
export JULIA_NUM_PRECOMPILE_TASKS=6
export JULIA_NUM_THREADS=6
export JULIA_CUDA_USE_COMPAT=false
export SEA_STARR_HOURS="$HOURS"
export SEA_STARR_NX="$NX"
export SEA_STARR_OUTPUT="$RUN"
export SEA_STARR_DATA="$C/inputs/mip_sources/seastarr_22241697"
export SEA_STARR_MEMBER="${SEA_STARR_MEMBER:-CTRL}"
export SEA_STARR_CHECKPOINT_HOURS="${SEA_STARR_CHECKPOINT_HOURS:-3}"

cd "$W"
{
  echo "job=$SLURM_JOB_ID node=$(hostname) gpu=$GPU_NAME uuid=$GPU_UUID"
  echo "commit=$(git rev-parse HEAD) branch=$(git rev-parse --abbrev-ref HEAD) dirty=$(git status --porcelain | wc -l)"
  echo "command: sbatch ... execution/sea_starr.sbatch $GPU_UUID $EXPERIMENT $HOURS $NX (member=$SEA_STARR_MEMBER)"
  echo "started=$(date -u +%FT%TZ)"
  sha256sum "$SEA_STARR_DATA"/SEA_STARR_${SEA_STARR_MEMBER}_SCM_driver.nc
  cp Manifest.toml "$RUN/Manifest.toml"
} | tee "$RUN/job.txt"

"$C/julia-1.12.7/bin/julia" --startup-file=no --project=. -e 'using Pkg; Pkg.precompile()' 2>&1 | tail -3
"$C/julia-1.12.7/bin/julia" --startup-file=no --project=. cases/sea_starr.jl
echo "finished=$(date -u +%FT%TZ)" >> "$RUN/job.txt"
touch "$RUN/COMPLETE"
