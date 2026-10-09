#!/bin/bash
#SBATCH --job-name=blab-seastarr
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=12
#SBATCH --mem=96G
#SBATCH --output=/shared/home/greg/breezelab-work/sea-starr/logs/slurm-%j.out
# SEA STARR launcher that selects a physically idle GPU at start time (COORDINATION.md procedure
# without a pre-chosen UUID): let Slurm place the job, then inside the allocation take the first
# A100/H100 on the node that has no compute process and whose shared lock is free, pin it, and run.
#   sbatch --partition=gpu-p4de-2c --time=<wall> execution/sea_starr_autogpu.sbatch <experiment> [hours] [Nx]
# Environment: SEA_STARR_MEMBER (CTRL), SEA_STARR_CHECKPOINT_HOURS (3), SEA_STARR_PICKUP (false).
set -euo pipefail
W=/shared/home/greg/breezelab-work/sea-starr
C=/shared/home/greg/breezelab-runs/20261006
EXPERIMENT=${1:?experiment name}
HOURS=${2:-66}
NX=${3:-192}
LOCKDIR=/shared/home/greg/breezelab-work/gpu-locks
mkdir -p "$LOCKDIR" "$W/logs" "$W/runs"

# Candidate GPUs: the Slurm-allocated one first (CUDA_VISIBLE_DEVICES may be an index or UUID), then every other GPU on the node.
allocated=${CUDA_VISIBLE_DEVICES:-}
mapfile -t ALL < <(nvidia-smi --query-gpu=uuid,name --format=csv,noheader | sed 's/, /,/')
candidates=()
for entry in "${ALL[@]}"; do
    uuid=${entry%%,*}; name=${entry#*,}
    case "$name" in *H100*|*A100*) ;; *) continue ;; esac
    idx=$(nvidia-smi --query-gpu=index --format=csv,noheader -i "$uuid")
    if [[ -n "$allocated" && ( "$allocated" == "$uuid" || "$allocated" == "$idx" ) ]]; then
        candidates=("$uuid" "${candidates[@]}")
    else
        candidates+=("$uuid")
    fi
done
GPU_UUID=""
for uuid in "${candidates[@]}"; do
    [ -z "$(nvidia-smi -i "$uuid" --query-compute-apps=pid --format=csv,noheader)" ] || { echo "busy: $uuid"; continue; }
    exec 9>"$LOCKDIR/$uuid.lock"
    if flock -n 9; then GPU_UUID=$uuid; break; else echo "locked: $uuid"; exec 9>&-; fi
done
[ -n "$GPU_UUID" ] || { echo "no physically idle, unlocked A100/H100 on $(hostname) at $(date -u +%FT%TZ)"; exit 1; }
export CUDA_VISIBLE_DEVICES="$GPU_UUID"
GPU_NAME=$(nvidia-smi -i "$GPU_UUID" --query-gpu=name --format=csv,noheader)

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
  echo "job=$SLURM_JOB_ID node=$(hostname) gpu=$GPU_NAME uuid=$GPU_UUID (auto-selected; slurm allocation: ${allocated:-none})"
  echo "commit=$(git rev-parse HEAD) branch=$(git rev-parse --abbrev-ref HEAD) dirty=$(git status --porcelain | wc -l)"
  echo "command: sbatch ... execution/sea_starr_autogpu.sbatch $EXPERIMENT $HOURS $NX (member=$SEA_STARR_MEMBER)"
  echo "started=$(date -u +%FT%TZ)"
  sha256sum "$SEA_STARR_DATA"/SEA_STARR_${SEA_STARR_MEMBER}_SCM_driver.nc
  cp Manifest.toml "$RUN/Manifest.toml"
} | tee "$RUN/job.txt"

"$C/julia-1.12.7/bin/julia" --startup-file=no --project=. -e 'using Pkg; Pkg.precompile()' 2>&1 | tail -3
"$C/julia-1.12.7/bin/julia" --startup-file=no --project=. cases/sea_starr.jl 2>&1 | tee "$RUN/run.log"
echo "finished=$(date -u +%FT%TZ)" >> "$RUN/job.txt"
touch "$RUN/COMPLETE"
