#!/bin/bash
# TRACER-MIP outer-domain launcher (wpcluster), following the shared UUID-pinning procedure:
#   sbatch --partition=gpu-p4de-2c --nodelist=<node> --gres=gpu:1 --cpus-per-task=12 --mem=96G --time=HH:MM:SS \
#          --export=ALL,EXPERIMENT=<name>,TRACER_MIP_NX=..,TRACER_MIP_NY=..,TRACER_MIP_STOP_HOURS=..,TRACER_MIP_PARENT=era5|synthetic \
#          execution/tracer_mip_outer.sbatch <GPU-UUID>
# Verify the GPU is physically idle (ssh <node> nvidia-smi) immediately before submitting.
#SBATCH --job-name=blab-tracer-mip
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=12
#SBATCH --mem=96G
#SBATCH --output=/shared/home/greg/breezelab-work/tracer-mip/logs/slurm-%j.out
set -euo pipefail
W=/shared/home/greg/breezelab-work/tracer-mip
C=/shared/home/greg/breezelab-runs/20261006
D=/shared/home/greg/breezelab-runs/docs-20261007-8dc824d
GPU_UUID=${1:?GPU UUID}
LOCKDIR=/shared/home/greg/breezelab-work/gpu-locks
mkdir -p "$LOCKDIR" "$W/logs" "$W/runs"
exec 9>"$LOCKDIR/$GPU_UUID.lock"
flock -n 9 || { echo "GPU lock $GPU_UUID is held by another job"; exit 1; }
[ -z "$(nvidia-smi -i "$GPU_UUID" --query-compute-apps=pid --format=csv,noheader)" ] || { echo 'Selected GPU is busy'; exit 1; }
export CUDA_VISIBLE_DEVICES="$GPU_UUID"
GPU_NAME=$(nvidia-smi -i "$CUDA_VISIBLE_DEVICES" --query-gpu=name --format=csv,noheader)
case "$GPU_NAME" in *H100*|*A100*) ;; *) echo "Unexpected GPU: $GPU_NAME"; exit 1 ;; esac

CPU_TAG=$(awk -F ': ' '/model name/ {print $2; exit}' /proc/cpuinfo | tr -c 'A-Za-z0-9' '_')
mkdir -p "$W/depot/$CPU_TAG"
export JULIA_DEPOT_PATH="$W/depot/$CPU_TAG:$D/depot-compute/$CPU_TAG:$D/depot-login:$C/depot:$HOME/.julia:"
export JULIA_PKG_PRECOMPILE_AUTO=0
export JULIA_NUM_PRECOMPILE_TASKS=6
export JULIA_NUM_THREADS=6
# CDS credentials for any in-run ERA5 access (values are never printed or logged).
set -a; source "$HOME/.bashrc" 2>/dev/null || true; set +a

EXPERIMENT=${EXPERIMENT:-outer_pilot}
RUN="$W/runs/${EXPERIMENT}_job${SLURM_JOB_ID}"
mkdir -p "$RUN"
export TRACER_MIP_OUTPUT_DIR="$RUN"
export TRACER_MIP_ERA5_DIR="${TRACER_MIP_ERA5_DIR:-$W/data/era5}"
cd "$W"
git rev-parse HEAD > "$RUN/COMMIT"
git status --porcelain > "$RUN/GIT_STATUS" || true
cp Manifest.toml "$RUN/Manifest.package.toml"
cp cases/tracer_mip/Manifest.toml "$RUN/Manifest.case.toml"
printf 'job=%s gpu=%s uuid=%s node=%s experiment=%s\nNX=%s NY=%s STOP_HOURS=%s PARENT=%s\n' \
  "$SLURM_JOB_ID" "$GPU_NAME" "$GPU_UUID" "$(hostname)" "$EXPERIMENT" \
  "${TRACER_MIP_NX:-protocol}" "${TRACER_MIP_NY:-protocol}" "${TRACER_MIP_STOP_HOURS:-protocol}" "${TRACER_MIP_PARENT:-era5}" > "$RUN/COMMAND"
( while true; do printf '%s,' "$(date -u +%FT%TZ)"; nvidia-smi -i "$GPU_UUID" --query-gpu=memory.used,memory.total,utilization.gpu --format=csv,noheader; sleep 30; done ) >> "$RUN/gpu_memory.csv" &
MONITOR=$!
trap 'kill $MONITOR 2>/dev/null || true' EXIT
echo "Starting TRACER-MIP $EXPERIMENT: job=$SLURM_JOB_ID gpu=$GPU_NAME commit=$(cat "$RUN/COMMIT")"
"$C/julia-1.12.7/bin/julia" --startup-file=no --project=cases/tracer_mip cases/tracer_mip.jl
touch "$RUN/COMPLETE"
