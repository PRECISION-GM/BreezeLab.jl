#!/bin/bash
# TRACER-MIP outer-domain launcher (wpcluster), following the shared UUID-pinning procedure:
#   sbatch --partition=gpu-p4de-2c [--nodelist=<node>] --gres=gpu:1 --cpus-per-task=12 --mem=96G --time=HH:MM:SS \
#          --export=ALL,EXPERIMENT=<name>,TRACER_MIP_NX=..,TRACER_MIP_NY=..,TRACER_MIP_STOP_HOURS=..,TRACER_MIP_PARENT=era5|synthetic \
#          execution/tracer_mip_outer.sbatch <GPU-UUID | auto>
# With a UUID, verify it is physically idle (ssh <node> nvidia-smi) immediately before submitting. With `auto` the job
# picks, at start, the first A100/H100 on its node (Slurm's own allocation first) that has no compute processes and
# whose lock in the shared lock directory it can take — for when UUIDs are unknown (e.g. replaced nodes).
# Requeue-safe: on a Slurm requeue (SLURM_RESTART_COUNT > 0) the job writes into <run>/segment_<n> and picks up from the
# latest checkpoint in <run>/checkpoints, so an outage never overwrites earlier output.
#SBATCH --job-name=blab-tracer-mip
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=12
#SBATCH --mem=96G
#SBATCH --output=/shared/home/greg/breezelab-work/tracer-mip/logs/slurm-%j.out
set -euo pipefail
W=/shared/home/greg/breezelab-work/tracer-mip          # depot root (compiled caches, shared sources)
CHECKOUT=${CHECKOUT:-$W}                                 # source tree to run (a worktree of another branch, if given)
C=/shared/home/greg/breezelab-runs/20261006
D=/shared/home/greg/breezelab-runs/docs-20261007-8dc824d
GPU_UUID=${1:?GPU UUID or auto}
LOCKDIR=/shared/home/greg/breezelab-work/gpu-locks
mkdir -p "$LOCKDIR" "$W/logs" "$W/runs"
claim_gpu() {   # $1 = UUID; succeeds only if the shared lock is free and the GPU has no compute processes
    exec 9>"$LOCKDIR/$1.lock"
    if flock -n 9 && [ -z "$(nvidia-smi -i "$1" --query-compute-apps=pid --format=csv,noheader)" ]; then return 0; fi
    exec 9>&-; return 1
}
if [ "$GPU_UUID" = auto ]; then
    SLURM_IDX=${CUDA_VISIBLE_DEVICES:-}; SLURM_IDX=${SLURM_IDX%%,*}
    SLURM_UUID=$( [ -n "$SLURM_IDX" ] && nvidia-smi -i "$SLURM_IDX" --query-gpu=uuid --format=csv,noheader 2>/dev/null || true)
    CANDIDATES=$( { [ -n "$SLURM_UUID" ] && echo "$SLURM_UUID"; nvidia-smi --query-gpu=uuid,name --format=csv,noheader | grep -E "A100|H100" | cut -d, -f1 || true; } | awk '!seen[$0]++' || true)
    GPU_UUID=""
    for U in $CANDIDATES; do claim_gpu "$U" && { GPU_UUID=$U; break; }; done
    [ -n "$GPU_UUID" ] || { echo "No physically idle A100/H100 on $(hostname)"; exit 1; }
    echo "Selected GPU $GPU_UUID on $(hostname) (Slurm allocation: ${SLURM_UUID:-none})"
else
    claim_gpu "$GPU_UUID" || { echo "GPU $GPU_UUID is locked or busy"; exit 1; }
fi
export CUDA_VISIBLE_DEVICES="$GPU_UUID"
GPU_NAME=$(nvidia-smi -i "$CUDA_VISIBLE_DEVICES" --query-gpu=name --format=csv,noheader)
case "$GPU_NAME" in *H100*|*A100*) ;; *) echo "Unexpected GPU: $GPU_NAME"; exit 1 ;; esac

CPU_TAG=$(awk -F ': ' '/model name/ {print $2; exit}' /proc/cpuinfo | tr -c 'A-Za-z0-9' '_')
mkdir -p "$W/depot/$CPU_TAG"
source "$W/execution/depot_layout.sh"
export JULIA_DEPOT_PATH="$W/depot/$CPU_TAG:$D/depot-compute/$CPU_TAG:$D/depot-login:$C/depot:$HOME/.julia:"
export JULIA_PKG_PRECOMPILE_AUTO=0
export JULIA_NUM_PRECOMPILE_TASKS=6
export JULIA_NUM_THREADS=6
# CDS credentials for any in-run ERA5 access (values are never printed or logged).
set -a; source "$HOME/.bashrc" 2>/dev/null || true; set +a

EXPERIMENT=${EXPERIMENT:-outer_pilot}
BASE="$W/runs/${EXPERIMENT}_job${SLURM_JOB_ID}"
RESTARTS=${SLURM_RESTART_COUNT:-0}
if [ "$RESTARTS" -gt 0 ]; then
    RUN="$BASE/segment_$RESTARTS"
    export TRACER_MIP_PICKUP=${TRACER_MIP_PICKUP:-latest}
else
    RUN="$BASE"
fi
mkdir -p "$RUN"
export TRACER_MIP_OUTPUT_DIR="$RUN"
export TRACER_MIP_CHECKPOINT_DIR=${TRACER_MIP_CHECKPOINT_DIR:-$BASE/checkpoints}
export TRACER_MIP_ERA5_DIR="${TRACER_MIP_ERA5_DIR:-$W/data/era5}"
cd "$CHECKOUT"
git rev-parse HEAD > "$RUN/COMMIT"
git status --porcelain > "$RUN/GIT_STATUS" || true
cp Manifest.toml "$RUN/Manifest.package.toml"
cp cases/tracer_mip/Manifest.toml "$RUN/Manifest.case.toml"
printf 'job=%s restart=%s gpu=%s uuid=%s node=%s experiment=%s\nNX=%s NY=%s STOP_HOURS=%s PARENT=%s CLOSURE=%s CFL=%s RELAX_WIDTH=%s SPONGE_DEPTH=%s FAST_MINUTES=%s CHECKPOINT_HOURS=%s PICKUP=%s\n' \
  "$SLURM_JOB_ID" "$RESTARTS" "$GPU_NAME" "$GPU_UUID" "$(hostname)" "$EXPERIMENT" \
  "${TRACER_MIP_NX:-protocol}" "${TRACER_MIP_NY:-protocol}" "${TRACER_MIP_STOP_HOURS:-protocol}" "${TRACER_MIP_PARENT:-era5}" \
  "${TRACER_MIP_CLOSURE:-tke}" "${TRACER_MIP_CFL:-0}" "${TRACER_MIP_RELAX_WIDTH:-5}" "${TRACER_MIP_SPONGE_DEPTH:-0}" \
  "${TRACER_MIP_FAST_MINUTES:-2.5}" "${TRACER_MIP_CHECKPOINT_HOURS:-3}" "${TRACER_MIP_PICKUP:-none}" > "$RUN/COMMAND"
( while true; do printf '%s,' "$(date -u +%FT%TZ)"; nvidia-smi -i "$GPU_UUID" --query-gpu=memory.used,memory.total,utilization.gpu --format=csv,noheader; sleep 30; done ) >> "$RUN/gpu_memory.csv" &
MONITOR=$!
trap 'kill $MONITOR 2>/dev/null || true' EXIT
echo "Starting TRACER-MIP $EXPERIMENT: job=$SLURM_JOB_ID gpu=$GPU_NAME commit=$(cat "$RUN/COMMIT")"
"$C/julia-1.12.7/bin/julia" --startup-file=no --project=cases/tracer_mip cases/tracer_mip.jl
touch "$RUN/COMPLETE"
