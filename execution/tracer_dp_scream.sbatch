#!/bin/bash
#SBATCH --job-name=blab-tracer-dps
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=12
#SBATCH --mem=96G
#SBATCH --output=/shared/home/greg/breezelab-work/tracer-dp-scream/logs/slurm-%j.out
# TRACER–DP-SCREAM launcher following breezelab-work/COORDINATION.md (UUID-pinned GPU, shared lock dir,
# physical-idleness check, per-CPU compiled depot overlay).
# Usage: sbatch --partition=<p> --nodelist=<node> --time=<wall> execution/tracer_dp_scream.sbatch \
#            <GPU UUID> <experiment> <stop_time seconds> [checkpoint interval s = 21600] [pickup_from dir]
# Environment knobs are passed through to cases/tracer_dp_scream.jl (Slurm exports the submitting environment):
#   TRACER_DP_SCREAM_NX (256), TRACER_DP_SCREAM_DX (200 m), TRACER_DP_SCREAM_SLICE_MINUTES (5), TRACER_DP_SCREAM_SLICE_HEIGHT (3000 m).
set -euo pipefail
W=/shared/home/greg/breezelab-work/tracer-dp-scream
C=/shared/home/greg/breezelab-runs/20261006
D=/shared/home/greg/breezelab-runs/docs-20261007-8dc824d
GPU_UUID=${1:?GPU UUID}
EXPERIMENT=${2:?experiment name}
STOP_TIME=${3:?stop time in seconds}
CHECKPOINT=${4:-21600}
PICKUP_FROM=${5:-}
LOCKS=/shared/home/greg/breezelab-work/gpu-locks
mkdir -p "$LOCKS" "$W/logs" "$W/runs"
exec 9>"$LOCKS/$GPU_UUID.lock"
flock -n 9 || { echo "GPU $GPU_UUID is locked by another BreezeLab job"; exit 1; }
[ -z "$(nvidia-smi -i "$GPU_UUID" --query-compute-apps=pid --format=csv,noheader)" ] || { echo 'Selected GPU is busy'; exit 1; }
export CUDA_VISIBLE_DEVICES="$GPU_UUID"
GPU_NAME=$(nvidia-smi -i "$GPU_UUID" --query-gpu=name --format=csv,noheader)
case "$GPU_NAME" in *H100*|*A100*) ;; *) echo "Unexpected GPU: $GPU_NAME"; exit 1 ;; esac
RUN="$W/runs/${EXPERIMENT}_job${SLURM_JOB_ID}"
mkdir -p "$RUN"
trap 'code=$?; if (( code != 0 )); then printf "exit=%s job=%s\n" "$code" "$SLURM_JOB_ID" > "$RUN/FAILED"; fi' EXIT
CPU_TAG=$(awk -F ': ' '/model name/ {print $2; exit}' /proc/cpuinfo | tr -c 'A-Za-z0-9' '_')
mkdir -p "$W/depot/$CPU_TAG"
export JULIA_DEPOT_PATH="$W/depot/$CPU_TAG:$D/depot-compute/$CPU_TAG:$D/depot-login:$C/depot:$HOME/.julia:"
export JULIA_PKG_PRECOMPILE_AUTO=0
export JULIA_NUM_PRECOMPILE_TASKS=6
export JULIA_NUM_THREADS=6
cd "$W"
{
  printf 'job=%s experiment=%s node=%s gpu=%s (%s) partition=%s started=%s\n' "$SLURM_JOB_ID" "$EXPERIMENT" "$(hostname)" "$GPU_UUID" "$GPU_NAME" "$SLURM_JOB_PARTITION" "$(date -u +%FT%TZ)"
  printf 'source_commit=%s branch=%s dirty=%s\n' "$(git rev-parse HEAD)" "$(git rev-parse --abbrev-ref HEAD)" "$(git status --porcelain | wc -l)"
  printf 'command=sbatch execution/tracer_dp_scream.sbatch %s %s %s %s %s\n' "$GPU_UUID" "$EXPERIMENT" "$STOP_TIME" "$CHECKPOINT" "$PICKUP_FROM"
  printf 'grid: NX=%s DX=%s m\n' "${TRACER_DP_SCREAM_NX:-256}" "${TRACER_DP_SCREAM_DX:-200}"
  printf 'julia=%s\n' "$C/julia-1.12.7/bin/julia"
} > "$RUN/RUN_INFO.txt"
sha256sum "$W/data/tracer_dp_scream/"*.nc > "$RUN/inputs.sha256"
cp "$W/Manifest.toml" "$W/Project.toml" "$RUN/"
export TRACER_DP_SCREAM_OUTPUT="$RUN" TRACER_DP_SCREAM_STOP_TIME="$STOP_TIME" TRACER_DP_SCREAM_CHECKPOINT="$CHECKPOINT" TRACER_DP_SCREAM_PICKUP_FROM="$PICKUP_FROM"
"$C/julia-1.12.7/bin/julia" --startup-file=no --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()' 2>&1 | tail -5
"$C/julia-1.12.7/bin/julia" --startup-file=no --project=. cases/tracer_dp_scream.jl
printf 'completed=%s\n' "$(date -u +%FT%TZ)" >> "$RUN/RUN_INFO.txt"
touch "$RUN/COMPLETE"
