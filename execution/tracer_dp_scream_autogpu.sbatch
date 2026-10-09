#!/bin/bash
#SBATCH --job-name=blab-tracer-dps
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=12
#SBATCH --mem=96G
#SBATCH --output=/shared/home/greg/breezelab-work/tracer-dp-scream/logs/slurm-%j.out
# TRACER–DP-SCREAM launcher that selects a physically idle GPU at start time (the COORDINATION.md procedure
# without a pre-chosen UUID, as execution/sea_starr_autogpu.sbatch): Slurm places the job; inside the
# allocation the first A100/H100 on the node with no compute process and a free shared lock is pinned.
# Submit from the checkout to run (its HEAD is recorded); outputs go to the case run directory.
#   cd <checkout at the commit to run>
#   TRACER_DP_SCREAM_NX=... TRACER_DP_SCREAM_DX=... TRACER_DP_SCREAM_SLICE_MINUTES=5 TRACER_DP_SCREAM_SLICE_HEIGHT=3000 \
#     sbatch --partition=gpu-p4de-2c --time=<wall> execution/tracer_dp_scream_autogpu.sbatch <experiment> <stop_time s> [checkpoint s = 21600] [pickup_from dir]
set -euo pipefail
W=${SLURM_SUBMIT_DIR:?submit from the checkout to run}
RUNS=/shared/home/greg/breezelab-work/tracer-dp-scream/runs
C=/shared/home/greg/breezelab-runs/20261006
D=/shared/home/greg/breezelab-runs/docs-20261007-8dc824d
EXPERIMENT=${1:?experiment name}
STOP_TIME=${2:?stop time in seconds}
CHECKPOINT=${3:-21600}
PICKUP_FROM=${4:-}
LOCKDIR=/shared/home/greg/breezelab-work/gpu-locks
mkdir -p "$LOCKDIR" "$RUNS"

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

RUN="$RUNS/${EXPERIMENT}_job${SLURM_JOB_ID}"
mkdir -p "$RUN"
trap 'code=$?; if (( code != 0 )); then printf "exit=%s job=%s\n" "$code" "$SLURM_JOB_ID" > "$RUN/FAILED"; fi' EXIT
CPU_TAG=$(awk -F ': ' '/model name/ {print $2; exit}' /proc/cpuinfo | tr -c 'A-Za-z0-9' '_')
DEPOT=/shared/home/greg/breezelab-work/tracer-dp-scream/depot
mkdir -p "$DEPOT/$CPU_TAG"
export JULIA_DEPOT_PATH="$DEPOT/$CPU_TAG:$D/depot-compute/$CPU_TAG:$D/depot-login:$C/depot:$HOME/.julia:"
export JULIA_PKG_PRECOMPILE_AUTO=0 JULIA_NUM_PRECOMPILE_TASKS=6 JULIA_NUM_THREADS=6 JULIA_CUDA_USE_COMPAT=false
cd "$W"
{
  printf 'job=%s experiment=%s node=%s gpu=%s (%s; auto-selected, slurm allocation: %s) partition=%s started=%s\n' \
         "$SLURM_JOB_ID" "$EXPERIMENT" "$(hostname)" "$GPU_UUID" "$GPU_NAME" "${allocated:-none}" "$SLURM_JOB_PARTITION" "$(date -u +%FT%TZ)"
  printf 'source=%s commit=%s branch=%s dirty=%s\n' "$W" "$(git rev-parse HEAD)" "$(git rev-parse --abbrev-ref HEAD)" "$(git status --porcelain --untracked-files=no | wc -l)"
  printf 'command=sbatch execution/tracer_dp_scream_autogpu.sbatch %s %s %s %s\n' "$EXPERIMENT" "$STOP_TIME" "$CHECKPOINT" "$PICKUP_FROM"
  printf 'grid: NX=%s DX=%s m; slices every %s min at %s m\n' "${TRACER_DP_SCREAM_NX:-256}" "${TRACER_DP_SCREAM_DX:-200}" \
         "${TRACER_DP_SCREAM_SLICE_MINUTES:-5}" "${TRACER_DP_SCREAM_SLICE_HEIGHT:-3000}"
  printf 'julia=%s\n' "$C/julia-1.12.7/bin/julia"
} > "$RUN/RUN_INFO.txt"
sha256sum "$W/data/tracer_dp_scream/"*.nc > "$RUN/inputs.sha256"
cp "$W/Manifest.toml" "$W/Project.toml" "$RUN/"
export TRACER_DP_SCREAM_OUTPUT="$RUN" TRACER_DP_SCREAM_STOP_TIME="$STOP_TIME" TRACER_DP_SCREAM_CHECKPOINT="$CHECKPOINT" TRACER_DP_SCREAM_PICKUP_FROM="$PICKUP_FROM"
"$C/julia-1.12.7/bin/julia" --startup-file=no --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()' 2>&1 | tail -5
"$C/julia-1.12.7/bin/julia" --startup-file=no --project=. cases/tracer_dp_scream.jl
printf 'completed=%s\n' "$(date -u +%FT%TZ)" >> "$RUN/RUN_INFO.txt"
touch "$RUN/COMPLETE"
