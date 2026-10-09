#!/bin/bash
#SBATCH --job-name=blab-seastarr-resume
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=12
#SBATCH --mem=96G
#SBATCH --output=/shared/home/greg/breezelab-work/sea-starr/logs/slurm-%j.out
# Resume a SEA STARR run from a checkpoint into a NEW segment directory, with the code the run started from:
#   sbatch --partition=gpu-p4de-2c --time=<wall> execution/sea_starr_resume.sbatch \
#          <code worktree> <experiment> <original run dir> <checkpoint file> [hours=66] [Nx=192] [member=CTRL] [twod_minutes]
# The GPU is selected at job start (physically idle, shared lock free, A100/H100 only), as in sea_starr_autogpu.sbatch;
# the checkpoint is copied into runs/<experiment>_job<id>/ and restored from there by execution/sea_starr_resume.jl.
set -euo pipefail
W=/shared/home/greg/breezelab-work/sea-starr
C=/shared/home/greg/breezelab-runs/20261006
HERE=/shared/home/greg/breezelab-work/sea-starr-wt/resume/execution
CODE=${1:?code worktree}; EXPERIMENT=${2:?experiment}; FROM=${3:?original run dir}; CKPT=${4:?checkpoint file}
HOURS=${5:-66}; NX=${6:-192}; MEMBER=${7:-CTRL}; TWOD=${8:-}
LOCKDIR=/shared/home/greg/breezelab-work/gpu-locks
mkdir -p "$LOCKDIR" "$W/logs" "$W/runs"

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
[ "$(realpath -m "$RUN")" != "$(realpath "$FROM")" ] || { echo "segment directory equals the original run directory"; exit 1; }
mkdir -p "$RUN"
trap 'code=$?; if (( code != 0 )); then printf "exit=%s job=%s\n" "$code" "$SLURM_JOB_ID" > "$RUN/FAILED"; fi' EXIT

CPU_TAG=$(awk -F ': ' '/model name/ {print $2; exit}' /proc/cpuinfo | tr -c 'A-Za-z0-9' '_')
mkdir -p "$W/depot/$CPU_TAG"
export JULIA_DEPOT_PATH="$W/depot/$CPU_TAG:/shared/home/greg/breezelab-runs/docs-20261007-8dc824d/depot-compute/$CPU_TAG:/shared/home/greg/breezelab-runs/docs-20261007-8dc824d/depot-login:$C/depot:$HOME/.julia:"
export JULIA_PKG_PRECOMPILE_AUTO=0 JULIA_NUM_PRECOMPILE_TASKS=6 JULIA_NUM_THREADS=6 JULIA_CUDA_USE_COMPAT=false
export SEA_STARR_ARCH=gpu SEA_STARR_HOURS="$HOURS" SEA_STARR_NX="$NX" SEA_STARR_MEMBER="$MEMBER" SEA_STARR_OUTPUT="$RUN"
export SEA_STARR_DATA="$C/inputs/mip_sources/seastarr_22241697" SEA_STARR_CHECKPOINT_HOURS="${SEA_STARR_CHECKPOINT_HOURS:-3}"
export SEA_STARR_RESUME_FROM="$FROM"
[ -n "$TWOD" ] && export SEA_STARR_TWOD_MINUTES="$TWOD"

cp "$HERE/sea_starr_resume.jl" "$RUN/"
cp "$CKPT" "$RUN/"
export SEA_STARR_RESUME_CHECKPOINT="$RUN/$(basename "$CKPT")"
{
  echo "job=$SLURM_JOB_ID node=$(hostname) gpu=$GPU_NAME uuid=$GPU_UUID (auto-selected; slurm allocation: ${allocated:-none})"
  echo "segment of $FROM; resumed from $CKPT (copied in)"
  echo "code=$CODE commit=$(git -C "$CODE" rev-parse HEAD) dirty=$(git -C "$CODE" status --porcelain | wc -l)"
  echo "driver=$HERE/sea_starr_resume.jl sha256=$(sha256sum "$HERE/sea_starr_resume.jl" | cut -c1-64)"
  echo "command: sbatch ... execution/sea_starr_resume.sbatch $*"
  echo "started=$(date -u +%FT%TZ)"
  cp "$CODE/Manifest.toml" "$RUN/Manifest.toml"
} | tee "$RUN/job.txt"

cd "$CODE"
"$C/julia-1.12.7/bin/julia" --startup-file=no --project="$CODE" -e 'using Pkg; Pkg.precompile()' 2>&1 | tail -3
"$C/julia-1.12.7/bin/julia" --startup-file=no --project="$CODE" "$RUN/sea_starr_resume.jl" 2>&1 | tee "$RUN/run.log"
echo "finished=$(date -u +%FT%TZ)" >> "$RUN/job.txt"
touch "$RUN/COMPLETE"
