#!/bin/bash
#SBATCH --job-name=blab-seastarr
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=12
#SBATCH --mem=96G
#SBATCH --requeue
#SBATCH --output=/shared/home/greg/breezelab-work/sea-starr/logs/slurm-%j.out
# Requeue-safe SEA STARR launcher (one output directory per segment; never reopens an earlier directory):
#   sbatch --partition=gpu-p4de-2c,gpu-prod --time=<wall> execution/sea_starr_segmented.sbatch <code worktree> <experiment> [hours=66] [Nx=192] [member=CTRL] [twod_minutes=3]
# First start: runs <code>/cases/sea_starr.jl into runs/<experiment>_job<id>/. After a requeue (SLURM_RESTART_COUNT > 0, or the first
# directory already holds output): writes runs/<experiment>_job<id>_r<k>/ and resumes from the highest-iteration checkpoint found in the
# earlier segments with <code>/execution/sea_starr_resume.jl (or starts fresh there if no checkpoint exists yet).
# GPU: physically idle, unlocked A100/H100 chosen at start (shared lock dir), as sea_starr_autogpu.sbatch. DRY_RUN=1 prints the plan only.
set -euo pipefail
W=/shared/home/greg/breezelab-work/sea-starr
C=/shared/home/greg/breezelab-runs/20261006
CODE=${1:?code worktree}; EXPERIMENT=${2:?experiment}; HOURS=${3:-66}; NX=${4:-192}; MEMBER=${5:-CTRL}; TWOD=${6:-3}
PREFIX="sea_starr_$(echo "$MEMBER" | tr 'A-Z' 'a-z')"
JOBID=${SLURM_JOB_ID:?}; RESTART=${SLURM_RESTART_COUNT:-0}
BASE="$W/runs/${EXPERIMENT}_job${JOBID}"
has_output() { compgen -G "$1/${PREFIX}_*.jld2" > /dev/null 2>&1; }
if [ "$RESTART" -eq 0 ] && ! has_output "$BASE"; then
    RUN="$BASE"; MODE=fresh; FROM=""; CKPT=""
else
    k=$RESTART; [ "$k" -lt 1 ] && k=1
    RUN="${BASE}_r${k}"; while [ -e "$RUN" ] && has_output "$RUN"; do k=$((k + 1)); RUN="${BASE}_r${k}"; done
    CKPT=""; best=-1
    for d in "$BASE" "${BASE}"_r*; do
        [ -d "$d" ] && [ "$d" != "$RUN" ] || continue
        for f in "$d"/${PREFIX}_checkpoint_iteration*.jld2; do
            [ -e "$f" ] || continue
            it=$(basename "$f" .jld2); it=${it##*_iteration}
            if [ "$it" -gt "$best" ]; then best=$it; CKPT=$f; fi
        done
    done
    if [ -n "$CKPT" ]; then MODE=resume; FROM=$(dirname "$CKPT"); else MODE=fresh; FROM=""; fi
fi
echo "plan: job=$JOBID restart=$RESTART mode=$MODE run=$RUN from=${FROM:-none} checkpoint=${CKPT:-none} code=$CODE"
[ "${DRY_RUN:-0}" = "1" ] && exit 0

LOCKDIR=/shared/home/greg/breezelab-work/gpu-locks
mkdir -p "$LOCKDIR" "$W/logs" "$W/runs"
allocated=${CUDA_VISIBLE_DEVICES:-}
mapfile -t ALL < <(nvidia-smi --query-gpu=uuid,name --format=csv,noheader | sed 's/, /,/')
candidates=()
for entry in "${ALL[@]}"; do
    uuid=${entry%%,*}; name=${entry#*,}
    case "$name" in *H100*|*A100*) ;; *) continue ;; esac
    idx=$(nvidia-smi --query-gpu=index --format=csv,noheader -i "$uuid")
    if [[ -n "$allocated" && ( "$allocated" == "$uuid" || "$allocated" == "$idx" ) ]]; then candidates=("$uuid" "${candidates[@]}"); else candidates+=("$uuid"); fi
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

mkdir -p "$RUN"
trap 'code=$?; if (( code != 0 )); then printf "exit=%s job=%s restart=%s\n" "$code" "$JOBID" "$RESTART" > "$RUN/FAILED"; fi' EXIT
CPU_TAG=$(awk -F ': ' '/model name/ {print $2; exit}' /proc/cpuinfo | tr -c 'A-Za-z0-9' '_')
mkdir -p "$W/depot/$CPU_TAG"
export JULIA_DEPOT_PATH="$W/depot/$CPU_TAG:/shared/home/greg/breezelab-runs/docs-20261007-8dc824d/depot-compute/$CPU_TAG:/shared/home/greg/breezelab-runs/docs-20261007-8dc824d/depot-login:$C/depot:$HOME/.julia:"
export JULIA_PKG_PRECOMPILE_AUTO=0 JULIA_NUM_PRECOMPILE_TASKS=6 JULIA_NUM_THREADS=6 JULIA_CUDA_USE_COMPAT=false
export SEA_STARR_ARCH=gpu SEA_STARR_HOURS="$HOURS" SEA_STARR_NX="$NX" SEA_STARR_MEMBER="$MEMBER" SEA_STARR_OUTPUT="$RUN"
export SEA_STARR_TWOD_MINUTES="$TWOD" SEA_STARR_CHECKPOINT_HOURS="${SEA_STARR_CHECKPOINT_HOURS:-3}"
export SEA_STARR_DATA="$C/inputs/mip_sources/seastarr_22241697"
{
  echo "job=$JOBID restart=$RESTART mode=$MODE node=$(hostname) gpu=$GPU_NAME uuid=$GPU_UUID (auto-selected; slurm allocation: ${allocated:-none})"
  echo "code=$CODE commit=$(git -C "$CODE" rev-parse HEAD) dirty=$(git -C "$CODE" status --porcelain | wc -l)"
  echo "launcher=$0 sha256=$(sha256sum "$0" | cut -c1-64)"
  echo "command: sbatch ... execution/sea_starr_segmented.sbatch $*"
  [ "$MODE" = resume ] && echo "segment of $FROM; resumed from $CKPT (copied in)"
  echo "started=$(date -u +%FT%TZ)"
  sha256sum "$SEA_STARR_DATA"/SEA_STARR_${MEMBER}_SCM_driver.nc
  cp "$CODE/Manifest.toml" "$RUN/Manifest.toml"
} | tee "$RUN/job.txt"
cd "$CODE"
"$C/julia-1.12.7/bin/julia" --startup-file=no --project="$CODE" -e 'using Pkg; Pkg.precompile()' 2>&1 | tail -3
if [ "$MODE" = resume ]; then
    cp "$CKPT" "$RUN/"
    export SEA_STARR_RESUME_FROM="$FROM" SEA_STARR_RESUME_CHECKPOINT="$RUN/$(basename "$CKPT")"
    "$C/julia-1.12.7/bin/julia" --startup-file=no --project="$CODE" "$CODE/execution/sea_starr_resume.jl" 2>&1 | tee "$RUN/run.log"
else
    "$C/julia-1.12.7/bin/julia" --startup-file=no --project="$CODE" "$CODE/cases/sea_starr.jl" 2>&1 | tee "$RUN/run.log"
fi
echo "finished=$(date -u +%FT%TZ)" >> "$RUN/job.txt"
touch "$RUN/COMPLETE"
