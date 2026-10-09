#!/bin/bash
#SBATCH --job-name=blab-ena-lasso
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=12
#SBATCH --mem=96G
#SBATCH --time=06:00:00
#SBATCH --output=/shared/home/greg/breezelab-work/ena-lasso/logs/slurm-%j.out
# Launch one official LASSO-ENA member on a physically idle H100/A100 with the shared
# UUID-pinning procedure (breezelab-work/COORDINATION.md). Usage from the checkout root:
#
#   sbatch --partition=<gpu-prod|gpu-p4de-2c> --nodelist=<node> --time=<wall> \
#       execution/submit_ena_lasso.sbatch <GPU UUID> <experiment> [julia args for cases/cli/run_case.jl ...]
#
# With no extra arguments the readable case script cases/ena_lasso.jl runs (member from
# ENA_LASSO_MEMBER, bundle from ENA_LASSO_BUNDLE); outputs land in runs/<experiment>_job<id>/.
set -euo pipefail
# BREEZELAB_CHECKOUT selects the source tree (e.g. a worktree at a merged main); runs/ and the
# launch record follow it. The per-CPU depot stays in the ena-lasso checkout so compile caches are shared.
W=${BREEZELAB_CHECKOUT:-/shared/home/greg/breezelab-work/ena-lasso}
D=/shared/home/greg/breezelab-work/ena-lasso
C=/shared/home/greg/breezelab-runs/20261006
GPU_UUID=${1:?GPU UUID or 'auto'}
EXPERIMENT=${2:?experiment name}
shift 2
mkdir -p /shared/home/greg/breezelab-work/gpu-locks "$W/logs" "$W/runs"
LOCKS=/shared/home/greg/breezelab-work/gpu-locks
# `auto`: take the first GPU of this node that has no compute process and no held lock (the
# campaigns share $LOCKS); the chosen UUID is pinned below exactly like an explicit one.
if [ "$GPU_UUID" = auto ]; then
    CHOSEN=""
    for uuid in $(nvidia-smi --query-gpu=uuid --format=csv,noheader); do
        [ -z "$(nvidia-smi -i "$uuid" --query-compute-apps=pid --format=csv,noheader)" ] || continue
        exec 9>"$LOCKS/$uuid.lock"
        if flock -n 9; then CHOSEN=$uuid; break; fi
        exec 9>&-
    done
    [ -n "$CHOSEN" ] || { echo "No physically idle GPU on $(hostname) at $(date -u +%FT%TZ); exiting without running"; exit 1; }
    GPU_UUID=$CHOSEN
else
    exec 9>"$LOCKS/$GPU_UUID.lock"
    flock -n 9 || { echo "GPU $GPU_UUID is locked by another campaign"; exit 1; }
fi
[ -z "$(nvidia-smi -i "$GPU_UUID" --query-compute-apps=pid --format=csv,noheader)" ] || { echo "Selected GPU is busy"; exit 1; }
export CUDA_VISIBLE_DEVICES="$GPU_UUID"
GPU_NAME=$(nvidia-smi -i "$CUDA_VISIBLE_DEVICES" --query-gpu=name --format=csv,noheader)
case "$GPU_NAME" in *H100*|*A100*) ;; *) echo "Unexpected GPU: $GPU_NAME"; exit 1 ;; esac
CPU_TAG=$(awk -F ': ' '/model name/ {print $2; exit}' /proc/cpuinfo | tr -c 'A-Za-z0-9' '_')
mkdir -p "$D/depot/$CPU_TAG"
export JULIA_DEPOT_PATH="$D/depot/$CPU_TAG:/shared/home/greg/breezelab-runs/docs-20261007-8dc824d/depot-compute/$CPU_TAG:/shared/home/greg/breezelab-runs/docs-20261007-8dc824d/depot-login:$C/depot:$HOME/.julia:"
export JULIA_PKG_PRECOMPILE_AUTO=0
export JULIA_NUM_PRECOMPILE_TASKS=6
export JULIA_NUM_THREADS=6
export JULIA_CUDA_USE_COMPAT=false
OUT="$W/runs/${EXPERIMENT}_job${SLURM_JOB_ID}"
mkdir -p "$OUT"
cd "$W"
printf 'job=%s gpu=%s (%s) node=%s source=%s output=%s\n' "$SLURM_JOB_ID" "$GPU_UUID" "$GPU_NAME" "$(hostname)" "$(git rev-parse HEAD)" "$OUT" | tee "$OUT/launch.txt"
printf 'command: %s\n' "$*" >> "$OUT/launch.txt"
if [ "$#" -eq 0 ]; then
    # Readable case script: redirect its output directory through the environment-free
    # convention by running it from a wrapper that sets pkgdir(...)/output → $OUT.
    ENA_LASSO_OUTPUT="$OUT" "$C/julia-1.12.7/bin/julia" --startup-file=no --project=. -e '
        using BreezeLab
        script = read(joinpath(pkgdir(BreezeLab), "cases", "ena_lasso.jl"), String)
        script = replace(script, "output_dir = joinpath(pkgdir(BreezeLab), \"output\", \"ena_lasso\", member)" => "output_dir = ENV[\"ENA_LASSO_OUTPUT\"]")
        include_string(Main, script, "ena_lasso.jl")'
else
    "$C/julia-1.12.7/bin/julia" --startup-file=no --project=. cases/cli/run_case.jl "$@" --output "$OUT"
fi
touch "$OUT/COMPLETE"
