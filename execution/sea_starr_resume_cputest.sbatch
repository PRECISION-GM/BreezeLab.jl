#!/bin/bash
#SBATCH --job-name=blab-seastarr-resume-cpu
#SBATCH --partition=cpu-c6a24xlarge
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=03:00:00
#SBATCH --output=/shared/home/greg/breezelab-work/sea-starr/logs/slurm-%j.out
# CPU pre-check of a resume (no GPU): restore the checkpoint read-only from the original run directory into a
# scratch segment directory, report continuity, and take a few steps. Same arguments as sea_starr_resume.sbatch.
set -euo pipefail
W=/shared/home/greg/breezelab-work/sea-starr
C=/shared/home/greg/breezelab-runs/20261006
HERE=/shared/home/greg/breezelab-work/sea-starr-wt/resume/execution
CODE=${1:?code worktree}; EXPERIMENT=${2:?experiment}; FROM=${3:?original run dir}; CKPT=${4:?checkpoint file}
HOURS=${5:-66}; NX=${6:-192}; MEMBER=${7:-CTRL}; TWOD=${8:-}
RUN="$W/runs/resume_cputest_${EXPERIMENT}_job${SLURM_JOB_ID}"
mkdir -p "$RUN"
CPU_TAG=$(awk -F ': ' '/model name/ {print $2; exit}' /proc/cpuinfo | tr -c 'A-Za-z0-9' '_')
mkdir -p "$W/depot/$CPU_TAG"
export JULIA_DEPOT_PATH="$W/depot/$CPU_TAG:/shared/home/greg/breezelab-runs/docs-20261007-8dc824d/depot-compute/$CPU_TAG:/shared/home/greg/breezelab-runs/docs-20261007-8dc824d/depot-login:$C/depot:$HOME/.julia:"
export JULIA_PKG_PRECOMPILE_AUTO=0 JULIA_NUM_PRECOMPILE_TASKS=16 JULIA_NUM_THREADS=16
export SEA_STARR_ARCH=cpu SEA_STARR_RESUME_STEPS="${SEA_STARR_RESUME_STEPS:-2}" SEA_STARR_HOURS="$HOURS" SEA_STARR_NX="$NX"
export SEA_STARR_MEMBER="$MEMBER" SEA_STARR_OUTPUT="$RUN" SEA_STARR_DATA="$C/inputs/mip_sources/seastarr_22241697"
export SEA_STARR_RESUME_FROM="$FROM" SEA_STARR_RESUME_CHECKPOINT="$CKPT"
[ -n "$TWOD" ] && export SEA_STARR_TWOD_MINUTES="$TWOD"
echo "job=$SLURM_JOB_ID node=$(hostname) cpu=$CPU_TAG code=$CODE commit=$(git -C "$CODE" rev-parse HEAD) from=$FROM ckpt=$CKPT started=$(date -u +%FT%TZ)"
cd "$CODE"
"$C/julia-1.12.7/bin/julia" --startup-file=no --project="$CODE" -e 'using Pkg; Pkg.precompile()' 2>&1 | tail -3
"$C/julia-1.12.7/bin/julia" --startup-file=no --project="$CODE" "$HERE/sea_starr_resume.jl" 2>&1
echo "finished=$(date -u +%FT%TZ)"
