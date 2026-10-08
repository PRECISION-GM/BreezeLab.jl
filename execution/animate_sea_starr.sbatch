#!/bin/bash
#SBATCH --job-name=blab-seastarr-anim
#SBATCH --partition=cpu-c6a24xlarge
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/shared/home/greg/breezelab-work/sea-starr/logs/slurm-%j.out
# Render the SEA STARR animations on the CPU partition from whatever output the 66 h runs
# have written so far (read-only; the running GPU jobs are not touched):
#   sbatch execution/animate_sea_starr.sbatch [extra arguments for analysis/animate_sea_starr.jl]
set -euo pipefail
W=/shared/home/greg/breezelab-work/sea-starr
C=/shared/home/greg/breezelab-runs/20261006
CPU_TAG=$(awk -F ': ' '/model name/ {print $2; exit}' /proc/cpuinfo | tr -c 'A-Za-z0-9' '_')
mkdir -p "$W/depot/$CPU_TAG" "$W/logs"
export JULIA_DEPOT_PATH="$W/depot/$CPU_TAG:/shared/home/greg/breezelab-runs/docs-20261007-8dc824d/depot-compute/$CPU_TAG:/shared/home/greg/breezelab-runs/docs-20261007-8dc824d/depot-login:$C/depot:$HOME/.julia:"
export JULIA_PKG_PRECOMPILE_AUTO=0
export JULIA_NUM_PRECOMPILE_TASKS=8
export JULIA_NUM_THREADS=8
cd "$W"
echo "job=$SLURM_JOB_ID node=$(hostname) cpu=$CPU_TAG commit=$(git rev-parse HEAD) started=$(date -u +%FT%TZ)"
"$C/julia-1.12.7/bin/julia" --startup-file=no --project=. -e 'using Pkg; Pkg.precompile()' 2>&1 | tail -3
"$C/julia-1.12.7/bin/julia" --startup-file=no --project=. analysis/animate_sea_starr.jl "$@" 2>&1
echo "finished=$(date -u +%FT%TZ)"
