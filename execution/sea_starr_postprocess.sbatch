#!/bin/bash
#SBATCH --job-name=blab-seastarr-post
#SBATCH --partition=cpu-c6a24xlarge
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/shared/home/greg/breezelab-work/sea-starr/logs/slurm-%j.out
# Post-process the SEA STARR members once their restart segments complete (submit with --dependency=afterok:<jobs>):
#   sbatch --dependency=afterok:<n100 seg2>:<n030 seg2> execution/sea_starr_postprocess.sbatch <ctrl dirs> <n100 dirs> <n030 dirs> <out dir>
# Each "dirs" argument is a comma-separated segment list. Runs the CTRL/N100/N030 analysis and re-renders the
# animations over the stitched segments, then copies figures/animations to campaign-figures.
set -euo pipefail
W=/shared/home/greg/breezelab-work/sea-starr
C=/shared/home/greg/breezelab-runs/20261006
WT=/shared/home/greg/breezelab-work/sea-starr-wt/resume
FIG=/shared/home/greg/breezelab-work/campaign-figures
CTRL=${1:?ctrl dirs}; N100=${2:?n100 dirs}; N030=${3:?n030 dirs}; OUT=${4:?output dir}
CPU_TAG=$(awk -F ': ' '/model name/ {print $2; exit}' /proc/cpuinfo | tr -c 'A-Za-z0-9' '_')
export JULIA_DEPOT_PATH="$W/depot/$CPU_TAG:/shared/home/greg/breezelab-runs/docs-20261007-8dc824d/depot-compute/$CPU_TAG:/shared/home/greg/breezelab-runs/docs-20261007-8dc824d/depot-login:$C/depot:$HOME/.julia:"
export JULIA_PKG_PRECOMPILE_AUTO=0 JULIA_NUM_PRECOMPILE_TASKS=8 JULIA_NUM_THREADS=8
cd "$WT"
echo "job=$SLURM_JOB_ID node=$(hostname) commit=$(git rev-parse HEAD) started=$(date -u +%FT%TZ)"
J="$C/julia-1.12.7/bin/julia --startup-file=no --project=$WT"
$J analysis/sea_starr_members.jl --ctrl "$CTRL" --n100 "$N100" --n030 "$N030" --out "$OUT"
$J analysis/animate_sea_starr.jl --ctrl "$CTRL" --n100 "$N100" --n030 "$N030" --copy-to "$FIG/animations"
mkdir -p "$FIG"
for f in members_timeseries members_state_space; do cp "$OUT/$f.png" "$FIG/seastarr_$f.png"; done
cp "$OUT/members_summary.md" "$FIG/seastarr_members_summary.md"
echo "finished=$(date -u +%FT%TZ)"
