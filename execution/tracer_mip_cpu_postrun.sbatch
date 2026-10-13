#!/bin/bash
# Full post-run analysis of a TRACER-MIP outer control (CPU partition). Chain it after the GPU job:
#   sbatch --partition=cpu-c6a24xlarge --dependency=afterany:<gpu job> --export=ALL,RUN=<run dir>,REFERENCE=<label=run dir> \
#          execution/tracer_mip_cpu_postrun.sbatch
# Runs, from this checkout: the ERA5 surface-wind check, the AMF1 site comparison, the extremes (stability) timeline,
# the multi-run comparison against ERA5 and AMF1 (RUN labelled by LABEL, plus REFERENCE), and the animations at the
# run's own cadence; copies figures/animations to the campaign-figures directory.
#SBATCH --job-name=blab-tmip-postrun
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=150G
#SBATCH --time=12:00:00
#SBATCH --output=/shared/home/greg/breezelab-work/tracer-mip/logs/slurm-%j.out
set -uo pipefail
CHECKOUT=${CHECKOUT:-$(cd "$(dirname "$0")/.." && pwd)}
W=${DEPOT_ROOT:-/shared/home/greg/breezelab-work/tracer-mip}
C=/shared/home/greg/breezelab-runs/20261006; D=/shared/home/greg/breezelab-runs/docs-20261007-8dc824d
CPU_TAG=$(awk -F ': ' '/model name/ {print $2; exit}' /proc/cpuinfo | tr -c 'A-Za-z0-9' '_')
source "$W/execution/depot_layout.sh"
export JULIA_DEPOT_PATH="$W/depot/$CPU_TAG:$D/depot-compute/$CPU_TAG:$D/depot-login:$C/depot:$HOME/.julia:"
export JULIA_PKG_PRECOMPILE_AUTO=0 JULIA_NUM_THREADS=8
J="$C/julia-1.12.7/bin/julia --startup-file=no --project=cases/tracer_mip"
RUN=${RUN:?run directory}; LABEL=${LABEL:-$(basename "$RUN")}
ERA5=${ERA5:-$W/data/era5}; ARM=${ARM:-$W/data/arm_tracer_20220807}
cd "$CHECKOUT"
echo "== $(date -u) node=$(hostname) commit=$(git rev-parse HEAD) run=$RUN status=$(ls "$RUN" | grep -c COMPLETE)"
step() { echo "== $(date -u) $1"; shift; "$@" 2>&1 | grep --line-buffered -v "^\s*$" | tail -40; }
step "surface wind vs ERA5"        $J analysis/tracer_mip_surface_wind_check.jl "$RUN" "$ERA5" "$RUN/diagnostics"
step "site comparison vs AMF1"     $J analysis/tracer_mip_site_comparison.jl "$RUN" "$ARM" "$RUN/analysis"
step "extremes timeline"           $J analysis/tracer_mip_blowup_diagnosis.jl "$RUN" "$RUN/analysis"
step "multi-run comparison"        $J analysis/tracer_mip_control_comparison.jl "$RUN/analysis/comparison" "$ERA5" "$ARM" "$LABEL=$RUN" ${REFERENCE:-}
# Playback speed in sim-min per video second (fps follows from the frame spacing); see analysis/animate_tracer_mip.jl
step "animations"                  env ANIMATION_MIN_PER_S=${ANIMATION_MIN_PER_S:-25} ANIMATION_SECTION_MIN_PER_S=${ANIMATION_SECTION_MIN_PER_S:-12.5} \
                                       ANIMATION_WINDOW_MIN_PER_S=${ANIMATION_WINDOW_MIN_PER_S:-12.5} ANIMATION_WINDOW_UTC=${ANIMATION_WINDOW_UTC:-17-01} \
                                       $J analysis/animate_tracer_mip.jl "$RUN" "$RUN/animations"
FIG=/shared/home/greg/breezelab-work/campaign-figures/tracer_mip_$(basename "$RUN"); mkdir -p "$FIG"
cp "$RUN"/animations/*.mp4 "$RUN"/animations/*.png "$RUN"/animations/animations.toml "$FIG"/ 2>/dev/null
cp "$RUN"/analysis/*.png "$RUN"/analysis/*.toml "$RUN"/analysis/comparison/* "$RUN"/diagnostics/* "$FIG"/ 2>/dev/null
cp analysis/README_animations.md "$FIG/README_animations.md" 2>/dev/null
printf 'Post-run analysis of %s (%s), commit %s, %s\n' "$RUN" "$(ls "$RUN" | grep -q COMPLETE && echo complete || echo INCOMPLETE)" "$(git rev-parse --short HEAD)" "$(date -u)" > "$FIG/PROVENANCE.txt"
echo "== $(date -u) copied to $FIG"
