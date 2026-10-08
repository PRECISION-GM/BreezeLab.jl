#!/bin/bash
# Re-render the TRACER-MIP animations once the control run has finished: snapshot its 10-min files,
# render, and copy the products to the campaign figures directory.
#   sbatch --partition=cpu-c6a24xlarge --dependency=afterany:<gpu job> --export=ALL,SRC=<run dir> execution/tracer_mip_cpu_animate_final.sbatch
#SBATCH --job-name=blab-tmip-anim-final
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=06:00:00
#SBATCH --output=/shared/home/greg/breezelab-work/tracer-mip/logs/slurm-%j.out
set -euo pipefail
SRC=${SRC:?run directory}
SNAP="$SRC/../anim_snapshot_final_$(basename "$SRC")"
mkdir -p "$SNAP/inner_region"
cp "$SRC"/outer_slices.jld2 "$SRC"/outer_surface.jld2 "$SRC"/provenance.toml "$SRC"/COMMIT "$SNAP"/
cp "$SRC"/inner_region/inner_region_state.jld2 "$SRC"/inner_region/README.txt "$SNAP"/inner_region/
date -u > "$SNAP/SNAPSHOT_UTC"; ls "$SRC" | grep -q COMPLETE && echo complete > "$SNAP/RUN_STATUS" || echo incomplete > "$SNAP/RUN_STATUS"
export RUN="$SNAP" OUT="$SRC/animations"
bash "$(dirname "$0")/tracer_mip_cpu_animate.sbatch"
FIG=/shared/home/greg/breezelab-work/campaign-figures/animations
mkdir -p "$FIG"
cp "$OUT"/*.mp4 "$OUT"/*.png "$OUT"/animations.toml "$FIG"/
cp "$(dirname "$0")/../analysis/README_animations.md" "$FIG/README.md"
printf '\nRendered %s from %s (run status: %s), commit %s\n' "$(date -u)" "$SNAP" "$(cat "$SNAP/RUN_STATUS")" "$(git -C "$(dirname "$0")/.." rev-parse HEAD)" >> "$FIG/README.md"
echo "== $(date -u) copied to $FIG"
