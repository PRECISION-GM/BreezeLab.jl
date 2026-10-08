#!/bin/bash
#SBATCH --job-name=blab-tracer-anim
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=03:00:00
#SBATCH --output=logs/slurm-%j.out
# CPU-partition render of analysis/animate_tracer_dp_scream.jl for one TRACER–DP-SCREAM run.
# Usage (from the repository root, `mkdir -p logs` first):
#   sbatch --partition=cpu-c6a24xlarge execution/animate_tracer_dp_scream.sbatch <run dir> <output dir> [<snapshot dir>]
# With a snapshot dir the run's *_slices/_timeseries/_profiles.jld2 and provenance.toml are copied there first so a
# still-running GPU job's files are never read while they are being appended; the snapshot is then rendered.
# Julia/depot layout follows breezelab-work/COORDINATION.md (campaign binary, per-CPU compiled depot under ./depot).
set -euo pipefail
RUN=${1:?run dir}; OUT=${2:?output dir}; SNAP=${3:-}
P=${SLURM_SUBMIT_DIR:-$(pwd)}
C=/shared/home/greg/breezelab-runs/20261006
D=/shared/home/greg/breezelab-runs/docs-20261007-8dc824d
CPU_TAG=$(awk -F ': ' '/model name/ {print $2; exit}' /proc/cpuinfo | tr -c 'A-Za-z0-9' '_')
mkdir -p "$P/depot/$CPU_TAG"
export JULIA_DEPOT_PATH="$P/depot/$CPU_TAG:$D/depot-compute/$CPU_TAG:$D/depot-login:$C/depot:$HOME/.julia:"
export JULIA_PKG_PRECOMPILE_AUTO=0 JULIA_NUM_THREADS=4 JULIA_NUM_PRECOMPILE_TASKS=8
cd "$P"; echo "commit=$(git rev-parse HEAD) node=$(hostname) started=$(date -u +%FT%TZ)"
J="$C/julia-1.12.7/bin/julia"
"$J" --startup-file=no --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()' 2>&1 | tail -3
SRC=$RUN
if [ -n "$SNAP" ]; then mkdir -p "$SNAP"; cp "$RUN"/*_{slices,timeseries,profiles}.jld2 "$RUN"/provenance.toml "$SNAP"/; SRC=$SNAP; fi
export DP_SCREAM_DATA=${DP_SCREAM_DATA:-$P/data/tracer_dp_scream}
time "$J" --startup-file=no --project=. analysis/animate_tracer_dp_scream.jl "$SRC" "$OUT"
ls -l "$OUT"; echo "finished=$(date -u +%FT%TZ)"
