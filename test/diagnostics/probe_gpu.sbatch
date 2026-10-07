#!/bin/bash
#SBATCH --job-name=p3probe
#SBATCH --partition=gpu-prod
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --mem=80G
#SBATCH --time=02:00:00
#SBATCH --output=output/slurm-probe-%j.out
# GPU stability probe (see test/diagnostics/p3_stability_probe.jl for the arguments), run from the
# submitting checkout so the ocmain worktree can be probed with the same script:
#   sbatch test/diagnostics/probe_gpu.sbatch Float32 0.5 45 p3_aer2
set -euo pipefail
export JULIA_DEPOT_PATH="$HOME/.julia-compute:$HOME/.julia:"
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/../..}"
mkdir -p output
echo "checkout: $(pwd) $(git rev-parse --short HEAD 2>/dev/null) | args: $*"
julia --startup-file=no --project=. test/diagnostics/p3_stability_probe.jl gpu "$@"
