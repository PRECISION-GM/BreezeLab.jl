#!/bin/bash
#SBATCH --job-name=breezelab
#SBATCH --partition=gpu-prod
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --mem=80G
#SBATCH --time=12:00:00
#SBATCH --output=output/slurm-%j.out
# Usage (from repository root):
#   sbatch --partition=gpu-prod execution/submit_gpu.sbatch cases/ena_covert.jl
#   sbatch --partition=gpu-p4de-2c execution/submit_gpu.sbatch cases/cli/run_case.jl --arch gpu ...
set -euo pipefail
if [ "$#" -eq 0 ]; then
    echo "Usage: sbatch execution/submit_gpu.sbatch <case.jl> [arguments...]" >&2
    exit 2
fi
# Separate compiled-code cache for the compute nodes (different CPU target than the login node);
# packages and artifacts are still read from the shared depot.
export JULIA_DEPOT_PATH="$HOME/.julia-compute:$HOME/.julia:"
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/..}"
mkdir -p output
export JULIA_CUDA_USE_COMPAT=false
julia --startup-file=no --project=. "$@"
