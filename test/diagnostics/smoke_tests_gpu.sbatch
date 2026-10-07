#!/bin/bash
#SBATCH --job-name=breezelab-smoke
#SBATCH --partition=gpu-prod
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=6
#SBATCH --mem=80G
#SBATCH --time=03:00:00
#SBATCH --output=output/slurm-smoke-%j.out
# GPU smoke tests before any long run: forcing-only, dry dynamics, cloudy P3 (30-60 min).
set -euo pipefail
# Separate compiled-code cache for the compute nodes (different CPU target than the login node);
# packages and artifacts are still read from the shared depot.
export JULIA_DEPOT_PATH="$HOME/.julia-compute:$HOME/.julia:"
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/../..}"
mkdir -p output
julia --startup-file=no --project=. test/diagnostics/smoke_tests.jl gpu
