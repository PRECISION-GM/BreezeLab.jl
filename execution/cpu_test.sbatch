#!/bin/bash
#SBATCH --job-name=blab-cpu-test
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/shared/home/greg/breezelab-work/tracer-dp-scream/logs/slurm-%j.out
# CPU-only precompile + Pkg.test() for the tracer-dp-scream checkout on a compute node whose CPU
# model matches the A100 nodes (AMD EPYC 7R13), so the compiled depot is reused by GPU jobs.
# Usage: sbatch --partition=cpu-c6a24xlarge execution/cpu_test.sbatch [fixtures]
set -euo pipefail
W=/shared/home/greg/breezelab-work/tracer-dp-scream
C=/shared/home/greg/breezelab-runs/20261006
D=/shared/home/greg/breezelab-runs/docs-20261007-8dc824d
CPU_TAG=$(awk -F ': ' '/model name/ {print $2; exit}' /proc/cpuinfo | tr -c 'A-Za-z0-9' '_')
echo "node=$(hostname) cpu_tag=$CPU_TAG nproc=$(nproc) started=$(date -u +%FT%TZ) commit=$(git -C $W rev-parse HEAD)"
mkdir -p "$W/depot/$CPU_TAG"
export JULIA_DEPOT_PATH="$W/depot/$CPU_TAG:$D/depot-compute/$CPU_TAG:$D/depot-login:$C/depot:$HOME/.julia:"
export JULIA_PKG_PRECOMPILE_AUTO=0
export JULIA_NUM_PRECOMPILE_TASKS=12
export JULIA_NUM_THREADS=8
export JULIA_CUDA_USE_COMPAT=false
cd "$W"
J="$C/julia-1.12.7/bin/julia"
time "$J" --startup-file=no --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
if [ "${1:-}" = "fixtures" ]; then
  "$J" --startup-file=no --project=. data_wrangling/make_tracer_fixtures.jl
  ls -l test/fixtures/*.nc
fi
time "$J" --startup-file=no --project=. -e 'using Pkg; Pkg.test(; allow_reresolve=false)'
echo "finished=$(date -u +%FT%TZ)"
