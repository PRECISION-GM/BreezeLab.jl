#!/bin/bash
# CPU-partition job: only the regional (NumericalEarth) software tests in the case sub-environment.
#SBATCH --job-name=blab-tmip-rtest
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=01:00:00
#SBATCH --output=/shared/home/greg/breezelab-work/tracer-mip/logs/slurm-%j.out
set -euo pipefail
W=/shared/home/greg/breezelab-work/tracer-mip
C=/shared/home/greg/breezelab-runs/20261006
D=/shared/home/greg/breezelab-runs/docs-20261007-8dc824d
CPU_TAG=$(awk -F ': ' '/model name/ {print $2; exit}' /proc/cpuinfo | tr -c 'A-Za-z0-9' '_')
source "$W/execution/depot_layout.sh"
export JULIA_DEPOT_PATH="$W/depot/$CPU_TAG:$D/depot-compute/$CPU_TAG:$D/depot-login:$C/depot:$HOME/.julia:"
export JULIA_PKG_PRECOMPILE_AUTO=0 JULIA_NUM_PRECOMPILE_TASKS=8 JULIA_NUM_THREADS=4
cd "$W"
echo "== $(date -u) node=$(hostname) commit=$(git rev-parse HEAD)"
"$C/julia-1.12.7/bin/julia" --startup-file=no --project=cases/tracer_mip -e 'using Pkg; Pkg.precompile("BreezeLab")
using Test, BreezeLab, Breeze, Oceananigans, Oceananigans.Units; using Dates: DateTime
test_grid(; Nx=8, Ny=8, Nz=24, Lz=6000) = RectilinearGrid(CPU(), Float64; size=(Nx, Ny, Nz), x=(0, 800), y=(0, 800), z=(0, Lz), halo=(5, 5, 5), topology=(Periodic, Periodic, Bounded))
@testset "regional" begin; include("test/tracer_mip_regional.jl"); end'
echo "== $(date -u) done"
