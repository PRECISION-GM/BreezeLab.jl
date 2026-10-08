#!/bin/bash
# CPU-partition job: instantiate + precompile the case sub-environment, stage ERA5 (credentials from
# ~/.bashrc, never printed), run the regional software tests. No GPU involved.
#SBATCH --job-name=blab-tmip-prep
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=80G
#SBATCH --time=04:00:00
#SBATCH --output=/shared/home/greg/breezelab-work/tracer-mip/logs/slurm-%j.out
set -euo pipefail
W=/shared/home/greg/breezelab-work/tracer-mip
C=/shared/home/greg/breezelab-runs/20261006
D=/shared/home/greg/breezelab-runs/docs-20261007-8dc824d
CPU_TAG=$(awk -F ': ' '/model name/ {print $2; exit}' /proc/cpuinfo | tr -c 'A-Za-z0-9' '_')
mkdir -p "$W/depot/$CPU_TAG" "$W/data/era5"
export JULIA_DEPOT_PATH="$W/depot/$CPU_TAG:$D/depot-compute/$CPU_TAG:$D/depot-login:$C/depot:$HOME/.julia:"
export JULIA_PKG_PRECOMPILE_AUTO=0 JULIA_NUM_PRECOMPILE_TASKS=20 JULIA_NUM_THREADS=8
set -a; source "$HOME/.bashrc" 2>/dev/null || true; set +a
J="$C/julia-1.12.7/bin/julia --startup-file=no"
cd "$W"
echo "== $(date -u) node=$(hostname) cpu_tag=$CPU_TAG commit=$(git rev-parse HEAD)"
echo "== instantiate case sub-environment"
$J --project=cases/tracer_mip -e 'using Pkg; Pkg.instantiate(); Pkg.precompile(; timing=true); Pkg.status()'
echo "== $(date -u) stage ERA5 aug07"
$J --project=cases/tracer_mip data_wrangling/fetch_era5_tracer_mip.jl aug07 "$W/data/era5" 6 2>&1 | grep --line-buffered -v -i "key:\|token" || true
echo "== $(date -u) regional software tests (synthetic parent)"
$J --project=cases/tracer_mip -e 'using Test, BreezeLab, Breeze, Oceananigans, Oceananigans.Units; using Dates: DateTime
test_grid(; Nx=8, Ny=8, Nz=24, Lz=6000) = RectilinearGrid(CPU(), Float64; size=(Nx, Ny, Nz), x=(0, 800), y=(0, 800), z=(0, Lz), halo=(5, 5, 5), topology=(Periodic, Periodic, Bounded))
@testset "regional" begin; include("test/tracer_mip_regional.jl"); end'
echo "== $(date -u) stage ERA5 jun17"
$J --project=cases/tracer_mip data_wrangling/fetch_era5_tracer_mip.jl jun17 "$W/data/era5" 6 2>&1 | grep --line-buffered -v -i "key:\|token" || true
echo "== $(date -u) done"
