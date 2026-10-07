#!/bin/bash
# CPU-partition job: precompile the package environment and run Pkg.test() (CPU suite).
#SBATCH --job-name=blab-tmip-test
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
mkdir -p "$W/depot/$CPU_TAG"
export JULIA_DEPOT_PATH="$W/depot/$CPU_TAG:$D/depot-compute/$CPU_TAG:$D/depot-login:$C/depot:$HOME/.julia:"
export JULIA_PKG_PRECOMPILE_AUTO=0 JULIA_NUM_PRECOMPILE_TASKS=20 JULIA_NUM_THREADS=8
J="$C/julia-1.12.7/bin/julia --startup-file=no"
cd "$W"
echo "== $(date -u) node=$(hostname) cpu_tag=$CPU_TAG commit=$(git rev-parse HEAD)"
[ -d data/covert2022_bin ] || cp -r /shared/home/greg/breezelab-runs/docs-20261007-8dc824d/source/data/covert2022_bin data/ 2>/dev/null || echo "no Covert inputs staged: ENA data tests will be skipped"
$J --project=. -e 'using Pkg; Pkg.precompile(; timing=true)'
echo "== $(date -u) Pkg.test"
$J --project=. -e 'using Pkg; Pkg.test(; allow_reresolve=false)'
echo "== $(date -u) done"
