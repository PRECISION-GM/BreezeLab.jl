# Source this on the wpcluster login node or inside a Slurm job (breezelab-work/COORDINATION.md Julia/depot
# setup for the ena-lasso checkout): `source execution/wpcluster_env.sh`; then `$JULIA --project ...`.
export JULIA=/shared/home/greg/breezelab-runs/20261006/julia-1.12.7/bin/julia
CPU_TAG=$(awk -F ': ' '/model name/ {print $2; exit}' /proc/cpuinfo | tr -c 'A-Za-z0-9' '_')
mkdir -p "/shared/home/greg/breezelab-work/ena-lasso/depot/$CPU_TAG"
export JULIA_DEPOT_PATH="/shared/home/greg/breezelab-work/ena-lasso/depot/$CPU_TAG:/shared/home/greg/breezelab-runs/docs-20261007-8dc824d/depot-compute/$CPU_TAG:/shared/home/greg/breezelab-runs/docs-20261007-8dc824d/depot-login:/shared/home/greg/breezelab-runs/20261006/depot:$HOME/.julia:"
export JULIA_PKG_PRECOMPILE_AUTO=0
export JULIA_NUM_PRECOMPILE_TASKS=2
export JULIA_NUM_THREADS=2
export JULIA_PKG_OFFLINE=true
