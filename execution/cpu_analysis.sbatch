#!/bin/bash
#SBATCH --job-name=blab-ena-lasso-cpu
#SBATCH --partition=cpu-c6a24xlarge
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --output=/shared/home/greg/breezelab-work/ena-lasso/logs/slurm-%j.out
# CPU-node Julia analysis for the ena-lasso checkout: sbatch execution/cpu_analysis.sbatch <script.jl> [args...]
set -euo pipefail
W=/shared/home/greg/breezelab-work/ena-lasso
cd "$W"
source execution/wpcluster_env.sh
export JULIA_NUM_THREADS=4 JULIA_NUM_PRECOMPILE_TASKS=8
echo "node $(hostname) cpu tag $CPU_TAG depot $JULIA_DEPOT_PATH"
"$JULIA" --startup-file=no --project=. "$@"
