#!/bin/bash
#SBATCH --job-name=blab-tmip-arm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=03:00:00
#SBATCH --output=/shared/home/greg/breezelab-work/tracer-mip/logs/slurm-%j.out
set -euo pipefail
set -a; source "$HOME/.bashrc" 2>/dev/null || true; set +a
cd "${CHECKOUT:-/shared/home/greg/breezelab-work/tracer-mip-ctrl}"
bash data_wrangling/fetch_arm_tracer_obs.sh "${OUT:?}" 2>&1 | grep --line-buffered -v -i "token\|user="
