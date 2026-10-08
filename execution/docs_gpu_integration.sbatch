#!/bin/bash
#SBATCH --job-name=blab-docs-gpu
#SBATCH --nodes=1 --ntasks=1 --gres=gpu:1 --cpus-per-task=12 --mem=96G --time=06:00:00
#SBATCH --output=/shared/home/greg/breezelab-work/integration/logs/slurm-%j.out
# Full GPU documentation build of the merged main (executes the ENA example on the GPU, renders all case pages).
# Usage: sbatch --partition=<p> --nodelist=<node> execution/docs_gpu_integration.sbatch <GPU UUID>
set -euo pipefail
W=/shared/home/greg/breezelab-work/integration; C=/shared/home/greg/breezelab-runs/20261006; D=/shared/home/greg/breezelab-runs/docs-20261007-8dc824d
GPU_UUID=${1:?GPU UUID}
LOCKS=/shared/home/greg/breezelab-work/gpu-locks; mkdir -p "$LOCKS"
exec 9>"$LOCKS/$GPU_UUID.lock"; flock -n 9 || { echo "GPU $GPU_UUID locked"; exit 1; }
[ -z "$(nvidia-smi -i "$GPU_UUID" --query-compute-apps=pid --format=csv,noheader)" ] || { echo 'Selected GPU is busy'; exit 1; }
export CUDA_VISIBLE_DEVICES="$GPU_UUID"
GPU_NAME=$(nvidia-smi -i "$GPU_UUID" --query-gpu=name --format=csv,noheader)
case "$GPU_NAME" in *H100*|*A100*) ;; *) echo "Unexpected GPU: $GPU_NAME"; exit 1 ;; esac
CPU_TAG=$(awk -F ': ' '/model name/ {print $2; exit}' /proc/cpuinfo | tr -c 'A-Za-z0-9' '_')
mkdir -p "$W/depot/$CPU_TAG"
export JULIA_DEPOT_PATH="$W/depot/$CPU_TAG:/shared/home/greg/breezelab-work/tracer-mip/depot/$CPU_TAG:/shared/home/greg/breezelab-work/tracer-dp-scream/depot/$CPU_TAG:/shared/home/greg/breezelab-work/sea-starr/depot/$CPU_TAG:$D/depot-compute/$CPU_TAG:$D/depot-login:$C/depot:$HOME/.julia:"
export JULIA_PKG_PRECOMPILE_AUTO=0 JULIA_NUM_PRECOMPILE_TASKS=6 JULIA_NUM_THREADS=6
cd $W; echo "commit=$(git rev-parse HEAD) node=$(hostname) gpu=$GPU_NAME started=$(date -u +%FT%TZ)"
[ -d data/covert2022_bin ] || "$C/julia-1.12.7/bin/julia" --startup-file=no --project=. data_wrangling/fetch_covert_inputs.jl 2>&1 | tail -2
"$C/julia-1.12.7/bin/julia" --startup-file=no docs/setup.jl 2>&1 | tail -2
"$C/julia-1.12.7/bin/julia" --startup-file=no --project=docs docs/make.jl
python3 - <<'PY'
import datetime, json, os, subprocess, pathlib
w = pathlib.Path('/shared/home/greg/breezelab-work/integration')
meta = {'source_commit': subprocess.check_output(['git','rev-parse','HEAD'], text=True, cwd=w).strip(),
        'job_id': os.environ['SLURM_JOB_ID'], 'gpu': os.environ.get('GPU_NAME', subprocess.check_output(['nvidia-smi','--query-gpu=name','--format=csv,noheader'], text=True).strip()),
        'partition': os.environ['SLURM_JOB_PARTITION'], 'completed_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'mode': 'full GPU example + all case pages', 'stop_time_seconds': 21600, 'grid': [256, 256, 192], 'microphysics': 'p3_n75'}
(w / 'gpu-build.json').write_text(json.dumps(meta, indent=2) + '\n')
PY
echo "finished=$(date -u +%FT%TZ)"; touch "$W/DOCS_COMPLETE"
