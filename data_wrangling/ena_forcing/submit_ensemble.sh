#!/bin/bash
# Submit every member of manifest.tsv (24 h, 128², sonde nudging + RRTMGP LW+SW) with run_ena.sbatch;
# appends member → job id to jobs.tsv.  Usage: ./submit_ensemble.sh [member-regex]
set -euo pipefail
root=$(cd "$(dirname "$0")" && pwd)
pattern=${1:-.}
tail -n +2 "$root/manifest.tsv" | while IFS=$'\t' read -r member micro rest; do
    [[ "$member" =~ $pattern ]] || continue
    input=$(echo "$rest" | awk -F'\t' '{print $(NF-1)}')
    output=$(echo "$rest" | awk -F'\t' '{print $NF}')
    mkdir -p "$output"
    cd "$output"
    id=$(BREEZELAB=$HOME/BreezeLab-ena-forcing DATA="$input" MICROPHYSICS="$micro" \
         sbatch --parsable --partition=gpua100largex4spot,gpua100largex4 --nice=10000 --job-name="ena-fe-${member%%_*}" \
         $HOME/ena-les/agentB/jobs/run_ena.sbatch r128 24 "$output" \
         --seed 1234 --closure none --scalar_advection_order 9 --dt 1 --cfl 0.8 --max_dt 5 \
         --nudging_tq 3600 --nudging_z 1200,1800 --nudging_uv 3600 --radiation rrtmgp)
    printf "%s\t%s\t%s\n" "$member" "$id" "$output" | tee -a "$root/jobs.tsv"
done
