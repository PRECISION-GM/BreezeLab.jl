#!/bin/bash
# The staged Covert-public-bin development benchmark runs (1M-control → P3-N75 → P3-aer2),
# full 256×256×192 grid at 35 m, 6 h (06-12 UTC 18 July 2017), one A100 each.
set -euo pipefail
cd "$(dirname "$0")/.."
# MODE=fixed (default): the preset's fixed dt = 0.5 s (SAM fidelity, ~2-3x more steps).
# MODE=adaptive: labelled override, CFL wizard up to --max_dt (recorded in provenance).
# FLOAT (default Float32, the preset precision). A single Float32/Float64 pair of chaotic LES
# probes differed late in the run; that is not evidence of a precision artefact, so both
# precisions are run as a controlled pair (FLOAT=Float64 for the second member).
# PARTITION=gpu-prod (H100) or gpu-p4de-2c (A100) on wpcluster; per-member overrides
# PARTITION_1M / PARTITION_N75 / PARTITION_AER2 select each partition.
# TIME is a conservative allocation ceiling, not a measured runtime prediction.
# RUNS selects members; each requests one GPU with CPUS and MEM.
# GRID=lasso runs on the 260-level LASSO grid (Δz = 25 m to 6 km) instead of the 192-level Covert
# grid of the prm/grd files (output suffix _lassogrid); MOMENTS=positive selects the positivity-only
# limiter for the number/volume moments (suffix _posmom); SLICE_INTERVAL=<seconds> overrides the
# 10-minute slice cadence (suffix _s<seconds>), e.g. 60 for smooth animations.
MODE=${MODE:-fixed}
FLOAT=${FLOAT:-Float32}
PARTITION=${PARTITION:-gpu-prod}
TIME=${TIME:-72:00:00}
RUNS=${RUNS:-"one_moment p3_n75 p3_aer2"}
common="--data data/covert2022_bin --preset covert_public_bin --Nx 256 --Ny 256 --float $FLOAT"
if [ "$MODE" = "adaptive" ]; then
    common="$common --max_dt ${MAX_DT:-1.5}"
    suffix="_adaptive"
else
    suffix=""
fi
[ "$FLOAT" = "Float64" ] && suffix="${suffix}_f64"
if [ "${GRID:-covert}" = "lasso" ]; then common="$common --lasso_grid true"; suffix="${suffix}_lassogrid"; fi
if [ "${MOMENTS:-plain}" = "positive" ]; then common="$common --moment_advection positive"; suffix="${suffix}_posmom"; fi
if [ -n "${SLICE_INTERVAL:-}" ]; then common="$common --slice_interval $SLICE_INTERVAL"; suffix="${suffix}_s${SLICE_INTERVAL}"; fi
# FORMULATION=StaticEnergy reproduces the pre-7-September runs (suffix _s); the default is
# LiquidIcePotentialTemperature. TAG=<text> appends a free suffix to the output directory.
if [ -n "${FORMULATION:-}" ]; then common="$common --formulation $FORMULATION"; [ "$FORMULATION" = "StaticEnergy" ] && suffix="${suffix}_s"; fi
if [ -n "${TAG:-}" ]; then suffix="${suffix}_${TAG}"; fi
submit() {  # submit <partition> <job name> <microphysics> [extra run_case options]
    local partition=$1 name=$2 microphysics=$3; shift 3
    sbatch --partition="$partition" --time="$TIME" --gres=gpu:1 --cpus-per-task="${CPUS:-4}" --mem="${MEM:-100G}" --job-name="$name" \
        execution/submit_gpu.sbatch cases/cli/run_case.jl --arch gpu $common --microphysics "$microphysics" "$@" --output "output/covert_public_bin_$microphysics$suffix"
}
for run in $RUNS; do
    case $run in
        one_moment) submit "${PARTITION_1M:-$PARTITION}"   lasso-1m   one_moment ;;
        p3_n75)     submit "${PARTITION_N75:-$PARTITION}"  lasso-n75  p3_n75 ;;
        p3_aer2)    submit "${PARTITION_AER2:-$PARTITION}" lasso-aer2 p3_aer2 --aerosol_replenishment diagnostic_ccn --aerosol_ss_cap 0.003 ;;
        p3_covert_n75) submit "${PARTITION_AER2:-$PARTITION}" covert-n75a p3_covert_n75 --aerosol_replenishment diagnostic_ccn ;;
        *) echo "unknown run $run" >&2; exit 1 ;;
    esac
done
