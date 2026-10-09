#!/bin/bash
# CPU-partition watchdog for a running outer control: every INTERVAL seconds snapshot the 10-min slices
# and append max|w| at 2 km, min/max T_sfc and max |u_sfc| of the newest records to monitor.csv.
#   sbatch --partition=cpu-c6a24xlarge --export=ALL,RUN=<run dir>,GPU_JOB=<id> execution/tracer_mip_monitor.sbatch
#SBATCH --job-name=blab-tmip-watch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=48:00:00
#SBATCH --output=/shared/home/greg/breezelab-work/tracer-mip/logs/slurm-%j.out
set -uo pipefail
CHECKOUT=${CHECKOUT:-/shared/home/greg/breezelab-work/tracer-mip-ctrl}
W=/shared/home/greg/breezelab-work/tracer-mip; C=/shared/home/greg/breezelab-runs/20261006; D=/shared/home/greg/breezelab-runs/docs-20261007-8dc824d
CPU_TAG=$(awk -F ': ' '/model name/ {print $2; exit}' /proc/cpuinfo | tr -c 'A-Za-z0-9' '_')
source "$W/execution/depot_layout.sh"
export JULIA_DEPOT_PATH="$W/depot/$CPU_TAG:$D/depot-compute/$CPU_TAG:$D/depot-login:$C/depot:$HOME/.julia:"
export JULIA_PKG_PRECOMPILE_AUTO=0 JULIA_NUM_THREADS=4
INTERVAL=${INTERVAL:-1800}; SNAP="$RUN/../watch_snapshot_$(basename "$RUN")"; mkdir -p "$SNAP"
CSV="$RUN/monitor.csv"; [ -f "$CSV" ] || echo "wall_utc,record,t_hours,max_abs_w_2km,max_abs_w_2km_interior,min_T_sfc,max_T_sfc,max_abs_u_sfc,inner_max_abs_w,inner_max_w_z,inner_min_T,inner_min_T_z,inner_max_ncl" > "$CSV"
cd "$CHECKOUT"
while /opt/slurm/bin/squeue -h -j "$GPU_JOB" -o %t | grep -q .; do
  if [ -f "$RUN/outer_slices.jld2" ]; then
    cp "$RUN/outer_slices.jld2" "$SNAP/outer_slices.jld2" 2>/dev/null && mkdir -p "$SNAP/inner_region" && cp "$RUN/inner_region/inner_region_state.jld2" "$SNAP/inner_region/" 2>/dev/null && \
    "$C/julia-1.12.7/bin/julia" --startup-file=no --project=cases/tracer_mip -e '
      using Oceananigans, Breeze, Printf; using Oceananigans.Fields: interior; using Dates
      snap = ARGS[1]; csv = ARGS[2]; last_done = parse(Int, ARGS[3])
      w = FieldTimeSeries(joinpath(snap, "outer_slices.jld2"), "w_2km"; backend = OnDisk()); T = FieldTimeSeries(joinpath(snap, "outer_slices.jld2"), "T_sfc"; backend = OnDisk()); u = FieldTimeSeries(joinpath(snap, "outer_slices.jld2"), "u_sfc"; backend = OnDisk())
      Nx, Ny, Nz = size(w.grid); rim = [i ≤ 10 || i > Nx - 10 || j ≤ 10 || j > Ny - 10 for i in 1:Nx, j in 1:Ny]
      using Oceananigans.Grids: znode, Center
      zc = [znode(Nx ÷ 2, Ny ÷ 2, k, w.grid, Center(), Center(), Center()) for k in 1:Nz]
      wi = FieldTimeSeries(joinpath(snap, "inner_region", "inner_region_state.jld2"), "w"; backend = OnDisk()); Ti = FieldTimeSeries(joinpath(snap, "inner_region", "inner_region_state.jld2"), "T"; backend = OnDisk())
      ni_ok = haskey(Dict(n => 1 for n in ("qᶜˡ",)), "qᶜˡ")
      open(csv, "a") do io
          for n in last_done + 1:min(length(w.times), length(wi.times))
              ww = abs.(dropdims(Array(interior(w[n])); dims = 3)); TT = dropdims(Array(interior(T[n])); dims = 3); uu = abs.(dropdims(Array(interior(u[n])); dims = 3))
              w3 = abs.(Array(interior(wi[n]))); T3 = Array(interior(Ti[n])); iw = argmax(w3); iT = argmin(T3)
              @printf(io, "%s,%d,%.3f,%.3f,%.3f,%.2f,%.2f,%.2f,%.3f,%.0f,%.2f,%.0f,%.3e\n", Dates.format(now(Dates.UTC), "yyyy-mm-ddTHH:MM"), n, w.times[n] / 3600, maximum(ww), maximum(ww[.!rim]), minimum(TT), maximum(TT), maximum(uu),
                      maximum(w3), zc[min(Tuple(iw)[3], Nz)], minimum(T3), zc[min(Tuple(iT)[3], Nz)], 0.0)
          end
      end' "$SNAP" "$CSV" "$(( $(wc -l < "$CSV") - 1 ))" 2>&1 | grep -v "^\s*$" | tail -2
    tail -1 "$CSV"
  fi
  sleep "$INTERVAL"
done
echo "== $(date -u) GPU job $GPU_JOB left the queue"
