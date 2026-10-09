#!/bin/bash
# Fetch TRACER AMF1 (La Porte, "hou" M1) observations for the August 7 2022 case through the ARM Live Data API.
#   set -a; source ~/.bashrc; set +a     # provides ARM_USERNAME / ARM_TOKEN (never printed)
#   bash data_wrangling/fetch_arm_tracer_obs.sh <out dir> [start=2022-08-07] [end=2022-08-08]
# Writes the files plus MANIFEST.txt (datastream, file, bytes, sha256). Default datastreams: surface met
# (houmetM1.b1), ceilometer (houceilM1.b1), ARSCL cloud boundaries (houarsclkazr1kolliasM1.c0), laser
# disdrometer (houldM1.b1), radiosondes (housondewnpnM1.b1). Per-file size cap 2 GB.
set -euo pipefail
OUT=${1:?output directory}; START=${2:-2022-08-07}; END=${3:-2022-08-08}
DATASTREAMS=${DATASTREAMS:-"houmetM1.b1 houceilM1.b1 houarsclkazr1kolliasM1.c0 houldM1.b1 housondewnpnM1.b1"}
[ -n "${ARM_USERNAME:-}" ] && [ -n "${ARM_TOKEN:-}" ] || { echo "ARM_USERNAME/ARM_TOKEN are not set (source your secret configuration)"; exit 2; }
mkdir -p "$OUT"; cd "$OUT"
MAN="$OUT/MANIFEST.txt"; echo "# ARM Live API, fetched $(date -u +%FT%TZ), window $START..$END" > "$MAN"
AUTH="$ARM_USERNAME:$ARM_TOKEN"
for DS in $DATASTREAMS; do
  LIST=$(curl -sS --max-time 120 "https://adc.arm.gov/armlive/livedata/query?user=$AUTH&ds=$DS&start=$START&end=$END&wt=json" | python3 -c 'import json,sys; d=json.load(sys.stdin); print("\n".join(d.get("files",[])))' 2>/dev/null || true)
  if [ -z "$LIST" ]; then echo "$DS: no files listed"; echo "$DS NONE" >> "$MAN"; continue; fi
  for F in $LIST; do
    [ -f "$F" ] || curl -sS --max-time 1800 --max-filesize 2000000000 -o "$F" "https://adc.arm.gov/armlive/livedata/saveData?user=$AUTH&file=$F" || { echo "$DS: download failed for $F"; rm -f "$F"; continue; }
    echo "$DS $F $(stat -c %s "$F") $(sha256sum "$F" | cut -d' ' -f1)" >> "$MAN"
  done
  echo "$DS: $(echo "$LIST" | wc -l) files"
done
echo "manifest: $MAN"
