#!/usr/bin/env bash
# Per-STEP time + peak-RSS for the gadgetize pipeline, using qsyn's `usage -t -m`
# after each command. `usage -t` prints that command's period CPU time (then resets);
# `usage -m` prints cumulative peak RSS (ru_maxrss is monotonic), so the delta between
# consecutive steps = the memory that step added to the peak.
#
# Capped (`ulimit -v`) + timed so it can never host-OOM (cf. 2026-06-23 incident).
#
# Usage: scripts/calib-steps.sh <circuit> [more...]
# Env: QSYN, VDIR, CAP_KB (default 8G), TMO (default 120s), OUTDIR (default bench2x2/calib-steps)
set -u
QSYN="${QSYN:-build/qsyn}"
VDIR="${VDIR:-/home/liao/quantum-circuit-optimization}"
CAP_KB="${CAP_KB:-8388608}"
TMO="${TMO:-120}"
OUTDIR="${OUTDIR:-bench2x2/calib-steps}"
mkdir -p "$OUTDIR/sani" "$OUTDIR/logs"
sanitize(){ sed -e 's/\bT\*/tdg/g' -e 's/\bS\*/sdg/g' -e 's/\bZ\*/z/g' -e 's/^cnot /tof /' "$1" > "$2"; }

# pipeline steps, in order (must match the dofile 1:1; phasepoly dropped when NO_PHASEPOLY=1)
if [ "${NO_PHASEPOLY:-0}" = 1 ]; then
  STEP_NAMES=(read convert_in tmerge hopt gadgetize convert_out)
else
  STEP_NAMES=(read convert_in tmerge hopt gadgetize phasepoly convert_out)
fi

for c in "$@"; do
  in="$VDIR/circuits/inputs/$c.qc"
  [ -f "$in" ] || { echo "$c: MISSING"; continue; }
  sani="$OUTDIR/sani/$c.qc"; sanitize "$in" "$sani"
  log="$OUTDIR/logs/$c.steps.log"

  # bare `usage` reports BOTH period time + peak memory (-t/-m are mutually exclusive).
  # single-line -c (multi-line strings don't parse). NO_PHASEPOLY=1 skips phasepoly+convert_out
  # (use for large circuits whose phasepoly TLEs, to isolate the gadgetize sub-steps).
  if [ "${NO_PHASEPOLY:-0}" = 1 ]; then
    cmd="qcir read $sani; usage; convert qcir tableau; usage; tableau optimize tmerge; usage; tableau optimize hopt; usage; tableau optimize gadgetize; usage; convert tableau qcir; usage; quit -f"
  else
    cmd="qcir read $sani; usage; convert qcir tableau; usage; tableau optimize tmerge; usage; tableau optimize hopt; usage; tableau optimize gadgetize; usage; tableau optimize phasepoly fasttodd; usage; convert tableau qcir; usage; quit -f"
  fi
  ( ulimit -v "$CAP_KB"; timeout -k5 "$TMO" "$QSYN" -q --no-version -c "$cmd" ) > "$log" 2>&1
  rc=$?

  # k-th 'Period time used' / 'Total memory used' = step k (positional, in command order)
  mapfile -t T < <(rg -o 'Period time used : [0-9.]+' "$log" | rg -o '[0-9.]+$')
  mapfile -t MEM < <(rg -o 'Total memory used: [0-9.]+' "$log" | rg -o '[0-9.]+')

  echo "### $c  (rc=$rc, cap=$((CAP_KB/1024/1024))G)"
  printf '%-12s %10s %12s %12s\n' step cpu_s peak_MiB add_MiB
  prev=0
  for i in "${!STEP_NAMES[@]}"; do
    t="${T[$i]:-NA}"; mem="${MEM[$i]:-NA}"
    if [ "$mem" != "NA" ]; then add=$(awk -v a="$mem" -v b="$prev" 'BEGIN{printf "%.2f", a-b}'); prev="$mem"; else add=NA; fi
    printf '%-12s %10s %12s %12s\n' "${STEP_NAMES[$i]}" "$t" "$mem" "$add"
  done
  echo
done
