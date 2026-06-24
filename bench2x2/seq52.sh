#!/usr/bin/env bash
# 52 GB sequential re-run of the 3 Vandaele MLE circuits (DEFAULT, hwb10, hwb11)
# on the AVX2 build, under /usr/bin/time -v, to see if more memory lets them
# complete and to capture peak RSS. 1800s time cap each (matches the sweep).
VDIR=/home/liao/quantum-circuit-optimization
BIN="$VDIR/target/release/quantum_circuit_optimization"
B=/home/liao/qsyn/bench2x2
for c in hwb10 hwb11 default; do
  log="$B/logs/${c}__52gb.log"
  echo "[52gb] start $c $(date '+%H:%M:%S')"
  ( ulimit -v 54525952; /usr/bin/time -v timeout -k 10 1800 env -C "$VDIR" "$BIN" "$VDIR/circuits/inputs/$c.qc" ) >"$log" 2>&1
  rc=$?
  if rg -qi 'allocation.*failed|dumped core' "$log"; then
    res="MLE(>52GB)"
  elif [ $rc -eq 124 ] || [ $rc -eq 137 ]; then
    res="TLE(1800s)"
  else
    t=$(rg -o 'T-count: [0-9]+' "$log" | tail -1)
    peak=$(rg -o 'Maximum resident set size .kbytes.: [0-9]+' "$log" | rg -o '[0-9]+$')
    res="$t  peakRSS=$((peak/1024/1024))GB"
  fi
  echo "[52gb] $c -> $res"
done
echo "[52gb] ALL DONE"
