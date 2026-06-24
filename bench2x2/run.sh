#!/usr/bin/env bash
# 2x2 FastTODD benchmark: {qsyn, Vandaele} x {pyzx, raw inputs}.
# Persistent + resumable: each (circuit,run) writes results/<c>__<run>.row;
# re-running skips existing rows. Per-process memory cap (ulimit -v) so a
# pathological circuit is killed (MLE) instead of OOM-ing the machine.
set -u
ROOT=/home/liao/qsyn
VDIR=/home/liao/quantum-circuit-optimization
BIN="$VDIR/target/release/quantum_circuit_optimization"
QSYN="$ROOT/build/qsyn"
B="$ROOT/bench2x2"
PAPER="$ROOT/benchmark/paper-todd.csv"
CAP="${CAP:-1800}"          # per-circuit wall cap (s)
MAXJOBS="${MAXJOBS:-4}"     # global concurrency (8 cores, 54GB)
MEMKB="${MEMKB:-12582912}"  # 12 GB virtual-memory cap per process
mkdir -p "$B/results" "$B/sani" "$B/logs"

sanitize(){ sed -e 's/\bT\*/tdg/g' -e 's/\bS\*/sdg/g' -e 's/\bZ\*/z/g' -e 's/^cnot /tof /' "$1" > "$2"; }
qubits_of(){ rg -m1 '^\.v ' "$1" 2>/dev/null | awk '{print NF-1}'; }

run_qsyn(){ # circuit run infile
  local c="$1" run="$2" in="$3" row="$B/results/$1__$2.row"
  [ -f "$row" ] && return
  local n sani log start sec rc st T= peak= note=
  n=$(qubits_of "$in"); sani="$B/sani/${c}__${run}.qc"; log="$B/logs/${c}__${run}.log"
  sanitize "$in" "$sani"; start=$(date +%s)
  ( ulimit -v "$MEMKB" 2>/dev/null; exec /usr/bin/time -v timeout -k 10 "$CAP" "$QSYN" -q --no-version -c \
      "qcir read $sani; convert qcir tableau; tableau optimize tmerge; tableau optimize hopt; tableau optimize phasepoly fasttodd; convert tableau qcir; qcir print --stat; quit -f" ) >"$log" 2>&1
  rc=$?; sec=$(( $(date +%s) - start ))
  peak=$(rg -o 'Maximum resident set size .kbytes.: [0-9]+' "$log" | rg -o '[0-9]+$' | tail -1)
  if [ $rc -eq 124 ] || [ $rc -eq 137 ]; then st=TLE
  elif rg -qi 'bad_alloc|out of memory|std::length_error' "$log"; then st=MLE
  elif [ $rc -ne 0 ] || rg -qi 'error' "$log"; then st=ERR; note="rc=$rc"
  else st=ok; T=$(rg -o 'T-family\s*:\s*[0-9]+' "$log" | rg -o '[0-9]+$' | tail -1); fi
  printf '%s,%s,%s,%s,%s,%s,%s,%s\n' "$c" "$run" "${n:-?}" "${T:-NA}" "$sec" "${peak:-NA}" "$st" "$note" >"$row"
  echo "[$run] $c -> ${T:-NA}T ${sec}s $st"
}

run_vand(){ # circuit run infile
  local c="$1" run="$2" in="$3" row="$B/results/$1__$2.row"
  [ -f "$row" ] && return
  local n log start sec rc st T= peak= note=
  n=$(qubits_of "$in"); log="$B/logs/${c}__${run}.log"; start=$(date +%s)
  ( ulimit -v "$MEMKB" 2>/dev/null; exec /usr/bin/time -v timeout -k 10 "$CAP" env -C "$VDIR" "$BIN" "$in" ) >"$log" 2>&1
  rc=$?; sec=$(( $(date +%s) - start ))
  peak=$(rg -o 'Maximum resident set size .kbytes.: [0-9]+' "$log" | rg -o '[0-9]+$' | tail -1)
  if rg -qi 'not implemented' "$log"; then st=NA; note=$(rg -i 'not implemented' "$log" | head -1)
  elif [ $rc -eq 124 ] || [ $rc -eq 137 ]; then st=TLE
  elif rg -qi 'bad_alloc|memory allocation' "$log"; then st=MLE
  elif [ $rc -ne 0 ]; then st=ERR; note="rc=$rc"
  else st=ok; T=$(rg -o 'T-count:\s*[0-9]+' "$log" | rg -o '[0-9]+$' | tail -1); fi
  printf '%s,%s,%s,%s,%s,%s,%s,%s\n' "$c" "$run" "${n:-?}" "${T:-NA}" "$sec" "${peak:-NA}" "$st" "$note" >"$row"
  echo "[$run] $c -> ${T:-NA}T ${sec}s $st"
}

throttle(){ while [ "$(jobs -rp | wc -l)" -ge "$MAXJOBS" ]; do wait -n; done; }

echo "[run] start $(date '+%H:%M:%S')  cap=${CAP}s maxjobs=$MAXJOBS memcap=$((MEMKB/1024/1024))GB"
while IFS=, read -r c _g _tier _n _hl _hp _ht _htt _vp _vt _vtt _vf file _notes; do
  [ "$c" = "circuit" ] && continue
  pyzx="$ROOT/$file"
  raw="$VDIR/circuits/inputs/$c.qc"; [ "$c" = "DEFAULT" ] && raw="$VDIR/circuits/inputs/default.qc"
  if [ -n "$file" ] && [ -f "$pyzx" ]; then
    run_qsyn "$c" qsyn_pyzx "$pyzx" & throttle
    run_vand "$c" vand_pyzx "$pyzx" & throttle
  fi
  if [ -f "$raw" ]; then
    run_qsyn "$c" qsyn_raw "$raw" & throttle
    run_vand "$c" vand_raw "$raw" & throttle
  fi
done < "$PAPER"
wait
{ echo "circuit,run,n,T,seconds,peak_kb,status,note"; cat "$B"/results/*.row 2>/dev/null | sort; } > "$B/all.csv"
echo "[run] DONE $(date '+%H:%M:%S') -> $B/all.csv"
