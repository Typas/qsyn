#!/usr/bin/env bash
# Phase-3 "after" sweep: re-run ONLY the qsyn columns (pyzx + raw) with the
# ported FastTODD binary, into a SEPARATE results dir so the original baseline
# in bench2x2/results/ (the slides' numbers) is left untouched. Resumable.
set -u
ROOT=/home/liao/qsyn
VDIR=/home/liao/quantum-circuit-optimization
QSYN="$ROOT/build/qsyn"
B="$ROOT/bench2x2"
RES="$B/results-new"          # <-- distinct from baseline results/
PAPER="$ROOT/benchmark/paper-todd.csv"
CAP="${CAP:-1800}"
MAXJOBS="${MAXJOBS:-4}"
MEMKB="${MEMKB:-12582912}"
mkdir -p "$RES" "$B/sani" "$B/logs-new"

sanitize(){ sed -e 's/\bT\*/tdg/g' -e 's/\bS\*/sdg/g' -e 's/\bZ\*/z/g' -e 's/^cnot /tof /' "$1" > "$2"; }
qubits_of(){ rg -m1 '^\.v ' "$1" 2>/dev/null | awk '{print NF-1}'; }

run_qsyn(){ # circuit run infile
  local c="$1" run="$2" in="$3" row="$RES/$1__$2.row"
  [ -f "$row" ] && return
  local n sani log start sec rc st T= note=
  n=$(qubits_of "$in"); sani="$B/sani/${c}__${run}.qc"; log="$B/logs-new/${c}__${run}.log"
  sanitize "$in" "$sani"; start=$(date +%s)
  ( ulimit -v "$MEMKB" 2>/dev/null; exec timeout -k 10 "$CAP" "$QSYN" -q --no-version -c \
      "qcir read $sani; convert qcir tableau; tableau optimize tmerge; tableau optimize hopt; tableau optimize phasepoly fasttodd; convert tableau qcir; qcir print --stat; quit -f" ) >"$log" 2>&1
  rc=$?; sec=$(( $(date +%s) - start ))
  if [ $rc -eq 124 ] || [ $rc -eq 137 ]; then st=TLE
  elif rg -qi 'bad_alloc|out of memory|std::length_error' "$log"; then st=MLE
  elif [ $rc -ne 0 ] || rg -qi 'error' "$log"; then st=ERR; note="rc=$rc"
  else st=ok; T=$(rg -o 'T-family\s*:\s*[0-9]+' "$log" | rg -o '[0-9]+$' | tail -1); fi
  printf '%s,%s,%s,%s,%s,%s,%s\n' "$c" "$run" "${n:-?}" "${T:-NA}" "$sec" "$st" "$note" >"$row"
  echo "[$run] $c -> ${T:-NA}T ${sec}s $st"
}

throttle(){ while [ "$(jobs -rp | wc -l)" -ge "$MAXJOBS" ]; do wait -n; done; }

echo "[rerun] start $(date '+%H:%M:%S')  cap=${CAP}s maxjobs=$MAXJOBS memcap=$((MEMKB/1024/1024))GB"
while IFS=, read -r c _g _tier _n _hl _hp _ht _htt _vp _vt _vtt _vf file _notes; do
  [ "$c" = "circuit" ] && continue
  pyzx="$ROOT/$file"
  raw="$VDIR/circuits/inputs/$c.qc"; [ "$c" = "DEFAULT" ] && raw="$VDIR/circuits/inputs/default.qc"
  [ -n "$file" ] && [ -f "$pyzx" ] && { run_qsyn "$c" qsyn_pyzx "$pyzx" & throttle; }
  [ -f "$raw" ] && { run_qsyn "$c" qsyn_raw "$raw" & throttle; }
done < "$PAPER"
wait
{ echo "circuit,run,n,T,seconds,status,note"; cat "$RES"/*.row 2>/dev/null | sort; } > "$B/all-new.csv"
echo "[rerun] DONE $(date '+%H:%M:%S') -> $B/all-new.csv"
