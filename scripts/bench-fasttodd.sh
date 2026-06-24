#!/usr/bin/env bash
# Benchmark phasepoly strategies (todd / tohpe / fasttodd) on the Amy/Maslov
# suite in benchmark/qc/optimized. Captures T-count + wall time per circuit per
# strategy into a CSV. Circuits run smallest-first (by file size) so partial
# results are useful early.
#
# Usage: scripts/bench-fasttodd.sh [-t SECONDS] [-o OUT.csv] [circuit-substr ...]
#   -t  per-run timeout in seconds (default 1800)
#   -o  output CSV (default docs/bench-fasttodd.csv)
#   trailing args filter circuits by substring (default: all)
set -uo pipefail

QSYN=${QSYN:-./build/qsyn}
TIMEOUT=1800
OUT=docs/bench-fasttodd.csv
CIRCDIR=benchmark/qc/optimized
STRATS=(todd tohpe fasttodd)

while getopts "t:o:" opt; do
  case $opt in
    t) TIMEOUT=$OPTARG ;;
    o) OUT=$OPTARG ;;
    *) echo "bad flag" >&2; exit 2 ;;
  esac
done
shift $((OPTIND - 1))
FILTERS=("$@")

[[ -x "$QSYN" ]] || { echo "qsyn not found/executable at $QSYN" >&2; exit 1; }

# smallest-first; only the *_pyzx.qc inputs (exclude *_pyzxtodd / *_tpar refs)
mapfile -t FILES < <(ls -Sr "$CIRCDIR"/*_pyzx.qc 2>/dev/null)

select_file() {  # keep file if it matches any filter (or no filters given)
  local f=$1
  [[ ${#FILTERS[@]} -eq 0 ]] && return 0
  local pat
  for pat in "${FILTERS[@]}"; do [[ $f == *"$pat"* ]] && return 0; done
  return 1
}

tmpd=$(mktemp -d -p "${TMPDIR:-/tmp}" benchfasttodd.XXXXXX)
trap 'rm -rf "$tmpd"' EXIT

echo "circuit,qubits,strategy,tcount,seconds,status" > "$OUT"
echo "[bench] qsyn=$QSYN timeout=${TIMEOUT}s out=$OUT" >&2

for f in "${FILES[@]}"; do
  select_file "$f" || continue
  circ=$(basename "$f" .qc); circ=${circ%_pyzx}
  for s in "${STRATS[@]}"; do
    dof="$tmpd/run.dof"
    cat > "$dof" <<EOF
qcir read $f
convert qcir tableau
tableau opt phasepoly $s
convert tableau qcir
qcir print --stat
EOF
    start=$SECONDS
    out=$(timeout "$TIMEOUT" "$QSYN" --no-version --qsynrc-path /dev/null "$dof" 2>&1)
    rc=$?
    secs=$((SECONDS - start))
    q=NA; t=NA; status=OK
    if [[ $rc -eq 124 ]]; then
      status=TIMEOUT
    elif [[ $rc -ne 0 ]]; then
      status=ERROR
    else
      line=$(printf '%s\n' "$out" | rg -o 'QCir \(.*depths\)' | tail -1)
      q=$(printf '%s' "$line" | rg -o '[0-9]+ qubits'  | rg -o '^[0-9]+')
      t=$(printf '%s' "$line" | rg -o '[0-9]+ T-gates' | rg -o '^[0-9]+')
      [[ -z $t ]] && status=PARSEFAIL
    fi
    printf '%s,%s,%s,%s,%s,%s\n' "$circ" "${q:-NA}" "$s" "${t:-NA}" "$secs" "$status" | tee -a "$OUT"
  done
done
echo "[bench] done -> $OUT" >&2
