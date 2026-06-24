#!/usr/bin/env bash
# Reproducible T-count / time / memory table for the gadgetize budget modes, for ONE circuit.
#
# Each mode runs the full chain  read -> convert -> tmerge -> hopt -> gadgetize[mode]
#   -> phasepoly fasttodd -> count T, under /usr/bin/time, and reports final T-count, wall time, peak RSS.
# Modes exceeding the time cap are marked TLE.
#
# Runs are SEQUENTIAL and isolated on purpose: parallel runs contend for CPU/RAM and make wall-time
# non-reproducible. Budgets are fixed integers so the table is deterministic. The `auto` row uses
# `--adaptive` (= detected available RAM) and is therefore MACHINE-DEPENDENT -- it is opt-in (RUN_AUTO=1)
# and labelled as not cross-machine reproducible.
#
# Portable: no absolute paths, circuit from this repo's benchmark/qc, work files in a mktemp dir.
# Needs /usr/bin/time and rg. Override binary with QSYN=/path/to/qsyn.
#
# Usage:   scripts/gadget-tcount-table.sh <circuit> <mode> [mode ...]
#   a mode is exactly one of:  none (no-gadgetize) | <size> e.g. 64M/1G (--memory-limit)
#                              | auto (--adaptive, MACHINE-DEPENDENT) | full (gadgetize all)
#   It runs ONLY the modes you list -- nothing is forced in.
# Examples:
#   scripts/gadget-tcount-table.sh hwb8 1G                       # just hwb8 @ 1G
#   scripts/gadget-tcount-table.sh hwb6 64M 128M 256M 512M auto  # hwb6 @ those budgets + adaptive
#   scripts/gadget-tcount-table.sh hwb8 none 64M 256M 1G full    # baseline + budgets + full
# Env: TMO (TLE cap seconds, default 1800)
set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
QSYN="${QSYN:-$REPO/build/qsyn}"; TMO="${TMO:-1800}"
WORK="$(mktemp -d)"; trap 'rm -rf "$WORK"' EXIT
[ -x "$QSYN" ] || { echo "qsyn not found at $QSYN (build it or set QSYN=)" >&2; exit 1; }

circuit_path() { case "$1" in
  tof3) echo benchmark/qc/optimized/tof_3_pyzx.qc;; hwb6) echo benchmark/qc/optimized/hwb6_pyzx.qc;;
  hwb8) echo benchmark/qc/optimized/hwb8_pyzx.qc;;
  ham15-low) echo benchmark/qc/optimized/ham15-low_pyzx.qc;; ham15-med) echo benchmark/qc/optimized/ham15-med_pyzx.qc;;
  ham15-high) echo benchmark/qc/optimized/ham15-high_pyzx.qc;;
  hwb10) echo benchmark/qc/vandaele2024/hwb10.qc;; hwb11) echo benchmark/qc/vandaele2024/hwb11.qc;;
  *) echo "";; esac; }
sanitize(){ sed -e 's/\bT\*/tdg/g' -e 's/\bS\*/sdg/g' -e 's/\bZ\*/z/g' -e 's/^cnot /tof /' "$1" > "$2"; }

CIRC="${1:-}"; shift || true
MODES=("$@")
[ -n "$CIRC" ] && [ ${#MODES[@]} -gt 0 ] || { echo "usage: $(basename "$0") <circuit> <mode> [mode ...]   (mode = none|<size>|auto|full)" >&2; exit 1; }
rel="$(circuit_path "$CIRC")"; [ -n "$rel" ] && [ -f "$REPO/$rel" ] || { echo "unknown/missing circuit: $CIRC" >&2; exit 1; }
QC="$WORK/$CIRC.qc"; sanitize "$REPO/$rel" "$QC"

# mode token -> (label, gadget-command). One circuit, many gadgetization modes.
mode_label(){ case "$1" in none) echo "no-gadgetize";; full) echo "full-gadgetize";;
  auto) echo "adaptive(auto)*MACHINE-DEP";; *) echo "adaptive ($1)";; esac; }
mode_gad(){ case "$1" in none) echo "";; full) echo "tableau optimize gadgetize; ";;
  auto) echo "tableau optimize gadgetize --adaptive; ";; *) echo "tableau optimize gadgetize --memory-limit $1; ";; esac; }

# label, gadget-command (empty = no-gadgetize). Runs once, prints one table row.
emit(){ local label="$1" gad="$2" log="$WORK/run.log"
  local cmd="qcir read $QC; convert qcir tableau; tableau optimize tmerge; tableau optimize hopt; ${gad}tableau optimize phasepoly fasttodd; convert tableau qcir; qcir print; quit -f"
  # /usr/bin/time must WRAP timeout (not the reverse): on TLE, timeout kills qsyn and time still
  # prints the peak + "Command terminated by signal". The reverse loses the report on TLE.
  /usr/bin/time -v timeout -k5 "$TMO" "$QSYN" -q --no-version -c "$cmd" > "$log" 2>&1
  local tg w p
  tg="$(rg -o '[0-9]+ T-gate' "$log" | rg -o '^[0-9]+')"
  w="$(rg -o 'Elapsed \(wall clock\).*' "$log" | rg -o '[0-9:.]+$')"
  p="$(rg -o 'Maximum resident set size \(kbytes\): [0-9]+' "$log" | rg -o '[0-9]+$')"
  awk -v l="$label" -v tg="${tg:-}" -v w="${w:-}" -v p="${p:-0}" 'BEGIN{
    printf "%-20s %-8s %-10s %.0f MiB\n", l, (tg==""?"TLE":tg), (w==""?"TLE":w), p/1024 }'
}

echo "circuit=$CIRC  cap=${TMO}s  (sequential, isolated runs)"
printf '%-20s %-8s %-10s %s\n' mode T-count wall peak
for m in "${MODES[@]}"; do emit "$(mode_label "$m")" "$(mode_gad "$m")"; done
