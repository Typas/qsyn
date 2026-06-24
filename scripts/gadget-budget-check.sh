#!/usr/bin/env bash
# Reproducible safety check for `tableau optimize gadgetize --memory-limit`.
#
# For each (circuit, budget) it runs the gadgetize pipeline under a cgroup v2 `memory.max` cap equal to
# the budget (the production RSS enforcer; NOT `ulimit -v`, which caps virtual address space, not RSS),
# and reports the planner's predicted peak vs the measured peak RSS and the return code. A correct
# predict-only model means: never OOM-killed (rc != 137) and predicted >= actual.
#
# Fully portable: no absolute paths, circuits resolved from this repo's benchmark/qc, work files in a
# mktemp dir that is removed on exit. Needs: systemd user instance with cgroup v2 (Linux/WSL2),
# /usr/bin/time, rg. Override the binary with QSYN=/path/to/qsyn.
#
# Usage:
#   scripts/gadget-budget-check.sh                       # default sweep (all 7 circuits)
#   scripts/gadget-budget-check.sh hwb10:48G hwb11:48G   # explicit circuit:budget pairs
set -uo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
QSYN="${QSYN:-$REPO/build/qsyn}"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

[ -x "$QSYN" ] || { echo "qsyn binary not found at $QSYN (build it, or set QSYN=...)" >&2; exit 1; }
command -v systemd-run >/dev/null || { echo "systemd-run not found (need cgroup v2 user instance)" >&2; exit 1; }

# circuit name -> repo-relative .qc path (in-repo; no external benchmark checkout needed)
circuit_path() {
  case "$1" in
    tof3)       echo "benchmark/qc/optimized/tof_3_pyzx.qc" ;;
    hwb6)       echo "benchmark/qc/optimized/hwb6_pyzx.qc" ;;
    ham15-low)  echo "benchmark/qc/optimized/ham15-low_pyzx.qc" ;;
    ham15-med)  echo "benchmark/qc/optimized/ham15-med_pyzx.qc" ;;
    ham15-high) echo "benchmark/qc/optimized/ham15-high_pyzx.qc" ;;
    hwb10)      echo "benchmark/qc/vandaele2024/hwb10.qc" ;;
    hwb11)      echo "benchmark/qc/vandaele2024/hwb11.qc" ;;
    *)          echo "" ;;
  esac
}

# normalize gate names qsyn's qc reader rejects from the pyzx/vandaele dumps
sanitize() { sed -e 's/\bT\*/tdg/g' -e 's/\bS\*/sdg/g' -e 's/\bZ\*/z/g' -e 's/^cnot /tof /' "$1" > "$2"; }

run_one() {
  local name="$1" budget="$2"
  local rel; rel="$(circuit_path "$name")"
  if [ -z "$rel" ] || [ ! -f "$REPO/$rel" ]; then printf '%-12s %-5s  MISSING (%s)\n' "$name" "$budget" "$rel"; return; fi
  local qc="$WORK/$name.qc"; sanitize "$REPO/$rel" "$qc"
  local log="$WORK/$name.$budget.log"
  systemd-run --user --scope -q -p MemoryMax="$budget" -p MemorySwapMax=0 \
    /usr/bin/time -v "$QSYN" -q --no-version -c \
    "logger debug; qcir read $qc; convert qcir tableau; tableau optimize tmerge; tableau optimize hopt; tableau optimize gadgetize --memory-limit $budget; quit -f" \
    > "$log" 2>&1
  local rc=$?
  local pred act
  pred="$(rg -o 'predicted peak [0-9]+ MiB' "$log" | rg -o '[0-9]+' | head -1)"
  act="$(rg -o 'Maximum resident set size \(kbytes\): [0-9]+' "$log" | rg -o '[0-9]+$')"
  # cap in MiB for a utilization figure (handles K/M/G/T suffixes)
  local capm; capm="$(awk -v s="$budget" 'BEGIN{u=substr(s,length(s),1); v=substr(s,1,length(s)-1)+0;
    if(u=="G")v*=1024; else if(u=="T")v*=1048576; else if(u=="K")v/=1024; print v}')"
  awk -v n="$name" -v b="$budget" -v rc="$rc" -v p="${pred:-NA}" -v a="${act:-NA}" -v capm="$capm" 'BEGIN{
    am = (a=="NA"?"NA":a/1024);
    safe = (rc!=137) ? "no-OOM" : "OOM!!";                 # the predict-only safety guarantee
    util = (a=="NA"||capm==0) ? "NA" : sprintf("%.0f%%", 100*(a/1024)/capm);   # actual/budget
    bound = (p=="NA"||a=="NA") ? "?" : (p+0 >= a/1024-1 ? "predict>=actual" : "UNDER!!"); # 1 MiB floor tolerance
    printf "%-12s cap=%-5s rc=%-3s predicted=%5s MiB  actual=%5.0f MiB  util=%-5s %s  %s\n", n,b,rc,p,am,util,safe,bound }'
}

pairs=("$@")
if [ ${#pairs[@]} -eq 0 ]; then
  pairs=(tof3:256M hwb6:256M ham15-low:256M ham15-med:256M ham15-high:256M hwb10:1G hwb11:1G)
fi
printf '%-12s %-9s %-6s %-22s %-18s %s\n' circuit budget rc predicted actual verdict
for pair in "${pairs[@]}"; do run_one "${pair%%:*}" "${pair##*:}"; done
