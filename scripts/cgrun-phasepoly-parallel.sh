#!/usr/bin/env bash
# L3 validation -- region-parallel phasepoly executor under a cgroup v2 memory.max cap.
#
# Checks the two guarantees the L3 plan (docs/plan-gadgetization-l3.md) mandates for the threaded
# phasepoly executor:
#   1. NO-OOM: with the cgroup memory.max cap set equal to the gadgetize/phasepoly budget, no thread
#      count OOM-kills (rc != 137). Predict-only admission means predicted peak >= actual RSS.
#   2. THREAD-COUNT-INVARIANT T-count: every thread count yields the same final T-count as --threads 1.
#
# Only the BUDGETED optimization (gadgetize + phasepoly) runs under the cap -- the tableau->qcir export
# is an unbudgeted transient (it rebuilds the full output circuit) and would dominate the peak, so it is
# excluded from the capped workload and the T-count is read in a separate uncapped run.
#
# Reproducible: fixed --memory-limit, never --adaptive. Portable: no absolute paths, work files in the
# scratch dir, results committed under bench2x2/l3-cgroup/. Needs a systemd user instance with cgroup v2
# (Linux/WSL2), /usr/bin/time, rg. Override the binary with QSYN=/path/to/qsyn; pin headroom with
# QSYN_HEADROOM=<pct>; pin an allocator with PRELOAD=/path/to/lib.so.
#
# Usage: scripts/cgrun-phasepoly-parallel.sh CIRCUIT BUDGET [THREADS...]
#   CIRCUIT  input .qc (e.g. bench2x2/sani/hwb10__qsyn_raw.qc)
#   BUDGET   fixed memory budget == cgroup cap (e.g. 256M)
#   THREADS  worker counts to sweep (default: 1 2 4 8)
set -uo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
QSYN="${QSYN:-$REPO/build/qsyn}"
CIRCUIT="$1"
BUDGET="$2"
shift 2
THREADS=("$@")
[ ${#THREADS[@]} -eq 0 ] && THREADS=(1 2 4 8)
OUTDIR="$REPO/bench2x2/l3-cgroup"
mkdir -p "$OUTDIR"

[ -x "$QSYN" ] || {
  echo "qsyn not found at $QSYN (build it, or set QSYN=...)" >&2
  exit 1
}
[ -f "$CIRCUIT" ] || {
  echo "circuit not found: $CIRCUIT" >&2
  exit 1
}
command -v systemd-run >/dev/null || {
  echo "systemd-run not found (need cgroup v2 user instance)" >&2
  exit 1
}

# BUDGET (K/M/G/T suffix) -> KiB, for the OOM-margin column.
cap_kib() { awk -v s="$1" 'BEGIN{u=substr(s,length(s),1); v=substr(s,1,length(s)-1)+0;
  if(u=="G")v*=1048576; else if(u=="M")v*=1024; else if(u=="T")v*=1073741824; else if(u=="K")v*=1; print int(v)}'; }
CAPKIB="$(cap_kib "$BUDGET")"

# shared pipeline up to (and including) gadgetize; phasepoly clause is appended per thread count.
PRE="logger debug; qcir read $CIRCUIT; convert qcir tableau; tableau optimize tmerge; tableau optimize hopt; tableau optimize gadgetize --memory-limit $BUDGET"
# env prefix carrying the optional allocator / headroom pins into qsyn's environment.
envwrap() {
  local -a e=(env)
  [ -n "${PRELOAD:-}" ] && e+=("LD_PRELOAD=$PRELOAD")
  [ -n "${QSYN_HEADROOM:-}" ] && e+=("QSYN_GADGET_HEADROOM=$QSYN_HEADROOM")
  printf '%s\n' "${e[@]}"
}
mapfile -t ENVV < <(envwrap)

pp_clause() { # thread count -> phasepoly command (>1 requires a budget)
  if [ "$1" -le 1 ]; then
    echo "tableau optimize phasepoly fasttodd"
  else echo "tableau optimize phasepoly fasttodd --threads $1 --memory-limit $BUDGET"; fi
}

csv="$OUTDIR/results.csv"
tbl="$OUTDIR/table.md"
echo "circuit,budget,allocator,headroom,threads,rc,peak_kib,predicted_mib,cap_kib,margin_pct,tcount,status" >"$csv"
alloc="${PRELOAD:+$(basename "$PRELOAD")}"
alloc="${alloc:-glibc}"
hr="${QSYN_HEADROOM:-default}"
base_tc=""
worst="ok"
printf '%-28s %-6s %-10s %-4s %-3s %9s %6s %8s %s\n' circuit budget alloc thr rc peak_kib pred_M margin status
for n in "${THREADS[@]}"; do
  pp="$(pp_clause "$n")"
  # (1) capped no-OOM run: gadgetize + phasepoly only, no export.
  out="$(systemd-run --user --scope -q -p MemoryMax="$BUDGET" -p MemorySwapMax=0 \
    "${ENVV[@]}" /usr/bin/time -v "$QSYN" -q --no-version -c "$PRE; $pp; quit -f" 2>&1)"
  rc=$?
  peak="$(printf '%s\n' "$out" | rg -o 'Maximum resident set size \(kbytes\): [0-9]+' | rg -o '[0-9]+$')"
  pred="$(printf '%s\n' "$out" | rg -o 'predicted peak [0-9]+ MiB' | rg -o '[0-9]+' | head -1)"
  # (2) uncapped T-count run: same optimization + export, to check thread-invariance.
  tco="$("${ENVV[@]}" "$QSYN" -q --no-version -c "$PRE; $pp; convert tableau qcir; qcir print --stat; quit -f" 2>&1)"
  tc="$(printf '%s\n' "$tco" | rg -o 'T-family +: [0-9]+' | rg -o '[0-9]+$')"
  [ -z "$base_tc" ] && base_tc="$tc"
  margin="$(awk -v p="${peak:-0}" -v c="$CAPKIB" 'BEGIN{ if(c>0) printf "%.1f", 100*(c-p)/c; else print "NA"}')"
  status="ok"
  [ "$rc" = "137" ] && {
    status="OOM"
    worst="OOM"
  }
  [ -n "$tc" ] && [ "$tc" != "$base_tc" ] && {
    status="T-DIVERGE($tc!=$base_tc)"
    worst="T-DIVERGE"
  }
  printf '%-28s %-6s %-10s %-4s %-3s %9s %6s %8s %s\n' "$(basename "$CIRCUIT")" "$BUDGET" "$alloc" "$n" "$rc" "${peak:-NA}" "${pred:-NA}" "$margin%" "$status"
  echo "$(basename "$CIRCUIT"),$BUDGET,$alloc,$hr,$n,$rc,${peak:-NA},${pred:-NA},$CAPKIB,$margin,${tc:-NA},$status" >>"$csv"
done

{
  echo "# L3 cgroup no-OOM + thread-invariance -- $(basename "$CIRCUIT") @ $BUDGET (allocator $alloc, headroom $hr)"
  echo
  echo "| threads | rc | peak_kib | predicted_mib | cap_kib | margin | tcount | status |"
  echo "|---------|----|----------|---------------|---------|--------|--------|--------|"
  tail -n +2 "$csv" | awk -F, '{printf "| %s | %s | %s | %s | %s | %s%% | %s | %s |\n",$5,$6,$7,$8,$9,$10,$11,$12}'
  echo
  echo "Cap == budget. status ok = no OOM (rc!=137) and T-count == --threads 1 baseline ($base_tc)."
} >"$tbl"
echo "-> $tbl ; $csv"
[ "$worst" = "ok" ] || {
  echo "FAIL: $worst" >&2
  exit 1
}
echo "PASS: no OOM, T-count thread-invariant across ${THREADS[*]}"
