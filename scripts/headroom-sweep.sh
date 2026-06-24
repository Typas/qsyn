#!/usr/bin/env bash
# L3 validation -- re-validate QSYN_GADGET_HEADROOM under threads, per allocator.
#
# The headroom percent (gadget_alloc_headroom_percent, adaptive_gadget.cpp:382; default 150, measured on
# glibc single-thread) is the one allocator-dependent knob in the predict-only memory model: real RSS =
# structural bytes x allocator overhead. Under the threaded phasepoly executor, k concurrent regions
# compound allocator fragmentation, and a different allocator returns freed gadgetize arenas differently,
# so the single-thread glibc 150 may no longer satisfy predicted >= actual. This sweep measures, per
# allocator and at the worst thread count, the headroom floor that keeps the run within the cgroup cap
# (no OOM) AND restores predicted >= actual.
#
# It drives scripts/cgrun-phasepoly-parallel.sh (the tested capped-run core) once per (allocator,
# headroom) cell at a fixed thread count, then aggregates peak vs predicted vs cap.
#
# Reproducible: fixed --memory-limit == cgroup cap, never --adaptive. Results committed under
# bench2x2/l3-headroom/. Needs the allocator .so libs at repo root (libjemalloc.so, libmimalloc.so.2/.3),
# systemd cgroup v2, rg. Override the binary with QSYN=/path/to/qsyn.
#
# Usage: scripts/headroom-sweep.sh [CIRCUIT] [BUDGET] [THREADS] [HEADROOMS...]
#   defaults: bench2x2/sani/hwb10__qsyn_raw.qc  256M  8  (headrooms 150 175 200)
set -uo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
QSYN="${QSYN:-$REPO/build/qsyn}"
CIRCUIT="${1:-$REPO/bench2x2/sani/hwb10__qsyn_raw.qc}"
BUDGET="${2:-256M}"
THREADS="${3:-8}" # worst thread count from the l3-cgroup sweep (fragmentation peaks here)
shift $(($# < 3 ? $# : 3)) || true
HEADROOMS=("$@")
[ ${#HEADROOMS[@]} -eq 0 ] && HEADROOMS=(150 175 200)

CGRUN="$REPO/scripts/cgrun-phasepoly-parallel.sh"
CGCSV="$REPO/bench2x2/l3-cgroup/results.csv" # cgrun writes its per-call row here
OUTDIR="$REPO/bench2x2/l3-headroom"
mkdir -p "$OUTDIR"
csv="$OUTDIR/results.csv"
tbl="$OUTDIR/table.md"

# allocator name -> LD_PRELOAD lib ("" = system glibc)
allocs=(glibc jemalloc mimalloc2 mimalloc3)
libfor() { case "$1" in
  glibc) echo "" ;;
  jemalloc) echo "$REPO/libjemalloc.so" ;;
  mimalloc2) echo "$REPO/libmimalloc.so.2" ;;
  mimalloc3) echo "$REPO/libmimalloc.so.3" ;;
  esac }

echo "allocator,budget,threads,headroom,rc,peak_kib,predicted_mib,cap_kib,predict_ge_actual,no_oom,verdict" >"$csv"
printf '%-10s %-4s %-3s %-3s %9s %6s %8s  %-13s %s\n' alloc thr hr rc peak_kib pred_M margin verdict
for a in "${allocs[@]}"; do
  lib="$(libfor "$a")"
  [ -z "$lib" ] || [ -f "$lib" ] || {
    echo "missing allocator lib for $a: $lib (skip)" >&2
    continue
  }
  for hr in "${HEADROOMS[@]}"; do
    PRELOAD="$lib" QSYN_HEADROOM="$hr" "$CGRUN" "$CIRCUIT" "$BUDGET" "$THREADS" >/dev/null 2>&1 || true
    # cgrun wrote exactly one data row (single thread count) -> read it.
    row="$(tail -n 1 "$CGCSV")"
    IFS=, read -r _c _b _al _hr _th rc peak pred capk margin _tc _st <<<"$row"
    pge="$(awk -v p="${pred:-0}" -v k="${peak:-0}" 'BEGIN{print (p*1024 >= k)?"yes":"no"}')"
    noo="$([ "$rc" = "137" ] && echo no || echo yes)"
    verdict="$([ "$pge" = yes ] && [ "$noo" = yes ] && echo SAFE || echo UNSAFE)"
    printf '%-10s %-4s %-3s %-3s %9s %6s %7s%%  %-13s %s\n' "$a" "$_th" "$hr" "$rc" "${peak:-NA}" "${pred:-NA}" "${margin:-NA}" "predge=$pge" "$verdict"
    echo "$a,$BUDGET,$_th,$hr,$rc,${peak:-NA},${pred:-NA},$capk,$pge,$noo,$verdict" >>"$csv"
  done
done

{
  echo "# L3 headroom re-validation under threads -- $(basename "$CIRCUIT") @ $BUDGET, threads=$THREADS"
  echo
  echo "Cap == budget. SAFE = predicted >= actual peak AND no OOM (rc!=137). The headroom floor per"
  echo "allocator is the smallest headroom whose row is SAFE."
  echo
  echo "| allocator | headroom | rc | peak_kib | predicted_mib | predict>=actual | no_oom | verdict |"
  echo "|-----------|----------|----|----------|---------------|-----------------|--------|---------|"
  tail -n +2 "$csv" | awk -F, '{printf "| %s | %s | %s | %s | %s | %s | %s | %s |\n",$1,$4,$5,$6,$7,$9,$10,$11}'
} >"$tbl"
echo "-> $tbl ; $csv"
