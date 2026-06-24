#!/usr/bin/env bash
# Measure the region-parallel phasepoly executor on one circuit: T-count, wall
# time, peak RSS. Reproducible -- the budget is a fixed --memory-limit, never
# --adaptive, so the result is identical across machines and thread counts.
#
# Usage: scripts/bench-phasepoly-parallel.sh CIRCUIT [THREADS] [GADGET] [BUDGET]
#   CIRCUIT  input circuit (.qc/.qasm)
#   THREADS  worker count (default 1 = sequential path)
#   GADGET   none | full | budget   (default none)
#   BUDGET   fixed memory budget (e.g. 8G); REQUIRED iff THREADS>1 or GADGET=budget,
#            and must be omitted otherwise (full gadgetize and the sequential path take no budget)
#   env QSYN  path to qsyn executable (default ./build/qsyn)
#
# GADGET=none is the only unitary path, so it is the only one whose output can be
# checked for equivalence against the input (`qcir equiv`); the gadgetized paths
# are measurement-based and report equiv=n/a.
set -euo pipefail

qsyn=${QSYN:-./build/qsyn}
timeout=${TIMEOUT:-1800} # per-run wall limit in seconds (default 30 min); env TIMEOUT overrides
circuit=$1
threads=${2:-1}
gadget=${3:-none}
budget=${4:-}

needs_budget=no
[[ $threads -gt 1 || $gadget == budget ]] && needs_budget=yes
if [[ $needs_budget == yes && -z $budget ]]; then
  echo "BUDGET required for THREADS>1 or GADGET=budget" >&2
  exit 1
fi
if [[ $needs_budget == no && -n $budget ]]; then
  echo "BUDGET must be omitted for sequential, no-budget runs" >&2
  exit 1
fi

dof=$(mktemp)
tf=$(mktemp)
trap 'rm -f "$dof" "$tf"' EXIT
{
  echo "qcir read $circuit"
  echo "convert qcir tableau"
  echo "tableau optimize tmerge"
  echo "tableau optimize hopt"
  case $gadget in
  none) ;;
  full) echo "tableau optimize gadgetize" ;;
  budget) echo "tableau optimize gadgetize --memory-limit $budget" ;;
  *)
    echo "GADGET must be none|full|budget" >&2
    exit 1
    ;;
  esac
  if [[ $threads -gt 1 ]]; then
    echo "tableau optimize phasepoly fasttodd --threads $threads --memory-limit $budget"
  else
    echo "tableau optimize phasepoly fasttodd"
  fi
  echo "convert tableau qcir"
  echo "qcir print --stat"
  echo "quit -f"
} >"$dof"

# Timed run: pipeline only, NO equiv (equiv is a heavy verification, not the workload being measured).
rc=0
out=$(timeout -k 30 "$timeout" /usr/bin/time -v "$qsyn" -f "$dof" 2>"$tf") || rc=$?
if [[ $rc -eq 124 ]]; then
  echo "circuit=$(basename "$circuit") gadget=$gadget threads=$threads budget=${budget:-none} status=TIMEOUT(${timeout}s)"
  exit 0
fi
tcount=$(printf '%s\n' "$out" | rg -o 'T-family +: [0-9]+' | rg -o '[0-9]+$')
wall=$(rg -o 'wall clock.*\): .*' "$tf" | rg -o '[0-9:.]+$')
rss_mb=$(($(rg -o 'Maximum resident set size \(kbytes\): [0-9]+' "$tf" | rg -o '[0-9]+$') / 1024))

# Separate untimed verification: only the unitary (no-gadget) path is equivalence-checkable.
# input qcir is id 0, the resynthesized one is id 1.
equiv=n/a
if [[ $gadget == none ]]; then
  ev=$(printf '%s\nqcir equiv 0 1\nquit -f\n' "$(sed '/^quit -f$/d' "$dof")" | timeout -k 30 "$timeout" "$qsyn" 2>&1) || true
  if printf '%s\n' "$ev" | rg -qi 'are equivalent'; then equiv=yes; else equiv=NO; fi
fi

echo "circuit=$(basename "$circuit") gadget=$gadget threads=$threads budget=${budget:-none} tcount=$tcount wall=$wall rss=${rss_mb}MB equiv=$equiv"
