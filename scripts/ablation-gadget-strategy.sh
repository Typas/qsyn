#!/usr/bin/env bash
# Gadgetize min-T merge-strategy ablation (docs/plan-gadgetize-mint.md).
#
# Holds the memory budget as a hard cap and makes T-count the objective: per circuit at its fixed
# gadgetize budget, sweep the four boundary-merge strategies and record final T-count (primary), peak
# RSS, and wall time. The lowest T-count wins; the winner is later folded into the planner and this
# scaffolding deleted.
#
# qsyn-only: this does NOT use the 2x2 Vandaele harness (no reference binary, no wasted runs). Every
# cell (circuit x strategy) is independent, small, and single-threaded, so all 8 run in PARALLEL; each
# is wrapped in `ulimit -v $CAP_KB` + `timeout` so a runaway region can never exhaust host memory.
# Phasepoly runs sequentially (no -j) for reproducibility.
#
# Reproducible: fixed per-circuit `--memory-limit` (never --adaptive, which is machine-dependent).
# Re-run with the same args to reproduce. Records go to $OUTDIR.
#
# Usage:
#   scripts/ablation-gadget-strategy.sh
# Env (with defaults):
#   QSYN=build/qsyn  CAP_KB=8388608 (8G ulimit -v safety net)  TMO=1800 (per-cell timeout, seconds)
#   OUTDIR=bench2x2/ablation
set -u

QSYN="${QSYN:-build/qsyn}"
CAP_KB="${CAP_KB:-8388608}"
TMO="${TMO:-1800}"
OUTDIR="${OUTDIR:-bench2x2/ablation}"

# Matrix: "label|sanitized-input|gadgetize-budget". Budgets sit in the multi-region regime, where the
# circuit splits into many bounded regions so the merge strategy matters (above it the circuit
# full-merges to one region and all strategies coincide).
MATRIX=(
  "hwb8|bench2x2/sani/hwb8__qsyn_raw.qc|256M"
  "ham15-high|bench2x2/sani/ham15-high__qsyn_raw.qc|64M"
)
STRATEGIES=(z-overlap max-merge max-terms balanced-fill)

mkdir -p "$OUTDIR/logs" "$OUTDIR/rows"
CSV="$OUTDIR/ablation.csv"
TABLE="$OUTDIR/table.md"

# Hard memory cap on every qsyn call. MALLOC_ARENA_MAX=1 + single-thread BLAS keep the virtual address
# space ~= RSS, so `ulimit -v` (which caps AS, not RSS) is a meaningful memory limit.
capped() { (
  ulimit -v "$CAP_KB"
  MALLOC_ARENA_MAX=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 "$@"
); }
peak_of() { rg -o 'Maximum resident set size .kbytes.: [0-9]+' "$1" | rg -o '[0-9]+$' | tail -1; }
status_of() { # $1=log $2=rc
  if [ "$2" -eq 124 ]; then
    echo TLE
  elif rg -qi 'bad_alloc|out of memory|length_error' "$1"; then
    echo MLE
  elif [ "$2" -ne 0 ]; then
    echo "ERR($2)"
  else echo ok; fi
}

# One cell: gadgetize with $strat under $budget, then phasepoly fasttodd, and record the row.
run_cell() {
  local label="$1" sani="$2" budget="$3" strat="$4"
  local log="$OUTDIR/logs/${label}__${strat}.log"
  capped /usr/bin/time -v timeout -k5 "$TMO" "$QSYN" -q --no-version -c \
    "qcir read $sani; convert qcir tableau; tableau optimize tmerge; tableau optimize hopt; tableau optimize gadgetize --strategy $strat -m $budget; tableau optimize phasepoly fasttodd; convert tableau qcir; qcir print --stat; quit -f" >"$log" 2>&1
  local rc=$? t peak wall st
  t=$(rg -o 'T-family\s*:\s*[0-9]+' "$log" | rg -o '[0-9]+$' | tail -1)
  peak=$(peak_of "$log")
  st=$(status_of "$log" "$rc")
  wall=$(rg -o 'wall clock.*: [0-9:.]+' "$log" | rg -o '[0-9:.]+$' | tail -1)
  echo "$label,$budget,$strat,${t:-},${peak:-},${wall:-},$rc,$st" >"$OUTDIR/rows/${label}__${strat}.row"
}

echo "[ablation] launching $((${#MATRIX[@]} * ${#STRATEGIES[@]})) cells in parallel (cap ${CAP_KB}KB, timeout ${TMO}s)"
for entry in "${MATRIX[@]}"; do
  IFS='|' read -r label sani budget <<<"$entry"
  if [ ! -f "$sani" ]; then
    echo "[ablation] MISSING INPUT: $sani" >&2
    continue
  fi
  for strat in "${STRATEGIES[@]}"; do
    run_cell "$label" "$sani" "$budget" "$strat" &
  done
done
wait
echo "[ablation] all cells done"

# Assemble CSV + a readable table, ordered by the matrix then strategy.
echo "circuit,budget,strategy,t_count,peak_kb,wall_s,rc,status" >"$CSV"
{
  echo "# Gadgetize merge-strategy ablation (T-count primary; lowest wins)"
  echo
  printf '| %-10s | %-7s | %-13s | %8s | %10s | %8s | %s |\n' circuit budget strategy t_count peak_kb wall status
  printf '|%s|%s|%s|%s|%s|%s|%s|\n' "------------" "---------" "---------------" "----------" "------------" "----------" "--------"
} >"$TABLE"
for entry in "${MATRIX[@]}"; do
  IFS='|' read -r label _ _ <<<"$entry"
  for strat in "${STRATEGIES[@]}"; do
    row="$OUTDIR/rows/${label}__${strat}.row"
    [ -f "$row" ] || continue
    cat "$row" >>"$CSV"
    IFS=',' read -r c b s t peak wall rc st <<<"$(cat "$row")"
    printf '| %-10s | %-7s | %-13s | %8s | %10s | %8s | %s |\n' "$c" "$b" "$s" "${t:-?}" "${peak:-?}" "${wall:-?}" "$st" >>"$TABLE"
  done
done
echo "[ablation] CSV   -> $CSV"
echo "[ablation] table -> $TABLE"
cat "$TABLE"
