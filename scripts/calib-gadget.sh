#!/usr/bin/env bash
# Adaptive-gadget L2 memory-model calibration sweep.
#
# For each circuit: (1) gadgetize-only -> merged region size (n,m) + GADGETIZE-STAGE peak RSS
#                       (isolates the widen/collapse cost -- the step that OOM'd at 56GB);
#                   (2) gadgetize + phasepoly fasttodd -> full-pipeline peak RSS, wall, status.
# EVERY qsyn invocation is wrapped in `ulimit -v $CAP_KB` + `timeout`, so a runaway
# region can never exhaust host memory (the 2026-06-23 hwb10 OOM/WSL-restart was an
# *uncapped* probe -- never again). Output is a CSV row per circuit.
#
# Reproducible: re-run with the same args to reproduce. Records go to $OUTDIR.
#
# Usage:
#   scripts/calib-gadget.sh [circuit ...]
# Env (with defaults):
#   QSYN=build/qsyn  VDIR=/home/liao/quantum-circuit-optimization
#   CAP_KB=8388608 (8G ulimit -v)  T_GAD=30  T_FULL=60 (timeouts, seconds)
#   OUTDIR=bench2x2/calib
# Do NOT pass hwb8/hwb10/hwb11 here: full-gadget on those TLEs/MLEs by design.
set -u

QSYN="${QSYN:-build/qsyn}"
VDIR="${VDIR:-/home/liao/quantum-circuit-optimization}"
CAP_KB="${CAP_KB:-8388608}"
T_GAD="${T_GAD:-30}"
T_FULL="${T_FULL:-60}"
OUTDIR="${OUTDIR:-bench2x2/calib}"

DEFAULT_CIRCUITS=(tof_4 tof_5 barenco_tof_4 barenco_tof_5 hwb6 ham15-low ham15-med
                  qcla_mod_7 mod_red_21 rc_adder_6 csla_mux_3 grover_5)
CIRCUITS=("$@"); [ ${#CIRCUITS[@]} -eq 0 ] && CIRCUITS=("${DEFAULT_CIRCUITS[@]}")

mkdir -p "$OUTDIR/sani" "$OUTDIR/logs"
CSV="$OUTDIR/calib.csv"

sanitize(){ sed -e 's/\bT\*/tdg/g' -e 's/\bS\*/sdg/g' -e 's/\bZ\*/z/g' -e 's/^cnot /tof /' "$1" > "$2"; }
# Hard memory cap on every qsyn call. MALLOC_ARENA_MAX=1 + single-thread BLAS keep the virtual
# address space ~= RSS, so `ulimit -v` (which caps AS, not RSS) is a meaningful memory limit --
# without it qsyn's arena VM alone exceeds ~1 GB even at tiny RSS (see worklog 2026-06-23).
capped(){ ( ulimit -v "$CAP_KB"; MALLOC_ARENA_MAX=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 "$@" ); }

# peak extractor for a /usr/bin/time -v log
peak_of(){ rg -o 'Maximum resident set size .kbytes.: [0-9]+' "$1" | rg -o '[0-9]+$' | tail -1; }
status_of(){ # $1=log $2=rc
  if   [ "$2" -eq 124 ]; then echo TLE
  elif rg -qi 'bad_alloc|out of memory|length_error' "$1"; then echo MLE
  elif [ "$2" -ne 0 ]; then echo "ERR($2)"
  else echo ok; fi; }

echo "circuit,n_gad,m_gad,gad_peak_kb,gad_st,full_peak_kb,full_wall_s,full_rc,full_st" > "$CSV"
printf '%-14s %6s %6s %11s %6s %11s %8s %s\n' circuit n_gad m_gad gad_peak gad_st full_peak full_wall full_st
for c in "${CIRCUITS[@]}"; do
  in="$VDIR/circuits/inputs/$c.qc"
  [ -f "$in" ] || { printf '%-14s %s\n' "$c" "MISSING-INPUT"; echo "$c,,,,,,,," >> "$CSV"; continue; }
  sani="$OUTDIR/sani/$c.qc"; sanitize "$in" "$sani"

  # STAGE 1: gadgetize-only, timed -> isolates widen/collapse peak
  golog="$OUTDIR/logs/$c.gadonly.log"
  capped /usr/bin/time -v timeout -k5 "$T_GAD" "$QSYN" -q --no-version -c \
    "qcir read $sani; convert qcir tableau; tableau optimize tmerge; tableau optimize hopt; tableau optimize gadgetize; convert tableau qcir; qcir print --stat; quit -f" > "$golog" 2>&1
  grc=$?
  n=$(rg -o '[0-9]+ qubits' "$golog" | rg -o '[0-9]+' | tail -1)
  m=$(rg -o 'T-family\s*:\s*[0-9]+' "$golog" | rg -o '[0-9]+$' | tail -1)
  gpeak=$(peak_of "$golog"); gst=$(status_of "$golog" "$grc")

  # STAGE 2: gadgetize + phasepoly, timed -> full-pipeline peak
  log="$OUTDIR/logs/$c.full.log"
  capped /usr/bin/time -v timeout -k5 "$T_FULL" "$QSYN" -q --no-version -c \
    "qcir read $sani; convert qcir tableau; tableau optimize tmerge; tableau optimize hopt; tableau optimize gadgetize; tableau optimize phasepoly fasttodd; convert tableau qcir; qcir print --stat; quit -f" > "$log" 2>&1
  rc=$?
  peak=$(peak_of "$log"); st=$(status_of "$log" "$rc")
  wall=$(rg -o 'wall clock.*: [0-9:.]+' "$log" | rg -o '[0-9:.]+$' | tail -1)

  printf '%-14s %6s %6s %11s %6s %11s %8s %s\n' "$c" "${n:-?}" "${m:-?}" "${gpeak:-?}" "$gst" "${peak:-?}" "${wall:-?}" "$st"
  echo "$c,${n:-},${m:-},${gpeak:-},$gst,${peak:-},${wall:-},$rc,$st" >> "$CSV"
done
echo "[calib] CSV -> $CSV"
