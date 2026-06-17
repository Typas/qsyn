#!/usr/bin/env bash
# Run the qsyn FastTODD pipeline (tmerge -> hopt -> phasepoly fasttodd, the
# no-ancilla Vandaele-comparable path) on one benchmark circuit and append a
# CSV row.
#
# Usage: todd-benchmark.sh <circuit-name> [timeout-seconds]
#
# Reads circuit metadata from benchmark/paper-todd.csv (column `file`),
# uses the sanitized copy in qsyn-bench/circuits/<circuit>.qc, and appends
#   circuit,n_qubits,pre_T,post_T,wall_seconds,status,equiv
# to qsyn-bench/results.csv (status: ok | DNF | error).
#
# equiv records whether the optimized circuit was proven equivalent to the
# original via `qcir equiv`:
#   pass       - proven equivalent (tableau, or tensor for <=7 qubits)
#   FAIL       - tensor contraction proved them NOT equivalent (real bug!)
#   unverified - could not be proven either way (qsyn warns this may be a
#                false negative; happens for >7-qubit circuits the tableau
#                method can't reduce and tensor contraction can't handle)
#   (empty)    - circuit did not reach the equivalence check (DNF/error)

set -u

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
QSYN="$REPO_ROOT/qsyn"
RUN_DIR="$REPO_ROOT/qsyn-bench"
CIRCUITS_DIR="$RUN_DIR/circuits"
RESULTS="$RUN_DIR/results.csv"
LOG_DIR="$RUN_DIR/logs"

CIRCUIT="${1:?usage: todd-benchmark.sh <circuit-name> [timeout-seconds]}"
TIMEOUT="${2:-1800}"

QC_FILE="$CIRCUITS_DIR/$CIRCUIT.qc"
mkdir -p "$LOG_DIR"
[ -f "$RESULTS" ] || echo "circuit,n_qubits,pre_T,post_T,wall_seconds,status,equiv" >"$RESULTS"

# Drop any previous row for this circuit so reruns overwrite.
if [ -f "$RESULTS" ]; then
	awk -F, -v c="$CIRCUIT" 'NR==1 || $1 != c' "$RESULTS" >"$RESULTS.tmp" && mv "$RESULTS.tmp" "$RESULTS"
fi

if [ ! -f "$QC_FILE" ]; then
	echo "$CIRCUIT,,,,0,error," >>"$RESULTS"
	echo "[$CIRCUIT] missing sanitized circuit file $QC_FILE (run \`just prepare\` first)" >&2
	exit 1
fi

parse_stat() { # stdin: qsyn output; prints "<qubits> <first T-family> <last T-family>"
	awk '
        /^QCir \(/    { gsub(/[(,]/, "", $2); if (!q) q = $2 }
        /^T-family/   { t = $3; if (pre == "") pre = t; post = t }
        END           { print q+0, pre+0, post+0 }
    '
}

# Fast pre-pass: parse + input stats only (no optimization).
PRE_OUT=$("$QSYN" -q --no-version -c "qcir read $QC_FILE; qcir print --stat; quit -f" 2>&1)
if echo "$PRE_OUT" | rg -q "error|Error"; then
	echo "$CIRCUIT,,,,0,error," >>"$RESULTS"
	echo "[$CIRCUIT] qsyn failed to read circuit:" >&2
	echo "$PRE_OUT" | head -5 >&2
	exit 1
fi
read -r N_QUBITS PRE_T _ <<<"$(echo "$PRE_OUT" | parse_stat)"

# Pipeline keeps the original circuit as QCir 0 and the optimized result as
# QCir 1 (created by `convert tableau qcir`), then proves equivalence between
# them with `qcir equiv 0 1`.
LOG="$LOG_DIR/$CIRCUIT.log"
START=$(date +%s)
# qsyn installs its own SIGTERM handler and does not abort mid-computation,
# so a plain `timeout` would hang forever once it sends SIGTERM. `-k 30`
# escalates to SIGKILL (uncatchable) 30s later, guaranteeing the cap holds.
timeout -k 30 "$TIMEOUT" "$QSYN" -q --no-version -c \
	"qcir read $QC_FILE; qcir print --stat; convert qcir tableau; tableau optimize tmerge; tableau optimize hopt; tableau optimize phasepoly fasttodd; convert tableau qcir; qcir print --stat; qcir equiv 0 1; quit -f" \
	>"$LOG" 2>&1
RC=$?
WALL=$(($(date +%s) - START))

# 124 = SIGTERM hit the cap; 137 = SIGKILL escalation (-k) had to force it.
if [ "$RC" -eq 124 ] || [ "$RC" -eq 137 ]; then
	echo "$CIRCUIT,$N_QUBITS,$PRE_T,,$WALL,DNF," >>"$RESULTS"
	echo "[$CIRCUIT] DNF after ${TIMEOUT}s (n=$N_QUBITS, pre T=$PRE_T)"
	exit 0
elif [ "$RC" -ne 0 ] || rg -q "error|Error" "$LOG"; then
	echo "$CIRCUIT,$N_QUBITS,$PRE_T,,$WALL,error," >>"$RESULTS"
	echo "[$CIRCUIT] pipeline error (rc=$RC), see $LOG" >&2
	exit 1
fi

# Classify the `qcir equiv` result. Note "are equivalent" is a substring of
# "are not equivalent", so test the negative form first.
if rg -q "are not equivalent" "$LOG"; then
	if rg -q "false negative|too large to check equivalence" "$LOG"; then
		EQUIV=unverified
	else
		EQUIV=FAIL
	fi
elif rg -q "are equivalent" "$LOG"; then
	EQUIV=pass
else
	EQUIV=unverified
fi

read -r _ _ POST_T <<<"$(parse_stat <"$LOG")"
echo "$CIRCUIT,$N_QUBITS,$PRE_T,$POST_T,$WALL,ok,$EQUIV" >>"$RESULTS"
echo "[$CIRCUIT] n=$N_QUBITS pre T=$PRE_T -> post T=$POST_T (${WALL}s) equiv=$EQUIV"
if [ "$EQUIV" = "FAIL" ]; then
	echo "[$CIRCUIT] WARNING: optimized circuit is NOT equivalent to the original!" >&2
fi
