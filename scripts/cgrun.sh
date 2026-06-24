#!/usr/bin/env bash
# Run a qsyn dofile under a cgroup v2 memory.max cap (RSS-based, arena/thread independent) and
# report peak RSS + return code. This is the production/calibration enforcer for the
# memory-budgeted gadgetize feature (cgroup memory.max, NOT ulimit -v which caps AS not RSS).
#
# Usage: scripts/cgrun.sh <cap, e.g. 12G|512M> <qsyn -c command string>
# Env:   QSYN (default build/qsyn)
# Prints: "rc=<code> peak=<KiB>" ; rc=137 means OOM-killed by the cgroup.
set -u
QSYN="${QSYN:-build/qsyn}"
CAP="$1"; shift
CMD="$*"
# /usr/bin/time -v gives the peak RSS even when the child is OOM-killed by the cgroup.
out=$(systemd-run --user --scope -q -p MemoryMax="$CAP" -p MemorySwapMax=0 \
        /usr/bin/time -v "$QSYN" -q --no-version -c "$CMD" 2>&1)
rc=$?
peak=$(printf '%s\n' "$out" | rg -o 'Maximum resident set size \(kbytes\): [0-9]+' | rg -o '[0-9]+$')
printf 'rc=%s peak=%s\n' "$rc" "${peak:-NA}"
printf '%s\n' "$out" | rg -v 'Maximum resident' | tail -3 >&2
