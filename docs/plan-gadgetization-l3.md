# L3 - parallel phasepoly executor

Executor shipped (see worklog). L1+L2 shipped. Forward-only.

## Rules
- No ASCII-code.
- No leak on internal processes.
- Place worklog, be brief.
- Do not pile trash up.
- Follow the coding style in the repository as found (Doxygen `/** @brief @param @return */` on function definitions; ASCII-only, TeX `\pi` not unicode).
- Reproducible only: committed scripts/tests, fixed `--memory-limit` (never `--adaptive`, which is machine-dependent); no throwaway ad-hoc measurements.
- No invented/magic defaults; parameters explicit.
- No hallucination: verify against the code, do not trust stale docs, state only what is checked.
- Decide by ablation/data, not by guess.

## Goal
Run the gadgetized `phasepoly` step on the independent regions in parallel, within budget M. Additive:
gadgetize stays sequential and unchanged; only phasepoly is threaded; no `fastTodd.cpp` changes;
`--threads 1` reproduces L2 exactly (per-region optimize is pure).

## Memory
`peak = max(gadgetize_peak, baseline + sum over the k concurrent regions of their phasepoly L-matrix)`.
- `kWidenCopies` is unchanged - it is gadgetize-stage and freed before phasepoly, so it never enters the
  threaded term.
- L2's phasepoly region-max becomes a k-region sum.

## Admission (mandatory; check + reserve must be ONE atomic step)
A bare atomic increment is NOT safe: two threads can both see headroom, both add, and overshoot M -> OOM
(check-then-act race). The add must commit only if it keeps `committed <= M`, indivisibly:
- mutex: `lock; if committed + predict(r) <= M { committed += predict(r); admit } else wait on CV; unlock`.
- or CAS loop: load `cur`; if `cur + predict(r) > M` then wait (do NOT add); else
  `CAS(committed: cur -> cur + predict(r))`, retry on failure.
Release on completion (notify the CV). Planner sets `k = max k with sum(k largest predicted) <= M`. With
predict >= actual per region, `sum(actual) <= sum(predict) <= M`, so no OOM. One OOM kills all threads.

## Allocator
- `MALLOC_ARENA_MAX=1` is out (serializes malloc -> no speedup). Enforce via cgroup `memory.max`.
- `QSYN_GADGET_HEADROOM` must be re-validated under threads, per allocator (k allocations compound
  fragmentation; mimalloc/jemalloc were +12-18% single-thread).
- Verify the per-thread footprint: disjoint tableau slice (1 + k L-matrices) vs private region copy.

## Speedup
Region-count-bound: per-region phasepoly wall ~5x per 2x region size, so parallelize across regions, not
within one - a single huge merged region does not parallelize.

## Tests
`--threads 1` byte-identical to L2; no-OOM under cgroup `memory.max` with k threads; results
thread-count-invariant.

## TODO -- parallel gadgetize (deferred)
Gadgetize is parallelizable across merge-groups (disjoint ranges), but L3 keeps it sequential. When
revisited: replace the running ancilla counter with a per-group prefix-sum base, write output by
index, extend admission control to gadgetize (its peak-halving move-out trick assumes sequential, so
concurrent groups sum their peaks), and confirm the per-group ops are free of global mutable state.
