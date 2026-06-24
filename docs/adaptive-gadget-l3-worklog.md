# L3 worklog -- parallel phasepoly executor

Plan: `docs/plan-gadgetization-l3.md`. Scope: seq gadgetize + parallel phasepoly, fastTodd untouched,
`--threads 1` byte-identical to L2.

## Progress
- Extracted `phasepoly_region_bytes(n, m_region)` from `predict_peak_bytes`. No behavior change.
- Threaded executor `optimize_phase_polynomial(Tableau&, strategy, num_threads, budget)`: builds
  independent Clifford-context units (disjoint subtableau indices), reserve-before-allocate budget
  gate (mutex+CV), std::thread pool. `num_threads<=1`/no-budget falls back to the sequential overload.
  Refuses up front if one region exceeds the budget; aborts all on failure/stop.
- Wired phasepoly subparser: `-j/--threads`, `-m/--memory-limit`, `--adaptive`. `--threads>1` requires
  a budget.
- Tests `tests/src/tableau/phasepoly_parallel.cpp`: byte-identical to sequential, thread-count-
  invariant (2 vs 8), budget-too-small refused with tableau untouched. Full suite green (55 cases).

## Validation (both done -- L3 closed)
Scripts: `scripts/cgrun-phasepoly-parallel.sh` (no-OOM + thread-invariance under a cgroup `memory.max`
cap == budget; caps the budgeted gadgetize+phasepoly only, reads T-count in a separate uncapped run
since the tableau->qcir export is an unbudgeted transient), `scripts/headroom-sweep.sh` (cross-allocator
headroom floor under threads). Results: `bench2x2/l3-cgroup/`, `bench2x2/l3-headroom/`.

1. No-OOM, k threads, large case. hwb10 @ 256M, threads {1,2,4,8}: every count `rc!=137` and T-count
   thread-invariant. Peak grows with k (fragmentation compounds): at the shipped 175 headroom, threads=8
   margin +9.0%. hwb11 @ 512M confirms at 2x scale.

2. Headroom re-validation under threads, per allocator (resolves the `adaptive_gadget.cpp` FIXME).
   glibc is the binding allocator: it retains freed gadgetize arenas, so the gadgetize and phasepoly
   stages SUM in RSS rather than `max()`. At the old default 150, threads=8 left predicted < actual
   (UNSAFE) on hwb10 @ 256M and was killed under the cap on hwb11 @ 512M. jemalloc / mimalloc2 / mimalloc3
   were already safe at 150 (mimalloc3 the *safest* -- the stale note claiming it was worst was wrong).

| allocator | safe at 150 (threads=8)? | floor |
|-----------|--------------------------|-------|
| glibc     | no (predicted<actual)    | 175   |
| jemalloc  | yes                      | 150   |
| mimalloc2 | yes                      | 150   |
| mimalloc3 | yes                      | 150   |

   Default raised 150 -> 175 (`gadget_alloc_headroom_percent`): glibc predicted>=actual with +8.8%
   (hwb10@256M) / +4.4% (hwb11@512M); others >=17%. The glibc margin compresses with scale, so very large
   budgets should pin a higher `QSYN_GADGET_HEADROOM`. The structural alternative to padding -- a
   `malloc_trim()` at the gadgetize->phasepoly boundary to return arenas and break the stage-sum -- is
   deferred (would let glibc drop back toward 150 and recover the slight T-count cost of under-merging).

## Deferred (out of L3 scope)
- fastTodd inner-threading (`docs/plan-multithread-fasttodd.md`); time-budgeted merge; parallel
  gadgetize (this plan's TODO). The min-T ablation scaffold (`GadgetizeMergeStrategy`, `--strategy`) is
  intentionally KEPT -- a time-budgeted merge variant may reuse it.
