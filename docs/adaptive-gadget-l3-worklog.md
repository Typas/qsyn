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

## Next
- No-OOM under cgroup `memory.max` with k threads on a large case (hwb10).
- Re-validate headroom percent under threads per allocator (fragmentation compounds).
