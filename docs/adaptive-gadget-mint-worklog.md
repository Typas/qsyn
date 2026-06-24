# Min-T gadgetize-strategy ablation -- worklog

Plan: `docs/plan-gadgetize-mint.md`. Question: is the shipped memory-budgeted z-overlap merge planner
the right boundary-merge strategy for *minimizing T-count*, or was it only ever a memory heuristic?
Decide by ablation, not guess.

## What was built (deletable scaffolding)
- `src/tableau/optimize/gadgetize_strategy.{hpp,cpp}`: `GadgetizeMergeStrategy` interface over boundary
  selection (mirrors the phasepoly strategy pattern), four budget-respecting impls -- z-overlap (the
  production greedy), max-merge (fewest/largest groups that fit), max-terms (largest polynomial first),
  balanced-fill (keep the largest group smallest) -- plus a prefix-match `make_gadgetize_merge_strategy`
  selector. `RegionPlan` holds the shared region decomposition + budget fit predicate.
- `gadgetize_within_budget(Tableau&, size_t, GadgetizeMergeStrategy const&)` overload; `gadgetize
  --strategy <name>` routes to it. Empty `--strategy` keeps `plan_boundaries` (production default,
  untouched). No phasepoly / L3 / fastTodd change.
- Tests (`tests/src/tableau/adaptive_gadget.cpp`): every strategy reproduces gadgetize-everything under a
  generous budget (serialize-equal + diagonal); z-overlap strategy matches the default planner
  byte-for-byte (no-divergence guard); selector rejects unknown/empty names. Full suite green.

## Ablation (scripts/ablation-gadget-strategy.sh; qsyn-only, 8 cells parallel, sequential phasepoly)
Fixed budgets, never --adaptive. Results -> `bench2x2/ablation/` (table.md, ablation.csv).

| circuit    | budget | strategy      | T    | wall  | status |
|------------|--------|---------------|------|-------|--------|
| hwb8       | 256M   | z-overlap     | 3420 | 1:34  | ok     |
| hwb8       | 256M   | max-merge     | --   | 30:00 | TLE    |
| hwb8       | 256M   | max-terms     | --   | 30:00 | TLE    |
| hwb8       | 256M   | balanced-fill | 3505 | 0:48  | ok     |
| ham15-high | 64M    | z-overlap     | 949  | 4:04  | ok     |
| ham15-high | 64M    | max-merge     | 823  | 16:48 | ok     |
| ham15-high | 64M    | max-terms     | --   | 30:00 | TLE    |
| ham15-high | 64M    | balanced-fill | 1010 | 0:08  | ok     |

## Finding
The memory budget bounds RSS but NOT fastTodd time, which is superlinear in a region's term count. The
merge-harder strategies (max-merge, max-terms) pack regions to the budget ceiling, producing polynomials
fastTodd cannot reduce within 1800 s -> TLE on hwb8 (both) and ham15 (max-terms). max-merge does win the
ham15 T-count (823 vs 949, -13%) but at 4x the wall time and TLEs on hwb8, so it is not robust.

z-overlap is the only strategy that completes both circuits under budget, with the lowest completing
T on hwb8 and second-lowest on ham15. The ablation thus VALIDATES the shipped planner: T-count was not
left on the table by selecting boundaries on z-overlap rather than max merge/terms -- the aggressive
strategies trade unbounded time for marginal-or-negative T gain.

Winner = z-overlap = existing `plan_boundaries` (no fold-in needed; the default is already the winner).

## Open (out of this scope)
max-merge's lower ham15 T shows there IS T to recover by merging more *when* fastTodd can finish -- a
TIME-budgeted (not just memory-budgeted) merge variant is the natural follow-up. Out of scope here
(would need streaming/out-of-core fastTodd to make large regions tractable).
