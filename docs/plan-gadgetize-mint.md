# Gadgetize min-T strategy -- ablation

Not started. Decides the gadgetize merge strategy by experiment, not by guess.

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

## Problem
Gadgetize exists to cut T-count: merging Hadamard-bounded regions feeds the phase-poly optimizer larger
polynomials with more cancellation. The shipped planner (`plan_boundaries` in `adaptive_gadget.cpp`)
instead targets memory -- greedy z-overlap merge under a budget -- so T-count only floats as a side
effect. No ablation was run to justify it.

## Goal
Fix the memory budget as a hard cap. Make T-count the objective. Choose the merge strategy by ablation:
run budget-respecting strategies under the *same* cap, measure final T-count, pick the lowest.

## Strategies (deletable scaffolding; each keeps every merged region's predicted peak <= M)
1. z-overlap -- current: merge the boundary whose regions share the most z-active qubits.
2. max-merge -- pack adjacent regions into the largest groups that each fit M (fewest regions).
3. max-terms -- merge to maximize Pauli-rotation terms per budgeted region.
4. balanced-fill -- fill regions evenly up to M, not greedily maxing one.

## Structure
- One isolated file: a merge-strategy interface (mirrors the phasepoly strategy pattern), the four
  impls, a name->strategy selector. Confined so teardown is a clean delete.
- `plan_boundaries` stays the default; new optional `gadgetize --strategy <name>` routes to a variant.
- No phasepoly / L3 / `fastTodd.cpp` change.

## Ablation
Committed reproducible script (fixed budget, never `--adaptive`). Per circuit at its fixed gadgetize
budget, sweep the four strategies, record T-count (primary), peak RSS, wall time. Phasepoly sequential.
Matrix:
- hwb8 (`bench2x2/sani/hwb8__qsyn_raw.qc`) @ 256M
- ham15-high (`bench2x2/sani/ham15-high__qsyn_raw.qc`) @ 64M

Budgets sit in the multi-region regime (gadgetize splits into many bounded regions, so strategy
matters). Above it the circuit full-merges to one region and all strategies coincide.

## Correctness
Gadget output is measurement-based -> not `qcir equiv` checkable. Reuse the existing adaptive-gadget
checks: diagonal-region invariant + reference T-count (`tests/src/tableau/adaptive_gadget.cpp`).

## Teardown
After the winner is chosen: delete the losing strategies, the `--strategy` flag, and the selector file;
fold the winner into `plan_boundaries`. Keep the ablation script + results table as the record.

## Out of scope
Out-of-core / streaming fasttodd (the other route to min-T without under-merging). Any L3 change.
