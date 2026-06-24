# Memory-budgeted partial Hadamard gadgetization + region threading - design

**Date:** 2026-06-22 (revised; supersedes the 2026-06-19 draft)
**Status:** L1 + substrate SHIPPED (commit 2a35608d); L2/L3 pre-implementation
**Branch:** feat/qsyn-adaptive-gadget (off fix/fasttodd-final)

> **L2 revision (2026-06-22, supersedes the design below where they conflict).** Locked with the
> user during L2 planning:
> 1. **No new `adaptive-gadget` command and no fused `adaptive_gadgetize_optimize` function.**
>    Instead, **extend the existing `gadgetize` command** with `--memory-limit`/`--adaptive`.
>    Gadgetization and `phasepoly` (FastTODD) stay **strictly separate, independently-runnable
>    steps** (the L1 ablation separation), so L3 cannot re-couple them.
> 2. **`--memory-limit` grammar:** absent -> FULL (current L1 all-internal, unchanged); `<size>`
>    (K/M/G) -> adaptive selective gadgetization with `<size>` as the **planner ceiling**; `auto`
>    (alias `--adaptive`) -> ceiling = available RAM.
> 3. **Enforcement is predict-only.** Exceeding `M` triggers the OOM killer -> whole process dead,
>    uncatchable; there is **no runtime cancel/monitor/fallback** ("try then put back" is
>    impossible). The planner's prediction is the sole safeguard; the OS (`ulimit -v` / cgroup
>    `memory.max`) is the external enforcer. Nothing in `fastTodd.cpp` / `optimize_phase_polynomial`
>    is touched. **macOS exception:** `RLIMIT_AS` is ignored there, so `--memory-limit` is
>    predict-only with no OS backstop on macOS.
> 4. **Memory model (corrected by the calibration gate, 2026-06-22; see
>    `docs/adaptive-gadget-worklog.md`).** The peak is **TOHPE**, not `fast_todd_move`:
>    `tohpe_once_iteration` holds **two** L-sized matrices at once (`l_matrix` + `l_matrix_transpose`,
>    each ~= `m*C(n,2)`) plus `augmented=identity(m)` and transient maps. The **dominant term is
>    `2*m*C(n,2)`, blowing up through `n`** (gadgetizing adds +1 qubit per internal H, so `C(n,2)`
>    explodes); the `m^2` term is secondary (high-m/low-n only). Real peak runs ~2-4x a single-matrix
>    estimate, and big full-gadget runs **TLE before MLE**, so exact fitting is blocked -> use the
>    analytic model + an empirical `FRAG` factor and a **conservative util target (~60-70%, NOT
>    99%)**. Coefficient a~=1 (one byte/bit) confirmed.
> 5. **L3 multithread:** the dispatcher parallelizes only the already-gadgetized `phasepoly` step
>    across independent regions, with **mandatory atomic reserve-before-allocate** admission (k
>    concurrent regions sum their peaks; one OOM kills all threads).

## 1. Problem & goal

Given a **fixed memory budget `M`**, answer two coupled questions for a Clifford+T
circuit: **how low can the T-count go, and how fast can it be produced?**

The two levers both spend the *same* resource - memory:

- **Gadgetization** (quality lever). Absorbing an internal Hadamard into an ancilla
  qubit merges two adjacent Hadamard-free regions into one larger phase polynomial,
  giving FastTODD a bigger 3-tensor-rank problem and more cancellation opportunity ->
  lower T. Cost: the merged region's qubit count `n` grows by +1 per absorbed
  Hadamard, and FastTODD's `L` matrix has `~C(n,2)` rows -> memory grows ~quadratically
  in `n`.
- **Threading** (speed lever). Hadamard-free regions are already independent
  (`strategy.optimize()` is a pure function of its `(clifford, polynomial)` - verified),
  so they can be optimized in parallel. Cost: each region in flight holds its own `L`
  matrix -> peak memory grows ~linearly in the number of concurrent regions.

The endpoints are known:

- **Gadgetize nothing** = qsyn today: memory-light (~5.5 GB even on big circuits), but
  higher T (optimizes each Hadamard-bounded region in isolation).
- **Gadgetize everything** = Vandaele's pipeline: T-count floor, but blows up to
  TLE/MLE on Hadamard-dense circuits (e.g. hwb11 ~14k Hadamards -> ~14k ancillas ->
  `L` matrix exceeds 52 GB). See `docs/r14944052-worklog.md` regime 3.

This feature lives **between** those endpoints: spend `M` where it most reduces T,
and use leftover headroom to go faster, never exceeding `M`.

### Decisions locked during brainstorming

| Decision               | Choice                                                                         |
|------------------------|--------------------------------------------------------------------------------|
| Deliverable            | `gadgetize --memory-limit` on the existing command (gadget-only)               |
| Budget enforcement     | **Predict-only**; external OOM (`ulimit -v`/cgroup) enforces; no runtime guard |
| Gadgetization ordering | **Benefit-ranked** (knapsack), not blind "absorb until full"                   |
| Benefit estimation     | **Static proxy** (no speculative trial runs during planning)                   |
| Thread scheduling      | (L3) atomic predict-based admission + LPT over region costs                    |
| Approach               | Separate steps (gadgetize || phasepoly), 3 independently-shippable layers       |

## 2. Architecture

### 2.1 Command interface

**Extend the existing `gadgetize` command** (`OptimizationMethod::gadgetize` in
`src/cmd/tableau_cmd.cpp`) with a memory-budget option -- no new command:

```
tableau optimize gadgetize [--memory-limit <size>|auto] [--adaptive]
```

- absent (`tableau optimize gadgetize`) -> **FULL** gadgetization = exactly L1's current
  all-internal behavior, **unchanged**.
- `--memory-limit <size>` (`K`/`M`/`G`) -> **adaptive**: the planner gadgetizes the most it can
  while keeping the predicted later-`phasepoly` peak under `<size>` (the **planner ceiling**).
- `--memory-limit auto`, alias `--adaptive` -> ceiling = available RAM (`available_memory_bytes()`).

`gadgetize` stays **gadget-only**; `phasepoly` (FastTODD) is run **separately** by the user, so the
two levers remain independently selectable for the ablation study (the L1 separation). There is **no
`adaptive-gadget` command and no fused `adaptive_gadgetize_optimize` function** -- gadgetization and
optimization are never combined into one call (so L3 cannot re-couple them). Threading is a future
`phasepoly`-side concern (sec 5), not a `gadgetize` flag.

**User-facing surface carries no internal jargon.** `--help`, default-level logs/errors, and the
procedure label use only: memory budget, T-count, qubits, ancillae. Never surfaced: knapsack,
benefit proxy, merge graph, Vandaele, FastTODD/TODD, L-matrix, the constant `K`, the layer numbers
L1/L2/L3. Those live in code comments and this doc only; planner internals sit in an anonymous
namespace.

### 2.2 Module

Existing `src/tableau/optimize/adaptive_gadget.{hpp,cpp}` (where L1 lives). L2 adds, additively, a
**subset-selective overload of the existing pure transform** -- *not* a combined optimize entry:

```cpp
// L1 (existing): gadgetize EVERY internal Hadamard (Full).
bool gadgetize_internal_hadamards(Tableau& tableau);

// L2 (new, additive): gadgetize only the chosen boundaries (the planner's subset).
bool gadgetize_internal_hadamards(Tableau& tableau,
                                  std::span<size_t const> selected_boundaries);
```

Both overloads are **pure `Tableau -> Tableau` transforms with no FastTODD inside.** The planner
(sec 3) that computes `selected_boundaries`, and the standalone reusable memory model it consumes,
live in the anonymous namespace of `adaptive_gadget.cpp`. `optimize_phase_polynomial` and the
`strategy.optimize(clifford, polynomial)` contract are **untouched**; `phasepoly` runs as its own
separate step.

### 2.3 Pipeline placement

`hopt` (`minimize_internal_hadamards`) runs first, then the (now budget-aware) `gadgetize`
transform, then `phasepoly` as a **separate** step: `tmerge -> hopt -> gadgetize [--memory-limit]
-> phasepoly fasttodd`. `hopt` first because it cheaply lowers the Hadamard count, so the budget is
spent on fewer, higher-value merge decisions.

### 2.4 Data flow (two separate steps, predict-only budget)

```
Tableau: C0 P0 C1 P1 C2 ... Pk      (Pi = diagonal phase-poly regions;
                                      Ci = Clifford blocks holding internal Hadamards)
   |
   v  gadgetize [--memory-limit M]  (single command)
   |    Planner (cheap, no optimization runs): choose selected_boundaries so the predicted
   |    later-phasepoly peak of every resulting region <= M*safety (predict-only; M never reached
   |    because exceeding it = OOM kill, uncatchable, no fallback).
   |    Then the pure subset-selective transform gadgetizes exactly those boundaries.
   v
Tableau' (some regions merged+gadgetized over n+k qubits, rest left non-gadget; ancillae |0>)
   |
   v  phasepoly fasttodd   (separate command; optimize_phase_polynomial, UNCHANGED)
   v
Tableau'' (optimized per region; peak stays under M because the planner sized the regions)
```

### 2.5 Three build layers (each shippable & testable)

- **L1 - Gadgetization core. SHIPPED (commit 2a35608d).** Substrate `extended_with_ancillae` + the
  pure all-internal gadgetize transform + bounded interior-merge + the structural/signature
  correctness net, via the `gadgetize` command (no planner). L2 adds a *selected-subset* overload.
- **L2 - Selection + planner + predict-only budget. SHIPPED + verified (commit pending, 2026-06-23).**
  Memory model SOLVED (profile-derived `predict_peak_bytes`, sec 3.2); enforcer = cgroup `memory.max`;
  hwb10 12G no longer OOMs (7.5 GB < 12 GB); predict >= actual on ham15-low/med/high, hwb6, tof3,
  hwb10, hwb11. Original text below.
- **L2 - Selection + planner + predict-only budget.** Add the **standalone reusable** memory model
  + static-proxy benefit ranking + knapsack selection + the subset-selective transform overload,
  wired into `gadgetize --memory-limit`. **Budget is predict-only:** the planner keeps every
  region's predicted peak under `M*safety`; exceeding `M` = OOM kill with **no recovery**, so there
  is **no runtime monitor, no cooperative cancel, no fallback** (impossible -- dead before you can
  put anything back). `M` is enforced externally (`ulimit -v` / cgroup); the margin is calibrated
  empirically (sec 6 calibration gate), not guessed; target ~99% utilization. `phasepoly` stays a
  separate, unchanged step. **After L1+L2 you can answer "how low can T go under `M`."**
- **L3 - Parallel `phasepoly` executor (additive).** A `RegionDispatcher` parallelizes only the
  **already-gadgetized** `phasepoly` step across the independent regions -- it must **not** re-fuse
  gadgetization with optimization. **Multithread memory hazard:** with k concurrent regions the
  peak is the *sum* of their peaks, and one OOM kills the whole process (all threads). So admission
  is **mandatory atomic reserve-before-allocate**: `committed_bytes += predicted(r)` only if it
  stays `<= M` (CAS/lock -- never check-then-start), released on completion; a worker waits on a CV
  when there is no headroom. L2's ~99% single-thread fill does **not** transfer (aggregate model
  error compounds kx, and summed peaks rarely coincide), so the planner sets `k = max k with
  sum (k largest predicted peaks) <= M`. The per-region optimize path has no shared mutable state
  (verified), so L3 needs **zero `fastTodd.cpp` changes**; `--threads 1` reproduces L2 exactly.

## 3. Planner (Phase 1)

### 3.1 Region & merge graph

After `hopt`/collapse the tableau is an alternating sequence of Clifford blocks `Ci`
and diagonal phase-poly regions `Pi`. Adjacent regions `Pi`, `Pi+1` are separated by
`Ci+1`, which contains some number `h(i)` of internal Hadamards. **Merging `Pi` and
`Pi+1`** means gadgetizing the Hadamards in `Ci+1`: each costs +1 ancilla and the two
regions union into one diagonal phase polynomial over the extended qubit set.

Candidate merges form a path graph over regions; the planner selects a subset of
boundary edges to contract.

### 3.2 Memory model -- the WHOLE gadgetized-region pipeline (calibrated)

**What is budgeted is the entire bundled gadgetize->phasepoly->convert-back cost of a region -- the
"monster" that distinguishes Vandaele from qsyn -- not FastTODD's matrices alone.** Root cause of
every sink is the same: **gadgetization inflates `n` (+1 qubit per internal H), and `n` appears,
often squared, in EVERY post-gadget stage.** Budgeting only FastTODD leaves the rest unbounded --
exactly how full-gadget hwb10 hit 56 GB -> OOM -> WSL2 restart in the *collapse* step, before FastTODD
even ran (incident 2026-06-23). The region's predicted peak is the **max over all post-gadget
stages** (same qsyn process; the merged-region tableau is resident throughout):

**SOLVED 2026-06-23 (profile-derived, supersedes the box below).** `predict_peak_bytes(n, s_clifford,
m_total, m_region)` in `adaptive_gadget.cpp` is a `max` over stages of EXACT sizeof/paper terms (every
constant is a `sizeof` or a paper dimension; the ONE flagged-empirical knob is `kFragPercent`):

```
floor = current_rss_bytes()            // LIVE plan-time RSS (QCir DAG ranges 13MB tof3 .. 140MB hwb11), not a constant
predict = floor + max(
  GADGETIZE = (3*(s_clifford*stab + m_total*rot) + 24*(2n)^2) * 150/100,  // ext+middle+out copies x glibc frag
  PHASEPOLY = 2*m_region*(n + C(n,2)) + m_region^2 + n*m_region )         // FastTODD L+Lt+augmented+A, 1 byte/bit
//  stab = 2n*(2n+1)/8 + 2n*32   (bit-packed dynamic_bitset + 32B PauliProduct struct)
//  rot  = (2n+1)/8 + 48 ;  s_clifford = #StabilizerTableau blocks ;  n = base_n + total selected ancillae (GLOBAL)
```
The old "40x unidentified allocation" was just the widened tableau in ~3 coexisting full-width copies;
the transform now `std::move`s them (byte-identical, ~halved peak). The historical multi-stage box
(below) was the right intuition; this is its measured form. Enforcer = cgroup `memory.max` (RSS).

```
predict_region_bytes(n, m) = BASELINE + FRAG * max(   // SUPERSEDED -- see solved model above
    WIDEN(n, m),          // extended_with_ancillae(k): stabilizer tableaux O((n)^2) + rotations O(m*n)
    COLLAPSE(n, m),       // collapse + extract_clifford_operators (Clifford synthesis) O(n^2) op-strings  <-- the 56GB hwb10 sink
    FASTTODD(n, m),       // = 2*m*(n + C(n,2)) + m*m + c*m*n   (TOHPE/FastTODD, calibrated below)
    SYNTH_OUT(n, m),      // convert tableau->qcir AFTER phasepoly: O(n^2) Cliffords + rotations -> QCir DAG
    RESIDENT(n, m)        // merged tableau held through phasepoly: ~m*n + O(n^2) boundaries
)
```
Per region: `n` = qubits incl. accumulated ancillas, `m` = terms, `C(n,2)=n(n-1)/2`. (`tmerge`/`hopt`
and input `convert` run *before* gadgetize, at small `n` -> not sinks.)
- **`FASTTODD` details (calibrated 2026-06-22).** Peak is **TOHPE** (`tohpe_once_iteration`,
  `fastTodd.cpp:311-400`), run before `fast_todd_move`: two L-sized one-byte-per-bit matrices at once
  -- `l_matrix` (`:322`) + `l_matrix_transpose` (`:323`), each ~= `m*C(n,2)` -- plus `augmented=
  identity(m)` (`:325`, m^2), `table`+`phase_poly_matrix` (~= `m*n`), transient maps (`:338-356`).
  FastTODD is **read-only here** -- modeled, not modified; no `--memory-limit` on it.
- **Dominant term `2*m*C(n,2)`, via `n`-inflation.** `m^2` is secondary (high-m/low-n only). `a~=1`
  (one byte/bit) confirmed; real peak runs ~2-4x a single-matrix estimate, and big full-gadget runs
  **TLE before MLE** so exact fitting is blocked -> `FRAG` (provisional ~=1.5) is the empirical knob,
  **util target conservative (~60-70%), NOT 99%** (overshoot = OOM death). `BASELINE` ~= 15-55 MB.
- Implement as ONE **standalone reusable module** (planner + L3 admission accountant both call it).
  The non-FastTODD stages (`WIDEN`/`COLLAPSE`/`SYNTH_OUT`/`RESIDENT`) need their allocation constants
  measured at implementation; refine all of `FRAG`/`BASELINE`/per-stage constants with the feature's
  own completing partial-gadget runs.

A **hard per-region constraint**: every region must satisfy `predict_region_bytes(region) <=
M*safety` (`safety` ~0.6-0.7), i.e. **the worst stage fits**. A merge that would violate it -- through
*any* stage (widen, collapse, FastTODD, output synthesis, or residency) -- is never selected. Because
exceeding `M` = OOM kill with no recovery, the prediction is the *only* safeguard.

### 3.3 Static benefit proxy

Estimate a merge's T-reduction *without* running FastTODD. Candidate proxy (tunable):

> **score(merge) = number of T-rotation parities that become combinable across the
> boundary** - i.e. count pairs `(p in Pi, q in Pi+1)` whose Pauli `z`-supports overlap
> (or are equal/complementary) once the gadget maps both regions into the common
> extended space.

Rationale: gadgetization helps precisely when T-rotations on the two sides share
parity structure that a single phase polynomial can cancel but two separate ones
cannot. The exact proxy is a calibration target (validate offline against trial-based
measurements on a handful of circuits); the architecture does not depend on the
specific formula.

### 3.4 Joint selection (knapsack + concurrency)

Coupled: more merges -> larger max region -> fewer regions can run concurrently within
`M`. Greedy heuristic (explainable; not provably optimal):

1. Compute base region costs; drop any merge whose merged region exceeds `M` (sec 3.2).
2. Sort remaining candidate merges by **benefit / delta-bytes** descending.
3. Apply merges in that order while the plan stays *feasible*.
   - **L2 feasibility (single-threaded):** just "max single region `<= M`."
   - **L3 feasibility:** there exists a thread count `k >= 1` with
     `sum(k largest region costs) <= M`.
4. **(L3 only)** set `k = min(max_threads, largest k with
   sum(k largest region costs) <= M)`.
5. **(L3 only)** emit schedule order = regions sorted largest-cost-first (LPT, for load
   balance).

Output: `MergePlan { merges, partition, per-region bytes, [k, order] }` (the `k`/`order`
fields are populated by L3; L2 leaves them empty and runs serially).

## 4. Gadgetization mechanism (Phase 2a)

Port Vandaele's `hadamard_gadgetization` (`quantum-circuit-optimization/src/circuit.rs:174`)
**selectively** - apply only to the Hadamards chosen by the planner, instead of all
internal Hadamards.

Per selected internal Hadamard on qubit `q` (semantic oracle = the Rust gadget):

- allocate a fresh ancilla `a`;
- emit the fixed gadget `S(a) . S(q) . CX(q,a) . S(a) . Z(a) . CX(a,q) . CX(q,a)`
  (no Hadamard inside the body);
- the data wire **stays on `q`** - subsequent operations are *not* remapped. The oracle
  only *records* `parent_ancilla[q] = a` into its ancilla->parent map for
  de-gadgetization bookkeeping; later non-H gates keep their original qubit indices
  (`circuit.rs:197-202`). L1 already does this ("data wire stays on q (no remap)");
- wrap the whole region: ancilla *prep* before the body, ancilla *uncompute* after,
  mirroring Vandaele's `c_out = anc + body + anc`.

Only Hadamards strictly between the first and last T are gadgetizable (leading/trailing
ones are not internal), matching the Rust `i < last && flag` guard.

**Integration note (main implementation risk).** qsyn's pipeline is tableau-native,
while the reference gadget is circuit-level. The likely implementation reconstructs the
affected span across the two regions, applies the gadget, and rebuilds a single
diagonal `vector<PauliRotation>` region over the extended qubit set. This is the
hardest correctness piece; pin it with golden equivalence tests (sec 6) using the Rust
function as the oracle.

### De-gadgetization & equivalence

The ancillae are **kept** (output is `n+k` qubits), each initialized |0>, exactly as the
oracle does (`c_out = anc + body + anc` keeps `anc`); the result is equivalent to the
input on the original qubits whenever the ancillae are |0>. They are not uncomputed back
to `n` qubits. **The gadgetized output is measurement-based, not unitary, so
`qcir equiv` cannot verify it** (see `adaptive-gadget-measurement-based.md`). The
gadgetized path is therefore gated by three qsyn-internal checks, not `qcir equiv`:

- **structural** - the merged region is diagonal (no trailing-H thread-through);
- **signature tensor** - the existing `has_same_signature_tensor` guard, applied per
  region in the extended space *before* de-gadgetization (FastTODD already runs this;
  gadgetized regions reuse it for free);
- **reference T-count** - `<=` the paper/Vandaele T on the golden circuits.

`qcir equiv 0 1` is still used, but only on the **non-gadget unitary path** (plain
`fasttodd`), where it is valid.

**The per-region default is plain non-gadget fasttodd.** The planner only gadgetizes the boundaries
it selects; every un-selected region stays the L1 non-gadget baseline (unitary, `qcir
equiv`-checkable). There is **no runtime fallback** -- see sec 5 (predict-only; OOM is fatal).

## 5. Budget enforcement (predict-only)

There is **no executor, no monitor, no `RegionDispatcher`, no runtime cancellation** in L2.
`gadgetize` runs its planner + the pure subset-selective transform and stops; `phasepoly` runs
afterward as the separate, unchanged `optimize_phase_polynomial`. Budget is guaranteed entirely by
**prediction**: the planner keeps every region's `predicted_bytes <= M*safety` (sec 3.2), so the
later `phasepoly` never reaches `M`.

### 5.1 Why predict-only (no runtime guard)

Exceeding `M` triggers the **OOM killer** -> the whole process dies, uncatchable. You cannot "try
then put back," cannot fall back to non-gadget at runtime (you are dead before you can). So the
prediction is the sole safeguard, and `fastTodd.cpp` / `optimize_phase_polynomial` are **untouched**
(no cooperative `stop_token`, contrary to the earlier draft). The real enforcer is the environment:
`ulimit -v` / cgroup `memory.max`. The margin `safety` is calibrated empirically (sec 6) just above
the model's measured worst-case residual -- tight (~99%), not padded.

**macOS:** `RLIMIT_AS` is ignored, so there is no OS backstop; `--memory-limit` is predict-only with
no hard guarantee on macOS (document in `--help`). `available_memory_bytes()` uses `sysctl
hw.memsize` + mach `vm_statistics64` there instead of `/proc/meminfo`.

### 5.2 Threading (L3, additive -- parallelizes `phasepoly` only)

A future `RegionDispatcher` parallelizes the **already-gadgetized** `phasepoly` step across the
independent regions; it never re-fuses gadgetization with optimization. The per-region optimize path
has no shared mutable state (only read-only globals: `stop_requested()`, the magic-static env cache
`fastTodd.cpp:38`), so it is re-entrant and L3 needs **zero `fastTodd.cpp` changes**.

- **Mandatory atomic admission.** With k concurrent regions the peak is the **sum** of their peaks,
  and a single OOM kills all threads. Before starting region `r` a worker must **atomically reserve**
  `committed_bytes += predicted_bytes(r)` only if it stays `<= M` (CAS/lock -- *never* check-then-
  start), releasing on completion; otherwise it waits on a CV. Reserve-before-allocate, so the live
  sum never exceeds `M`.
- **Utilization does not transfer.** L2's ~99% single-thread fill cannot hold under k threads
  (aggregate model error compounds kx; summed peaks rarely coincide), so the planner sets
  `k = max k with sum (k largest predicted peaks) <= M` and trades fill for parallelism.
- **Determinism.** Results are thread-count-invariant (per-region optimize is pure); `--threads 1`
  reproduces the L2 sequential reference exactly.

## 6. Testing

- **Calibration gate (do FIRST, on the current L1 binary).** Before building the planner, validate
  the sec-3.2 model: on `qcla_mod_7` / `ham15-med` run no-gadget vs full-`gadgetize` + `phasepoly`
  under `/usr/bin/time -v`; extract merged `(n,m)`; **bracket the true peak by ratcheting `ulimit
  -v` down** to the smallest value that still completes; fit the multi-term coefficients; the
  production `safety` margin = just above the measured worst-case residual. Record in the worklog.
- **Unit (planner).** On synthetic region-cost lists: every selected region fits `M*safety`; higher
  benefit/byte merges picked first; infeasible merges rejected; selection fills toward the ceiling.
- **Gadgetized correctness.** `tof_3`, `gf2_5`, `barenco_tof_3` via `gadgetize` (selective + full)
  then `phasepoly fasttodd` -> T <= plain non-gadget `fasttodd`, gated by structural-diagonal +
  `has_same_signature_tensor` + reference-T (**not** `qcir equiv` - measurement-based). Tiny budget
  -> gadgetizes fewer merges; a budget too small for even one merge -> **error** (never a silent
  zero-gadget result; no-gadget = omit the step). Full mode byte-identical to L1.
- **Non-gadget correctness.** `phasepoly fasttodd` without `gadgetize` is unitary -> `qcir equiv 0
  1` passes (the verifiable baseline; bare `gadgetize` = FULL is measurement-based, not this).
- **Memory (headline).** Under `ulimit -v M` / cgroup `memory.max = M`: `gadgetize --memory-limit M`
  then `phasepoly fasttodd` on a large circuit completes **without OOM** and `peak/M >= safety`.
- **Stress.** Hadamard-dense `hwb8`/`ham15-high` at a fixed budget complete **without OOM** - the
  headline win over Vandaele's gadgetize-all - recording achieved T, peak RSS, utilization.
- **Original code intact.** `phasepoly fasttodd` byte-identical to pre-branch on the golden dofiles.

## 7. Risks & open questions

1. **Benefit-proxy quality.** A static proxy may mis-rank merges. Mitigation: calibrate
   the formula offline against trial-based measurements on a few circuits; the
   architecture is proxy-agnostic.
2. **Memory-model accuracy (now safety-critical).** Predict-only enforcement means an
   under-prediction = OOM death, not graceful degradation, so the multi-term model (incl. the m^2
   `augmented` term, sec 3.2) must be accurate to ~1% on large circuits. Validated by the sec-6
   calibration gate; re-calibrate if `Row` becomes bit-packed.
3. **No core change to FastTODD.** Enforcement is predict-only; the earlier "defaulted
   `jthread::stop_token` cooperative cancel" is **dropped** -- `fastTodd.cpp` /
   `optimize_phase_polynomial` are untouched. The L3 multithread hazard (sum-of-peaks OOM) is
   handled by atomic reserve-before-allocate admission in the dispatcher (sec 5.2), not in FastTODD.
4. **Gadgetization correctness** - porting the circuit-level gadget into the
   tableau-native pipeline; oracle = Vandaele `circuit.rs:174`. Because the output is
   measurement-based, it is pinned by the structural-diagonal + `has_same_signature_tensor`
   + reference-T net, **not** `qcir equiv` (which is valid only on the non-gadget path).
5. **CMakeLists merge note (no SIMD in L1).** L1 has no SIMD - no `-march`/AVX anywhere.
   Its only `CMakeLists.txt` change is four link deps (sul::dynamic_bitset, xtl, xtensor,
   xtensor-blas) on the unit-test target, needed because the new tableau/qcir gtests pull
   headers that use them. SIMD lives on `feat/qsyn-simd-build-flags`. The two branches
   touch the same `target_link_libraries` block, so a textual merge conflict is expected
   *when they are eventually merged* - resolve it then; nothing to do until the merge.
6. **Orthogonal to bit-packing.** This feature *manages* memory; it does not shrink the
   large per-region constant factor from the one-byte-per-bit `BooleanMatrix::Row`
   (worklog Phase-2/3 future work). Both are independently worthwhile.

## 8. Out of scope (YAGNI)

- Bit-packing / SIMD for `BooleanMatrix::Row` (separate future work).
- Distributed/multi-process execution (single-process thread pool only).
- Auto-tuning the safety margin or proxy weights at runtime (fixed constants first).
- Speculative/trial-based benefit estimation (explicitly rejected: static proxy chosen).
- Making the `RegionDispatcher` permanent or general - it is explicitly temporary and
  swappable, isolating the L3 threading so it can be replaced without touching L1/L2.
- Uncomputing ancillae back to `n` qubits - the design keeps them (initialized |0>),
  matching the oracle; an `n`-qubit output is not a goal.
