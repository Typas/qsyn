# Adaptive-gadget L2 -- worklog

Branch: `feat/qsyn-adaptive-gadget`. Raw, code-grounded source of truth for L2 (budget-aware
selective gadgetization), in the style of `docs/r14944052-worklog.md`. Dated entries; exact
`file:line` refs; **VERIFIED** vs **ASSUMED** marked explicitly; later entries supersede earlier
ones where they conflict. The polished design is `docs/plan-gadgetization.md` (L2 revision banner).

All line numbers refer to the state at the dated entry; re-verify before relying on them.

---

## 2026-06-22 -- L2 design locked (supersedes the `docs/plan-gadgetization.md` body)

Decisions, with the reasoning, as agreed with the user during planning:

1. **Extend the existing `gadgetize` command** with `--memory-limit`/`--adaptive`; **no** new
   `adaptive-gadget` command, **no** fused `adaptive_gadgetize_optimize` function. Gadgetization and
   `phasepoly` (FastTODD) stay strictly separate, independently-runnable steps -- the L1 ablation
   separation (`adaptive_gadget.hpp:45-50`). Rationale: a fused command/function silently re-couples
   the two levers; keeping them separate is what makes L3 an additive change, not a rewrite.
2. **`--memory-limit` grammar:** absent -> FULL (current L1 all-internal, byte-identical, unchanged);
   `<size>` (K/M/G) -> adaptive selective gadgetization, `<size>` = planner ceiling; `auto` (alias
   `--adaptive`) -> ceiling = available RAM.
3. **Predict-only enforcement.** VERIFIED with the user: exceeding `M` triggers the OOM killer ->
   whole process dead, uncatchable. "Try then put back" is impossible -- you are dead before you can
   fall back. So there is **no runtime monitor, no cooperative cancel, no fallback**; the planner's
   prediction is the only safeguard. The OS (`ulimit -v` / cgroup `memory.max`) is the external
   enforcer. Consequence: **nothing in `fastTodd.cpp` / `optimize_phase_polynomial` is touched** --
   the earlier "defaulted `jthread::stop_token`" plan (old doc sec 5.1 / risk #3) is **dropped**.
4. **Tight, not conservative.** User wants >=99% utilization on large circuits, not a 75-80% waste.
   Since overshoot is fatal, this forces an *accurate* model, and a margin calibrated empirically
   (see calibration gate below), never guessed. ~~99% is a target pending the gate.~~ (resolved
   2026-06-23, see the "MEMORY MODEL SOLVED" entry)
5. **Multithread-friendly so L3 is additive.** L3 parallelizes only the already-gadgetized
   `phasepoly` step; mandatory **atomic reserve-before-allocate** admission (k concurrent regions
   sum their peaks; one OOM kills all threads). Planner sets `k = max k with sum (k largest predicted
   peaks) <= M`. `--threads 1` = exact L2.

### VERIFIED facts from reading the code (2026-06-22)

- **m^2 `augmented` finding (key).** `fast_todd_move` (`src/tableau/optimize/fastTodd.cpp:484-512`)
  allocates, persistently and one-byte-per-bit:
  - `matrix` (L^T): `n_terms x width`, `width = n + C(n,2)` (`fastTodd.cpp:490-497`);
  - **`augmented = dvlab::identity(n_terms)`: `m x m`** (`fastTodd.cpp:498`) -- **dominates for large
    `m`**, and is the term the old `K*C(n,2)*m` model omits. Likely the real MLE driver.
  - `table` (transposed phase-poly): `m x n` (`fastTodd.cpp:484`);
  - `term_index` hashmap: `m` Rows of `n` bits + bucket overhead (`fastTodd.cpp:508-512`).
  So model: `predicted ~= a*m*(n+C(n,2)) + b*m^2 + c*m*n + d*m + e`. To be implemented as a
  **standalone reusable function** (L3 admission reuses it).
- **No shared mutable state in the optimize path** (so L3 is re-entrant, needs zero fastTodd edits):
  only read-only globals -- `extern stop_requested()` (`fastTodd.cpp:29`) and the magic-static env
  cache `tohpe_use_iteration()` (`fastTodd.cpp:37-44`, init-once/read-only). VERIFIED by scan.
- **`stop_token` IS vendored** (`vendor/jthread/{jthread.hpp,stop_token.hpp}`) and used by the CLI
  (`src/cli/cli.hpp:147,232`, SIGINT -> `cli.cpp:229` `request_stop()`). Noted only because the old
  plan leaned on it; predict-only no longer needs it.
- **L1 `gadgetize` has no no-gadget mode**: it always does FULL (`adaptive_gadget.cpp:150-224`);
  no-gadget today = don't run the command. **No dofile/golden calls the bare `gadgetize` CLI**
  (gtests call the C++ funcs directly), so changing command defaults is golden-safe.
- **Portability:** `/proc/meminfo`, cgroup, and `ulimit -v` (RLIMIT_AS) enforcement are Linux/WSL
  only. **macOS ignores RLIMIT_AS** -> predict-only there, no OS backstop; `available_memory_bytes()`
  must use `sysctl hw.memsize` + mach `vm_statistics64`. `getrusage ru_maxrss` already portable
  (`src/util/usage.cpp` has the `__APPLE__` branch).

### TODO (next): calibration gate -- *test first*, before building the planner
On the current L1 binary, `qcla_mod_7` / `ham15-med` (from the benchmark circuit inputs; chosen because
full-gadget ~=10s/14s, no-gadget ~=0s, memory diverges hard):
- no-gadget: `... tmerge; hopt; phasepoly fasttodd` under `/usr/bin/time -v` -> `peak_kb`, secs.
- full-gadget: `... tmerge; hopt; gadgetize; phasepoly fasttodd` -> `peak_kb`, secs; extract merged
  `(n,m)` (`m` from the `num_terms before FastTODD` debug line `fastTodd.cpp:736`).
- Bracket the true peak by ratcheting `ulimit -v` down to the smallest value that still completes.
- Fit the multi-term model; report residual %; set production `safety` from the measured worst case.
- **Record results here before writing any planner code.**

---

## 2026-06-22 -- Calibration gate, run 1 (current L1 binary `build/qsyn` v0.8.1)

Pipeline: `logger debug; qcir read; convert qcir tableau; tableau optimize tmerge; hopt;
[gadgetize;] phasepoly fasttodd; convert tableau qcir; qcir print --stat`. `/usr/bin/time -v`.
`m` = last `num_terms before FastTODD` (`fastTodd.cpp:736`). Circuits sanitized
(`T*->tdg, S*->sdg, Z*->z, cnot->tof`).

| circuit    | mode     | T   |  m  | n(after) | peak_kb | wall  |
|------------|----------|----:|----:|---------:|--------:|------:|
| qcla_mod_7 | no-gad   | 236 |  --  |    --     |  15028  | 0.03s |
| qcla_mod_7 | full-gad | 163 | 237 |    84    |  21968  | 25s   |
| ham15-med  | no-gad   | 212 |  --  |    --     |  14508  | 0.02s |
| ham15-med  | full-gad | 141 | 212 |    71    |  18660  | 8.2s  |

**Finding (supersedes the assumption that these circuits stress memory):** at `m~=237, n~=84` the
FastTODD matrices are tiny -- `matrix = m*(n+C(84,2)) ~= 237*3570 ~= 0.81 MB`, `augmented = m^2 ~= 56 KB`
-- **sub-MB**, dwarfed by a **~14.5-15 MB fixed process baseline** (`e`). So:
- **Time and memory do NOT correlate here.** The 25 s for qcla_mod_7 is FastTODD *compute*
  (O(m^2) candidate pairs x Gaussian elimination over `width~=C(n,2)`), not allocation. Picking
  circuits by "Vandaele ~=10 s" gives compute-heavy, **memory-light** cases.
- These runs pin only the **baseline `e ~= 14.7 MB`**; they cannot calibrate `a` (matrix) or `b`
  (m^2) because those terms are < 1 MB here. **Need large-`m` circuits** (thousands of terms) where
  the matrices reach 100s of MB-GBs. -> moving to Hadamard-dense `hwb8` (and up) for the variable-
  term calibration. `ulimit -v` bracketing is only meaningful once the variable part dominates.
- T-count win confirmed and large: 236->163 (qcla_mod_7), 212->141 (ham15-med). Gadgetization is
  clearly worth it; the budget question is real only at the high-`m` end.

## 2026-06-22 -- Calibration run 2: gadgetized region sizes + the REAL dominant term

Ran `... hopt; gadgetize; convert tableau qcir; qcir print --stat` (gadgetize-only, fast -- no
phasepoly) to get the merged region's `n` (qubits) and `m` (~= T-family) deterministically:

| circuit    | n (after gadget) |  m   | matrix `m*C(n,2)` | augmented `m^2` |
|------------|-----------------:|-----:|------------------:|---------------:|
| hwb6       |               27 |   75 |             28 KB |          6 KB  |
| ham15-med  |               71 |  212 |           0.54 MB |         45 KB  |
| qcla_mod_7 |               84 |  237 |           0.83 MB |         56 KB  |
| hwb8       |             1115 | 3517 |           2.19 GB |         12 MB  |

**CORRECTION (supersedes the 2026-06-22 "m^2 is the dominant omitted term" claim).** For the
*gadgetization* use case the dominant term is **`matrix = m*width`, `width = n + C(n,2)`**, and it
blows up through **`n`**, not `m`: gadgetizing adds +1 qubit per internal Hadamard, so hwb8's
~1100 internal H push `n` to **1115** -> `C(n,2)=621k` -> matrix ~= **2.19 GB**, while `augmented = m^2`
is only **12 MB**. So the doc's *original* `K*C(n,2)*m` instinct was directionally right; my L2
"m^2 dominates" over-correction is **wrong for this regime**. m^2 only dominates when `m >> C(n,2)`
(many terms, few qubits) -- the opposite of what gadgetization produces.
- **Model:** keep both terms -- `a*m*(n+C(n,2)) + b*m^2 + c*m*n + d*m + e` -- but note the matrix term
  `a*m*C(n,2)` is the GB-scale driver for Hadamard-dense gadgetization; m^2 matters only in the
  high-m/low-n corner.
- **Coefficient `a` ~= 1 confirmed (order of magnitude).** hwb8 full-gadget TLE'd (rc=124, 300 s)
  with peak **1.3 GB**; predicted matrix at `m=3517,n=1115` = 2.19 GB (1.3 GB ~= matrix at the
  TOHPE-reduced `m`~=2100). This rules out an 8x bit-packing error (would be 17 GB) and a /8 error
  (274 MB). The one-byte-per-bit `BooleanMatrix::Row` => memory is deterministic and predictable.
- **New constraint discovered: TIME, not memory, binds first for big gadgetized regions.** hwb8
  full-gadget TLE'd before MLE -- `fast_todd_move` is O(m^2 x reduce-over-width). So full gadgetize on
  Hadamard-dense circuits is compute-bound; the memory budget will often select *below* the
  memory-MLE point anyway because the compute is intractable. (Time is L3's "how fast"; for L2 the
  planner still budgets memory, but note the user may hit TLE before the memory cap. Out of L2
  scope to optimize, but record it.)
- **Gate status:** the central question -- *is accurate, tight, predict-only budgeting feasible?* --
  is **YES**: peak is dominated by the analytic, one-byte-per-bit `matrix = m*width` term (a~=1).
  Precise `e`/overhead fit + the `ulimit -v` tight-bracket validation still want one *completing*
  mid-size run (hundreds of MB); attempting `ham15-high` next.

## 2026-06-22 -- Calibration run 3: ham15-high -- the overhead gap (blocks naive 99%)

`ham15-high` gadgetized -> **n=351, m=1019**. Full-gadget + phasepoly: **TLE again** (rc=124, 180 s),
peak **255,696 KB ~= 256 MB**.

- Predicted `matrix = m*width = 1019*(351+C(351,2)) = 1019*61776 ~= 63 MB`; `augmented = m^2 ~= 1 MB`.
- **Observed peak ~= 256 MB ~= 4x the matrix-only estimate.** The extra is transients + TOHPE +
  process baseline + allocator fragmentation:
  - per-candidate-pair transient in `fast_todd_move` (`fastTodd.cpp:540-579`): `r_mat` = `(n+1)`
    Rows of `width` ~= `(n+1)/m x matrix` ~= 0.3-0.35x matrix (hwb8 and ham15-high both ~0.32x) --
    allocated/freed every `(i,j)` pair, a fragmentation amplifier;
  - `fasttodd_once` runs the **full TOHPE fixpoint first** (`fastTodd.cpp:639-651`) with its own
    structures, before the `matrix` is ever built;
  - process baseline (~=15-56 MB depending on circuit size).
- **Confound: both large gadget runs TLE before completing**, so the 256 MB / 1.3 GB peaks are at an
  unknown intermediate `m` (TOHPE-reduced), not the clean `m` we measured. Clean coefficient fitting
  from completing large runs is therefore **blocked** -- `fast_todd_move` is O(m^2*width), so any
  region big enough to stress memory also TLEs.

### Gate verdict (honest)
- **Feasible:** memory is deterministic and predictable in SHAPE; the matrix term dominates; a~=1.
- **But 99% tight is NOT safe as-is:** the real peak runs ~2-4x over the naive matrix term, the
  overhead is variable, and TLE blocks clean fitting. A tight 99% margin against the analytic value
  would OOM (fatal). Realistic safe utilization with the current (uninstrumented) model is more like
  **~25-40% of M** (multiplier ~2.5-4x), or we must model the transient term explicitly + a
  calibrated fragmentation factor.
- **Options to recover tightness** (need user decision, see below): (a) accept the conservative
  multiplier (low utilization, safe); (b) model transient `(n+1)*width` + TOHPE + a measured
  fragmentation factor, calibrated via (c); (c) add a **one-line `spdlog::debug` in `fast_todd_move`**
  logging `m,n` at matrix-allocation time so peak attributes to a known `m` -- exact calibration, but
  it touches `fastTodd.cpp` (log only, byte-identical behavior) which the user asked to keep
  untouched; (d) accept TLE as the practical binder and budget loosely.

## 2026-06-22 -- DECISION + final model: "model overhead analytically" (user choice)

User chose **option (b): model the overhead analytically, no `fastTodd.cpp` edit, ~60-75% util.**

**TOHPE is the real peak (VERIFIED by reading `tohpe_once_iteration`, `fastTodd.cpp:311-400`).**
`fasttodd_once` runs the TOHPE fixpoint *before* `fast_todd_move`, and TOHPE allocates **two**
L-sized one-byte-per-bit matrices at once -- `l_matrix` (`:322`) **and** `l_matrix_transpose`
(`:323`), each ~= `m*C(n,2)` -- plus `augmented_matrix = identity(m)` (`:325`, m^2), `table`+
`phase_poly_matrix` (~= `m*n` each), and a transient `std::unordered_map` keyed by `table[i]+table[j]`
(`:338-356`) holding up to O(|y|*(m-|y|)) Rows of `n` bytes. So TOHPE's footprint ~= **2*m*C(n,2) +
m^2 + maps**, which is *larger* than `fast_todd_move`'s (one matrix + augmented). This is most of the
~4x gap over the single-matrix estimate (~=2x from the two L-matrices + the maps/transients +
baseline + fragmentation).

### Final standalone memory model (to implement as a reusable function)
```
predict_region_bytes(n, m) =
    FRAG * (  2 * m * (n + C(n,2))     // TOHPE l_matrix + l_matrix_transpose  (DOMINANT, n-inflation)
            +     m * m                // augmented = identity(m)
            +  c * m * n               // table + phase_poly_matrix + term maps (c ~ 3-4)
           )
  + BASELINE                            // process baseline, ~15-55 MB (grows with circuit size)
```
- `C(n,2) = n(n-1)/2`. The `2*m*C(n,2)` term dominates for gadgetized (high-n) regions; `m^2` only
  for the high-m/low-n corner.
- **`FRAG`** (fragmentation + unmodeled-map slack) is the empirical knob; provisional **FRAG ~= 1.5**
  from the partial-peak cross-checks (ham15-high: analytic core `2*1019*61776 + 1019^2 + 4*1019*351`
  ~= 126+1+1.4 = **~128 MB**; observed TLE-peak **256 MB** => ratio ~2x incl. baseline+maps; hwb8
  analytic core `2*3517*622170` ~= **4.4 GB** vs 1.3 GB observed-at-reduced-m, consistent). **Util
  target ~60-70%** (i.e. require `predict <= 0.65*M`), refined once partial-gadget runs (the L2
  feature itself) yield completing mid-size data.
- **a~=1 (one byte per bit) holds** -- no 8x packing error. Memory IS predictable; the conservatism
  is in `FRAG`/util, deliberately, because overshoot = OOM death.

## 2026-06-23 -- INCIDENT: uncapped hwb10 gadgetize OOM'd -> WSL2 restart

Kernel log: `Out of memory: Killed process 20225 (qsyn) total-vm:71777416kB,
anon-rss:56775868kB` -- qsyn reached **56.7 GB RSS / 71.7 GB VM** on a 54 GB WSL2 VM -> Linux
OOM-killer fired -> WSL2 utility VM destabilized/restarted (all shells exit 1). The Hyper-V vPCI
**Event 33101 x4** and the **DefenderApiLogger** "max size" warning were **downstream of the
restart**, not the cause.

**Root cause (operator error):** the `gadgetize`-only probe command was run **without** the
`ulimit -v 12G` wrapper that the phasepoly runs had. On hwb10 (raw 16 qubits, no-gadget m=15840),
full gadgetization inflates `n` to the thousands; the `collapse`/tableau build for that merged
region **alone** ballooned to 56 GB. Confirmed the cap was absent (71 GB VM >> 12 GB).

**Lessons:**
1. **Operational:** wrap *every* qsyn invocation in `ulimit -v` -- including "cheap" probes. Never
   run uncapped on hwb-class circuits. (And `ulimit -v` over-budget -> `bad_alloc`->SIGABRT->WSL
   CaptureCrash, which is itself disruptive -- so prefer caps **comfortably above** a *completing*
   workload, and never ratchet-to-crash on WSL2.)
2. **Design (important):** the **gadgetize + `collapse` step is NOT memory-free** -- `collapse` on a
   huge merged region (n in the thousands) explodes on its own, before FastTODD. The earlier plan
   called the transform "cheap"; that is **false at scale**. The planner must size regions so that
   *both* the gadgetize/collapse build *and* the later FastTODD peak fit `M`. In the real feature
   this is automatic (planner selects the subset before gadgetizing), so the merged region stays
   bounded; full-gadget hwb10 = unbounded = the Vandaele MLE the feature prevents. Incident
   **validates the premise.**
3. **Model scope:** `predict_region_bytes(n,m)` (or a sibling) should also bound the collapse-build
   cost, not only the TOHPE/FastTODD peak. Investigate collapse's allocation when implementing.

## 2026-06-23 -- Calibration run 4: capped COMPLETING sweep (reproducible via scripts/calib-gadget.sh)

Every qsyn call `ulimit -v 8G` + `timeout`; CSV -> `bench2x2/calib/calib.csv`. All rc=0
(full gadgetize + phasepoly fasttodd, completing):

| circuit       |  n |   m | peak_kb | wall  |
|---------------|---:|----:|--------:|------:|
| tof_4         | 11 |  23 |   12936 | 0.00s |
| barenco_tof_4 | 14 |  28 |   12916 | 0.00s |
| tof_5         | 15 |  31 |   12928 | 0.00s |
| barenco_tof_5 | 20 |  40 |   12972 | 0.01s |
| csla_mux_3    | 21 |  62 |   13224 | 0.02s |
| rc_adder_6    | 24 |  47 |   13280 | 0.02s |
| hwb6          | 27 |  75 |   13800 | 0.05s |
| mod_red_21    | 28 |  73 |   13732 | 0.06s |
| ham15-low     | 46 |  97 |   14968 | 0.46s |
| grover_5      | 78 | 166 |   17568 | 5.64s |
| ham15-med     | 71 | 212 |   18856 | 8.04s |
| qcla_mod_7    | 84 | 237 |   22220 | 25.5s |

Fit:
- **BASELINE ~= 12.9 MB** (floor at n=11,m=23 where matrix term is negligible). Solid.
- Small-circuit variable part (peak - BASELINE) ~ 4.5-10x the analytic `2*m*C(n,2)` core but only
  single-digit MB absolute (vector-capacity rounding, transient convert/synthesis buffers). **Not
  budget-relevant** -- these fit any realistic `M`.
- Budgeting regime = large-`n` gadgetized (GBs): `2*m*C(n,2)` dominates; partial-peak cross-checks
  give **FRAG ~= 2** (ham15-high analytic ~128 MB vs observed 256 MB; hwb8 consistent -- two
  L-matrices alone are 2x).

**Production model (provisional; refine with the feature's own completing partial runs):**
`predict_region_bytes(n,m) = 12.9MB + 2.0*( 2*m*C(n,2) + m*m + 4*m*n )`, used with `safety ~= 0.7`
(require predict <= 0.7*M). Conservative by construction. Reproduce: `scripts/calib-gadget.sh`.

## 2026-06-23 -- Calibration run 5: gadgetize-STAGE peak isolated (was missing)

Run 4 measured only the combined peak. Added `/usr/bin/time -v` to the gadgetize-only call in
`scripts/calib-gadget.sh` to isolate the widen/collapse/convert-back stage (`gad_peak`) vs the full
pipeline (`full_peak`). Capped 8G, all ok:

| circuit    |  n |   m | gad_peak | full_peak | gad-base | full-base |
|------------|---:|----:|---------:|----------:|---------:|----------:|
| tof_4      | 11 |  23 |    12724 |     12796 |    ~0    |    ~0     |
| hwb6       | 27 |  75 |    13464 |     13800 |  0.7 MB  |  1.1 MB   |
| ham15-low  | 46 |  97 |    14636 |     15112 |  1.9 MB  |  2.4 MB   |
| grover_5   | 78 | 166 |    16496 |     17732 |  3.8 MB  |  5.0 MB   |
| ham15-med  | 71 | 212 |    16336 |     18764 |  3.6 MB  |  6.0 MB   |
| qcla_mod_7 | 84 | 237 |    17472 |     21964 |  4.8 MB  |  9.3 MB   |

(baseline ~12.7 MB.)

**Findings (confirms the multi-stage model; closes the "gadgetization measurement missing" gap):**
- The **gadgetize stage has real, measurable cost** -- ~half the full peak's variable part at these
  sizes (qcla_mod_7: 4.8 of 9.3 MB). It is **not** "cheap". The stage cost is mostly the
  `convert tableau->qcir` synthesis of the gadgetized circuit (QCir gate DAG over n+k qubits) plus
  widen/collapse, since the tableau itself (O(n^2) + m*n bits) is sub-MB here.
- **`full_peak = max(gad_stage, fasttodd_stage)`** -- sequential stages, same process. At small `n`
  FastTODD dominates (full > gad). **At large `n` the gadgetize stage dominates**: hwb10
  gadgetize-only OOM'd at 56 GB *before* FastTODD ran (collapse/convert is O(n^2)-ish and `n`
  explodes). So both `COLLAPSE`/`SYNTH_OUT` and `FASTTODD` terms are required; budgeting only
  FastTODD would miss the binding stage for Hadamard-dense circuits.
- Large-`n` gadgetize-stage isolation (hwb8/hwb10 gad-only peak) is **deferred** -- those run near
  the host ceiling and risk SIGABRT->CaptureCrash; the `COLLAPSE`/`SYNTH_OUT` constants are taken
  analytically (conservative) for now and refined with the feature's own completing partial runs.

## 2026-06-23 -- Large-regime reality: TIME binds before MEMORY; phasepoly peak unmeasurable at scale

- **hwb8 large-n anchor (gadgetize-only, completes):** n=1115, m=3517 -> **gad-stage peak = 1.33 GB**
  (`bench2x2/calib-large/calib.csv`). This is the only clean GB-scale measurement; it is the
  GADGETIZE stage (collapse + Clifford synthesis + convert tableau->qcir), NOT phasepoly.
- **The phasepoly/FastTODD stage peak at large m is effectively UNMEASURABLE.** Vandaele-full
  >1800 s; FastTODD is O(m^2 * width). It TLEs long before its memory peak is reached, so the
  large-regime `FASTTODD` term stays **analytic** (`2*m*C(n,2)`) -- no completing run can anchor it.
  hwb8 "full_peak 1.32 GB" was the gad stage caught at the 5 s cap, NOT phasepoly (correction: an
  earlier note implied ~4 GB phasepoly was "shown" -- it was the *formula* value, never observed).
- **TIME binds before MEMORY for Hadamard-dense circuits.** hwb8/hwb10 full-gadget TLE on compute;
  hwb10 even MLEs in the GADGETIZE stage (collapse/convert at n in the thousands) before phasepoly
  runs. So the binding constraint is often compute time, not the memory cap.
- **Design consequence (in scope to note, not to build):** selective gadgetization bounds *both*.
  The same `2*m*C(n,2)` / region-size lever that caps memory also caps phasepoly time, so a tighter
  `M` -> smaller selected regions -> faster FastTODD. L2 budgets memory only; explicit time
  budgeting is out of scope (the "how fast" / L3 axis), but `M` doubles as a rough time knob.
  The planner's per-region feasibility already covers the gadgetize stage (incl. the hwb10 MLE
  point) via the multi-stage `max(...)` model.

## 2026-06-23 -- Calibration run 6: per-STEP breakdown (scripts/calib-steps.sh, `usage` after each step)

`usage` (bare) prints period CPU time + cumulative peak RSS (ru_maxrss, monotonic) -> the add_MiB
delta = what each step added to the peak. hwb6 (n=27) vs hwb8 (n=1115), gadgetize-only
(NO_PHASEPOLY=1). NB: the running instance pre-dated the label fix, so hwb8's "phasepoly" slot is
actually `convert_out` (relabeled below).

**hwb8 (n=1115, m=3517), gadgetize-only:**

| step        | cpu_s  | peak_MiB | add_MiB |
|-------------|-------:|---------:|--------:|
| read        |   0.00 |     5.8  |   5.8   |
| convert_in  |   0.10 |     7.3  |   1.5   |
| tmerge      |   0.19 |     7.3  |   0.0   |
| hopt        |   0.06 |    13.2  |   5.9   |
| **gadgetize** | **590.95** | **1319** | **+1306** |
| convert_out |   4.17 |    1319  |  +0.0   |

(hwb6 n=27: every step add < 0.7 MiB -- baseline regime, no signal.)

**DECISIVE FINDING:** the entire gadgetize-stage cost -- **1.3 GB and 591 s CPU** -- is in the
**`gadgetize` step itself (collapse / `extract_clifford_operators` Clifford synthesis)**.
`convert tableau->qcir` (`convert_out`) is **negligible** (+0 MiB, 4 s). Therefore:
- The hwb10 56 GB MLE is **collapse**, specifically -- not the convert-back, not phasepoly.
- `collapse` is the gadgetize-stage monster in **both memory and compute**: gadgetize-only is
  already 591 s CPU at n=1115, and would be slower + >52 GB at hwb10 scale.
- **Model:** `COLLAPSE` is THE gadgetize-stage term to get right; `SYNTH_OUT`/`WIDEN`/`RESIDENT` are
  negligible by comparison. Anchor: collapse(n=1115, m=3517) ~= 1.3 GB. (Two points -- hwb6 tiny,
  hwb8 1.3 GB -- don't pin the exponent precisely; order is between n^2 and n^3 / ~ m*n^2. Refine
  with the feature's own partial runs. Conservatively bound collapse by the same `2*m*C(n,2)`-scale
  term used for FastTODD, which is >= the observed collapse here.)
- **NB collapse is also a TIME sink** (591 s) -- gadgetize-only can TLE at large n, independent of
  phasepoly. Reinforces that selective gadgetization (smaller n) is needed for tractability, not
  just memory.

## 2026-06-23 -- L2 IMPLEMENTED + verified

Shipped (6 files, +540/-37; `fastTodd.cpp` and `tableau_optimization.*` UNTOUCHED -- original
optimizer intact, confirmed by `git diff --name-only`):
- `dvlab_string.hpp`: `parse_memory_size` (K/M/G + bare bytes). `sysdep.hpp`:
  `available_memory_bytes` (Linux /proc/meminfo + cgroup; macOS sysctl+mach; Windows).
- `adaptive_gadget.{hpp,cpp}`: refactored L1's body into a shared `gadgetize_span` helper; added the
  selective overload `gadgetize_internal_hadamards(t, selected_boundaries)`, the `predict_region_bytes`
  model, the greedy `plan_boundaries` (anonymous ns; benefit = shared z-active qubits; fills to 70%
  ceiling), and `gadgetize_within_budget`.
- `tableau_cmd.cpp`: `gadgetize --memory-limit <size> | --adaptive` (mutually exclusive); no flag =
  Full (unchanged); too-small budget / invalid value -> error (never a silent no-gadget result).

Verification:
- **Unit:** all `[adaptive-gadget]` (60 assertions, 7 cases) incl. the **golden `select-all ==
  no-arg` byte-identical invariant**, subset (boundary-0 merged, boundary-1 separate), empty/OOB ->
  false, `predict_region_bytes` conservative-upper-bound, `gadgetize_within_budget` end to end.
- **Full suite:** 712,940 assertions / 52 cases PASS -> no regressions; FastTODD goldens intact.
- **End-to-end (qcla_mod_7, base n=26):** bare `gadgetize` -> full (84 q, T=163); `--memory-limit
  48M` -> **partial 82 q** (between non-gadget 26 and full 84) -- selective path proven;
  `60M/100M` -> full; `40M` -> "too small" error; `foo` -> "invalid"; `--memory-limit X --adaptive`
  -> mutual-exclusion error; `--adaptive` -> runs (available RAM).

Remaining / known limits (not blockers):
- The 32 MiB `BASELINE` makes the *partial* window narrow on tiny circuits (any merge must clear
  ~46 MiB). Fine -- those fit any real budget; refine constants with bigger completing partial runs.
- Stress bench (hwb8/ham15-high under `ulimit -v`, peak/M) is limited by **phasepoly TLE** (not OOM)
  at that scale, so a clean wall-time "win" demo isn't reachable for the largest circuits; the
  memory behavior is already characterized in runs 4-6.

## 2026-06-23 -- Adaptive end-to-end on hwb8/hwb10 (under ulimit) + the ulimit-AS finding

**CRITICAL ulimit finding.** `ulimit -v` limits *virtual* address space (AS), and qsyn reserves
huge AS in glibc/OpenBLAS per-thread arenas (RSS tiny). Bare `qcir read hwb8` under `ulimit -v 1G`
**bad_allocs** (rc=134) at ~129 MB RSS -- qsyn's *baseline* AS alone exceeds 1 G. So a `--memory-limit`
budget enforced by `ulimit -v` is unusable unless arenas are capped: run with
**`MALLOC_ARENA_MAX=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1`** so AS ~= RSS. (Document this for
users; or prefer cgroup `memory.max`, which limits RSS not AS.)

**hwb8 `gadgetize --memory-limit 1G`, `ulimit -v 1G`, arena-capped:** rc=0 (NOT killed); peak
**949 MB < 1 G**; result `1113 qubits, 462077 gates, H-gate 2217 (2 internal)`, T-family 3517.
=> **2 internal Hadamards remain (< 5)** -- meets the bar; the planner left 2 boundaries ungadgetized
to fit 1 G. Wall ~2 min. NB peak 949 MB exceeded the 0.7x1G=700 MB predicted ceiling (model
under-predicted ~1.35x) but stayed under the 1 G cap -- safe here, but the conservative margin should
account for this (lower the util target or raise FRAG).

## 2026-06-23 -- hwb10 12G FAILS: model under-predicts the collapse stage (real gap)

**hwb10 `gadgetize --memory-limit 12G`, `ulimit -v 12G`, arena-capped: rc=134 (bad_alloc at ~11.9 GB
≈ 12 G), wall ~75 s.** Same result WITH and WITHOUT the trailing `convert tableau qcir` -> it is the
`gadgetize` command itself (planner + selective `collapse`), not a convert artifact.

**Root cause:** `predict_region_bytes` was simplified to the FastTODD term `2*m*C(n,2)` on the
assumption (from hwb8, run 6) that the gadgetize/collapse stage is <= the FastTODD term. That holds
for hwb8 but **not for hwb10**: the collapse stage of the planner's selected large-`n` groups exceeds
the prediction, so the planner over-selects and actual collapse blows past 12 G. The multi-stage
model (doc sec 3.2: `peak = max(COLLAPSE, FASTTODD, ...)`) was right; collapsing it to FastTODD-only
was wrong for large circuits.

**Status of the user's stress table:**
- hwb8 `--memory-limit 1G`: PASS (not killed; peak 949 MB; 2 internal H <5).
- hwb10 `--memory-limit 12G`: **FAIL (killed)** -- predict-only guarantee violated; needs a
  collapse term.
- `--adaptive`: NOT run uncapped -- with the model under-predicting, it would over-select and risk a
  host OOM / WSL restart (the 2026-06-23 incident). Blocked until the model is fixed.

**Fix (next):** add a calibrated COLLAPSE term to `predict_region_bytes` so
`predict = BASELINE + max(FRAG*fasttodd_var, c_coll*collapse(n,m))`, calibrated from capped runs
(hwb8 collapse 1.33 GB @ n=1115 is one point; need >=1 more at larger n). Then the planner stops
selecting before collapse exceeds M. Re-run the hwb8/hwb10/--adaptive table after.

## 2026-06-23 -- Enforcement: cgroup memory.max, NOT ulimit -v (L3-critical)

`MALLOC_ARENA_MAX=1` made `ulimit -v` usable in L2 by collapsing glibc's per-thread arenas (each
reserves ~64 MB of AS) to one, so AS ~= RSS. **But it is an L2-only band-aid and BREAKS L3:** a
single malloc arena serializes every thread on one allocator lock, destroying the parallel
`phasepoly` speedup that is L3's entire purpose. Root cause: **`ulimit -v` caps virtual address
space (AS), not RSS**, and AS is inflated by arenas/BLAS buffers/mmap reservations that are never
resident.

**Decision:** the production / L3 memory enforcer is **cgroup `memory.max`** (limits RSS, arena- and
thread-count independent), not `ulimit -v`. The model predicts **RSS** (`ru_maxrss`/`peak_kb`), so
cgroup is the natural match; `ulimit -v` never was. `--memory-limit auto` already reads cgroup
(`available_memory_bytes`). `ulimit -v` + `MALLOC_ARENA_MAX=1` are kept ONLY as a single-thread
*calibration* convenience. Calibration caveat: arena=1 fragments less than multi-arena production, so
multi-arena RSS runs slightly higher -> keep the FRAG margin (and re-validate under cgroup at L3).
Doc `plan-gadgetization.md` to be updated (enforcer = cgroup; ulimit -v demoted to calibration only).

## 2026-06-23 -- ~~MEMORY MODEL IS UNSOLVED~~. Anti-overfitting rules (READ BEFORE TOUCHING THE MODEL)

(Resolved 2026-06-23, see the "MEMORY MODEL SOLVED" entry below. Kept for the anti-overfitting lessons.)

The selective-gadgetize transform + planner skeleton work and are tested. **The memory-cost model
that the planner uses is NOT trustworthy yet** -- do not believe any committed constant. I burned a
long session overfitting; these are the hard-won rules so the next session does not repeat it.

### DON'Ts (each is a mistake I actually made this session)
1. **Do NOT fit a constant to one circuit.** `collapse ~= 1072 * n^2` was fitted to hwb8 alone.
   Measured `C = collapse/n^2`: **ham15-high (n=351,m=1019) -> 463**, **hwb8 (n=1115,m=3517) -> 1072**.
   2.3x spread over two circuits => `C` is NOT constant => the 1072 constant is hwb8 overfit and
   under-predicts higher-m circuits (hwb10, m=15840) -> OOM. Any single-circuit fit will fail.
2. **Do NOT calibrate on small circuits.** tof_4/ham15-low/qcla/grover are baseline-dominated
   (variable part sub-MB on a ~13 MB floor) -> they tell you NOTHING about the asymptotic term that
   causes OOM at GB scale. Fitting on them and extrapolating is guaranteed wrong.
3. **Do NOT fixate on hwb8.** The feature must be generic across hwb8/hwb10/hwb11/etc. Tuning until
   "hwb8 passes" is meaningless.
4. **Do NOT conflate the two baselines.** RSS floor ~= 13 MB (small circuits). AS floor (what
   `ulimit -v` caps) is much larger (glibc arenas + OpenBLAS + stacks reserve virtual memory that is
   never resident). I wrongly inflated the model's RSS baseline to 96 MB to mask an AS-vs-RSS gap.
   The model predicts **RSS**; enforce with **cgroup memory.max (RSS)**, not `ulimit -v` (AS).
5. **Do NOT drop a term because it is small in one regime.** A term negligible on small/low-m
   circuits dominates on large/high-m ones. Collapse depends on BOTH n and m; an `n^2`-only or
   `n^3`-only model is wrong.
6. **Do NOT guess coefficients to match a measurement.** My code-counted allocations are ~10x BELOW
   the measured peak (see below) -> I am missing the dominant allocation. Guessing a coefficient to
   bridge that gap is fitting noise.

### What is actually known (measurements, gadgetize stage; arena-capped, /usr/bin/time -v RSS)
| circuit    |  n_gad |   m   | gadgetize peak | notes |
|------------|-------:|------:|---------------:|-------|
| ham15-high |    351 |  1019 |   ~57 MB       | completes |
| hwb8 full  |   1115 |  3517 |   ~1.33 GB (no-arena) / ~0.97 GB (arena) | completes gad-only |
| hwb10 12G  |  large | 15840 |   OOM >12 GB   | aborts in gadgetize/collapse |

- collapse peak is roughly invariant to the planner's safety/baseline knobs in my runs -> the
  planner was not actually gating selection (per-group predict let many small groups all pass while
  the true cost tracks the GLOBAL n_total). Unresolved: per-group vs global n_total.
- **Code-count mismatch (the key blocker):** `AGSynthesisStrategy::synthesize`
  (`stabilizer_tableau.cpp:230`) builds `clifford_ops` = O(n^2) `CliffordOperator` (~24 B each) ->
  ~120-240 MB at n=1115. Resident rotations `m*n` bits, tableaux `n^2` bits = MB. **Sum is ~10x below
  the measured ~1 GB.** So a dominant allocation is UNIDENTIFIED by reading the code.
- **Likely confound:** the calib "gadgetize-only" commands INCLUDED `convert tableau qcir`; the
  output QCir DAG (hwb8: 462k gates) may be the real GB sink, not collapse. Must measure
  gadgetize-only (quit, NO convert) vs +convert separately, under a NON-capping ulimit.

### The ONLY correct next step: PROFILE, don't fit
Run hwb8 gadgetize (and gadgetize+convert separately) under `valgrind --tool=massif` (or RSS-sample
per step), identify the actual dominant allocation and its scaling in (n,m) from the profile, THEN
derive the term from its code dimensions (generic, no per-circuit fit). Model predicts RSS; only the
~13 MB floor + one fragmentation factor are empirical. Enforce via cgroup, not ulimit -v.
`MALLOC_ARENA_MAX=1` is an L2-single-thread calibration crutch -- it serializes malloc and BREAKS L3.

### Current code state (to fix next session)
`predict_region_bytes` in `adaptive_gadget.cpp` currently has `baseline=16MB`, `1072*n^2 + m*n`
(the hwb8 overfit) with planner per-group n and safety 0.9. **Treat all those numbers as wrong
placeholders.** The selective transform, planner structure, `gadgetize --memory-limit/--adaptive`
wiring, and tests are sound; only the cost model needs the profile-derived rewrite.

### GATE: PASSED (with the conservative-util caveat)
Tight *predict-only* budgeting is feasible; "~99%" is **not** safe (real peak 2-4x the naive term,
TLE blocks exact fitting), so production targets **~60-70% util** with the analytic model above +
`FRAG`. The model + `FRAG` are refined later with the feature's own completing partial runs. No
`fastTodd.cpp` edit. -> proceed to implementation (parsing -> memory model fn -> selective transform +
planner -> command wiring).

## 2026-06-23 -- MEMORY MODEL SOLVED (profile-derived, no magic) + L2 safety verified

Profiled the gadgetize sink with in-code `/proc/self/statm` RSS probes (valgrind absent; hwb8 gadgetize
is 13 min so massif was impractical). **The 40x "unidentified allocation" was not exotic** -- it is the
widened tableau held in multiple coexisting copies:

- `extended_with_ancillae(k)` materializes the WHOLE input tableau at the global width n (= base_n +
  total selected ancillae): `s_clifford` StabilizerTableaux, each `2n*(2n+1)/8` bits PLUS `2n *
  sizeof(PauliProduct=32B)` struct overhead (the overhead dominates at small n -- the earlier
  small-circuit under-prediction). hwb8: `ext` jumped 24 MB -> 702 MB at n=1115, s_clifford~=1145.
- the selective path then held 2-3 more full-width copies: the per-group `interior`, the `middle`
  accumulator, and `out` (which **copied** `middle`). Total hwb8 gadgetize peak ~1.42 GB; hwb10 ~1.15 GB.

**Two fixes (both reduce memory AND simplify the model):**
1. Transform: in `gadgetize_span`, `interior.push_back(std::move(ext[i]))` (consume ext piecewise -- the
   group's subtableaux are not read elsewhere), and in the selective assembly `out.push_back(std::move(st))`
   (don't copy `middle` into `out`). Result is byte-identical (golden test still passes); peak ~halved.
2. `predict_peak_bytes(n, s_clifford, m_total, m_region)` = max over stages of EXACT sizeof/paper terms:
   - GADGETIZE = `kWidenCopies(=3) * (s_clifford*stab_bytes + m_total*rot_bytes) + 24*(2n)^2`, x `kFragPercent(=150)`.
     `stab_bytes = 2n*(2n+1)/8 + 2n*32`; `rot_bytes = (2n+1)/8 + 48`.
   - PHASEPOLY = `2*m_region*(n + C(n,2)) + m_region^2 + n*m_region` (FastTODD L + Lt + augmented + A; 1
     byte/bit BooleanMatrix::Row; cols `n+C(n,2)` per Heyfron2018/Vandaele2024).
   Every constant is a `sizeof` or a paper dimension. The ONLY non-sizeof factor is `kFragPercent=150`
   (glibc fragmentation, FLAGGED empirical in the code). The floor is **live RSS** read at plan time
   (`current_rss_bytes`), NOT a constant -- the resident input (QCir DAG) ranges 13 MB (tof3) -> 140 MB
   (hwb11), so a fixed floor under-counts big inputs (this was the residual hwb11 miss).
3. Planner: `plan_boundaries` predicts at the **global** width `base_n + total_selected_ancillae` with
   `m = largest region` (NOT the old per-group width, which is what OOM'd hwb10). `s_clifford = tableau.size()-R`.

**Validation (cgroup `memory.max`, the production enforcer; `ulimit -v` retired):** predict >= actual on
every measured point.

| circuit    | budget | predict | actual | within budget |
|------------|-------:|--------:|-------:|:-------------:|
| ham15-high |   30M  |  29 MiB | 24 MiB | yes           |
| ham15-high |   60M  |  59     | 38     | yes           |
| hwb10      |  512M  | 511     | 456    | yes           |
| hwb11      |  512M  | 509     | 433    | yes           |
| hwb10      |   1G   | 1023    | 797    | yes (cgroup rc=0) |
| hwb11      |   1G   | 1023    | 873    | yes           |
| **hwb10**  | **12G**| **12269** | **7714** | **rc=0, NO OOM** |

**Headline regression fixed:** hwb10 `--memory-limit 12G` previously OOM'd (rc=134 at ~11.9 GB); now
completes at 7.5 GB < 12 GB cap (predict 12.0 GB >= actual 7.5 GB), wall ~59 s. Monotonicity
`peak(full) > peak(adaptive) > peak(no-gadget)` holds wherever above the ~13 MB RSS noise floor
(ham15-high 17 -> 23 -> 72 MB; hwb10/hwb11 even more) and is strict in the prediction always. Tiny
circuits (tof3/hwb6/ham15-low/med) are floor-clamped (gadget growth < allocator noise) and fit any real
budget. Full suite 438,639 assertions / 52 cases PASS; FastTODD/`tableau_optimization`/`stabilizer_tableau`
UNTOUCHED; RSS probes reverted.

## 2026-06-23 -- model is SAFE but LOOSE at large budgets (FIXME); reproducible tooling; T-table

Post-"solved" review exposed that the cost model, while safe (never under-predicts -> never OOMs), is
conservative and gets looser as n grows. Documented as a `FIXME` at the `kWidenCopies`/`kFragPercent`
constants in `adaptive_gadget.cpp`.

- **Utilization drops with budget** (hwb10, `--memory-limit` vs measured peak): 512M->89%, 1G->78%,
  12G->63%, 48G->51% (predict/actual = 1.12, 1.28, 1.59, 1.95). So at a big budget on a big circuit the
  planner gadgetizes FEWER Hadamards than would fit, leaving T on the table. Root cause: the gadgetize
  multiplier is a FIXED ~4.5x envelope (`kWidenCopies=3` x `kFragPercent=150`), but the true peak/
  structural ratio is not constant -- it falls as n grows (the `output` block count shrinks as
  gadgetization merges regions; `ext` is `std::move`-consumed). hwb8 FULL path (one group, no `middle`
  accumulation) sits at the structural floor ~2.08x of `s_clifford*stab(n)` -> floor copies = 2.
- **The two multipliers are NOT paper/sizeof.** `kWidenCopies` is code-structural (countable copies);
  `kFragPercent` is glibc-allocator-specific (meaningless under a different allocator). A tight fix needs
  an n-dependent multiplier; the exact form is undetermined -- log/power fits over the few feasible points
  were under-constrained, and large budgets (>~50 GB actual) are UNMEASURABLE on this 57 GB host, so the
  asymptotic floor where util levels off cannot be observed here. Any clean form must be a positive,
  bounded decay (a `const - b*ln n` line is unphysical -> goes negative). Left as FIXME, not "fitted".
- **phasepoly time ~5x per 2x budget** (hwb8 full chain, T / wall): no-gad 3510/0.6s; 64M 3487/2.9s;
  128M 3471/15s; 256M 3425/76s; 512M 3372/355s. T falls monotonically; 1G ~= 5x512M ~= 1775s sits right
  on the 1800 s TLE boundary; >=2G is hopeless. So on hwb8 the completing adaptive curve tops out ~1G,
  and `--adaptive`(auto) == full (detects ~53 GB, gadgetizes fully) -> TLE in phasepoly.
- **Reproducible tooling (committed):** `scripts/gadget-budget-check.sh` (cgroup-capped predict-vs-actual
  + util/no-OOM, per (circuit,budget)) and `scripts/gadget-tcount-table.sh` (single circuit x chosen
  modes -> T / wall / peak; SEQUENTIAL so timing is reproducible; modes = none|<size>|auto|full). Both are
  portable: repo-relative, in-repo circuits, mktemp work dir, no absolute paths. Bug found+fixed: must
  nest `/usr/bin/time -v timeout ... qsyn` (not the reverse), else a TLE kills `time` before it prints and
  the peak is lost.
- **Inspectability:** added a debug-level diagnostic in `plan_boundaries` ("gadgetize budget X MiB:
  selected N of M internal Hadamards (n A -> B, s_clifford, m_region); predicted peak ...") -- silent at the
  default user level (no jargon leak), shown with `logger debug`. The floor is the LIVE plan-time RSS
  (`current_rss_bytes`, /proc/self/statm), not a constant (input QCir DAG resident ranges 13->140 MB).

## 2026-06-24 -- allocator knob + cross-allocator sweep + hwb8 full-chain T-table

- **Allocator factor is now user-defined (decision A).** `kFragPercent` -> env `QSYN_GADGET_HEADROOM`
  (a percent; default 150 measured on glibc), read once. Documented as THE one allocator-dependent,
  non-sizeof knob; re-validate per allocator (jemalloc/tcmalloc/mimalloc/musl). Verified live: 75/150/300
  -> 173/105/58 Hadamards selected on ham15-high @ 40M.
- **hwb8 full optimization chain (read -> ... -> gadgetize[mode] -> phasepoly fasttodd), HEADROOM=150, glibc:**

  | mode          | final T | wall   | peak RSS |
  |---------------|--------:|-------:|---------:|
  | no-gadgetize  |    3510 |   0.6s |   60 MiB |
  | adaptive 64M  |    3487 |   2.9s |   68 MiB |
  | adaptive 128M |    3471 |  15.2s |  107 MiB |
  | adaptive 256M |    3425 |  75.8s |  172 MiB |
  | adaptive 512M |    3372 | 355s   |  288 MiB |
  | adaptive 1G   |     TLE |    TLE |        - |
  | adaptive auto |     TLE |    TLE |        - |
  | full          |     TLE |    TLE |        - |

  T falls monotonically with budget; cost rises. **phasepoly wall ~5x per 2x budget** -> 1G ~1775s sits on
  the 1800s cap (TLE is the 5x law, NOT parallel contention -- the law held cleanly across contention
  levels). Completing curve tops out ~1G; `auto`==`full` (detects ~53GB). TLE peaks not captured (the
  time/timeout nesting bug, since fixed in the script).

- **Cross-allocator (HEADROOM=100 = pure structural predict; LD_PRELOAD; /usr/bin/time peak only):**

  | hwb8 256M      | glibc | arena=1 | mimalloc2 | mimalloc3 | jemalloc |
  |----------------|------:|--------:|----------:|----------:|---------:|
  | wall (m:ss)    |  3:15 |    3:17 |      3:21 |      3:20 |     3:30 |
  | T-count        |  3409 |    3409 |      3409 |      3409 |     3408 |
  | peak RSS (MiB) |   231 |     230 |       236 |       236 |      244 |

  | hwb8 512M      | glibc | arena=1 | mimalloc2 | mimalloc3 | jemalloc |
  |----------------|------:|--------:|----------:|----------:|---------:|
  | wall (m:ss)    | 16:47 |   16:50 |     16:14 |     17:15 |    17:15 |
  | T-count        |  3346 |    3346 |      3346 |      3346 |     3340 |
  | peak RSS (MiB) |   397 |     397 |       444 |       467 |      445 |

  Findings: **glibc is the LOWEST-RSS allocator** (mimalloc3 worst, +18% @512M); spread GROWS with scale
  (~6% @256M -> ~18% @512M); `arena=1` == glibc (bloat is NOT multi-arena frag -- likely the phasepoly
  BooleanMatrix); T-count allocator-invariant. **Implication: the glibc-calibrated default 150 is the
  OPTIMISTIC case** -- under mimalloc3 a glibc prediction under-covers ~18% -> OOM risk (open: keep 150 +
  doc "raise for non-glibc", or bump default to the worst allocator). Methodology: never read live `ps`
  RSS -- it is mid-run, always <= the `ru_maxrss` peak; only /usr/bin/time peak counts (the 728M `auto` and
  "mimalloc uses half" were both live-transient misreads).
