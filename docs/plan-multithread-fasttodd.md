# Multithreading fastTodd internally -- direct-edit plan

Inner parallelism: edit `src/tableau/optimize/fastTodd.cpp` to parallelize its own algorithm.
Distinct from L3. L3 threads *across* disjoint regions and leaves fastTodd untouched; this plan threads
*inside* a single fastTodd call.

## Rules
- No ASCII-code.
- No leak on internal processes.
- Place worklog, be brief.
- Do not pile trash up.

## Scope
- Direct edits to `fastTodd.cpp` are allowed here. ("fastTodd immutable" holds only under L3 -- there it
  must be threaded around, not edited.)
- Target: fastTodd's internal hot loops (per worklog/profiling, not assumed here).
- Complementary to L3: L3 = outer region parallelism; this = inner per-call parallelism.

## Relationship to L3
- L3 and this can compose: outer k_regions x inner k_threads. Nested threading needs a thread/core
  budget so the product does not oversubscribe cores or compound memory. Resolve before combining.

## Safety
- Editing internals can introduce shared mutable state across inner threads -- the exact thing L3 relied
  on being absent. Inner threads must keep per-thread data disjoint or synchronized.
- Determinism: parallel reductions must stay byte-identical to single-threaded, or document any
  tolerated nondeterminism.

## Status
Not started. Needs fastTodd internal hotspot profiling before design.
