# Plan: Opt-in SIMD via build flags (arm64 / x86-64)

## Goal

Let users build qsyn with architecture-targeted compiler flags so the compiler
auto-vectorizes the existing hot loops — primarily the GF(2) row-XOR in the
FastTODD / TOHPE / TODD phase-polynomial optimization path. No source changes,
no intrinsics, no data-layout changes.

## Decisions (settled)

- **Scope:** build flags only. The compiler does the vectorization; we just tell
  it which ISA to target. (Bit-packing rows into 64-bit words + hand-written
  AVX2/NEON kernels would yield more, but is explicitly *not* this project.)
- **Dispatch:** build-time `-march`. No runtime CPU detection / kernel
  multi-versioning.
- **Default:** `off` — today's behavior. SIMD is strictly opt-in, so CI, Docker,
  and any distributed binaries are unaffected unless they ask for it.
- **Levels:** `off` / `generic` / `modern` / `native` (intent-based, arch-neutral
  names — see table).

## Background (why this is worth doing)

- The build currently sets `-O3` with **no `-march`/`-mtune`**
  (`CMakeLists.txt:13`), so codegen targets the baseline x86-64 / armv8-a ISA and
  the compiler cannot emit AVX2/FMA or newer NEON-class instructions.
- The hot path is `BooleanMatrix::Row::operator+=` in
  `src/util/boolean_matrix.cpp` — a `for` loop XOR-ing two
  `std::vector<unsigned char>` rows. It is already instrumented with the
  `QSYN_TOHPE_PROFILE` op/byte counters, confirming it as the inner loop of the
  TODD/FastTODD/TOHPE work (Vandaele 2024).
- The Vandaele Rust reference is already built with
  `RUSTFLAGS="-C target-cpu=native"` (`justfile:113`). So for an apples-to-apples
  comparison, qsyn should be buildable with an equivalent `native` tuning.

> **Note on ceiling.** Rows store one GF(2) bit per `unsigned char` (one byte per
> bit). Auto-vectorization still helps (wider byte-wise XOR per instruction), but
> the larger win — packing 64 bits per word — is deliberately out of scope here.
> Capture before/after numbers so we know whether the follow-up bit-packing
> project is worth it.

## Level → flag mapping

| Level             | Intent                   | x86-64 flags                  | arm64 (aarch64) flags | Portable?             |
|-------------------|--------------------------|-------------------------------|-----------------------|-----------------------|
| `off` *(default)* | unchanged                | *(none)*                      | *(none)*              | yes — baseline ISA    |
| `generic`         | portable but vectorizing | `-march=x86-64-v2`            | `-march=armv8-a`      | yes — any modern CPU  |
| `modern`          | assume AVX2/FMA-class HW | `-march=x86-64-v3`            | `-march=armv8.2-a`    | most CPUs since ~2015 |
| `native`          | this exact build host    | `-march=native -mtune=native` | `-mcpu=native`        | **no** — host only    |

Notes:
- `modern` maps to AVX2+FMA on x86 (`x86-64-v3`) and to `armv8.2-a` on arm64.
  These are analogous "one generation past baseline" targets, not identical
  feature sets — documented in the option help so arm64 users aren't misled.
- On Apple Silicon, `-mcpu=native` / `-march=native` interacts with the existing
  GCC-on-arm64 arch-clearing logic in `CMakeLists.txt` (~lines 23–27); the
  flag-validation step below makes any unsupported flag a no-op rather than a
  build break.

## Implementation

All changes are in the build system. No `src/` edits.

### Step 1 — Add the `QSYN_SIMD` cache option (CMakeLists.txt)

Near the top of `CMakeLists.txt`, after `CheckCXXCompilerFlag` is included
(it already is — `include(CheckCXXCompilerFlag)`):

```cmake
set(QSYN_SIMD "off" CACHE STRING
    "SIMD/arch tuning level: off | generic | modern | native")
set_property(CACHE QSYN_SIMD PROPERTY STRINGS off generic modern native)
```

### Step 2 — Resolve level + arch into a candidate flag list

Define a function that maps `(QSYN_SIMD, CMAKE_SYSTEM_PROCESSOR)` to candidate
flags, validates each with `check_cxx_compiler_flag`, and collects the supported
ones into a cache/internal variable `QSYN_SIMD_FLAGS`.

```cmake
function(qsyn_resolve_simd_flags out_var)
  set(_flags "")
  if(QSYN_SIMD STREQUAL "off")
    set(${out_var} "" PARENT_SCOPE)
    return()
  endif()

  # arch family
  string(TOLOWER "${CMAKE_SYSTEM_PROCESSOR}" _arch)
  if(_arch MATCHES "x86_64|amd64")
    set(_is_x86 TRUE)
  elseif(_arch MATCHES "aarch64|arm64")
    set(_is_arm TRUE)
  else()
    message(WARNING "QSYN_SIMD=${QSYN_SIMD}: unrecognized arch '${_arch}', ignoring.")
    set(${out_var} "" PARENT_SCOPE)
    return()
  endif()

  set(_candidates "")
  if(_is_x86)
    if(QSYN_SIMD STREQUAL "generic")
      list(APPEND _candidates "-march=x86-64-v2")
    elseif(QSYN_SIMD STREQUAL "modern")
      list(APPEND _candidates "-march=x86-64-v3")
    elseif(QSYN_SIMD STREQUAL "native")
      list(APPEND _candidates "-march=native" "-mtune=native")
    endif()
  elseif(_is_arm)
    if(QSYN_SIMD STREQUAL "generic")
      list(APPEND _candidates "-march=armv8-a")
    elseif(QSYN_SIMD STREQUAL "modern")
      list(APPEND _candidates "-march=armv8.2-a")
    elseif(QSYN_SIMD STREQUAL "native")
      list(APPEND _candidates "-mcpu=native")
    endif()
  endif()

  # validate each flag; drop unsupported ones with a warning
  foreach(_f IN LISTS _candidates)
    string(MAKE_C_IDENTIFIER "QSYN_SIMD_SUPPORTS_${_f}" _cache_key)
    check_cxx_compiler_flag("${_f}" ${_cache_key})
    if(${_cache_key})
      list(APPEND _flags "${_f}")
    else()
      message(WARNING "QSYN_SIMD: compiler rejects '${_f}', skipping it.")
    endif()
  endforeach()

  set(${out_var} "${_flags}" PARENT_SCOPE)
endfunction()

qsyn_resolve_simd_flags(QSYN_SIMD_FLAGS)
message(STATUS "QSYN_SIMD=${QSYN_SIMD} -> flags: ${QSYN_SIMD_FLAGS}")
```

### Step 3 — Apply flags to the relevant targets

Add the resolved flags as `PRIVATE` compile options on:

- `${QSYN_LIB_NAME}` — **most important**; `boolean_matrix.cpp` and the tableau
  optimization sources compile here.
- `${CMAKE_PROJECT_NAME}` (the `qsyn` executable).
- `${UNIT_TEST_NAME}` — so tests exercise the same codegen.

```cmake
if(QSYN_SIMD_FLAGS)
  target_compile_options(${QSYN_LIB_NAME}     PRIVATE ${QSYN_SIMD_FLAGS})
  target_compile_options(${CMAKE_PROJECT_NAME} PRIVATE ${QSYN_SIMD_FLAGS})
  target_compile_options(${UNIT_TEST_NAME}    PRIVATE ${QSYN_SIMD_FLAGS})
endif()
```

(Place each call after its target is defined.) These are compile *options*, not
warning flags, so they don't interact with the existing `-Werror`.

### Step 4 — Thread the option through the Makefile

In the `configure` / `configure-debug` targets in `Makefile`, forward an optional
`QSYN_SIMD` variable so `make release QSYN_SIMD=native` works:

```make
QSYN_SIMD ?= off
# ...in the cmake invocation, add:
	-DQSYN_SIMD=$(QSYN_SIMD) \
```

### Step 5 — Document usage

- Add a short section to `README.md` (build options): the four levels, the
  default-`off` rationale, and the portability warning for `native`/`modern`.
- For Docker/release images, keep `off` (or `generic` if a vectorizing-but-
  portable floor is wanted) — never `native`, since the build host ≠ run host.

## Verification

1. **Flags actually apply.** Configure with `-DQSYN_SIMD=native`, then inspect
   `compile_commands.json` for `boolean_matrix.cpp` — confirm the `-march` flag is
   present.
2. **Vectorization landed.** Compile `boolean_matrix.cpp` with `-fopt-info-vec`
   (GCC) / `-Rpass=loop-vectorize` (Clang) and confirm the `operator+=` loop
   vectorizes; or `objdump -d` the object and look for vector mnemonics
   (`vpxor`/`ymm` on x86, NEON `eor v…` on arm64) absent in the `off` build.
3. **Correctness unchanged.** `make test` passes identically for `off`,
   `generic`, `modern`, and (on the build host) `native`.
4. **Speedup measured.** Run `scripts/bench-fasttodd.sh` for `off` vs `native`
   and record T-count-identical results with wall-time deltas. Optionally use
   `QSYN_TOHPE_PROFILE=1` to confirm the row-add op count is unchanged (only
   per-op cost should drop).

## Risks & mitigations

| Risk                                                           | Mitigation                                                                                  |
|----------------------------------------------------------------|---------------------------------------------------------------------------------------------|
| `native`/`modern` binary `SIGILL`s on an older CPU             | Default is `off`; portability warning documented; distribution builds stay `off`/`generic`. |
| Compiler doesn't know a flag (old GCC/Clang, Apple GCC quirks) | Every flag gated by `check_cxx_compiler_flag`; unsupported → warn + skip, never fail.       |
| Unknown `CMAKE_SYSTEM_PROCESSOR` (e.g. cross-compile)          | Function warns and emits no flags.                                                          |
| Small/no speedup because of byte-per-bit layout                | Expected; benchmark numbers decide whether the follow-up bit-packing project is justified.  |

## Out of scope (explicit non-goals)

- Bit-packing GF(2) rows into 64-bit words.
- Hand-written AVX2 / AVX-512 / NEON intrinsics or kernels.
- Runtime CPU-feature detection / function multi-versioning.
- `-ffast-math` and other FP-semantics-altering flags. The hot path is integer
  GF(2) work and gains nothing; meanwhile fast-math could perturb the
  tensor/extractor floating-point paths. Deliberately excluded.

## Follow-up (if benchmarks justify it)

If Step-4 numbers show the byte-per-bit layout is the bottleneck, open a separate
project for bit-packed rows (`uint64_t` words) — at which point auto-vectorization
of the word loop, and optionally explicit SIMD, become the high-value next step.
