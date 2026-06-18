# TODD benchmark automation. See docs/benchmark-comparison.md for results.

run_dir := "qsyn-bench"
circuits_dir := run_dir / "circuits"
paper_csv := "benchmark/paper-todd.csv"
python := `command -v uv >/dev/null 2>&1 && echo "uv run" || echo "python3"`
vandaele_dir := "../quantum-circuit-optimization"
verify_circuits := "tof_3 barenco_tof_3 gf2^4_mult gf2^5_mult"

# List available recipes.
default:
    @just --list

# Build qsyn (release).
build:
    make release

# Sanitize paper benchmark circuits into {{circuits_dir}} (T* -> tdg, S* -> sdg) and smoke-parse each.
prepare:
    #!/usr/bin/env bash
    set -euo pipefail
    mkdir -p "{{circuits_dir}}"
    skipped=()
    while IFS=, read -r circuit _group _tier _n _label _hp _ht _htt _vp _vt _vtt _vf file _notes; do
        [ "$circuit" = "circuit" ] && continue
        if [ -z "$file" ]; then
            skipped+=("$circuit (no circuit file)")
            continue
        fi
        out="{{circuits_dir}}/$circuit.qc"
        sed -e 's/\bT\*/tdg/g' -e 's/\bS\*/sdg/g' -e 's/\bZ\*/z/g' -e 's/^cnot /tof /' "$file" > "$out"
        if ! ./qsyn -q --no-version -c "qcir read $out; quit -f" 2>&1 | rg -q "error|Error"; then
            echo "ok    $circuit  <- $file"
        else
            skipped+=("$circuit (qsyn cannot parse $file)")
            rm -f "$out"
        fi
    done < "{{paper_csv}}"
    if [ "${#skipped[@]}" -gt 0 ]; then
        printf '%s\n' "${skipped[@]}" > "{{run_dir}}/skipped.txt"
        echo "--- skipped (see {{run_dir}}/skipped.txt):"
        printf '  %s\n' "${skipped[@]}"
    fi

# Run the TODD pipeline on one circuit; appends to {{run_dir}}/results.csv.
bench circuit timeout="1800":
    ./scripts/todd-benchmark.sh {{circuit}} {{timeout}}

# Run all fast-tier circuits (small/medium; minutes overall).
bench-all timeout="1800":
    #!/usr/bin/env bash
    set -euo pipefail
    while IFS=, read -r circuit _group tier rest; do
        [ "$tier" = "fast" ] || continue
        ./scripts/todd-benchmark.sh "$circuit" "{{timeout}}" || true
    done < "{{paper_csv}}"

# Run the slow-tier circuits (large T-counts; may take hours, DNF rows on timeout).
bench-all-slow timeout="7200":
    #!/usr/bin/env bash
    set -euo pipefail
    while IFS=, read -r circuit _group tier rest; do
        [ "$tier" = "slow" ] || continue
        ./scripts/todd-benchmark.sh "$circuit" "{{timeout}}" || true
    done < "{{paper_csv}}"

# Generate {{run_dir}}/benchmark-comparison.md; review, then copy to docs/.
report:
    {{python}} scripts/todd-report.py

# Remove generated benchmark outputs (keeps paper-todd.csv).
clean-bench:
    rm -rf "{{circuits_dir}}" "{{run_dir}}/logs" "{{run_dir}}/results.csv" "{{run_dir}}/skipped.txt"

# Phase-2: equivalence + T-count for todd/tohpe/fasttodd on small circuits -> verify-small.csv
verify-small:
    #!/usr/bin/env bash
    set -euo pipefail
    mkdir -p "{{run_dir}}"
    out="{{run_dir}}/verify-small.csv"
    echo "circuit,strategy,tcount,equiv" > "$out"
    for c in {{verify_circuits}}; do
        qc="{{circuits_dir}}/$c.qc"
        for s in todd tohpe fasttodd; do
            log=$(./build/qsyn -q --no-version -c "qcir read $qc; convert qcir tableau; tableau optimize tmerge; tableau optimize hopt; tableau optimize phasepoly $s; convert tableau qcir; qcir print --stat; qcir equiv 0 1; quit -f" 2>&1)
            t=$(echo "$log" | rg -oP 'T-family\s*:\s*\K\d+' | tail -1)
            eq=$(echo "$log" | rg -qi "are equivalent" && echo pass || echo CHECK)
            echo "$c,$s,${t:-NA},$eq" | tee -a "$out"
        done
    done

# Phase-2: prove the FastTODD win — gf2 circuits finish fast (not TLE) -> verify-fast.csv
verify-fast timeout="600":
    #!/usr/bin/env bash
    set -euo pipefail
    mkdir -p "{{run_dir}}"
    out="{{run_dir}}/verify-fast.csv"
    echo "circuit,tcount,equiv,seconds,status" > "$out"
    for c in gf2^9_mult gf2^16_mult; do
        qc="{{circuits_dir}}/$c.qc"
        start=$(date +%s)
        if log=$(timeout -k 10 {{timeout}} ./build/qsyn -q --no-version -c "qcir read $qc; convert qcir tableau; tableau optimize tmerge; tableau optimize hopt; tableau optimize phasepoly fasttodd; convert tableau qcir; qcir print --stat; qcir equiv 0 1; quit -f" 2>&1); then st=ok; else st=TLE; fi
        sec=$(( $(date +%s) - start ))
        t=$(echo "$log" | rg -oP 'T-family\s*:\s*\K\d+' | tail -1)
        eq=$(echo "$log" | rg -qi "are equivalent" && echo pass || echo CHECK)
        echo "$c,${t:-NA},$eq,$sec,$st" | tee -a "$out"
    done

# Phase-2: cross-check one circuit's fasttodd T-count against the Vandaele reference binary.
verify-vandaele circuit:
    #!/usr/bin/env bash
    set -euo pipefail
    RUSTFLAGS="-C target-cpu=native" cargo build -r --manifest-path "{{vandaele_dir}}/Cargo.toml"
    echo "--- Vandaele FastTODD ({{circuit}}):"
    "{{vandaele_dir}}/target/release/quantum_circuit_optimization" FastTODD "{{vandaele_dir}}/circuits/inputs/{{circuit}}.qc" | rg "T-count"

# Phase-2 gate: build, small-circuit equiv/T, gf2 fast-win, unit tests.
verify: build verify-small verify-fast
    make test
