# TODD benchmark automation. See docs/benchmark-comparison.md for results.

run_dir := "qsyn-bench"
circuits_dir := run_dir / "circuits"
paper_csv := "benchmark/paper-todd.csv"
python := `command -v uv >/dev/null 2>&1 && echo "uv run" || echo "python3"`

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
