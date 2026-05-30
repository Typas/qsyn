/**
 * @file
 * @brief implementation of the fastTodd phase polynomial optimization
 * @copyright Copyright(c) 2024 DVLab, GIEE, NTU, Taiwan
 */

#include <spdlog/spdlog.h>

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdlib>
#include <cstdint>
#include <iterator>
#include <optional>
#include <ranges>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include "../tableau_optimization.hpp"
#include "./todd.hpp"
#include "tableau/pauli_rotation.hpp"
#include "tableau/stabilizer_tableau.hpp"
#include "util/boolean_matrix.hpp"
#include "util/ordered_hashmap.hpp"
#include "util/phase.hpp"
#include "util/util.hpp"

#include "fmt/core.h"

extern bool stop_requested();

namespace qsyn::tableau {
using namespace todd;
using namespace signature;
using namespace qsyn::tableau;

namespace {
using Clock = std::chrono::steady_clock;

bool tohpe_trace_enabled() {
    static bool const enabled = [] {
        auto const* raw = std::getenv("QSYN_TOHPE_TRACE");
        return raw != nullptr && raw[0] == '1';
    }();
    return enabled;
}

bool tohpe_trace_compare_enabled() {
    static bool const enabled = [] {
        auto const* raw = std::getenv("QSYN_TOHPE_TRACE_COMPARE");
        return raw != nullptr && raw[0] == '1';
    }();
    return enabled;
}

bool tohpe_trace_stop_on_divergence() {
    static bool const enabled = [] {
        auto const* raw = std::getenv("QSYN_TOHPE_TRACE_STOP_ON_DIVERGENCE");
        return raw != nullptr && raw[0] == '1';
    }();
    return enabled;
}

bool tohpe_use_fast_port() {
    static bool const enabled = [] {
        auto const* raw = std::getenv("QSYN_TOHPE_USE_FAST_PORT");
        if (raw == nullptr) return true;
        return raw[0] != '0';
    }();
    return enabled;
}

std::string row_to_bitstring(dvlab::BooleanMatrix::Row const& row) {
    std::string bits;
    bits.reserve(row.size());
    for (auto const bit : row) {
        bits.push_back(bit ? '1' : '0');
    }
    return bits;
}

std::string bytes_to_bitstring(std::vector<unsigned char> const& bits) {
    std::string ret;
    ret.reserve(bits.size());
    for (auto const bit : bits) {
        ret.push_back(bit ? '1' : '0');
    }
    return ret;
}

std::vector<std::string> polynomial_term_bitstrings(Polynomial const& polynomial) {
    if (polynomial.empty()) {
        return {};
    }
    auto const table = dvlab::transpose(load_phase_poly_matrix(polynomial));
    auto ret         = std::vector<std::string>{};
    ret.reserve(table.num_rows());
    for (auto const& row : table) {
        ret.emplace_back(row_to_bitstring(row));
    }
    std::sort(ret.begin(), ret.end());
    return ret;
}

std::string summarize_term_bitstrings(std::vector<std::string> const& terms) {
    return fmt::format("[{}]", fmt::join(terms, ","));
}

std::string summarize_indices(std::vector<size_t> const& indices) {
    return fmt::format("[{}]", fmt::join(indices, ","));
}

void tohpe_trace(std::string_view impl, std::string_view stage, std::string const& message) {
    if (!tohpe_trace_enabled()) return;
    fmt::println("[TOHPE_TRACE] impl={} stage={} {}", impl, stage, message);
}

struct TohpeProfiler {
    bool const enabled = [] {
        auto const* raw = std::getenv("QSYN_TOHPE_PROFILE");
        return raw != nullptr && raw[0] == '1';
    }();

    uint64_t tohpe_once_calls = 0;

    uint64_t kernel_rows_processed = 0;
    uint64_t kernel_entry_pivots   = 0;
    uint64_t kernel_row_adds       = 0;
    uint64_t clear_column_calls    = 0;
    uint64_t clear_column_row_adds = 0;
    uint64_t score_seed_inserts    = 0;
    uint64_t score_pair_count      = 0;
    uint64_t score_map_updates     = 0;
    uint64_t removed_rows          = 0;

    /// `get_row_products` inside post-properization rebuild only (see `tohpe_once`).
    uint64_t ns_rebuild_row_products    = 0;
    /// `get_l_matrix` + transpose inside post-properization rebuild only.
    uint64_t ns_rebuild_get_l_transpose = 0;

    ~TohpeProfiler() {
        if (!enabled) return;
        if (tohpe_once_calls == 0) return;
        auto const to_ms = [](uint64_t ns) -> double { return static_cast<double>(ns) / 1'000'000.0; };
        auto const sum_ns = ns_rebuild_row_products + ns_rebuild_get_l_transpose;
        auto const sum_ms = to_ms(sum_ns);
        auto const pct_of_pair = [sum_ms, &to_ms](uint64_t ns) {
            if (sum_ms == 0.0) return 0.0;
            return (to_ms(ns) * 100.0) / sum_ms;
        };
        // Use fmt::println so output is visible even when spdlog level is "off" (see examples/t-opt.dof).
        fmt::println("[TOHPE_PROFILE] row_products_ms={:.3f} ({:.1f}%) l_matrix_ms={:.3f} ({:.1f}%)",
                     to_ms(ns_rebuild_row_products), pct_of_pair(ns_rebuild_row_products),
                     to_ms(ns_rebuild_get_l_transpose), pct_of_pair(ns_rebuild_get_l_transpose));
        fmt::println("[TOHPE_PROFILE] row_add_ops={} row_add_bytes={}",
                     dvlab::profiled_row_add_ops(),
                     dvlab::profiled_row_add_bytes());
    }
};

TohpeProfiler g_tohpe_profiler;

struct ScopedNsTimer {
    explicit ScopedNsTimer(uint64_t& target) : _target(target), _begin(Clock::now()) {}
    ~ScopedNsTimer() {
        _target += static_cast<uint64_t>(std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - _begin).count());
    }
    uint64_t& _target;
    Clock::time_point _begin;
};

/**
 * @brief Transform matrix to reduced-row-echelon form, store row operations in augmented_matrix, and get it pivots.
 *
 * @param matrix, augmented_matrix, pivots
 * @return null
 */
dvlab::BooleanMatrix::Row kernel(dvlab::BooleanMatrix& matrix, dvlab::BooleanMatrix& augmented_matrix, std::unordered_map<size_t, size_t>& pivots) {
    spdlog::trace("kernel start");
    g_tohpe_profiler.kernel_rows_processed += matrix.num_rows();
    g_tohpe_profiler.kernel_entry_pivots += pivots.size();
    for (auto i : std::views::iota(0ul, matrix.num_rows())) {
        if (pivots.contains(i)) continue;
        for (auto const& [key, value] : pivots) {
            if (matrix[i][value]) {
                ++g_tohpe_profiler.kernel_row_adds;
                matrix[i] += matrix[key];
                augmented_matrix[i] += augmented_matrix[key];
            }
        }

        if (matrix[i].is_zeros()) {
            spdlog::debug("matrix[{}] is zero, return:", i);
            augmented_matrix[i].print_row();
            return augmented_matrix[i];
        }
        size_t const first_one_idx = std::distance(matrix[i].begin(), std::find(matrix[i].begin(), matrix[i].end(), true));

        auto pivot           = matrix[i];
        auto augmented_pivot = augmented_matrix[i];
        for (auto const& [key, _] : pivots) {
            if (matrix[key][first_one_idx]) {
                ++g_tohpe_profiler.kernel_row_adds;
                matrix[key] += pivot;
                augmented_matrix[key] += augmented_pivot;
            }
        }
        pivots[i] = first_one_idx;
    }
    spdlog::trace("kernel end");
    return dvlab::BooleanMatrix::Row(0);
}

/**
 * @brief Transform matrix to reduced-row-echelon form, store row operations in augmented_matrix, and get it pivots.
 *
 * @param matrix, augmented_matrix, pivots
 * @return null
 */
dvlab::BooleanMatrix get_l_matrix(dvlab::BooleanMatrix const& phase_poly_matrix, dvlab::BooleanMatrix const& row_products) {
    auto l_matrix       = phase_poly_matrix;  // copy
    auto const num_rows = phase_poly_matrix.num_rows();

    // [Note] the index here is different from the index of original repository impl
    auto const get_row_product_idx = [&num_rows](size_t a, size_t b) {
        if (a > b) {
            std::swap(a, b);
        }
        return (a * num_rows) - (a * (a + 1) / 2) + b - a - 1;
    };

    auto const id_vec = std::views::iota(0ul, num_rows) | tl::to<std::vector>();

    for (auto const& [a, b] : dvlab::combinations<2>(id_vec)) {
        auto const new_row = row_products[get_row_product_idx(a, b)];
        l_matrix.push_row(new_row);
    }

    return l_matrix;
}

[[maybe_unused]] int calculate_score(dvlab::BooleanMatrix::Row const& y, std::vector<std::pair<int, int>> const& s_matrix) {
    const int abs_y = static_cast<int>(y.sum() % 2);
    int ret         = -1 * abs_y;

    for (auto& indexes : s_matrix) {
        if (indexes.first == indexes.second) {
            ret += y[indexes.first] + 2 * (y[indexes.first] == 0) * (abs_y);
        } else {
            ret += 2 * (y[indexes.first] ^ y[indexes.second]);
        }
    }
    return ret;
}

bool compare_row_value(dvlab::BooleanMatrix::Row const& a, dvlab::BooleanMatrix::Row const& b) {
    DVLAB_ASSERT(a.size() == b.size(), "Row compare should have same size.");
    for (auto i : std::views::iota(0ul, a.size())) {
        if (!a[i] ^ b[i]) continue;
        return a[i] < b[i];
    }
    // a == b
    return false;
}

void clear_column(size_t idx, dvlab::BooleanMatrix& matrix, dvlab::BooleanMatrix& augmented_matrix, std::unordered_map<size_t, size_t>& pivots) {
    ++g_tohpe_profiler.clear_column_calls;
    if (!pivots.contains(idx)) return;
    auto val = pivots[idx];
    pivots.erase(idx);

    if (!augmented_matrix[idx][idx]) {
        for (auto j : std::views::iota(0ul, matrix.num_rows())) {
            if (!augmented_matrix[j][idx]) {
                continue;
            }
            pivots[j]             = val;
            auto col              = matrix[j];
            auto augmented_col    = augmented_matrix[j];
            matrix[j]             = matrix[idx];
            augmented_matrix[j]   = augmented_matrix[idx];
            matrix[idx]           = col;
            augmented_matrix[idx] = augmented_col;
            break;
        }
    }

    auto col           = matrix[idx];
    auto augmented_col = augmented_matrix[idx];
    for (auto j : std::views::iota(0ul, matrix.num_rows())) {
        if (augmented_matrix[j][idx] && idx != j) {
            ++g_tohpe_profiler.clear_column_row_adds;
            matrix[j] += col;
            augmented_matrix[j] += augmented_col;
        }
    }
}

/** Build one row of transpose(L) from a term, matching C++ get_row_products/get_l_matrix column order. */
dvlab::BooleanMatrix::Row build_l_transpose_row_from_term(dvlab::BooleanMatrix::Row const& term_z, size_t n_qubits) {
    DVLAB_ASSERT(term_z.size() >= n_qubits, "term row shorter than n_qubits");
    std::vector<unsigned char> vec;
    vec.reserve(n_qubits + (n_qubits * (n_qubits - 1) / 2));
    for (size_t q = 0; q < n_qubits; ++q) {
        vec.emplace_back(term_z[q]);
    }
    for (size_t a = 0; a < n_qubits; ++a) {
        for (size_t b = a + 1; b < n_qubits; ++b) {
            vec.emplace_back(static_cast<unsigned char>(term_z[a] & term_z[b]));
        }
    }
    return dvlab::BooleanMatrix::Row(vec);
}

dvlab::BooleanMatrix get_z_matrix(dvlab::BooleanMatrix const& phase_poly_matrix) {
    auto z_matrix_transposed = dvlab::transpose(phase_poly_matrix);
    auto const num_terms     = phase_poly_matrix.num_cols();
    auto const n_qubits      = phase_poly_matrix.num_rows();

    // Enumerate pairs of existing phase-polynomial terms. Using the pair count here
    // produces synthetic indices past the last real column and can segfault.
    auto const id_vec = std::views::iota(0ul, num_terms) | tl::to<std::vector>();
    auto seen_z       = std::unordered_set<dvlab::BooleanMatrix::Row, dvlab::BooleanMatrixRowHash>();

    if (num_terms < 2) {
        return z_matrix_transposed;
    }

    for (auto const& [a, b] : dvlab::combinations<2>(id_vec)) {
        if (stop_requested()) {
            return phase_poly_matrix;
        }
        dvlab::BooleanMatrix::Row z(n_qubits);
        for (size_t k = 0; k < n_qubits; ++k) {
            z[k] = phase_poly_matrix[k][a] ^ phase_poly_matrix[k][b];
        }

        if (seen_z.contains(z)) {
            continue;
        }
        seen_z.insert(z);

        z_matrix_transposed.push_row(z);
    }
    return z_matrix_transposed;
}

std::vector<std::vector<std::pair<int, int>>> get_s_matrices(dvlab::BooleanMatrix const& phase_poly_matrix, dvlab::BooleanMatrix const& z_matrix) {
    std::vector<std::vector<std::pair<int, int>>> s;
    auto phase_poly_matrix_transposed = dvlab::transpose(phase_poly_matrix);
    auto const n_qubits               = phase_poly_matrix.num_rows();
    auto const num_terms              = phase_poly_matrix.num_cols();
    auto const id_vec                 = std::views::iota(0ul, num_terms) | tl::to<std::vector>();

    for (auto& z : z_matrix.get_matrix()) {
        std::vector<std::pair<int, int>> s_z;
        // a != b
        for (auto const& [a, b] : dvlab::combinations<2>(id_vec)) {
            dvlab::BooleanMatrix::Row tmp(n_qubits);
            for (size_t k = 0; k < n_qubits; ++k) {
                tmp[k] = phase_poly_matrix[k][a] ^ phase_poly_matrix[k][b];
            }
            if (tmp == z) {
                s_z.push_back(std::make_pair(a, b));
            }
        }
        // a == b
        for (auto a : std::views::iota(0ul, num_terms)) {
            if (phase_poly_matrix_transposed.get_row(a) == z) {
                s_z.push_back(std::make_pair(a, a));
            }
        }
        s.push_back(s_z);
    }

    return s;
}

Polynomial tohpe_once(Polynomial const& polynomial) {
    if (polynomial.empty()) {
        return polynomial;
    }

    // auto const n_qubits = polynomial.front().n_qubits();

    // Each column represents a term in the phase polynomial, and each row represents a qubit.
    auto const phase_poly_matrix = load_phase_poly_matrix(polynomial);

    auto const idx_vec = std::views::iota(0ul, polynomial.size()) | tl::to<std::vector>();

    auto const row_products = get_row_products(phase_poly_matrix);

    auto const l_matrix             = get_l_matrix(phase_poly_matrix, row_products);
    auto const z_matrix             = get_z_matrix(phase_poly_matrix);
    auto const s_matrices           = get_s_matrices(phase_poly_matrix, z_matrix);
    auto const nullspace_transposed = get_nullspace_transposed(l_matrix);

    tohpe_trace("slow", "input", fmt::format("terms={}", summarize_term_bitstrings(polynomial_term_bitstrings(polynomial))));

    // check if y satisfies tohpe condition
    if (nullspace_transposed.is_empty()) {
        tohpe_trace("slow", "no-nullspace", "nullspace is empty");
        return polynomial;
    }

    auto phase_poly_matrix_copy = phase_poly_matrix;
    for (auto const& y : nullspace_transposed) {
        if (y.is_zeros()) {
            continue;
        } else if (y.sum() != y.size() && y.sum() % 2 == 0) {
            // y is candidate

            int max_score    = std::numeric_limits<int>::min();
            size_t max_index = 0;
            for (auto a : std::views::iota(0ul, s_matrices.size())) {
                auto const& s_matrix = s_matrices[a];
                auto score           = calculate_score(y, s_matrix);
                if (score > max_score) {
                    max_score = score;
                    max_index = a;
                }
            }

            auto& chosen_z = z_matrix[max_index];

            spdlog::debug("Found a TOHPE move");
            spdlog::debug("- z: {}", fmt::join(chosen_z, ""));
            spdlog::debug("- y: {}", fmt::join(y, ""));

            for (auto const i : std::views::iota(0ul, phase_poly_matrix_copy.num_rows())) {
                if (chosen_z[i] == 1) {
                    phase_poly_matrix_copy[i] += y;
                }
            }

            phase_poly_matrix_copy = dvlab::transpose(phase_poly_matrix_copy);
            if (y.sum() % 2 == 1) {
                phase_poly_matrix_copy.push_row(chosen_z);
            }
            auto const output = from_boolean_matrix(phase_poly_matrix_copy);
            tohpe_trace("slow", "move", fmt::format("y={} z={} output={}",
                                                     row_to_bitstring(y), row_to_bitstring(chosen_z),
                                                     summarize_term_bitstrings(polynomial_term_bitstrings(output))));
            return output;
        }
    }
    // no candidate, return same matrix
    auto const output = from_boolean_matrix(dvlab::transpose(phase_poly_matrix_copy));
    tohpe_trace("slow", "no-candidate", fmt::format("output={}", summarize_term_bitstrings(polynomial_term_bitstrings(output))));
    return output;
}

Polynomial tohpe_once_fast_port(Polynomial const& polynomial, size_t max_moves = static_cast<size_t>(-1)) {
    if (polynomial.empty()) {
        return polynomial;
    }
    if (g_tohpe_profiler.enabled) {
        ++g_tohpe_profiler.tohpe_once_calls;
    }
    auto const n_qubits = polynomial.front().n_qubits();

    auto phase_poly_matrix = load_phase_poly_matrix(polynomial);
    auto table             = dvlab::transpose(phase_poly_matrix);

    auto row_products = get_row_products(phase_poly_matrix);

    auto l_matrix           = get_l_matrix(phase_poly_matrix, row_products);
    auto l_matrix_transpose = dvlab::transpose(l_matrix);

    dvlab::BooleanMatrix augmented_matrix = dvlab::identity(table.num_rows());
    std::unordered_map<size_t, size_t> pivots;

    tohpe_trace("fast_port", "input", fmt::format("terms={}", summarize_term_bitstrings(polynomial_term_bitstrings(polynomial))));

    size_t move_idx = 0;
    while (true) {
        spdlog::trace("loop {}", move_idx);
        auto n_terms                  = table.num_rows();
        dvlab::BooleanMatrix::Row y   = kernel(l_matrix_transpose, augmented_matrix, pivots);
        auto const pivots_before_move = pivots.size();
        if (y.is_zeros()) {
            spdlog::trace("No y is found, end tohpe.");
            tohpe_trace("fast_port", "no-y", fmt::format("move={} pivots={}", move_idx, pivots_before_move));
            break;
        }

        std::unordered_map<dvlab::BooleanMatrix::Row, int, dvlab::BooleanMatrixRowHash> map;
        auto const y_parity = y.sum() % 2;
        for (auto i : std::views::iota(0ul, n_terms)) {
            if (y_parity && !y[i]) {
                map[table[i]] = 1;
                ++g_tohpe_profiler.score_seed_inserts;
            } else if (!y_parity && y[i]) {
                map[table[i]] = 1;
                ++g_tohpe_profiler.score_seed_inserts;
            }
        }

        for (auto i : std::views::iota(0ul, n_terms)) {
            if (!y[i]) continue;
            for (auto j : std::views::iota(0ul, n_terms)) {
                if (y[j]) continue;
                ++g_tohpe_profiler.score_pair_count;
                auto const z = table[i] + table[j];
                if (!map.contains(z)) {
                    map[z] = 0;
                }
                map[z] += 2;
                ++g_tohpe_profiler.score_map_updates;
            }
        }

        auto max_score  = int{0};
        auto max_z      = dvlab::BooleanMatrix::Row(n_qubits);
        bool have_max_z = false;

        for (auto const& [z, score] : map) {
            if (score > max_score || (have_max_z && score == max_score && compare_row_value(z, max_z))) {
                max_score  = score;
                max_z      = z;
                have_max_z = true;
            }
        }

        if (max_score <= 0) {
            tohpe_trace("fast_port", "nonpositive-score", fmt::format("move={} y={} best_score={}", move_idx, row_to_bitstring(y), max_score));
            break;
        }

        std::vector<unsigned char> to_update(y.begin(), y.begin() + static_cast<long>(n_terms));
        if (y_parity) {
            table.push_zeros_row();
            l_matrix_transpose.push_zeros_row();
            augmented_matrix.push_zeros_column();
            augmented_matrix.push_zeros_row();
            augmented_matrix[augmented_matrix.num_rows() - 1][augmented_matrix.num_cols() - 1] = 1;
            to_update.push_back(1);
        }
        for (auto i : std::views::iota(0ul, to_update.size())) {
            if (to_update[i] == 0) continue;
            table[i] += max_z;
        }

        std::unordered_map<dvlab::BooleanMatrix::Row, size_t, dvlab::BooleanMatrixRowHash> hashmap;
        std::vector<size_t> to_remove;
        n_terms = table.num_rows();
        for (auto i : std::views::iota(0ul, n_terms)) {
            if (table[i].is_zeros()) {
                to_remove.push_back(i);
                continue;
            } else if (hashmap.contains(table[i])) {
                to_remove.push_back(hashmap[table[i]]);
                to_remove.push_back(i);
                hashmap.erase(table[i]);
            } else {
                hashmap[table[i]] = i;
            }
        }
        std::sort(to_remove.begin(), to_remove.end(), std::greater<size_t>());
        g_tohpe_profiler.removed_rows += to_remove.size();

        for (auto const& i : to_remove) {
            clear_column(i, l_matrix_transpose, augmented_matrix, pivots);
            table[i] = table[table.num_rows() - 1];
            table.erase_row(table.num_rows() - 1);
            l_matrix_transpose[i] = l_matrix_transpose[l_matrix_transpose.num_rows() - 1];
            l_matrix_transpose.erase_row(l_matrix_transpose.num_rows() - 1);
            augmented_matrix[i] = augmented_matrix[augmented_matrix.num_rows() - 1];
            augmented_matrix.erase_row(augmented_matrix.num_rows() - 1);
            to_update[i] = to_update.back();
            to_update.pop_back();
            if (pivots.contains(table.num_rows())) {
                pivots[i] = pivots[table.num_rows()];
                pivots.erase(table.num_rows());
            }
            auto const row_count = table.num_rows();
            for (auto j : std::views::iota(0ul, augmented_matrix.num_rows())) {
                auto& row = augmented_matrix[j];
                if (row.size() <= row_count) continue;

                if (i < row.size()) {
                    if (row[i] != row[row_count]) {
                        row[i] ^= 1;
                    }
                }
                if (row[row_count]) {
                    row[row_count] ^= 1;
                }
            }
        }

        auto const aug_width = table.num_rows();
        for (auto i : std::views::iota(0ul, aug_width)) {
            if (augmented_matrix[i].size() > aug_width) {
                std::vector<unsigned char> trimmed(augmented_matrix[i].begin(), augmented_matrix[i].begin() + static_cast<ptrdiff_t>(aug_width));
                augmented_matrix[i].set_row(std::move(trimmed));
            }
        }

        std::vector<size_t> update_indices;
        update_indices.reserve(to_update.size());
        for (auto i : std::views::iota(0ul, to_update.size())) {
            if (to_update[i]) {
                update_indices.push_back(i);
            }
        }
        for (auto const i : update_indices) {
            clear_column(i, l_matrix_transpose, augmented_matrix, pivots);
            {
                ScopedNsTimer sub(g_tohpe_profiler.ns_rebuild_get_l_transpose);
                l_matrix_transpose[i] = build_l_transpose_row_from_term(table[i], n_qubits);
            }
            std::vector<unsigned char> bv(table.num_rows(), 0);
            bv[i] = 1;
            augmented_matrix[i].set_row(bv);
        }

        auto const output = from_boolean_matrix(table);
        tohpe_trace("fast_port", "move", fmt::format("move={} y={} z={} to_update={} to_remove={} pivots_before={} pivots_after={} output={}",
                                                      move_idx, row_to_bitstring(y), row_to_bitstring(max_z),
                                                      bytes_to_bitstring(to_update), summarize_indices(to_remove),
                                                      pivots_before_move, pivots.size(),
                                                      summarize_term_bitstrings(polynomial_term_bitstrings(output))));
        ++move_idx;
        if (move_idx >= max_moves) {
            return output;
        }
    }
    auto const output = from_boolean_matrix(table);
    tohpe_trace("fast_port", "output", fmt::format("terms={}", summarize_term_bitstrings(polynomial_term_bitstrings(output))));
    return output;
}
}  // namespace

std::pair<StabilizerTableau, Polynomial> TohpePhasePolynomialOptimizationStrategy::optimize(StabilizerTableau const& clifford, Polynomial const& polynomial) const {
    if (polynomial.empty()) {
        fmt::println("Polynomial is empty, returning the input Clifford and polynomial");
        return {clifford, polynomial};
    }

    if (std::ranges::any_of(polynomial, [](PauliRotation const& rotation) { return 4 % rotation.phase().denominator() != 0; })) {
        spdlog::error("Failed to perform TODD optimization: the polynomial contains a non-4th-root-of-unity phase!!");
        return {clifford, polynomial};
    }

    auto ret_clifford   = clifford;
    auto ret_polynomial = polynomial;

    properize(ret_clifford, ret_polynomial);

    auto multi_linear_polynomial = MultiLinearPolynomial();
    multi_linear_polynomial.add_rotations(ret_polynomial, false);

    spdlog::debug("num_terms before TOHPE: {}", ret_polynomial.size());
    spdlog::trace("Polynomial before TOHPE:\n{}", fmt::join(ret_polynomial, "\n"));

    size_t move_idx = 0;
    while (true) {
        auto const num_terms = ret_polynomial.size();
        auto const input_summary = summarize_term_bitstrings(polynomial_term_bitstrings(ret_polynomial));
        std::optional<Polynomial> fast_port_polynomial;
        if (tohpe_trace_compare_enabled()) {
            fast_port_polynomial = tohpe_once_fast_port(ret_polynomial, 1);
        }
        auto const slow_polynomial = tohpe_use_fast_port() ? tohpe_once_fast_port(ret_polynomial) : tohpe_once(ret_polynomial);
        if (tohpe_trace_compare_enabled()) {
            auto const slow_summary = summarize_term_bitstrings(polynomial_term_bitstrings(slow_polynomial));
            auto const fast_summary = summarize_term_bitstrings(polynomial_term_bitstrings(*fast_port_polynomial));
            auto const same_output  = slow_summary == fast_summary;
            tohpe_trace("compare", same_output ? "agree" : "diverge",
                        fmt::format("move={} input={} slow={} fast={}", move_idx, input_summary, slow_summary, fast_summary));
            if (!same_output && tohpe_trace_stop_on_divergence()) {
                ret_polynomial = slow_polynomial;
                break;
            }
        }
        ret_polynomial = slow_polynomial;
        ++move_idx;
        if (ret_polynomial.empty() || ret_polynomial.size() == num_terms) {
            break;
        }
        spdlog::debug("num_terms after TOHPE: {}", ret_polynomial.size());
        spdlog::trace("Polynomial after TOHPE:\n{}", fmt::join(ret_polynomial, "\n"));
    }

    auto optimized_multi_linear_polynomial = MultiLinearPolynomial();
    optimized_multi_linear_polynomial.add_rotations(ret_polynomial, false);
    if (!multi_linear_polynomial.has_same_signature_tensor(optimized_multi_linear_polynomial)) {
        spdlog::error("Failed to perform TOHPE optimization: the post-optimization polynomial does not have the same signature tensor as the pre-optimization polynomial!!");
        return {clifford, polynomial};
    }

    multi_linear_polynomial.add_rotations(ret_polynomial, true);

    if (auto clifford_ops = multi_linear_polynomial.extract_clifford_operators(); clifford_ops.has_value()) {
        ret_clifford.apply(*clifford_ops);
    } else {
        spdlog::error("Failed to perform TOHPE optimization: the post-optimization polynomial does not have the same signature as the pre-optimization polynomial!!");
        return {clifford, polynomial};
    }

    return {ret_clifford, ret_polynomial};
}

std::pair<StabilizerTableau, Polynomial> FastToddPhasePolynomialOptimizationStrategy::optimize(StabilizerTableau const& clifford, Polynomial const& polynomial) const {
    return {clifford, polynomial};
}

}  // namespace qsyn::tableau
